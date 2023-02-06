# 1 "gem_main.f90"
!       3D Flux Tube Toroidal Electromagnetic GK Delta-f Code
!   global variables...
program gem_main

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none
  integer :: n,i,j,k,ip,ns
  real :: dum1,dum2,dum3,dum4
  
  call initialize
  call mpi_barrier(mpi_comm_world,ierr)

  starttm=MPI_WTIME()  
  do  timestep=ncurr,nm
     tcurr = tcurr+dt
     
     poisson_start_tm=poisson_start_tm+MPI_WTIME()
     call poisson(timestep-1,0)
     poisson_end_tm=poisson_end_tm+MPI_WTIME()

     call reporter(timestep-1)

     if(mod(timestep,1000)==0)then
        do i=0,last 
           !                 if(myid==i)write(*,*)myid,mm(1),mme
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do
     end if

  end do
  lasttm=MPI_WTIME()
  
  call mpi_barrier(mpi_comm_world,ierr)

  tottm=lasttm-starttm

  call MPI_REDUCE(tottm,total_tm,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)  
  
  total_tm=total_tm/numprocs
  
  dum1 = poisson_end_tm-poisson_start_tm
  call MPI_REDUCE(dum1,poisson_tm,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)  
  poisson_tm = poisson_tm/numprocs

   if(myid==0)write(*,*)'right caculation'
   if(myid==0)write(*,*)'total_tm',total_tm, 'poisson_tm',poisson_tm

  !  write(*,*)'ps time=',pstm,'tot time=',tottm
  do i=0,last 
     !            if(myid==i)write(*,*)myid,mm(1),mme
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end do

100 call MPI_FINALIZE(ierr)

end program gem_main

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine hybinit

  use gem_com
  use gem_equil
  implicit none
  GCLR = int(MyId/ntube)
  GLST = numprocs/ntube-1
  TCLR = mod(MyId,ntube)
  TLST = ntube-1

  !***MPI variables

  !      if (GCLR.eq.GLST) 
  !     %     mykm=km-GLST*mykm
  mykm = 1
  rngbr = GCLR+1
  if (GCLR.eq.GLST)rngbr = 0

  lngbr = GCLR-1
  if (GCLR.eq.0) lngbr = GLST

  idnxt = TCLR+1
  if (TCLR.eq.TLST)idnxt = 0

  idprv = TCLR-1
  if (TCLR.eq.0) idprv = TLST
  !     
  !      return
end subroutine hybinit

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine init

  use gem_com
  use gem_equil
  use omp_lib
  implicit none

!include 'fftw3.f03'
  character(len=62) dumchar
  INTEGER :: i,j,k,m,n,ns,idum,i1,j1,k1,nthreads,iam,i_err
  INTEGER :: mm1,mm2,lr1
  real :: mims1,tets1,q1,kappan,kappat,r,qr,th,cost,dum,zdum
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp
  real :: grp,gxdgyp,jacp,jfnp,gn0ep,gt0ep,gt0ip,grdgtp,gthp,gsfp,gthfp
  real,DIMENSION (1:5) :: gn0sp,gt0sp
  real :: wx0,wx1,wz0,wz1,b

!$omp parallel private(iam)
nthreads = omp_get_num_threads()
iam = omp_get_thread_num()
!write(*,*)'myid, nthreads, iam = ',myid, nthreads, iam
!$omp end parallel
n_omp = nthreads


  !jycheng
  namelist /primary_parameters/ itube,mimp,mcmp,chgi,chgc,imx,jmx,kmx,nsmx,ntube,lxa, &
                                lymult,jcnt,dt,nm,nsm,amp,r0a,ifluid,ipbm,ias,isg,amie,rneui, &
                                betai,nonlin1,nonlin2,nonline,iflut,iexb,ipara,vcut,wecut,micell,mecell, &
                                qbeam,mbeam,isonew
  namelist /control_parameters/ iperi,iperidf,delra,delri,delre,delrn,nlow,xshape, &
                                yshape,zshape,iput,iget,igetmx,ision,isiap,peritr,llk,mlk,onemd,izonal, &
                                adiabatic_electron,ineq0,nzcrt,npze,npzi,npzc,npzb,iphbf,iapbf, &
                                idpbf,cut,kxcut,kycut,bcut,c4,vexbsw,vparsw,mach,gamma_E,isuni,isunie, &
                                lr1,iflr,iorb,nxsrc,nesrc,nzsrc,gammah,ghzon,gamgtc,gamtoy,gamgyro,nzgyro
  namelist /diagnostics_parameters/ icrs_sec,ipg,isphi,nplot,xnplt,isft,mynf,frmax,ifskp,idg
  namelist /fluxtube/ lxmult,Rovera,elon0,selon0,tria0,stria0,rmaj0p,q0,shat0,teti,tcti,rhoia,Rovlni,Rovlti, &
                       Rovlne,Rovlte,Rovlnc,Rovltc,ncne,nuacs
  namelist /others/ nrst,eprs,tor,vpp,vt0,yd0,nmx

  IU=cmplx(0.,1.)
  pi=4.0*atan(1.0)
  pi2 = pi*2.

  call ppinit_mpi(myid,numprocs)
  last=numprocs-1
  !the initial timestep index
  timestep=0
  !the initial timestep 
  tcurr=0.
  open(unit=115,file='gem.in',status='old',action='read')
  read(115,nml=primary_parameters)
  read(115,nml=control_parameters)
  read(115,nml=diagnostics_parameters)
  read(115,nml=fluxtube)
  read(115,nml=others)
  close(115)
  if(idg==1)then
     do i = 0,last
        if(myid==i)then
           open(9, file='plot',status='unknown',position='append')
           write(9,*)'pass read namelists',myid
           close(9)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     end do
  end if


  nonlin(1)=nonlin1
  nonlin(2)=nonlin2


# 168
!$acc enter data copyin(nonlin)





# 174
  ntube=int(numprocs/kmx)
  if((mod(numprocs,kmx) .ne. 0) .and. myid==0)then
    write(*,*)'WARNING: MPI number is a multiple of kmx, ntube is not an integer'
    stop
  endif

  !mm1=int(imx,8)*int(jmx,8)*int(kmx,8)*int(micell,8)
  mm1=imx*jmx*kmx*micell
  !mm2=int(imx,8)*int(jmx,8)*int(kmx,8)*int(mecell,8)
  mm2=imx*jmx*kmx*mecell
  mmx=int(real(mm1/int(kmx*ntube))*1.5)
  mmxe=int(real(mm2/int(kmx*ntube))*1.5)
  call ppinit_decomp(myid,numprocs,ntube,tube_comm,grid_comm)
  call hybinit
  call mpi_barrier(mpi_comm_world,ierr)
  if(idg==1)then
     do i = 0,last
        if(myid==i)then
           open(9, file='plot',status='unknown',position='append')
           write(9,*)'pass hybinit',myid
           close(9)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     end do
  end if
  
  im=imx;jm=jmx;km=kmx
  cgpfacx = 4
  cgpfacy = 2
  icgp = imx/cgpfacx
  jcgp = jmx/cgpfacy

  if(isft==1) iget = 1
  nfreq = 400
  call new_gem_com()
  if(idg==1)then
     do i = 0,last
        if(myid==i)then
           open(9, file='plot',status='unknown',position='append')
           write(9,*)'pass new_equil',myid
           close(9)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     end do
  end if
  !      rin = r0-lx/2
  !      rout = r0+lx/2
  rina = r0a-lxa/2.
  routa = r0a+lxa/2.

  call new_equil()
  lx=lxa*a
  ly=2.*pi*r0/q0abs/lymult
  !beta has been used all the time, so the definition is fixed, declared in gem_com.f. betai is declared in equil.f
  !in gem.in betai is intended to be GYRO's beta_e (within a factor of 2). It is defined to be (mu0 ne Te)/Bunit**2
  beta=betai 
  br0 = rmaj0
  lr0 = r0
  qp = q0p
  lz = pi2*q0abs*rmaj0
  delz = lz/ntheta
  rneu = nuacs
  kycut=pi2/ly*(jcnt-1)/2*1.1
  if(itube==1)then
     kxcut=pi2*shat0*kycut*4
  end if

  if(myid.eq.master)then
     open(19, file='eqdat', status='unknown')
     write(19,*) 'r0a,lxa,beta,a,rmaj0=', r0a,lxa,beta,a,rmaj0
     write(19,*)'elon0,tria0,rmaj0p,stria0,selon0,q0,shat=', elon0,tria0,rmaj0p,stria0,selon0,q0,r0*q0p/q0
     write(19,*)'xn0s(1,nr/2),xn0s(2,nr/2),t0s(1,nr/2),t0s(2,nr/2)=',xn0s(1,nr/2),xn0s(2,nr/2),t0s(1,nr/2),t0s(2,nr/2)
     write(19,*)'capns(1,nr/2),capns(2,nr/2),capts(1,nr/2),capts(2,nr/2)=', capns(1,nr/2),capns(2,nr/2),capts(1,nr/2),capts(2,nr/2)
     write(19,*)'capne(nr/2),capte(nr/2)=', capne(nr/2),capte(nr/2)
     write(19,*)'t0i0,t0e0,ni0,ne0=', t0i(nr/2),t0e(nr/2),xn0i(nr/2),xn0e(nr/2)
     write(19,*)'a/cs=', a/sqrt(t0e(nr/2)/mimp)
     write(19,*)'i,   sf,   ni,   ne,   nc,   ti,   tc,   capni,   capnc, captc'
     do i = 0,nr
        write(19,99)i,sf(i),xn0e(i),xn0i(i),xn0c(i),t0e(i),t0i(i),t0c(i),capne(i),capni(i),capnc(i),capte(i),capti(i),captc(i)
     end do
     write(19,99)i,cn0e,cn0s(1),cn0s(2),cn0s(3)
     close(19)
  end if
98 format(7(1x,e14.7))
99 format(1x,i3,2x,15(1x,e10.3))
  ! equilibrium data for transport coefficients calculation by mablab with flux data
  if(myid.eq.master)then
     open(912, file='eqflux', status='unknown')
     do j = 1, nsubd
        i = int(nr*(2*j-1)/(2*nsubd))
        write(912,9991) xn0e(i)*cn0e,capne(i),xn0c(i)*cn0c,capnc(i),t0i(i),capti(i)
     enddo
     close(912)
  endif
9991 format(6(2x,e12.4)) 
  !
  !      write(*,*)'br0,lr0,q0,qp = ', br0,lr0,q0,qp

  dte = dt
  iadi = 0
  if(isg.gt.0.)fradi = isg
  if(ifluid.eq.0)then
     iadi = 1
     fradi = 1.0
  end if

  !     begin reading species info, ns=1,nsm...
  ns = 1
  mims(3)=mbeam
  q(3)=qbeam
  mims(1)=mimp;mims(2)=mcmp;q(1)=chgi;q(2)=chgc


# 287
  !$acc update device(mims,q)





# 293
  rhoia = mims(1)*sqrt(T0e(nr/2)/mims(1)) /(q(1)*1.) /a
  tmm(1)=mm1
  tmm(2)=mm2
  mm(:)=int(mm1/numprocs)
  mme = int(mm2/numprocs)
  if (MyId.eq.Last) mm(ns)=mm1-Last*mm(ns)
  !     write(*,*)'in init  ',Myid,mm(ns)
  tets(1)=1
  lr(1)=lr1
  lr(2)=lr1
  lr(3)=lr1


# 306
  !$acc update device(lr)




  
# 312
  pzcrite = abs(psi(nr)-psi(0))/br0/npze
  encrit = 0.2*t0e(nr/2)
  pzcrit(1) = q(1)*abs(psi(nr)-psi(0))/br0/npzi
  pzcrit(2) = q(2)*abs(psi(nr)-psi(0))/br0/npzc
  pzcrit(3) = q(3)*abs(psi(nr)-psi(0))/br0/npzb

  kapn(ns)=kappan
  kapt(ns)=kappat

  emass = 1./amie
  qel = -1.

  vsphere = sqrt(vcut*1.1) * sqrt(tge/emass)
  cv = 4.0/3.0*pi * vsphere**3

  ecutsrc = vcut
  
  if(iget.eq.1) amp=0.
  !     totvol is the square for now...
  dx=lx/float(im)
  dy=ly/float(jm)
  dz=lz/float(km)
  !      totvol=lx*ly*lz
  dxcgp = lx/icgp
  dycgp = ly/jcgp
  dxsrc = lx/nxsrc
  desrc = ecutsrc/nesrc
  
  e0=lr0/q0/br0
  !     
  do i=0,nxpp
     xg(i)=i*dx !dx*(tclr*nxpp+i)
  enddo
  do j=0,jm
     yg(j)=dy*float(j)
  enddo
  kcnt=1

  do k=0,mykm
     n=GCLR*kcnt+k 
     zg(k)=dz*float(n)
  enddo

  !      jcnt = 3  !jmx/ntube                                                                                                                                                                         
  mstart = 0
  ntor0 = mstart+1
  do m = 0,jcnt-1
     isgnft(m) = 1
     j1 = mstart+int((float(m)+1.0)/2)
     jft(m) = j1
     if(m==0)then
        isgnft(m) = 1
        jft(m) = 0
     end if
     if(m>0.and.mod(m,2)==0)then
        isgnft(m) = -1
        jft(m) = jmx-j1
     end if
  end do

  !     initialize bfld   

  zfnth(0) = 0.
  do j = 1,ntheta
     zfnth(j) = zfnth(j-1)+dth*q0*br0*(1./jfn(j-1)+1./jfn(j))/2
  end do
  if(q0<0)zfnth = zfnth+lz

  thfnz(0) = -pi
  thfnz(ntheta) = pi
  if(q0<0.)then
     thfnz(0) = pi
     thfnz(ntheta) = -pi
  end if

  if(q0>0)then
     k = 0
     do j = 1,ntheta-1
        zdum = j*lz/ntheta
        do i = k,ntheta-1
           if(zfnth(i)<=zdum.and.zfnth(i+1)>zdum)then
              k = i
              dum = (zdum-zfnth(i))*dth/(zfnth(i+1)-zfnth(i))
              thfnz(j) = i*dth-pi+dum
              go to 127
           end if
        end do
127     continue
     end do
  end if

  if(q0<0)then
     k = 0
     do j = 1,ntheta-1
        zdum = lz-j*lz/ntheta
        do i = k,ntheta-1
           if(zfnth(i)>=zdum.and.zfnth(i+1)<zdum)then
              k = i
              dum = (zdum-zfnth(i))*dth/(zfnth(i+1)-zfnth(i))
              thfnz(ntheta-j) = i*dth-pi+dum
              go to 128
           end if
        end do
128     continue
     end do
  end if

  do i1 = 0,nxpp
     r = xg(i1)-0.5*lx+lr0
     do k1 = 0,mykm
        k = int(zg(k1)/delz)
        k = min(k,ntheta-1)
        wz0 = ((k+1)*delz-zg(k1))/delz
        wz1 = 1-wz0
        th = wz0*thfnz(k)+wz1*thfnz(k+1)
        i = int((r-rin)/dr)
        i = min(i,nr-1)
        wx0 = (rin+(i+1)*dr-r)/dr
        wx1 = 1.-wx0
        k = int((th+pi)/dth)
        k = min(k,ntheta-1)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1.-wz0
        dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
             +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
        dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
             +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
        grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
             +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 

        grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
             +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 
        gthp = wx0*wz0*gth(i,k)+wx0*wz1*gth(i,k+1) &
             +wx1*wz0*gth(i+1,k)+wx1*wz1*gth(i+1,k+1) 

        bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
        radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
             +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
        dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
             +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
        qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
             +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
        grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
             +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
        jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
             +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1) 
        gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
             +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
        gthfp = wx0*wz0*thflx(i,k)+wx0*wz1*thflx(i,k+1) &
             +wx1*wz0*thflx(i+1,k)+wx1*wz1*thflx(i+1,k+1) 

        fp = wx0*f(i)+wx1*f(i+1)
        gsfp = wx0*sf(i)+wx1*sf(i+1)                
        jfnp = wz0*jfn(k)+wz1*jfn(k+1)
        psipp = wx0*psip(i)+wx1*psip(i+1)        
        gn0ep = wx0*xn0e(i)+wx1*xn0e(i+1)        
        gt0ep = wx0*t0e(i)+wx1*t0e(i+1)        
        do ns = 1, nsm
           gn0sp(ns) = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)
           gt0sp(ns) = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)           
           gn0s(ns,i1) = gn0sp(ns)
           gt0s(ns,i1) = gt0sp(ns)
        enddo
        gt0ip = wx0*t0s(1,i)+wx1*t0s(1,i+1)        
        b=1.-tor+tor*bfldp
        cfx(i1,k1) = br0/b**3*fp/radiusp*dbdtp*grcgtp
        cfy(i1,k1) = br0/b**3*fp/radiusp* &
             (dydrp*dbdtp-lr0/q0*qhatp*dbdrp)*grcgtp
        bdgxcgy(i1,k1) = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp             
        bmag(i1,k1) = b
        jac(i1,k1) = jacp*jfnp
        bdgrzn(i1,k1) = q0*br0/radiusp/b*psipp*grcgtp/jfnp
        gn0e(i1) = gn0ep
        gt0e(i1) = gt0ep
        gt0i(i1) = gt0ip

        ggxdgy(i1,k1) = gxdgyp
        ggy2(i1,k1) = dydrp**2*grp**2 + (r0/q0*qhatp*gthp)**2 + 2*dydrp*r0/q0*qhatp*grdgtp
        ggx(i1,k1) = grp
        gsf(i1) = gsfp
        gthf(i1,k1) = gthfp
     end do
  end do

  iseed = -(1777+myid*13)
  idum = ran2(iseed)
  phi = 0.
  camp = 0.

!  !$acc update device(phi)
!  !$acc update device(apar)
  
  if(myid.eq.master)then
     write(*,*)zfnth(ntheta),thfnz(ntheta/2),thfnz(ntheta/2+1)
     if(myid.eq.master)open(9, file='plot', &
          status='unknown',position='append')
     write(9,*)'dt,beta= ',dt, beta
     write(9,*)'amp = ',amp
     write(9,*)'peritr,ifluid= ',peritr,ifluid
     write(9,*)'tor,nonlin(1:2),nonline= ',tor,nonlin(1),nonlin(2),nonline
     write(9,*)'isuni= ',isuni, 'amie= ',amie
     write(9,*)'kxcut,kycut,bcut,wecut= ',kxcut,kycut,bcut,wecut
     write(9,*)'fradi,isg= ',fradi,isg
     write(9,*)'llk,mlk,onemd =',llk,mlk,onemd
     write(9,*)'vcut= ', vcut
     write(9,*)'rneu= ', rneu
     write(9,*)'V-ExB switch= ', vexbsw
     write(9,*)'V-parallel switch= ', vparsw
     write(9,*)'tot mmi,tot mme= ',mm1,mm2
     write(9,*)'pzcrite,encrit = ',pzcrite,encrit
     write(9,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
     write(9,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
     write(9,*) 't0,kyrhoi_local=',t0i(nr/2),2*pi*sqrt(mims(1))*sqrt(t0i(nr/2))/ly
     write(9,*) 'isonew,nzcrt,nzsrc=',isonew,nzcrt,nzsrc      
     close(9)
  end if

  if(myid.eq.master)then
     write(*,*)'dt,beta= ',dt, beta
     write(*,*)'amp,vpp,yd0 = ',amp,vpp,yd0
     write(*,*)'peritr,ifluid= ',peritr,ifluid
     write(*,*)'tor,nonlin= ',tor,nonlin(1),nonlin(2)
     write(*,*)'isuni= ',isuni, 'amie= ',amie
     write(*,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
     write(*,*)'fradi,isg= ',fradi,isg
     write(*,*)'llk,mlk,onemd =',llk,mlk,onemd
     write(*,*)'vcut= ', vcut
     write(*,*)'rneu= ', rneu
     write(*,*)'V-ExB switch= ', vexbsw
     write(*,*)'V-parallel switch= ', vparsw
     write(*,*)'mm1= ',mm1
     write(*,*)'pzcrite,encrit = ',pzcrite,encrit
     write(*,*)'nue0 = ',nue0(1),nue0(nr/2),nue0(nr-1)
     write(*,*)'xn0e(1),xnir0 = ',xn0e(1),xnir0
     write(*,*)'frequ, eru = ', frequ, eru
     write(*,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
     write(*,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
     write(*,*) 't0,kyrhoi_local=',t0i(nr/2),2*pi*sqrt(mims(1))*sqrt(t0i(nr/2))/ly
     write(*,*) 'coefu = ', xu**2*frequ
     write(*,*) 'ktheta*rhos = ',2*pi*sqrt(mims(1))*sqrt(t0e(nr/2))/ly
     write(*,*) 'cs/a, q0, q0p, s^hat = ',sqrt(t0e(nr/2)/2.)/a, q0, q0p, q0p/q0*r0
     write(*,*) 'rho* = rhos/a = ', sqrt(mims(1))*sqrt(t0e(nr/2))/a
     write(*,*) 'f0p,psip(nr/2),Bunit,candyf0p = ',f0p,psip(nr/2),bunit,candyf0p
     write(*,*) 'lxa min = ', ly*q0/(2*pi*r0*q0p)/a
     write(*,*) 't0i(nr/2)= ', t0i(nr/2)
     write(*,*) 'Gyrokrs = ', 2*pi*sqrt(mims(1))*sqrt(t0e(nr/2))/ly/bunit
  end if
  close(115)


# 563
!$acc update device( bfld,qhat,radius,gr,gth,grdgt,grcgt)
!$acc update device( gxdgy,dydr,dbdr,dbdth,dqhdr,jacob)
!$acc update device( yfn,hght,thflx) 
!$acc update device( rmaj,rmajp,elon,selon,tria,stria, psi)
!$acc update device( f,psip,sf,jacoba,jfn,zfnth,thfnz)
!$acc update device( t0i,t0e,t0b,t0c,t0ip,t0ep,t0bp,t0cp)
!$acc update device( xn0i,xn0e,xn0c,xn0b,xn0ip,xn0ep,xn0bp)
!$acc update device( xn0cp,vpari,vparc,vparb)
!$acc update device( vparip,vparcp,vparbp)
!$acc update device( capti,capte,captb,captc,capni,capne)
!$acc update device( capnb,capnc,zeff,nue0,phinc,phincp)
!$acc update device( er,upari,dipdr)
!$acc update device( e1gx,e1gy,e2gx,e2gy,bdge1gx,bdge1gy,bdge2gx,bdge2gy)
!$acc update device( psip2)
!$acc update device( curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz,bdcrvb)
!$acc update device( t0s,xn0s,capts,capns,vpars,vparsp)
!$acc update device( cn0s,n0smax,tgis)




  !      return
# 585
end subroutine init

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine gkps(nstep,ip)   

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none

  INTEGER :: ns
  real :: b,gam0,gam1,delyz,th,bf,dum,r,qr,shat
  real,DIMENSION(1:5) :: b1,b2,ter
  real :: kx1,kx2,ky
  real,dimension(:),allocatable :: akx,aky
  real :: sgnx,sgny,sz,myfe
  INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx
  COMPLEX :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
  COMPLEX :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
  COMPLEX :: sl(1:imx-1,0:jcnt-1,0:1)
  COMPLEX :: aphik(0:nxpp-1),myaphik(0:nxpp-1)
  real :: myaph(0:nxpp),aph(0:nxpp),u(0:imx,0:jmx,0:1)
  complex :: cdum
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp,grdgtp,gthp
  real :: wx0,wx1,wz0,wz1
  character(len=70) fname
  character(len=5) holdmyid

  save ifirst,akx,aky
  if(idg==1)write(*,*)'enter gkps'
  !   now do field solve...

  !      phi = 0.
  !      return
  !      temp3dxy = 0.
  aphik = 0.
  myaphik = 0.

  !  find rho(kx,ky)
  do i = 0,nxpp
     do j = 0,jm
        do k = 0,1
           u(i,j,k) = ran2(iseed)*nstep
        end do
     end do
  end do
  call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)
  if(idg==1)write(*,*)'pass first fft'

  temp3dxy = 0.

!$omp parallel do private(k,j,myj)
  do i = 0,jcnt*2-1
     k = int(i/jcnt)
     j = i-k*jcnt
        myj=jft(j)
        sl(1:imx-1,j,k) = v(1:imx-1,j,k)
        call ZGETRS('N',imx-1,1,mxg(:,:,j,k),imx-1,ipivg(:,j,k), &
             sl(:,j,k),imx-1,INFO) 
        temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
        temp3dxy(0,myj,k) = 0.
  end do

  !  from rho(kx,ky) to phi(kx,ky)
  do k = 0,1
     do j = 0,jmx-1
        do i = 0,imx-1
           temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx !*formphi(i,j,k)
        end do
     end do
  end do

  !  from phi(kx,ky) to phi(x,y)
  do k=0,mykm
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3dxy(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

100 continue
  do i = 0,nxpp-1
     do j = 0,jm-1
        do k = 0,mykm
           phi(i,j,k) = temp3dxy(i,j,k)
        end do
     end do
  end do

  !    x-y boundary points 

  if(idg==1)write(*,*)'pass enfz', myid

  !      return
end subroutine gkps

!      End of gkps....

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine spec(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: i,j,k,l,m,n
  real :: pf,efe,efi,pfi,efc,pfc,efb,pfb,pflxgb,eflxgb,x
  real :: pf_em,efe_em,efi_em,pfi_em,efc_em,pfc_em,efb_em,pfb_em
  real :: tdum,v(0:imx),denz(0:imx)
  real :: dtmp(0:imx,0:jmx,0:1),titmp(0:imx,0:jmx,0:1),tetmp(0:imx,0:jmx,0:1)


  tdum = tcurr-dt
  if(myid.eq.master)then
     open(9, file='plot', status='unknown',position='append')

     write(*,10)tdum,rmsphi(n)

10   format(1x,f10.1,2(2x,e15.8))
     write(9,10)tdum,rmsphi(n)
  endif

end subroutine spec

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

real function ran2(idum)

  parameter( IM1=2147483563,  &
       IM2=2147483399, &
       AM=1.0/IM1,&
       IMM1=IM1-1,&
       IA1=40014,&
       IA2=40692,&
       IQ1=53668,&
       IQ2=52774,&
       IR1=12211,&
       IR2=3791,&
       NTAB=32,&
       NDIV=1+IMM1/NTAB,&
       EPS=1.2e-7,&
       RNMX=1.0-EPS &
       )

  integer :: j,k,idum2=123456789,iy=0,iv(0:NTAB-1)
  real :: temp

  save idum2, iy,iv
  !      write(*,*)'idum2,iy  ',idum2,iy
  if(idum.le.0)then
     if(-idum.lt.1)then
        idum=1
     else
        idum = -idum
     end if
     idum2 = idum
     do j = NTAB+7,0,-1
        k = idum/IQ1
        idum = IA1*(idum-k*IQ1)-k*IR1
        if(idum.lt.0)idum = idum+IM1
        if(j.lt.NTAB)iv(j) = idum
     end do
     iy = iv(0)
  end if

  k = idum/IQ1
  idum = IA1*(idum-k*IQ1)-k*IR1
  if(idum.lt.0)idum = idum+IM1
  k = idum2/IQ2
  idum2 = IA2*(idum2-k*IQ2)-k*IR2
  if(idum2.lt.0)idum2 = idum2+IM2
  j = iy/NDIV
  iy = iv(j)-idum2
  iv(j) = idum
  if(iy<1)iy = iy+IMM1
  temp = AM*iy
  if(temp>RNMX)then
     ran2 = RNMX
  else
     ran2 = temp
  end if
  return
end function ran2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine initialize
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none

!  include "fftw3.f03"
  real :: dum,dum1,dum2,jacp,xndum,r,wx0,wx1
  !        complex(8),dimension(0:1) :: x,y
  real,dimension(0:1) :: x,y
  integer :: n,i,j,k,ip,iret

  call init
  if(idg==1)then
     do i = 0,last
        if(myid==i)then
           open(9, file='plot',status='unknown',position='append')
           write(9,*)'pass init',myid
           close(9)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     end do
  end if

  call mpi_barrier(mpi_comm_world,ierr)
  
  dum = 0.
  dum1 = 0.
  do i = 0,im-1
     r = xg(i)-0.5*lx+lr0
     j = int((r-rin)/dr)
     j = min(j,nr-1)
     wx0 = (rin+(j+1)*dr-r)/dr
     wx1 = 1.-wx0
     xndum = wx0*xn0e(j)+wx1*xn0e(j+1)
     dum = dum+(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
     dum1 = dum1+xndum*(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
  end do
  call MPI_ALLREDUCE(dum,jacp,1,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)

  call MPI_ALLREDUCE(dum1,dum2,1,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)

  totvol = dx*ly*dz*jacp    
  n0=float(tmm(1))/totvol
  n0e=mme*numprocs/totvol
  if(idg==1)then
     do i = 0,last
        if(myid==i)then
           open(9, file='plot',status='unknown',position='append')
           write(9,*)'totvol,jacp,dum2=',totvol,jacp,dum2,myid
           close(9)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     end do
  end if

  !  calculate the volume of each radial subdomain
  do k = 1,nsubd
     dum = 0.
     do i = (k-1)*im/nsubd,k*im/nsubd-1
        r = xg(i)-0.5*lx+lr0
        j = int((r-rin)/dr)
        j = min(j,nr-1)
        wx0 = (rin+(j+1)*dr-r)/dr
        wx1 = 1.-wx0
        dum = dum+(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
     end do
     call MPI_ALLREDUCE(dum,jacp,1,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
     vol(k) = dx*ly*dz*jacp    
  end do

  !     initialize particle quantities...
  if( cut.eq.0.) cut=1.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call ccfft('x',0,imx,0.0,tmpx,coefx,workx,0)
  call ccfft('y',0,jmx,0.0,tmpy,coefy,worky,0)
  call ccfft('z',0,kmx,0.0,tmpz,coefz,workz,0)
  call dsinf(1,x,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)

  ncurr = 1

  if(ifluid==1 .and. iperi==0)call gkps_init
end subroutine initialize
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine poisson(n,ip)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,i1,j,k,ip,it,iter=1
  real :: myrmsphi,rmp(20),myavap(0:imx-1)
  real :: myjaca(0:imx-1),jaca(0:imx-1)

  call gkps(n,ip)

  myrmsphi=0.
  do k=0,mykm-1
     do j=0,jm-1
        do i1=0,im-1
           myrmsphi=myrmsphi+phi(i1,j,k)*phi(i1,j,k)
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(myrmsphi,rmsphi(n),1, &
       MPI_REAL8,                               &
       MPI_SUM,TUBE_COMM,ierr)
  rmsphi(n)=sqrt(rmsphi(n)/(im*jm*km))

  if(idg.eq.1)write(*,*)'pass poisson'
!  !$acc update device(phi)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine poisson
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine reporter(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip

  if(mod(n,xnplt).eq.0) then
     call spec(n)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end subroutine reporter
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dcmpy(u,v)
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx,id
  INTEGER :: recvcnt(0:ntube-1)
  real :: u(0:imx-1,0:jmx-1,0:1)
  complex :: v(0:imx-1,0:jcnt-1,0:1),myv(0:imx-1,0:jcnt-1,0:1)
  complex :: sbuf(0:imx*jmx*2-1),rbuf(0:imx*jcnt*2-1)
  complex :: temp3d(0:imx-1,0:jmx-1,0:1)
  real :: kx,ky,kx0,th,shat,sgny

  do k=0,mykm
     do j=0,jm-1
        do i=0,imx-1
           temp3d(i,j,k)=u(i,j,k)
        enddo
     enddo
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)
        end do
     end do
  enddo

  do m = 0,jcnt-1
     myv(:,m,:) = temp3d(:,jft(m),:)
  end do

  cnt = 2*jcnt*imx
  call mpi_allreduce(myv,v,cnt,MPI_DOUBLE_COMPLEX,mpi_sum, &
       grid_comm,ierr)

  !      return
end subroutine dcmpy
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gkps_init

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none

  INTEGER :: ns
  real :: b,gam0,gam1,delyz,th,bf,dum,r,qr,shat
  real,DIMENSION(1:5) :: b1,b2,ter
  real :: kx1,kx2,ky
  real,dimension(:),allocatable :: akx,aky
  real,dimension(:,:,:,:,:),allocatable:: gamb1,gamb2
  real :: sgnx,sgny,sz,myfe
  INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx
  complex :: cdum
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp,grdgtp,gthp
  real :: wx0,wx1,wz0,wz1
  character(len=70) fname
  character(len=5) holdmyid

  if(idg==1)then
     do i = 0,last
        if(myid==i)then
           open(9, file='plot',status='unknown',position='append')
           write(9,*)'enter gkps_init',myid
           close(9)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     end do
  end if

  write(holdmyid,'(I5.5)') MyId
  fname='./matrix/'//'mx_phi_'//holdmyid
  !     form factors....
  if (igetmx==0) then
     allocate(akx(0:imx-1),aky(0:jcnt-1), &
          gamb1(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1), &
          gamb2(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1))
     if(idg==1)then
        do i = 0,last
           if(myid==i)then
              open(9, file='plot',status='unknown',position='append')
              write(9,*)'after allocate gkps_init',myid
              close(9)
           end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do
     end if


     do l=0,im-1
        do m=0,jcnt-1
           j = tclr*jcnt+m
           do i=0,im-1 
              r = lr0-lx/2+i*dx
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              do n = 0,1
                 k = gclr*kcnt+n
                 k1 = int(k*dz/delz)
                 k1 = min(k1,ntheta-1)
                 wz0 = ((k1+1)*delz-dz*k)/delz
                 wz1 = 1-wz0
                 th = wz0*thfnz(k1)+wz1*thfnz(k1+1)
                 k = int((th+pi)/dth)
                 k = min(k,ntheta-1)
                 wz0 = (-pi+(k+1)*dth-th)/dth
                 wz1 = 1.-wz0
                 bfldp = wx0*wz0*bfld(i1,k)+wx0*wz1*bfld(i1,k+1) &
                      +wx1*wz0*bfld(i1+1,k)+wx1*wz1*bfld(i1+1,k+1) 
                 dydrp = wx0*wz0*dydr(i1,k)+wx0*wz1*dydr(i1,k+1) &
                      +wx1*wz0*dydr(i1+1,k)+wx1*wz1*dydr(i1+1,k+1) 
                 qhatp = wx0*wz0*qhat(i1,k)+wx0*wz1*qhat(i1,k+1) &
                      +wx1*wz0*qhat(i1+1,k)+wx1*wz1*qhat(i1+1,k+1) 
                 grp = wx0*wz0*gr(i1,k)+wx0*wz1*gr(i1,k+1) &
                      +wx1*wz0*gr(i1+1,k)+wx1*wz1*gr(i1+1,k+1) 
                 gthp = wx0*wz0*gth(i1,k)+wx0*wz1*gth(i1,k+1) &
                      +wx1*wz0*gth(i1+1,k)+wx1*wz1*gth(i1+1,k+1) 
                 gxdgyp = wx0*wz0*gxdgy(i1,k)+wx0*wz1*gxdgy(i1,k+1) &
                      +wx1*wz0*gxdgy(i1+1,k)+wx1*wz1*gxdgy(i1+1,k+1) 
                 grdgtp = wx0*wz0*grdgt(i1,k)+wx0*wz1*grdgt(i1,k+1) &
                      +wx1*wz0*grdgt(i1+1,k)+wx1*wz1*grdgt(i1+1,k+1) 

                 m1 = mstart+int((float(m)+1.0)/2)
                 if(m==0)m1=0
                 sgny = isgnft(m)

                 ky=sgny*2.*pi*float(m1)/ly
                 kx1=pi*float(l)/lx
                 kx2=-pi*float(l)/lx
                 bf=bfldp

                 do ns = 1, nsm
                    b1(ns)=mims(ns)*(kx1*kx1*grp**2 + &
                         ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                         +2*dydrp*lr0/q0*qhatp*grdgtp) &
                         +2*kx1*ky*gxdgyp)/(bf*bf)*ter(ns)/(q(ns)*q(ns))

                    b2(ns)=mims(ns)*(kx2*kx2*grp**2 + &
                         ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                         +2*dydrp*lr0/q0*qhatp*grdgtp) &
                         +2*kx2*ky*gxdgyp)/(bf*bf)*ter(ns)/(q(ns)*q(ns))

                    call MPI_BARRIER(MPI_COMM_WORLD,ierr)                    

                    call srcbes(b1(ns),gam0,gam1)
                    gamb1(ns,l,m,i,n)=gam0
                    call srcbes(b2(ns),gam0,gam1)
                    gamb2(ns,l,m,i,n)=gam0
                 enddo

              enddo
           enddo
        enddo
     enddo
     
     if(idg==1)then
        do i = 0,last
           if(myid==i)then
              open(9, file='plot',status='unknown',position='append')
              write(9,*)'before assembling mxg',myid
              close(9)
           end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do
     end if

     do k = 0,1
        do j = 0,jcnt-1
           do i = 1,imx-1
              r = lr0-lx/2+i*dx
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              do ix = 1,imx-1
                 mxg(i,ix,j,k) = 0.0
                 if(i==ix)mxg(i,ix,j,k) = fradi*gn0e(i)*cn0e/gt0e(i)
                 do ikx = 0,imx-1
                    do ns = 1, nsm
                       mxg(i,ix,j,k) = mxg(i,ix,j,k)+q(ns)*sin(ix*ikx*pi/imx)* &
                            ((1-gamb1(ns,ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                            (1-gamb2(ns,ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                            /ter(ns)*cn0s(ns)*gn0s(ns,i)/(IU*imx)
                    end do
                 end do
              end do
           end do
        end do
     end do
     if(idg==1)then
        do i = 0,last
           if(myid==i)then
              open(9, file='plot',status='unknown',position='append')
              write(9,*)'before zgetrf',myid
              close(9)
           end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do
     end if

!$omp parallel do private(k,j)
     do i = 0,jcnt*2-1
        k = int(i/jcnt)
        j = i-k*jcnt
        call ZGETRF(imx-1,imx-1,mxg(:,:,j,k),imx-1,ipivg(:,j,k),INFO )
     end do
     if(idg==1)then
        do i = 0,last
           if(myid==i)then
              open(9, file='plot',status='unknown',position='append')
              write(9,*)'after zgetrf',myid
              close(9)
           end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do
     end if

 open(10000+MyId,file=fname,form='unformatted',status='unknown')
 do i = 1,imx-1
    do m = 0,jcnt-1
       do k = 0,1
          write(10000+MyId)(mxg(i,j,m,k),j=1,imx-1),ipivg(i,m,k)
       end do
    end do
 end do
 close(10000+myid)

endif

if(igetmx.eq.1) then
 open(10000+MyId,file=fname,form='unformatted',status='old')
 do i = 1,imx-1
    do m = 0,jcnt-1
       do k = 0,1
          read(10000+MyId)(mxg(i,j,m,k),j=1,imx-1),ipivg(i,m,k)
       end do
    end do
 end do
 close(10000+myid)
end if

  if(idg==1)write(*,*)'pass form factors'
end subroutine gkps_init

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
