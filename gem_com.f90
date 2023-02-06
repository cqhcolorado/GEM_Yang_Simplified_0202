module gem_com

  !common data used for gem

  use mpi
  use, intrinsic :: iso_c_binding
  use gem_pputil

  implicit none

  INTERFACE
     real function revers(num,n)
     end function revers

     real function ran2(i)
     end function ran2

     real function en3(s)
       real :: s
     end function en3
  END INTERFACE

  REAL :: poisson_start_tm=0,poisson_end_tm=0,poisson_tot_tm=0
    
  real :: poisson0_start_tm=0,poisson0_end_tm=0,poisson0_tot_tm=0
  real :: total_tm=0.0
  real :: startinit=0.0, lastinit=0.0, totinit=0.0, init_tm=0.0
  real :: poisson_tm
  
  integer :: imx,jmx,kmx,mmx,mmxe,nmx,nsmx,nsubd=8,&
       ntube,nxpp,ngdx=5,nb=6, &
       negrd=16,nlgrd=16

  character(len=70) outname
  REAL :: endtm,begtm,pstm
  REAL :: starttm,lasttm,tottm
  real :: aux1(50000),aux2(20000)
  real,dimension(:),allocatable :: workx,worky,workz,xsinin,xsinout
  complex,dimension(:),allocatable :: tmpx,xin,xout,yin,yout,zin,zout
  complex,dimension(:),allocatable :: tmpy
  complex,dimension(:),allocatable :: tmpz
  complex,dimension(:,:),allocatable :: tmpxyin,tmpxyout
  type(C_PTR) :: planx,iplanx,plany,iplany,planz,iplanz,plansinx

  integer :: icgp,jcgp,cgpfacx,cgpfacy,nxsrc,nesrc
  integer :: mme,mmb
  REAL, dimension(:,:),allocatable :: rwx,rwy
  INTEGER,dimension(:),allocatable :: mm,tmm,lr
  integer :: micell,mecell !jycheng
  integer :: nonlin1,nonlin2 !jycheng
  REAL,dimension(:),allocatable :: tets,mims,q
  REAL,dimension(:),allocatable :: kapn, kapt
  INTEGER :: timestep,im,jm,km,mykm,iseed,nrst,nfreq,isft,mynf,ifskp,iphbf,iapbf,idpbf
  real,dimension(:),allocatable :: time
  REAL :: dx,dy,dxcgp,dycgp,dz,dxsrc,desrc,ecutsrc,pi,pi2,dt,dte,totvol,n0,n0e,tcurr,rmpp,rmaa,eprs
  REAL :: lx,ly,lz,xshape,yshape,zshape,pzcrit(5),pzcrite,encrit
  INTEGER :: nm,nsm,kcnt,jcnt,ncurr,llk,mlk,onemd,iflr,iorb
  integer :: izonal,adiabatic_electron,ineq0,iflut,nlow,ntor0,mstart,iexb
  REAL :: cut,amp,tor,amie,isg,rneu,rneui,emass,qel,mbeam,qbeam,teth,vexbsw,vparsw,gammah,gamgtc,gamtoy,ghzon,gamgyro
  REAL :: c4,fradi,kxcut,kycut,bcut,wecut,ftrap,adwn,adwe,adwp,frmax
  INTEGER :: iput,iget,igetmx,idg,kzlook,ision,isiap,peritr,iadi,ipred,icorr,jpred,jcorr
  REAL,DIMENSION(:,:),allocatable :: yyamp,yyre,yyim
  complex,dimension(:,:),allocatable :: camp,campf
  REAL :: br0,lr0,qp,e0,vcut,vpp,vt0,yd0,cv,vsphere
  integer :: nonlin(5),nonline,ipara,isuni,isunie,ifluid,nopz,nopi(5),nowi(5),noen,nowe,novpar,isonew
  integer :: ipbm,ias 
  complex :: IU
  real,dimension(:),allocatable :: coefx,coefy,coefz
  complex,dimension(1:8) :: apk,ptk,dpdtk
  integer,dimension(1:8) :: lapa,mapa,napa
  real :: mrtio(0:1),aven,avptch
  integer :: icrs_sec,ipg,isphi
  integer,dimension(0:255) :: isgnft,jft

  REAL,DIMENSION(:,:,:,:),allocatable :: den
  real,dimension(:,:,:),allocatable :: phi
  REAL,DIMENSION(:),allocatable :: xg,yg,zg
  real,dimension(:,:),allocatable :: cfx,cfy,jac,bmag,bdgxcgy,bdgrzn,ggxdgy,ggy2,ggx,gthf
  real,dimension(:),allocatable :: gn0e,gt0e,gt0i,avap,dtez,gsf
  real,dimension(:,:),allocatable :: gn0s,gt0s,dtiz

  
  !              Various diagnostic arrays and scalars
  !    plotting constants

  INTEGER :: nplot,xnplt,imovie=1000000000,nzcrt,npze,npzi,npzc,npzb,nzsrc,nzgyro
  REAL :: contu,wmax

  !    energy diagnostic arrays

  REAL,DIMENSION(:,:),allocatable :: ke
  REAL,DIMENSION(:),allocatable :: fe,te
  REAL,DIMENSION(:),allocatable :: rmsphi,rmsapa,avewe
  REAL,DIMENSION(:,:),allocatable :: nos,avewi

  !    flux diagnostics
  REAL,DIMENSION(:),allocatable :: vol

  !   kr, ktheta spectrum plots
  REAL,DIMENSION(:,:),allocatable :: phik

  !     weighty variables
  INTEGER,dimension(:),allocatable :: deljp,deljm
  INTEGER,dimension(:,:),allocatable :: jpl
  INTEGER,dimension(:,:),allocatable :: jpn
  INTEGER,dimension(:,:),allocatable :: jmi
  INTEGER,dimension(:,:),allocatable :: jmn
  REAL,DIMENSION(:),allocatable :: weightp,weightm
  REAL,DIMENSION(:),allocatable :: weightpn,weightmn


  complex,dimension(:,:,:,:),allocatable :: mxg,mxa,mxd
  integer,dimension(:,:,:),allocatable :: ipivg,ipiva,ipivd

  !      MPI variables
  !  include '/usr/include/mpif.h'

  integer,parameter :: Master=0
  integer :: numprocs,n_omp
  INTEGER :: MyId,Last,cnt,ierr
  INTEGER :: GRID_COMM,TUBE_COMM
  INTEGER :: GCLR,TCLR,GLST,TLST
  INTEGER :: stat(MPI_STATUS_SIZE)
  INTEGER :: lngbr,rngbr,idprv,idnxt

  character(len=*) directory
  parameter(directory='./dump/')

  character(len=*) outdir
  parameter(outdir='./out/')

  !real :: ran2,revers
  !integer :: mod
  !real :: amod
  save

contains

  subroutine new_gem_com()
    nxpp = imx !/ntube
    allocate(workx(4*imx),worky(4*jmx),workz(4*kmx),xsinin(imx),xsinout(imx))
    allocate(tmpx(0:imx-1),xin(imx),xout(imx),yin(jmx),yout(jmx),zin(kmx),zout(kmx))
    allocate(tmpy(0:jmx-1))
    allocate(tmpz(0:kmx-1))
    allocate(tmpxyin(1:imx,1:jmx),tmpxyout(1:imx,1:jmx))

    allocate(rwx(5,4),rwy(5,4))
    allocate(mm(5),tmm(5),lr(5))
    allocate(tets(5),mims(5),q(5))
    
#ifdef OPENACC
!$acc enter data create(lr,mims,q)
#endif

#ifdef OPENMP
#endif

    allocate(kapn(5),kapt(5))
    allocate(time(0:nmx))
    allocate(yyamp(0:jmx,0:4),yyre(0:jmx,0:4),yyim(0:jmx,0:4),camp(0:6,0:50000),campf(0:6,0:nfreq-1))
    allocate(coefx(100+8*imx),coefy(100+8*jmx),coefz(100+8*kmx))

    ALLOCATE( den(nsmx,0:nxpp,0:jmx,0:1))
    allocate( phi(0:nxpp,0:jmx,0:1))

!!$acc enter data create(phi)

    ALLOCATE( xg(0:nxpp),yg(0:jmx),zg(0:1))

    allocate( cfx(0:nxpp,0:1),cfy(0:nxpp,0:1),jac(0:nxpp,0:1))
    allocate( bmag(0:nxpp,0:1),bdgxcgy(0:nxpp,0:1),bdgrzn(0:nxpp,0:1),gthf(0:nxpp,0:1),gsf(0:nxpp))
    allocate( ggxdgy(0:nxpp,0:1),ggy2(0:nxpp,0:1),ggx(0:nxpp,0:1))
    allocate (gn0e(0:nxpp),gt0e(0:nxpp),gt0i(0:nxpp),avap(0:nxpp),dtez(0:imx))
    allocate (gn0s(1:5,0:nxpp),gt0s(1:5,0:nxpp),dtiz(1:5,0:imx))

    
!    Various diagnostic arrays and scalars
!    plotting constants

    ALLOCATE( rmsphi(0:nmx),rmsapa(0:nmx))
    ALLOCATE( nos(nsmx,0:nmx))

    !    flux diagnostics
    ALLOCATE( vol(1:nsubd))


    !     weighty variables
    ALLOCATE( deljp(0:nxpp),deljm(0:nxpp))
    ALLOCATE( jpl(0:nxpp,0:jmx))
    ALLOCATE( jpn(0:nxpp,0:jmx))
    ALLOCATE( jmi(0:nxpp,0:jmx))
    ALLOCATE( jmn(0:nxpp,0:jmx))
    ALLOCATE( weightp(0:nxpp),weightm(0:nxpp))
    ALLOCATE( weightpn(0:nxpp),weightmn(0:nxpp))

    allocate(mxg(imx-1,imx-1,0:jcnt-1,0:1),mxa(imx-1,imx-1,0:jcnt-1,0:1),mxd(imx-1,imx-1,0:jcnt-1,0:1), &
             ipivg(imx-1,0:jcnt-1,0:1),ipiva(imx-1,0:jcnt-1,0:1),ipivd(imx-1,0:jcnt-1,0:1))
  end subroutine new_gem_com

end module gem_com










