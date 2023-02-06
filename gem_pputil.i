# 1 "gem_pputil.f90"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE gem_pputil
!
!  use fft_wrapper

# 7
  use openacc

# 9
  use mpi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ppinit_mpi,ppinit_decomp,ppexit


  INTEGER, SAVE :: me, nvp,npp,GCLR,TCLR
  INTEGER, SAVE :: TUBE_COMM,GRID_COMM
!
CONTAINS

!======================================================================
!mpi initialization
subroutine ppinit_mpi(idproc,nproc)
  use mpi
  integer, intent(out) :: idproc,nproc
  integer :: ierr,npp

  call mpi_init(ierr)
  if (ierr.ne.MPI_SUCCESS) then
    write(*,*) 'problem with mpi_init: ierr=',ierr
    stop
  endif
  !set the size of the MPI in the mpi_comm_world MPI space
  call mpi_comm_size(mpi_comm_world,npp,ierr)
  if (ierr.ne.MPI_SUCCESS) then
    write(*,*) 'problem with mpi_comm_size: ierr=',ierr
    stop
  endif
  !set the rank of the MPI in the mpi_comm_world MPI space
  call mpi_comm_rank(mpi_comm_world,me,ierr)
  if (ierr.ne.MPI_SUCCESS) then
    write(*,*) 'problem with mpi_comm_rank: ierr=',ierr
    stop
  endif
  nproc=npp
  idproc=me 
end subroutine ppinit_mpi
!particle decomp
subroutine ppinit_decomp(mype,nproc,ntube,com1,com2)
  use mpi
  integer, intent(in) :: mype,nproc,ntube
  integer, intent(out) :: com1,com2
  integer :: ierr

  !the rank of the particle decomposition
  gclr=int(mype/ntube)
  !the rank in the particle decomposition
  tclr=mod(mype,ntube)

  !set the particle MPI space
  call mpi_comm_split(mpi_comm_world,gclr,tclr,grid_comm,ierr)
  call mpi_comm_split(mpi_comm_world,tclr,gclr,tube_comm,ierr)

  com1=tube_comm
  com2=grid_comm
  nvp=nproc/ntube
end subroutine ppinit_decomp
!===========================================================================
  SUBROUTINE ppexit
    INTEGER :: ierr
    CALL MPI_FINALIZE(ierr)
    STOP
  END SUBROUTINE ppexit
!
!===========================================================================

END MODULE gem_pputil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
