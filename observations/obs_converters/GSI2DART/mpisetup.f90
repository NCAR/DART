module mpisetup

! Adapted from .../src/enkf/mpisetup.f90

!$$$  module documentation block
!
! module: mpisetup                     initialize MPI.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract: Initialize and finalize MPI, create variables numproc (total # of MPI tasks,
!  same for each task) and nproc (MPI task #, different for each task). MPI
!  subroutine names and constants are imported via the 'mpif.h' include file.
!
! Public Subroutines
!   mpi_initialize: Initialize MPI.
!   mpi_cleanup: Wait for all tasks, then finalize MPI.
!   all the MPI routines (included via the mpif.h include file).
!
! Public Variables:
!   nproc:  MPI task # (root = 0).
!   numproc:  total number of MPI tasks.
!   mpi_status:  MPI status
!   mpi_realkind: double/single precision. not explicitly used
!
! program history log:
!   2009-02-23  Initial version.
!
! attributes:
!   language: f95
!
!$$$

use kinds, only: r_kind, r_single, r_double
use mpi_utilities_mod, only : task_count, my_task_id, initialize_mpi_utilities, finalize_mpi_utilities
implicit none
! mpi definitions.
! include 'mpif.h'
integer, public :: numproc, nproc
! integer, public :: mpi_status(mpi_status_size)
! integer, public :: mpi_realkind

contains

subroutine mpi_initialize()
integer :: ierr
! call mpi_init(ierr)

call initialize_mpi_utilities('gsi_to_dart')
! nproc is process number, numproc is total number of processes.
! call mpi_comm_rank(mpi_comm_world,nproc,ierr)
! call mpi_comm_size(mpi_comm_world,numproc,ierr)
numproc = task_count()
nproc = my_task_id()

if (nproc == 0) print *,'running on ',numproc,' processors ...'
if (r_kind == r_single) then
   ! mpi_realkind = mpi_real4
   ! if (nproc == 0) print *,'mpi_realkind set to mpi_real4',mpi_real4
else if (r_kind == r_double) then
   ! mpi_realkind = mpi_real8
   ! if (nproc == 0) print *,'mpi_realkind set to mpi_real8',mpi_real8
else
   print *,'ERROR: illegal r_kind (must be single or double)'
   call mpi_cleanup()
   ! call finalize_mpi_utilities()
endif

end subroutine mpi_initialize

subroutine mpi_cleanup()
integer :: ierr
flush(6,err=10)
flush(0,err=10)
10 continue
!call mpi_barrier(mpi_comm_world,ierr)
if (nproc == 0) write(6,*) 'all done!'
!call mpi_finalize(ierr)
call finalize_mpi_utilities()
! if (ierr /= 0) then
!  print *, 'MPI_Finalize error status = ',ierr
! end if
end subroutine mpi_cleanup

end module mpisetup
