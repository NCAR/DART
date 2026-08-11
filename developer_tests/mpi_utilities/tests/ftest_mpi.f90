! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program ftest_mpi

! simple MPI fortran program.  use to test running interactively
! with MPI parallel communication libraries.  warning -- this program
! may compile without obvious errors, but at runtime, unless MPI_Init()
! returns 0 as the error code, there is a good chance the compile and
! link phase did not succeed.

! Most fortran MPI implementations provide either a fortran 90 module
! which defines the interfaces to the MPI library routines, or an include
! file which defines constants.  Try to use the module if it is available.

use mpi

implicit none

!include "mpif.h"

! integer variables
integer :: ierror, myrank, totalprocs, rc

   print *, "program start"

   ierror = -999
   call MPI_Init(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Init() did not succeed, error code = ", ierror
      print *, "If error code is -999, the most likely problem is that"
      print *, "the right MPI libraries were not found at compile time."
      stop
   endif

   print *, "MPI initialized successfully"

   myrank = -1
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Comm_rank() did not succeed, error code = ", ierror
      stop
   endif
   print *, "My MPI rank is: ", myrank

   totalprocs = -1
   call MPI_Comm_size(MPI_COMM_WORLD, totalprocs, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Comm_size() did not succeed, error code = ", ierror
      stop
   endif
   print *, "Total MPI tasks: ", totalprocs

   ! This is not really an MPI test, but we do use execute_command_line to
   ! start model advances in async=2 and async=4 modes
   call do_execute_command_line("echo hello world", rc)
   if (rc /= 0) print *, 'call to execute_command_line() returned error'

   ierror = -999
   call MPI_Finalize(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Finalize() did not succeed, error code = ", ierror
      stop
   endif

   print *, "All MPI calls succeeded, test passed."
   print *, "program end"

contains

!> wrapper so you only have to make this work in a single place

subroutine do_execute_command_line(execute, rc)

character(len=*), intent(in)  :: execute
integer,          intent(out) :: rc

integer :: exitstat

    call execute_command_line(trim(execute), exitstat=exitstat)

    rc = exitstat

end subroutine do_execute_command_line


end program ftest_mpi
