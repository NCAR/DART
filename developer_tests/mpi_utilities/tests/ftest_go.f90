! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ftest_go

! test using named pipes to synchronize two separate executables.
! used in the async=4 mode of running a parallel filter and running
! a parallel model advance, where filter does not exit during the
! model advances.

! run this in background and then run ftest_go to wake it up.

! The following 2 build tips are the 2 places where different installations
! of MPI seem to vary the most.  Some systems have an include file, some
! have a F90 module.  Some require an interface block to use the system()
! function, some give an error if it is here.   You can use this program
! to figure out which combinations work on your system.  Then go into the
! $DART/mpi_utilities and make the same two changes in mpi_utilities_mod.f90,
! and just the system() change (if needed) in null_mpi_utilities_mod.f90.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BUILD TIP 1:
! Most fortran MPI implementations provide either a fortran 90 module
! which defines the interfaces to the MPI library routines, or an include
! file which defines constants.  Try to use the module if it is available.

use mpi

implicit none

!include "mpif.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BUILD TIP 2:
! Some systems require this interface block in order to use the system()
! function.  However, some other systems complain if this is here... 
! If this is a problem your program will not link and most likely give 
! you an error about an undefined symbol (something like '_system_').  
! Comment this block in or out as needed.

!  ! interface block for getting return code back from system() routine
!  interface
!   function system(string)
!    character(len=*) :: string
!    integer :: system
!   end function system
!  end interface
!  ! end block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! integer variables
integer :: ierror, myrank, totalprocs

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

   call restart_task()

   ierror = -999
   call MPI_Finalize(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Finalize() did not succeed, error code = ", ierror
      stop
   endif

   print *, "All MPI calls succeeded, test passed."
   print *, "program end"

contains

!-----------------------------------------------------------------------------
subroutine block_task()

! block by reading a named pipe file until some other task
! writes a string into it.  this ensures the task is not
! spinning and using CPU cycles, but is asleep waiting in
! the kernel.   one subtlety with this approach is that even
! though named pipes are created in the filesystem, they are
! implemented in the kernel, so on a multiprocessor machine
! the write into the pipe file must occur on the same PE as
! the reader is waiting.  see the 'wakeup_filter' program for
! the MPI job which spreads out on all the PEs for this job
! and writes into the file from the correct PE.

character(len = 32) :: fifo_name
integer :: rc
logical :: verbose

   verbose = .true.


   ! the i5.5 format below will not handle task counts larger than this.
   if (totalprocs > 99999) then
      print *, 'cannot handle task counts > 99999'
      call exit
   endif
   
   if (verbose) write(*,*) 'putting to sleep task id ', myrank


   ! if you change this in any way, change the corresponding string in 
   ! restart_task() below.

   write(fifo_name, '(a, i5.5)') "filter_lock", myrank
   
   if (verbose) write(*,*) 'removing any previous lock file: '//trim(fifo_name)
   rc = system('rm -f '//trim(fifo_name)//' '//char(0))

   if (verbose) write(*,*) 'made fifo, named: '//trim(fifo_name)
   rc = system('mkfifo '//trim(fifo_name)//' '//char(0))
   
   if (verbose) write(*,*) 'ready to read from lock file: '//trim(fifo_name)
   rc = system('cat < '//trim(fifo_name)//' '//char(0))
   
   if (verbose) write(*,*) 'got response, removing lock file: '//trim(fifo_name)
   rc = system('rm -f '//trim(fifo_name)//' '//char(0))

end subroutine block_task

!-----------------------------------------------------------------------------
subroutine restart_task()


character(len = 32) :: fifo_name
integer :: rc
logical :: verbose

   verbose = .true.

   ! the i5.5 format below will not handle task counts larger than this.
   if (totalprocs > 99999) then
      print *, 'cannot handle task counts > 99999'
      call exit
   endif

   if (verbose) write(*,*) 'waking up task id ', myrank

   write(fifo_name,"(a,i5.5)") "filter_lock", myrank

   if (verbose) write(*,*) 'ready to write to lock file: '//trim(fifo_name)
   rc = system('echo restart > '//trim(fifo_name)//' '//char(0))
   
   if (verbose) write(*,*) 'response was read from lock file: '//trim(fifo_name)
   

end subroutine restart_task

end program ftest_go

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
