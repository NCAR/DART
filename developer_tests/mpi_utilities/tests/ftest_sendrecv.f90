! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program ftest_sendrecv

! MPI fortran program that uses parts of the DART library to test
! the send and receive functions in the mpi_utilities_mod.f90 file.
! THIS IS NOT A STANDALONE PROGRAM!

! The following 2 build tips are the 2 places where different installations
! of MPI seem to vary the most.  Some systems have an include file, some
! have a F90 module.  Some require an interface block to use the system()
! function, some give an error if it is here.   You can use this program
! to figure out which combinations work on your system.  Then go into 
! $DART/assimilation_code/utilities/mpi_utilities_mod.f90 and make the 
! same two changes there. also possibly (if needed) null_mpi_utilities_mod.f90.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BUILD TIP 1:
! Most fortran MPI implementations provide either a fortran 90 module
! which defines the interfaces to the MPI library routines, or an include
! file which defines constants.  Try to use the module if it is available.

use mpi

! these are DART modules which must be built before running this test.
use types_mod
use utilities_mod
use time_manager_mod
use mpi_utilities_mod

implicit none

!include "mpif.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BUILD TIP 2:
! Some systems require this interface block in order to use the system()
! function.  However, some other systems complain if this is here.
! The 'fixsystem' script in the assimilation_code/modules/utilities
! directory will try to fix this for you.  If it does not work,
! your program will not link and most likely give you an error 
! about an undefined symbol (something like '_system_').  
! Comment this block in or out as needed.

! interface block for getting return code back from system() routine

! !!SYSTEM_BLOCK_EDIT START COMMENTED_OUT
!   interface
!    function system(string)
!     character(len=*) :: string
!     integer :: system
!    end function system
!   end interface
! !!SYSTEM_BLOCK_EDIT END COMMENTED_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: BSIZE = 101

real(r8) :: buf(BSIZE)
integer :: i

! integer variables
integer :: myrank, totalprocs

   call initialize_mpi_utilities('ftest_sendrecv')

   myrank = my_task_id()
   print *, "My MPI rank is: ", myrank

   totalprocs = task_count()
   if (myrank == 0) print *, "Total MPI tasks: ", totalprocs

   do i = 1, BSIZE-1, 1
      if (myrank == 0) print *, "Exchanging ", i, " items between tasks 0 and 1"
      if (myrank == 0) then
         buf(:) = 1
         call send_to(1, buf(1:i))
      else if (myrank == 1) then
         buf(:) = 0
         call receive_from(0, buf(1:i))
         if (any(buf(1:i) == 0)) then
            print *, 'NOT EVERYTHING WAS SENT'
            stop
         endif
      endif
   enddo

   call finalize_mpi_utilities()

end program ftest_sendrecv

