! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program ftest_sendrecv_big

! MPI fortran program that uses parts of the DART library to test
! the send and receive functions in the mpi_utilities_mod.f90 file.
! THIS IS NOT A STANDALONE PROGRAM!


! this module should be on your system if MPI is installed
use mpi

! these are DART modules which must be built before running this test.
use types_mod
use utilities_mod
use time_manager_mod
use mpi_utilities_mod

implicit none

! some older installations only installed a header file instead of
! a module (module is better).  only if there is no other option,
! comment out 'use mpi' above and comment in the include file here.
!include "mpif.h"


integer(i8), parameter :: BSIZE = 3350074489_i8

real(r8), allocatable :: buf(:)

! integer variables
integer :: myrank, totalprocs

   call initialize_mpi_utilities('ftest_sendrecv_big')

   allocate(buf(BSIZE))
   myrank = my_task_id()
   print *, "My MPI rank is: ", myrank

   totalprocs = task_count()
   if (myrank == 0) print *, "Total MPI tasks: ", totalprocs
   if (myrank == 0) print *, "Exchanging buffer of size ", size(buf, kind=i8), "between tasks 0 and 1"
   if (myrank == 0) then
      buf(:) = 1
      call send_to(1, buf)
   else if (myrank == 1) then
      buf(:) = 0
      call receive_from(0, buf)
      if (any(buf == 0)) then
         print *, 'NOT EVERYTHING WAS SENT'
         stop
      endif
   endif

   call finalize_mpi_utilities()

end program ftest_sendrecv_big

