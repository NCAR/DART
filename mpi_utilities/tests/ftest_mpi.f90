! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ftest_mpi

! simple MPI fortran program.  use to test running interactively
! with MPI parallel communication libraries.  warning -- this program
! may compile without obvious errors, but at runtime, unless MPI_Init()
! returns 0 as the error code, there is a good chance the compile and
! link phase did not succeed.

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

 ! interface block for getting return code back from system() routine
 interface
  function system(string)
   character(len=*) :: string
   integer :: system
  end function system
 end interface
 ! end block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

   ! This is not really an MPI test, but we do use the system() function to
   ! start model advances in async=2 and async=4 modes, and to get this
   ! program to compile may involve editing the mpi_utilities module in dart.
   rc = system("echo hello world " // char(0))

   ierror = -999
   call MPI_Finalize(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Finalize() did not succeed, error code = ", ierror
      stop
   endif

   print *, "All MPI calls succeeded, test passed."
   print *, "program end"

end program ftest_mpi

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
