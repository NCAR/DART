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

! the NAG compiler needs these special definitions enabled.
! the 'fixsystem' script in the assimilation_code/modules/utilities dir
! should fix this for you.  please leave the BLOCK comment lines unchanged.

! !!NAG_BLOCK_EDIT START COMMENTED_OUT
! !#ifdef __NAG__
!
! use F90_unix_proc, only : sleep, system, exit
!
! !! these are the calling sequences for NAG compiler
! !  PURE SUBROUTINE SLEEP(SECONDS,SECLEFT)
! !    INTEGER,INTENT(IN) :: SECONDS
! !    INTEGER,OPTIONAL,INTENT(OUT) :: SECLEFT
! !
! !  SUBROUTINE SYSTEM(STRING,STATUS,ERRNO)
! !    CHARACTER*(*),INTENT(IN) :: STRING
! !    INTEGER,OPTIONAL,INTENT(OUT) :: STATUS,ERRNO
! !
! !!also used in exit_all outside this module
! !  SUBROUTINE EXIT(STATUS)
! !    INTEGER,OPTIONAL :: STATUS
! !! end block
!
! !#endif
! !!NAG_BLOCK_EDIT END COMMENTED_OUT

implicit none

!include "mpif.h"



! interface block for getting return code back from system() routine
! the 'fixsystem' script in the assimilation_code/modules/utilities dir
! should fix this for you.  please leave the BLOCK comment lines unchanged.


! !!SYSTEM_BLOCK_EDIT START COMMENTED_OUT
! !#if .not. defined (__GFORTRAN__) .and. .not. defined(__NAG__)
! ! interface block for getting return code back from system() routine
! interface
!  function system(string)
!   character(len=*) :: string
!   integer :: system
!  end function system
! end interface
! ! end block
! !#endif
! !!SYSTEM_BLOCK_EDIT END COMMENTED_OUT



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
   ! program to compile may involve updating the fixsystem script and/or
   ! the mpi_utilities_mod.f90 for new compilers.
   call do_system("echo hello world", rc)
   if (rc /= 0) print *, 'call to system() returned error'

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
!> 'shell_name' is a namelist item and normally is the null string.
!> on at least on cray system, the compute nodes only had one type
!> of shell and you had to specify it.

subroutine do_system(execute, rc)

character(len=*), intent(in)  :: execute
integer,          intent(out) :: rc

! fill this in if an explicit shell name is required by your system
character(len=32) :: shell_name = ''

! !!NAG_BLOCK_EDIT START COMMENTED_OUT
!  call system(trim(shell_name)//' '//trim(execute)//' '//char(0), errno=rc)
! !!NAG_BLOCK_EDIT END COMMENTED_OUT
! !!OTHER_BLOCK_EDIT START COMMENTED_IN
    rc = system(trim(shell_name)//' '//trim(execute)//' '//char(0))
! !!OTHER_BLOCK_EDIT END COMMENTED_IN

end subroutine do_system



end program ftest_mpi

