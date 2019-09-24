! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: parse_args_test.f90 11289 2017-03-10 21:56:06Z hendric@ucar.edu $

! this test reads stdin.  to automate it, run it with:
!  cat bob | ./parse_args_test
! or
!  ./parse_args_test < bob
!
! where the file bob contains one line per test with arguments to be parsed.

program parse_args_test

use utilities_mod,     only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, &
                              E_ERR, E_MSG
use parse_args_mod,    only : get_args_from_string

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/preprocess/developer_tests/utilities/parse_args_test.f90 $"
character(len=*), parameter :: revision = "$Revision: 11289 $"
character(len=*), parameter :: revdate  = "$Date: 2017-03-10 14:56:06 -0700 (Fri, 10 Mar 2017) $"

integer :: iunit = 5
integer :: ierr, nargs, i
character(len=512) :: nextline
character(len=64) :: args(100)

!----------------------------------------------------------------------

! main code here
 
! initialize the dart libs
call initialize_utilities('parse_args_test')

print *, 'enter a line of text:'

GETMORE : do
   read(iunit, '(A)', iostat=ierr) nextline
   if (ierr /= 0) exit GETMORE

   call get_args_from_string(nextline, nargs, args)

   print *, 'got ', nargs, ' arguments from the input line'
   do i=1, nargs
      print *, 'arg ', i, 'is: ', '"'//trim(args(i))//'"'
   enddo

enddo GETMORE

! finalize parse_args_test
call error_handler(E_MSG,'parse_args_test','Finished successfully.',source,revision,revdate)
call finalize_utilities()

! end of main code

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/preprocess/developer_tests/utilities/parse_args_test.f90 $
! $Id: parse_args_test.f90 11289 2017-03-10 21:56:06Z hendric@ucar.edu $
! $Revision: 11289 $
! $Date: 2017-03-10 14:56:06 -0700 (Fri, 10 Mar 2017) $
