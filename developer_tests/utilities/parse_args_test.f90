! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! this test reads stdin.  to automate it, run it with:
!  cat bob | ./parse_args_test
! or
!  ./parse_args_test < bob
!
! where the file bob contains one line per test with arguments to be parsed.

program parse_args_test

use utilities_mod,     only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_ERR, E_MSG
use parse_args_mod,    only : get_args_from_string, get_name_val_pairs_from_string

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "parse_args_test.f90"
character(len=*), parameter :: revision = ""
character(len=*), parameter :: revdate  = ""

integer :: iunit = 5
integer :: ierr, nargs, i
character(len=512) :: nextline
character(len=64) :: args(100)
character(len=64) :: names(100)
character(len=64) :: vals(100)
logical :: continue_line

!----------------------------------------------------------------------

! main code here
 
! initialize the dart libs
call initialize_utilities('parse_args_test')

print *, 'enter lines with words separated by blanks.'
print *, '(enter control-D to end)'


GETMORE : do
   read(iunit, '(A)', iostat=ierr) nextline
   if (ierr /= 0) exit GETMORE
   if (nextline == 'NEXT') exit GETMORE

   call get_args_from_string(nextline, nargs, args)

   print *, 'input line: "'//trim(nextline)//'"'
   print *, 'got ', nargs, ' arguments from the input line'
   do i=1, nargs
      print *, 'arg ', i, 'is: ', '"'//trim(args(i))//'"'
   enddo

enddo GETMORE

print *, 'enter lines with name=value separated by blanks.'
print *, 'enter control-D to end'

!rewind(iunit)
GETMORE2 : do
   read(iunit, '(A)', iostat=ierr) nextline
   if (ierr /= 0) exit GETMORE2
   if (nextline == 'NEXT') exit GETMORE2

   call get_name_val_pairs_from_string(nextline, nargs, names, vals, continue_line)

   print *, 'input line: "'//trim(nextline)//'"'
   print *, 'got ', nargs, ' arguments from the input line'
   do i=1, nargs
      print *, 'name  ', i, 'is: ', '"'//trim(names(i))//'"'
      print *, 'value ', i, 'is: ', '"'//trim(vals(i))//'"'
   enddo
   print *, 'continue line = ', continue_line

enddo GETMORE2

! finalize parse_args_test
call error_handler(E_MSG,'parse_args_test','Finished successfully.',source,revision,revdate)
call finalize_utilities()

! end of main code

!----------------------------------------------------------------------

end program

