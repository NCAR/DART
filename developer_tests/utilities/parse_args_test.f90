! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program parse_args_test

use utilities_mod,     only : initialize_utilities, finalize_utilities
use parse_args_mod,    only : get_args_from_string, get_name_val_pairs_from_string

use test  ! fortran-testanything

implicit none

integer, parameter :: MAXVALS = 100
integer :: nargs
character(len=512) :: nextline
character(len=64) :: args(MAXVALS)
character(len=64) :: names(MAXVALS)
character(len=64) :: vals(MAXVALS)
logical :: continue_line

character(len=20) :: testnum

! the ascii code for the backslash character is 92.
character(len=1), parameter :: BACKSLASH = ACHAR(92)

!----------------------------------------------------------------------

! main code here
 
! initialize the dart libs
call initialize_utilities('parse_args_test', standalone_program=.true.)

! this count should match the number of tests expected to be 
! executed if the test program runs to completion.
call plan(57)

! parse tests - divide the blank-separated words on a line
! into individual words.  note no interpretation of the
! contents is expected - the returns are always strings.

testnum = 'test 1'
call resetvals()
nextline = 'one two three'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 3),         trim(testnum)//', nargs')
call ok((args(1) == 'one'),   trim(testnum)//', arg 1')
call ok((args(2) == 'two'),   trim(testnum)//', arg 2')
call ok((args(3) == 'three'), trim(testnum)//', arg 3')
call ok((args(4) == ''),      trim(testnum)//', arg 4')


testnum = 'test 2'
call resetvals()
nextline = ' one      two three  '
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 3),         trim(testnum)//', nargs')
call ok((args(1) == 'one'),   trim(testnum)//', arg 1')
call ok((args(2) == 'two'),   trim(testnum)//', arg 2')
call ok((args(3) == 'three'), trim(testnum)//', arg 3')
call ok((args(4) == ''),      trim(testnum)//', arg 4')


testnum = 'test 3'
call resetvals()
nextline = '     one        two        three'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 3),         trim(testnum)//', nargs')
call ok((args(1) == 'one'),   trim(testnum)//', arg 1')
call ok((args(2) == 'two'),   trim(testnum)//', arg 2')
call ok((args(3) == 'three'), trim(testnum)//', arg 3')
call ok((args(4) == ''),      trim(testnum)//', arg 4')


testnum = 'test 4'
call resetvals()
nextline = 'one       two th/ree'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 3),          trim(testnum)//', nargs')
call ok((args(1) == 'one'),    trim(testnum)//', arg 1')
call ok((args(2) == 'two'),    trim(testnum)//', arg 2')
call ok((args(3) == 'th/ree'), trim(testnum)//', arg 3')
call ok((args(4) == ''),       trim(testnum)//', arg 4')


testnum = 'test 5'
call resetvals()
nextline = '"quoted string"'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 1),                 trim(testnum)//', nargs')
call ok((args(1) == 'quoted string'), trim(testnum)//', arg 1')
call ok((args(2) == ''),              trim(testnum)//', arg 2')


testnum = 'test 6'
call resetvals()
nextline = '  "quoted string"  "second string  with blanks"  '
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 2),                              trim(testnum)//', nargs')
call ok((args(1) == 'quoted string'),              trim(testnum)//', arg 1')
call ok((args(2) == 'second string  with blanks'), trim(testnum)//', arg 2')
call ok((args(3) == ''),                           trim(testnum)//', arg 3')


testnum = 'test 7'
call resetvals()
nextline = '"quoted" "string"'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 2),           trim(testnum)//', nargs')
call ok((args(1) == 'quoted'),  trim(testnum)//', arg 1')
call ok((args(2) == 'string'),  trim(testnum)//', arg 2')
call ok((args(3) == ''),        trim(testnum)//', arg 3')


testnum = 'test 8'
call resetvals()
nextline = 'string'//BACKSLASH//'abc'
call get_args_from_string(nextline, nargs, args)
 
call ok((nargs == 1),             trim(testnum)//', nargs')
call ok((args(1) == 'stringabc'), trim(testnum)//', arg 1')
call ok((args(2) == ''),          trim(testnum)//', arg 2')


testnum = 'test 9'
call resetvals()
nextline = 'directory'//BACKSLASH//BACKSLASH//'abc'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 1),                               trim(testnum)//', nargs')
call ok((args(1) == 'directory'//BACKSLASH//'abc'), trim(testnum)//', arg 1')
call ok((args(2) == ''),                            trim(testnum)//', arg 2')


! start of name/value pair tests

testnum = 'test 10'
call resetvals()
nextline = 'a=1 b=2 c=3'
call get_name_val_pairs_from_string(nextline, nargs, names, vals, continue_line)

call ok((nargs == 3),           trim(testnum)//', nargs')
call ok((names(1) == 'a'),      trim(testnum)//', name 1')
call ok((vals(1)  == '1'),      trim(testnum)//', val 1')
call ok((names(2) == 'b'),      trim(testnum)//', name 2')
call ok((vals(2)  == '2'),      trim(testnum)//', val 2')
call ok((names(3) == 'c'),      trim(testnum)//', name 3')
call ok((vals(3)  == '3'),      trim(testnum)//', val 3')
call ok((names(4) == ''),       trim(testnum)//', name 4')
call ok((vals(4)  == ''),       trim(testnum)//', val 4')
call ok((.not. continue_line),  trim(testnum)//', no continue_line')


testnum = 'test 11'
call resetvals()
nextline = 'a = 1.0 b   =  2.0  c=3.0  &'
call get_name_val_pairs_from_string(nextline, nargs, names, vals, continue_line)

call ok((nargs == 3),           trim(testnum)//', nargs')
call ok((names(1) == 'a'),      trim(testnum)//', name 1')
call ok((vals(1)  == '1.0'),    trim(testnum)//', val 1')
call ok((names(2) == 'b'),      trim(testnum)//', name 2')
call ok((vals(2)  == '2.0'),    trim(testnum)//', val 2')
call ok((names(3) == 'c'),      trim(testnum)//', name 3')
call ok((vals(3)  == '3.0'),    trim(testnum)//', val 3')
call ok((names(4) == ''),       trim(testnum)//', name 4')
call ok((vals(4)  == ''),       trim(testnum)//', val 4')
call ok(continue_line,          trim(testnum)//', has continue_line')


! finalize parse_args_test
call finalize_utilities()

! end of main code
!----------------------------------------------------------------------

contains

subroutine resetvals()

integer :: i

nargs = -888
do i=1, MAXVALS
   args(i) = ''
   names(i) = ''
   vals(i) = ''
enddo

end subroutine resetvals


end program

