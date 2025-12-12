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

!----------------------------------------------------------------------

! main code here
 
! initialize the dart libs
call initialize_utilities('parse_args_test')

call resetvals()
nextline = 'one two three'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 3),         'testline 1, nargs')
call ok((args(1) == 'one'),   'testline 1, arg 1')
call ok((args(2) == 'two'),   'testline 1, arg 2')
call ok((args(3) == 'three'), 'testline 1, arg 3')
call ok((args(4) == ''),      'testline 1, arg 4')


call resetvals()
nextline = '"quoted string"'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 1),                 'testline 2, nargs')
call ok((args(1) == 'quoted string'), 'testline 2, arg 1')
call ok((args(2) == ''),              'testline 2, arg 2')


call resetvals()
nextline = '1 20 300 4000'
call get_args_from_string(nextline, nargs, args)

call ok((nargs == 4),        'testline 3, nargs')
call ok((args(1) == '1'),    'testline 3, arg 1')
call ok((args(2) == '20'),   'testline 3, arg 2')
call ok((args(3) == '300'),  'testline 3, arg 3')
call ok((args(4) == '4000'), 'testline 3, arg 4')
call ok((args(5) == ''),     'testline 3, arg 5')


call resetvals()
nextline = 'a=1 b=2 c=3'
call get_name_val_pairs_from_string(nextline, nargs, names, vals, continue_line)

call ok((nargs == 3),           'testline 4, nargs')
call ok((names(1) == 'a'),      'testline 4, name 1')
call ok((vals(1)  == '1'),      'testline 4, val 1')
call ok((names(2) == 'b'),      'testline 4, name 2')
call ok((vals(2)  == '2'),      'testline 4, val 2')
call ok((names(3) == 'c'),      'testline 4, name 3')
call ok((vals(3)  == '3'),      'testline 4, val 3')
call ok((names(4) == ''),       'testline 4, name 4')
call ok((vals(4)  == ''),       'testline 4, val 4')
call ok((.not. continue_line),  'testline 4, continue_line')


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

