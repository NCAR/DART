! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! this test reads stdin.  to automate it, run it with:
!  cat bob | ./parse_words_test
! or
!  ./parse_words_test < bob
!
! where the file bob contains one line per test with words to be parsed.

program parse_words_test

use utilities_mod,     only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_ERR, E_MSG
use parse_args_mod,    only : get_csv_words_from_string

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "parse_words_test.f90"
character(len=*), parameter :: revision = ""
character(len=*), parameter :: revdate  = ""

integer :: iunit = 5
integer :: ierr, nwords, i
character(len=512) :: nextline
character(len=64) :: words(100)
character(len=64) :: names(100)
character(len=64) :: vals(100)
logical :: continue_line

!----------------------------------------------------------------------

! main code here
 
! initialize the dart libs
call initialize_utilities('parse_words_test')

print *, 'enter lines with words separated by commas.'
print *, '(enter control-D to end)'


GETMORE : do
   read(iunit, '(A)', iostat=ierr) nextline
   if (ierr /= 0) exit GETMORE
   if (nextline == 'NEXT') exit GETMORE

   call get_csv_words_from_string(nextline, ',', nwords, words)

   print *, 'input line: "'//trim(nextline)//'"'
   print *, 'got ', nwords, ' words from the input line'
   do i=1, nwords
      print *, 'word ', i, 'is: ', '"'//trim(words(i))//'"'
   enddo

enddo GETMORE

print *, 'enter lines with words separated by semicolons.'
print *, '(enter control-D to end)'


GETMORE2 : do
   read(iunit, '(A)', iostat=ierr) nextline
   if (ierr /= 0) exit GETMORE2
   if (nextline == 'NEXT') exit GETMORE2

   call get_csv_words_from_string(nextline, ';', nwords, words)

   print *, 'input line: "'//trim(nextline)//'"'
   print *, 'got ', nwords, ' words from the input line'
   do i=1, nwords
      print *, 'word ', i, 'is: ', '"'//trim(words(i))//'"'
   enddo

enddo GETMORE2

! finalize parse_words_test
call error_handler(E_MSG,'parse_words_test','Finished successfully.',source,revision,revdate)
call finalize_utilities()

! end of main code

!----------------------------------------------------------------------

end program

