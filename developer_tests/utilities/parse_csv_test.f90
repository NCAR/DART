! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program parse_csv_test

use utilities_mod,        only : initialize_utilities, finalize_utilities, &
                                 open_file, close_file
use parse_args_mod,       only : get_csv_words_from_string

use test   ! fortran-testanything

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "parse_csv_test.f90"

integer :: nwords, i
character(len=512) :: nextline
character(len=64) :: words(100)

!----------------------------------------------------------------------

! main code here
 
! initialize the dart libs
call initialize_utilities('parse_csv_test', standalone_program=.true.)
call plan(19)

nextline = "these,are,comma,separated,words,on,a,line"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 8), "line 1, word count")
call ok((words(1) == "these"), "line 1, word 1")
call ok((words(8) == "line"), "line 1, word 8")

nextline = "one,two,three"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 3), "line 2, word count")
call ok((words(1) == "one"), "line 2, word 1")
call ok((words(3) == "three"), "line 2, word 3")

nextline = '"quoted string1"'//",nonquote,'quote two','embedded_"//'"'//"_quote'"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 4), "line 3, word count")
call ok((words(1) == "quoted string1"), "line 3, word 1")
call ok((words(2) == "nonquote"), "line 3, word 2")
call ok((words(3) == "quote two"), "line 3, word 3")
call ok((words(4) == 'embedded_"_quote'), "line 3, word 4")

nextline = "1;2;3"
call get_csv_words_from_string(nextline, ';', nwords, words)
call ok((nwords == 3), "line 4, word count")
call ok((words(1) == "1"), "line 4, word 1")
call ok((words(3) == "3"), "line 4, word 3")

nextline = "single_string"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 1), "line 5, word count")
call ok((words(1) == "single_string"), "line 5, word 1")

nextline = "comma separated strings,with embedded blanks"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 2), "line 6, word count")
call ok((words(1) == "comma separated strings"), "line 6, word 1")
call ok((words(2) == "with embedded blanks"), "line 6, word 2")

call finalize_utilities()

! end of main code

!----------------------------------------------------------------------

end program

