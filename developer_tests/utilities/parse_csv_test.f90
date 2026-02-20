! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! test parsing comma separated values (CSV) from a single input string.
! tests include missing columns, embedded blanks, quoted values,
! and alternative separators.

program parse_csv_test

use utilities_mod,        only : initialize_utilities, finalize_utilities, &
                                 open_file, close_file
use read_csv_mod,         only : get_csv_words_from_string

use test   ! fortran-testanything

implicit none

integer :: nwords, i
character(len=512) :: nextline
character(len=64) :: words(100)
character(len=64) :: testname

!----------------------------------------------------------------------

! main code here
 
! initialize the dart libs
call initialize_utilities('parse_csv_test', standalone_program=.true.)
call plan(20)

testname = "basic case"
nextline = "these,are,comma,separated,words,on,a,line"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 8),         trim(testname)//", word count")
call ok((words(1) == "these"), trim(testname)//", word 1")
call ok((words(8) == "line"),  trim(testname)//", word 8")

testname = "empty column"
nextline = "one,,three,four"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 4),         trim(testname)//", word count")
call ok((words(1) == "one"),   trim(testname)//", word 1")
call ok((words(2) == ""),      trim(testname)//", word 2")
call ok((words(3) == "three"), trim(testname)//", word 3")

testname = "quotes"
nextline = '"quoted string1"'//",nonquote,'quote two','embedded_"//'"'//"_quote'"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 4),                    trim(testname)//", word count")
call ok((words(1) == "quoted string1"),   trim(testname)//", word 1")
call ok((words(2) == "nonquote"),         trim(testname)//", word 2")
call ok((words(3) == "quote two"),        trim(testname)//", word 3")
call ok((words(4) == 'embedded_"_quote'), trim(testname)//", word 4")

testname = "different separator"
nextline = "1;2;3"
call get_csv_words_from_string(nextline, ';', nwords, words)
call ok((nwords == 3),     trim(testname)//", word count")
call ok((words(1) == "1"), trim(testname)//", word 1")
call ok((words(3) == "3"), trim(testname)//", word 3")

testname = "single word"
nextline = "single_string"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 1),                 trim(testname)//", word count")
call ok((words(1) == "single_string"), trim(testname)//", word 1")

testname = "embedded blanks"
nextline = "comma separated strings,with embedded blanks"
call get_csv_words_from_string(nextline, ',', nwords, words)
call ok((nwords == 2),                           trim(testname)//", word count")
call ok((words(1) == "comma separated strings"), trim(testname)//", word 1")
call ok((words(2) == "with embedded blanks"),    trim(testname)//", word 2")

call finalize_utilities()

! end of main code

!----------------------------------------------------------------------

end program

