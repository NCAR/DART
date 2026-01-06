! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! test cases for CSV read routines.
!
! creates test files, runs test, and deletes file.
!
 

program csv_read_test

use types_mod,            only : r8, MISSING_R8, MISSING_I
use utilities_mod,        only : initialize_utilities, finalize_utilities, &
                                 open_file, close_file
use mpi_utilities_mod,    only : shell_execute
use read_csv_mod,         only : csv_file_type, csv_get_field, &
                                 csv_open, csv_close, csv_get_nrows

use test    ! fortran-testanything

implicit none


character(len=512) :: string1, string2
character(len=64)  :: testname

! File variables
integer            :: io, iunit

! if you want to keep the generated input files for further testing
! set this to true and the file deletes will not happen.
logical :: keep_testfiles = .false.

character(len=512), allocatable :: dat(:)
character, parameter :: BACKSLASH = ACHAR(92)

! csv file handle
type(csv_file_type) :: cf

!------------------------------------------------------------------------

! start test
call initialize_utilities('csv_read_test', standalone_program=.true.)
call plan(44)   ! should be +1 with bad_column() and +5 with short_array() tests

testname = "basic case"
call basic_test("test1", testname)

testname = "empty fields"
call empty_test("test2", testname)

testname = "semicolon separated"
call diff_separator("test3", testname)

! this one provokes a (correct) fatal error
!testname = "bad column"
!call bad_column("test4", testname)

testname = "blank test"
call blank_test("test5", testname)

! this one provokes a (correct) fatal error
!testname = "short arrays"
!call short_array("test6", testname)

testname = "bad num columns"
call bad_cols("test7", testname)

testname = "autodetect separators"
call sep_test("test8", testname)

testname = "type mismatch"
call bad_type("test9", testname)

testname = "mixed numerics"
call mixed("test10", testname)


! end test
call finalize_utilities()


contains


!-----------------------------------------------
! Basic test:
! Create a test file and close it.
! Open csv file and get data.
! 3 columns of data, nothing strange.

subroutine basic_test(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(2)
integer :: tot(2)

! create test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,Date,Total"
write(iunit, "(A)") "Alice,12/1/22,200"
write(iunit, "(A)") "Bob,1/1/25,400"
call close_file(iunit)


! start test
! assume comma as separator
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 2),         trim(testname)//", number of rows")

! Read the data
! the return data type is set by the type of the third argument.
! this one is a string.
call csv_get_field(cf, 'Name', nam, testname)
call ok((nam(1) == 'Alice'),  trim(testname)//", read name 1")
call ok((nam(2) == 'Bob'),    trim(testname)//", read name 2")

! this one is an integer
call csv_get_field(cf, 'Total', tot, testname)
call ok((tot(1) == 200),      trim(testname)//", read total 1")
call ok((tot(2) == 400),      trim(testname)//", read total 2")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine basic_test

!-----------------------------------------------
! test for blank column, blanks in column headers

subroutine empty_test(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(4)
integer :: tot(4)


! create the test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,,Date,Total Days"
write(iunit, "(A)") "Alice,,12/1/2022,200"
write(iunit, "(A)") "Bob,,1/1/2025,400"
write(iunit, "(A)") "Carl,,12/1/2020,60"
write(iunit, "(A)") "David,,4/1/2021,160"
call close_file(iunit)


! run the test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 4), trim(testname)//", number of rows")

! Read the data
call csv_get_field(cf, 'Name', nam, fname)
call ok((nam(2) == "Bob"),   trim(testname)//",  read name 2")
call ok((nam(4) == "David"), trim(testname)//",  read name 4")

call csv_get_field(cf, 'Total Days', tot, fname)
call ok((tot(1) == 200),     trim(testname)//",  read total 1")
call ok((tot(4) == 160),     trim(testname)//",  read total 4")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine empty_test

!-----------------------------------------------
! test different separator character
! and reals

subroutine diff_separator(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(3), idents(3)
real(r8) :: tot(3)

! create test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name;;Id;Total"
write(iunit, "(A)") "Alice;;12.1.2022;200.2"
write(iunit, "(A)") "Bob;;1,1,2,0,2,5;4.0e3"
write(iunit, "(A)") "Carl;;12 4 2020;60"
call close_file(iunit)


! run test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 3) , trim(testname)//", number of rows")

! Read the data
call csv_get_field(cf, 'Name', nam, fname)
call ok((nam(1) == "Alice"), trim(testname)//", read name 1")
call ok((nam(2) == "Bob"),   trim(testname)//", read name 2")
call ok((nam(3) == "Carl"),  trim(testname)//", read name 3")

call csv_get_field(cf, 'Total', tot, fname)
call ok((tot(1) == 200.2_r8),   trim(testname)//", read total 1")
call ok((tot(2) == 4000_r8),    trim(testname)//", read total 2")
call ok((tot(3) == 60._r8),     trim(testname)//", read total 3")

call csv_get_field(cf, 'Id', idents, fname)
call ok((idents(1) == '12.1.2022'),   trim(testname)//", read ident 1")
call ok((idents(2) == '1,1,2,0,2,5'), trim(testname)//", read ident 2")
call ok((idents(3) == '12 4 2020'),   trim(testname)//", read ident 3")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine diff_separator

!-----------------------------------------------
! try to retrieve a bad column  name

subroutine bad_column(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(2)
integer :: tot(2)

! create test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,Date,Total"
write(iunit, "(A)") "Alice,12/1/22,200"
write(iunit, "(A)") "Bob,1/1/25,400"
call close_file(iunit)


! start test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 2),         trim(testname)//", number of rows")

! try to read a non-existant column
call csv_get_field(cf, 'Age', nam, testname)

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine bad_column

!-----------------------------------------------
! different separator and embedded blanks?

subroutine blank_test(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit

character(len=256) :: nam(2)
real(r8) :: tot(2)


! create file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name: Date: Total"
write(iunit, "(A)") "Bob: 1/1/25: 400.04"
write(iunit, "(A)") "Alice: 12/1/22: 202.8"
call close_file(iunit)


! run test
call csv_open(fname, cf, forced_delim=':', context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 2),         trim(testname)//", number of rows")

! Read the data
call csv_get_field(cf, 'Name', nam, fname)
call ok((nam(1) == 'Bob'),    trim(testname)//", read name 1")
call ok((nam(2) == 'Alice'),  trim(testname)//", read name 2")

call csv_get_field(cf, 'Total', tot, fname)
call ok((tot(1) == 400.04_r8),   trim(testname)//", read total 1")
call ok((tot(2) == 202.8_r8),    trim(testname)//", read total 2")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine blank_test


!-----------------------------------------------
! test too short an array as input

subroutine short_array(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(2)
integer :: tot(2)

! create test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,Total"
write(iunit, "(A)") "Alice,200"
write(iunit, "(A)") "Bob,400"
write(iunit, "(A)") "Carl,60"
write(iunit, "(A)") "David,160"
call close_file(iunit)


! start test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 4),         trim(testname)//", number of rows")

! try to read the data into an array that is too short.
! string
call csv_get_field(cf, 'Name', nam, testname)
call ok((nam(1) == 'Bob'),    trim(testname)//", read name 1")
call ok((nam(2) == 'Alice'),  trim(testname)//", read name 2")

! integer
call csv_get_field(cf, 'Total', tot, testname)
call ok((tot(1) == 400),      trim(testname)//", read total 1")
call ok((tot(2) == 200),      trim(testname)//", read total 2")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine short_array


!-----------------------------------------------
! test badly formed input file - wrong number of columns
! in some of the lines

subroutine bad_cols(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(4)
integer :: tot(4)

! create test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,Weekly,Total"
write(iunit, "(A)") "Alice,200,1000"  ! right
write(iunit, "(A)") "Bob,400"         ! no third number
write(iunit, "(A)") "60,Carl,,60"     ! too many numbers
write(iunit, "(A)") "David"           ! single value
call close_file(iunit)


! start test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 4),         trim(testname)//", number of rows")

! try to read the data into an array that is too short.
! string
call csv_get_field(cf, 'Name', nam, testname)
call ok((nam(1) == 'Alice'),  trim(testname)//", read name 1")
call ok((nam(2) == 'Bob'),    trim(testname)//", read name 2")

! integer
call csv_get_field(cf, 'Total', tot, testname)
call ok((tot(1) == 1000),        trim(testname)//", read total 1")
call ok((tot(2) == MISSING_R8),  trim(testname)//", read total 2")
call ok((tot(3) == MISSING_R8),  trim(testname)//", read total 3")
call ok((tot(4) == MISSING_R8),  trim(testname)//", read total 4")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine bad_cols


!-----------------------------------------------
! test wrong type passed in for a column

subroutine bad_type(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(2)
integer :: tot(2)

! create test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,Total"
write(iunit, "(A)") "Alice,abc"
write(iunit, "(A)") "Bob,#$%"
call close_file(iunit)


! start test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 2),    trim(testname)//", number of rows")

! try to read a non-numeric string into an integer array
! integer
call csv_get_field(cf, 'Total', tot, testname)
call ok((tot(1) == MISSING_I),  trim(testname)//", read total 1")
call ok((tot(2) == MISSING_I),  trim(testname)//", read total 2")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine bad_type


!-----------------------------------------------
! test automatic separator test

subroutine sep_test(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(2)
integer :: tot(2)

! create test file
! this header has 1 more comma than semicolon, so should
! detect comma is the separator.
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name;X,Day;,Week;,Month,Total"
write(iunit, "(A)") "Alice;T,1,2,3,6"
write(iunit, "(A)") "Bob;R,2,4,6,12"
call close_file(iunit)


! start test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 2),    trim(testname)//", number of rows")

! make sure columns are being counted correctly
call csv_get_field(cf, 'Total', tot, testname)
call ok((tot(1) == 6),  trim(testname)//", read total 1")
call ok((tot(2) == 12), trim(testname)//", read total 2")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine sep_test


!-----------------------------------------------
! test column with mixed types

subroutine mixed(fname, testname)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: testname

integer :: i, nrows
integer :: iunit
character(len=256) :: nam(5)
real(r8) :: tot(5)

! create test file
iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,Total"
write(iunit, "(A)") "Alice,123.45"
write(iunit, "(A)") "Bob,abc"
write(iunit, "(A)") "Charles,333.33"
write(iunit, "(A)") "David,x1r2t3"
write(iunit, "(A)") "Edgar,123"
call close_file(iunit)


! start test
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 5),    trim(testname)//", number of rows")

! try to read a mix of numbers and non-numeric strings into a real array

call csv_get_field(cf, 'Total', tot, testname)
call ok((tot(1) == 123.45_r8),  trim(testname)//", read total 1")
call ok((tot(2) == MISSING_R8), trim(testname)//", read total 2")
call ok((tot(3) == 333.33_r8),  trim(testname)//", read total 3")
call ok((tot(4) == MISSING_R8), trim(testname)//", read total 4")
call ok((tot(5) == 123.0_r8),   trim(testname)//", read total 5")

call csv_close(cf)


! delete the test file
call delete_file(fname)

end subroutine mixed

!-----------------------------------------------
! clean up

subroutine delete_file(fname)
character(len=*), intent(in) :: fname

integer :: rc

! do not delete test files if true
if (keep_testfiles) return

rc = shell_execute(BACKSLASH//"rm "//trim(fname))

end subroutine delete_file


end program csv_read_test
