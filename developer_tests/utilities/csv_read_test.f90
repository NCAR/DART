! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! test cases for CSV read/parse routines.

program csv_read_test

use types_mod,            only : r8, MISSING_R8
use utilities_mod,        only : initialize_utilities, finalize_utilities, &
                                 file_exist, E_ERR, open_file, close_file
use mpi_utilities_mod,    only : shell_execute
use read_csv_mod,         only : csv_file_type, csv_get_field, &
                                 csv_open, csv_close, csv_get_nrows

use test    ! fortran-testanything

implicit none

character(len=*), parameter :: source = 'csv_read_test'

character(len=512) :: string1, string2

! File variables
character(len=256) :: fname
integer            :: io, iunit, filenum

character(len=512), allocatable :: dat(:)

! csv obs file 
type(csv_file_type) :: cf

!------------------------------------------------------------------------

! start test
call initialize_utilities(source, standalone_program=.true.)
call plan(24)

fname = "test_1"
call create_data_1(fname)
call test_1(fname)
call delete_file(fname)

fname = "test_2"
call create_data_2(fname)
call test_2(fname)
call delete_file(fname)

fname = "test_3"
call create_data_3(fname)
call test_3(fname)
call delete_file(fname)

fname = "test_4"
call create_data_4(fname)
call test_4(fname)
call delete_file(fname)

! end test
call finalize_utilities()


contains


!-----------------------------------------------
! Create test data file 
subroutine create_data_1(fname)
character(len=*), intent(in) :: fname

integer :: iunit, rc

iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name, Date, Total"
write(iunit, "(A)") "Bob, 1/1/25, 400"
write(iunit, "(A)") "Alice, 12/1/22, 200"
call close_file(iunit)

end subroutine create_data_1

!-----------------------------------------------
! Open csv file and get data
subroutine test_1(fname)
character(len=*), intent(in) :: fname

integer :: i, nrows

character(len=256) :: nam(2)
integer :: tot(2)

! Open csv file and get dims
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 2),         trim(fname)//", number of rows")

! Read the data
call csv_get_field(cf, 'Name', nam, fname)
call ok((nam(1) == 'Bob'),    trim(fname)//", read name 1")
call ok((nam(2) == 'Alice'),  trim(fname)//", read name 2")

call csv_get_field(cf, 'Total', tot, fname)
call ok((tot(1) == 400),      trim(fname)//", read total 1")
call ok((tot(2) == 200),      trim(fname)//", read total 2")

call csv_close(cf)

end subroutine test_1

!-----------------------------------------------
! Create test data file 
subroutine create_data_2(fname)
character(len=*), intent(in) :: fname

integer :: iunit, rc

iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name,,Date,Total"
write(iunit, "(A)") "Bob,,1/1/2025,400"
write(iunit, "(A)") "Alice,,12/1/2022,200"
write(iunit, "(A)") "Carl,,12/1/2020,60"
call close_file(iunit)

end subroutine create_data_2

!-----------------------------------------------
! Open csv file and get data
subroutine test_2(fname)
character(len=*), intent(in) :: fname

integer :: i, nrows

character(len=256) :: nam(3)
integer :: tot(3)

! Open csv file and get dims
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 3), trim(fname)//", number of rows")

! Read the data
call csv_get_field(cf, 'Name', nam, fname)
call ok((nam(1) == "Bob"),   trim(fname)//",  read name 1")
call ok((nam(2) == "Alice"), trim(fname)//",  read name 2")
call ok((nam(3) == "Carl"),  trim(fname)//",  read name 3")

call csv_get_field(cf, 'Total', tot, fname)
call ok((tot(1) == 400),     trim(fname)//",  read total 1")
call ok((tot(2) == 200),     trim(fname)//",  read total 2")
call ok((tot(3) == 60),      trim(fname)//",  read total 3")

call csv_close(cf)

end subroutine test_2

!-----------------------------------------------
! Create test data file 
subroutine create_data_3(fname)
character(len=*), intent(in) :: fname

integer :: iunit, rc

iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name;;Date;Total"
write(iunit, "(A)") "Bob;;1,1,2,0,2,5;40000"
write(iunit, "(A)") "Alice;;12.1.2022;200"
write(iunit, "(A)") "Carl;;12 4 2020;60"
call close_file(iunit)


end subroutine create_data_3

!-----------------------------------------------
! Open csv file and get data
subroutine test_3(fname)
character(len=*), intent(in) :: fname

integer :: i, nrows

character(len=256) :: nam(3)
integer :: tot(3)

! Open csv file and get dims
call csv_open(fname, cf, context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 3) , trim(fname)//", number of rows")


! Read the data
call csv_get_field(cf, 'Name', nam, fname)
call ok((nam(1) == "Bob"),   trim(fname)//", read name 1")
call ok((nam(2) == "Alice"), trim(fname)//", read name 2")
call ok((nam(3) == "Carl"),  trim(fname)//", read name 3")

call csv_get_field(cf, 'Total', tot, fname)
call ok((tot(1) == 40000),   trim(fname)//", read total 1")
call ok((tot(2) == 200),     trim(fname)//", read total 2")
call ok((tot(3) == 60),      trim(fname)//", read total 3")

call csv_close(cf)

end subroutine test_3

!-----------------------------------------------
! Create test data file 
subroutine create_data_4(fname)
character(len=*), intent(in) :: fname

integer :: iunit, rc

iunit = open_file(fname, action="write")
write(iunit, "(A)") "Name: Date: Total"
write(iunit, "(A)") "Bob: 1/1/25: 400"
write(iunit, "(A)") "Alice: 12/1/22: 200"
call close_file(iunit)

end subroutine create_data_4

!-----------------------------------------------
! Open csv file and get data
subroutine test_4(fname)
character(len=*), intent(in) :: fname

integer :: i, nrows

character(len=256) :: nam(2)
integer :: tot(2)

! Open csv file and get dims
call csv_open(fname, cf, forced_delim=':', context=fname)

nrows = csv_get_nrows(cf)
call ok((nrows == 2),         trim(fname)//", number of rows")

! Read the data
call csv_get_field(cf, 'Name', nam, fname)
call ok((nam(1) == 'Bob'),    trim(fname)//", read name 1")
call ok((nam(2) == 'Alice'),  trim(fname)//", read name 2")

call csv_get_field(cf, 'Total', tot, fname)
call ok((tot(1) == 400),      trim(fname)//", read total 1")
call ok((tot(2) == 200),      trim(fname)//", read total 2")

call csv_close(cf)

end subroutine test_4

!-----------------------------------------------
! clean up
subroutine delete_file(fname)
character(len=*), intent(in) :: fname

integer :: rc

rc = shell_execute("rm "//trim(fname))

end subroutine delete_file


end program csv_read_test
