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
use parse_args_mod,       only : csv_file_type, csv_get_obs_num, csv_get_field, &
                                 csv_open, csv_close, csv_print_header

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
!  Declare namelist parameters
logical            :: debug             = .true.


namelist /csv_read_test_nml/ debug

! start test
call initialize_utilities()

fname = "test1"
call create_data_1(fname)
call test_1(fname)
call delete_file(fname)

fname = "test2"
call create_data_2(fname)
call test_2(fname)
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

if (debug) rc = shell_execute("cat "//trim(fname))

end subroutine create_data_1

!-----------------------------------------------
! Open csv file and get data
subroutine test_1(fname)
character(len=*), intent(in) :: fname

integer :: i, nrows

character(len=*), parameter :: routine = "test_1"
character(len=256) :: nam(2)
integer :: tot(2)

! Open csv file and get dims
call csv_open(fname, cf, routine)
nrows = cf%nrows
if (debug) print *, "number of rows found = ", nrows

if (nrows /= 2) print *, "TEST 1 FAIL: BAD NUMBER OF ROWS"

if (debug) call csv_print_header(cf)

! Read the data
call csv_get_field(cf, 'Name', nam, routine)
if (nam(1) /= "Bob") print *, "TEST 1 FAIL:  READ BAD NAME 1 DATA"
if (nam(2) /= "Alice") print *, "TEST 1 FAIL:  READ BAD NAME 2 DATA"

call csv_get_field(cf, 'Total', tot, routine)
if (tot(1) /= 400) print *, "TEST 1 FAIL:  READ BAD TOTAL 1 DATA"
if (tot(2) /= 200) print *, "TEST 1 FAIL:  READ BAD TOTAL 2 DATA"

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

if (debug) rc = shell_execute("cat "//trim(fname))

end subroutine create_data_2

!-----------------------------------------------
! Open csv file and get data
subroutine test_2(fname)
character(len=*), intent(in) :: fname

integer :: i, nrows

character(len=*), parameter :: routine = "test_1"
character(len=256) :: nam(3)
integer :: tot(3)

! Open csv file and get dims
call csv_open(fname, cf, routine)
nrows = cf%nrows
if (debug) print *, "number of rows found = ", nrows

if (nrows /= 3) print *, "TEST 2 FAIL: BAD NUMBER OF ROWS"

if (debug) call csv_print_header(cf)

! Read the data
call csv_get_field(cf, 'Name', nam, routine)
if (nam(1) /= "Bob") print *, "TEST 2 FAIL:  READ BAD NAME 1 DATA"
if (nam(2) /= "Alice") print *, "TEST 2 FAIL:  READ BAD NAME 2 DATA"
if (nam(3) /= "Carl") print *, "TEST 2 FAIL:  READ BAD NAME 3 DATA"

call csv_get_field(cf, 'Total', tot, routine)
if (tot(1) /= 400) print *, "TEST 2 FAIL:  READ BAD TOTAL 1 DATA"
if (tot(2) /= 200) print *, "TEST 2 FAIL:  READ BAD TOTAL 2 DATA"
if (tot(3) /= 60) print *, "TEST 2 FAIL:  READ BAD TOTAL 3 DATA"

call csv_close(cf)

end subroutine test_2

!-----------------------------------------------
! clean up
subroutine delete_file(fname)
character(len=*), intent(in) :: fname

integer :: rc

rc = shell_execute("rm "//trim(fname))

end subroutine delete_file


end program csv_read_test
