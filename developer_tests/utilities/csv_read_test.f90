! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! test cases for CSV read/parse routines.

program csv_read_test

use types_mod,            only : r8, MISSING_R8
use utilities_mod,        only : initialize_utilities, finalize_utilities, &
                                 file_exist, E_ERR, open_file, close_file
use parse_args_mod,       only : csv_file_type, csv_get_obs_num, csv_get_field,       &
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

! call test 1
fname = "test1"
call create_data_1(fname)
call test_1(fname)
call delete_data_1(fname)

! end test
call finalize_utilities()


contains


!-----------------------------------------------
! Create test data file 
subroutine create_data_1(fname)
character(len=*), intent(in) :: fname

integer :: iunit, nrows
character(len=*), parameter :: routine = "create_data_1"

iunit = open_file(fname, action="write")
write(iunit, *) "Name, Date, Total"
write(iunit, *) "Bob, 1/1/25, 400"
write(iunit, *) "Alice, 12/1/22, 200"
call close_file(iunit)

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

if (nrows /= 2) print *, "TEST 1 FAILED"

! Read the data
call csv_get_field(cf, 'Name', nam, routine)
call csv_get_field(cf, 'Total', tot, routine)

call csv_close(cf)

end subroutine test_1

!-----------------------------------------------
! clean up
subroutine delete_data_1(fname)
character(len=*), intent(in) :: fname

integer :: iunit, rc
character(len=*), parameter :: routine = "delete_data_1"

character(len=256) :: filedel

write(filedel, *) "echo rm ", fname
call system(filedel, rc)

end subroutine delete_data_1


end program csv_read_test
