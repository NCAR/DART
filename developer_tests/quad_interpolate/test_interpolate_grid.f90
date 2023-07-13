! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


program test_interpolate_grid

use regular_grid_mod, only : create_grid0, create_field, set_data_function
use types_mod, only : r8, MISSING_R8, metadatalength
use functions_mod, only : f_sine, f_sum, f_row, f_col, f_rand
use grid_mod, only : grid_type, dump_grid
use target_grid_mod, only : create_gridT
use utilities_mod, only : initialize_utilities, finalize_utilities, &
                          array_dump, register_module, do_nml_file, do_nml_term, &
                          nmlfileunit, find_namelist_in_file, &
                          check_namelist_read
use write_grid_mod, only : write_grid

use quad_interp_mod, only : do_interp, set_quad_grid_opts


implicit none

! irregular grids have lon(:,:) lat(:,:)
! so the split between reg and irreg has to
! be somewhere.  do we need separate progs
! for interp to reg and interp to irreg ?


! maybe grid0 should also be a grid type?

real(r8), allocatable :: lon(:)
real(r8), allocatable :: lat(:)
real(r8), allocatable :: field0(:,:), fieldT(:,:), field1(:,:), field2(:,:)

real(r8) :: resolution = 1
real(r8) :: grid_origin_longitude = 0.0_r8
real(r8) :: grid_end_longitude = 360.0_r8  ! minus resolution if global
real(r8) :: grid_origin_latitude = -90.0_r8
real(r8) :: grid_end_latitude = 90.0_r8


real(r8) :: grid_corners(2,2)   ! (min, max lon), (min, max lat)
real(r8) :: globe_corners(2,2) 
data globe_corners(:,1)  /   0.0_r8, 360.0_r8 /
data globe_corners(:,2)  / -90.0_r8,  90.0_r8 /
integer  :: ni, nj
integer  :: i,j ! loop variables

integer :: case = 1
logical :: is_regular = .true.
integer :: debug = 0

logical :: grid_global = .true.
logical :: grid_spans_lon_zero = .true.
logical :: grid_pole_wrap = .true.

character(len=256) :: data_function = "sine"
character(len=256) :: target_filename
character(len=metadatalength) :: lon_name, lat_name
type(grid_type) :: grid0, gridT, grid1, grid2

namelist /test_interpolate_grid_nml/ is_regular, target_filename, &
         case, debug, resolution, lon_name, lat_name, resolution, &
         grid_global, grid_spans_lon_zero, grid_pole_wrap, data_function, &
         grid_origin_longitude, grid_origin_latitude, grid_end_longitude, &
         grid_end_latitude


call initialize_utilities("test_interpolate_grid")
call read_namelist()
call set_quad_grid_opts(grid_global, grid_spans_lon_zero, grid_pole_wrap)
call set_data_function(data_function)


! make source grid and give it data

call set_grid_corners(grid_origin_longitude, grid_end_longitude, grid_origin_latitude, grid_end_latitude)
call create_grid0(resolution, grid_corners, grid0, ni, nj)

allocate(field0(ni,nj))
call create_field(grid0, field0)
call write_grid(grid0, field0, "field0.nc")

! make target grid, no data

call create_gridT(target_filename, lon_name, lat_name, is_regular, gridT)
allocate(fieldT(gridT%nlon, gridT%nlat))


! call quad utils to move data from src to dst

call do_interp(grid0, field0, gridT, fieldT)
call write_grid(gridT, fieldT, "fieldT.nc")


! create a copy of the src grid with no data

call create_grid0(resolution, grid_corners, grid1, ni, nj)
allocate(field1(ni,nj))
field1(:,:) = MISSING_R8

! call quad utils to move data back from dst to src2

call do_interp(gridT, fieldT, grid1, field1)
call write_grid(grid1, field1, "field1.nc")


! compare field0 and field1 and evaluate success or failure
! call evaluate_results(grid0, field0, grid1, field1)

! create a denser copy of the src grid with no data

call create_grid0(resolution/2.0, grid_corners, grid2, ni, nj)
allocate(field2(ni,nj))
field2(:,:) = MISSING_R8

! call quad utils to move data back from dst to src2

call do_interp(gridT, fieldT, grid2, field2)
call write_grid(grid2, field2, "field2.nc")


! compare field0 and field1 and evaluate success or failure
! call evaluate_results(grid0, field0, grid1, field1)

deallocate(field0, fieldT, field1, field2)

call finalize_utilities()


contains

!------------------------------------------------------------------
subroutine read_namelist()          
                                 
integer :: iunit, io

! Read the namelist 
call find_namelist_in_file("input.nml", "test_interpolate_grid_nml", iunit)
read(iunit, nml = test_interpolate_grid_nml, iostat = io)
call check_namelist_read(iunit, io, "test_interpolate_grid_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=test_interpolate_grid_nml)
if (do_nml_term()) write(     *     , nml=test_interpolate_grid_nml)

end subroutine read_namelist


!------------------------------------------------------------------
subroutine set_grid_corners(xmin, xmax, ymin, ymax)

real(r8), intent(in) :: xmin, xmax, ymin, ymax

grid_corners(1,1) = xmin
grid_corners(2,1) = xmax
grid_corners(1,2) = ymin
grid_corners(2,2) = ymax

end subroutine set_grid_corners

end program test_interpolate_grid
