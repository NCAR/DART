! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! write a grid and data to a netcdf file for plotting
! and diagnosis

module write_grid_mod

use types_mod, only : r8, metadatalength, MISSING_R8

! fix this to use only once we know what
! routines are needed
use grid_mod

use  netcdf_utilities_mod, only :  nc_create_file, nc_close_file, &
                                   nc_define_dimension, nc_define_double_variable, &
                                   nc_end_define_mode, nc_add_global_creation_time, &
                                   nc_add_attribute_to_variable, nc_put_variable

implicit none

private

interface get_range
 module procedure get_range_1d
 module procedure get_range_2d
end interface

public :: write_grid

contains


!---------------------------------------------
subroutine write_grid(grid, field, output_file)

type(grid_type), intent(in) :: grid
real(r8), intent(in) :: field(:,:)
character(len=*), intent(in) :: output_file

logical :: is_regular

is_regular = is_grid_type_regular(grid)

if (is_regular) then
   call write_reg_grid(grid, field, output_file)
else
   call write_irreg_grid(grid, field, output_file)
endif

end subroutine write_grid

!---------------------------------------------
subroutine write_reg_grid(grid, field, output_file)

type(grid_type), intent(in) :: grid
real(r8), intent(in) :: field(:,:)
character(len=*), intent(in) :: output_file

integer :: i, ncid
real(r8) :: data_range(2)
integer :: data_size(2)
character(len=*), parameter :: routine = "write_reg_grid"


ncid = nc_create_file(output_file, routine)

call nc_add_global_creation_time(ncid, routine)

call nc_define_dimension(ncid, "lon",  grid%nlon, routine)
call nc_define_dimension(ncid, "lat",  grid%nlat, routine)

call get_range(grid%irlon, data_range)
call nc_define_double_variable(ncid, "lon", "lon", routine)
call nc_add_attribute_to_variable(ncid, "lon", "range", data_range, routine)
call nc_add_attribute_to_variable(ncid, "lon", "cartesian_axis", "X", routine)

call get_range(grid%irlat, data_range)
call nc_define_double_variable(ncid, "lat", "lat", routine)
call nc_add_attribute_to_variable(ncid, "lat", "range", data_range, routine)
call nc_add_attribute_to_variable(ncid, "lat", "cartesian_axis", "Y", routine)

call get_range(field, data_range)
call nc_define_double_variable(ncid, "field", (/ "lon", "lat" /) , routine)
call nc_add_attribute_to_variable(ncid, "field", "long_name", "data", routine)
call nc_add_attribute_to_variable(ncid, "field", "range", data_range, routine)
call nc_add_attribute_to_variable(ncid, "field", "_FillValue", MISSING_R8, routine)
call nc_add_attribute_to_variable(ncid, "field", "missing_value", MISSING_R8, routine)

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)


! Fill the dimension variables, then the data field
call nc_put_variable(ncid, "lon", grid%irlon, routine)
call nc_put_variable(ncid, "lat", grid%irlat, routine)

call nc_put_variable(ncid, "field", field, routine)

! close and return
call nc_close_file(ncid, routine)

end subroutine write_reg_grid

!---------------------------------------------
subroutine write_irreg_grid(grid, field, output_file)

type(grid_type), intent(in) :: grid
real(r8), intent(in) :: field(:,:)
character(len=*), intent(in) :: output_file

integer :: i, ncid
real(r8) :: data_range(2)
integer :: data_size(2)
character(len=*), parameter :: routine = "write_irreg_grid"


ncid = nc_create_file(output_file, routine)

call nc_add_global_creation_time(ncid, routine)

call nc_define_dimension(ncid, "i", grid%nlon, routine)
call nc_define_dimension(ncid, "j", grid%nlat, routine)

call get_range(grid%iilon, data_range)
call nc_define_double_variable(ncid, "lon", (/ "i", "j" /), routine)
call nc_add_attribute_to_variable(ncid, "lon", "range", data_range, routine)

call get_range(grid%iilat, data_range)
call nc_define_double_variable(ncid, "lat", (/ "i", "j" /), routine)
call nc_add_attribute_to_variable(ncid, "lat", "range", data_range, routine)

call get_range(field, data_range)
call nc_define_double_variable(ncid, "field", (/ "i", "j" /), routine)
call nc_add_attribute_to_variable(ncid, "field", "long_name", "data", routine)
call nc_add_attribute_to_variable(ncid, "field", "range", data_range, routine)
call nc_add_attribute_to_variable(ncid, "field", "_FillValue", MISSING_R8, routine)
call nc_add_attribute_to_variable(ncid, "field", "missing_value", MISSING_R8, routine)

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)


! Fill the dimension variables, then the data field
call nc_put_variable(ncid, "lon", grid%iilon, routine)
call nc_put_variable(ncid, "lat", grid%iilat, routine)

call nc_put_variable(ncid, "field", field, routine)

! close and return
call nc_close_file(ncid, routine)

end subroutine write_irreg_grid

!---------------------------------------------
subroutine get_range_1d(data_array, data_range)

real(r8), intent(in) :: data_array(:)
real(r8), intent(out) :: data_range(2)

data_range(1) = minval(data_array)
data_range(2) = maxval(data_array)

end subroutine get_range_1d

!---------------------------------------------
subroutine get_range_2d(data_array, data_range)

real(r8), intent(in) :: data_array(:,:)
real(r8), intent(out) :: data_range(2)

data_range(1) = minval(data_array)
data_range(2) = maxval(data_array)

end subroutine get_range_2d
!---------------------------------------------

!---------------------------------------------
!---------------------------------------------

end module write_grid_mod
