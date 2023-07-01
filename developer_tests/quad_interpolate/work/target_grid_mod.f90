! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! read in an empty target grid from a netcdf file

module target_grid_mod

use types_mod, only : r8, metadatalength

! fix this to use only once we know what
! routines are needed
use grid_mod
use read_grid_mod
use netcdf_utilities_mod

implicit none

private

public :: create_gridT

contains


!---------------------------------------------
subroutine create_gridT(input_file, lon_name, lat_name, is_regular, gridT)

character(len=*), intent(in) :: input_file
character(len=*), intent(in) :: lon_name, lat_name
logical, intent(in) :: is_regular
type(grid_type), intent(out) :: gridT

integer :: i, ncid

ncid = nc_open_file_readonly(input_file, "create_gridT")

! if these need info from the file, pass in ncid 
call set_lon_lat_names(gridT, lon_name, lat_name)
call set_grid_type(gridT, is_regular)

call read_grid(ncid, gridT)

call nc_close_file(ncid)


end subroutine create_gridT

!---------------------------------------------
subroutine set_lon_lat_names(gridT, lon_name, lat_name)

type(grid_type), intent(inout) :: gridT
character(len=metadatalength) :: lon_name, lat_name

call set_grid_names(gridT, lon_name, lat_name)

end subroutine set_lon_lat_names

!---------------------------------------------
subroutine set_grid_type(gridT, is_regular)

type(grid_type), intent(inout) :: gridT
logical, intent(in) :: is_regular

if (is_regular) then
   call set_grid_type_regular(gridT)
else
   call set_grid_type_irregular(gridT)
endif

end subroutine set_grid_type

!---------------------------------------------
!---------------------------------------------

end module target_grid_mod
