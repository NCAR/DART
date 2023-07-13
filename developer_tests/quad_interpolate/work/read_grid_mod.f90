! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! read in an empty target grid from a netcdf file

module read_grid_mod

use types_mod, only : r8, metadatalength
use utilities_mod, only : error_handler, E_ERR

! fix these to use only once we know what
! routines are needed
use netcdf_utilities_mod 
use grid_mod

implicit none

private

character(len=*), parameter :: source = 'read_grid_mod'
character(len=256) :: msgstring

public :: read_grid

contains


!---------------------------------------------
subroutine read_grid(ncid, gridT)

integer, intent(in) :: ncid
type(grid_type), intent(inout) :: gridT

if (is_grid_type_regular(gridT)) then
   call read_reggrid(ncid, gridT)
else
   call read_irreggrid(ncid, gridT)
endif

end subroutine read_grid

!---------------------------------------------

!---------------------------------------------
subroutine read_reggrid(ncid, gridT)
 
integer :: ncid
type(grid_type), intent(inout) :: gridT

integer :: nlon, nlat
character(len=metadatalength) :: lon_name, lat_name

call get_grid_names(gridT, lon_name, lat_name)

call get_dim_size(ncid, lon_name, nlon)
call get_dim_size(ncid, lat_name, nlat)

call set_grid_sizes(gridT, nlon, nlat)

call allocate_grid_space(gridT)

! add an interface to get data arrays in grid mod

call read_1d_array(ncid, lon_name, gridT%irlon)
call read_1d_array(ncid, lat_name, gridT%irlat)

end subroutine read_reggrid

!---------------------------------------------
subroutine read_irreggrid(ncid, gridT)
 
integer :: ncid
type(grid_type), intent(inout) :: gridT

integer :: nlon, nlat
character(len=metadatalength) :: lon_name, lat_name

call get_grid_names(gridT, lon_name, lat_name)

call get_both_dim_sizes(ncid, lon_name, nlon, nlat)

call set_grid_sizes(gridT, nlon, nlat)

!print *, 'read irreggrid:'
!print *, 'lon, lat: ', trim(lon_name), ', ', trim(lat_name)
!print *, 'nlon, nlat: ', nlon, nlat

call allocate_grid_space(gridT)

! add an interface to get data arrays in grid mod

call read_2d_array(ncid, lon_name, gridT%iilon)
call read_2d_array(ncid, lat_name, gridT%iilat)

! call dump_grid(gridT)

end subroutine read_irreggrid

!---------------------------------------------

!---------------------------------------------
subroutine  get_dim_size(ncid, varname, dimsize)

integer :: ncid
character(len=*), intent(in) :: varname
integer, intent(out) :: dimsize

integer :: ndim

call nc_get_variable_num_dimensions(ncid, varname, ndim)
if (ndim /= 1) then
   write(msgstring, *) 'expected '//trim(varname)//' to be 1d, is ', ndim
   call error_handler(E_ERR, source, msgstring)
endif

call nc_get_variable_size(ncid, varname, dimsize)

end subroutine  get_dim_size

!---------------------------------------------
subroutine  get_both_dim_sizes(ncid, varname, dimsize1, dimsize2)

integer :: ncid
character(len=*), intent(in) :: varname
integer, intent(out) :: dimsize1, dimsize2

integer :: ndim
integer :: dimsizes(2)

call nc_get_variable_num_dimensions(ncid, varname, ndim)
if (ndim /= 2) then
   write(msgstring, *) 'expected '//trim(varname)//' to be 2d, is ', ndim
   call error_handler(E_ERR, source, msgstring)
endif

call nc_get_variable_size(ncid, varname, dimsizes)

dimsize1 = dimsizes(1)
dimsize2 = dimsizes(2)

end subroutine  get_both_dim_sizes

!---------------------------------------------
subroutine read_1d_array(ncid, varname, vardata)

integer,  intent(in)         :: ncid
character(len=*), intent(in) :: varname
real(r8), intent(out)        :: vardata(:)

call nc_get_variable(ncid, varname, vardata)

end subroutine read_1d_array

!---------------------------------------------
subroutine read_2d_array(ncid, varname, vardata)

integer,  intent(in)         :: ncid
character(len=*), intent(in) :: varname
real(r8), intent(out)        :: vardata(:, :)

call nc_get_variable(ncid, varname, vardata)

end subroutine read_2d_array

!---------------------------------------------

end module read_grid_mod
