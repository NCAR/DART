! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! define a grid

module grid_mod

use types_mod, only : r8, metadatalength

implicit none

private

integer, parameter :: GRID_UNKNOWN = 0
integer, parameter :: GRID_REGULAR = 1
integer, parameter :: GRID_IRREGULAR = 2

type grid_type
 integer :: type = GRID_UNKNOWN
 integer :: nlon = 0
 integer :: nlat = 0
 character(len=metadatalength) :: lon_name = 'unknown'
 character(len=metadatalength) :: lat_name = 'unknown'
 real(r8), allocatable :: rrlon(:),   rrlat(:)
 real(r8), allocatable :: irlon(:,:), irlat(:,:)
end type

public :: grid_type, set_grid_type_regular, &
          is_grid_type_regular, &
          set_grid_type_irregular, &
          is_grid_type_irregular, &
          set_grid_sizes, get_grid_sizes, &
          set_grid_names, get_grid_names, &
          allocate_grid_space

contains

!---------------------------------------------
subroutine set_grid_type_regular(grid)

type(grid_type), intent(inout) :: grid

grid%type = GRID_REGULAR

end subroutine

!---------------------------------------------
subroutine set_grid_type_irregular(grid)

type(grid_type), intent(inout) :: grid

grid%type = GRID_IRREGULAR

end subroutine

!---------------------------------------------
function is_grid_type_regular(grid)

type(grid_type), intent(in) :: grid
logical :: is_grid_type_regular


is_grid_type_regular = (grid%type == GRID_REGULAR)

end function

!---------------------------------------------
function is_grid_type_irregular(grid)

type(grid_type), intent(in) :: grid
logical :: is_grid_type_irregular

is_grid_type_irregular = (grid%type == GRID_IRREGULAR)

end function

!---------------------------------------------
subroutine set_grid_sizes(grid, nlon, nlat)

type(grid_type), intent(inout) :: grid
integer, intent(in) :: nlon, nlat

grid%nlon = nlon
grid%nlat = nlat

end subroutine

!---------------------------------------------
subroutine get_grid_sizes(grid, nlon, nlat)

type(grid_type), intent(inout) :: grid
integer, intent(out) :: nlon, nlat

nlon = grid%nlon
nlat = grid%nlat

end subroutine

!---------------------------------------------
subroutine set_grid_names(grid, lon_name, lat_name)

type(grid_type), intent(inout) :: grid
character(len=*), intent(in) :: lon_name, lat_name

grid%lon_name = lon_name
grid%lat_name = lat_name

end subroutine

!---------------------------------------------
subroutine get_grid_names(grid, lon_name, lat_name)

type(grid_type), intent(inout) :: grid
character(len=*), intent(out) :: lon_name, lat_name

lon_name = grid%lon_name
lat_name = grid%lat_name

end subroutine


!---------------------------------------------
! call after type and sizes are set
subroutine allocate_grid_space(grid)

type(grid_type), intent(inout) :: grid

! add checks here

if (grid%type == GRID_REGULAR) then
   allocate(grid%rrlon(grid%nlon), grid%rrlat(grid%nlat))
else
   allocate(grid%irlon(grid%nlon, grid%nlat))
   allocate(grid%irlat(grid%nlon, grid%nlat))
endif

end subroutine

!---------------------------------------------

end module grid_mod
