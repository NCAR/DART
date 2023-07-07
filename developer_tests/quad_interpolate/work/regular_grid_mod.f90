! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module regular_grid_mod

use types_mod, only : r8
use grid_mod

implicit none

private

public :: create_grid0, create_field

contains


!---------------------------------------------
subroutine create_grid0(resolution, grid, ni, nj)

real(r8), intent(in) :: resolution ! decimal degree 1, 0.1
type(grid_type), intent(inout) :: grid
integer, intent(out) :: ni, nj


integer :: i, j

! regular lat/lon grid goes from [0-360), [-90,90]

ni = floor(360.d0 / resolution)
nj = floor(180.d0 / resolution) 

call set_grid_type_regular(grid)
call set_grid_sizes(grid, ni, nj)
call set_grid_names(grid, "lon", "lat")
call allocate_grid_space(grid)

! make a cover routine for setting lon/lat arrays

grid%irlon(1) = 0.0_r8

do i = 2, ni
   grid%irlon(i) = grid%irlon(i-1) + resolution
enddo


grid%irlat(1) = -90.0_r8

do j = 2, nj
   grid%irlat(j) = grid%irlat(j-1) + resolution
enddo

end subroutine create_grid0


!---------------------------------------------
subroutine create_field(grid, field, func)

type(grid_type), intent(inout) :: grid
real(r8), intent(inout) :: field(:,:)

interface 
   pure function func(x,y)
      use types_mod, only: r8
      real(r8), intent(in) :: x,y
      real(r8) :: func
   end function
end interface

integer :: i,j,ni,nj

call get_grid_sizes(grid, ni, nj)

do j = 1, nj
   do i = 1, ni
      field(i,j) = func(grid%irlon(i), grid%irlat(j))
   enddo
enddo

end subroutine create_field

end module regular_grid_mod
