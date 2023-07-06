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
subroutine create_grid0(resolution, grid, n)

real(r8), intent(in) :: resolution ! decimal degree 1, 0.1
type(grid_type), intent(inout) :: grid
integer, intent(out) :: n


integer :: i

n = floor(360.d0 / resolution)

! ni = floor(360.d0 / resolution)
! nj = floor(180.d0 / resolution) 

print *, 'starting create grid, n = ', n
print *, 'passed in resolution = ', resolution

call set_grid_type_regular(grid)
call set_grid_sizes(grid, n, n)
call set_grid_names(grid, "lon", "lat")
call allocate_grid_space(grid)

! make a cover routine for setting lon/lat arrays

grid%irlon(1) = 0.0_r8

do i = 2, n
   grid%irlon(i) = grid%irlon(i-1) + resolution
enddo


grid%irlat(1) = -90.0_r8

do i = 2, n
   grid%irlat(i) = grid%irlat(i-1) + (resolution/2.0_r8)
enddo

print *, 'leaving create grid, n = ', n

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
print *, 'got grid sizes of ', ni, nj

do i = 1, ni
   do j = 1, nj
print *, i, j, field(i, j),  grid%irlon(i), grid%irlat(j)
      field(i,j) = func(grid%irlon(i), grid%irlat(j))
   enddo
enddo

end subroutine create_field

end module regular_grid_mod
