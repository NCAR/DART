! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module regular_grid_mod

use types_mod, only : r8
use functions_mod, only : f_sine, f_sum, f_row, f_col, f_rand, f_long, f_lati
use grid_mod

implicit none

private

public :: create_grid0, create_field, set_data_function


interface 
   pure function func(x,y)
      use types_mod, only: r8
      real(r8), intent(in) :: x,y
      real(r8) :: func
   end function
end interface
 
procedure(func), pointer :: func_ptr => NULL()


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
subroutine create_field(grid, field)

type(grid_type), intent(inout) :: grid
real(r8), intent(inout) :: field(:,:)

integer :: i,j,ni,nj

call get_grid_sizes(grid, ni, nj)

do j = 1, nj
   do i = 1, ni
      field(i,j) = func_ptr(grid%irlon(i), grid%irlat(j))
print *, "i, j, f =", i, j, field(i,j)
   enddo
enddo

end subroutine create_field

!---------------------------------------------
subroutine set_data_function(func_name)

character(len=*), intent(in) :: func_name

!interface 
!   pure function the_func(x,y)
!      use types_mod, only: r8
!      real(r8), intent(in) :: x,y
!      real(r8) :: func
!   end function
!end interface

select case ( func_name )

case ("sine")
   func_ptr => f_sine

case ("sum")
   func_ptr => f_sum

case ("row")
   func_ptr => f_row

case ("column")
   func_ptr => f_col

case ("random")
   func_ptr => f_rand

case ("longitude")
   func_ptr => f_long

case ("latitude")
   func_ptr => f_lati

case default 
   print *, "data function must be one of:"
   print *, "sine, sum, row, column, random, longitude, latitude"
   call exit(-1)

end select


end subroutine set_data_function

!---------------------------------------------

end module regular_grid_mod
