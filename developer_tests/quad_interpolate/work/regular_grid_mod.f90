! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module regular_grid_mod

use types_mod, only : r8

implicit none

private

public :: create_grid0, create_field

contains


!---------------------------------------------
subroutine create_grid0(resolution, lon, lat, n)

real(r8), intent(in) :: resolution ! decimal degree 1, 0.1
real(r8), allocatable, intent(out) :: lon(:)
real(r8), allocatable, intent(out) :: lat(:)
integer, intent(out) :: n

integer :: i

n = floor(360.d0 / resolution)

allocate(lon(n), lat(n))

lon(1) = 0

do i = 2, n
   lon(i) = lon(i-1) + resolution
enddo

lat = lon

end subroutine create_grid0


!---------------------------------------------
subroutine create_field(n, lon, lat, field, func)

real(r8), intent(in) :: lon(n)
real(r8), intent(in) :: lat(n)
real(r8), intent(out) :: field(n,n)

interface 
   pure function func(x,y)
      use types_mod, only: r8
      real(r8), intent(in) :: x,y
      real(r8) :: func
   end function
end interface

integer :: i,j,n

do i = 1, n
   do j = 1, n
      field(i,j) = func(lon(i), lat(j))
   enddo
enddo

end subroutine create_field

end module regular_grid_mod
