! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


program test_interpolate_grid

use regular_grid_mod, only : create_grid0, create_field
use types_mod, only : r8
use functions_mod, only : sine
use target_grid_mod, only : create_gridT, target_grid

implicit none

! irregular grids have lon(:,:) lat(:,:)
! so the split between reg and irreg has to
! be somewhere.  do we need separate progs
! for interp to reg and interp to irreg ?


! maybe grid0 should also be a grid type?

real(r8), allocatable :: lon(:)
real(r8), allocatable :: lat(:)
real(r8), allocatable :: field(:,:)

real(r8) :: resolution
integer  :: n
integer  :: i,j ! loop variables

resolution = 1
call create_grid0(resolution, lon, lat, n)
allocate(field(n,n))
call create_field(n, lon, lat, field, sine)

! export FORT_FMT_RECL=10240 for easy print with ifort
do i = 1, n
  print*, (field(i,j), j = 1, n)
enddo

end program test_interpolate_grid
