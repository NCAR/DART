! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module quad_interp_mod

! use the quad interpolate code on a fully regular grid
! and on an irregular grid.

use        types_mod, only : r8, i8, MISSING_R8

use grid_mod

use quad_utils_mod, only : quad_interp_handle, init_quad_interp, finalize_quad_interp, set_quad_coords,       &
                           quad_lon_lat_locate, quad_lon_lat_evaluate, GRID_QUAD_FULLY_REGULAR,               &
                           GRID_QUAD_IRREG_SPACED_REGULAR, GRID_QUAD_FULLY_IRREGULAR, GRID_QUAD_UNKNOWN_TYPE, &
                           QUAD_LOCATED_UNKNOWN, QUAD_LOCATED_CELL_CENTERS, QUAD_LOCATED_LON_EDGES,           &
                           QUAD_LOCATED_LAT_EDGES, QUAD_LOCATED_CELL_CORNERS


implicit none
private


public :: do_reg_interp

type(quad_interp_handle) :: h

integer  :: i, j, k
integer  :: nx, ny
integer  :: nrx, nry
integer  :: four_lons(4), four_lats(4)

real(r8) :: data_del_lon, data_del_lat, sample_del_lon, sample_del_lat
integer  :: lon_indices(4), lat_indices(4)
real(r8) :: lon_fract, lat_fract
integer  :: istatus
real(r8) :: invals(4), outval

contains

subroutine do_reg_interp(from, fromfield, to, tofield)

type(grid_type), intent(in) :: from
real(r8), intent(in) :: fromfield(:,:)
type(grid_type), intent(inout) :: to
real(r8), intent(out) :: tofield(:,:)



call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, from%nlon, from%nlat, QUAD_LOCATED_CELL_CENTERS, .false., .false., .false., h)
call set_quad_coords(h, from%irlon, from%irlat)


do j=1, to%nlat
   do i=1, to%nlon

      call quad_lon_lat_locate(h, to%irlon(i), to%irlat(j), lon_indices, lat_indices, &
                               lon_fract, lat_fract, istatus)
      if (istatus /= 0) then
         tofield(i, j) = MISSING_R8 
         cycle
      endif

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      do k=1, 4
         invals(k) = fromfield(lon_indices(k), lat_indices(k))
      enddo

      if (any(invals == MISSING_R8)) then
         tofield(i, j) = MISSING_R8
         cycle
      endif

      call quad_lon_lat_evaluate(h, lon_fract, lat_fract, invals, outval, istatus)

      tofield(i, j) = outval

   enddo
enddo

call finalize_quad_interp(h)

end subroutine


!%! !------------------------------------------------------------
!%! 
!%! subroutine do_irreg_interp()
!%! 
!%! call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, QUAD_LOCATED_CELL_CENTERS, .false., .false., .false., h)
!%! call set_quad_coords(h, data_lons, data_lats)
!%! 
!%! ! for each location in the sampling grid, interpolate a data value
!%! do i=1, nrx
!%!    do j=1, nry
!%!       call quad_lon_lat_locate(h, sample_lons(i), sample_lats(j), lon_indices, lat_indices, istatus)
!%! 
!%!       if (istatus /= 0) then
!%!          interp_data(i, j) = MISSING_R8 
!%!          cycle
!%!       endif
!%! 
!%!       ! get values of data at lon/lat bot/top indices, counterclockwise around quad
!%!       do k=1, 4
!%!          invals(k) = grid_data(lon_indices(k), lat_indices(k))
!%!       enddo
!%! 
!%!       call quad_lon_lat_evaluate(h, sample_lons(i), sample_lats(j), four_lons, four_lats, &
!%!                                  invals, outval, istatus)
!%! 
!%!       interp_data(i, j) = outval
!%! 
!%!    enddo
!%! enddo
!%! 
!%! call finalize_quad_interp(h)
!%! 
!%! end subroutine
!%! 
!------------------------------------------------------------

end module quad_interp_mod

