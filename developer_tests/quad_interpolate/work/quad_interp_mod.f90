! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module quad_interp_mod

! this needs 4 variants of the routines:  
!   from regular to regular
!   from regular to irregular
!   from irregular to regular
!   from irregular to irregular
! 
! the differences are in how you initialize the quad code
! based on the 'from' grid, and how you access the destination
! locations on the 'to' grid.,
!

use        types_mod, only : r8, i8, MISSING_R8

use grid_mod

use quad_utils_mod, only : quad_interp_handle, init_quad_interp, finalize_quad_interp, set_quad_coords,       &
                           quad_lon_lat_locate, quad_lon_lat_evaluate, GRID_QUAD_FULLY_REGULAR,               &
                           GRID_QUAD_IRREG_SPACED_REGULAR, GRID_QUAD_FULLY_IRREGULAR, GRID_QUAD_UNKNOWN_TYPE, &
                           QUAD_LOCATED_UNKNOWN, QUAD_LOCATED_CELL_CENTERS, QUAD_LOCATED_LON_EDGES,           &
                           QUAD_LOCATED_LAT_EDGES, QUAD_LOCATED_CELL_CORNERS


implicit none
private


public :: do_interp, set_quad_grid_opts

type(quad_interp_handle) :: h

integer  :: i, j, k
integer  :: nx, ny
integer  :: nrx, nry
integer  :: four_lons(4), four_lats(4)

!real(r8) :: data_del_lon, data_del_lat, sample_del_lon, sample_del_lat
integer  :: four_lon_inds(4), four_lat_inds(4)
real(r8) :: lon_fract, lat_fract
integer  :: istatus
real(r8) :: invals(4), outval

logical :: grid_global = .true.
logical :: grid_spans_lon_zero = .true.
logical :: grid_pole_wrap = .true.

contains

!------------------------------------------------------------------
subroutine set_quad_grid_opts(global, spans_lon_zero, pole_wrap)

logical :: global
logical :: spans_lon_zero
logical :: pole_wrap
logical :: grid_rotate

grid_global = global
grid_spans_lon_zero = spans_lon_zero
grid_pole_wrap = pole_wrap

end subroutine set_quad_grid_opts

!------------------------------------------------------------------
subroutine do_interp(from, fromfield, to, tofield)

type(grid_type), intent(in) :: from
real(r8), intent(in) :: fromfield(:,:)
type(grid_type), intent(inout) :: to
real(r8), intent(out) :: tofield(:,:)

logical :: from_reg, to_reg
 
from_reg = is_grid_type_regular(from)
to_reg   = is_grid_type_regular(to)
 
!print *, 'in do_interp'
!print *, 'is_reg for from, to = ', from_reg, to_reg
!call dump_grid(from, label = 'from')
!call dump_grid(to, label = 'to')

if (from_reg) then
   if (to_reg) then
      call do_reg_reg_interp(from, fromfield, to, tofield)
   else
      call do_reg_irreg_interp(from, fromfield, to, tofield)
   endif
else
   if (to_reg) then
      call do_irreg_reg_interp(from, fromfield, to, tofield)
   else
      call do_irreg_irreg_interp(from, fromfield, to, tofield)
   endif
endif

end subroutine do_interp

!------------------------------------------------------------------

subroutine do_reg_reg_interp(from, fromfield, to, tofield)

type(grid_type), intent(in) :: from
real(r8), intent(in) :: fromfield(:,:)
type(grid_type), intent(inout) :: to
real(r8), intent(out) :: tofield(:,:)


call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, from%nlon, from%nlat, &
                      QUAD_LOCATED_CELL_CENTERS, grid_global, grid_spans_lon_zero, grid_pole_wrap, h)
call set_quad_coords(h, from%irlon, from%irlat)

do j=1, to%nlat
   do i=1, to%nlon

      call quad_lon_lat_locate(h, to%irlon(i), to%irlat(j), four_lon_inds, four_lat_inds, &
                               lon_fract, lat_fract, istatus)
      if (istatus /= 0) then
         tofield(i, j) = MISSING_R8 
         cycle
      endif

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      do k=1, 4
         invals(k) = fromfield(four_lon_inds(k), four_lat_inds(k))
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

!------------------------------------------------------------

subroutine do_reg_irreg_interp(from, fromfield, to, tofield)

type(grid_type), intent(in) :: from
real(r8), intent(in) :: fromfield(:,:)
type(grid_type), intent(inout) :: to
real(r8), intent(out) :: tofield(:,:)


call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, from%nlon, from%nlat, &
                      QUAD_LOCATED_CELL_CENTERS, grid_global, grid_spans_lon_zero, grid_pole_wrap, h)
call set_quad_coords(h, from%irlon, from%irlat)


do j=1, to%nlat
   do i=1, to%nlon

      call quad_lon_lat_locate(h, to%iilon(i,j), to%iilat(i,j), four_lon_inds, four_lat_inds, &
                               lon_fract, lat_fract, istatus)
      if (istatus /= 0) then
         tofield(i, j) = MISSING_R8 
         cycle
      endif

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      do k=1, 4
         invals(k) = fromfield(four_lon_inds(k), four_lat_inds(k))
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

!------------------------------------------------------------

subroutine do_irreg_reg_interp(from, fromfield, to, tofield)

type(grid_type), intent(in) :: from
real(r8), intent(in) :: fromfield(:,:)
type(grid_type), intent(inout) :: to
real(r8), intent(out) :: tofield(:,:)


call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, from%nlon, from%nlat, &
                      QUAD_LOCATED_CELL_CENTERS, grid_global, grid_spans_lon_zero, grid_pole_wrap, h)
call set_quad_coords(h, from%iilon, from%iilat)


do j=1, to%nlat
   do i=1, to%nlon

      call quad_lon_lat_locate(h, to%irlon(i), to%irlat(j), four_lon_inds, four_lat_inds, istatus)
      if (istatus /= 0) then
         tofield(i, j) = MISSING_R8 
         cycle
      endif

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      do k=1, 4
         invals(k) = fromfield(four_lon_inds(k), four_lat_inds(k))
      enddo

      if (any(invals == MISSING_R8)) then
         tofield(i, j) = MISSING_R8
         cycle
      endif

      call quad_lon_lat_evaluate(h, to%irlon(i), to%irlat(j), four_lon_inds, four_lat_inds, &
                                  invals, outval, istatus)

      tofield(i, j) = outval

   enddo
enddo

call finalize_quad_interp(h)

end subroutine

!------------------------------------------------------------

subroutine do_irreg_irreg_interp(from, fromfield, to, tofield)

type(grid_type), intent(in) :: from
real(r8), intent(in) :: fromfield(:,:)
type(grid_type), intent(inout) :: to
real(r8), intent(out) :: tofield(:,:)


call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, from%nlon, from%nlat, &
                      QUAD_LOCATED_CELL_CENTERS, grid_global, grid_spans_lon_zero, grid_pole_wrap, h)
call set_quad_coords(h, from%iilon, from%iilat)


do j=1, to%nlat
   do i=1, to%nlon

      call quad_lon_lat_locate(h, to%iilon(i,j), to%iilat(i,j), four_lon_inds, four_lat_inds, istatus)
      if (istatus /= 0) then
         tofield(i, j) = MISSING_R8 
         cycle
      endif

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      do k=1, 4
         invals(k) = fromfield(four_lon_inds(k), four_lat_inds(k))
      enddo

      if (any(invals == MISSING_R8)) then
         tofield(i, j) = MISSING_R8
         cycle
      endif

      call quad_lon_lat_evaluate(h, to%iilon(i,j), to%iilat(i,j), four_lon_inds, four_lat_inds, &
                                  invals, outval, istatus)

      tofield(i, j) = outval

   enddo
enddo

call finalize_quad_interp(h)

end subroutine

!------------------------------------------------------------

end module quad_interp_mod

