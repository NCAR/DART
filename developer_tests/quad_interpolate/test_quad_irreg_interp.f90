! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program test_quad_irreg_interp

! intended to show how the state structure and quad code can be used
! together.  start with a simple regional grid and work out from there.

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, i8, MISSING_R8, deg2rad, rad2deg
use    utilities_mod, only : error_handler, initialize_utilities, finalize_utilities
use   random_seq_mod, only : init_random_seq, random_seq_type, &
                             random_uniform, random_gaussian

use quad_utils_mod, only : quad_interp_handle, init_quad_interp, finalize_quad_interp, set_quad_coords,       &
                           quad_lon_lat_locate, quad_lon_lat_evaluate, GRID_QUAD_FULLY_REGULAR,               &
                           GRID_QUAD_IRREG_SPACED_REGULAR, GRID_QUAD_FULLY_IRREGULAR, GRID_QUAD_UNKNOWN_TYPE, &
                           QUAD_LOCATED_UNKNOWN, QUAD_LOCATED_CELL_CENTERS, QUAD_LOCATED_LON_EDGES,           &
                           QUAD_LOCATED_LAT_EDGES, QUAD_LOCATED_CELL_CORNERS


implicit none

integer :: debug = 0

type(quad_interp_handle) :: h

! data grid size
integer, parameter :: nx = 9
integer, parameter :: ny = 5

! locations of data grid corners
real(r8) :: data_lons(nx, ny) = MISSING_R8
real(r8) :: data_lats(nx, ny) = MISSING_R8

! extents of the data grid (these mimic a regional model's grid)
real(r8) :: start_lon = 100.0_r8
real(r8) :: end_lon   = 150.5_r8
real(r8) :: start_lat = -11.4_r8
real(r8) :: end_lat   =  34.1_r8

! angle to rotate data grid in degrees
! positive is counterclockwise; will rotate
! around lower left grid point (start lon/lat).
 real(r8) :: angle = 10.0_r8
!real(r8) :: angle = 45.0_r8
!real(r8) :: angle = 30.0_r8
!real(r8) :: angle =  90.0_r8
!real(r8) :: angle = -30.0_r8
!real(r8) :: angle = -10.0_r8
!real(r8) :: angle = 0.0_r8

! deform grid by this fraction of the deltas
real(r8) :: lon_def = 0.25_r8
real(r8) :: lat_def = 0.25_r8

! data values on the grid
real(r8) :: grid_data(nx, ny) = MISSING_R8
integer  :: data_choice = 0   ! see code for selection values

! percent of data values that should be marked 'missing data'
!real(r8) :: miss_percent =   0.0_r8    ! none
 real(r8) :: miss_percent =   3.0_r8    ! 3%
!real(r8) :: miss_percent = 100.0_r8    ! all

! sampling grid size
integer, parameter :: nrx = 210
integer, parameter :: nry = 150

! locations of sampling grid
real(r8) :: sample_lons(nrx) = MISSING_R8
real(r8) :: sample_lats(nry) = MISSING_R8

! extents of the sampling grid
real(r8) :: sample_start_lon = 110.0_r8
real(r8) :: sample_end_lon   = 140.0_r8
real(r8) :: sample_start_lat = -20.0_r8
real(r8) :: sample_end_lat   =  30.0_r8

! where interpolated values are stored on reg grid
real(r8) :: interp_data(nrx, nry) = MISSING_R8


type(random_seq_type) :: ran

integer  :: i, j, k
real(r8) :: lon_del, lat_del, sample_lon_del, sample_lat_del
integer  :: lon_bot, lat_bot, lon_top, lat_top
integer  :: four_lons(4), four_lats(4)
integer  :: istatus
real(r8) :: invals(4), outval
integer  :: iunit_orig, iunit_interp

call initialize_utilities('test_quad_irreg_interp')
call init_random_seq(ran)

lon_del = (end_lon - start_lon) / (nx-1)
lat_del = (end_lat - start_lat) / (ny-1)

! "data grid" corners and data vals
do i=1, nx
   do j=1, ny
      ! generate locations of the corners of all the quads
      data_lons(i, j) = start_lon + (i-1)*lon_del + deform(lon_del, lon_def, ran)
      data_lats(i, j) = start_lat + (j-1)*lat_del + deform(lat_del, lat_def, ran)
      if (angle /= 0.0_r8) &
         call rotate(data_lons(i, j), data_lats(i, j), angle, start_lon, start_lat)

      ! generate the data values on the corners.  pick one:
      select case (data_choice) 
      case (1)
         ! increasing monotonically 
         grid_data(i, j) = (j-1)*nx + i
      case (2)
         ! constant by row
         grid_data(i, j) = j
      case (3)
         ! constant by column
         grid_data(i, j) = i
      case (4)
         ! based on lon only
         grid_data(i, j) = data_lons(i, j)
      case (5) 
         ! based on lat only
         grid_data(i, j) = data_lats(i, j)
      case (6)
         ! random between (0-10)
         grid_data(i, j) = random_uniform(ran) * 10.0_r8
      case default
         ! gaussian with mean 0 and stddev 1
         grid_data(i, j) = random_gaussian(ran, 0.0_r8, 1.0_r8)
      end select

      if (miss_percent > 0.0_r8) then
        if (random_uniform(ran) * 100.0_r8 < miss_percent) grid_data(i, j) = MISSING_R8 
      endif

   enddo
enddo

sample_lon_del = (sample_end_lon - sample_start_lon) / nrx
sample_lat_del = (sample_end_lat - sample_start_lat) / nry

! "sampled grid" spacing along each axis
do i=1, nrx
   sample_lons(i) = sample_start_lon + (i-1)*sample_lon_del
enddo
do j=1, nry
   sample_lats(j) = sample_start_lat + (j-1)*sample_lat_del
enddo

! end of data setup - now call interp routines

call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, QUAD_LOCATED_CELL_CENTERS, .false., .false., .false., h)
call set_quad_coords(h, data_lons, data_lats)

! for each location in the sampling grid, interpolate a data value
do i=1, nrx
   do j=1, nry
      call quad_lon_lat_locate(h, sample_lons(i), sample_lats(j), four_lons, four_lats, istatus)
      if (istatus /= 0) then
         !print *, 'location outside of grid: ', sample_lons(i), sample_lats(j)
         interp_data(i, j) = MISSING_R8 
         cycle
      endif
      if (debug > 0) print *, i, j, four_lons(1), four_lons(3), four_lats(1), four_lats(3), sample_lons(i), sample_lats(j)

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      do k=1, 4
         invals(k) = grid_data(four_lons(k), four_lats(k))
      enddo

      call quad_lon_lat_evaluate(h, sample_lons(i), sample_lats(j), four_lons, four_lats, &
                                 invals, outval, istatus)

      interp_data(i, j) = outval

   enddo
enddo

! this program doesn't currently have any missing locations - but i'll test that next.
if (debug > 0) print *, 'number of missing values in  input data: ', count(grid_data(:,:) == MISSING_R8)
if (debug > 0) print *, 'number of missing values in output data: ', count(interp_data(:,:) == MISSING_R8)

call writeit_2d('data_lons_2d_irreg_test.txt', nx, ny, data_lons)
call writeit_2d('data_lats_2d_irreg_test.txt', nx, ny, data_lats)
call writeit_2d('data_data_2d_irreg_test.txt', nx, ny, grid_data)

call writeit_1d('sample_lons_1d_irreg_test.txt', nrx, sample_lons)
call writeit_1d('sample_lats_1d_irreg_test.txt', nry, sample_lats)
call writeit_2d('sample_data_2d_irreg_test.txt', nrx, nry, interp_data)

call finalize_quad_interp(h)
call finalize_utilities('test_quad_irreg_interp')

if (debug > 0) print *, 'closed files and finalized interp handle'

contains

!------------------------------------------------------------
! rotate vector a counterclockwise by angle theta, relative
! to the given origin point.

subroutine rotate(x, y, theta, x0, y0)
 real(r8), intent(inout) :: x, y
 real(r8), intent(in)    :: theta
 real(r8), intent(in)    :: x0, y0

real(r8) :: a(2), b(2)
real(r8) :: r(2,2)
real(r8) :: rads

a(1) = x - x0
a(2) = y - y0

rads = theta * deg2rad

r(1,1) = cos(rads)
r(1,2) = sin(rads)
r(2,1) = sin(-rads)
r(2,2) = cos(rads)

b(1) = r(1,1)*a(1) + r(1,2)*a(2)
b(2) = r(2,1)*a(1) + r(2,2)*a(2)

x = b(1) + x0
y = b(2) + y0

end subroutine rotate

!------------------------------------------------------------
! compute +/- a random value based on a width and percentage
! of that width

function deform(width, fraction, seq)

use random_seq_mod

 real(r8), intent(in) :: width
 real(r8), intent(in) :: fraction
 type(random_seq_type), intent(inout) :: seq
 real(r8)             :: deform

real(r8) :: val

! random val between -1 and 1
val = (random_uniform(seq) * 2.0_r8) - 1.0_r8

deform = val * width * fraction

end function deform

!------------------------------------------------------------

subroutine writeit_1d(fname, nx, dataarray)

use    utilities_mod

 character(len=*), intent(in) :: fname
 integer, intent(in) :: nx
 real(r8), intent(in) :: dataarray(:)

integer :: i, j, iunit

iunit = open_file(fname, action='write')

do i=1, nx
   write(iunit, *) dataarray(i)
enddo

call close_file(iunit)

end subroutine writeit_1d

!------------------------------------------------------------

subroutine writeit_2d(fname, nx, ny, dataarray)

use    utilities_mod

 character(len=*), intent(in) :: fname
 integer, intent(in) :: nx, ny
 real(r8), intent(in) :: dataarray(nx, ny)

integer :: i, j, iunit

iunit = open_file(fname, action='write')

do i=1, nx
   do j=1, ny
      write(iunit, *) dataarray(i, j)
   enddo
enddo

call close_file(iunit)

end subroutine writeit_2d

!------------------------------------------------------------

end program test_quad_irreg_interp

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
