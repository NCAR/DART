! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program test_quad_reg_interp

! example of how to use the quad interpolate code on a fully regular grid.

!>@todo FIXME include the state struct or fold this into a version of model_mod_check

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

! data grid.  these are the values we will interpolate from.

! data grid size 
! (we compute delta lon, delta lat based on these vals)
integer, parameter :: ndx = 9
integer, parameter :: ndy = 5

! extents of the data grid (these mimic a regional model's grid)
real(r8) :: data_start_lon = 100.0_r8
real(r8) :: data_end_lon   = 150.5_r8
real(r8) :: data_start_lat = -11.4_r8
real(r8) :: data_end_lat   =  34.1_r8

! these aren't needed for the interpolation, but they're written
! out for ease of plotting the results
real(r8) :: grid_lons(ndx)
real(r8) :: grid_lats(ndy)

! data values on the grid
real(r8) :: grid_data(ndx, ndy) = MISSING_R8
integer  :: data_choice = 0   ! see code for selection values

! percent of data values that should be marked 'missing data'
!real(r8) :: miss_percent =   0.0_r8    ! none
 real(r8) :: miss_percent =   3.0_r8    ! 3%
!real(r8) :: miss_percent = 100.0_r8    ! all


! interpolation test grid.  we construct a different grid
! and call the interpolation code on each corner of this
! other grid.  called 'sampling grid' to differentiate it
! from the 'data grid'.  usually much denser so we can look
! for discontinuties or errors in the interp code.

! sampling grid size
integer, parameter :: nsx = 210
integer, parameter :: nsy = 150

! locations of sampling grid
real(r8) :: sample_lons(nsx, nsy) = MISSING_R8
real(r8) :: sample_lats(nsx, nsy) = MISSING_R8

! extents of the sampling grid
real(r8) :: sample_start_lon = 110.0_r8
real(r8) :: sample_end_lon   = 140.0_r8
real(r8) :: sample_start_lat = -20.0_r8
real(r8) :: sample_end_lat   =  30.0_r8

! angle to rotate sampling grid in degrees
! positive is counterclockwise; will rotate
! around lower left grid point (start lon/lat).
!real(r8) :: angle =  10.0_r8
!real(r8) :: angle =  45.0_r8
 real(r8) :: angle = -65.0_r8
!real(r8) :: angle =  30.0_r8
!real(r8) :: angle =  90.0_r8
!real(r8) :: angle = -30.0_r8
!real(r8) :: angle = -10.0_r8
!real(r8) :: angle =   0.0_r8

! deform grid by this fraction of the deltas
!real(r8) :: lon_def = 0.25_r8
!real(r8) :: lat_def = 0.25_r8
 real(r8) :: lon_def = 0.01_r8
 real(r8) :: lat_def = 0.01_r8
!real(r8) :: lon_def = 0.0_r8
!real(r8) :: lat_def = 0.0_r8

! where interpolated values are stored on reg grid
real(r8) :: interp_data(nsx, nsy) = MISSING_R8


type(random_seq_type) :: ran

integer  :: i, j
real(r8) :: data_del_lon, data_del_lat, sample_del_lon, sample_del_lat
integer  :: lon_indices(4), lat_indices(4)
real(r8) :: lon_fract, lat_fract
integer  :: istatus
real(r8) :: invals(4), outval
integer  :: iunit_orig, iunit_interp

call initialize_utilities('test_quad_reg_interp')
call init_random_seq(ran)


! "data grid" corners and data vals

data_del_lon = (data_end_lon - data_start_lon) / ndx
data_del_lat = (data_end_lat - data_start_lat) / ndy

do i=1, ndx
   do j=1, ndy
      ! generate the data values on the corners.  pick one:
      select case (data_choice) 
      case (1)
         ! increasing monotonically 
         grid_data(i, j) = (j-1)*ndx + i
      case (2)
         ! constant by row
         grid_data(i, j) = j
      case (3)
         ! constant by column
         grid_data(i, j) = i
      case (4)
         ! based on lon only
         grid_data(i, j) = data_start_lon + (data_del_lon * (i-1))
      case (5) 
         ! based on lat only
         grid_data(i, j) = data_start_lat + (data_del_lat * (j-1))
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

sample_del_lon = (sample_end_lon - sample_start_lon) / nsx
sample_del_lat = (sample_end_lat - sample_start_lat) / nsy

! "sampled grid" spacing along each axis
do i=1, nsx
   do j=1, nsy
      sample_lons(i, j) = sample_start_lon + (i-1)*sample_del_lon + deform(sample_del_lon, lon_def, ran)
      sample_lats(i, j) = sample_start_lat + (j-1)*sample_del_lat + deform(sample_del_lat, lat_def, ran)

      ! generate locations of the corners of all the quads
      if (angle /= 0.0_r8) &
         call rotate(sample_lons(i, j), sample_lats(i, j), angle, sample_start_lon, sample_start_lat)
enddo
enddo

! end of data setup - now call interp routines


call init_quad_interp(GRID_QUAD_FULLY_REGULAR, ndx, ndy, QUAD_LOCATED_CELL_CENTERS, .false., .false., .false., h)
call set_quad_coords(h, data_start_lon, data_del_lon, data_start_lat, data_del_lat)

do i=1, nsx
   do j=1, nsy

      !>this interface now returns an array of 4 index combinations 
      !>so the calling code can do loops from 1 to 4 instead of making 
      !>combinations of lat/lon bot/top in the right order for eval. 

      call quad_lon_lat_locate(h, sample_lons(i,j), sample_lats(i,j), lon_indices, lat_indices, &
                               lon_fract, lat_fract, istatus)
      if (istatus /= 0) then
         !print *, 'location outside of grid: ', sample_lons(i,j), sample_lats(i,j)
         interp_data(i, j) = MISSING_R8 
         cycle
      endif
      if(debug > 0)print *, i, j, lon_indices, lat_indices, lon_fract, lat_fract, sample_lons(i,j), sample_lats(i,j)

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      invals(1) = grid_data(lon_indices(1), lat_indices(1))
      invals(2) = grid_data(lon_indices(2), lat_indices(2))
      invals(3) = grid_data(lon_indices(3), lat_indices(3))
      invals(4) = grid_data(lon_indices(4), lat_indices(4))

      ! where does this go?  inside quad_lon_lat_evaluate() or here?
      if (any(invals == MISSING_R8)) then
         interp_data(i, j) = MISSING_R8
         cycle
      endif

      call quad_lon_lat_evaluate(h, lon_fract, lat_fract, invals, outval, istatus)

      if (istatus == 0) then
         interp_data(i, j) = outval
      else
         interp_data(i, j) = MISSING_R8
      endif

   enddo
enddo

! this program doesn't currently have any missing locations - but i'll test that next.
if (debug > 0) print *, 'number of missing values in  input data: ', count(grid_data(:,:) == MISSING_R8)
if (debug > 0) print *, 'number of missing values in output data: ', count(interp_data(:,:) == MISSING_R8)

! generate these only for output for plotting the results
! these aren't needed for the interpolation.
do i=1,ndx
   grid_lons(i) = data_start_lon + (data_del_lon * (i-1))
enddo
do j=1,ndy
   grid_lats(j) = data_start_lat + (data_del_lat * (j-1))
enddo

call writeit_1d('data_lons_1d_reg_test.txt', ndx, grid_lons)
call writeit_1d('data_lats_1d_reg_test.txt', ndy, grid_lats)
call writeit_2d('data_data_2d_reg_test.txt', ndx, ndy, grid_data)

call writeit_2d('sample_lons_2d_reg_test.txt', nsx, nsy, sample_lons)
call writeit_2d('sample_lats_2d_reg_test.txt', nsx, nsy, sample_lats)
call writeit_2d('sample_data_2d_reg_test.txt', nsx, nsy, interp_data)

call finalize_quad_interp(h)
if (debug > 0) print *, 'closed files and finalized interp handle'

call finalize_utilities('test_quad_reg_interp')


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

end program test_quad_reg_interp

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
