! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength, DEG2RAD, RAD2DEG, &
                             vtablenamelength

use time_manager_mod, only : time_type, set_time, get_time

use     location_mod, only : location_type, get_close_type, get_dist, &
                             loc_get_close_obs => get_close_obs,      &
                             loc_get_close_state => get_close_state,  &
                             set_location, query_location,            &
                             get_location, VERTISHEIGHT, VERTISUNDEF, &
                             VERTISLEVEL

use    utilities_mod, only : error_handler, E_MSG, E_ALLMSG,  &
                             nmlfileunit, do_nml_file, do_nml_term,                &
                             find_namelist_in_file, check_namelist_read,           &
                             find_enclosing_indices

use obs_kind_mod,     only : QTY_GEOMETRIC_HEIGHT

use netcdf_utilities_mod,  only : nc_add_global_attribute, nc_synchronize_file,           &
                                 nc_add_global_creation_time, nc_begin_define_mode,       & 
                                 nc_end_define_mode, nc_open_file_readonly, nc_close_file, &
                                 nc_get_dimension_size, nc_get_variable

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, get_dart_vector_index, get_domain_size, &
                                get_model_variable_indices, get_varid_from_kind,      &
                                get_variable_name

use ensemble_manager_mod,  only : ensemble_type

use cube_sphere_grid_tools_mod, only : col_index_to_lat_lon,    &
                                       get_bounding_box,  get_grid_delta

! These routines are passed through from default_model_mod.
use default_model_mod, only : pert_model_copies, read_model_time, write_model_time,    &
                              init_time => fail_init_time,                             &
                              init_conditions => fail_init_conditions,                 &
                              convert_vertical_obs, convert_vertical_state, adv_1step, &
                              parse_variables_clamp, MAX_STATE_VARIABLE_FIELDS_CLAMP

implicit none
private

! routines required by DART code - will be called from filter and other DART executables. 
public :: get_model_size,                      &
          get_state_meta_data,                 &
          model_interpolate,                   &
          end_model,                           &
          static_init_model,                   &
          nc_write_model_atts,                 &
          get_close_obs,                       &
          get_close_state,                     &
          pert_model_copies,                   &
          convert_vertical_obs,                &
          convert_vertical_state,              &
          read_model_time,                     &
          adv_1step,                           &
          init_time,                           &
          init_conditions,                     &
          shortest_time_between_assimilations, &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "aether_cube_sphere/model_mod"

! Error codes
integer, parameter :: INVALID_VERT_COORD_ERROR_CODE   = 15
integer, parameter :: INVALID_ALTITUDE_VAL_ERROR_CODE = 17
integer, parameter :: INVALID_MODEL_LEVEL_ERROR_CODE  = 18
integer, parameter :: UNKNOWN_OBS_QTY_ERROR_CODE      = 20

! Error message strings                
character(len=512) :: string1
character(len=512) :: string2

integer  :: iunit, io

logical :: module_initialized = .false.
integer :: dom_id ! used to access the state structure
type(time_type) :: assimilation_time_step

! Geometry variables that are used throughout the module; read from a template file
integer               :: np                ! Number of grid rows across a face
real(r8)              :: del, half_del     ! Grid row spacing and half of that
integer               :: ncenter_altitudes ! The number of altitudes and the altitudes
real(r8), allocatable :: center_altitude(:)

! The basic geometry variables need to be accessed by the test_grid_box routine
! Also need get_state_index
public :: np, ncenter_altitudes, get_state_index

! Current model time needed for computing location for scalar F10.7
type(time_type) :: state_time

! Horizontal column dimension rather than being direct functions of latitude and longitude.
integer               :: no_third_dimension = -99

! Namelist for options to be set at runtime.
character(len=256)              :: template_file           = 'filter_input_0001.nc'
integer                         :: time_step_days          = 0
integer                         :: time_step_seconds       = 3600
integer, parameter              :: MAX_STATE_VARIABLES     = 100
integer, parameter              :: NUM_STATE_TABLE_COLUMNS = 5
character(len=vtablenamelength) ::                &
          variables(MAX_STATE_VARIABLE_FIELDS_CLAMP) = ' '  ! Table of state variables and associated metadata

namelist /model_nml/ template_file, time_step_days, time_step_seconds, variables

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. 

subroutine static_init_model()

module_initialized = .true.

! Read the namelist contents
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated.
assimilation_time_step = set_time(time_step_seconds, time_step_days)

! Define which variables are in the model state
dom_id = add_domain(template_file, parse_variables_clamp(variables))

! Get the altitudes and the number of grid rows
call read_template_file()

! Get the grid spacing; these quantities are in module storage
call get_grid_delta(np, del, half_del)

! NOTE TO AETHER MODELERS
! Need a way to get time from Aether for scalar F10.7 location definition
! Just set to arbitrary time for now
state_time = set_time(0, 1)

end subroutine static_init_model

!------------------------------------------------------------------

! Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = get_domain_size(dom_id)

end function get_model_size

!------------------------------------------------------------------

! Given a state handle, a location, and a state quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case a positive istatus should be returned.

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: qty
real(r8),            intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,             intent(out) :: istatus(ens_size)


! Vertical interpolation variables
integer     :: below_index, above_index, enclosing_status, which_vertical, level
integer     :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4), num_bound_points, var_id, i
integer(i8) :: bounding_state_index(4, 2)

real(r8)    :: lon_lat_alt(3), fract
real(r8)    :: grid_pt_lat(4), grid_pt_lon(4), pt_lat, pt_lon, bounding_value(4, 2, ens_size)
real(r8)    :: below_values(ens_size), above_values(ens_size)

! Initialize module if not already done
if(.not. module_initialized ) call static_init_model

! Set all obs to MISSING_R8 initially
expected_obs(:) = MISSING_R8
! Successful return status default
istatus(:) = 0

! Determine the vertical location type
lon_lat_alt    = get_location(location)
pt_lat         = lon_lat_alt(2) * DEG2RAD
pt_lon         = lon_lat_alt(1) * DEG2RAD
which_vertical = nint(query_location(location, 'WHICH_VERT'))

! Only heights currently supported for observations; fail if other is selected
if (.not. (which_vertical == VERTISHEIGHT .or. which_vertical == VERTISLEVEL )) then
   istatus = INVALID_VERT_COORD_ERROR_CODE
   return
endif

! Geometric height is a special case
if(qty == QTY_GEOMETRIC_HEIGHT) then
   level = nint(lon_lat_alt(3))
   if(level < 1 .or. level > ncenter_altitudes) then
      istatus = INVALID_MODEL_LEVEL_ERROR_CODE
   else
      expected_obs = center_altitude(nint(lon_lat_alt(3))) 
   endif
   return
endif

! See if the state contains the obs quantity
var_id = get_varid_from_kind(dom_id, qty)
if(var_id <= 0) then
   istatus = UNKNOWN_OBS_QTY_ERROR_CODE
   return
endif

! Find the bounding vertical levels and the fractional distance between
if(which_vertical == VERTISHEIGHT) then
   call find_enclosing_indices(ncenter_altitudes, center_altitude, lon_lat_alt(3), & 
      below_index, above_index, fract, enclosing_status)
   if (enclosing_status /= 0) then
      istatus = INVALID_ALTITUDE_VAL_ERROR_CODE
      return
   endif
else
   ! If VERTISLEVEL, bounds are both the vertical level
   ! This does twice the needed computation, could get rid of this
   above_index = abs(nint(lon_lat_alt(3)))
   ! Fail if level is outside of model bounds
   if(above_index < 1 .or. above_index > ncenter_altitudes) then
      istatus = INVALID_MODEL_LEVEL_ERROR_CODE
      return
   endif
   below_index = above_index
endif

! If the vertical location is acceptable, then do the horizontal interpolation
! Find the enclosing triangle or quad
call get_bounding_box(pt_lat, pt_lon, del, half_del, np, grid_face, &
   grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

! Map grid_face, latitude index and longitude index to 1D index used in the state vector
! Then get the state values
do i = 1, num_bound_points
   bounding_state_index(i, 1) =  get_state_index(grid_face(i), grid_lat_ind(i), &
      grid_lon_ind(i), below_index, var_id)
   bounding_value(i, 1, :) = get_state(bounding_state_index(i, 1), state_handle)
   bounding_state_index(i, 2) =  get_state_index(grid_face(i), grid_lat_ind(i), &
      grid_lon_ind(i), above_index, var_id)
   bounding_value(i, 2, :) = get_state(bounding_state_index(i, 2), state_handle)
enddo

! Do inverse distance weighted horizontal interpolation on both levels
below_values = idw_interp(ens_size, RAD2DEG*pt_lat, RAD2DEG*pt_lon, &
   RAD2DEG*grid_pt_lat, RAD2DEG*grid_pt_lon, bounding_value(:, 1, :), num_bound_points)
above_values = idw_interp(ens_size, RAD2DEG*pt_lat, RAD2DEG*pt_lon, &
   RAD2DEG*grid_pt_lat, RAD2DEG*grid_pt_lon, bounding_value(:, 2, :), num_bound_points)

! Do the vertical interpolation, linear in height to get final
call vert_interp(ens_size, below_values, above_values, fract, expected_obs)

end subroutine model_interpolate

!------------------------------------------------------------------

! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if(.not. module_initialized) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations

!------------------------------------------------------------------

! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: qty

! Local variables

integer  :: lev_index, col_index, my_var_id, my_qty, dom_id
real(r8) :: lat, lon
integer  :: seconds, days ! for f10.7 location
real(r8) :: longitude     ! for f10.7 location

if(.not. module_initialized) call static_init_model

call get_model_variable_indices(index_in, col_index, lev_index, &
   no_third_dimension, var_id=my_var_id, dom_id=dom_id, kind_index=my_qty)

! F10.7 scalar location does not have a column
if(trim(get_variable_name(dom_id, my_var_id)) == 'SCALAR_F10.7') then
   ! AETHER MODELERS SHOULD REFINE THIS AS NEEDED
   ! Set the location as per TIEGCM example for now
   ! f10_7 is most accurately located at local noon at equator.
   ! 360.0 degrees in 86400 seconds, 43200 secs == 12:00 UTC == longitude 0.0
   call get_time(state_time, seconds, days)
   longitude = 360.0_r8 * real(seconds,r8) / 86400.0_r8 - 180.0_r8
   if (longitude < 0.0_r8) longitude = longitude + 30.0_r8  
   write(string1,*)'Longitude assigned for F10.7 state variable is', longitude
   call error_handler(E_ALLMSG, 'get_state_meta_data', string1, source)
   location = set_location(longitude, 0.0_r8,  400000.0_r8, VERTISUNDEF)
   return                    
end if      
   
! Get the latitude and longitude of this columm; These are in radians
call col_index_to_lat_lon(col_index, np, del, half_del, lat, lon)

! Set the location type, lat and lon converted to degrees
location = set_location(RAD2DEG*lon, RAD2DEG*lat, &
   center_altitude(lev_index), VERTISHEIGHT)

! Set the physical quantity, e.g. QTY_TEMPERATURE
if(present(qty)) qty = my_qty

end subroutine get_state_meta_data

!------------------------------------------------------------------

! Any model specific distance calcualtion can be done here
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs

!------------------------------------------------------------------

! Any model specific distance calcualtion can be done here
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     ! location of interest
integer,                       intent(in)    :: base_type    ! observation TYPE
type(location_type),           intent(inout) :: locs(:)      ! state locations
integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
integer,                       intent(out)   :: num_close    ! how many are close
integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'

call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)

end subroutine get_close_state

!------------------------------------------------------------------

! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

end subroutine end_model

!------------------------------------------------------------------

! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if(.not. module_initialized) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!-----------------------------------------------------------------------

subroutine read_template_file()

integer               :: number_of_columns, ncid

! Gets altitudes and number of points per face row from an Aether template file
ncid = nc_open_file_readonly(trim(template_file))

! Get the number of vertical levels
ncenter_altitudes =  nc_get_dimension_size(ncid, 'z')

! Allocate space for vertical levels
allocate(center_altitude(ncenter_altitudes))

! Get the vertical levels
call nc_get_variable(ncid, 'alt', center_altitude)

! Get the number of columns
number_of_columns = nc_get_dimension_size(ncid, 'col')

call nc_close_file(ncid)

! Compute the number of grid rows across a face
np = nint(sqrt(number_of_columns / 6.0_r8))

end subroutine read_template_file

!-----------------------------------------------------------------------

function get_state_index(face, lat_ind, lon_ind, lev_ind, var_ind)

integer             :: get_state_index
integer, intent(in) :: face, lat_ind, lon_ind, lev_ind, var_ind

! Given the cube face, latitude (first) index, longitude (second) index on the face,
! the level index and the variable index, returns the state index for use
! by get_state_meta_data. Needs to know the number of lat and lon points across the face, np.

! This function makes the explicit assumption that the state is mapped to 
! the state index in the following fashion:
! From fastest to slowest varying index:
! 1. longitude index : 1 to np
! 2. latitude index  : 1 to np
! 3. face            : 0 to 5 consistent with Aether defs
! 4. level index     : 1 to n_lev
! 4. variable index  : 1 to number of variables

integer :: column

! Get the index of the column in DART storage
column = lon_ind + np * ((lat_ind - 1) + np * face)
get_state_index = get_dart_vector_index(column, lev_ind, no_third_dimension, dom_id, var_ind)

end function get_state_index

!-----------------------------------------------------------------------

! Performs IDW interpolation using great-circle distances for a triangle or quad
! Modified quad_idw_interp imported from the cice model_mod
! This should eventually be in its own module

function idw_interp(ens_size, lat, lon, y_corners, x_corners, p, num_corners)
      
integer,   intent(in) :: ens_size
real(r8)              :: idw_interp(ens_size) ! Interpolated value at (lat, lon).
real(r8),  intent(in) :: lat, lon ! Interpolation point (latitude, longitude) in degrees
real(r8),  intent(in) :: y_corners(:), x_corners(:) ! corner points (latitude, longitude) in degrees
real(r8),  intent(in) :: p(:, :) ! Values at the quadrilaterals corner points, second dimension is ens_size
integer,   intent(in) :: num_corners

! Set the power for the inverse distances
real(r8), parameter :: power = 2.0_r8 ! Power for IDW (squared distance)

! This value of epsilon radians is a distance of approximately 1 mm on the earth's surface
real(r8), parameter :: epsilon_radians = 1.56e-11_r8

type(location_type) :: corner(num_corners), point
real(r8)            :: distances(num_corners), inv_power_dist(num_corners)
                 
integer :: i, n

! Compute the distances from the point to each corner
point = set_location(lon, lat, MISSING_R8, VERTISUNDEF)
do i = 1, num_corners
   corner(i) = set_location(x_corners(i), y_corners(i), MISSING_R8, VERTISUNDEF)
   distances(i) = get_dist(point, corner(i), no_vert=.true.)
end do

if(minval(distances) < epsilon_radians) then
   ! To avoid any round off issues, if smallest distance is less than epsilon radians
   ! just assign the value at the closest gridpoint to the interpolant
   idw_interp = p(minloc(distances, 1), :)
   return
else
   ! Get the inverse distances raised to the power
   inv_power_dist = 1.0_r8 / (distances ** power)

   ! Calculate the weights for each grid point and sum up weighted values
   do n = 1, ens_size
      idw_interp(n) = sum(inv_power_dist(1:num_corners) * p(1:num_corners, n)) / sum(inv_power_dist)
   end do
endif

! Round-off can lead to result being outside of range of gridpoints
! Test for now and fix if this happens
do n = 1, ens_size
   ! If all vertices have the same value, just return that value
   ! This avoids some issues with roundoff leading to interpolated being outside of range
   if(all(p(2:num_corners, n) == p(1, n))) then
      idw_interp(n) = p(1, n)
   elseif(idw_interp(n) < minval(p(:, n)) .or. idw_interp(n) > maxval(p(:, n))) then
      write(string1,*)'IDW interpolation result is outside of range of grid point values'
      write(string2, *) 'Interpolated value, min and max are: ', &
         idw_interp(n), minval(p(:, n)), maxval(p(:, n))
      call error_handler(E_MSG, 'idw_interp', string1, source, text2=string2)
   
      ! Fixing out of range
      idw_interp(n) = max(idw_interp(n), minval(p(:, n)))
      idw_interp(n) = min(idw_interp(n), maxval(p(:, n)))
   endif

end do

end function idw_interp

!-----------------------------------------------------------------------
! interpolate in the vertical between 2 arrays of items.

! vert_fracts: 0 is 100% of the first level and 
!              1 is 100% of the second level

subroutine vert_interp(nitems, levs1, levs2, vert_fract, out_vals)

integer,  intent(in)  :: nitems
real(r8), intent(in)  :: levs1(nitems), levs2(nitems), vert_fract
real(r8), intent(out) :: out_vals(nitems)
   
out_vals(:) = levs1(:) * (1.0_r8 - vert_fract) +  levs2(:) * vert_fract
   
end subroutine vert_interp

!-----------------------------------------------------------------------

!===================================================================
! End of model_mod
!===================================================================
end module model_mod
