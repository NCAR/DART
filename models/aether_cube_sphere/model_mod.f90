! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

use           netcdf

use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength, DEG2RAD, RAD2DEG, PI

use time_manager_mod, only : time_type, set_time

use     location_mod, only : location_type, get_close_type, get_dist, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, query_location, &
                             get_location, VERTISHEIGHT, VERTISUNDEF

use    utilities_mod, only : register_module, error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, to_upper, &
                             find_enclosing_indices

use obs_kind_mod,     only : get_index_for_quantity

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, nc_open_file_readonly, & 
                                 nc_close_file
use distributed_state_mod, only : get_state
use state_structure_mod, only : add_domain, get_dart_vector_index, get_domain_size, &
                                get_model_variable_indices, get_varid_from_kind

use transform_state_mod, only : file_type

use ensemble_manager_mod, only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, read_model_time, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step

implicit none
private

! routines required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          adv_1step,              &
          init_time,              &
          init_conditions,        &
          shortest_time_between_assimilations, &
          write_model_time

! Routine for comprehensive test of interpolation
public :: test_grid_box

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! message strings                
character(len=512) :: string1
character(len=512) :: string2

integer  :: iunit, io

logical :: module_initialized = .false.
integer :: dom_id ! used to access the state structure
type(time_type) :: assimilation_time_step

! Parameter controlling tolerance for ???
real(r8), parameter :: roundoff = 1.0e-12_r8

! Geometry variables that are used throughout the module; read from a template file
integer               :: np                ! Number of grid rows across a face
real(r8)              :: del, half_del     ! Grid row spacing and half of that
integer               :: ncenter_altitudes ! The number of altitudes and the altitudes
real(r8), allocatable :: center_altitude(:)

! Horizontal column dimension rather than being direct functions of latitude and longitude.
integer                          :: no_third_dimension = -99

! Error codes
integer, parameter :: INVALID_VERT_COORD_ERROR_CODE = 15
integer, parameter :: INVALID_ALTITUDE_VAL_ERROR_CODE = 17
integer, parameter :: UNKNOWN_OBS_QTY_ERROR_CODE = 20

! Use the namelist for options to be set at runtime.
character(len=256) :: template_file = 'model_restart.nc'
integer  :: time_step_days      = 0
integer  :: time_step_seconds   = 3600

integer, parameter              :: MAX_STATE_VARIABLES     = 100
integer, parameter              :: NUM_STATE_TABLE_COLUMNS = 5
character(len=vtablenamelength) :: variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ''

type :: var_type
    integer :: count
    character(len=64), allocatable :: names(:)
    integer,           allocatable :: qtys(:)
    real(r8),          allocatable :: clamp_values(:, :)
    logical,           allocatable :: updates(:)
end type var_type

namelist /model_nml/ template_file, time_step_days, time_step_seconds, variables

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. 

subroutine static_init_model()

type(var_type) :: var
real(r8)       :: cube_side

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source)

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
assimilation_time_step = set_time(time_step_seconds, &
                                  time_step_days)

var = assign_var(variables, MAX_STATE_VARIABLES)

! Get the altitudes and the number of grid rows
call read_template_file()

! Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
! two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2.0_r8 * sqrt(1.0_r8 / 3.0_r8)
! These grid spacings are in module storage since they are used repeatedly in many routines
del = cube_side / np
half_del = del / 2.0_r8

! Define which variables are in the model state
dom_id = add_domain(template_file, var%count, var%names, var%qtys, &
                    var%clamp_values, var%updates)

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
!

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

! Location values stored in a vector
real(r8), dimension(3)                           :: lon_lat_alt

! Vertical interpolation variables
integer     :: below_index, above_index, enclosing_status, which_vertical
real(r8)    :: fraction

real(r8)    :: grid_pt_lat(4), grid_pt_lon(4), pt_lat, pt_lon, bounding_value(4, 2, ens_size)
real(r8)    :: below_values(ens_size), above_values(ens_size)
integer(i8) :: bounding_state_index(4, 2)
integer     :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4), num_bound_points
integer     :: var_id, n_lev, i

write(*, *) 'ENTERING MODEL_INTERPOLATE'

! Initialize module if not already done
if ( .not. module_initialized ) call static_init_model

! Set all obs to MISSING_R8 initially
expected_obs(:) = MISSING_R8
! Successful default
istatus(:) = 0

! Determine the vertical location type
lon_lat_alt = get_location(location)
pt_lat = lon_lat_alt(2) * DEG2RAD
pt_lon = lon_lat_alt(1) * DEG2RAD
which_vertical = nint(query_location(location, 'WHICH_VERT'))

! Only height currently supported for observations; fail if other is selected
if (.not. which_vertical == VERTISHEIGHT ) then
   istatus(:) = INVALID_VERT_COORD_ERROR_CODE
   return
endif

! See if the state contains the obs quantity
var_id = get_varid_from_kind(dom_id, qty)
write(*, *) 'GETTING VARID ', var_id, qty
if(var_id <= 0) then
   istatus = UNKNOWN_OBS_QTY_ERROR_CODE
   return
endif

! Find the bounding vertical levels and the fractional distance between
call find_enclosing_indices(ncenter_altitudes, center_altitude, lon_lat_alt(3), below_index, &
   above_index, fraction, enclosing_status)
if (enclosing_status /= 0) then
   istatus(:) = INVALID_ALTITUDE_VAL_ERROR_CODE
   return
endif

! If the vertical location is acceptable, then do the horizontal interpolation
! Find the enclosing triangle or quad
call get_bounding_box(pt_lat, pt_lon, &
   grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

! Map the grid_face, latitude index and longitude index to the one dimensional index used in the state vector
! Then get the state values
do i = 1, num_bound_points
   bounding_state_index(i, 1) =  get_state_index(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
      below_index, var_id)
   bounding_value(i, 1, :) = get_state(bounding_state_index(i, 1), state_handle)
   bounding_state_index(i, 2) =  get_state_index(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
      above_index, var_id)
   bounding_value(i, 2, :) = get_state(bounding_state_index(i, 2), state_handle)
enddo

! Do inverse distance weighted horizontal interpolation on both levels
below_values =  idw_interp(ens_size, RAD2DEG*pt_lat, &
   RAD2DEG*pt_lon, RAD2DEG*grid_pt_lat, RAD2DEG*grid_pt_lon, bounding_value(:, 1, :), num_bound_points)
above_values =  idw_interp(ens_size, RAD2DEG*pt_lat, &
   RAD2DEG*pt_lon, RAD2DEG*grid_pt_lat, RAD2DEG*grid_pt_lon, bounding_value(:, 2, :), num_bound_points)

! Do the vertical interpolation, linear in height to get final
call vert_interp(ens_size, below_values, above_values, fraction, expected_obs)

write(*, *) 'EXITING MODEL_INTERPOLATE ', expected_obs

end subroutine model_interpolate


!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty

! Local variables

integer :: lev_index, col_index, my_var_id, my_qty

real(r8) :: lat, lon

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, col_index, lev_index, no_third_dimension, &
                                var_id=my_var_id, kind_index=my_qty)

! Get the latitude and longitude of this columm; These lats and lons are in radians
call col_index_to_lat_lon(col_index, lat, lon)

! Set the location type
location = set_location(RAD2DEG*lon, RAD2DEG*lat, center_altitude(lev_index), VERTISHEIGHT)

! should be set to the physical quantity, e.g. QTY_TEMPERATURE
if (present(qty)) qty = my_qty

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

if ( .not. module_initialized ) call static_init_model

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
! Parse the table of variables characteristics into arrays for easier access.

function assign_var(variables, MAX_STATE_VARIABLES) result(var)

character(len=vtablenamelength), intent(in) :: variables(:, :)
integer, intent(in)                         :: MAX_STATE_VARIABLES

type(var_type) :: var
integer        :: ivar
character(len=vtablenamelength) :: table_entry

!-----------------------------------------------------------------------
! Codes for interpreting the NUM_STATE_TABLE_COLUMNS of the variables table
integer, parameter :: NAME_INDEX      = 1 ! ... variable name
integer, parameter :: QTY_INDEX       = 2 ! ... DART qty
integer, parameter :: MIN_VAL_INDEX   = 3 ! ... minimum value if any
integer, parameter :: MAX_VAL_INDEX   = 4 ! ... maximum value if any
integer, parameter :: UPDATE_INDEX    = 5 ! ... update (state) or not

! Loop through the variables array to get the actual count of the number of variables
do ivar = 1, MAX_STATE_VARIABLES
   ! If the element is an empty string, the loop has exceeded the extent of the variables
   if (variables(1, ivar) == '') then
      var%count = ivar-1
      exit
   endif 
enddo

! Allocate the arrays in the var derived type
allocate(var%names(var%count), var%qtys(var%count), var%clamp_values(var%count, 2), var%updates(var%count))

do ivar = 1, var%count
   
   var%names(ivar) = trim(variables(NAME_INDEX, ivar))

   table_entry = variables(QTY_INDEX, ivar)
   call to_upper(table_entry)

   var%qtys(ivar) = get_index_for_quantity(table_entry)

   if (variables(MIN_VAL_INDEX, ivar) /= 'NA') then
      read(variables(MIN_VAL_INDEX, ivar), '(d16.8)') var%clamp_values(ivar,1)
   else
      var%clamp_values(ivar,1) = MISSING_R8
   endif

   if (variables(MAX_VAL_INDEX, ivar) /= 'NA') then
      read(variables(MAX_VAL_INDEX, ivar), '(d16.8)') var%clamp_values(ivar,2)
   else
      var%clamp_values(ivar,2) = MISSING_R8
   endif

   table_entry = variables(UPDATE_INDEX, ivar)
   call to_upper(table_entry)

   if (table_entry == 'UPDATE') then
      var%updates(ivar) = .true.
   else
      var%updates(ivar) = .false.
   endif

enddo
   
end function assign_var

!-----------------------------------------------------------------------

subroutine read_template_file()

integer               :: dimid, varid, number_of_columns
character(len=256)    :: name
type(file_type)       :: templatefile

! Getting the altitudes and the number of points per face row from
! This should be getting this from one of the retart files; is that info available at the model_mod level?
! Working with template file for now which has to be named in the model_mod_nml
templatefile%file_path = trim(template_file)
templatefile%ncid = nc_open_file_readonly(templatefile%file_path)

! Get the number of vertical levels
templatefile%ncstatus = nf90_inq_dimid(templatefile%ncid, 'z', dimid)
templatefile%ncstatus = nf90_inquire_dimension(templatefile%ncid, dimid, name, ncenter_altitudes)

! Allocate space for vertical levels
allocate(center_altitude(ncenter_altitudes))

! Get the vertical levels
templatefile%ncstatus = nf90_inq_varid(templatefile%ncid, 'alt', varid)
templatefile%ncstatus = nf90_get_var(templatefile%ncid, varid, center_altitude)

! Get the number of columns
templatefile%ncstatus = nf90_inq_dimid(templatefile%ncid, 'col', dimid)
templatefile%ncstatus = nf90_inquire_dimension(templatefile%ncid, dimid, name, number_of_columns)

call nc_close_file(templatefile%ncid)

! Compute the number of grid rows across a face
np = nint(sqrt(number_of_columns / 6.0_r8))

end subroutine read_template_file



!-----------------------------------------------------------------------
! Determines if the projection of a point p onto the plane of a triangle with vertices
! v1, v2 and v3 is inside the triangle or not. Computes the areas of each of the triangles
! between p and a pair of vertices. These should sum to the area of the triangle if p
! is inside and be larger than that if p is outside. 

function is_point_in_triangle(v1, v2, v3, p)

logical              :: is_point_in_triangle
real(r8), intent(in) :: v1(3), v2(3), v3(3), p(3)    

real(r8) :: a(3), b(3), perp(3), unit_perp(3), p_proj(3)
real(r8) :: offset, len_s1, len_s2, len_s3, len_p1, len_p2, len_p3, at, at1, at2, at3
real(r8) :: area_dif, dif_frac, threshold

! Get the projection of the point p onto the plane containing the triangle
! Start by getting perpendicular vector to plane by cross product
a = v1 - v2
b = v2 - v3
perp(1) = a(2) * b(3) - a(3) * b(2)
perp(2) = a(3) * b(1) - a(1) * b(3)
perp(3) = a(1) * b(2) - a(2) * b(1)
! Get unit vector in direction of perp
unit_perp = perp / sqrt(dot_product(perp, perp))
! Projection of vector from v1 to p on the unit perp vector is how much to move to get to plane
offset = dot_product((p-v1), unit_perp)
p_proj = p - offset * unit_perp

! Compute lengths of the sides
len_s1 = sqrt(dot_product(v1-v2, v1-v2))
len_s2 = sqrt(dot_product(v3-v2, v3-v2))
len_s3 = sqrt(dot_product(v1-v3, v1-v3))

! Compute the lengths from the point p
len_p1 = sqrt(dot_product(p_proj-v1, p_proj-v1))
len_p2 = sqrt(dot_product(p_proj-v2, p_proj-v2))
len_p3 = sqrt(dot_product(p_proj-v3, p_proj-v3))

! Area of triangle
at = heron(len_s1, len_s2, len_s3)

! Compute areas of sub triangles
at1 = heron(len_p1, len_p2, len_s1)
at2 = heron(len_p2, len_p3, len_s2)
at3 = heron(len_p3, len_p1, len_s3)

! Difference between sub triangles and the triangle area
area_dif = at1 + at2 + at3 - at

! Quadrilaterals on the interior of the cube sphere sides are really spherical quads,
! There sides are great circles. This routine assumes that the triangles composing the quads
! have straight sides in regular space. The algorithm finds points that are inside the 
! spherical quads. These quads actually 'bulge' out compared to the regular sides, so it is possible
! to have points that are inside the spherical quad but just barely outside of the regular
! quads. This threshold is tuned so that these points still show as inside. The tuning is for
! np = 18 (number of points along a grid face is 18). Fewer points might require a larger
! threshold while more points might be okay with a smaller one.
threshold = 0.002_r8

dif_frac = area_dif / at
is_point_in_triangle = abs(dif_frac) < threshold

end function is_point_in_triangle
!-----------------------------------------------------------------------

function lat_lon_to_xyz(lat, lon)

real(r8)             :: lat_lon_to_xyz(3)
real(r8), intent(in) :: lat, lon

lat_lon_to_xyz(1) = cos(lat) * cos(lon)
lat_lon_to_xyz(2) = cos(lat) * sin(lon)
lat_lon_to_xyz(3) = sin(lat)

end function lat_lon_to_xyz

!-----------------------------------------------------------------------

! Given the face (from 0 to 5) the number of lons/lats on a face (np) and
! the indices of the i and j grid point, returns the latitude and longitude of the point

subroutine grid_to_lat_lon(face, lat_ind, lon_ind, lat, lon)

integer,  intent(in)  :: face, lat_ind, lon_ind
real(r8), intent(out) :: lat, lon

real(r8) :: x, y, blon, blat, rot_angle
real(r8) :: vect(3), RZ(3, 3), rot_vect(3)


x = sqrt(1.0_r8/3.0_r8) - (half_del + del * (lon_ind - 1))
if(face == 5) then 
   y =  sqrt(1.0_r8/3.0_r8) - (half_del + del * (lat_ind - 1))
else
   y = -sqrt(1.0_r8/3.0_r8) + (half_del + del * (lat_ind - 1))
endif

if(face < 4) then
   blon = atan2(sqrt(1.0_r8 / 3.0_r8), x)
   blat = atan2(y, sqrt(1.0_r8/3.0_r8 + x**2.0_r8))
   blon = blon - PI/4.0_r8

   ! Above is for face 0; add PI/2 for each additional face tangent to equator
   lon = blon + PI/2.0_r8 * face; 
   lat = blat;
elseif(face == 4 .or. face == 5) then
   ! Face 4 is tangent to south pole
   lon = atan2(y, x)
   lat = atan2(sqrt(1.0_r8/3.0_r8), sqrt(x**2.0_r8 + y**2.0_r8))

   if(face == 4) lat = -lat

   !  Get ready for rotation
   vect = lat_lon_to_xyz(lat, lon)

   ! Then rotate 45 degrees around Z
   rot_angle = -PI/4.0_r8;

   ! Create the rotation matrix
   RZ(1, 1:3) = [cos(rot_angle),  sin(rot_angle), 0.0_r8]
   RZ(2, 1:3) = [-sin(rot_angle), cos(rot_angle), 0.0_r8]
   RZ(3, 1:3) = [0.0_r8,          0.0_r8,         1.0_r8]
   rot_vect = matmul(RZ, vect)

   lat = asin(rot_vect(3))
   lon = atan2(rot_vect(2), rot_vect(1))
   ! Note that there are inconsistent treatments of the value near longitude
   ! 0 in the grid files for Aether. Some points have a value near or just less
   ! than 2PI, other points have values just greater than 0. This code 
   ! avoids values near to 2PI and have 0 instead.
   if(lon < 0.0_r8) lon =  lon + 2.0_r8*PI
   if(lon >= 2.0_r8*PI) lon = 0.0_r8

endif

end subroutine grid_to_lat_lon

!-----------------------------------------------------------------------

subroutine lat_lon_to_grid(lat, lon, face, lat_ind, lon_ind)

real(r8), intent(in)  :: lat, lon
integer,  intent(out) :: face, lat_ind, lon_ind

real(r8) :: len(2)

! Get the face and the length along the two imbedded cube faces for the point
call get_face(lat, lon, face, len);

! Figure out which interval this is in along each cube face; This gives 0 to np grid indices
lon_ind = nint((len(1) + half_del) / del)
lat_ind = nint((len(2) + half_del) / del)

end subroutine lat_lon_to_grid

!-----------------------------------------------------------------------

subroutine fix_face(face, lat_grid, lon_grid, f_face, f_lat_grid, f_lon_grid, edge, corner)

integer, intent(in)  :: face, lat_grid, lon_grid
integer, intent(out) :: f_face, f_lat_grid, f_lon_grid
logical, intent(out) :: edge, corner

integer :: left_neighbor(6),   right_neighbor(6)
integer :: bottom_neighbor(6), top_neighbor(6)
integer :: left_lon_grid(6),   right_lon_grid(6)
integer :: left_lat_grid(6),   right_lat_grid(6)
integer :: bottom_lon_grid(6), top_lon_grid(6)
integer :: bottom_lat_grid(6), top_lat_grid(6)
! For points past the edge of a face, finds the corresponding points on the adjacent face
! Need to do something more special for corner points

! Default is not a corner or an edge
corner = .false.
edge   = .false.

! Just return if no edge
if(lon_grid > 0 .and. lon_grid < np + 1 .and. lat_grid > 0 .and. lat_grid < np + 1) then
   f_face = face
   f_lon_grid = lon_grid
   f_lat_grid = lat_grid
   return
endif

if((lat_grid == 0 .or. lat_grid == np + 1) .and. (lon_grid == 0 .or. lon_grid == np + 1)) then
   corner = .true.
   f_face     = -99
   f_lon_grid = -99
   f_lat_grid = -99
   return
endif

! Otherwise, on an edge
edge = .true.
! Deal with each side of faces separately
if(lon_grid == 0) then
   ! On left edge not corner
   left_neighbor = [3, 0, 1, 2, 0, 0]
   f_face = left_neighbor(face + 1)
   left_lon_grid = [np, np, np, np, lat_grid, np+1-lat_grid]
   f_lon_grid = left_lon_grid(face + 1)
   left_lat_grid = [lat_grid, lat_grid, lat_grid, lat_grid, 1, np]
   f_lat_grid = left_lat_grid(face + 1)
elseif(lon_grid == np + 1) then
   ! On right edge not corner
   right_neighbor = [1, 2, 3, 0, 2, 2]
   f_face = right_neighbor(face + 1)
   right_lon_grid = [1, 1, 1, 1, np+1-lat_grid, lat_grid]
   f_lon_grid = right_lon_grid(face + 1)
   right_lat_grid = [lat_grid, lat_grid, lat_grid, lat_grid, 1, np]
   f_lat_grid = right_lat_grid(face + 1)
elseif(lat_grid == 0) then
   ! On bottom edge not corner
   bottom_neighbor = [4, 4, 4, 4, 3, 1]
   f_face = bottom_neighbor(face + 1)
   bottom_lon_grid = [1, lon_grid, np, np+1-lon_grid, np+1-lon_grid, lon_grid]
   f_lon_grid = bottom_lon_grid(face + 1)
   bottom_lat_grid = [lon_grid, np, np+1-lon_grid, 1, 1, np]
   f_lat_grid = bottom_lat_grid(face + 1)
elseif(lat_grid == np + 1) then
   ! On top edge not corner
   top_neighbor = [5, 5, 5, 5, 1, 3]
   f_face = top_neighbor(face + 1)
   top_lon_grid = [1, lon_grid, np, np+1-lon_grid, lon_grid, np+1-lon_grid]
   f_lon_grid = top_lon_grid(face + 1)
   top_lat_grid = [np+1-lon_grid, 1, lon_grid, np, 1, np]
   f_lat_grid = top_lat_grid(face + 1)
endif

end subroutine fix_face

!-----------------------------------------------------------------------

subroutine get_face(lat, lon_in, face, len)

real(r8), intent(in)  :: lat, lon_in
integer,  intent(out) :: face
real(r8), intent(out) :: len(2)

integer :: side, rside, rside2
real(r8) :: inv_sqrt_2, rlon, rlon2, gama, gamb, lon
real(r8) :: vec(3), rot_vec(3), rot_vec2(3), rot(3, 3), rot2(3, 3), lon_grid(2), lon_grid_m(2)
! Returns which face contains (lat, lon) and the length from the edge of the point
! along each of the great circle axes.

! Range adjustment
lon = lon_in
if(lon >= 2.0_r8*PI) lon = 0.0_r8

! Convert lat lon to x y z on unit sphere
vec = lat_lon_to_xyz(lat, lon)

! Get the longitudes for this point in the two rotated spaces

!====================================================================
! Following code shows individual rotations;
! Single matrix at the end multiplies them together off-line for single rotation
! Can collapse these to a single rotation vector for efficiency
! Rotation 90 degrees around y to put pole on equator
!RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
! Rotate 45 degrees around x
!RX = [1 0 0; 0 cosd(45) sind(45); 0 -sind(45) cosd(45)];
! Then rotate 45 degrees around Z
!RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
! Get longitude in the two rotated spaces
!rot_vec = RZ * RX * RY * vec;
!====================================================================
inv_sqrt_2 = 1.0_r8 / sqrt(2.0_r8);
rot(1, 1:3) = [0.5_r8,     0.5_r8,      -inv_sqrt_2]
rot(2, 1:3) = [0.5_r8,     0.5_r8,      inv_sqrt_2]
rot(3, 1:3) = [inv_sqrt_2, -inv_sqrt_2, 0.0_r8]
rot_vec = matmul(rot, vec)
! Compute the longitude in the rotated space
rlon = atan2(rot_vec(2), rot_vec(1))
if(rlon < 0.0_r8) rlon = rlon + 2.0_r8*PI

!====================================================================
! Can collapse these to a single rotation vector for efficiency
! Rotation 90 degrees around y to put pole on equator
!RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
! Rotate -45 degrees around x
!RX2 = [1 0 0; 0 cosd(-45) sind(-45); 0 -sind(-45) cosd(-45)];
! Then rotate 45 degrees around Z
!RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
! Get longitude in the two rotated spaces
!rot_vec2 = RZ * RX2 * RY * vec;
!====================================================================
rot2(1, 1:3) = [-0.5_r8,     0.5_r8,        -inv_sqrt_2]
rot2(2, 1:3) = [-0.5_r8,     0.5_r8,        inv_sqrt_2]
rot2(3, 1:3) = [inv_sqrt_2,  inv_sqrt_2,    0.0_r8]
rot_vec2 = matmul(rot2, vec)
! Compute the longitude in the rotated space
rlon2 = atan2(rot_vec2(2), rot_vec2(1))
if(rlon2 < 0.0_r8) rlon2 = rlon2 + 2.0_r8*PI

! Which non-polar side could we be on, 1 to 4
side = floor(lon / (PI/2.0_r8)) + 1.0_r8
! Which rotated 1 side are we on 
rside = floor(rlon / (PI/2.0_r8)) + 1.0_r8
! Which rotated 2 side
rside2 = floor(rlon2 / (PI/2.0_r8)) + 1.0_r8

! Figure out the face from here (0 to 5, 4 is south, 5 is north)
! These are consistent with the numbering on Aether grid files for the cubed sphere
if    ( side == 1 .and. rside  == 1) then
   face = 0; lon_grid(1) = lon;  lon_grid(2) = rlon
elseif( side == 2 .and. rside2 == 1) then
   face = 1; lon_grid(1) = lon;  lon_grid(2) = rlon2
elseif( side == 3 .and. rside  == 3) then 
   face = 2; lon_grid(1) = lon;  lon_grid(2) = rlon
elseif( side == 4 .and. rside2 == 3) then
   face = 3; lon_grid(1) = lon;  lon_grid(2) = rlon2
elseif(rside == 4 .and. rside2 == 4) then
   face = 4; lon_grid(1) = rlon; lon_grid(2) = rlon2
elseif(rside == 2 .and. rside2 == 2) then
   face = 5; lon_grid(1) = rlon; lon_grid(2) = rlon2
endif

! Can also use the fact that the projection is equidistant on the imbedded cube to get what fraction 
! across the imbedded rectangle we are
! Take the longitudes and turn them into a number between -sqrt(1/3) and sqrt(1/3)
lon_grid_m = mod(lon_grid, PI/2.0_r8)

! Use law of sines to go from lon back to position along edge of imbedded cube 
! The triangle of interest has a side of length 2sqrt(1/3) (1/2 of the planar diagonal of the imbedded cube)
! The angles adjacent to this side are the longitude and 45 degrees
! The angle opposite the side of length 2sqrt(1/3) is PI - (longitude + PI/4)
! The side opposite the longitude is how far along the side of the cube
! The cube side is 2sqrt(1/3), so the length along the side is between zero and this value
gama = PI - (PI/4.0_r8 + lon_grid_m(1))
len(1) = sqrt(2.0_r8/3.0_r8) * sin(lon_grid_m(1)) / sin(gama)

gamb = PI - (PI/4 + lon_grid_m(2))
len(2) = sqrt(2.0_r8/3.0_r8) * sin(lon_grid_m(2)) / sin(gamb)

! If we are on sides 2 or 3, the lengths need to be modified because the grid storage
! for Aether goes from smallest latitude to largest and the longitudes of the shifted
! poles are going the opposite way
if(face == 2 .or. face == 3) len(2) = 2.0_r8 * sqrt(1.0_r8/3.0_r8) - len(2)

! Same for face 4 (the bottom) but it's the other coordinate that's reversed
if(face == 4) len(1) = 2.0_r8*sqrt(1.0_r8/3.0_r8) - len(1)

end subroutine get_face

!-----------------------------------------------------------------------
subroutine get_corners(face, lat_grid, lon_grid, lat, lon, &
   f_face, f_lat_grid, f_lon_grid, num_bound_points)
   
integer, intent(in) :: face, lat_grid, lon_grid
real(r8), intent(in) :: lat, lon
integer, intent(out) :: f_face(4), f_lat_grid(4), f_lon_grid(4), num_bound_points

integer :: corner, quad, i
integer :: quad_lon_grid(3, 4), quad_lat_grid(3, 4), quad_face(3, 4)
real(r8) :: pxyz(3), qxyz(4, 3), grid_pt_lat, grid_pt_lon

! Checks to see if the point under consideration is at a corner
! If it is, return the face, lat_index, and lon_index for each of the three bounding points

! Default is to find a triangle
num_bound_points = 3

if(face == 0) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 1
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 2
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 5
   else
      corner = 6;
   endif
elseif(face == 1) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 2
   elseif(lat_grid == 0    .and. lon_grid == np+1) then 
      corner = 3
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 6
   else
      corner = 7
   endif
elseif(face == 2) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 3
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 4
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 7
   else
      corner = 8
   endif
elseif(face == 3) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 4
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 1
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 8
   else
      corner = 5;
   endif
elseif(face == 4) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 1
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 4
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 2
   else
      corner = 3
   endif
elseif(face == 5) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 6
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 7
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 5
   else
      corner = 8;
   endif
endif

! Harvest the information on the grid points bounding the appropriate corner
! Arrays of info for adjacent quads for bulges (three of them, first index)
quad_lon_grid(1:3, 1:4)  = -99
quad_lat_grid(1:3, 1:4)  = -99
quad_face(1:3, 1:4) = -99

if(corner == 1) then
   f_face(1:3) =     [3,  0, 4]
   f_lon_grid(1:3) = [np, 1, 1]
   f_lat_grid(1:3) = [1,  1, 1]
   quad_face(1, 1:4) = [3,  0, 0, 3]
   quad_face(2, 1:4) = [0,  0, 4, 4]
   quad_face(3, 1:4) = [3,  3, 4, 4]
   quad_lat_grid(1, 1:4) = [1,  1, 2, 2]
   quad_lat_grid(2, 1:4) = [1,  1, 1, 2]
   quad_lat_grid(3, 1:4) = [1,  1, 1, 1]
   quad_lon_grid(1, 1:4) = [np,   1,  1, np]
   quad_lon_grid(2, 1:4) = [1,    2,  1, 1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1, 2 ]
elseif(corner == 2) then
   f_face(1:3) =     [0,  1, 4 ]
   f_lon_grid(1:3) = [np, 1, 1 ]
   f_lat_grid(1:3) = [1,  1, np]
   quad_face(1, 1:4) = [0, 1, 1, 0]
   quad_face(2, 1:4) = [1, 1, 4, 4]
   quad_face(3, 1:4) = [0, 0, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2,  2   ]
   quad_lat_grid(2, 1:4) = [1, 1, np, np  ]
   quad_lat_grid(3, 1:4) = [1, 1, np, np-1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  2,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  1 ]
elseif(corner == 3) then
   f_face(1:3) =     [1,  2, 4 ]
   f_lon_grid(1:3) = [np, 1, np]
   f_lat_grid(1:3) = [1,  1, np]
   quad_face(1, 1:4) = [1, 2, 2, 1]
   quad_face(2, 1:4) = [2, 2, 4, 4]
   quad_face(3, 1:4) = [1, 1, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2,    2 ]
   quad_lat_grid(2, 1:4) = [1, 1, np-1, np]
   quad_lat_grid(3, 1:4) = [1, 1, np,   np]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np  ]
   quad_lon_grid(2, 1:4) = [1,    2,  np, np  ]
   quad_lon_grid(3, 1:4) = [np-1, np, np, np-1]
elseif(corner == 4) then
   f_face(1:3) =     [2,  3, 4 ]
   f_lon_grid(1:3) = [np, 1, np]
   f_lat_grid(1:3) = [1,  1, 1 ]
   quad_face(1, 1:4) = [2, 3, 3,2]
   quad_face(2, 1:4) = [3, 3, 4, 4]
   quad_face(3, 1:4) = [2, 2, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2, 2]
   quad_lat_grid(2, 1:4) = [1, 1, 1, 1]
   quad_lat_grid(3, 1:4) = [1, 1, 1, 2]
   quad_lon_grid(1, 1:4) = [np,   1,  1,    np]
   quad_lon_grid(2, 1:4) = [1,    2,  np-1, np]
   quad_lon_grid(3, 1:4) = [np-1, np, np,   np]
elseif(corner == 5) then
   f_face(1:3) =     [3,  0,  5 ]
   f_lon_grid(1:3) = [np, 1,  1 ]
   f_lat_grid(1:3) = [np, np, np]
   quad_face(1, 1:4) = [3, 0, 0, 3]
   quad_face(2, 1:4) = [0, 0, 5, 5]
   quad_face(3, 1:4) = [3, 3, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np  ]
   quad_lat_grid(2, 1:4) = [np,   np,   np, np-1]
   quad_lat_grid(3, 1:4) = [np,   np,   np, np  ]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  1,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  2 ]
elseif(corner == 6) then
   f_face(1:3) =     [0,  1, 5 ]
   f_lon_grid(1:3) = [np, 1,  1]
   f_lat_grid(1:3) = [np, np, 1]
   quad_face(1, 1:4) = [0, 1, 1, 0]
   quad_face(2, 1:4) = [1, 1, 5, 5]
   quad_face(3, 1:4) = [0, 0, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np]
   quad_lat_grid(2, 1:4) = [np,   np,   1,  1 ]
   quad_lat_grid(3, 1:4) = [np,   np,   1,  2 ]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  2,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  1 ]
elseif(corner == 7) then
   f_face(1:3) =     [1,  2,  5 ]
   f_lon_grid(1:3) = [np, 1,  np]
   f_lat_grid(1:3) = [np, np, 1 ]
   quad_face(1, 1:4) = [1, 2, 2, 1]
   quad_face(2, 1:4) = [2, 2, 5, 5]
   quad_face(3, 1:4) = [1, 1, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np]
   quad_lat_grid(2, 1:4) = [np,   np,   1,   2]
   quad_lat_grid(3, 1:4) = [np,   np,   1,   1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np  ]
   quad_lon_grid(2, 1:4) = [1,    2,  np, np  ]
   quad_lon_grid(3, 1:4) = [np-1, np, np, np-1]
elseif(corner == 8) then
   f_face(1:3) =     [2,  3,  5 ]
   f_lon_grid(1:3) = [np, 1,  np]
   f_lat_grid(1:3) = [np, np, np]
   quad_face(1, 1:4) = [2, 3, 3, 2]
   quad_face(2, 1:4) = [3, 3, 5, 5]
   quad_face(3, 1:4) = [2, 2, 5, 5];
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np  ]
   quad_lat_grid(2, 1:4) = [np,   np,   np, np  ]
   quad_lat_grid(3, 1:4) = [np,   np,   np, np-1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,    np]
   quad_lon_grid(2, 1:4) = [1,    2,  np-1, np]
   quad_lon_grid(3, 1:4) = [np-1, np, np,   np]
endif

! Load up the array for the point
pxyz = lat_lon_to_xyz(lat, lon)

! Get lats and lons of the triangle vertices
do i = 1, 3
   call grid_to_lat_lon(f_face(i), f_lat_grid(i), f_lon_grid(i), grid_pt_lat, grid_pt_lon)
   ! Convert to x, y, z coords to check for whether points are in tris/quads
   qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat, grid_pt_lon)
enddo

! See if the point is in the triangle; if so, all is good
if(is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz)) return

! If it's not in the triangle, have to check the three adjacent quads at the corner
num_bound_points = 4

do quad = 1, 3
   ! Compute lat and lon for a quad
   do i = 1, 4
         call grid_to_lat_lon(quad_face(quad, i), quad_lat_grid(quad, i), quad_lon_grid(quad, i), &
            grid_pt_lat, grid_pt_lon)
      ! Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat, grid_pt_lon)
   enddo

   ! See if the point is inside this quad
   if(is_point_in_quad(qxyz, pxyz)) then
      f_face = quad_face(quad, 1:4)
      f_lat_grid = quad_lat_grid(quad, 1:4)
      f_lon_grid = quad_lon_grid(quad, 1:4)
      return
   endif
enddo

! Falling of the end is not happy
!!!fprintf('UNEXPECTED FAILURE IN GET_CORNERS.M\n');
!!!stop

end subroutine get_corners

!-----------------------------------------------------------------------

! Given the latitude and longitude of a point, returns the face, array indices, latitude
! and longitude of the bounding three or four grid points along the number of points
! (3 triangle; 4 quad). np is the number of grid points across each face of the cube sphere.

subroutine get_bounding_box(lat, lon, &
   grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

real(r8), intent(in)  :: lat, lon
real(r8), intent(out) :: grid_pt_lat(4), grid_pt_lon(4)
integer,  intent(out) :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4), num_bound_points

real(r8) :: len(2), qxyz(4, 3), pxyz(3)
integer  :: face, low_grid(2), hi_grid(2), i, my_pt, corner_index
integer  :: lat_grid(4), lon_grid(4), face1_pts(2), face2_pts(2), face1_count, face2_count
logical  :: on_edge, edge, corner

! Get the face and the length along the two imbedded cube faces for the point
call get_face(lat, lon, face, len);

! Figure out which interval this is in along each cube face; This gives 0 to np grid indices
low_grid(1) = floor((len(1) + half_del) / del)
low_grid(2) = floor((len(2) + half_del) / del)
hi_grid = low_grid + 1

! Get the indices for the lat and lon directions: Points go counterclockwise starting from lower left 
! For now assume this is a quad, but will correct below if it is a triangle
lat_grid(1) = low_grid(2); lat_grid(2) = hi_grid(2); lat_grid(3) = lat_grid(2); lat_grid(4) = lat_grid(1)
lon_grid(1) = low_grid(1); lon_grid(2) = lon_grid(1); lon_grid(3) = hi_grid(1); lon_grid(4) = lon_grid(3)

! If points are on the edge map to adjacent faces
on_edge = .false.
do i = 1, 4
   call fix_face(face, lat_grid(i), lon_grid(i), &
      grid_face(i), grid_lat_ind(i), grid_lon_ind(i), edge, corner)
   ! If any point is on an edge, on_edge is true
   if(edge) on_edge = .true.
   if(corner) then
      corner_index = i
      exit
   endif
enddo

! If it's at a corner, need to find the triangles in a different fashion
! It is possible that the point initially looks like it is in a corner due to the fact 
! that the edges of the grid are on great circles from the corresponding faces, but the
! grid point at the edge are not connected by these great circles.
if(corner) then
   call get_corners(face, lat_grid(corner_index), lon_grid(corner_index), lat, lon, &
      grid_face, grid_lat_ind, grid_lon_ind, num_bound_points)
   if(num_bound_points == 4) corner = .false.
else
   ! If not initially at a corner it's definitely in a quad
   num_bound_points = 4
endif

! Compute the lat and lon corresponding to these point
do i = 1, num_bound_points
   call grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
      grid_pt_lat(i), grid_pt_lon(i))
enddo

! Make on_edge true only if we are on an edge but not at a corner
on_edge = (on_edge .and.  .not. corner)

if(on_edge) then
   ! If this is an edge, may need to revise box selection
   ! See if the point is in the box (approximately)
   ! Load up the arrays for the vertex points
   do i = 1, num_bound_points   
      ! Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat(i), grid_pt_lon(i))
   enddo

   ! Convert point to xyz
   pxyz = lat_lon_to_xyz(lat, lon);

   if(.not. is_point_in_quad(qxyz, pxyz)) then
      ! Not in this box, need to move 'equatorward'
      ! Find indices (from 1 to 4) of points on the same face
      face1_pts(1:2) = 0; face2_pts(1:2) = 0;
      face1_count = 0;    face2_count = 0;
      do i = 1, 4
         if(grid_face(i) == grid_face(1)) then
            face1_count = face1_count + 1;    face1_pts(face1_count) = i
         else
            face2_count = face2_count + 1;    face2_pts(face2_count) = i
         endif
      enddo

      ! First process points of the first face
      ! Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face1_pts(1)) == grid_lon_ind(face1_pts(2))) then
         ! Adjust the face1 latitudes
         do i = 1, 2
            my_pt = face1_pts(i)
            if(grid_lat_ind(my_pt) > np/2) then
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1
            endif
         enddo
      else
         ! Adjust the face1 longitudes
         do i = 1, 2
            my_pt = face1_pts(i)
            if(grid_lon_ind(my_pt) > np/2) then
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1
            endif
         enddo
      endif

      ! Do the same thing for face2 
      ! Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face2_pts(1)) == grid_lon_ind(face2_pts(2))) then
         ! Adjust the face1 latitudes
         do i = 1, 2
            my_pt = face2_pts(i);
            if(grid_lat_ind(my_pt) > np/2) then
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1
            endif
         enddo
      else
         ! Adjust the face1 longitudes
         do i = 1, 2
            my_pt = face2_pts(i);
            if(grid_lon_ind(my_pt) > np/2) then
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1
            endif
         enddo
      endif

      ! Compute the lat and lon corresponding to these point
      do i = 1, num_bound_points
         call grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
            grid_pt_lat(i), grid_pt_lon(i))
      enddo
   endif

endif

end subroutine get_bounding_box

!-----------------------------------------------------------------------

subroutine test_grid_box

integer  :: i, j, num_bound_points
integer  :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4)
real(r8) :: pt_lon_d, pt_lat_d, pt_lon, pt_lat 
real(r8) :: qxyz(4, 3), pxyz(3), grid_pt_lat(4), grid_pt_lon(4)
logical  :: inside

type(location_type) :: location
integer             :: qty, lon_count, lat_count, my_face, my_level, my_qty, my_lon_ind, my_lat_ind
integer(i8)         :: state_index, state_index2
real(r8)            :: lon_lat_hgt(3), my_lat, my_lon, base_dist, dist_sum
integer             :: test_face, test_lat_ind, test_lon_ind, col_index, test_col_index


! Temporary test that grid_to_lat_lon and lat_lon_to_grid are inverses of each other
do my_face = 0, 5
   do my_lat_ind = 1, np
      do my_lon_ind = 1, np
         call grid_to_lat_lon(my_face, my_lat_ind, my_lon_ind, pt_lat, pt_lon)
         call lat_lon_to_grid(pt_lat, pt_lon, test_face, test_lat_ind, test_lon_ind)
         if(my_face .ne. test_face .or. my_lat_ind .ne. test_lat_ind .or. my_lon_ind .ne. test_lon_ind) then
            write(*, *) 'lat_lon_to_grid is not inverse of grid_to_lat_lon'
            write(*, *) my_face, test_face, my_lat_ind, test_lat_ind, my_lon_ind, test_lon_ind
            stop
         endif

         col_index = my_lon_ind + (my_lat_ind - 1) * np + my_face * np*np
         call col_index_to_lat_lon(col_index, pt_lat, pt_lon)
         test_col_index = lat_lon_to_col_index(pt_lat, pt_lon)
write(*, *) col_index, test_col_index
         if(col_index .ne. test_col_index) then
            write(*, *) 'lat_lon_to_col_index is not inverse of col_index_to_lat_lon'
            write(*, *) my_face, my_lat_ind, my_lon_ind, col_index, test_col_index
            stop
         endif

      enddo
   enddo
enddo

! Test points for the following:
! 1. Does the bounding box found contain the observed point?
! 2. Are the computed vertex latitudes and longitudes the same as those in the Aether grid files?

! Largest edges are on the quads in the center of a face
! Get distance along side of center quad
do i = 1, 2
   ! Traverse half of the rows (each np across), plus halfway across the next row
   state_index = (np/2)*np + np/2 + i - 1
   call get_state_meta_data(state_index, location, qty)
   lon_lat_hgt = get_location(location)       
   qxyz(i, 1:3) = lat_lon_to_xyz(DEG2RAD*lon_lat_hgt(2), DEG2RAD*lon_lat_hgt(1))
enddo
base_dist = sqrt(sum((qxyz(1, :) - qxyz(2, :))**2))

do lon_count = 0, 3600
   pt_lon_d = lon_count / 10.0_r8
   write(*, *) pt_lon_d
   do lat_count = -900, 900
   pt_lat_d = lat_count / 10.0_r8


! Convert to radians
pt_lon = DEG2RAD * pt_lon_d
pt_lat = DEG2RAD * pt_lat_d

! Get the x, y, z coords for this point
pxyz = lat_lon_to_xyz(pt_lat, pt_lon);

call get_bounding_box(pt_lat, pt_lon, &
   grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

do i = 1, num_bound_points
   ! Convert to x, y, z coords to check for whether points are in tris/quads
   qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat(i), grid_pt_lon(i));
enddo


! Get the latitude and longitude of the bounding points from get_state_meta_data as a confirmation test
do i = 1, num_bound_points
   state_index = get_state_index(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
      1, 1)
   call get_state_meta_data(state_index, location, qty)
   lon_lat_hgt = get_location(location)       
   ! Deal with Aether file round off
   if(abs(lon_lat_hgt(1) - 360.0_r8) < 0.0001) lon_lat_hgt(1) = 0.0_r8 
   if(abs(RAD2DEG*grid_pt_lat(i) - lon_lat_hgt(2)) > 0.0001_r8 .or. &
      abs(RAD2DEG*grid_pt_lon(i) - lon_lat_hgt(1)) > 0.0001_r8) then
      write(*, *) grid_pt_lat(i), grid_pt_lon(i), lon_lat_hgt(2), lon_lat_hgt(1)
      stop
   endif
enddo

if(num_bound_points == 3) then
   ! See if the point is inside a local approximately tangent triangle
   inside = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz);
else
   ! Or quadrilateral
   inside = is_point_in_quad(qxyz, pxyz);
endif

if(.not. inside) then
   write(*, *) 'ERROR: inside is false'
   write(*, *)pt_lat, pt_lon, num_bound_points
   stop
endif

! Also check on distance to vertices; this greatly reduces the possibility that
! bounding boxes that are bigger than they should be are being found
dist_sum = 0.0_r8
do i = 1, num_bound_points
   ! Compute sum of distances between point and each of the vertices
   dist_sum = dist_sum + sqrt(sum((qxyz(i, :) - pxyz)**2))
end do

if(num_bound_points == 4) then
   ! For quad, sum should be less than 3.5 times the baseline
   if(dist_sum / base_dist > 3.5_r8) then
      write(*, *) 'ratio of sum of distances to vertices is too large'
      do i = 1, num_bound_points
         write(*, *) 'grid ', i, grid_pt_lat(i), grid_pt_lon(i)
         write(*, *) 'grid xyz ', i, qxyz(i, :) 
      enddo
      write(*, *) 'point ', pt_lat, pt_lon
      write(*, *) 'point xyz ', pxyz
      stop
   endif
elseif(num_bound_points == 3) then
   ! For quad, sum should be less than 3 times the baseline
   if(dist_sum / base_dist > 3.0_r8) then
      write(*, *) 'ratio of sum of distances to vertices is too large'
      stop
   endif
endif

enddo
enddo

!-------------------
! Block that loops through all state variables and confirms that the algorithms for mapping
! state vector (face/lon/lat) indices and get_state_meta_data correcty match up.
do my_qty = 1, 2
   do my_level = 1, 44
      do my_face = 0, 5
         do my_lat_ind = 1, np
            do my_lon_ind = 1, np
               state_index = get_state_index(my_face, my_lat_ind, my_lon_ind, &
                  my_level, my_qty)
               call get_state_meta_data(state_index, location, qty)
               lon_lat_hgt = get_location(location)       

               ! Want to compare the lat lon directly from code to that from get_state_meta_data
               call grid_to_lat_lon(my_face, my_lat_ind, my_lon_ind, my_lat, my_lon)
                
               ! ROUNDOFF FROM AETHER, NEED TO FIX somewhere.
               if(abs(lon_lat_hgt(1) - 360.0_r8) < 0.0001) lon_lat_hgt(1) = 0.0_r8 

               ! Check that things are consistent            
               if(abs(RAD2DEG*my_lat - lon_lat_hgt(2)) > 0.0001_r8 .or. &
                  abs(RAD2DEG*my_lon - lon_lat_hgt(1)) > 0.0001_r8) then
                  write(*, *) 'ERROR: grid points not appropriately mapping'
                  write(*, *) my_face, my_qty, my_level, my_lat_ind, my_lon_ind
                  stop
               endif
   
            enddo
         enddo
      enddo
   enddo
enddo
!-------------------------------

end subroutine test_grid_box

!-----------------------------------------------------------------------

function is_point_in_quad(v, p)

logical              :: is_point_in_quad
real(r8), intent(in) :: v(4, 3), p(3)

logical :: inside_t(4)

integer :: i

! See if the point is inside this quad; it's inside if it's in one or more contained triangles
inside_t(1) = is_point_in_triangle(v(1, :), v(2, :), v(3, :), p)
inside_t(2) = is_point_in_triangle(v(1, :), v(2, :), v(4, :), p)
inside_t(3) = is_point_in_triangle(v(1, :), v(3, :), v(4, :), p)
inside_t(4) = is_point_in_triangle(v(2, :), v(3, :), v(4, :), p)

is_point_in_quad = any(inside_t)

end function is_point_in_quad
   

!-----------------------------------------------------------------------

! Computes Herons formula to get area of triangle from lenghts of sides
! Super accuracy is not needed in the area calculation here

function heron(a, b, c)

real(r8)             :: heron
real(r8), intent(in) :: a, b, c

real(r8) :: s, arg

s = (a + b + c) /2
arg = (s * (s - a) * (s - b) * (s - c))

! Make sure we don't roundoff to a negative
if(arg <= 0.0_r8) then
   heron = 0.0_r8
else
   heron = sqrt(arg)
endif

end function heron

!-----------------------------------------------------------------------

function get_state_index(face, lat_ind, lon_ind, lev_ind, var_ind)

integer             :: get_state_index
integer, intent(in) :: face, lat_ind, lon_ind, lev_ind, var_ind

! Given the cube face, latitude (first) index, longitude (second) index on the face,
! the level index and the variable index, returns the state index for use
! by get_state_meta_data. Needs to know (are these in global storage) the number
! of lat and lon points across the face (np)

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
!!!write(*, *) 'in get_state_index calling get_dart_vector_index'
!!!write(*, *) 'column, lev_index, no_third, var_ind', column, lev_ind, no_third_dimension, var_ind
get_state_index = get_dart_vector_index(column, lev_ind, no_third_dimension, dom_id, var_ind)

end function get_state_index

!-----------------------------------------------------------------------

subroutine col_index_to_lat_lon(col_index, lat, lon)

integer,  intent(in)  :: col_index
real(r8), intent(out) :: lat, lon

integer :: face, resid, lat_ind, lon_ind

! Given the index of a horizontal column, returns the latitude and longitude in radians

! Which face are we on? np**2 points per face
face = (col_index - 1) / (np**2)
resid = col_index - face * np**2

! Get latitude index
lat_ind = (resid - 1) / np + 1

lon_ind = resid - (lat_ind - 1) * np

! Get the corresponding latitude and longitude
call grid_to_lat_lon(face, lat_ind, lon_ind, lat, lon)

end subroutine col_index_to_lat_lon

!-----------------------------------------------------------------------

function lat_lon_to_col_index(lat, lon)

integer              :: lat_lon_to_col_index
real(r8), intent(in) :: lat, lon

integer :: face, lat_ind, lon_ind

! Get the face, lat_ind and lon_ind
call lat_lon_to_grid(lat, lon, face, lat_ind, lon_ind)

! Get column index using first level, first variable
lat_lon_to_col_index = get_state_index(face, lat_ind, lon_ind, lev_ind = 1, var_ind = 1)

end function lat_lon_to_col_index

!-----------------------------------------------------------------------

! Temporary use of modified quad_idw_interp imported from the cice model_mod
! This may/should be in its own module

function idw_interp(ens_size, lat, lon, y_corners, x_corners, p, num_corners)

! Performs IDW interpolation using great-circle distances for a triangle or quad
      
real(r8)              :: idw_interp(ens_size) ! Interpolated value at (lon, lat).
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: lon, lat ! Interpolation point (longitude, latitude) in degrees
real(r8),  intent(in) :: y_corners(:), x_corners(:) ! corner points (latitude, longitude) in degrees
real(r8),  intent(in) :: p(:, :) ! Values at the quadrilaterals corner points, second dimension is ens_size
integer,   intent(in) :: num_corners

! Set the power for the inverse distances
real(r8), parameter :: power = 2.0_r8 ! Power for IDW (squared distance)

! This value of epsilon radians is a distance of approximately 1 mm
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

write(*, *) 'point ', lat, lon
do i = 1, num_corners
   write(*, *) 'bound ', i, y_corners(i), x_corners(i)
enddo
write(*, *) 'distances ', distances
 
if(minval(distances) < epsilon_radians) then
   ! To avoid any round off issues, if smallest distance is less than epsilon radians
   ! just assign the value at the closest gridpoint to the interpolant
   write(*, *) 'in the close distance case'
   idw_interp = p(minloc(distances, 1), :)
else
   ! Get the inverse distances raised to the power
   inv_power_dist = 1.0_r8 / (distances ** power)

   ! Calculate the weights for each grid point and sum up weighted values
   do n = 1, ens_size
      idw_interp(n) = sum(inv_power_dist(1:num_corners) * p(1:num_corners, n)) / sum(inv_power_dist)
   end do
endif

do n = 1, ens_size
   write(*, *) 'idw_interp ', n, idw_interp(n), p(1:num_corners, n)
enddo


! Unclear if round-off could ever lead to result being outside of range of gridpoints
! Test for now and terminate if this happens 
do n = 1, ens_size
   if(idw_interp(n) < minval(p(:, n)) .or. idw_interp(n) > maxval(p(:, n))) then
      write(string1,*)'IDW interpolation result is outside of range of grid point values'
write(*, *) 'n ', n, idw_interp(n), minval(p(:, n)), maxval(p(:, n))
      write(string2, *) 'Interpolated value, min and max are: ', &
         idw_interp(n), minval(p(:, n)), maxval(p(:, n))
write(*, *) 'string2 ', string2
      call error_handler(E_ERR, 'idw_interp', string1, &
         source, revision, revdate, text2=string2)
   endif

   ! Fixing out of range; this will not happen with current error check 
   idw_interp(n) = max(idw_interp(n), minval(p(:, n)))
   idw_interp(n) = min(idw_interp(n), maxval(p(:, n)))
end do

end function idw_interp


!-----------------------------------------------------------------------
! interpolate in the vertical between 2 arrays of items.

! vert_fracts: 0 is 100% of the first level and 
!              1 is 100% of the second level

subroutine vert_interp(nitems, levs1, levs2, vert_fract, out_vals)

   integer,  intent(in)  :: nitems
   real(r8), intent(in)  :: levs1(nitems)
   real(r8), intent(in)  :: levs2(nitems)
   real(r8), intent(in)  :: vert_fract
   real(r8), intent(out) :: out_vals(nitems)
   
   out_vals(:) = (levs1(:) * (1.0_r8 - vert_fract)) + &
                 (levs2(:) *           vert_fract )
   
end subroutine vert_interp

!-----------------------------------------------------------------------

!===================================================================
! End of model_mod
!===================================================================
end module model_mod
