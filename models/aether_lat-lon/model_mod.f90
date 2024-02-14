! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

!-----------------------------------------------------------------------
!
! Interface for Aether
!
!-----------------------------------------------------------------------

use types_mod, only : &
    r8, i8, MISSING_R8, vtablenamelength

use time_manager_mod, only : &
    time_type, set_time, set_calendar_type

use location_mod, only : &
    location_type, get_close_type, &
    get_close_obs, get_close_state, &
    is_vertical, set_location, &
    VERTISHEIGHT, query_location, get_location

use utilities_mod, only : &
    open_file, close_file, &
    error_handler, E_ERR, E_MSG, E_WARN, &
    nmlfileunit, do_nml_file, do_nml_term,  &
    find_namelist_in_file, check_namelist_read, to_upper, &
    find_enclosing_indices

use obs_kind_mod, only : get_index_for_quantity

use netcdf_utilities_mod, only : &
    nc_add_global_attribute, nc_synchronize_file, &
    nc_add_global_creation_time, &
    nc_begin_define_mode, nc_end_define_mode, &
    nc_open_file_readonly, nc_get_dimension_size, nc_create_file, &
    nc_get_variable

use quad_utils_mod, only : &
    quad_interp_handle, init_quad_interp, set_quad_coords, &
    quad_lon_lat_locate, quad_lon_lat_evaluate, &
    GRID_QUAD_FULLY_REGULAR, QUAD_LOCATED_CELL_CENTERS

use state_structure_mod, only : &
    add_domain, get_dart_vector_index, get_domain_size, &
    get_model_variable_indices, get_varid_from_kind, &
    state_structure_info

use distributed_state_mod, only : get_state

use ensemble_manager_mod, only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : &
    pert_model_copies, read_model_time, write_model_time, &
    init_time => fail_init_time, &
    init_conditions => fail_init_conditions, &
    convert_vertical_obs, convert_vertical_state, adv_1step

implicit none
private

! routines required by DART code - will be called from filter and other DART executables.
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

character(len=256), parameter :: source   = 'aether_lat-lon/model_mod.f90'

logical :: module_initialized = .false.
integer :: dom_id ! used to access the state structure
type(time_type) :: assimilation_time_step 

!-----------------------------------------------------------------------
! Default values for namelist
character(len=256) :: template_file = 'filter_input_0001.nc'
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

!-----------------------------------------------------------------------
! Dimensions

character(len=4), parameter :: LEV_DIM_NAME = 'alt'
character(len=4), parameter :: LAT_DIM_NAME = 'lat'
character(len=4), parameter :: LON_DIM_NAME = 'lon'
character(len=4), parameter :: TIME_DIM_NAME = 'time'

character(len=4), parameter :: LEV_VAR_NAME = 'alt'
character(len=4), parameter :: LAT_VAR_NAME = 'lat'
character(len=4), parameter :: LON_VAR_NAME = 'lon'
character(len=4), parameter :: TIME_VAR_NAME = 'time'

! Filter
! To be assigned in assign_dimensions (for filter) 
! or get_grid_from_blocks (aether_to_dart, dart_to_aether).
real(r8), allocatable  :: levs(:), lats(:), lons(:)

integer  :: nlev, nlat, nlon
real(r8) :: lon_start, lon_delta, lat_start, lat_delta, lat_end

!-----------------------------------------------------------------------
! to be assigned in the verify_variables subroutine

type(quad_interp_handle) :: quad_interp

integer, parameter :: GENERAL_ERROR_CODE = 99
integer, parameter :: INVALID_VERT_COORD_ERROR_CODE = 15
integer, parameter :: INVALID_LATLON_VAL_ERROR_CODE = 16
integer, parameter :: INVALID_ALTITUDE_VAL_ERROR_CODE = 17
integer, parameter :: UNKNOWN_OBS_QTY_ERROR_CODE = 20

type(time_type)    :: state_time ! module-storage declaration of current model time
character(len=512) :: error_string_1, error_string_2

contains

!-----------------------------------------------------------------------
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

integer  :: iunit, io
type(var_type) :: var

module_initialized = .true.

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

call set_calendar_type('GREGORIAN')

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call assign_dimensions()

! Dimension start and deltas needed for set_quad_coords
lon_start = lons(1)
lon_delta = lons(2) - lons(1)
lat_start = lats(1)
lat_delta = lats(2) - lats(1)

var = assign_var(variables, MAX_STATE_VARIABLES)

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime 
! feel free to hardcode it and remove from the namelist.
assimilation_time_step = set_time(time_step_seconds, time_step_days)

! Define which variables are in the model state
! This is using add_domain_from_file (arg list matches)
dom_id = add_domain(template_file, var%count, var%names, var%qtys, &
                    var%clamp_values, var%updates)

call state_structure_info(dom_id)


call init_quad_interp(GRID_QUAD_FULLY_REGULAR, nlon, nlat,                    &
                      QUAD_LOCATED_CELL_CENTERS,                              &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=quad_interp)

call set_quad_coords(quad_interp, lon_start, lon_delta, lat_start, lat_delta)

end subroutine static_init_model

!-----------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = get_domain_size(dom_id)

end function get_model_size

!-----------------------------------------------------------------------
! Use quad_utils_mod to interpalate the ensemble to the ob location.

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: qty
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! Local storage

character(len=*), parameter :: routine = 'model_interpolate'

real(r8) :: loc_array(3), llon, llat, lvert, lon_fract, lat_fract
integer  :: four_lons(4), four_lats(4)
integer  :: status1, which_vert, varid
real(r8) :: quad_vals(4, ens_size)

if ( .not. module_initialized ) call static_init_model

! Assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, expected_obs will be set to the
! good values, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

expected_obs = MISSING_R8          ! the DART bad value flag
istatus = GENERAL_ERROR_CODE ! unknown error

! Get the individual locations values

loc_array  = get_location(location)
llon       = loc_array(1)
llat       = loc_array(2)
lvert      = loc_array(3)
which_vert = nint(query_location(location))

! Only height and level for vertical location type is supported at this point
if (.not. is_vertical(location, "HEIGHT") .and. .not. is_vertical(location, "LEVEL")) THEN
     istatus = INVALID_VERT_COORD_ERROR_CODE
     return
endif

! See if the state contains the obs quantity
varid = get_varid_from_kind(dom_id, qty)

if (varid > 0) then
    istatus = 0
else
    istatus = UNKNOWN_OBS_QTY_ERROR_CODE
endif

! get the indices for the 4 corners of the quad in the horizontal, plus
! the fraction across the quad for the obs location
call quad_lon_lat_locate(quad_interp, llon, llat, & 
                         four_lons, four_lats, lon_fract, lat_fract, status1)
if (status1 /= 0) then
   istatus(:) = INVALID_LATLON_VAL_ERROR_CODE  ! cannot locate enclosing horizontal quad
   return
endif

call get_quad_vals(state_handle, ens_size, varid, four_lons, four_lats, &
                   loc_array, which_vert, quad_vals, istatus)
if (any(istatus /= 0)) return

! do the horizontal interpolation for each ensemble member
call quad_lon_lat_evaluate(quad_interp, lon_fract, lat_fract, ens_size, &
                           quad_vals, expected_obs, istatus)

! All good.
istatus(:) = 0
    
end subroutine model_interpolate

!-----------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations



!-----------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional ,  intent(out) :: qty

! Local variables

integer :: lat_index, lon_index, lev_index
integer :: my_var_id, my_qty

if ( .not. module_initialized ) call static_init_model

! Restart data is ordered (lev,lat,lon) (translated from C to fortran).
call get_model_variable_indices(index_in, lev_index, lat_index, lon_index, &
                                var_id=my_var_id, kind_index=my_qty)

! should be set to the actual location using set_location()
location = set_location(lons(lon_index), lats(lat_index), levs(lev_index), VERTISHEIGHT)

! should be set to the physical quantity, e.g. QTY_TEMPERATURE
if (present(qty)) qty = my_qty

end subroutine get_state_meta_data

!-----------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

end subroutine end_model

!-----------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

character(len=*), parameter :: routine = 'nc_write_model_atts'

if ( .not. module_initialized ) call static_init_model

! It is already in define mode from nc_create_file.
! OR NOT, if called by create_and_open_state_output
call nc_begin_define_mode(ncid)

! nc_write_model_atts is called by create_and_open_state_output,
!   which calls nf90_enddef before it.
call nc_add_global_creation_time(ncid, routine)

call nc_add_global_attribute(ncid, "model_source", source, routine)
call nc_add_global_attribute(ncid, "model", "aether", routine)

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
! Read dimension information from the template file and use 
! it to assign values to variables.

subroutine assign_dimensions()

integer  :: ncid
character(len=24), parameter :: routine = 'assign_dimensions'

call error_handler(E_MSG, routine, 'reading filter input ['//trim(template_file)//']')

ncid = nc_open_file_readonly(template_file, routine)

! levels
nlev = nc_get_dimension_size(ncid, trim(LEV_DIM_NAME), routine)
allocate(levs(nlev))
call nc_get_variable(ncid, trim(LEV_VAR_NAME), levs, routine)

! latitiude
nlat = nc_get_dimension_size(ncid, trim(LAT_DIM_NAME), routine)
allocate(lats(nlat))
call nc_get_variable(ncid, trim(LAT_VAR_NAME), lats, routine)

! longitude
nlon = nc_get_dimension_size(ncid, trim(LON_DIM_NAME), routine)
allocate(lons(nlon))
call nc_get_variable(ncid, trim(LON_VAR_NAME), lons, routine)

end subroutine assign_dimensions

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
! Extract state values needed by the interpolation from all ensemble members.

subroutine get_quad_vals(state_handle, ens_size, varid, four_lons, four_lats, &
                         lon_lat_vert, which_vert, quad_vals, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: varid
integer,             intent(in)  :: four_lons(4), four_lats(4)
real(r8),            intent(in)  :: lon_lat_vert(3)
integer,             intent(in)  :: which_vert
real(r8),            intent(out) :: quad_vals(4, ens_size)
integer,             intent(out) :: istatus(ens_size)

integer  :: lev1, lev2, stat
real(r8) :: vert_val, vert_fract
character(len=512) :: error_string_1

character(len=*), parameter :: routine = 'get_quad_vals'

quad_vals(:,:) = MISSING_R8
istatus(:) = GENERAL_ERROR_CODE

vert_val = lon_lat_vert(3)

if ( which_vert == VERTISHEIGHT ) then
    call find_enclosing_indices(nlev, levs(:), vert_val, lev1, lev2, &
                                vert_fract, stat, log_scale = .false.)

    if (stat /= 0) then
        istatus     = INVALID_ALTITUDE_VAL_ERROR_CODE
    end if
else
    istatus(:) = INVALID_VERT_COORD_ERROR_CODE
    write(error_string_1, *) 'unsupported vertical type: ', which_vert
    call error_handler(E_ERR, routine, error_string_1, source)
endif

! we have all the indices and fractions we could ever want.
! now get the data values at the bottom levels, the top levels, 
! and do vertical interpolation to get the 4 values in the columns.
! the final horizontal interpolation will happen later.

if (varid > 0) then

    call get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                               lev1, lev2, vert_fract, varid, quad_vals, istatus)
else 
    write(error_string_1, *) 'unsupported variable: ', varid
    call error_handler(E_ERR, routine, error_string_1, source)
endif

if (any(istatus /= 0)) return

! when you get here, istatus() was set either by passing it to a
! subroutine, or setting it explicitly here.
end subroutine get_quad_vals

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
! Extract the state values at the corners of the 2 quads used for interpolation.

subroutine get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                 lev1, lev2, vert_fract, varid, quad_vals, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_lons(4), four_lats(4)
integer,             intent(in) :: lev1, lev2
real(r8),            intent(in) :: vert_fract
integer,             intent(in) :: varid
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer     :: icorner
integer(i8) :: state_indx
real(r8)    :: vals1(ens_size), vals2(ens_size)
real(r8)    :: qvals(ens_size)

character(len=*), parameter :: routine = 'get_four_state_values:'

do icorner = 1, 4

   ! Most rapidly varying dim must be first
   state_indx = get_dart_vector_index(lev1, four_lats(icorner), &
                                      four_lons(icorner), dom_id, varid)

   if (state_indx < 0) then
      write(error_string_1,'(A)') 'Could not find dart state index from '
      write(error_string_2,'(A,3F15.4)') 'lon, lat, and lev1 index :', &
           four_lons(icorner), four_lats(icorner), lev1
      call error_handler(E_ERR, routine, error_string_1, source, &
                         text2=error_string_2)
      return
   endif

   vals1(:) = get_state(state_indx, state_handle)    ! all the ensemble members for level (i)

   state_indx = get_dart_vector_index(lev2, four_lats(icorner), &
                                      four_lons(icorner), dom_id, varid)

   if (state_indx < 0) then
      write(error_string_1,'(A)') 'Could not find dart state index from '
      write(error_string_2,'(A,3F15.4)') 'lon, lat, and lev2 index :', &
           four_lons(icorner), four_lats(icorner), lev2
      call error_handler(E_ERR, routine, error_string_1, source, &
                         text2=error_string_2)
      return
   endif

   vals2(:) = get_state(state_indx, state_handle)    ! all the ensemble members for level (i)

   ! if directly using quad_vals here, it would create a temporary array and give a warning
   call vert_interp(ens_size, vals1, vals2, vert_fract, qvals)
   quad_vals(icorner, :) = qvals
enddo

istatus = 0

end subroutine get_four_state_values

!-----------------------------------------------------------------------
! End of model_mod
!-----------------------------------------------------------------------
end module model_mod

