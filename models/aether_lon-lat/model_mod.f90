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
    r4, r8, i8, MISSING_R4, MISSING_R8, vtablenamelength, MISSING_I, RAD2DEG

use time_manager_mod, only : &
    time_type, set_calendar_type, set_time, get_time, set_date, &
    print_date, print_time

use location_mod, only : &
    location_type, get_close_type, &
    loc_get_close_obs => get_close_obs, &
    loc_get_close_state => get_close_state, &
    is_vertical, set_location, &
    VERTISHEIGHT, query_location, get_location

use utilities_mod, only : &
    open_file, close_file, file_exist, register_module, &
    error_handler, E_ERR, E_MSG, E_WARN, &
    nmlfileunit, do_output, do_nml_file, do_nml_term,  &
    find_namelist_in_file, check_namelist_read, to_upper, &
    find_enclosing_indices

use obs_kind_mod, only : QTY_GEOMETRIC_HEIGHT

use netcdf_utilities_mod, only : &
    nc_add_global_attribute, nc_synchronize_file, &
    nc_add_global_creation_time, &
    nc_begin_define_mode, nc_end_define_mode, &
    nc_open_file_readonly, nc_get_dimension_size, nc_create_file, &
    nc_close_file, nc_get_variable, nc_define_dimension, &
    nc_define_real_variable, nc_define_real_scalar, nc_open_file_readwrite, &
    nc_add_attribute_to_variable, nc_put_variable, &
    nc_get_attribute_from_variable, NF90_FILL_REAL

use quad_utils_mod, only : &
    quad_interp_handle, init_quad_interp, set_quad_coords, &
    quad_lon_lat_locate, quad_lon_lat_evaluate, &
    GRID_QUAD_FULLY_REGULAR, QUAD_LOCATED_CELL_CENTERS

use obs_kind_mod, only : get_index_for_quantity

use state_structure_mod, only : &
    add_domain, get_dart_vector_index, get_domain_size, &
    get_model_variable_indices, get_varid_from_kind

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
! TODO: Is nc_write_model_vars no longer mandatory?
!       Tiegcm has it listed, but it's just a pass-through to-from default_model_mod
!       which has a do-nothing version, and a note "currently unused".
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

public :: restart_files_to_netcdf, &
          netcdf_to_restart_files, &
          block_file_name

character(len=256), parameter :: source   = 'aether_lon-lat/model_mod.f90'
character(len=32 ), parameter :: revision = ''
character(len=32 ), parameter :: revdate  = ''

logical :: module_initialized = .false.
integer :: dom_id ! used to access the state structure
type(time_type) :: assimilation_time_step 

!-----------------------------------------------------------------------
! Default values for namelist
! TODO: replace model_nml:filter_io_filename with filter_io_root, 
!            so that namelist doesn't need to be changed for each member
character(len=256) :: filter_io_filename = 'filter_input_0001.nc'
integer  :: time_step_days      = 0
integer  :: time_step_seconds   = 3600
integer  :: debug = 0

! TODO: Should this be defined here, or does it come from netcdf_utilities_mod.f90?
!     It's a public parameter from that module, which gets it from the netcdf module
!     https://docs.unidata.ucar.edu/netcdf-fortran/current/f90-variables.html#f90-variables-introduction
!     integer, parameter :: NF90_MAX_NAME = 256
!     This module uses vtablenamelength instead (which is shorter = less white space output 
!     to diagnostics).
integer, parameter              :: MAX_STATE_VARIABLES     = 100
integer, parameter              :: NUM_STATE_TABLE_COLUMNS = 6
character(len=vtablenamelength) :: variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ' ' 

namelist /model_nml/ filter_io_filename, time_step_days, time_step_seconds, debug, variables

!-----------------------------------------------------------------------
! aether_to_dart namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: aether_restart_dirname    = 'none'
! An ensemble of file names is created using this root and $member in it,
character (len = vtablenamelength) :: filter_io_root = 'filter_input'

namelist /aether_to_dart_nml/ aether_restart_dirname, filter_io_root, variables, debug

! dart_to_aether namelist parameters with default values.
!-----------------------------------------------------------------------

namelist /dart_to_aether_nml/ aether_restart_dirname, filter_io_root, variables, debug

!-----------------------------------------------------------------------
! Codes for interpreting the NUM_STATE_TABLE_COLUMNS of the variables table
! VT_ORIGININDX is used differently from the usual domains context.
! It does not provide full path+filenames.  Here it is used by aether_to_dart and
! dart_to_aether to identify whether a variable comes from a neutrals or ions block file.
! It is not used by filter's definition of a domain.
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_ORIGININDX   = 5 ! file of origin
integer, parameter :: VT_STATEINDX    = 6 ! ... update (state) or not

!-----------------------------------------------------------------------
! Dimensions

! TODO: using length * causes(?) a problem when calling nc_define_var_real_Nd
! with the list of dim_names in this order.  nc_define also uses size * 
! and apparently looks at the first one, sees that it's size 3, and assumes that for all.
!   routine: nc_define_var_real_Nd
!   message:  "Temperature" inquire dimension id for dim "tim": 
!   errcode      -46= NetCDF: Invalid dimension ID or name
character(len=4), parameter :: LEV_DIM_NAME = 'alt'
character(len=4), parameter :: LAT_DIM_NAME = 'lat'
character(len=4), parameter :: LON_DIM_NAME = 'lon'
character(len=4), parameter :: TIME_DIM_NAME = 'time'

character(len=4), parameter :: LEV_VAR_NAME = 'alt'
character(len=4), parameter :: LAT_VAR_NAME = 'lat'
character(len=4), parameter :: LON_VAR_NAME = 'lon'
character(len=4), parameter :: TIME_VAR_NAME = 'time'

! Aether
! number of blocks along each dim
integer :: nblocks_lon=MISSING_I, nblocks_lat=MISSING_I, nblocks_lev=MISSING_I
integer :: nx_per_block, ny_per_block, nz_per_block

! TODO: should nghost be read from the namelist?
!       Aaron; not in the foreseeable future.
integer, parameter :: nghost = 2   ! number of ghost cells on all edges

! Filter
! To be assigned in assign_dimensions (for filter) 
! or get_grid_from_blocks (aether_to_dart, dart_to_aether).
real(r8), allocatable  :: levs(:), lats(:), lons(:)
! TODO: Sort out the precision of levs... in filter_*.nc versus Aether restarts.
! Can't just change this to r4.
! I'll need to read the dims from filter_input_0001.nc into r4 temp array,
! then convert to these r8 vars.

integer  :: nlev, nlat, nlon
real(r8) :: lon_start, lon_delta, lat_start, lat_delta, lat_end

!-----------------------------------------------------------------------
! Day 0 in Aether's calendar is (+/1 a day) -4710/11/24 0 UTC
! integer               :: aether_ref_day = 2451545.0_r8  ! cJULIAN2000 in Aether = day of date 2000/01/01.
character(len=32)     :: calendar = 'GREGORIAN'

! But what we care about is the ref time for the times in the files, which is 1965-1-1 00:00
integer, dimension(:) :: aether_ref_date(5) = (/1965,1,1,0,0/)  ! y,mo,d,h,m (secs assumed 0)
type(time_type)       :: aether_ref_time
integer               :: aether_ref_ndays, aether_ref_nsecs

!-----------------------------------------------------------------------
! to be assigned in the verify_variables subroutine
integer  :: nvar, nvar_neutral, nvar_ion

character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES)
real(r8)                        :: var_ranges(MAX_STATE_VARIABLES,2)
logical                         :: var_update(MAX_STATE_VARIABLES)
integer                         :: var_qtys(MAX_STATE_VARIABLES)

type(quad_interp_handle) :: quad_interp

integer, parameter :: GENERAL_ERROR_CODE = 99
integer, parameter :: INVALID_VERT_COORD_ERROR_CODE = 15
integer, parameter :: INVALID_LATLON_VAL_ERROR_CODE = 16
integer, parameter :: INVALID_ALTITUDE_VAL_ERROR_CODE = 17
integer, parameter :: UNKNOWN_OBS_QTY_ERROR_CODE = 20

type(time_type) :: state_time ! module-storage declaration of current model time

character(len=512) :: error_string_1, error_string_2

contains

!-----------------------------------------------------------------------
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

integer  :: iunit, io

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source)

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type(calendar)

! Debug global att creation time
! This filter_io_filename comes from the namelist (filter_input_0001.nc)
! Somehow filter is creating 'filter_output_0001.nc' when it dies.
call assign_dimensions(filter_io_filename, levs, lats, lons, nlev, nlat, nlon)

! Dimension start and deltas needed for set_quad_coords
lon_start = lons(1)
lon_delta = lons(2) - lons(1)
lat_start = lats(1)
lat_delta = lats(2) - lats(1)

call verify_variables(variables, filter_io_filename, nvar, var_names, var_qtys, var_ranges, var_update)

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime 
! feel free to hardcode it and remove from the namelist.
assimilation_time_step = set_time(time_step_seconds, time_step_days)

! Define which variables are in the model state
! This is using add_domain_from_file (arg list matches)
dom_id = add_domain(filter_io_filename, nvar, var_names, var_qtys, var_ranges, var_update)

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

IF (debug > 85) then
   write(error_string_1,'(A,3F15.4)')  'requesting interpolation at ', llon, llat, lvert
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
end if

! Only height and level for vertical location type is supported at this point
if (.not. is_vertical(location, "HEIGHT") .and. .not. is_vertical(location, "LEVEL")) THEN
     istatus = INVALID_VERT_COORD_ERROR_CODE
     return
endif

if (qty == QTY_GEOMETRIC_HEIGHT .and. is_vertical(location, "LEVEL")) then
   if (nint(lvert) < 1 .or. nint(lvert) > size(levs,1)) then
      expected_obs = MISSING_R8
      istatus = 1
   else
      expected_obs = levs(nint(lvert))
      istatus = 0
   endif
   return ! Early Return
endif

! do we know how to interpolate this quantity?
call ok_to_interpolate(qty, varid, status1)

if (status1 /= 0) then
   if(debug > 12) then
      write(error_string_1,'(A,I5,A)') 'Did not find observation quantity ', qty, &
           ' in the state vector'
      call error_handler(E_WARN, routine, error_string_1, source, revision, revdate)
   endif
   istatus(:) = status1   ! this quantity not in the state vector
   return
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
! character(len=*), parameter :: routine = 'get_state_meta_data'

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

! character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                       num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs

!-----------------------------------------------------------------------
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

! character(len=*), parameter :: routine = 'get_close_state'


call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_ind, dist, ens_handle)

end subroutine get_close_state

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

! Debug global att creation time; This requires being in define mode.
! nc_write_model_atts is called by create_and_open_state_output,
!   which calls nf90_enddef before it.
call nc_add_global_creation_time(ncid, routine)

call nc_add_global_attribute(ncid, "model_source", source, routine)
call nc_add_global_attribute(ncid, "model", "aether", routine)
! TODO Shouldn't the calendar type be defined here?
!      It's defined in the time variable = good enough for write_model_time.

! call nc_end_define_mode(ncid)

end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
! Add dimension variable contents to the file.

subroutine def_fill_dimvars(ncid)

integer, intent(in) :: ncid

character(len=*), parameter :: routine = 'def_fill_dimvars'

! call nc_begin_define_mode(ncid)

! Global atts for aether_to_dart and dart_to_aether.
call nc_add_global_creation_time(ncid, routine)
call nc_add_global_attribute(ncid, "model_source", source, routine)
call nc_add_global_attribute(ncid, "model", "aether", routine)

! define grid dimensions
call nc_define_dimension(ncid, trim(LEV_DIM_NAME),  nlev,  routine)
call nc_define_dimension(ncid, trim(LAT_DIM_NAME),  nlat,  routine)
call nc_define_dimension(ncid, trim(LON_DIM_NAME),  nlon,  routine)

! define grid variables
! z
call nc_define_real_variable(     ncid, trim(LEV_VAR_NAME), (/ trim(LEV_DIM_NAME) /), routine)
call nc_add_attribute_to_variable(ncid, trim(LEV_VAR_NAME), 'units',     'm', routine)
call nc_add_attribute_to_variable(ncid, trim(LEV_VAR_NAME), 'long_name', 'height above mean sea level', routine)

! latitude
call nc_define_real_variable(     ncid, trim(LAT_VAR_NAME), (/ trim(LAT_DIM_NAME) /),  routine)
call nc_add_attribute_to_variable(ncid, trim(LAT_VAR_NAME), 'units',     'degrees_north', routine)
call nc_add_attribute_to_variable(ncid, trim(LAT_VAR_NAME), 'long_name', 'latitude', routine)

! longitude
call nc_define_real_variable(     ncid, trim(LON_VAR_NAME), (/ trim(LON_VAR_NAME) /), routine)
call nc_add_attribute_to_variable(ncid, trim(LON_VAR_NAME), 'units', 'degrees_east', routine)
call nc_add_attribute_to_variable(ncid, trim(LON_VAR_NAME), 'long_name', 'longitude',  routine)

! Dimension 'time' will no longer be created by write_model_time,
! or by nc_define_unlimited_dimension.  It will be a scalar variable.
! time
call nc_define_real_scalar(       ncid, trim(TIME_VAR_NAME), routine)
call nc_add_attribute_to_variable(ncid, trim(TIME_VAR_NAME), 'calendar', 'gregorian', routine)
call nc_add_attribute_to_variable(ncid, trim(TIME_VAR_NAME), 'units', 'days since 1601-01-01 00:00:00', routine)
call nc_add_attribute_to_variable(ncid, trim(TIME_VAR_NAME), 'long_name', 'gregorian_days',  routine)

call nc_end_define_mode(ncid)

call nc_put_variable(ncid, trim(LEV_VAR_NAME),  levs,  routine)
call nc_put_variable(ncid, trim(LAT_VAR_NAME),  lats,  routine)
call nc_put_variable(ncid, trim(LON_VAR_NAME),  lons,  routine)
! time will be written elsewhere.

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine def_fill_dimvars

!-----------------------------------------------------------------------
! Read dimension information from the template file and use 
! it to assign values to variables.

subroutine assign_dimensions(filter_io_filename, levs, lats, lons, nlev, nlat, nlon)

character(len=*),      intent(in)  :: filter_io_filename
real(r8), allocatable, intent(out) :: levs(:), lats(:), lons(:)
integer,               intent(out) :: nlev, nlat, nlon

integer  :: ncid
character(len=24), parameter :: routine = 'assign_dimensions'

call error_handler(E_MSG, routine, 'reading filter input ['//trim(filter_io_filename)//']')

ncid = nc_open_file_readonly(filter_io_filename, routine)

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
! Parse the table of variables' characteristics into arrays for easier access.

subroutine verify_variables(variables, file, nvar, &
                            var_names, var_qtys, var_ranges, var_update)

character(len=*), intent(in)    :: variables(:,:)
character(len=*), intent(inout) :: file
integer,          intent(out)   :: nvar
character(len=*), intent(out)   :: var_names(:)
real(r8),         intent(out)   :: var_ranges(:,:)
logical,          intent(out)   :: var_update(:)
integer,          intent(out)   :: var_qtys(:)

character(len=*), parameter :: routine = 'verify_variables'

integer  :: io, i, quantity
real(r8) :: minvalue, maxvalue

character(len=vtablenamelength) :: varname
character(len=vtablenamelength) :: dartstr
character(len=vtablenamelength) :: minvalstring
character(len=vtablenamelength) :: maxvalstring
character(len=vtablenamelength) :: state_or_aux

nvar = 0
MY_LOOP : do i = 1, size(variables,2)

! TODO Why define these intermediate strings?  Is the code clearer or faster?
   varname      = variables(VT_VARNAMEINDX,i)
   dartstr      = variables(VT_KINDINDX,i)
   minvalstring = variables(VT_MINVALINDX,i)
   maxvalstring = variables(VT_MAXVALINDX,i)
   state_or_aux = variables(VT_STATEINDX,i)

   if ( varname == ' ' .and. dartstr == ' ' ) exit MY_LOOP ! Found end of list.

   if ( varname == ' ' .or.  dartstr == ' ' ) then
      error_string_1 = 'model_nml: variable list not fully specified'
      error_string_2 = 'reading from "'//trim(filter_io_filename)//'"'
      call error_handler(E_ERR, routine, error_string_1, &
                         source, revision, revdate, text2=error_string_2)
   endif

   ! The internal DART routines check if the variable name is valid.

   ! Make sure DART kind is valid
   quantity = get_index_for_quantity(dartstr)
   if( quantity < 0 ) then
      write(error_string_1,'(''there is no obs_kind "'',a,''" in obs_kind_mod.f90'')') &
           trim(dartstr)
      call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
   endif

   ! All good to here - fill the output variables

   nvar = nvar + 1
   if (variables(VT_ORIGININDX,i) == 'neutrals') nvar_neutral = nvar_neutral + 1
   if (variables(VT_ORIGININDX,i) == 'ions') nvar_ion = nvar_ion + 1
   var_names( nvar)   = varname
   var_qtys(  nvar)   = quantity
   var_ranges(nvar,:) = (/ MISSING_R8, MISSING_R8 /)
   var_update(nvar)   = .false.   ! at least initially

   ! convert the [min,max]valstrings to numeric values if possible
   read(minvalstring,*,iostat=io) minvalue
   if (io == 0) var_ranges(nvar,1) = minvalue

   read(maxvalstring,*,iostat=io) maxvalue
   if (io == 0) var_ranges(nvar,2) = maxvalue

   call to_upper(state_or_aux)
   if (state_or_aux == 'UPDATE') var_update(nvar) = .true.

enddo MY_LOOP

if (nvar == MAX_STATE_VARIABLES) then
   error_string_1 = 'WARNING: you may need to increase "MAX_STATE_VARIABLES"'
   write(error_string_2,'(''you have specified at least '',i4,'' perhaps more.'')') nvar
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate, text2=error_string_2)
endif

end subroutine verify_variables

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
    call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
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
    call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
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
      call error_handler(E_ERR, routine, error_string_1, source, revision, revdate, &
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
      call error_handler(E_ERR, routine, error_string_1, source, revision, revdate, &
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
! return 0 (ok) if we know how to interpolate this quantity.
! if it is a field in the state, return the variable id from
! the state structure.  if not in the state, varid will return -1

subroutine ok_to_interpolate(qty, varid, istatus)

integer, intent(in)  :: qty
integer, intent(out) :: varid
integer, intent(out) :: istatus

! See if the state contains the obs quantity
varid = get_varid_from_kind(dom_id, qty)

! in the state vector
if (varid > 0) then
   istatus = 0
   return
endif

! add any quantities that can be interpolated to this list if they
! are not in the state vector.
select case (qty)
   case (QTY_GEOMETRIC_HEIGHT)
      istatus = 0
   case default
      istatus = UNKNOWN_OBS_QTY_ERROR_CODE
end select
    
end subroutine ok_to_interpolate

!-----------------------------------------------------------------------
! Converts Aether restart files to a netCDF file
!
! This routine needs:
!
! 1.  A base dirname for the restart files (aether_restart_dirname).
! The filenames have the format 'dirname/{neutrals,ions}_mMMMM_gBBBB.rst'  
! where BBBB is the block number, MMMM is the member number,
! and they have leading 0s.   Blocks start in the
! southwest corner of the lat/lon grid and go east first, 
! then to the west end of the next row north and end in the northeast corner. 
! 
! In the process, the routine will find:
!
! 1. The number of blocks in Lon and Lat (nblocks_lon, nblocks_lat)
!
! 2. The number of lons and lats in a single grid block  (nx_per_block, ny_per_block, nz_per_block)
!
! 3. The overall grid size, {nlon,nlat,nalt} when you've read in all the blocks. 
!
! 4. The number of neutral species (and probably a mapping between
!    the species number and the variable name)  (nvar_neutral)
!
! 5. The number of ion species (ditto - numbers <-> names) (nvar_ion)
!
! In addition to reading in the state data, it fills Longitude, Latitude, and Altitude arrays.  
! This grid is orthogonal and rectangular but can have irregular spacing along
! any of the three dimensions.

subroutine restart_files_to_netcdf(member)

integer, intent(in) :: member

integer :: ncid

character(len=*), parameter :: routine = 'restart_files_to_netcdf'

if (module_initialized ) then
    write(error_string_1,'(3A)')'The aether static_init_model was already initialized but ', &
         trim(routine), ' uses a separate initialization procedure'
    call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
end if

call static_init_blocks("aether_to_dart_nml")

write(filter_io_filename,'(2A, I0.4, A3)') trim(filter_io_root),'_', member + 1,'.nc'
! nc_create_file does not leave define mode
ncid = nc_create_file(filter_io_filename)

call error_handler(E_MSG, '', '')
write(error_string_1,'(3A)') 'converting Aether restart files in directory ', &
                 "'"//trim(aether_restart_dirname)//"'"
write(error_string_2,'(3A)') ' to the NetCDF file ', "'"//trim(filter_io_filename)//"'"
call error_handler(E_MSG, routine, error_string_1, text2=error_string_2)
call error_handler(E_MSG, '', '')

! TODO: we haven't settled on the mechanism for identifying the state vector field names and source.
!       (defined type, arrays, named indices,...)
! TODO: def_fill_dimvars functionality was in nc_write_model_atts but shouldn't have been.  
!       I separated nc_write_model_atts into to parts and this is one of them.
!       Is this the best place for the call?  It's in the "define" section for the filter_input file.
!       It works.
call def_fill_dimvars(ncid)

! Write_model_time will make a time variable, if needed, which it is not.
call write_model_time(ncid, state_time)

! Define (non-time) variables
call restarts_to_filter(aether_restart_dirname, ncid, member, define=.true.)

! Read and convert (non-time) variables
call restarts_to_filter(aether_restart_dirname, ncid, member, define=.false.)
! subr. called by this routine closes the file only if define = .true.
call nc_close_file(ncid)

call error_handler(E_MSG, '', '')
write(error_string_1,'(3A)') 'Successfully converted the Aether restart files to ', &
                 "'"//trim(filter_io_filename)//"'"
call error_handler(E_MSG, routine, error_string_1)
call error_handler(E_MSG, '', '')
    
end subroutine restart_files_to_netcdf
    
!-----------------------------------------------------------------------
! Writes the state variables from a dart state vector (1d array)
! into Aether netcdf restart file sets.

subroutine netcdf_to_restart_files(member)
    
integer, intent(in) :: member

integer :: ncid
character(len=*), parameter :: routine = 'netcdf_to_restart_files:'
    
! write out the state vector data.  
! when this routine returns all the data has been written.
    
if (module_initialized ) then
    write(error_string_1,'(3A)')'The aether mod was already initialized but ', &
      trim(routine), ' uses a separate initialization procedure'
    call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
end if

call static_init_blocks("dart_to_aether_nml")

write(filter_io_filename,'(2A,I0.4,A3)') trim(filter_io_root),'_',member + 1,'.nc'

call error_handler(E_MSG, routine, '', '', revision, revdate)
write(error_string_1,'(3A)') 'Extracting fields from DART file ', "'"//trim(filter_io_filename)//"'"
write(error_string_2,'(3A)') 'into Aether restart files in directory ', &
     "'"//trim(aether_restart_dirname)//"'"
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate, text2=error_string_2)

ncid = nc_open_file_readonly(filter_io_filename, routine)

call filter_to_restarts(ncid, member)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------
call error_handler(E_MSG, routine,'','', revision, revdate)
write(error_string_1,'(3A)') 'Successfully converted to the Aether restart files in directory'
write(error_string_2,'(3A)') "'"//trim(aether_restart_dirname)//"'"
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate, text2=error_string_2)

call nc_close_file(ncid)
    
end subroutine netcdf_to_restart_files

!-----------------------------------------------------------------------
! ? Will this need to open the grid_{below,corners,down,left} filetypes?
!   This code can handle it; a longer filetype passed in, and no member.
! ? Aether output files?

function block_file_name(filetype, memnum, blocknum)

character(len=*), intent(in)  :: filetype  ! one of {grid,ions,neutrals}
integer,          intent(in)  :: blocknum
integer,          intent(in)  :: memnum
character(len=128)            :: block_file_name

character(len=*), parameter :: routine = 'block_file_name'

block_file_name = trim(filetype)
if (memnum   >= 0) write(block_file_name, '(A,A2,I0.4)') trim(block_file_name), '_m', memnum
if (blocknum >= 0) write(block_file_name, '(A,A2,I0.4)') trim(block_file_name), '_g', blocknum
block_file_name = trim(block_file_name)//'.nc'
if ( debug > 0 ) then
   write(error_string_1,'("filename, memnum, blocknum = ",A,2(1x,i5))') &
        trim(block_file_name), memnum, blocknum
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
endif
    
end function block_file_name

!-----------------------------------------------------------------------
! Like static_init_model, but for aether_to_dart and dart_to_aether.
!   Read the namelist, 
!   parse the 'variables' table,
!   get the Aether grid information
!   convert the Aether time into a DART time.

subroutine static_init_blocks(nml)

character(len=*), intent(in) :: nml

character(len=128) :: aether_filter_io_filename
integer            :: iunit, io

character(len=*), parameter :: routine = 'static_init_blocks'

if (module_initialized) return ! only need to do this once

! This prevents subroutines called from here from calling static_init_mod.
module_initialized = .true.

!------------------
! Read the namelist

call find_namelist_in_file("input.nml", trim(nml), iunit)
if (trim(nml) == 'aether_to_dart_nml') then
   read(iunit, nml = aether_to_dart_nml, iostat = io)
   ! Record the namelist values used for the run
   if (do_nml_file()) write(nmlfileunit, nml=aether_to_dart_nml)
   if (do_nml_term()) write(     *     , nml=aether_to_dart_nml)
else if (trim(nml) == 'dart_to_aether_nml') then
   read(iunit, nml = dart_to_aether_nml, iostat = io)
   ! Record the namelist values used for the run
   if (do_nml_file()) write(nmlfileunit, nml=dart_to_aether_nml)
   if (do_nml_term()) write(     *     , nml=dart_to_aether_nml)
endif
call check_namelist_read(iunit, io, trim(nml)) ! closes, too.


! error-check, convert namelist input to arrays.
! 'variables' comes from the namelist in input.nml
call verify_variables(variables, filter_io_filename, nvar, var_names, var_qtys, var_ranges, var_update)

!--------------------------------
! TODO:  Set the time step 
! Ensures model_advance_time is multiple of 'dynamics_timestep'

! Aether uses Julian time internally, andor a Julian calendar
! (days from the start of the calendar), depending on the context)
call set_calendar_type( calendar )   

!--------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the block restart files, could be stretched ...

call get_grid_info_from_blocks(aether_restart_dirname, nlon, nlat, nlev, nblocks_lon, &
               nblocks_lat, nblocks_lev, lat_start, lat_end, lon_start)

if( debug  > 0 ) then
    write(error_string_1,'(A,3I5)') 'grid dims are ', nlon, nlat, nlev
    call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
endif

! Opens and closes the grid block file, but not the filter netcdf file.
call get_grid_from_blocks(aether_restart_dirname, nblocks_lon, nblocks_lat, nblocks_lev, &
   nx_per_block, ny_per_block, nz_per_block)
! , lons, lats, levs )

! Convert the Aether reference date (not calendar day = 0 date)
! to the days and seconds of the calendar set in model_mod_nml.
aether_ref_time = set_date(aether_ref_date(1), aether_ref_date(2), aether_ref_date(3), &
                     aether_ref_date(4), aether_ref_date(5))
call get_time(aether_ref_time, aether_ref_nsecs, aether_ref_ndays)

! Get the model time from a restart file.
aether_filter_io_filename = block_file_name(variables(VT_ORIGININDX,1), 0, 0)
state_time = read_aether_time(trim(aether_restart_dirname)//'/'//trim(aether_filter_io_filename))

if ( debug > 0 ) then
  write(error_string_1,'("grid: nlon, nlat, nlev =",3(1x,i5))') nlon, nlat, nlev
  call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
endif
    
end subroutine static_init_blocks

!-----------------------------------------------------------------------
! Read Aether's info file (UAM.in) to get a description of the restart file blocks' grids.

subroutine get_grid_info_from_blocks(restart_dirname, nlon, nlat, nlev,     &
                                     nblocks_lon, nblocks_lat, nblocks_lev, &
                                     lat_start, lat_end, lon_start)

character(len=*), intent(in)  :: restart_dirname
integer,          intent(out) :: nlon   ! Number of Longitude centers
integer,          intent(out) :: nlat   ! Number of Latitude  centers
integer,          intent(out) :: nlev   ! Number of Vertical grid centers
integer,          intent(out) :: nblocks_lon, nblocks_lat, nblocks_lev
real(r8),         intent(out) :: lat_start, lat_end, lon_start

character(len=100) :: c_line
character(len=256) :: file_loc
integer            :: i, iunit, ios

! TODO: get the grid info from a namelist (98 variables), instead of Aether's UAM.in.  
!       Then remove functions read_in_*.
!       The rest of the UAM.in contents are for running Aether.
!       Can wait until aether_to_dart push is done.
character(len=*), parameter :: filename = 'UAM.in'
character(len=*), parameter :: routine = 'get_grid_info_from_blocks'

! get the ball rolling ...

nblocks_lon = 0
nblocks_lat = 0
nblocks_lev = 0
lat_start   = 0.0_r8
lat_end     = 0.0_r8
lon_start   = 0.0_r8

write(file_loc,'(a,''/'',a)') trim(restart_dirname), trim(filename)

if (debug > 4) then
write(error_string_1,'(3A)') 'Now opening Aether UAM file: ', trim(file_loc)
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
end if


iunit = open_file(trim(file_loc), action='read')

UAMREAD : do i = 1, 1000000

read(iunit,'(a)',iostat=ios) c_line

if (ios /= 0) then
! If we get to the end of the file or hit a read error without
! finding what we need, die.
write(error_string_1,'(3A)') 'cannot find #GRID in ', trim(file_loc)
call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
endif

if (c_line(1:5) .ne. "#GRID") cycle UAMREAD

nblocks_lon = read_in_int( iunit,'nblocks_lon', trim(file_loc))
nblocks_lat = read_in_int( iunit,'nblocks_lat', trim(file_loc))
nblocks_lev = read_in_int( iunit,'nblocks_lev', trim(file_loc))
lat_start   = read_in_real(iunit,'lat_start',  trim(file_loc))
lat_end     = read_in_real(iunit,'lat_end',    trim(file_loc))
lon_start   = read_in_real(iunit,'lon_start',  trim(file_loc))

exit UAMREAD

enddo UAMREAD

if (debug > 4) then
write(error_string_1,'(3A)') 'Successfully read Aether UAM grid file:', trim(file_loc)
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,I5)') '   nblocks_lon:', nblocks_lon
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,I5)') '   nblocks_lat:', nblocks_lat
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,I5)') '   nblocks_lev:', nblocks_lev
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,F15.4)') '   lat_start:', lat_start
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,F15.4)') '   lat_end:', lat_end
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,F15.4)') '   lon_start:', lon_start
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
end if

call close_file(iunit)

end subroutine get_grid_info_from_blocks

!-----------------------------------------------------------------------
! Read block grid values (2D arrays) from a grid NetCDF file.
! Allocate and fill the full-domain 1-D dimension arrays (lon, lat, levs)

subroutine get_grid_from_blocks(dirname, nblocks_lon, nblocks_lat, nblocks_lev, &
                                nx_per_block, ny_per_block, nz_per_block)
! ,       &
!                                 lons, lats, levs )

character(len=*), intent(in)  :: dirname
integer,          intent(in)  :: nblocks_lon ! Number of Longitude blocks
integer,          intent(in)  :: nblocks_lat ! Number of Latitude  blocks
integer,          intent(in)  :: nblocks_lev ! Number of Altitude  blocks
integer,          intent(out) :: nx_per_block ! Number of non-halo Longitude centers per block
integer,          intent(out) :: ny_per_block ! Number of non-halo Latitude  centers per block
integer,          intent(out) :: nz_per_block ! Number of Vertical grid centers
! real(r8), allocatable, intent(inout) :: lons(:), lats(:), levs(:)

integer               :: nb, offset, ncid, nboff
integer               :: starts(3), ends(3), xcount, ycount, zcount
character(len=128)    :: filename
real(r4), allocatable :: temp(:,:,:)

character(len=*), parameter :: routine = 'get_grid_from_blocks'

! TODO: Here it needs to read the x,y,z  from a NetCDF block file(s),
!       in order to calculate the n[xyz]_per_block dimensions. 
!       grid_g0000.nc looks like a worthy candidate, but a restart could be used.
write (filename,'(2A)')  trim(dirname),'/grid_g0000.nc'
ncid = nc_open_file_readonly(filename, routine)

! The grid (and restart) file variables have halos, so strip them off
! to get the number of actual data values in each dimension of the block.
nx_per_block = nc_get_dimension_size(ncid, 'x', routine) - (2 * nghost)
ny_per_block = nc_get_dimension_size(ncid, 'y', routine) - (2 * nghost)
nz_per_block = nc_get_dimension_size(ncid, 'z', routine)

nlon = nblocks_lon * nx_per_block
nlat = nblocks_lat * ny_per_block
nlev = nblocks_lev * nz_per_block     

write(error_string_1,'(A,I5)')  'nlon = ', nlon
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,I5)')  'nlat = ', nlat
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
write(error_string_1,'(A,I5)')  'nlev = ', nlev
call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)

! TODO; do these need to be deallocated somewhere?
!       Probably not; this is only done once, and these arrays are needed
!       through most of the a2d and d2a programs.
allocate( lons( nlon ))
allocate( lats( nlat ))
allocate( levs( nlev ))

if (debug > 4) then
   write(error_string_1,'(2A)') 'Successfully read Aether grid file:', trim(filename)
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   write(error_string_1,'(A,I5)') '   nx_per_block:', nx_per_block
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   write(error_string_1,'(A,I5)') '   ny_per_block:', ny_per_block
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   write(error_string_1,'(A,I5)') '   nz_per_block:', nz_per_block
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
endif

! A temp array large enough to hold any of the 3D
! Lon, Lat or Alt arrays from a block plus ghost cells.
! The restart files have C-indexing (fastest changing dim is the last).
allocate(temp( 1:nz_per_block, &
               1-nghost:ny_per_block+nghost, &
               1-nghost:nx_per_block+nghost))
temp = MISSING_R4 

starts(1) = 1 - nghost
starts(2) = 1 - nghost
starts(3) = 1
ends(1)   = nx_per_block + nghost
ends(2)   = ny_per_block + nghost
ends(3)   = nz_per_block
xcount = nx_per_block + (2 * nghost)
ycount = ny_per_block + (2 * nghost)
zcount = nz_per_block
if ( debug > 0 ) then
   write(error_string_1,'(2(A,3i5),A,3(1X,i5))') &
        'starts = ',starts, 'ends = ',ends, '[xyz]counts = ',xcount,ycount,zcount
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
endif

! go across the south-most block row picking up all longitudes
do nb = 1, nblocks_lon

   ! filename is trimmed by passage to open_block_file + "len=*" there.
   filename = block_file_name('grid', -1, nb-1)
   ncid = open_block_file(filename, 'read')

   ! Read 3D array and extract the longitudes of the non-halo data of this block.
   ! The restart files have C-indexing (fastest changing dim is the last),
   ! So invert the dimension bounds.
   call nc_get_variable(ncid, 'Longitude', &
        temp(starts(3):ends(3), starts(2):ends(2), starts(1):ends(1)), &
        context=routine, &
        nc_count=(/ zcount,ycount,xcount /))

   offset = (nx_per_block * (nb - 1))
   lons(offset+1:offset+nx_per_block) = temp(1,1,1:nx_per_block)

   call nc_close_file(ncid)
enddo

! go up west-most block row picking up all latitudes
do nb = 1, nblocks_lat

   ! Aether's block name counter start with 0, but the lat values can come from 
   ! any lon=const column of blocks. 
   nboff = ((nb - 1) * nblocks_lon)
   filename = block_file_name('grid', -1, nboff)
   ncid = open_block_file(filename, 'read')

   call nc_get_variable(ncid, 'Latitude', &
        temp(starts(3):ends(3), starts(2):ends(2), starts(1):ends(1)), &
        context=routine, nc_count=(/zcount,ycount,xcount/))


   offset = (ny_per_block * (nb - 1))
   lats(offset+1:offset+ny_per_block) = temp(1,1:ny_per_block,1)

   call nc_close_file(ncid)
enddo


! this code assumes UseTopography is false - that all columns share
! the same altitude array, so we can read it from the first block.
! if this is not the case, this code has to change.

filename = block_file_name('grid', -1, 0)
ncid = open_block_file(filename, 'read')

temp = MISSING_R8
call nc_get_variable(ncid, 'Altitude', &
     temp(starts(3):ends(3), starts(2):ends(2), starts(1):ends(1)), &
     context=routine, nc_count=(/zcount,ycount,xcount/))

levs(1:nz_per_block) = temp(1:nz_per_block,1,1)

call nc_close_file(ncid)

deallocate(temp)

! convert from radians into degrees
lons = lons * RAD2DEG
lats = lats * RAD2DEG

if (debug > 4) then
   print *, routine, 'All lons ', lons
   print *, routine, 'All lats ', lats
   print *, routine, 'All levs ', levs
endif

if ( debug > 1 ) then ! Check dimension limits
   write(error_string_1,'(A,2F15.4)') 'LON range ', minval(lons), maxval(lons)
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   write(error_string_1,'(A,2F15.4)') 'LAT range ', minval(lats), maxval(lats)
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   write(error_string_1,'(A,2F15.4)') 'ALT range ', minval(levs), maxval(levs)
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
endif

end subroutine get_grid_from_blocks

!-----------------------------------------------------------------------
! Read the Aether restart file time and convert to a DART time.

function read_aether_time(filename)
type(time_type)              :: read_aether_time
character(len=*), intent(in) :: filename

integer  :: ncid
integer  :: tsimulation   ! the time read from a restart file; seconds from aether_ref_date.
integer  :: ndays, nsecs

character(len=*), parameter :: routine = 'read_aether_time'

tsimulation = MISSING_I

ncid = open_block_file(filename, 'read')
call nc_get_variable(ncid, 'time', tsimulation, context=routine)
call nc_close_file(ncid, routine, filename)

! Calculate the DART time of the file time.
! TODO: review calculation of ndays in read_aether_time
ndays     = tsimulation / 86400
nsecs     = tsimulation - (ndays * 86400)
! The ref day is not finished, but don't need to subtract 1 because 
! that was accounted for in the integer calculation of ndays.
ndays     = aether_ref_ndays + ndays
read_aether_time = set_time(nsecs, ndays)

if (do_output()) &
    call print_time(read_aether_time, routine//': time in restart file '//filename)
if (do_output()) &
    call print_date(read_aether_time, routine//': date in restart file '//filename)

if (debug > 8) then
   write(error_string_1,'(A,I5)')'tsimulation ', tsimulation
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   write(error_string_1,'(A,I5)')'ndays       ', ndays
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   write(error_string_1,'(A,I5)')'nsecs       ', nsecs
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)

   call print_date(aether_ref_time, routine//':model base date')
   call print_time(aether_ref_time, routine//':model base time')
endif

end function read_aether_time

!-----------------------------------------------------------------------
! Convert Aether's non-CF-compliant names into CF-compliant names for filter.
! For the ions, it moves the name of the ion from the end of the variable names
! to the beginning.

function aether_name_to_dart(varname)

character(len=vtablenamelength), intent(in) :: varname

character(len=vtablenamelength) :: aether_name_to_dart, aether
character(len=64)               :: parts(8), var_root
integer                         :: char_num, first, i_parts, aether_len, end_str

aether     = trim(varname)
aether_len = len_trim(varname)
parts = ''

! Look for the last ' '.  The characters after that are the species.
! If there's no ' ', the whole string is the species.
char_num = 0
char_num = scan(trim(aether),' ', back=.true.)
var_root = aether(char_num+1:aether_len)
! purge_chars removes unwanted [()\] 
parts(1) = purge_chars( trim(var_root),')(\', plus_minus=.true.)
! TODO: keep aether_name_to_dart diagnostic?  Then add routine, error_handler.
! print*,'var_root, parts(1) = ', var_root, parts(1) 
end_str  = char_num

! Tranform remaining pieces of varname into DART versions.
char_num = MISSING_I
first = 1
i_parts = 2
P_LOOP: do
   ! This returns the position of the first blank *within the substring* passed in.
   char_num = scan(aether(first:end_str),' ', back=.false.)
   if (char_num > 0 .and. first < aether_len) then
      parts(i_parts) = purge_chars(aether(first:first+char_num-1), '.)(\', plus_minus=.true.)

      first   = first + char_num 
      i_parts = i_parts + 1
   else
      exit P_LOOP
   endif
enddo P_LOOP

! Construct the DART field name from the parts
aether_name_to_dart = trim(parts(1))
i_parts = 2
Build : do
   if (trim(parts(i_parts)) /= '') then
      aether_name_to_dart = trim(aether_name_to_dart)//'_'//trim(parts(i_parts))
      i_parts = i_parts + 1
   else
      exit Build
   endif
enddo Build
   
end function aether_name_to_dart
   
!-----------------------------------------------------------------------
! Replace undesirable characters with better.
   
function purge_chars(ugly_string, chars, plus_minus)
   
character (len=*), intent(in) :: ugly_string, chars
logical,           intent(in) :: plus_minus
character (len=64)            :: purge_chars

character (len=256) :: temp_str

integer :: char_num, end_str, pm_num

! Trim is not needed here
temp_str = ugly_string
end_str  = len_trim(temp_str)
char_num = MISSING_I
Squeeze : do 
   ! Returns 0 if chars are not found
   char_num = scan(temp_str, chars)
   ! Need to change it to a char that won't be found by scan in the next iteration,
   ! and can be easily removed.
   if (char_num > 0) then
      ! Squeeze out the character
      temp_str(char_num:end_str-1) = temp_str(char_num+1:end_str)
      temp_str(end_str:end_str) = ''
!       temp_str(char_num:char_num) = ' '
   else
      exit Squeeze
   endif
enddo Squeeze

! Replace + and - with pos and neg.  Assume there's only 1.
temp_str = trim(adjustl(temp_str))
end_str  = len_trim(temp_str)
pm_num   = scan(trim(temp_str),'+-', back=.false.)
if (pm_num == 0 .or. .not. plus_minus) then
   purge_chars = trim(temp_str)
else
   if (temp_str(pm_num:pm_num) == '+') then
      purge_chars = temp_str(1:pm_num-1)//'pos'
   else if (temp_str(pm_num:pm_num) == '-') then
      purge_chars = temp_str(1:pm_num-1)//'neg'
   endif
   if (pm_num + 1 <= end_str) &
      purge_chars = trim(purge_chars)//temp_str(pm_num+1:end_str)
endif
   
end function purge_chars

!-----------------------------------------------------------------------
! Open an Aether restart block file (neutral, ion, ...?)

function open_block_file(filename, rw)

! filename is trimmed by this definition
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: rw   ! 'read' or 'readwrite'
integer                      :: open_block_file

character(len=*), parameter :: routine = 'open_block_file'

if ( .not. file_exist(filename) ) then
   write(error_string_1,'(4A)') 'cannot open file ', filename,' for ', rw
   call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
endif

if (debug > 0) then
   write(error_string_1,'(4A)') 'Opening file ', trim(filename), ' for ', rw
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
end if


if (rw == 'read') then
   open_block_file = nc_open_file_readonly(filename, routine)
else if (rw == 'readwrite') then
   open_block_file = nc_open_file_readwrite(filename, routine)
else
   error_string_1 = ': must be called with rw={read,readwrite}, not '//rw
   call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
endif


if (debug > 80) then
   write(error_string_1,'(4A)') 'Returned file descriptor is ', open_block_file
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
end if
    
end function open_block_file

!-----------------------------------------------------------------------
! Open all restart files (blocks x {neutrals,ions}) for 1 member 
! and transfer the requested variable contents to the filter input file.
! This is called with 'define' = 
! .true.  define variables in the file or 
! .false. transfer the data from restart files to a filter_inpu.nc file.

subroutine restarts_to_filter(dirname, ncid_output, member, define)

character(len=*), intent(in)  :: dirname
integer,          intent(in)  :: ncid_output, member
logical,          intent(in)  :: define

integer :: ib, jb, ib_loop, jb_loop

if (define) then
   ! if define, run one block.
   ! the block_to_filter_io call defines the variables in the whole domain netCDF file.
   ib_loop = 1
   jb_loop = 1
   ! nc_write_model_atts puts it in define, and takes it out.
   call nc_begin_define_mode(ncid_output)
else
   ! if not define, and run all blocks.
   ! the block_to_filter_io call adds the (ib,jb) block to a netCDF variable 
   ! in order to make a file containing the data for all the blocks.
   ib_loop = nblocks_lon
   jb_loop = nblocks_lat
end if

do jb = 1, jb_loop
   do ib = 1, ib_loop

      call block_to_filter_io(ncid_output, dirname, ib, jb, member, define)

   enddo
enddo

if (define) then
   call nc_end_define_mode(ncid_output)
endif
    
end subroutine restarts_to_filter

!-----------------------------------------------------------------------
! Read in a real number from the UAM.in file.
! TODO: the file name should not be filter_io_filename.

function read_in_real(iunit, varname, filter_io_filename)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname, filter_io_filename
real(r8)                     :: read_in_real

character(len=100) :: c_line
integer            :: i, ios
character(len=*), parameter :: routine = 'read_in_real'

! Read a line 
read(iunit,'(a)',iostat=ios) c_line
if (ios /= 0) then
   write(error_string_1,'(4A)') 'cannot find '//trim(varname)//' in '//trim(filter_io_filename)
   call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
endif

! Remove anything after a space or TAB
i=index(c_line,' ');     if( i > 0 ) c_line(i:len(c_line))=' '
i=index(c_line,char(9)); if( i > 0 ) c_line(i:len(c_line))=' '

! Now that we have a line with nothing else ... parse it
read(c_line,*,iostat=ios) read_in_real

if(ios /= 0) then
   write(error_string_1,'(4A)')'unable to read '//trim(varname)//' in '//trim(filter_io_filename)
   call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
endif
        
end function read_in_real

!-----------------------------------------------------------------------
! Read in an integer from the UAM.in file.
! TODO: the file name should not be filter_io_filename.

function read_in_int(iunit, varname, filter_io_filename)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname, filter_io_filename
integer                      :: read_in_int

character(len=100) :: c_line
integer            :: i, ios
character(len=*), parameter :: routine = 'read_in_int'

! Read a line 
read(iunit,'(a)',iostat=ios) c_line
if (ios /= 0) then
   write(error_string_1,'(4A)') 'cannot find '//trim(varname)//' in '//trim(filter_io_filename)
   call error_handler(E_ERR,'get_grid_dims', error_string_1, source, revision, revdate)
endif

! Remove anything after a space or TAB
i=index(c_line,' ');     if( i > 0 ) c_line(i:len(c_line))=' '
i=index(c_line,char(9)); if( i > 0 ) c_line(i:len(c_line))=' '

read(c_line,*,iostat=ios) read_in_int

if(ios /= 0) then
   write(error_string_1,'(4A)')'unable to read '//trim(varname)//' in '//trim(filter_io_filename)
   call error_handler(E_ERR, routine, error_string_1, source, revision, revdate, &
                      text2=c_line)
endif
    
end function read_in_int

!-----------------------------------------------------------------------
! Open all restart files (neutrals,ions) for a block and read in the requested data items.
! The write_filter_io calls will write the data to the filter_input.nc.

subroutine write_filter_io(data3d, varname, ib, jb, ncid)

real(r4), intent(in) :: data3d(1:nz_per_block, &
                               1-nghost:ny_per_block+nghost, &
                               1-nghost:nx_per_block+nghost)

character(len=vtablenamelength), intent(in) :: varname
integer,  intent(in) :: ib, jb
integer,  intent(in) :: ncid

integer :: starts(3)

character(len=*), parameter :: routine = 'write_filter_io'

! write(varname,'(A)') trim(variables(VT_VARNAMEINDX,ivar))

! to compute the start, consider (ib-1)*nx_per_block+1
starts(1) = 1
starts(2) = (jb-1) * ny_per_block + 1
starts(3) = (ib-1) * nx_per_block + 1
! TODO: convert to error_msg
! print*, routine,'; starts = ', starts
! print*, routine,'; counts = ', nz_per_block, ny_per_block, nx_per_block,1

call nc_put_variable(ncid, varname, &
                     data3d(1:nz_per_block,1:ny_per_block,1:nx_per_block), &
                     context=routine, nc_start=starts, &
                     nc_count=(/nz_per_block,ny_per_block,nx_per_block/))
! TODO: convert to error_msg
! print*, routine,': filled varname = ', varname 
   
end subroutine write_filter_io

!-----------------------------------------------------------------------
! Transfer variable data from a block restart file to the filter_input.nc file. 
! It's called with 2 modes:
! define = .true.  define the NC variables in the filter_input.nc 
! define = .false. write the data from a block to the NC file using write_filter_io.

subroutine block_to_filter_io(ncid_output, dirname, ib, jb, member, define)

integer,          intent(in) :: ncid_output
character(len=*), intent(in) :: dirname
integer,          intent(in) :: ib, jb
integer,          intent(in) :: member
logical,          intent(in) :: define

real(r4), allocatable :: temp1d(:), temp2d(:,:), temp3d(:,:,:)
real(r4), allocatable :: alt1d(:), density_ion_e(:,:,:)
integer               :: ivar, nb, ncid_input
! TEC? integer               :: maxsize
!      logical               :: no_idensity
!      real(r4)              :: temp0d 
character(len=32)     :: att_val 
character(len=128)    :: file_root 
character(len=256)    :: filename
character(len=vtablenamelength) :: varname, dart_varname

character(len=*), parameter :: routine = 'block_to_filter_io'

! The block number, as counted in Aether.
! Lower left is 0, increase to the East, then 1 row farther north, West to East.
nb = (jb - 1) * nblocks_lon + ib - 1

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nghost:max(nx_per_block, ny_per_block, nz_per_block) + nghost))

! treat alt specially since we want to derive TEC here
! TODO: See density_ion_e too.
allocate( alt1d(1-nghost:max(nx_per_block, ny_per_block, nz_per_block) + nghost))

! temp array large enough to hold any 2D field 
allocate(temp2d(1-nghost:ny_per_block+nghost, &
                1-nghost:nx_per_block+nghost))

! TODO: We need all altitudes, but there might be vertical blocks in the future.
!       But there would be no vertical halos.
!       Make nzcount adapt to whether there are blocks.
!       And temp needs to have C-ordering, which is what the restart files have.
! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1:nz_per_block, &
                1-nghost:ny_per_block+nghost, &
                1-nghost:nx_per_block+nghost))

! TODO: Waiting for e- guidance from Aaron.
! save density_ion_e to compute TEC
allocate(density_ion_e(1:nz_per_block, &
                       1-nghost:ny_per_block+nghost, &
                       1-nghost:nx_per_block+nghost))

! TODO: Aether gives a unique name to each (of 6) velocity components.
!       Do we want to use a temp4d array to handle them?  
!       They are independent variables in the block files (and state).
! ! temp array large enough to hold velocity vect, etc
! maxsize = max(3, nvar_ion)
! allocate(temp4d(1-nghost:nx_per_block+nghost, &
!                 1-nghost:ny_per_block+nghost, &
!                 1-nghost:nz_per_block+nghost, maxsize))


! TODO; Does Aether need a replacement for these Density fields?  Yes.
!       But they are probably read by the loops below.
!       Don't need to fetch index because Aether has NetCDF restarts,
!       so just loop over the field names to read.
! 
! ! assume we could not find the electron density for VTEC calculations
! no_idensity = .true.
! 
! if (inum > 0) then
!    ! one or more items in the state vector need to replace the
!    ! data in the output file.  loop over the index list in order.
!    j = 1
! ! TODO:   electron density is not in the restart files, but it's needed for TEC
!           In Aether they will be from an ions file, but now only from an output file (2023-10-30).
!           Can that be handled like the neutrals and ions files, using variables(VT_ORIGININDX,:)
!           to build an output file name?  Are outputs in block form?
!                ! save the electron density for TEC computation
!                density_ion_e(:,:,:) = temp3d(:,:,:)

! Handle the 2 restart file types (ions and neutrals).
! Each field has a file type associated with it: variables(VT_ORIGININDX,f_index)
! TODO: for now require that all neutrals are listed in variables before the ions.

file_root = variables(VT_ORIGININDX,1)
filename = block_file_name(file_root, member, nb)
ncid_input = open_block_file(filename, 'read')

! TODO: prints > ERR_MSG?
if (debug >= 100 .and. do_output()) print*,'block_to_filter_io: nvar_neutral = ', nvar_neutral
do ivar = 1, nvar_neutral
   ! The nf90 functions cannot read the variable names with the '\'s in them.
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   if (debug >= 100 .and. do_output()) print*, routine,'varname = ', varname
   ! Translate the Aether field name into a DART field name.
   dart_varname = aether_name_to_dart(varname)

   ! TODO: Given the subroutine name, perhaps these definition sections should be 
   !       one call higher up, with the same loop around it.
   if (define) then
   ! Define the variable in the filter_input.nc file (the output from this program).
   ! The calling routine entered define mode.

      if (debug > 10 .and. do_output()) then
         write(error_string_1,'(A,I0,2A)') 'Defining ivar = ', ivar,':', dart_varname
         call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
      end if
   
      call nc_define_real_variable(ncid_output, dart_varname, &
           (/ LEV_DIM_NAME, LAT_DIM_NAME, LON_DIM_NAME/) )
      call nc_get_attribute_from_variable(ncid_input, varname, 'units', att_val, routine)
      call nc_add_attribute_to_variable(ncid_output, dart_varname, 'units', att_val, routine)

   else if (file_root == 'neutrals') then
   ! Read 3D array and extract the non-halo data of this block.
! TODO: There are no 2D or 1D fields in ions or neutrals, but there could be; different temp array.
      call nc_get_variable(ncid_input, varname, temp3d, context=routine)
      if (debug >= 100 .and. do_output()) then
         ! TODO convert to error_handler?  Or diagnostics are no longer useful?
         print*,'block_to_filter_io: temp3d = ', temp3d(1,1,1), temp3d(15,15,15), varname
         print*,'block_to_filter_io: define = ', define
      endif
      call write_filter_io(temp3d, dart_varname, ib, jb, ncid_output)
   else
      write(error_string_1,'(A,I3,A)') 'Trying to read neutrals, but variables(', &
           VT_ORIGININDX,ivar , ') /= "neutrals"'
      call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
   endif

enddo
call nc_close_file(ncid_input)

file_root = variables(VT_ORIGININDX,nvar_neutral+1)
filename = block_file_name(file_root, member, nb)
ncid_input = open_block_file(filename, 'read')

do ivar = nvar_neutral +1, nvar_neutral + nvar_ion
   ! Purging \ from aether name.
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   dart_varname = aether_name_to_dart(varname)

   if (define) then

      if (debug > 10 .and. do_output()) then
         write(error_string_1,'(A,I0,2A)') 'Defining ivar = ', ivar,':', dart_varname
         call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
      end if
   
      call nc_define_real_variable(ncid_output, dart_varname, &
                                   (/ LEV_DIM_NAME, LAT_DIM_NAME, LON_DIM_NAME/) )
      call nc_get_attribute_from_variable(ncid_input, varname, 'units', att_val, routine)
      call nc_add_attribute_to_variable(ncid_output, dart_varname, 'units', att_val, routine)
      print*, routine,': defined ivar, dart_varname, att = ', ivar, dart_varname, att_val

   else if (file_root == 'ions') then
      call nc_get_variable(ncid_input, varname, temp3d, context=routine)
      call write_filter_io(temp3d, dart_varname, ib, jb, ncid_output)
   else
      write(error_string_1,'(A,I3,A)') 'Trying to read ions, but variables(', &
           VT_ORIGININDX,ivar , ') /= "ions"'
      call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
   endif

enddo

! Leave file open if fields were just added (define = .false.),
! so that time can be added.
if (define) call nc_close_file(ncid_input)

! TODO: Does Aether need TEC to be calculated? Yes
! ! add the VTEC as an extended-state variable
! ! NOTE: This variable will *not* be written out to the Aether restart files 
! 
! if (no_idensity) then
!    write(error_string_1,*) 'Cannot compute the VTEC without the electron density'
!    call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
! end if
! 
!       temp2d = 0._r8
!       ! compute the TEC integral
!       do i =1,nz_per_block-1 ! approximate the integral over the altitude as a sum of trapezoids
!          ! area of a trapezoid: A = (h2-h1) * (f2+f1)/2
!          temp2d(:,:) = temp2d(:,:) + ( alt1d(i+1)-alt1d(i) )  * &
!                        ( density_ion_e(:,:,i+1)+density_ion_e(:,:,i) ) /2.0_r8
!       end do  
!       ! convert temp2d to TEC units
!       temp2d = temp2d/1e16_r8
!    call write_block_to_filter2d(temp2d, ivals(1), block, ncid, define) 

! TODO: Does Aether need f10_7 to be calculated or processed? Yes
! !gitm_index = get_index_start(domain_id, 'VerticalVelocity')
! call get_index_from_gitm_varname('f107', inum, ivals)
! if (inum > 0) then
!   call write_block_to_filter0d(temp0d, ivals(1), ncid, define) !see comments in the body of the subroutine
! endif
! 

deallocate(temp1d, temp2d, temp3d)
deallocate(alt1d, density_ion_e)
   
end subroutine block_to_filter_io

!-----------------------------------------------------------------------
! Extract (updated) variables from a filter_output.nc file
! and write to existing block restart files.

subroutine filter_to_restarts(ncid, member)

integer, intent(in) :: member, ncid

real(r4), allocatable           :: fulldom1d(:), fulldom3d(:,:,:)
character(len=256)              :: file_root
integer                         :: ivar
character(len=vtablenamelength) :: varname, dart_varname

character(len=*), parameter :: routine = 'filter_to_restarts'

! Space for full domain field (read from filter_output.nc)
! and halo around the full domain
allocate(fulldom3d(1:nlev, &
                   1-nghost:nlat+nghost, &
                   1-nghost:nlon+nghost))

! get the dirname, construct the filenames inside open_block_file

! >>> TODO: Not all fields have halos suitable for calculating gradients.  
!     These do (2023-11-8): neutrals; temperature, O, O2, N2, and the horizontal winds. 
!                           ions; none.
!     The current model_mod will fill all neutral halos anyway, 
!     since that's simpler and won't break the model.
!     TODO: add an attribute to the variables (?) to denote whether a field 
!           should have its halo filled.
do ivar = 1, nvar_neutral
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   if (debug >= 0 .and. do_output()) then
      write(error_string_1,'("varname = ",A)') trim(varname)
      call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   endif
   dart_varname = aether_name_to_dart(varname)

   file_root = trim(variables(VT_ORIGININDX,ivar))
   if (file_root == 'neutrals') then
      ! This parameter is available through the `use netcdf` command.
      fulldom3d = NF90_FILL_REAL
      
      call nc_get_variable(ncid, dart_varname, fulldom3d(1:nlev,1:nlat,1:nlon), &
                           nc_count=(/ nlev,nlat,nlon,1 /), context=routine)
      ! TODO: ncount not needed?  Reading the whole field.

      ! Copy updated field values to full domain halo.
      ! Block domains+halos will be easily read from this.
      call add_halo_fulldom3d(fulldom3d)

      call filter_io_to_blocks(fulldom3d, varname, file_root, member)
   else
      ! TODO: error; varname is inconsistent with VT_ORIGININDX
   endif

enddo

do ivar = nvar_neutral + 1, nvar_neutral + nvar_ion
   varname = purge_chars(trim(variables(VT_VARNAMEINDX,ivar)), '\', plus_minus=.false.)
   dart_varname = aether_name_to_dart(varname)

   file_root = trim(variables(VT_ORIGININDX,ivar))
   if (debug >= 0 .and. do_output()) then
      write(error_string_1,'("varname, dart_varname, file_root = ",3(2x,A))') &
             trim(varname), trim(dart_varname), file_root 
      call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   endif

   if (file_root == 'ions') then
      fulldom3d = NF90_FILL_REAL
      call nc_get_variable(ncid, dart_varname, fulldom3d(1:nlev,1:nlat,1:nlon), &
                           nc_count=(/ nlev,nlat,nlon,1 /), context=routine)
      !? ncount not needed?  Reading the whole field.

      ! 2023-11: ions do not have real or used data in their halos.
      !          Make this clear by leaving the halos filled with MISSING_R4
      !          TODO: Will this be translated into NetCDF missing_value?
      ! call add_halo_fulldom3d(fulldom3d)

      call filter_io_to_blocks(fulldom3d, varname, file_root, member)

   else
      ! TODO: error; varname is inconsistent with VT_ORIGININDX
   endif
enddo

deallocate(fulldom3d)
!, fulldom1d
   
end subroutine filter_to_restarts

!-----------------------------------------------------------------------
! Copy updated data from the full domain into the halo regions,
! in preparation for extracting haloed blocks into the block restart files.
! First the halos past the East and West edges are taken from the wrap-around points.
! The halos beyond the edge latitudes in the North and South 
! are taken by reaching over the pole to a longitude that's half way around the globe.
! This is independent of the number of blocks.

subroutine add_halo_fulldom3d(fulldom3d)

! Space for full domain field (read from filter_output.nc)
! and halo around the full domain
real(r4), intent(inout) :: fulldom3d(1:nz_per_block,       &
                                     1-nghost:nlat+nghost, &
                                     1-nghost:nlon+nghost)  

integer               :: g, i, j, haflat, haflon
real(r4), allocatable :: normed(:,:)
character(len=16)     :: debug_format

character(len=*), parameter :: routine = 'add_halo_fulldom3d'

! An array for debugging by renormalizing an altitude of fulldom3d.
allocate(normed(1-nghost:nlat+nghost, &
                1-nghost:nlon+nghost))

haflat = nlat / 2
haflon = nlon / 2

do g = 1,nghost
   ! left; reach around the date line.
   !         There's no data at the ends of the halos for this copy.
   fulldom3d  (:,1:nlat,     1-g) &
   = fulldom3d(:,1:nlat,nlon+1-g)

   ! right
   fulldom3d  (:,1:nlat,nlon+g) &
   = fulldom3d(:,1:nlat,g)

   ! bottom; reach over the S Pole for halo values.
   !         There is data at the ends of the halos for these.)

   fulldom3d  (:, 1-g ,1-nghost       :haflon) &
   = fulldom3d(:,   g ,1-nghost+haflon:nlon)
   fulldom3d  (:, 1-g ,haflon+1:nlon) &
   = fulldom3d(:,   g ,       1:haflon)
   ! Last 2 (halo) points on the right edge (at the bottom)
   fulldom3d  (:, 1-g ,  nlon+1:  nlon+nghost) &
   = fulldom3d(:,   g ,haflon+1:haflon+nghost)

   ! top
   fulldom3d  (:, nlat  +g, 1-nghost       :haflon) &
   = fulldom3d(:, nlat+1-g, 1-nghost+haflon:nlon)
   fulldom3d  (:, nlat  +g, haflon+1:nlon) &
   = fulldom3d(:, nlat+1-g,        1:haflon)
   ! Last 2 (halo) points on the right edge (at the top)
   fulldom3d  (:, nlat  +g,   nlon+1:  nlon+nghost) &
   = fulldom3d(:, nlat+1-g, haflon+1:haflon+nghost)
enddo

if (any(fulldom3d == MISSING_R4)) then
   error_string_1 = 'ERROR: some fulldom3d contain MISSING_R4 after halos'
   call error_handler(E_ERR, routine, error_string_1, source, revision, revdate)
endif

! TODO: Keep halo corners check for future use?
!       Add more robust rescaling.
! Debug; print the 4x4 arrays (corners & middle) 
! to see whether values are copied correctly
! Level 44 values range from 800-eps to 805.  I don't want to see the 80.
! For O+ range from 0 to 7e+11, but are close to 1.1082e+10 near the corners.
! 2023-12-20; Aaron sent new files with 54 levels.
if (debug >= 100 .and. do_output()) then
   if (fulldom3d(54,10,10) > 1.e+10) then
      normed = fulldom3d(54,:,:) - 1.1092e+10
      debug_format = '(3(4E10.4,2X))'
   else if (fulldom3d(54,10,10) < 1000._r4) then
      normed = fulldom3d(54,:,:) - 800._r4
      debug_format = '(3(4F10.5,2X))'
   endif
   
   ! Debug HDF5 
   write(error_string_1,'("normed_field(10,nlat+1,nlon+2) = ",3(1x,i5))') normed(nlat+1,nlon+2)
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
   
   ! 17 format debug_format
   print*,'top'
   do j = nlat+2, nlat-1, -1
      write(*,debug_format) (normed(j,i), i=      -1,       2), &
                            (normed(j,i), i=haflon-1,haflon+2), &
                            (normed(j,i), i=  nlon-1,  nlon+2)
   enddo
   print*,'middle'
   do j = haflat+2, haflat-1 , -1
      write(*,debug_format) (normed(j,i), i=      -1,       2), &
                            (normed(j,i), i=haflon-1,haflon+2), &
                            (normed(j,i), i=  nlon-1,  nlon+2)
   enddo
   print*,'bottom'
   do j = 2,-1, -1
      write(*,debug_format) (normed(j,i), i=      -1,       2), &
                            (normed(j,i), i=haflon-1,haflon+2), &
                            (normed(j,i), i=  nlon-1,  nlon+2)
   enddo
endif

deallocate(normed)

end subroutine add_halo_fulldom3d

!-----------------------------------------------------------------------
! Transfer part of the full field into a block restart file.

subroutine filter_io_to_blocks(fulldom3d, varname, file_root, member)

real(r4), intent(in) :: fulldom3d(1:nz_per_block,       &
                                  1-nghost:nlat+nghost, &
                                  1-nghost:nlon+nghost  )
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: file_root
integer,          intent(in) :: member

! Don't collect velocity components (6 of them)
!   real(r4) :: temp0d 
! , temp1d(:)   ?
integer            :: ncid_output
integer            :: ib, jb, nb
integer            :: starts(3), ends(3), xcount, ycount, zcount
character(len=256) :: block_file

character(len=*), parameter :: routine = 'filter_io_to_blocks'

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
! allocate(temp1d(1-nghost:max(nx_per_block,ny_per_block,nz_per_block)+nghost))


zcount = nz_per_block
ycount = ny_per_block + (2 * nghost)
xcount = nx_per_block + (2 * nghost)


if (debug > 0 .and. do_output()) then
   write(error_string_1,'(A,I0,A,I0,A)') 'Now putting the data for ', nblocks_lon, &
        ' blocks lon by ',nblocks_lat,' blocks lat'
   call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
end if

starts(1) = 1
ends(1)   = nz_per_block

do jb = 1, nblocks_lat
   starts(2) = (jb - 1) * ny_per_block - nghost + 1
   ends(2)   =  jb      * ny_per_block + nghost

   do ib = 1, nblocks_lon
      starts(3) = (ib - 1) * nx_per_block - nghost + 1
      ends(3)   =  ib      * nx_per_block + nghost

      nb = (jb - 1) * nblocks_lon + ib - 1

      block_file = block_file_name(trim(file_root), member, nb)
      ncid_output = open_block_file(block_file, 'readwrite')
   
      ! TODO: error checking; does the block file have the field in it?
      if ( debug > 0 .and. do_output()) then
        write(error_string_1,'(A,3(2X,i5))') "block, ib, jb = ", nb, ib, jb
        call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
        write(error_string_1,'(3(A,3i5))') &
             'starts = ',starts, 'ends = ',ends, '[xyz]counts = ',xcount,ycount,zcount
        call error_handler(E_MSG, routine, error_string_1, source, revision, revdate)
      endif      

      call nc_put_variable(ncid_output, trim(varname), &
           fulldom3d(starts(1):ends(1), starts(2):ends(2), starts(3):ends(3)), &
           context=routine, nc_count=(/ zcount,ycount,xcount /) )

      call nc_close_file(ncid_output)

   enddo
enddo

! 
! TODO: ? Add f107 and Rho to the restart files
!    call read_filter_io_block0d(ncid, ivals(1), data0d)
!    if (data0d < 0.0_r8) data0d = 60.0_r8 !alex
!    write(ounit) data0d

end subroutine filter_io_to_blocks

!-----------------------------------------------------------------------
! End of model_mod
!-----------------------------------------------------------------------
end module model_mod

