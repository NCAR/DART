! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

!-------------------------------------------------------------------------------
!
! Interface for HAO-TIEGCM
!
!-------------------------------------------------------------------------------

use        types_mod, only : r4, r8, i8, MISSING_R8, MISSING_R4, PI, &
                             earth_radius, gravity, obstypelength, MISSING_I

use time_manager_mod, only : time_type, set_calendar_type, set_time_missing,        &
                             set_time, get_time, print_time,                        &
                             set_date, get_date, print_date,                        &
                             operator(*),  operator(+), operator(-),                &
                             operator(>),  operator(<), operator(/),                &
                             operator(/=), operator(<=)

use     location_mod, only : location_type,                                         &
                             loc_get_close_obs => get_close_obs,                    &
                             loc_get_close_state => get_close_state,                &
                             set_location, get_location,                            &
                             get_dist, query_location,                              &
                             get_close_type, VERTISUNDEF,                           &
                             VERTISPRESSURE, VERTISHEIGHT, VERTISLEVEL,             &
                             vertical_localization_on, set_vertical

use    utilities_mod, only : file_exist, open_file, close_file, logfileunit,        &
                             error_handler, E_ERR, E_MSG, E_WARN, nmlfileunit,      &
                             do_output, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, nc_check, register_module,   &
                             file_to_text, find_textfile_dims, to_upper

use     obs_kind_mod, only : QTY_U_WIND_COMPONENT,           &
                             QTY_V_WIND_COMPONENT,           &
                             QTY_TEMPERATURE,                &! neutral temperature obs
                             QTY_PRESSURE,                   &! neutral pressure obs
                             QTY_MOLEC_OXYGEN_MIXING_RATIO,  &! neutral composition obs
                             QTY_1D_PARAMETER,               &
                             QTY_GEOPOTENTIAL_HEIGHT,        &
                             QTY_GEOMETRIC_HEIGHT,           &
                             QTY_VERTICAL_TEC,               &! total electron content
                             get_index_for_quantity

use mpi_utilities_mod,only : my_task_id

use default_model_mod, only : adv_1step,                                &
                              init_conditions => fail_init_conditions,  &
                              init_time => fail_init_time,              &
                              nc_write_model_vars,                      &
                              pert_model_copies

use state_structure_mod, only : add_domain, get_dart_vector_index, add_dimension_to_variable, &
                                finished_adding_domain, state_structure_info, &
                                get_domain_size, set_parameter_value, get_model_variable_indices, &
                                get_num_dims, get_dim_name, get_variable_name, &
                                get_varid_from_varname, get_num_varids_from_kind, &
                                get_varid_from_kind, get_varids_from_kind, &
                                hyperslice_domain, get_num_domains

use distributed_state_mod, only : get_state, get_state_array

use ensemble_manager_mod, only : ensemble_type

use netcdf_utilities_mod, only : nc_synchronize_file, nc_add_global_attribute, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_define_dimension, nc_end_define_mode, &
                                 nc_put_variable,nc_add_attribute_to_variable, &
                                 nc_define_real_variable, nc_define_character_variable

use typesizes  !HK do we need these with netcdf_utilities_mod?
use netcdf

implicit none
private

!DART mandatory public interfaces
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          get_close_obs,          &
          get_close_state,        &
          shortest_time_between_assimilations, &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          write_model_time

!DART pass through interfaces
public :: adv_1step,              &
          init_conditions,        &
          init_time,              &
          pert_model_copies

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = 'tiegcm/model_mod.f90'
character(len=32 ), parameter :: revision = ''
character(len=128), parameter :: revdate  = ''

!-------------------------------------------------------------------------------
! namelist with default values

! IMPORTANT: Change output file names in tiegcm.nml to match these names
! i.e.  OUTPUT='tiegcm_restart_p.nc'
!       SECOUT='tiegcm_s.nc'
character(len=256) :: tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
character(len=256) :: tiegcm_secondary_file_name = 'tiegcm_s.nc'
character(len=256) :: tiegcm_namelist_file_name  = 'tiegcm.nml'
integer            :: debug = 0
logical            :: estimate_f10_7 = .false.
logical            :: initialize_f10_7 = .false.
integer            :: assimilation_period_seconds = 3600

integer, parameter :: MAX_NUM_VARIABLES = 30
integer, parameter :: MAX_NUM_COLUMNS = 6
character(len=NF90_MAX_NAME) :: variables(MAX_NUM_VARIABLES * MAX_NUM_COLUMNS) = ' '

namelist /model_nml/ tiegcm_restart_file_name, &
                     tiegcm_secondary_file_name, tiegcm_namelist_file_name, &
                     variables, debug, estimate_f10_7, initialize_f10_7, &
                     assimilation_period_seconds
                     

!-------------------------------------------------------------------------------
! define model parameters

integer                               :: nilev, nlev, nlon, nlat
real(r8),dimension(:),    allocatable :: lons, lats, levs, ilevs, plevs, pilevs
! HK are plevs, pilves per ensemble member?
real(r8)                              :: TIEGCM_reference_pressure
integer                               :: time_step_seconds
integer                               :: time_step_days
type(time_type)                       :: time_step

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_ORIGININDX   = 5 ! ... file of origin
integer, parameter :: VT_STATEINDX    = 6 ! ... update (state) or not

character(len=obstypelength) :: variable_table(MAX_NUM_VARIABLES, MAX_NUM_COLUMNS)

! include_vTEC = .true.  vTEC must be calculated from other vars
! include_vTEC = .false. just ignore vTEC altogether

!HK why are there two of these?
logical  :: include_vTEC = .true.
logical  :: include_vTEC_in_state = .false.

! IMPORTANT: 1 D model parameters (e.g., F107) are read in from "tiegcm.nml"
! (note "estimate_f10_7" option is still under
! development by Tomoko Matsuo as of June 24, 2011)

real(r8)        :: f10_7
type(time_type) :: state_time ! module-storage declaration of current model time

integer(i8)           :: model_size ! the state vector length
integer :: nfields  ! number of tiegcm variables in DART state
! global domain id to be used by routines in state_structure_mod
integer :: domain_id(3) ! restart, secondary, calculate
integer, parameter :: RESTART_DOM = 1
integer, parameter :: SECONDARY_DOM = 2
integer, parameter :: CONSTRUCT_DOM = 3

! lon and lat grid specs. 2.5 degree or 5 degree grid
real(r8)  :: bot_lon        = MISSING_R8
real(r8)  :: top_lon        = MISSING_R8
real(r8)  :: delta_lon      = MISSING_R8
real(r8)  :: bot_lat        = MISSING_R8
real(r8)  :: top_lat        = MISSING_R8
real(r8)  :: delta_lat      = MISSING_R8
integer   :: zero_lon_index = MISSING_I


! FOR NOW OBS LOCATIONS ARE EXPECTED GIVEN IN HEIGHT [m],
! AND SO VERTICAL LOCALIZATION COORDINATE IS *always* HEIGHT
! (note that gravity adjusted geopotential height (ZG)
!  read in from "tiegcm_s.nc" *WARNING* ZG is 'cm', DART is mks)
integer               :: ivarZG

character(len=512)    :: string1, string2, string3
logical, save         :: module_initialized = .false.

!===============================================================================
contains
!===============================================================================

subroutine static_init_model()
!-------------------------------------------------------------------------------
!

integer :: iunit, io
real(r8) :: total_steps

if (module_initialized) return ! only need to do this once

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the namelist entry for model_mod from input.nml
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

if (do_output()) then
   write(     *     ,*)'static_init_model: debug level is ',debug
   write(logfileunit,*)'static_init_model: debug level is ',debug
endif

! Read in TIEGCM namelist input file (just for definition)
! Read in TIEGCM grid definition etc from TIEGCM restart file
! Read in TIEGCM auxiliary variables from TIEGCM 'secondary' file

call read_TIEGCM_namelist(tiegcm_namelist_file_name)
call read_TIEGCM_definition(tiegcm_restart_file_name)

if ( estimate_f10_7 ) then
  if (initialize_f10_7) then
   write(string1, '(f10.5)') f10_7
   call error_handler(E_MSG, 'initalizing f10_7 as', string1 )
  else
   call error_handler(E_MSG, 'reading f10_7 estimate', 'from netcdf file')
  endif
endif

! error-check, convert namelist input to variable_table, and build the
! state structure
call verify_variables()

call set_calendar_type('Gregorian')

! Convert the last year/day/hour/minute to a dart time.
state_time = read_model_time(tiegcm_restart_file_name)

! Ensure assimilation_period is a multiple of the dynamical timestep
! The time is communicated to TIEGCM through their "STOP" variable,
! which is an array of length 3 corresponding to day-of-year, hour, minute
! SO - there is some combination of 'STEP' and assimilation_period_seconds
! that must be an integer number of minutes.

time_step_seconds = time_step_seconds + time_step_days*86400

if (assimilation_period_seconds < time_step_seconds) then
   write(string1,*)'assimilation_period_seconds must be >= STEP'
   write(string2,*)' input.nml: assimilation_period_seconds ',assimilation_period_seconds
   write(string3,*)'tiegcm.nml: STEP ',time_step_seconds
   call error_handler(E_ERR,'static_init_model',string1, &
              source, revision, revdate, text2=string2,text3=string3)
endif

total_steps = real(assimilation_period_seconds,r8)/real(time_step_seconds,r8)

if ( time_step_seconds*nint(total_steps) /= assimilation_period_seconds) then
   write(string1,*)'assimilation_period_seconds must be an integer number of tiegcm "STEP"s'
   write(string2,*)' input.nml: assimilation_period_seconds ',assimilation_period_seconds
   write(string3,*)'tiegcm.nml: STEP ',time_step_seconds
   call error_handler(E_ERR,'static_init_model',string1, &
              source, revision, revdate, text2=string2,text3=string3)
endif

if ( mod(assimilation_period_seconds,60) /= 0 ) then
   write(string1,*)'assimilation_period_seconds must be an integer number of tiegcm "STEP"s'
   write(string2,*)'assimilation_period_seconds=',assimilation_period_seconds, &
                   ' STEP=',time_step_seconds
   write(string3,*)'AND must be an integer number of minutes because of tiegcm "STOP"'
   call error_handler(E_ERR,'static_init_model',string1, &
              source, revision, revdate, text2=string2,text3=string3)
endif

time_step = set_time(assimilation_period_seconds, 0)

end subroutine static_init_model


!-------------------------------------------------------------------------------

function get_model_size()
! Returns the size of the model as an integer. Required for all
! applications.

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



!-------------------------------------------------------------------------------


subroutine model_interpolate(state_handle, ens_size, location, iqty, obs_val, istatus)
! Given  a location, and a model state variable qty,
! interpolates the state variable field to that location and returns
! the value in obs_val for each ensemble member.
! The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The iqty variable
! is a integer that specifies the qty of field (for
! instance temperature, zonal wind component, etc.).

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: iqty
real(r8),           intent(out) :: obs_val(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer  :: which_vert
integer  :: lat_below, lat_above, lon_below, lon_above ! these are indices
integer  :: zero_lon_index
real(r8) :: lon_fract, lat_fract
real(r8) :: lon, lat, lon_lat_lev(3)
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real(r8), dimension(ens_size) :: val11, val12, val21, val22
real(r8) :: height
integer  :: level, bogus_level
integer  :: dom_id, var_id

if ( .not. module_initialized ) call static_init_model

! Default for failure return
istatus(:) = 1
obs_val(:) = MISSING_R8

! Failure codes
! 11 QTY_GEOPOTENTIAL_HEIGHT is unsupported
! 22 unsupported veritcal coordinate
! 33 level given < or > model levels
! 44 quantity not part of the state
! 55 outside state (can not extrapolate above or below)
! 66 unknown vertical stagger

! GITM uses a vtec routine in obs_def_upper_atm_mod:get_expected_gnd_gps_vtec()
! TIEGCM has its own vtec routine, so we should use it. This next block ensures that.
! The get_expected_gnd_gps_vtec() tries to interpolate QTY_GEOPOTENTIAL_HEIGHT
! when it does, this will kill it. 

if ( iqty == QTY_GEOPOTENTIAL_HEIGHT ) then
   istatus(:) = 11
   write(string1,*)'QTY_GEOPOTENTIAL_HEIGHT currently unsupported'
   call error_handler(E_ERR,'model_interpolate',string1,source, revision, revdate)
endif


! Get the position
lon_lat_lev = get_location(location)
lon         = lon_lat_lev(1) ! degree
lat         = lon_lat_lev(2) ! degree
height      = lon_lat_lev(3) ! level (int) or height (real)
level       = int(lon_lat_lev(3))  ! HK can I just convert lon_lat_lev to int always?


which_vert = nint(query_location(location))

call compute_bracketing_lat_indices(lat, lat_below, lat_above, lat_fract)
call compute_bracketing_lon_indices(lon, lon_below, lon_above, lon_fract)

! Pressure is not part of the state vector
! pressure is static data on plevs/pilevs
if ( iqty == QTY_PRESSURE) then
   if (which_vert == VERTISLEVEL) then
      ! Some variables need plevs, some need pilevs
      ! We only need the height (aka level)
      ! the obs_def_upper_atm_mod.f90:get_expected_O_N2_ratio routines queries
      ! for the pressure at the model levels - EXACTLY - so ...
      ! FIXME ... at present ... the only time model_interpolate
      ! gets called with QTY_PRESSURE is to calculate density, which
      ! requires other variables that only live on the midpoints.
      ! I cannot figure out how to generically decide when to
      ! use plevs vs. pilevs

      ! Check to make sure vertical level is possible.
      if ((level < 1) .or. (level > nlev)) then
         istatus(:) = 33
         return
      else
         obs_val(:) = plevs(level)
         istatus(:) = 0
      endif
   elseif (which_vert == VERTISHEIGHT) then

      ! If PRESSURE; calculate the pressure from several variables. HK I don't understand this.
      ! vert_interp() interpolates the state column to
      ! the same vertical height as the observation.
      ! THEN, it is the same as the 2D case.

      ! FIXME ... is it possible to try to get a pressure with which_vert == undefined
      ! At present, vert_interp will simply fail because height is a negative number.
      ! HK what are you supposed to do for pressure with VERTISUNDEF? level 1?

      call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_below, height, iqty, val11, istatus)
      if (any(istatus /= 0)) return  !HK bail at the first failure
      call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_above, height, iqty, val12, istatus)
      if (any(istatus /= 0)) return
      call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_below, height, iqty, val21, istatus)
      if (any(istatus /= 0)) return
      call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_above, height, iqty, val22, istatus)
      obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
   else

      write(string1,*) 'vertical coordinate type:',which_vert,' cannot be handled'
      call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)

   endif ! which vert

   return

endif ! end of QTY_PRESSURE

! check if qty is in the state vector
call find_qty_in_state(iqty, dom_id, var_id)
if (dom_id < 0 ) then
   istatus(:) = 44
   return
endif

if( which_vert == VERTISHEIGHT ) then

   ! vert_interp() interpolates the state column to
   ! the same vertical height as the observation.
   ! THEN, it is the same as the 2D case.

   call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_below, height, iqty, val11, istatus)
   if (any(istatus /= 0)) return  !HK bail at the first failure
   call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_above, height, iqty, val12, istatus)
   if (any(istatus /= 0)) return
   call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_below, height, iqty, val21, istatus)
   if (any(istatus /= 0)) return
   call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_above, height, iqty, val22, istatus)
   obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
   istatus = 0
elseif( which_vert == VERTISLEVEL) then
   ! Check to make sure vertical level is possible.
   if ((level < 1) .or. (level > nilev)) then
     istatus(:) = 33
     return
   endif

   ! one use of model_interpolate is to allow other modules/routines
   ! the ability to 'count' the model levels. To do this, create observations
   ! with locations on model levels and 'interpolate' for QTY_GEOMETRIC_HEIGHT.
   ! When the interpolation fails, you've gone one level too far.
   ! HK why does it have to be QTY_GEOMETRIC_HEIGHT?

   val11(:) = get_state(get_dart_vector_index(lon_below, lat_below, level, domain_id(dom_id), var_id ), state_handle)
   val12(:) = get_state(get_dart_vector_index(lon_below, lat_above, level, domain_id(dom_id), var_id ), state_handle)
   val21(:) = get_state(get_dart_vector_index(lon_above, lat_below, level, domain_id(dom_id), var_id ), state_handle)
   val22(:) = get_state(get_dart_vector_index(lon_above, lat_above, level, domain_id(dom_id), var_id ), state_handle)
   obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
   istatus = 0

elseif( which_vert == VERTISUNDEF) then
   bogus_level  = 1  !HK what should this be?  Do only 2D fields have VERTISUNDEF?
   val11(:) = get_state(get_dart_vector_index(lon_below, lat_below, bogus_level, domain_id(dom_id), var_id), state_handle)
   val12(:) = get_state(get_dart_vector_index(lon_below, lat_above, bogus_level, domain_id(dom_id), var_id), state_handle)
   val21(:) = get_state(get_dart_vector_index(lon_above, lat_below, bogus_level, domain_id(dom_id), var_id), state_handle)
   val22(:) = get_state(get_dart_vector_index(lon_above, lat_above, bogus_level, domain_id(dom_id), var_id), state_handle)
   obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
   istatus(:) = 0

else

   write(string1,*) 'vertical coordinate type:',which_vert,' cannot be handled'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)

endif

end subroutine model_interpolate

!-------------------------------------------------------------------------------
function shortest_time_between_assimilations()
type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations

!-------------------------------------------------------------------------------


subroutine get_state_meta_data(index_in, location, var_qty)
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_qty

integer  :: remainder
integer  :: relindx, absindx
integer  :: lon_index, lat_index, lev_index
integer  :: local_qty, var_id, dom_id
integer  :: seconds, days ! for f10.7 location
real(r8) :: longitude ! for f10.7 location
character(len=NF90_MAX_NAME) :: dim_name
integer  :: d

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, lev_index, var_id=var_id, dom_id=dom_id, kind_index=local_qty)

!print*, 'HK index_in, lon_index, lat_index, lev_index', index_in, lon_index, lat_index, lev_index, var_id

if(present(var_qty)) var_qty = local_qty

!HK check for f10.7 by varname?
if (get_variable_name(dom_id, var_id) == 'f10_7') then
   ! f10_7 is most accurately located at local noon at equator.
   ! 360.0 degrees in 86400 seconds, 43200 secs == 12:00 UTC == longitude 0.0

   call get_time(state_time, seconds, days)
   longitude = 360.0_r8 * real(seconds,r8) / 86400.0_r8 - 180.0_r8
   if (longitude < 0.0_r8) longitude = longitude + 360.0_r8
   location = set_location(longitude, 0.0_r8,  400000.0_r8, VERTISUNDEF)
   return
end if

! search for either ilev or lev
dim_name = ilev_or_lev(dom_id, var_id)

select case (trim(dim_name))
   case ('ilev')
      location  = set_location(lons(lon_index), lats(lat_index), ilevs(lev_index), VERTISLEVEL)
   case ('lev') ! height on midpoint
      location  = set_location(lons(lon_index), lats(lat_index), levs(lev_index), VERTISLEVEL)
   case default
    call error_handler(E_ERR, 'get_state_meta_data', 'expecting ilev or ilat dimension')
    ! HK TODO 2D variables.
end select

end subroutine get_state_meta_data


!-------------------------------------------------------------------------------


subroutine end_model()
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

end subroutine end_model


!-------------------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file.
! What should go in here?

subroutine nc_write_model_atts( ncid, dom_id ) 

integer, intent(in)  :: ncid      ! netCDF file identifier
integer, intent(in)  :: dom_id

!-------------------------------------------------------------------------------
! variables for the namelist output
!-------------------------------------------------------------------------------

character(len=70), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_tiegcm_namelist

character(len=*), parameter :: routine = 'nc_write_model_atts'

if ( .not. module_initialized ) call static_init_model

! Write Global Attributes
call nc_begin_define_mode(ncid, routine)

call nc_add_global_creation_time(ncid, routine)

call nc_add_global_attribute(ncid, "model_source", source, routine)
call nc_add_global_attribute(ncid, "model", "TIEGCM", routine)


! define grid dimensions
call nc_define_dimension(ncid, 'lon',  nlon,  routine)
call nc_define_dimension(ncid, 'lat',  nlat,  routine)
call nc_define_dimension(ncid, 'lev',  nlev,  routine)
call nc_define_dimension(ncid, 'ilev', nilev, routine)

! define grid variables
! longitude
call nc_define_real_variable(     ncid, 'lon', (/ 'lon' /), routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'long_name', 'geographic longitude (-west, +east)',  routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'units', 'degrees_east', routine)

! latitude
call nc_define_real_variable(     ncid, 'lat', (/ 'lat' /),  routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'long_name', 'geographic latitude (-south, +north)', routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'units',     'degrees_north', routine)

! levs
call nc_define_real_variable(     ncid, 'lev', (/ 'lev' /), routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'long_name',      'midpoint levels', routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'short name',     'ln(p0/p)', routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'positive',       'up', routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'standard_name',  'atmosphere_ln_pressure_coordinate', routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'formula_terms',  'p0: p0 lev: lev', routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'formula',  'p(k) = p0 * exp(-lev(k))', routine)


! ilevs
call nc_define_real_variable(     ncid, 'ilev', (/ 'ilev' /), routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'long_name',      'interface levels', routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'short name',     'ln(p0/p)', routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'positive',       'up', routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'standard_name',  'atmosphere_ln_pressure_coordinate', routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'formula_terms',  'p0: p0 lev: ilev', routine)
call nc_add_attribute_to_variable(ncid, 'lev',  'formula',         'p(k) = p0 * exp(-ilev(k))', routine)


!-------------------------------------------------------------------------------
! Determine shape of namelist.
! long lines are truncated when read into textblock
!-------------------------------------------------------------------------------

call find_textfile_dims(tiegcm_namelist_file_name, nlines, linelen)
if (nlines > 0) then
   has_tiegcm_namelist = .true.
   allocate(textblock(nlines))
   textblock = ''
   call nc_define_dimension(ncid, 'tiegcmNMLnlines',  nlines,  routine)
   call nc_define_dimension(ncid, 'linelen',  len(textblock(1)),  routine)
   call nc_define_character_variable(ncid, 'tiegcm_nml', (/ 'linelen        ', 'tiegcmNMLnlines' /), routine)
   call nc_add_attribute_to_variable(ncid, 'tiegcm_nml', 'long_name', &
         'contents of '//trim(tiegcm_namelist_file_name), routine)
else
   has_tiegcm_namelist = .false.
endif

call nc_end_define_mode(ncid, routine)

!-------------------------------------------------------------------------------
! Write variables
!-------------------------------------------------------------------------------

! Fill in the coordinate variables

where (lons >= 180.0_r8) lons = lons - 360.0_r8
call nc_put_variable(ncid, 'lon',  lons,  routine)
call nc_put_variable(ncid, 'lat',  lats,  routine)
call nc_put_variable(ncid, 'lev',  levs,   routine)
call nc_put_variable(ncid, 'ilev', ilevs,  routine)

! Fill tiegcm in namelist variable
if (has_tiegcm_namelist) then
   call file_to_text(tiegcm_namelist_file_name, textblock)
   call nc_put_variable(ncid, 'tiegcm_nml', textblock, routine)
endif

! flush any pending i/o to disk
call nc_synchronize_file(ncid, routine)
if (has_tiegcm_namelist) deallocate(textblock)

end subroutine nc_write_model_atts


!-------------------------------------------------------------------------------
! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate.
! FOR NOW VERTICAL LOCALIZATION IS DONE ONLY IN HEIGHT (ZG)
! OBS VERTICAL LOCATION IS GIVEN IN HEIGHT (model_interpolate)
! STATE VERTICAL LOCATION IS GIVEN IN HEIGHT
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:)
integer(i8),                   intent(in)     :: loc_indx(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

integer :: k, q_ind
integer :: n
integer :: istatus

n = size(locs)

if (vertical_localization_on()) then ! need to get height
  call convert_vertical_state(state_handle, n, locs, loc_qtys, loc_indx, VERTISHEIGHT, istatus)  ! HK Do we care about istatus?
endif

call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                       num_close, close_ind, dist)

! Make the ZG part of the state vector far from everything so it does not get updated.
! Scroll through all the obs_loc(:) and obs_kind(:) elements

do k = 1,num_close
   q_ind  = close_ind(k)
   if (loc_qtys(q_ind) == QTY_GEOMETRIC_HEIGHT) then
      if (do_output() .and. (debug > 99)) then
         write(     *     ,*)'get_close_obs ZG distance is ', &
                     dist(k),' changing to ',10.0_r8 * PI
         write(logfileunit,*)'get_close_obs ZG distance is ', &
                     dist(k),' changing to ',10.0_r8 * PI
      endif
      dist(k) = 10.0_r8 * PI
   endif
enddo


if (estimate_f10_7) then
! f10_7 is given a location of latitude 0.0 and the longitude
! of local noon. By decreasing the distance from the observation
! to the dynamic f10_7 location we are allowing the already close
! observations to have a larger impact in the parameter estimation.
! 0.25 is heuristic. The 'close' observations have already been
! determined by the cutoff. Changing the distance here does not
! allow more observations to impact anything.
   do k = 1, num_close
      q_ind  = close_ind(k)
      if  (loc_qtys(q_ind) == QTY_1D_PARAMETER) then
         dist(k) = dist(k)*0.25_r8
      endif
   enddo
endif


end subroutine get_close_state

!-------------------------------------------------------------------------------


subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, state_handle)
! Given a DART ob (referred to as "base") and a set of obs priors or
! state variables returns the subset of close ones to the "base" ob, their
! indices, and their distances to the "base" ob...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate.
! FOR NOW VERTICAL LOCALIZATION IS DONE ONLY IN HEIGHT (ZG)
! OBS VERTICAL LOCATION IS GIVEN IN HEIGHT (model_interpolate)
! STATE VERTICAL LOCATION IS GIVEN IN HEIGHT

! Note that both base_obs_loc and obs_loc are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling
! routine. The calling routine is always filter_assim and these arrays are local
! arrays within filter_assim. In other words, these modifications will only
! matter within filter_assim, but will not propagate backwards to filter.

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                       num_close, close_ind, dist)

end subroutine get_close_obs

!-------------------------------------------------------------------------------

subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

integer  :: current_vert_type, i
real(r8) :: height(1)
integer  :: local_status(1)

character(len=*), parameter :: routine = 'convert_vertical_obs'

if ( which_vert == VERTISHEIGHT .or. which_vert == VERTISUNDEF) then
  status(:) = 0
  return
endif

! HK can I just call model_interpolate here?
do i = 1, num
   current_vert_type = nint(query_location(locs(i)))
   if (( current_vert_type == which_vert ) .or. &
       ( current_vert_type == VERTISUNDEF)) then
      status(i) = 0
      cycle
   endif

  call model_interpolate(state_handle, 1, locs(i), QTY_GEOMETRIC_HEIGHT, height, local_status )
  
  if (local_status(1) == 0) then
     call set_vertical(locs(i), height(1), VERTISHEIGHT)
  endif
  status(i) = local_status(1)

enddo

end subroutine convert_vertical_obs

!-------------------------------------------------------------------------------
! HK what can you localize in? height only?
subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

integer :: var_id, dom_id, lon_index, lat_index, lev_index
integer :: i, d
real(r8) :: height(1), height1(1), height2(1)
character(len=NF90_MAX_NAME) :: dim_name
integer(i8) :: height_idx


if  ( which_vert /= VERTISHEIGHT ) then
   call error_handler(E_ERR,'convert_vertical_state', 'only supports VERTISHEIGHT')
endif

istatus = 0 !HK what are you doing with this?

do i = 1, num

   call get_model_variable_indices(loc_indx(i), lon_index, lat_index, lev_index, var_id=var_id, dom_id=dom_id)
   
   ! search for either ilev or lev
   dim_name = ilev_or_lev(dom_id, var_id)
   
   select case (trim(dim_name))
      case ('ilev')
         height_idx = get_dart_vector_index(lon_index, lat_index, lev_index, &
                                            domain_id(SECONDARY_DOM), ivarZG)
         height = get_state(height_idx, state_handle)/100.0_r8
   
      case ('lev') ! height on midpoint
        height_idx = get_dart_vector_index(lon_index, lat_index, lev_index, &
                                   domain_id(SECONDARY_DOM), ivarZG)
        height1 = get_state(height_idx, state_handle)/100.0_r8
        height_idx = get_dart_vector_index(lon_index, lat_index, lev_index+1, &
                           domain_id(SECONDARY_DOM), ivarZG)
        height2 = get_state(height_idx, state_handle)/100.0_r8
        height = (height1 + height2) / 2.0_r8
   
      case default
       call error_handler(E_ERR, 'convert_vertical_state', 'expecting ilev or ilat dimension')
   end select
   
   locs(i) = set_location(lons(lon_index), lats(lat_index), height(1), VERTISHEIGHT)

end do

end subroutine convert_vertical_state

!-------------------------------------------------------------------------------

function read_model_time(filename, ncfileid, lasttime)
! Gets the latest time in the netCDF file.
character(len=*),  intent(in) :: filename
integer, optional, intent(in) :: ncfileid
integer, optional, intent(in) :: lasttime
type(time_type) :: read_model_time

integer :: ncid, TimeDimID, time_dimlen, DimID, dimlen, VarID

integer,  parameter         :: nmtime = 3
integer,  dimension(nmtime) :: mtime  ! day, hour, minute
integer                     :: year, doy, utsec
integer, allocatable, dimension(:,:) :: mtimetmp
integer, allocatable, dimension(:)   :: yeartmp

if ( present(ncfileid) ) then
   ncid = ncfileid
else
   call nc_check(nf90_open(filename, NF90_NOWRITE, ncid), &
              'read_model_time','open '//trim(filename))
endif

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
        'read_model_time', 'inquire id of time')
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
        'read_model_time', 'inquire_dimension time')
call nc_check(nf90_inq_dimid(ncid, 'mtimedim', DimID), &
        'read_model_time', 'inq_dimid mtimedim')
call nc_check(nf90_inquire_dimension(ncid,     DimID, len=dimlen), &
        'read_model_time', 'inquire_dimension mtimedim')

if (present(lasttime)) then
   if (lasttime /= time_dimlen) then
      write(string1, *) trim(filename), ' last time index is ', &
                        time_dimlen, ' desired ', lasttime
      call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
   endif
endif

if (dimlen /= nmtime) then
   write(string1, *) trim(filename), ' mtimedim = ',dimlen, ' DART expects ', nmtime
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

allocate(mtimetmp(dimlen, time_dimlen), yeartmp(time_dimlen))

!... get mtime
call nc_check(nf90_inq_varid(ncid, 'mtime', VarID), &
        'read_model_time', 'inquire id of time')
call nc_check(nf90_get_var(ncid, VarID, values=mtimetmp), &
        'read_model_time', 'get_var mtime')

!... get year
call nc_check(nf90_inq_varid(ncid, 'year', VarID), 'read_model_time', 'inq_varid year')
call nc_check(nf90_get_var(ncid, VarID, values=yeartmp), 'read_model_time', 'get_var year')

! pick off the latest/last
mtime = mtimetmp(:,time_dimlen)
year  = yeartmp(   time_dimlen)

deallocate(mtimetmp,yeartmp)

doy   =  mtime(1)
utsec = (mtime(2)*60 + mtime(3))*60
read_model_time = set_time(utsec, doy-1) + set_date(year, 1, 1)  ! Jan 1 of whatever year.

if (do_output()) then
   write(*,*) trim(filename)//':read_model_time: tiegcm [year, doy, hour, minute]', &
            year, mtime
   call print_date(read_model_time, str=trim(filename)//':read_model_time: date ')
   call print_time(read_model_time, str=trim(filename)//':read_model_time: time ')
endif

if ( .not. present(ncfileid) ) then
   call nc_check(nf90_close(ncid), 'read_model_time', 'close '//trim(filename))
endif

end function read_model_time

!-------------------------------------------------------------------------------
subroutine write_model_time(ncid, dart_time)

use typeSizes
use netcdf

integer,         intent(in) :: ncid
type(time_type), intent(in) :: dart_time

integer :: dim_ids(2), var_id, ret
integer :: year, month, day, hour, minute, second
character(len=19) :: timestring

end subroutine write_model_time

!===============================================================================
! Routines below here are private to the module
!===============================================================================


subroutine read_TIEGCM_namelist(file_name)
! Under certain situations, the value of f10.7 is a parameter to be estimated
! and needs to be added to state vector

character(len=*), intent(in) :: file_name
integer  :: iunit, io
integer  :: daysec = 86400

!-------------------------------------------------------------------------------
! 1/3/2011, the namelist definition taken from $TGCMROOT/tiegcm1.93/src/input.F
!           the following parameter values are from params.F
!           modify the namelist definition for future tiegcm updates

integer,parameter :: mxind_time = 500 ! max number of time-dependent solar index points
integer,parameter :: mxhvols = 100    ! max number of output history file
integer,parameter :: mxseries = 10    ! max number of time series for primary histories
integer,parameter :: mxfsech = 100    ! max number of fields on secondary histories

! Namelist user input variables:

character(len=80):: &
     &  label,           &! optional generic text label for this run
     &  tempdir,         &! temporary directory
     &  magvol,          &! file name or mss path to magnetic data file
     &  amievol           ! file or mss path of amie data file (optional)

! date and calday are no longer supported, and are replaced by start_day,
! start_year, and calendar_advance. Date and calday are retained here so
! error usage statements can be issued if user sets one of them.

integer :: &
     &  start_day,       &! starting day of year (integer 0->365)
     &  start_year,      &! starting year (4-digit integer yyyy)
     &  calendar_advance,&! if > 0, advance calendar day from start_day
     &  date(3),         &! old: model starting year, day ( 2 ints yyyy,dd)
     &  calday,          &! old: starting calendar day (0-mxday)
     &  mxday,           &! calendar day (0-mxday)
     &  step,            &! model time step (integer seconds)
     &  dispose,         &! dispose output files to mss if dispose==1 or 2
     &  eddy_dif,        &! 0/1 flag for DOY dependent eddy diffusion (difk, dift, xmue)
     &  dynamo,          &! 0/1 flag for dynamo
     &  tideann,         &! 0/1 flag for annual tide (deprecated as of May 2008)
     &  aurora,          &! 0/1 flag for aurora
     &  ntask_lat,       &! number of tasks in latitude  dimension
     &  ntask_lon         ! number of tasks in longitude dimension
real :: &
     &  tide(10),        &! semidiurnal tide amplitudes and phases
     &  tide2(2),        &! diurnal tide amplitude and phase
     &  tide3m3(2),      &! 2-day wave amplitude and phase
     &  f107 = MISSING_R4,            &! 10.7 cm daily solar flux
     &  f107a = MISSING_R4,           &! 10.7 cm average (81-day) solar flux
     &  colfac,          &! collision factor
     &  amie_ibkg         ! AMIE_IBKG (not sure...)
!
! Input parameters that can be either constant or time-dependent:
real :: &
     &  power,           &! hemispheric power (gw) (hpower on histories)
     &  ctpoten,         &! cross-cap potential (volts)
     &  bximf,           &! BX component of IMF
     &  byimf,           &! BY component of IMF
     &  bzimf,           &! BZ component of IMF in nT
     &  swvel,           &! Solar wind velocity in km/s
     &  swden,           &! Solar wind density in #/cm3
     &  al,              &! AL lower magnetic auroral activity index in nT
     &  kp                ! Kp index
real,dimension(4,mxind_time) :: power_time,ctpoten_time,           &
     &  bximf_time,byimf_time,bzimf_time,swvel_time,swden_time,al_time,  &
     &  kp_time
integer :: &
     &  ntimes_ctpoten,ntimes_power,ntimes_bximf,ntimes_byimf,           &
     &  ntimes_bzimf,ntimes_swden,ntimes_swvel,ntimes_al,ntimes_kp
logical :: aluse    ! logical to use AL in Weimer 2001 model or not

! Parameters as read from namelist:
real :: rd_power,rd_ctpoten,rd_f107,rd_f107a,rd_bximf,rd_byimf,    &
     &  rd_bzimf,rd_swvel,rd_swden,rd_kp
!
! If indices_interp==1, time-dependent indices (power_time, ctpoten_time, etc)
! will be interpolated to model time, otherwise they will change only
! when the given values change. This has no effect on indices given as constants.

integer :: indices_interp=1

! Import data file names:

integer,parameter :: mxlen_filename=80
character(len=mxlen_filename) ::                                   &
!
! 4/2/08 btf: Introducing Weimer 2005 model (wei05sc.F).
!             Retain ability to call either the 2001 or 2005 weimer models
!             for now, to facilitate comparison runs, so potential_model
!             can be either WEIMER01 or WEIMER05.
!
     &  potential_model,   &! electric potential model used
                            ! Values can be 'HEELIS', 'WEIMER', or 'NONE'
                            ! If absent, the default value is set to 'HEELIS'
     &  weimer_ncfile,     &! path to netcdf weimer01 coefficients file
     &  wei05sc_ncfile,    &! path to netcdf data files for weimer05 model
     &  gpi_ncfile,        &! mss path or file path to netcdf gpi data file
     &  ncep_ncfile,       &! ncep data file (time-gcm only)
     &  see_ncfile,        &! mss path or file path to netcdf SEE flux data file
     &  imf_ncfile,        &! mss path or disk file path to netcdf IMF data file
     &  gswm_mi_di_ncfile, &! gswm migrating diurnal data file
     &  gswm_mi_sdi_ncfile,&! gswm migrating semi-diurnal data file
     &  gswm_nm_di_ncfile, &! gswm non-migrating diurnal data file
     &  gswm_nm_sdi_ncfile,&! gswm non-migrating semi-diurnal data file
     &  saber_ncfile,      &! SABER data (T,Z)
     &  tidi_ncfile,       &! TIDI data (U,V)
     &  seeflux,           &! SEE measured solar flux spectrum
     &  amienh,            &! Northern hemisphere AMIE input
     &  amiesh              ! Southern hemisphere AMIE input
!
!     integer,parameter :: ngpivars = 4
!     real :: gpi_vars(ngpivars) ! f107,f107a,power,ctpoten
!     character(len=16) ::
!    |  gpi_names(ngpivars)      ! names of gpi_vars

! Primary history user input (dimension parameters are in params.h):
character(len=80) :: &
     &  source,            &! file containing source history (optional)
        output(mxhvols)     ! output file(s) (required)
integer ::           &
     &  source_start(3),   &! source history model time
     &  start(3,mxseries), &! primary history model start time(s)
     &  stop(3,mxseries),  &! primary history model stop time(s)
     &  hist(3,mxseries),  &! primary history disk write frequency
     &  save(3,mxseries),  &! primary history file save frequency
     &  mxhist_prim,       &! max number of histories per primary file
     &  msreten,           &! retention period for history files
     &  noutput             ! number of output files given
!
! Secondary history user input (dimension parameters are in params.h):
character(len=80) ::   &
     &  secsource,           &! file containing source sec_history (for mhd)
     &  secout(mxhvols)       ! secondary history output file(s)
character(len=16) ::   &
     &  secflds(mxfsech)      ! secondary history output fields
integer ::             &
     &  secstart(3,mxseries),&! secondary history model start time(s)
     &  secstop(3,mxseries), &! secondary history model stop time(s)
     &  sechist(3,mxseries), &! secondary history disk write frequency
     &  secsave(3,mxseries), &! secondary history file save frequency
     &  mxhist_sech,         &! max number of histories per secondary file
     &  sech_nbyte            ! 4 or 8: write real or double values to secondary file
!
! Namelist for read:
namelist/tgcm_input/                                        &
     &  label,tempdir,magvol,amievol,date,calday,step,dispose,    &
     &  source,source_start,output,start,stop,hist,save,          &
     &  secout,secstart,secstop,sechist,secsave,secflds,          &
     &  potential_model,eddy_dif,dynamo,tide,tide2,tide3m3,       &
     &  f107,f107a,power,ctpoten,bximf,byimf,bzimf,swvel,swden,al,&
     &  kp,colfac,tideann,aurora,gpi_ncfile,gswm_mi_di_ncfile,    &
     &  gswm_mi_sdi_ncfile,gswm_nm_di_ncfile,gswm_nm_sdi_ncfile,  &
     &  mxhist_prim,mxhist_sech,msreten,ntask_lat,ntask_lon,      &
     &  start_day,start_year,calendar_advance,see_ncfile,         &
     &  ctpoten_time,power_time,bximf_time,byimf_time,bzimf_time, &
     &  kp_time,al_time,swden_time,swvel_time,indices_interp,     &
     &  imf_ncfile,saber_ncfile,tidi_ncfile,sech_nbyte, amie_ibkg, &
     &  seeflux, amienh, amiesh 


!-------------------------------------------------------------------------------

if( .not. file_exist(file_name)) then
   write(string1,*) trim(file_name),' not available.'
   call error_handler(E_ERR,'read_TIEGCM_namelist',string1,source,revision,revdate)
endif

call error_handler(E_MSG,'read_TIEGCM_namelist:','reading namelist ['//trim(file_name)//']')

! Read the namelist entry tgcm_input from tiegcm.nml
! TJH It would be nice to read the namelist and skip all the ';' in column 1.
! Basically, we are just getting the value of f10.7 and saving it.

call find_namelist_in_file('tiegcm.nml', 'tgcm_input', iunit)
read(iunit, nml = tgcm_input, iostat = io)
call check_namelist_read(iunit, io, 'tgcm_input')

if (step >= daysec) then
    time_step_days    = int(step/daysec)
    time_step_seconds = mod(step,daysec)
else
    time_step_days    = 0
    time_step_seconds = step
endif

if (do_output() .and. (debug > 1)) then
   write(string1,*) '..  tiegcm time_step_days    is ',time_step_days
   write(string2,*) 'tiegcm time_step_seconds is ',time_step_seconds
   call error_handler(E_MSG,'read_TIEGCM_namelist:',string1,text2=string2)
endif

f10_7 = f107  ! save this in module storage

if (do_output() .and. (debug > 1)) then
   write(string1,*) '..  f107 from tiegcm.nml     is ',f107
   call error_handler(E_MSG,'read_TIEGCM_namelist:',string1)
endif

end subroutine read_TIEGCM_namelist


!-------------------------------------------------------------------------------


subroutine read_TIEGCM_definition(file_name)
! Read TIEGCM grid definition and Geopotential from a tiegcm restart file
! fills metadata storage variables:
! lons(:), nlon
! lats(:), nlat
! lev(:),  nlev
! ilev(:), nilev
! plevs(:)
! pilevs(:)
! Converts the tiegcm longitude (-+180) to (0 360)
! Sets the grid specs

character(len=*), intent(in) :: file_name
integer  :: ncid, VarID, DimID, TimeDimID
real(r8) :: p0

if( .not. file_exist(file_name)) then
   write(string1,*) trim(file_name),' not available.'
   call error_handler(E_ERR,'read_TIEGCM_definition',string1,source,revision,revdate)
endif

call error_handler(E_MSG,'read_TIEGCM_definition:','reading restart ['//trim(file_name)//']')

call nc_check(nf90_open(file_name, NF90_NOWRITE, ncid), &
       'read_TIEGCM_definition','open '//trim(file_name))

! Make sure time is the unlimited dimension

call nc_check(nf90_inquire(ncid, unlimitedDimId = TimeDimID), &
       'read_TIEGCM_definition', 'inquire id of unlimited dimension time')
call nc_check(nf90_inq_dimid(ncid, 'time', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid time')

if( TimeDimID /= DimID ) then
   write(string1,*) trim(file_name),' does not have the "time" dimension as unlimited.'
   write(string2,*) 'This is a fundamental requirement for DART/TIEGCM'
   call error_handler(E_ERR,'read_TIEGCM_definition',string1,source,revision,revdate)
endif

! longitude - TIEGCM uses values +/- 180, DART uses values [0,360]
! HK check this against "Compute bracketing lon indices"

call nc_check(nf90_inq_dimid(ncid, 'lon', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid lon')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nlon), 'read_TIEGCM_definition', &
                  'inquire_dimension lon')
allocate(lons(nlon))
call nc_check(nf90_inq_varid(ncid, 'lon', VarID), 'read_TIEGCM_definition', &
                  'inq_varid lon')
call nc_check(nf90_get_var(ncid, VarID, values=lons), 'read_TIEGCM_definition', &
                  'get_var lon')

where (lons < 0.0_r8) lons = lons + 360.0_r8

! latitiude

call nc_check(nf90_inq_dimid(ncid, 'lat', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid lat')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nlat), 'read_TIEGCM_definition', &
                  'inquire_dimension lat')
allocate(lats(nlat))
call nc_check(nf90_inq_varid(ncid, 'lat', VarID), 'read_TIEGCM_definition', &
                  'inq_varid lat')
call nc_check(nf90_get_var(ncid, VarID, values=lats), 'read_TIEGCM_definition', &
                  'get_var lat')

! pressure

call nc_check(nf90_inq_varid(ncid, 'p0', VarID), 'read_TIEGCM_definition', &
                  'inq_varid p0')
call nc_check(nf90_get_var(ncid, VarID, values=p0), 'read_TIEGCM_definition', &
                  'get_var p0')

TIEGCM_reference_pressure = p0

call nc_check(nf90_inq_dimid(ncid, 'lev', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid lev')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nlev), 'read_TIEGCM_definition', &
                  'inquire_dimension lev')
! top level is not viable. The lower boundary condition is stored in the top level
nlev = nlev - 1

allocate(levs(nlev), plevs(nlev))

call nc_check(nf90_inq_varid(ncid, 'lev', VarID), 'read_TIEGCM_definition', &
                  'inq_varid lev')
call nc_check(nf90_get_var(ncid, VarID, values=levs), 'read_TIEGCM_definition', &
                  'get_var lev')

plevs = p0 * exp(-levs) * 100.0_r8 ![Pa] = 100* [millibars] = 100* [hPa]

call nc_check(nf90_inq_dimid(ncid, 'ilev', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid ilev')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nilev), 'read_TIEGCM_definition', &
                  'inquire_dimension ilev')
! top level is not viable. The lower boundary condition is stored in the top level
nilev = nilev - 1
allocate(ilevs(nilev), pilevs(nilev))

call nc_check(nf90_inq_varid(ncid, 'ilev', VarID), 'read_TIEGCM_definition', &
                  'inq_varid ilev')
call nc_check(nf90_get_var(ncid, VarID, values=ilevs), 'read_TIEGCM_definition', &
                  'get_var ilev')

pilevs = p0 * exp(-ilevs) * 100.0_r8 ! [Pa] = 100* [millibars] = 100* [hPa]

if ((nlev+1) .ne. nilev) then
   write(string1,*) 'number of midpoints should be 1 less than number of interfaces.'
   write(string2,*) 'number of midpoints  is nlev  = ',nlev
   write(string3,*) 'number of interfaces is nilev = ',nilev
   call error_handler(E_MSG,'read_TIEGCM_definition', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif


! Get lon and lat grid specs !HK these are static
bot_lon        = lons(1)                         ! 180.
delta_lon      = abs((lons(1)-lons(2)))          ! 5. or 2.5
zero_lon_index = int(bot_lon/delta_lon) + 1      ! 37 or 73
top_lon        = lons(nlon)                      ! 175. or 177.5
bot_lat        = lats(1)                         !
top_lat        = lats(nlat)                      !
delta_lat      = abs((lats(1)-lats(2)))          !

end subroutine read_TIEGCM_definition

!-------------------------------------------------------------------------------
! Fill up the variable_table from the namelist item 'variables'
! The namelist item variables is where a user specifies
! which variables they want in the DART state:
! variable name, dart qty, clamping min, clamping max, origin file, update or not

subroutine verify_variables()

integer :: nfields_restart       ! number of variables from restart file
integer :: nfields_secondary     ! number of variables from secondary file
integer :: nfields_constructed   ! number of constructed state variables

integer  :: i, nrows, ncols
integer  :: index1, indexN, varsize
integer  :: ncid, ncid1, ncid2, ncerr, VarID, dimlen

!HK obstypelength vs NF90_MAX_NAME?
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: filename
character(len=NF90_MAX_NAME) :: state_or_aux

nrows = size(variable_table,1) !HK these are MAX_NUM_VARIABLES, MAX_NUM_COLUMNS
ncols = size(variable_table,2)

! Convert the (input) 1D array "variables" into a table with six columns.
! The number of rows in the table correspond to the number of variables in the
! DART state vector.
! Column 1 is the netCDF variable name.
! Column 2 is the corresponding DART kind.
! Column 3 is the minimum value ("NA" if there is none) Not Applicable
! Column 4 is the maximum value ("NA" if there is none) Not Applicable
! Column 5 is the file of origin 'restart' or 'secondary' or 'calculate'
! Column 6 is whether or not the variable should be updated in the restart file.

nfields = 0
nfields_restart = 0
nfields_secondary = 0
nfields_constructed = 0

ROWLOOP : do i = 1, nrows

   varname      = trim(variables(ncols*i - 5))
   dartstr      = trim(variables(ncols*i - 4))
   minvalstring = trim(variables(ncols*i - 3))
   maxvalstring = trim(variables(ncols*i - 2))
   filename     = trim(variables(ncols*i - 1))
   state_or_aux = trim(variables(ncols*i    ))

   call to_upper(filename)
   call to_upper(state_or_aux) ! update or not

   variable_table(i,VT_VARNAMEINDX) = trim(varname)
   variable_table(i,VT_KINDINDX)    = trim(dartstr)
   variable_table(i,VT_MINVALINDX)  = trim(minvalstring)
   variable_table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   variable_table(i,VT_ORIGININDX)  = trim(filename)
   variable_table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ((variable_table(i,1) == ' ') ) exit ROWLOOP

   ! Any other condition is an error.
   if ( any(variable_table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:variables not fully specified.'
      string2 = 'Must be 6 entries per variable, last known variable name is'
      string3 = trim(variable_table(i,1))
      call error_handler(E_ERR,'get_variables_in_domain',string1, &
          source,revision,revdate,text2=string2,text3=string3)
   endif

   nfields=nfields+1
   if (variable_table(i,VT_ORIGININDX) == 'RESTART') nfields_restart = nfields_restart+1
   if (variable_table(i,VT_ORIGININDX) == 'SECONDARY') nfields_secondary = nfields_secondary+1
   if (variable_table(i,VT_ORIGININDX) == 'CALCULATE') nfields_constructed = nfields_constructed+1

enddo ROWLOOP

! Make sure variable exists in netCDF file, as long as it is not vTEC.
! vTEC gets constructed.

if (varname == 'VTEC') include_vTEC_in_state = .true.


! Do we need to augment the state vector with the parameter to estimate?
if ( estimate_f10_7 ) then

   nfields = nfields + 1
   nfields_constructed = nfields_constructed+1
   variable_table(nfields,VT_VARNAMEINDX) = 'f10_7'
   variable_table(nfields,VT_KINDINDX)    = 'QTY_1D_PARAMETER'
   variable_table(nfields,VT_MINVALINDX)  = 'NA'
   variable_table(nfields,VT_MAXVALINDX)  = 'NA'
   variable_table(nfields,VT_ORIGININDX)  = 'CALCULATE'
   variable_table(nfields,VT_STATEINDX)   = 'UPDATE'

endif

if (include_vTEC_in_state) then
   ! FIXME ... check to make sure all required variables are part of the DART
   ! state vector. These should be part of the DART state vector :
   ! If this is the case, then we _could_ use a more standard obs_def approach
   ! and call model_interpolate to return the VTEC on demand. What we have now
   ! is basically 'cached' the forward obs operator and created VTEC. HOWEVER,
   ! if we are doing prior state-space inflation, this has an impact on the
   ! VTEC in the state vector. Also ... same for ZG ... this is not great.

   ! NE
   ! TI
   ! TE
   ! OP

endif

! Record the contents of the DART state vector
if (do_output() .and. (debug > 99)) then
   do i = 1,nfields
      write(*,'(''variable'',i4,'' is '',a12,1x,a32,4(1x,a20))') i, &
             trim(variable_table(i,1)), &
             trim(variable_table(i,2)), &
             trim(variable_table(i,3)), &
             trim(variable_table(i,4)), &
             trim(variable_table(i,5)), &
             trim(variable_table(i,6))
      write(logfileunit,'(''variable'',i4,'' is '',a12,1x,a32,4(1x,a20))') i, &
             trim(variable_table(i,1)), &
             trim(variable_table(i,2)), &
             trim(variable_table(i,3)), &
             trim(variable_table(i,4)), &
             trim(variable_table(i,5)), &
             trim(variable_table(i,6))
   enddo
endif


call load_up_state_structure_from_file('tiegcm_restart_p.nc', nfields_restart, 'RESTART', RESTART_DOM)
call load_up_state_structure_from_file('tiegcm_s.nc', nfields_secondary, 'SECONDARY', SECONDARY_DOM)
if (estimate_f10_7) then
   call load_up_calculated_variables(nfields_constructed, 'CALCULATE', CONSTRUCT_DOM)
   model_size = get_domain_size(RESTART_DOM) + get_domain_size(SECONDARY_DOM) &
                          + get_domain_size(CONSTRUCT_DOM)
else
   model_size = get_domain_size(RESTART_DOM) + get_domain_size(SECONDARY_DOM)
endif

! set ivar. ZG is in the secondary domain
ivarZG = get_varid_from_varname(domain_id(SECONDARY_DOM), 'ZG')

end subroutine verify_variables

!-------------------------------------------------------------------------------
! Adds a domain to the state structure from a netcdf file
! Called from verify_variables
subroutine load_up_state_structure_from_file(filename, nvar, domain_name, domain_num)

character(len=*), intent(in) :: filename ! filename to read from
integer,          intent(in) :: nvar ! number of variables in domain
character(len=*), intent(in) :: domain_name ! restart, secondary
integer,          intent(in) :: domain_num

integer :: i,j

character(len=NF90_MAX_NAME), allocatable :: var_names(:)
real(r8), allocatable :: clamp_vals(:,:)
integer, allocatable :: kind_list(:)
logical, allocatable :: update_list(:)


allocate(var_names(nvar), kind_list(nvar), &
     clamp_vals(nvar,2), update_list(nvar))

update_list(:) = .true. ! default to update state variable
clamp_vals(:,:) = MISSING_R8 ! default to no clamping

j = 0
do i = 1, nfields
   if (variable_table(i,VT_ORIGININDX) == trim(domain_name)) then
      j = j+1
      var_names(j) = variable_table(i, VT_VARNAMEINDX)
      kind_list(j) = get_index_for_quantity(variable_table(i, VT_KINDINDX))
      if (variable_table(i, VT_MINVALINDX) /= 'NA') then
         read(variable_table(i, VT_MINVALINDX), '(d16.8)') clamp_vals(j,1)
      endif
      if (variable_table(i, VT_MAXVALINDX) /= 'NA') then
        read(variable_table(i, VT_MAXVALINDX), '(d16.8)') clamp_vals(j,2)
      endif
      if (variable_table(i, VT_STATEINDX) == 'NO_COPY_BACK') then
         update_list(j) = .false.
      endif
   endif
enddo

domain_id(domain_num) = add_domain(filename, nvar, &
                          var_names, kind_list, clamp_vals, update_list)

! remove top level from all variables
call hyperslice_domain(domain_id(domain_num), 'lev', nlev)
call hyperslice_domain(domain_id(domain_num), 'ilev', nlev)

deallocate(var_names, kind_list, clamp_vals, update_list)

end subroutine load_up_state_structure_from_file
!-------------------------------------------------------------------------------
! calcualted variables do not have a netcdf file
! HK or do they? What happens for multiple assimilation cycles?
subroutine load_up_calculated_variables(nvar, domain_name, domain_num)

integer,          intent(in) :: nvar ! number of variables in domain
character(len=*), intent(in) :: domain_name ! restart, secondary
integer,          intent(in) :: domain_num

integer :: i,j

character(len=NF90_MAX_NAME), allocatable :: var_names(:)
real(r8), allocatable :: clamp_vals(:,:)
integer, allocatable :: kind_list(:)
logical, allocatable :: update_list(:)


allocate(var_names(nvar), kind_list(nvar), &
     clamp_vals(nvar,2), update_list(nvar))

update_list(:) = .true. ! default to update state variable
clamp_vals(:,:) = MISSING_R8 ! default to no clamping

j = 0
do i = 1, nfields
   if (variable_table(i,VT_ORIGININDX) == domain_name) then
      j = j+1
      var_names(j) = variable_table(i, VT_VARNAMEINDX)
      kind_list(j) = get_index_for_quantity(variable_table(i, VT_KINDINDX))
      if (variable_table(i, VT_MINVALINDX) /= 'NA') then
         read(variable_table(i, VT_MINVALINDX), '(d16.8)') clamp_vals(j,1)
      endif
      if (variable_table(i, VT_MAXVALINDX) /= 'NA') then
        read(variable_table(i, VT_MAXVALINDX), '(d16.8)') clamp_vals(j,2)
      endif
      if (variable_table(i, VT_STATEINDX) == 'NO_COPY_BACK') then
         update_list(j) = .false.
      endif
   endif
enddo

domain_id(domain_num) = add_domain(nvar, var_names, &
       kind_list, clamp_vals, update_list, init_parameter_estimate=initialize_f10_7)

do i = 1, nvar ! HK f10_7 and VTEC and anything else?
  if (var_names(i) == 'f10_7') then
    ! dimensions to match what? Lanai has single value f10_7
    call add_dimension_to_variable(domain_id(domain_num), i, 'parameter', 1)
    ! initialize to namelist value, or read from file
    if (initialize_f10_7) call set_parameter_value(domain_id(domain_num), i, f10_7)
  endif

  !HK I don't understand why VTEC is a state variable vs. calculating it on the fly.
  if(var_names(i) == 'VTEC') then
     call error_handler(E_ERR,'load_up_calculated_variables', 'VTEC not sure yet.')
  endif
enddo

call finished_adding_domain(domain_id(domain_num))

end subroutine load_up_calculated_variables

!-------------------------------------------------------------------------------


!HK subroutine create_vtec( ncid, last_time, vTEC)
!HK !
!HK ! Create the vTEC from constituents in the netCDF file.
!HK !
!HK 
!HK integer,                  intent(in)  :: ncid
!HK integer,                  intent(in)  :: last_time
!HK real(r8), dimension(:,:), intent(out) :: vTEC
!HK 
!HK real(r8), allocatable, dimension(:,:,:) :: NE, TI, TE
!HK real(r8), allocatable, dimension(:,:,:) :: NEm_extended, ZG_extended
!HK real(r8), allocatable, dimension(:,:)   :: GRAVITYtop, Tplasma, Hplasma
!HK real(r8), allocatable, dimension(:)     :: delta_ZG, NE_middle
!HK 
!HK real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! m^2 * kg / s^2 / K
!HK real(r8), PARAMETER :: omass      = 2.678e-26_r8 ! mass of atomic oxgen kg
!HK 
!HK real(r8) :: earth_radiusm
!HK integer  :: VarID, nlev10, j, k
!HK 
!HK allocate( NE(nlon,nlat,nilev), NEm_extended(nlon,nlat,nilev+10), &
!HK           ZG_extended(nlon,nlat,nilev+10))
!HK allocate( TI(nlon,nlat,nlev), TE(nlon,nlat,nlev) )
!HK allocate( GRAVITYtop(nlon,nlat), Tplasma(nlon,nlat), Hplasma(nlon,nlat) )
!HK allocate( delta_ZG(nlev+9), NE_middle(nlev+9) )
!HK 
!HK !... NE (interfaces)
!HK call nc_check(nf90_inq_varid(ncid, 'NE', VarID), 'create_vtec', 'inq_varid NE')
!HK call nc_check(nf90_get_var(ncid, VarID, values=NE,     &
!HK                    start = (/    1,    1,     1, last_time /),   &
!HK                    count = (/ nlon, nlat, nilev,         1 /)),&
!HK                    'create_vtec', 'get_var NE')
!HK 
!HK !... ZG (interfaces) already read into the module and converted to metres
!HK 
!HK !... TI (midpoints)
!HK call nc_check(nf90_inq_varid(ncid, 'TI', VarID), 'create_vtec', 'inq_varid TI')
!HK call nc_check(nf90_get_var(ncid, VarID, values=TI,     &
!HK                    start = (/    1,    1,    1, last_time /),   &
!HK                    count = (/ nlon, nlat, nlev,         1 /)), &
!HK                    'create_vtec', 'get_var TI')
!HK 
!HK !... TE (midpoints)
!HK call nc_check(nf90_inq_varid(ncid, 'TE', VarID), 'create_vtec', 'inq_varid TE')
!HK call nc_check(nf90_get_var(ncid, VarID, values=TE,     &
!HK                    start = (/    1,    1,    1, last_time /),   &
!HK                    count = (/ nlon, nlat, nlev,         1 /)), &
!HK                    'create_vtec', 'get_var TE')
!HK 
!HK ! Construct vTEC given the parts
!HK 
!HK earth_radiusm = earth_radius * 1000.0_r8 ! Convert earth_radius in km to m
!HK NE            = NE * 1.0e+6_r8           ! Convert NE in #/cm^3 to #/m^3
!HK 
!HK ! Gravity at the top layer
!HK GRAVITYtop(:,:) = gravity * (earth_radiusm / (earth_radiusm + ZG(:,:,nilev))) ** 2
!HK 
!HK ! Plasma Temperature
!HK Tplasma(:,:) = (TI(:,:,nlev-1) + TE(:,:,nlev-1)) / 2.0_r8
!HK 
!HK ! Compute plasma scale height
!HK Hplasma = (2.0_r8 * k_constant / omass ) * Tplasma / GRAVITYtop
!HK 
!HK ! NE is extrapolated to 10 more layers
!HK nlev10  = nlev + 10
!HK 
!HK  ZG_extended(:,:,1:nilev) = ZG
!HK NEm_extended(:,:,1:nilev) = NE
!HK 
!HK do j = nlev, nlev10
!HK    NEm_extended(:,:,j) = NEm_extended(:,:,j-1) * exp(-0.5_r8)
!HK     ZG_extended(:,:,j) =  ZG_extended(:,:,j-1) + Hplasma(:,:) / 2.0_r8
!HK enddo
!HK 
!HK ! finally calculate vTEC - one gridcell at a time.
!HK 
!HK do j = 1, nlat
!HK do k = 1, nlon
!HK     delta_ZG(1:(nlev10-1)) =  ZG_extended(k,j,2:nlev10) -  ZG_extended(k,j,1:(nlev10-1))
!HK    NE_middle(1:(nlev10-1)) = (NEm_extended(k,j,2:nlev10) + NEm_extended(k,j,1:(nlev10-1))) / 2.0_r8
!HK    vTEC(k,j) = sum(NE_middle * delta_ZG) * 1.0e-16_r8 ! Convert to TECU (1.0e+16 #/m^2)
!HK enddo
!HK enddo
!HK 
!HK deallocate( NE, NEm_extended, ZG_extended)
!HK deallocate( TI, TE )
!HK deallocate( GRAVITYtop, Tplasma, Hplasma )
!HK deallocate( delta_ZG, NE_middle )
!HK 
!HK end subroutine create_vtec


!-------------------------------------------------------------------------------

subroutine vert_interp(state_handle, n, dom_id, var_id, lon_index, lat_index, height, iqty, &
                       val, istatus)
! returns the value at an arbitrary height on an existing horizontal grid location.
! istatus == 0 is a 'good' return code.

type(ensemble_type), intent(in) :: state_handle
integer,          intent(in)  :: n ! ensemble_size
integer,          intent(in)  :: dom_id
integer,          intent(in)  :: var_id
integer,          intent(in)  :: lon_index
integer,          intent(in)  :: lat_index
real(r8),         intent(in)  :: height
integer,          intent(in)  :: iqty
real(r8),         intent(out) :: val(n)
integer,          intent(out) :: istatus(n)

logical :: is_pressure
character(len=NF90_MAX_NAME) :: vertstagger

! Presume the worst. Failure.
istatus    = 1
val        = MISSING_R8

is_pressure = (iqty == QTY_PRESSURE)
if (is_pressure) then
   vertstagger = 'ilev'  !HK what should this be? Does it matter?
else
   vertstagger = ilev_or_lev(dom_id, var_id)
endif

if (vertstagger == 'ilev') then
  call vert_interp_ilev(state_handle, height, n, lon_index, lat_index, is_pressure, &
                          dom_id, var_id, val, istatus)
elseif (vertstagger == 'lev') then
  call vert_interp_lev(state_handle, height, n, lon_index, lat_index, is_pressure, &
                          dom_id, var_id, val, istatus)
endif

end subroutine vert_interp

!-------------------------------------------------------------------------------
subroutine find_qty_in_state(iqty, which_dom, var_id)
! Finds the first variable of the appropriate DART KIND
!
! Note from FindVar_by_kind:
! FIXME There is some confusion about using the T-minus-1 variables
! in this construct. Both TN and TN_NM have the same dart_kind,
! so we use the first one ... but it is not guaranteed that TN
! must preceed TN_NM, for example.
! HK we can force X rather than X_MN

integer, intent(in)  :: iqty
integer, intent(out) :: which_dom
integer, intent(out) :: var_id

integer :: num_same_kind, id, k
integer, allocatable :: multiple_kinds(:), n
character(NF90_MAX_NAME) :: varname

which_dom = -1
var_id = -1

!HK can you ever have the same variable in restart and secondary?
!   I don't think so
do id = 1, get_num_domains() ! RESTART_DOM, SECONDARY_DOM, CONSTRUCT_DOM

   num_same_kind = get_num_varids_from_kind(domain_id(id), iqty)
   if (num_same_kind == 0 ) cycle
   if (num_same_kind  > 1 ) then ! need to pick which one you want
     which_dom = id
     allocate(multiple_kinds(num_same_kind))
     call get_varids_from_kind(domain_id(id), iqty, multiple_kinds)
     do k = 1, num_same_kind
       varname = adjustl(get_variable_name(domain_id(id), multiple_kinds(k)))
       n = len(trim(varname))
       if (n <= 2) then ! variable name can not be X_MN
          var_id = multiple_kinds(k)
          exit
       elseif (trim(varname(n-2:n)) == '_NM') then ! variable name is _MN
          cycle ! assuming we want the X, not the X_MN
       else
         var_id = multiple_kinds(k)
         exit
       endif
     enddo
     deallocate(multiple_kinds)
   else !
      which_dom = id
      var_id = get_varid_from_kind(domain_id(id), iqty)
   endif
enddo

end subroutine find_qty_in_state

!-------------------------------------------------------------------------------
! find enclosing lon indices
! Compute bracketing lon indices:
! TIEGCM [-180 175]  DART [180, 185, ..., 355, 0, 5, ..., 175]
! HK in read_TIEGCM_definition the longitudes are converted to 0 to 360
subroutine compute_bracketing_lon_indices(lon, idx_below, idx_above, fraction)

real(r8), intent(in)  :: lon ! longitude
integer,  intent(out) :: idx_below, idx_above ! index in lons()
real(r8), intent(out) :: fraction ! fraction to use for interpolation

if(lon > top_lon .and. lon < bot_lon) then     ! at wraparound point [175 < lon < 180]
   idx_below = nlon
   idx_above = 1
   fraction = (lon - top_lon) / delta_lon
elseif (lon >= bot_lon) then                  ! [180 <= lon <= 360]
   idx_below = int((lon - bot_lon) / delta_lon) + 1
   idx_above = idx_below + 1
   fraction = (lon - lons(idx_below)) / delta_lon
else                                           ! [0 <= lon <= 175 ]
   idx_below = int((lon - 0.0_r8) / delta_lon) + zero_lon_index
   idx_above = idx_below + 1
   fraction = (lon - lons(idx_below)) / delta_lon
endif


end subroutine compute_bracketing_lon_indices

!-------------------------------------------------------------------------------
! on ilev
subroutine vert_interp_ilev(state_handle, height, n, lon_index, lat_index, is_pressure, &
                                    dom_id, var_id, val, istatus)

type(ensemble_type), intent(in) :: state_handle
real(r8), intent(in)  :: height
integer,  intent(in)  :: n ! ensemble size
integer,  intent(in)  :: lon_index
integer,  intent(in)  :: lat_index
logical,  intent(in)  :: is_pressure
integer,  intent(in)  :: dom_id, var_id
real(r8), intent(out) :: val(n) ! interpolated value
integer,  intent(out) :: istatus(n)

integer :: lev_bottom(n)
integer :: lev_top(n)
real(r8) :: frac_lev(n)
integer  :: k, i
real(r8) :: zgrid(n), delta_z(n), z2(n), zgrid_top(n), zgrid_bottom(n)
logical  :: found(n) ! track which ensemble members have been located
real(r8) :: val_top(n), val_bottom(n)
integer(i8)  :: indx_top(n), indx_bottom(n) ! state vector indice

istatus    = 1
frac_lev   = MISSING_R8
lev_top    = -1
lev_bottom = -1
found = .false.

   zgrid_bottom(:) = get_state(get_dart_vector_index(lon_index,lat_index,1, &
                               domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8
   zgrid_top(:) = get_state(get_dart_vector_index(lon_index,lat_index, nilev, &
                               domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8

   ! cannot extrapolate below bottom or beyond top
   do i = 1, n
      if ((zgrid_bottom(i) > height) .or. (zgrid_top(i) < height)) then
        istatus(i) = 55
      endif
   enddo
   if (any(istatus == 55)) return ! fail if any ensemble member fails

   ! Figure out what level is above/below, and by how much
   h_loop_interface : do k = 2, nilev

      zgrid(:) = get_state(get_dart_vector_index(lon_index,lat_index,k, &
                            domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8

      ! per ensemble member
      do i = 1, n
         if (found(i)) cycle
         if (height <= zgrid(i)) then
           found(i) = .true.
           lev_top(i)    = k
           lev_bottom(i) = lev_top(i) - 1
           if (all(found)) exit h_loop_interface
         endif
      enddo

   enddo h_loop_interface

   do i = 1, n
     indx_top(i) = get_dart_vector_index(lon_index,lat_index,lev_top(i), domain_id(SECONDARY_DOM), ivarZG)
     indx_bottom(i) = get_dart_vector_index(lon_index,lat_index,lev_bottom(i), domain_id(SECONDARY_DOM), ivarZG)
   enddo

   call get_state_array(zgrid(:), indx_top(:), state_handle)

   call get_state_array(z2(:), indx_bottom(:), state_handle)

   delta_z(:) = zgrid(:) - z2(:)
   frac_lev(:) = (zgrid(:) - height)/delta_z(:)


   if (is_pressure) then ! get fom plevs (pilevs?) array TODO Lanai is always plves

      val_top(:)    = plevs(lev_top(:))     !pressure at midpoint [Pa]
      val_bottom(:) = plevs(lev_bottom(:))  !pressure at midpoint [Pa]
      val(:)        = exp(frac_lev(:) * log(val_bottom(:)) + (1.0 - frac_lev(:)) * log(val_top(:)))
  
   else  ! get from state vector

      do i = 1, n
        indx_top(:) = get_dart_vector_index(lon_index,lat_index,lev_top(i), dom_id, var_id)
        indx_bottom(i) = get_dart_vector_index(lon_index,lat_index,lev_bottom(i), dom_id, var_id)
      enddo

      call get_state_array(val_top, indx_top(:), state_handle)
      call get_state_array(val_bottom, indx_bottom(:), state_handle)

      val(:) = frac_lev(:) * val_bottom(:)  + (1.0 - frac_lev(:)) * val_top(:)

   endif

  istatus(:) = 0

end subroutine vert_interp_ilev

!-------------------------------------------------------------------------------
! on lev
subroutine vert_interp_lev(state_handle, height, n, lon_index, lat_index, is_pressure, &
                                     dom_id, var_id, val, istatus)

type(ensemble_type), intent(in) :: state_handle
real(r8), intent(in)  :: height
integer,  intent(in)  :: n ! ensemble size
integer,  intent(in)  :: lon_index
integer,  intent(in)  :: lat_index
logical,  intent(in)  :: is_pressure
integer,  intent(in)  :: dom_id, var_id
real(r8), intent(out) :: val(n)  ! interpolated value
integer,  intent(out) :: istatus(n)

integer :: lev_bottom(n)
integer :: lev_top(n)
real(r8) :: frac_lev(n)

integer  :: k, i
real(r8) :: delta_z(n)
real(r8) :: zgrid_upper(n), zgrid_lower(n)
integer(i8)  :: indx_top(n), indx_bottom(n) ! state vector indices
logical  :: found(n) ! track which ensemble members have been located
real(r8) :: val_top(n), val_bottom(n)

istatus    = 1
frac_lev   = MISSING_R8
lev_top    = -1
lev_bottom = -1
found = .false.

   ! Variable is on level midpoints, not ilevels.
   ! Get height as the average of the ilevels.

   ! ilev index    1      2      3      4    ...  27    28    29
   ! ilev value  -7.00, -6.50, -6.00, -5.50, ... 6.00, 6.50, 7.00 ;
   !  lev value     -6.75, -6.25, -5.75, -5.25, ... 6.25, 6.75
   !  lev index        1      2      3      4    ...  27    28

   !mid_level 1
   zgrid_lower(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,1, &
                          domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8)     +  &
                    (get_state(get_dart_vector_index(lon_index,lat_index,2, &
                          domain_id(SECONDARY_DOM), ivarZG), state_handle) /100.0_r8)  ) / 2.0_r8

   !mid_level nlev
   zgrid_upper(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,nilev-1, &
                             domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8)     +  &
                     (get_state(get_dart_vector_index(lon_index,lat_index,nilev, &
                          domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8) ) / 2.0_r8

   ! cannot extrapolate below bottom or beyond top
   do i = 1, n
      if ((zgrid_lower(i) > height) .or. (zgrid_upper(i) < height)) then
        istatus(i) = 55
      endif
   enddo
   if (any(istatus == 55)) return ! ! fail if any ensemble member fails

   ! Figure out what level is above/below, and by how much
   h_loop_midpoint: do k = 2, nilev-1

    zgrid_upper(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,k, &
                          domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8 )     +  &
                       (get_state(get_dart_vector_index(lon_index,lat_index,k+1, &
                          domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8) ) / 2.0_r8

      ! per ensemble member
      do i = 1, n
         if (found(i)) cycle
         if (height <= zgrid_upper(i)) then
            found(i) = .true.
            lev_top(i)    = k
            lev_bottom(i) = lev_top(i) - 1
            if (all(found)) exit h_loop_midpoint
         endif
      enddo

   enddo h_loop_midpoint

   zgrid_lower(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,k-1, &
                         domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8)     +  &
                      (get_state(get_dart_vector_index(lon_index,lat_index,k, &
                         domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8) ) / 2.0_r8


   zgrid_upper(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,k, &
                         domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8)     +  &
                      (get_state(get_dart_vector_index(lon_index,lat_index,k+1, &
                         domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8) ) / 2.0_r8

   delta_z(:) = zgrid_upper(:) - zgrid_lower(:)

   where (zgrid_upper == zgrid_lower)  ! avoid divide by zero ! HK does this happen?
      frac_lev = 0
   elsewhere
      frac_lev = (zgrid_upper - height)/delta_z
   endwhere

if (is_pressure) then ! get fom plevs (pilevs?) array TODO Lanai is always plves

   val_top(:)    = plevs(lev_top(:))     !pressure at midpoint [Pa]
   val_bottom(:) = plevs(lev_bottom(:))  !pressure at midpoint [Pa]
   val(:)        = exp(frac_lev(:) * log(val_bottom(:)) + (1.0 - frac_lev(:)) * log(val_top(:)))

else ! get from state vector

   do i = 1, n
     indx_top(:) = get_dart_vector_index(lon_index,lat_index,lev_top(i), dom_id, var_id)
     indx_bottom(i) = get_dart_vector_index(lon_index,lat_index,lev_bottom(i), dom_id, var_id)
   enddo

   call get_state_array(val_top, indx_top(:), state_handle)
   call get_state_array(val_bottom, indx_bottom(:), state_handle)

   val(:) = frac_lev(:) * val_bottom(:)  + (1.0 - frac_lev(:)) * val_top(:)

endif

istatus(:) = 0

end subroutine vert_interp_lev

!-------------------------------------------------------------------------------
! Compute neighboring lat rows: TIEGCM [-87.5, 87.5] DART [-90, 90]
! HK note from model_interpolate: What should be done?
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE IS NOT GREAT!
subroutine compute_bracketing_lat_indices(lat, idx_below, idx_above, fraction)

real(r8), intent(in)  :: lat ! latitude
integer,  intent(out) :: idx_below, idx_above ! index in lats()
real(r8), intent(out) :: fraction ! fraction to use for interpolation

if(lat >= bot_lat .and. lat <= top_lat) then ! -87.5 <= lat <= 87.5
   idx_below = int((lat - bot_lat) / delta_lat) + 1
   idx_above = idx_below + 1
   fraction = (lat - lats(idx_below) ) / delta_lat
else if(lat < bot_lat) then ! South of bottom lat
   idx_below = 1
   idx_above = 1
   fraction = 1.0_r8
else                        ! North of top lat
   idx_below = nlat
   idx_above = nlat
   fraction = 1.0_r8
endif

end subroutine compute_bracketing_lat_indices

!-------------------------------------------------------------------------------
function interpolate(n, lon_fract, lat_fract, val11, val12, val21, val22) result(obs_val)

integer,  intent(in) :: n ! number of ensemble members
real(r8), intent(in) :: lon_fract, lat_fract
real(r8), dimension(n), intent(in) :: val11, val12, val21, val22
real(r8), dimension(n) :: obs_val

real(r8) :: a(n, 2)

a(:, 1) = lon_fract * val21(:) + (1.0_r8 - lon_fract) * val11(:)
a(:, 2) = lon_fract * val22(:) + (1.0_r8 - lon_fract) * val12(:)

obs_val(:) = lat_fract * a(:,2) + (1.0_r8 - lat_fract) * a(:,1)

end function interpolate

!-------------------------------------------------------------------------------
function ilev_or_lev(dom_id, var_id) result(dim_name)

integer, intent(in) :: dom_id
integer, intent(in) :: var_id
character(len=NF90_MAX_NAME) :: dim_name

integer :: d
! search for either ilev or lev
dim_name = 'null'
do d = 1, get_num_dims(dom_id, var_id)
   dim_name = get_dim_name(dom_id, var_id, d)
   if (dim_name == 'ilev' .or. dim_name == 'lev') exit
enddo

end function ilev_or_lev
!===============================================================================
! End of model_mod
!===============================================================================
end module model_mod
