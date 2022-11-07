! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! TODOs
!  - Nick Dietrich fix_mmr. When do to this?
!  - model_time
!  - get_state_meta_data 2D variables
!  - test vtec

module model_mod

!-------------------------------------------------------------------------------
!
! Interface for HAO-TIEGCM 2.0
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
                             get_close_obs,                                         &
                             loc_get_close_state => get_close_state,                &
                             set_location, get_location,                            &
                             get_dist, query_location,                              &
                             get_close_type, VERTISUNDEF,                           &
                             VERTISPRESSURE, VERTISHEIGHT, VERTISLEVEL,             &
                             vertical_localization_on, set_vertical

use    utilities_mod, only : open_file, close_file, logfileunit,                    &
                             error_handler, E_ERR, E_MSG, E_WARN, nmlfileunit,      &
                             do_output, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, register_module,             &
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
                                get_domain_size, get_model_variable_indices, &
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
                                 nc_define_real_variable, &
                                 nc_check, nc_open_file_readonly, nc_get_dimension_size, &
                                 nc_close_file, nc_get_variable

use dart_time_io_mod,     only : write_model_time

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


character(len=256) :: tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
character(len=256) :: tiegcm_secondary_file_name = 'tiegcm_s.nc'
integer            :: debug = 0
logical            :: estimate_f10_7 = .false.
character(len=256) :: f10_7_file_name = 'f10_7.nc'
integer            :: assimilation_period_seconds = 3600
real(r8)           :: model_res = 5.0_r8

integer, parameter :: MAX_NUM_VARIABLES = 30
integer, parameter :: MAX_NUM_COLUMNS = 6
character(len=NF90_MAX_NAME) :: variables(MAX_NUM_VARIABLES * MAX_NUM_COLUMNS) = ' '

namelist /model_nml/ tiegcm_restart_file_name, &
                     tiegcm_secondary_file_name, &
                     variables, debug, estimate_f10_7, &
                     f10_7_file_name, &
                     assimilation_period_seconds, model_res
                     

!-------------------------------------------------------------------------------
! define model parameters

! nilev is number of interaface levels
! nlev is number of midpoint levels
integer                               :: nilev, nlev, nlon, nlat
real(r8),dimension(:),    allocatable :: lons, lats, levs, ilevs, plevs, pilevs
! levels + top level boundary condition for nlev.
integer                               :: all_nlev
real(r8),dimension(:),    allocatable :: all_levs
! HK are plevs, pilves per ensemble member?
real(r8)                              :: TIEGCM_reference_pressure
integer                               :: time_step_seconds
integer                               :: time_step_days
type(time_type)                       :: time_step

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! variable name
integer, parameter :: VT_KINDINDX     = 2 ! DART quantity
integer, parameter :: VT_MINVALINDX   = 3 ! minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! maximum value if any
integer, parameter :: VT_ORIGININDX   = 5 ! file of origin
integer, parameter :: VT_STATEINDX    = 6 ! update (state) or not

character(len=obstypelength) :: variable_table(MAX_NUM_VARIABLES, MAX_NUM_COLUMNS)

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


! Obs locations are expected to be given in height [m] or level,
! and so vertical localization coordinate is *always* height.
! Note that gravity adjusted geopotential height (ZG) is read in
! "tiegcm_s.nc". ZG is 'cm', dart is mks
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

if (module_initialized) return ! only need to do this once

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

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

! Read in TIEGCM grid definition from TIEGCM restart file
call read_TIEGCM_definition(tiegcm_restart_file_name)

if ( estimate_f10_7 ) then
   call error_handler(E_MSG, 'f10_7 part of DART state', source)
endif

! error-check, convert namelist input to variable_table, and build the
! state structure
call verify_variables()

call set_calendar_type('Gregorian')

! Convert the last year/day/hour/minute to a dart time.
state_time = read_model_time(tiegcm_restart_file_name)

! Assumes assimilation_period is a multiple of the dynamical timestep
! TIEGCM namelist has variable "STOP"
! which is an array of length 3 corresponding to day-of-year, hour, minute
time_step = set_time(assimilation_period_seconds, 0)

end subroutine static_init_model


!-------------------------------------------------------------------------------

function get_model_size()
! Returns the size of the model as an integer.

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



!-------------------------------------------------------------------------------


subroutine model_interpolate(state_handle, ens_size, location, iqty, obs_val, istatus)
! Given a location, and a model state variable qty,
! interpolates the state variable field to that location.
! obs_val is the interpolated value for each ensemble member
! istatus is the success (0) or failure of the interpolation

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: iqty
real(r8),           intent(out) :: obs_val(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer  :: which_vert
integer  :: lat_below, lat_above, lon_below, lon_above ! these are indices
real(r8) :: lon_fract, lat_fract
real(r8) :: lon, lat, lon_lat_lev(3)
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
level       = int(lon_lat_lev(3))


which_vert = nint(query_location(location))

call compute_bracketing_lat_indices(lat, lat_below, lat_above, lat_fract)
call compute_bracketing_lon_indices(lon, lon_below, lon_above, lon_fract)

! Pressure is not part of the state vector
! pressure is static data on plevs/pilevs
if ( iqty == QTY_PRESSURE) then
   if (which_vert == VERTISLEVEL) then
      ! @todo from Lanai code:
      !   Some variables need plevs, some need pilevs
      !   We only need the height (aka level)
      !   the obs_def_upper_atm_mod.f90:get_expected_O_N2_ratio routines queries
      !   for the pressure at the model levels - EXACTLY - so ...
      !   FIXME ... at present ... the only time model_interpolate
      !   gets called with QTY_PRESSURE is to calculate density, which
      !   requires other variables that only live on the midpoints.
      !   I cannot figure out how to generically decide when to
      !   use plevs vs. pilevs

      ! Check to make sure vertical level is possible.
      if ((level < 1) .or. (level > nlev)) then
         istatus(:) = 33
         return
      else
         obs_val(:) = plevs(level)
         istatus(:) = 0
      endif
   elseif (which_vert == VERTISHEIGHT) then

      ! @todo from Lanai code:
      !   FIXME ... is it possible to try to get a pressure with which_vert == undefined
      !   At present, vert_interp will simply fail because height is a negative number.
      !   @todo HK what are you supposed to do for pressure with VERTISUNDEF? level 1?

      call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_below, height, iqty, val11, istatus)
      if (any(istatus /= 0)) return  ! bail at the first failure
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


if ( iqty == QTY_VERTICAL_TEC ) then ! extrapolate vtec

   call extrapolate_vtec(state_handle, ens_size, lon_below, lat_below, val11)
   call extrapolate_vtec(state_handle, ens_size, lon_below, lat_above, val11)
   call extrapolate_vtec(state_handle, ens_size, lon_above, lat_below, val11)
   call extrapolate_vtec(state_handle, ens_size, lon_above, lat_above, val11)
   obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
   istatus(:) = 0

   return
endif

! check if qty is in the state vector
call find_qty_in_state(iqty, dom_id, var_id)
if (dom_id < 0 ) then
   istatus(:) = 44
   return
endif

if( which_vert == VERTISHEIGHT ) then

   call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_below, height, iqty, val11, istatus)
   if (any(istatus /= 0)) return  ! bail at the first failure
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
! Given an integer index into the state vector, returns the
! associated location and optionally the variable quantity.

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_qty

integer  :: lon_index, lat_index, lev_index
integer  :: local_qty, var_id, dom_id
integer  :: seconds, days ! for f10.7 location
real(r8) :: longitude     ! for f10.7 location
character(len=NF90_MAX_NAME) :: dim_name

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, lev_index, var_id=var_id, dom_id=dom_id, kind_index=local_qty)

if(present(var_qty)) var_qty = local_qty

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
   case ('lev')
      location  = set_location(lons(lon_index), lats(lat_index), levs(lev_index), VERTISLEVEL)
   case default
    call error_handler(E_ERR, 'get_state_meta_data', 'expecting ilev or ilat dimension')
    ! HK @todo 2D variables.
end select

end subroutine get_state_meta_data


!-------------------------------------------------------------------------------


subroutine end_model()
! Does any shutdown and clean-up needed for model.

end subroutine end_model


!-------------------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file.
subroutine nc_write_model_atts( ncid, dom_id )

integer, intent(in)  :: ncid      ! netCDF file identifier
integer, intent(in)  :: dom_id

real(r8), allocatable :: temp_lons(:)
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
call nc_define_dimension(ncid, 'lev',  all_nlev,  routine)
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


call nc_end_define_mode(ncid, routine)

!-------------------------------------------------------------------------------
! Write variables
!-------------------------------------------------------------------------------

! Fill in the coordinate variables

! longitude - TIEGCM uses values +/- 180, DART uses values [0,360]
allocate(temp_lons(nlon))
temp_lons = lons
where (temp_lons >= 180.0_r8) temp_lons = temp_lons - 360.0_r8
call nc_put_variable(ncid, 'lon',  temp_lons,  routine)
call nc_put_variable(ncid, 'lat',  lats,  routine)
call nc_put_variable(ncid, 'lev',  all_levs,   routine)
call nc_put_variable(ncid, 'ilev', ilevs,  routine)
deallocate(temp_lons)

! flush any pending i/o to disk
call nc_synchronize_file(ncid, routine)

end subroutine nc_write_model_atts


!-------------------------------------------------------------------------------
! Vertical localization is done only in height (ZG).
! obs vertical location is given in height (model_interpolate).
! state vertical location is given in height.
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
! HK Note if you have inflation on ZG has been inflated.
! Scroll through all the obs_loc(:) and obs_kind(:) elements

do k = 1,num_close
   q_ind  = close_ind(k)
   if (loc_qtys(q_ind) == QTY_GEOMETRIC_HEIGHT) then
      if (do_output() .and. (debug > 99)) then
         write(     *     ,*)'get_close_state ZG distance is ', &
                     dist(k),' changing to ',10.0_r8 * PI
         write(logfileunit,*)'get_close_state ZG distance is ', &
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

subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus(:)

integer  :: current_vert_type, i
real(r8) :: height(1)
integer  :: local_status(1)

character(len=*), parameter :: routine = 'convert_vertical_obs'

if ( which_vert == VERTISHEIGHT .or. which_vert == VERTISUNDEF) then
  istatus(:) = 0
  return
endif

do i = 1, num
   current_vert_type = nint(query_location(locs(i)))
   if (( current_vert_type == which_vert ) .or. &
       ( current_vert_type == VERTISUNDEF)) then
      istatus(i) = 0
      cycle
   endif

  call model_interpolate(state_handle, 1, locs(i), QTY_GEOMETRIC_HEIGHT, height, local_status )
  
  if (local_status(1) == 0) call set_vertical(locs(i), height(1), VERTISHEIGHT)
  istatus(i) = local_status(1)

enddo

end subroutine convert_vertical_obs

!-------------------------------------------------------------------------------
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
integer :: i
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

function read_model_time(filename)
character(len=*),  intent(in) :: filename
type(time_type) :: read_model_time

integer :: ncid, time_dimlen, dimlen

integer,  parameter         :: nmtime = 3
integer,  dimension(nmtime) :: mtime  ! day, hour, minute
integer                     :: year, doy, utsec
integer, allocatable, dimension(:,:) :: mtimetmp
integer, allocatable, dimension(:)   :: yeartmp

character(len=*), parameter :: routine = 'read_model_time'

ncid = nc_open_file_readonly(filename, routine)

time_dimlen = nc_get_dimension_size(ncid, 'time', routine)
dimlen = nc_get_dimension_size(ncid, 'mtimedim', routine)

if (dimlen /= nmtime) then
   write(string1, *) trim(filename), ' mtimedim = ',dimlen, ' DART expects ', nmtime
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

allocate(mtimetmp(dimlen, time_dimlen), yeartmp(time_dimlen))

call nc_get_variable(ncid, 'mtime', mtimetmp, routine)
call nc_get_variable(ncid, 'year', yeartmp, routine)

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

call nc_close_file(ncid, routine, filename)

end function read_model_time


!===============================================================================
! Routines below here are private to the module
!===============================================================================

subroutine read_TIEGCM_definition(file_name)
! Read TIEGCM grid definition from a tiegcm restart file
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
integer  :: ncid, DimID, TimeDimID
real(r8) :: p0

character(len=*), parameter :: routine = 'read_TIEGCM_definition'

call error_handler(E_MSG,routine,'reading restart ['//trim(file_name)//']')

ncid = nc_open_file_readonly(file_name, routine)

! longitude - TIEGCM uses values +/- 180, DART uses values [0,360]
nlon = nc_get_dimension_size(ncid, 'lon', routine)
allocate(lons(nlon))
call nc_get_variable(ncid, 'lon', lons, routine)
where (lons < 0.0_r8) lons = lons + 360.0_r8

! latitiude
nlat = nc_get_dimension_size(ncid, 'lat', routine)
allocate(lats(nlat))
call nc_get_variable(ncid, 'lat', lats, routine)

! pressure
call nc_get_variable(ncid, 'p0', p0, routine)
TIEGCM_reference_pressure = p0

! level
all_nlev = nc_get_dimension_size(ncid, 'lev', routine)
! top level is not viable. The lower boundary condition is stored in the top level
nlev = all_nlev - 1
allocate(all_levs(all_nlev),levs(nlev), plevs(nlev))
call nc_get_variable(ncid, 'lev', all_levs, routine)

levs=all_levs(1:nlev)
plevs = p0 * exp(-levs) * 100.0_r8 ![Pa] = 100* [millibars] = 100* [hPa]

! ilevel
nilev = nc_get_dimension_size(ncid, 'ilev', routine)
allocate(ilevs(nilev), pilevs(nilev))
call nc_get_variable(ncid, 'ilev', ilevs, routine)

pilevs = p0 * exp(-ilevs) * 100.0_r8 ! [Pa] = 100* [millibars] = 100* [hPa]

if ((nlev+1) .ne. nilev) then
   write(string1,*) 'number of midpoints should be 1 less than number of interfaces.' !HK is the top level for nilev not a boundary condition?
   write(string2,*) 'number of midpoints  is nlev  = ',nlev
   write(string3,*) 'number of interfaces is nilev = ',nilev
   call error_handler(E_MSG,'read_TIEGCM_definition', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif


! Get lon and lat grid specs
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

character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: filename
character(len=NF90_MAX_NAME) :: state_or_aux

nrows = size(variable_table,1) ! these are MAX_NUM_VARIABLES, MAX_NUM_COLUMNS
ncols = size(variable_table,2)

! Convert the (input) 1D array "variables" into a table with six columns.
! The number of rows in the table correspond to the number of variables in the
! DART state vector.
! Column 1 is the netCDF variable name.
! Column 2 is the corresponding DART kind.
! Column 3 is the minimum value ("NA" if there is none) Not Applicable
! Column 4 is the maximum value ("NA" if there is none) Not Applicable
! Column 5 is the file of origin tiegcm 'restart' or 'secondary'
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
   if (variable_table(i,VT_ORIGININDX) == 'CALCULATE') nfields_constructed = nfields_constructed + 1

enddo ROWLOOP

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

if (nfields_secondary == 0) call error_handler(E_ERR, 'ZG is required in &model_nml::variables', source)

call load_up_state_structure_from_file(tiegcm_restart_file_name, nfields_restart, 'RESTART', RESTART_DOM)
call load_up_state_structure_from_file(tiegcm_secondary_file_name, nfields_secondary, 'SECONDARY', SECONDARY_DOM)

if (estimate_f10_7) then
   if (nfields_constructed == 0) then
      call error_handler(E_ERR, 'expecting f10.7 in &model_nml::variables', source)
   endif
   call load_up_state_structure_from_file(f10_7_file_name, nfields_constructed, 'CALCULATE', CONSTRUCT_DOM)
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

! remove top level from all lev variables - this is the boundary condition
call hyperslice_domain(domain_id(domain_num), 'lev', nlev)

deallocate(var_names, kind_list, clamp_vals, update_list)

end subroutine load_up_state_structure_from_file
!-------------------------------------------------------------------------------

subroutine extrapolate_vtec(state_handle, ens_size, lon_index, lat_index, vTEC)
!
! Create the vTEC from constituents in state.
!

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index, lat_index
real(r8),            intent(out) :: vTEC(ens_size)

! n(i)levs x ensmeble size
real(r8), allocatable, dimension(:,:) :: NE, ZG
real(r8), allocatable, dimension(:,:) :: TI, TE
real(r8), allocatable, dimension(:,:) :: NEm_extended, ZG_extended
real(r8), allocatable, dimension(:,:)     :: delta_ZG, NE_middle
real(r8), dimension(ens_size)   :: GRAVITYtop, Tplasma, Hplasma

real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! m^2 * kg / s^2 / K
real(r8), PARAMETER :: omass      = 2.678e-26_r8 ! mass of atomic oxgen kg

real(r8) :: earth_radiusm
integer  :: nlevX, nilevX, j, i, var_id
integer(i8) :: idx

! NE,ZG are extrapolated
!  20 more layers for 2.5 degree resolution
!  10 more layers for 5 degree resolution
if (model_res == 2.5) then
  nlevX  = nlev + 20
  nilevX  = nilev + 20
else
  nlevX = nlev + 10
  nilevX = nilev + 10
endif


allocate( NE(nilev, ens_size), NEm_extended(nilevX, ens_size), &
          ZG(nilev, ens_size), ZG_extended(nilevX, ens_size))
allocate( TI(nlev, ens_size), TE(nlev, ens_size) )
allocate( delta_ZG(nlevX-1, ens_size), NE_middle(nlevX-1, ens_size) )

! NE (interfaces)
var_id = get_varid_from_varname(domain_id(RESTART_DOM), 'NE')
do i = 1, nilev
   idx = get_dart_vector_index(lon_index,lat_index, i, &
                            domain_id(RESTART_DOM), var_id)
   NE(i, :) = get_state(idx, state_handle)
enddo

! ZG (interfaces)
do i = 1, nilev
  idx = get_dart_vector_index(lon_index,lat_index, i, &
                         domain_id(RESTART_DOM), var_id)
  ZG(i, :) = get_state(idx, state_handle)
enddo

! TI (midpoints)
var_id = get_varid_from_varname(domain_id(RESTART_DOM), 'TI')
do i = 1, nlev
   idx = get_dart_vector_index(lon_index,lat_index, i, &
                          domain_id(RESTART_DOM), var_id)
   TI(i, :) = get_state(idx, state_handle)
enddo

! TE (midpoints)
var_id = get_varid_from_varname(domain_id(RESTART_DOM), 'TE')
do i = 1, nlev
   idx = get_dart_vector_index(lon_index,lat_index, i, &
                          domain_id(RESTART_DOM), var_id)
   TE(i, :) = get_state(idx, state_handle)
enddo

! Construct vTEC given the parts

earth_radiusm = earth_radius * 1000.0_r8 ! Convert earth_radius in km to m
NE            = NE * 1.0e+6_r8           ! Convert NE in #/cm^3 to #/m^3

! Gravity at the top layer
GRAVITYtop(:) = gravity * (earth_radiusm / (earth_radiusm + ZG(nilev,:))) ** 2

! Plasma Temperature
Tplasma(:) = (TI(nlev-1,:) + TE(nlev-1,:)) / 2.0_r8

! Compute plasma scale height
Hplasma(:) = (2.0_r8 * k_constant / omass ) * Tplasma(:) / GRAVITYtop(:)

ZG_extended(1:nilev,:) = ZG
NEm_extended(1:nilev,:) = NE

do j = nlev, nlevX
   NEm_extended(j,:) = NEm_extended(j-1,:) * exp(-0.5_r8)
    ZG_extended(j,:) =  ZG_extended(j-1,:) + Hplasma(:) / 2.0_r8
enddo

delta_ZG(1:(nlevX-1),:) =  ZG_extended(2:nlevX,:) -  ZG_extended(1:(nlevX-1),:)
NE_middle(1:(nlevX-1),:) = (NEm_extended(2:nlevX,:) + NEm_extended(1:(nlevX-1),:)) / 2.0_r8

do i = 1, ens_size
   vTEC(i) = sum(NE_middle(:,i) * delta_ZG(:,i)) * 1.0e-16_r8 ! Convert to TECU (1.0e+16 #/m^2)
enddo

deallocate( NE, NEm_extended, ZG, ZG_extended)
deallocate( TI, TE )
deallocate( delta_ZG, NE_middle )

end subroutine extrapolate_vtec


!-------------------------------------------------------------------------------

subroutine vert_interp(state_handle, n, dom_id, var_id, lon_index, lat_index, height, iqty, &
                       val, istatus)
! returns the value at an arbitrary height on an existing horizontal grid location.
! istatus == 0 is success.

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
   vertstagger = 'ilev'
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
! Returns the variable id for a given DART qty
! Will return X rather than X_MN variable.

integer, intent(in)  :: iqty
integer, intent(out) :: which_dom
integer, intent(out) :: var_id

integer :: num_same_kind, id, k
integer, allocatable :: multiple_kinds(:), n
character(NF90_MAX_NAME) :: varname

which_dom = -1
var_id = -1

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
subroutine compute_bracketing_lon_indices(lon, idx_below, idx_above, fraction)

real(r8), intent(in)  :: lon ! longitude
integer,  intent(out) :: idx_below, idx_above ! index in lons()
real(r8), intent(out) :: fraction ! fraction to use for interpolation

if(lon >= top_lon .and. lon < bot_lon) then     ! at wraparound point [175 <= lon < 180]
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

   where (zgrid == z2)  ! avoid divide by zero
      frac_lev = 0.0_r8
      delta_z = 0.0_r8
   elsewhere
      delta_z = (zgrid - z2)/100.0_r8
      frac_lev = (zgrid/100.0_r8 - height)/delta_z
   endwhere

   if (is_pressure) then ! get fom plevs (pilevs?) array @todo HK Lanai is always plves

      val_top(:)    = plevs(lev_top(:))     !pressure at midpoint [Pa]
      val_bottom(:) = plevs(lev_bottom(:))  !pressure at midpoint [Pa]
      val(:)        = exp(frac_lev(:) * log(val_bottom(:)) + (1.0 - frac_lev(:)) * log(val_top(:)))
  
   else  ! get from state vector

      do i = 1, n
        indx_top(i) = get_dart_vector_index(lon_index,lat_index,lev_top(i), dom_id, var_id)
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

integer :: lev(n), lev_minus_one(n), lev_plus_one(n)
real(r8) :: frac_lev(n)

integer  :: k, i
real(r8) :: delta_z(n)
real(r8) :: zgrid_upper(n), zgrid_lower(n) ! ZG on midpoints
real(r8) :: z_k(n), z_k_minus_one(n), z_k_plus_one(n)  ! ZG on ilves
integer(i8)  :: indx_top(n), indx_bottom(n) ! state vector indices for qty
integer(i8)  :: indx(n), indx_minus_one(n), indx_plus_one(n) ! state vector indices for ZG
logical  :: found(n) ! track which ensemble members have been located
real(r8) :: val_top(n), val_bottom(n)

istatus    = 1
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
            lev(i) = k
            lev_minus_one(i) = lev(i) - 1
            lev_plus_one(i) = lev(i) + 1
            if (all(found)) exit h_loop_midpoint
         endif
      enddo

   enddo h_loop_midpoint

   do i = 1, n
     indx(i) = get_dart_vector_index(lon_index,lat_index,lev(i), domain_id(SECONDARY_DOM), ivarZG)
     indx_minus_one(i) = get_dart_vector_index(lon_index,lat_index,lev_minus_one(i), domain_id(SECONDARY_DOM), ivarZG)
     indx_plus_one(i) = get_dart_vector_index(lon_index,lat_index,lev_plus_one(i), domain_id(SECONDARY_DOM), ivarZG)
   enddo

   call get_state_array(z_k(:),indx(:), state_handle)
   call get_state_array(z_k_minus_one, indx_minus_one(:), state_handle)  
   call get_state_array(z_k_plus_one, indx_plus_one(:), state_handle)


   !lower midpoint   
   zgrid_lower(:) = ( z_k(:) + z_k_minus_one ) / 2.0_r8 / 100.0_r8
   
   ! upper midpoint
   zgrid_upper(:) = ( z_k(:) + z_k_plus_one ) / 2.0_r8 / 100.0_r8

   where (zgrid_upper == zgrid_lower)  ! avoid divide by zero
      frac_lev = 0.0_r8
      delta_z = 0.0_r8
   elsewhere
      delta_z = zgrid_upper - zgrid_lower
      frac_lev = (zgrid_upper - height)/delta_z
   endwhere

if (is_pressure) then ! get fom plevs 

   val_top(:)    = plevs(lev(:))     !pressure at midpoint [Pa]
   val_bottom(:) = plevs(lev_minus_one(:))  !pressure at midpoint [Pa]
   val(:)        = exp(frac_lev(:) * log(val_bottom(:)) + (1.0 - frac_lev(:)) * log(val_top(:)))

else ! get from state vector

   do i = 1, n
     indx_top(i) = get_dart_vector_index(lon_index,lat_index,lev(i), dom_id, var_id)
     indx_bottom(i) = get_dart_vector_index(lon_index,lat_index,lev_minus_one(i), dom_id, var_id)
   enddo

   call get_state_array(val_top, indx_top(:), state_handle)
   call get_state_array(val_bottom, indx_bottom(:), state_handle)

   val(:) = frac_lev(:) * val_bottom(:)  + (1.0 - frac_lev(:)) * val_top(:)

endif

istatus(:) = 0

end subroutine vert_interp_lev

!-------------------------------------------------------------------------------
! Compute neighboring lat rows: TIEGCM [-87.5, 87.5] DART [-90, 90]
! Poles >|87.5| set to |87.5| 
subroutine compute_bracketing_lat_indices(lat, idx_below, idx_above, fraction)

real(r8), intent(in)  :: lat ! latitude
integer,  intent(out) :: idx_below, idx_above ! index in lats()
real(r8), intent(out) :: fraction ! fraction to use for interpolation

if(lat >= bot_lat .and. lat < top_lat) then ! -87.5 <= lat < 87.5
   idx_below = int((lat - bot_lat) / delta_lat) + 1
   idx_above = idx_below + 1
   fraction = (lat - lats(idx_below) ) / delta_lat
else if(lat < bot_lat) then ! South of bottom lat
   idx_below = 1
   idx_above = 1
   fraction = 1.0_r8
else                        ! On or North of top lat
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
