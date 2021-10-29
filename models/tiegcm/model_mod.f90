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
                             earth_radius, gravity, obstypelength

use time_manager_mod, only : time_type, set_calendar_type, set_time_missing,        &
                             set_time, get_time, print_time,                        &
                             set_date, get_date, print_date,                        &
                             operator(*),  operator(+), operator(-),                &
                             operator(>),  operator(<), operator(/),                &
                             operator(/=), operator(<=)

use     location_mod, only : location_type,                                         &
                             loc_get_close_obs => get_close_obs,                    &
                             set_location, get_location,                            &
                             get_dist, query_location,                              &
                             get_close_type, VERTISUNDEF,                           &
                             VERTISPRESSURE, VERTISHEIGHT, VERTISLEVEL

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

use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use mpi_utilities_mod,only : my_task_id

use default_model_mod, only : adv_1step,                                &
                              init_conditions => fail_init_conditions,  &
                              init_time => fail_init_time,              &
                              nc_write_model_vars,                      &
                              pert_model_copies

use state_structure_mod, only : add_domain, get_dart_vector_index, add_dimension_to_variable, &
                                finished_adding_domain, state_structure_info, &
                                get_domain_size, set_parameter_value, get_model_variable_indices, &
                                get_num_dims, get_dim_name

use distributed_state_mod, only : get_state

use ensemble_manager_mod, only : ensemble_type

use netcdf_utilities_mod, only : nc_synchronize_file

use typesizes  !HK do we needs these with netcdf_utilities_mod?
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

!TIEGCM specific routines
public :: get_f107_value

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = 'tiegcm/model_mod.f90'
character(len=32 ), parameter :: revision = ''
character(len=128), parameter :: revdate  = ''

!-------------------------------------------------------------------------------
! namelist with default values
! output_state_vector = .true.  results in a "state-vector"   netCDF file
! output_state_vector = .false. results in a "prognostic-var" netCDF file

! IMPORTANT: Change output file names in tiegcm.nml to match these names
! i.e.  OUTPUT='tiegcm_restart_p.nc'
!       SECOUT='tiegcm_s.nc'
character(len=256) :: tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
character(len=256) :: tiegcm_secondary_file_name = 'tiegcm_s.nc'
character(len=256) :: tiegcm_namelist_file_name  = 'tiegcm.nml'
logical            :: output_state_vector = .false.
integer            :: debug = 0
logical            :: estimate_f10_7 = .false.
integer            :: assimilation_period_seconds = 3600

integer, parameter :: MAX_NUM_VARIABLES = 30
integer, parameter :: MAX_NUM_COLUMNS = 6
character(len=NF90_MAX_NAME) :: variables(MAX_NUM_VARIABLES * MAX_NUM_COLUMNS) = ' '

namelist /model_nml/ tiegcm_restart_file_name, &
                     tiegcm_secondary_file_name, tiegcm_namelist_file_name, &
                     variables, debug, estimate_f10_7, assimilation_period_seconds

!-------------------------------------------------------------------------------
! define model parameters

integer                               :: nilev, nlev, nlon, nlat
real(r8),dimension(:),    allocatable :: lons, lats, levs, ilevs, plevs, pilevs
! HK are plevs, pilves per ensemble member?
real(r8)                              :: TIEGCM_missing_value !! global attribute
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

! FOR NOW OBS LOCATIONS ARE EXPECTED GIVEN IN HEIGHT [m],
! AND SO VERTICAL LOCALIZATION COORDINATE IS *always* HEIGHT
! (note that gravity adjusted geopotential height (ZG)
!  read in from "tiegcm_s.nc" *WARNING* ZG is 'cm', DART is mks)

real(r8), allocatable, dimension(:,:,:) :: ZG
integer               :: ivarZG

logical               :: first_pert_call = .true.
type(random_seq_type) :: random_seq
character(len=512)    :: string1, string2, string3
logical, save         :: module_initialized = .false.

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum


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
!HK call read_TIEGCM_secondary(tiegcm_secondary_file_name)

! error-check, convert namelist input to variable_table, and build the
! state structure
call verify_variables()

!HK There are two files and a calculated variable VTEC:
! tiegcm_restart_file_name
! tiegcm_secondary_file_name
! calcualted variable VTEC

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
integer,            intent(out) :: istatus(ens_size)
real(r8),           intent(out) :: obs_val(ens_size) !< array of interpolated values

integer  :: ivar
integer  :: i, which_vert
integer  :: lat_below, lat_above, lon_below, lon_above
integer  :: zero_lon_index
real(r8) :: lon_fract, lat_fract
real(r8) :: lon, lat, lon_lat_lev(3)
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real(r8) :: val(2,2), a(2)
real(r8) :: height
integer  :: level

if ( .not. module_initialized ) call static_init_model

! Default for failure return
istatus = 1
obs_val = MISSING_R8

! GITM uses a vtec routine in obs_def_upper_atm_mod:get_expected_gnd_gps_vtec()
! TIEGCM has its own vtec routine, so we should use it. This next block ensures that.
! The get_expected_gnd_gps_vtec() tries to interpolate QTY_GEOPOTENTIAL_HEIGHT
! when it does, this will kill it. 

if ( iqty == QTY_GEOPOTENTIAL_HEIGHT ) then
   write(string1,*)'QTY_GEOPOTENTIAL_HEIGHT currently unsupported'
   call error_handler(E_ERR,'model_interpolate',string1,source, revision, revdate)
endif

! Get the position

lon_lat_lev = get_location(location)
lon         = lon_lat_lev(1) ! degree
lat         = lon_lat_lev(2) ! degree

which_vert = nint(query_location(location))

!HK overdoing these if statments - I don't think the vertical location
! changes thoughout this routine, check once and be done.
if( which_vert == VERTISHEIGHT ) then
   height = lon_lat_lev(3)
elseif( which_vert == VERTISLEVEL ) then
   level  = int(lon_lat_lev(3))
elseif( which_vert == VERTISUNDEF) then
   height = MISSING_R8
   level  = -1
else
   which_vert = nint(query_location(location))
   write(string1,*) 'vertical coordinate type:',which_vert,' cannot be handled'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
endif

! Check to make sure vertical level is possible.
if (which_vert == VERTISLEVEL) then
   if ((level < 1) .or. (level > nilev)) return
endif

if ((iqty == QTY_PRESSURE) .and. (which_vert == VERTISLEVEL)) then
   ! Some variables need plevs, some need pilevs
   ! We only need the height (aka level)
   ! the obs_def_upper_atm_mod.f90:get_expected_O_N2_ratio routines queries
   ! for the pressure at the model levels - EXACTLY - so ...
   ! FIXME ... at present ... the only time model_interpolate
   ! gets called with QTY_PRESSURE is to calculate density, which
   ! requires other variables that only live on the midpoints.
   ! I cannot figure out how to generically decide when to
   ! use plevs vs. pilevs 
   if (level <= nlev) then
      obs_val(:) = plevs(level)
      istatus(:) = 0
   endif
   return
endif

! Get lon and lat grid specs
bot_lon        = lons(1)                         ! 180.
delta_lon      = abs((lons(1)-lons(2)))          ! 5. or 2.5
zero_lon_index = int(bot_lon/delta_lon) + 1      ! 37 or 73
top_lon        = lons(nlon)                      ! 175. or 177.5
bot_lat        = lats(1)                         !
top_lat        = lats(nlat)                      !
delta_lat      = abs((lats(1)-lats(2)))          !

! Compute bracketing lon indices:
! TIEGCM [-180 175]  DART [180, 185, ..., 355, 0, 5, ..., 175]
if(lon > top_lon .and. lon < bot_lon) then     ! at wraparound point [175 < lon < 180]
   lon_below = nlon
   lon_above = 1
   lon_fract = (lon - top_lon) / delta_lon
elseif (lon >= bot_lon) then                  ! [180 <= lon <= 360]
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - lons(lon_below)) / delta_lon
else                                           ! [0 <= lon <= 175 ]
   lon_below = int((lon - 0.0_r8) / delta_lon) + zero_lon_index
   lon_above = lon_below + 1
   lon_fract = (lon - lons(lon_below)) / delta_lon
endif

! Compute neighboring lat rows: TIEGCM [-87.5, 87.5] DART [-90, 90]
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE IS NOT GREAT!
if(lat >= bot_lat .and. lat <= top_lat) then ! -87.5 <= lat <= 87.5
   lat_below = int((lat - bot_lat) / delta_lat) + 1
   lat_above = lat_below + 1
   lat_fract = (lat - lats(lat_below) ) / delta_lat
else if(lat < bot_lat) then ! South of bottom lat
   lat_below = 1
   lat_above = 1
   lat_fract = 1.0_r8
else                        ! North of top lat
   lat_below = nlat
   lat_above = nlat
   lat_fract = 1.0_r8
endif

! If PRESSURE; calculate the pressure from several variables.
! vert_interp() interpolates the state column to
! the same vertical height as the observation.
! THEN, it is the same as the 2D case.

! FIXME ... is it possible to try to get a pressure with which_vert == undefined
! At present, vert_interp will simply fail because height is a negative number.
if (iqty == QTY_PRESSURE) then

!HK  enesmble size vs column
   !call vert_interp(obs_val(:),lon_below,lat_below,height,iqty,'ilev',-1,val(1,1),istatus(1))
   !if (istatus(1) == 0) &
   !call vert_interp(obs_val(:),lon_below,lat_above,height,iqty,'ilev',-1,val(1,2),istatus(1))
   !if (istatus(1) == 0) &
   !call vert_interp(obs_val(:),lon_above,lat_below,height,iqty,'ilev',-1,val(2,1),istatus(1))
   !if (istatus(1) == 0) &
   !call vert_interp(obs_val(:),lon_above,lat_above,height,iqty,'ilev',-1,val(2,2),istatus(1))

   ! Now that we have the four surrounding points and their relative weights,
   ! actually perform the bilinear horizontal interpolation to get the value at the
   ! (arbitrary) desired location.

   if (istatus(1) == 0) then
      do i = 1, 2
         a(i) = lon_fract * val(2, i) + (1.0_r8 - lon_fract) * val(1, i)
      end do
      obs_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)
   endif

   return
endif

! If it is not pressure ...
! FindVar_by_kind would fail with iqty == QTY_PRESSURE, so we have
! to calculate the pressure separately before this part.

!HK
!ivar = FindVar_by_kind(iqty)
!if (ivar < 0) return ! as a failure

! Now, need to find the values for the four corners

if ( which_vert == VERTISUNDEF ) then ! (time, lat, lon)
   ! time is always a singleton dimension

   !val(1,1) = x(get_dart_vector_index(ivar, indx1=lon_below, indx2=lat_below, indx3=1))
   !val(1,2) = x(get_dart_vector_index(ivar, indx1=lon_below, indx2=lat_above, indx3=1))
   !val(2,1) = x(get_dart_vector_index(ivar, indx1=lon_above, indx2=lat_below, indx3=1))
   !val(2,2) = x(get_dart_vector_index(ivar, indx1=lon_above, indx2=lat_above, indx3=1))
   !istatus  = 0

elseif (which_vert == VERTISLEVEL) then

   ! one use of model_interpolate is to allow other modules/routines
   ! the ability to 'count' the model levels. To do this, create observations
   ! with locations on model levels and 'interpolate' for QTY_GEOMETRIC_HEIGHT.
   ! When the interpolation fails, you've gone one level too far. 
   ! The only geometric height variable we have is ZG, and its a 4D variable.

   !val(1,1) = x(get_dart_vector_index(ivar, indx1=lon_below, indx2=lat_below, indx3=level))
   !val(1,2) = x(get_dart_vector_index(ivar, indx1=lon_below, indx2=lat_above, indx3=level))
   !val(2,1) = x(get_dart_vector_index(ivar, indx1=lon_above, indx2=lat_below, indx3=level))
   !val(2,2) = x(get_dart_vector_index(ivar, indx1=lon_above, indx2=lat_above, indx3=level))
   !istatus = 0

elseif (which_vert == VERTISHEIGHT) then

   ! vert_interp() interpolates the state column to
   ! the same vertical height as the observation.
   ! THEN, it is the same as the 2D case.
!HK Yo ensemble size vs column
   !call vert_interp(obs_val(:), lon_below, lat_below, height, iqty, &
   !          progvar(ivar)%verticalvar, ivar, val(1,1), istatus(1))
   !if (istatus(1) == 0) &
   !call vert_interp(obs_val(:), lon_below, lat_above, height, iqty, &
   !          progvar(ivar)%verticalvar, ivar, val(1,2), istatus(1))
   !if (istatus(1) == 0) &
   !call vert_interp(obs_val(:), lon_above, lat_below, height, iqty, &
   !          progvar(ivar)%verticalvar, ivar, val(2,1), istatus(1))
   !if (istatus(1) == 0) &
   !call vert_interp(obs_val(:), lon_above, lat_above, height, iqty, &
   !          progvar(ivar)%verticalvar, ivar, val(2,2), istatus(1))

else

   write(string1,*)'has failed'
   write(string2,*)'for some reason'
   call error_handler(E_ERR,'model_interpolate', string1, &
              source, revision, revdate, text2=string2)

endif

! Now that we have the four surrounding points and their relative weights,
! actually perform the bilinear interpolation to get the value at the
! (arbitrary) desired location.

if (istatus(1) == 0) then
   do i = 1, 2
      a(i) = lon_fract * val(2, i) + (1.0_r8 - lon_fract) * val(1, i)
   end do
   obs_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)
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
integer  :: ivar, seconds, days
real(r8) :: height, longitude

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, lev_index, var_id=var_id, dom_id=dom_id, kind_index=local_qty)

!HK check for f10.7 by varname?

!location  = set_location(lons(lon_index), lats(lat_index), height, VERTISHEIGHT)

if(present(var_qty)) var_qty = local_qty

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

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_tiegcm_namelist


if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! Determine shape of namelist.
! long lines are truncated when read into textblock
!-------------------------------------------------------------------------------

call find_textfile_dims(tiegcm_namelist_file_name, nlines, linelen)
if (nlines > 0) then
  has_tiegcm_namelist = .true.
else
  has_tiegcm_namelist = .false.
endif

! HK the netcdf calls need updating to use netcdf utilities
!if (has_tiegcm_namelist) then
!   allocate(textblock(nlines))
!   textblock = ''
!
!   call nc_check(nf90_def_dim(ncid=ncid, name='tiegcmNMLnlines', &
!          len = nlines, dimid = nlinesDimID), &
!          'nc_write_model_atts', 'def_dim tiegcmNMLnlines')
!
!   call nc_check(nf90_def_var(ncid,name='tiegcm_nml', xtype=nf90_char, &
!          dimids = (/ linelenDimID, nlinesDimID /), varid=nmlVarID), &
!          'nc_write_model_atts', 'def_var tiegcm_namelist')
!
!   call nc_check(nf90_put_att(ncid, nmlVarID, 'long_name', &
!          'contents of '//trim(tiegcm_namelist_file_name)), &
!          'nc_write_model_atts', 'put_att tiegcm_namelist')
!
!endif
!
!
!!-------------------------------------------------------------------------------
!! Fill the variables we can
!!-------------------------------------------------------------------------------
!
!if (has_tiegcm_namelist) then
!   call file_to_text(tiegcm_namelist_file_name, textblock)
!   call nc_check(nf90_put_var(ncid, nmlVarID, textblock ), &
!                 'nc_write_model_atts', 'put_var nmlVarID')
!   deallocate(textblock)
!endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!-------------------------------------------------------------------------------
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:)
integer(i8),                   intent(in)     :: loc_indx(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

call error_handler(E_ERR, 'watch out', 'not done yetc')

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
! STATE VERTICAL LOCATION IS GIVEN IN HEIGHT (get_state_meta_data)

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

integer                              :: k, t_ind

!HK! Finds all the locations or observations that are close.
!HKcall loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
!HK                                           num_close, close_ind)
!HK
!HK! Make the ZG part of the state vector far from everything so it does not get updated.
!HK! Scroll through all the obs_loc(:) and obs_kind(:) elements
!HK
!HKdo k = 1,num_close
!HK   t_ind  = close_ind(k)
!HK   if (obs_kind(t_ind) == QTY_GEOMETRIC_HEIGHT) then
!HK      if (do_output() .and. (debug > 99)) then
!HK         write(     *     ,*)'get_close_obs ZG distance is ', &
!HK                     dist(k),' changing to ',10.0_r8 * PI
!HK         write(logfileunit,*)'get_close_obs ZG distance is ', &
!HK                     dist(k),' changing to ',10.0_r8 * PI
!HK      endif
!HK      dist(k) = 10.0_r8 * PI
!HK   endif
!HKenddo
!HK
!HK
!HKif (estimate_f10_7) then
!HK   do k = 1, num_close
!HK
!HK      t_ind  = close_ind(k)
!HK
!HK      ! This was part of an experimental setup - if the distance is LARGE,
!HK      ! the state variables will not be updated. By increasing the distance,
!HK      ! the values in the DART vector will remain unchanged. If you allow
!HK      ! the DART vector to be updated, the posterior observation operators
!HK      ! will be impacted - which is usually the desire. You can then avoid
!HK      ! impacting the tiegcm forecast through the input.nml 'NO_COPY_BACK' feature.
!HK
!HK      if (    (obs_kind(t_ind) == QTY_MOLEC_OXYGEN_MIXING_RATIO) &
!HK         .or. (obs_kind(t_ind) == QTY_U_WIND_COMPONENT) &
!HK         .or. (obs_kind(t_ind) == QTY_V_WIND_COMPONENT) &
!HK         .or. (obs_kind(t_ind) == QTY_TEMPERATURE) ) then
!HK      !  dist(k) = 10.0_r8 * PI
!HK
!HK      elseif  (obs_kind(t_ind) == QTY_1D_PARAMETER) then
!HK         ! f10_7 is given a location of latitude 0.0 and the longitude
!HK         ! of local noon. By decreasing the distance from the observation
!HK         ! to the dynamic f10_7 location we are allowing the already close
!HK         ! observations to have a larger impact in the parameter estimation.
!HK         ! 0.25 is heuristic. The 'close' observations have already been 
!HK         ! determined by the cutoff. Changing the distance here does not
!HK         ! allow more observations to impact anything.
!HK         dist(k) = dist(k)*0.25_r8
!HK      endif
!HK
!HK   enddo
!HKendif

end subroutine get_close_obs


!-------------------------------------------------------------------------------

!===============================================================================
! TIEGCM public routines
!===============================================================================


function get_f107_value(x)
real(r8), dimension(:), intent(in) :: x
real(r8)                           :: get_f107_value

integer :: ivar

if ( .not. module_initialized ) call static_init_model

! If the f10_7 is part of the DART state, return that value.
! If it is not part of the DART state, just return the value from
! the TIEGCM namelist.

if (estimate_f10_7) then

 !HK get_state of f10_7

   call error_handler(E_ERR,'get_f107_value', 'no f10_7 in DART state', &
        source, revision, revdate)

else
   get_f107_value = f10_7
endif

end function get_f107_value



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
nlev = nlev - 1 ! top level is not viable
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
   call error_handler(E_ERR,'read_TIEGCM_definition', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

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
integer  :: ivar, index1, indexN, varsize
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

   string1 = 'Estimating f10_7 is not supported.'
   string2 = 'This feature is under development.'
   string3 = 'If you want to experiment with this, change E_ERR to E_MSG in "get_variables_in_domain".'
   call error_handler(E_MSG, 'get_variables_in_domain:', string1, &
                      source, revision, revdate, text2=string2, text3=string3)

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
call load_up_calculated_variables(nfields_constructed, 'CALCULATE', CONSTRUCT_DOM)

model_size = get_domain_size(RESTART_DOM) + get_domain_size(SECONDARY_DOM) &
                          + get_domain_size(CONSTRUCT_DOM)

print*, 'restart size', get_domain_size(RESTART_DOM)
print*, 'secondary size', get_domain_size(SECONDARY_DOM)
print*, 'constructed size', get_domain_size(CONSTRUCT_DOM)

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

domain_id(domain_num) = add_domain(nvar, var_names, kind_list, clamp_vals, update_list, .true.)

do i = 1, nvar ! HK f10_7 and VTEC and anything else?
  if (var_names(i) == 'f10_7') then
    ! dimensions to match what? Lanai has single value f10_7
    call add_dimension_to_variable(domain_id(domain_num), i, 'parameter', 1)
    call set_parameter_value(domain_id(domain_num), i, f10_7)
  endif

  !HK I don't understand why VTEC is a state variable vs. calculating it on the fly.
  if(var_names(i) == 'VTEC') then
     call error_handler(E_ERR,'load_up_calculated_variables', 'VTEC not sure yet.')
  endif
enddo

call finished_adding_domain(domain_id(domain_num))

end subroutine load_up_calculated_variables

!-------------------------------------------------------------------------------


subroutine create_vtec( ncid, last_time, vTEC)
!
! Create the vTEC from constituents in the netCDF file.
!

integer,                  intent(in)  :: ncid
integer,                  intent(in)  :: last_time
real(r8), dimension(:,:), intent(out) :: vTEC

real(r8), allocatable, dimension(:,:,:) :: NE, TI, TE
real(r8), allocatable, dimension(:,:,:) :: NEm_extended, ZG_extended
real(r8), allocatable, dimension(:,:)   :: GRAVITYtop, Tplasma, Hplasma
real(r8), allocatable, dimension(:)     :: delta_ZG, NE_middle

real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! m^2 * kg / s^2 / K
real(r8), PARAMETER :: omass      = 2.678e-26_r8 ! mass of atomic oxgen kg

real(r8) :: earth_radiusm
integer  :: VarID, nlev10, j, k

allocate( NE(nlon,nlat,nilev), NEm_extended(nlon,nlat,nilev+10), &
          ZG_extended(nlon,nlat,nilev+10))
allocate( TI(nlon,nlat,nlev), TE(nlon,nlat,nlev) )
allocate( GRAVITYtop(nlon,nlat), Tplasma(nlon,nlat), Hplasma(nlon,nlat) )
allocate( delta_ZG(nlev+9), NE_middle(nlev+9) )

!... NE (interfaces)
call nc_check(nf90_inq_varid(ncid, 'NE', VarID), 'create_vtec', 'inq_varid NE')
call nc_check(nf90_get_var(ncid, VarID, values=NE,     &
                   start = (/    1,    1,     1, last_time /),   &
                   count = (/ nlon, nlat, nilev,         1 /)),&
                   'create_vtec', 'get_var NE')

!... ZG (interfaces) already read into the module and converted to metres

!... TI (midpoints)
call nc_check(nf90_inq_varid(ncid, 'TI', VarID), 'create_vtec', 'inq_varid TI')
call nc_check(nf90_get_var(ncid, VarID, values=TI,     &
                   start = (/    1,    1,    1, last_time /),   &
                   count = (/ nlon, nlat, nlev,         1 /)), &
                   'create_vtec', 'get_var TI')

!... TE (midpoints)
call nc_check(nf90_inq_varid(ncid, 'TE', VarID), 'create_vtec', 'inq_varid TE')
call nc_check(nf90_get_var(ncid, VarID, values=TE,     &
                   start = (/    1,    1,    1, last_time /),   &
                   count = (/ nlon, nlat, nlev,         1 /)), &
                   'create_vtec', 'get_var TE')

! Construct vTEC given the parts

earth_radiusm = earth_radius * 1000.0_r8 ! Convert earth_radius in km to m
NE            = NE * 1.0e+6_r8           ! Convert NE in #/cm^3 to #/m^3

! Gravity at the top layer
GRAVITYtop(:,:) = gravity * (earth_radiusm / (earth_radiusm + ZG(:,:,nilev))) ** 2

! Plasma Temperature
Tplasma(:,:) = (TI(:,:,nlev-1) + TE(:,:,nlev-1)) / 2.0_r8

! Compute plasma scale height
Hplasma = (2.0_r8 * k_constant / omass ) * Tplasma / GRAVITYtop

! NE is extrapolated to 10 more layers
nlev10  = nlev + 10

 ZG_extended(:,:,1:nilev) = ZG
NEm_extended(:,:,1:nilev) = NE

do j = nlev, nlev10
   NEm_extended(:,:,j) = NEm_extended(:,:,j-1) * exp(-0.5_r8)
    ZG_extended(:,:,j) =  ZG_extended(:,:,j-1) + Hplasma(:,:) / 2.0_r8
enddo

! finally calculate vTEC - one gridcell at a time.

do j = 1, nlat
do k = 1, nlon
    delta_ZG(1:(nlev10-1)) =  ZG_extended(k,j,2:nlev10) -  ZG_extended(k,j,1:(nlev10-1))
   NE_middle(1:(nlev10-1)) = (NEm_extended(k,j,2:nlev10) + NEm_extended(k,j,1:(nlev10-1))) / 2.0_r8
   vTEC(k,j) = sum(NE_middle * delta_ZG) * 1.0e-16_r8 ! Convert to TECU (1.0e+16 #/m^2)
enddo
enddo

deallocate( NE, NEm_extended, ZG_extended)
deallocate( TI, TE )
deallocate( GRAVITYtop, Tplasma, Hplasma )
deallocate( delta_ZG, NE_middle )

end subroutine create_vtec

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

istatus = 0

do i = 1, num

   call get_model_variable_indices(loc_indx(i), lon_index, lat_index, lev_index, var_id=var_id, dom_id=dom_id)
   
   ! search for either ilev or lev
   dim_name = 'null'
   do d = 1, get_num_dims(dom_id, var_id)
      dim_name = get_dim_name(dom_id, var_id, d)
      if (dim_name == 'ilev' .or. dim_name == 'lev') exit
   enddo
   
   select case (trim(dim_name))
      case ('ilev')
         height = get_state(loc_indx(i), state_handle)
   
      case ('lev') ! height on midpoint
        height1 = get_state(loc_indx(i), state_handle)
        height2 = get_state(get_dart_vector_index(lon_index, lat_index, loc_indx(i)+1, dom_id, var_id), state_handle)
        height = (height1 + height2) / 2.0_r8
   
      case default
       call error_handler(E_ERR, 'convert_vertical_state', 'expecting ilev or ilat dimension')
   
   end select
   
   locs(i) = set_location(lons(lon_index), lats(lat_index), height(1), VERTISHEIGHT)

end do


end subroutine convert_vertical_state

!-------------------------------------------------------------------------------

subroutine vert_interp(x, lon_index, lat_index, height, iqty, vertstagger, &
                       ivar, val, istatus)
! returns the value at an arbitrary height on an existing horizontal grid location.
! istatus == 0 is a 'good' return code.

real(r8),         intent(in)  :: x(:)
integer,          intent(in)  :: lon_index
integer,          intent(in)  :: lat_index
real(r8),         intent(in)  :: height
integer,          intent(in)  :: iqty
character(len=*), intent(in)  :: vertstagger
integer,          intent(in)  :: ivar
real(r8),         intent(out) :: val
integer,          intent(out) :: istatus

integer  :: k, lev_top, lev_bottom
real(r8) :: zgrid, delta_z, zgrid_top, zgrid_bottom
real(r8) :: zgrid_upper, zgrid_lower
real(r8) :: val_top, val_bottom, frac_lev

! Presume the worst. Failure.
istatus    = 1
val        = MISSING_R8
delta_z    = MISSING_R8
frac_lev   = MISSING_R8
lev_top    = 0
lev_bottom = 0

!HK if ( vertstagger == 'ilev') then
!HK 
!HK    zgrid_bottom = x(get_dart_vector_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=1    ))
!HK    zgrid_top    = x(get_dart_vector_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=nilev))
!HK 
!HK    ! cannot extrapolate below bottom or beyond top ... fail ...
!HK    if ((zgrid_bottom > height) .or. (zgrid_top < height)) return
!HK 
!HK    ! Figure out what level is above/below, and by how much
!HK    h_loop_interface : do k = 2, nilev
!HK 
!HK       zgrid = x(get_dart_vector_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=k))
!HK 
!HK       if (height <= zgrid) then
!HK          lev_top    = k
!HK          lev_bottom = lev_top - 1
!HK          delta_z    = zgrid - x(get_dart_vector_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=lev_bottom))
!HK          frac_lev   = (zgrid - height)/delta_z
!HK          exit h_loop_interface
!HK       endif
!HK 
!HK    enddo h_loop_interface
!HK 
!HK elseif ( vertstagger == 'lev') then
!HK    ! Variable is on level midpoints, not ilevels.
!HK    ! Get height as the average of the ilevels.
!HK 
!HK    ! ilev index    1      2      3      4    ...  27    28    29
!HK    ! ilev value  -7.00, -6.50, -6.00, -5.50, ... 6.00, 6.50, 7.00 ;
!HK    !  lev value     -6.75, -6.25, -5.75, -5.25, ... 6.25, 6.75
!HK    !  lev index        1      2      3      4    ...  27    28
!HK 
!HK    !mid_level 1
!HK    zgrid_bottom = (x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=1)) + &
!HK                    x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=2))) / 2.0_r8
!HK 
!HK    !mid_level nlev
!HK    zgrid_top    = (x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=nilev-1)) + &
!HK                    x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=nilev))) / 2.0_r8
!HK 
!HK    ! cannot extrapolate below bottom or beyond top ... fail ...
!HK    if ((zgrid_bottom > height) .or. (zgrid_top < height)) return
!HK 
!HK    ! Figure out what level is above/below, and by how much
!HK    h_loop_midpoint: do k = 2, nilev-1
!HK 
!HK      lev_bottom = k-1
!HK      lev_top    = k
!HK 
!HK      zgrid_lower = &
!HK              (x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k-1 )) + &
!HK               x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k   ))) / 2.0_r8
!HK 
!HK      zgrid_upper = &
!HK              (x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k  )) + &
!HK               x(get_dart_vector_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k+1))) / 2.0_r8
!HK 
!HK      if (height  <= zgrid_upper) then
!HK         if (zgrid_upper == zgrid_lower) then ! avoid divide by zero
!HK            frac_lev = 0.0_r8  ! the fraction does not matter ...
!HK         else
!HK            delta_z  = zgrid_upper - zgrid_lower
!HK            frac_lev = (zgrid_upper - height)/delta_z
!HK         endif
!HK         exit h_loop_midpoint
!HK      endif
!HK 
!HK    enddo h_loop_midpoint
!HK 
!HK else
!HK    write(string1,*)'Unknown vertical stagger ',trim(vertstagger)
!HK    call error_handler(E_MSG,'vert_interp:', string1, source, revision, revdate )
!HK endif
!HK 
!HK ! Check to make sure we didn't fall through the h_loop ... unlikely (impossible?)
!HK if ( (frac_lev == MISSING_R8) .or. (lev_top == 0) .or. (lev_bottom == 0) ) then
!HK    write(string1,*)'Should not be here ... fell through loop.'
!HK    call error_handler(E_MSG,'vert_interp:', string1, source, revision, revdate )
!HK    return
!HK endif
!HK 
!HK istatus = 0 ! If we made it this far, it worked.
!HK 
!HK if (iqty == QTY_PRESSURE) then ! log-linear interpolation in height
!HK 
!HK    val_top    = plevs(lev_top)     !pressure at midpoint [Pa]
!HK    val_bottom = plevs(lev_bottom)  !pressure at midpoint [Pa]
!HK    val        = exp(frac_lev * log(val_bottom) + (1.0 - frac_lev) * log(val_top))
!HK 
!HK else ! simple linear interpolation in height
!HK 
!HK    val_top    = x(get_dart_vector_index(ivar, indx1=lon_index, indx2=lat_index, indx3=lev_top))
!HK    val_bottom = x(get_dart_vector_index(ivar, indx1=lon_index, indx2=lat_index, indx3=lev_bottom))
!HK    val        =     frac_lev *     val_bottom  + (1.0 - frac_lev) *     val_top
!HK endif

end subroutine vert_interp

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
! End of model_mod
!===============================================================================
end module model_mod
