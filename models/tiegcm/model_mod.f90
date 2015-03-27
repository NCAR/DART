! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

!-------------------------------------------------------------------------------
!
! Interface for HAO-TIEGCM
!
!-------------------------------------------------------------------------------

use        types_mod, only : r4, r8, MISSING_R8, MISSING_R4, PI, &
                             earth_radius, gravity, obstypelength

use time_manager_mod, only : time_type, set_calendar_type, set_time_missing,        &
                             set_time, get_time, print_time,                        &
                             set_date, get_date, print_date,                        &
                             operator(*),  operator(+), operator(-),                &
                             operator(>),  operator(<), operator(/),                &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_close_maxdist_init,                 &
                             get_close_obs_init, loc_get_close_obs => get_close_obs,&
                             set_location, get_location, query_location,            &
                             get_dist, vert_is_height, horiz_dist_only,             &
                             get_close_type, vert_is_undef, VERTISUNDEF,            &
                             VERTISPRESSURE, VERTISHEIGHT, vert_is_level

use    utilities_mod, only : file_exist, open_file, close_file, logfileunit,        &
                             error_handler, E_ERR, E_MSG, E_WARN, nmlfileunit,      &
                             do_output, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, nc_check, register_module,   &
                             file_to_text, find_textfile_dims, to_upper

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,           &
                             KIND_V_WIND_COMPONENT,           &
                             KIND_TEMPERATURE,                &! neutral temperature obs
                             KIND_PRESSURE,                   &! neutral pressure obs
                             KIND_MOLEC_OXYGEN_MIXING_RATIO,  &! neutral composition obs
                             KIND_1D_PARAMETER,               &
                             KIND_GEOPOTENTIAL_HEIGHT,        &
                             KIND_GEOMETRIC_HEIGHT,           &
                             KIND_VERTICAL_TEC,               &! total electron content
                             get_raw_obs_kind_index, get_raw_obs_kind_name

use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use mpi_utilities_mod,only : my_task_id

use typesizes
use netcdf

implicit none
private

!DART mandatory public interfaces
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

!TIEGCM specific routines
public ::   read_TIEGCM_restart, &
          update_TIEGCM_restart, &
          get_restart_file_name, &
          get_f107_value,        &
          test_interpolate

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   '$URL$'
character(len=32 ), parameter :: revision = '$Revision$'
character(len=128), parameter :: revdate  = '$Date$'

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

integer, parameter :: max_num_variables = 30
integer, parameter :: max_num_columns = 6
character(len=NF90_MAX_NAME) :: variables(max_num_variables * max_num_columns) = ' '

namelist /model_nml/ output_state_vector, tiegcm_restart_file_name, &
                     tiegcm_secondary_file_name, tiegcm_namelist_file_name, &
                     variables, debug, estimate_f10_7, assimilation_period_seconds

!-------------------------------------------------------------------------------
! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens ! ntime, [nlev,] nlat, nlon
   integer  :: rank
   integer  :: varsize     ! prod(dimlens(1:rank))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   character(len=obstypelength) :: verticalvar
   character(len=obstypelength) :: kind_string
   integer  :: dart_kind
   integer  :: xtype
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   real(r8) :: missingR8
   real(r4) :: missingR4
   logical  :: has_missing_value
   logical  :: update  ! is this a state variable (updateable)
   character(len=NF90_MAX_NAME) :: origin
end type progvartype

type(progvartype), dimension(max_num_variables) :: progvar
integer :: nfields  ! number of tiegcm variables in DART state

!-------------------------------------------------------------------------------
! define model parameters

integer                               :: nilev, nlev, nlon, nlat
real(r8),dimension(:),    allocatable :: lons, lats, levs, ilevs, plevs, pilevs
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

character(len=obstypelength) :: variable_table(max_num_variables, max_num_columns)

! include_vTEC = .true.  vTEC must be calculated from other vars
! include_vTEC = .false. just ignore vTEC altogether

logical  :: include_vTEC = .true.
logical  :: include_vTEC_in_state = .false.

! IMPORTANT: 1 D model parameters (e.g., F107) are read in from "tiegcm.nml"
! (note "estimate_f10_7" option is still under
! development by Tomoko Matsuo as of June 24, 2011)

real(r8)        :: f10_7
type(time_type) :: state_time ! module-storage declaration of current model time

integer                :: model_size
real(r8), allocatable  :: ens_mean(:)

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

interface prog_var_to_vector
   module procedure var1d_to_vector
   module procedure var2d_to_vector
   module procedure var3d_to_vector
   module procedure var4d_to_vector
end interface

interface apply_attributes
   module procedure apply_attributes_1D
   module procedure apply_attributes_2D
   module procedure apply_attributes_3D
   module procedure apply_attributes_4D
end interface

!===============================================================================
contains
!===============================================================================



subroutine static_init_model()
!-------------------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

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
call read_TIEGCM_secondary(tiegcm_secondary_file_name)

! error-check and convert namelist input to variable_table
! and the 'progvar' database in the scope of the module
call verify_variables(variables, nfields)

! Compute overall model size

model_size = progvar(nfields)%indexN

if (do_output() .and. (debug > 0)) then
   write(*,*) 'nlon  = ', nlon
   write(*,*) 'nlat  = ', nlat
   write(*,*) 'nlev  = ', nlev
   write(*,*) 'nilev = ', nilev
   write(*,*) 'model_size = ', model_size
   if (estimate_f10_7) &
   write(*,*) 'estimating f10.7 ... init value ',f10_7
endif

allocate (ens_mean(model_size))
ens_mean = MISSING_R8

! Might as well use the Gregorian Calendar
call set_calendar_type('Gregorian')

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


subroutine init_conditions(x)
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'no good way to specify initial conditions'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

x(:) = MISSING_R8 ! just to silence compiler messages

end subroutine init_conditions


!-------------------------------------------------------------------------------


subroutine adv_1step(x, time)
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a
! NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

write(string1,*) 'TIEGCM cannot be advanced as a subroutine from within DART.'
write(string2,*) 'check your input.nml setting of "async" and try again.'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

! just to silence compiler messages

x(:) = MISSING_R8
call print_time(time)

end subroutine adv_1step


!-------------------------------------------------------------------------------


function get_model_size()
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!-------------------------------------------------------------------------------


subroutine init_time(time)
! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'no good way to specify initial time'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate)

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time


!-------------------------------------------------------------------------------


subroutine model_interpolate(x, location, ikind, obs_val, istatus)
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The ikind variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.).

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: ikind
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

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
! The get_expected_gnd_gps_vtec() tries to interpolate KIND_GEOPOTENTIAL_HEIGHT
! when it does, this will kill it. 

if ( ikind == KIND_GEOPOTENTIAL_HEIGHT ) then
   write(string1,*)'KIND_GEOPOTENTIAL_HEIGHT currently unsupported'
   call error_handler(E_ERR,'model_interpolate',string1,source, revision, revdate)
endif

! Get the position

lon_lat_lev = get_location(location)
lon         = lon_lat_lev(1) ! degree
lat         = lon_lat_lev(2) ! degree

if(     vert_is_height(location) ) then
   height = lon_lat_lev(3)
elseif( vert_is_level(location) ) then
   level  = int(lon_lat_lev(3))
elseif( vert_is_undef(location) ) then
   height = MISSING_R8
   level  = -1
else
   which_vert = nint(query_location(location))
   write(string1,*) 'vertical coordinate type:',which_vert,' cannot be handled'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
endif

! Check to make sure vertical level is possible.
if (vert_is_level(location)) then
   if ((level < 1) .or. (level > nilev)) return
endif

if ((ikind == KIND_PRESSURE) .and. (vert_is_level(location))) then
   ! Some variables need plevs, some need pilevs
   ! We only need the height (aka level)
   ! the obs_def_upper_atm_mod.f90:get_expected_O_N2_ratio routines queries
   ! for the pressure at the model levels - EXACTLY - so ...
   ! FIXME ... at present ... the only time model_interpolate
   ! gets called with KIND_PRESSURE is to calculate density, which
   ! requires other variables that only live on the midpoints.
   ! I cannot figure out how to generically decide when to
   ! use plevs vs. pilevs 
   if (level <= nlev) then
      obs_val = plevs(level)
      istatus = 0
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
if (ikind == KIND_PRESSURE) then

   call vert_interp(x,lon_below,lat_below,height,ikind,'ilev',-1,val(1,1),istatus)
   if (istatus == 0) &
   call vert_interp(x,lon_below,lat_above,height,ikind,'ilev',-1,val(1,2),istatus)
   if (istatus == 0) &
   call vert_interp(x,lon_above,lat_below,height,ikind,'ilev',-1,val(2,1),istatus)
   if (istatus == 0) &
   call vert_interp(x,lon_above,lat_above,height,ikind,'ilev',-1,val(2,2),istatus)

   ! Now that we have the four surrounding points and their relative weights,
   ! actually perform the bilinear horizontal interpolation to get the value at the
   ! (arbitrary) desired location.

   if (istatus == 0) then
      do i = 1, 2
         a(i) = lon_fract * val(2, i) + (1.0_r8 - lon_fract) * val(1, i)
      end do
      obs_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)
   endif

   return
endif

! If it is not pressure ...
! FindVar_by_kind would fail with ikind == KIND_PRESSURE, so we have
! to calculate the pressure separately before this part.

ivar = FindVar_by_kind(ikind)
if (ivar < 0) return ! as a failure

! Now, need to find the values for the four corners

if ((progvar(ivar)%rank == 3) .or. (vert_is_undef(location))) then ! (time, lat, lon)
   ! time is always a singleton dimension

   val(1,1) = x(get_index(ivar, indx1=lon_below, indx2=lat_below, indx3=1))
   val(1,2) = x(get_index(ivar, indx1=lon_below, indx2=lat_above, indx3=1))
   val(2,1) = x(get_index(ivar, indx1=lon_above, indx2=lat_below, indx3=1))
   val(2,2) = x(get_index(ivar, indx1=lon_above, indx2=lat_above, indx3=1))
   istatus  = 0

elseif ((progvar(ivar)%rank == 4) .and. (vert_is_level(location))) then

   ! one use of model_interpolate is to allow other modules/routines
   ! the ability to 'count' the model levels. To do this, create observations
   ! with locations on model levels and 'interpolate' for KIND_GEOMETRIC_HEIGHT.
   ! When the interpolation fails, you've gone one level too far. 
   ! The only geometric height variable we have is ZG, and its a 4D variable.

   val(1,1) = x(get_index(ivar, indx1=lon_below, indx2=lat_below, indx3=level))
   val(1,2) = x(get_index(ivar, indx1=lon_below, indx2=lat_above, indx3=level))
   val(2,1) = x(get_index(ivar, indx1=lon_above, indx2=lat_below, indx3=level))
   val(2,2) = x(get_index(ivar, indx1=lon_above, indx2=lat_above, indx3=level))
   istatus = 0

elseif ((progvar(ivar)%rank == 4) .and. (vert_is_height(location))) then

   ! vert_interp() interpolates the state column to
   ! the same vertical height as the observation.
   ! THEN, it is the same as the 2D case.

   call vert_interp(x, lon_below, lat_below, height, ikind, &
             progvar(ivar)%verticalvar, ivar, val(1,1), istatus)
   if (istatus == 0) &
   call vert_interp(x, lon_below, lat_above, height, ikind, &
             progvar(ivar)%verticalvar, ivar, val(1,2), istatus)
   if (istatus == 0) &
   call vert_interp(x, lon_above, lat_below, height, ikind, &
             progvar(ivar)%verticalvar, ivar, val(2,1), istatus)
   if (istatus == 0) &
   call vert_interp(x, lon_above, lat_above, height, ikind, &
             progvar(ivar)%verticalvar, ivar, val(2,2), istatus)

else

   write(string1,*)trim(progvar(ivar)%varname),'has unsupported rank',progvar(ivar)%rank
   write(string2,*)trim(progvar(ivar)%varname),'or unsupported vertical system ',which_vert,&
                   'for that rank.'
   call error_handler(E_ERR,'model_interpolate', string1, &
              source, revision, revdate, text2=string2)

endif

! Now that we have the four surrounding points and their relative weights,
! actually perform the bilinear interpolation to get the value at the
! (arbitrary) desired location.

if (istatus == 0) then
   do i = 1, 2
      a(i) = lon_fract * val(2, i) + (1.0_r8 - lon_fract) * val(1, i)
   end do
   obs_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)
endif

end subroutine model_interpolate


!-------------------------------------------------------------------------------


function get_model_time_step()
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = time_step

end function get_model_time_step


!-------------------------------------------------------------------------------


subroutine get_state_meta_data(index_in, location, var_kind)
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_kind

integer  :: remainder
integer  :: relindx, absindx
integer  :: lon_index, lat_index, lev_index
integer  :: ivar, seconds, days
real(r8) :: height, longitude

if ( .not. module_initialized ) call static_init_model

ivar   = Find_Variable_by_index(index_in,'get_state_meta_data')
relindx = index_in - progvar(ivar)%index1 + 1

if     (progvar(ivar)%rank == 0) then  ! scalars ... no location
   ! f10_7 is most accurately located at local noon at equator.
   ! 360.0 degrees in 86400 seconds, 43200 secs == 12:00 UTC == longitude 0.0

   call get_time(state_time, seconds, days)
   longitude = 360.0_r8 * real(seconds,r8) / 86400.0_r8 - 180.0_r8
   if (longitude < 0.0_r8) longitude = longitude + 360.0_r8
   location = set_location(longitude, 0.0_r8,  400000.0_r8, VERTISUNDEF)

elseif (progvar(ivar)%rank == 1) then  ! time?
   write(string1,*)trim(progvar(ivar)%varname),'has unsupported shape (1D)'
   write(string2,*)'dimension ('//trim(progvar(ivar)%dimnames(1))// &
                   ') ... unknown location'
   call error_handler(E_ERR,'get_state_meta_data', string1, &
              source, revision, revdate, text2=string2)

elseif (progvar(ivar)%rank == 2) then  ! something, and time? lat & lon but no time?
   write(string1,*)trim(progvar(ivar)%varname), &
                   'has unsupported shape (2D) ... unknown location'
   write(string2,*)'dimension 1 = ('//trim(progvar(ivar)%dimnames(1))//')'
   write(string3,*)'dimension 2 = ('//trim(progvar(ivar)%dimnames(2))//')'
   call error_handler(E_ERR,'get_state_meta_data', string1, &
              source, revision, revdate, text2=string2, text3=string3)

elseif (progvar(ivar)%rank == 3) then  ! [Fortran ordering] = lon, lat, time
   ! The time dimension is always length 1, so it really doesn't matter.

   lat_index = 1 + (relindx - 1) / nlon ! relies on integer arithmetic
   lon_index = relindx - (lat_index-1) * nlon

   if (do_output() .and. (debug > 3)) then
      absindx = get_index(ivar, indx1=lon_index, indx2=lat_index, knownindex=index_in)
   endif

   if (trim(progvar(ivar)%varname) == 'VTEC') then
      ! assign arbitrary height to allow for localization
      location = set_location(lons(lon_index), lats(lat_index), 300000.0_r8, VERTISHEIGHT)
   else
      location = set_location(lons(lon_index), lats(lat_index), 0.0_r8, VERTISUNDEF)
   endif

elseif (progvar(ivar)%rank == 4) then  ! [Fortran ordering] = lon, lat, lev, time

   lev_index = 1 + (relindx - 1) / (nlon * nlat)
   remainder = relindx - (lev_index-1) * nlon * nlat
   lat_index = 1 + (remainder - 1) / nlon
   lon_index = remainder - (lat_index-1) * nlon
   height    = get_height(ivar, lon_index, lat_index, lev_index)

   if (do_output() .and. (debug > 3)) then
      absindx = get_index( ivar, indx1=lon_index, indx2=lat_index, indx3=lev_index, knownindex=index_in)
   endif

   location  = set_location(lons(lon_index), lats(lat_index), height, VERTISHEIGHT)

else
   write(string1,*)'Problem with DART variable ',trim(progvar(ivar)%varname)
   write(string2,*)'has unsupported number (',progvar(ivar)%rank,') of dimensions.'
   call error_handler(E_ERR,'get_state_meta_data', string1, &
                      source, revision, revdate, text2=string2)
endif

! If the type is wanted, return it
if(present(var_kind)) var_kind = progvar(ivar)%dart_kind

end subroutine get_state_meta_data


!-------------------------------------------------------------------------------


subroutine end_model()
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if ( .not. module_initialized ) call static_init_model


end subroutine end_model


!-------------------------------------------------------------------------------


function nc_write_model_atts( ncid ) result (ierr)
! TJH 20 Dec 2013 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
!
! TJH 20 Dec 2013 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...

integer, intent(in)  :: ncid      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

integer :: myndims
integer :: ivar, VarID
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
character(len=NF90_MAX_NAME) :: varname

integer :: lonDimID, latDimID, levDimID, ilevDimID
integer :: lonVarID, latVarID, levVarID, ilevVarID

!-------------------------------------------------------------------------------
! variables for the namelist output
!-------------------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_tiegcm_namelist

!-------------------------------------------------------------------------------
! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.
!-------------------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

real(r8), allocatable :: temp_lon(:)

integer :: i

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncid refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_Inquire(ncid, nDimensions, nVariables, &
              nAttributes, unlimitedDimID), 'nc_write_model_atts','inquire')
call nc_check(nf90_Redef(ncid),'nc_write_model_atts','redef')

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncid, name='NMLlinelen', dimid = linelenDimID), &
       'nc_write_model_atts', 'inq_dimid NMLlinelen')
call nc_check(nf90_inq_dimid(ncid=ncid, name='copy', dimid=MemberDimID),&
       'nc_write_model_atts', 'inq_dimid copy')
call nc_check(nf90_inq_dimid(ncid=ncid, name='time', dimid=  TimeDimID),&
       'nc_write_model_atts', 'inq_dimid time')

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
                     ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncid, name='StateVariable', &
                        len=model_size, dimid = StateVarDimID),&
       'nc_write_model_atts', 'state def_dim')

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,str1    ),&
       'nc_write_model_atts', 'creation put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model_source'  ,source  ),&
       'nc_write_model_atts', 'source put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model_revision',revision),&
       'nc_write_model_atts', 'revision put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model_revdate' ,revdate ),&
       'nc_write_model_atts', 'revdate put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model','TIEGCM'         ),&
       'nc_write_model_atts', 'model put')

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

if (has_tiegcm_namelist) then
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncid, name='tiegcmNMLnlines', &
          len = nlines, dimid = nlinesDimID), &
          'nc_write_model_atts', 'def_dim tiegcmNMLnlines')

   call nc_check(nf90_def_var(ncid,name='tiegcm_nml', xtype=nf90_char, &
          dimids = (/ linelenDimID, nlinesDimID /), varid=nmlVarID), &
          'nc_write_model_atts', 'def_var tiegcm_namelist')

   call nc_check(nf90_put_att(ncid, nmlVarID, 'long_name', &
          'contents of '//trim(tiegcm_namelist_file_name)), &
          'nc_write_model_atts', 'put_att tiegcm_namelist')

endif

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncid,name='StateVariable', xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID), &
              'nc_write_model_atts', 'statevariable def_var')
   call nc_check(nf90_put_att(ncid, StateVarVarID, 'long_name', 'State Variable ID'), &
              'nc_write_model_atts', 'statevariable long_name')
   call nc_check(nf90_put_att(ncid, StateVarVarID, 'units',     'indexical'), &
              'nc_write_model_atts', 'statevariable units')
   call nc_check(nf90_put_att(ncid, StateVarVarID, 'valid_range', (/ 1, model_size /)), &
              'nc_write_model_atts', 'statevariable valid_range')

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncid, name='state', xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid=StateVarID), 'nc_write_model_atts', 'state def_var')
   call nc_check(nf90_put_att(ncid, StateVarID, 'long_name', 'model state or fcopy'), &
              'nc_write_model_atts', 'state long_name')

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncid), 'nc_write_model_atts', 'state enddef')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncid, StateVarVarID, (/ (i,i=1,model_size) /) ), &
              'nc_write_model_atts', 'state put_var')

else

   !----------------------------------------------------------------------------
   ! We need to process the reshaped variables.
   !----------------------------------------------------------------------------
   ! Define the dimensions IDs
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid= ncid ,  name='lon', &
             & len =  nlon,  dimid= lonDimID), 'nc_write_model_atts lon')
   call nc_check(nf90_def_dim(ncid= ncid ,  name='lat', &
             & len =  nlat,  dimid= latDimID), 'nc_write_model_atts lat')
   call nc_check(nf90_def_dim(ncid= ncid ,  name='lev', &
             & len =  nlev,  dimid= levDimID), 'nc_write_model_atts lev')
   call nc_check(nf90_def_dim(ncid= ncid ,  name='ilev', &
             & len = nilev,  dimid=ilevDimID), 'nc_write_model_atts ilev')

   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid, name='lon', &
             & xtype=nf90_double, dimids=lonDimID, varid=lonVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, lonVarID, &
             & 'long_name', 'geographic longitude (-west, +east)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, lonVarID, 'units', 'degrees_east'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, lonVarID, 'valid_range', &
             & (/ -180.0_r8, 180.0_r8 /)),'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid, name='lat', &
             & xtype=nf90_double, dimids=latDimID, varid=latVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, latVarID, &
             & 'long_name', 'geographic latitude (-south +north)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, latVarID, 'units', 'degrees_north'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, latVarID, 'valid_range', &
             & (/ -90.0_r8, 90.0_r8 /)),'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid, name='lev', &
             & xtype=nf90_double, dimids=levDimID, varid=levVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'long_name', 'midpoint levels'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'short_name', 'ln(p0/p)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'units', ''),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'positive', 'up'),&
             'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid, name='ilev', &
             & xtype=nf90_double, dimids=ilevDimID, varid=ilevVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'long_name', 'interface levels'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'short_name', 'ln(p0/p)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'units', ''),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'positive', 'up'),&
             'nc_write_model_atts')

   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and their Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncid, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncid,name=trim(varname),xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(varname)//' def_var' )

      call nc_check(nf90_put_att(ncid,VarID,'long_name',trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(varname)//' put_att long_name' )
      call nc_check(nf90_put_att(ncid,VarID,'units',    trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(varname)//' put_att units' )
      call nc_check(nf90_put_att(ncid,VarID,'DART_kind',trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(varname)//' put_att dart_kind' )
   enddo

   call nc_check(nf90_enddef(ncid), 'nc_write_model_atts', 'reshaped enddef')

   !----------------------------------------------------------------------------
   ! Fill the variables
   !----------------------------------------------------------------------------
   allocate(temp_lon(size(lons)))  ! Convert longitudes back to a monotonic array
   temp_lon = lons                 ! instead of [180 ... 355,0 ... 175]
   where(temp_lon >= 180.0_r8) temp_lon = temp_lon - 360.0_r8

   call nc_check(nf90_put_var(ncid, lonVarID, temp_lon),'nc_write_model_atts','put_var lons')
   call nc_check(nf90_put_var(ncid, latVarID, lats),'nc_write_model_atts','put_var lats')
   call nc_check(nf90_put_var(ncid, levVarID, levs),'nc_write_model_atts','put_var levs')
   call nc_check(nf90_put_var(ncid,ilevVarID,ilevs),'nc_write_model_atts','put_var ilevs')

   deallocate(temp_lon)

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_tiegcm_namelist) then
   call file_to_text(tiegcm_namelist_file_name, textblock)
   call nc_check(nf90_put_var(ncid, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var nmlVarID')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncid), 'nc_write_model_atts', 'sync')
if (do_output() .and. (debug > 1)) write (*,*) 'nc_write_model_atts: netCDF file ', ncid, ' is synched '

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!-------------------------------------------------------------------------------


function nc_write_model_vars( ncid, statevec, copyindex, timeindex ) result (ierr)
! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...

integer,                intent(in) :: ncid
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncid refers to an open netCDF file,
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_Inquire(ncid, nDimensions, nVariables, &
                  nAttributes, unlimitedDimID), 'nc_write_model_vars', 'inquire')

call nc_check(nf90_inq_dimid(ncid, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy ')

call nc_check(nf90_inq_dimid(ncid, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time ')

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncid, 'state', StateVarID), &
          'nc_write_model_vars', 'state inq_varid' )
   call nc_check(NF90_put_var(ncid, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)), &
          'nc_write_model_vars', 'state put_var')

else

   !----------------------------------------------------------------------------
   ! We need to process the reshaped variables.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields

      varname = trim(progvar(ivar)%varname)
      string2 = trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried.
      ! This requires that Time is the unlimited dimension (the last one in Fortran),
      ! and that 'copy' is the second-to-last. The variables declared in the DART
      ! diagnostic files are required to have the same shape as in the source
      ! restart file. If Time is present there, it must also be the 'last' one.

      call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      ! This block creates the netCDF control information for the hyperslab indexing
      mystart(:) = 1
      mycount(:) = 1
      DimCheck : do i = 1,ncNdims

         write(string1,'(A,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncid,dimIDs(i),name=dimnames(i),len=dimlen), &
               'nc_write_model_vars', trim(string1))

         ! Dont care about the length of these two, setting to 1 later.
         if (dimIDs(i) == CopyDimID) cycle DimCheck
         if (dimIDs(i) == TimeDimID) cycle DimCheck

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*)trim(string2),' dim/dimlen ',i,dimlen, &
                                     ' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         mycount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if (do_output() .and. (debug > 3)) then
         write(string1,*)'nc_write_model_vars',trim(progvar(ivar)%varname),' file id ',ncid
         write(string2,*)'nc_write_model_vars '//trim(varname)//' start is ',mystart(1:ncNdims)
         write(string3,*)'nc_write_model_vars '//trim(varname)//' count is ',mycount(1:ncNdims)
         call error_handler(E_MSG,'nc_write_model_vars',string1,text2=string2,text3=string3)
      endif

      ! Now that we have the hyperslab indices ... it is easy.

      call vector_to_prog_var(statevec, ivar, ncid, VarID, ncNdims, &
                              mystart, mycount, limit=.false.)

   enddo

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

if (do_output() .and. (debug > 1)) &
   write (*,*) 'nc_write_model_vars: Finished filling variables '

call nc_check(nf90_sync(ncid), 'nc_write_model_vars', 'sync')

if (do_output() .and. (debug > 1)) &
   write (*,*) 'nc_write_model_vars: netCDF file is synched '

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars


!-------------------------------------------------------------------------------


subroutine pert_model_state(state, pert_state, interf_provided)
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)    :: state(:)
real(r8), intent(out)   :: pert_state(:)
logical,  intent(out)   :: interf_provided

integer                 :: i, variable_type
type(location_type)     :: temp_loc

if ( .not. module_initialized ) call static_init_model

! An interface is provided
interf_provided = .true.

! If first call initialize random sequence
! CAUTION: my_task_id is NOT ensemble member number
! For example, my_task_id will be in [0,N-1]
! if a single instance of the model using N MPI tasks.

if(first_pert_call) then
   call init_random_seq(random_seq,my_task_id())
   first_pert_call = .false.
endif

do i = 1, get_model_size()
   call get_state_meta_data(i, temp_loc, variable_type)
   if(variable_type == KIND_1D_PARAMETER) then
      pert_state(i) = random_gaussian(random_seq,state(i),20.0_r8)
   else
      ! FIXME ... is this really what you want to do - no variability in the states.
      pert_state(i) = state(i)
   endif
end do

end subroutine pert_model_state


!-------------------------------------------------------------------------------


subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                            num_close, close_ind, dist)
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

type(get_close_type), intent(in)     :: gc
type(location_type),  intent(inout)  :: base_obs_loc, obs_loc(:)
integer,              intent(in)     :: base_obs_kind, obs_kind(:)
integer,              intent(out)    :: num_close, close_ind(:)
real(r8),             intent(out)    :: dist(:)

integer                              :: k, t_ind

! Finds all the locations or observations that are close.
call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                                           num_close, close_ind, dist)

! Make the ZG part of the state vector far from everything so it does not get updated.
! Scroll through all the obs_loc(:) and obs_kind(:) elements

do k = 1,num_close
   t_ind  = close_ind(k)
   if (obs_kind(t_ind) == KIND_GEOMETRIC_HEIGHT) then
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
   do k = 1, num_close

      t_ind  = close_ind(k)

      ! This was part of an experimental setup - if the distance is LARGE,
      ! the state variables will not be updated. By increasing the distance,
      ! the values in the DART vector will remain unchanged. If you allow
      ! the DART vector to be updated, the posterior observation operators
      ! will be impacted - which is usually the desire. You can then avoid
      ! impacting the tiegcm forecast through the input.nml 'NO_COPY_BACK' feature.

      if (    (obs_kind(t_ind) == KIND_MOLEC_OXYGEN_MIXING_RATIO) &
         .or. (obs_kind(t_ind) == KIND_U_WIND_COMPONENT) &
         .or. (obs_kind(t_ind) == KIND_V_WIND_COMPONENT) &
         .or. (obs_kind(t_ind) == KIND_TEMPERATURE) ) then
      !  dist(k) = 10.0_r8 * PI

      elseif  (obs_kind(t_ind) == KIND_1D_PARAMETER) then
         ! f10_7 is given a location of latitude 0.0 and the longitude
         ! of local noon. By decreasing the distance from the observation
         ! to the dynamic f10_7 location we are allowing the already close
         ! observations to have a larger impact in the parameter estimation.
         ! 0.25 is heuristic. The 'close' observations have already been 
         ! determined by the cutoff. Changing the distance here does not
         ! allow more observations to impact anything.
         dist(k) = dist(k)*0.25_r8
      endif

   enddo
endif

end subroutine get_close_obs


!-------------------------------------------------------------------------------


subroutine ens_mean_for_model(filter_ens_mean)
! Stores provided ensemble mean within the module for later use

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

ens_mean = filter_ens_mean

end subroutine ens_mean_for_model


!===============================================================================
! TIEGCM public routines
!===============================================================================


subroutine read_TIEGCM_restart(file_name, statevec, model_time)
!
! Read TIEGCM restart file fields and pack it into a DART vector
!
character(len=*),       intent(in)  :: file_name
real(r8), dimension(:), intent(out) :: statevec
type(time_type),        intent(out) :: model_time

integer :: ncid, VarID
integer :: time_dimlen, dimlen, ncNdims
integer :: LonDimID, LatDimID, LevDimID, IlevDimID, TimeDimID
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount

real(r8), allocatable, dimension(:)       :: temp1D
real(r8), allocatable, dimension(:,:)     :: temp2D
real(r8), allocatable, dimension(:,:,:)   :: temp3D
real(r8), allocatable, dimension(:,:,:,:) :: temp4D

integer :: i, ivar

if ( .not. module_initialized ) call static_init_model

if( .not. file_exist(file_name)) then
   write(string1,*)trim(file_name)//' does not exist.'
   call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call error_handler(E_MSG,'read_TIEGCM_restart:','reading restart ['//trim(file_name)//']')

! Open the netCDF file and then read all the static information.

call nc_check( nf90_open( file_name, NF90_NOWRITE, ncid ), &
                              'read_TIEGCM_restart', 'open')

!... check for matching dimensions
call nc_check( nf90_inq_dimid(ncid, 'lon', LonDimID), &
         'read_TIEGCM_restart', 'inq_dimid lon')
call nc_check( nf90_inquire_dimension(ncid, LonDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension lon')
if (dimlen .ne. nlon) then
  write(string1, *) trim(file_name), ' dim_lon = ',dimlen, ' DART expects ',nlon
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'lat', LatDimID), &
         'read_TIEGCM_restart', 'inq_dimid lat')
call nc_check( nf90_inquire_dimension(ncid, LatDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension lat')
if (dimlen .ne. nlat) then
  write(string1, *) trim(file_name), ' dim_lat = ',dimlen, ' DART expects ',nlat
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'lev', LevDimID), &
         'read_TIEGCM_restart', 'inq_dimid lev')
call nc_check( nf90_inquire_dimension(ncid, LevDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension lev')
if (dimlen .ne. (nlev+1)) then
  write(string1, *) trim(file_name), ' dim_lev = ',dimlen, ' DART expects ',nlev+1
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'ilev', IlevDimID), &
         'read_TIEGCM_restart', 'inq_dimid ilev')
call nc_check( nf90_inquire_dimension(ncid, IlevDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension ilev')
if (dimlen .ne. nilev) then
  write(string1, *) trim(file_name), ' dim_ilev = ',dimlen, ' DART expects ',nilev
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'time', TimeDimID), &
         'read_TIEGCM_restart', 'inq_dimid time')
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'read_TIEGCM_restart', 'inquire_dimension time')

! Loop over all variables in DART state.
! Make sure the shape of the variable matches what we expect.
! Get the hyperslab with the most current (last) time.

ReadVariable: do ivar = 1,nfields

   string2 = trim(file_name)//' '//trim(progvar(ivar)%varname)

   if (trim(progvar(ivar)%varname) == 'f10_7') then
      statevec( progvar(ivar)%index1 ) = f10_7
      cycle ReadVariable
   endif

   ! ZG is already a module variable that is required.
   if (trim(progvar(ivar)%varname) == 'ZG') then
      call prog_var_to_vector(ivar, ZG, statevec)
      ! deallocate(ZG)
      ! FIXME After the create_vtec() rewrite, ZG should be deallocated.
      ! It can be accessed directly from the DART state
      ! maybe should also be removed from the read_tiegcm_seconday() routine
      cycle ReadVariable
   endif

   if (trim(progvar(ivar)%varname) == 'VTEC') then
      allocate(temp2D(nlon,nlat))
      call create_vtec(ncid, time_dimlen, temp2D)
      call prog_var_to_vector(ivar, temp2D, statevec)
      deallocate(temp2D)
      cycle ReadVariable
   endif

   call nc_check(nf90_inq_varid(ncid, trim(progvar(ivar)%varname), VarID), &
          'read_TIEGCM_restart', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'read_TIEGCM_restart', 'inquire '//trim(string2))

   if (ncNdims /= progvar(ivar)%rank) then
      write(string1,*)trim(string2),'has ',ncNDims,'dimensions.'
      write(string3,*)'same thing in DART has ',progvar(ivar)%rank,' dimensions.'
      call error_handler(E_ERR, 'read_TIEGCM_restart', string1, &
                      source, revision, revdate, text2=string3)
   endif

   mystart(:) = 1
   mycount(:) = 1
   DimCheck : do i = 1,progvar(ivar)%rank

      write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
               'read_TIEGCM_restart', trim(string1))

      ! Check the shape of the variable
      if ( dimIDs(i) == LevDimID ) dimlen = dimlen - 1

      if ( dimIDs(i) == TimeDimID ) then
         dimlen = 1
      elseif ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' has dim/dimlen ',i,dimlen
         write(string3,*)' expected dimlen of ', progvar(ivar)%dimlens(i)
         call error_handler(E_ERR, 'read_TIEGCM_restart', string1, &
                         source, revision, revdate, text2=string3)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = time_dimlen
   where(dimIDs == TimeDimID) mycount = 1
   where(dimIDs ==  LonDimID) mycount = nlon
   where(dimIDs ==  LatDimID) mycount = nlat
   where(dimIDs ==  LevDimID) mycount = nlev
   where(dimIDs == iLevDimID) mycount = nilev

   if (do_output() .and. (debug > 3)) then
      write(*,*)'read_TIEGCM_restart:',trim(string2)
      write(*,*)'read_TIEGCM_restart: start ',mystart(1:ncNdims)
      write(*,*)'read_TIEGCM_restart: count ',mycount(1:ncNdims)
      write(*,*)
   endif

   if     (progvar(ivar)%rank == 1) then
      allocate(temp1D(mycount(1)))
      call nc_check(nf90_get_var(ncid, VarID, values=temp1D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call apply_attributes(ivar, ncid, VarID, temp1D)
      call prog_var_to_vector(ivar, temp1D, statevec)
      deallocate(temp1D)

   elseif (progvar(ivar)%rank == 2) then
      allocate(temp2D(mycount(1), mycount(2)))
      call nc_check(nf90_get_var(ncid, VarID, values=temp2D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call apply_attributes(ivar, ncid, VarID, temp2D)
      call prog_var_to_vector(ivar, temp2D, statevec)
      deallocate(temp2D)

   elseif (progvar(ivar)%rank == 3) then
      allocate(temp3D(mycount(1), mycount(2), mycount(3)))
      call nc_check(nf90_get_var(ncid, VarID, values=temp3D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call apply_attributes(ivar, ncid, VarID, temp3D)
      call prog_var_to_vector(ivar, temp3D, statevec)
      deallocate(temp3D)

   elseif (progvar(ivar)%rank == 4) then
      allocate(temp4D(mycount(1), mycount(2), mycount(3), mycount(4)))
      call nc_check(nf90_get_var(ncid, VarID, values=temp4D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call apply_attributes(ivar, ncid, VarID, temp4D)
      call prog_var_to_vector(ivar, temp4D, statevec)
      deallocate(temp4D)
   else
      write(string1,*) trim(string2),' has unsupported number of dims (', &
                                     progvar(ivar)%rank,')'
      call error_handler(E_ERR, 'read_TIEGCM_restart', string1, &
                         source, revision, revdate)
   endif

enddo ReadVariable

! Convert the last year/doy/hour/minute to a dart time.
model_time = get_state_time(trim(file_name), ncid)
state_time = model_time   ! state_time is scoped for entire module

call nc_check( nf90_close(ncid), 'read_TIEGCM_restart', 'close')

end subroutine read_TIEGCM_restart


!-------------------------------------------------------------------------------


subroutine update_TIEGCM_restart(statevec, filename, dart_time)
! Updates TIEGCM restart file fields with the reshaped variables
! in the DART state vector. The variable MISSING attribute has to be
! reinstated ... etc.

real(r8), dimension(:), intent(in) :: statevec
character(len=*),       intent(in) :: filename
type(time_type),        intent(in) :: dart_time

integer :: i, ivar
integer :: ncid, VarID, TimeDimID, LevDimID
integer :: dimlen, ncNdims
integer :: timeindex
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimnames

type(time_type) :: file_time

if ( .not. module_initialized ) call static_init_model

if( .not. file_exist(filename)) then
   write(string1,*) trim(filename),' not available.'
   call error_handler(E_ERR,'update_TIEGCM_restart',string1,source,revision,revdate)
endif

if (do_output()) &
   call error_handler(E_MSG,'update_TIEGCM_restart','opening',text2=filename)

call nc_check(nf90_open(trim(filename),NF90_WRITE,ncid),'update_TIEGCM_restart','open')

call SanityCheck(trim(filename), ncid, TimeDimID=TimeDimID, ntimes=timeindex)

call nc_check(nf90_inq_dimid(ncid, 'lev', LevDimID), &
                   'update_TIEGCM_restart', 'inq_dimid lev'//trim(filename))

file_time = get_state_time(trim(filename), ncid, timeindex)

if ( file_time /= dart_time ) then
   call print_time(dart_time,'DART   current time',logfileunit)
   call print_time(file_time,'TIEGCM current time',logfileunit)
   call print_time(dart_time,'DART   current time')
   call print_time(file_time,'TIEGCM current time')
   write(string1,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'update_TIEGCM_restart',string1,source,revision,revdate)
endif

if (do_output() .and. (debug > 0)) then
   call print_date(file_time,'update_TIEGCM_restart: date in restart file '//trim(filename))
   call print_time(file_time,'update_TIEGCM_restart: time in restart file '//trim(filename))
   call print_date(file_time,'update_TIEGCM_restart: date in restart file '//trim(filename),logfileunit)
   call print_time(file_time,'update_TIEGCM_restart: time in restart file '//trim(filename),logfileunit)
endif

UPDATE : do ivar = 1,nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   ! Skip the variables that are not supposed to be updated.
   if ( .not. progvar(ivar)%update ) cycle UPDATE

   ! Ensure netCDF variable is conformable with DART progvar quantity.
   ! At the same time, create the mystart(:) and mycount(:) arrays.

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
                 'update_TIEGCM_restart', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
                 'update_TIEGCM_restart', 'inquire '//trim(string2))

   mystart(:) = 1
   mycount(:) = 1
   DimCheck : do i = 1,ncNdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid,dimIDs(i),name=dimnames(i),len=dimlen), &
                    'update_TIEGCM_restart', string1)

      ! Dont care about the length of Time Dimension, setting to 1 later.
      if (dimIDs(i) == TimeDimID) cycle DimCheck
      if (dimIDs(i) == LevDimID) dimlen = dimlen - 1

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen
         write(string2,*)'expected it to be ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR, 'update_TIEGCM_restart', string1, &
                      source, revision, revdate, text2=string2)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   if (dimIDs(ncNdims) /= TimeDimID) then
      write(string1,*) trim(string2),&
         ' required to have "Time" as the last/unlimited dimension'
      write(string2,*)' last dimension is ',trim(dimnames(ncNdims))
      call error_handler(E_ERR, 'update_TIEGCM_restart', string1, &
                   source, revision, revdate, text2=string2)
   endif

   where(dimIDs == TimeDimID) mystart = timeindex
   where(dimIDs == TimeDimID) mycount = 1

   if (do_output() .and. (debug > 3)) then
      write(*,*)'update_TIEGCM_restart '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(*,*)'update_TIEGCM_restart '//trim(varname)//' count is ',mycount(1:ncNdims)
      do i = 1,NcNdims
         write(*,*)'update_TIEGCM_restart '//trim(varname)//' dim ',&
                   i,' is ',trim(dimnames(i))
      enddo
   endif

   ! When limit=.true. the top layer of the variables with 'lev' as the 
   ! vertical coordinate are not updated - and - the range of values for
   ! each variable is limited as per the namelist settings.
   call vector_to_prog_var(statevec, ivar, ncid, VarID, ncNdims, &
                           mystart, mycount, limit=.true.)

enddo UPDATE

call nc_check(nf90_sync( ncid), 'update_TIEGCM_restart', 'sync')
call nc_check(nf90_close(ncid), 'update_TIEGCM_restart', 'close')

end subroutine update_TIEGCM_restart


!-------------------------------------------------------------------------------


function get_restart_file_name()
character(len=256) :: get_restart_file_name

if ( .not. module_initialized ) call static_init_model
get_restart_file_name = adjustl(tiegcm_restart_file_name)

end function get_restart_file_name


!-------------------------------------------------------------------------------


function get_f107_value(x)
real(r8), dimension(:), intent(in) :: x
real(r8)                           :: get_f107_value

integer :: ivar

if ( .not. module_initialized ) call static_init_model

! If the f10_7 is part of the DART state, return that value.
! If it is not part of the DART state, just return the value from
! the TIEGCM namelist.

if (estimate_f10_7) then

   VARLOOP : do ivar = 1,nfields
      if (progvar(ivar)%varname == 'f10_7') then
         get_f107_value = x(progvar(ivar)%index1)
         return
      endif
   enddo VARLOOP

   call error_handler(E_ERR,'get_f107_value', 'no f10_7 in DART state', &
        source, revision, revdate)

else
   get_f107_value = f10_7
endif

end function get_f107_value


!-------------------------------------------------------------------------------


subroutine test_interpolate(x, locarray)
! This is used to debug other functions. Used by check_model_mod - ONLY.
! Not use in any other DART program.

real(r8), dimension(:), intent(in) :: x
real(r8), dimension(3), intent(in) :: locarray

type(location_type) :: location
real(r8) :: obs_val
integer  :: istatus

location = set_location(locarray(1), locarray(2), locarray(3), VERTISHEIGHT)

call model_interpolate(x, location, KIND_PRESSURE, obs_val, istatus)

write(*,*)'test_interpolate: PRESSURE value at ',locarray, &
          ' is ',obs_val,' status is ',istatus,' (0 is good)'

call model_interpolate(x, location, KIND_VERTICAL_TEC, obs_val, istatus)

write(*,*)'test_interpolate: VERTICAL_TEC value at ',locarray, &
          ' is ',obs_val,' status is ',istatus,' (0 is good)'

end subroutine test_interpolate


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


subroutine read_TIEGCM_secondary(file_name)
! Read TIEGCM geometric height (ZG) from a tiegcm secondary output file

character(len=*), intent(in):: file_name

integer :: ncid
integer :: TimeDimID, time_dimlen, VarID

real(r8) :: spvalR8, spvalR4

if( .not. file_exist(file_name)) then
  write(string1,*) trim(file_name),' not available.'
  call error_handler(E_ERR,'read_TIEGCM_secondary',string1,source,revision,revdate)
endif

call error_handler(E_MSG,'read_TIEGCM_secondary:', &
           'reading secondary ['//trim(file_name)//']')

call nc_check(nf90_open(file_name, NF90_NOWRITE, ncid), 'read_TIEGCM_secondary', 'open')

call SanityCheck(file_name, ncid, TimeDimID=TimeDimID, ntimes=time_dimlen)

allocate(ZG(nlon,nlat,nilev)) ! comes from module storage

!... actually read the target variable
call nc_check(nf90_inq_varid(ncid, 'ZG', VarID), 'read_TIEGCM_secondary', 'inq_varid ZG')
call nc_check(nf90_get_var(ncid, VarID, values=ZG,  &
                         start = (/ 1, 1, 1, time_dimlen /),    &
                         count = (/ nlon, nlat, nilev, 1 /)),    &
                         'read_TIEGCM_secondary', 'get_var ZG')

if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
   where(ZG == spvalR8) ZG = MISSING_R8
endif

! ... check units and convert them to meters if need be.
if (nf90_get_att(ncid, VarID, 'units' , string1) == NF90_NOERR) then
   if(trim(string1) == 'cm') then
      call error_handler(E_MSG,'read_TIEGCM_secondary:', &
          'Converting ZG from cm to meters.')
      where(ZG /= MISSING_R8) ZG = ZG/100.0_r8
   elseif(trim(string1) == 'm') then
      call error_handler(E_MSG,'read_TIEGCM_secondary:', &
          'ZG already in meters.')
   else
      call error_handler(E_ERR,'read_TIEGCM_secondary', &
          'ZG has unknown units',source, revision, revdate,text2=string1)
   endif
endif

call nc_check(nf90_close(ncid),'read_TIEGCM_secondary', 'close')

end subroutine read_TIEGCM_secondary


!-------------------------------------------------------------------------------


subroutine verify_variables( variables, ngood )
! This routine checks the user input against the variables available in the
! input netcdf file to see if it is possible to construct the DART state vector
! specified by the input.nml:model_nml:variables  variable.
!
! There is an additional complication that if the user requests VTEC to be part
! of the state, there are good scientific reasons to include the variables
! that are used to derive VTEC to be part of the DART state. So, if VTEC is
! included, then there are other variables that get included even if not
! explicitly mentioned in model_nml:variables

character(len=*), intent(in)  :: variables(:)
integer,          intent(out) :: ngood

integer  :: i, nrows, ncols
integer  :: ivar, index1, indexN, varsize
integer  :: ncid, ncid1, ncid2, ncerr, VarID, dimlen
real(r8) :: spvalR8, spvalR4

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: filename
character(len=NF90_MAX_NAME) :: state_or_aux
character(len=obstypelength) :: dimname

nrows = size(variable_table,1)
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

ngood = 0
MyLoop : do i = 1, nrows

   varname      = trim(variables(ncols*i - 5))
   dartstr      = trim(variables(ncols*i - 4))
   minvalstring = trim(variables(ncols*i - 3))
   maxvalstring = trim(variables(ncols*i - 2))
   filename     = trim(variables(ncols*i - 1))
   state_or_aux = trim(variables(ncols*i    ))

   call to_upper(filename)
   call to_upper(state_or_aux)

   variable_table(i,VT_VARNAMEINDX) = trim(varname)
   variable_table(i,VT_KINDINDX)    = trim(dartstr)
   variable_table(i,VT_MINVALINDX)  = trim(minvalstring)
   variable_table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   variable_table(i,VT_ORIGININDX)  = trim(filename)
   variable_table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ((variable_table(i,1) == ' ') ) exit MyLoop

   ! Any other condition is an error.
   if ( any(variable_table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:variables not fully specified.'
      string2 = 'Must be 6 entries per variable, last known variable name is'
      string3 = trim(variable_table(i,1))
      call error_handler(E_ERR,'verify_variables',string1, &
          source,revision,revdate,text2=string2,text3=string3)
   endif

   ! Make sure variable exists in netCDF file, as long as it is not vTEC.
   ! vTEC gets constructed.

   if (varname == 'VTEC') include_vTEC_in_state = .true.

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''No obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_variables',string1,source,revision,revdate)
   endif

   ngood = ngood + 1

enddo MyLoop

! Do we need to augment the state vector with the parameter to estimate?
if ( estimate_f10_7 ) then

   string1 = 'Estimating f10_7 is not supported.'
   string2 = 'This feature is under development.'
   string3 = 'If you want to experiment with this, change E_ERR to E_MSG in "verify_variables".'
   call error_handler(E_ERR, 'verify_variables:', string1, &
                      source, revision, revdate, text2=string2, text3=string3)

   ngood = ngood + 1
   variable_table(ngood,VT_VARNAMEINDX) = 'f10_7'
   variable_table(ngood,VT_KINDINDX)    = 'KIND_1D_PARAMETER'
   variable_table(ngood,VT_MINVALINDX)  = 'NA'
   variable_table(ngood,VT_MAXVALINDX)  = 'NA'
   variable_table(ngood,VT_ORIGININDX)  = 'CALCULATE'
   variable_table(ngood,VT_STATEINDX)   = 'UPDATE'

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
   do i = 1,ngood
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

!-------------------------------------------------------------------------------
! Now that we know how many variables, etc., fill metadata structure 'progvar'
! read_TIEGCM_restart() uses this structure to figure out what to put in the
! DART state vector.
!-------------------------------------------------------------------------------

call nc_check(nf90_open(tiegcm_restart_file_name, NF90_NOWRITE, ncid1), &
              'verify_variables','open '//trim(tiegcm_restart_file_name))

call nc_check(nf90_open(tiegcm_secondary_file_name, NF90_NOWRITE, ncid2), &
              'verify_variables','open '//trim(tiegcm_secondary_file_name))

index1  = 1;
indexN  = 0;
FillLoop : do ivar = 1, ngood

   varname = trim(variable_table(ivar,VT_VARNAMEINDX))

   ! Check to see which file contains the variable.
   if     (variable_table(ivar,VT_ORIGININDX) == 'RESTART') then
      ncid     = ncid1
      filename = tiegcm_restart_file_name
   elseif (variable_table(ivar,VT_ORIGININDX) == 'SECONDARY') then
      ncid     = ncid2
      filename = tiegcm_secondary_file_name
   elseif (variable_table(ivar,VT_ORIGININDX) == 'CALCULATE') then
      ncid     = -1
      filename = 'calculate'
   else
      write(string1,'(''unknown option ['',a,''] for variable '',a)') &
                             trim(variable_table(ivar,VT_ORIGININDX)), trim(varname)
      string2 = 'valid options are "RESTART", "SECONDARY", or "CALCULATE"'
      call error_handler(E_ERR,'verify_variables',string1,source,revision,revdate,text2=string2)
   endif

   progvar(ivar)%varname           = varname
   progvar(ivar)%long_name         = 'blank'
   progvar(ivar)%units             = 'blank'
   progvar(ivar)%dimnames          = 'blank'
   progvar(ivar)%dimlens           = -1
   progvar(ivar)%rank              = -1
   progvar(ivar)%varsize           = -1
   progvar(ivar)%index1            = -1
   progvar(ivar)%indexN            = -1
   progvar(ivar)%verticalvar       = 'undefined'
   progvar(ivar)%kind_string       = trim(variable_table(ivar,VT_KINDINDX))
   progvar(ivar)%dart_kind         = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%xtype             = -1
   progvar(ivar)%missingR8         = MISSING_R8
   progvar(ivar)%missingR4         = MISSING_R4
   progvar(ivar)%rangeRestricted   = BOUNDED_NONE
   progvar(ivar)%minvalue          = MISSING_R8
   progvar(ivar)%maxvalue          = MISSING_R8
   progvar(ivar)%has_missing_value = .false.
   progvar(ivar)%update            = .false.
   progvar(ivar)%origin            = trim(filename)

   if ((variable_table(ivar,VT_STATEINDX)  == 'UPDATE') .and. &
       (variable_table(ivar,VT_ORIGININDX) == 'RESTART')) progvar(ivar)%update = .true.

   ! Convert the information in columns 3 and 4 into the appropriate variables
   call SetVariableLimits(ivar)

   if (trim(varname) == 'f10_7') then
      ! parameter to estimate, value comes from tiegcm.nml
      varsize = 1
      progvar(ivar)%long_name         = '10.7 cm daily solar flux'
      progvar(ivar)%units             = 'none'
      progvar(ivar)%dimnames(1)       = 'parameter'
      progvar(ivar)%dimlens(1)        = 0
      progvar(ivar)%rank              = 0
      progvar(ivar)%varsize           = varsize
      progvar(ivar)%index1            = index1
      progvar(ivar)%indexN            = index1 + varsize - 1
      progvar(ivar)%xtype             = NF90_DOUBLE
      index1                          = index1 + varsize ! sets up for next variable

      cycle FillLoop
   endif

   if (trim(varname) == 'VTEC') then
      ! 2D variable not in netCDF file, but useful to be part of state.
      varsize = nlon*nlat*1
      progvar(ivar)%long_name         = 'Total Electron Content'
      progvar(ivar)%units             = 'TECU'
      progvar(ivar)%dimnames(1:3)     = (/'lon ','lat ','time'/)
      progvar(ivar)%dimlens(1:3)      = (/nlon,  nlat,      1 /)
      progvar(ivar)%rank              = 3
      progvar(ivar)%varsize           = varsize
      progvar(ivar)%index1            = index1
      progvar(ivar)%indexN            = index1 + varsize - 1
      progvar(ivar)%xtype             = NF90_DOUBLE
      index1                          = index1 + varsize ! sets up for next variable

      cycle FillLoop
   endif

   if (ncid < 0 ) then
      write(string1,'(a,'' is supposed to be calculated.'')') trim(varname)
      string2 = 'Not supported at this time.'
      call error_handler(E_ERR, 'verify_variables', string1, &
                 source,revision,revdate, text2=string2)
   endif

   if (trim(varname) == 'ZG') ivarZG = ivar

   ! Now that the "special" cases are handled, just read the information
   ! from the input netCDF file and insert into our structure.

   string2 = trim(filename)//' '//trim(varname)

   ncerr = NF90_inq_varid(ncid, trim(varname), VarID)
   if (ncerr /= NF90_NOERR) then
      write(string1,'(''No variable '',a,'' in '',a)') trim(varname), trim(filename)
      string2 = 'If variable name is unexpected, check against input namelist.'
      string3 = 'Make sure each variable has a DART kind, min, and max values.'
      call error_handler(E_ERR, 'verify_variables', string1, &
      source,revision,revdate, text2=string2, text3=string3)
   endif

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%rank, xtype=progvar(ivar)%xtype), &
            'verify_variables', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'verify_variables', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if ( trim(varname) == 'ZG' )  then
      progvar(ivar)%units = 'meters'  ! converted as necessary in read_TIEGCM_secondary
   elseif( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'verify_variables', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Saving any FillValue so I can use it when I read and write ...

   if (progvar(ivar)%xtype == NF90_DOUBLE) then
      if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
         progvar(ivar)%missingR8         = spvalR8
         progvar(ivar)%has_missing_value = .true.
      endif
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR4) == NF90_NOERR) then
         progvar(ivar)%missingR4         = spvalR4
         progvar(ivar)%has_missing_value = .true.
      endif
   else
      ! The missing_value and _FillValue attributes, specifically ...
      write(string1,*)'do not support variables of type ',progvar(ivar)%xtype
      call error_handler(E_ERR,'verify_variables',string1,source,revision,revdate)
   endif

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%rank

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                          'verify_variables', string1)
      ! read_TIEGCM_restart only reads the latest time step, so no matter how many
      ! time steps are defined, the time dimension is only 1
      if (trim(dimname) == 'time') dimlen = 1

      ! record what kind of vertical coordinate system is used for this variable
      if (trim(dimname) ==  'lev') progvar(ivar)%verticalvar =  'lev'
      if (trim(dimname) == 'ilev') progvar(ivar)%verticalvar = 'ilev'

      ! TN,UN,VN,O1,O2 actually have useless 'top' levels.
      if (trim(dimname) ==  'lev') dimlen = dimlen - 1

      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname

      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize = varsize
   progvar(ivar)%index1  = index1
   progvar(ivar)%indexN  = index1 + varsize - 1
   index1                = index1 + varsize      ! sets up for next variable

enddo FillLoop

ReportLoop : do ivar = 1, ngood
if (do_output() .and. (debug > 2)) then
   write(logfileunit,*)
   write(logfileunit,*)trim(progvar(ivar)%varname),' variable number ',ivar
   write(logfileunit,*)'long_name         ',trim(progvar(ivar)%long_name)
   write(logfileunit,*)'units             ',trim(progvar(ivar)%units)
   write(logfileunit,*)'dimnames          ',progvar(ivar)%dimnames(1:progvar(ivar)%rank)
   write(logfileunit,*)'dimlens           ',progvar(ivar)%dimlens( 1:progvar(ivar)%rank)
   write(logfileunit,*)'rank              ',progvar(ivar)%rank
   write(logfileunit,*)'varsize           ',progvar(ivar)%varsize
   write(logfileunit,*)'index1            ',progvar(ivar)%index1
   write(logfileunit,*)'indexN            ',progvar(ivar)%indexN
   write(logfileunit,*)'dart_kind         ',progvar(ivar)%dart_kind
   write(logfileunit,*)'kind_string       ',trim(progvar(ivar)%kind_string)
   write(logfileunit,*)'verticalvar       ',trim(progvar(ivar)%verticalvar)
   write(logfileunit,*)'xtype             ',progvar(ivar)%xtype
   write(logfileunit,*)'rangeRestricted   ',progvar(ivar)%rangeRestricted
   write(logfileunit,*)'minvalue          ',progvar(ivar)%minvalue
   write(logfileunit,*)'maxvalue          ',progvar(ivar)%maxvalue
   write(logfileunit,*)'missingR8         ',progvar(ivar)%missingR8
   write(logfileunit,*)'missingR4         ',progvar(ivar)%missingR4
   write(logfileunit,*)'has_missing_value ',progvar(ivar)%has_missing_value
   write(logfileunit,*)'update            ',progvar(ivar)%update
   write(logfileunit,*)'origin            ',trim(progvar(ivar)%origin)

   write(    *      ,*)
   write(    *      ,*)trim(progvar(ivar)%varname),' variable number ',ivar
   write(    *      ,*)'long_name         ',trim(progvar(ivar)%long_name)
   write(    *      ,*)'units             ',trim(progvar(ivar)%units)
   write(    *      ,*)'dimnames          ',progvar(ivar)%dimnames(1:progvar(ivar)%rank)
   write(    *      ,*)'dimlens           ',progvar(ivar)%dimlens( 1:progvar(ivar)%rank)
   write(    *      ,*)'rank              ',progvar(ivar)%rank
   write(    *      ,*)'varsize           ',progvar(ivar)%varsize
   write(    *      ,*)'index1            ',progvar(ivar)%index1
   write(    *      ,*)'indexN            ',progvar(ivar)%indexN
   write(    *      ,*)'dart_kind         ',progvar(ivar)%dart_kind
   write(    *      ,*)'kind_string       ',trim(progvar(ivar)%kind_string)
   write(    *      ,*)'verticalvar       ',trim(progvar(ivar)%verticalvar)
   write(    *      ,*)'xtype             ',progvar(ivar)%xtype
   write(    *      ,*)'rangeRestricted   ',progvar(ivar)%rangeRestricted
   write(    *      ,*)'minvalue          ',progvar(ivar)%minvalue
   write(    *      ,*)'maxvalue          ',progvar(ivar)%maxvalue
   write(    *      ,*)'missingR8         ',progvar(ivar)%missingR8
   write(    *      ,*)'missingR4         ',progvar(ivar)%missingR4
   write(    *      ,*)'has_missing_value ',progvar(ivar)%has_missing_value
   write(    *      ,*)'update            ',progvar(ivar)%update
   write(    *      ,*)'origin            ',trim(progvar(ivar)%origin)
endif
enddo ReportLoop

call nc_check(nf90_close(ncid1), 'verify_variables', &
        'close '//trim(tiegcm_restart_file_name))
call nc_check(nf90_close(ncid2), 'verify_variables', &
        'close '//trim(tiegcm_secondary_file_name))

end subroutine verify_variables


!-------------------------------------------------------------------------------


subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)
! I am trying to preserve the original shape of the variable as much as possible.
!
! the netCDF declarations look like : variable(time,       level, lat, lon) becomes
!                                     variable(time, copy, level, lat, lon)
!
! Since 'time' or 'Time' is the unlimited dimension in both ... I can skip it
! in the DEFDIM loop.

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i, mydimid

ndims = 0

DEFDIM : do i = 1,progvar(ivar)%rank

   if ((trim(progvar(ivar)%dimnames(i)) == 'Time') .or. &
       (trim(progvar(ivar)%dimnames(i)) == 'time')) cycle DEFDIM

   call nc_check(nf90_inq_dimid(ncid=ncid,name=progvar(ivar)%dimnames(i),dimid=mydimid), &
          'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid

enddo DEFDIM

ndims = ndims + 1               ! The next-to-last dimension is 'copy'
dimids(ndims) = memberdimid
ndims = ndims + 1               ! The last dimension is unlimited == time
dimids(ndims) = unlimitedDimid

if (do_output() .and. (debug > 99)) then
   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%rank)
   write(logfileunit,*)' thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%rank)
   write(     *     ,*)' thus dimids ',dimids(1:ndims)
endif

return
end subroutine define_var_dims


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


subroutine vert_interp(x, lon_index, lat_index, height, ikind, vertstagger, ivar, val, istatus)
! returns the value at an arbitrary height on an existing horizontal grid location.
! istatus == 0 is a 'good' return code.

real(r8),         intent(in)  :: x(:)
integer,          intent(in)  :: lon_index
integer,          intent(in)  :: lat_index
real(r8),         intent(in)  :: height
integer,          intent(in)  :: ikind
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

if ( vertstagger == 'ilev') then

   zgrid_bottom = x(get_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=1    ))
   zgrid_top    = x(get_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=nilev))

   ! cannot extrapolate below bottom or beyond top ... fail ...
   if ((zgrid_bottom > height) .or. (zgrid_top < height)) return

   ! Figure out what level is above/below, and by how much
   h_loop_interface : do k = 2, nilev

      zgrid = x(get_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=k))

      if (height <= zgrid) then
         lev_top    = k
         lev_bottom = lev_top - 1
         delta_z    = zgrid - x(get_index(ivarZG,indx1=lon_index,indx2=lat_index,indx3=lev_bottom))
         frac_lev   = (zgrid - height)/delta_z
         exit h_loop_interface
      endif

   enddo h_loop_interface

elseif ( vertstagger == 'lev') then
   ! Variable is on level midpoints, not ilevels.
   ! Get height as the average of the ilevels.

   ! ilev index    1      2      3      4    ...  27    28    29
   ! ilev value  -7.00, -6.50, -6.00, -5.50, ... 6.00, 6.50, 7.00 ;
   !  lev value     -6.75, -6.25, -5.75, -5.25, ... 6.25, 6.75
   !  lev index        1      2      3      4    ...  27    28

   !mid_level 1
   zgrid_bottom = (x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=1)) + &
                   x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=2))) / 2.0_r8

   !mid_level nlev
   zgrid_top    = (x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=nilev-1)) + &
                   x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=nilev))) / 2.0_r8

   ! cannot extrapolate below bottom or beyond top ... fail ...
   if ((zgrid_bottom > height) .or. (zgrid_top < height)) return

   ! Figure out what level is above/below, and by how much
   h_loop_midpoint: do k = 2, nilev-1

     lev_bottom = k-1
     lev_top    = k

     zgrid_lower = &
             (x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k-1 )) + &
              x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k   ))) / 2.0_r8

     zgrid_upper = &
             (x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k  )) + &
              x(get_index(ivarZG, indx1=lon_index, indx2=lat_index, indx3=k+1))) / 2.0_r8

     if (height  <= zgrid_upper) then
        if (zgrid_upper == zgrid_lower) then ! avoid divide by zero
           frac_lev = 0.0_r8  ! the fraction does not matter ...
        else
           delta_z  = zgrid_upper - zgrid_lower
           frac_lev = (zgrid_upper - height)/delta_z
        endif
        exit h_loop_midpoint
     endif

   enddo h_loop_midpoint

else
   write(string1,*)'Unknown vertical stagger ',trim(vertstagger)
   call error_handler(E_MSG,'vert_interp:', string1, source, revision, revdate )
endif

! Check to make sure we didn't fall through the h_loop ... unlikely (impossible?)
if ( (frac_lev == MISSING_R8) .or. (lev_top == 0) .or. (lev_bottom == 0) ) then
   write(string1,*)'Should not be here ... fell through loop.'
   call error_handler(E_MSG,'vert_interp:', string1, source, revision, revdate )
   return
endif

istatus = 0 ! If we made it this far, it worked.

if (ikind == KIND_PRESSURE) then ! log-linear interpolation in height

   val_top    = plevs(lev_top)     !pressure at midpoint [Pa]
   val_bottom = plevs(lev_bottom)  !pressure at midpoint [Pa]
   val        = exp(frac_lev * log(val_bottom) + (1.0 - frac_lev) * log(val_top))

else ! simple linear interpolation in height

   val_top    = x(get_index(ivar, indx1=lon_index, indx2=lat_index, indx3=lev_top))
   val_bottom = x(get_index(ivar, indx1=lon_index, indx2=lat_index, indx3=lev_bottom))
   val        =     frac_lev *     val_bottom  + (1.0 - frac_lev) *     val_top
endif

end subroutine vert_interp


!-------------------------------------------------------------------------------


function get_height(ivar, lonindex, latindex, levindex)
! TIEGCM's 'natural' vertical coordinate is pressure, DART needs it in height.
!
! Need the ensemble mean value of ZG for this location.
!
! ZG exists on ilev coordinates ... "interface levels"
! Variables that have the same vertical coordinate system are easy.
! All the variables on "midpoint levels" must be computed.
!
! ilev:long_name = "interface levels" ;
!  lev:long_name = "midpoint levels" ;
!
! For example ... IF:
! ilev index    1      2      3      4    ...  27    28    29
! ilev value  -7.00, -6.50, -6.00, -5.50, ... 6.00, 6.50, 7.00 ;
!  lev value      -6.75, -6.25, -5.75, -5.25, ... 6.25, 6.75;
!  lev index        1      2      3      4    ...  27    28

integer, intent(in) :: ivar
integer, intent(in) :: lonindex
integer, intent(in) :: latindex
integer, intent(in) :: levindex
real(r8)            :: get_height

integer  ::  index1,  index2

if (trim(progvar(ivar)%verticalvar) == 'ilev') then

   if (levindex > nilev) then
      write(string1,*)'requesting out-of-bounds level [',levindex,']'
      write(string2,*)'for variable ',trim(progvar(ivar)%varname)
      call error_handler(E_ERR,'get_height', string1, &
                         source, revision, revdate, text2=string2)
   endif

     index1 = get_index(ivarZG, indx1=lonindex, indx2=latindex, indx3=levindex)
     get_height = ens_mean(index1)

elseif (trim(progvar(ivar)%verticalvar) == 'lev') then

   if (levindex > nlev) then
      write(string1,*)'requesting out-of-bounds level [',levindex,']'
      write(string2,*)'for variable ',trim(progvar(ivar)%varname)
      call error_handler(E_ERR,'get_height', string1, &
                         source, revision, revdate, text2=string2)
   endif

   ! Since ZG is defined for level interfaces, requests for heights on 
   ! midpoints must be calculated.
   !
   ! incoming index  1 should be an average of ZG ilev 1+2
   ! incoming index  2 should be an average of ZG ilev 2+3
   ! ...
   ! incoming index 28 should be an average of ZG ilev 28+29
   ! incoming index 29 is not possible

   index1     = get_index(ivarZG, indx1=lonindex, indx2=latindex, indx3=levindex   )
   index2     = get_index(ivarZG, indx1=lonindex, indx2=latindex, indx3=levindex+1 )
   get_height = (ens_mean(index1) + ens_mean(index2)) / 2.0_r8

else
   write(string1,*)'unknown vertical coordinate system <', &
                                  trim(progvar(ivar)%verticalvar),'>'
   write(string2,*)'on variable ',trim(progvar(ivar)%varname)
   call error_handler(E_ERR,'get_height', string1, &
                      source, revision, revdate, text2=string2)
endif

end function get_height


!-------------------------------------------------------------------------------


function get_index(ivar, indx1, indx2, indx3, indx4, knownindex)
!
! returns the statevector index for a lon,lat[,lev[,time]]
!
! This module preserves the number of dimension in the TIEGCM netCDF file.
! NE(time, ilev, lat, lon) becomes a fortran variable
! NE(lon, lat, ilev,    1)
!
! Because DART preserves "time"  as a singleton dimension.
! There is no point requesting the "time" index.
! Put another way, you can call this routine with two indices
! for a rank 3 variable, as long as the last dimension is 1

integer,           intent(in) :: ivar
integer, optional, intent(in) :: indx1
integer, optional, intent(in) :: indx2
integer, optional, intent(in) :: indx3
integer, optional, intent(in) :: indx4
integer, optional, intent(in) :: knownindex
integer                       :: get_index

integer :: ndim1, ndim2, ndim3, ndim4
integer :: ix1, ix2, ix3, ix4

real(r8) :: height

! Set default values, then use input value

ix1 = 1
ix2 = 1
ix3 = 1
ix4 = 1

if ( present(indx1) ) ix1 = indx1
if ( present(indx2) ) ix2 = indx2
if ( present(indx3) ) ix3 = indx3
if ( present(indx4) ) ix4 = indx4

if     (progvar(ivar)%rank == 0) then ! scalars

   get_index = progvar(ivar)%index1

   if ( present(indx1) .or. present(indx2) .or. present(indx3) .or. present(indx4) ) then
      write(string1,*) trim(progvar(ivar)%varname),'called with extra optional arguments.'
      write(string2,*) ' dim1 ', present(indx1), ' dim2 ', present(indx2), &
                       ' dim3 ', present(indx3), ' dim4 ', present(indx4)
      call error_handler(E_ERR, 'get_index', string1, &
            source, revision, revdate, text2=string2)
   endif

elseif (progvar(ivar)%rank == 1) then

   get_index = progvar(ivar)%index1 - 1 + ix1

   if ( present(indx2) .or. present(indx3) .or. present(indx4) ) then
      write(string1,*) trim(progvar(ivar)%varname),'called with extra optional arguments.'
      write(string2,*) ' dim2 ', present(indx2), ' dim3 ', present(indx3), &
                       ' dim4 ', present(indx4)
      call error_handler(E_ERR, 'get_index', string1, &
            source, revision, revdate, text2=string2)
   endif

elseif (progvar(ivar)%rank == 2) then

   ndim1     = progvar(ivar)%dimlens(1)
   get_index = progvar(ivar)%index1 - 1 + ix1 + &
                                         (ix2-1) *  ndim1

   if ( present(indx3) .or. present(indx4) ) then
      write(string1,*) trim(progvar(ivar)%varname),'called with extra optional arguments.'
      write(string2,*) 'dim3 ',present(indx3), ' dim4 ', present(indx4)
      call error_handler(E_ERR, 'get_index', string1, &
            source, revision, revdate, text2=string2)
   endif

elseif (progvar(ivar)%rank == 3) then

   ndim1     = progvar(ivar)%dimlens(1)
   ndim2     = progvar(ivar)%dimlens(2)
   get_index = progvar(ivar)%index1 - 1 + ix1 + &
                                         (ix2-1) *  ndim1 + &
                                         (ix3-1) * (ndim1*ndim2)

   if ( present(indx4) ) then
      write(string1,*) trim(progvar(ivar)%varname),'called with an extra optional argument.'
      write(string2,*) 'dim4 ', present(indx4)
      call error_handler(E_ERR, 'get_index', string1, &
            source, revision, revdate, text2=string2)
   endif

elseif (progvar(ivar)%rank == 4) then ! Fortran shape is lon.lat.lev.time

   ndim1     = progvar(ivar)%dimlens(1)
   ndim2     = progvar(ivar)%dimlens(2)
   ndim3     = progvar(ivar)%dimlens(3)
   get_index = progvar(ivar)%index1 - 1 + ix1 + &
                                         (ix2-1) *  ndim1 + &
                                         (ix3-1) * (ndim1*ndim2) + &
                                         (ix4-1) * (ndim1*ndim2*ndim3)

else
   write(string1,*) trim(progvar(ivar)%varname)//' has unsupported shape.'
   call error_handler(E_ERR, 'get_index', string1, source, revision, revdate )
endif

! This is useful when testing
if (present(knownindex)) then
   if (knownindex /= get_index) then

      if (present(indx1) .and. present(indx2) .and. present(indx3) .and. present(indx4)) then
         height = get_height(ivar, indx1, indx2, indx3)
         write(*,*)'index_in, loni, lati, levi, timei', &
             knownindex, indx1, indx2, indx3, indx4, lons(indx1), lats(indx2), height
      else if (present(indx1) .and. present(indx2) .and. present(indx3)) then
         height = get_height(ivar, indx1, indx2, indx3)
         write(*,*)'index_in, loni, lati, levi', &
             knownindex, indx1, indx2, indx3, lons(indx1), lats(indx2), height
      else if (present(indx1) .and. present(indx2)) then
         write(*,*)'index_in, loni, lati', &
             knownindex, indx1, indx2, lons(indx1), lats(indx2)
      else if (present(indx1)) then
         write(*,*)'index_in, loni', &
             knownindex, indx1, lons(indx1)
      else 
         write(*,*)'index_in', &
             knownindex
      endif

      write(string1,*)'FAILURE ... original index /= get_index() value'
      write(string2,*)'original index = ',knownindex
      write(string3,*)'get_index()    = ',get_index
      call error_handler(E_ERR,'get_index', string1, &
              source, revision, revdate, text2=string2, text3=string3)
   endif
endif

end function get_index


!-------------------------------------------------------------------------------


subroutine vector_to_prog_var(statevec, ivar, ncid, VarID, numdims, &
                              mystart, mycount, limit)
! Unpacks a DART state vector into appropriate variable shape and then writes it
! to the netCDF file. The shape of the variable is specified by the 'mycount' array.
! If 'limit' has a value of TRUE, the min/max values of the variable are
! restricted to those values specified in the 'progvar' array.
! Typically, the DART files (True_State.nc, Prior_Diag.nc, Posterior_Diag.nc)
! are not restricted, but the files for TIEGCM may be.

real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: ivar
integer,                intent(in) :: ncid
integer,                intent(in) :: VarID
integer,                intent(in) :: numdims
integer,  dimension(:), intent(in) :: mystart
integer,  dimension(:), intent(inout) :: mycount
logical,                intent(in) :: limit

real(r8)                                    :: data_0d
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:)   :: data_4d_array
real(r8), allocatable, dimension(:,:,:,:,:) :: data_5d_array

integer :: icount, idim1, idim2, idim3, idim4, idim5

if ( progvar(ivar)%rank == 0 ) then  ! handle scalars

   data_0d = statevec(progvar(ivar)%index1)

   ! We only want to apply the limits when converting for TIEGCM
   ! In general, we do not apply limits for DART diagnostics.
   if ( limit ) then

      if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         if ((data_0d /= MISSING_R8) .and. &
             (data_0d > progvar(ivar)%maxvalue)) &
              data_0d = progvar(ivar)%maxvalue
      endif

      if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         if ((data_0d /= MISSING_R8) .and. &
             (data_0d < progvar(ivar)%minvalue)) &
              data_0d = progvar(ivar)%minvalue
      endif

      if ( progvar(ivar)%has_missing_value ) then
         if (data_0d == MISSING_R8) data_0d = progvar(ivar)%missingR8
      endif

   endif

   call nc_check(nf90_put_var(ncid, VarID, (/ data_0d /), &
           start = mystart(1:numdims), count=mycount(1:numdims)), &
           'vector_to_prog_var', 'put_var '//trim(progvar(ivar)%varname))

elseif ( numdims == 1 ) then

   allocate(data_1d_array( mycount(1) ))

   icount = progvar(ivar)%index1

   do idim1 = 1, mycount(1)
      data_1d_array(idim1) = statevec(icount)
      icount = icount + 1
   enddo

   if ( (icount-1) /= progvar(ivar)%indexN ) then
      write(string1, *) trim(progvar(ivar)%varname), 'supposed to end at ', &
                             progvar(ivar)%indexN
      write(string2, *) 'it ended at ',icount-1
      call error_handler(E_ERR, 'vector_to_prog_var', string1, &
                source, revision, revdate, text2=string2)
   endif

   if ( limit ) then

      if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_1d_array /= MISSING_R8) .and. &
                (data_1d_array > progvar(ivar)%maxvalue)) &
                 data_1d_array = progvar(ivar)%maxvalue
      endif

      if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_1d_array /= MISSING_R8) .and. &
                (data_1d_array < progvar(ivar)%minvalue)) &
                 data_1d_array = progvar(ivar)%minvalue
      endif

      if ( progvar(ivar)%has_missing_value ) then
         where (data_1d_array == MISSING_R8) data_1d_array = progvar(ivar)%missingR8
      endif

   endif

   call nc_check(nf90_put_var(ncid, VarID, data_1d_array, &
           start = mystart(1:numdims), count=mycount(1:numdims)), &
           'vector_to_prog_var', 'put_var '//trim(progvar(ivar)%varname))
   deallocate(data_1d_array)

elseif ( numdims == 2 ) then

   allocate(data_2d_array(mycount(1), mycount(2)))

   icount = progvar(ivar)%index1

   do idim2 = 1, mycount(2)
   do idim1 = 1, mycount(1)
      data_2d_array(idim1,idim2) = statevec(icount)
      icount = icount + 1
   enddo
   enddo

   if ( (icount-1) /= progvar(ivar)%indexN ) then
      write(string1, *) trim(progvar(ivar)%varname), 'supposed to end at ', &
                             progvar(ivar)%indexN
      write(string2, *) 'it ended at ',icount-1
      call error_handler(E_ERR, 'vector_to_prog_var', string1, &
                source, revision, revdate, text2=string2)
   endif

   if ( limit ) then

      if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_2d_array /= MISSING_R8) .and. &
                (data_2d_array > progvar(ivar)%maxvalue)) &
                 data_2d_array = progvar(ivar)%maxvalue
      endif

      if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_2d_array /= MISSING_R8) .and. &
                (data_2d_array < progvar(ivar)%minvalue)) &
                 data_2d_array = progvar(ivar)%minvalue
      endif

      if ( progvar(ivar)%has_missing_value ) then
         where (data_2d_array == MISSING_R8) data_2d_array = progvar(ivar)%missingR8
      endif
   endif

   call nc_check(nf90_put_var(ncid, VarID, data_2d_array, &
           start = mystart(1:numdims), count=mycount(1:numdims)), &
           'vector_to_prog_var', 'put_var '//trim(progvar(ivar)%varname))
   deallocate(data_2d_array)

elseif ( numdims == 3 ) then

   allocate(data_3d_array( mycount(1), mycount(2), mycount(3) ))

   icount = progvar(ivar)%index1

   do idim3 = 1, mycount(3)
   do idim2 = 1, mycount(2)
   do idim1 = 1, mycount(1)
      data_3d_array(idim1,idim2,idim3) = statevec(icount)
      icount = icount + 1
   enddo
   enddo
   enddo

   if ( (icount-1) /= progvar(ivar)%indexN ) then
      write(string1, *) trim(progvar(ivar)%varname), 'supposed to end at ', &
                             progvar(ivar)%indexN
      write(string2, *) 'it ended at ',icount-1
      call error_handler(E_ERR, 'vector_to_prog_var', string1, &
                source, revision, revdate, text2=string2)
   endif

   if ( limit ) then

      if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_3d_array /= MISSING_R8) .and. &
                (data_3d_array > progvar(ivar)%maxvalue)) &
                 data_3d_array = progvar(ivar)%maxvalue
      endif

      if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_3d_array /= MISSING_R8) .and. &
                (data_3d_array < progvar(ivar)%minvalue)) &
                 data_3d_array = progvar(ivar)%minvalue
      endif

      if ( progvar(ivar)%has_missing_value ) then
         where (data_3d_array == MISSING_R8) data_3d_array = progvar(ivar)%missingR8
      endif

   endif

   if (do_output() .and. (debug > 3)) then
      write(string2,*)trim(progvar(ivar)%varname)//' start is ',mystart(1:numdims)
      write(string3,*)trim(progvar(ivar)%varname)//' count is ',mycount(1:numdims)
      call error_handler(E_MSG,'vector_to_prog_var',' ',text2=string2,text3=string3)
   endif

   call nc_check(nf90_put_var(ncid, VarID, data_3d_array, &
           start = mystart(1:numdims), count=mycount(1:numdims)), &
           'vector_to_prog_var', 'put_var '//trim(progvar(ivar)%varname))
   deallocate(data_3d_array)

elseif ( numdims == 4 ) then

   allocate(data_4d_array( mycount(1), mycount(2), mycount(3), mycount(4) ))

   icount = progvar(ivar)%index1

   do idim4 = 1, mycount(4)
   do idim3 = 1, mycount(3)
   do idim2 = 1, mycount(2)
   do idim1 = 1, mycount(1)
      data_4d_array(idim1,idim2,idim3,idim4) = statevec(icount)
      icount = icount + 1
   enddo
   enddo
   enddo
   enddo

   if ( (icount-1) /= progvar(ivar)%indexN ) then
      write(string1, *) trim(progvar(ivar)%varname), 'supposed to end at ', &
                             progvar(ivar)%indexN
      write(string2, *) 'it ended at ',icount-1
      call error_handler(E_ERR, 'vector_to_prog_var', string1, &
                source, revision, revdate, text2=string2)
   endif

   if ( limit ) then

      if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_4d_array /= MISSING_R8) .and. &
                (data_4d_array > progvar(ivar)%maxvalue)) &
                 data_4d_array = progvar(ivar)%maxvalue
      endif

      if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_4d_array /= MISSING_R8) .and. &
                (data_4d_array < progvar(ivar)%minvalue)) &
                 data_4d_array = progvar(ivar)%minvalue
      endif

      if ( progvar(ivar)%has_missing_value ) then
         where (data_4d_array == MISSING_R8) data_4d_array = progvar(ivar)%missingR8
      endif

   endif

   if (do_output() .and. (debug > 3)) then
      write(string2,*)trim(progvar(ivar)%varname)//' start is ',mystart(1:numdims)
      write(string3,*)trim(progvar(ivar)%varname)//' count is ',mycount(1:numdims)
      call error_handler(E_MSG,'vector_to_prog_var',' ',text2=string2,text3=string3)
   endif

   call nc_check(nf90_put_var(ncid, VarID, data_4d_array(:,:,1:mycount(3),:), &
           start = mystart(1:numdims), count=mycount(1:numdims)), &
           'vector_to_prog_var', 'put_var '//trim(progvar(ivar)%varname))

   deallocate(data_4d_array)

elseif ( numdims == 5 ) then

   allocate(data_5d_array( mycount(1), mycount(2), mycount(3), mycount(4), mycount(5) ))

   icount = progvar(ivar)%index1

   do idim5 = 1, mycount(5)
   do idim4 = 1, mycount(4)
   do idim3 = 1, mycount(3)
   do idim2 = 1, mycount(2)
   do idim1 = 1, mycount(1)
      data_5d_array(idim1,idim2,idim3,idim4,idim5) = statevec(icount)
      icount = icount + 1
   enddo
   enddo
   enddo
   enddo
   enddo

   if ( (icount-1) /= progvar(ivar)%indexN ) then
      write(string1, *) trim(progvar(ivar)%varname), 'supposed to end at ', &
                             progvar(ivar)%indexN
      write(string2, *) 'it ended at ',icount-1
      call error_handler(E_ERR, 'vector_to_prog_var', string1, &
                source, revision, revdate, text2=string2)
   endif

   if ( limit ) then

      if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_5d_array /= MISSING_R8) .and. &
                (data_5d_array > progvar(ivar)%maxvalue)) &
                 data_5d_array = progvar(ivar)%maxvalue
      endif

      if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
          (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
         where ((data_5d_array /= MISSING_R8) .and. &
                (data_5d_array < progvar(ivar)%minvalue)) &
                 data_5d_array = progvar(ivar)%minvalue
      endif

      if ( progvar(ivar)%has_missing_value ) then
         where (data_5d_array == MISSING_R8) data_5d_array = progvar(ivar)%missingR8
      endif

   endif

   if (do_output() .and. (debug > 3)) then
      write(string2,*)trim(progvar(ivar)%varname)//' start is ',mystart(1:numdims)
      write(string3,*)trim(progvar(ivar)%varname)//' count is ',mycount(1:numdims)
      call error_handler(E_MSG,'vector_to_prog_var',' ',text2=string2,text3=string3)
   endif

   call nc_check(nf90_put_var(ncid, VarID, data_5d_array, &
           start = mystart(1:numdims), count=mycount(1:numdims)), &
           'vector_to_prog_var', 'put_var '//trim(progvar(ivar)%varname))
   deallocate(data_5d_array)

else

   write(string1,*)'cannot handle tiegcm variables with ',numdims,' dimensions.'
   write(string2,*)trim(progvar(ivar)%varname),' has rank ', &
                        progvar(ivar)%dimnames(1:progvar(ivar)%rank)
   call error_handler(E_ERR, 'vector_to_prog_var', string1, &
             source, revision, revdate, text2=string2)

endif

end subroutine vector_to_prog_var

!===============================================================================
! These routines are overloaded to prog_var_to_vector()
!===============================================================================

subroutine var1d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:)
real(r8), intent(out) ::   x(:)

integer :: icount, idim1

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:1D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim1 = 1,size(var,1)
!  write(8,*)'var ',ifield, icount, idim1
   x(icount) = var(idim1)
   icount = icount + 1
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:1D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var1d_to_vector

!-------------------------------------------------------------------------------

subroutine var2d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:,:)
real(r8), intent(out) :: x(:)

integer :: icount,idim1,idim2

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:2D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
!  write(8,*)'var ',ifield, icount, idim1, idim2
   x(icount) = var(idim1,idim2)
   icount = icount + 1
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:2D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var2d_to_vector

!-------------------------------------------------------------------------------

subroutine var3d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:,:,:)
real(r8), intent(out) :: x(:)

integer :: icount,idim1,idim2,idim3

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:3D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim3 = 1,size(var,3)
do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
!  write(8,*)'var ',ifield, icount, idim1, idim2, idim3
   x(icount) = var(idim1,idim2,idim3)
   icount = icount + 1
enddo
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:3D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var3d_to_vector

!-------------------------------------------------------------------------------

subroutine var4d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:,:,:,:)
real(r8), intent(out) :: x(:)

integer :: icount,idim1,idim2,idim3,idim4

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:4D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim4 = 1,size(var,4)
do idim3 = 1,size(var,3)
do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
   x(icount) = var(idim1,idim2,idim3,idim4)
   icount = icount + 1
enddo
enddo
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:4D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var4d_to_vector

!-------------------------------------------------------------------------------


function FindVar_by_kind(ikind)
! Finds the first variable of the appropriate DART KIND
!
! FIXME There is some confusion about using the T-minus-1 variables
! in this construct. Both TN and TN_NM have the same dart_kind,
! so we use the first one ... but it is not guaranteed that TN
! must preceed TN_NM, for example.

integer, intent(in) :: ikind
integer             :: FindVar_by_kind

integer :: ivar 

FindVar_by_kind = -1

VARLOOP : do ivar = 1,nfields
   if (progvar(ivar)%dart_kind == ikind) then
      FindVar_by_kind = ivar
      return
   endif
enddo VARLOOP

if (do_output() .and. (debug > 99)) then
   write(string1, *) 'unable to find a variable with a DART kind of ',ikind
   call error_handler(E_MSG, 'FindVar_by_kind', string1, source, revision, revdate )
endif

end function FindVar_by_kind


!-------------------------------------------------------------------------------


function Find_Variable_by_index(myindx, msgstring)
! Given an index into the DART state vector, return the index of metadata
! variable 'progvar' responsible for this portion of the state vector
integer,          intent(in) :: myindx
character(len=*), intent(in) :: msgstring
integer                      :: Find_Variable_by_index

integer :: ivar

Find_Variable_by_index = -1

FindIndex : do ivar = 1,nfields
   if ((myindx >= progvar(ivar)%index1)  .and. &
       (myindx <= progvar(ivar)%indexN)) then
      Find_Variable_by_index = ivar
      exit FindIndex
   endif
enddo FindIndex

if (Find_Variable_by_index < 0) then
   write(string1,*)'index ',myindx,' is out of range of all variables.'
   write(string2,*)'model size is ',model_size
   call error_handler(E_ERR, 'Find_Variable_by_index'//trim(msgstring), string1, &
                      source, revision, revdate, text2=string2 )
endif

end function Find_Variable_by_index


!-------------------------------------------------------------------------------


subroutine SanityCheck(filename, ncid, LonDimID, LatDimID, LevDimID, iLevDimId, &
                       TimeDimID, ntimes)
! Just make sure the dimensions of the netCDF file match those
! that were used for static_init_model()

character(len=*),  intent(in)  :: filename
integer,           intent(in)  :: ncid
integer, optional, intent(out) :: LonDimID
integer, optional, intent(out) :: LatDimID
integer, optional, intent(out) :: LevDimID
integer, optional, intent(out) :: iLevDimID
integer, optional, intent(out) :: TimeDimID
integer, optional, intent(out) :: ntimes

integer :: dimlen, DimID, DimID2

call nc_check(nf90_inq_dimid(ncid, 'lon', DimID), 'SanityCheck', 'inq_dimid lon')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
   'SanityCheck', 'inquire_dimension lon')
if (dimlen .ne. nlon) then
   write(string1, *) trim(filename), ' dim_lon = ',dimlen, ' DART expects ',nlon
   call error_handler(E_ERR,'SanityCheck',string1,source,revision,revdate)
endif
if (present(LonDimID)) LonDimID = dimlen


call nc_check(nf90_inq_dimid(ncid, 'lat', DimID), 'SanityCheck', 'inq_dimid lat')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
   'SanityCheck', 'inquire_dimension lat')
if (dimlen .ne. nlat) then
   write(string1, *) trim(filename), ' dim_lat = ',dimlen, ' DART expects ',nlat
   call error_handler(E_ERR,'SanityCheck',string1,source,revision,revdate)
endif
if (present(LatDimID)) LatDimID = dimlen

! TIEGCM has a useless top level for everything defined on 'lev'
! 'nlev' had a local value (i.e. in DART) that is 1 less than 'nilev'

call nc_check(nf90_inq_dimid(ncid, 'lev', DimID), 'SanityCheck', 'inq_dimid lev')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
   'SanityCheck', 'inquire_dimension lev')
if (dimlen .ne. (nlev+1)) then
   write(string1, *) trim(filename), ' dim_lev = ',dimlen, ' DART expects ',nlev
   call error_handler(E_ERR,'SanityCheck',string1,source,revision,revdate)
endif
if (present(LevDimID)) LevDimID = dimlen


call nc_check(nf90_inq_dimid(ncid, 'ilev', DimID), 'SanityCheck', 'inq_dimid ilev')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
   'SanityCheck', 'inquire_dimension ilev')
if (dimlen .ne. nilev) then
   write(string1, *) trim(filename), ' dim_ilev = ',dimlen, ' DART expects ',nilev
   call error_handler(E_ERR,'SanityCheck',string1,source,revision,revdate)
endif
if (present(iLevDimID)) iLevDimID = dimlen


call nc_check(nf90_inq_dimid(ncid, 'time', DimID), 'SanityCheck', 'inq_dimid time')
call nc_check(nf90_inquire(ncid, unlimitedDimId = DimID2), &
   'SanityCheck', 'inquire id of unlimited dimension')
if (DimID .ne. DimID2) then
   write(string1, *) trim(filename), 'has an unlimited dimension ID of ',DimID2
   write(string2, *) 'and a time dimension ID of ',DimID, &
                     '. DART requires them to be the same.'
   call error_handler(E_ERR, 'SanityCheck', string1, &
              source, revision, revdate, text2=string2)
endif
if (present(TimeDimID)) TimeDimID = DimID


call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
   'SanityCheck', 'inquire_dimension time')
if (present(ntimes)) ntimes = dimlen

end subroutine SanityCheck



function get_state_time(filename, ncid, lasttime)
! Gets the latest time in the netCDF file.
character(len=*),  intent(in) :: filename
integer,           intent(in) :: ncid
integer, optional, intent(in) :: lasttime
type(time_type) :: get_state_time

integer :: TimeDimID, time_dimlen, DimID, dimlen, VarID

integer,  parameter         :: nmtime = 3
integer,  dimension(nmtime) :: mtime  ! day, hour, minute
integer                     :: year, doy, utsec
integer, allocatable, dimension(:,:) :: mtimetmp
integer, allocatable, dimension(:)   :: yeartmp

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
        'get_state_time', 'inquire id of time')
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
        'get_state_time', 'inquire_dimension time')
call nc_check(nf90_inq_dimid(ncid, 'mtimedim', DimID), &
        'get_state_time', 'inq_dimid mtimedim')
call nc_check(nf90_inquire_dimension(ncid,     DimID, len=dimlen), &
        'get_state_time', 'inquire_dimension mtimedim')

if (present(lasttime)) then
   if (lasttime /= time_dimlen) then
      write(string1, *) trim(filename), ' last time index is ', &
                        time_dimlen, ' desired ', lasttime
      call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
   endif
endif

if (dimlen /= nmtime) then
   write(string1, *) trim(filename), ' mtimedim = ',dimlen, ' DART expects ', nmtime
   call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
endif

allocate(mtimetmp(dimlen, time_dimlen), yeartmp(time_dimlen))

!... get mtime
call nc_check(nf90_inq_varid(ncid, 'mtime', VarID), &
        'get_state_time', 'inquire id of time')
call nc_check(nf90_get_var(ncid, VarID, values=mtimetmp), &
        'get_state_time', 'get_var mtime')

!... get year
call nc_check(nf90_inq_varid(ncid, 'year', VarID), 'get_state_time', 'inq_varid year')
call nc_check(nf90_get_var(ncid, VarID, values=yeartmp), 'get_state_time', 'get_var year')

! pick off the latest
mtime = mtimetmp(:,time_dimlen)
year  = yeartmp(   time_dimlen)

deallocate(mtimetmp,yeartmp)

doy   =  mtime(1)
utsec = (mtime(2)*60 + mtime(3))*60
get_state_time = set_time(utsec, doy-1) + set_date(year, 1, 1)  ! Jan 1 of whatever year.

if (do_output()) then
   write(*,*) trim(filename)//':get_state_time: tiegcm [year, doy, hour, minute]', &
            year, mtime
   call print_date(get_state_time, str=trim(filename)//':get_state_time: date ')
   call print_time(get_state_time, str=trim(filename)//':get_state_time: time ')
endif

end function get_state_time




subroutine SetVariableLimits(ivar)

! Convert the information in the variable_table into appropriate
! numerical limits for each variable.
! If the limit does not apply, it is set to MISSING_R8, even if
! it is the maximum that does not apply.

integer, intent(in) :: ivar

integer  :: ios
real(r8) :: minvalue, maxvalue

! set the default values

minvalue = MISSING_R8
maxvalue = MISSING_R8
progvar(ivar)%minvalue = MISSING_R8
progvar(ivar)%maxvalue = MISSING_R8

read(variable_table(ivar,VT_MINVALINDX),*,iostat=ios) minvalue
if (ios == 0) progvar(ivar)%minvalue = minvalue

read(variable_table(ivar,VT_MAXVALINDX),*,iostat=ios) maxvalue
if (ios == 0) progvar(ivar)%maxvalue = maxvalue

! rangeRestricted == BOUNDED_NONE  == 0 ... unlimited range
! rangeRestricted == BOUNDED_BELOW == 1 ... minimum, but no maximum
! rangeRestricted == BOUNDED_ABOVE == 2 ... maximum, but no minimum
! rangeRestricted == BOUNDED_BOTH  == 3 ... minimum and maximum

if (   (progvar(ivar)%minvalue /= MISSING_R8) .and. &
       (progvar(ivar)%maxvalue /= MISSING_R8) ) then
   progvar(ivar)%rangeRestricted = BOUNDED_BOTH

elseif (progvar(ivar)%maxvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_ABOVE

elseif (progvar(ivar)%minvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_BELOW

else
   progvar(ivar)%rangeRestricted = BOUNDED_NONE

endif

if (do_output() .and. (debug > 99)) then
   if (progvar(ivar)%rangeRestricted == BOUNDED_NONE) then
      write(*,*)'TJH ',trim(progvar(ivar)%varname),' has no limits.'
   elseif (progvar(ivar)%rangeRestricted == BOUNDED_BELOW) then
      write(*,*)'TJH ',trim(progvar(ivar)%varname),' has min ',progvar(ivar)%minvalue
   elseif (progvar(ivar)%rangeRestricted == BOUNDED_ABOVE) then
      write(*,*)'TJH ',trim(progvar(ivar)%varname),' has max ',progvar(ivar)%maxvalue
   elseif (progvar(ivar)%rangeRestricted == BOUNDED_BOTH) then
      write(*,*)'TJH ',trim(progvar(ivar)%varname),' has min ',progvar(ivar)%minvalue, &
                                                   ' and max ',progvar(ivar)%maxvalue
   else
      write(string1,*)'Unable to determine if ',trim(progvar(ivar)%varname), &
                      ' is range-restricted.'
      write(string2,*)'minimum value string is ',trim(variable_table(ivar,3))
      write(string3,*)'maximum value string is ',trim(variable_table(ivar,4))
      call error_handler(E_ERR,'SetVariableLimits',string1, &
         source,revision,revdate,text2=string2,text3=string3)
   endif
endif

! Check to make sure min is less than max if both are specified.

if ( progvar(ivar)%rangeRestricted == BOUNDED_BOTH ) then
   if (maxvalue < minvalue) then
      write(string1,*)'&model_nml state_variable input error for ',trim(progvar(ivar)%varname)
      write(string2,*)'minimum value (',minvalue,') must be less than '
      write(string3,*)'maximum value (',maxvalue,')'
      call error_handler(E_ERR,'SetVariableLimits',string1, &
         source,revision,revdate,text2=string2,text3=string3)
   endif
endif

end subroutine SetVariableLimits



subroutine apply_attributes_1D(ivar, ncid, VarID, temp1D)
! The missing value attributes have been determined by verify_variables()

integer,  intent(in)    :: ivar
integer,  intent(in)    :: ncid
integer,  intent(in)    :: VarID
real(r8), intent(inout) :: temp1D(:)

integer  :: nfrc
real(r8) :: missr8
real(r4) :: missr4

if (progvar(ivar)%xtype == NF90_FLOAT) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr4)
   if (nfrc == NF90_NOERR) where (temp1D == missr4) temp1D = MISSING_R8
elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr8)
   if (nfrc == NF90_NOERR) where (temp1D == missr8) temp1D = MISSING_R8
endif

end subroutine apply_attributes_1D


subroutine apply_attributes_2D(ivar, ncid, VarID, temp2D)
! The missing value attributes have been determined by verify_variables()

integer,  intent(in)    :: ivar
integer,  intent(in)    :: ncid
integer,  intent(in)    :: VarID
real(r8), intent(inout) :: temp2D(:,:)

integer  :: nfrc
real(r8) :: missr8
real(r4) :: missr4

if (progvar(ivar)%xtype == NF90_FLOAT) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr4)
   if (nfrc == NF90_NOERR) where (temp2D == missr4) temp2D = MISSING_R8
elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr8)
   if (nfrc == NF90_NOERR) where (temp2D == missr8) temp2D = MISSING_R8
endif

end subroutine apply_attributes_2D


subroutine apply_attributes_3D(ivar, ncid, VarID, temp3D)
! The missing value attributes have been determined by verify_variables()

integer,  intent(in)    :: ivar
integer,  intent(in)    :: ncid
integer,  intent(in)    :: VarID
real(r8), intent(inout) :: temp3D(:,:,:)

integer  :: nfrc
real(r8) :: missr8
real(r4) :: missr4

if (progvar(ivar)%xtype == NF90_FLOAT) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr4)
   if (nfrc == NF90_NOERR) where (temp3D == missr4) temp3D = MISSING_R8
elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr8)
   if (nfrc == NF90_NOERR) where (temp3D == missr8) temp3D = MISSING_R8
endif

end subroutine apply_attributes_3D


subroutine apply_attributes_4D(ivar, ncid, VarID, temp4D)
! The missing value attributes have been determined by verify_variables()

integer,  intent(in)    :: ivar
integer,  intent(in)    :: ncid
integer,  intent(in)    :: VarID
real(r8), intent(inout) :: temp4D(:,:,:,:)

integer  :: nfrc
real(r8) :: missr8
real(r4) :: missr4

if (progvar(ivar)%xtype == NF90_FLOAT) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr4)
   if (nfrc == NF90_NOERR) where (temp4D == missr4) temp4D = MISSING_R8
elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
   nfrc = nf90_get_att(ncid, VarID, 'missing_value', missr8)
   if (nfrc == NF90_NOERR) where (temp4D == missr8) temp4D = MISSING_R8
endif

end subroutine apply_attributes_4D


!===============================================================================
! End of model_mod
!===============================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
