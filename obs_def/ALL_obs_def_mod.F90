! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use        types_mod, only : r8, missing_i, missing_r8, RAD2DEG
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, file_exist, &
                             open_file, check_nml_error, logfileunit, close_file
use     location_mod, only : location_type, read_location, write_location, set_location, &
                             interactive_location, set_location_missing
!WRF use     location_mod, only : query_location
use time_manager_mod, only : time_type, read_time, write_time, set_time, set_time_missing, &
                             interactive_time
use  assim_model_mod, only : get_state_meta_data, interpolate


! Need to include use statements for whatever observation kind detail modules are being used
#IFDEF raw_state_1d_integral
   use obs_def_raw_state_mod, only : write_1d_integral, read_1d_integral, &
                                        interactive_1d_integral, get_expected_1d_integral
#ENDIF

implicit none
private

interface assignment(=)
   module procedure copy_obs_def
end interface

public :: init_obs_def, get_obs_def_location, get_obs_kind, get_obs_def_time, &
   get_obs_def_error_variance, set_obs_def_location, set_obs_def_kind, set_obs_def_time, &
   set_obs_def_error_variance, interactive_obs_def, write_obs_def, read_obs_def, &
   obs_def_type, get_expected_obs_from_def, &
   set_radar_obs_def, destroy_obs_def, copy_obs_def, assignment(=)
!WRF public set_obs_def_platform, get_obs_def_platform

! Public for obs kinds
public :: KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, KIND_P, KIND_W, KIND_QR, KIND_TD, &
          KIND_RHO, KIND_VR, KIND_REF, KIND_U10, KIND_V10, KIND_T2, KIND_Q2, KIND_TD2

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type obs_def_type
! In revision, obs_kind module is responsible for taking care of identity obs kinds, too
   private
   type(location_type)   :: location      ! center of mass, so to speak
   integer               :: kind
   type(time_type)       :: time
   real(r8)              :: error_variance
   integer               :: key             ! Used by specialized observation types
end type obs_def_type

logical, save :: module_initialized = .false.




! ADD A LONG TABLE OF DEFINED BUFR INDICES, ETC.
! Definition of observation kind types:

! KIND_U   = zonal wind component
! KIND_V   = meridional wind component
! KIND_PS  = Surface pressure
! KIND_T   = Temperature
! KIND_QV  = Specific humidity (mixing ratio)
! KIND_P   = Pressure
! KIND_W   = Vertical velocity
! KIND_QR  = Rainwater mixing ratio
! KIND_TD  = Dew point temperature
! KIND_RHO = Density
! KIND_VR  = Doppler radar radial velocity
! KIND_REF = Radar reflectivity
! KIND_U10 = zonal windERR component at 10 m AGL
! KIND_V10 = meridional wind component at 10 m AGL
! KIND_T2  = Temperature at 2 m AGL
! KIND_Q2  = Specific humidity (mixing ratio) at 2 m AGL
! KIND_TD2 = Dew point temperature at 2 m AGL

integer, parameter :: KIND_U = 1, KIND_V = 2, KIND_PS = 3, KIND_T = 4,   &
                      KIND_QV = 5, KIND_P = 6, KIND_W = 7, KIND_QR = 8, KIND_TD = 10, &
                      KIND_RHO = 11, &
                      KIND_VR = 100, KIND_REF = 101, &
                      KIND_U10 = 200, KIND_V10 = 201, KIND_T2 = 202, KIND_Q2 = 203, &
                      KIND_TD2 = 204

integer, parameter :: RADIOSONDE_U_WIND_COMPONENT                          = 1, &
                      RADIOSONDE_V_WIND_COMPONENT                          = 2, &
                      SURFACE_PRESSURE                                     = 3, &
                      RADIOSONDE_TEMPERATURE                               = 4, &
                      RADIOSONDE_SPECIFIC_HUMIDITY                         = 5, &
                      RADIOSONDE_PRESSURE                                  = 6, &
                      VERTICAL_VELOCITY                                    = 7, &
                      RAINWATER_MIXING_RATIO                               = 8, &
                      RAW_STATE_VARIABLE                                   = 9, &
                      DEW_POINT_TEMPERATURE                                = 10, &
                      DENSITY                                              = 11, &
                      DOPPLER_RADIAL_VELOCITY                              = 12, &
                      RADAR_REFLECTIVITY                                   = 13, &
                      U_10_METER_WIND                                      = 14, &
                      V_10_METER_WIND                                      = 15, &
                      TEMPERATURE_2_METER                                  = 16, &
                      SPECIFIC_HUMIDITY_2_METER                            = 17, &
                      DEW_POINT_2_METER                                    = 18, &
                      RAW_STATE_1D_INTEGRAL                                = 19



integer, parameter :: max_obs_kinds = 19
integer :: num_kind_assimilate, num_kind_evaluate

type obs_kind_type
   integer              :: index
   character(len = 32) :: name
   logical              :: assimilate
   logical              :: evaluate
end type obs_kind_type

type(obs_kind_type) :: obs_kind_info(max_obs_kinds) = &
   (/obs_kind_type(RADIOSONDE_U_WIND_COMPONENT, 'radiosonde_u_wind_component',   .false., .false.), &
     obs_kind_type(RADIOSONDE_V_WIND_COMPONENT, 'radiosonde_v_wind_component',   .false., .false.), &
     obs_kind_type(SURFACE_PRESSURE,  'surface_pressure',                        .false., .false.), &
     obs_kind_type(RADIOSONDE_TEMPERATURE,  'radiosonde_temperature',            .false., .false.), &
     obs_kind_type(RADIOSONDE_SPECIFIC_HUMIDITY, 'radiosonde_specific_humidity', .false., .false.), &
     obs_kind_type(RADIOSONDE_PRESSURE, 'radiosonde_pressure',                   .false., .false.), &
     obs_kind_type(VERTICAL_VELOCITY,  'vertical_velocity',                      .false., .false.), &
     obs_kind_type(RAINWATER_MIXING_RATIO, 'rainwater_mixing_ratio',             .false., .false.), &
     obs_kind_type(RAW_STATE_VARIABLE,  'raw_state_variable',                    .false., .false.), &
     obs_kind_type(DEW_POINT_TEMPERATURE,  'dew_point_temperature',              .false., .false.), &
     obs_kind_type(DENSITY, 'density',                                           .false., .false.), &
     obs_kind_type(DOPPLER_RADIAL_VELOCITY, 'doppler_radial_velocity',           .false., .false.), &
     obs_kind_type(RADAR_REFLECTIVITY, 'radar_reflectivity',                     .false., .false.), &
     obs_kind_type(U_10_METER_WIND, 'u_10_meter_wind',                           .false., .false.), &
     obs_kind_type(V_10_METER_WIND, 'v_10_meter_wind',                           .false., .false.), &
     obs_kind_type(TEMPERATURE_2_METER, 'temperature_2_meter',                   .false., .false.), &
     obs_kind_type(SPECIFIC_HUMIDITY_2_METER, 'specific_humidity_2_meter',       .false., .false.), &
     obs_kind_type(DEW_POINT_2_METER, 'dew_point_2_meter',                       .false., .false.), &
     obs_kind_type(RAW_STATE_1D_INTEGRAL, 'raw_state_1d_integral',               .false., .false.) /)
! Namelist array to turn on any requested observation types
character(len = 129) :: assimilate_these_obs_types(max_obs_kinds) = 'null'
character(len = 129) :: evaluate_these_obs_types(max_obs_kinds) = 'null'

namelist /obs_def_nml/ assimilate_these_obs_types, evaluate_these_obs_types

contains


  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

integer :: iunit, ierr, io, i, j
character(len = 169) :: err_string

call register_module(source, revision, revdate)
module_initialized = .true.

if(file_exist('input.nml')) then
   iunit = open_file(fname = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = obs_def_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'obs_def_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

call error_handler(E_MSG,'initialize_module','obs_def_nml values are',' ',' ',' ')
write(logfileunit, *) 'assimilate_these_obs_types'
write(*, *) 'ASSIMILATE_THESE_OBS_TYPES'
do i = 1, max_obs_kinds
   if(assimilate_these_obs_types(i) == 'null') goto 22
   write(logfileunit, *) trim(assimilate_these_obs_types(i))
   write(*, *) trim(assimilate_these_obs_types(i))
   num_kind_assimilate = i
end do
22 write(logfileunit, *) 'evaluate_these_obs_types'
write(*, *) 'EVALUATE_THESE_OBS_TYPES'
do i = 1, max_obs_kinds
   if(evaluate_these_obs_types(i) == 'null') goto 33
   write(logfileunit, *) trim(evaluate_these_obs_types(i))
   write(*, *) trim(evaluate_these_obs_types(i))
   num_kind_evaluate = i
end do
33 continue

! Figure out which kinds are being used, look for errors
! Start by loading up kinds to assimilate
do i = 1, num_kind_assimilate
   ! Search for the matching string
   do j = 1, max_obs_kinds
      if(assimilate_these_obs_types(i) == obs_kind_info(j)%name) then
         obs_kind_info(j)%assimilate = .true.
         goto 44
      endif
   end do
   ! Falling off the end is an error
   write(err_string, *) trim(assimilate_these_obs_types(i)), &
      ' from obs_def_nml is not a legal observation kind'
   call error_handler(E_ERR, 'initialize_module', err_string, source, revision, revdate)
   44 continue
end do

! Now look for kinds to evaluate
do i = 1, num_kind_evaluate
   ! Search for the matching string
   do j = 1, max_obs_kinds
      if(evaluate_these_obs_types(i) == obs_kind_info(j)%name) then
         obs_kind_info(j)%evaluate = .true.
         goto 55
      endif
   end do
   ! Falling off the end is an error
   write(err_string, *) trim(evaluate_these_obs_types(i)), &
      ' from obs_def_nml is not a legal observation kind'
   call error_handler(E_ERR, 'initialize_module', err_string, source, revision, revdate)
   55 continue
end do

! Make it an error to ask to assimilate AND evaluate the same obs kind
do i = 1, max_obs_kinds
   if(obs_kind_info(i)%evaluate .and. obs_kind_info(i)%assimilate) then
      write(err_string, *) 'Illegal to evaluate and assimilate same kind ', trim(obs_kind_info(i)%name)
      call error_handler(E_ERR, 'initialize_module', err_string, source, revision, revdate)
   endif
end do

end subroutine initialize_module


!----------------------------------------------------------------------------

subroutine init_obs_def(obs_def, location, kind, time, error_variance)
! Need to add additional component arguments as optionals as needed

! Constructor for an obs_def

type(obs_def_type), intent(out) :: obs_def
type(location_type), intent(in) :: location
integer,             intent(in) :: kind
type(time_type),     intent(in) :: time
real(r8),            intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

obs_def%location = location
obs_def%kind = kind
obs_def%time = time
obs_def%error_variance = error_variance
! No key assigned for standard observation defs
obs_def%key = -1

end subroutine init_obs_def

!---------------------------------------------------------------------

subroutine copy_obs_def(obs_def1, obs_def2)

! Copy function to be overloaded with '='

type(obs_def_type), intent(out) :: obs_def1
type(obs_def_type), intent(in) :: obs_def2

if ( .not. module_initialized ) call initialize_module

obs_def1%location = obs_def2%location
obs_def1%kind = obs_def2%kind
obs_def1%time = obs_def2%time
obs_def1%error_variance = obs_def2%error_variance
obs_def1%key = obs_def2%key
!WRF obs_def1%platform = obs_def2%platform
!deallocate(obs_def1%platform_qc)
!allocate(obs_def1%platform_qc(size(obs_def2%platform_qc))
! Should this be pointer assignment or regular
!obs_def1%platform_qc >= or == obs_def2%platform_qc
!obs_def1%aperture = obs_def2%aperture

end subroutine copy_obs_def

!----------------------------------------------------------------------------

function get_obs_def_error_variance(obs_def)

type(obs_def_type), intent(in) :: obs_def
real(r8)                       :: get_obs_def_error_variance

if ( .not. module_initialized ) call initialize_module

get_obs_def_error_variance = obs_def%error_variance

end function get_obs_def_error_variance

!----------------------------------------------------------------------------

function get_obs_def_location(obs_def)

! Returns observation location.

type(location_type)            :: get_obs_def_location
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_location = obs_def%location

end function get_obs_def_location

!----------------------------------------------------------------------------

function get_obs_kind(obs_def)

! Returns observation kind

integer                        :: get_obs_kind
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_kind = obs_def%kind

end function get_obs_kind

!----------------------------------------------------------------------------

function get_obs_def_time(obs_def)

! Returns observation time

type(time_type)                :: get_obs_def_time
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_time = obs_def%time

end function get_obs_def_time

!----------------------------------------------------------------------------

!WRF function get_obs_def_platform(obs_def)

! Returns the platform of an obs_def

!WRF type(platform_type)            :: get_obs_def_platform
!WRF type(obs_def_type), intent(in) :: obs_def

!WRF if ( .not. module_initialized ) call initialize_module

!WRF get_obs_def_platform = obs_def%platform

!WRF end function get_obs_def_platform

!----------------------------------------------------------------------------

subroutine set_obs_def_location(obs_def, location)

! Sets the location of an obs_def

type(obs_def_type), intent(inout) :: obs_def
type(location_type),   intent(in) :: location

if ( .not. module_initialized ) call initialize_module

obs_def%location = location

end subroutine set_obs_def_location

!----------------------------------------------------------------------------

subroutine set_obs_def_error_variance(obs_def, error_variance)

! Sets the error variance of an obs_def

type(obs_def_type), intent(inout) :: obs_def
real(r8), intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

obs_def%error_variance = error_variance

end subroutine set_obs_def_error_variance

!----------------------------------------------------------------------------

subroutine set_obs_def_kind(obs_def, kind)

! Sets the kind of an obs_def

type(obs_def_type), intent(inout) :: obs_def
integer,               intent(in) :: kind

if ( .not. module_initialized ) call initialize_module

obs_def%kind = kind

end subroutine set_obs_def_kind

!----------------------------------------------------------------------------

subroutine set_obs_def_time(obs_def, time)

! Sets the time of an obs_def

type(obs_def_type), intent(inout) :: obs_def
type(time_type), intent(in) :: time

if ( .not. module_initialized ) call initialize_module

obs_def%time = time

end subroutine set_obs_def_time

!----------------------------------------------------------------------------

!WRF subroutine set_obs_def_platform(obs_def, platform)

! Sets the platform of an obs_def

!WRF type(obs_def_type),  intent(inout) :: obs_def
!WRF type(platform_type), intent(in)    :: platform

!WRF if ( .not. module_initialized ) call initialize_module

!WRF obs_def%platform = platform

!WRF end subroutine set_obs_def_platform


!----------------------------------------------------------------------------

subroutine get_expected_obs_from_def(key, obs_def, obs_kind_ind, state, obs_val, &
   istatus, assimilate_this_ob, evaluate_this_ob)

! Compute forward operator for a particular obs_def
integer, intent(in) :: key
type(obs_def_type), intent(in) :: obs_def
integer, intent(in) :: obs_kind_ind
real(r8), intent(in) :: state(:)
real(r8), intent(out) :: obs_val
integer, intent(out) :: istatus
logical, intent(out) :: assimilate_this_ob, evaluate_this_ob

type(location_type) :: location

! Load up the assimilate and evaluate status for this observation kind
assimilate_this_ob = obs_kind_info(obs_kind_ind)%assimilate
evaluate_this_ob = obs_kind_info(obs_kind_ind)%evaluate

! If not being assimilated or evaluated return with missing_r8 and istatus -99???
if(obs_kind_info(obs_kind_ind)%assimilate .or. obs_kind_info(obs_kind_ind)%evaluate) then
   location = get_obs_def_location(obs_def)
   ! Compute the forward operator;
   select case(obs_kind_ind)
      #IFDEF raw_state_variable
         case(RAW_STATE_VARIABLE)
         call interpolate(state, location, 1, obs_val, istatus)
      #ENDIF
      #IFDEF raw_state_1d_integral
         case(RAW_STATE_1D_INTEGRAL)
            call get_expected_1d_integral(state, location, obs_def%key, obs_val, istatus)
      #ENDIF
      #IFDEF radiosonde_u_wind_component
         case(RADIOSONDE_U_WIND_COMPONENT)
         !!!call interpolate(state, location, TYPE_U, obs_val, istatus)
         call interpolate(state, location, 1, obs_val, istatus)
      #ENDIF
      #IFDEF radiosonde_v_wind_component
         case(RADIOSONDE_V_WIND_COMPONENT)
         !!!call interpolate(state, location, TYPE_V, obs_val, istatus)
         call interpolate(state, location, 2, obs_val, istatus)
      #ENDIF
      #IFDEF radiosonde_temperature
         case(SURFACE_PRESSURE)
         !!!call interpolate(state, location, TYPE_PS, obs_val, istatus)
         call interpolate(state, location, 3, obs_val, istatus)
      #ENDIF
      #IFDEF radiosonde_temperature
         case(RADIOSONDE_TEMPERATURE)
         !!!call interpolate(state, location, TYPE_T, obs_val, istatus)
         call interpolate(state, location, 4, obs_val, istatus)
      #ENDIF
   end select
else
   ! Not computing forward operator
   obs_val = missing_r8
   istatus = -99
endif

end subroutine get_expected_obs_from_def

!----------------------------------------------------------------------------

subroutine read_obs_def(ifile, obs_def, key, fform)

! Reads an obs_def from file which is just an integer unit number in the
! current preliminary implementation.

type(obs_def_type),      intent(inout) :: obs_def
integer,                 intent(in)    :: ifile
integer,                 intent(in)    :: key
character(len=*), intent(in), optional :: fform

character(len=5)  :: header
character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Begin by reading five character ascii header, then location, kind, error variance, index

! Need to add additional error checks on read
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      continue
   CASE DEFAULT
      read(ifile, 11) header
11    Format(a5)
      if(header /= 'obdef') then
   call error_handler(E_ERR,'read_obs_def', 'Expected location header "obdef" in input file', &
                      source, revision, revdate)
      endif
END SELECT

! Read the location, kind, time and error variance
obs_def%location = read_location(ifile, fileformat)
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) obs_def%kind
   CASE DEFAULT
      read(ifile, '(a5)' ) header
      if(header /= 'kind ') then
         call error_handler(E_ERR,'read_kind', 'Expected kind header "kind " in input file', &
                            source, revision, revdate)
      endif
      read(ifile, *) obs_def%kind
END SELECT

! This kind may have its own module that needs to read more
select case(obs_def%kind)
   case(RAW_STATE_VARIABLE)
   case(RADIOSONDE_U_WIND_COMPONENT)
   case(RADIOSONDE_V_WIND_COMPONENT)
   case(RADIOSONDE_TEMPERATURE)
   case(SURFACE_PRESSURE)
   #IFDEF raw_state_1d_integral
      case(RAW_STATE_1D_INTEGRAL)
         call read_1d_integral(obs_def%key, ifile, fileformat)
   #ENDIF
end select

! Read the time for the observation
obs_def%time = read_time(ifile, fileformat)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) obs_def%error_variance
   CASE DEFAULT
      read(ifile, *) obs_def%error_variance
END SELECT

!WRF call read_platform(ifile, obs_def%platform, fform)

end subroutine read_obs_def

!----------------------------------------------------------------------------

subroutine write_obs_def(ifile, obs_def, key, fform)

! Writes an obs_def to file.

integer,                    intent(in) :: ifile
type(obs_def_type),         intent(in) :: obs_def
integer,                    intent(in) :: key
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Write the 5 character identifier for verbose formatted output
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      continue
   CASE DEFAULT
      write(ifile, 11)
11    format('obdef')
END SELECT

! Write out the location, kind and error variance
call write_location(ifile, obs_def%location, fileformat)
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) obs_def%kind
   CASE DEFAULT
      write(ifile, '(''kind'')' )
      write(ifile, *) obs_def%kind
END SELECT

! This kind may have its own module that needs to write more
select case(obs_def%kind)
   case(RAW_STATE_VARIABLE)
   case(RADIOSONDE_U_WIND_COMPONENT)
   case(RADIOSONDE_V_WIND_COMPONENT)
   case(RADIOSONDE_TEMPERATURE)
   case(SURFACE_PRESSURE)
   #IFDEF raw_state_1d_integral
      case(RAW_STATE_1D_INTEGRAL)
         call write_1d_integral(obs_def%key, ifile, fileformat)
   #ENDIF
end select

call write_time(ifile, obs_def%time, fileformat)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) obs_def%error_variance
   CASE DEFAULT
      write(ifile, *) obs_def%error_variance
END SELECT

end subroutine write_obs_def


subroutine interactive_obs_def(obs_def, key)
!---------------------------------------------------------------------------
!
! Allows interactive creation of an observation

type(obs_def_type), intent(inout) :: obs_def
integer,               intent(in) :: key

integer :: i

if ( .not. module_initialized ) call initialize_module

! Get the observation kind WANT A STRING OPTION, TOO?
write(*, *) 'Input the index for the observation kind'
write(*, *) '-1 * state variable index for identity obs OR '
do i = 1, max_obs_kinds
   if(obs_kind_info(i)%assimilate .or. obs_kind_info(i)%evaluate) &
      write(*, *) obs_kind_info(i)%index, trim(obs_kind_info(i)%name)
end do
read(*, *) obs_def%kind


! Input any special stuff for this kind
select case(obs_def%kind)
   #IFDEF doppler_radial_velocity
      call interactive_doppler_radial_velocity()
   #ENDIF
   #IFDEFF radar_reflectivity
   #ENDIF radar_reflectivity
   #IFDEF raw_state_1d_integral
      case(RAW_STATE_1D_INTEGRAL)
         call interactive_1d_integral(obs_def%key)
   #ENDIF
end select



! If the kind is an identity observation, don't need to call location
! Get location from state meta_data
if(obs_def%kind < 0) then
! Get the location of this from model
   call get_state_meta_data(-1 * obs_def%kind, obs_def%location)
else! Get the location
   call interactive_location(obs_def%location)
endif

! Get the time
call interactive_time(obs_def%time)

write(*, *) 'Input error variance for this observation definition '
read(*, *) obs_def%error_variance

! TJH -- might want to do some sort of error checking (i.e. for positive values)

end subroutine interactive_obs_def


subroutine set_radar_obs_def(rad_loc,rgate,raz,elev_rad,ae,dir,var,obs_def)
!---------------------------------------------------------------------------
!
! Allows creation of a radar observation

type(location_type),    intent(in)    :: rad_loc
real(r8),               intent(in)    :: rgate,raz,elev_rad,ae,var
real(r8),               intent(in)    :: dir(3)
type(obs_def_type),     intent(inout) :: obs_def

real(r8) :: h, spath, x, y, rad_lon, rad_lat, obs_lon, obs_lat, vloc

if ( .not. module_initialized ) call initialize_module

! Set the observation kind
obs_def%kind = KIND_VR

! Doviak & Zrniv, 1993: Doppler radar and weather observations, eq. 2.28b-c

h = sqrt( rgate*rgate + ae*ae + 2.0_r8 *rgate*ae*sin(elev_rad) ) - ae
spath = ae * asin(rgate * cos(elev_rad) / (ae + h))

x = spath*sin(raz)
y = spath*cos(raz)

!WRF rad_lon = query_location(rad_loc, 'lon')
!WRF rad_lat = query_location(rad_loc, 'lat')

obs_lat = y/ae + rad_lat
obs_lon = x/(ae*cos(rad_lat + y/(2.0_r8*ae))) + rad_lon

!WRF vloc = query_location(rad_loc, 'vloc')

vloc = vloc + h

obs_lon = obs_lon*RAD2DEG
obs_lat = obs_lat*RAD2DEG

!WRF obs_def%location = set_location(obs_lon, obs_lat, vloc, 3)

! Set the time
obs_def%time = set_time(0, 0)

obs_def%error_variance = var

!WRF call set_platform_location(obs_def%platform, rad_loc)
!WRF call set_platform_orientation(obs_def%platform, dir)

!WRF call set_obs_def_platform(obs_def, obs_def%platform)

end subroutine set_radar_obs_def


!----------------------------------------------------------------


subroutine destroy_obs_def(obs_def)
! TECHNICALLY NEED TO CALL DESTRUCTORS FOR ALL SUBCOMPONENTS, NO ALLOCATED STORAGE YET
! obs_def_type has the following components:
! type(location_type) :: location      ! center of mass, so to speak
! type(obs_kind_type) :: kind          ! keyword, BUFR values for now
! type(time_type)     :: time
! real(r8)            :: error_variance

type(obs_def_type), intent(inout) :: obs_def

if ( .not. module_initialized ) call initialize_module

call set_obs_def_location(       obs_def, set_location_missing() )
obs_def%kind = missing_i
call set_obs_def_time(           obs_def, set_time_missing() )
call set_obs_def_error_variance( obs_def, missing_r8)

end subroutine destroy_obs_def

!-------------------------------------------------------

end module obs_def_mod
