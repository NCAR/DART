! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! BEGIN DART PREPROCESS TYPE DEFINITIONS
!   MPD_ABSOLUTE_HUMIDITY,    QTY_ABSOLUTE_HUMIDITY
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_abs_humidity_mod, only : get_expected_absolute_humidity
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(MPD_ABSOLUTE_HUMIDITY)
!            call get_expected_absolute_humidity(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(MPD_ABSOLUTE_HUMIDITY)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(MPD_ABSOLUTE_HUMIDITY)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(MPD_ABSOLUTE_HUMIDITY)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE

module obs_def_abs_humidity_mod

! Module contributed by Michael Ying. Thank You Michael! 

use        types_mod, only : r8, missing_r8, gas_constant, gas_constant_v
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, E_ALLMSG
use     location_mod, only : location_type, set_location, get_location, write_location, &
                             read_location, is_vertical
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_TEMPERATURE, QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO

use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: get_expected_absolute_humidity

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_abs_humidity_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save       :: module_initialized   = .false.
logical, save       :: first_time_warn_low  = .true.
logical, save       :: first_time_warn_high = .true.
character(len=512)  :: string1, string2
real(r8), parameter :: MIN_VALUE = 0.0
real(r8), parameter :: MAX_VALUE = 0.1

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------------

subroutine get_expected_absolute_humidity(state_handle, ens_size, location, ah, istatus)

type(ensemble_type), intent(in)     :: state_handle
integer,             intent(in)     :: ens_size
type(location_type), intent(in)     :: location
real(r8),            intent(out)    :: ah(ens_size)    ! absolute humidity
integer,             intent(out)    :: istatus(ens_size)

real(r8), dimension(ens_size) :: qvap, tmpk, pres, rho_dry, tmpv
integer,  dimension(ens_size) :: qvap_istatus, tmpk_istatus, pres_istatus
real(r8)                      :: xyz(3)
logical :: return_now

if ( .not. module_initialized ) call initialize_module

istatus = 0

!  interpolate the mixing ratio to the location
call interpolate(state_handle, ens_size, location, QTY_VAPOR_MIXING_RATIO, qvap, qvap_istatus)
call track_status(ens_size, qvap_istatus, ah, istatus, return_now)

where (istatus == 0 .and. qvap < 0.0_r8) qvap = epsilon(0.0_r8)  ! clamping?
where (istatus /= 0) istatus = 99
if(return_now) return

!  interpolate the temperature to the desired location
call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, tmpk, tmpk_istatus)
call track_status(ens_size, tmpk_istatus, ah, istatus, return_now)

where (tmpk <= 0.0_r8)
   ah = missing_r8
   istatus = 99
endwhere

if(all(istatus /= 0)) return

!  interpolate the pressure, if observation location is not pressure
if ( is_vertical(location, "PRESSURE") ) then
   xyz = get_location(location)
   pres = xyz(3)
else
   ! pressure comes back in pascals (not hPa or mb)
   call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pres, pres_istatus)
   call track_status(ens_size, pres_istatus, ah, istatus, return_now)

   where (pres <= 0.0_r8 .or. pres >= 120000.0_r8)
      ah = missing_r8
      istatus = 99
   end where
   if(all(istatus /= 0)) return
endif

!  Compute the ah
where (istatus == 0)
   tmpv = tmpk * (1.0_r8 + qvap * gas_constant_v / gas_constant)
   rho_dry = pres / (gas_constant * tmpv)
   ah = qvap * rho_dry
end where

! Warnings - first time only
if (first_time_warn_low) then
   if (any(ah < MIN_VALUE .and. istatus == 0)) then
      write(string1, *) 'checking absolute humidity value computed by the forward operator'
      write(string2,*)'all values lower than ',MIN_VALUE,' will be set to ',MIN_VALUE
      call error_handler(E_MSG,'get_expected_absolute_humidity', string1, &
                         text2=string2, text3='this message will only print once')
      first_time_warn_low = .false.
   endif
endif
if (first_time_warn_high) then
   if (any(ah > MAX_VALUE .and. istatus == 0)) then
      write(string1, *) 'checking absolute humidity value computed by the forward operator'
      write(string2,*)'all values larger than ',MAX_VALUE,' will be set to ',MAX_VALUE
      call error_handler(E_MSG,'get_expected_absolute_humidity', string1, &
                         text2=string2, text3='this message will only print once')
      first_time_warn_high = .false.
   endif
endif

! Actually force the absolute humidity to be between MIN_VALUE and MAX_VALUE
where(ah < MIN_VALUE .and. istatus == 0) ah = MIN_VALUE
where(ah > MAX_VALUE .and. istatus == 0) ah = MAX_VALUE

return
end subroutine get_expected_absolute_humidity

!----------------------------------------------------------------------------

end module obs_def_abs_humidity_mod
! END DART PREPROCESS MODULE CODE
