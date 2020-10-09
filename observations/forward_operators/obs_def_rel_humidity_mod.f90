! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! BEGIN DART PREPROCESS TYPE DEFINITIONS
!RADIOSONDE_RELATIVE_HUMIDITY,    QTY_RELATIVE_HUMIDITY
!DROPSONDE_RELATIVE_HUMIDITY,     QTY_RELATIVE_HUMIDITY
!AIRCRAFT_RELATIVE_HUMIDITY,      QTY_RELATIVE_HUMIDITY
!ACARS_RELATIVE_HUMIDITY,         QTY_RELATIVE_HUMIDITY
!MARINE_SFC_RELATIVE_HUMIDITY,    QTY_RELATIVE_HUMIDITY
!LAND_SFC_RELATIVE_HUMIDITY,      QTY_RELATIVE_HUMIDITY
!METAR_RELATIVE_HUMIDITY_2_METER, QTY_RELATIVE_HUMIDITY
!AIRS_RELATIVE_HUMIDITY,          QTY_RELATIVE_HUMIDITY
!MESONET_RELATIVE_HUMIDITY,       QTY_RELATIVE_HUMIDITY
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_rel_humidity_mod, only : get_expected_relative_humidity
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY, MESONET_RELATIVE_HUMIDITY)
!            call get_expected_relative_humidity(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY, MESONET_RELATIVE_HUMIDITY)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY, MESONET_RELATIVE_HUMIDITY)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY, MESONET_RELATIVE_HUMIDITY)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_rel_humidity_mod

use             types_mod, only : r8, missing_r8, L_over_Rv
use         utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, E_ALLMSG
use          location_mod, only : location_type, set_location, get_location, write_location, &
                                  read_location, is_vertical
use       assim_model_mod, only : interpolate
use          obs_kind_mod, only : QTY_TEMPERATURE, QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO
use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: get_expected_relative_humidity

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_rel_humidity_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save       :: module_initialized   = .false.
logical, save       :: first_time_warn_low  = .true.
logical, save       :: first_time_warn_high = .true.
character(len=512)  :: msgstring
real(r8), parameter :: MIN_VALUE = 1.0e-9
real(r8), parameter :: MAX_VALUE = 1.1

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------------

subroutine get_expected_relative_humidity(state_handle, ens_size, location, rh, istatus)

type(ensemble_type), intent(in)     :: state_handle
integer,             intent(in)     :: ens_size
type(location_type), intent(in)     :: location
real(r8),            intent(out)    :: rh(ens_size)    ! relative humidity (fraction)
integer,             intent(out)    :: istatus(ens_size)

real(r8), dimension(ens_size) :: qvap, tmpk, pres, es, qsat
integer,  dimension(ens_size) :: qvap_istatus, tmpk_istatus, pres_istatus
real(r8)                      :: xyz(3)
logical :: return_now

if ( .not. module_initialized ) call initialize_module

istatus = 0

!  interpolate the mixing ratio to the location
call interpolate(state_handle, ens_size, location, QTY_VAPOR_MIXING_RATIO, qvap, qvap_istatus)
call track_status(ens_size, qvap_istatus, rh, istatus, return_now)

where (istatus == 0 .and. qvap < 0.0_r8) qvap = epsilon(0.0_r8)  ! clamping?
where (istatus /= 0) istatus = 99
if(return_now) return

!  interpolate the temperature to the desired location
call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, tmpk, tmpk_istatus)
call track_status(ens_size, tmpk_istatus, rh, istatus, return_now)

where (tmpk <= 0.0_r8)
   rh = missing_r8
   istatus = 99
endwhere

if(all(istatus /= 0)) return  ! not using return_now because where clause modifies istatus.

!  interpolate the pressure, if observation location is not pressure
if ( is_vertical(location, "PRESSURE") ) then
   xyz = get_location(location)
   pres = xyz(3)
else
   ! pressure comes back in pascals (not hPa or mb)
   call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pres, pres_istatus)
   call track_status(ens_size, pres_istatus, rh, istatus, return_now)

   where (pres <= 0.0_r8 .or. pres >= 120000.0_r8)
      rh = missing_r8
      istatus = 99
   end where
   if(all(istatus /= 0)) return  ! not using return_now because where clause modifies istatus.

endif

!  Compute the rh
where (istatus == 0)
   es   = 611.0_r8 * exp(L_over_Rv * (1.0_r8 / 273.0_r8 - 1.0_r8 / tmpk))
   qsat = 0.622_r8 * es / (pres - es)
   rh   = qvap / qsat
end where

! Warnings - first time only

! nsc - my dilemma: each task does these obs independently.  if you use E_ALLMSG you'll
! get a message from every processor which runs into this case.  if you use E_MSG you'll
! only see those on task 0.  for now i've removed the test and you'll get the message (once)
! even if there aren't any values outside the valid range.  i can't see a good way to
! enforce that you get only one message when this could happen on any task and probably
! will happen on many tasks.

if (first_time_warn_low) then
   !if (any(rh < MIN_VALUE .and. istatus == 0)) then
   !   write(msgstring, '(A,F12.6)') 'values lower than low limit detected, e.g.', minval(rh, istatus==0)
      write(msgstring, *) 'checking relative humidity value computed by the forward operator'
      call error_handler(E_MSG,'get_expected_relative_humidity', msgstring,      &
                         text2='all values lower than 1e-9 will be set to 1e-9', &
                         text3='this message will only print once')
      first_time_warn_low = .false.
   !endif
endif
if (first_time_warn_high) then
   !if (any(rh > MAX_VALUE .and. istatus == 0)) then
   !   write(msgstring, '(A,F12.6)') 'values higher than high limit detected, e.g.', maxval(rh, istatus==0)
      write(msgstring, *) 'checking relative humidity value computed by the forward operator'
      call error_handler(E_MSG,'get_expected_relative_humidity', msgstring,      &
                         text2='all values larger than 1.1 will be set to 1.1', &
                         text3='this message will only print once')
      first_time_warn_high = .false.
   !endif
endif

! Actually force the relative humidity to be between MIN_VALUE and MAX_VALUE
where(rh < MIN_VALUE .and. istatus == 0) rh = MIN_VALUE
where(rh > MAX_VALUE .and. istatus == 0) rh = MAX_VALUE

return
end subroutine get_expected_relative_humidity

!----------------------------------------------------------------------------

end module obs_def_rel_humidity_mod
! END DART PREPROCESS MODULE CODE

