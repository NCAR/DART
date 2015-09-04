! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
!RADIOSONDE_RELATIVE_HUMIDITY,    KIND_RELATIVE_HUMIDITY
!DROPSONDE_RELATIVE_HUMIDITY,     KIND_RELATIVE_HUMIDITY
!AIRCRAFT_RELATIVE_HUMIDITY,      KIND_RELATIVE_HUMIDITY
!ACARS_RELATIVE_HUMIDITY,         KIND_RELATIVE_HUMIDITY
!MARINE_SFC_RELATIVE_HUMIDITY,    KIND_RELATIVE_HUMIDITY
!LAND_SFC_RELATIVE_HUMIDITY,      KIND_RELATIVE_HUMIDITY
!METAR_RELATIVE_HUMIDITY_2_METER, KIND_RELATIVE_HUMIDITY
!AIRS_RELATIVE_HUMIDITY,          KIND_RELATIVE_HUMIDITY
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_rel_humidity_mod, only : get_expected_relative_humidity_distrib
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY)
!            call get_expected_relative_humidity_distrib(state_ens_handle, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(RADIOSONDE_RELATIVE_HUMIDITY, DROPSONDE_RELATIVE_HUMIDITY, &
!              AIRCRAFT_RELATIVE_HUMIDITY,   ACARS_RELATIVE_HUMIDITY,     &
!              MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_RELATIVE_HUMIDITY,  &
!              METAR_RELATIVE_HUMIDITY_2_METER, AIRS_RELATIVE_HUMIDITY)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_rel_humidity_mod

use        types_mod, only : r8, missing_r8, L_over_Rv
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, write_location, &
                             read_location, vert_is_pressure
use  assim_model_mod, only : interpolate_distrib
use     obs_kind_mod, only : KIND_TEMPERATURE, KIND_PRESSURE, KIND_VAPOR_MIXING_RATIO

use ensemble_manager_mod, only : ensemble_type, copies_in_window

implicit none
private

public :: get_expected_relative_humidity_distrib

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save       :: module_initialized   = .false.
logical, save       :: first_time_warn_low  = .true.
logical, save       :: first_time_warn_high = .true.
character(len=64)   :: msgstring
real(r8), parameter :: MIN_VALUE = 1.0e-9
real(r8), parameter :: MAX_VALUE = 1.1

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------------

subroutine get_expected_relative_humidity_distrib(state_ens_handle, location, rh, istatus)

type(ensemble_type), intent(in)     :: state_ens_handle
type(location_type), intent(in)     :: location
real(r8),            intent(out)    :: rh(:)    ! relative humidity (fraction)
integer,             intent(out)    :: istatus(:)

real(r8), allocatable :: qvap(:), tmpk(:), pres(:), es(:), qsat(:)
real(r8)              :: xyz(3)
integer,  allocatable :: track_status(:)
integer              :: e, ens_size

if ( .not. module_initialized ) call initialize_module

ens_size = copies_in_window(state_ens_handle)

allocate(qvap(ens_size), tmpk(ens_size), pres(ens_size), es(ens_size), qsat(ens_size))
allocate(track_status(ens_size))

!  interpolate the mixing ratio to the location
call interpolate_distrib(location, KIND_VAPOR_MIXING_RATIO, istatus, qvap, state_ens_handle)
do e = 1, ens_size
   if (istatus(e) /= 0 .or. qvap(e) < 0.0_r8) then
      if (istatus(e) == 0) then
         qvap(e) = epsilon(0.0_r8)
      else
         rh(e) = missing_r8
         if (istatus(e) == 0) istatus(e) = 99 !HK early return?
      endif
   endif
enddo

track_status = istatus

!  interpolate the temperature to the desired location
call interpolate_distrib(location, KIND_TEMPERATURE, istatus, tmpk, state_ens_handle)
do e = 1, ens_size
   if (istatus(e) /= 0 .or. tmpk(e) <= 0.0_r8) then
      rh(e) = missing_r8
      if (istatus(e) == 0) istatus(e) = 99
   endif
enddo

do e = 1, ens_size
   if (istatus(e) /= 0 ) track_status(e) = istatus(e)
enddo

!  interpolate the pressure, if observation location is not pressure
if ( vert_is_pressure(location) ) then
   xyz = get_location(location)
   pres = xyz(3)
else
   ! pressure comes back in pascals (not hPa or mb)
   call interpolate_distrib(location, KIND_PRESSURE, istatus, pres, state_ens_handle)
   do e = 1, ens_size
      if (istatus(e) /= 0 .or. pres(e) <= 0.0_r8 .or. pres(e) >= 120000.0_r8)  then
         rh(e) = missing_r8
         if (istatus(e) == 0) istatus(e) = 99
      endif
   enddo
endif

do e = 1, ens_size
   if (istatus(e) /= 0 ) track_status(e) = istatus(e)
enddo

!  Compute the rh - HK is having missing_r8 going to cause problems
es   = 611.0_r8 * exp(L_over_Rv * (1.0_r8 / 273.0_r8 - 1.0_r8 / tmpk))
qsat = 0.622_r8 * es / (pres - es)
rh   = qvap / qsat

do e = 1, ens_size

   if (track_status(e) == 0 ) then
      !  Check for unreasonable values and return limits
      if (rh(e) < MIN_VALUE) then
         if (first_time_warn_low) then
            write(msgstring, '(A,F12.6)') 'values lower than low limit detected, e.g.', rh
            call error_handler(E_MSG,'get_expected_relative_humidity', msgstring,      &
                               text2='all values lower than 1e-9 will be set to 1e-9', &
                               text3='this message will only print once')
            first_time_warn_low = .false.
         endif
         rh(e) = MIN_VALUE
      endif

      if (rh(e) > MAX_VALUE) then
         if (first_time_warn_high) then
            write(msgstring, '(A,F12.6)') 'values higher than high limit detected, e.g.', rh
            call error_handler(E_MSG,'get_expected_relative_humidity', msgstring,      &
                            text2='all values larger than 1.1 will be set to 1.1', &
                               text3='this message will only print once')
            first_time_warn_high = .false.
         endif
         rh(e) = MAX_VALUE
      endif
   else
      rh(e) = missing_r8 ! put missing values back in
   endif
enddo

deallocate(qvap, tmpk, pres, es, qsat)

return
end subroutine get_expected_relative_humidity_distrib

!----------------------------------------------------------------------------

end module obs_def_rel_humidity_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
