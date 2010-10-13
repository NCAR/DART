! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! BEGIN DART PREPROCESS KIND LIST
! DEWPOINT,                KIND_DEWPOINT
! DEWPOINT_2_METER,        KIND_DEWPOINT
! BUOY_DEWPOINT,           KIND_DEWPOINT
! SHIP_DEWPOINT,           KIND_DEWPOINT
! SYNOP_DEWPOINT,          KIND_DEWPOINT
! AIREP_DEWPOINT,          KIND_DEWPOINT
! AMDAR_DEWPOINT,          KIND_DEWPOINT
! PILOT_DEWPOINT,          KIND_DEWPOINT
! BOGUS_DEWPOINT,          KIND_DEWPOINT
! METAR_DEWPOINT_2_METER,  KIND_DEWPOINT
! RADIOSONDE_DEWPOINT,     KIND_DEWPOINT
! DROPSONDE_DEWPOINT,      KIND_DEWPOINT
! AIRCRAFT_DEWPOINT,       KIND_DEWPOINT
! ACARS_DEWPOINT,          KIND_DEWPOINT
! MARINE_SFC_DEWPOINT,     KIND_DEWPOINT
! LAND_SFC_DEWPOINT,       KIND_DEWPOINT
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_dew_point_mod, only : get_expected_dew_point
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(DEWPOINT)
!            call get_expected_dew_point(state, location, 1, obs_val, istatus)
!         case(AIREP_DEWPOINT, AMDAR_DEWPOINT, PILOT_DEWPOINT, BOGUS_DEWPOINT)
!            call get_expected_dew_point(state, location, 1, obs_val, istatus)
!         case(RADIOSONDE_DEWPOINT, AIRCRAFT_DEWPOINT, ACARS_DEWPOINT, DROPSONDE_DEWPOINT)
!            call get_expected_dew_point(state, location, 1, obs_val, istatus)
!         case(DEWPOINT_2_METER)
!            call get_expected_dew_point(state, location, 2, obs_val, istatus)
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT)
!            call get_expected_dew_point(state, location, 2, obs_val, istatus)
!         case(MARINE_SFC_DEWPOINT, LAND_SFC_DEWPOINT)
!            call get_expected_dew_point(state, location, 2, obs_val, istatus)
!         case(METAR_DEWPOINT_2_METER)
!            call get_expected_dew_point(state, location, 2, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(DEWPOINT, DEWPOINT_2_METER)
!            continue
!         case(METAR_DEWPOINT_2_METER)
!            continue
!         case(AIREP_DEWPOINT, AMDAR_DEWPOINT, PILOT_DEWPOINT, BOGUS_DEWPOINT)
!            continue
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT)
!            continue
!         case(RADIOSONDE_DEWPOINT, AIRCRAFT_DEWPOINT, ACARS_DEWPOINT, DROPSONDE_DEWPOINT)
!            continue
!         case(MARINE_SFC_DEWPOINT, LAND_SFC_DEWPOINT)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(DEWPOINT, DEWPOINT_2_METER)
!            continue
!         case(METAR_DEWPOINT_2_METER)
!            continue
!         case(AIREP_DEWPOINT, AMDAR_DEWPOINT, PILOT_DEWPOINT, BOGUS_DEWPOINT)
!            continue
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT)
!            continue
!         case(RADIOSONDE_DEWPOINT, AIRCRAFT_DEWPOINT, ACARS_DEWPOINT, DROPSONDE_DEWPOINT)
!            continue
!         case(MARINE_SFC_DEWPOINT, LAND_SFC_DEWPOINT)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(DEWPOINT, DEWPOINT_2_METER)
!            continue
!         case(METAR_DEWPOINT_2_METER)
!            continue
!         case(AIREP_DEWPOINT, AMDAR_DEWPOINT, PILOT_DEWPOINT, BOGUS_DEWPOINT)
!            continue
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT)
!            continue
!         case(RADIOSONDE_DEWPOINT, AIRCRAFT_DEWPOINT, ACARS_DEWPOINT, DROPSONDE_DEWPOINT)
!            continue
!         case(MARINE_SFC_DEWPOINT, LAND_SFC_DEWPOINT)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_dew_point_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, missing_r8, t_kelvin
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location , write_location, &
                             read_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_SURFACE_PRESSURE, KIND_VAPOR_MIXING_RATIO, KIND_PRESSURE

implicit none
private

public :: get_expected_dew_point

! version controlled file description for error handling, do not edit
character(len=128) :: &
source   = "$URL$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module



subroutine get_expected_dew_point(state_vector, location, key, td, istatus)

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: td               ! dewpoint (K)
integer,             intent(out) :: istatus

integer  :: ipres
real(r8) :: qv                            ! water vapor mixing ratio (kg/kg)
real(r8) :: e_mb                          ! water vapor pressure (mb)
real(r8), PARAMETER :: e_min = 0.001_r8   ! threshold for minimum vapor pressure (mb),
                                          !   to avoid problems near zero in Bolton's equation
real(r8) :: p_Pa                          ! pressure (Pa)
real(r8) :: p_mb                          ! pressure (mb)

character(len=129) :: errstring

if ( .not. module_initialized ) call initialize_module

if(key == 1) then
   ipres = KIND_PRESSURE
elseif(key == 2) then
   ipres = KIND_SURFACE_PRESSURE
else
   write(errstring,*)'key has to be 1 (upper levels) or 2 (2-meter), got ',key
   call error_handler(E_ERR,'get_expected_dew_point', errstring, &
        source, revision, revdate)
endif

call interpolate(state_vector, location, ipres, p_Pa, istatus)
if (istatus /= 0) then
   td = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_VAPOR_MIXING_RATIO, qv, istatus)
if (istatus /= 0) then
   td = missing_r8
   return
endif
if (qv < 0.0_r8 .or. qv >= 1.0_r8) then
   td = missing_r8
   if (istatus == 0) istatus = 1
   return
endif

!------------------------------------------------------------------------------
!  Compute water vapor pressure.
!------------------------------------------------------------------------------

p_mb = p_Pa * 0.01_r8

e_mb = qv * p_mb / (0.622_r8 + qv)
e_mb = max(e_mb, e_min)

!------------------------------------------------------------------------------
!  Use Bolton's approximation to compute dewpoint.
!------------------------------------------------------------------------------

td = t_kelvin + (243.5_r8 / ((17.67_r8 / log(e_mb/6.112_r8)) - 1.0_r8) )

end subroutine get_expected_dew_point

!----------------------------------------------------------------------------

end module obs_def_dew_point_mod
! END DART PREPROCESS MODULE CODE
