! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! DEWPOINT,                QTY_DEWPOINT
! DEWPOINT_2_METER,        QTY_DEWPOINT
! BUOY_DEWPOINT,           QTY_DEWPOINT
! SHIP_DEWPOINT,           QTY_DEWPOINT
! SYNOP_DEWPOINT,          QTY_DEWPOINT
! AIREP_DEWPOINT,          QTY_DEWPOINT
! AMDAR_DEWPOINT,          QTY_DEWPOINT
! PILOT_DEWPOINT,          QTY_DEWPOINT
! BOGUS_DEWPOINT,          QTY_DEWPOINT
! AIRS_DEWPOINT,           QTY_DEWPOINT
! METAR_DEWPOINT_2_METER,  QTY_DEWPOINT
! RADIOSONDE_DEWPOINT,     QTY_DEWPOINT
! DROPSONDE_DEWPOINT,      QTY_DEWPOINT
! AIRCRAFT_DEWPOINT,       QTY_DEWPOINT
! ACARS_DEWPOINT,          QTY_DEWPOINT
! MARINE_SFC_DEWPOINT,     QTY_DEWPOINT
! LAND_SFC_DEWPOINT,       QTY_DEWPOINT
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_dew_point_mod, only : get_expected_dew_point
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(DEWPOINT)
!            call get_expected_dew_point(state_handle,  ens_size, location, 1, expected_obs, istatus)
!         case(AIREP_DEWPOINT, AMDAR_DEWPOINT, PILOT_DEWPOINT, BOGUS_DEWPOINT, AIRS_DEWPOINT)
!            call get_expected_dew_point(state_handle, ens_size,  location, 1, expected_obs, istatus)
!         case(RADIOSONDE_DEWPOINT, AIRCRAFT_DEWPOINT, ACARS_DEWPOINT, DROPSONDE_DEWPOINT)
!            call get_expected_dew_point(state_handle,  ens_size, location, 1, expected_obs, istatus)
!
!         case(DEWPOINT_2_METER)
!            call get_expected_dew_point(state_handle,  ens_size, location, 2, expected_obs, istatus)
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT)
!            call get_expected_dew_point(state_handle,  ens_size, location, 2, expected_obs, istatus)
!         case(MARINE_SFC_DEWPOINT, LAND_SFC_DEWPOINT)
!            call get_expected_dew_point(state_handle,  ens_size, location, 2, expected_obs, istatus)
!         case(METAR_DEWPOINT_2_METER)
!            call get_expected_dew_point(state_handle,  ens_size, location, 2, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(DEWPOINT, DEWPOINT_2_METER)
!            continue
!         case(METAR_DEWPOINT_2_METER)
!            continue
!         case(AIREP_DEWPOINT, AMDAR_DEWPOINT, PILOT_DEWPOINT, BOGUS_DEWPOINT)
!            continue
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT, AIRS_DEWPOINT)
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
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT, AIRS_DEWPOINT)
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
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT, AIRS_DEWPOINT)
!            continue
!         case(RADIOSONDE_DEWPOINT, AIRCRAFT_DEWPOINT, ACARS_DEWPOINT, DROPSONDE_DEWPOINT)
!            continue
!         case(MARINE_SFC_DEWPOINT, LAND_SFC_DEWPOINT)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_dew_point_mod

use        types_mod, only : r8, missing_r8, t_kelvin
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location , write_location, &
                             read_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_SURFACE_PRESSURE, QTY_VAPOR_MIXING_RATIO, QTY_PRESSURE

use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: get_expected_dew_point

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------
subroutine get_expected_dew_point(state_handle, ens_size, location, key, td, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(in)    :: location
integer,             intent(in)    :: key
real(r8),            intent(out)   :: td(ens_size)   ! dewpoint (K)
integer,             intent(out)   :: istatus(ens_size)

integer  :: ipres
real(r8) :: qv(ens_size)            ! water vapor mixing ratio (kg/kg)
real(r8) :: e_mb(ens_size)          ! water vapor pressure (mb)
real(r8),   PARAMETER :: e_min = 0.001_r8 ! threshold for minimum vapor pressure (mb),
                                          !   to avoid problems near zero vapor pressure in Bolton's equation
real(r8) :: p_Pa(ens_size)          ! pressure (Pa)
real(r8) :: p_mb(ens_size)          ! pressure (mb)
real(r8) :: log_term(ens_size)      ! Intermediate term in computation of dewpoint temperature
!> @todo make strings longer
character(len=129) :: errstring
logical  :: return_now
integer  :: qv_istatus(ens_size)
integer  :: ipres_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

istatus = 0

if(key == 1) then
   ipres = QTY_PRESSURE
elseif(key == 2) then
   ipres = QTY_SURFACE_PRESSURE
else
   write(errstring,*)'key has to be 1 (upper levels) or 2 (2-meter), got ',key
   call error_handler(E_ERR,'get_expected_dew_point', errstring, &
        source, revision, revdate)
endif

call interpolate(state_handle, ens_size, location, ipres, p_Pa, ipres_istatus)
call track_status(ens_size, ipres_istatus, td, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size,location, QTY_VAPOR_MIXING_RATIO, qv, qv_istatus)
call track_status(ens_size, qv_istatus, td, istatus, return_now)
if (return_now) return


where (qv < 0.0_r8 .or. qv >= 1.0_r8)
   istatus = 1
   td = missing_r8
end where



!------------------------------------------------------------------------------
!  Compute water vapor pressure.
!------------------------------------------------------------------------------
where (istatus == 0 ) ! To avoid possible FPE with missing_r8

   p_mb = p_Pa * 0.01_r8

   e_mb = qv * p_mb / (0.622_r8 + qv)
   e_mb = max(e_mb, e_min)

   !------------------------------------------------------------------------------
   !  Use Bolton's approximation to compute dewpoint.
   !  Bolton, David, 1980: The Computation of Equivalent Potential Temperature. 
   !  Monthly Weather Review. 108 (7): 1046-1053.
   !  Bolton does not explicitly have this formula, but the pieces are there.
   !  He uses the three constant values used here.
   !  Other authors, and the Wikipedia page on dewpoint, provide this formula
   !  explicitly and suggest 6.1121 instead of 6.112.
   !------------------------------------------------------------------------------

   ! The following expression can fail numerically for dewpoints very close to 0 C
   !td = t_kelvin + (243.5_r8 / ((17.67_r8 / log(e_mb/6.112_r8)) - 1.0_r8) )

   ! A numerically robust formula that avoids the failure near dewpoints of 0 C
   log_term = log(e_mb / 6.112_r8)
   td = t_kelvin + 243.5_r8 * log_term / (17.67_r8 - log_term)

elsewhere
   td = missing_r8
end where

end subroutine get_expected_dew_point

!----------------------------------------------------------------------------

end module obs_def_dew_point_mod
! END DART PREPROCESS MODULE CODE

