! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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
! AIRS_DEWPOINT,           KIND_DEWPOINT
! METAR_DEWPOINT_2_METER,  KIND_DEWPOINT
! RADIOSONDE_DEWPOINT,     KIND_DEWPOINT
! DROPSONDE_DEWPOINT,      KIND_DEWPOINT
! AIRCRAFT_DEWPOINT,       KIND_DEWPOINT
! ACARS_DEWPOINT,          KIND_DEWPOINT
! MARINE_SFC_DEWPOINT,     KIND_DEWPOINT
! LAND_SFC_DEWPOINT,       KIND_DEWPOINT
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_dew_point_mod, only : get_expected_dew_point_distrib
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(DEWPOINT)
!            call get_expected_dew_point_distrib(state_ens_handle,  location, 1, expected_obs, istatus)
!         case(AIREP_DEWPOINT, AMDAR_DEWPOINT, PILOT_DEWPOINT, BOGUS_DEWPOINT, AIRS_DEWPOINT)
!            call get_expected_dew_point_distrib(state_ens_handle,  location, 1, expected_obs, istatus)
!         case(RADIOSONDE_DEWPOINT, AIRCRAFT_DEWPOINT, ACARS_DEWPOINT, DROPSONDE_DEWPOINT)
!            call get_expected_dew_point_distrib(state_ens_handle,  location, 1, expected_obs, istatus)
!
!         case(DEWPOINT_2_METER)
!            call get_expected_dew_point_distrib(state_ens_handle,  location, 2, expected_obs, istatus)
!         case(BUOY_DEWPOINT, SHIP_DEWPOINT, SYNOP_DEWPOINT)
!            call get_expected_dew_point_distrib(state_ens_handle,  location, 2, expected_obs, istatus)
!         case(MARINE_SFC_DEWPOINT, LAND_SFC_DEWPOINT)
!            call get_expected_dew_point_distrib(state_ens_handle,  location, 2, expected_obs, istatus)
!         case(METAR_DEWPOINT_2_METER)
!            call get_expected_dew_point_distrib(state_ens_handle,  location, 2, expected_obs, istatus)
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
use  assim_model_mod, only : interpolate_distrib
use     obs_kind_mod, only : KIND_SURFACE_PRESSURE, KIND_VAPOR_MIXING_RATIO, KIND_PRESSURE

use data_structure_mod, only : ensemble_type, copies_in_window

implicit none
private

public :: get_expected_dew_point_distrib

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



subroutine get_expected_dew_point_distrib(state_ens_handle,  location, key, td, istatus)

type(ensemble_type)                :: state_ens_handle
type(location_type), intent(in)    :: location
integer,             intent(in)    :: key
real(r8),            intent(out)   :: td(:)              ! dewpoint (K)
integer,             intent(out)   :: istatus(:)

integer  :: ipres
real(r8), allocatable :: qv(:)            ! water vapor mixing ratio (kg/kg)
real(r8), allocatable :: e_mb(:)          ! water vapor pressure (mb)
real(r8),   PARAMETER :: e_min = 0.001_r8 ! threshold for minimum vapor pressure (mb),
                                          !   to avoid problems near zero in Bolton's equation
real(r8), allocatable :: p_Pa(:)          ! pressure (Pa)
real(r8), allocatable :: p_mb(:)          ! pressure (mb)
integer  :: e, ens_size
real(r8), allocatable :: track_status(:)
character(len=129) :: errstring

if ( .not. module_initialized ) call initialize_module

ens_size = copies_in_window(state_ens_handle)
allocate(qv(ens_size), p_Pa(ens_size), p_mb(ens_size), e_mb(ens_size), track_status(ens_size))

if(key == 1) then
   ipres = KIND_PRESSURE
elseif(key == 2) then
   ipres = KIND_SURFACE_PRESSURE
else
   write(errstring,*)'key has to be 1 (upper levels) or 2 (2-meter), got ',key
   call error_handler(E_ERR,'get_expected_dew_point', errstring, &
        source, revision, revdate)
endif

call interpolate_distrib(location, ipres, istatus, p_Pa, state_ens_handle)
if ( all(istatus /= 0) ) then
   td(:) = missing_r8
   return
endif

track_status = istatus

call interpolate_distrib(location, KIND_VAPOR_MIXING_RATIO, istatus, qv, state_ens_handle)
if ( all(istatus /= 0) ) then
   td(:) = missing_r8
   return
endif

do e = 1, ens_size
   if (qv(e) < 0.0_r8 .or. qv(e) >= 1.0_r8) then
      if (istatus(e) == 0) istatus(e) = 1
   endif
enddo
do e = 1, ens_size
   if (istatus(e) /= 0 ) then
      track_status(e) = istatus(e)
      td(e) = missing_r8
   endif
enddo

!------------------------------------------------------------------------------
!  Compute water vapor pressure.
!------------------------------------------------------------------------------
!HK can you get weird problems if some values are missing_r8?
p_mb = p_Pa * 0.01_r8

e_mb = qv * p_mb / (0.622_r8 + qv)
e_mb = max(e_mb, e_min)

!------------------------------------------------------------------------------
!  Use Bolton's approximation to compute dewpoint.
!------------------------------------------------------------------------------

td = t_kelvin + (243.5_r8 / ((17.67_r8 / log(e_mb/6.112_r8)) - 1.0_r8) )

do e = 1, ens_size
  if (track_status(e) /= 0 ) td(e) = missing_r8
enddo

istatus = track_status

end subroutine get_expected_dew_point_distrib

!----------------------------------------------------------------------------

end module obs_def_dew_point_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
