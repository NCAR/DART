! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!>@todo  make a namelist-defined grid (lower left corner, nx, ny, dx, dy)
!> allocate, call interp for t, u, v, p and fill it.  store the location
!> of the obs being interpolated.  compute the values.  wait for subsequent
!> calls for this same location.
!>

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! VORTEX_LAT, QTY_VORTEX_LAT
! VORTEX_LON, QTY_VORTEX_LON
! VORTEX_PMIN, QTY_VORTEX_PMIN
! VORTEX_WMAX, QTY_VORTEX_WMAX
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_vortex_mod, only : get_expected_vortex_info
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(VORTEX_LAT)
!            call get_expected_vortex_info(state_handle, ens_size, location, 'lat', expected_obs, istatus)
!         case(VORTEX_LON)
!            call get_expected_vortex_info(state_handle, ens_size, location, 'lon', expected_obs, istatus)
!         case(VORTEX_PMIN)
!            call get_expected_vortex_info(state_handle, ens_size, location, 'pmi', expected_obs, istatus)
!         case(VORTEX_WMAX)
!            call get_expected_vortex_info(state_handle, ens_size, location, 'wma', expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(VORTEX_LAT)
!         continue
!      case(VORTEX_LON)
!         continue
!      case(VORTEX_PMIN)
!         continue
!      case(VORTEX_WMAX)
!         continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(VORTEX_LAT)
!         continue
!      case(VORTEX_LON)
!         continue
!      case(VORTEX_PMIN)
!         continue
!      case(VORTEX_WMAX)
!         continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(VORTEX_LAT)
!         continue
!      case(VORTEX_LON)
!         continue
!      case(VORTEX_PMIN)
!         continue
!      case(VORTEX_WMAX)
!         continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_vortex_mod

use        types_mod, only : r8, missing_r8, ps0, PI, gravity
use    utilities_mod, only : register_module, error_handler, E_ERR
use     location_mod, only : location_type
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, &
                             QTY_TEMPERATURE, QTY_VERTICAL_VELOCITY, &
                             QTY_RAINWATER_MIXING_RATIO, QTY_DENSITY, &
                             QTY_VORTEX_LAT, QTY_VORTEX_LON, QTY_VORTEX_PMIN, &
                             QTY_VORTEX_WMAX

use ensemble_manager_mod,  only : ensemble_type

implicit none
private

public :: get_expected_vortex_info

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

! Storage for the special information required for observations of this type
integer, parameter  :: max_vortex_obs = 3000

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine get_expected_vortex_info(state_handle, ens_size, location, whichinfo, vinfo, istatus)

!
! Return vortex info according to whichinfo
! whichinfo='lat', vinfo = center lat
! whichinfo='lon', vinfo = center lon
! whichinfo='pmi', vinfo =  minimum sea level pressure
!

type(ensemble_type), intent(in)     :: state_handle
integer,             intent(in)     :: ens_size
type(location_type), intent(in)     :: location
character(len=3),    intent(in)     :: whichinfo
real(r8),            intent(out)    :: vinfo(ens_size)
integer,             intent(out)    :: istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

if (whichinfo == 'lat') then
   call interpolate(state_handle, ens_size, location, QTY_VORTEX_LAT, vinfo, istatus)
elseif (whichinfo == 'lon') then
   call interpolate(state_handle, ens_size, location, QTY_VORTEX_LON, vinfo, istatus)
elseif (whichinfo == 'pmi') then
   call interpolate(state_handle, ens_size, location, QTY_VORTEX_PMIN, vinfo, istatus)
elseif (whichinfo == 'wma') then
   call interpolate(state_handle, ens_size, location, QTY_VORTEX_WMAX, vinfo, istatus)
endif

end subroutine get_expected_vortex_info

end module obs_def_vortex_mod
! END DART PREPROCESS MODULE CODE

