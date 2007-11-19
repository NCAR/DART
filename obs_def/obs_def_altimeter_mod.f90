module obs_def_altimeter_mod

! BEGIN DART PREPROCESS KIND LIST
! RADIOSONDE_SURFACE_ALTIMETER, KIND_SURFACE_PRESSURE
! MARINE_SFC_ALTIMETER, KIND_SURFACE_PRESSURE
! LAND_SFC_ALTIMETER, KIND_SURFACE_PRESSURE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_altimeter_mod, only : get_expected_altimeter, compute_altimeter
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, LAND_SFC_ALTIMETER)
!            call get_expected_altimeter(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, LAND_SFC_ALTIMETER)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, LAND_SFC_ALTIMETER)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, LAND_SFC_ALTIMETER)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF


use        types_mod, only : r8, missing_r8, t_kelvin
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location , write_location, &
                             read_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_SURFACE_PRESSURE, KIND_SURFACE_ELEVATION

implicit none
private

public :: get_expected_altimeter, compute_altimeter

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_altimeter_mod.f90,v $", &
revision = "$Revision: 3093 $", &
revdate  = "$Date: 2007-09-01 16:01:48 -0600 (Sat, 01 Sep 2007) $"

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module


subroutine get_expected_altimeter(state_vector, location, altimeter_setting, istatus)

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: altimeter_setting     ! altimeter (hPa)
integer,             intent(out) :: istatus

real(r8) :: psfc                ! surface pressure value   (Pa)
real(r8) :: hsfc                ! surface elevation level  (m above SL)

if ( .not. module_initialized ) call initialize_module

!  interpolate the surface pressure to the desired location
call interpolate(state_vector, location, KIND_SURFACE_PRESSURE, psfc, istatus)
if (istatus /= 0) then
   altimeter_setting = missing_r8
   return
endif

!  interpolate the surface elevation to the desired location
call interpolate(state_vector, location, KIND_SURFACE_ELEVATION, hsfc, istatus)
if (istatus /= 0) then
   altimeter_setting = missing_r8
   return
endif

!  Compute the altimeter setting given surface pressure and height, altimeter is hPa
altimeter_setting = compute_altimeter(psfc * 0.01, hsfc)

if (altimeter_setting < 880.0_r8 .or. altimeter_setting >= 1100.0_r8) then
   altimeter_setting = missing_r8
   if (istatus == 0) istatus = 1
   return
endif

return
end subroutine get_expected_altimeter


function compute_altimeter(psfc, hsfc)

real(r8), parameter :: k1 = 0.190284_r8
real(r8), parameter :: k2 = 8.4228807E-5_r8

real(r8), intent(in) :: psfc  !  (hPa)
real(r8), intent(in) :: hsfc  !  (m above MSL)

real(r8) :: compute_altimeter !  (hPa)

compute_altimeter = ((psfc - 0.3) ** k1 + k2 * hsfc) ** (1 / k1)

return
end function compute_altimeter

!----------------------------------------------------------------------------

end module obs_def_altimeter_mod
