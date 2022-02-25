! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


! BEGIN DART PREPROCESS TYPE DEFINITIONS
! RADIOSONDE_SURFACE_ALTIMETER, QTY_SURFACE_PRESSURE
! DROPSONDE_SURFACE_ALTIMETER,  QTY_SURFACE_PRESSURE
! MARINE_SFC_ALTIMETER,         QTY_SURFACE_PRESSURE
! LAND_SFC_ALTIMETER,           QTY_SURFACE_PRESSURE
! METAR_ALTIMETER,              QTY_SURFACE_PRESSURE
! MESONET_SURFACE_ALTIMETER,    QTY_SURFACE_PRESSURE
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_altimeter_mod, only : get_expected_altimeter, compute_altimeter
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER, MESONET_SURFACE_ALTIMETER)
!            call get_expected_altimeter(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER, MESONET_SURFACE_ALTIMETER)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER, MESONET_SURFACE_ALTIMETER)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER, MESONET_SURFACE_ALTIMETER)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_altimeter_mod

use             types_mod, only : r8, missing_r8
use         utilities_mod, only : register_module
use          location_mod, only : location_type
use       assim_model_mod, only : interpolate
use          obs_kind_mod, only : QTY_SURFACE_PRESSURE, QTY_SURFACE_ELEVATION
use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: get_expected_altimeter, compute_altimeter

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_altimeter_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine get_expected_altimeter(state_handle, ens_size, location, altimeter_setting, istatus)

type(ensemble_type), intent(in)     :: state_handle
integer,             intent(in)     :: ens_size
type(location_type), intent(in)     :: location
real(r8),            intent(out)    :: altimeter_setting(ens_size)     ! altimeter (hPa)
integer,             intent(out)    :: istatus(ens_size)

real(r8) :: psfc(ens_size)                ! surface pressure value   (Pa)
real(r8) :: hsfc(ens_size)                ! surface elevation level  (m above SL)
integer :: psfc_istatus(ens_size)
integer :: hsfc_istatus(ens_size)
logical :: return_now
integer :: imem

real(r8), parameter :: MIN_REALISTIC_ALTIMETER = 880.0_r8      ! (hPa)
real(r8), parameter :: MAX_REALISTIC_ALTIMETER = 1100.0_r8     ! (hPa)

if ( .not. module_initialized ) call initialize_module

istatus = 0 ! Need to initialize this to zero for track_status.
altimeter_setting = missing_r8

! interpolate the surface pressure to the desired location. 
! The implicit units of pressure here are Pascals.
call interpolate(state_handle, ens_size, location, QTY_SURFACE_PRESSURE, psfc, psfc_istatus)
call track_status(ens_size, psfc_istatus, altimeter_setting, istatus, return_now)
if (return_now) return

!  interpolate the surface elevation to the desired location
call interpolate(state_handle, ens_size, location, QTY_SURFACE_ELEVATION, hsfc, hsfc_istatus)
call track_status(ens_size, hsfc_istatus, altimeter_setting, istatus, return_now)
if (return_now) return

!  Compute the altimeter setting given surface pressure (hPa) and height (m), altimeter is hPa
do imem = 1, ens_size
   if(istatus(imem) == 0) altimeter_setting(imem) = compute_altimeter(psfc(imem) * 0.01_r8, hsfc(imem))
enddo

where (altimeter_setting <  MIN_REALISTIC_ALTIMETER .or. &
       altimeter_setting >= MAX_REALISTIC_ALTIMETER)
   altimeter_setting = MISSING_R8
   istatus = 1
endwhere

end subroutine get_expected_altimeter

!----------------------------------------------------------------------

function compute_altimeter(psfc, hsfc)

real(r8), parameter :: k1 = 0.190284_r8
real(r8), parameter :: k2 = 8.4228807E-5_r8

real(r8), intent(in) :: psfc  !  (hPa)
real(r8), intent(in) :: hsfc  !  (m above MSL)

real(r8) :: compute_altimeter !  (hPa)

compute_altimeter = ((psfc - 0.3_r8) ** k1 + k2 * hsfc) ** (1.0_r8 / k1)

return
end function compute_altimeter

!----------------------------------------------------------------------------

end module obs_def_altimeter_mod

! END DART PREPROCESS MODULE CODE

