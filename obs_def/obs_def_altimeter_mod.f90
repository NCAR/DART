! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! RADIOSONDE_SURFACE_ALTIMETER, KIND_SURFACE_PRESSURE
! DROPSONDE_SURFACE_ALTIMETER,  KIND_SURFACE_PRESSURE
! MARINE_SFC_ALTIMETER,         KIND_SURFACE_PRESSURE
! LAND_SFC_ALTIMETER,           KIND_SURFACE_PRESSURE
! METAR_ALTIMETER,              KIND_SURFACE_PRESSURE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_altimeter_mod, only : get_expected_altimeter_distrib, compute_altimeter_distrib
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER)
!            call get_expected_altimeter_distrib(state_ens_handle, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(RADIOSONDE_SURFACE_ALTIMETER, DROPSONDE_SURFACE_ALTIMETER, MARINE_SFC_ALTIMETER, &
!              LAND_SFC_ALTIMETER, METAR_ALTIMETER)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_altimeter_mod

use           types_mod, only : r8, missing_r8, t_kelvin
use       utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use        location_mod, only : location_type, set_location, get_location, write_location, &
                                read_location
use     assim_model_mod, only : interpolate_distrib
use        obs_kind_mod, only : KIND_SURFACE_PRESSURE, KIND_SURFACE_ELEVATION
use  data_structure_mod, only : ensemble_type, copies_in_window

implicit none
private

public :: get_expected_altimeter_distrib, compute_altimeter_distrib

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

subroutine get_expected_altimeter_distrib(state_ens_handle, location, altimeter_setting, istatus)

type(ensemble_type), intent(in)     :: state_ens_handle
type(location_type), intent(in)     :: location
real(r8),            intent(out)    :: altimeter_setting(:)     ! altimeter (hPa)
integer,             intent(out)    :: istatus(:)

real(r8), allocatable :: psfc(:)                ! surface pressure value   (Pa)
real(r8), allocatable :: hsfc(:)                ! surface elevation level  (m above SL)
integer,  allocatable :: track_status(:)
integer               :: e, ens_size

ens_size = copies_in_window(state_ens_handle)

if ( .not. module_initialized ) call initialize_module

allocate(psfc(ens_size), hsfc(ens_size), track_status(ens_size))

track_status = 0
altimeter_setting = missing_r8

!  interpolate the surface pressure to the desired location
call interpolate_distrib(location, KIND_SURFACE_PRESSURE, istatus, psfc, state_ens_handle)
do e = 1, ens_size
   if (istatus(e) /= 0) then
      altimeter_setting(e) = missing_r8
   endif
enddo
track_status = istatus

!  interpolate the surface elevation to the desired location
call interpolate_distrib(location, KIND_SURFACE_ELEVATION, istatus, hsfc, state_ens_handle)
do e = 1, ens_size
   if (istatus(e) /= 0) then
      altimeter_setting(e) = missing_r8
      track_status(e) = istatus(e)
   endif
enddo

!  Compute the altimeter setting given surface pressure and height, altimeter is hPa
do e = 1, ens_size ! looping to avoid FPE with missing_r8
   if(track_status(e) == 0) altimeter_setting(e:e) = compute_altimeter_distrib(psfc(e:e) * 0.01_r8, hsfc(e:e), 1)
enddo

istatus = track_status

do e = 1, ens_size
   if (altimeter_setting(e) < 880.0_r8 .or. altimeter_setting(e) >= 1100.0_r8) then
      altimeter_setting(e) = missing_r8
      if (istatus(e) == 0) istatus(e) = 1
   endif
enddo

deallocate(psfc, hsfc, track_status)

return
end subroutine get_expected_altimeter_distrib

!----------------------------------------------------------------------

function compute_altimeter_distrib(psfc, hsfc, n)

real(r8), parameter :: k1 = 0.190284_r8
real(r8), parameter :: k2 = 8.4228807E-5_r8

integer,  intent(in) :: n
real(r8), intent(in) :: psfc(n)  !  (hPa)
real(r8), intent(in) :: hsfc(n)  !  (m above MSL)

real(r8) :: compute_altimeter_distrib(n) !  (hPa)

compute_altimeter_distrib = ((psfc - 0.3_r8) ** k1 + k2 * hsfc) ** (1.0_r8 / k1)

return
end function compute_altimeter_distrib

!----------------------------------------------------------------------------

end module obs_def_altimeter_mod

! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
