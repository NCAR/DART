! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
!BUOY_U_WIND_COMPONENT,        KIND_U_WIND_COMPONENT,      COMMON_CODE
!BUOY_V_WIND_COMPONENT,        KIND_V_WIND_COMPONENT,      COMMON_CODE
!BUOY_SURFACE_PRESSURE,        KIND_SURFACE_PRESSURE,      COMMON_CODE
!BUOY_TEMPERATURE,             KIND_TEMPERATURE,           COMMON_CODE
!SHIP_U_WIND_COMPONENT,        KIND_U_WIND_COMPONENT,      COMMON_CODE
!SHIP_V_WIND_COMPONENT,        KIND_V_WIND_COMPONENT,      COMMON_CODE
!SHIP_SURFACE_PRESSURE,        KIND_SURFACE_PRESSURE,      COMMON_CODE
!SHIP_TEMPERATURE,             KIND_TEMPERATURE,           COMMON_CODE
!SYNOP_U_WIND_COMPONENT,       KIND_U_WIND_COMPONENT,      COMMON_CODE
!SYNOP_V_WIND_COMPONENT,       KIND_V_WIND_COMPONENT,      COMMON_CODE
!SYNOP_SURFACE_PRESSURE,       KIND_SURFACE_PRESSURE,      COMMON_CODE
!SYNOP_SPECIFIC_HUMIDITY,      KIND_SPECIFIC_HUMIDITY,     COMMON_CODE
!SYNOP_TEMPERATURE,            KIND_TEMPERATURE,           COMMON_CODE
!AIREP_U_WIND_COMPONENT,       KIND_U_WIND_COMPONENT,      COMMON_CODE
!AIREP_V_WIND_COMPONENT,       KIND_V_WIND_COMPONENT,      COMMON_CODE
!AIREP_PRESSURE,               KIND_PRESSURE,              COMMON_CODE
!AIREP_TEMPERATURE,            KIND_TEMPERATURE,           COMMON_CODE
!AMDAR_U_WIND_COMPONENT,       KIND_U_WIND_COMPONENT,      COMMON_CODE
!AMDAR_V_WIND_COMPONENT,       KIND_V_WIND_COMPONENT,      COMMON_CODE
!AMDAR_PRESSURE,               KIND_PRESSURE,              COMMON_CODE
!AMDAR_TEMPERATURE,            KIND_TEMPERATURE,           COMMON_CODE
!PILOT_U_WIND_COMPONENT,       KIND_U_WIND_COMPONENT,      COMMON_CODE
!PILOT_V_WIND_COMPONENT,       KIND_V_WIND_COMPONENT,      COMMON_CODE
!PILOT_PRESSURE,               KIND_PRESSURE,              COMMON_CODE
!PILOT_TEMPERATURE,            KIND_TEMPERATURE,           COMMON_CODE
!BOGUS_U_WIND_COMPONENT,       KIND_U_WIND_COMPONENT,      COMMON_CODE
!BOGUS_V_WIND_COMPONENT,       KIND_V_WIND_COMPONENT,      COMMON_CODE
!BOGUS_PRESSURE,               KIND_PRESSURE,              COMMON_CODE
!BOGUS_TEMPERATURE,            KIND_TEMPERATURE,           COMMON_CODE
!PROFILER_U_WIND_COMPONENT,    KIND_U_WIND_COMPONENT,      COMMON_CODE
!PROFILER_V_WIND_COMPONENT,    KIND_V_WIND_COMPONENT,      COMMON_CODE
!PROFILER_PRESSURE,            KIND_PRESSURE,              COMMON_CODE
!SATEM_THICKNESS,              KIND_TEMPERATURE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_gts_mod, only : get_expected_thickness_distrib
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(SATEM_THICKNESS)
!            call get_expected_thickness_distrib(state_ens_handle, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! case(SATEM_THICKNESS)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(SATEM_THICKNESS)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(SATEM_THICKNESS)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_gts_mod

use        types_mod, only : r8, missing_r8, gravity, gas_constant, gas_constant_v
use    utilities_mod, only : register_module
use     location_mod, only : location_type, set_location, get_location , &
                             VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use  assim_model_mod, only : interpolate_distrib
use     obs_kind_mod, only : KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, KIND_PRESSURE

use data_structure_mod, only : ensemble_type

implicit none
private

public :: get_expected_thickness_distrib

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

contains



subroutine initialize_module
!-----------------------------------------------------------------------------
call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module



subroutine get_expected_thickness_distrib(state_ens_handle, location, thickness, istatus)
!-----------------------------------------------------------------------------
! inputs:
!    state_vector:    DART state vector
!
! output parameters:
!    thickness: modeled satem thickness between the specified layers. (unit: m)
!    istatus:  =0 normal; =1 outside of domain.
!---------------------------------------------
! Hui Liu  NCAR/IMAGE  April 9, 2008
!---------------------------------------------
implicit none

type(ensemble_type)              :: state_ens_handle
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: thickness(:)
integer,             intent(out) :: istatus(:)

! local variables

integer :: k
integer, parameter :: nlevel = 9, num_press_int = 10
real(r8) :: lon, lat, pressure, obsloc(3), obs_level(nlevel), press(num_press_int+1)
real(r8) :: press_top, press_bot, press_int

integer,  allocatable :: istatus0(:), istatus2(:), track_status(:)

data  obs_level/85000.0, 70000.0, 50000.0, 40000.0, 30000.0, 25000.0, 20000.0, 15000.0, 10000.0/

real(r8), parameter::  rdrv1 = gas_constant_v/gas_constant - 1.0_r8

real(r8) :: lon2, p
type(location_type) :: location2
integer :: which_vert
real(r8), allocatable :: t(:), q(:), tv(:)

integer :: e, ens_size

if ( .not. module_initialized ) call initialize_module

ens_size = state_ens_handle%num_copies -5

allocate(istatus0(ens_size), istatus2(ens_size), track_status(ens_size))

obsloc   = get_location(location)
lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90
pressure = obsloc(3)                       ! (Pa)

press_bot = -999.9_r8
press_top = -999.9_r8

! find the bottom and top of the layer -
! The 'bottom' is defined to be the incoming observation value.
! The 'top' is the next higher observation level.
! If the bottom is the highest observation level;
! there is no top, and it is an ERROR.

if(pressure > obs_level(1) ) then
   press_bot = pressure
   press_top = obs_level(1)
endif

do k = 1, nlevel-1
   if( abs(pressure - obs_level(k)) < 0.01_r8 ) then
      press_bot = obs_level(k)
      press_top = obs_level(k+1)
      exit
   endif
end do

!  define the vertical intervals for the integration of Tv. Use 10 sub_layers
!  for the vertical integration.

press_int = (press_bot - press_top)/num_press_int

press(1) = press_bot
do k=2, num_press_int+1
   press(k) =  press(1) - press_int * (k-1)
end do

!   define the location of the interpolated points.

lon2 = lon
if(lon > 360.0_r8 ) lon2 = lon - 360.0_r8
if(lon <   0.0_r8 ) lon2 = lon + 360.0_r8

which_vert = VERTISPRESSURE
istatus0  = 0
istatus2  = 0
thickness = 0.0_r8

track_status = 0

do k=2, num_press_int+1, 2
   p = press(k)
   location2 = set_location(lon2, lat, p,  which_vert)

   call interpolate_distrib(location2,  KIND_TEMPERATURE,       istatus0, t, state_ens_handle)
   call interpolate_distrib(location2,  KIND_SPECIFIC_HUMIDITY, istatus2, q, state_ens_handle)

   if(all(istatus0 > 0) .or. all(istatus2 > 0)  ) then
      istatus = 1
      thickness = missing_r8
   endif

   ! t :  Kelvin, from top to bottom
   ! q :  kg/kg, from top to bottom

   tv    = t * (1.0_r8 + rdrv1 * q )         ! virtual temperature
   thickness = thickness + gas_constant/gravity * tv * log(press(k-1)/press(k+1))

   do e = 1, ens_size
      ! expected values are around 1 KM; this is a check for
      ! wildly out of range values.
      if( abs(thickness(e)) > 10000.0_r8 ) then
         istatus(e) = 2
         thickness(e) = missing_r8
      endif
      if (istatus(e) /= 0) track_status(e) = istatus(e)

   enddo

end do

istatus = track_status

return
end subroutine get_expected_thickness_distrib

end module obs_def_gts_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
