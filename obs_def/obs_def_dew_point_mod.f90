! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_dew_point_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_dew_point_mod.f90,v $
! $Revision$
! $Date$
! $Author$
! $Name:  $

! BEGIN DART PREPROCESS KIND LIST
! DEW_POINT_TEMPERATURE, KIND_DEW_POINT_TEMPERATURE
! DEW_POINT_2_METER, KIND_DEW_POINT_TEMPERATURE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_dew_point_mod, only : get_expected_dew_point
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(DEW_POINT_TEMPERATURE)
!            call get_expected_dew_point(state, location, 1, obs_val, istatus)
!         case(DEW_POINT_2_METER)
!            call get_expected_dew_point(state, location, 2, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(DEW_POINT_TEMPERATURE, DEW_POINT_2_METER)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(DEW_POINT_TEMPERATURE, DEW_POINT_2_METER)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(DEW_POINT_TEMPERATURE, DEW_POINT_2_METER)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

use        types_mod, only : r8, missing_r8, t_kelvin
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location , write_location, &
                             read_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_SURFACE_PRESSURE, KIND_VAPOR_MIXING_RATIO, KIND_PRESSURE

implicit none
private

public :: get_expected_dew_point

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_dew_point_mod.f90,v $", &
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
