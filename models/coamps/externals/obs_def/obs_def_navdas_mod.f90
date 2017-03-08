! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! GEOPOTENTIAL_HEIGHT,   QTY_GEOPOTENTIAL_HEIGHT
! LOG_SPECIFIC_HUMIDITY, QTY_SPECIFIC_HUMIDITY
! AIR_TEMPERATURE,       QTY_TEMPERATURE
! U_WIND_COMPONENT,      QTY_U_WIND_COMPONENT
! V_WIND_COMPONENT,      QTY_V_WIND_COMPONENT
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_navdas_mod, only : get_expected_geopotential_height,&
!                                  get_expected_log_specific_humidity,&
!                                  get_expected_air_temperature,&
!                                  get_expected_wind_component
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(GEOPOTENTIAL_HEIGHT)
!            call get_expected_geopotential_height(state, location, &
!                                                  obs_val, istatus) 
!         case(LOG_SPECIFIC_HUMIDITY)
!            call get_expected_log_specific_humidity(state, location,&
!                                                    obs_val, istatus) 
!         case(AIR_TEMPERATURE)
!            call get_expected_air_temperature(state, location,&
!                                              obs_val,istatus) 
!         case(U_WIND_COMPONENT)
!            call get_expected_wind_component(state, location,1,&
!                                             obs_val,istatus) 
!         case(V_WIND_COMPONENT)
!            call get_expected_wind_component(state, location,2,&
!                                             obs_val,istatus) 
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(GEOPOTENTIAL_HEIGHT)
!            continue
!         case(LOG_SPECIFIC_HUMIDITY)
!            continue
!         case(AIR_TEMPERATURE)
!            continue
!         case(U_WIND_COMPONENT)
!            continue
!         case(V_WIND_COMPONENT)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(GEOPOTENTIAL_HEIGHT)
!            continue
!         case(LOG_SPECIFIC_HUMIDITY)
!            continue
!         case(AIR_TEMPERATURE)
!            continue
!         case(U_WIND_COMPONENT)
!            continue
!         case(V_WIND_COMPONENT)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(GEOPOTENTIAL_HEIGHT)
!            continue
!         case(LOG_SPECIFIC_HUMIDITY)
!            continue
!         case(AIR_TEMPERATURE)
!            continue
!         case(U_WIND_COMPONENT)
!            continue
!         case(V_WIND_COMPONENT)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_navdas_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

  use        types_mod, only : r8, missing_r8, t_kelvin, gas_constant,ps0
  use    utilities_mod, only : register_module
  use     location_mod, only : location_type, vert_is_pressure,&
                               query_location 
  use  assim_model_mod, only : interpolate
  use     obs_kind_mod, only : QTY_POTENTIAL_TEMPERATURE, &
                               QTY_SPECIFIC_HUMIDITY,     &
                               QTY_U_WIND_COMPONENT,      &
                               QTY_V_WIND_COMPONENT,      &
                               QTY_GEOPOTENTIAL_HEIGHT

  implicit none
  private

  public :: get_expected_geopotential_height
  public :: get_expected_log_specific_humidity
  public :: get_expected_air_temperature
  public :: get_expected_wind_component

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  logical, save :: module_initialized = .false.

  ! For some reason specific heat is not included in types_mod
  real(r8), parameter :: Cp = 1004.0_r8

  ! Switches for the wind component function
  integer, parameter  :: U = 1
  integer, parameter  :: V = 2

contains

  ! initialize_module
  ! -----------------
  ! Handles any NAVDAS module initialization tasks
  subroutine initialize_module
    call register_module(source, revision, revdate)
    module_initialized = .true.
  end subroutine initialize_module

  ! get_expected_geopotential_height
  ! --------------------------------
  ! Forward operator for Z
  subroutine get_expected_geopotential_height(state_vector, location,&
       geop_height, istatus)  
    real(r8),            intent(in)  :: state_vector(:)
    type(location_type), intent(in)  :: location
    real(r8),            intent(out) :: geop_height
    integer,             intent(out) :: istatus

    if ( .not. module_initialized ) call initialize_module
    ! Specific humidity at this location
    call interpolate(state_vector, location, QTY_GEOPOTENTIAL_HEIGHT,&
                     geop_height, istatus)

  end subroutine get_expected_geopotential_height

  ! get_expected_log_specific_humidity
  ! ----------------------------------
  ! Forward operator for ln(q_v)
  subroutine get_expected_log_specific_humidity(state_vector,&
                                                location, ln_qv,&
                                                istatus) 
    real(r8),            intent(in)  :: state_vector(:)
    type(location_type), intent(in)  :: location
    real(r8),            intent(out) :: ln_qv
    integer,             intent(out) :: istatus

    real(r8) :: qv

    if ( .not. module_initialized ) call initialize_module

    ! Specific humidity at this location
    call interpolate(state_vector, location, QTY_SPECIFIC_HUMIDITY,&
                     qv, istatus)
    if (istatus /= 0) then
       ln_qv = missing_r8
       return
    endif
    
    ! We have the specific humidity, now calculate natural log
    ln_qv = log(qv)

  end subroutine get_expected_log_specific_humidity

  ! get_expected_air_temperature
  ! ----------------------------
  ! Forward operator for T (in Kelvin) - note that this currently
  ! assumes that the location is 3D in pressure coordinates
  subroutine get_expected_air_temperature(state_vector, location, t,&
                                          istatus)  
    real(r8),            intent(in)  :: state_vector(:)
    type(location_type), intent(in)  :: location
    real(r8),            intent(out) :: t
    integer,             intent(out) :: istatus

    real(r8) :: pot_t   ! Potential temperature
    real(r8) :: pres    ! Working pressure
    real(r8) :: rocp    ! R over Cp

    if ( .not. module_initialized ) call initialize_module

    ! Potential temperature at this location
    call interpolate(state_vector, location,&
                     QTY_POTENTIAL_TEMPERATURE,pot_t, istatus)
    if (istatus /= 0) then
       t = missing_r8
       return
    endif
    
    ! Pressure - for now, assume that we are working in pressure
    ! coordinates only.  Once we add other capabilities, change this
    ! to call model_interpolate for QTY_PRESSURE
    if (.not. vert_is_pressure(location)) then
       t = missing_r8
       return
    end if
    pres = query_location(location,'VLOC')

    ! Potential temperature is T * (P00/P) ** (R/Cp)
    rocp = gas_constant / Cp
    t = pot_t * (pres / ps0) ** rocp

  end subroutine get_expected_air_temperature

  ! get_expected_wind_component
  ! ---------------------------
  ! Forward operator for U and V (m/s)
  subroutine get_expected_wind_component(state_vector, location, key,&
                                         wind, istatus) 
    real(r8),            intent(in)  :: state_vector(:)
    type(location_type), intent(in)  :: location
    integer,             intent(in)  :: key
    real(r8),            intent(out) :: wind
    integer,             intent(out) :: istatus

    if ( .not. module_initialized ) call initialize_module

    if (key .eq. U) then
       call interpolate(state_vector, location, QTY_U_WIND_COMPONENT,&
                        wind, istatus)
       if (istatus /= 0) then
          wind = missing_r8
          return
       endif
    else if (key .eq. V) then
       call interpolate(state_vector, location, QTY_V_WIND_COMPONENT,&
                        wind, istatus)
       if (istatus /= 0) then
          wind = missing_r8
          return
       endif
    else
       wind = missing_r8
       return
    endif

  end subroutine get_expected_wind_component


end module obs_def_navdas_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
