! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! RADIOSONDE_GEOPOTENTIAL_HEIGHT,   QTY_GEOPOTENTIAL_HEIGHT,    COMMON_CODE
! RADIOSONDE_LOG_SPECIFIC_HUMIDITY, QTY_SPECIFIC_HUMIDITY,      COMMON_CODE
! RADIOSONDE_POT_TEMPERATURE,       QTY_POTENTIAL_TEMPERATURE,  COMMON_CODE
! RADIOSONDE_U_WIND_COMPONENT,      QTY_U_WIND_COMPONENT,       COMMON_CODE
! RADIOSONDE_V_WIND_COMPONENT,      QTY_V_WIND_COMPONENT,       COMMON_CODE
! RADIOSONDE_AIR_TEMPERATURE,       QTY_TEMPERATURE
! EVAL_GEOPOTENTIAL_HEIGHT,         QTY_GEOPOTENTIAL_HEIGHT,    COMMON_CODE
! EVAL_LOG_SPECIFIC_HUMIDITY,       QTY_SPECIFIC_HUMIDITY,      COMMON_CODE
! EVAL_POT_TEMPERATURE,             QTY_POTENTIAL_TEMPERATURE,  COMMON_CODE
! EVAL_U_WIND_COMPONENT,            QTY_U_WIND_COMPONENT,       COMMON_CODE
! EVAL_V_WIND_COMPONENT,            QTY_V_WIND_COMPONENT,       COMMON_CODE
! EVAL_AIR_TEMPERATURE,             QTY_TEMPERATURE
! Z_TC_SYNTH,  QTY_GEOPOTENTIAL_HEIGHT,    COMMON_CODE
! U_TC_SYNTH,  QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_TC_SYNTH,  QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_TC_SYNTH,  QTY_TEMPERATURE
! Z_RAOB,      QTY_GEOPOTENTIAL_HEIGHT,    COMMON_CODE
! U_RAOB,      QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_RAOB,      QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_RAOB,      QTY_TEMPERATURE
! U_DROPSONDE, QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_DROPSONDE, QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_DROPSONDE, QTY_TEMPERATURE
! Z_PIBAL,     QTY_GEOPOTENTIAL_HEIGHT,    COMMON_CODE
! U_PIBAL,     QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_PIBAL,     QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_PIBAL,     QTY_TEMPERATURE
! U_AIREP,     QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_AIREP,     QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_AIREP,     QTY_TEMPERATURE
! U_AMDAR,     QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_AMDAR,     QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_AMDAR,     QTY_TEMPERATURE
! U_ACARS,     QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_ACARS,     QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_ACARS,     QTY_TEMPERATURE
! U_MDCRS,     QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_MDCRS,     QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_MDCRS,     QTY_TEMPERATURE
! U_CLD_WNDS1, QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_CLD_WNDS1, QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_CLD_WNDS2, QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_CLD_WNDS2, QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_METEO,     QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_METEO,     QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_GOES,      QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_GOES,      QTY_V_WIND_COMPONENT,       COMMON_CODE
! T_TOVS_T,    QTY_TEMPERATURE
! P_SFC_LAND,  QTY_SURFACE_PRESSURE      
! U_SFC_LAND,  QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_SFC_LAND,  QTY_V_WIND_COMPONENT,       COMMON_CODE
! P_SFC_SHIP,  QTY_SURFACE_PRESSURE      
! U_SFC_SHIP,  QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_SFC_SHIP,  QTY_V_WIND_COMPONENT,       COMMON_CODE
! VORTEX_LON,  QTY_VORTEX_LON,             COMMON_CODE
! VORTEX_LAT,  QTY_VORTEX_LAT,             COMMON_CODE
! U_WSAT_WIND, QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_WSAT_WIND, QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_SCAT_WINDS,QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_SCAT_WINDS,QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_QSCAT_WIND,QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_QSCAT_WIND,QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_ASCAT_WIND,QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_ASCAT_WIND,QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_TERRA_WIND,QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_TERRA_WIND,QTY_V_WIND_COMPONENT,       COMMON_CODE
! U_AQUA_WIND, QTY_U_WIND_COMPONENT,       COMMON_CODE
! V_AQUA_WIND, QTY_V_WIND_COMPONENT,       COMMON_CODE
! WS_SSMI_FF1, QTY_VELOCITY
! WS_SSMI_FF2, QTY_VELOCITY
! TPPW_SSMI,   QTY_TOTAL_PRECIPITABLE_WATER, COMMON_CODE
! TPPW_WSAT,   QTY_TOTAL_PRECIPITABLE_WATER, COMMON_CODE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_navdas_mod, only : &
!                                  get_expected_altimeter,&
!                                  get_expected_air_temperature,&
!                                  get_expected_windspeed
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(RADIOSONDE_AIR_TEMPERATURE,EVAL_AIR_TEMPERATURE, &
!        T_RAOB,T_PIBAL,T_AIREP,T_AMDAR,T_ACARS,T_MDCRS, &
!        T_TOVS_T, T_TC_SYNTH, T_DROPSONDE)
!        call get_expected_air_temperature(state, location, obs_val, istatus) 
!   case(P_SFC_LAND, P_SFC_SHIP)
!        call get_expected_altimeter(state, location, obs_val, istatus) 
!   case(WS_SSMI_FF1, WS_SSMI_FF2)
!        call get_expected_windspeed(state, location, obs_val, istatus) 
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(RADIOSONDE_AIR_TEMPERATURE,EVAL_AIR_TEMPERATURE, &
!        T_RAOB,T_PIBAL,T_AIREP,T_AMDAR,T_ACARS,T_MDCRS, &
!        T_TOVS_T, T_TC_SYNTH, T_DROPSONDE)
!     continue
!   case(P_SFC_LAND, P_SFC_SHIP)
!     continue
!   case(WS_SSMI_FF1, WS_SSMI_FF2)
!     continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(RADIOSONDE_AIR_TEMPERATURE,EVAL_AIR_TEMPERATURE, &
!        T_RAOB,T_PIBAL,T_AIREP,T_AMDAR,T_ACARS,T_MDCRS, &
!        T_TOVS_T, T_TC_SYNTH, T_DROPSONDE)
!     continue
!   case(P_SFC_LAND, P_SFC_SHIP)
!     continue
!   case(WS_SSMI_FF1, WS_SSMI_FF2)
!     continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(RADIOSONDE_AIR_TEMPERATURE,EVAL_AIR_TEMPERATURE, &
!        T_RAOB,T_PIBAL,T_AIREP,T_AMDAR,T_ACARS,T_MDCRS, &
!        T_TOVS_T, T_TC_SYNTH, T_DROPSONDE)
!     continue
!   case(P_SFC_LAND, P_SFC_SHIP)
!     continue
!   case(WS_SSMI_FF1, WS_SSMI_FF2)
!     continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_navdas_mod

  use        types_mod, only : r8, missing_r8, t_kelvin, gas_constant,ps0

  use    utilities_mod, only : register_module

  use     location_mod, only : location_type, vert_is_pressure,&
                               query_location 

  use  assim_model_mod, only : interpolate

  use     obs_kind_mod, only : QTY_POTENTIAL_TEMPERATURE, &
                               QTY_SURFACE_PRESSURE,      &
                               QTY_U_WIND_COMPONENT,      &
                               QTY_V_WIND_COMPONENT,      &
                               QTY_SURFACE_ELEVATION

  use  coamps_intrinsic_mod, only : compute_altimeter

  implicit none
  private

  public :: get_expected_windspeed
  public :: get_expected_air_temperature
  public :: get_expected_altimeter

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  logical, save :: module_initialized = .false.

  ! For some reason specific heat is not included in types_mod
  real(r8), parameter :: Cp = 1004.0_r8

contains

  ! initialize_module
  ! -----------------
  ! Handles any NAVDAS module initialization tasks
  subroutine initialize_module
    call register_module(source, revision, revdate)
    module_initialized = .true.
  end subroutine initialize_module

  ! get_expected_windspeed
  ! ----------------------------
  ! Forward operator for windspeed
  subroutine get_expected_windspeed(state_vector, location, wspd, istatus)  
    real(r8),            intent(in)  :: state_vector(:)
    type(location_type), intent(in)  :: location
    real(r8),            intent(out) :: wspd
    integer,             intent(out) :: istatus

    real(r8) :: uwind   ! zonal wind component
    real(r8) :: vwind   ! meridional wind component

    if ( .not. module_initialized ) call initialize_module

    ! Zonal wind at this location
    call interpolate(state_vector, location, QTY_U_WIND_COMPONENT, uwind, istatus)
    if (istatus /= 0) then
       wspd = missing_r8
       return
    endif

    ! Meridional wind at this location
    call interpolate(state_vector, location, QTY_V_WIND_COMPONENT, vwind, istatus)
    if (istatus /= 0) then
       wspd = missing_r8
       return
    endif

    wspd = sqrt(uwind**2 + vwind**2)
    
  end subroutine get_expected_windspeed

  ! get_expected_air_temperature
  ! ----------------------------
  ! Forward operator for T (in Kelvin) - note that this currently
  ! assumes that the location is 3D in pressure coordinates
  subroutine get_expected_air_temperature(state_vector, location, t, istatus)  
    real(r8),            intent(in)  :: state_vector(:)
    type(location_type), intent(in)  :: location
    real(r8),            intent(out) :: t
    integer,             intent(out) :: istatus

    real(r8) :: pot_t   ! Potential temperature
    real(r8) :: pres    ! Working pressure
    real(r8) :: rocp    ! R over Cp

    if ( .not. module_initialized ) call initialize_module

    ! Potential temperature at this location
    call interpolate(state_vector, location, QTY_POTENTIAL_TEMPERATURE,pot_t, istatus)
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

  ! get_expected_altimeter
  ! ----------------------------
  ! Forward operator for altimeter
  subroutine get_expected_altimeter(state_vector, location, altimeter, istatus)
    real(r8),            intent(in)  :: state_vector(:)
    type(location_type), intent(in)  :: location
    real(r8),            intent(out) :: altimeter
    integer,             intent(out) :: istatus

    real(r8)                         :: sfc_pres
    real(r8)                         :: sfc_elev

    if ( .not. module_initialized ) call initialize_module

    call interpolate(state_vector, location, QTY_SURFACE_PRESSURE, sfc_pres, istatus)
    if (istatus /= 0) then
       altimeter = missing_r8
       return
    endif

    call interpolate(state_vector, location, QTY_SURFACE_ELEVATION, sfc_elev, istatus)
    if (istatus /= 0) then
       altimeter = missing_r8
       return
    endif

    altimeter = compute_altimeter(sfc_pres, sfc_elev)

    if (altimeter < 88000.0_r8 .or. altimeter >= 110000.0_r8) then
       istatus = 1
       altimeter = missing_r8
       return
    endif

    return
  end subroutine get_expected_altimeter

end module obs_def_navdas_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
