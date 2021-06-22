! Additions by Luke Madaus, Nov 2015, to generate simple
! Surface observations for CM1 OSSE experiments

! !!! Note about Specific Humidity observations:
! !!! UNITS in original BUFR are g/kg; This is converted to kg/kg by
! !!! the BUFR to obs_sequence conversion programs making it unnecessary
! !!! to multiply by 1000 at assimilation time.
! !!! PLEASE pay attention to units for specific humidity in models.

! BEGIN DART PREPROCESS TYPE DEFINITIONS
!TEMPERATURE_2M,             QTY_2M_TEMPERATURE,       COMMON_CODE
!U_WIND_10M,                 QTY_10M_U_WIND_COMPONENT, COMMON_CODE
!V_WIND_10M,                 QTY_10M_V_WIND_COMPONENT, COMMON_CODE
!SURFACE_PRESSURE,           QTY_SURFACE_PRESSURE,     COMMON_CODE
!SPECIFIC_HUMIDITY_2M,       QTY_SPECIFIC_HUMIDITY,    COMMON_CODE
! END DART PREPROCESS TYPE DEFINITIONS



