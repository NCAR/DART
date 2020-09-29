! Additions by Luke Madaus, Nov 2015, to generate simple
! Surface observations for CM1 OSSE experiments

! BEGIN DART PREPROCESS KIND LIST
!TEMPERATURE_2M,             QTY_2M_TEMPERATURE,       COMMON_CODE
!U_WIND_10M,                 QTY_10M_U_WIND_COMPONENT, COMMON_CODE
!V_WIND_10M,                 QTY_10M_V_WIND_COMPONENT, COMMON_CODE
!SURFACE_PRESSURE,           QTY_SURFACE_PRESSURE,     COMMON_CODE
!SPECIFIC_HUMIDITY_2M,       QTY_SPECIFIC_HUMIDITY,    COMMON_CODE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
! END DART PREPROCESS MODULE CODE





! !!! Note about Specific Humidity observations:
! !!! UNITS in original BUFR are g/kg; This is converted to kg/kg by
! !!! the BUFR to obs_sequence conversion programs making it unnecessary
! !!! to multiply by 1000 at assimilation time.
! !!! PLEASE pay attention to units for specific humidity in models.



! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
