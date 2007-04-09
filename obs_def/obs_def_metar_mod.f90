! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! BEGIN DART PREPROCESS KIND LIST
! METAR_U_10_METER_WIND, KIND_U_WIND_COMPONENT
! METAR_V_10_METER_WIND, KIND_V_WIND_COMPONENT
! METAR_TEMPERATURE_2_METER, KIND_TEMPERATURE
! METAR_SPECIFIC_HUMIDITY_2_METER, KIND_SPECIFIC_HUMIDITY
! METAR_SURFACE_PRESSURE, KIND_SURFACE_PRESSURE
! METAR_POT_TEMP_2_METER, KIND_POTENTIAL_TEMPERATURE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the metar obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(METAR_U_10_METER_WIND)
!            call interpolate(state, location, KIND_U_WIND_COMPONENT, obs_val, istatus)
!         case(METAR_V_10_METER_WIND)
!            call interpolate(state, location, KIND_V_WIND_COMPONENT, obs_val, istatus)
!         case(METAR_TEMPERATURE_2_METER)
!            call interpolate(state, location, KIND_TEMPERATURE, obs_val, istatus)
!         case(METAR_SPECIFIC_HUMIDITY_2_METER)
!            call interpolate(state, location, KIND_SPECIFIC_HUMIDITY, obs_val, istatus)
!         case(METAR_SURFACE_PRESSURE)
!            call interpolate(state, location, KIND_SURFACE_PRESSURE, obs_val, istatus)
!         case(METAR_POT_TEMP_2_METER)
!            call interpolate(state, location, KIND_POTENTIAL_TEMPERATURE, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, METAR_TEMPERATURE_2_METER, &
!              METAR_SPECIFIC_HUMIDITY_2_METER, METAR_SURFACE_PRESSURE, &
!              METAR_POT_TEMP_2_METER)
!            continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, METAR_TEMPERATURE_2_METER, &
!              METAR_SPECIFIC_HUMIDITY_2_METER, METAR_SURFACE_PRESSURE, &
!              METAR_POT_TEMP_2_METER)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, METAR_TEMPERATURE_2_METER, &
!              METAR_SPECIFIC_HUMIDITY_2_METER, METAR_SURFACE_PRESSURE, &
!              METAR_POT_TEMP_2_METER)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
