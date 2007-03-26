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
!EVAL_U_WIND_COMPONENT,  KIND_U_WIND_COMPONENT
!EVAL_V_WIND_COMPONENT,  KIND_V_WIND_COMPONENT
!EVAL_SURFACE_PRESSURE,  KIND_SURFACE_PRESSURE
!EVAL_TEMPERATURE,       KIND_TEMPERATURE
!EVAL_SPECIFIC_HUMIDITY, KIND_SPECIFIC_HUMIDITY
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the reanalysis bufr obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(EVAL_U_WIND_COMPONENT)
!            call interpolate(state, location, KIND_U_WIND_COMPONENT, obs_val, istatus)         
!         case(EVAL_V_WIND_COMPONENT)
!            call interpolate(state, location, KIND_V_WIND_COMPONENT, obs_val, istatus)         
!         case(EVAL_TEMPERATURE)
!            call interpolate(state, location, KIND_TEMPERATURE, obs_val, istatus)
!
!         case(EVAL_SPECIFIC_HUMIDITY)
!            call interpolate(state, location, KIND_SPECIFIC_HUMIDITY, obs_val, istatus)
!         case(EVAL_SURFACE_PRESSURE)
!            call interpolate(state, location, KIND_SURFACE_PRESSURE, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!case(EVAL_U_WIND_COMPONENT, EVAL_V_WIND_COMPONENT, EVAL_SURFACE_PRESSURE, &
!   EVAL_TEMPERATURE, EVAL_SPECIFIC_HUMIDITY)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!case(EVAL_U_WIND_COMPONENT, EVAL_V_WIND_COMPONENT, EVAL_SURFACE_PRESSURE, &
!   EVAL_TEMPERATURE, EVAL_SPECIFIC_HUMIDITY)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!case(EVAL_U_WIND_COMPONENT, EVAL_V_WIND_COMPONENT, EVAL_SURFACE_PRESSURE, &
!   EVAL_TEMPERATURE, EVAL_SPECIFIC_HUMIDITY)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

