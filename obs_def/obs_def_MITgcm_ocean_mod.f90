! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! <next few lines under version control, do not edit>
! $URL$
! $Id: $
! $Revision$
! $Date: $

! BEGIN DART PREPROCESS KIND LIST
!TEMPERATURE,            KIND_TEMPERATURE
!SALINITY,               KIND_SALINITY
!U_CURRENT_COMPONENT,    KIND_U_WIND_COMPONENT
!V_CURRENT_COMPONENT,    KIND_V_WIND_COMPONENT
!SEA_SURFACE_HEIGHT,     KIND_SEA_SURFACE_HEIGHT
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the MITgcm_ocean obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(TEMPERATURE)
!            call interpolate(state, location, KIND_TEMPERATURE, obs_val, istatus)         
!         case(SALINITY)
!            call interpolate(state, location, KIND_SALINITY, obs_val, istatus)         
!         case(U_CURRENT_COMPONENT)
!            call interpolate(state, location, KIND_U_CURRENT_COMPONENT, obs_val, istatus)         
!         case(V_CURRENT_COMPONENT)
!            call interpolate(state, location, KIND_V_CURRENT_COMPONENT, obs_val, istatus)         
!         case(SEA_SURFACE_HEIGHT)
!            call interpolate(state, location, KIND_SEA_SURFACE_HEIGHT, obs_val, istatus)
!
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!case(TEMPERATURE, SALINITY,  U_CURRENT_COMPONENT, V_CURRENT_COMPONENT, SEA_SURFACE_HEIGHT)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!case(TEMPERATURE, SALINITY,  U_CURRENT_COMPONENT, V_CURRENT_COMPONENT, SEA_SURFACE_HEIGHT)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!case(TEMPERATURE, SALINITY,  U_CURRENT_COMPONENT, V_CURRENT_COMPONENT, SEA_SURFACE_HEIGHT)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

