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
!U_WIND_COMPONENT,    KIND_U_WIND_COMPONENT
!V_WIND_COMPONENT,    KIND_V_WIND_COMPONENT
!GEOPOTENTIAL_HEIGHT, KIND_GEOPOTENTIAL_HEIGHT
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the pe2lyr obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(U_WIND_COMPONENT)
!            call interpolate(state, location, KIND_U_WIND_COMPONENT, obs_val, istatus)         
!         case(V_WIND_COMPONENT)
!            call interpolate(state, location, KIND_V_WIND_COMPONENT, obs_val, istatus)         
!         case(GEOPOTENTIAL_HEIGHT)
!            call interpolate(state, location, KIND_GEOPOTENTIAL_HEIGHT, obs_val, istatus)
!
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!case(U_WIND_COMPONENT, V_WIND_COMPONENT, GEOPOTENTIAL_HEIGHT)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!case(U_WIND_COMPONENT, V_WIND_COMPONENT, GEOPOTENTIAL_HEIGHT)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!case(U_WIND_COMPONENT, V_WIND_COMPONENT, GEOPOTENTIAL_HEIGHT)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

