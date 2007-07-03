! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/obs_def/obs_def_simple_advection_mod.f90 $
! $Id: obs_def_simple_advection_mod.f90 2713 2007-03-26 04:09:04Z thoar $
! $Revision: 2713 $
! $Date: 2007-03-25 22:09:04 -0600 (Sun, 25 Mar 2007) $

! BEGIN DART PREPROCESS KIND LIST
!VELOCITY,                     KIND_VELOCITY
!TRACER_CONCENTRATION,         KIND_TRACER_CONCENTRATION
!TRACER_SOURCE,                KIND_TRACER_SOURCE
!MEAN_SOURCE,                  KIND_MEAN_SOURCE
!SOURCE_PHASE,                 KIND_SOURCE_PHASE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the simple_advection obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(VELOCITY)
!            call interpolate(state, location, KIND_VELOCITY, obs_val, istatus)         
!         case(TRACER_CONCENTRATION)
!            call interpolate(state, location, KIND_TRACER_CONCENTRATION, obs_val, istatus)         
!         case(TRACER_SOURCE)
!            call interpolate(state, location, KIND_TRACER_SOURCE, obs_val, istatus)
!         case(MEAN_SOURCE)
!            call interpolate(state, location, KIND_MEAN_SOURCE, obs_val, istatus)
!         case(SOURCE_PHASE)
!            call interpolate(state, location, KIND_SOURCE_PHASE, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!case(VELOCITY, TRACER_CONCENTRATION, TRACER_SOURCE, MEAN_SOURCE, SOURCE_PHASE)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!case(VELOCITY, TRACER_CONCENTRATION, TRACER_SOURCE, MEAN_SOURCE, SOURCE_PHASE)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!case(VELOCITY, TRACER_CONCENTRATION, TRACER_SOURCE, MEAN_SOURCE, SOURCE_PHASE)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
