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
! EFGWORO, KIND_GRAV_WAVE_DRAG_EFFIC
! FRACLDV, KIND_GRAV_WAVE_STRESS_FRACTION
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the GWD obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(EFGWORO)
!      ! call get_expected_EFGWORO(state, location, 1, obs_val, istatus)
!   case(FRACLDV)
!      ! call get_expected_FRACLDV(state, location, 2, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(EFGWORO,FRACLDV)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(EFGWORO,FRACLDV)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(EFGWORO,FRACLDV)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
