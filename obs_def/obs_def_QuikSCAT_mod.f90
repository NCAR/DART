! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section,
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! <next few lines automatically updated by version control software, do not edit>
! $Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_QuikSCAT_mod.f90,v $
! $URL$
! $Revision$
! $Date$
! $Author$
! $Id$

! This module supports the observation types from the SeaWinds instrument
! on the QuiKSCAT satellite.
! http://winds.jpl.nasa.gov/missions/quikscat/index.cfm
! Since the data already have an NCEP BUFR 'code' of QKSWND,
! we will drag that along.

! While nothing specific needs to be done to use QuikSCAT winds, declaring
! a specific type for them allows for the ability to differentiate them
! from other wind observation types, allowing for impact assessment, for example.

! BEGIN DART PREPROCESS KIND LIST
! QKSWND_U_WIND_COMPONENT,  KIND_U_WIND_COMPONENT
! QKSWND_V_WIND_COMPONENT,  KIND_V_WIND_COMPONENT
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!! No use statements are required for the QuikSCAT scatterometer wind obs_def module
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(QKSWND_U_WIND_COMPONENT)
!            call interpolate(state, location, KIND_U_WIND_COMPONENT, obs_val, istatus)         
!         case(QKSWND_V_WIND_COMPONENT)
!            call interpolate(state, location, KIND_V_WIND_COMPONENT, obs_val, istatus)         
!
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!case(QKSWND_U_WIND_COMPONENT, QKSWND_V_WIND_COMPONENT)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!case(QKSWND_U_WIND_COMPONENT, QKSWND_V_WIND_COMPONENT)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!case(QKSWND_U_WIND_COMPONENT, QKSWND_V_WIND_COMPONENT)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

