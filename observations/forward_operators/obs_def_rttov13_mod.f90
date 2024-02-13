! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!----------------------------------------------------------------------
! This module provides support for computing forward operators for
! radiance observations by calling the rttov radiative transfer model.
!
! the additional metadata in each observation includes:
!
!  OBS            X
! rttov
!    sat az/el
!    sun az/el
!    platform
!    sat_id
!    instrument
!    channel
!    <anything else useful>
!
! The following is a list of all the radiances types DART can assimilate
! using RTTOV with the following convention:
!
! (PLATFORM)_(SATELLITE)_(SENSOR)_RADIANCE, QTY_RADIANCE
! 
! where   PLATFORM    is the satellite series (e.g. NOAA or DMSP),
!         SATELLITE   is the satellite number (e.g. 18 for NOAA 18),
!         SENSOR      is the satellite sensor name (e.g. AIRS for Aqua/EOS 2)
!
! This list is sorted in numerical order of RTTOV platform, satellite, 
! sensor ID. The list is sorted first on platform id, then in case of
! duplicates on satellite id, then in case of duplicates on sensor ID.
! Adding a new sensor should follow this convention.
!
! Below the list, the first and last radiance will be used in the select 
! case statements as follows:
!
!   case(NOAA_1_VTPR1_RADIANCE:CLOUDSAT_1_CPR_TB)
!
! This ensures that all radiance types between these two numbers will be
! handled the same way. This is a Fortran construct for a case statement 
! between two numbers, which both must be parameters. Preprocess ensures 
! that all types listed in descending order in observation def files are 
! given consecutive, ascending parameter values. 
!
! N.B. that if another observation type is added before or after the start/end
! (currently NOAA_1_VTPR1_RADIANCE and CLOUDSAT_1_CPR_TB, respectively), 
! all case statements in this obs def file must be modified appropriately to 
! ensure correct processing.  
!
! The units for the microwave brightness temperatures (any observations 
! that have a quantity of QTY_BRIGHTNESS_TEMPERATURE) is degrees kelvin.
!
! The units for all observations with a quantity of QTY_RADANCE are as
! described in the RTTOV v12 user guide V1.3 (p54): "mW/cm-1/sr/sq.m"
! Doc ID: NWPSAF-MO-UD-037, Date: 05/03/2019
!
! The observation converters are responsible for providing these 
! observations with the correct units.
!----------------------------------------------------------------------

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! NOAA_1_VTPR1_RADIANCE,        QTY_RADIANCE
! NOAA_2_VTPR1_RADIANCE,        QTY_RADIANCE
! NOAA_3_VTPR1_RADIANCE,        QTY_RADIANCE
! NOAA_4_VTPR1_RADIANCE,        QTY_RADIANCE
! NOAA_5_HIRS_RADIANCE,         QTY_RADIANCE
! NOAA_5_MSU_TB,                QTY_BRIGHTNESS_TEMPERATURE
! NOAA_5_AVHRR_RADIANCE,        QTY_RADIANCE
! NOAA_6_HIRS_RADIANCE,         QTY_RADIANCE
! NOAA_6_MSU_TB,                QTY_BRIGHTNESS_TEMPERATURE
! NOAA_6_AVHRR_RADIANCE,        QTY_RADIANCE
! NOAA_7_HIRS_RADIANCE,         QTY_RADIANCE
! NOAA_7_MSU_TB,                QTY_BRIGHTNESS_TEMPERATURE
! NOAA_7_AVHRR_RADIANCE,        QTY_RADIANCE
! NOAA_8_HIRS_RADIANCE,         QTY_RADIANCE
! NOAA_8_MSU_TB,                QTY_BRIGHTNESS_TEMPERATURE
! NOAA_8_AVHRR_RADIANCE,        QTY_RADIANCE
! NOAA_9_HIRS_RADIANCE,         QTY_RADIANCE
! NOAA_9_MSU_TB,                QTY_BRIGHTNESS_TEMPERATURE
! NOAA_9_AVHRR_RADIANCE,        QTY_RADIANCE
! NOAA_10_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_10_MSU_TB,               QTY_BRIGHTNESS_TEMPERATURE
! NOAA_10_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_11_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_11_MSU_TB,               QTY_BRIGHTNESS_TEMPERATURE
! NOAA_11_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_12_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_12_MSU_TB,               QTY_BRIGHTNESS_TEMPERATURE
! NOAA_12_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_13_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_14_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_14_MSU_TB,               QTY_BRIGHTNESS_TEMPERATURE
! NOAA_14_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_15_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_15_AMSUA_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_15_AMSUB_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_15_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_16_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_16_AMSUA_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_16_AMSUB_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_16_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_17_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_17_AMSUA_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_17_AMSUB_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_17_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_18_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_18_AMSUA_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_18_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_18_MHS_TB,               QTY_BRIGHTNESS_TEMPERATURE
! NOAA_19_HIRS_RADIANCE,        QTY_RADIANCE
! NOAA_19_AMSUA_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NOAA_19_AVHRR_RADIANCE,       QTY_RADIANCE
! NOAA_19_MHS_TB,               QTY_BRIGHTNESS_TEMPERATURE
! NOAA_20_ATMS_TB,              QTY_BRIGHTNESS_TEMPERATURE
! NOAA_20_VIIRS_RADIANCE,       QTY_RADIANCE
! DMSP_8_SSMI_TB,               QTY_BRIGHTNESS_TEMPERATURE
! DMSP_9_SSMI_TB,               QTY_BRIGHTNESS_TEMPERATURE
! DMSP_10_SSMI_TB,              QTY_BRIGHTNESS_TEMPERATURE
! DMSP_11_SSMI_TB,              QTY_BRIGHTNESS_TEMPERATURE
! DMSP_11_SSMT2_TB,             QTY_BRIGHTNESS_TEMPERATURE
! DMSP_12_SSMI_TB,              QTY_BRIGHTNESS_TEMPERATURE
! DMSP_12_SSMT2_TB,             QTY_BRIGHTNESS_TEMPERATURE
! DMSP_13_SSMI_TB,              QTY_BRIGHTNESS_TEMPERATURE
! DMSP_14_SSMI_TB,              QTY_BRIGHTNESS_TEMPERATURE
! DMSP_14_SSMT2_TB,             QTY_BRIGHTNESS_TEMPERATURE
! DMSP_15_SSMI_TB,              QTY_BRIGHTNESS_TEMPERATURE
! DMSP_15_SSMT2_TB,             QTY_BRIGHTNESS_TEMPERATURE
! DMSP_16_SSMIS_TB,             QTY_BRIGHTNESS_TEMPERATURE
! DMSP_17_SSMIS_TB,             QTY_BRIGHTNESS_TEMPERATURE
! DMSP_18_SSMIS_TB,             QTY_BRIGHTNESS_TEMPERATURE
! DMSP_19_SSMIS_TB,             QTY_BRIGHTNESS_TEMPERATURE
! METEOSAT_1_MVIRI_RADIANCE,    QTY_RADIANCE
! METEOSAT_2_MVIRI_RADIANCE,    QTY_RADIANCE
! METEOSAT_3_MVIRI_RADIANCE,    QTY_RADIANCE
! METEOSAT_4_MVIRI_RADIANCE,    QTY_RADIANCE
! METEOSAT_5_MVIRI_RADIANCE,    QTY_RADIANCE
! METEOSAT_6_MVIRI_RADIANCE,    QTY_RADIANCE
! METEOSAT_7_MVIRI_RADIANCE,    QTY_RADIANCE
! GOES_4_SOUNDER_RADIANCE,      QTY_RADIANCE
! GOES_5_SOUNDER_RADIANCE,      QTY_RADIANCE
! GOES_6_SOUNDER_RADIANCE,      QTY_RADIANCE
! GOES_7_SOUNDER_RADIANCE,      QTY_RADIANCE
! GOES_8_IMAGER_RADIANCE,       QTY_RADIANCE
! GOES_8_SOUNDER_RADIANCE,      QTY_RADIANCE
! GOES_9_IMAGER_RADIANCE,       QTY_RADIANCE
! GOES_9_SOUNDER_RADIANCE,      QTY_RADIANCE
! GOES_10_IMAGER_RADIANCE,      QTY_RADIANCE
! GOES_10_SOUNDER_RADIANCE,     QTY_RADIANCE
! GOES_11_IMAGER_RADIANCE,      QTY_RADIANCE
! GOES_11_SOUNDER_RADIANCE,     QTY_RADIANCE
! GOES_12_IMAGER_RADIANCE,      QTY_RADIANCE
! GOES_12_SOUNDER_RADIANCE,     QTY_RADIANCE
! GOES_13_IMAGER_RADIANCE,      QTY_RADIANCE
! GOES_13_SOUNDER_RADIANCE,     QTY_RADIANCE
! GOES_14_IMAGER_RADIANCE,      QTY_RADIANCE
! GOES_14_SOUNDER_RADIANCE,     QTY_RADIANCE
! GOES_15_IMAGER_RADIANCE,      QTY_RADIANCE
! GOES_15_SOUNDER_RADIANCE,     QTY_RADIANCE
! GOES_16_ABI_RADIANCE,         QTY_RADIANCE
! GOES_17_ABI_RADIANCE,         QTY_RADIANCE
! GOES_18_ABI_RADIANCE,         QTY_RADIANCE
! GOES_19_ABI_RADIANCE,         QTY_RADIANCE
! GMS_1_IMAGER_RADIANCE,        QTY_RADIANCE
! GMS_2_IMAGER_RADIANCE,        QTY_RADIANCE
! GMS_3_IMAGER_RADIANCE,        QTY_RADIANCE
! GMS_4_IMAGER_RADIANCE,        QTY_RADIANCE
! GMS_5_IMAGER_RADIANCE,        QTY_RADIANCE
! FY2_2_VISSR_RADIANCE,         QTY_RADIANCE
! FY2_3_VISSR_RADIANCE,         QTY_RADIANCE
! FY2_4_VISSR_RADIANCE,         QTY_RADIANCE
! FY2_5_VISSR_RADIANCE,         QTY_RADIANCE
! FY2_7_VISSR_RADIANCE,         QTY_RADIANCE
! TRMM_1_TMI_TB,                QTY_BRIGHTNESS_TEMPERATURE
! ERS_1_ATSR_RADIANCE,          QTY_RADIANCE
! ERS_1_MWR_TB,                 QTY_BRIGHTNESS_TEMPERATURE
! ERS_2_ATSR_RADIANCE,          QTY_RADIANCE
! ERS_2_MWR_TB,                 QTY_BRIGHTNESS_TEMPERATURE
! EOS_1_MODIS_RADIANCE,         QTY_RADIANCE
! EOS_1_ASTER_RADIANCE,         QTY_RADIANCE
! EOS_2_AMSUA_TB,               QTY_BRIGHTNESS_TEMPERATURE
! EOS_2_AIRS_RADIANCE,          QTY_RADIANCE
! EOS_2_HSB_TB,                 QTY_BRIGHTNESS_TEMPERATURE
! EOS_2_MODIS_RADIANCE,         QTY_RADIANCE
! EOS_2_AMSRE_TB,               QTY_BRIGHTNESS_TEMPERATURE
! METOP_1_HIRS_RADIANCE,        QTY_RADIANCE
! METOP_1_AMSUA_TB,             QTY_BRIGHTNESS_TEMPERATURE
! METOP_1_AVHRR_RADIANCE,       QTY_RADIANCE
! METOP_1_MHS_TB,               QTY_BRIGHTNESS_TEMPERATURE
! METOP_2_HIRS_RADIANCE,        QTY_RADIANCE
! METOP_2_AMSUA_TB,             QTY_BRIGHTNESS_TEMPERATURE
! METOP_2_AVHRR_RADIANCE,       QTY_RADIANCE
! METOP_2_MHS_TB,               QTY_BRIGHTNESS_TEMPERATURE
! METOP_3_AVHRR_RADIANCE,       QTY_RADIANCE
! ENVISAT_1_ATSR_RADIANCE,      QTY_RADIANCE
! ENVISAT_1_MWR_TB,             QTY_BRIGHTNESS_TEMPERATURE
! MSG_1_SEVIRI_RADIANCE,        QTY_RADIANCE
! MSG_2_SEVIRI_RADIANCE,        QTY_RADIANCE
! MSG_3_SEVIRI_RADIANCE,        QTY_RADIANCE
! MSG_4_SEVIRI_RADIANCE,        QTY_RADIANCE
! MSG_4_SEVIRI_TB,              QTY_BRIGHTNESS_TEMPERATURE
! MSG_4_SEVIRI_BDRF,            QTY_BI_DIRECTIONAL_REFLECTANCE
! FY1_3_MVISR_RADIANCE,         QTY_RADIANCE
! FY1_4_MVISR_RADIANCE,         QTY_RADIANCE
! MTSAT_1_IMAGER_RADIANCE,      QTY_RADIANCE
! MTSAT_2_IMAGER_RADIANCE,      QTY_RADIANCE
! CORIOLIS_1_WINDSAT_TB,        QTY_BRIGHTNESS_TEMPERATURE
! JPSS_0_ATMS_TB,               QTY_BRIGHTNESS_TEMPERATURE
! JPSS_0_VIIRS_RADIANCE,        QTY_RADIANCE
! SENTINEL3_1_SLSTR_RADIANCE,   QTY_RADIANCE
! SENTINEL3_2_SLSTR_RADIANCE,   QTY_RADIANCE
! MEGHATR_1_SAPHIR_TB,          QTY_BRIGHTNESS_TEMPERATURE
! MEGHATR_1_MADRAS_TB,          QTY_BRIGHTNESS_TEMPERATURE
! FY3_1_MWTS_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_1_MWHS_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_1_IRAS_RADIANCE,          QTY_RADIANCE
! FY3_1_MWRI_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_2_MWTS_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_2_MWHS_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_2_MWRI_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_3_MWRI_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_3_MWTS2_TB,               QTY_BRIGHTNESS_TEMPERATURE
! FY3_3_MWHS2_TB,               QTY_BRIGHTNESS_TEMPERATURE
! FY3_3_MERSI1_RADIANCE,        QTY_RADIANCE
! FY3_4_MWRI_TB,                QTY_BRIGHTNESS_TEMPERATURE
! FY3_4_MWTS2_TB,               QTY_BRIGHTNESS_TEMPERATURE
! FY3_4_MWHS2_TB,               QTY_BRIGHTNESS_TEMPERATURE
! FY3_4_MERSI2_RADIANCE,        QTY_RADIANCE
! COMS_1_MI_RADIANCE,           QTY_RADIANCE
! METEOR_M_1_MSUMR_RADIANCE,    QTY_RADIANCE
! METEOR_M_2_MSUMR_RADIANCE,    QTY_RADIANCE
! METEOR_M_2_MTVZAGY_TB,        QTY_BRIGHTNESS_TEMPERATURE
! CALIPSO_1_IIR_RADIANCE,       QTY_RADIANCE
! GCOM_W_1_AMSR2_TB,            QTY_BRIGHTNESS_TEMPERATURE
! NIMBUS_3_MRIR_RADIANCE,       QTY_RADIANCE
! NIMBUS_4_THIR_RADIANCE,       QTY_RADIANCE
! NIMBUS_5_THIR_RADIANCE,       QTY_RADIANCE
! NIMBUS_6_HIRS_RADIANCE,       QTY_RADIANCE
! NIMBUS_6_SCAMS_TB,            QTY_BRIGHTNESS_TEMPERATURE
! NIMBUS_6_THIR_RADIANCE,       QTY_RADIANCE
! NIMBUS_7_SMMR_TB,             QTY_BRIGHTNESS_TEMPERATURE
! NIMBUS_7_THIR_RADIANCE,       QTY_RADIANCE
! HIMAWARI_8_AHI_RADIANCE,      QTY_RADIANCE
! HIMAWARI_9_AHI_RADIANCE,      QTY_RADIANCE
! MTG_1_FCI_RADIANCE,           QTY_RADIANCE
! SARAL_1_ALTIKA_TB,            QTY_BRIGHTNESS_TEMPERATURE
! METOPSG_1_ICI_TB,             QTY_BRIGHTNESS_TEMPERATURE
! METOPSG_1_METIMAGE_RADIANCE,  QTY_RADIANCE
! METOPSG_1_MWS_TB,             QTY_BRIGHTNESS_TEMPERATURE
! METOPSG_1_MWI_TB,             QTY_BRIGHTNESS_TEMPERATURE
! LANDSAT_4_TM_RADIANCE,        QTY_RADIANCE
! LANDSAT_5_TM_RADIANCE,        QTY_RADIANCE
! LANDSAT_7_TM_RADIANCE,        QTY_RADIANCE
! LANDSAT_8_TIRS_RADIANCE,      QTY_RADIANCE
! JASON_2_AMR_TB,               QTY_BRIGHTNESS_TEMPERATURE
! GPM_1_GMI_TB,                 QTY_BRIGHTNESS_TEMPERATURE
! GPM_1_DPR_TB,                 QTY_BRIGHTNESS_TEMPERATURE
! INSAT3_4_IMAGER_RADIANCE,     QTY_RADIANCE
! INSAT3_4_SOUNDER_RADIANCE,    QTY_RADIANCE
! INSAT3_5_IMAGER_RADIANCE,     QTY_RADIANCE
! INSAT3_5_SOUNDER_RADIANCE,    QTY_RADIANCE
! TICFIRE_1_MBFIRI_RADIANCE,    QTY_RADIANCE
! ISS_1_ECOSTRES_RADIANCE,      QTY_RADIANCE
! HJ1_2_IRMSS_RADIANCE,         QTY_RADIANCE
! GKOMPSAT2_1_AMI_RADIANCE,     QTY_RADIANCE
! GCOM_C_1_SGLI_RADIANCE,       QTY_RADIANCE
! SMOS_1_MIRAS_TB,              QTY_BRIGHTNESS_TEMPERATURE
! ORS_6_COWVR_TB,               QTY_BRIGHTNESS_TEMPERATURE
! FY4_1_AGRI_RADIANCE,          QTY_RADIANCE
! TROPICS_0_TROPICS_TB,         QTY_BRIGHTNESS_TEMPERATURE
! GF5_1_VIMS_RADIANCE,          QTY_RADIANCE
! HY2_1_MWRI_TB,                QTY_BRIGHTNESS_TEMPERATURE
! CLOUDSAT_1_CPR_TB,            QTY_BRIGHTNESS_TEMPERATURE 
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_rttov_mod, only : read_rttov_metadata, &
!                                write_rttov_metadata, &
!                          interactive_rttov_metadata, &
!                                get_expected_radiance
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(NOAA_1_VTPR1_RADIANCE:CLOUDSAT_1_CPR_TB)
!         call get_expected_radiance(obs_kind_ind, state_handle, ens_size, location, obs_def%key, get_obs_def_type_of_obs(obs_def), expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(NOAA_1_VTPR1_RADIANCE:CLOUDSAT_1_CPR_TB)
!      call read_rttov_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(NOAA_1_VTPR1_RADIANCE:CLOUDSAT_1_CPR_TB)
!      call write_rttov_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(NOAA_1_VTPR1_RADIANCE:CLOUDSAT_1_CPR_TB)
!      call interactive_rttov_metadata(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_rttov_mod

use        types_mod, only : r8, MISSING_R8, MISSING_I, obstypelength

use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, E_ALLMSG, &
                             ascii_file_format, nmlfileunit, do_nml_file, &
                             do_nml_term, check_namelist_read, find_namelist_in_file, &
                             interactive_r, interactive_i, open_file, file_exist, &
                             close_file

use     location_mod, only : location_type, set_location, get_location, &
                             VERTISUNDEF, VERTISHEIGHT, VERTISLEVEL, VERTISSURFACE

use  assim_model_mod, only : interpolate

use obs_def_utilities_mod, only : track_status
use  ensemble_manager_mod, only : ensemble_type
use         utilities_mod, only : to_upper

! This code contains the DART to RTTOV interface. It uses a data structure of
! platforms/satellites/sensors to store runtime information for each sensor.
! This runtime information is initialized on the fly as needed to compute the
! RTTOV forward operator. 
!
! Copyright:
!    RTTOV was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    RTTOV is copyright 2017, EUMETSAT, All Rights Reserved.

! RTTOV module containing useful RTTOV constants
use rttov_const, only :     &
      errorstatus_success, &
      errorstatus_fatal,   &
      platform_name,       &
      inst_name,           &
      pmax,                &
      pmin

! RTTOV module rttov_types contains definitions of all RTTOV data types
use rttov_types, only :     &
      rttov_options,       &
      rttov_options_scatt, &
      rttov_scatt_coef,    &
      rttov_coefs,         &
      rttov_profile,       &
      rttov_profile_cloud, &
      rttov_transmission,  &
      rttov_radiance,      &
      rttov_chanprof,      &
      rttov_emissivity,    &
      rttov_reflectance

! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
use parkind1, only : jpim, jprb, jplm


! There are so many radiance observation types that it is impractical to list them all.
use     obs_kind_mod

implicit none

private

! FIXME: have a separate module only for the metadata
! routines that will be "use"d by this module and the
! observation converter.

!FIXME: the converters only need the metadata routines.
! filter needs them all.  but with a complicated forward
! operator, the convert now needs to be built with the
! entire rttov library even though it will never call it.
! how can we exclude the get_expected code from the
! converter?  preprocessor options?  split this file
! into two separate ones?  forward operator vs obs metadata?  
! then the converters just use the metadata and filter
! uses both...  or something.

public ::         set_visir_metadata, &
                     set_mw_metadata, &
                  get_visir_metadata, & 
                     get_mw_metadata, &
                 read_rttov_metadata, &
                write_rttov_metadata, &
          interactive_rttov_metadata, &
               get_expected_radiance, &
            get_rttov_option_logical, &
                         get_channel

! The rttov_test.f90 program uses these, but no one else should.

public ::     atmos_profile_type, &
          trace_gas_profile_type, &
              cloud_profile_type, &
            aerosol_profile_type, &
               rttov_sensor_type, &
             visir_metadata_type, &
                mw_metadata_type, &
            sensor_runtime_setup, &
             atmos_profile_setup, &
         trace_gas_profile_setup, &
           aerosol_profile_setup, &
             cloud_profile_setup, &
             read_sensor_db_file, &
                get_rttov_sensor, &
                do_forward_model

! routines for unit testing
public :: test_unit_setup,       &
          test_set_metadata,     &
          test_unit_teardown,    &
          test_metadata,         &
          test_key_within_range, &
          test_subkey_within_range

! Metadata for rttov observations.
type visir_metadata_type
   real(jprb)    :: sat_az      ! azimuth of satellite position (degrees)
   real(jprb)    :: sat_ze      ! zenith of satellite position (degrees)
   real(jprb)    :: sun_az      ! azimuth of solar position (degrees, only used with add_solar)
   real(jprb)    :: sun_ze      ! zenith of solar position  (degrees, only used with add_solar)
   integer(jpim) :: platform_id ! see rttov user guide, table 2
   integer(jpim) :: sat_id      ! see rttov user guide, table 2
   integer(jpim) :: sensor_id   ! see rttov user guide, table 3
   integer(jpim) :: channel     ! each channel is a different obs
   real(jprb)    :: specularity ! specularity (0-1, only used with do_lambertian)
end type visir_metadata_type

type mw_metadata_type
   real(jprb)    :: sat_az      ! azimuth of satellite position (degrees)
   real(jprb)    :: sat_ze      ! zenith of satellite position (degrees)
   integer(jpim) :: platform_id ! see rttov user guide, table 2
   integer(jpim) :: sat_id      ! see rttov user guide, table 2
   integer(jpim) :: sensor_id   ! see rttov user guide, table 3
   integer(jpim) :: channel     !  each channel is a different obs
   real(jprb)    :: mag_field   ! strength of mag_field (Gauss, )
   real(jprb)    :: cosbk       ! cosine of angle between mag field and viewing angle
   real(jprb)    :: fastem_p1 ! FASTEM land/sea ice parameter 1
   real(jprb)    :: fastem_p2 ! FASTEM land/sea ice parameter 2
   real(jprb)    :: fastem_p3 ! FASTEM land/sea ice parameter 3
   real(jprb)    :: fastem_p4 ! FASTEM land/sea ice parameter 4
   real(jprb)    :: fastem_p5 ! FASTEM land/sea ice parameter 5
end type mw_metadata_type

! DART container type to hold the essential atmosphere and surface fields.
! For 2D fields, the order is (ens_size, numlevels), while 1D fields are 
! (ens_size).
type atmos_profile_type
   real(r8), allocatable :: temperature(:,:)    ! mandatory, level temperature (K)
   real(r8), allocatable :: pressure(:,:)       ! mandatory, level pressure (Pa)
   real(r8), allocatable :: moisture(:,:)       ! mandatory, level water vapor (kg/kg)
   real(r8), allocatable :: sfc_p(:)            ! mandatory, surface pressure (Pa)
   real(r8), allocatable :: s2m_t(:)            ! mandatory, 2 meter temp (K)
   real(r8), allocatable :: skin_temp(:)        ! mandatory, surface skin temp (K)
   real(r8), allocatable :: sfc_elev(:)         ! mandatory, surface elevation (m)
   real(r8), allocatable :: surftype(:)         ! mandatory, surface type (land=0, water=1, seaice = 2) 
   real(r8), allocatable :: s2m_q(:)            ! optional, 2 meter wator vapor (kg/kg) (used if add_q2m)
   real(r8), allocatable :: s10m_u(:)           ! optional, 10 meter u wind (m/s) (used if add_uv10m)
   real(r8), allocatable :: s10m_v(:)           ! optional, 10 meter v wind (m/s) (used if add_uv10m)
   real(r8), allocatable :: wfetch(:)           ! optional, wind fetch (m) (used if use_wfetch)
   real(r8), allocatable :: water_type(:)       ! optional, water type (fresh=0, ocean=1) (used if use_water_type)
   real(r8), allocatable :: sfc_salinity(:)     ! optional, ocean salinity (practial salinity unit) (used if use_salinity)
   real(r8), allocatable :: sfc_foam_frac(:)    ! optional, foam fraction (0-1) (used if supply_foam_fraction)
   real(r8), allocatable :: sfc_snow_frac(:)    ! optional, snow cover (0-1) (used if use_sfc_snow_frac)
end type atmos_profile_type

! container type for trace gasses. Note these values are on levels, so the 
! size of the arrays are (ens_size, numlevels)
type trace_gas_profile_type
   real(r8), allocatable :: ozone(:,:)         ! ozone concentration (kg/kg) (used if add_ozone)
   real(r8), allocatable :: co2(:,:)           ! CO2 concentration   (kg/kg) (used if add_co2)
   real(r8), allocatable :: n2o(:,:)           ! N2O concentration   (kg/kg) (used if add_n2o)
   real(r8), allocatable :: ch4(:,:)           ! CH4 concentration   (kg/kg) (used if add_ch4)
   real(r8), allocatable :: co(:,:)            ! CO concentration    (kg/kg) (used if add_co)
   real(r8), allocatable :: so2(:,:)           ! SO2 concentration   (kg/kg) (used if add_so2)
end type trace_gas_profile_type

! Container type for aerosols. Note these values are on layers, so the 
! size of the arrays are (ens_size, numlevels-1)
type aerosol_profile_type
   real(r8), allocatable :: insoluble(:,:)                  ! INSO (kg/kg), OPAC only
   real(r8), allocatable :: water_soluble(:,:)              ! WASO, OPAC
   real(r8), allocatable :: soot(:,:)                       ! SOOT, OPAC
   real(r8), allocatable :: sea_salt_accum(:,:)             ! SSAM, OPAC
   real(r8), allocatable :: sea_salt_coarse(:,:)            ! SSCM, OPAC
   real(r8), allocatable :: mineral_nucleus(:,:)            ! MINM, OPAC
   real(r8), allocatable :: mineral_accum(:,:)              ! MIAM, OPAC
   real(r8), allocatable :: mineral_coarse(:,:)             ! MICM, OPAC
   real(r8), allocatable :: mineral_transport(:,:)          ! MITR, OPAC
   real(r8), allocatable :: sulphated_droplets(:,:)         ! SUSO, OPAC
   real(r8), allocatable :: volcanic_ash(:,:)               ! VOLA, OPAC
   real(r8), allocatable :: new_volcanic_ash(:,:)           ! VAPO, OPAC
   real(r8), allocatable :: asian_dust(:,:)                 ! ASDU, OPAC
   real(r8), allocatable :: black_carbon(:,:)               ! BCAR, CAMS
   real(r8), allocatable :: dust_bin1(:,:)                  ! DUS1, CAMS
   real(r8), allocatable :: dust_bin2(:,:)                  ! DUS2, CAMS
   real(r8), allocatable :: dust_bin3(:,:)                  ! DUS3, CAMS
   real(r8), allocatable :: ammonium_sulphate(:,:)          ! SULP, CAMS
   real(r8), allocatable :: sea_salt_bin1(:,:)              ! SSA1, CAMS
   real(r8), allocatable :: sea_salt_bin2(:,:)              ! SSA2, CAMS
   real(r8), allocatable :: sea_salt_bin3(:,:)              ! SSA3, CAMS
   real(r8), allocatable :: hydrophilic_organic_matter(:,:) ! OMAT, CAMS
end type aerosol_profile_type

! container type for clouds - note RTTOV uses different cloud fields depending on scheme, frequency, type, etc. Note these values are on layers, not levels, so the 
! size of the arrays are (ens_size, numlevels-1)
type cloud_profile_type
   real(r8), allocatable :: hydro(:,:,:)         ! (MW) hydrometeor concentrations, dimension (imem, nlevels, nhydro)
   real(r8), allocatable :: hydro_frac(:)    ! (MW) cloud fraction (0-1) dimension (imem, nlevels, nhydro_frac=1)
   real(r8), allocatable :: cfrac(:,:)         ! (VIS/IR) cloud fractional cover (0-1) (imem, nlevels)
   real(r8), allocatable :: simple_cfrac(:)    ! (VIS/IR)    cloud fraction for simple cloud (0-1) 
   real(r8), allocatable :: ctp(:)             ! (VIS/IR)    cloud top pressure for simple cloud (hPa) 
   real(r8), allocatable :: w(:,:)             ! (VIS/IR)    vertical velocity (used for classification of cumulus vs. stratus)
   real(r8), allocatable :: clw(:,:)           ! (MW) cloud non-precipitating water; "clear air MW radiance", outside RTTOV-SCATT
   real(r8), allocatable :: rain(:,:)          ! (VIS/IR) cloud precipitating water
   real(r8), allocatable :: ciw(:,:)           ! (VIS/IR) cloud non-precipitating ice concentration
   real(r8), allocatable :: snow(:,:)          ! (VIS/IR) cloud precipitating ice (fluffy)
   real(r8), allocatable :: graupel(:,:)       ! (VIS/IR) cloud precipitating ice (soft hail / snow pellets)
   real(r8), allocatable :: hail(:,:)          ! (VIS/IR) cloud precipitating ice (hard hail)
   real(r8), allocatable :: clwde(:,:)         ! (VIS/IR)    cloud liquid effective diameter if clw_scheme = 2
   real(r8), allocatable :: icede(:,:)         ! (VIS/IR)    cloud liquid effective diameter if ice_scheme = 1
end type cloud_profile_type

! Container for the DART/rttov run-time structures to be used (per sensor)
! This is setup per sensor in accordance with Figure 1 of the RTTOV user guide,
! which indicates the allocation is per coefficient file, which is
! per sensor. RTTOV_alloc_direct also takes in the coefficient file. It's
! not clear if some of these structures could be reused across different
! sensors - to be safe, they are gathered together here. There will be 
! one runtime per sensor, but the runtime will only be allocated when the 
! sensor is used by a forward operator.
type rttov_sensor_runtime_type
   type(rttov_options),       pointer :: opts               => null() ! Options for RTTOV-DIRECT 
   type(rttov_options_scatt), pointer :: opts_scatt         => null() ! Options for RTTOV-SCATT
   type(rttov_coefs)                  :: coefs                        ! Sensor coefficients structure
   type(rttov_scatt_coef),    pointer :: coefs_scatt        => null() ! Sensor coefficients structure
   type(rttov_chanprof),      pointer :: chanprof(:)        => null() ! Input channel/profile list
   integer(kind=jpim),        pointer :: frequencies(:)     => null() ! Indexes into Mietable lookup
   integer(kind=jpim),        pointer :: frequencies_all(:) => null() ! Indexes into Mietable lookup
   logical(kind=jplm),        pointer :: calcemis(:)        => null() ! Flag to indicate calculation of emissivity within RTTOV
   type(rttov_emissivity),    pointer :: emissivity(:)      => null() ! Input/output surface emissivity
   logical(kind=jplm),        pointer :: calcrefl(:)        => null() ! Flag to indicate calculation of BRDF within RTTOV
   type(rttov_reflectance),   pointer :: reflectance(:)     => null() ! Input/output surface BRDF
   type(rttov_profile),       pointer :: profiles(:)        => null() ! Input profiles
   type(rttov_profile_cloud), pointer :: cld_profiles(:)    => null() ! Input profiles
   type(rttov_transmission)           :: transmission                 ! Output transmittances
   type(rttov_radiance)               :: radiance                     ! Output radiances
end type rttov_sensor_runtime_type

! a simple fortran type containing the generic sensor information 
type rttov_sensor_type
   integer              :: platform_id
   integer              :: satellite_id
   integer              :: sensor_id
   character(len=512)   :: base_obs_type
   character(len=3)     :: sensor_type_name ! mw, ir, or vis
   integer              :: sensor_type_num  ! vis = 1, ir = 2, mw = 3
   character(len=512)   :: coefficient_file
   integer, allocatable :: channels(:)      ! list of channels to use

   ! Runtime per-sensor RTTOV structures used to execute the forward operator
   type(rttov_sensor_runtime_type), pointer :: runtime => null()

   ! Linked-list of rttov_sensor_types, sorted (asc) by sensor_id 
   type(rttov_sensor_type),         pointer :: next_sensor  => null()
end type rttov_sensor_type

type rttov_satellite_type
   ! this is a singly-linked sorted (asc by sensor_id) list of sensors on the satellite
   type(rttov_sensor_type), pointer :: head
end type rttov_satellite_type

type rttov_platform_type
   ! this is a trie - each array index corresponds to the satellite number
   type(rttov_satellite_type), pointer :: sats(:)
end type rttov_platform_type

! this is a trie - each array index corresponds to the platform number
type(rttov_platform_type), pointer :: platforms(:)

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_rttov13_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

character(len=512) :: string1, string2
logical            :: module_initialized = .false.
logical            :: arrays_prealloced  = .false.
integer            :: iunit, rc

integer,                   allocatable :: obstype_metadata(:,:) ! key & subkey
type(visir_metadata_type), pointer     :: visir_obs_metadata(:)
type(visir_metadata_type)              :: missing_visir_metadata
type(mw_metadata_type),    pointer     :: mw_obs_metadata(:)
type(mw_metadata_type)                 :: missing_mw_metadata

character(len=5), parameter :: VISIR_STRING = 'visir'
character(len=5), parameter :: MW_STRING    = 'mw   '

integer, parameter :: NO_OBS = -1
integer, parameter :: VISIR = 1
integer, parameter :: MW = 2

! row in obstype_metadata(:,:)
integer, parameter :: SUBTYPE = 1 
integer, parameter :: SUBKEY = 2

logical :: debug = .false.
integer :: MAXrttovkey = 100000  !FIXME - some initial number of obs
integer ::    rttovkey = 0       ! useful length of metadata arrays
integer ::    visirnum = 0
integer ::       mwnum = 0

character(len=512)   :: rttov_sensor_db_file = 'unspecified'

! -----------------------------------------------------------------------------
! DART/RTTOV options in the input.nml namelist.
! 
! DART exposes all of the RTTOV 13 options available and passes them to 
! RTTOV with little to no additional checking for consistency. The default in 
! most cases can be used and need not be specified in the namelist. 
!
! See the RTTOV user guide Annex O for more information on these options.
!
! Note a new option must be also be transfered in 
! initialize_rttov_sensor_runtime
! -----------------------------------------------------------------------------
logical              :: first_lvl_is_sfc     = .true.   ! is level 1 the surface (true) or top of atmosphere (false)?
logical              :: mw_clear_sky_only    = .false.  ! only use clear-sky for MW (plus clw emission if clw_data is true) or full RTTOV-SCATT (false)?
integer              :: interp_mode          = 1        ! Interpolation mode: Rochon on OD (1), Log-linear (2), Rochon on log-linear OD (3), Rochon on WF (4), Rochon on log-linear WF (5)
logical              :: do_checkinput        = .true.   ! check if profiles are within absolute and regression limits
logical              :: apply_reg_limits     = .false.  ! clamp to min/max values
logical              :: verbose              = .true.   ! if false, only fatal errors output 
logical              :: fix_hgpl             = .true.   ! surface elevation assigned to 2m pressure (true) or surface pressure (true)
logical              :: do_lambertian        = .false.  ! treat surface as Lambertian instead of specular? (all)
logical              :: lambertian_fixed_angle = .true. ! use fixed angle for Lambertian calculations? (all, do_lambertian only)
logical              :: rad_down_lin_tau     = .true.   ! use linear-in-tau approximation? (all)
real(r8)             :: max_zenith_angle     = 75.      ! maximum zenith angle to accept (in degrees) (all)
logical              :: use_q2m              = .false.  ! use surface humidity? (all)
logical              :: use_uv10m            = .false.  ! use u and v 10 meters? (all, used in sea surface emissivity and BRDF models)
logical              :: use_wfetch           = .false.  ! use wind fetch (length of water wind has blown over in m)  (all, used in sea surface BRDF models)
logical              :: use_water_type       = .false.  ! use water type (0 = fresh, ocean = 1) (all, used in surface BRDF atlas and models)
logical              :: addrefrac            = .true.   ! enable atmospheric refraction (all) 
logical              :: plane_parallel       = .false.  ! treat atmosphere as strictly plane-parallel? (all)
logical              :: use_salinity         = .false.  ! use ocean salinity (in practical salinity units) (MW, FASTEM 4-6 and TESSEM2)
logical              :: cfrac_data           = .false.  ! specify cloud fraction? (VIS/IR/MW)
logical              :: clw_data             = .false.  ! specify non-precip cloud liquid water? (VIS/IR/MW)
logical              :: rain_data            = .false.  ! specify precip cloud liquid water? (VIS/IR/MW)
logical              :: ciw_data             = .false.  ! specify non-precip cloud ice? (VIS/IR)
logical              :: snow_data            = .false.  ! specify precip cloud fluffy ice? (VIS/IR/MW)
logical              :: graupel_data         = .false.  ! specify precip cloud soft-hail? (VIS/IR/MW)
logical              :: hail_data            = .false.  ! specify precip cloud hard-hail? (VIS/IR/MW)
logical              :: w_data               = .false.  ! specify vertical velocity (used for classifying clouds as cumulus versus stratus)? (VIS/IR)
integer              :: clw_scheme           = 2        ! Liebe (1) or Rosenkranz (2) or TKC (3) (MW, clear-sky only)
real(r8)             :: clw_cloud_top        = 322.0_r8   ! lower hPa limit for clw calculations; clw at lower pressures is ignored (MW, clear-sky only)
integer              :: fastem_version       = 6        ! MW sea-surface emissivity model to use (0-6). 1-6: FASTEM version 1-6, 0: TESSEM2 (MW)
logical              :: supply_foam_fraction = .false.  ! include foam fraction in skin%foam_fraction? FASTEM only. (MW)
logical              :: use_totalice         = .false.  ! Specify totalice instead of precip/non-precip ice (MW, RTTOV-SCATT only)
logical              :: use_zeeman           = .false.  ! Simulate Zeeman effect (MW)
real(r8)             :: cc_threshold         = 0.001_r8   ! if effective cloud fraction below this value, treat simulation as clear-sky (MW, 0-1, RTTOV-SCATT only)
logical              :: ozone_data           = .false.  ! specify ozone profiles? (VIS/IR)
logical              :: co2_data             = .false.  ! specify CO2 profiles? (VIS/IR)
logical              :: n2o_data             = .false.  ! specify N2O profiles? (VIS/IR)
logical              :: co_data              = .false.  ! specify CO profiles? (VIS/IR)
logical              :: ch4_data             = .false.  ! specify CH4 profiles? (VIS/IR)
logical              :: so2_data             = .false.  ! specify SO2 profiles? (VIS/IR)
logical              :: addsolar             = .false.  ! include solar calculations (VIS/IR)
logical              :: rayleigh_single_scatt = .true.  ! if false, disable Rayleigh (VIS, addsolar only)
logical              :: do_nlte_correction   = .false.  ! if true include non-LTE bias correction for hires sounders (VIS/IR)
integer              :: solar_sea_brdf_model = 2        ! JONSWAP (1) or Elfouhaily (2) (VIS)
integer              :: ir_sea_emis_model    = 2        ! ISEM (1) or IREMIS (2) (IR)
logical              :: use_sfc_snow_frac    = .false.  ! use sfc snow cover (0-1) (IR, used in emis atlas)
logical              :: add_aerosl           = .false.  ! enable aerosol scattering (VIS/IR)
integer              :: aerosl_type          = 1        ! OPAC (1) or CAMS (2) (VIS/IR, add_aerosl only)
logical              :: add_clouds           = .true.   ! enable cloud scattering (VIS/IR)
integer              :: ice_scheme           = 1        ! SSEC (1) or Baran 2014 (2) or Baran 2018 (3) (VIS/IR, add_clouds only)
logical              :: use_icede            = .false.  ! use ice effective diameter (IR, add_clouds, ice_scheme = 1) 
integer              :: idg_scheme           = 2        ! Ou and Liou (1), Wyser (2), Boudala (3), McFarquar (2003) (VIS/IR, add_clouds only, ice_scheme = 1)
logical              :: user_aer_opt_param   = .false.  ! specify aerosol scattering properties (VIS/IR, add_clouds only)
logical              :: user_cld_opt_param   = .false.  ! specify cloud scattering properties (VIS/IR, add_clouds only)
logical              :: grid_box_avg_cloud   = .true.   ! cloud concentrations are grid box averages. False = concentrations for cloudy layer only. (VIS/IR, add_clouds and not user_cld_opt_param only)
real(r8)             :: cldcol_threshold     = -1.0_r8    ! threshold for cloud stream weights for scattering (VIS/IR, add_clouds only)
integer              :: cloud_overlap        = 1        ! default: 1 (max/random overlap)
real(r8)             :: cc_low_cloud_top     = 750.0_r8   ! cloud fraction maximum in layers from ToA down to specified hPa (VIS/IR, cloud_overlap only)
integer              :: ir_scatt_model       = 2        ! DOM (1) or Chou-scaling (2) (IR, add_clouds or add_aerosl only)
integer              :: vis_scatt_model      = 1        ! DOM (1), single scat (2), or MFASIS (3) (VIS, addsolar and add_clouds or add_aerosl only)
integer              :: dom_nstreams         = 8        ! number of streams to use with DOM (VIS/IR, add_clouds or add_aerosl and DOM model only, must be >= 2 and even)
real(r8)             :: dom_accuracy         = 0.0_r8     ! convergence criteria for DOM (VIS/IR, add_clouds or addaerosol and DOM model only)
real(r8)             :: dom_opdep_threshold  = 0.0_r8     ! DOM ignores layers below this optical depth (VIS/IR, add_clouds or addaerosol and DOM model only)
logical              :: addpc                = .false.  ! do principal component calculations? (VIS/IR)
integer              :: npcscores            = -1       ! number of PC scores to use (VIS/IR, addpc only)
logical              :: addradrec            = .false.  ! reconstruct the radiances (VIS/IR, addpc only)
integer              :: ipcreg               = 1        ! number of predictors, see Table 29 of user guide (VIS/IR, addpc only)
logical              :: use_htfrtc           = .false.  ! use HTFRTC of Havemann 2018  
integer              :: htfrtc_n_pc          = -1       ! number of PCs to use (HTFRTC only, max 300)
logical              :: htfrtc_simple_cloud  = .false.  ! use simple-cloud scattering (HTFRTC only)
logical              :: htfrtc_overcast      = .false.  ! calculate overcast radiances (HTFRTC only)
real(r8)             :: wfetc_value          = 100000.0_r8 ! Real wfetc Wind fetch (m) (length of water over which the wind has blown, typical
                                                           ! value 100000m for open ocean). Used if wfetc not provided by model.

namelist / obs_def_rttov_nml/ rttov_sensor_db_file,   &
                              first_lvl_is_sfc,       &
                              mw_clear_sky_only,      &
                              interp_mode,            &
                              do_checkinput,          &
                              apply_reg_limits,       &
                              verbose,                &
                              fix_hgpl,               &
                              do_lambertian,          &
                              lambertian_fixed_angle, &
                              rad_down_lin_tau,       &
                              max_zenith_angle,       &
                              use_q2m,                &
                              use_uv10m,              &
                              use_wfetch,             &
                              use_water_type,         &
                              addrefrac,              &
                              plane_parallel,         &
                              use_salinity,           &
                              cfrac_data,             &
                              clw_data,               &
                              rain_data,              &
                              ciw_data,               &
                              snow_data,              &
                              graupel_data,           &
                              hail_data,              &
                              w_data,                 &
                              clw_scheme,             &
                              clw_cloud_top,          &
                              fastem_version,         &
                              supply_foam_fraction,   &
                              use_totalice,           &
                              use_zeeman,             &
                              cc_threshold,           &
                              ozone_data,             &
                              co2_data,               &
                              n2o_data,               &
                              co_data,                &
                              ch4_data,               &
                              so2_data,               &
                              addsolar,              &
                              rayleigh_single_scatt,  &
                              do_nlte_correction,     &
                              solar_sea_brdf_model,   &
                              ir_sea_emis_model,      &
                              use_sfc_snow_frac,      &
                              add_aerosl,             &
                              aerosl_type,            &
                              add_clouds,             &
                              ice_scheme,             &
                              use_icede,              &
                              idg_scheme,             &
                              user_aer_opt_param,     &
                              user_cld_opt_param,     &
                              grid_box_avg_cloud,     &
                              cldcol_threshold,       &
                              cloud_overlap,          &
                              cc_low_cloud_top,       &
                              ir_scatt_model,         &
                              vis_scatt_model,        &
                              dom_nstreams,           &
                              dom_accuracy,           &
                              dom_opdep_threshold,    &
                              addpc,                  &
                              npcscores,              &
                              addradrec,              &
                              ipcreg,                 &
                              use_htfrtc,             &
                              htfrtc_n_pc,            &
                              htfrtc_simple_cloud,    &
                              htfrtc_overcast,        &
                              wfetc_value

type(atmos_profile_type)     :: atmos
type(trace_gas_profile_type) :: trace_gas
type(aerosol_profile_type)   :: aerosols
type(cloud_profile_type)     :: clouds
  
! include the interface files as per the RTTOV standard
include "rttov_direct.interface"
include "rttov_read_coefs.interface"
include "rttov_dealloc_coefs.interface"
include "rttov_read_scattcoeffs.interface"
include "rttov_dealloc_scattcoeffs.interface"
include "rttov_scatt_setupindex.interface"
include "rttov_scatt.interface"
include "rttov_alloc_direct.interface"
include "rttov_alloc_scatt_prof.interface"
include "rttov_user_options_checkinput.interface"
include "rttov_print_opts.interface"
include "rttov_print_opts_scatt.interface"
include "rttov_print_profile.interface"
include "rttov_print_cld_profile.interface"
include "rttov_skipcommentline.interface"

!--------------------------
!
integer(kind=jpim), parameter :: ioout = 0    ! stdout for now

integer, parameter :: NUM_PLATFORMS_INITIAL  = 60 ! first guess of # of platforms
integer, parameter :: NUM_SATELLITES_INITIAL = 25 ! first guess of # of satellites per platform

! arrays of length nlevels, initialized in do_forward_model
integer,  allocatable :: lvlidx(:)
integer,  allocatable :: ly1idx(:)
integer,  allocatable :: ly2idx(:)
real(jprb), allocatable :: totalwater(:)
real(jprb), allocatable :: totalice(:)

contains

!----------------------------------------------------------------------
! Add a sensor to the platform/sat/sensor data structures with the 
! given information. Note that this will create the platform, satellite,
! and/or sensors as necessary.
! 
subroutine add_sensor(platform_id, satellite_id, sensor_id, sensor_type_name, &
   base_obs_type, coef_file, channels)

   integer,          intent(in) :: platform_id        ! the RTTOV platform id
   integer,          intent(in) :: satellite_id       ! the RTTOV satellite id
   integer,          intent(in) :: sensor_id          ! the RTTOV sensor id
   character(len=*), intent(in) :: sensor_type_name   ! the type of sensor (vis, ir, mw)
   character(len=*), intent(in) :: base_obs_type      ! the DART observation type without _RADIANCE, _TB, or _BDRF
   character(len=*), intent(in) :: coef_file          ! the RTTOV coefficient file
   integer,          intent(in) :: channels(:)        ! the channels to simulate (or 0-length for all)

   integer :: i, new_size

   type(rttov_platform_type),  pointer :: tmp_platforms(:)   ! temporary platform array, for dynamic growing
   type(rttov_satellite_type), pointer :: tmp_satellites(:)  ! temporary satellite array, for dynamic growing

   type(rttov_platform_type),  pointer :: platform          
   type(rttov_satellite_type), pointer :: satellite         
   type(rttov_sensor_type),    pointer :: sensor           
   type(rttov_sensor_type),    pointer :: new_sensor
   type(rttov_sensor_type),    pointer :: last_sensor

   character(len=*), parameter :: routine = 'add_sensor'

   character(len=512) :: string1
   character(len=512) :: string2

   string2 = 'name: ' // trim(base_obs_type)

   ! check 
   if (platform_id <= 0 .or. satellite_id < 0 .or. sensor_id < 0) then
      write(string1,*)'Invalid platform/satellite/sensor id:',&
         platform_id, ',',satellite_id,',',sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if

   ! allocate the platforms trie if necessary
   if (.not. associated(platforms)) then
      allocate(platforms(NUM_PLATFORMS_INITIAL))

      do i=1,size(platforms)
         ! by default, leave the contained sats trie null
         platforms(i) % sats => null()
      end do

      if (debug) then
         write(string1,*)'allocated platforms:',size(platforms)
         call error_handler(E_MSG, routine, string1, source, revision, revdate)
      end if
   end if

   ! check if we need to expand platforms array
   ! can't combine this with the above as platform_id may be too large at different times
   if (platform_id > size(platforms)) then

      new_size = 2*size(platforms)

      if (platform_id > new_size) then
         write(string1,*) 'Error: platform ID is too large: ',platform_id,'vs.',new_size
         call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
      else
         write(string1,*)'Resizing platform size from ',size(platforms),'to',&
            new_size
         call error_handler(E_WARN, routine, string1, source, revision, revdate)

         ! increase the size of platforms array
         allocate(tmp_platforms(new_size))

         do i=1,size(platforms)
            ! copy the pointer to the sats array
            tmp_platforms(i) % sats => platforms(i) % sats 
         end do

         do i=size(platforms)+1,new_size
            tmp_platforms(i) % sats => null()
         end do

         ! make the switch. deallocate the platform array
         deallocate(platforms)
         platforms => tmp_platforms
      end if 
   end if

   platform => platforms(platform_id)

   ! allocate the sats trie if necessary
   if (.not. associated(platform % sats)) then
      ! there are satellites with id 0, so start from 0 index
      allocate(platform % sats(0:NUM_SATELLITES_INITIAL-1))

      ! initialize the sensor linked list to null
      do i=0,size(platform % sats)-1
         platform % sats(i) % head => null()
      end do

      if (debug) then
         write(string1,*) 'allocated platform sats:',platform_id,size(platform % sats)
         call error_handler(E_MSG, routine, string1, source, revision, revdate)
      end if
   end if

   ! check if we would need to expand the sats array
   if (satellite_id > size(platform % sats)-1) then

      new_size = 2*size(platform % sats)

      if (satellite_id > new_size-1) then
         write(string1,*) 'Error: satellite ID is too large: ',satellite_id,'vs.',new_size
         call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
      else
         write(string1,*) 'Resizing platform % sats size from ',size(platform % sats),'to',&
            new_size
         call error_handler(E_WARN, routine, string1, source, revision, revdate)

         allocate(tmp_satellites(0:new_size))

         ! copy the pointers for the linked list over
         do i=0,size(platform % sats)-1
            tmp_satellites(i) % head => platform % sats(i) % head
         end do

         do i=size(platform % sats),new_size-1
            tmp_satellites(i) % head => null()
         end do

         ! make the switch
         deallocate(platform % sats)
         platform % sats => tmp_satellites
      end if 
   end if

   ! get the satellite pointer from the sats trie
   satellite => platform % sats(satellite_id)

   ! setup our new sensor structure to be added
   allocate(new_sensor)
   new_sensor % platform_id   = platform_id
   new_sensor % satellite_id  = satellite_id
   new_sensor % sensor_id     = sensor_id
   new_sensor % base_obs_type = base_obs_type
   new_sensor % sensor_type_name = sensor_type_name

   if (trim(sensor_type_name) == 'vis') then
      new_sensor % sensor_type_num = 1
   elseif (trim(sensor_type_name) == 'ir') then
      new_sensor % sensor_type_num = 2
   elseif (trim(sensor_type_name) == 'mw') then
      new_sensor % sensor_type_num = 3
   elseif (trim(sensor_type_name) == 'po') then
      new_sensor % sensor_type_num = 4
   elseif (trim(sensor_type_name) == 'hi') then
      new_sensor % sensor_type_num = 5
   else
      write(string1,*) 'Error: unknown sensor_type_name: ', trim(sensor_type_name),':',&
         platform_id,satellite_id,sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if

   new_sensor % coefficient_file = coef_file
   allocate(new_sensor % channels(size(channels)))
   new_sensor % channels(:) = channels(:)
   new_sensor % next_sensor => null()

   ! traverse the sensor list. Insert the node so the list is sorted.
   last_sensor => null()
   sensor => satellite % head

   linkedlist: do
      if (.not. associated(sensor)) then
         ! reached the end of the list. Add new_sensor to the end
         if (.not. associated(last_sensor)) then
            ! update head => new_sensor
            satellite % head => new_sensor
         else
            ! update last_sensor % next => sensor
            last_sensor % next_sensor => new_sensor
         end if

         exit ! the linked list loop
      end if

      ! check for duplicates
      if (sensor_id == sensor % sensor_id) then
         write(string1,*) 'Error: tried to add the same sensor to a satellite twice:',&
            platform_id,satellite_id,sensor_id
         call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
      end if

      ! if the sensor id to add is greater than the current one, we found our spot:
      ! we need to put the new sensor in between sensor and last_sensor
      if (sensor_id > sensor % sensor_id) then

         ! update new_sensor % next => sensor
         new_sensor % next_sensor => sensor

         if (associated(last_sensor)) then
            ! update last_sensor % next => new_sensor % next => sensor
            last_sensor % next_sensor => new_sensor
         else
            ! update head => new_sensor % next => sensor
            satellite % head => new_sensor
         end if

         exit ! the linked list loop
      else
         ! traverse forward in the list and keep looking
         last_sensor => sensor
         sensor => sensor % next_sensor
      end if
   end do linkedlist
end subroutine add_sensor

!----------------------------------------------------------------------
! Read the sensor DB file and setup the relevant data structures.
! This function will read through the DB file twice: once to 
! discover how many lines are in the file and the maximum line
! length, and a second time to actually parse the contents of the
! file. The contents of the file are assumed to be an "unformatted"
! (aka text) CSV file with the following contents:
! 
! <base_obs_type>,<platform_id>,<satellite_id>,<sensor_id>,<sensor_type>,<coef_file>,[channel list]
!
! where: 
!     base_obs_type is the DART observation quantity minus _RADIANCE, _TB, or _BDRF
!     platform_id   is the RTTOV platform id
!     satellite_id  is the RTTOV satellite id
!     sensor_id     is the RTTOV sensor id
!     sensor_type   is vis, ir, or mw (visible, infrared, or microwave)
!     coef_file     is the RTTOV coefficient file
!     channel_list  is an optional list of channels (can be zero-length, meaning all available channels should be used) 
! 
subroutine read_sensor_db_file(sensor_db_file)

character(len=*), intent(in) :: sensor_db_file

! a buffer to use for reading the file piece by piece
integer, parameter :: buffer_len = 80
character(len=buffer_len) :: buffer
character(len=2048) :: next_line

integer :: nlines, io, size_read, line_len, maxline_len
integer :: i
character(len=100) :: imsg

character(len=256) :: base_obs_type
integer :: dbUnit
integer :: platform_id
integer :: satellite_id
integer :: sensor_id
integer :: ntok, toknum, lastind
integer, allocatable :: channels(:)
character(len=16)   :: sensor_type
character(len=512)  :: coef_file
character(len=512)  :: token
! assume the file is delimited by commas
character(len=*), parameter :: delimiter = ','

character(len=*), parameter :: routine = 'read_sensor_db_file'

character(len=512) :: string1, string2

if ( file_exist(sensor_db_file) ) then
   dbUnit = open_file(sensor_db_file, action='read')
else
   write(string1,*)'sensor_db_file: "'//adjustl(trim(sensor_db_file))//'"'
   write(string2,*)'does not exist.'
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

! throw away the header row
read(dbUnit,'(A)') buffer

! Now count the number of lines in the file.
! While we are at it, find the size of a string large enough to hold
! the longest line for error checking
nlines = 0
maxline_len = -1

testlineloop: do 
   line_len = 0
   ! the next loop finds the length of each line
   bufferloop: do 
      read(dbUnit,'(A)',iostat=io,iomsg=imsg,advance='no',size=size_read) buffer

      if (io == 0) then
         ! 0 = no error, read the full buffer, continue
         line_len = line_len + buffer_len
      elseif (io == -1) then
         ! -1 = end of file reached. Exit both loops.
         exit ! the buffer loop
      elseif (io <= 0) then
         ! Different compilers have different end of record iostat values
         ! All return a negative value not equal to -1, however.
         line_len = line_len + size_read
         nlines = nlines + 1
         maxline_len = max(maxline_len,line_len)
         exit ! the buffer loop
      else
         ! for any other error, quit
         write(string1,*)'fatal error:',io,imsg
         call error_handler(E_ERR, routine, string1, source, revision, revdate)
      end if
   end do bufferloop

   if (io == -1) then
      ! all compilers have -1 = end of file
      exit ! the line loop
   end if
end do testlineloop

if (debug) then
   write(string1,*)'found nlines:',nlines,'max len:',maxline_len
   call error_handler(E_MSG, routine, string1, source, revision, revdate)
end if

if (maxline_len > len(next_line)) then
   write(string1,*)'The length of a line in the ' // trim(sensor_db_file) &
      //' file exceeds the maximum: ',maxline_len,'vs.',len(next_line)
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

rewind(dbUnit)

! throw away the header again
read(dbUnit,'(A)') buffer

lineloop: do
   read(dbUnit,'(A)', iostat=io) next_line
   if (io == 0) then
      if (debug) then
         write(string1,*)'read the line:',trim(next_line)
         call error_handler(E_MSG, routine, string1, source, revision, revdate)
      end if
   else
      exit ! the line loop
   end if

   toknum = 0
   lastind = 1

   ntok = 0

   ! first pass: count the number of tokens
   do i=1,len(next_line)
      if (next_line(i:i) == delimiter .or. i==len(next_line)) then
         ntok = ntok + 1
      end if
   end do

   if (debug) then
      write(string1,*)'found ntokens:',ntok
      call error_handler(E_MSG, routine, string1, source, revision, revdate)
   end if

   allocate(channels(ntok-6))

   ! second pass: parse the tokens
   do i=1,len(next_line)
      if (next_line(i:i) == delimiter .or. i == len(next_line)) then
         toknum = toknum + 1
         token = adjustl(trim(next_line(lastind:i-1)))
         if (debug) then
            write(string1,*)'found the token :',toknum,trim(token)
            call error_handler(E_MSG, routine, string1, source, revision, revdate)
         end if
         select case (toknum)
            case (1)
               base_obs_type = token
            case (2)
               platform_id = str2int(token)
            case (3)
               satellite_id = str2int(token)
            case (4)
               sensor_id = str2int(token)
            case (5)
               sensor_type = token
            case (6)
               coef_file = token
            case default
               channels(toknum-6) = str2int(token)
         end select

         lastind = i + 1
      end if
   end do

   if (debug) then
      write(string1,*)'values:',trim(base_obs_type),platform_id,satellite_id,sensor_id,&
         trim(sensor_type),trim(coef_file),channels
      call error_handler(E_MSG, routine, string1, source, revision, revdate)
   end if

   call add_sensor(platform_id, satellite_id, sensor_id, sensor_type, base_obs_type, &
      coef_file, channels)

   deallocate(channels)
end do lineloop

call close_file(dbUnit)

end subroutine read_sensor_db_file

!----------------------------------------------------------------------
! Function to return the sensor associated with the given platform/sat/
! sensor id. If no sensor is found or if invalid id combinations are 
! passed in, an error will be generated.
! 
function get_rttov_sensor(instrument_id) result(sensor)

integer, intent(in) :: instrument_id(3)

integer :: platform_id
integer :: satellite_id
integer :: sensor_id

type(rttov_platform_type),  pointer :: platform
type(rttov_satellite_type), pointer :: satellite
type(rttov_sensor_type),    pointer :: sensor

character(len=*), parameter :: routine = 'get_rttov_sensor'

character(len=512) :: string1

platform_id  = instrument_id(1)
satellite_id = instrument_id(2)
sensor_id    = instrument_id(3)

if (platform_id <= 0 .or. platform_id > size(platforms) .or. &
    satellite_id < 0 .or. sensor_id < 0) then
   write(string1,*)'Invalid platform/satellite/sensor id:', platform_id, ',',&
      satellite_id,',',sensor_id
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

platform => platforms(platform_id)

if (satellite_id > size(platform % sats)-1) then
   write(string1,*)'Invalid satellite id:', satellite_id,size(platform % sats)-1
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

satellite => platform % sats(satellite_id)

sensor => satellite % head

linkedlist: do
   ! reached the end of the list, but didn't find what we were looking for.
   if (.not. associated(sensor)) then
      write(string1,*)'Could not find the sensor with the platform/satellite/sensor id:',&
         platform_id,satellite_id,sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if

   ! check for our id
   if (sensor_id == sensor % sensor_id) then
      ! found it - return without error (function result is sensor)
      return
   end if

   ! if the sensor id is greater than the current one, we are out of luck
   if (sensor_id > sensor % sensor_id) then
      write(string1,*)'Could not find the sensor with the platform/satellite/sensor id:',&
         platform_id,satellite_id,sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   else
      ! moving on to the next sensor in the list
      sensor => sensor % next_sensor 
   end if
end do linkedlist 

end function get_rttov_sensor

!----------------------------------------------------------------------
! Clean up (deallocate and nullify) the platforms/satellite/sensors data 
! structures.
! 
subroutine clean_up_sensors()

type(rttov_platform_type),  pointer :: platform
type(rttov_satellite_type), pointer :: satellite
type(rttov_sensor_type),    pointer :: sensor
type(rttov_sensor_type),    pointer :: next_sensor

integer :: i, j, error_status

do i=1,size(platforms)
   platform => platforms(i)

   do j=1,size(platform % sats)
      satellite => platform % sats(j)

      sensor => satellite % head

      do while (associated(sensor))
         next_sensor => sensor % next_sensor

         if (allocated(sensor % channels)) then
            deallocate(sensor % channels)
         end if

         if (associated(sensor % runtime)) then
            call sensor_runtime_takedown(sensor % runtime, error_status) 
            deallocate(sensor % runtime)
         end if

         deallocate(sensor)

         sensor => next_sensor
      end do

      nullify(satellite % head)
   end do

   deallocate(platform % sats)
   nullify(platform % sats)
end do

deallocate(platforms)
nullify(platforms)
end subroutine clean_up_sensors

!----------------------------------------------------------------------
! A helper function to convert a character(len=*) into an integer
!
function str2int(str) result(intv)

character(len=*), intent(in) :: str

integer :: strlength
integer :: intv

strlength = len_trim(adjustl(str))
if (strlength > 9) then
   print *,'Error: integer string length is greater than 9 digits long.'
   stop
end if 

read(str,*) intv

end function str2int

!----------------------------------------------------------------------
! Setup the sensor runtime structures for use with RTTOV in DART. Note
! that this loads the RTTOV coefficient file, allocates memory, sets the
! RTTOV options, and provides some basic-level of error checking.
! 
subroutine sensor_runtime_setup(sensor, ens_size, nlevs, opts, opts_scatt, &
   use_totalice)

type(rttov_sensor_type),             intent(inout)  :: sensor
integer,                             intent(in)     :: ens_size
integer,                             intent(in)     :: nlevs
type(rttov_options),                 pointer        :: opts
type(rttov_options_scatt), optional, pointer        :: opts_scatt
logical,                   optional, intent(in)     :: use_totalice ! only relevant if opts_scatt is present

character(len=*), parameter :: routine = 'sensor_runtime_setup'

character(len=512) :: string1
character(len=512) :: string2 ! used to hold the sensor % base_obs_type, don't overwrite

type(rttov_sensor_runtime_type), pointer :: runtime

integer :: nchannels

! Return error status of RTTOV subroutine calls
integer(kind=jpim)               :: errorstatus 

logical :: do_totalice

integer :: instrument(3) ! platform_id, satellite_id, sensor_id
integer :: alloc_status
integer :: i, ich, nch

logical(kind=jplm), pointer :: use_chan(:,:)  => null() ! flags to specify channels to simulate

! used below, don't overwrite
string2 = 'name: ' // trim(sensor % base_obs_type)

instrument(1) = sensor % platform_id
instrument(2) = sensor % satellite_id
instrument(3) = sensor % sensor_id

if (associated(sensor % runtime)) then
   write(string1,*)'Runtime information already initialized for platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

if (debug) then
   write(string1,*) 'Now initializing the platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_MSG, routine, string1, source, revision, revdate, text2=string2)
end if

allocate(runtime)
sensor % runtime => runtime

! --------------------------------------------------------------------------
! 1. Set options 
! --------------------------------------------------------------------------

runtime % opts => opts

if (present(opts_scatt)) then
   runtime % opts_scatt => opts_scatt
end if

if (present(use_totalice)) then
   do_totalice = use_totalice
else
   do_totalice = .false.
end if

! --------------------------------------------------------------------------
! 2. Read coefficients
! --------------------------------------------------------------------------

if (debug) then
   write(string1,*)'The coefficient file is:',trim(sensor % coefficient_file)
   call error_handler(E_MSG, routine, string1, source, revision, revdate)
end if

call rttov_read_coefs(errorstatus, runtime % coefs, runtime % opts, instrument=instrument)

if (errorstatus /= errorstatus_success) then
   write(string1,*)'Error reading coefs for platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
endif

if (present(opts_scatt)) then
   allocate(runtime % coefs_scatt)
   call rttov_read_scattcoeffs(errorstatus, runtime % opts_scatt, runtime % coefs, &
      runtime % coefs_scatt)
   if (errorstatus /= errorstatus_success) then
      write(string1,*) 'Error reading scatt coefs for platform/sat/sensor id combination:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   endif
end if

if (size(sensor % channels) /= 0) then
   nchannels = size(sensor % channels)
else
   nchannels = runtime % coefs % coef % fmv_chn
end if

! Ensure input number of channels is not higher than number stored in coefficient file
IF (nchannels > runtime % coefs % coef % fmv_chn) THEN
  write(string1,*) 'Number of requested channels too large: ',nchannels,&
      ' vs. ',runtime % coefs % coef % fmv_chn
  call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif

! Ensure the options and coefficients are consistent
if (runtime % opts % config % do_checkinput) then
   if (debug) then
      write(string1,*) 'Now running rttov_user_options_checkinput'
      call error_handler(E_MSG, routine, string1, source, revision, revdate)
   end if

   call rttov_user_options_checkinput(errorstatus, runtime % opts, runtime % coefs)

   if (errorstatus /= errorstatus_success) then
     write(string1,*) 'error in rttov_user_options_checkinput'
     call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if
end if

! --------------------------------------------------------------------------
! 3. Allocate RTTOV input and output structures
! --------------------------------------------------------------------------

! In DART RTTOV_Direct, we only simulate a single channel at a time for all ens members
allocate(runtime % profiles(ens_size))

! Allocate structures for rttov_direct
CALL rttov_alloc_direct(                  &
      errorstatus,                        &
      1_jpim,                             &  ! 1 => allocate
      ens_size,                           &
      ens_size,                           &
      nlevs,                              &
      runtime % chanprof,                 &
      runtime % opts,                     &
      runtime % profiles,                 &
      runtime % coefs,                    &
      runtime % transmission,             &
      runtime % radiance,                 &
      calcemis    =runtime % calcemis,    &
      emissivity  =runtime % emissivity,  &
      calcrefl    =runtime % calcrefl,    &
      reflectance =runtime % reflectance, &
      init=.TRUE._jplm)

IF (errorstatus /= errorstatus_success) THEN
   write(string1,*) 'allocation error for rttov_direct structures'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif

if (debug) then
   write(string1,*) 'Successfully initialized rttov_direct platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_MSG, routine, string1, source, revision, revdate, text2=string2)
end if

if (present(opts_scatt)) then
   ! Allocate cld_profile for rttov_scatt
   allocate(runtime % cld_profiles(ens_size), stat=alloc_status)
   if (alloc_status /= 0) then
      write(string1,*) 'allocation error for cld_profiles array'
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if

   CALL rttov_alloc_scatt_prof(               &
         err=errorstatus,                     &
         nprof=ens_size,                      &
         cld_profiles=runtime % cld_profiles, &
         nlev=nlevs,                          &
         nhydro=5,                            & ! default 
         nhydro_frac=1,                       & ! 1 cfrac profile
         ! use_totalice=do_totalice,            &
         asw=1,                               & ! 1 = allocate
         init=.TRUE._jplm) !,                    &
         ! mmr_snowrain=.TRUE._jplm)              ! true = kg/kg units for clouds

   if (errorstatus /= errorstatus_success) then
      write(string1,*) 'allocation error for rttov_direct structures'
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   endif

   allocate(use_chan(ens_size,runtime % coefs % coef % fmv_chn))
   allocate(runtime % frequencies(ens_size))

   ! only the channels to simulate will be set
   allocate(runtime % frequencies_all(runtime % coefs % coef % fmv_chn))

   if (size(sensor % channels) /= 0) then
      nch = size(sensor % channels)
   else
      nch = runtime % coefs % coef % fmv_chn
   end if

   do i=1,ens_size
      runtime % chanprof(i) % prof = i
   end do

   do i=1,nch
      use_chan(:,:) = .FALSE._jplm

      if (size(sensor % channels) /= 0) then
         ich = sensor % channels(i)
      else
         ich = i
      end if

      ! Set use_chan to .TRUE. only for the ith required channel
      use_chan(:,ich) = .TRUE._jplm

      runtime % chanprof(:) % chan = ich

      ! Populate chanprof and frequencies arrays for this one channel
      call rttov_scatt_setupindex (                    &
            errorstatus=errorstatus,                   &
            nprofiles=ens_size,                        &
            n_chan=runtime % coefs % coef % fmv_chn,   &
            coef_rttov=runtime % coefs,                &
            coef_scatt=runtime % coefs_scatt,          &
            nchannels=ens_size,                        &
            chanprof=runtime % chanprof,               &
            frequencies=runtime % frequencies,         &
            lchannel_subset=use_chan )

      runtime % frequencies_all(ich) = runtime % frequencies(1)
   end do

   if (debug) then
      write(string1,*) 'Successfully initialized RTTOV-scatt cloud profiles for platform/sat/sensor id combination:',&
         instrument
      call error_handler(E_MSG, routine, string1, source, revision, revdate, text2=string2)
   end if
end if

end subroutine sensor_runtime_setup

subroutine atmos_profile_setup(atmos, ens_size, numlevels, use_q2m, &
      use_uv10m, use_wfetch, use_water_type, use_salinity,          &
      supply_foam_fraction,  use_sfc_snow_frac)

type(atmos_profile_type), intent(inout) :: atmos
integer,                  intent(in)    :: ens_size
integer,                  intent(in)    :: numlevels
logical,                  intent(in)    :: use_q2m
logical,                  intent(in)    :: use_uv10m
logical,                  intent(in)    :: use_wfetch
logical,                  intent(in)    :: use_water_type
logical,                  intent(in)    :: use_salinity
logical,                  intent(in)    :: supply_foam_fraction
logical,                  intent(in)    :: use_sfc_snow_frac

allocate(atmos%temperature(ens_size, numlevels), &
         atmos%   moisture(ens_size, numlevels), &
         atmos%   pressure(ens_size, numlevels), &
         atmos%      sfc_p(ens_size),          &
         atmos%      s2m_t(ens_size),          &
         atmos%  skin_temp(ens_size),          &
         atmos%   sfc_elev(ens_size),          &
         atmos%   surftype(ens_size))

! zero the arrays as well
atmos%temperature = 0.0_jprb
atmos%   moisture = 0.0_jprb
atmos%   pressure = 0.0_jprb
atmos%      sfc_p = 0.0_jprb
atmos%      s2m_t = 0.0_jprb
atmos%  skin_temp = 0.0_jprb
atmos%   sfc_elev = 0.0_jprb
atmos%   surftype = 0.0_jprb

if (use_q2m) then
   allocate(atmos%s2m_q(ens_size))
   atmos%s2m_q = 0.0_jprb
end if

if (use_uv10m) then
   allocate(atmos%s10m_u(ens_size))
   allocate(atmos%s10m_v(ens_size))

   atmos%s10m_u = 0.0_jprb
   atmos%s10m_v = 0.0_jprb
end if

if (use_wfetch) then
   allocate(atmos%wfetch(ens_size))

   atmos%wfetch = 0.0_jprb
end if

if (use_water_type) then
   allocate(atmos%water_type(ens_size))

   atmos%water_type = -1_jpim
end if

if (use_salinity) then
   allocate(atmos%sfc_salinity(ens_size))

   atmos%sfc_salinity = 0.0_jprb
end if

if (supply_foam_fraction) then
   allocate(atmos%sfc_foam_frac(ens_size))
   atmos%sfc_foam_frac = 0.0_jprb
end if

if (use_sfc_snow_frac) then
   allocate(atmos%sfc_snow_frac(ens_size))
   atmos%sfc_snow_frac = 0.0_jprb
end if

end subroutine atmos_profile_setup

subroutine trace_gas_profile_setup(trace_gas, ens_size, numlevels,  &
      ozone_data, co2_data, n2o_data, ch4_data, co_data, so2_data)

type(trace_gas_profile_type), intent(inout) :: trace_gas
integer,                      intent(in)    :: ens_size
integer,                      intent(in)    :: numlevels
logical,                      intent(in)    :: ozone_data
logical,                      intent(in)    :: co2_data
logical,                      intent(in)    :: n2o_data
logical,                      intent(in)    :: ch4_data
logical,                      intent(in)    :: co_data
logical,                      intent(in)    :: so2_data

if (ozone_data) then
   allocate(trace_gas%ozone(ens_size, numlevels))
   trace_gas%ozone = 0.0_jprb
end if

if (co2_data) then
   allocate(trace_gas%co2(ens_size, numlevels))
   trace_gas%co2 = 0.0_jprb
end if

if (n2o_data) then
   allocate(trace_gas%n2o(ens_size, numlevels))
   trace_gas%n2o = 0.0_jprb
end if

if (ch4_data) then
   allocate(trace_gas%ch4(ens_size, numlevels))
   trace_gas%ch4 = 0.0_jprb
end if

if (co_data) then
   allocate(trace_gas%co(ens_size, numlevels))
   trace_gas%co = 0.0_jprb
end if

if (so2_data) then
   allocate(trace_gas%so2(ens_size, numlevels))
   trace_gas%co = 0.0_jprb
end if

end subroutine trace_gas_profile_setup


subroutine aerosol_profile_setup(aerosols, ens_size, numlevels,  &
      aerosl_type)

type(aerosol_profile_type), intent(inout) :: aerosols
integer,                    intent(in)    :: ens_size
integer,                    intent(in)    :: numlevels
integer,                    intent(in)    :: aerosl_type

character(len=512) :: string1
character(len=*), parameter :: routine = 'aerosol_profile_setup'

if (aerosl_type == 1) then
   ! OPAC
   allocate(aerosols%insoluble(ens_size,numlevels)) 
   allocate(aerosols%water_soluble(ens_size,numlevels)) 
   allocate(aerosols%soot(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_accum(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_coarse(ens_size,numlevels)) 
   allocate(aerosols%mineral_nucleus(ens_size,numlevels)) 
   allocate(aerosols%mineral_accum(ens_size,numlevels)) 
   allocate(aerosols%mineral_coarse(ens_size,numlevels)) 
   allocate(aerosols%mineral_transport(ens_size,numlevels)) 
   allocate(aerosols%sulphated_droplets(ens_size,numlevels)) 
   allocate(aerosols%volcanic_ash(ens_size,numlevels)) 
   allocate(aerosols%new_volcanic_ash(ens_size,numlevels)) 
   allocate(aerosols%asian_dust(ens_size,numlevels)) 

   aerosols%insoluble = 0.0_jprb
   aerosols%water_soluble = 0.0_jprb
   aerosols%soot = 0.0_jprb
   aerosols%sea_salt_accum = 0.0_jprb
   aerosols%sea_salt_coarse = 0.0_jprb
   aerosols%mineral_nucleus = 0.0_jprb
   aerosols%mineral_accum = 0.0_jprb
   aerosols%mineral_coarse = 0.0_jprb
   aerosols%mineral_transport = 0.0_jprb
   aerosols%sulphated_droplets = 0.0_jprb
   aerosols%volcanic_ash = 0.0_jprb
   aerosols%new_volcanic_ash = 0.0_jprb
   aerosols%asian_dust = 0.0_jprb
elseif (aerosl_type == 2) then
   ! CAMS
   allocate(aerosols%black_carbon(ens_size,numlevels)) 
   allocate(aerosols%dust_bin1(ens_size,numlevels)) 
   allocate(aerosols%dust_bin2(ens_size,numlevels)) 
   allocate(aerosols%dust_bin3(ens_size,numlevels)) 
   allocate(aerosols%ammonium_sulphate(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_bin1(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_bin2(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_bin3(ens_size,numlevels)) 
   allocate(aerosols%hydrophilic_organic_matter(ens_size,numlevels)) 

   aerosols%black_carbon = 0.0_jprb
   aerosols%dust_bin1 = 0.0_jprb
   aerosols%dust_bin2 = 0.0_jprb
   aerosols%dust_bin3 = 0.0_jprb
   aerosols%ammonium_sulphate = 0.0_jprb
   aerosols%sea_salt_bin1 = 0.0_jprb
   aerosols%sea_salt_bin2 = 0.0_jprb
   aerosols%sea_salt_bin3 = 0.0_jprb
   aerosols%hydrophilic_organic_matter = 0.0_jprb
else
   ! error
   write(string1,*)"Unknown aerosl_type:",aerosl_type
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

end subroutine aerosol_profile_setup

subroutine cloud_profile_setup(clouds, ens_size, numlevels,  &
      cfrac_data, clw_data, clw_scheme, rain_data, ciw_data, &
      ice_scheme, use_icede, snow_data, graupel_data,        &
      hail_data, w_data, htfrtc_simple_cloud)

type(cloud_profile_type), intent(inout) :: clouds
integer,                  intent(in)    :: ens_size
integer,                  intent(in)    :: numlevels
logical,                  intent(in)    :: cfrac_data
logical,                  intent(in)    :: clw_data
integer,                  intent(in)    :: clw_scheme
logical,                  intent(in)    :: rain_data
logical,                  intent(in)    :: ciw_data
integer,                  intent(in)    :: ice_scheme
logical,                  intent(in)    :: use_icede
logical,                  intent(in)    :: snow_data
logical,                  intent(in)    :: graupel_data
logical,                  intent(in)    :: hail_data
logical,                  intent(in)    :: w_data
logical,                  intent(in)    :: htfrtc_simple_cloud

! RTTOV wants layers, but models probably prefer levels
if (cfrac_data) then
   allocate(clouds%cfrac(ens_size, numlevels))
   clouds%cfrac = 0.0_jprb
end if

if (clw_data) then
   allocate(clouds%clw(ens_size, numlevels))
   clouds%clw = 0.0_jprb

   if (clw_scheme == 2) then
      allocate(clouds%clwde(ens_size, numlevels)) 
      clouds%clwde = 20.0_jprb  ! lkugler default value
   end if
end if

if (rain_data) then
   allocate(clouds%rain(ens_size, numlevels))
   clouds%rain = 0.0_jprb
end if

if (ciw_data) then
   allocate(clouds%ciw(ens_size, numlevels)) 
   clouds%ciw = 0.0_jprb
   if (ice_scheme == 1 .and. use_icede) then
      allocate(clouds%icede(ens_size, numlevels))
      clouds%icede = 60.0_jprb  ! lkugler default value
   end if
end if

if (snow_data) then
   allocate(clouds%snow(ens_size, numlevels))
   clouds%snow = 0.0_jprb
end if

if (graupel_data) then
   allocate(clouds%graupel(ens_size, numlevels))
   clouds%graupel = 0.0_jprb
end if

if (hail_data) then
   allocate(clouds%hail(ens_size, numlevels))
   clouds%hail = 0.0_jprb
end if

if (w_data) then
   allocate(clouds%w(ens_size, numlevels))
   clouds%w = 0.0_jprb
end if

if (htfrtc_simple_cloud) then
   allocate(clouds%simple_cfrac(ens_size))
   allocate(clouds%ctp(ens_size))
   clouds%simple_cfrac = 0.0_jprb
   clouds%ctp = 0.0_jprb
end if


end subroutine cloud_profile_setup

!----------------------------------------------------------------------
! Run the forward model for the given sensor, preallocating arrays and
! initializing the sensor if necessary (which will be on the first call
! to this operator and the first time using the sensor, respectively).
!
! istatus: Returns  0 if everything is OK, >0 if error occured.
!                 101 = zenith angle was greater than max_zenith_angle
!                 102 = non-monotonic pressure detected

subroutine do_forward_model(ens_size, nlevels, flavor, location, &
   atmos, trace_gas, clouds, aerosols, sensor, channel,          &
   first_lvl_is_sfc, mw_clear_sky_only, clw_scheme, ice_scheme,  &
   idg_scheme, aerosl_type, do_lambertian, use_totalice,         &
   use_zeeman, radiances, error_status, visir_md, mw_md)

integer,                                intent(in)  :: ens_size
integer,                                intent(in)  :: nlevels
integer,                                intent(in)  :: flavor
type(location_type),                    intent(in)  :: location
type(atmos_profile_type),               intent(in)  :: atmos
type(trace_gas_profile_type),           intent(in)  :: trace_gas
type(cloud_profile_type),               intent(in)  :: clouds
type(aerosol_profile_type),             intent(in)  :: aerosols
type(rttov_sensor_type),                pointer     :: sensor
integer,                                intent(in)  :: channel
logical,                                intent(in)  :: first_lvl_is_sfc
logical,                                intent(in)  :: mw_clear_sky_only
integer,                                intent(in)  :: clw_scheme
integer,                                intent(in)  :: ice_scheme
integer,                                intent(in)  :: idg_scheme
integer,                                intent(in)  :: aerosl_type
logical,                                intent(in)  :: do_lambertian
logical,                                intent(in)  :: use_totalice
logical,                                intent(in)  :: use_zeeman
real(r8),                               intent(out) :: radiances(ens_size)
integer,                                intent(out) :: error_status(ens_size)
type(visir_metadata_type),     pointer, intent(in)  :: visir_md
type(mw_metadata_type),        pointer, intent(in)  :: mw_md

character(len=obstypelength) :: obs_qty_string
integer                      :: obs_type_num

integer :: ilvl, imem 

! observation location variables
real(r8) :: lon, lat, obsloc(3)

integer(kind=jpim) :: j, nch

! Return error status of RTTOV subroutine calls
integer(kind=jpim)               :: errorstatus 

type(rttov_sensor_runtime_type), pointer :: runtime

character(len=512) :: string1
character(len=512) :: string2

character(len=*), parameter :: routine = 'do_forward_model'

real(jprb) :: maxw

logical :: is_visir
logical :: is_mw
logical :: is_cumulus
integer :: instrument(3)
integer :: surftype

if (.not. associated(sensor)) then
   write(string1,*)'Passed an unassociated sensor'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

error_status(:) = 0 ! 0 is success

string2 = 'name: ' // trim(sensor % base_obs_type)

instrument(1) = sensor % platform_id
instrument(2) = sensor % satellite_id
instrument(3) = sensor % sensor_id

is_visir = associated(visir_md)
is_mw    = associated(mw_md)

if (.not. is_visir .and. .not. is_mw) then
   write(string1,*)'Neither vis/ir nor mw metadata were present for platform/sat/sensor id combination:',&
      sensor % platform_id,sensor % satellite_id,sensor % sensor_id
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

if (is_visir .and. is_mw) then
   write(string1,*)'Both vis/ir and mw metadata were present (only one can be specified) for platform/sat/sensor id combination:',&
      sensor % platform_id,sensor % satellite_id,sensor % sensor_id
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

runtime => sensor % runtime

if (.not. associated(runtime)) then
   write(string1,*)'Runtime information not initialized for platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

obsloc   = get_location(location)

lon      = obsloc(1) ! degree: 0 to 360
lat      = obsloc(2) ! degree: -90 to 90
! ignore the 3rd dimension of the observation location for now
!height   = obsloc(3) ! (m)

! --------------------------------------------------------------------------
! 4. Build the list of profile/channel indices in chanprof
! --------------------------------------------------------------------------

nch = 0_jpim
do j = 1, ens_size
    nch = nch + 1_jpim
    runtime % chanprof(nch) % prof = j
    runtime % chanprof(nch) % chan = channel
end do

! We would like a level index array to allow either surface first or surface last order

! One would assume the number of levels would not change between calls, but check
if (allocated(lvlidx) .and. size(lvlidx) /= nlevels) then
   write(string1,*)'Number of levels changed from ',size(lvlidx),' to ',nlevels,&
      ' for platform/sat/sensor id combination:',instrument
   call error_handler(E_WARN, routine, string1, source, revision, revdate, text2=string2)
   
   ! deallocate so we can try again here with the right number of levels
   deallocate(lvlidx)
   deallocate(ly1idx)
   deallocate(ly2idx)
   deallocate(totalwater)
   deallocate(totalice)
end if

! (re)allocate if necessary
if (.not. allocated(lvlidx)) then
   allocate(lvlidx(nlevels))
   allocate(ly1idx(nlevels-1))
   allocate(ly2idx(nlevels-1))
   allocate(totalwater(nlevels))
   allocate(totalice(nlevels))
end if

! finally set the array to the correct order
if (first_lvl_is_sfc) then
   do ilvl=1,nlevels
      lvlidx(ilvl) = nlevels - ilvl + 1
   end do
else
   do ilvl=nlevels,1,-1
      lvlidx(ilvl) = ilvl
   end do
end if

! used for averaging levels to layers
ly1idx = lvlidx(1:nlevels-1)
ly2idx = lvlidx(2:nlevels)

! Loop over all of the ensemble members
! There is one profile per ensemble member
DO imem = 1, ens_size

   ! check that the pressure is monotonically increasing from TOA down
   do ilvl = 1, nlevels-1
      if (atmos % pressure(imem,lvlidx(ilvl)) >= atmos % pressure(imem,lvlidx(ilvl+1))) then
         if (debug) then
            write(string1,*) 'For ens #',imem,', pressure ',lvlidx(ilvl),' was greater than or equal to pressure ',&
                lvlidx(ilvl+1),':',atmos % pressure(imem,lvlidx(ilvl)),' >= ',atmos%pressure(imem,lvlidx(ilvl+1))
            call error_handler(E_ALLMSG,routine,string1,source,revision,revdate)
         end if

         radiances(:) = MISSING_R8
         error_status(:) = 102
         return
      end if
   end do

   surftype   = nint(atmos%surftype(imem))

   runtime % profiles(imem) % nlevels = nlevels
   runtime % profiles(imem) % nlayers = nlevels - 1
   runtime % profiles(imem) % gas_units = 1 ! 1 = kg/kg, 2 = ppmv
   runtime % profiles(imem) % mmr_cldaer = .true. ! kg/kg
   runtime % profiles(imem) % clw_scheme = clw_scheme
   runtime % profiles(imem) % ice_scheme = ice_scheme
   runtime % profiles(imem) % icede_param = idg_scheme
   runtime % profiles(imem) % clwde_param = 1  ! default CLW Deff parameterisation

   ! set the required profile variables 
   runtime % profiles(imem) % p(:) = atmos % pressure(imem,lvlidx)/100.0_jprb  ! Pa -> hPa
   runtime % profiles(imem) % t(:) = atmos % temperature(imem,lvlidx) 
   runtime % profiles(imem) % q(:) = max(atmos % moisture(imem,lvlidx),1e-8_jprb) 

   ! set trace gases if opts present and individual flags are used
   ! the arrays are assumed to have been allocated - this was already checked in obs_def
   if (is_visir .or. mw_clear_sky_only) then
      if (runtime % opts % rt_all % ozone_data) then
         runtime % profiles(imem) % o3(:) = max(trace_gas % ozone(imem,lvlidx),0.0_r8)
      end if

      if (runtime % opts % rt_all % co2_data) then
         runtime % profiles(imem) % co2(:) = max(trace_gas % co2(imem,lvlidx),0.0_r8)
      end if

      if (runtime % opts % rt_all % n2o_data) then
         runtime % profiles(imem) % n2o(:) = max(trace_gas % n2o(imem,lvlidx),0.0_r8)
      end if

      if (runtime % opts % rt_all % ch4_data) then
         runtime % profiles(imem) % ch4(:) = max(trace_gas % ch4(imem,lvlidx),0.0_r8)
      end if

      if (runtime % opts % rt_all % co_data) then
         runtime % profiles(imem) % co(:)  = max(trace_gas % co(imem,lvlidx),0.0_r8)
      end if

      if (runtime % opts % rt_all % so2_data) then
         runtime % profiles(imem) % so2(:) = max(trace_gas % so2(imem,lvlidx),0.0_r8)
      end if

      ! "clear-sky" cloud-liquid water data for microwave (not RTTOV-SCATT)
      if (is_mw .and. runtime % opts % rt_mw % clw_data) then
         runtime % profiles(imem) % clw(:) = max(clouds % clw(imem,lvlidx),0.0_r8)
      end if

      ! add aerosols
      if (runtime % opts % rt_ir % addaerosl) then
         if (aerosl_type == 1) then 
            ! OPAC, convert from levels to layers
            runtime % profiles(imem) % aerosols(1,:) = max(0.5_jprb*(aerosols % insoluble(imem,ly1idx) +      &
                                                              aerosols % insoluble(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(2,:) = max(0.5_jprb*(aerosols % water_soluble(imem,ly1idx) +  &
                                                              aerosols % water_soluble(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(3,:) = max(0.5_jprb*(aerosols % soot(imem,ly1idx) +           &
                                                              aerosols % soot(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(4,:) = max(0.5_jprb*(aerosols % sea_salt_accum(imem,ly1idx) +      &
                                                              aerosols % sea_salt_accum(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(5,:) = max(0.5_jprb*(aerosols % sea_salt_coarse(imem,ly1idx) +     &
                                                              aerosols % sea_salt_coarse(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(6,:) = max(0.5_jprb*(aerosols % mineral_nucleus(imem,ly1idx) +     &
                                                              aerosols % mineral_nucleus(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(7,:) = max(0.5_jprb*(aerosols % mineral_accum(imem,ly1idx) +       &
                                                              aerosols % mineral_accum(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(8,:) = max(0.5_jprb*(aerosols % mineral_coarse(imem,ly1idx) +      &
                                                              aerosols % mineral_coarse(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(9,:) = max(0.5_jprb*(aerosols % mineral_transport(imem,ly1idx) +   &
                                                              aerosols % mineral_transport(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(10,:) = max(0.5_jprb*(aerosols % sulphated_droplets(imem,ly1idx) + &
                                                               aerosols % sulphated_droplets(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(11,:) = max(0.5_jprb*(aerosols % volcanic_ash(imem,ly1idx) +       &
                                                               aerosols % volcanic_ash(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(12,:) = max(0.5_jprb*(aerosols % new_volcanic_ash(imem,ly1idx) +   &
                                                               aerosols % new_volcanic_ash(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(13,:) = max(0.5_jprb*(aerosols % asian_dust(imem,ly1idx) +         &
                                                               aerosols % asian_dust(imem,ly2idx)),0.0_r8)
         elseif (aerosl_type == 2) then
            ! CAMS, convert from levels to levels
            runtime % profiles(imem) % aerosols(1,:) = max(0.5_jprb*(aerosols % black_carbon(imem,ly1idx) +               &
                                                              aerosols % black_carbon(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(2,:) = max(0.5_jprb*(aerosols % dust_bin1(imem,ly1idx) +                  &
                                                              aerosols % dust_bin1(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(3,:) = max(0.5_jprb*(aerosols % dust_bin2(imem,ly1idx) +                  &
                                                              aerosols % dust_bin2(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(4,:) = max(0.5_jprb*(aerosols % dust_bin3(imem,ly1idx) +                  &
                                                              aerosols % dust_bin3(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(5,:) = max(0.5_jprb*(aerosols % ammonium_sulphate(imem,ly1idx) +          &
                                                              aerosols % ammonium_sulphate(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(6,:) = max(0.5_jprb*(aerosols % sea_salt_bin1(imem,ly1idx) +              &
                                                              aerosols % sea_salt_bin1(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(7,:) = max(0.5_jprb*(aerosols % sea_salt_bin2(imem,ly1idx) +              &
                                                              aerosols % sea_salt_bin2(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(8,:) = max(0.5_jprb*(aerosols % sea_salt_bin3(imem,ly1idx) +              &
                                                              aerosols % sea_salt_bin3(imem,ly2idx)),0.0_r8)
            runtime % profiles(imem) % aerosols(9,:) = max(0.5_jprb*(aerosols % hydrophilic_organic_matter(imem,ly1idx) + &
                                                              aerosols % hydrophilic_organic_matter(imem,ly2idx)),0.0_r8)
         else
            write(string1,*)'Unknown aerosol type ',aerosl_type,&
               ' for platform/sat/sensor id combination:',instrument
            call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
         end if
      end if ! add aerosols

      ! add IR clouds
      if (runtime % opts % rt_ir % addclouds) then
         ! first find the total water and ice components

         totalwater(:) = 0.0_jprb

         if (allocated(clouds % clw)) then
            totalwater(:) = totalwater(:) + max(clouds % clw(imem,:),0.0_r8)
         end if

         if (allocated(clouds % rain)) then
            totalwater(:) = totalwater(:) + max(clouds % rain(imem,:),0.0_r8)
         end if

         totalice(:) = 0.0_r8

         if (allocated(clouds % ciw)) then
            totalice(:) = totalice(:) + max(clouds % ciw(imem,:),0.0_r8)
         end if 

         if (allocated(clouds % snow)) then
            totalice(:) = totalice(:) + max(clouds % snow(imem,:),0.0_r8)
         end if 

         if (allocated(clouds % graupel)) then
            totalice(:) = totalice(:) + max(clouds % graupel(imem,:),0.0_r8)
         end if 

         if (allocated(clouds % hail)) then
            totalice(:) = totalice(:) + max(clouds % hail(imem,:),0.0_r8)
         end if 

         ! Classify liquid water cloud type. If the maximum absolute w is > 0.5 m/s, classify as a cumulus cloud.
         ! FIXME: consider adding 0.5 as a namelist parameter
         if (allocated(clouds % w)) then
            maxw = maxval(abs(clouds % w(imem,:)))
            is_cumulus = maxw > 0.5_jprb ! m/s
         else
            ! assume cumulus if w not provided
            is_cumulus = .true. 
         end if

         ! depending on the vertical velocity and land type, classify clouds the way RTTOV wants 
         runtime % profiles(imem) % cloud(:,:) = 0.0_jprb
         if (.not. is_cumulus) then
            ! stratus
            if (surftype == 0) then
               ! 1, Stratus Continental, STCO (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(1,:) = 0.5_jprb*(totalwater(ly1idx) + &
                                                                 totalwater(ly2idx))
            else
               ! 2, Stratus Maritime, STMA (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(2,:) = 0.5_jprb*(totalwater(ly1idx) + &
                                                              totalwater(ly2idx))
            end if
         else
            ! cumulus
            if (surftype == 0) then
               ! 3, Cumulus Continental Clean, CUCC (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(3,:) = 0.5_jprb*(totalwater(ly1idx) + &
                                                              totalwater(ly2idx))
               ! FIXME: give user (or obs) control over CUCC versus CUCP
               ! 4, Cumulus Continental Polluted, CUCP (kg/kg), convert levels to layers
               ! runtime % profiles(imem) % cloud(4,:) = 0.5_jprb*(totalwater(ly1idx) + totalwater(ly2idx))
            else
               ! 5, Cumulus Maritime, CUMA (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(5,:) = 0.5_jprb*(totalwater(ly1idx) + &
                                                              totalwater(ly2idx))
            end if
         end if

         ! allow specification of cloud water effective diameter
         if (allocated(clouds % clwde)) then
            ! Liquid water effective diameter (microns), convert levels to layers
            runtime % profiles(imem)% clwde(:) = 0.5_jprb*(clouds % clwde(imem,ly1idx) + &
                                                           clouds % clwde(imem,ly2idx))
         end if

         ! all types of ice cloud (kg/kg) go in 6, Ice Cloud (CIRR, but not just cirrus)
         ! convert levels to layers
         runtime % profiles(imem) % cloud(6,:) = 0.5_jprb*(totalice(ly1idx) + &
                                                           totalice(ly2idx))
         
         ! allow specification of ice effective diameter
         if (allocated(clouds % icede)) then
            ! Ice effective diameter (microns), convert levels to layers
            runtime % profiles(imem) % icede(:) = 0.5_jprb*(clouds % icede(imem,ly1idx) + &
                                                            clouds % icede(imem,ly2idx))
         end if

         if (allocated(clouds % cfrac)) then
            ! Cloud fraction (0-1), convert levels to layers
            runtime % profiles(imem) % cfrac(:) = min(max(0.5_jprb*(clouds % cfrac(imem,ly1idx) + &
                                                         clouds % cfrac(imem,ly2idx)),0.0_r8),1.0_r8)
         else
            ! Assume cloud fraction is 1 everywhere. No harm if no clouds.
            runtime % profiles(imem) % cfrac(:) = 1.0_jprb
         end if
      end if ! add IR clouds
   else if (is_mw) then
      ! RTTOV-SCATT, add MW clouds

      ! Changes v12->v13: 
      ! variable `clouds % cc` removed
      ! Instead of named variables (e.g. rain, snow, hail, ...),
      ! a number of hydrometeor types can be specified in the hydro(nlevels,nhydro) member array (see guide)

      ! Removed variables: cc, clw, ciw, totalice, rain, sp
      ! New variables: hydro_frac(nlevels, nhydro_frac), hydro(nlevels, nhydro)

      ! nhydro = number of hydrometeor types (e.g. rain, snow, ...)
      ! nhydro_frac = 1 or nhydro

      if (allocated(clouds % cfrac) .and. runtime % opts_scatt % lusercfrac) then
         ! Use custom cfrac values
         ! TODO: specify cfrac (scalar?!)
         ! runtime % cld_profiles(imem) % cfrac = ? not implemented
      else
         ! normally calculated internally in RTTOV-SCATT
         runtime % cld_profiles(imem) % cfrac = -1
      end if

      ! cloud fraction per hydrometeor type 
      ! TODO: How do we get this from model data? From the 3D rain field?
      runtime % cld_profiles(imem) % hydro_frac(:,:) = 1.0_jprb

      ! code proposed, depends on the hydrotables of RTTOV?
      ! TODO: adapt to hydrotable, change indices of hydro(:,X) <--- here
      runtime % cld_profiles(imem) % hydro = 0.0_jprb
      if (allocated(clouds % clw)) then
         runtime % cld_profiles(imem) % hydro(:,0) = max(clouds % clw(imem,:),0.0_r8)
      endif
      if (allocated(clouds % rain)) then
         runtime % cld_profiles(imem) % hydro(:,1) = max(clouds % rain(imem,:),0.0_r8)
      endif
      if (allocated(clouds % ciw)) then
         runtime % cld_profiles(imem) % hydro(:,2) = max(clouds % ciw(imem,:),0.0_r8)
      endif
      if (allocated(clouds % snow)) then
         runtime % cld_profiles(imem) % hydro(:,3) = max(clouds % snow(imem,:),0.0_r8)
      endif
      if (allocated(clouds % hail)) then
         runtime % cld_profiles(imem) % hydro(:,4) = max(clouds % hail(imem,:),0.0_r8)
      endif

      ! also add "half-level pressures" as requested by RTTOV-Scatt
      runtime % cld_profiles(imem) % ph(2:nlevels) = 0.5_jprb*(atmos % pressure(imem,ly1idx)+atmos % pressure(imem,ly2idx))/100.0_jprb
      runtime % cld_profiles(imem) % ph(nlevels+1) = atmos % sfc_p(imem)/100.0_jprb
      runtime % cld_profiles(imem) % ph(1) = 0.0_jprb
   else
      write(string1,*)'Neither opts or opts_scatter were available for platform/sat/sensor id combination:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if ! has runtime % opts or runtime % opts_scatt

   ! 2m above surface pressure (hPa)
   runtime % profiles(imem) % s2m % p  = atmos % sfc_p(imem)/100.0_jprb ! convert from Pa to hPa
   ! 2m above surface temperature (K)
   runtime % profiles(imem) % s2m % t  = atmos % s2m_t(imem)

   ! set the optional variables
   if (allocated(atmos % s2m_q)) then
      ! 2m above surface water vapor (kg/kg)
      runtime % profiles(imem) % s2m % q  = atmos % s2m_q(imem)
   end if
 
   if (allocated(atmos % s10m_u) .and. allocated(atmos % s10m_v)) then
      ! 10 m above surface u and v wind (m/s)
      runtime % profiles(imem) % s2m % u  = atmos % s10m_u(imem)
      runtime % profiles(imem) % s2m % v  = atmos % s10m_v(imem)
   end if
 
   if (allocated(atmos % wfetch)) then
      ! Wind fetch over the ocean (m)
      runtime % profiles(imem) % s2m % wfetc = atmos % wfetch(imem)  
   else
      runtime % profiles(imem) % s2m % wfetc = wfetc_value
   end if
   
   ! Surface type (0=land, 1=sea, 2=sea-ice)
   runtime % profiles(imem) % skin % surftype  = surftype

   if (allocated(atmos % water_type)) then 
      ! Water type (0=fresh, 1=ocean)
      runtime % profiles(imem) % skin % watertype = nint(atmos % water_type(imem))
   end if

   ! Surface skin temperature (K) 
   runtime % profiles(imem) % skin % t = atmos % skin_temp(imem)

   if (allocated(atmos % sfc_salinity)) then
      ! Surface ocean salinity (practical salinity unit)  
      runtime % profiles(imem) % skin % salinity = atmos % sfc_salinity(imem)
   end if

   if (allocated(atmos % sfc_foam_frac)) then
      ! Surface foam fraction (0-1)
      runtime % profiles(imem) % skin % foam_fraction = min(max(atmos % sfc_foam_frac(imem),0.0_r8),1.0_r8)
   end if

   if (allocated(atmos % sfc_snow_frac)) then
      ! Surface snow fraction (0-1)
      runtime % profiles(imem) % skin % snow_fraction = min(max(atmos % sfc_snow_frac(imem),0.0_r8),1.0_r8)
   end if

   if (is_mw) then
      ! FASTEM parameters, see RTTOV user guide e.g. Table 21
      runtime % profiles(imem) % skin % fastem(1)   = mw_md%fastem_p1
      runtime % profiles(imem) % skin % fastem(2)   = mw_md%fastem_p2
      runtime % profiles(imem) % skin % fastem(3)   = mw_md%fastem_p3
      runtime % profiles(imem) % skin % fastem(4)   = mw_md%fastem_p4
      runtime % profiles(imem) % skin % fastem(5)   = mw_md%fastem_p5
   end if

   if (is_visir .and. do_lambertian) then
      ! Surface specularity (0-1)
      runtime % emissivity(imem) % specularity = min(max(visir_md % specularity,0.0_jprb),1.0_jprb)
   end if

   if (allocated(clouds % simple_cfrac)) then
      ! Simple (column) cloud fraction, 0-1
      runtime % profiles(imem) % cfraction = min(max(clouds % simple_cfrac(imem),0.0_r8),1.0_r8)
   end if

   if (allocated(clouds % ctp)) then
      ! Simple (column) cloud top pressure, hPa
      runtime % profiles(imem) % ctp = clouds % ctp(imem)
   end if

   if (is_visir) then
      ! Sat. zenith and azimuth angles (degrees)
      if (visir_md % sat_ze > max_zenith_angle) then
        radiances(:) = MISSING_R8
        error_status(:) = 101
        return
      end if
      runtime % profiles(imem) % zenangle    = max(visir_md % sat_ze,0.0_jprb)
      runtime % profiles(imem) % azangle     = visir_md % sat_az

      ! Solar zenith and azimuth angles (degrees), only relevant if use_solar
      runtime % profiles(imem) % sunzenangle = max(visir_md % sun_ze,0.0_jprb)
      runtime % profiles(imem) % sunazangle  = visir_md % sun_az
   else
      if (mw_md % sat_ze > max_zenith_angle) then
        radiances(:) = MISSING_R8
        error_status(:) = 101
        return
      end if
      ! Sat. zenith and azimuth angles (degrees)
      runtime % profiles(imem) % zenangle    = max(mw_md % sat_ze,0.0_jprb)
      runtime % profiles(imem) % azangle     = mw_md % sat_az
   end if

   ! Elevation (km), latitude and longitude (degrees)
   runtime % profiles(imem) % latitude  = lat
   runtime % profiles(imem) % longitude = lon
   runtime % profiles(imem) % elevation = atmos % sfc_elev(imem)/1000.0_jprb ! m -> km

   if (use_zeeman .and. is_mw) then
      runtime % profiles(imem) % Be    = mw_md % mag_field ! mag field strength, Gauss
      runtime % profiles(imem) % cosbk = mw_md % cosbk     ! cosine of angle between mag field and viewing angle
   end if
end do ! profile data

! --------------------------------------------------------------------------
! 6. Specify surface emissivity and reflectance
! --------------------------------------------------------------------------

! Here we assume zero specularity
runtime % emissivity(:) % specularity = 0._jprb

! Calculate emissivity within RTTOV where the input emissivity value is
! zero or less (all channels in this case)
runtime % calcemis(:) = (runtime % emissivity(:) % emis_in <= 0._jprb)

! In this example we have no values for input reflectances
runtime % reflectance(:) % refl_in = 0._jprb

! Calculate BRDF within RTTOV where the input BRDF value is zero or less
! (all channels in this case)
runtime % calcrefl(:) = (runtime % reflectance(:) % refl_in <= 0._jprb)

! Use default cloud top BRDF for simple cloud in VIS/NIR channels
runtime % reflectance(:) % refl_cloud_top = 0._jprb

if (debug) then
   print *,'is_visir:',is_visir
   call rttov_print_opts(runtime % opts, lu=ioout)
   call rttov_print_profile(runtime % profiles(1), lu=ioout)
   if (associated(runtime % opts_scatt)) then
      call rttov_print_opts_scatt(runtime % opts_scatt, lu=ioout)
   end if
end if


if (is_visir .or. mw_clear_sky_only) then
   ! Call RTTOV forward model
   call rttov_direct (                         &
           errorstatus,                        & ! out   error flag
           runtime % chanprof,                 & ! in    channel and profile index structure
           runtime % opts,                     & ! in    options structure
           runtime % profiles,                 & ! in    profile array
           runtime % coefs,                    & ! in    coefficients structure
           runtime % transmission,             & ! inout computed transmittances
           runtime % radiance,                 & ! inout computed radiances
           calcemis    = runtime % calcemis,   & ! in    flag for internal emissivity calcs
           emissivity  = runtime % emissivity, & ! inout input/output emissivities per channel
           calcrefl    = runtime % calcrefl,   & ! in    flag for internal BRDF calcs
           reflectance = runtime % reflectance ) ! inout input/output BRDFs per channel

   if (errorstatus /= errorstatus_success) then
      write(string1,*)'RTTOV direct error for platform/satellite/sensor id:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if
else
   ! set the frequencies based on which channel is being called
   runtime % frequencies(:) = runtime % frequencies_all(channel)

   if (debug) then
      call rttov_print_opts_scatt(runtime % opts_scatt, lu=ioout)
      call rttov_print_cld_profile(runtime % cld_profiles(1), lu=ioout)
   end if

   call rttov_scatt (          &
      errorstatus,            & ! out   error flag
      runtime % opts_scatt,   & ! in    RTTOV-SCATT options structure
      nlevels,                & ! in    number of profile levels
      runtime % chanprof,     & ! in    channel and profile index structure
      runtime % frequencies,  & ! in    channel indexes for Mietable lookup
      runtime % profiles,     & ! in    profile array
      runtime % cld_profiles, & ! in    cloud/hydrometeor profile array
      runtime % coefs,        & ! in    coefficients structure
      runtime % coefs_scatt,  & ! in    Mietable structure
      runtime % calcemis,     & ! in    flag for internal emissivity calcs
      runtime % emissivity,   & ! inout input/output emissivities per channel
      runtime % radiance )      ! inout computed radiances
   
   if (errorstatus /= errorstatus_success) then
      write(string1,*)'RTTOV scatt error for platform/satellite/sensor id:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if
end if

! utility from obs_kind_mod for getting string of obs qty
obs_qty_string = get_name_for_type_of_obs(flavor)
obs_type_num   = get_quantity_for_type_of_obs(flavor)

if (obs_type_num == QTY_RADIANCE) then
   do imem = 1, ens_size
      radiances(imem) = runtime % radiance % total(imem)
   end do
   if (debug) then
      print*, 'RADIANCE % TOTAL for ',trim(obs_qty_string),'= ', radiances(:)
   end if
elseif (obs_type_num == QTY_BRIGHTNESS_TEMPERATURE) then
   do imem = 1, ens_size
      radiances(imem) = runtime % radiance % bt(imem)
   end do
   if (debug) then
      print*, 'RADIANCE % BT for ',trim(obs_qty_string),'= ', radiances(:)
   end if
elseif (obs_type_num == QTY_BI_DIRECTIONAL_REFLECTANCE) then
   do imem = 1, ens_size
      radiances(imem) = runtime % radiance % refl(imem)
   end do
   if (debug) then
      print*, 'RADIANCE % REFL for ',trim(obs_qty_string),'= ', radiances(:)
   end if
else
   call error_handler(E_ERR, 'unknown observation quantity for ' // trim(obs_qty_string), &
      source, revision, revdate)
end if

IF (errorstatus /= errorstatus_success) THEN
  WRITE (*,*) 'rttov_direct error'
  !CALL rttov_exit(errorstatus)
ENDIF

end subroutine do_forward_model

!--------------------------------------------------------------------------

subroutine sensor_runtime_takedown(runtime, error_status)

type(rttov_sensor_runtime_type), pointer :: runtime
integer, intent(out) :: error_status

integer :: nchanprof
integer :: nlevels

! Return error status of RTTOV subroutine calls
integer(kind=jpim)               :: errorstatus 

!FIXME - in the forward operator code we won't be deallocating
! this structure.

! initialize error status tp success
error_status = errorstatus_success

! Deallocate all RTTOV arrays and structures

nchanprof = size(runtime % chanprof)
nlevels = size(runtime % profiles(1) % p)

call rttov_alloc_direct(               &
      errorstatus,                     &
      0_jpim,                          &  ! 0 => deallocate
      nchanprof,                       &  ! DART only does one channel at a time, so nchanprof = nprof
      nchanprof,                       &
      nlevels,                         &
      runtime % chanprof,              &
      runtime % opts,                  &
      runtime % profiles,              &
      runtime % coefs,                 &
      runtime % transmission,          &
      runtime % radiance,              &
      calcemis=runtime % calcemis,     &
      emissivity=runtime % emissivity, &
      calcrefl=runtime % calcrefl,     &
      reflectance=runtime % reflectance)

IF (errorstatus /= errorstatus_success) THEN
  WRITE(*,*) 'deallocation error for rttov_direct structures'
  error_status = errorstatus
  return
ENDIF

CALL rttov_dealloc_coefs(errorstatus, runtime % coefs)
IF (errorstatus /= errorstatus_success) THEN
  WRITE(*,*) 'coefs deallocation error'
ENDIF

end subroutine sensor_runtime_takedown

!----------------------------------------------------------------------------

subroutine initialize_module

character(len=*), parameter :: routine = 'initialize_module'

call register_module(source, revision, revdate)

module_initialized = .true.

missing_visir_metadata%sat_az       = MISSING_R8
missing_visir_metadata%sat_ze       = MISSING_R8
missing_visir_metadata%sun_az       = MISSING_R8
missing_visir_metadata%sun_ze       = MISSING_R8
missing_visir_metadata%platform_id  = MISSING_I
missing_visir_metadata%sat_id       = MISSING_I
missing_visir_metadata%sensor_id    = MISSING_I
missing_visir_metadata%channel      = MISSING_I
missing_visir_metadata%specularity  = MISSING_R8

missing_mw_metadata%sat_az       = MISSING_R8
missing_mw_metadata%sat_ze       = MISSING_R8
missing_mw_metadata%platform_id  = MISSING_I
missing_mw_metadata%sat_id       = MISSING_I
missing_mw_metadata%sensor_id    = MISSING_I
missing_mw_metadata%channel      = MISSING_I
missing_mw_metadata%mag_field    = MISSING_R8
missing_mw_metadata%cosbk        = MISSING_R8
missing_mw_metadata%fastem_p1 = MISSING_R8
missing_mw_metadata%fastem_p2 = MISSING_R8
missing_mw_metadata%fastem_p3 = MISSING_R8
missing_mw_metadata%fastem_p4 = MISSING_R8
missing_mw_metadata%fastem_p5 = MISSING_R8

allocate(obstype_metadata(2, MAXrttovkey))
allocate(visir_obs_metadata(MAXrttovkey))
allocate(mw_obs_metadata(MAXrttovkey))

obstype_metadata(:,:) = NO_OBS
visir_obs_metadata(:) = missing_visir_metadata
mw_obs_metadata(:)    = missing_mw_metadata

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_rttov_nml", iunit)
read(iunit, nml = obs_def_rttov_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_rttov_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_rttov_nml)
if (do_nml_term()) write(     *     , nml=obs_def_rttov_nml)

if (debug) then
   call error_handler(E_MSG,routine,'',source,revision,revdate)
   
   write(string1,*)'rttov_sensor_db_file   - ',trim(rttov_sensor_db_file)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'first_lvl_is_sfc       - ',first_lvl_is_sfc
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'mw_clear_sky_only      - ',mw_clear_sky_only
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'interp_mode            - ',interp_mode
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'do_checkinput          - ',do_checkinput
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'apply_reg_limits       - ',apply_reg_limits
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'verbose                - ',verbose
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'fix_hgpl               - ',fix_hgpl
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'do_lambertian          - ',do_lambertian
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'lambertian_fixed_angle - ',lambertian_fixed_angle
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'rad_down_lin_tau       - ',rad_down_lin_tau
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_q2m                - ',use_q2m
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_uv10m              - ',use_uv10m
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_wfetch             - ',use_wfetch
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_water_type         - ',use_water_type
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'addrefrac              - ',addrefrac
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'plane_parallel         - ',plane_parallel
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_salinity           - ',use_salinity
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cfrac_data             - ',cfrac_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'clw_data               - ',clw_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'rain_data              - ',rain_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ciw_data               - ',ciw_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'snow_data              - ',snow_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'graupel_data           - ',graupel_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'hail_data              - ',hail_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'w_data                 - ',w_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'clw_scheme             - ',clw_scheme
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'clw_cloud_top          - ',clw_cloud_top
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'fastem_version         - ',fastem_version
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'supply_foam_fraction   - ',supply_foam_fraction
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_totalice           - ',use_totalice
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_zeeman             - ',use_zeeman
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cc_threshold           - ',cc_threshold
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ozone_data             - ',ozone_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'co2_data               - ',co2_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'n2o_data               - ',n2o_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'co_data                - ',co_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ch4_data               - ',ch4_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'so2_data               - ',so2_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'addsolar               - ',addsolar
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'rayleigh_single_scatt  - ',rayleigh_single_scatt
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'do_nlte_correction     - ',do_nlte_correction
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'solar_sea_brdf_model   - ',solar_sea_brdf_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ir_sea_emis_model      - ',ir_sea_emis_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_sfc_snow_frac      - ',use_sfc_snow_frac
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'add_aerosl             - ',add_aerosl
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'aerosl_type            - ',aerosl_type
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'add_clouds             - ',add_clouds
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ice_scheme             - ',ice_scheme
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_icede              - ',use_icede
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'idg_scheme             - ',idg_scheme
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'user_aer_opt_param     - ',user_aer_opt_param
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'user_cld_opt_param     - ',user_cld_opt_param
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'grid_box_avg_cloud     - ',grid_box_avg_cloud
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cldcol_threshold       - ',cldcol_threshold
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cloud_overlap          - ',cloud_overlap
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cc_low_cloud_top   - ',cc_low_cloud_top
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ir_scatt_model         - ',ir_scatt_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'vis_scatt_model        - ',vis_scatt_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'dom_nstreams           - ',dom_nstreams
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'dom_accuracy           - ',dom_accuracy
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'dom_opdep_threshold    - ',dom_opdep_threshold
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'addpc                  - ',addpc
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'npcscores              - ',npcscores
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'addradrec              - ',addradrec
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ipcreg                 - ',ipcreg
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_htfrtc             - ',use_htfrtc
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'htfrtc_n_pc            - ',htfrtc_n_pc
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'htfrtc_simple_cloud    - ',htfrtc_simple_cloud
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'htfrtc_overcast        - ',htfrtc_overcast
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if

call read_sensor_db_file(rttov_sensor_db_file)

end subroutine initialize_module

!----------------------------------------------------------------------
! Initialize a RTTOV sensor runtime. A rttov_sensor_type instance contains 
! information such as options and coefficients that are initialized in a "lazy"
! fashion only when it will be used for the first time.

subroutine initialize_rttov_sensor_runtime(sensor,ens_size,nlevels)

type(rttov_sensor_type), pointer    :: sensor
integer,                 intent(in) :: ens_size
integer,                 intent(in) :: nlevels

type(rttov_options),       pointer  :: opts
type(rttov_options_scatt), pointer  :: opts_scatt
logical :: is_mw
logical :: is_vis
logical :: is_ir
logical :: is_visir

if (.not. module_initialized) then
   call initialize_module()
end if

if (.not. associated(sensor)) then
   write(string1,*) "The sensor must be initialized before initializing a rttov sensor runtime"
   call error_handler(E_ERR, 'initialize_rttov', string1, source, revision, revdate)
end if

! some logicals to ease in setting the options
is_mw  = sensor%sensor_type_name == 'mw'
is_vis = sensor%sensor_type_name == 'vis'
is_ir  = sensor%sensor_type_name == 'ir'
is_visir = is_vis .or. is_ir

allocate(opts)
! these options have descriptions above, except for where defaults are set outside of user-control
opts % interpolation % addinterp        = .true.      ! Allow interpolation of input profile
opts % interpolation % reg_limit_extrap = .true.      ! Extrapolate beyond top of model intelligently 
opts % interpolation % interp_mode      = interp_mode    
opts % interpolation % lgradp           = .false.     ! Do not allow TL/AD/K of user pressure levels
opts % interpolation % spacetop         = .true.      ! Treat model top as space boundary (not recommended to be false)

opts % config % do_checkinput    = do_checkinput    
opts % config % apply_reg_limits = apply_reg_limits 
opts % config % verbose          = verbose         
opts % config % fix_hgpl         = fix_hgpl

opts % rt_all % switchrad              = is_mw                 ! switch to BT for AD/K if microwave
opts % rt_all % do_lambertian          = do_lambertian
opts % rt_all % lambertian_fixed_angle = lambertian_fixed_angle
opts % rt_all % rad_down_lin_tau       = rad_down_lin_tau
opts % rt_all % use_q2m                = use_q2m
opts % rt_all % addrefrac              = addrefrac
opts % rt_all % plane_parallel         = plane_parallel

if (is_mw) then
   ! mw options could be used with direct or scatt, so set them here
   opts % rt_mw % clw_data              = clw_data
   opts % rt_mw % clw_scheme            = clw_scheme
   opts % rt_mw % clw_cloud_top         = clw_cloud_top
   opts % rt_mw % fastem_version        = fastem_version
   opts % rt_mw % supply_foam_fraction  = supply_foam_fraction
end if

if (is_mw .and. .not. mw_clear_sky_only) then
   ! use RTTOV-SCATT 
   allocate(opts_scatt)

   opts_scatt % config % do_checkinput    = do_checkinput
   opts_scatt % config % apply_reg_limits = apply_reg_limits
   opts_scatt % config % verbose          = verbose
   opts_scatt % config % fix_hgpl         = fix_hgpl

   opts_scatt % rad_down_lin_tau      = rad_down_lin_tau
   opts_scatt % interp_mode           = interp_mode
   opts_scatt % lgradp                = .false.               ! Do not allow TL/AD/K of user pressure levels
   opts_scatt % reg_limit_extrap      = .true.                ! intelligently extend beyond the model top
   opts_scatt % fastem_version        = fastem_version
   opts_scatt % supply_foam_fraction  = supply_foam_fraction
   opts_scatt % use_q2m               = use_q2m
   opts_scatt % lusercfrac            = .false.               ! Have RTTOV-SCATT calculate the effective cloud fraction 
   opts_scatt % cc_threshold          = cc_threshold

   call sensor_runtime_setup(sensor,                &
                             ens_size=ens_size,     &
                             nlevs=nlevels,         &
                             opts=opts,             &
                             opts_scatt=opts_scatt, &
                             use_totalice=use_totalice)
else
   opts % rt_all % ozone_data            = ozone_data
   opts % rt_all % co2_data              = co2_data
   opts % rt_all % n2o_data              = n2o_data
   opts % rt_all % co_data               = co_data
   opts % rt_all % ch4_data              = ch4_data
   opts % rt_all % so2_data              = so2_data
   opts % rt_ir % addsolar              = addsolar
   opts % rt_ir % rayleigh_single_scatt = rayleigh_single_scatt
   opts % rt_ir % do_nlte_correction    = do_nlte_correction
   opts % rt_ir % solar_sea_brdf_model  = solar_sea_brdf_model
   opts % rt_ir % ir_sea_emis_model     = ir_sea_emis_model
   opts % rt_ir % addaerosl             = add_aerosl
   opts % rt_ir % addclouds             = add_clouds
   opts % rt_ir % user_aer_opt_param    = user_aer_opt_param
   opts % rt_ir % user_cld_opt_param    = user_cld_opt_param
   opts % rt_ir % grid_box_avg_cloud    = grid_box_avg_cloud
   opts % rt_ir % cldcol_threshold      = cldcol_threshold
   opts % rt_ir % cloud_overlap         = cloud_overlap
   opts % rt_ir % cc_low_cloud_top      = cc_low_cloud_top
   opts % rt_ir % ir_scatt_model        = ir_scatt_model
   opts % rt_ir % vis_scatt_model       = vis_scatt_model
   opts % rt_ir % dom_nstreams          = dom_nstreams
   opts % rt_ir % dom_accuracy          = dom_accuracy
   opts % rt_ir % dom_opdep_threshold   = dom_opdep_threshold

   opts % rt_ir % pc % addpc     = addpc
   opts % rt_ir % pc % npcscores = npcscores
   opts % rt_ir % pc % addradrec = addradrec
   opts % rt_ir % pc % ipcbnd    = 1           ! default, currently the only valid value
   opts % rt_ir % pc % ipcreg    = ipcreg

   opts % htfrtc_opts % htfrtc        = use_htfrtc
   opts % htfrtc_opts % n_pc_in       = htfrtc_n_pc
   opts % htfrtc_opts % reconstruct   = .true. ! By default reconstruct radiances
   opts % htfrtc_opts % simple_cloud  = htfrtc_simple_cloud
   opts % htfrtc_opts % overcast      = htfrtc_overcast

   call sensor_runtime_setup(sensor,                &
                             ens_size=ens_size,     &
                             nlevs=nlevels,         &
                             opts=opts)
end if



end subroutine initialize_rttov_sensor_runtime

!----------------------------------------------------------------------
! Fill the module storage metadata for a particular visible/infrared 
! observation.
!
! Visible / infrared observations have several auxillary metadata variables.
! Other than the key, which is standard DART fare, the RTTOV satellite azimuth 
! and satellite zenith angle must be specified. See the RTTOV user guide for 
! more information (in particular, see figure 4). If the *addsolar*
! namelist value is set to true, then the solar azimuth and solar zenith angles must
! be specified - again see the RTTOV user guide. In addition to the platform/satellite/
! sensor ID numbers, which are the RTTOV unique identifiers, the channel specifies
! the chanenl number in the RTTOV coefficient file. Finally, if *do_lambertian* is 
! true, specularity must be specified here. Again, see the RTTOV user guide for more 
! information.

subroutine set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
   platform_id, sat_id, sensor_id, channel, specularity)

integer,  intent(out) :: key
real(r8), intent(in)  :: sat_az, sat_ze
real(r8), intent(in)  :: sun_az, sun_ze ! only relevant if addsolar
integer,  intent(in)  :: platform_id, sat_id, sensor_id, channel
real(r8), intent(in)  :: specularity    ! only relevant if do_lambertian

if ( .not. module_initialized ) call initialize_module

visirnum = visirnum + 1  ! increase module storage used counters
rttovkey = rttovkey + 1

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(rttovkey,'set_visir_metadata', VISIR)

key = rttovkey ! now that we know its legal
obstype_metadata(SUBTYPE, key) = VISIR 
obstype_metadata(SUBKEY, key) = visirnum

visir_obs_metadata(visirnum) % sat_az      = sat_az
visir_obs_metadata(visirnum) % sat_ze      = sat_ze
visir_obs_metadata(visirnum) % sun_az      = sun_az
visir_obs_metadata(visirnum) % sun_ze      = sun_ze
visir_obs_metadata(visirnum) % platform_id = platform_id
visir_obs_metadata(visirnum) % sat_id      = sat_id
visir_obs_metadata(visirnum) % sensor_id   = sensor_id
visir_obs_metadata(visirnum) % channel     = channel
visir_obs_metadata(visirnum) % specularity = specularity

end subroutine set_visir_metadata

!----------------------------------------------------------------------
! Fill the module storage metadata for a particular microwave 
! observation.
! Microwave observations have several auxillary metadata variables.
! Other than the key, which is standard DART fare, the RTTOV satellite azimuth 
! and satellite zenith angle must be specified. See the RTTOV user guide for 
! more information (in particular, see figure 4). In addition to the platform/satellite/
! sensor ID numbers, which are the RTTOV unique identifiers, the channel specifies
! the chanenl number in the RTTOV coefficient file. In addition, if 
! use_zeeman is true, the magnetic field and cosine of the angle
! between the magnetic field and angle of propagation must be specified. 
! See the RTTOV user guide for more information. Finally, the fastem parameters for
! land must be specified here. This may be difficult for observations to set, so 
! default values (see table 21 in the RTTOV user guide) can be used until a better 
! solution is devised. 

subroutine set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, &
   channel, mag_field, cosbk, fastem_p1, fastem_p2, fastem_p3,         &
   fastem_p4, fastem_p5)

integer,  intent(out) :: key
real(r8), intent(in)  :: sat_az, sat_ze
integer,  intent(in)  :: platform_id, sat_id, sensor_id, channel
real(r8), intent(in)  :: mag_field, cosbk                            ! only relevant with use_zeeman
real(r8), intent(in)  :: fastem_p1, fastem_p2, fastem_p3, & ! 
                         fastem_p4, fastem_p5 !5 FASTEM parameters describing the land/sea ice surface, see Table 21 (MW, used in FASTEM and for land/sea ice for TESSM2)

if ( .not. module_initialized ) call initialize_module

mwnum    = mwnum + 1    ! increase module storage used counters
rttovkey = rttovkey + 1 

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(rttovkey,'set_mw_metadata', MW)

key = rttovkey ! now that we know its legal
obstype_metadata(SUBTYPE, key) = MW    
obstype_metadata(SUBKEY, key) = mwnum


mw_obs_metadata(mwnum) % sat_az       = sat_az
mw_obs_metadata(mwnum) % sat_ze       = sat_ze
mw_obs_metadata(mwnum) % platform_id  = platform_id
mw_obs_metadata(mwnum) % sat_id       = sat_id
mw_obs_metadata(mwnum) % sensor_id    = sensor_id
mw_obs_metadata(mwnum) % channel      = channel
mw_obs_metadata(mwnum) % mag_field    = mag_field
mw_obs_metadata(mwnum) % cosbk        = cosbk
mw_obs_metadata(mwnum) % fastem_p1    = fastem_p1
mw_obs_metadata(mwnum) % fastem_p2    = fastem_p2
mw_obs_metadata(mwnum) % fastem_p3    = fastem_p3
mw_obs_metadata(mwnum) % fastem_p4    = fastem_p4
mw_obs_metadata(mwnum) % fastem_p5    = fastem_p5

end subroutine set_mw_metadata


!----------------------------------------------------------------------
! Query the metadata in module storage for a particular observation.
! See set_visir_metadata for more information on these fields.

subroutine get_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
   platform_id, sat_id, sensor_id, channel, specularity)

integer,  intent(in)  :: key
real(r8), intent(out) :: sat_az, sat_ze, sun_az, sun_ze
integer,  intent(out) :: platform_id, sat_id, sensor_id, channel
real(r8), intent(out) :: specularity

integer :: mykey

character(len=*), parameter :: routine = 'get_visir_metadata'

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the useful length of the metadata arrays.
call key_within_range(key,routine)

if (obstype_metadata(SUBTYPE,key) /= VISIR) then
   write(string1,*)'The obstype metadata for key ',key,'is not visir as expected'
   call error_handler(E_ERR, routine, string1, source, &
      revision, revdate)
end if

mykey = obstype_metadata(SUBKEY,key)

if (.not. is_valid_subkey(mykey, VISIR)) then
   write(string1,*)'The visir-specific key ',mykey,'is invalid.'
   write(string2,*)'Size of visir_obs_metadata:', visirnum
   call error_handler(E_ERR, routine, string1, source, &
      revision, revdate, text2=string2)
end if

sat_az      = visir_obs_metadata(mykey) % sat_az
sat_ze      = visir_obs_metadata(mykey) % sat_ze
sun_az      = visir_obs_metadata(mykey) % sun_az
sun_ze      = visir_obs_metadata(mykey) % sun_ze
platform_id = visir_obs_metadata(mykey) % platform_id
sat_id      = visir_obs_metadata(mykey) % sat_id
sensor_id   = visir_obs_metadata(mykey) % sensor_id
channel     = visir_obs_metadata(mykey) % channel
specularity = visir_obs_metadata(mykey) % specularity

end subroutine get_visir_metadata

!----------------------------------------------------------------------
! Query the metadata in module storage for a particular observation.
! See set_mw_metadata for more information on these fields.

subroutine get_mw_metadata(key, sat_az, sat_ze,        &
   platform_id, sat_id, sensor_id, channel, mag_field, &
   cosbk, fastem_p1, fastem_p2, fastem_p3,    &
   fastem_p4, fastem_p5)

integer,  intent(in)  :: key
real(r8), intent(out) :: sat_az, sat_ze
integer,  intent(out) :: platform_id, sat_id, sensor_id, channel
real(r8), intent(out) :: mag_field, cosbk
real(r8), intent(out) :: fastem_p1, fastem_p2, fastem_p3, &
                         fastem_p4, fastem_p5

character(len=*), parameter :: routine = 'get_mw_metadata'

integer :: mykey

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the useful length of the metadata arrays.
call key_within_range(key,routine)

if (obstype_metadata(SUBTYPE,key) /= MW) then
   write(string1,*)'The obstype metadata for key ',key,'is not MW as expected'
   call error_handler(E_ERR, routine, string1, source, &
      revision, revdate)
end if

mykey = obstype_metadata(SUBKEY, key)

if (.not. is_valid_subkey(mykey, MW)) then
   write(string1,*)'The mw-specific key ',mykey,'is invalid.'
   write(string2,*)'Size of mw_obs_metadata:',mwnum 
   call error_handler(E_ERR, routine, string1, source, &
      revision, revdate, text2=string2)
end if

sat_az       = mw_obs_metadata(mykey) % sat_az
sat_ze       = mw_obs_metadata(mykey) % sat_ze
platform_id  = mw_obs_metadata(mykey) % platform_id
sat_id       = mw_obs_metadata(mykey) % sat_id
sensor_id    = mw_obs_metadata(mykey) % sensor_id
channel      = mw_obs_metadata(mykey) % channel
mag_field    = mw_obs_metadata(mykey) % mag_field
cosbk        = mw_obs_metadata(mykey) % cosbk
fastem_p1    = mw_obs_metadata(mykey) % fastem_p1
fastem_p2    = mw_obs_metadata(mykey) % fastem_p2
fastem_p3    = mw_obs_metadata(mykey) % fastem_p3
fastem_p4    = mw_obs_metadata(mykey) % fastem_p4
fastem_p5    = mw_obs_metadata(mykey) % fastem_p5

end subroutine get_mw_metadata


!----------------------------------------------------------------------
! This routine reads the metadata for the visir/mw radiance/tb obs
subroutine read_rttov_metadata(key, obsID, ifile, fform)
integer,          intent(out)          :: key    ! index into local metadata
integer,          intent(in)           :: obsID
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
logical           :: is_asciifile
integer           :: ierr
character(len=5)  :: header
integer           :: oldkey
real(r8)          :: sat_az, sat_ze, sun_az, sun_ze
integer           :: platform_id, sat_id, sensor_id, channel
real(r8)          :: specularity
real(r8)          :: mag_field, cosbk
real(r8)          :: fastem_p1, fastem_p2, fastem_p3, &
                     fastem_p4, fastem_p5

logical :: is_visir

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

write(string2,*)'observation #',obsID

if ( is_asciifile ) then
   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_rttov_metadata','header',string2)
   if (trim(header) == trim(VISIR_STRING)) then
      ! Load the visible/IR data from the file
      read(ifile, *, iostat=ierr) sat_az, sat_ze, sun_az, sun_ze
      call check_iostat(ierr,'read_rttov_metadata','sat,sun az/ze',string2)
      read(ifile, *, iostat=ierr) platform_id, sat_id, sensor_id, channel
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, *, iostat=ierr) specularity
      call check_iostat(ierr,'read_rttov_metadata','specularity',string2)
      read(ifile, *, iostat=ierr) oldkey
      call check_iostat(ierr,'read_rttov_metadata','oldkey',string2)
      is_visir = .true.
   elseif (trim(header) == trim(MW_STRING)) then
      ! Load the visible/IR data from the file
      read(ifile, *, iostat=ierr) sat_az, sat_ze
      call check_iostat(ierr,'read_rttov_metadata','sat az/ze',string2)
      read(ifile, *, iostat=ierr) platform_id, sat_id, sensor_id, channel
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, *, iostat=ierr) mag_field, cosbk
      call check_iostat(ierr,'read_rttov_metadata','mag_field/cosbk',string2)
      read(ifile, *, iostat=ierr) fastem_p1, fastem_p2, fastem_p3, &
                                  fastem_p4, fastem_p5
      call check_iostat(ierr,'read_rttov_metadata','fastem_p1-5',string2)
      read(ifile, *, iostat=ierr) oldkey
      call check_iostat(ierr,'read_rttov_metadata','oldkey',string2)
      is_visir = .false.
   else
      write(string1,*)"Expected radiance header ["//visir_STRING//"] or ["//mw_STRING//"]" // &
           " in input file, got ["//header//"]"
      call error_handler(E_ERR, 'read_rttov_metadata', string1, source, revision, revdate, text2=string2)
   endif
else
   read(ifile, iostat=ierr) header
   call  check_iostat(ierr,'read_rttov_metadata','header',string2)

   if (trim(header) == trim(VISIR_STRING)) then
      ! Load the visible/IR data from the file
      read(ifile, iostat=ierr) sat_az, sat_ze, sun_az, sun_ze
      call check_iostat(ierr,'read_rttov_metadata','sat,sun az/ze',string2)
      read(ifile, iostat=ierr) platform_id, sat_id, sensor_id, channel
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, iostat=ierr) specularity
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, iostat=ierr) oldkey
      call check_iostat(ierr,'read_rttov_metadata','oldkey',string2)
      is_visir = .true.
   elseif (trim(header) == trim(MW_STRING)) then
      ! Load the visible/IR data from the file
      read(ifile, iostat=ierr) sat_az, sat_ze
      call check_iostat(ierr,'read_rttov_metadata','sat az/ze',string2)
      read(ifile, iostat=ierr) platform_id, sat_id, sensor_id, channel
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, iostat=ierr) mag_field, cosbk
      call check_iostat(ierr,'read_rttov_metadata','mag_field/cosbk',string2)
      read(ifile, iostat=ierr) fastem_p1, fastem_p2, fastem_p3, &
                                  fastem_p4, fastem_p5
      call check_iostat(ierr,'read_rttov_metadata','fastem_p1-5',string2)
      read(ifile, iostat=ierr) oldkey
      call check_iostat(ierr,'read_rttov_metadata','oldkey',string2)
      is_visir = .false.
   else
      write(string1,*)"Expected radiance header ["//visir_STRING//"] or ["//mw_STRING//"]" // &
         " in input file, got ["//header//"]"
      call error_handler(E_ERR, 'read_rttov_metadata', string1, &
         source, revision, revdate, text2=string2)
   end if
endif

! The oldkey is thrown away.

! Store the metadata in module storage. The new key is returned.
if (is_visir) then
   call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
      platform_id, sat_id, sensor_id, channel, specularity)
else
   call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, &
      channel, mag_field, cosbk, fastem_p1, fastem_p2, fastem_p3,   &
      fastem_p4, fastem_p5)
end if

end subroutine read_rttov_metadata


!----------------------------------------------------------------------
! writes the metadata for radiance observations.

subroutine write_rttov_metadata(key, ifile, fform)

integer,           intent(in)           :: key
integer,           intent(in)           :: ifile
character(len=*),  intent(in), optional :: fform

logical  :: is_asciifile
real(r8) :: sat_az, sat_ze, sun_az, sun_ze
integer  :: platform_id, sat_id, sensor_id, channel
real(r8) :: specularity
real(r8) :: mag_field, cosbk
real(r8) :: fastem_p1, fastem_p2, fastem_p3, &
            fastem_p4, fastem_p5

character(len=5)  :: header
character(len=*), parameter :: routine = 'write_rttov_metadata'

if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation.

is_asciifile = ascii_file_format(fform)

select case (obstype_metadata(SUBTYPE,key))

   case (VISIR)
      call get_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
         platform_id, sat_id, sensor_id, channel, specularity)
   
      header = VISIR_STRING
   
      if (is_asciifile) then
         write(ifile, *) header
         write(ifile, *) sat_az, sat_ze, sun_az, sun_ze
         write(ifile, *) platform_id, sat_id, sensor_id, channel
         write(ifile, *) specularity
         write(ifile, *) key
      else
         write(ifile   ) header
         write(ifile   ) sat_az, sat_ze, sun_az, sun_ze
         write(ifile   ) platform_id, sat_id, sensor_id, channel
         write(ifile   ) specularity
         write(ifile   ) key
      endif
   case (MW)
      call get_mw_metadata(key, sat_az, sat_ze,                  &
         platform_id, sat_id, sensor_id, channel, mag_field, cosbk, &
         fastem_p1, fastem_p2, fastem_p3, fastem_p4,    &
         fastem_p5)
   
      header = MW_STRING
   
      if (is_asciifile) then
         write(ifile, *) header
         write(ifile, *) sat_az, sat_ze
         write(ifile, *) platform_id, sat_id, sensor_id, channel
         write(ifile, *) mag_field, cosbk
         write(ifile, *) fastem_p1, fastem_p2, fastem_p3, &
                         fastem_p4, fastem_p5 
         write(ifile, *) key
      else
         write(ifile   ) header
         write(ifile   ) sat_az, sat_ze
         write(ifile   ) platform_id, sat_id, sensor_id, channel
         write(ifile   ) mag_field, cosbk
         write(ifile   ) fastem_p1, fastem_p2, fastem_p3, &
                         fastem_p4, fastem_p5 
         write(ifile   ) key
      end if
   case default
      call error_handler(E_ERR, routine, 'unknown metadata type', source, revision, revdate)
end select

end subroutine write_rttov_metadata


!----------------------------------------------------------------------
subroutine interactive_rttov_metadata(key)
integer, intent(out) :: key

real(r8)          :: sat_az, sat_ze, sun_az, sun_ze
integer           :: platform_id, sat_id, sensor_id, channel
real(r8)          :: specularity
real(r8)          :: mag_field, cosbk
real(r8)          :: fastem_p1, fastem_p2, fastem_p3, &
                     fastem_p4, fastem_p5

integer :: visir_int
logical :: is_visir

if ( .not. module_initialized ) call initialize_module

! Prompt for input for the required metadata
visir_int = interactive_i('freq type    Visible/IR or MW observation? [1 = vis/ir, 0 = microwave]', minvalue=0, maxvalue=1)

if (visir_int == 0) then
   is_visir = .false.
else
   is_visir = .true.
end if

sat_az   = interactive_r('sat_az        satellite azimuth [degrees]', minvalue = 0.0_r8, maxvalue = 360.0_r8)
sat_ze   = interactive_r('sat_ze        satellite zenith [degrees]',  minvalue = 0.0_r8, maxvalue = 90.0_r8)
if (is_visir .and. addsolar) then
   sun_az   = interactive_r('sun_az        solar azimuth [degrees]',     minvalue = 0.0_r8, maxvalue = 360.0_r8)
   sun_ze   = interactive_r('sun_ze        solar zenith [degrees]',      minvalue = 0.0_r8, maxvalue = 90.0_r8)
else
   sun_az = MISSING_R8
   sun_ze = MISSING_R8
end if

platform_id = interactive_i('platform id   RTTOV Platform number [see docs]',     minvalue = 1)
sat_id      = interactive_i('satellite id  RTTOV Satellite ID number [see docs]', minvalue = 0)
sensor_id   = interactive_i('sensor id     RTTOV Sensor number [see docs]',       minvalue = 1)
channel     = interactive_i('channel       Instrument channel number [see docs]', minvalue = 1)

if (is_visir .and. do_lambertian) then
   specularity = interactive_r('specularity   (0-1)',     minvalue = 0.0_r8, maxvalue = 1.0_r8)
else
   specularity = MISSING_R8
end if

if (.not. is_visir .and. use_zeeman) then
   mag_field   = interactive_r('magnetic fld  Earth magnetic field strength (Gauss)',                    minvalue = 0.2_r8,  maxvalue = 0.7_r8)
   cosbk       = interactive_r('cosbk         Cosine of angle between magnetic field and viewing angle', minvalue = -1.0_r8, maxvalue = 1.0_r8)
else
   mag_field   = MISSING_R8
   cosbk       = MISSING_R8
end if

if (.not. is_visir) then
   fastem_p1 = interactive_r('fastem land1  1st fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_p2 = interactive_r('fastem land2  2nd fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_p3 = interactive_r('fastem land3  3rd fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_p4 = interactive_r('fastem land4  4th fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_p5 = interactive_r('fastem land5  5th fastem land/sea ice parameter', minvalue = 0.0_r8)
else
   fastem_p1 = MISSING_R8
   fastem_p2 = MISSING_R8
   fastem_p3 = MISSING_R8
   fastem_p4 = MISSING_R8
   fastem_p5 = MISSING_R8
end if

if (is_visir) then
   call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
      platform_id, sat_id, sensor_id, channel, specularity)
else
   call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, &
      channel, mag_field, cosbk, fastem_p1, fastem_p2, fastem_p3,   &
      fastem_p4, fastem_p5)
end if

end subroutine interactive_rttov_metadata


!----------------------------------------------------------------------
subroutine get_expected_radiance(obs_kind_ind, state_handle, ens_size, location, key, flavor, val, istatus)

integer,             intent(in)  :: obs_kind_ind
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location          ! location of obs
integer,             intent(in)  :: key               ! key into module metadata
integer,             intent(in)  :: flavor            ! flavor of obs
real(r8),            intent(out) :: val(ens_size)     ! value of obs
integer,             intent(out) :: istatus(ens_size) ! status of the calculation

integer  :: platform_id, sat_id, sensor_id, channel
integer  :: instrument(3)

integer :: this_istatus(ens_size)

integer  :: i
real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat
real(r8) :: loc_value(ens_size)
type(location_type) :: loc
integer :: maxlevels, numlevels

type(location_type) :: loc_undef        ! surface location
type(location_type) :: loc_sea          ! surface location
type(location_type) :: loc_sfc          ! surface location
type(location_type) :: loc_2m           ! surface location
type(location_type) :: loc_10m          ! surface location

type(rttov_sensor_type), pointer :: sensor

character(len=*), parameter :: routine = 'get_expected_radiance'

type(visir_metadata_type), pointer :: visir_md
type(mw_metadata_type),    pointer :: mw_md

logical :: return_now

!=================================================================================

if ( .not. module_initialized ) call initialize_module

val = 0.0_r8 ! set return value early

! Make sure the desired key is within the useful length of the metadata arrays.
call key_within_range(key, routine)

select case (obstype_metadata(SUBTYPE, key))

case (VISIR)
   visir_md => visir_obs_metadata(obstype_metadata(SUBKEY,key))
   mw_md    => null()

   platform_id = visir_md % platform_id
   sat_id      = visir_md % sat_id
   sensor_id   = visir_md % sensor_id
   channel     = visir_md % channel
case (MW)
   visir_md => null()
   mw_md    => mw_obs_metadata(obstype_metadata(SUBKEY,key))

   platform_id = mw_md % platform_id
   sat_id      = mw_md % sat_id
   sensor_id   = mw_md % sensor_id
   channel     = mw_md % channel
case default
   call error_handler(E_ERR, routine, 'unknown metadata type', source, revision, revdate)
end select

!=================================================================================
! Determine the number of model levels 
! using only the standard DART interfaces to model
!=================================================================================

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

! if necessary, preallocate the column arrays
if ( .not. arrays_prealloced) then
   numlevels = 0

   maxlevels = 10000   ! something larger than we think will exist
   COUNTLEVELS : do i = 1,maxlevels
      loc = set_location(loc_lon, loc_lat, real(i,r8), VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc, QTY_PRESSURE, loc_value, this_istatus)
      if ( any(this_istatus /= 0 ) ) exit COUNTLEVELS
      numlevels = numlevels + 1
   enddo COUNTLEVELS

   if ((numlevels == maxlevels) .or. (numlevels == 0)) then
      write(string1,'(A,I0)') 'FAILED to determine number of levels in model:', &
         numlevels
         
      if (debug) call error_handler(E_ALLMSG,routine,string1,source,revision,revdate)
      istatus = 1
      val     = MISSING_R8
      return
   else
       if (debug) write(*,*) routine // ' we have ',numlevels,' model levels'
   endif

   ! allocate the necessary fields
   call atmos_profile_setup(atmos, ens_size, numlevels, use_q2m, & 
      use_uv10m, use_wfetch, use_water_type, use_salinity,       &
      supply_foam_fraction, use_sfc_snow_frac)

   call trace_gas_profile_setup(trace_gas, ens_size, numlevels,  &
      ozone_data, co2_data, n2o_data, ch4_data, co_data, so2_data)

   if (add_aerosl) then
      call aerosol_profile_setup(aerosols, ens_size, numlevels,  &
         aerosl_type)
   end if

   call cloud_profile_setup(clouds, ens_size, numlevels,     &
      cfrac_data, clw_data, clw_scheme, rain_data, ciw_data, &
      ice_scheme, use_icede, snow_data, graupel_data,        &
      hail_data, w_data, htfrtc_simple_cloud)

   arrays_prealloced = .true.
else
   ! load the number of levels from the temperature array, assumed to have been initialized
   numlevels = size(atmos%temperature,2)
end if 

instrument(1) = platform_id
instrument(2) = sat_id
instrument(3) = sensor_id

! also check if the sensor runtime needs to be initialized
sensor => get_rttov_sensor(instrument)

if (.not. associated(sensor)) then
   write(string1,*)'Could not find the platform/sat/sensor id combination:',&
      platform_id,sat_id,sensor_id
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

if (.not. associated(sensor%runtime)) then
   call initialize_rttov_sensor_runtime(sensor,ens_size,numlevels)
end if

! Set all of the istatuses back to zero for track_status
istatus = 0

GETLEVELDATA : do i = 1,numlevels
   loc = set_location(loc_lon, loc_lat, real(i,r8), VERTISLEVEL)

   call interpolate(state_handle, ens_size, loc, QTY_PRESSURE, atmos%pressure(:,i), this_istatus)
   call check_status('QTY_PRESSURE', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .true., return_now)
   if (return_now) return

   call interpolate(state_handle, ens_size, loc, QTY_TEMPERATURE, atmos%temperature(:, i), this_istatus)
   call check_status('QTY_TEMPERATURE', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .true., return_now)
   if (return_now) return

   call interpolate(state_handle, ens_size, loc, QTY_VAPOR_MIXING_RATIO, atmos%moisture(:, i), this_istatus)
   call check_status('QTY_VAPOR_MIXING_RATIO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .true., return_now)
   if (return_now) return

   if (ozone_data) then
      call interpolate(state_handle, ens_size, loc, QTY_O3, trace_gas%ozone(:, i), this_istatus)
      call check_status('QTY_O3', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (co2_data) then
      call interpolate(state_handle, ens_size, loc, QTY_CO2, trace_gas%co2(:, i), this_istatus)
      call check_status('QTY_CO2', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (n2o_data) then
      call interpolate(state_handle, ens_size, loc, QTY_N2O, trace_gas%n2o(:, i), this_istatus)
      call check_status('QTY_N2O', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (ch4_data) then
      call interpolate(state_handle, ens_size, loc, QTY_CH4, trace_gas%ch4(:, i), this_istatus)
      call check_status('QTY_CH4', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (co_data) then
      call interpolate(state_handle, ens_size, loc, QTY_CO, trace_gas%co(:, i), this_istatus)
      call check_status('QTY_CO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (so2_data) then
      call interpolate(state_handle, ens_size, loc, QTY_SO2, trace_gas%so2(:, i), this_istatus)
      call check_status('QTY_SO2', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (add_aerosl) then
      if (aerosl_type == 1) then
         ! OPAC
         call interpolate(state_handle, ens_size, loc, QTY_INSOLUBLE_AER, aerosols%insoluble(:, i), this_istatus)
         call check_status('QTY_INSOLUBLE_AER', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_H2O_SOLUBLE_AER, aerosols%water_soluble(:, i), this_istatus)
         call check_status('QTY_H2O_SOLUBLE_AER', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_SOOT, aerosols%soot(:, i), this_istatus)
         call check_status('QTY_SOOT', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_SEASALT_ACCUM, aerosols%sea_salt_accum(:, i), this_istatus)
         call check_status('QTY_SEASALT_ACCUM', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_SEASALT_COARSE, aerosols%sea_salt_coarse(:, i), this_istatus)
         call check_status('QTY_SEASALT_COARSE', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_NUCLEUS, aerosols%mineral_nucleus(:, i), this_istatus)
         call check_status('QTY_MINERAL_NUCLEUS', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_ACCUM, aerosols%mineral_accum(:, i), this_istatus)
         call check_status('QTY_MINERAL_ACCUM', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_COARSE, aerosols%mineral_coarse(:, i), this_istatus)
         call check_status('QTY_MINERAL_COARSE', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_TRANSPORTED, aerosols%mineral_transport(:, i), this_istatus)
         call check_status('QTY_MINERAL_TRANSPORT', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_SULPHATED_DROPS, aerosols%sulphated_droplets(:, i), this_istatus)
         call check_status('QTY_SULPHATED_DROPS', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_VOLCANIC_ASH, aerosols%volcanic_ash(:, i), this_istatus)
         call check_status('QTY_VOLCANIC_ASH', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_NEW_VOLCANIC_ASH, aerosols%new_volcanic_ash(:, i), this_istatus)
         call check_status('QTY_NEW_VOLCANIC_ASH', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_ASIAN_DUST, aerosols%asian_dust(:, i), this_istatus)
         call check_status('QTY_ASIAN_DUST', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return
      elseif (aerosl_type == 2) then
         ! CAMS
         call interpolate(state_handle, ens_size, loc, QTY_BLACK_CARBON, aerosols%black_carbon(:, i), this_istatus)
         call check_status('QTY_BLACK_CARBON', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_DUST_BIN1, aerosols%dust_bin1(:, i), this_istatus)
         call check_status('QTY_DUST_BIN1', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_DUST_BIN2, aerosols%dust_bin2(:, i), this_istatus)
         call check_status('QTY_DUST_BIN2', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_DUST_BIN3, aerosols%dust_bin3(:, i), this_istatus)
         call check_status('QTY_DUST_BIN3', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_AMMONIUM_SULPHATE, aerosols%ammonium_sulphate(:, i), this_istatus)
         call check_status('QTY_AMMONIUM_SULPHATE', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_SEA_SALT_BIN1, aerosols%sea_salt_bin1(:, i), this_istatus)
         call check_status('QTY_SEA_SALT_BIN1', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_SEA_SALT_BIN2, aerosols%sea_salt_bin2(:, i), this_istatus)
         call check_status('QTY_SEA_SALT_BIN2', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_SEA_SALT_BIN3, aerosols%sea_salt_bin3(:, i), this_istatus)
         call check_status('QTY_SEA_SALT_BIN3', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return

         call interpolate(state_handle, ens_size, loc, QTY_HYDROPHILIC_ORGANIC_MATTER, aerosols%hydrophilic_organic_matter(:, i), this_istatus)
         call check_status('QTY_HYDROPHILIC_ORGANIC_MATTER', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return
      end if
   end if

   if (w_data) then
      ! specify vertical velocity
      call interpolate(state_handle, ens_size, loc, QTY_VERTICAL_VELOCITY, clouds%w(:, i), this_istatus)
      call check_status('QTY_VERTICAL_VELOCITY', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (cfrac_data) then
      ! specify cloud fraction
      call interpolate(state_handle, ens_size, loc, QTY_CLOUD_FRACTION, clouds%cfrac(:, i), this_istatus)
      call check_status('QTY_CLOUD_FRACTION', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   ! clwde should be specified when clw_scheme == 2 (takes particle diameter from model); clw_scheme = 1 would parametrize diameters depending on cloud type
   if (clw_scheme == 2) then
      ! The effective diameter must also be specified with clw_scheme 2
      ! call interpolate(state_handle, ens_size, loc, QTY_CLOUDWATER_DE, clouds%clwde(:, i), this_istatus)
      !clouds%clwde(:, i) = 2*1e6*clouds%clwde(:, i)  ! convert from WRF variable radius in m to DART diameter in micrometer
      call check_status('QTY_CLOUDWATER_DE', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (clw_data) then
      ! specify non-precip cloud liquid water ! for MW (but not used by RTTOV-SCATT anymore since v13)
      call interpolate(state_handle, ens_size, loc, QTY_CLOUDWATER_MIXING_RATIO, clouds%clw(:, i), this_istatus)
      call check_status('QTY_CLOUDWATER_MIXING_RATIO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (rain_data) then
      ! specify precip cloud liquid water (i.e. rain)
      call interpolate(state_handle, ens_size, loc, QTY_RAINWATER_MIXING_RATIO, clouds%rain(:, i), this_istatus)
      call check_status('QTY_RAINWATER_MIXING_RATIO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (ciw_data) then
      ! specify non-precip cloud ice
      call interpolate(state_handle, ens_size, loc, QTY_ICE_MIXING_RATIO, clouds%ciw(:, i), this_istatus)
      call check_status('QTY_ICE_MIXING_RATIO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return

      if (ice_scheme == 1 .and. use_icede) then
         ! if use_icede with ice_scheme 1, must also specify ice effective diameter
         !call interpolate(state_handle, ens_size, loc, QTY_CLOUD_ICE_DE, clouds%icede(:, i), this_istatus)
         !clouds%icede(:, i) = 2*1e6*clouds%icede(:, i)  ! convert from WRF variable radius in m to DART diameter in micrometer
         call check_status('QTY_CLOUD_ICE_DE', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
         if (return_now) return
      end if
   end if

   if (snow_data) then
      ! specify precip fluffy ice (i.e. snow)
      call interpolate(state_handle, ens_size, loc, QTY_SNOW_MIXING_RATIO, clouds%snow(:, i), this_istatus)
      call check_status('QTY_SNOW_MIXING_RATIO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (graupel_data) then
      ! specify precip soft-hail (i.e. graupel)
      call interpolate(state_handle, ens_size, loc, QTY_GRAUPEL_MIXING_RATIO, clouds%graupel(:, i), this_istatus)
      call check_status('QTY_GRAUPEL_MIXING_RATIO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

   if (hail_data) then
      ! specify precip hard-hail (i.e. hail)
      call interpolate(state_handle, ens_size, loc, QTY_HAIL_MIXING_RATIO, clouds%hail(:, i), this_istatus)
      call check_status('QTY_HAIL_MIXING_RATIO', ens_size, this_istatus, val, loc, istatus, routine, source, revision, revdate, .false., return_now)
      if (return_now) return
   end if

end do GETLEVELDATA

! Set the surface fields
! Some models check that the difference between 'surface' locations and 
! the surface of the model is within some namelist-specified tolerance.
! To ensure that these interpolations return useful values, get the model
! definiiton of the surface and add the appropriate height, if any.

! Simply check if we can get the model surface elevation (loc_sea is not used)
loc_sea = set_location(loc_lon, loc_lat, 0.0_r8, VERTISSURFACE )
call interpolate(state_handle, ens_size, loc_sea, QTY_SURFACE_ELEVATION, atmos%sfc_elev(:), this_istatus)
call check_status('QTY_SURFACE_ELEVATION', ens_size, this_istatus, val, loc_sea, istatus, routine, source, revision, revdate, .true., return_now)
if (return_now) return

!>@todo FIXME If the model does not check for surface consistency,
!>            should we continue ...

! These are the locations of interest. 
loc_undef = set_location(loc_lon, loc_lat,                MISSING_R8, VERTISUNDEF )
loc_sfc   = set_location(loc_lon, loc_lat,  0.0_r8+atmos%sfc_elev(1), VERTISSURFACE )
loc_2m    = set_location(loc_lon, loc_lat,  2.0_r8+atmos%sfc_elev(1), VERTISSURFACE )
loc_10m   = set_location(loc_lon, loc_lat, 10.0_r8+atmos%sfc_elev(1), VERTISSURFACE )

call interpolate(state_handle, ens_size, loc_sfc, QTY_SURFACE_PRESSURE, atmos%sfc_p(:), this_istatus)
call check_status('QTY_SURFACE_PRESSURE', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .true., return_now)
if (return_now) return

call interpolate(state_handle, ens_size, loc_2m, QTY_2M_TEMPERATURE, atmos%s2m_t(:), this_istatus)
call check_status('QTY_2M_TEMPERATURE', ens_size, this_istatus, val, loc_2m, istatus, routine, source, revision, revdate, .true., return_now)
if (return_now) return

call interpolate(state_handle, ens_size, loc_sfc, QTY_SKIN_TEMPERATURE, atmos%skin_temp(:), this_istatus)
call check_status('QTY_SKIN_TEMPERATURE', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .true., return_now)
if (return_now) return

! set to 2m if an error

call interpolate(state_handle, ens_size, loc_sfc, QTY_SURFACE_TYPE, atmos%surftype(:), this_istatus)
call check_status('QTY_SURFACE_TYPE', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .true., return_now)
if (return_now) return

! if not available, lookup by lat/lon?

if (use_q2m) then
   call interpolate(state_handle, ens_size, loc_2m, QTY_2M_SPECIFIC_HUMIDITY, atmos%s2m_q(:), this_istatus)
   call check_status('QTY_2M_SPECIFIC_HUMIDITY', ens_size, this_istatus, val, loc_2m, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (use_uv10m) then
   call interpolate(state_handle, ens_size, loc_10m, QTY_10M_U_WIND_COMPONENT, atmos%s10m_u(:), this_istatus)
   call check_status('QTY_10M_U_WIND_COMPONENT', ens_size, this_istatus, val, loc_10m, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
   call interpolate(state_handle, ens_size, loc_10m, QTY_10M_V_WIND_COMPONENT, atmos%s10m_v(:), this_istatus)
   call check_status('QTY_10M_V_WIND_COMPONENT', ens_size, this_istatus, val, loc_10m, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (use_wfetch) then
   call interpolate(state_handle, ens_size, loc_sfc, QTY_WIND_FETCH, atmos%wfetch(:), this_istatus)
   call check_status('QTY_WIND_FETCH', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (use_water_type) then
   call interpolate(state_handle, ens_size, loc_sfc, QTY_WATER_TYPE, atmos%water_type(:), this_istatus)
   call check_status('QTY_WATER_TYPE', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (use_salinity) then
   call interpolate(state_handle, ens_size, loc_sfc, QTY_SALINITY, atmos%sfc_salinity(:), this_istatus)
   call check_status('QTY_SALINITY', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (supply_foam_fraction) then
   call interpolate(state_handle, ens_size, loc_sfc, QTY_FOAM_FRAC, atmos%sfc_foam_frac(:), this_istatus)
   call check_status('QTY_FOAM_FRAC', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (use_sfc_snow_frac) then
   call interpolate(state_handle, ens_size, loc_sfc, QTY_SNOWCOVER_FRAC, atmos%sfc_snow_frac(:), this_istatus)
   call check_status('QTY_SNOWCOVER_FRAC', ens_size, this_istatus, val, loc_sfc, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (add_clouds .and. htfrtc_simple_cloud) then
   ! specify simple cloud information - per column
   call interpolate(state_handle, ens_size, loc_undef, QTY_COLUMN_CLOUD_FRAC, clouds%simple_cfrac(:), this_istatus)
   call check_status('QTY_COLUMN_CLOUD_FRAC', ens_size, this_istatus, val, loc_undef, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return

   call interpolate(state_handle, ens_size, loc_undef, QTY_CLOUD_TOP_PRESSURE, clouds%ctp(:), this_istatus)
   call check_status('QTY_CLOUD_TOP_PRESSURE', ens_size, this_istatus, val, loc_undef, istatus, routine, source, revision, revdate, .false., return_now)
   if (return_now) return
end if

if (debug) then
   print*, 'interpolate pressure    = ', atmos%pressure(1,1),    '...', atmos%pressure(1,numlevels)
   print*, 'interpolate temperature = ', atmos%temperature(1,1), '...', atmos%temperature(1,numlevels)
   print*, 'interpolate moisture    = ', atmos%moisture(1,1),    '...', atmos%moisture(1,numlevels)
end if

call do_forward_model(ens_size=ens_size,                    &
                      nlevels=numlevels,                    &
                      flavor=flavor,                        &
                      location=location,                    &
                      atmos=atmos,                          &
                      trace_gas=trace_gas,                  &
                      clouds=clouds,                        &
                      aerosols=aerosols,                    &
                      sensor=sensor,                        &
                      channel=channel,                      &
                      first_lvl_is_sfc=first_lvl_is_sfc,    &
                      mw_clear_sky_only=mw_clear_sky_only,  &
                      clw_scheme=clw_scheme,                &
                      ice_scheme=ice_scheme,                &
                      idg_scheme=idg_scheme,                &
                      aerosl_type=aerosl_type,              &
                      do_lambertian=do_lambertian,          &
                      use_totalice=use_totalice,            &
                      use_zeeman=use_zeeman,                &
                      radiances=val,                        &
                      error_status=this_istatus,            &
                      visir_md=visir_md,                    &
                      mw_md=mw_md) 

! copy the status from this_istatus to istatus, set missing if error
call track_status(ens_size, this_istatus, val, istatus, return_now)

if (debug) then
   print*, 'istatus  = ', istatus 
   print*, 'radiance = ', val
end if

end subroutine get_expected_radiance


!----------------------------------------------------------------------
! this is a fatal error routine.  use with care - if a forward
! operator fails it should return a bad error code, not die.
! useful in read/write routines where it is fatal to fail.

subroutine check_iostat(istat, routine, varname, msgstring)
integer,          intent(in) :: istat
character(len=*), intent(in) :: routine
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: msgstring

if ( istat /= 0 ) then
   write(string1,*)'istat should be 0 but is ',istat,' for '//varname
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=msgstring)
end if

end subroutine check_iostat


!----------------------------------------------------------------------
! Make sure we are addressing within the metadata arrays

subroutine key_within_range(key, routine)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

if (is_valid_key(key)) then
   ! we are still within limits
   return
else
   ! Bad news. Tell the user.
   write(string1, *) 'key (',key,') not within known range ( 1,', rttovkey,')'
   call error_handler(E_ERR,routine,string1,source,revision,revdate) 
end if

end subroutine key_within_range

!----------------------------------------------------------------------
! Make sure the key is within the useful range of the metdata arrays
function is_valid_key(key)

integer, intent(in) :: key
logical :: is_valid_key

is_valid_key = (key >0 .and. key <= rttovkey)

end function is_valid_key
!-----------------------------------------------------------------------
! Make sure the subkey is within the useful range of the metadata arrays
function is_valid_subkey(skey, stype)

integer, intent(in) :: skey  ! subkey 
integer, intent(in) :: stype ! subype (visir or mw)
logical :: is_valid_subkey

if (skey <=0 ) then
  is_valid_subkey = .false.
  return
endif

select case (stype)
case (VISIR)
  is_valid_subkey = (skey <= visirnum)
case (MW)
  is_valid_subkey = (skey <= mwnum)
case default
  is_valid_subkey = .false.
end select 

end function is_valid_subkey
 
!----------------------------------------------------------------------
! If the allocatable metadata arrays are not big enough ... try again

subroutine grow_metadata(key, routine, obstype)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine
integer,          intent(in) :: obstype

integer :: orglength, newlength
type(visir_metadata_type), allocatable :: safe_visir_metadata(:)
type(mw_metadata_type),    allocatable :: safe_mw_metadata(:)
integer,                   allocatable :: safe_obstype_metadata(:,:)

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key > 2*size(obstype_metadata(1,:))) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

select case (obstype)
   case (VISIR)
      if (visirnum > size(visir_obs_metadata)) then
         ! we need to grow visir_obs_metadata
         orglength = size(visir_obs_metadata)
         newlength = 2 * orglength
   
         ! News. Tell the user we are increasing storage.
         write(string1, *) 'Warning: subkey (',visirnum,') exceeds visir_obs_metadata length (',orglength,')'
         write(string2, *) 'Increasing visir_obs_metadata to length ',newlength
         call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
   
         allocate(safe_visir_metadata(orglength))
         safe_visir_metadata(:) = visir_obs_metadata(:)
   
         deallocate(visir_obs_metadata)
         allocate(visir_obs_metadata(newlength))
   
         visir_obs_metadata = missing_visir_metadata
         visir_obs_metadata(1:orglength) = safe_visir_metadata(:)
         deallocate(safe_visir_metadata)
     endif
   case (MW)
      ! duplicate the above since we can't use any object-oriented magic
      if (mwnum > size(mw_obs_metadata)) then
         ! we need to grow mw_obs_metadata
         orglength = size(mw_obs_metadata)
         newlength = 2 * orglength
   
         ! News. Tell the user we are increasing storage.
         write(string1, *) 'Warning: subkey (',mwnum,') exceeds mw_obs_metadata length (',orglength,')'
         write(string2, *) 'Increasing mw_obs_metadata to length ',newlength
         call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
   
         allocate(safe_mw_metadata(orglength))
         safe_mw_metadata(:) = mw_obs_metadata(:)
   
         deallocate(mw_obs_metadata)
         allocate(mw_obs_metadata(newlength))
   
         mw_obs_metadata = missing_mw_metadata
         mw_obs_metadata(1:orglength) = safe_mw_metadata(:)
         deallocate(safe_mw_metadata)
      end if
   case default
      call error_handler(E_ERR, routine, 'unknown metadata type', source, revision, revdate)
end select


if ( key > size(obstype_metadata(1,:)) ) then
   ! we need to grow obstype_metadata
   orglength = size(obstype_metadata(1,:))
   newlength = 2 * orglength

   ! News. Tell the user we are increasing storage.
   write(string1, *) 'Warning: key (',key,') exceeds obs_metadata length (',orglength,')'
   write(string2, *) 'Increasing obs_metadata to length ',newlength
   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
  
   allocate(safe_obstype_metadata(2,orglength))
   safe_obstype_metadata(:,:) = obstype_metadata(:,:)
   deallocate(obstype_metadata)
   allocate(obstype_metadata(2, newlength)) 

   obstype_metadata(:,:) = NO_OBS
   obstype_metadata(:,1:orglength) = safe_obstype_metadata(:,1:orglength)
   deallocate(safe_obstype_metadata) 
endif

end subroutine grow_metadata

!----------------------------------------------------------------------
!> override track_status behavior to always error and print an error message if the field
!> cannot be found

subroutine check_status(field_name, ens_size, this_istatus, val, loc, istatus, &
   routine, source, revision, revdate, required, return_now)

character(len=*),    intent(in)  :: field_name
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: this_istatus(ens_size)
real(r8),            intent(out) :: val(ens_size)
type(location_type), intent(in)  :: loc
integer,             intent(out) :: istatus(ens_size) ! status of the calculation
character(len=*),    intent(in)  :: routine
character(len=*),    intent(in)  :: source
character(len=*),    intent(in)  :: revision
character(len=*),    intent(in)  :: revdate
logical,             intent(in)  :: required
logical,             intent(out) :: return_now

real(r8) :: locv(3)

call track_status(ens_size, this_istatus, val, istatus, return_now)

if (return_now) then
   locv = get_location(loc)
   if (required) then
      write(string1,*) 'Could not find required field ' // trim(field_name), ' istatus:',istatus,&
         'location:',locv(1),'/',locv(2),'/',locv(3)
      call error_handler(E_WARN,routine,string1,source,revision,revdate)
   else if (debug) then
      write(string1,*) 'Could not find requested field ' // trim(field_name), ' istatus:',istatus,&
         'location:',locv(1),'/',locv(2),'/',locv(3)
      call error_handler(E_ALLMSG,routine,string1,source,revision,revdate)
   end if
end if

end subroutine check_status

!----------------------------------------------------------------------
! Return the logical value of the RTTOV parameter associated with the field_name.

function get_rttov_option_logical(field_name) result(p)
   character(*), intent(in) :: field_name

   logical :: p

   integer,          parameter   :: duc = ichar('A') - ichar('a')
   character(len=:), allocatable :: fname
   
   integer   :: slen


   ! copy the string over to an appropriate size
   slen = len_trim(field_name)
   allocate(character(len=slen) :: fname)
   fname = trim(adjustl(field_name))
   call to_upper(fname)

   select case (fname)
      case('FIRST_LVL_IS_SFC')
         p = FIRST_LVL_IS_SFC
      case('MW_CLEAR_SKY_ONLY')
         p = MW_CLEAR_SKY_ONLY
      case('DO_CHECKINPUT')
         p = DO_CHECKINPUT
      case('APPLY_REG_LIMITS')
         p = APPLY_REG_LIMITS
      case('VERBOSE')
         p = VERBOSE
      case('FIX_HGPL')
         p = FIX_HGPL
      case('DO_LAMBERTIAN')
         p = DO_LAMBERTIAN
      case('LAMBERTIAN_FIXED_ANGLE')
         p = LAMBERTIAN_FIXED_ANGLE
      case('RAD_DOWN_LIN_TAU')
         p = RAD_DOWN_LIN_TAU
      case('USE_Q2M')
         p = USE_Q2M
      case('USE_UV10M')
         p = USE_UV10M
      case('USE_WFETCH')
         p = USE_WFETCH
      case('USE_WATER_TYPE')
         p = USE_WATER_TYPE
      case('ADDREFRAC')
         p = ADDREFRAC
      case('PLANE_PARALLEL')
         p = PLANE_PARALLEL
      case('USE_SALINITY')
         p = USE_SALINITY
      case('CFRAC_DATA')
         p = CFRAC_DATA
      case('CLW_DATA')
         p = CLW_DATA
      case('RAIN_DATA')
         p = RAIN_DATA
      case('CIW_DATA')
         p = CIW_DATA
      case('SNOW_DATA')
         p = SNOW_DATA
      case('GRAUPEL_DATA')
         p = GRAUPEL_DATA
      case('HAIL_DATA')
         p = HAIL_DATA
      case('W_DATA')
         p = W_DATA
      case('SUPPLY_FOAM_FRACTION')
         p = SUPPLY_FOAM_FRACTION
      case('USE_TOTALICE')
         p = USE_TOTALICE
      case('USE_ZEEMAN')
         p = USE_ZEEMAN
      case('OZONE_DATA')
         p = OZONE_DATA
      case('CO2_DATA')
         p = CO2_DATA
      case('N2O_DATA')
         p = N2O_DATA
      case('CO_DATA')
         p = CO_DATA
      case('CH4_DATA')
         p = CH4_DATA
      case('SO2_DATA')
         p = SO2_DATA
      case('ADDSOLAR')
         p = ADDSOLAR
      case('RAYLEIGH_SINGLE_SCATT')
         p = RAYLEIGH_SINGLE_SCATT
      case('DO_NLTE_CORRECTION')
         p = DO_NLTE_CORRECTION
      case('USE_SFC_SNOW_FRAC')
         p = USE_SFC_SNOW_FRAC
      case('ADD_AEROSL')
         p = ADD_AEROSL
      case('ADD_CLOUDS')
         p = ADD_CLOUDS
      case('USE_ICEDE')
         p = USE_ICEDE
      case('USER_AER_OPT_PARAM')
         p = USER_AER_OPT_PARAM
      case('USER_CLD_OPT_PARAM')
         p = USER_CLD_OPT_PARAM
      case('GRID_BOX_AVG_CLOUD')
         p = GRID_BOX_AVG_CLOUD
      case('ADDPC')
         p = ADDPC
      case('ADDRADREC')
         p = ADDRADREC
      case('USE_HTFRTC')
         p = USE_HTFRTC
      case('HTFRTC_SIMPLE_CLOUD')
         p = HTFRTC_SIMPLE_CLOUD
      case('HTFRTC_OVERCAST')
         p = HTFRTC_OVERCAST
      case default
         write(string1,*) "Unknown logical field ",fname
         call error_handler(E_ERR, 'get_rttov_option_logical', string1, source, revision, revdate)
   end select

end function get_rttov_option_logical

!-----------------------------------------------------------------------
! A function to return the key associated with a VISIR/MW metadata.
! Note that this function (and the module storage of VISIR/MW metadata)
! should be rethought to allow for the notion of a "channel" to
! live more harmoniously with other types of observations.
!
! FIXME: reintegrate radiance_obs_seq_to_netcdf and obs_seq_to_netcdf

function get_channel(flavor, key) result(channel)

   integer, intent(in) :: flavor
   integer, intent(in) :: key
   integer             :: channel

   real(r8) :: sat_az, sat_ze, sun_az, sun_ze
   integer  :: platform_id, sat_id, sensor_id
   real(r8) :: mag_field, cosbk, specularity
   real(r8) :: fastem_p1, fastem_p2, fastem_p3, fastem_p4, fastem_p5

   character(len=*), parameter :: routine = 'get_channel'
   channel = MISSING_R8

   ! If the observation is not supported by this module, there is no channel
   ! This is delicate in that all types supported by this module are consecutively
   ! numbered. If new types are added, this will need to change.
   if (flavor < NOAA_1_VTPR1_RADIANCE .or. flavor > CLOUDSAT_1_CPR_TB) return

   ! Retrieve channel from different metadata types
   ! All the other arguments are mandatory, but not needed here.

   select case (obstype_metadata(SUBTYPE,key))

   case ( VISIR )
      call get_visir_metadata(key, &
                           sat_az, sat_ze, sun_az, sun_ze, &
                           platform_id, sat_id, sensor_id, channel, &
                           specularity)
   case ( MW )
      call get_mw_metadata(key, &
                        sat_az, sat_ze, &
                        platform_id, sat_id, sensor_id, channel, &
                        mag_field, cosbk, &
                        fastem_p1, fastem_p2, fastem_p3, fastem_p4, fastem_p5)
   case default
      call error_handler(E_ERR, routine, 'unknown metadata type', source, revision, revdate)
   end select

   if (debug) write(*,*)'get_channel: key,channel ',key,channel

end function get_channel


!-----------------------------------------------------------------------
! Test functions below this point
!-----------------------------------------------------------------------
function test_unit_setup(metadata_start)

integer, intent(in) :: metadata_start
logical :: test_unit_setup

! Test requires module not initialized to start
if ( module_initialized ) then
  test_unit_setup = .false.
else
  test_unit_setup = .true.
endif

! Overwrite module globals for testing
! start at 1 to watch the metadata grow 
MAXrttovkey = metadata_start
call initialize_module

end function test_unit_setup

!-----------------------------------------------------------------------
! Using a teardown for unit tests as there is no end_module
subroutine test_unit_teardown()

module_initialized = .false.
MAXrttovkey = 0
rttovkey = 0
visirnum = 0
mwnum = 0

deallocate(obstype_metadata, visir_obs_metadata, mw_obs_metadata)

end subroutine test_unit_teardown

!-----------------------------------------------------------------------
! inputs: n_ir : number of visir obs
!         n_mw : number of microwave obs
! outputs metadata_size[3] :  size(1,obstype_metadata)
!                             size(visir_obs_metadata)
!                             size(mw_obs_metadata)
!  Loop to have a look at how grow_metadata works
!  The actul data in the obs is junk
function test_set_metadata(n_ir, n_mw) result(metadata_size)

integer, intent(in)  :: n_ir, n_mw
integer  :: metadata_size(3) 

integer  :: key 
real(r8) :: sat_az, sat_ze, sun_az, sun_ze
integer  :: platform_id, sat_id, sensor_id, channel
real(r8) :: specularity 

real(r8) ::  mag_field, cosbk     
real(r8) :: fastem_p1, fastem_p2, fastem_p3, fastem_p4, fastem_p5 

integer :: ii

do ii = 1, n_ir
  call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
        platform_id, sat_id, sensor_id, channel, specularity)
enddo

do ii = 1, n_mw
  call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, &
     channel, mag_field, cosbk, fastem_p1, fastem_p2, fastem_p3,         &
     fastem_p4, fastem_p5)
enddo

metadata_size(1) = size(obstype_metadata,2)
metadata_size(2) = size(visir_obs_metadata)
metadata_size(3) = size(mw_obs_metadata)

end function test_set_metadata

!-----------------------------------------------------------------------
! passes metadata out to external test routine
subroutine test_metadata(metadata)

integer, intent(out) :: metadata(:,:)
metadata = obstype_metadata

end subroutine test_metadata

!-----------------------------------------------------------------------
! call is_valid_key for a bunch of values
subroutine test_key_within_range(keysin,inrange)

integer, intent(in)  :: keysin(:)
logical, intent(out) :: inrange(:)

integer :: ii, key

do ii = 1, size(keysin)
   key = keysin(ii)
   inrange(ii) = is_valid_key(key)
enddo

end subroutine test_key_within_range
!-----------------------------------------------------------------------
subroutine test_subkey_within_range(subkeysin,stype,inrange)

integer, intent(in)  :: subkeysin(:)
integer, intent(in)  :: stype
logical, intent(out) :: inrange(:)

integer :: ii, key

do ii = 1, size(subkeysin)
   key = subkeysin(ii)
   inrange(ii) = is_valid_subkey(key,stype)
enddo

end subroutine test_subkey_within_range
!-----------------------------------------------------------------------
end module obs_def_rttov_mod

! END DART PREPROCESS MODULE CODE
