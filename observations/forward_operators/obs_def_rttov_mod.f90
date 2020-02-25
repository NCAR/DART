! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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
!----------------------------------------------------------------------

! BEGIN DART PREPROCESS KIND LIST
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
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_rttov_mod, only : read_rttov_metadata, &
!                                write_rttov_metadata, &
!                          interactive_rttov_metadata, &
!                                get_expected_radiance
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(NOAA_1_VTPR1_RADIANCE:CLOUDSAT_1_CPR_TB)
!         call get_expected_radiance(obs_kind_ind, state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
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

use        types_mod, only : r8, PI, metadatalength, MISSING_R8, MISSING_I
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             logfileunit, get_unit, open_file, close_file, nc_check, &
                             file_exist, ascii_file_format, nmlfileunit, do_nml_file, &
                             do_nml_term, check_namelist_read, find_namelist_in_file
use     location_mod, only : location_type, set_location, get_location, VERTISUNDEF, &
                             VERTISHEIGHT, VERTISLEVEL, set_location_missing
use     obs_kind_mod

use  assim_model_mod, only : interpolate

use obs_def_utilities_mod, only : track_status
use ensemble_manager_mod,  only : ensemble_type

use rttov_interface_mod,   only : visir_metadata_type,             &
                                  mw_metadata_type,                &
                                  rttov_sensor_type,               &
                                  rttov_sensor_runtime_type,       &
                                  atmos_profile_type,              &
                                  aerosol_profile_type,            &
                                  cloud_profile_type,              &
                                  trace_gas_profile_type,          &
                                  get_rttov_sensor,                &
                                  read_sensor_db_file,             &
                                  sensor_runtime_setup,            &
                                  do_forward_model,                &
                                  sensor_runtime_takedown            ! unused at present?

use parkind1,              only : jpim, jprb, jplm

use rttov_types, only :    &
      rttov_options,       &
      rttov_options_scatt

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

public ::            set_visir_metadata, &
                     set_mw_metadata,    &
                     get_visir_metadata, &
                     get_mw_metadata,    &
                    read_rttov_metadata, &
                   write_rttov_metadata, &
             interactive_rttov_metadata, &
                  get_expected_radiance

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2
logical            :: module_initialized = .false.
logical            :: arrays_prealloced  = .false.
integer            :: iunit, rc

logical,                   allocatable :: obstype_metadata(:)
integer,                   allocatable :: obstype_subkey(:)
type(visir_metadata_type), pointer     :: visir_obs_metadata(:)
type(visir_metadata_type)              :: missing_visir_metadata
type(mw_metadata_type),    pointer     :: mw_obs_metadata(:)
type(mw_metadata_type)                 :: missing_mw_metadata

character(len=5), parameter :: VISIR_STRING = 'visir_rttov'
character(len=5), parameter :: MW_STRING    = 'mw_rttov'

logical :: debug = .false.
integer :: MAXrttovkey = 100000  !FIXME - some initial number of obs
integer ::    rttovkey = 0       ! useful length of metadata arrays
integer ::    visirnum = 0
integer ::       mwnum = 0

integer(kind=jpim), parameter :: maxchannels = 10000

character(len=512)   :: rttov_sensor_db_file = 'unspecified'

! -----------------------------------------------------------------------------
! DART/RTTOV options in the input.nml namelist.
! 
! DART exposes all of the RTTOV 12.3 options available and passes them to 
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
logical              :: fix_hgpl             = .false.  ! surface elevation assigned to 2m pressure (true) or surface pressure (false)
logical              :: do_lambertian        = .false.  ! treat surface as Lambertian instead of specular? (all)
logical              :: lambertian_fixed_angle = .true. ! use fixed angle for Lambertian calculations? (all, do_lambertian only)
logical              :: rad_down_lin_tau     = .true.   ! use linear-in-tau approximation? (all)
logical              :: use_q2m              = .false.  ! use surface humidity? (all)
logical              :: use_uv10m            = .false.  ! use u and v 10 meters? (all, used in sea surface emissivity and BRDF models)
logical              :: use_wfetch           = .false.  ! use wind fetch (length of water wind has blown over in m)  (all, used in sea surface BRDF models)
logical              :: use_water_type       = .false.  ! use water type (0 = fresh, ocean = 1) (all, used in surface BRDF atlas and models)
logical              :: addrefrac            = .false.  ! enable atmospheric refraction (all) 
logical              :: plane_parallel       = .false.  ! treat atmosphere as strictly plane-parallel? (all)
logical              :: use_salinity         = .false.  ! use ocean salinity (in practical salinity units) (MW, FASTEM 4-6 and TESSEM2)
logical              :: use_fastem_params    = .false.  ! use 5 FASTEM parameters describing the land/sea ice surface, see Table 21 (MW, used in FASTEM)
logical              :: use_specularity      = .false.  ! use specularity (VIS/IR, only if do_lambertian is true)
logical              :: apply_band_correction= .true.   ! apply band correction from coef file? (MW)
logical              :: cfrac_data           = .false.  ! specify cloud fraction? (VIS/IR/MW)
logical              :: clw_data             = .false.  ! specify non-precip cloud liquid water? (VIS/IR/MW)
logical              :: rain_data            = .false.  ! specify precip cloud liquid water? (VIS/IR/MW)
logical              :: ciw_data             = .false.  ! specify non-precip cloud ice? (VIS/IR/MW)
logical              :: snow_data            = .false.  ! specify precip cloud fluffy ice? (VIS/IR/MW)
logical              :: graupel_data         = .false.  ! specify precip cloud soft-hail? (VIS/IR/MW)
logical              :: hail_data            = .false.  ! specify precip cloud hard-hail? (VIS/IR/MW)
logical              :: w_data               = .false.  ! specify vertical velocity (used for classifying clouds as cumulus versus stratus)? (VIS/IR)
integer              :: clw_scheme           = 2        ! Liebe (1) or Rosenkranz (2) or TKC (3) (MW, clear-sky only)
real(8)              :: clw_cloud_top        = 322.d0   ! lower hPa limit for clw calculations; clw at lower pressures is ignored (MW, clear-sky only)
integer              :: fastem_version       = 6        ! MW sea-surface emissivity model to use (0-6). 1-6: FASTEM version 1-6, 0: TESSEM2 (MW)
logical              :: supply_foam_fraction = .false.  ! include foam fraction in skin%foam_fraction? FASTEM only. (MW)
logical              :: use_totalice         = .false.  ! Specify totalice instead of precip/non-precip ice (MW, RTTOV-SCATT only)
logical              :: use_zeeman           = .false.  ! Simulate Zeeman effect (MW)
real(r8)             :: cc_threshold         = 0.05d0   ! if effective cloud fraction below this value, treat simulation as clear-sky (MW, 0-1, RTTOV-SCATT only)
logical              :: ozone_data           = .false.  ! specify ozone profiles? (VIS/IR)
logical              :: co2_data             = .false.  ! specify CO2 profiles? (VIS/IR)
logical              :: n2o_data             = .false.  ! specify N2O profiles? (VIS/IR)
logical              :: co_data              = .false.  ! specify CO profiles? (VIS/IR)
logical              :: ch4_data             = .false.  ! specify CH4 profiles? (VIS/IR)
logical              :: so2_data             = .false.  ! specify SO2 profiles? (VIS/IR)
logical              :: add_solar            = .false.  ! include solar calculations (VIS/IR)
logical              :: rayleigh_single_scatt = .true.  ! if false, disable Rayleigh (VIS, add_solar only)
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
real(r8)             :: cldstr_threshold     = -1.d0    ! threshold for cloud stream weights for scattering (VIS/IR, add_clouds only)
logical              :: cldstr_simple        = .false.  ! If true, one clear column, one cloudy column (VIS/IR, add_clouds only)
real(r8)             :: cldstr_low_cloud_top = 750.d0   ! cloud fraction maximum in layers from ToA down to specified hPa (VIS/IR, cldstr_simple only)
integer              :: ir_scatt_model       = 2        ! DOM (1) or Chou-scaling (2) (IR, add_clouds or add_aerosl only)
integer              :: vis_scatt_model      = 1        ! DOM (1), single scat (2), or MFASIS (3) (VIS, add_solar and add_clouds or add_aerosl only)
integer              :: dom_nstreams         = 8        ! number of streams to use with DOM (VIS/IR, add_clouds or add_aerosl and DOM model only, must be >= 2 and even)
real(r8)             :: dom_accuracy         = 0.d0     ! convergence criteria for DOM (VIS/IR, add_clouds or addaerosol and DOM model only)
real(r8)             :: dom_opdep_threshold  = 0.d0     ! DOM ignores layers below this optical depth (VIS/IR, add_clouds or addaerosol and DOM model only)
logical              :: addpc                = .false.  ! do principal component calculations? (VIS/IR)
integer              :: npcscores            = -1       ! number of PC scores to use (VIS/IR, addpc only)
logical              :: addradrec            = .false.  ! reconstruct the radiances (VIS/IR, addpc only)
integer              :: ipcreg               = 1        ! number of predictors, see Table 29 of user guide (VIS/IR, addpc only)
logical              :: use_htfrtc           = .false.  ! use HTFRTC of Havemann 2018  
integer              :: htfrtc_n_pc          = -1       ! number of PCs to use (HTFRTC only, max 300)
logical              :: htfrtc_simple_cloud  = .false.  ! use simple-cloud scattering (HTFRTC only)
logical              :: htfrtc_overcast      = .false.  ! calculate overcast radiances (HTFRTC only)

namelist / obs_def_rttov_nml/ rttov_sensor_db_file,   &
                              mw_clear_sky_only,      &
                              interp_mode,            &
                              do_checkinput,          &
                              apply_reg_limits,       &
                              verbose,                &
                              fix_hgpl,               &
                              do_lambertian,          &
                              lambertian_fixed_angle, &
                              rad_down_lin_tau,       &
                              use_q2m,                &
                              use_uv10m,              &
                              use_wfetch,             &
                              use_water_type,         &
                              addrefrac,              &
                              plane_parallel,         &
                              use_salinity,           &
                              use_fastem_params,      &
                              use_specularity,        &
                              apply_band_correction,  &
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
                              add_solar,              &
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
                              cldstr_threshold,       &
                              cldstr_simple,          &
                              cldstr_low_cloud_top,   &
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
                              htfrtc_overcast

type(atmos_profile_type)     :: atmos
type(trace_gas_profile_type) :: trace_gas
type(aerosol_profile_type)   :: aerosols
type(cloud_profile_type)     :: clouds

integer               :: year
integer               :: month
integer               :: day

contains

!----------------------------------------------------------------------------

subroutine initialize_module

integer :: istatus
integer :: ios

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
missing_mw_metadata%fastem_land1 = MISSING_R8
missing_mw_metadata%fastem_land2 = MISSING_R8
missing_mw_metadata%fastem_land3 = MISSING_R8
missing_mw_metadata%fastem_land4 = MISSING_R8
missing_mw_metadata%fastem_land5 = MISSING_R8

allocate(obstype_metadata(2*MAXrttovkey))
allocate(obstype_subkey(2*MAXrttovkey))
allocate(visir_obs_metadata(MAXrttovkey))
allocate(mw_obs_metadata(MAXrttovkey))

obstype_subkey(:)    = -1
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
   write(string1,*)'mw_clear_sky_only      -',mw_clear_sky_only
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'interp_mode            -',interp_mode
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'do_checkinput          -',do_checkinput
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'apply_reg_limits       -',apply_reg_limits
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'verbose                -',verbose
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'fix_hgpl               -',fix_hgpl
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'do_lambertian          -',do_lambertian
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'lambertian_fixed_angle -',lambertian_fixed_angle
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'rad_down_lin_tau       -',rad_down_lin_tau
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_q2m                -',use_q2m
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_uv10m              -',use_uv10m
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_wfetch             -',use_wfetch
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_water_type         -',use_water_type
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'addrefrac              -',addrefrac
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'plane_parallel         -',plane_parallel
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_salinity           -',use_salinity
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_fastem_params      -',use_fastem_params
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_specularity        -',use_specularity
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'apply_band_correction  -',apply_band_correction
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cfrac_data             -',cfrac_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'clw_data               -',clw_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'rain_data              -',rain_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ciw_data               -',ciw_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'snow_data              -',snow_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'graupel_data           -',graupel_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'hail_data              -',hail_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'w_data                 -',w_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'clw_scheme             -',clw_scheme
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'clw_cloud_top          -',clw_cloud_top
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'fastem_version         -',fastem_version
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'supply_foam_fraction   -',supply_foam_fraction
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_totalice           -',use_totalice
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_zeeman             -',use_zeeman
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cc_threshold           -',cc_threshold
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ozone_data             -',ozone_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'co2_data               -',co2_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'n2o_data               -',n2o_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'co_data                -',co_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ch4_data               -',ch4_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'so2_data               -',so2_data
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'add_solar              -',add_solar
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'rayleigh_single_scatt  -',rayleigh_single_scatt
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'do_nlte_correction     -',do_nlte_correction
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'solar_sea_brdf_model   -',solar_sea_brdf_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ir_sea_emis_model      -',ir_sea_emis_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_sfc_snow_frac      -',use_sfc_snow_frac
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'add_aerosl             -',add_aerosl
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'aerosl_type            -',aerosl_type
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'add_clouds             -',add_clouds
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ice_scheme             -',ice_scheme
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_icede              -',use_icede
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'idg_scheme             -',idg_scheme
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'user_aer_opt_param     -',user_aer_opt_param
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'user_cld_opt_param     -',user_cld_opt_param
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'grid_box_avg_cloud     -',grid_box_avg_cloud
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cldstr_threshold       -',cldstr_threshold
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cldstr_simple          -',cldstr_simple
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'cldstr_low_cloud_top   -',cldstr_low_cloud_top
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ir_scatt_model         -',ir_scatt_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'vis_scatt_model        -',vis_scatt_model
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'dom_nstreams           -',dom_nstreams
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'dom_accuracy           -',dom_accuracy
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'dom_opdep_threshold    -',dom_opdep_threshold
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'addpc                  -',addpc
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'npcscores              -',npcscores
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'addradrec              -',addradrec
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ipcreg                 -',ipcreg
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'use_htfrtc             -',use_htfrtc
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'htfrtc_n_pc            -',htfrtc_n_pc
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'htfrtc_simple_cloud           -',htfrtc_simple_cloud
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'htfrtc_overcast               -',htfrtc_overcast
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if

call read_sensor_db_file(rttov_sensor_db_file)

end subroutine initialize_module

subroutine initialize_rttov_sensor_runtime(sensor,ens_size,nlevels)

type(rttov_sensor_type), pointer    :: sensor
integer,                 intent(in) :: ens_size
integer,                 intent(in) :: nlevels

integer :: istatus
type(rttov_options),       pointer  :: opts
type(rttov_options_scatt), pointer  :: opts_scatt
logical :: is_mw
logical :: is_vis
logical :: is_ir
logical :: is_visir

if (.not. module_initialized) then
   write(string1,*)"The module must be initialized before initializing a rttov sensor runtime"
   call error_handler(E_ERR, 'initialize_rttov', string1, source, revision, revdate)
end if

if (.not. associated(sensor)) then
   write(string1,*)"The sensor must be initialized before initializing a rttov sensor runtime"
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

! mw options could be used with direct or scatt, so set them here
opts % rt_mw % apply_band_correction = apply_band_correction
opts % rt_mw % clw_data              = clw_data
opts % rt_mw % clw_calc_on_coef_lev  = .false.                 ! calculate CLW optical depths on input levels
opts % rt_mw % clw_scheme            = clw_scheme
opts % rt_mw % clw_cloud_top         = clw_cloud_top
opts % rt_mw % fastem_version        = fastem_version
opts % rt_mw % supply_foam_fraction  = supply_foam_fraction

if (is_mw .and. .not. mw_clear_sky_only) then
   ! use RTTOV-SCATT 
   allocate(opts_scatt)

   opts_scatt % config % do_checkinput    = do_checkinput
   opts_scatt % config % apply_reg_limits = apply_reg_limits
   opts_scatt % config % verbose          = verbose
   opts_scatt % config % fix_hgpl         = fix_hgpl

   opts_scatt % rad_down_lin_tau      = rad_down_lin_tau
   opts_scatt % apply_band_correction = apply_band_correction
   opts_scatt % interp_mode           = interp_mode
   opts_scatt % lgradp                = .false.               ! Do not allow TL/AD/K of user pressure levels
   opts_scatt % reg_limit_extrap      = .true.                ! intelligently extend beyond the model top
   opts_scatt % fastem_version        = fastem_version
   opts_scatt % supply_foam_fraction  = supply_foam_fraction
   opts_scatt % use_q2m               = use_q2m
   opts_scatt % lradiance             = .false.               ! Calculate Brightness Temperatures 
   opts_scatt % lusercfrac            = .false.               ! Have RTTOV-SCATT calculate the effective cloud fraction 
   opts_scatt % cc_threshold          = cc_threshold

   call sensor_runtime_setup(sensor,                &
                             ens_size=ens_size,     &
                             nlevs=nlevels,         &
                             opts=opts,             &
                             opts_scatt=opts_scatt)
else
   opts % rt_ir % ozone_data            = ozone_data
   opts % rt_ir % co2_data              = co2_data
   opts % rt_ir % n2o_data              = n2o_data
   opts % rt_ir % co_data               = co_data
   opts % rt_ir % ch4_data              = ch4_data
   opts % rt_ir % so2_data              = so2_data
   opts % rt_ir % addsolar              = add_solar
   opts % rt_ir % rayleigh_single_scatt = rayleigh_single_scatt
   opts % rt_ir % do_nlte_correction    = do_nlte_correction
   opts % rt_ir % solar_sea_brdf_model  = solar_sea_brdf_model
   opts % rt_ir % ir_sea_emis_model     = ir_sea_emis_model
   opts % rt_ir % addaerosl             = add_aerosl
   opts % rt_ir % addclouds             = add_clouds
   opts % rt_ir % user_aer_opt_param    = user_aer_opt_param
   opts % rt_ir % user_cld_opt_param    = user_cld_opt_param
   opts % rt_ir % grid_box_avg_cloud    = grid_box_avg_cloud
   opts % rt_ir % cldstr_threshold      = cldstr_threshold
   opts % rt_ir % cldstr_simple         = cldstr_simple
   opts % rt_ir % cldstr_low_cloud_top  = cldstr_low_cloud_top
   opts % rt_ir % ir_scatt_model        = ir_scatt_model
   opts % rt_ir % vis_scatt_model       = vis_scatt_model
   opts % rt_ir % dom_nstreams          = dom_nstreams
   opts % rt_ir % dom_accuracy          = dom_accuracy
   opts % rt_ir % dom_opdep_threshold   = dom_opdep_threshold

   opts % rt_ir % pc % addpc     = addpc
   opts % rt_ir % pc % npcscores = npcscores
   opts % rt_ir % pc % addradrec = addradrec
   opts % rt_ir % pc % ipcbnd    = 1           ! This is the correct default for RTTOV 12.3 
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

subroutine set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
   platform_id, sat_id, sensor_id, channel, specularity)

integer,  intent(out) :: key
real(r8), intent(in)  :: sat_az, sat_ze
real(r8), intent(in)  :: sun_az, sun_ze ! only relevant if add_solar
integer,  intent(in)  :: platform_id, sat_id, sensor_id, channel
real(r8), intent(in)  :: specularity    ! only relevant if do_lambertian

if ( .not. module_initialized ) call initialize_module

visirnum = visirnum + 1  ! increase module storage used counters
rttovkey = rttovkey + 1
obstype_metadata(rttovkey) = .true. ! .true. for visir
  obstype_subkey(rttovkey) = visirnum

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(rttovkey,'set_visir_metadata', .true.)

key = rttovkey ! now that we know its legal

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

subroutine set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, &
   channel, mag_field, cosbk, fastem_land1, fastem_land2, fastem_land3,         &
   fastem_land4, fastem_land5)

integer,  intent(out) :: key
real(r8), intent(in)  :: sat_az, sat_ze
integer,  intent(in)  :: platform_id, sat_id, sensor_id, channel
real(r8), intent(in)  :: mag_field, cosbk                            ! only relevant with use_zeeman
real(r8), intent(in)  :: fastem_land1, fastem_land2, fastem_land3, & ! 
                         fastem_land4, fastem_land5                  ! only relevant with use_fastem_params

if ( .not. module_initialized ) call initialize_module

mwnum    = mwnum + 1    ! increase module storage used counters
rttovkey = rttovkey + 1 
obstype_metadata(rttovkey) = .false. ! .false. for MW 
  obstype_subkey(rttovkey) = mwnum

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(rttovkey,'set_mw_metadata', .false.)

key = rttovkey ! now that we know its legal

mw_obs_metadata(mwnum) % sat_az       = sat_az
mw_obs_metadata(mwnum) % sat_ze       = sat_ze
mw_obs_metadata(mwnum) % platform_id  = platform_id
mw_obs_metadata(mwnum) % sat_id       = sat_id
mw_obs_metadata(mwnum) % sensor_id    = sensor_id
mw_obs_metadata(mwnum) % channel      = channel
mw_obs_metadata(mwnum) % mag_field    = mag_field
mw_obs_metadata(mwnum) % cosbk        = cosbk
mw_obs_metadata(mwnum) % fastem_land1 = fastem_land1
mw_obs_metadata(mwnum) % fastem_land2 = fastem_land2
mw_obs_metadata(mwnum) % fastem_land3 = fastem_land3
mw_obs_metadata(mwnum) % fastem_land4 = fastem_land4
mw_obs_metadata(mwnum) % fastem_land5 = fastem_land5

end subroutine set_mw_metadata


!----------------------------------------------------------------------
! Query the metadata in module storage for a particular observation.

subroutine get_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
   platform_id, sat_id, sensor_id, channel, specularity)

integer,  intent(in)  :: key
real(r8), intent(out) :: sat_az, sat_ze, sun_az, sun_ze
integer,  intent(out) :: platform_id, sat_id, sensor_id, channel
real(r8), intent(out) :: specularity

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key,'get_visir_metadata')

sat_az      = visir_obs_metadata(key) % sat_az
sat_ze      = visir_obs_metadata(key) % sat_ze
sun_az      = visir_obs_metadata(key) % sun_az
sun_ze      = visir_obs_metadata(key) % sun_ze
platform_id = visir_obs_metadata(key) % platform_id
sat_id      = visir_obs_metadata(key) % sat_id
sensor_id   = visir_obs_metadata(key) % sensor_id
channel     = visir_obs_metadata(key) % channel
specularity = visir_obs_metadata(key) % specularity

end subroutine get_visir_metadata


subroutine get_mw_metadata(key, sat_az, sat_ze,        &
   platform_id, sat_id, sensor_id, channel, mag_field, &
   cosbk, fastem_land1, fastem_land2, fastem_land3,    &
   fastem_land4, fastem_land5)

integer,  intent(in)  :: key
real(r8), intent(out) :: sat_az, sat_ze
integer,  intent(out) :: platform_id, sat_id, sensor_id, channel
real(r8), intent(out) :: mag_field, cosbk
real(r8), intent(out) :: fastem_land1, fastem_land2, fastem_land3, &
                         fastem_land4, fastem_land5

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key,'get_mw_metadata')

sat_az       = mw_obs_metadata(key) % sat_az
sat_ze       = mw_obs_metadata(key) % sat_ze
platform_id  = mw_obs_metadata(key) % platform_id
sat_id       = mw_obs_metadata(key) % sat_id
sensor_id    = mw_obs_metadata(key) % sensor_id
channel      = mw_obs_metadata(key) % channel
mag_field    = mw_obs_metadata(key) % mag_field
cosbk        = mw_obs_metadata(key) % cosbk
fastem_land1 = mw_obs_metadata(key) % fastem_land1
fastem_land2 = mw_obs_metadata(key) % fastem_land2
fastem_land3 = mw_obs_metadata(key) % fastem_land3
fastem_land4 = mw_obs_metadata(key) % fastem_land4
fastem_land5 = mw_obs_metadata(key) % fastem_land5

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
real(r8)          :: fastem_land1, fastem_land2, fastem_land3, &
                     fastem_land4, fastem_land5

logical :: is_visir

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

write(string2,*)'observation #',obsID

if ( is_asciifile ) then
   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_rttov_metadata','header',string2)
   if (trim(header) == trim(visir_STRING)) then
      ! Load the visible/IR data from the file
      read(ifile, *, iostat=ierr) sat_az, sat_ze, sun_az, sun_ze
      call check_iostat(ierr,'read_rttov_metadata','sat,sun az/ze',string2)
      read(ifile, *, iostat=ierr) platform_id, sat_id, sensor_id, channel
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, *, iostat=ierr) specularity
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, *, iostat=ierr) oldkey
      call check_iostat(ierr,'read_rttov_metadata','oldkey',string2)
      is_visir = .true.
   elseif (trim(header) == trim(mw_STRING)) then
      ! Load the visible/IR data from the file
      read(ifile, *, iostat=ierr) sat_az, sat_ze
      call check_iostat(ierr,'read_rttov_metadata','sat az/ze',string2)
      read(ifile, *, iostat=ierr) platform_id, sat_id, sensor_id, channel
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, *, iostat=ierr) mag_field, cosbk
      call check_iostat(ierr,'read_rttov_metadata','mag_field/cosbk',string2)
      read(ifile, *, iostat=ierr) fastem_land1, fastem_land2, fastem_land3, &
                                  fastem_land4, fastem_land5
      call check_iostat(ierr,'read_rttov_metadata','fastem_land1-5',string2)
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

   if (trim(header) == trim(visir_STRING)) then
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
   elseif (trim(header) == trim(mw_STRING)) then
      ! Load the visible/IR data from the file
      read(ifile, iostat=ierr) sat_az, sat_ze
      call check_iostat(ierr,'read_rttov_metadata','sat az/ze',string2)
      read(ifile, iostat=ierr) platform_id, sat_id, sensor_id, channel
      call check_iostat(ierr,'read_rttov_metadata','platform/sat_id/sensor/channel',string2)
      read(ifile, iostat=ierr) mag_field, cosbk
      call check_iostat(ierr,'read_rttov_metadata','mag_field/cosbk',string2)
      read(ifile, iostat=ierr) fastem_land1, fastem_land2, fastem_land3, &
                                  fastem_land4, fastem_land5
      call check_iostat(ierr,'read_rttov_metadata','fastem_land1-5',string2)
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
      channel, mag_field, cosbk, fastem_land1, fastem_land2, fastem_land3,   &
      fastem_land4, fastem_land5)
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
real(r8) :: fastem_land1, fastem_land2, fastem_land3, &
            fastem_land4, fastem_land5

logical  :: is_visir

is_visir = obstype_metadata(key)

if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation.

is_asciifile = ascii_file_format(fform)

if (is_visir) then
   call get_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
      platform_id, sat_id, sensor_id, channel, specularity)

   if (is_asciifile) then
      write(ifile, *) trim(visir_STRING)
      write(ifile, *) sat_az, sat_ze, sun_az, sun_ze
      write(ifile, *) platform_id, sat_id, sensor_id, channel
      write(ifile, *) specularity
      write(ifile, *) key
   else
      write(ifile   ) trim(visir_STRING)
      write(ifile   ) sat_az, sat_ze, sun_az, sun_ze
      write(ifile   ) platform_id, sat_id, sensor_id, channel
      write(ifile   ) specularity
      write(ifile   ) key
   endif
else
   call get_mw_metadata(key, sat_az, sat_ze,                  &
      platform_id, sat_id, sensor_id, channel, mag_field, cosbk, &
      fastem_land1, fastem_land2, fastem_land3, fastem_land4,    &
      fastem_land5)

   if (is_asciifile) then
      write(ifile, *) trim(mw_STRING)
      write(ifile, *) sat_az, sat_ze
      write(ifile, *) platform_id, sat_id, sensor_id, channel
      write(ifile, *) mag_field, cosbk
      write(ifile, *) fastem_land1, fastem_land2, fastem_land3, &
                      fastem_land4, fastem_land5 
      write(ifile, *) key
   else
      write(ifile   ) trim(mw_STRING)
      write(ifile   ) sat_az, sat_ze
      write(ifile   ) platform_id, sat_id, sensor_id, channel
      write(ifile   ) mag_field, cosbk
      write(ifile   ) fastem_land1, fastem_land2, fastem_land3, &
                      fastem_land4, fastem_land5 
      write(ifile   ) key
   end if
end if


end subroutine write_rttov_metadata


!----------------------------------------------------------------------
subroutine interactive_rttov_metadata(key)
integer, intent(out) :: key

real(r8)          :: sat_az, sat_ze, sun_az, sun_ze
integer           :: platform_id, sat_id, sensor_id, channel
real(r8)          :: specularity
real(r8)          :: mag_field, cosbk
real(r8)          :: fastem_land1, fastem_land2, fastem_land3, &
                     fastem_land4, fastem_land5

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
if (is_visir .and. add_solar) then
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

if (.not. is_visir .and. use_fastem_params) then
   fastem_land1 = interactive_r('fastem land1  1st fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_land2 = interactive_r('fastem land2  2nd fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_land3 = interactive_r('fastem land3  3rd fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_land4 = interactive_r('fastem land4  4th fastem land/sea ice parameter', minvalue = 0.0_r8)
   fastem_land5 = interactive_r('fastem land5  5th fastem land/sea ice parameter', minvalue = 0.0_r8)
else
   fastem_land1 = MISSING_R8
   fastem_land2 = MISSING_R8
   fastem_land3 = MISSING_R8
   fastem_land4 = MISSING_R8
   fastem_land5 = MISSING_R8
end if

if (is_visir) then
   call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
      platform_id, sat_id, sensor_id, channel, specularity)
else
   call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, &
      channel, mag_field, cosbk, fastem_land1, fastem_land2, fastem_land3,   &
      fastem_land4, fastem_land5)
end if

end subroutine interactive_rttov_metadata


!----------------------------------------------------------------------
! prompt for a real value, optionally setting min and/or max limits
! loops until valid value input.

function interactive_r(str1,minvalue,maxvalue)
real(r8)                       :: interactive_r
character(len=*),   intent(in) :: str1
real(r8), optional, intent(in) :: minvalue
real(r8), optional, intent(in) :: maxvalue


! Prompt and ensure value is in range if limits are specified

if (present(minvalue) .and. present(maxvalue)) then

   interactive_r = minvalue - 1.0_r8
   MINMAXLOOP : do while ((interactive_r < minvalue) .or. (interactive_r > maxvalue))
      write(*, *) 'Enter '//str1
      read( *, *) interactive_r
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive_r = minvalue - 1.0_r8
   MINLOOP : do while (interactive_r < minvalue)
      write(*, *) 'Enter '//str1
      read( *, *) interactive_r
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive_r = maxvalue + 1.0_r8
   MAXLOOP : do while (interactive_r > maxvalue) 
      write(*, *) 'Enter '//str1
      read( *, *) interactive_r
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive_r
endif

end function interactive_r


!----------------------------------------------------------------------
! prompt for an integer value, optionally setting min and/or max limits
! loops until valid value input.

function interactive_i(str1,minvalue,maxvalue)
integer                        :: interactive_i
character(len=*),   intent(in) :: str1
integer,  optional, intent(in) :: minvalue
integer,  optional, intent(in) :: maxvalue

! Prompt with a minimum amount of error checking

if (present(minvalue) .and. present(maxvalue)) then

   interactive_i = minvalue - 1
   MINMAXLOOP : do while ((interactive_i < minvalue) .and. (interactive_i > maxvalue))
      write(*, *) 'Enter '//str1
      read( *, *) interactive_i
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive_i = minvalue - 1
   MINLOOP : do while (interactive_i < minvalue)
      write(*, *) 'Enter '//str1
      read( *, *) interactive_i
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive_i = maxvalue + 1
   MAXLOOP : do while (interactive_i > maxvalue)
      write(*, *) 'Enter '//str1
      read( *, *) interactive_i
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive_i
endif

end function interactive_i


!----------------------------------------------------------------------

subroutine get_expected_radiance(obs_kind_ind, state_handle, ens_size, location, key, val, istatus)

integer,             intent(in)  :: obs_kind_ind
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location          ! location of obs
integer,             intent(in)  :: key               ! key into module metadata
real(r8),            intent(out) :: val(ens_size)     ! value of obs
integer,             intent(out) :: istatus(ens_size) ! status of the calculation

real(r8) :: sat_az, sat_ze, sun_az, sun_ze
integer  :: platform_id, sat_id, sensor_id, channel

integer :: this_istatus(ens_size)

integer  :: i, zi
real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat
real(r8) :: loc_value(ens_size), radiance(ens_size)
type(location_type) :: loc
integer :: imem, maxlevels, numlevels
integer :: error_status
logical :: return_now

type(rttov_sensor_type), pointer :: sensor

character(len=*), parameter :: routine = 'get_expected_radiance'

type(visir_metadata_type), pointer :: visir_md
type(mw_metadata_type),    pointer :: mw_md

logical :: is_visir

!=================================================================================

if ( .not. module_initialized ) call initialize_module

val = 0.0_r8 ! set return value early

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key, routine)

is_visir = obstype_metadata(key)

if (is_visir) then
   visir_md =  visir_obs_metadata(obstype_subkey(key))
   mw_md    => null()
else
   visir_md => null()
   mw_md    =  mw_obs_metadata(obstype_subkey(key))
end if

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

   if (debug) then
      print *,'istatus is:',this_istatus
   end if

   if ((numlevels == maxlevels) .or. (numlevels == 0)) then
      write(string1,'(A,I0)') 'FAILED to determine number of levels in model:', &
         numlevels
         
      if (debug) call error_handler(E_MSG,routine,string1,source,revision,revdate)
      istatus = 1
      val     = MISSING_R8
      return
   else
       if (debug) write(*,*) routine // ' we have ',numlevels,' model levels'
   endif

   ! allocate the necessary fields
   allocate(atmos%temperature(ens_size, numlevels), &
            atmos%   moisture(ens_size, numlevels), &
            atmos%   pressure(ens_size, numlevels), &
            atmos%      sfc_p(ens_size),            &
            atmos%      s2m_t(ens_size),            &
            atmos%  skin_temp(ens_size),            &
            atmos%   sfc_elev(ens_size),            &
            atmos%   surftype(ens_size))

   if (use_q2m) then
      allocate(atmos%s2m_q(ens_size))
   end if
  
   if (use_uv10m) then
      allocate(atmos%s10m_u(ens_size))
      allocate(atmos%s10m_v(ens_size))
   end if

   if (use_wfetch) then
      allocate(atmos%wfetch(ens_size))
   end if

   if (use_water_type) then
      allocate(atmos%water_type(ens_size))
   end if

   if (use_salinity) then
      allocate(atmos%sfc_salinity(ens_size))
   end if

   if (supply_foam_fraction) then
      allocate(atmos%sfc_foam_frac(ens_size))
   end if

   if (use_sfc_snow_frac) then
      allocate(atmos%sfc_snow_frac(ens_size))
   end if

   if (ozone_data) then
      allocate(trace_gas%ozone(ens_size, numlevels))
   end if

   if (co2_data) then
      allocate(trace_gas%co2(ens_size, numlevels))
   end if

   if (n2o_data) then
      allocate(trace_gas%n2o(ens_size, numlevels))
   end if

   if (ch4_data) then
      allocate(trace_gas%ch4(ens_size, numlevels))
   end if

   if (co_data) then
      allocate(trace_gas%co(ens_size, numlevels))
   end if

   if (add_aerosl) then
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
      else
         ! error
         write(string1,*)"Unknown aerosl_type:",aerosl_type
         call error_handler(E_ERR, routine, string1, source, revision, revdate)
      end if
   end if

   ! RTTOV wants layers, but models probably prefer levels
   if (cfrac_data) then
      allocate(clouds%cfrac(ens_size, numlevels))
   end if

   if (clw_data) then
      allocate(clouds%clw(ens_size, numlevels))

      if (clw_scheme == 2) then
         allocate(clouds%clwde(ens_size, numlevels)) 
      end if
   end if

   if (rain_data) then
      allocate(clouds%rain(ens_size, numlevels))
   end if

   if (ciw_data) then
      allocate(clouds%ciw(ens_size, numlevels))
      if (ice_scheme == 1 .and. use_icede) then
         allocate(clouds%icede(ens_size, numlevels))
      end if
   end if

   if (snow_data) then
      allocate(clouds%snow(ens_size, numlevels))
   end if

   if (graupel_data) then
      allocate(clouds%graupel(ens_size, numlevels))
   end if

   if (hail_data) then
      allocate(clouds%hail(ens_size, numlevels))
   end if

   if (w_data) then
      allocate(clouds%w(ens_size, numlevels))
   end if

   if (htfrtc_simple_cloud) then
      allocate(clouds%simple_cfrac(ens_size))
      allocate(clouds%ctp(ens_size))
   end if
else
   ! load the number of levels from the temperature array, assumed to have been initialized
   numlevels = size(atmos%temperature,2)
end if 

! also check if the sensor runtime needs to be initialized
sensor => get_rttov_sensor(platform_id, sat_id, sensor_id)

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
   call check_status('QTY_PRESSURE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

   call interpolate(state_handle, ens_size, loc, QTY_TEMPERATURE, atmos%temperature(:, i), this_istatus)
   call check_status('QTY_TEMPERATURE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

   call interpolate(state_handle, ens_size, loc, QTY_SPECIFIC_HUMIDITY, atmos%moisture(:, i), this_istatus)
   call check_status('QTY_SPECIFIC_HUMIDITY', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

   if (ozone_data) then
      call interpolate(state_handle, ens_size, loc, QTY_O3, trace_gas%ozone(:, i), this_istatus)
      call check_status('QTY_O3', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (co2_data) then
      call interpolate(state_handle, ens_size, loc, QTY_CO2, trace_gas%co2(:, i), this_istatus)
      call check_status('QTY_CO2', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (n2o_data) then
      call interpolate(state_handle, ens_size, loc, QTY_N2O, trace_gas%n2o(:, i), this_istatus)
      call check_status('QTY_N2O', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (ch4_data) then
      call interpolate(state_handle, ens_size, loc, QTY_CH4, trace_gas%ch4(:, i), this_istatus)
      call check_status('QTY_CH4', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (co_data) then
      call interpolate(state_handle, ens_size, loc, QTY_CO, trace_gas%co(:, i), this_istatus)
      call check_status('QTY_CO', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (add_aerosl) then
      if (aerosl_type == 1) then
         ! OPAC
         call interpolate(state_handle, ens_size, loc, QTY_INSOLUBLE_AER, aerosols%insoluble(:, i), this_istatus)
         call check_status('QTY_INSOLUBLE_AER', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_H2O_SOLUBLE_AER, aerosols%water_soluble(:, i), this_istatus)
         call check_status('QTY_H2O_SOLUBLE_AER', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_SOOT, aerosols%soot(:, i), this_istatus)
         call check_status('QTY_SOOT', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_SEASALT_ACCUM, aerosols%sea_salt_accum(:, i), this_istatus)
         call check_status('QTY_SEASALT_ACCUM', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_SEASALT_COARSE, aerosols%sea_salt_coarse(:, i), this_istatus)
         call check_status('QTY_SEASALT_COARSE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_NUCLEUS, aerosols%mineral_nucleus(:, i), this_istatus)
         call check_status('QTY_MINERAL_NUCLEUS', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_ACCUM, aerosols%mineral_accum(:, i), this_istatus)
         call check_status('QTY_MINERAL_ACCUM', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_COARSE, aerosols%mineral_coarse(:, i), this_istatus)
         call check_status('QTY_MINERAL_COARSE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_MINERAL_TRANSPORTED, aerosols%mineral_transport(:, i), this_istatus)
         call check_status('QTY_MINERAL_TRANSPORT', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_SULPHATED_DROPS, aerosols%sulphated_droplets(:, i), this_istatus)
         call check_status('QTY_SULPHATED_DROPS', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_VOLCANIC_ASH, aerosols%volcanic_ash(:, i), this_istatus)
         call check_status('QTY_VOLCANIC_ASH', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_NEW_VOLCANIC_ASH, aerosols%new_volcanic_ash(:, i), this_istatus)
         call check_status('QTY_NEW_VOLCANIC_ASH', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_ASIAN_DUST, aerosols%asian_dust(:, i), this_istatus)
         call check_status('QTY_ASIAN_DUST', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
      elseif (aerosl_type == 2) then
         ! CAMS
         call interpolate(state_handle, ens_size, loc, QTY_BLACK_CARBON, aerosols%black_carbon(:, i), this_istatus)
         call check_status('QTY_BLACK_CARBON', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_DUST_BIN1, aerosols%dust_bin1(:, i), this_istatus)
         call check_status('QTY_DUST_BIN1', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_DUST_BIN2, aerosols%dust_bin2(:, i), this_istatus)
         call check_status('QTY_DUST_BIN2', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_DUST_BIN3, aerosols%dust_bin3(:, i), this_istatus)
         call check_status('QTY_DUST_BIN3', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_AMMONIUM_SULPHATE, aerosols%ammonium_sulphate(:, i), this_istatus)
         call check_status('QTY_AMMONIUM_SULPHATE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_SEA_SALT_BIN1, aerosols%sea_salt_bin1(:, i), this_istatus)
         call check_status('QTY_SEA_SALT_BIN1', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_SEA_SALT_BIN2, aerosols%sea_salt_bin2(:, i), this_istatus)
         call check_status('QTY_SEA_SALT_BIN2', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_SEA_SALT_BIN3, aerosols%sea_salt_bin3(:, i), this_istatus)
         call check_status('QTY_SEA_SALT_BIN3', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

         call interpolate(state_handle, ens_size, loc, QTY_HYDROPHILIC_ORGANIC_MATTER, aerosols%hydrophilic_organic_matter(:, i), this_istatus)
         call check_status('QTY_HYDROPHILIC_ORGANIC_MATTER', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
      end if
   end if

   if (cfrac_data) then
      ! specify cloud fraction
      call interpolate(state_handle, ens_size, loc, QTY_CLOUD_FRACTION, clouds%cfrac(:, i), this_istatus)
      call check_status('QTY_CLOUD_FRACTION', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (clw_data) then
      ! specify non-precip cloud liquid water
      call interpolate(state_handle, ens_size, loc, QTY_CLOUDWATER_MIXING_RATIO, clouds%clw(:, i), this_istatus)
      call check_status('QTY_CLOUDWATER_MIXING_RATIO', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

      if (clw_scheme == 2) then
         ! The effective diameter must also be specified with clw_scheme 2
         call interpolate(state_handle, ens_size, loc, QTY_CLOUDWATER_DE, clouds%clwde(:, i), this_istatus)
         call check_status('QTY_CLOUDWATER_DE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
      end if
   end if

   if (rain_data) then
      ! specify precip cloud liquid water (i.e. rain)
      call interpolate(state_handle, ens_size, loc, QTY_RAINWATER_MIXING_RATIO, clouds%rain(:, i), this_istatus)
      call check_status('QTY_RAINWATER_MIXING_RATIO', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (ciw_data) then
      ! specify non-precip cloud ice
      call interpolate(state_handle, ens_size, loc, QTY_ICE_MIXING_RATIO, clouds%ciw(:, i), this_istatus)
      call check_status('QTY_ICE_MIXING_RATIO', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

      if (ice_scheme == 1 .and. use_icede) then
         ! if use_icede with ice_scheme 1, must also specify ice effective diameter
         call interpolate(state_handle, ens_size, loc, QTY_CLOUD_ICE_DE, clouds%icede(:, i), this_istatus)
         call check_status('QTY_CLOUD_ICE_DE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
      end if
   end if

   if (snow_data) then
      ! specify precip fluffy ice (i.e. snow)
      call interpolate(state_handle, ens_size, loc, QTY_SNOW_MIXING_RATIO, clouds%snow(:, i), this_istatus)
      call check_status('QTY_SNOW_MIXING_RATIO', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (graupel_data) then
      ! specify precip soft-hail (i.e. graupel)
      call interpolate(state_handle, ens_size, loc, QTY_GRAUPEL_MIXING_RATIO, clouds%graupel(:, i), this_istatus)
      call check_status('QTY_GRAUPEL_MIXING_RATIO', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (hail_data) then
      ! specify precip hard-hail (i.e. hail)
      call interpolate(state_handle, ens_size, loc, QTY_HAIL_MIXING_RATIO, clouds%hail(:, i), this_istatus)
      call check_status('QTY_HAIL_MIXING_RATIO', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

   if (w_data) then
      ! specify vertical velocity
      call interpolate(state_handle, ens_size, loc, QTY_VERTICAL_VELOCITY, clouds%w(:, i), this_istatus)
      call check_status('QTY_VERTICAL_VELOCITY', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   end if

end do GETLEVELDATA

loc = set_location(loc_lon, loc_lat, 1.0d0, VERTISLEVEL )

! set the surface fields
call interpolate(state_handle, ens_size, loc, QTY_SURFACE_PRESSURE, atmos%sfc_p(:), this_istatus)
call check_status('QTY_SURFACE_PRESSURE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

call interpolate(state_handle, ens_size, loc, QTY_SURFACE_ELEVATION, atmos%sfc_elev(:), this_istatus)
call check_status('QTY_SURFACE_ELEVATION', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

call interpolate(state_handle, ens_size, loc, QTY_2M_TEMPERATURE, atmos%s2m_t(:), this_istatus)
call check_status('QTY_2M_TEMPERATURE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

call interpolate(state_handle, ens_size, loc, QTY_SKIN_TEMPERATURE, atmos%skin_temp(:), this_istatus)
call check_status('QTY_SKIN_TEMPERATURE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

! set to 2m if an error

call interpolate(state_handle, ens_size, loc, QTY_SURFACE_TYPE, atmos%surftype(:), this_istatus)
call check_status('QTY_SURFACE_TYPE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .true.)

! if not available, lookup by lat/lon?

if (use_q2m) then
   call interpolate(state_handle, ens_size, loc, QTY_2M_SPECIFIC_HUMIDITY, atmos%s2m_q(:), this_istatus)
   call check_status('QTY_2M_SPECIFIC_HUMIDITY', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (use_uv10m) then
   call interpolate(state_handle, ens_size, loc, QTY_10M_U_WIND_COMPONENT, atmos%s10m_u(:), this_istatus)
   call check_status('QTY_10M_U_WIND_COMPONENT', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
   call interpolate(state_handle, ens_size, loc, QTY_10M_V_WIND_COMPONENT, atmos%s10m_v(:), this_istatus)
   call check_status('QTY_10M_V_WIND_COMPONENT', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (use_wfetch) then
   call interpolate(state_handle, ens_size, loc, QTY_WIND_FETCH, atmos%wfetch(:), this_istatus)
   call check_status('QTY_WIND_FETCH', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (use_water_type) then
   call interpolate(state_handle, ens_size, loc, QTY_WATER_TYPE, atmos%water_type(:), this_istatus)
   call check_status('QTY_WATER_TYPE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (use_salinity) then
   call interpolate(state_handle, ens_size, loc, QTY_SALINITY, atmos%sfc_salinity(:), this_istatus)
   call check_status('QTY_SALINITY', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (supply_foam_fraction) then
   call interpolate(state_handle, ens_size, loc, QTY_FOAM_FRAC, atmos%sfc_foam_frac(:), this_istatus)
   call check_status('QTY_FOAM_FRAC', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (use_sfc_snow_frac) then
   call interpolate(state_handle, ens_size, loc, QTY_SNOWCOVER_FRAC, atmos%sfc_snow_frac(:), this_istatus)
   call check_status('QTY_SNOWCOVER_FRAC', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (add_clouds .and. htfrtc_simple_cloud) then
   ! specify simple cloud information - per column
   call interpolate(state_handle, ens_size, loc, QTY_COLUMN_CLOUD_FRAC, clouds%simple_cfrac(:), this_istatus)
   call check_status('QTY_COLUMN_CLOUD_FRAC', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)

   call interpolate(state_handle, ens_size, loc, QTY_CLOUD_TOP_PRESSURE, clouds%ctp(:), this_istatus)
   call check_status('QTY_CLOUD_TOP_PRESSURE', ens_size, this_istatus, val, istatus, routine, source, revision, revdate, .false.)
end if

if (debug) then
   print*, 'interpolate pressure    = ', atmos%pressure(1,1),    '...', atmos%pressure(1,numlevels)
   print*, 'interpolate temperature = ', atmos%temperature(1,1), '...', atmos%temperature(1,numlevels)
   print*, 'interpolate moisture    = ', atmos%moisture(1,1),    '...', atmos%moisture(1,numlevels)
end if

call do_forward_model(ens_size=ens_size,                    &
                      nlevels=numlevels,                    &
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
                      use_fastem_params=use_fastem_params,  &
                      radiances=radiance,                   &
                      error_status=this_istatus,            &
                      visir_md=visir_md,                    &
                      mw_md=mw_md) 

if (debug) then
   print*, 'istatus  = ', istatus 
   print*, 'radiance = ', radiance
end if

!=================================================================================
! ... and finally set the return the radiance forward operator value

where (istatus == 0) val = radiance
where (istatus /= 0) val = missing_r8

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

integer :: maxkey 

maxkey = size(obstype_metadata)

if ((key > 0) .and. (key <= maxkey)) then
   ! we are still within limits
   return
else
   ! Bad news. Tell the user.
   write(string1, *) 'key (',key,') not within known range ( 1,', maxkey,')'
   call error_handler(E_ERR,routine,string1,source,revision,revdate) 
end if

end subroutine key_within_range

!----------------------------------------------------------------------
! If the allocatable metadata arrays are not big enough ... try again

subroutine grow_metadata(key, routine, is_visir)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine
logical,          intent(in) :: is_visir

integer :: orglength, newlength
type(visir_metadata_type), allocatable :: safe_visir_metadata(:)
type(mw_metadata_type),    allocatable :: safe_mw_metadata(:)
logical,                   allocatable :: safe_obstype_metadata(:)
integer,                   allocatable :: safe_obstype_subkey(:)

integer :: current_obstype_length
integer :: new_obstype_length

current_obstype_length = size(obstype_metadata)
new_obstype_length      = current_obstype_length

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key >= 2*current_obstype_length) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

if (is_visir) then
   if ((key > 0) .and. (key <= size(visir_obs_metadata))) then
      ! we are still within limits
      return
   else
      ! we need to grow visir_obs_metadata
      orglength = size(visir_obs_metadata)
      newlength = 2 * orglength

      ! News. Tell the user we are increasing storage.
      write(string1, *) 'key (',key,') exceeds visir_obs_metadata length (',orglength,')'
      write(string2, *) 'Increasing visir_obs_metadata to length ',newlength
      call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

      allocate(safe_visir_metadata(orglength))
      safe_visir_metadata(:) = visir_obs_metadata(:)

      deallocate(visir_obs_metadata)
        allocate(visir_obs_metadata(newlength))

      visir_obs_metadata(1:orglength)           = safe_visir_metadata(:)
      visir_obs_metadata(orglength+1:newlength) = missing_visir_metadata

      deallocate(safe_visir_metadata)
      new_obstype_length = current_obstype_length + (newlength-orglength) 
   end if
else
   ! duplicate the above since we can't use any object-oriented magic
   if ((key > 0) .and. (key <= size(mw_obs_metadata))) then
      ! we are still within limits
      return
   else
      ! we need to grow mw_obs_metadata
      orglength = size(mw_obs_metadata)
      newlength = 2 * orglength

      ! News. Tell the user we are increasing storage.
      write(string1, *) 'key (',key,') exceeds mw_obs_metadata length (',orglength,')'
      write(string2, *) 'Increasing mw_obs_metadata to length ',newlength
      call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

      allocate(safe_mw_metadata(orglength))
      safe_mw_metadata(:) = mw_obs_metadata(:)

      deallocate(mw_obs_metadata)
        allocate(mw_obs_metadata(newlength))

      mw_obs_metadata(1:orglength)           = safe_mw_metadata(:)
      mw_obs_metadata(orglength+1:newlength) = missing_mw_metadata

      deallocate(safe_mw_metadata)
      new_obstype_length = current_obstype_length + (newlength-orglength) 
   end if
end if

if (current_obstype_length /= new_obstype_length) then
   ! expand the obstype metadata as well, copying over as above
   allocate(safe_obstype_metadata(current_obstype_length))
   safe_obstype_metadata(:) = obstype_metadata(:)

   deallocate(obstype_metadata)
     allocate(obstype_metadata(new_obstype_length))

   obstype_metadata(1:current_obstype_length) = safe_obstype_metadata(:)
   ! obstype_metadata(current_obstype_length+1:new_obstype_length) = -1
   ! no default for obstype_metadata as it is a logical that is true for visir

   allocate(safe_obstype_subkey(current_obstype_length))
   safe_obstype_subkey(:) = obstype_subkey(:)

   deallocate(obstype_subkey)
     allocate(obstype_subkey(new_obstype_length))

   obstype_subkey(1:current_obstype_length) = safe_obstype_subkey(:)
   obstype_subkey(current_obstype_length+1:new_obstype_length) = -1
end if

end subroutine grow_metadata

!----------------------------------------------------------------------

subroutine check_status(field_name, ens_size, this_istatus, val, istatus, routine, source, revision, revdate, required)

   character(len=*), intent(in)  :: field_name
   integer,          intent(in)  :: ens_size
   integer,          intent(in)  :: this_istatus(ens_size)
   real(r8),         intent(out) :: val(ens_size)
   integer,          intent(out) :: istatus(ens_size) ! status of the calculation
   character(len=*), intent(in)  :: routine
   character(len=*), intent(in)  :: source
   character(len=*), intent(in)  :: revision
   character(len=*), intent(in)  :: revdate
   logical,          intent(in)  :: required

   logical :: error

   ! override track_status behavior to always error and print an error message if the field
   ! cannot be found

   call track_status(ens_size, this_istatus, val, istatus, error)

   if (error) then
      if (required) then
         call error_handler(E_ERR,routine,'Could not find required field ' // trim(field_name),source,revision,revdate)
      else
         call error_handler(E_ERR,routine,'Could not find requested field ' // trim(field_name),source,revision,revdate)
      end if
   end if

end subroutine check_status

end module obs_def_rttov_mod

! END DART PREPROCESS MODULE CODE
