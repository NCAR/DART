! This code is not protected by the DART copyright agreement.

! adapted from example AIRS readers provided by 
! http://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/AIRS/3.8_ScienceDataSoftwareTools/V6_FORTRAN_C_READERS.tar.gz
!
! Sep 2020, tjh

module amsua_support_mod

! the contents of this file are an amalgam of:
!  amsua_bt_typ.inc
!  amsua_bt_struct.inc
!  amsua_bt_rdr.f
! although in several cases they were modified to use
! fortran 90 derived types, free format continuation lines, 
! and to avoid using the deprecated syntax 'double precision'.

use         types_mod, only : r8, deg2rad, rad2deg, PI, MISSING_R8

use     utilities_mod, only : error_handler, E_MSG, E_ERR, &
                              is_longitude_between, register_module

use  time_manager_mod, only : time_type, get_date, set_date,            &
                              get_time, set_time, set_calendar_type,    &
                              GREGORIAN, print_date, print_time

use  obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                              set_obs_values, obs_sequence_type,              &
                              obs_type, set_copy_meta_data, set_qc_meta_data, &
                              print_obs_seq_summary

use      location_mod, only : location_type, VERTISUNDEF, &
                              set_location, get_location

use      obs_kind_mod, only : get_index_for_type_of_obs

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use obs_def_rttov_mod, only : set_mw_metadata, &
                              get_rttov_option_logical

use hdf5

implicit none
private

public :: amsua_bt_gran_type, AMSUA_BT_CHANNEL, read_amsua_bt_granule , make_obs_sequence

interface read_attribute
   module procedure read_char_attribute
   module procedure read_integer_attribute
   module procedure read_short_integer_attribute
   module procedure read_byte_attribute
   module procedure read_real_attribute
   module procedure read_real8_attribute
end interface

! Using the HDF-EOS5 functions simply required
! putting an 'he5_' in front of the function. As far as I can tell from
! https://hdfeos.org/examples/fort_he5_swath.php, the calling structure is
! the same.
!
! integer :: he5_swopen         ! open a swath file
! integer :: he5_swattach       ! attatch to a swath object
! integer :: he5_swinqswath     ! retieves number and names of swaths in file
! integer :: he5_swrdfld        ! read data from a data field
! integer :: he5_swrdattr       ! read swath attribute
! integer :: he5_swdetach       ! detatching from the swath object
! integer :: he5_swclose        ! closing the file

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'amsua_support_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

character(len=512) :: string1

! ----------------------------------------------------------------------
! START of amsua_bt_typ.inc - (modified)
! ----------------------------------------------------------------------

integer, parameter :: AMSUA_BT_GEOXTRACK    = 30  ! cross-track footprints
integer, parameter :: AMSUA_BT_GEOTRACK     = 45  ! along-track scanlines
integer, parameter :: AMSUA_BT_CHANNEL      = 15  ! number of microwave channels
integer, parameter :: AMSUA_BT_CALXTRACK    =  4 
integer, parameter :: AMSUA_BT_SPACEXTRACK  =  2 
integer, parameter :: AMSUA_BT_BBXTRACK     =  2 
integer, parameter :: AMSUA_BT_WARMPRTA11   =  5 
integer, parameter :: AMSUA_BT_WARMPRTA12   =  5 
integer, parameter :: AMSUA_BT_WARMPRTA2    =  7 

! Type definitions for each record type

      TYPE  amsua_bt_QA_bb_PRT_a11_t

! Attributes

        real           min
        real           max
        real           mean
        real           dev
        integer      num_in
        integer      num_lo
        integer      num_hi
        integer      num_bad
        real           range_min
        real           range_max
        byte           missing
        integer      max_track
        integer      max_xtrack
        integer      min_track
        integer      min_xtrack

      END TYPE


      TYPE  amsua_bt_QA_bb_PRT_a12_t

! Attributes

        real           min
        real           max
        real           mean
        real           dev
        integer      num_in
        integer      num_lo
        integer      num_hi
        integer      num_bad
        real           range_min
        real           range_max
        byte           missing
        integer      max_track
        integer      max_xtrack
        integer      min_track
        integer      min_xtrack
      END TYPE


      TYPE  amsua_bt_QA_bb_PRT_a2_t

! Attributes

        real           min
        real           max
        real           mean
        real           dev
        integer      num_in
        integer      num_lo
        integer      num_hi
        integer      num_bad
        real           range_min
        real           range_max
        byte           missing
        integer      max_track
        integer      max_xtrack
        integer      min_track
        integer      min_xtrack
      END TYPE


      TYPE  amsua_bt_QA_rec_PRT_a11_t

! Attributes

        real           min
        real           max
        real           mean
        real           dev
        integer      num_in
        integer      num_lo
        integer      num_hi
        integer      num_bad
        real           range_min
        real           range_max
        byte           missing
        integer      max_track
        integer      max_xtrack
        integer      min_track
        integer      min_xtrack
      END TYPE


      TYPE  amsua_bt_QA_rec_PRT_a12_t

! Attributes

        real           min
        real           max
        real           mean
        real           dev
        integer      num_in
        integer      num_lo
        integer      num_hi
        integer      num_bad
        real           range_min
        real           range_max
        byte           missing
        integer      max_track
        integer      max_xtrack
        integer      min_track
        integer      min_xtrack
      END TYPE


      TYPE  amsua_bt_QA_rec_PRT_a2_t

! Attributes

        real           min
        real           max
        real           mean
        real           dev
        integer      num_in
        integer      num_lo
        integer      num_hi
        integer      num_bad
        real           range_min
        real           range_max
        byte           missing
        integer      max_track
        integer      max_xtrack
        integer      min_track
        integer      min_xtrack
      END TYPE


      TYPE  amsua_bt_bb_signals_t

! Attributes

        real           min( AMSUA_BT_CHANNEL, 2)
        real           max( AMSUA_BT_CHANNEL, 2)
        real           mean( AMSUA_BT_CHANNEL, 2)
        real           dev( AMSUA_BT_CHANNEL, 2)
        integer        num( AMSUA_BT_CHANNEL, 2)
        integer        num_bad( AMSUA_BT_CHANNEL, 2)
        integer        max_track( AMSUA_BT_CHANNEL, 2)
        integer        max_xtrack( AMSUA_BT_CHANNEL, 2)
        integer        min_track( AMSUA_BT_CHANNEL, 2)
        integer        min_xtrack( AMSUA_BT_CHANNEL, 2)
      END TYPE


      TYPE  amsua_bt_space_signals_t

! Attributes

        real           min( AMSUA_BT_CHANNEL, 2)
        real           max( AMSUA_BT_CHANNEL, 2)
        real           mean( AMSUA_BT_CHANNEL, 2)
        real           dev( AMSUA_BT_CHANNEL, 2)
        integer        num( AMSUA_BT_CHANNEL, 2)
        integer        num_bad( AMSUA_BT_CHANNEL, 2)
        integer        max_track( AMSUA_BT_CHANNEL, 2)
        integer        max_xtrack( AMSUA_BT_CHANNEL, 2)
        integer        min_track( AMSUA_BT_CHANNEL, 2)
        integer        min_xtrack( AMSUA_BT_CHANNEL, 2)
      END TYPE


      TYPE  amsua_bt_gain_stats_t

! Attributes

        real           min( AMSUA_BT_CHANNEL)
        real           max( AMSUA_BT_CHANNEL)
        real           mean( AMSUA_BT_CHANNEL)
        real           dev( AMSUA_BT_CHANNEL)
        integer        num( AMSUA_BT_CHANNEL)
        integer        num_bad( AMSUA_BT_CHANNEL)
        integer        max_track( AMSUA_BT_CHANNEL)
        integer        max_xtrack( AMSUA_BT_CHANNEL)
        integer        min_track( AMSUA_BT_CHANNEL)
        integer        min_xtrack( AMSUA_BT_CHANNEL)
      END TYPE


      TYPE  amsua_bt_QA_unfiltered_scene_count_t

! Attributes

        real           min( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        real           max( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        real           mean( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        real           dev( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        integer        num( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        integer        num_bad( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        integer        max_track( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        integer        max_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        integer        min_track( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
        integer        min_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
      END TYPE


      TYPE  amsua_bt_QA_unfiltered_BB_count_t

! Attributes

        real           min( AMSUA_BT_CHANNEL, 2)
        real           max( AMSUA_BT_CHANNEL, 2)
        real           mean( AMSUA_BT_CHANNEL, 2)
        real           dev( AMSUA_BT_CHANNEL, 2)
        integer        num( AMSUA_BT_CHANNEL, 2)
        integer        num_bad( AMSUA_BT_CHANNEL, 2)
        integer        max_track( AMSUA_BT_CHANNEL, 2)
        integer        max_xtrack( AMSUA_BT_CHANNEL, 2)
        integer        min_track( AMSUA_BT_CHANNEL, 2)
        integer        min_xtrack( AMSUA_BT_CHANNEL, 2)
      END TYPE


      TYPE  amsua_bt_QA_unfiltered_space_count_t

! Attributes

        real           min( AMSUA_BT_CHANNEL, 2)
        real           max( AMSUA_BT_CHANNEL, 2)
        real           mean( AMSUA_BT_CHANNEL, 2)
        real           dev( AMSUA_BT_CHANNEL, 2)
        integer        num( AMSUA_BT_CHANNEL, 2)
        integer        num_bad( AMSUA_BT_CHANNEL, 2)
        integer        max_track( AMSUA_BT_CHANNEL, 2)
        integer        max_xtrack( AMSUA_BT_CHANNEL, 2)
        integer        min_track( AMSUA_BT_CHANNEL, 2)
        integer        min_xtrack( AMSUA_BT_CHANNEL, 2)
      END TYPE


      TYPE  amsua_bt_QA_cal_coef_a0_t

! Attributes

        real           min( AMSUA_BT_CHANNEL)
        real           max( AMSUA_BT_CHANNEL)
        real           mean( AMSUA_BT_CHANNEL)
        real           dev( AMSUA_BT_CHANNEL)
        integer        num( AMSUA_BT_CHANNEL)
        integer        num_bad( AMSUA_BT_CHANNEL)
        integer        max_track( AMSUA_BT_CHANNEL)
        integer        max_xtrack( AMSUA_BT_CHANNEL)
        integer        min_track( AMSUA_BT_CHANNEL)
        integer        min_xtrack( AMSUA_BT_CHANNEL)
      END TYPE


      TYPE  amsua_bt_QA_cal_coef_a1_t

! Attributes

        real           min( AMSUA_BT_CHANNEL)
        real           max( AMSUA_BT_CHANNEL)
        real           mean( AMSUA_BT_CHANNEL)
        real           dev( AMSUA_BT_CHANNEL)
        integer        num( AMSUA_BT_CHANNEL)
        integer        num_bad( AMSUA_BT_CHANNEL)
        integer        max_track( AMSUA_BT_CHANNEL)
        integer        max_xtrack( AMSUA_BT_CHANNEL)
        integer        min_track( AMSUA_BT_CHANNEL)
        integer        min_xtrack( AMSUA_BT_CHANNEL)
      END TYPE


      TYPE  amsua_bt_QA_cal_coef_a2_t

! Attributes

        real           min( AMSUA_BT_CHANNEL)
        real           max( AMSUA_BT_CHANNEL)
        real           mean( AMSUA_BT_CHANNEL)
        real           dev( AMSUA_BT_CHANNEL)
        integer        num( AMSUA_BT_CHANNEL)
        integer        num_bad( AMSUA_BT_CHANNEL)
        integer        max_track( AMSUA_BT_CHANNEL)
        integer        max_xtrack( AMSUA_BT_CHANNEL)
        integer        min_track( AMSUA_BT_CHANNEL)
        integer        min_xtrack( AMSUA_BT_CHANNEL)
      END TYPE


      TYPE  amsua_bt_QA_bb_raw_noise_counts_t

! Attributes

        real           min( AMSUA_BT_CHANNEL)
        real           max( AMSUA_BT_CHANNEL)
        real           mean( AMSUA_BT_CHANNEL)
        real           dev( AMSUA_BT_CHANNEL)
        integer        num( AMSUA_BT_CHANNEL)
        integer        num_bad( AMSUA_BT_CHANNEL)
        integer        max_track( AMSUA_BT_CHANNEL)
        integer        max_xtrack( AMSUA_BT_CHANNEL)
        integer        min_track( AMSUA_BT_CHANNEL)
        integer        min_xtrack( AMSUA_BT_CHANNEL)
      END TYPE


      TYPE  amsua_bt_QA_sv_raw_noise_counts_t

! Attributes

        real           min( AMSUA_BT_CHANNEL)
        real           max( AMSUA_BT_CHANNEL)
        real           mean( AMSUA_BT_CHANNEL)
        real           dev( AMSUA_BT_CHANNEL)
        integer        num( AMSUA_BT_CHANNEL)
        integer        num_bad( AMSUA_BT_CHANNEL)
        integer        max_track( AMSUA_BT_CHANNEL)
        integer        max_xtrack( AMSUA_BT_CHANNEL)
        integer        min_track( AMSUA_BT_CHANNEL)
        integer        min_xtrack( AMSUA_BT_CHANNEL)
      END TYPE

! ----------------------------------------------------------------------
! END of amsua_bt_typ.inc - (modified)
! ----------------------------------------------------------------------
! Content from 'Section 1.5 Data Disclaimer'
!
! *Invalid Values
! Fields in Level 1B and Level 2 data products may contain an invalid value:
! -9999 for floating-point and 16-bit and 32-bit integers
! -1 or 255 for 8-bit fields.
!
! Content from 'Section 2.4 Key data fields'
!
! * qa_scanline
! Bit field for each scanline (bit 0 set if sun glint in scanline; bit 1 set if costal
! crossing in scanline, bit 2 set if some channels had excessive NeDT estimated)
!
! * qa_channel
! Bit field by channel for each scanline (bit 0 set if all space view counts bad; bit 1 set
! if space view counts marginal; bit 2 set if space view counts could not be smoothed;
! bit 3 set if all blackbody counts bad; bit 4 set if blackbody counts marginal; bit 5 set
! if blackbody counts could not be smoothed; bit 6 set if unable to calculate
! calibration coefficients; bit 7 set if excessive NeDT estimated)
! 
! * antenna_temp
! calibrated, geolocated channel-by-channel AMSU observed raw antenna temperature (K)
! 
! * brightness_temp
! calibrated, geolocated channel-by-channel AMSU sidelobe-corrected antenna temperature (K)
! 
! * brightness_temp_err
! error estimate for brightness_temp (K)
! 
! * landFrac
! fraction of AMSU footprint that is land (0.0 -> 1.0)
! 
! * landFrac_err
! error estimate for landFrac
! 
! ----------------------------------------------------------------------
! START of amsua_bt_struct.inc - (greatly modified)
! ----------------------------------------------------------------------

! Record holds an entire granule of amsua_bt
TYPE  amsua_bt_gran_type

   ! Attributes
   character*256  processing_level
   character*256  instrument
   character*256  DayNightFlag
   character*256  AutomaticQAFlag
   integer      NumTotalData
   integer      NumProcessData
   integer      NumSpecialData
   integer      NumBadData
   integer      NumMissingData
   integer      NumLandSurface
   integer      NumOceanSurface
   character*256  node_type
   integer      start_year
   integer      start_month
   integer      start_day
   integer      start_hour
   integer      start_minute
   real         start_sec
   integer      start_orbit
   integer      end_orbit
   integer      orbit_path
   integer      start_orbit_row
   integer      end_orbit_row
   integer      granule_number
   integer      num_scansets
   integer      num_scanlines
   real*8       start_Latitude
   real*8       start_Longitude
   real*8       start_Time
   real*8       end_Latitude
   real*8       end_Longitude
   real*8       end_Time
   real         eq_x_longitude
   real*8       eq_x_tai
   integer      orbitgeoqa
   integer*2    num_satgeoqa
   integer*2    num_glintgeoqa
   integer*2    num_moongeoqa
   integer*2    num_ftptgeoqa
   integer*2    num_zengeoqa
   integer*2    num_demgeoqa
   integer*2    num_fpe
   integer*2    LonGranuleCen
   integer*2    LatGranuleCen
   integer*2    LocTimeGranuleCen
   integer      num_scanlines_not_norm_mode_a1
   integer      num_scanlines_not_norm_mode_a2
   integer      num_missing_scanlines_a1
   integer      num_missing_scanlines_a2
   integer      num_data_gaps_a1
   integer      num_data_gaps_a2
   integer      num_instr_mode_changes_a1
   integer      num_instr_mode_changes_a2
   integer      num_scanlines_rec_cal_prob_a11
   integer      num_scanlines_rec_cal_prob_a12
   integer      num_scanlines_rec_cal_prob_a2
   integer      num_scanlines_sig_coast_xing
   integer      num_scanlines_sig_sun_glint
   integer      MoonInViewMWCount

   TYPE(amsua_bt_QA_bb_PRT_a11_t)  QA_bb_PRT_a11
   TYPE(amsua_bt_QA_bb_PRT_a12_t)  QA_bb_PRT_a12
   TYPE(amsua_bt_QA_bb_PRT_a2_t)   QA_bb_PRT_a2
   TYPE(amsua_bt_QA_rec_PRT_a11_t) QA_rec_PRT_a11
   TYPE(amsua_bt_QA_rec_PRT_a12_t) QA_rec_PRT_a12
   TYPE(amsua_bt_QA_rec_PRT_a2_t)  QA_rec_PRT_a2

   character*256  granules_present

! Geolocation fields
! Latitude ... geodetic latitude, degrees N [-90,90]
! Longitude .. geodetic longitude, degrees E [-180,180]
! Time     ... 'shutter' TAI Time, floating point elapsed seconds since 1 Jan 1993

   real*8  Latitude(AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real*8 Longitude(AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real*8      Time(AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)

! Data Fields

   real       scanang(    AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       satheight(  AMSUA_BT_GEOTRACK)
   real       satroll(    AMSUA_BT_GEOTRACK)
   real       satpitch(   AMSUA_BT_GEOTRACK)
   real       satyaw(     AMSUA_BT_GEOTRACK)
   integer    satgeoqa(   AMSUA_BT_GEOTRACK)
   integer*2  glintgeoqa( AMSUA_BT_GEOTRACK)
   integer*2  moongeoqa(  AMSUA_BT_GEOTRACK)
   integer    ftptgeoqa(  AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   integer*2  zengeoqa(   AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   integer*2  demgeoqa(   AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real*8     nadirTAI(                                AMSUA_BT_GEOTRACK)
   real*8     sat_lat(                                 AMSUA_BT_GEOTRACK)
   real*8     sat_lon(                                 AMSUA_BT_GEOTRACK)
   byte       scan_node_type(                          AMSUA_BT_GEOTRACK)
   real       satzen(              AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       satazi(              AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       solzen(              AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       solazi(              AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       glintlat(                                AMSUA_BT_GEOTRACK)
   real       glintlon(                                AMSUA_BT_GEOTRACK)
   integer*2  sun_glint_distance(  AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       topog(               AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       topog_err(           AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       landFrac(            AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       landFrac_err(        AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       antenna_temp(        AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       brightness_temp(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
   real       brightness_temp_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)

   real       center_freq(                   AMSUA_BT_CHANNEL) ! center channel frequency (GHz)
   real       IF_offset_1(                   AMSUA_BT_CHANNEL) ! offset of first intermediate frequency (MHz)
   real       IF_offset_2(                   AMSUA_BT_CHANNEL) ! offset of second    ..       frequency (MHz)
   real       bandwidth(                     AMSUA_BT_CHANNEL) ! bandwith of sum of 1,2, or 4 channels (MHz)
   integer    num_calibrated_scanlines(      AMSUA_BT_CHANNEL)
   integer    num_scanlines_ch_cal_problems( AMSUA_BT_CHANNEL)
   real       NeDT(                          AMSUA_BT_CHANNEL) ! instrument noise level estimated from warm count scatter

   TYPE(amsua_bt_bb_signals_t)    bb_signals
   TYPE(amsua_bt_space_signals_t) space_signals
   TYPE(amsua_bt_gain_stats_t)    gain_stats
   TYPE(amsua_bt_QA_unfiltered_scene_count_t) QA_unfiltered_scene_count
   TYPE(amsua_bt_QA_unfiltered_BB_count_t)    QA_unfiltered_BB_count
   TYPE(amsua_bt_QA_unfiltered_space_count_t) QA_unfiltered_space_count
   TYPE(amsua_bt_QA_cal_coef_a0_t)            QA_cal_coef_a0
   TYPE(amsua_bt_QA_cal_coef_a1_t)            QA_cal_coef_a1
   TYPE(amsua_bt_QA_cal_coef_a2_t)            QA_cal_coef_a2
   TYPE(amsua_bt_QA_bb_raw_noise_counts_t)    QA_bb_raw_noise_counts
   TYPE(amsua_bt_QA_sv_raw_noise_counts_t)    QA_sv_raw_noise_counts
   integer    state1( AMSUA_BT_GEOTRACK)
   integer    state2( AMSUA_BT_GEOTRACK)
   real       cal_coef_a0(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
   real       cal_coef_a0_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
   real       cal_coef_a1(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
   real       cal_coef_a1_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
   real       cal_coef_a2(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
   real       cal_coef_a2_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
   byte       a1_ColdCalPstion(      AMSUA_BT_GEOTRACK)
   byte       a2_ColdCalPstion(      AMSUA_BT_GEOTRACK)
   byte       a1_PLO_Redundncy(      AMSUA_BT_GEOTRACK)
   byte       a11_mux_temp_used(     AMSUA_BT_GEOTRACK)
   real       a11_receiver_temp(     AMSUA_BT_GEOTRACK)
   real       a11_target_temp(       AMSUA_BT_GEOTRACK)
   byte       a12_mux_temp_used(     AMSUA_BT_GEOTRACK)
   real       a12_receiver_temp(     AMSUA_BT_GEOTRACK)
   real       a12_target_temp(       AMSUA_BT_GEOTRACK)
   byte       a2_diplexer_temp_used( AMSUA_BT_GEOTRACK)
   real       a2_receiver_temp(      AMSUA_BT_GEOTRACK)
   real       a2_target_temp(        AMSUA_BT_GEOTRACK)
   byte       qa_scanline(           AMSUA_BT_GEOTRACK)
   byte       qa_receiver_a11(       AMSUA_BT_GEOTRACK)
   byte       qa_receiver_a12(       AMSUA_BT_GEOTRACK)
   byte       qa_receiver_a2(        AMSUA_BT_GEOTRACK)
   byte       qa_channel(            AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
END TYPE

! ----------------------------------------------------------------------
! END of amsua_bt_struct.inc
! ----------------------------------------------------------------------

contains

! ----------------------------------------------------------------------
! START of amsua_bt_rdr.f - (heavily modified)
! ----------------------------------------------------------------------

! This function was autogenerated by the mkezio program to read an
! AIRS swath (AKA 'granule') of type "L1B_AMSU" from file given by the
! file_name argument into a buffer pointed to by the amsua_bt_gran
! argument.  The caller owns the buffer.  The entire granule
! is read -- every attribute and field, the whole lat/lon/time
! extent.
!
! Errors opening the file, etc. are fatal and cause STOP.
! Problems reading individual attributes or fields are reported to
! the console but do not interrupt program flow.
!
! Modified by TJH Sep 2020: 
!      use F90 syntax
!      call a function to read the attributes
!      renamed the subroutine from amsua_bt_rdr 
!      note: 'swath' seems historical, 'granule' is in the current documentation.

subroutine read_amsua_bt_granule(file_name, amsua_bt_gran)

IMPLICIT NONE

character(len=*),         intent(in)  :: file_name
TYPE(amsua_bt_gran_type), intent(out) :: amsua_bt_gran

integer :: statn                ! HDF-EOS status. 0 for success
integer :: fid                  ! HDF-EOS file ID
integer :: swid                 ! HDF-EOS swath ID
integer :: nchar                ! Number of characters
character(len=256) :: swathname ! Name of swath
integer :: nswath               ! Number of swaths
integer :: start(10) = 0        ! start of each dimensions for Swath I/O
                                ! 0 => start with first element
integer :: stride(10) = 1       ! stride of each dimensions for Swath I/O
                                ! 1 => use every element
integer :: edge(10)             ! size of each dimension for swath I/O
                                ! will be set for each individual read

integer :: he5_swopen           ! open a swath file
integer :: he5_swinqswath       ! retrieves number and names of swaths in file
integer :: he5_swattach         ! attatch to a swath object
integer :: he5_swrdfld          ! read data from a data field
integer :: he5_swdetach         ! detatching from the swath object
integer :: he5_swclose          ! closes the file

if ( .not. module_initialized ) call initialize_module

fid = he5_swopen(file_name, 101)
if (fid .eq. -1) then
  print *, "Error ", fid, " opening file ", file_name
  stop
end if

! Get name of swath(s)
nswath = he5_swinqswath(file_name, swathname, nchar)
if (nswath .ne. 1) then
  print *, "he5_swinqswath found ", nswath, " swaths for file ", file_name, " Need exactly 1"
  stop
end if

! There's exactly one swath.  Make sure it is the right one.
if (swathname .ne. 'L1B_AMSU') then
  print *, "Error: bad swath name ", swathname, " in file ", file_name
  print *, "Expected L1B_AMSU"
  stop
end if

! Attach to (open) the one swath.
swid = he5_swattach(fid, swathname)
if (swid .eq. -1) then
  print *, "Failed to attach to swath ", swathname, " in file ", file_name
  stop
end if

! Attributes
call read_attribute(swid, "processing_level", amsua_bt_gran%processing_level)
call read_attribute(swid, "instrument", amsua_bt_gran%instrument)
call read_attribute(swid, "DayNightFlag", amsua_bt_gran%DayNightFlag)
call read_attribute(swid, "AutomaticQAFlag", amsua_bt_gran%AutomaticQAFlag)
call read_attribute(swid, "NumTotalData", amsua_bt_gran%NumTotalData)
call read_attribute(swid, "NumProcessData", amsua_bt_gran%NumProcessData)
call read_attribute(swid, "NumSpecialData", amsua_bt_gran%NumSpecialData)
call read_attribute(swid, "NumBadData", amsua_bt_gran%NumBadData)
call read_attribute(swid, "NumMissingData",amsua_bt_gran%NumMissingData)
call read_attribute(swid, "NumLandSurface", amsua_bt_gran%NumLandSurface)
call read_attribute(swid, "NumOceanSurface", amsua_bt_gran%NumOceanSurface)
call read_attribute(swid, "node_type", amsua_bt_gran%node_type)
call read_attribute(swid, "start_year", amsua_bt_gran%start_year)
call read_attribute(swid, "start_month", amsua_bt_gran%start_month)
call read_attribute(swid, "start_day", amsua_bt_gran%start_day)
call read_attribute(swid, "start_hour", amsua_bt_gran%start_hour)
call read_attribute(swid, "start_minute", amsua_bt_gran%start_minute)
call read_attribute(swid, "start_sec", amsua_bt_gran%start_sec)
call read_attribute(swid, "start_orbit", amsua_bt_gran%start_orbit)
call read_attribute(swid, "end_orbit", amsua_bt_gran%end_orbit)
call read_attribute(swid, "orbit_path", amsua_bt_gran%orbit_path)
call read_attribute(swid, "start_orbit_row", amsua_bt_gran%start_orbit_row)
call read_attribute(swid, "end_orbit_row", amsua_bt_gran%end_orbit_row)
call read_attribute(swid, "granule_number", amsua_bt_gran%granule_number)
call read_attribute(swid, "num_scansets", amsua_bt_gran%num_scansets)
call read_attribute(swid, "num_scanlines", amsua_bt_gran%num_scanlines)
call read_attribute(swid, "start_Latitude", amsua_bt_gran%start_Latitude)
call read_attribute(swid, "start_Longitude", amsua_bt_gran%start_Longitude)
call read_attribute(swid, "start_Time", amsua_bt_gran%start_Time)
call read_attribute(swid, "end_Latitude", amsua_bt_gran%end_Latitude)
call read_attribute(swid, "end_Longitude", amsua_bt_gran%end_Longitude)
call read_attribute(swid, "end_Time", amsua_bt_gran%end_Time)
call read_attribute(swid, "eq_x_longitude", amsua_bt_gran%eq_x_longitude)
call read_attribute(swid, "eq_x_tai", amsua_bt_gran%eq_x_tai)
call read_attribute(swid, "orbitgeoqa", amsua_bt_gran%orbitgeoqa)
call read_attribute(swid, "num_satgeoqa", amsua_bt_gran%num_satgeoqa)
call read_attribute(swid, "num_glintgeoqa", amsua_bt_gran%num_glintgeoqa)
call read_attribute(swid, "num_moongeoqa", amsua_bt_gran%num_moongeoqa)
call read_attribute(swid, "num_ftptgeoqa", amsua_bt_gran%num_ftptgeoqa)
call read_attribute(swid, "num_zengeoqa", amsua_bt_gran%num_zengeoqa)
call read_attribute(swid, "num_demgeoqa", amsua_bt_gran%num_demgeoqa)
call read_attribute(swid, "num_fpe", amsua_bt_gran%num_fpe)
call read_attribute(swid, "LonGranuleCen", amsua_bt_gran%LonGranuleCen)
call read_attribute(swid, "LatGranuleCen", amsua_bt_gran%LatGranuleCen)
call read_attribute(swid, "LocTimeGranuleCen", amsua_bt_gran%LocTimeGranuleCen)
call read_attribute(swid, "num_scanlines_not_norm_mode_a1", amsua_bt_gran%num_scanlines_not_norm_mode_a1)
call read_attribute(swid, "num_scanlines_not_norm_mode_a2", amsua_bt_gran%num_scanlines_not_norm_mode_a2)
call read_attribute(swid, "num_missing_scanlines_a1", amsua_bt_gran%num_missing_scanlines_a1)
call read_attribute(swid, "num_missing_scanlines_a2", amsua_bt_gran%num_missing_scanlines_a2)
call read_attribute(swid, "num_data_gaps_a1", amsua_bt_gran%num_data_gaps_a1)
call read_attribute(swid, "num_data_gaps_a2", amsua_bt_gran%num_data_gaps_a2)
call read_attribute(swid, "num_instr_mode_changes_a1", amsua_bt_gran%num_instr_mode_changes_a1)
call read_attribute(swid, "num_instr_mode_changes_a2", amsua_bt_gran%num_instr_mode_changes_a2)
call read_attribute(swid, "num_scanlines_rec_cal_prob_a11", amsua_bt_gran%num_scanlines_rec_cal_prob_a11)
call read_attribute(swid, "num_scanlines_rec_cal_prob_a12", amsua_bt_gran%num_scanlines_rec_cal_prob_a12)
call read_attribute(swid, "num_scanlines_rec_cal_prob_a2", amsua_bt_gran%num_scanlines_rec_cal_prob_a2)
call read_attribute(swid, "num_scanlines_sig_coast_xing", amsua_bt_gran%num_scanlines_sig_coast_xing)
call read_attribute(swid, "num_scanlines_sig_sun_glint", amsua_bt_gran%num_scanlines_sig_sun_glint)
call read_attribute(swid, "MoonInViewMWCount", amsua_bt_gran%MoonInViewMWCount)
call read_attribute(swid, "QA_bb_PRT_a11.min", amsua_bt_gran%QA_bb_PRT_a11%min)
call read_attribute(swid, "QA_bb_PRT_a11.max", amsua_bt_gran%QA_bb_PRT_a11%max)
call read_attribute(swid, "QA_bb_PRT_a11.mean", amsua_bt_gran%QA_bb_PRT_a11%mean)
call read_attribute(swid, "QA_bb_PRT_a11.dev", amsua_bt_gran%QA_bb_PRT_a11%dev)
call read_attribute(swid, "QA_bb_PRT_a11.num_in", amsua_bt_gran%QA_bb_PRT_a11%num_in)
call read_attribute(swid, "QA_bb_PRT_a11.num_lo", amsua_bt_gran%QA_bb_PRT_a11%num_lo)
call read_attribute(swid, "QA_bb_PRT_a11.num_hi", amsua_bt_gran%QA_bb_PRT_a11%num_hi)
call read_attribute(swid, "QA_bb_PRT_a11.num_bad", amsua_bt_gran%QA_bb_PRT_a11%num_bad)
call read_attribute(swid, "QA_bb_PRT_a11.range_min", amsua_bt_gran%QA_bb_PRT_a11%range_min)
call read_attribute(swid, "QA_bb_PRT_a11.range_max", amsua_bt_gran%QA_bb_PRT_a11%range_max)
call read_attribute(swid, "QA_bb_PRT_a11.missing", amsua_bt_gran%QA_bb_PRT_a11%missing)
call read_attribute(swid, "QA_bb_PRT_a11.max_track", amsua_bt_gran%QA_bb_PRT_a11%max_track)
call read_attribute(swid, "QA_bb_PRT_a11.max_xtrack", amsua_bt_gran%QA_bb_PRT_a11%max_xtrack)
call read_attribute(swid, "QA_bb_PRT_a11.min_track", amsua_bt_gran%QA_bb_PRT_a11%min_track)
call read_attribute(swid, "QA_bb_PRT_a11.min_xtrack", amsua_bt_gran%QA_bb_PRT_a11%min_xtrack)
call read_attribute(swid, "QA_bb_PRT_a12.min", amsua_bt_gran%QA_bb_PRT_a12%min)
call read_attribute(swid, "QA_bb_PRT_a12.max", amsua_bt_gran%QA_bb_PRT_a12%max)
call read_attribute(swid, "QA_bb_PRT_a12.mean", amsua_bt_gran%QA_bb_PRT_a12%mean)
call read_attribute(swid, "QA_bb_PRT_a12.dev", amsua_bt_gran%QA_bb_PRT_a12%dev)
call read_attribute(swid, "QA_bb_PRT_a12.num_in", amsua_bt_gran%QA_bb_PRT_a12%num_in)
call read_attribute(swid, "QA_bb_PRT_a12.num_lo", amsua_bt_gran%QA_bb_PRT_a12%num_lo)
call read_attribute(swid, "QA_bb_PRT_a12.num_hi", amsua_bt_gran%QA_bb_PRT_a12%num_hi)
call read_attribute(swid, "QA_bb_PRT_a12.num_bad", amsua_bt_gran%QA_bb_PRT_a12%num_bad)
call read_attribute(swid, "QA_bb_PRT_a12.range_min", amsua_bt_gran%QA_bb_PRT_a12%range_min)
call read_attribute(swid, "QA_bb_PRT_a12.range_max", amsua_bt_gran%QA_bb_PRT_a12%range_max)
call read_attribute(swid, "QA_bb_PRT_a12.missing", amsua_bt_gran%QA_bb_PRT_a12%missing)
call read_attribute(swid, "QA_bb_PRT_a12.max_track", amsua_bt_gran%QA_bb_PRT_a12%max_track)
call read_attribute(swid, "QA_bb_PRT_a12.max_xtrack", amsua_bt_gran%QA_bb_PRT_a12%max_xtrack)
call read_attribute(swid, "QA_bb_PRT_a12.min_track", amsua_bt_gran%QA_bb_PRT_a12%min_track)
call read_attribute(swid, "QA_bb_PRT_a12.min_xtrack", amsua_bt_gran%QA_bb_PRT_a12%min_xtrack)
call read_attribute(swid, "QA_bb_PRT_a2.min", amsua_bt_gran%QA_bb_PRT_a2%min)
call read_attribute(swid, "QA_bb_PRT_a2.max", amsua_bt_gran%QA_bb_PRT_a2%max)
call read_attribute(swid, "QA_bb_PRT_a2.mean", amsua_bt_gran%QA_bb_PRT_a2%mean)
call read_attribute(swid, "QA_bb_PRT_a2.dev", amsua_bt_gran%QA_bb_PRT_a2%dev)
call read_attribute(swid, "QA_bb_PRT_a2.num_in", amsua_bt_gran%QA_bb_PRT_a2%num_in)
call read_attribute(swid, "QA_bb_PRT_a2.num_lo", amsua_bt_gran%QA_bb_PRT_a2%num_lo)
call read_attribute(swid, "QA_bb_PRT_a2.num_hi", amsua_bt_gran%QA_bb_PRT_a2%num_hi)
call read_attribute(swid, "QA_bb_PRT_a2.num_bad", amsua_bt_gran%QA_bb_PRT_a2%num_bad)
call read_attribute(swid, "QA_bb_PRT_a2.range_min", amsua_bt_gran%QA_bb_PRT_a2%range_min)
call read_attribute(swid, "QA_bb_PRT_a2.range_max", amsua_bt_gran%QA_bb_PRT_a2%range_max)
call read_attribute(swid, "QA_bb_PRT_a2.missing", amsua_bt_gran%QA_bb_PRT_a2%missing)
call read_attribute(swid, "QA_bb_PRT_a2.max_track", amsua_bt_gran%QA_bb_PRT_a2%max_track)
call read_attribute(swid, "QA_bb_PRT_a2.max_xtrack", amsua_bt_gran%QA_bb_PRT_a2%max_xtrack)
call read_attribute(swid, "QA_bb_PRT_a2.min_track", amsua_bt_gran%QA_bb_PRT_a2%min_track)
call read_attribute(swid, "QA_bb_PRT_a2.min_xtrack", amsua_bt_gran%QA_bb_PRT_a2%min_xtrack)
call read_attribute(swid, "QA_rec_PRT_a11.min", amsua_bt_gran%QA_rec_PRT_a11%min)
call read_attribute(swid, "QA_rec_PRT_a11.max", amsua_bt_gran%QA_rec_PRT_a11%max)
call read_attribute(swid, "QA_rec_PRT_a11.mean", amsua_bt_gran%QA_rec_PRT_a11%mean)
call read_attribute(swid, "QA_rec_PRT_a11.dev", amsua_bt_gran%QA_rec_PRT_a11%dev)
call read_attribute(swid, "QA_rec_PRT_a11.num_in", amsua_bt_gran%QA_rec_PRT_a11%num_in)
call read_attribute(swid, "QA_rec_PRT_a11.num_lo", amsua_bt_gran%QA_rec_PRT_a11%num_lo)
call read_attribute(swid, "QA_rec_PRT_a11.num_hi", amsua_bt_gran%QA_rec_PRT_a11%num_hi)
call read_attribute(swid, "QA_rec_PRT_a11.num_bad", amsua_bt_gran%QA_rec_PRT_a11%num_bad)
call read_attribute(swid, "QA_rec_PRT_a11.range_min", amsua_bt_gran%QA_rec_PRT_a11%range_min)
call read_attribute(swid, "QA_rec_PRT_a11.range_max", amsua_bt_gran%QA_rec_PRT_a11%range_max)
call read_attribute(swid, "QA_rec_PRT_a11.missing", amsua_bt_gran%QA_rec_PRT_a11%missing)
call read_attribute(swid, "QA_rec_PRT_a11.max_track", amsua_bt_gran%QA_rec_PRT_a11%max_track)
call read_attribute(swid, "QA_rec_PRT_a11.max_xtrack", amsua_bt_gran%QA_rec_PRT_a11%max_xtrack)
call read_attribute(swid, "QA_rec_PRT_a11.min_track", amsua_bt_gran%QA_rec_PRT_a11%min_track)
call read_attribute(swid, "QA_rec_PRT_a11.min_xtrack", amsua_bt_gran%QA_rec_PRT_a11%min_xtrack)
call read_attribute(swid, "QA_rec_PRT_a12.min", amsua_bt_gran%QA_rec_PRT_a12%min)
call read_attribute(swid, "QA_rec_PRT_a12.max", amsua_bt_gran%QA_rec_PRT_a12%max)
call read_attribute(swid, "QA_rec_PRT_a12.mean", amsua_bt_gran%QA_rec_PRT_a12%mean)
call read_attribute(swid, "QA_rec_PRT_a12.dev", amsua_bt_gran%QA_rec_PRT_a12%dev)
call read_attribute(swid, "QA_rec_PRT_a12.num_in", amsua_bt_gran%QA_rec_PRT_a12%num_in)
call read_attribute(swid, "QA_rec_PRT_a12.num_lo", amsua_bt_gran%QA_rec_PRT_a12%num_lo)
call read_attribute(swid, "QA_rec_PRT_a12.num_hi", amsua_bt_gran%QA_rec_PRT_a12%num_hi)
call read_attribute(swid, "QA_rec_PRT_a12.num_bad", amsua_bt_gran%QA_rec_PRT_a12%num_bad)
call read_attribute(swid, "QA_rec_PRT_a12.range_min", amsua_bt_gran%QA_rec_PRT_a12%range_min)
call read_attribute(swid, "QA_rec_PRT_a12.range_max", amsua_bt_gran%QA_rec_PRT_a12%range_max)
call read_attribute(swid, "QA_rec_PRT_a12.missing", amsua_bt_gran%QA_rec_PRT_a12%missing)
call read_attribute(swid, "QA_rec_PRT_a12.max_track", amsua_bt_gran%QA_rec_PRT_a12%max_track)
call read_attribute(swid, "QA_rec_PRT_a12.max_xtrack", amsua_bt_gran%QA_rec_PRT_a12%max_xtrack)
call read_attribute(swid, "QA_rec_PRT_a12.min_track", amsua_bt_gran%QA_rec_PRT_a12%min_track)
call read_attribute(swid, "QA_rec_PRT_a12.min_xtrack", amsua_bt_gran%QA_rec_PRT_a12%min_xtrack)
call read_attribute(swid, "QA_rec_PRT_a2.min", amsua_bt_gran%QA_rec_PRT_a2%min)
call read_attribute(swid, "QA_rec_PRT_a2.max", amsua_bt_gran%QA_rec_PRT_a2%max)
call read_attribute(swid, "QA_rec_PRT_a2.mean", amsua_bt_gran%QA_rec_PRT_a2%mean)
call read_attribute(swid, "QA_rec_PRT_a2.dev", amsua_bt_gran%QA_rec_PRT_a2%dev)
call read_attribute(swid, "QA_rec_PRT_a2.num_in", amsua_bt_gran%QA_rec_PRT_a2%num_in)
call read_attribute(swid, "QA_rec_PRT_a2.num_lo", amsua_bt_gran%QA_rec_PRT_a2%num_lo)
call read_attribute(swid, "QA_rec_PRT_a2.num_hi", amsua_bt_gran%QA_rec_PRT_a2%num_hi)
call read_attribute(swid, "QA_rec_PRT_a2.num_bad", amsua_bt_gran%QA_rec_PRT_a2%num_bad)
call read_attribute(swid, "QA_rec_PRT_a2.range_min", amsua_bt_gran%QA_rec_PRT_a2%range_min)
call read_attribute(swid, "QA_rec_PRT_a2.range_max", amsua_bt_gran%QA_rec_PRT_a2%range_max)
call read_attribute(swid, "QA_rec_PRT_a2.missing", amsua_bt_gran%QA_rec_PRT_a2%missing)
call read_attribute(swid, "QA_rec_PRT_a2.max_track", amsua_bt_gran%QA_rec_PRT_a2%max_track)
call read_attribute(swid, "QA_rec_PRT_a2.max_xtrack", amsua_bt_gran%QA_rec_PRT_a2%max_xtrack)
call read_attribute(swid, "QA_rec_PRT_a2.min_track", amsua_bt_gran%QA_rec_PRT_a2%min_track)
call read_attribute(swid, "QA_rec_PRT_a2.min_xtrack", amsua_bt_gran%QA_rec_PRT_a2%min_xtrack)
call read_attribute(swid, "granules_present", amsua_bt_gran%granules_present)

call dump_attributes(amsua_bt_gran)

! Geolocation fields
edge(1) = AMSUA_BT_GEOXTRACK
edge(2) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "Latitude", start, stride, edge, amsua_bt_gran%Latitude)

statn = he5_swrdfld(swid, "Longitude", start, stride, edge, amsua_bt_gran%Longitude)

statn = he5_swrdfld(swid, "Time", start, stride, edge, amsua_bt_gran%Time)


! Data Fields
edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "scanang", start, stride, edge, amsua_bt_gran%scanang)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "satheight", start, stride, edge, amsua_bt_gran%satheight)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "satroll", start, stride, edge, amsua_bt_gran%satroll)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "satpitch", start, stride, edge, amsua_bt_gran%satpitch)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "satyaw", start, stride, edge, amsua_bt_gran%satyaw)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "satgeoqa", start, stride, edge, amsua_bt_gran%satgeoqa)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "glintgeoqa", start, stride, edge, amsua_bt_gran%glintgeoqa)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "moongeoqa", start, stride, edge, amsua_bt_gran%moongeoqa)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "ftptgeoqa", start, stride, edge, amsua_bt_gran%ftptgeoqa)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "zengeoqa", start, stride, edge, amsua_bt_gran%zengeoqa)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "demgeoqa", start, stride, edge, amsua_bt_gran%demgeoqa)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "nadirTAI", start, stride, edge, amsua_bt_gran%nadirTAI)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "sat_lat", start, stride, edge, amsua_bt_gran%sat_lat)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "sat_lon", start, stride, edge, amsua_bt_gran%sat_lon)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "scan_node_type", start, stride, edge, amsua_bt_gran%scan_node_type)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "satzen", start, stride, edge, amsua_bt_gran%satzen)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "satazi", start, stride, edge, amsua_bt_gran%satazi)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "solzen", start, stride, edge, amsua_bt_gran%solzen)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "solazi", start, stride, edge, amsua_bt_gran%solazi)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "glintlat", start, stride, edge, amsua_bt_gran%glintlat)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "glintlon", start, stride, edge, amsua_bt_gran%glintlon)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "sun_glint_distance", start, stride, edge, amsua_bt_gran%sun_glint_distance)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "topog", start, stride, edge, amsua_bt_gran%topog)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "topog_err", start, stride, edge, amsua_bt_gran%topog_err)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "landFrac", start, stride, edge, amsua_bt_gran%landFrac)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = he5_swrdfld(swid, "landFrac_err", start, stride, edge, amsua_bt_gran%landFrac_err)

edge(3) = AMSUA_BT_GEOTRACK
edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "antenna_temp", start, stride, edge, amsua_bt_gran%antenna_temp)

edge(3) = AMSUA_BT_GEOTRACK
edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "brightness_temp", start, stride, edge, amsua_bt_gran%brightness_temp)

edge(3) = AMSUA_BT_GEOTRACK
edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "brightness_temp_err", start, stride, edge, amsua_bt_gran%brightness_temp_err)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "center_freq", start, stride, edge, amsua_bt_gran%center_freq)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "IF_offset_1", start, stride, edge, amsua_bt_gran%IF_offset_1)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "IF_offset_2", start, stride, edge, amsua_bt_gran%IF_offset_2)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bandwidth", start, stride, edge, amsua_bt_gran%bandwidth)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "num_calibrated_scanlines", start, stride, edge, amsua_bt_gran%num_calibrated_scanlines)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "num_scanlines_ch_cal_problems", start, stride, edge, amsua_bt_gran%num_scanlines_ch_cal_problems)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.min", start, stride, edge, amsua_bt_gran%bb_signals%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.max", start, stride, edge, amsua_bt_gran%bb_signals%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.mean", start, stride, edge, amsua_bt_gran%bb_signals%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.dev", start, stride, edge, amsua_bt_gran%bb_signals%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.num", start, stride, edge, amsua_bt_gran%bb_signals%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.num_bad", start, stride, edge, amsua_bt_gran%bb_signals%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.max_track", start, stride, edge, amsua_bt_gran%bb_signals%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.max_xtrack", start, stride, edge, amsua_bt_gran%bb_signals%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.min_track", start, stride, edge, amsua_bt_gran%bb_signals%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "bb_signals.min_xtrack", start, stride, edge, amsua_bt_gran%bb_signals%min_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.min", start, stride, edge, amsua_bt_gran%space_signals%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.max", start, stride, edge, amsua_bt_gran%space_signals%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.mean", start, stride, edge, amsua_bt_gran%space_signals%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.dev", start, stride, edge, amsua_bt_gran%space_signals%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.num", start, stride, edge, amsua_bt_gran%space_signals%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.num_bad", start, stride, edge, amsua_bt_gran%space_signals%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.max_track", start, stride, edge, amsua_bt_gran%space_signals%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.max_xtrack", start, stride, edge, amsua_bt_gran%space_signals%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.min_track", start, stride, edge, amsua_bt_gran%space_signals%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "space_signals.min_xtrack", start, stride, edge, amsua_bt_gran%space_signals%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.min", start, stride, edge, amsua_bt_gran%gain_stats%min)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.max", start, stride, edge, amsua_bt_gran%gain_stats%max)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.mean", start, stride, edge, amsua_bt_gran%gain_stats%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.dev", start, stride, edge, amsua_bt_gran%gain_stats%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.num", start, stride, edge, amsua_bt_gran%gain_stats%num)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.num_bad", start, stride, edge, amsua_bt_gran%gain_stats%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.max_track", start, stride, edge, amsua_bt_gran%gain_stats%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.max_xtrack", start, stride, edge, amsua_bt_gran%gain_stats%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.min_track", start, stride, edge, amsua_bt_gran%gain_stats%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "gain_stats.min_xtrack", start, stride, edge, amsua_bt_gran%gain_stats%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "NeDT", start, stride, edge, amsua_bt_gran%NeDT)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.min", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%min)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.max", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%max)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.mean", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%mean)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.dev", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%dev)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.num", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%num)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.num_bad", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%num_bad)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.max_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%max_track)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.max_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%max_xtrack)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.min_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%min_track)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_scene_count.min_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%min_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.min", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.max", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.mean", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.dev", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.num", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.num_bad", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.max_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.max_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.min_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_BB_count.min_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%min_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.min", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.max", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.mean", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.dev", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.num", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.num_bad", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.max_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.max_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.min_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_unfiltered_space_count.min_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.min", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%min)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.max", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%max)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.mean", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.dev", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.num", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%num)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.num_bad", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.max_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.max_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.min_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a0.min_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.min", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%min)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.max", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%max)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.mean", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.dev", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.num", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%num)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.num_bad", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.max_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.max_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.min_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a1.min_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.min", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%min)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.max", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%max)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.mean", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.dev", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.num", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%num)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.num_bad", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.max_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.max_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.min_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_cal_coef_a2.min_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.min", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%min)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.max", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%max)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.mean", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.dev", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.num", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%num)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.num_bad", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.max_track", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.max_xtrack", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.min_track", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_bb_raw_noise_counts.min_xtrack", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.min", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%min)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.max", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%max)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.mean", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.dev", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.num", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%num)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.num_bad", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.max_track", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.max_xtrack", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.min_track", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "QA_sv_raw_noise_counts.min_xtrack", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%min_xtrack)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "state1", start, stride, edge, amsua_bt_gran%state1)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "state2", start, stride, edge, amsua_bt_gran%state2)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "cal_coef_a0", start, stride, edge, amsua_bt_gran%cal_coef_a0)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "cal_coef_a0_err", start, stride, edge, amsua_bt_gran%cal_coef_a0_err)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "cal_coef_a1", start, stride, edge, amsua_bt_gran%cal_coef_a1)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "cal_coef_a1_err", start, stride, edge, amsua_bt_gran%cal_coef_a1_err)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "cal_coef_a2", start, stride, edge, amsua_bt_gran%cal_coef_a2)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "cal_coef_a2_err", start, stride, edge, amsua_bt_gran%cal_coef_a2_err)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a1_ColdCalPstion", start, stride, edge, amsua_bt_gran%a1_ColdCalPstion)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a2_ColdCalPstion", start, stride, edge, amsua_bt_gran%a2_ColdCalPstion)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a1_PLO_Redundncy", start, stride, edge, amsua_bt_gran%a1_PLO_Redundncy)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a11_mux_temp_used", start, stride, edge, amsua_bt_gran%a11_mux_temp_used)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a11_receiver_temp", start, stride, edge, amsua_bt_gran%a11_receiver_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a11_target_temp", start, stride, edge, amsua_bt_gran%a11_target_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a12_mux_temp_used", start, stride, edge, amsua_bt_gran%a12_mux_temp_used)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a12_receiver_temp", start, stride, edge, amsua_bt_gran%a12_receiver_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a12_target_temp", start, stride, edge, amsua_bt_gran%a12_target_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a2_diplexer_temp_used", start, stride, edge, amsua_bt_gran%a2_diplexer_temp_used)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a2_receiver_temp", start, stride, edge, amsua_bt_gran%a2_receiver_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "a2_target_temp", start, stride, edge, amsua_bt_gran%a2_target_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "qa_scanline", start, stride, edge, amsua_bt_gran%qa_scanline)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "qa_receiver_a11", start, stride, edge, amsua_bt_gran%qa_receiver_a11)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "qa_receiver_a12", start, stride, edge, amsua_bt_gran%qa_receiver_a12)

edge(1) = AMSUA_BT_GEOTRACK
statn = he5_swrdfld(swid, "qa_receiver_a2", start, stride, edge, amsua_bt_gran%qa_receiver_a2)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = he5_swrdfld(swid, "qa_channel", start, stride, edge, amsua_bt_gran%qa_channel)


! Final clean-up
statn = he5_swdetach(swid)
if (statn /= 0 ) print *, "Error detaching from input file ", file_name
statn = he5_swclose(fid)
if (statn /= 0 ) print *, "Error closing input file ", file_name


return
end subroutine read_amsua_bt_granule

! ----------------------------------------------------------------------
! END of amsua_bt_rdr.f
! ----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module


subroutine read_char_attribute(swid,variable_name,variable,filename)
integer,                    intent(in)    :: swid
character(len=*),           intent(in)    :: variable_name
character(len=*),           intent(inout) :: variable
character(len=*), optional, intent(in)    :: filename

integer :: he5_swrdattr
integer :: statn

statn = he5_swrdattr(swid, variable_name, variable)
if (statn /= 0) then
   print *,'ERROR reading "'//trim(variable_name)//'"'
   if (present(filename)) print *,'         from "'//trim(filename)//'"'
   stop
endif

end subroutine read_char_attribute


subroutine read_integer_attribute(swid,variable_name,variable,filename)
integer,                    intent(in)    :: swid
character(len=*),           intent(in)    :: variable_name
integer,                    intent(inout) :: variable
character(len=*), optional, intent(in)    :: filename

integer :: he5_swrdattr
integer :: statn

statn = he5_swrdattr(swid, variable_name, variable)
if (statn /= 0) then
   print *,'ERROR reading "'//trim(variable_name)//'" from'
   if (present(filename)) print *,'         from "'//trim(filename)//'"'
   stop
endif

end subroutine read_integer_attribute


subroutine read_short_integer_attribute(swid,variable_name,variable,filename)
integer,                    intent(in)    :: swid
character(len=*),           intent(in)    :: variable_name
integer*2,                  intent(inout) :: variable
character(len=*), optional, intent(in)    :: filename

integer :: he5_swrdattr
integer :: statn

statn = he5_swrdattr(swid, variable_name, variable)
if (statn /= 0) then
   print *,'ERROR reading "'//trim(variable_name)//'" from'
   if (present(filename)) print *,'         from "'//trim(filename)//'"'
   stop
endif

end subroutine read_short_integer_attribute


subroutine read_byte_attribute(swid,variable_name,variable,filename)
integer,                    intent(in)    :: swid
character(len=*),           intent(in)    :: variable_name
byte,                       intent(inout) :: variable
character(len=*), optional, intent(in)    :: filename

integer :: he5_swrdattr
integer :: statn

statn = he5_swrdattr(swid, variable_name, variable)
if (statn /= 0) then
   print *,'ERROR reading "'//trim(variable_name)//'" from'
   if (present(filename)) print *,'         from "'//trim(filename)//'"'
   stop
endif

end subroutine read_byte_attribute


subroutine read_real_attribute(swid,variable_name,variable,filename)
integer,                    intent(in)    :: swid
character(len=*),           intent(in)    :: variable_name
real,                       intent(inout) :: variable
character(len=*), optional, intent(in)    :: filename

integer :: he5_swrdattr
integer :: statn

statn = he5_swrdattr(swid, variable_name, variable)
if (statn /= 0) then
   print *,'ERROR reading "'//trim(variable_name)//'" from'
   if (present(filename)) print *,'         from "'//trim(filename)//'"'
   stop
endif

end subroutine read_real_attribute


subroutine read_real8_attribute(swid,variable_name,variable,filename)
integer,                    intent(in)    :: swid
character(len=*),           intent(in)    :: variable_name
real*8,                     intent(inout) :: variable
character(len=*), optional, intent(in)    :: filename

integer :: he5_swrdattr
integer :: statn

statn = he5_swrdattr(swid, variable_name, variable)
if (statn /= 0) then
   print *,'ERROR reading "'//trim(variable_name)//'" from'
   if (present(filename)) print *,'         from "'//trim(filename)//'"'
   stop
endif

end subroutine read_real8_attribute


!------------------------------------------------------------------------------
!  extract the requested AMSUA channel observations from a granule
!  and convert to DART observation format.  allow caller to specify
!  a bounding box and only extract data within that region.

subroutine make_obs_sequence (seq, granule, lon1, lon2, lat1, lat2, &
                             use_channels, scan_thin, pix_thin)

type(obs_sequence_type),    intent(inout) :: seq
type(amsua_bt_gran_type),   intent(in)    :: granule
real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
logical,                    intent(in)    :: use_channels(:)
integer,                    intent(in)    :: scan_thin, pix_thin

type(obs_type) :: obs, prev_obs

integer :: num_copies, num_qc
! max possible obs from this one granule. in practice if the
! real number of processed channels is very much smaller, make
! another parameter so we don't allocate all these unused obs
! (takes time & space) and then delete them at the end.
integer :: max_num

logical :: is_first_obs
type(time_type) :: pre_time

character(len=*), parameter :: routine = 'make_obs_sequence'

if ( .not. module_initialized ) then
    write(string1,*) 'The data has not been read yet.' 
    call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

! one observation data value and one quality control value
! per obs.  if you change these you have to set additional
! metadata for them below.
num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence
max_num = AMSUA_BT_GEOXTRACK * AMSUA_BT_GEOTRACK * AMSUA_BT_CHANNEL

call init_obs_sequence(seq, num_copies, num_qc, max_num)

! set meta data of obs_seq
call set_copy_meta_data(seq, 1, 'observation')
call set_qc_meta_data(seq, 1, 'QC')

! Initialize the obs variables
call init_obs(     obs, 1, 1)
call init_obs(prev_obs, 1, 1)

is_first_obs = .true.

call add_granule_observations(seq,lon1,lon2,lat1,lat2,use_channels,&
                            scan_thin, pix_thin, obs, prev_obs,  &
                            is_first_obs, pre_time, granule)

end subroutine make_obs_sequence

!-----------------------------------------------------------------------------

subroutine add_granule_observations(seq, lon1, lon2, lat1, lat2, use_channels,  &
                            scan_thin, pix_thin, obs, prev_obs, is_first_obs, &
                            pre_time, swath)

type(obs_sequence_type),    intent(inout) :: seq
real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
logical,                    intent(in)    :: use_channels(:)
integer,                    intent(in)    :: scan_thin, pix_thin
type(obs_type),             intent(inout) :: obs, prev_obs
logical,                    intent(inout) :: is_first_obs
type(time_type),            intent(inout) :: pre_time
type(amsua_bt_gran_type),   intent(in)    :: swath

integer :: ipix, iscan, ichan
integer :: days, seconds
integer :: obs_num, key
integer :: which_vert 

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_err
real(r8) :: rqc

real(r8) :: lam1, lam2, phi1, phi2, x, y

real(8) :: sat_az, sat_ze
integer :: platform_id, sat_id, sensor_id
real(8) :: mag_field, cosbk
real(8) :: fastem_p1
real(8) :: fastem_p2
real(8) :: fastem_p3
real(8) :: fastem_p4
real(8) :: fastem_p5

type(time_type) :: obs_time
character(len=*), parameter :: routine = 'add_granule_observations:'

integer :: robstype

! things known (and constant) about the input data and rttov
! consistent with values in rttov_sensor_db.csv

platform_id = 9  ! EOS
sat_id      = 2  ! satellite
sensor_id   = 3  ! AMSUA 

!------------------------------------------------------------------------------
!  loop over all observations within the file

obs_num = 1
which_vert = VERTISUNDEF

! assign each observation the correct observation type
robstype = get_index_for_type_of_obs('EOS_2_AMSUA_TB')
if (robstype < 1) then
   string1 = 'unknown observation type EOS_2_AMSUA_TB'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

!        Latitude(                  AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
!       Longitude(                  AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
! brightness_temp(AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
! rows are along-track, stepping in the direction the satellite is moving
!TJH scanloop:  do iscan=1,swath%nscans
!TJH 
!TJH    ! if we're going to subset rows, we will cycle here
!TJH    if (scan_thin > 0) then
!TJH       if (modulo(iscan, scan_thin) /= 0) cycle scanloop
!TJH    endif
!TJH 
!TJH    ! columns are across-track, varying faster than rows.
!TJH    pixloop:  do ipix=1,swath%npixels
!TJH 
!TJH       ! if we're going to subset columns, ditto
!TJH       if (pix_thin > 0) then
!TJH          if (modulo(ipix, pix_thin) /= 0) cycle pixloop
!TJH       endif
!TJH 
!TJH       ! check channel quality control
!TJH       rqc = swath%Quality(ipix, iscan)
!TJH       ! reject bad scans here - cycle pixloop
!TJH       if (rqc /= 0) cycle pixloop
!TJH 
!TJH       ! observation lat, lon:
!TJH       olat  = swath%Latitude (ipix,iscan) ! valid range [ -90.00,  90.00]
!TJH       olon  = swath%Longitude(ipix,iscan) ! valid range [-180.00, 180.00]
!TJH 
!TJH       ! verify the location is not outside valid limits.  AIRS  uses -180/180
!TJH       if((olon > 180.0_r8) .or. (olon < -180.0_r8) .or.  &
!TJH          (olat >  90.0_r8) .or. (olat <  -90.0_r8)) then
!TJH          write(*,*)'WARNING : invalid location.  col,row,lon,lat = ', ipix,iscan,olon,olat
!TJH          cycle pixloop
!TJH       endif
!TJH 
!TJH       ! reject observations outside the bounding box (allowing wrapping)
!TJH       if(( olat < lat1) .or. ( olat > lat2 ) .or. &
!TJH          (.not. is_longitude_between(olon, lon1, lon2))) cycle pixloop
!TJH 
!TJH       ! make sure lon is between 0 and 360
!TJH       if (olon < 0.0_r8) olon = olon + 360.0_r8
!TJH 
!TJH       ! set the zenith angle (aka earth incidence angle)
!TJH       sat_ze = swath%incidenceAngle(1,ipix,iscan)
!TJH 
!TJH       lam1 = deg2rad*swath%longitude(ipix,iscan)
!TJH       lam2 = deg2rad*swath%SClongitude(iscan)
!TJH       phi1 = deg2rad*swath%latitude(ipix,iscan)
!TJH       phi2 = deg2rad*swath%SClatitude(iscan)
!TJH 
!TJH       ! calculate the bearing between the obs lat/lon and the SClat/lon
!TJH       y = sin(lam2-lam1)*cos(phi2)
!TJH       x = cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lam2-lam1)
!TJH       sat_az = rad2deg*atan2(y,x)
!TJH 
!TJH       obs_time = set_date(int(swath%year(iscan)),int(swath%month(iscan)), &
!TJH          int(swath%DayOfMonth(iscan)), int(swath%Hour(iscan)),            &
!TJH          int(swath%Minute(iscan)), int(swath%Second(iscan)))
!TJH 
!TJH       call get_time(obs_time, seconds, days)
!TJH 
!TJH       channel_loop: do ichan=1, swath%nchannels
!TJH 
!TJH          if (.not. use_channels(ichan+swath%offset)) cycle channel_loop
!TJH 
!TJH          ! create the radiance obs for this observation, add to sequence
!TJH 
!TJH          ! apparently -9999 is missing data, outside of qc mechanism
!TJH          obs_value = swath%Tc(ichan, ipix, iscan)
!TJH          if (obs_value < 0.0_r8) cycle channel_loop
!TJH 
!TJH          obs_err = swath%noiseLevel(ichan)
!TJH 
!TJH          ! column integrated value, so no vertical location
!TJH          vloc = 0.0_r8
!TJH 
!TJH          ! We don't yet have specularity data to add to the observations.
!TJH          if (get_rttov_option_logical('use_zeeman')) then
!TJH             write(string1,*) 'GMI observations do not yet support the Zeeman effect'
!TJH             call error_handler(E_ERR,routine,string1,source,revision,revdate)
!TJH          else
!TJH             mag_field = MISSING_R8
!TJH             cosbk = MISSING_R8
!TJH          end if
!TJH 
!TJH          ! FIXME: consider adding an atlas for looking up FASTEM parameters
!TJH          ! From the RTTOV User guide:
!TJH          ! Surface type
!TJH          ! Typical RTTOV default for land: 3.0, 5.0, 15.0, 0.1, 0.3
!TJH          ! Summer land surface:
!TJH          !        Forest: 1.7, 1.0, 163.0, 0.0, 0.5
!TJH          !        Open grass: 2.2, 1.3, 138.0, 0.0, 0.42
!TJH          !        Bare soil: 2.3, 1.9, 21.8, 0.0, 0.5
!TJH          ! Winter surface type:
!TJH          !        Forest and snow: 2.9, 3.4, 27.0, 0.0, 0.0
!TJH          !        Deep dry snow: 3.0, 24.0, 60.0, 0.1, 0.15
!TJH          !        Frozen soil: 117.8, 2.0, 0.19, 0.2 ,0.35
!TJH          ! Sea ice
!TJH          !        Grease ice: 23.7, 7.7, 17.3, 0.0, 0.15
!TJH          !        Baltic nilas: 1.6, 3.3, 2.2, 0.0, 0.0
!TJH          !        New ice (no snow): 2.9, 3.4, 27.0, 0.0, 0.0
!TJH          !        New ice (snow): 2.2, 3.7, 122.0, 0.0, 0.15
!TJH          !        Brash ice: 3.0, 5.5, 183.0, 0.0, 0.0
!TJH          !        Compact pack ice: 2.0, 1700000.0, 49000000.0, 0.0, 0.0
!TJH          !        Fast ice: 1.5, 77.8, 703.0, 0.1, 0.35
!TJH          !        Lake ice + snow: 1.8, 67.1, 534.0, 0.1, 0.15
!TJH          !        Multi-year ice: 1.5, 85000.0, 4700000.0, 0.0, 0.0
!TJH          fastem_p1 = 3.0d0
!TJH          fastem_p2 = 5.0d0
!TJH          fastem_p3 = 15.0d0
!TJH          fastem_p4 = 0.1d0
!TJH          fastem_p5 = 0.3d0
!TJH 
!TJH          ! add additional metadata for this obs type.  returns key to use in create call
!TJH          call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, & 
!TJH             ichan+swath%offset, mag_field, cosbk, fastem_p1, fastem_p2, fastem_p3, &
!TJH             fastem_p4, fastem_p5)
!TJH 
!TJH          call create_3d_obs(olat, olon, vloc, which_vert, obs_value, robstype, &
!TJH                             obs_err, days, seconds, rqc, obs, key)
!TJH 
!TJH          call add_obs_to_seq(seq, obs, obs_time, prev_obs, pre_time, is_first_obs)
!TJH 
!TJH          obs_num = obs_num + 1
!TJH 
!TJH       enddo channel_loop
!TJH    enddo pixloop
!TJH enddo scanloop
!TJH 
!TJH ! Print a little summary
!TJH call print_obs_seq_summary(seq)
!TJH 
!TJH write(string1,*) 'Converted ',obs_num-1,' obs for swath ',swath%dset_prefix, &
!TJH                    '; total AMSUA obs = ',key
!TJH call error_handler(E_MSG, routine, string1, source, revision, revdate)

end subroutine add_granule_observations



subroutine dump_attributes(granule)
type (amsua_bt_gran_type), intent(in) :: granule

print *,' AutomaticQAFlag ', granule%AutomaticQAFlag
print *,' DayNightFlag ', granule%DayNightFlag
print *,' LatGranuleCen ', granule%LatGranuleCen
print *,' LocTimeGranuleCen ', granule%LocTimeGranuleCen
print *,' LonGranuleCen ', granule%LonGranuleCen
print *,' MoonInViewMWCount ', granule%MoonInViewMWCount
print *,' NumBadData ', granule%NumBadData
print *,' NumLandSurface ', granule%NumLandSurface
print *,' NumMissingData ', granule%NumMissingData
print *,' NumOceanSurface ', granule%NumOceanSurface
print *,' NumProcessData ', granule%NumProcessData
print *,' NumSpecialData ', granule%NumSpecialData
print *,' NumTotalData ', granule%NumTotalData
print *,' eq_x_longitude ', granule%eq_x_longitude
print *,' eq_x_tai ', granule%eq_x_tai
print *,' granule_number', granule%granule_number
print *,' instrument ', granule%instrument
print *,' node_type ', granule%node_type
print *,' num_data_gaps_a1 ', granule%num_data_gaps_a1
print *,' num_data_gaps_a2 ', granule%num_data_gaps_a2
print *,' num_demgeoqa ', granule%num_demgeoqa
print *,' num_fpe ', granule%num_fpe
print *,' num_ftptgeoqa ', granule%num_ftptgeoqa
print *,' num_glintgeoqa ', granule%num_glintgeoqa
print *,' num_instr_mode_changes_a1 ', granule%num_instr_mode_changes_a1
print *,' num_instr_mode_changes_a2 ', granule%num_instr_mode_changes_a2
print *,' num_missing_scanlines_a1 ', granule%num_missing_scanlines_a1
print *,' num_missing_scanlines_a2 ', granule%num_missing_scanlines_a2
print *,' num_moongeoqa ', granule%num_moongeoqa
print *,' num_satgeoqa ', granule%num_satgeoqa
print *,' num_scanlines ', granule%num_scanlines
print *,' num_scanlines_not_norm_mode_a1 ', granule%num_scanlines_not_norm_mode_a1
print *,' num_scanlines_not_norm_mode_a2 ', granule%num_scanlines_not_norm_mode_a2
print *,' num_scanlines_rec_cal_prob_a11 ', granule%num_scanlines_rec_cal_prob_a11
print *,' num_scanlines_rec_cal_prob_a12 ', granule%num_scanlines_rec_cal_prob_a12
print *,' num_scanlines_rec_cal_prob_a2 ', granule%num_scanlines_rec_cal_prob_a2
print *,' num_scanlines_sig_coast_xing ', granule%num_scanlines_sig_coast_xing
print *,' num_scanlines_sig_sun_glint ', granule%num_scanlines_sig_sun_glint
print *,' num_scansets ', granule%num_scansets
print *,' num_zengeoqa ', granule%num_zengeoqa
print *,' orbit_path ', granule%orbit_path
print *,' orbitgeoqa ', granule%orbitgeoqa
print *,' processing_level ', granule%processing_level
print *,' start_Latitude ', granule%start_Latitude
print *,' end_Latitude ', granule%end_Latitude
print *,' start_Longitude ', granule%start_Longitude
print *,' end_Longitude ', granule%end_Longitude
print *,' start_Time ', granule%start_Time
print *,' end_Time ', granule%end_Time
print *,' start_year ', granule%start_year
print *,' start_month ', granule%start_month
print *,' start_day ', granule%start_day
print *,' start_hour ', granule%start_hour
print *,' start_minute ', granule%start_minute
print *,' start_sec ', granule%start_sec
print *,' start_orbit ', granule%start_orbit
print *,' end_orbit ', granule%end_orbit
print *,' start_orbit_row ', granule%start_orbit_row
print *,' end_orbit_row ', granule%end_orbit_row

end subroutine dump_attributes

end module

