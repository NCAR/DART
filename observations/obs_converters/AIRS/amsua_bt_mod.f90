
module msua_bt_mod

! the contents of this file are an amalgam of:
!    amsua_bt_typ.inc
!    amsua_bt_struct.inc
! modified to use fortran 90

use types_mod, only : r4, r8, digits12

!-------------------------------------------------------------------------------
! start of amsua_bt_typ.inc
!-------------------------------------------------------------------------------

integer, parameter :: AMSUA_BT_GEOXTRACK    = 30
integer, parameter :: AMSUA_BT_GEOTRACK     = 45
integer, parameter :: AMSUA_BT_CHANNEL      = 15
integer, parameter :: AMSUA_BT_CALXTRACK    =  4
integer, parameter :: AMSUA_BT_SPACEXTRACK  =  2
integer, parameter :: AMSUA_BT_BBXTRACK     =  2
integer, parameter :: AMSUA_BT_WARMPRTA11   =  5
integer, parameter :: AMSUA_BT_WARMPRTA12   =  5
integer, parameter :: AMSUA_BT_WARMPRTA2    =  7

! Type definitions for each record type

TYPE  amsua_bt_QA_bb_PRT_a11_t
  real         min
  real         max
  real         mean
  real         dev
  integer      num_in
  integer      num_lo
  integer      num_hi
  integer      num_bad
  real         range_min
  real         range_max
  byte         missing
  integer      max_track
  integer      max_xtrack
  integer      min_track
  integer      min_xtrack
END TYPE


TYPE  amsua_bt_QA_bb_PRT_a12_t
  real         min
  real         max
  real         mean
  real         dev
  integer      num_in
  integer      num_lo
  integer      num_hi
  integer      num_bad
  real         range_min
  real         range_max
  byte         missing
  integer      max_track
  integer      max_xtrack
  integer      min_track
  integer      min_xtrack
END TYPE


TYPE  amsua_bt_QA_bb_PRT_a2_t
  real         min
  real         max
  real         mean
  real         dev
  integer      num_in
  integer      num_lo
  integer      num_hi
  integer      num_bad
  real         range_min
  real         range_max
  byte         missing
  integer      max_track
  integer      max_xtrack
  integer      min_track
  integer      min_xtrack
END TYPE


TYPE  amsua_bt_QA_rec_PRT_a11_t
  real         min
  real         max
  real         mean
  real         dev
  integer      num_in
  integer      num_lo
  integer      num_hi
  integer      num_bad
  real         range_min
  real         range_max
  byte         missing
  integer      max_track
  integer      max_xtrack
  integer      min_track
  integer      min_xtrack
END TYPE


TYPE  amsua_bt_QA_rec_PRT_a12_t
  real         min
  real         max
  real         mean
  real         dev
  integer      num_in
  integer      num_lo
  integer      num_hi
  integer      num_bad
  real         range_min
  real         range_max
  byte         missing
  integer      max_track
  integer      max_xtrack
  integer      min_track
  integer      min_xtrack
END TYPE


TYPE  amsua_bt_QA_rec_PRT_a2_t

  real         min
  real         max
  real         mean
  real         dev
  integer      num_in
  integer      num_lo
  integer      num_hi
  integer      num_bad
  real         range_min
  real         range_max
  byte         missing
  integer      max_track
  integer      max_xtrack
  integer      min_track
  integer      min_xtrack
END TYPE


TYPE  amsua_bt_bb_signals_t
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

!-------------------------------------------------------------------------------
! end of amsua_bt_typ.inc
!-------------------------------------------------------------------------------
! start of of amsua_bt_struct.inc
!-------------------------------------------------------------------------------

! Record holds an entire granule of amsua_bt
TYPE  amsua_bt_gran_t

  character(len=256) :: instrument
  character(len=256) :: DayNightFlag
  character(len=256) :: AutomaticQAFlag
  integer      NumTotalData
  integer      NumProcessData
  integer      NumSpecialData
  integer      NumBadData
  integer      NumMissingData
  integer      NumLandSurface
  integer      NumOceanSurface
  character(len=256) :: node_type
  integer      start_year
  integer      start_month
  integer      start_day
  integer      start_hour
  integer      start_minute
  real   start_sec
  integer      start_orbit
  integer      end_orbit
  integer      orbit_path
  integer      start_orbit_row
  integer      end_orbit_row
  integer      granule_number
  integer      num_scansets
  integer      num_scanlines
  real(r8) start_Latitude
  real(r8) start_Longitude
  real(r8) start_Time
  real(r8) end_Latitude
  real(r8) end_Longitude
  real(r8) end_Time
  real   eq_x_longitude
  real(r8) eq_x_tai
  integer      orbitgeoqa
  integer*2  num_satgeoqa
  integer*2  num_glintgeoqa
  integer*2  num_moongeoqa
  integer*2  num_ftptgeoqa
  integer*2  num_zengeoqa
  integer*2  num_demgeoqa
  integer*2  num_fpe
  integer*2  LonGranuleCen
  integer*2  LatGranuleCen
  integer*2  LocTimeGranuleCen
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

TYPE(amsua_bt_QA_bb_PRT_a11_t) QA_bb_PRT_a11

TYPE(amsua_bt_QA_bb_PRT_a12_t) QA_bb_PRT_a12

TYPE(amsua_bt_QA_bb_PRT_a2_t) QA_bb_PRT_a2

TYPE(amsua_bt_QA_rec_PRT_a11_t) QA_rec_PRT_a11

TYPE(amsua_bt_QA_rec_PRT_a12_t) QA_rec_PRT_a12

TYPE(amsua_bt_QA_rec_PRT_a2_t) QA_rec_PRT_a2

  character(len=256) ::   granules_present

  ! Geolocation fields
  real(r8) Latitude(AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r8) Longitude(AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r8) Time(AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)

  ! Data Fields
  real   scanang( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   satheight( AMSUA_BT_GEOTRACK)
  real   satroll( AMSUA_BT_GEOTRACK)
  real   satpitch( AMSUA_BT_GEOTRACK)
  real   satyaw( AMSUA_BT_GEOTRACK)
  integer  satgeoqa( AMSUA_BT_GEOTRACK)
  integer*2  glintgeoqa( AMSUA_BT_GEOTRACK)
  integer*2  moongeoqa( AMSUA_BT_GEOTRACK)
  integer  ftptgeoqa( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  integer*2  zengeoqa( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  integer*2  demgeoqa( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r8)         nadirTAI( AMSUA_BT_GEOTRACK)
  real(r8)         sat_lat( AMSUA_BT_GEOTRACK)
  real(r8)         sat_lon( AMSUA_BT_GEOTRACK)
  byte   scan_node_type( AMSUA_BT_GEOTRACK)
  real   satzen( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   satazi( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   solzen( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   solazi( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   glintlat( AMSUA_BT_GEOTRACK)
  real   glintlon( AMSUA_BT_GEOTRACK)
  integer*2  sun_glint_distance( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   topog( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   topog_err( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   landFrac( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   landFrac_err( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   antenna_temp( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   brightness_temp( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   brightness_temp_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real   center_freq( AMSUA_BT_CHANNEL)
  real   IF_offset_1( AMSUA_BT_CHANNEL)
  real   IF_offset_2( AMSUA_BT_CHANNEL)
  real   bandwidth( AMSUA_BT_CHANNEL)
  integer  num_calibrated_scanlines( AMSUA_BT_CHANNEL)
  integer  num_scanlines_ch_cal_problems( AMSUA_BT_CHANNEL)
TYPE(amsua_bt_bb_signals_t) bb_signals

TYPE(amsua_bt_space_signals_t) space_signals

TYPE(amsua_bt_gain_stats_t) gain_stats

  real   NeDT( AMSUA_BT_CHANNEL)
TYPE(amsua_bt_QA_unfiltered_scene_count_t) QA_unfiltered_scene_count

TYPE(amsua_bt_QA_unfiltered_BB_count_t) QA_unfiltered_BB_count

TYPE(amsua_bt_QA_unfiltered_space_count_t) QA_unfiltered_space_count

TYPE(amsua_bt_QA_cal_coef_a0_t) QA_cal_coef_a0

TYPE(amsua_bt_QA_cal_coef_a1_t) QA_cal_coef_a1

TYPE(amsua_bt_QA_cal_coef_a2_t) QA_cal_coef_a2

TYPE(amsua_bt_QA_bb_raw_noise_counts_t) QA_bb_raw_noise_counts

TYPE(amsua_bt_QA_sv_raw_noise_counts_t) QA_sv_raw_noise_counts
  integer  state1( AMSUA_BT_GEOTRACK)
  integer  state2( AMSUA_BT_GEOTRACK)
  real   cal_coef_a0( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real   cal_coef_a0_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real   cal_coef_a1( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real   cal_coef_a1_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real   cal_coef_a2( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real   cal_coef_a2_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  byte   a1_ColdCalPstion( AMSUA_BT_GEOTRACK)
  byte   a2_ColdCalPstion( AMSUA_BT_GEOTRACK)
  byte   a1_PLO_Redundncy( AMSUA_BT_GEOTRACK)
  byte   a11_mux_temp_used( AMSUA_BT_GEOTRACK)
  real   a11_receiver_temp( AMSUA_BT_GEOTRACK)
  real   a11_target_temp( AMSUA_BT_GEOTRACK)
  byte   a12_mux_temp_used( AMSUA_BT_GEOTRACK)
  real   a12_receiver_temp( AMSUA_BT_GEOTRACK)
  real   a12_target_temp( AMSUA_BT_GEOTRACK)
  byte   a2_diplexer_temp_used( AMSUA_BT_GEOTRACK)
  real   a2_receiver_temp( AMSUA_BT_GEOTRACK)
  real   a2_target_temp( AMSUA_BT_GEOTRACK)
  byte   qa_scanline( AMSUA_BT_GEOTRACK)
  byte   qa_receiver_a11( AMSUA_BT_GEOTRACK)
  byte   qa_receiver_a12( AMSUA_BT_GEOTRACK)
  byte   qa_receiver_a2( AMSUA_BT_GEOTRACK)
  byte   qa_channel( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
END TYPE

!-------------------------------------------------------------------------------
! end of of amsua_bt_struct.inc
!-------------------------------------------------------------------------------

end module amsua_bt_mod
