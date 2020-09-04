      integer   AMSUA_BT_GEOXTRACK                      
      parameter(AMSUA_BT_GEOXTRACK                      =    30)
      integer   AMSUA_BT_GEOTRACK                       
      parameter(AMSUA_BT_GEOTRACK                       =    45)
      integer   AMSUA_BT_CHANNEL                        
      parameter(AMSUA_BT_CHANNEL                        =    15)
      integer   AMSUA_BT_CALXTRACK                      
      parameter(AMSUA_BT_CALXTRACK                      =     4)
      integer   AMSUA_BT_SPACEXTRACK                    
      parameter(AMSUA_BT_SPACEXTRACK                    =     2)
      integer   AMSUA_BT_BBXTRACK                       
      parameter(AMSUA_BT_BBXTRACK                       =     2)
      integer   AMSUA_BT_WARMPRTA11                     
      parameter(AMSUA_BT_WARMPRTA11                     =     5)
      integer   AMSUA_BT_WARMPRTA12                     
      parameter(AMSUA_BT_WARMPRTA12                     =     5)
      integer   AMSUA_BT_WARMPRTA2                      
      parameter(AMSUA_BT_WARMPRTA2                      =     7)

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

        real           min( 15, 2)
        real           max( 15, 2)
        real           mean( 15, 2)
        real           dev( 15, 2)
        integer        num( 15, 2)
        integer        num_bad( 15, 2)
        integer        max_track( 15, 2)
        integer        max_xtrack( 15, 2)
        integer        min_track( 15, 2)
        integer        min_xtrack( 15, 2)
      END TYPE


      TYPE  amsua_bt_space_signals_t

! Attributes

        real           min( 15, 2)
        real           max( 15, 2)
        real           mean( 15, 2)
        real           dev( 15, 2)
        integer        num( 15, 2)
        integer        num_bad( 15, 2)
        integer        max_track( 15, 2)
        integer        max_xtrack( 15, 2)
        integer        min_track( 15, 2)
        integer        min_xtrack( 15, 2)
      END TYPE


      TYPE  amsua_bt_gain_stats_t

! Attributes

        real           min( 15)
        real           max( 15)
        real           mean( 15)
        real           dev( 15)
        integer        num( 15)
        integer        num_bad( 15)
        integer        max_track( 15)
        integer        max_xtrack( 15)
        integer        min_track( 15)
        integer        min_xtrack( 15)
      END TYPE


      TYPE  amsua_bt_QA_unfiltered_scene_count_t

! Attributes

        real           min( 15, 30)
        real           max( 15, 30)
        real           mean( 15, 30)
        real           dev( 15, 30)
        integer        num( 15, 30)
        integer        num_bad( 15, 30)
        integer        max_track( 15, 30)
        integer        max_xtrack( 15, 30)
        integer        min_track( 15, 30)
        integer        min_xtrack( 15, 30)
      END TYPE


      TYPE  amsua_bt_QA_unfiltered_BB_count_t

! Attributes

        real           min( 15, 2)
        real           max( 15, 2)
        real           mean( 15, 2)
        real           dev( 15, 2)
        integer        num( 15, 2)
        integer        num_bad( 15, 2)
        integer        max_track( 15, 2)
        integer        max_xtrack( 15, 2)
        integer        min_track( 15, 2)
        integer        min_xtrack( 15, 2)
      END TYPE


      TYPE  amsua_bt_QA_unfiltered_space_count_t

! Attributes

        real           min( 15, 2)
        real           max( 15, 2)
        real           mean( 15, 2)
        real           dev( 15, 2)
        integer        num( 15, 2)
        integer        num_bad( 15, 2)
        integer        max_track( 15, 2)
        integer        max_xtrack( 15, 2)
        integer        min_track( 15, 2)
        integer        min_xtrack( 15, 2)
      END TYPE


      TYPE  amsua_bt_QA_cal_coef_a0_t

! Attributes

        real           min( 15)
        real           max( 15)
        real           mean( 15)
        real           dev( 15)
        integer        num( 15)
        integer        num_bad( 15)
        integer        max_track( 15)
        integer        max_xtrack( 15)
        integer        min_track( 15)
        integer        min_xtrack( 15)
      END TYPE


      TYPE  amsua_bt_QA_cal_coef_a1_t

! Attributes

        real           min( 15)
        real           max( 15)
        real           mean( 15)
        real           dev( 15)
        integer        num( 15)
        integer        num_bad( 15)
        integer        max_track( 15)
        integer        max_xtrack( 15)
        integer        min_track( 15)
        integer        min_xtrack( 15)
      END TYPE


      TYPE  amsua_bt_QA_cal_coef_a2_t

! Attributes

        real           min( 15)
        real           max( 15)
        real           mean( 15)
        real           dev( 15)
        integer        num( 15)
        integer        num_bad( 15)
        integer        max_track( 15)
        integer        max_xtrack( 15)
        integer        min_track( 15)
        integer        min_xtrack( 15)
      END TYPE


      TYPE  amsua_bt_QA_bb_raw_noise_counts_t

! Attributes

        real           min( 15)
        real           max( 15)
        real           mean( 15)
        real           dev( 15)
        integer        num( 15)
        integer        num_bad( 15)
        integer        max_track( 15)
        integer        max_xtrack( 15)
        integer        min_track( 15)
        integer        min_xtrack( 15)
      END TYPE


      TYPE  amsua_bt_QA_sv_raw_noise_counts_t

! Attributes

        real           min( 15)
        real           max( 15)
        real           mean( 15)
        real           dev( 15)
        integer        num( 15)
        integer        num_bad( 15)
        integer        max_track( 15)
        integer        max_xtrack( 15)
        integer        min_track( 15)
        integer        min_xtrack( 15)
      END TYPE



! Record holds an entire granule of amsua_bt
      TYPE  amsua_bt_gran_t

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
        real   start_sec
        integer      start_orbit
        integer      end_orbit
        integer      orbit_path
        integer      start_orbit_row
        integer      end_orbit_row
        integer      granule_number
        integer      num_scansets
        integer      num_scanlines
        double precision start_Latitude
        double precision start_Longitude
        double precision start_Time
        double precision end_Latitude
        double precision end_Longitude
        double precision end_Time
        real   eq_x_longitude
        double precision eq_x_tai
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

        character*256  granules_present

! Geolocation fields
        double precision Latitude(AMSUA_BT_GEOXTRACK,
     &                          AMSUA_BT_GEOTRACK)
        double precision Longitude(AMSUA_BT_GEOXTRACK,
     &                           AMSUA_BT_GEOTRACK)
        double precision Time(AMSUA_BT_GEOXTRACK,
     &                      AMSUA_BT_GEOTRACK)

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
        double precision         nadirTAI( AMSUA_BT_GEOTRACK)
        double precision         sat_lat( AMSUA_BT_GEOTRACK)
        double precision         sat_lon( AMSUA_BT_GEOTRACK)
        byte   scan_node_type( AMSUA_BT_GEOTRACK)
        real   satzen( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   satazi( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   solzen( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   solazi( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   glintlat( AMSUA_BT_GEOTRACK)
        real   glintlon( AMSUA_BT_GEOTRACK)
        integer*2  sun_glint_distance( AMSUA_BT_GEOXTRACK, 
     &	AMSUA_BT_GEOTRACK)
        real   topog( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   topog_err( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   landFrac( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   landFrac_err( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
        real   antenna_temp( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, 
     &	AMSUA_BT_GEOTRACK)
        real   brightness_temp( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, 
     &	AMSUA_BT_GEOTRACK)
        real   brightness_temp_err( AMSUA_BT_CHANNEL, 
     &	AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
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
      TYPE(amsua_bt_QA_unfiltered_scene_count_t) 
     &      QA_unfiltered_scene_count

      TYPE(amsua_bt_QA_unfiltered_BB_count_t) QA_unfiltered_BB_count

      TYPE(amsua_bt_QA_unfiltered_space_count_t) 
     &      QA_unfiltered_space_count

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

