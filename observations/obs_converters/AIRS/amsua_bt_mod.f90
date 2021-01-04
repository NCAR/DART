
module amsua_bt_mod

! the contents of this file are a consolidation of:
!    amsua_bt_typ.inc
!    amsua_bt_struct.inc
! modified to use fortran 90, removed redundant type declarations

use types_mod, only : i2, r4, digits12

implicit none
private

public :: amsua_bt_granule, amsua_bt_rdr, &
          AMSUA_BT_GEOXTRACK,  AMSUA_BT_GEOTRACK,    AMSUA_BT_CHANNEL, &
          AMSUA_BT_CALXTRACK,  AMSUA_BT_SPACEXTRACK, AMSUA_BT_BBXTRACK, &
          AMSUA_BT_WARMPRTA11, AMSUA_BT_WARMPRTA12,  AMSUA_BT_WARMPRTA2

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

TYPE  amsua_bt_record
   real(r4) :: min
   real(r4) :: max
   real(r4) :: mean
   real(r4) :: dev
   integer  :: num_in
   integer  :: num_lo
   integer  :: num_hi
   integer  :: num_bad
   real(r4) :: range_min
   real(r4) :: range_max
   byte     :: missing
   integer  :: max_track
   integer  :: max_xtrack
   integer  :: min_track
   integer  :: min_xtrack
END TYPE

TYPE  amsua_bt_bb_signals
   real(r4) :: min(        AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   real(r4) :: max(        AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   real(r4) :: mean(       AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   real(r4) :: dev(        AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   integer  :: num(        AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   integer  :: num_bad(    AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   integer  :: max_track(  AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   integer  :: max_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   integer  :: min_track(  AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
   integer  :: min_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_BBXTRACK)
END TYPE

TYPE  amsua_bt_space_signals
   real(r4) :: min(        AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   real(r4) :: max(        AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   real(r4) :: mean(       AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   real(r4) :: dev(        AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   integer  :: num(        AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   integer  :: num_bad(    AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   integer  :: max_track(  AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   integer  :: max_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   integer  :: min_track(  AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
   integer  :: min_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_SPACEXTRACK)
END TYPE

TYPE  amsua_bt_scene_count
   real(r4) :: min(        AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   real(r4) :: max(        AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   real(r4) :: mean(       AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   real(r4) :: dev(        AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   integer  :: num(        AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   integer  :: num_bad(    AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   integer  :: max_track(  AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   integer  :: max_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   integer  :: min_track(  AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
   integer  :: min_xtrack( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK)
END TYPE

TYPE  amsua_bt_stats
   real(r4) :: min(        AMSUA_BT_CHANNEL)
   real(r4) :: max(        AMSUA_BT_CHANNEL)
   real(r4) :: mean(       AMSUA_BT_CHANNEL)
   real(r4) :: dev(        AMSUA_BT_CHANNEL)
   integer  :: num(        AMSUA_BT_CHANNEL)
   integer  :: num_bad(    AMSUA_BT_CHANNEL)
   integer  :: max_track(  AMSUA_BT_CHANNEL)
   integer  :: max_xtrack( AMSUA_BT_CHANNEL)
   integer  :: min_track(  AMSUA_BT_CHANNEL)
   integer  :: min_xtrack( AMSUA_BT_CHANNEL)
END TYPE

!-------------------------------------------------------------------------------
! end of amsua_bt_typ.inc
!-------------------------------------------------------------------------------
! start of of amsua_bt_struct.inc
!-------------------------------------------------------------------------------

! Record holds an entire granule of amsua_bt
TYPE  amsua_bt_granule

   ! Attributes
  character(len=256) :: processing_level
  character(len=256) :: instrument
  character(len=256) :: DayNightFlag
  character(len=256) :: AutomaticQAFlag
  integer            :: NumTotalData
  integer            :: NumProcessData
  integer            :: NumSpecialData
  integer            :: NumBadData
  integer            :: NumMissingData
  integer            :: NumLandSurface
  integer            :: NumOceanSurface
  character(len=256) :: node_type
  integer            :: start_year
  integer            :: start_month
  integer            :: start_day
  integer            :: start_hour
  integer            :: start_minute
  real(r4)           :: start_sec
  integer            :: start_orbit
  integer            :: end_orbit
  integer            :: orbit_path
  integer            :: start_orbit_row
  integer            :: end_orbit_row
  integer            :: granule_number
  integer            :: num_scansets
  integer            :: num_scanlines
  real(digits12)     :: start_Latitude
  real(digits12)     :: start_Longitude
  real(digits12)     :: start_Time
  real(digits12)     :: end_Latitude
  real(digits12)     :: end_Longitude
  real(digits12)     :: end_Time
  real(r4)           :: eq_x_longitude
  real(digits12)     :: eq_x_tai
  integer            :: orbitgeoqa
  integer(i2)        :: num_satgeoqa
  integer(i2)        :: num_glintgeoqa
  integer(i2)        :: num_moongeoqa
  integer(i2)        :: num_ftptgeoqa
  integer(i2)        :: num_zengeoqa
  integer(i2)        :: num_demgeoqa
  integer(i2)        :: num_fpe
  integer(i2)        :: LonGranuleCen
  integer(i2)        :: LatGranuleCen
  integer(i2)        :: LocTimeGranuleCen
  integer            :: num_scanlines_not_norm_mode_a1
  integer            :: num_scanlines_not_norm_mode_a2
  integer            :: num_missing_scanlines_a1
  integer            :: num_missing_scanlines_a2
  integer            :: num_data_gaps_a1
  integer            :: num_data_gaps_a2
  integer            :: num_instr_mode_changes_a1
  integer            :: num_instr_mode_changes_a2
  integer            :: num_scanlines_rec_cal_prob_a11
  integer            :: num_scanlines_rec_cal_prob_a12
  integer            :: num_scanlines_rec_cal_prob_a2
  integer            :: num_scanlines_sig_coast_xing
  integer            :: num_scanlines_sig_sun_glint
  integer            :: MoonInViewMWCount

  TYPE(amsua_bt_record) :: QA_bb_PRT_a11
  TYPE(amsua_bt_record) :: QA_bb_PRT_a12
  TYPE(amsua_bt_record) :: QA_bb_PRT_a2
  TYPE(amsua_bt_record) :: QA_rec_PRT_a11
  TYPE(amsua_bt_record) :: QA_rec_PRT_a12
  TYPE(amsua_bt_record) :: QA_rec_PRT_a2

  character(len=256) :: granules_present

  ! Geolocation fields
  ! Latitude ... geodetic latitude, degrees N [-90,90]
  ! Longitude .. geodetic longitude, degrees E [-180,180]
  ! Time     ... 'shutter' TAI Time, floating point elapsed seconds since 1 Jan 1993

  real(digits12) :: latitude( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(digits12) :: longitude(AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(digits12) :: Time(     AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)

  ! Data Fields

  real(r4)       :: scanang(            AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: satheight(          AMSUA_BT_GEOTRACK)
  real(r4)       :: satroll(            AMSUA_BT_GEOTRACK)
  real(r4)       :: satpitch(           AMSUA_BT_GEOTRACK)
  real(r4)       :: satyaw(             AMSUA_BT_GEOTRACK)
  integer        :: satgeoqa(           AMSUA_BT_GEOTRACK)
  integer(i2)    :: glintgeoqa(         AMSUA_BT_GEOTRACK)
  integer(i2)    :: moongeoqa(          AMSUA_BT_GEOTRACK)
  integer        :: ftptgeoqa(          AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  integer(i2)    :: zengeoqa(           AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  integer(i2)    :: demgeoqa(           AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(digits12) :: nadirTAI(           AMSUA_BT_GEOTRACK)
  real(digits12) :: sat_lat(            AMSUA_BT_GEOTRACK)
  real(digits12) :: sat_lon(            AMSUA_BT_GEOTRACK)
  byte           :: scan_node_type(     AMSUA_BT_GEOTRACK)
  real(r4)       :: satzen(             AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: satazi(             AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: solzen(             AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: solazi(             AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: glintlat(           AMSUA_BT_GEOTRACK)
  real(r4)       :: glintlon(           AMSUA_BT_GEOTRACK)
  integer(i2)    :: sun_glint_distance( AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: topog(              AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: topog_err(          AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: landFrac(           AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: landFrac_err(       AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: antenna_temp(        AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: brightness_temp(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
  real(r4)       :: brightness_temp_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)

  real(r4)       :: center_freq(                   AMSUA_BT_CHANNEL) ! center channel frequency (GHz)
  real(r4)       :: IF_offset_1(                   AMSUA_BT_CHANNEL) ! offset of first intermediate frequency (MHz)
  real(r4)       :: IF_offset_2(                   AMSUA_BT_CHANNEL) ! offset of second    ..       frequency (MHz)
  real(r4)       :: bandwidth(                     AMSUA_BT_CHANNEL) ! bandwith of sum of 1,3, or 4 channels (MHz)
  integer        :: num_calibrated_scanlines(      AMSUA_BT_CHANNEL)
  integer        :: num_scanlines_ch_cal_problems( AMSUA_BT_CHANNEL)
  real(r4)       :: NeDT(AMSUA_BT_CHANNEL) ! instrument noise level estimated from warm count scatter

  TYPE(amsua_bt_bb_signals)    :: bb_signals
  TYPE(amsua_bt_space_signals) :: space_signals
  TYPE(amsua_bt_stats)         :: gain_stats
  TYPE(amsua_bt_scene_count)   :: QA_unfiltered_scene_count
  TYPE(amsua_bt_bb_signals)    :: QA_unfiltered_BB_count
  TYPE(amsua_bt_space_signals) :: QA_unfiltered_space_count
  TYPE(amsua_bt_stats)         :: QA_cal_coef_a0
  TYPE(amsua_bt_stats)         :: QA_cal_coef_a1
  TYPE(amsua_bt_stats)         :: QA_cal_coef_a2
  TYPE(amsua_bt_stats)         :: QA_bb_raw_noise_counts
  TYPE(amsua_bt_stats)         :: QA_sv_raw_noise_counts

  integer  :: state1( AMSUA_BT_GEOTRACK)
  integer  :: state2( AMSUA_BT_GEOTRACK)
  real(r4) :: cal_coef_a0(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real(r4) :: cal_coef_a0_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real(r4) :: cal_coef_a1(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real(r4) :: cal_coef_a1_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real(r4) :: cal_coef_a2(     AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  real(r4) :: cal_coef_a2_err( AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
  byte     :: a1_ColdCalPstion(      AMSUA_BT_GEOTRACK)
  byte     :: a2_ColdCalPstion(      AMSUA_BT_GEOTRACK)
  byte     :: a1_PLO_Redundncy(      AMSUA_BT_GEOTRACK)
  byte     :: a11_mux_temp_used(     AMSUA_BT_GEOTRACK)
  real(r4) :: a11_receiver_temp(     AMSUA_BT_GEOTRACK)
  real(r4) :: a11_target_temp(       AMSUA_BT_GEOTRACK)
  byte     :: a12_mux_temp_used(     AMSUA_BT_GEOTRACK)
  real(r4) :: a12_receiver_temp(     AMSUA_BT_GEOTRACK)
  real(r4) :: a12_target_temp(       AMSUA_BT_GEOTRACK)
  byte     :: a2_diplexer_temp_used( AMSUA_BT_GEOTRACK)
  real(r4) :: a2_receiver_temp(      AMSUA_BT_GEOTRACK)
  real(r4) :: a2_target_temp(        AMSUA_BT_GEOTRACK)
  byte     :: qa_scanline(           AMSUA_BT_GEOTRACK)
  byte     :: qa_receiver_a11(       AMSUA_BT_GEOTRACK)
  byte     :: qa_receiver_a12(       AMSUA_BT_GEOTRACK)
  byte     :: qa_receiver_a2(        AMSUA_BT_GEOTRACK)
  byte     :: qa_channel(            AMSUA_BT_CHANNEL, AMSUA_BT_GEOTRACK)
END TYPE

!-------------------------------------------------------------------------------
! end of amsua_bt_struct.inc
!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------
! start of amsua_bt_rdr.f90
!-------------------------------------------------------------------------------

! This function was based on the output of the mkezio program to read
! an AIRS swath of type "L1B_AMSU" from file given by the
! file_name argument into a buffer pointed to by the amsua_bt_gran
! argument.  The caller owns the buffer.  The entire granule
! is read -- every attribute and field, the whole lat/lon/time
! extent.
!
! Errors opening the file, etc. are fatal and cause STOP.
! Problems reading individual attributes or fields are reported to
! the console but do not interrupt program flow.

subroutine amsua_bt_rdr(file_name, amsua_bt_gran)

character(len=*),       intent(in)  :: file_name
type(amsua_bt_granule), intent(out) :: amsua_bt_gran

integer :: statn                  ! HDF-EOS status. 0 for success
integer :: fid                    ! HDF-EOS file ID
integer :: swid                   ! HDF-EOS swath ID
integer :: nchar                  ! Number of characters
character(len=256) :: swathname   ! Name of swath
integer :: nswath                 ! Number of swaths
integer :: start(10) = 0          ! start of each dimensions for Swath I/O
                                  ! 0 => start with first element
integer :: stride(10) = 1         ! stride of each dimensions for Swath I/O
                                  ! 1 => use every element
integer :: edge(10)               ! size of each dimension for swath I/O
                                  ! will be set for each individual read
integer :: swopen, swinqswath, swattach
integer :: swrdfld, swrdattr
integer :: swdetach, swclose

fid = swopen(file_name, 1)
if (fid .eq. -1) then
  print *, "Error ", fid, " opening file ", file_name
  stop
end if

! Get name of swath(s)
nswath = swinqswath(file_name, swathname, nchar)
if (nswath .ne. 1) then
  print *, "swinqswath found ", nswath, " swaths for file ", file_name, " Need exactly 1"
  stop
end if

! There's exactly one swath.  Make sure it is the right one.
if (swathname .ne. 'L1B_AMSU') then
  print *, "Error: bad swath name ", swathname, " in file ", file_name
  print *, "Expected L1B_AMSU"
  stop
end if

! Attach to (open) the one swath.
swid = swattach(fid, swathname)
if (swid .eq. -1) then
  print *, "Failed to attach to swath ", swathname, " in file ", file_name
  stop
end if

! Attributes
statn = swrdattr(swid, "processing_level", amsua_bt_gran%processing_level)
statn = swrdattr(swid, "instrument", amsua_bt_gran%instrument)
statn = swrdattr(swid, "DayNightFlag", amsua_bt_gran%DayNightFlag)
statn = swrdattr(swid, "AutomaticQAFlag", amsua_bt_gran%AutomaticQAFlag)
statn = swrdattr(swid, "NumTotalData", amsua_bt_gran%NumTotalData)
statn = swrdattr(swid, "NumProcessData", amsua_bt_gran%NumProcessData)
statn = swrdattr(swid, "NumSpecialData", amsua_bt_gran%NumSpecialData)
statn = swrdattr(swid, "NumBadData", amsua_bt_gran%NumBadData)
statn = swrdattr(swid, "NumMissingData",amsua_bt_gran%NumMissingData)
statn = swrdattr(swid, "NumLandSurface", amsua_bt_gran%NumLandSurface)
statn = swrdattr(swid, "NumOceanSurface", amsua_bt_gran%NumOceanSurface)
statn = swrdattr(swid, "node_type", amsua_bt_gran%node_type)
statn = swrdattr(swid, "start_year", amsua_bt_gran%start_year)
statn = swrdattr(swid, "start_month", amsua_bt_gran%start_month)
statn = swrdattr(swid, "start_day", amsua_bt_gran%start_day)
statn = swrdattr(swid, "start_hour", amsua_bt_gran%start_hour)
statn = swrdattr(swid, "start_minute", amsua_bt_gran%start_minute)
statn = swrdattr(swid, "start_sec", amsua_bt_gran%start_sec)
statn = swrdattr(swid, "start_orbit", amsua_bt_gran%start_orbit)
statn = swrdattr(swid, "end_orbit", amsua_bt_gran%end_orbit)
statn = swrdattr(swid, "orbit_path", amsua_bt_gran%orbit_path)
statn = swrdattr(swid, "start_orbit_row", amsua_bt_gran%start_orbit_row)
statn = swrdattr(swid, "end_orbit_row", amsua_bt_gran%end_orbit_row)
statn = swrdattr(swid, "granule_number", amsua_bt_gran%granule_number)
statn = swrdattr(swid, "num_scansets", amsua_bt_gran%num_scansets)
statn = swrdattr(swid, "num_scanlines", amsua_bt_gran%num_scanlines)
statn = swrdattr(swid, "start_Latitude", amsua_bt_gran%start_Latitude)
statn = swrdattr(swid, "start_Longitude", amsua_bt_gran%start_Longitude)
statn = swrdattr(swid, "start_Time", amsua_bt_gran%start_Time)
statn = swrdattr(swid, "end_Latitude", amsua_bt_gran%end_Latitude)
statn = swrdattr(swid, "end_Longitude", amsua_bt_gran%end_Longitude)
statn = swrdattr(swid, "end_Time", amsua_bt_gran%end_Time)
statn = swrdattr(swid, "eq_x_longitude", amsua_bt_gran%eq_x_longitude)
statn = swrdattr(swid, "eq_x_tai", amsua_bt_gran%eq_x_tai)
statn = swrdattr(swid, "orbitgeoqa", amsua_bt_gran%orbitgeoqa)
statn = swrdattr(swid, "num_satgeoqa", amsua_bt_gran%num_satgeoqa)
statn = swrdattr(swid, "num_glintgeoqa", amsua_bt_gran%num_glintgeoqa)
statn = swrdattr(swid, "num_moongeoqa", amsua_bt_gran%num_moongeoqa)
statn = swrdattr(swid, "num_ftptgeoqa", amsua_bt_gran%num_ftptgeoqa)
statn = swrdattr(swid, "num_zengeoqa", amsua_bt_gran%num_zengeoqa)
statn = swrdattr(swid, "num_demgeoqa", amsua_bt_gran%num_demgeoqa)
statn = swrdattr(swid, "num_fpe", amsua_bt_gran%num_fpe)
statn = swrdattr(swid, "LonGranuleCen", amsua_bt_gran%LonGranuleCen)
statn = swrdattr(swid, "LatGranuleCen", amsua_bt_gran%LatGranuleCen)
statn = swrdattr(swid, "LocTimeGranuleCen", amsua_bt_gran%LocTimeGranuleCen)
statn = swrdattr(swid, "num_scanlines_not_norm_mode_a1", amsua_bt_gran%num_scanlines_not_norm_mode_a1)
statn = swrdattr(swid, "num_scanlines_not_norm_mode_a2", amsua_bt_gran%num_scanlines_not_norm_mode_a2)
statn = swrdattr(swid, "num_missing_scanlines_a1", amsua_bt_gran%num_missing_scanlines_a1)
statn = swrdattr(swid, "num_missing_scanlines_a2", amsua_bt_gran%num_missing_scanlines_a2)
statn = swrdattr(swid, "num_data_gaps_a1", amsua_bt_gran%num_data_gaps_a1)
statn = swrdattr(swid, "num_data_gaps_a2", amsua_bt_gran%num_data_gaps_a2)
statn = swrdattr(swid, "num_instr_mode_changes_a1", amsua_bt_gran%num_instr_mode_changes_a1)
statn = swrdattr(swid, "num_instr_mode_changes_a2", amsua_bt_gran%num_instr_mode_changes_a2)
statn = swrdattr(swid, "num_scanlines_rec_cal_prob_a11", amsua_bt_gran%num_scanlines_rec_cal_prob_a11)
statn = swrdattr(swid, "num_scanlines_rec_cal_prob_a12", amsua_bt_gran%num_scanlines_rec_cal_prob_a12)
statn = swrdattr(swid, "num_scanlines_rec_cal_prob_a2", amsua_bt_gran%num_scanlines_rec_cal_prob_a2)
statn = swrdattr(swid, "num_scanlines_sig_coast_xing", amsua_bt_gran%num_scanlines_sig_coast_xing)
statn = swrdattr(swid, "num_scanlines_sig_sun_glint", amsua_bt_gran%num_scanlines_sig_sun_glint)
statn = swrdattr(swid, "MoonInViewMWCount", amsua_bt_gran%MoonInViewMWCount)
statn = swrdattr(swid, "QA_bb_PRT_a11.min", amsua_bt_gran%QA_bb_PRT_a11%min)
statn = swrdattr(swid, "QA_bb_PRT_a11.max", amsua_bt_gran%QA_bb_PRT_a11%max)
statn = swrdattr(swid, "QA_bb_PRT_a11.mean", amsua_bt_gran%QA_bb_PRT_a11%mean)
statn = swrdattr(swid, "QA_bb_PRT_a11.dev", amsua_bt_gran%QA_bb_PRT_a11%dev)
statn = swrdattr(swid, "QA_bb_PRT_a11.num_in", amsua_bt_gran%QA_bb_PRT_a11%num_in)
statn = swrdattr(swid, "QA_bb_PRT_a11.num_lo", amsua_bt_gran%QA_bb_PRT_a11%num_lo)
statn = swrdattr(swid, "QA_bb_PRT_a11.num_hi", amsua_bt_gran%QA_bb_PRT_a11%num_hi)
statn = swrdattr(swid, "QA_bb_PRT_a11.num_bad", amsua_bt_gran%QA_bb_PRT_a11%num_bad)
statn = swrdattr(swid, "QA_bb_PRT_a11.range_min", amsua_bt_gran%QA_bb_PRT_a11%range_min)
statn = swrdattr(swid, "QA_bb_PRT_a11.range_max", amsua_bt_gran%QA_bb_PRT_a11%range_max)
statn = swrdattr(swid, "QA_bb_PRT_a11.missing", amsua_bt_gran%QA_bb_PRT_a11%missing)
statn = swrdattr(swid, "QA_bb_PRT_a11.max_track", amsua_bt_gran%QA_bb_PRT_a11%max_track)
statn = swrdattr(swid, "QA_bb_PRT_a11.max_xtrack", amsua_bt_gran%QA_bb_PRT_a11%max_xtrack)
statn = swrdattr(swid, "QA_bb_PRT_a11.min_track", amsua_bt_gran%QA_bb_PRT_a11%min_track)
statn = swrdattr(swid, "QA_bb_PRT_a11.min_xtrack", amsua_bt_gran%QA_bb_PRT_a11%min_xtrack)
statn = swrdattr(swid, "QA_bb_PRT_a12.min", amsua_bt_gran%QA_bb_PRT_a12%min)
statn = swrdattr(swid, "QA_bb_PRT_a12.max", amsua_bt_gran%QA_bb_PRT_a12%max)
statn = swrdattr(swid, "QA_bb_PRT_a12.mean", amsua_bt_gran%QA_bb_PRT_a12%mean)
statn = swrdattr(swid, "QA_bb_PRT_a12.dev", amsua_bt_gran%QA_bb_PRT_a12%dev)
statn = swrdattr(swid, "QA_bb_PRT_a12.num_in", amsua_bt_gran%QA_bb_PRT_a12%num_in)
statn = swrdattr(swid, "QA_bb_PRT_a12.num_lo", amsua_bt_gran%QA_bb_PRT_a12%num_lo)
statn = swrdattr(swid, "QA_bb_PRT_a12.num_hi", amsua_bt_gran%QA_bb_PRT_a12%num_hi)
statn = swrdattr(swid, "QA_bb_PRT_a12.num_bad", amsua_bt_gran%QA_bb_PRT_a12%num_bad)
statn = swrdattr(swid, "QA_bb_PRT_a12.range_min", amsua_bt_gran%QA_bb_PRT_a12%range_min)
statn = swrdattr(swid, "QA_bb_PRT_a12.range_max", amsua_bt_gran%QA_bb_PRT_a12%range_max)
statn = swrdattr(swid, "QA_bb_PRT_a12.missing", amsua_bt_gran%QA_bb_PRT_a12%missing)
statn = swrdattr(swid, "QA_bb_PRT_a12.max_track", amsua_bt_gran%QA_bb_PRT_a12%max_track)
statn = swrdattr(swid, "QA_bb_PRT_a12.max_xtrack", amsua_bt_gran%QA_bb_PRT_a12%max_xtrack)
statn = swrdattr(swid, "QA_bb_PRT_a12.min_track", amsua_bt_gran%QA_bb_PRT_a12%min_track)
statn = swrdattr(swid, "QA_bb_PRT_a12.min_xtrack", amsua_bt_gran%QA_bb_PRT_a12%min_xtrack)
statn = swrdattr(swid, "QA_bb_PRT_a2.min", amsua_bt_gran%QA_bb_PRT_a2%min)
statn = swrdattr(swid, "QA_bb_PRT_a2.max", amsua_bt_gran%QA_bb_PRT_a2%max)
statn = swrdattr(swid, "QA_bb_PRT_a2.mean", amsua_bt_gran%QA_bb_PRT_a2%mean)
statn = swrdattr(swid, "QA_bb_PRT_a2.dev", amsua_bt_gran%QA_bb_PRT_a2%dev)
statn = swrdattr(swid, "QA_bb_PRT_a2.num_in", amsua_bt_gran%QA_bb_PRT_a2%num_in)
statn = swrdattr(swid, "QA_bb_PRT_a2.num_lo", amsua_bt_gran%QA_bb_PRT_a2%num_lo)
statn = swrdattr(swid, "QA_bb_PRT_a2.num_hi", amsua_bt_gran%QA_bb_PRT_a2%num_hi)
statn = swrdattr(swid, "QA_bb_PRT_a2.num_bad", amsua_bt_gran%QA_bb_PRT_a2%num_bad)
statn = swrdattr(swid, "QA_bb_PRT_a2.range_min", amsua_bt_gran%QA_bb_PRT_a2%range_min)
statn = swrdattr(swid, "QA_bb_PRT_a2.range_max", amsua_bt_gran%QA_bb_PRT_a2%range_max)
statn = swrdattr(swid, "QA_bb_PRT_a2.missing", amsua_bt_gran%QA_bb_PRT_a2%missing)
statn = swrdattr(swid, "QA_bb_PRT_a2.max_track", amsua_bt_gran%QA_bb_PRT_a2%max_track)
statn = swrdattr(swid, "QA_bb_PRT_a2.max_xtrack", amsua_bt_gran%QA_bb_PRT_a2%max_xtrack)
statn = swrdattr(swid, "QA_bb_PRT_a2.min_track", amsua_bt_gran%QA_bb_PRT_a2%min_track)
statn = swrdattr(swid, "QA_bb_PRT_a2.min_xtrack", amsua_bt_gran%QA_bb_PRT_a2%min_xtrack)
statn = swrdattr(swid, "QA_rec_PRT_a11.min", amsua_bt_gran%QA_rec_PRT_a11%min)
statn = swrdattr(swid, "QA_rec_PRT_a11.max", amsua_bt_gran%QA_rec_PRT_a11%max)
statn = swrdattr(swid, "QA_rec_PRT_a11.mean", amsua_bt_gran%QA_rec_PRT_a11%mean)
statn = swrdattr(swid, "QA_rec_PRT_a11.dev", amsua_bt_gran%QA_rec_PRT_a11%dev)
statn = swrdattr(swid, "QA_rec_PRT_a11.num_in", amsua_bt_gran%QA_rec_PRT_a11%num_in)
statn = swrdattr(swid, "QA_rec_PRT_a11.num_lo", amsua_bt_gran%QA_rec_PRT_a11%num_lo)
statn = swrdattr(swid, "QA_rec_PRT_a11.num_hi", amsua_bt_gran%QA_rec_PRT_a11%num_hi)
statn = swrdattr(swid, "QA_rec_PRT_a11.num_bad", amsua_bt_gran%QA_rec_PRT_a11%num_bad)
statn = swrdattr(swid, "QA_rec_PRT_a11.range_min", amsua_bt_gran%QA_rec_PRT_a11%range_min)
statn = swrdattr(swid, "QA_rec_PRT_a11.range_max", amsua_bt_gran%QA_rec_PRT_a11%range_max)
statn = swrdattr(swid, "QA_rec_PRT_a11.missing", amsua_bt_gran%QA_rec_PRT_a11%missing)
statn = swrdattr(swid, "QA_rec_PRT_a11.max_track", amsua_bt_gran%QA_rec_PRT_a11%max_track)
statn = swrdattr(swid, "QA_rec_PRT_a11.max_xtrack", amsua_bt_gran%QA_rec_PRT_a11%max_xtrack)
statn = swrdattr(swid, "QA_rec_PRT_a11.min_track", amsua_bt_gran%QA_rec_PRT_a11%min_track)
statn = swrdattr(swid, "QA_rec_PRT_a11.min_xtrack", amsua_bt_gran%QA_rec_PRT_a11%min_xtrack)
statn = swrdattr(swid, "QA_rec_PRT_a12.min", amsua_bt_gran%QA_rec_PRT_a12%min)
statn = swrdattr(swid, "QA_rec_PRT_a12.max", amsua_bt_gran%QA_rec_PRT_a12%max)
statn = swrdattr(swid, "QA_rec_PRT_a12.mean", amsua_bt_gran%QA_rec_PRT_a12%mean)
statn = swrdattr(swid, "QA_rec_PRT_a12.dev", amsua_bt_gran%QA_rec_PRT_a12%dev)
statn = swrdattr(swid, "QA_rec_PRT_a12.num_in", amsua_bt_gran%QA_rec_PRT_a12%num_in)
statn = swrdattr(swid, "QA_rec_PRT_a12.num_lo", amsua_bt_gran%QA_rec_PRT_a12%num_lo)
statn = swrdattr(swid, "QA_rec_PRT_a12.num_hi", amsua_bt_gran%QA_rec_PRT_a12%num_hi)
statn = swrdattr(swid, "QA_rec_PRT_a12.num_bad", amsua_bt_gran%QA_rec_PRT_a12%num_bad)
statn = swrdattr(swid, "QA_rec_PRT_a12.range_min", amsua_bt_gran%QA_rec_PRT_a12%range_min)
statn = swrdattr(swid, "QA_rec_PRT_a12.range_max", amsua_bt_gran%QA_rec_PRT_a12%range_max)
statn = swrdattr(swid, "QA_rec_PRT_a12.missing", amsua_bt_gran%QA_rec_PRT_a12%missing)
statn = swrdattr(swid, "QA_rec_PRT_a12.max_track", amsua_bt_gran%QA_rec_PRT_a12%max_track)
statn = swrdattr(swid, "QA_rec_PRT_a12.max_xtrack", amsua_bt_gran%QA_rec_PRT_a12%max_xtrack)
statn = swrdattr(swid, "QA_rec_PRT_a12.min_track", amsua_bt_gran%QA_rec_PRT_a12%min_track)
statn = swrdattr(swid, "QA_rec_PRT_a12.min_xtrack", amsua_bt_gran%QA_rec_PRT_a12%min_xtrack)
statn = swrdattr(swid, "QA_rec_PRT_a2.min", amsua_bt_gran%QA_rec_PRT_a2%min)
statn = swrdattr(swid, "QA_rec_PRT_a2.max", amsua_bt_gran%QA_rec_PRT_a2%max)
statn = swrdattr(swid, "QA_rec_PRT_a2.mean", amsua_bt_gran%QA_rec_PRT_a2%mean)
statn = swrdattr(swid, "QA_rec_PRT_a2.dev", amsua_bt_gran%QA_rec_PRT_a2%dev)
statn = swrdattr(swid, "QA_rec_PRT_a2.num_in", amsua_bt_gran%QA_rec_PRT_a2%num_in)
statn = swrdattr(swid, "QA_rec_PRT_a2.num_lo", amsua_bt_gran%QA_rec_PRT_a2%num_lo)
statn = swrdattr(swid, "QA_rec_PRT_a2.num_hi", amsua_bt_gran%QA_rec_PRT_a2%num_hi)
statn = swrdattr(swid, "QA_rec_PRT_a2.num_bad", amsua_bt_gran%QA_rec_PRT_a2%num_bad)
statn = swrdattr(swid, "QA_rec_PRT_a2.range_min", amsua_bt_gran%QA_rec_PRT_a2%range_min)
statn = swrdattr(swid, "QA_rec_PRT_a2.range_max", amsua_bt_gran%QA_rec_PRT_a2%range_max)
statn = swrdattr(swid, "QA_rec_PRT_a2.missing", amsua_bt_gran%QA_rec_PRT_a2%missing)
statn = swrdattr(swid, "QA_rec_PRT_a2.max_track", amsua_bt_gran%QA_rec_PRT_a2%max_track)
statn = swrdattr(swid, "QA_rec_PRT_a2.max_xtrack", amsua_bt_gran%QA_rec_PRT_a2%max_xtrack)
statn = swrdattr(swid, "QA_rec_PRT_a2.min_track", amsua_bt_gran%QA_rec_PRT_a2%min_track)
statn = swrdattr(swid, "QA_rec_PRT_a2.min_xtrack", amsua_bt_gran%QA_rec_PRT_a2%min_xtrack)
statn = swrdattr(swid, "granules_present", amsua_bt_gran%granules_present)

! Geolocation fields
edge(1) = AMSUA_BT_GEOXTRACK
edge(2) = AMSUA_BT_GEOTRACK
statn = swrdfld(swid, "Latitude",  start, stride, edge, amsua_bt_gran%Latitude)
statn = swrdfld(swid, "Longitude", start, stride, edge, amsua_bt_gran%Longitude)
statn = swrdfld(swid, "Time",      start, stride, edge, amsua_bt_gran%Time)


! Data Fields
edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "scanang", start, stride, edge, amsua_bt_gran%scanang)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "satheight", start, stride, edge, amsua_bt_gran%satheight)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "satroll", start, stride, edge, amsua_bt_gran%satroll)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "satpitch", start, stride, edge, amsua_bt_gran%satpitch)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "satyaw", start, stride, edge, amsua_bt_gran%satyaw)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "satgeoqa", start, stride, edge, amsua_bt_gran%satgeoqa)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "glintgeoqa", start, stride, edge, amsua_bt_gran%glintgeoqa)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "moongeoqa", start, stride, edge, amsua_bt_gran%moongeoqa)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "ftptgeoqa", start, stride, edge, amsua_bt_gran%ftptgeoqa)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "zengeoqa", start, stride, edge, amsua_bt_gran%zengeoqa)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "demgeoqa", start, stride, edge, amsua_bt_gran%demgeoqa)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "nadirTAI", start, stride, edge, amsua_bt_gran%nadirTAI)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "sat_lat", start, stride, edge, amsua_bt_gran%sat_lat)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "sat_lon", start, stride, edge, amsua_bt_gran%sat_lon)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "scan_node_type", start, stride, edge, amsua_bt_gran%scan_node_type)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "satzen", start, stride, edge, amsua_bt_gran%satzen)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "satazi", start, stride, edge, amsua_bt_gran%satazi)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "solzen", start, stride, edge, amsua_bt_gran%solzen)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "solazi", start, stride, edge, amsua_bt_gran%solazi)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "glintlat", start, stride, edge, amsua_bt_gran%glintlat)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "glintlon", start, stride, edge, amsua_bt_gran%glintlon)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "sun_glint_distance", start, stride, edge, amsua_bt_gran%sun_glint_distance)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "topog", start, stride, edge, amsua_bt_gran%topog)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "topog_err", start, stride, edge, amsua_bt_gran%topog_err)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "landFrac", start, stride, edge, amsua_bt_gran%landFrac)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_GEOXTRACK
statn = SWrdfld(swid, "landFrac_err", start, stride, edge, amsua_bt_gran%landFrac_err)

edge(3) = AMSUA_BT_GEOTRACK
edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "antenna_temp", start, stride, edge, amsua_bt_gran%antenna_temp)

edge(3) = AMSUA_BT_GEOTRACK
edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "brightness_temp", start, stride, edge, amsua_bt_gran%brightness_temp)

edge(3) = AMSUA_BT_GEOTRACK
edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "brightness_temp_err", start, stride, edge, amsua_bt_gran%brightness_temp_err)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "center_freq", start, stride, edge, amsua_bt_gran%center_freq)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "IF_offset_1", start, stride, edge, amsua_bt_gran%IF_offset_1)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "IF_offset_2", start, stride, edge, amsua_bt_gran%IF_offset_2)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bandwidth", start, stride, edge, amsua_bt_gran%bandwidth)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "num_calibrated_scanlines", start, stride, edge, amsua_bt_gran%num_calibrated_scanlines)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "num_scanlines_ch_cal_problems", start, stride, edge, amsua_bt_gran%num_scanlines_ch_cal_problems)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.min", start, stride, edge, amsua_bt_gran%bb_signals%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.max", start, stride, edge, amsua_bt_gran%bb_signals%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.mean", start, stride, edge, amsua_bt_gran%bb_signals%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.dev", start, stride, edge, amsua_bt_gran%bb_signals%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.num", start, stride, edge, amsua_bt_gran%bb_signals%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.num_bad", start, stride, edge, amsua_bt_gran%bb_signals%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.max_track", start, stride, edge, amsua_bt_gran%bb_signals%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.max_xtrack", start, stride, edge, amsua_bt_gran%bb_signals%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.min_track", start, stride, edge, amsua_bt_gran%bb_signals%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "bb_signals.min_xtrack", start, stride, edge, amsua_bt_gran%bb_signals%min_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.min", start, stride, edge, amsua_bt_gran%space_signals%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.max", start, stride, edge, amsua_bt_gran%space_signals%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.mean", start, stride, edge, amsua_bt_gran%space_signals%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.dev", start, stride, edge, amsua_bt_gran%space_signals%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.num", start, stride, edge, amsua_bt_gran%space_signals%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.num_bad", start, stride, edge, amsua_bt_gran%space_signals%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.max_track", start, stride, edge, amsua_bt_gran%space_signals%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.max_xtrack", start, stride, edge, amsua_bt_gran%space_signals%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.min_track", start, stride, edge, amsua_bt_gran%space_signals%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "space_signals.min_xtrack", start, stride, edge, amsua_bt_gran%space_signals%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.min", start, stride, edge, amsua_bt_gran%gain_stats%min)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.max", start, stride, edge, amsua_bt_gran%gain_stats%max)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.mean", start, stride, edge, amsua_bt_gran%gain_stats%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.dev", start, stride, edge, amsua_bt_gran%gain_stats%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.num", start, stride, edge, amsua_bt_gran%gain_stats%num)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.num_bad", start, stride, edge, amsua_bt_gran%gain_stats%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.max_track", start, stride, edge, amsua_bt_gran%gain_stats%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.max_xtrack", start, stride, edge, amsua_bt_gran%gain_stats%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.min_track", start, stride, edge, amsua_bt_gran%gain_stats%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "gain_stats.min_xtrack", start, stride, edge, amsua_bt_gran%gain_stats%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "NeDT", start, stride, edge, amsua_bt_gran%NeDT)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.min", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%min)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.max", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%max)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.mean", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%mean)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.dev", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%dev)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.num", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%num)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.num_bad", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%num_bad)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.max_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%max_track)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.max_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%max_xtrack)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.min_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%min_track)

edge(2) = AMSUA_BT_GEOXTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_scene_count.min_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_scene_count%min_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.min", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.max", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.mean", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.dev", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.num", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.num_bad", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.max_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.max_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.min_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_BB_count.min_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_BB_count%min_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.min", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%min)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.max", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%max)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.mean", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%mean)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.dev", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%dev)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.num", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%num)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.num_bad", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%num_bad)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.max_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%max_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.max_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%max_xtrack)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.min_track", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%min_track)

edge(2) = 2
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_unfiltered_space_count.min_xtrack", start, stride, edge, amsua_bt_gran%QA_unfiltered_space_count%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.min", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%min)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.max", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%max)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.mean", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.dev", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.num", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%num)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.num_bad", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.max_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.max_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.min_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a0.min_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a0%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.min", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%min)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.max", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%max)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.mean", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.dev", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.num", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%num)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.num_bad", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.max_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.max_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.min_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a1.min_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a1%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.min", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%min)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.max", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%max)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.mean", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.dev", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.num", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%num)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.num_bad", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.max_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.max_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.min_track", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_cal_coef_a2.min_xtrack", start, stride, edge, amsua_bt_gran%QA_cal_coef_a2%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.min", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%min)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.max", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%max)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.mean", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.dev", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.num", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%num)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.num_bad", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.max_track", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.max_xtrack", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.min_track", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_bb_raw_noise_counts.min_xtrack", start, stride, edge, amsua_bt_gran%QA_bb_raw_noise_counts%min_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.min", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%min)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.max", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%max)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.mean", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%mean)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.dev", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%dev)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.num", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%num)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.num_bad", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%num_bad)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.max_track", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%max_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.max_xtrack", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%max_xtrack)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.min_track", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%min_track)

edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "QA_sv_raw_noise_counts.min_xtrack", start, stride, edge, amsua_bt_gran%QA_sv_raw_noise_counts%min_xtrack)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "state1", start, stride, edge, amsua_bt_gran%state1)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "state2", start, stride, edge, amsua_bt_gran%state2)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "cal_coef_a0", start, stride, edge, amsua_bt_gran%cal_coef_a0)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "cal_coef_a0_err", start, stride, edge, amsua_bt_gran%cal_coef_a0_err)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "cal_coef_a1", start, stride, edge, amsua_bt_gran%cal_coef_a1)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "cal_coef_a1_err", start, stride, edge, amsua_bt_gran%cal_coef_a1_err)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "cal_coef_a2", start, stride, edge, amsua_bt_gran%cal_coef_a2)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "cal_coef_a2_err", start, stride, edge, amsua_bt_gran%cal_coef_a2_err)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a1_ColdCalPstion", start, stride, edge, amsua_bt_gran%a1_ColdCalPstion)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a2_ColdCalPstion", start, stride, edge, amsua_bt_gran%a2_ColdCalPstion)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a1_PLO_Redundncy", start, stride, edge, amsua_bt_gran%a1_PLO_Redundncy)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a11_mux_temp_used", start, stride, edge, amsua_bt_gran%a11_mux_temp_used)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a11_receiver_temp", start, stride, edge, amsua_bt_gran%a11_receiver_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a11_target_temp", start, stride, edge, amsua_bt_gran%a11_target_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a12_mux_temp_used", start, stride, edge, amsua_bt_gran%a12_mux_temp_used)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a12_receiver_temp", start, stride, edge, amsua_bt_gran%a12_receiver_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a12_target_temp", start, stride, edge, amsua_bt_gran%a12_target_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a2_diplexer_temp_used", start, stride, edge, amsua_bt_gran%a2_diplexer_temp_used)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a2_receiver_temp", start, stride, edge, amsua_bt_gran%a2_receiver_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "a2_target_temp", start, stride, edge, amsua_bt_gran%a2_target_temp)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "qa_scanline", start, stride, edge, amsua_bt_gran%qa_scanline)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "qa_receiver_a11", start, stride, edge, amsua_bt_gran%qa_receiver_a11)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "qa_receiver_a12", start, stride, edge, amsua_bt_gran%qa_receiver_a12)

edge(1) = AMSUA_BT_GEOTRACK
statn = SWrdfld(swid, "qa_receiver_a2", start, stride, edge, amsua_bt_gran%qa_receiver_a2)

edge(2) = AMSUA_BT_GEOTRACK
edge(1) = AMSUA_BT_CHANNEL
statn = SWrdfld(swid, "qa_channel", start, stride, edge, amsua_bt_gran%qa_channel)

! Final clean-up
statn = swdetach(swid)
statn = swclose(fid)

return
end subroutine amsua_bt_rdr



end module amsua_bt_mod
