! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download 

module amsua_netCDF_support_mod

use         types_mod, only : r4, r8, digits12, deg2rad, rad2deg, PI, MISSING_R8

use     utilities_mod, only : error_handler, E_MSG, E_ERR, &
                              is_longitude_between

use  time_manager_mod, only : time_type, get_date, set_date,            &
                              get_time, set_time, set_calendar_type,    &
                              GREGORIAN, print_date, print_time,        &
                              operator(+)

use  obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                              set_obs_values, obs_sequence_type, get_num_obs, &
                              obs_type, set_copy_meta_data, set_qc_meta_data, &
                              print_obs_seq_summary, get_first_obs, get_next_obs, &
                              read_obs_seq

use      location_mod, only : location_type, VERTISUNDEF, &
                              set_location, get_location

use      obs_kind_mod, only : get_index_for_type_of_obs

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use netcdf_utilities_mod, only : nc_open_file_readonly, &
                                 nc_get_dimension_size, &
                                 nc_get_global_attribute, &
                                 nc_get_variable, &
                                 nc_define_dimension, &
                                 nc_define_unlimited_dimension, &
                                 nc_define_double_variable, &
                                 nc_define_real_variable, &
                                 nc_define_integer_variable, &
                                 nc_add_attribute_to_variable, &
                                 nc_put_variable, &
                                 nc_end_define_mode, &
                                 NF90_MAX_NAME, NF90_MAX_VAR_DIMS, &
                                 nc_close_file

use obs_def_rttov_mod, only : set_mw_metadata, &
                              get_rttov_option_logical

use amsua_bt_mod, only : amsua_bt_granule, &
                         AMSUA_BT_GEOXTRACK,  AMSUA_BT_GEOTRACK,    AMSUA_BT_CHANNEL, &
                         AMSUA_BT_CALXTRACK,  AMSUA_BT_SPACEXTRACK, AMSUA_BT_BBXTRACK, &
                         AMSUA_BT_WARMPRTA11, AMSUA_BT_WARMPRTA12,  AMSUA_BT_WARMPRTA2

implicit none
private

public :: initialize_amsua_netcdf, &
          channel_list_to_indices, &
          amsua_bt_granule, &
          read_amsua_bt_netCDF_granule, &
          make_obs_sequence, &
          define_amsua_dimensions, &
          define_amsua_variables, &
          fill_amsua_variables, &
          max_possible_obs, &
          append_or_create, &
          combine_sequences

character(len=*), parameter :: source = 'amsua_netCDF_support_mod.f90'

logical, save :: module_initialized = .false.

character(len=512) :: string1, string2, string3

! flag to specify how much run-time output to create.
! 0 == none, higher is more.
integer :: verbosity = 0

contains


!------------------------------------------------------------------------------
!>

subroutine initialize_module

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module


!------------------------------------------------------------------------------
!>

subroutine initialize_amsua_netcdf(filename, verbose)

character(len=*), intent(in) :: filename
integer,          intent(in) :: verbose

integer :: ncid
integer :: GEOXTRACK,  GEOTRACK,    CHANNEL, &
           CALXTRACK,  SPACEXTRACK, BBXTRACK, &
           WARMPRTA11, WARMPRTA12,  WARMPRTA2
logical :: mismatch

character(len=*), parameter :: context = 'initialize_amsua_netcdf'

if ( .not. module_initialized ) call initialize_module

! set verbosity for entire module
verbosity = verbose

ncid = nc_open_file_readonly(filename,context)

GEOXTRACK   = nc_get_dimension_size(ncid, 'GeoXTrack',   context)
GEOTRACK    = nc_get_dimension_size(ncid, 'GeoTrack',    context)
CHANNEL     = nc_get_dimension_size(ncid, 'Channel',     context)
CALXTRACK   = nc_get_dimension_size(ncid, 'CalXTrack',   context)
SPACEXTRACK = nc_get_dimension_size(ncid, 'SpaceXTrack', context)
BBXTRACK    = nc_get_dimension_size(ncid, 'BBXTrack',    context)
WARMPRTA11  = nc_get_dimension_size(ncid, 'WarmPRTA11',  context)
WARMPRTA12  = nc_get_dimension_size(ncid, 'WarmPRTA12',  context)
WARMPRTA2   = nc_get_dimension_size(ncid, 'WarmPRTA2',   context)

call nc_close_file(ncid,context)

! compare to values from amsua_bt_mod.f90 just to make sure
mismatch = .false.
if ( AMSUA_BT_GEOXTRACK   /= GEOXTRACK   ) mismatch = .true.
if ( AMSUA_BT_GEOTRACK    /= GEOTRACK    ) mismatch = .true.
if ( AMSUA_BT_CHANNEL     /= CHANNEL     ) mismatch = .true.
if ( AMSUA_BT_CALXTRACK   /= CALXTRACK   ) mismatch = .true.
if ( AMSUA_BT_SPACEXTRACK /= SPACEXTRACK ) mismatch = .true.
if ( AMSUA_BT_BBXTRACK    /= BBXTRACK    ) mismatch = .true.
if ( AMSUA_BT_WARMPRTA11  /= WARMPRTA11  ) mismatch = .true.
if ( AMSUA_BT_WARMPRTA12  /= WARMPRTA12  ) mismatch = .true.
if ( AMSUA_BT_WARMPRTA2   /= WARMPRTA2   ) mismatch = .true.

if (mismatch) then
   call error_handler(E_ERR,'initialize_amsua_netcdf','assumptions wrong',source,&
              text2='dimensions in file do not match expected dimensions')
endif

end subroutine initialize_amsua_netcdf


!------------------------------------------------------------------------------
!> Convert character representation of desired channels to the channel index
!> The desired input can be either an alphanumeric or numeric representation.
!>
!> from the documentation
!> A2-1, A2-2 | A1-1, A1-2, A1-3,   A1-4, A1-5    A1-6, A1-7 ...  A1-13
!> 23.8, 31.4 | 50.3, 52.8, 53.596, 54.4, 54.940, 55.5, 57.29034, 89

subroutine channel_list_to_indices(channel_list,use_channels)

character(len=*), intent(in)  :: channel_list(:)
logical,          intent(out) :: use_channels(:)

! from the dump of the converted netCDF file (chs 9-14 ?)
! center_freq = 23.8,     31.4,     50.3,     52.8,     53.596, ...
!               54.4,     54.94,    55.5,     57.29034, 57.29034, ...
!               57.29034, 57.29034, 57.29034, 57.29034, 89 ;

integer :: ichannel

if ( .not. module_initialized ) call initialize_module

use_channels = .false.

CHANNELS: do ichannel = 1,size(channel_list)
   select case (channel_list(ichannel))
      case('1','A2-1','23.8')
         use_channels(1) = .true.
      case('2','A2-2','31.4')
         use_channels(2) = .true.
      case('3','A1-1','50.3')
         use_channels(3) = .true.
      case('4','A1-2','52.8')
         use_channels(4) = .true.
      case('5','A1-3','53.596')
         use_channels(5) = .true.
      case('6','A1-4','54.4')
         use_channels(6) = .true.
      case('7','A1-5','54.94')
         use_channels(7) = .true.
      case('8','A1-6','55.5')
         use_channels(8) = .true.
      case('9','A1-7')
         use_channels(9) = .true.
      case('10','A1-8')
         use_channels(10) = .true.
      case('11','A1-9')
         use_channels(11) = .true.
      case('12','A1-10')
         use_channels(12) = .true.
      case('13','A1-11')
         use_channels(13) = .true.
      case('14','A1-12')
         use_channels(14) = .true.
      case('15','A1-13','89')
         use_channels(15) = .true.
      case ('null')
         exit CHANNELS
      case default
         write(string1,*)'<'//trim(channel_list(ichannel))//'> is not a valid channel declaration'
         call error_handler(E_ERR,'channel_list_to_indices','unknown channel input', &
                    source,text2=trim(channel_list(ichannel)))
    end select
enddo CHANNELS

if (all(use_channels .eqv. .false.)) then
   call error_handler(E_ERR,'channel_list_to_indices','no valid input channels',source) 
endif

! FIXME ... print a summary of the channels and frequences being used

end subroutine channel_list_to_indices


!------------------------------------------------------------------------------
!> Routine to ingest an entire netCDF-format granule as converted by h4tonccf_nc4.
!>
!> Based on the 'amsua_bt_rdr' routine which was part of the NASA-provided software.
!> Modified by TJH Sep 2020: 
!> note: 'swath' seems historical, 'granule' is in the current documentation.

subroutine read_amsua_bt_netCDF_granule(file_name, granule)

character(len=*),       intent(in)  :: file_name
TYPE(amsua_bt_granule), intent(out) :: granule

character(len=*), parameter :: context = 'read_amsua_bt_netCDF_granule'
character(len=*), parameter :: SWATHNAME = 'L1B_AMSU'

integer :: swid                   ! netCDF file identifier
integer :: nchar
character(len=256) :: swath_title ! Name of swath

if ( .not. module_initialized ) call initialize_module

swid = nc_open_file_readonly(file_name, context)

! Get name of swath, make sure it is the right one
! turns out there are non-printable characters in what you get from the netCDF file

call nc_get_global_attribute(swid, 'SWATH_TITLE__', swath_title, context)

nchar = len_trim(SWATHNAME)

if (swath_title(1:nchar) .ne. SWATHNAME) then
   write(string1,*) 'ERROR: bad SWATH_TITLE__ global attribute in file "'//trim(file_name)//'"'
   write(string2,*) 'Expected "'//SWATHNAME//'"'
   write(string3,*) 'got      "'//swath_title(1:nchar)//'"'
   call error_handler(E_ERR, context, string1, &
              source, text2=string2, text3=string3)
end if

! FIXME ... unsupported types are commented out for now.
! Support for them must be provided in netcdf_utilities_mod.f90

! Attributes
call nc_get_global_attribute(swid, "processing_level", granule%processing_level, context)
call nc_get_global_attribute(swid, "instrument", granule%instrument, context)
call nc_get_global_attribute(swid, "DayNightFlag", granule%DayNightFlag, context)
call nc_get_global_attribute(swid, "AutomaticQAFlag", granule%AutomaticQAFlag, context)
call nc_get_global_attribute(swid, "NumTotalData", granule%NumTotalData, context)
call nc_get_global_attribute(swid, "NumProcessData", granule%NumProcessData, context)
call nc_get_global_attribute(swid, "NumSpecialData", granule%NumSpecialData, context)
call nc_get_global_attribute(swid, "NumBadData", granule%NumBadData, context)
call nc_get_global_attribute(swid, "NumMissingData",granule%NumMissingData, context)
call nc_get_global_attribute(swid, "NumLandSurface", granule%NumLandSurface, context)
call nc_get_global_attribute(swid, "NumOceanSurface", granule%NumOceanSurface, context)
call nc_get_global_attribute(swid, "node_type", granule%node_type, context)
call nc_get_global_attribute(swid, "start_year", granule%start_year, context)
call nc_get_global_attribute(swid, "start_month", granule%start_month, context)
call nc_get_global_attribute(swid, "start_day", granule%start_day, context)
call nc_get_global_attribute(swid, "start_hour", granule%start_hour, context)
call nc_get_global_attribute(swid, "start_minute", granule%start_minute, context)
call nc_get_global_attribute(swid, "start_sec", granule%start_sec, context)
call nc_get_global_attribute(swid, "start_orbit", granule%start_orbit, context)
call nc_get_global_attribute(swid, "end_orbit", granule%end_orbit, context)
call nc_get_global_attribute(swid, "orbit_path", granule%orbit_path, context)
call nc_get_global_attribute(swid, "start_orbit_row", granule%start_orbit_row, context)
call nc_get_global_attribute(swid, "end_orbit_row", granule%end_orbit_row, context)
call nc_get_global_attribute(swid, "granule_number", granule%granule_number, context)
call nc_get_global_attribute(swid, "num_scansets", granule%num_scansets, context)
call nc_get_global_attribute(swid, "num_scanlines", granule%num_scanlines, context)
call nc_get_global_attribute(swid, "start_Latitude", granule%start_Latitude, context)
call nc_get_global_attribute(swid, "start_Longitude", granule%start_Longitude, context)
call nc_get_global_attribute(swid, "start_Time", granule%start_Time, context)
call nc_get_global_attribute(swid, "end_Latitude", granule%end_Latitude, context)
call nc_get_global_attribute(swid, "end_Longitude", granule%end_Longitude, context)
call nc_get_global_attribute(swid, "end_Time", granule%end_Time, context)
call nc_get_global_attribute(swid, "eq_x_longitude", granule%eq_x_longitude, context)
call nc_get_global_attribute(swid, "eq_x_tai", granule%eq_x_tai, context)
call nc_get_global_attribute(swid, "orbitgeoqa", granule%orbitgeoqa, context)
! call nc_get_global_attribute(swid, "num_satgeoqa", granule%num_satgeoqa, context)
! call nc_get_global_attribute(swid, "num_glintgeoqa", granule%num_glintgeoqa, context)
! call nc_get_global_attribute(swid, "num_moongeoqa", granule%num_moongeoqa, context)
! call nc_get_global_attribute(swid, "num_ftptgeoqa", granule%num_ftptgeoqa, context)
! call nc_get_global_attribute(swid, "num_zengeoqa", granule%num_zengeoqa, context)
! call nc_get_global_attribute(swid, "num_demgeoqa", granule%num_demgeoqa, context)
! call nc_get_global_attribute(swid, "num_fpe", granule%num_fpe, context)
! call nc_get_global_attribute(swid, "LonGranuleCen", granule%LonGranuleCen, context)
! call nc_get_global_attribute(swid, "LatGranuleCen", granule%LatGranuleCen, context)
! call nc_get_global_attribute(swid, "LocTimeGranuleCen", granule%LocTimeGranuleCen, context)
call nc_get_global_attribute(swid, "num_scanlines_not_norm_mode_a1", granule%num_scanlines_not_norm_mode_a1, context)
call nc_get_global_attribute(swid, "num_scanlines_not_norm_mode_a2", granule%num_scanlines_not_norm_mode_a2, context)
call nc_get_global_attribute(swid, "num_missing_scanlines_a1", granule%num_missing_scanlines_a1, context)
call nc_get_global_attribute(swid, "num_missing_scanlines_a2", granule%num_missing_scanlines_a2, context)
call nc_get_global_attribute(swid, "num_data_gaps_a1", granule%num_data_gaps_a1, context)
call nc_get_global_attribute(swid, "num_data_gaps_a2", granule%num_data_gaps_a2, context)
call nc_get_global_attribute(swid, "num_instr_mode_changes_a1", granule%num_instr_mode_changes_a1, context)
call nc_get_global_attribute(swid, "num_instr_mode_changes_a2", granule%num_instr_mode_changes_a2, context)
call nc_get_global_attribute(swid, "num_scanlines_rec_cal_prob_a11", granule%num_scanlines_rec_cal_prob_a11, context)
call nc_get_global_attribute(swid, "num_scanlines_rec_cal_prob_a12", granule%num_scanlines_rec_cal_prob_a12, context)
call nc_get_global_attribute(swid, "num_scanlines_rec_cal_prob_a2", granule%num_scanlines_rec_cal_prob_a2, context)
call nc_get_global_attribute(swid, "num_scanlines_sig_coast_xing", granule%num_scanlines_sig_coast_xing, context)
call nc_get_global_attribute(swid, "num_scanlines_sig_sun_glint", granule%num_scanlines_sig_sun_glint, context)
call nc_get_global_attribute(swid, "MoonInViewMWCount", granule%MoonInViewMWCount, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_min", granule%QA_bb_PRT_a11%min, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_max", granule%QA_bb_PRT_a11%max, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_mean", granule%QA_bb_PRT_a11%mean, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_dev", granule%QA_bb_PRT_a11%dev, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_num_in", granule%QA_bb_PRT_a11%num_in, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_num_lo", granule%QA_bb_PRT_a11%num_lo, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_num_hi", granule%QA_bb_PRT_a11%num_hi, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_num_bad", granule%QA_bb_PRT_a11%num_bad, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_range_min", granule%QA_bb_PRT_a11%range_min, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_range_max", granule%QA_bb_PRT_a11%range_max, context)
! call nc_get_global_attribute(swid, "QA_bb_PRT_a11_missing", granule%QA_bb_PRT_a11%missing, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_max_track", granule%QA_bb_PRT_a11%max_track, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_max_xtrack", granule%QA_bb_PRT_a11%max_xtrack, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_min_track", granule%QA_bb_PRT_a11%min_track, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a11_min_xtrack", granule%QA_bb_PRT_a11%min_xtrack, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_min", granule%QA_bb_PRT_a12%min, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_max", granule%QA_bb_PRT_a12%max, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_mean", granule%QA_bb_PRT_a12%mean, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_dev", granule%QA_bb_PRT_a12%dev, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_num_in", granule%QA_bb_PRT_a12%num_in, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_num_lo", granule%QA_bb_PRT_a12%num_lo, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_num_hi", granule%QA_bb_PRT_a12%num_hi, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_num_bad", granule%QA_bb_PRT_a12%num_bad, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_range_min", granule%QA_bb_PRT_a12%range_min, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_range_max", granule%QA_bb_PRT_a12%range_max, context)
! call nc_get_global_attribute(swid, "QA_bb_PRT_a12_missing", granule%QA_bb_PRT_a12%missing, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_max_track", granule%QA_bb_PRT_a12%max_track, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_max_xtrack", granule%QA_bb_PRT_a12%max_xtrack, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_min_track", granule%QA_bb_PRT_a12%min_track, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a12_min_xtrack", granule%QA_bb_PRT_a12%min_xtrack, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_min", granule%QA_bb_PRT_a2%min, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_max", granule%QA_bb_PRT_a2%max, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_mean", granule%QA_bb_PRT_a2%mean, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_dev", granule%QA_bb_PRT_a2%dev, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_num_in", granule%QA_bb_PRT_a2%num_in, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_num_lo", granule%QA_bb_PRT_a2%num_lo, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_num_hi", granule%QA_bb_PRT_a2%num_hi, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_num_bad", granule%QA_bb_PRT_a2%num_bad, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_range_min", granule%QA_bb_PRT_a2%range_min, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_range_max", granule%QA_bb_PRT_a2%range_max, context)
! call nc_get_global_attribute(swid, "QA_bb_PRT_a2_missing", granule%QA_bb_PRT_a2%missing, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_max_track", granule%QA_bb_PRT_a2%max_track, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_max_xtrack", granule%QA_bb_PRT_a2%max_xtrack, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_min_track", granule%QA_bb_PRT_a2%min_track, context)
call nc_get_global_attribute(swid, "QA_bb_PRT_a2_min_xtrack", granule%QA_bb_PRT_a2%min_xtrack, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_min", granule%QA_rec_PRT_a11%min, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_max", granule%QA_rec_PRT_a11%max, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_mean", granule%QA_rec_PRT_a11%mean, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_dev", granule%QA_rec_PRT_a11%dev, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_num_in", granule%QA_rec_PRT_a11%num_in, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_num_lo", granule%QA_rec_PRT_a11%num_lo, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_num_hi", granule%QA_rec_PRT_a11%num_hi, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_num_bad", granule%QA_rec_PRT_a11%num_bad, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_range_min", granule%QA_rec_PRT_a11%range_min, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_range_max", granule%QA_rec_PRT_a11%range_max, context)
! call nc_get_global_attribute(swid, "QA_rec_PRT_a11.missing", granule%QA_rec_PRT_a11%missing, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_max_track", granule%QA_rec_PRT_a11%max_track, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_max_xtrack", granule%QA_rec_PRT_a11%max_xtrack, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_min_track", granule%QA_rec_PRT_a11%min_track, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a11_min_xtrack", granule%QA_rec_PRT_a11%min_xtrack, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_min", granule%QA_rec_PRT_a12%min, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_max", granule%QA_rec_PRT_a12%max, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_mean", granule%QA_rec_PRT_a12%mean, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_dev", granule%QA_rec_PRT_a12%dev, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_num_in", granule%QA_rec_PRT_a12%num_in, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_num_lo", granule%QA_rec_PRT_a12%num_lo, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_num_hi", granule%QA_rec_PRT_a12%num_hi, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_num_bad", granule%QA_rec_PRT_a12%num_bad, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_range_min", granule%QA_rec_PRT_a12%range_min, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_range_max", granule%QA_rec_PRT_a12%range_max, context)
! call nc_get_global_attribute(swid, "QA_rec_PRT_a12_missing", granule%QA_rec_PRT_a12%missing, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_max_track", granule%QA_rec_PRT_a12%max_track, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_max_xtrack", granule%QA_rec_PRT_a12%max_xtrack, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_min_track", granule%QA_rec_PRT_a12%min_track, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a12_min_xtrack", granule%QA_rec_PRT_a12%min_xtrack, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_min", granule%QA_rec_PRT_a2%min, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_max", granule%QA_rec_PRT_a2%max, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_mean", granule%QA_rec_PRT_a2%mean, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_dev", granule%QA_rec_PRT_a2%dev, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_num_in", granule%QA_rec_PRT_a2%num_in, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_num_lo", granule%QA_rec_PRT_a2%num_lo, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_num_hi", granule%QA_rec_PRT_a2%num_hi, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_num_bad", granule%QA_rec_PRT_a2%num_bad, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_range_min", granule%QA_rec_PRT_a2%range_min, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_range_max", granule%QA_rec_PRT_a2%range_max, context)
! call nc_get_global_attribute(swid, "QA_rec_PRT_a2_missing", granule%QA_rec_PRT_a2%missing, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_max_track", granule%QA_rec_PRT_a2%max_track, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_max_xtrack", granule%QA_rec_PRT_a2%max_xtrack, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_min_track", granule%QA_rec_PRT_a2%min_track, context)
call nc_get_global_attribute(swid, "QA_rec_PRT_a2_min_xtrack", granule%QA_rec_PRT_a2%min_xtrack, context)
call nc_get_global_attribute(swid, "granules_present", granule%granules_present, context)

if (verbosity > 2) call dump_attributes(granule,file_name)

! Geolocation fields
call nc_get_variable(swid, "latitude",  granule%latitude, context)
call nc_get_variable(swid, "longitude", granule%longitude, context)
call nc_get_variable(swid, "Time",      granule%Time, context)

! Data Fields
call nc_get_variable(swid, "scanang",   granule%scanang, context)
call nc_get_variable(swid, "satheight", granule%satheight, context)
call nc_get_variable(swid, "satroll",   granule%satroll, context)
call nc_get_variable(swid, "satpitch",  granule%satpitch, context)
call nc_get_variable(swid, "satyaw",    granule%satyaw, context)

call nc_get_variable(swid, "satgeoqa",   granule%satgeoqa, context)
call nc_get_variable(swid, "glintgeoqa", granule%glintgeoqa, context)
call nc_get_variable(swid, "moongeoqa",  granule%moongeoqa, context)
call nc_get_variable(swid, "ftptgeoqa",  granule%ftptgeoqa, context)
call nc_get_variable(swid, "zengeoqa",   granule%zengeoqa, context)
call nc_get_variable(swid, "demgeoqa",   granule%demgeoqa, context)

call nc_get_variable(swid, "nadirTAI", granule%nadirTAI, context)
call nc_get_variable(swid, "sat_lat", granule%sat_lat, context)
call nc_get_variable(swid, "sat_lon", granule%sat_lon, context)
! call nc_get_variable(swid, "scan_node_type", granule%scan_node_type, context)

call nc_get_variable(swid, "satzen", granule%satzen, context)
call nc_get_variable(swid, "satazi", granule%satazi, context)
call nc_get_variable(swid, "solzen", granule%solzen, context)
call nc_get_variable(swid, "solazi", granule%solazi, context)
call nc_get_variable(swid, "glintlat", granule%glintlat, context)
call nc_get_variable(swid, "glintlon", granule%glintlon, context)
call nc_get_variable(swid, "sun_glint_distance", granule%sun_glint_distance, context)
call nc_get_variable(swid, "topog", granule%topog, context)
call nc_get_variable(swid, "topog_err", granule%topog_err, context)
call nc_get_variable(swid, "landFrac", granule%landFrac, context)
call nc_get_variable(swid, "landFrac_err", granule%landFrac_err, context)
call nc_get_variable(swid, "antenna_temp", granule%antenna_temp, context)
call nc_get_variable(swid, "brightness_temp", granule%brightness_temp, context)
call nc_get_variable(swid, "brightness_temp_err", granule%brightness_temp_err, context)
call nc_get_variable(swid, "center_freq", granule%center_freq, context)
call nc_get_variable(swid, "IF_offset_1", granule%IF_offset_1, context)
call nc_get_variable(swid, "IF_offset_2", granule%IF_offset_2, context)
call nc_get_variable(swid, "bandwidth", granule%bandwidth, context)
call nc_get_variable(swid, "num_calibrated_scanlines", granule%num_calibrated_scanlines, context)
call nc_get_variable(swid, "num_scanlines_ch_cal_problems", granule%num_scanlines_ch_cal_problems, context)
call nc_get_variable(swid, "bb_signals_min", granule%bb_signals%min, context)
call nc_get_variable(swid, "bb_signals_max", granule%bb_signals%max, context)
call nc_get_variable(swid, "bb_signals_mean", granule%bb_signals%mean, context)
call nc_get_variable(swid, "bb_signals_dev", granule%bb_signals%dev, context)
call nc_get_variable(swid, "bb_signals_num", granule%bb_signals%num, context)
call nc_get_variable(swid, "bb_signals_num_bad", granule%bb_signals%num_bad, context)
call nc_get_variable(swid, "bb_signals_max_track", granule%bb_signals%max_track, context)
call nc_get_variable(swid, "bb_signals_max_xtrack", granule%bb_signals%max_xtrack, context)
call nc_get_variable(swid, "bb_signals_min_track", granule%bb_signals%min_track, context)
call nc_get_variable(swid, "bb_signals_min_xtrack", granule%bb_signals%min_xtrack, context)
call nc_get_variable(swid, "space_signals_min", granule%space_signals%min, context)
call nc_get_variable(swid, "space_signals_max", granule%space_signals%max, context)
call nc_get_variable(swid, "space_signals_mean", granule%space_signals%mean, context)
call nc_get_variable(swid, "space_signals_dev", granule%space_signals%dev, context)
call nc_get_variable(swid, "space_signals_num", granule%space_signals%num, context)
call nc_get_variable(swid, "space_signals_num_bad", granule%space_signals%num_bad, context)
call nc_get_variable(swid, "space_signals_max_track", granule%space_signals%max_track, context)
call nc_get_variable(swid, "space_signals_max_xtrack", granule%space_signals%max_xtrack, context)
call nc_get_variable(swid, "space_signals_min_track", granule%space_signals%min_track, context)
call nc_get_variable(swid, "space_signals_min_xtrack", granule%space_signals%min_xtrack, context)
call nc_get_variable(swid, "gain_stats_min", granule%gain_stats%min, context)
call nc_get_variable(swid, "gain_stats_max", granule%gain_stats%max, context)
call nc_get_variable(swid, "gain_stats_mean", granule%gain_stats%mean, context)
call nc_get_variable(swid, "gain_stats_dev", granule%gain_stats%dev, context)
call nc_get_variable(swid, "gain_stats_num", granule%gain_stats%num, context)
call nc_get_variable(swid, "gain_stats_num_bad", granule%gain_stats%num_bad, context)
call nc_get_variable(swid, "gain_stats_max_track", granule%gain_stats%max_track, context)
call nc_get_variable(swid, "gain_stats_max_xtrack", granule%gain_stats%max_xtrack, context)
call nc_get_variable(swid, "gain_stats_min_track", granule%gain_stats%min_track, context)
call nc_get_variable(swid, "gain_stats_min_xtrack", granule%gain_stats%min_xtrack, context)
call nc_get_variable(swid, "NeDT", granule%NeDT, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_min", granule%QA_unfiltered_scene_count%min, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_max", granule%QA_unfiltered_scene_count%max, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_mean", granule%QA_unfiltered_scene_count%mean, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_dev", granule%QA_unfiltered_scene_count%dev, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_num", granule%QA_unfiltered_scene_count%num, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_num_bad", granule%QA_unfiltered_scene_count%num_bad, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_max_track", granule%QA_unfiltered_scene_count%max_track, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_max_xtrack", granule%QA_unfiltered_scene_count%max_xtrack, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_min_track", granule%QA_unfiltered_scene_count%min_track, context)
call nc_get_variable(swid, "QA_unfiltered_scene_count_min_xtrack", granule%QA_unfiltered_scene_count%min_xtrack, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_min", granule%QA_unfiltered_BB_count%min, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_max", granule%QA_unfiltered_BB_count%max, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_mean", granule%QA_unfiltered_BB_count%mean, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_dev", granule%QA_unfiltered_BB_count%dev, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_num", granule%QA_unfiltered_BB_count%num, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_num_bad", granule%QA_unfiltered_BB_count%num_bad, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_max_track", granule%QA_unfiltered_BB_count%max_track, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_max_xtrack", granule%QA_unfiltered_BB_count%max_xtrack, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_min_track", granule%QA_unfiltered_BB_count%min_track, context)
call nc_get_variable(swid, "QA_unfiltered_BB_count_min_xtrack", granule%QA_unfiltered_BB_count%min_xtrack, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_min", granule%QA_unfiltered_space_count%min, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_max", granule%QA_unfiltered_space_count%max, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_mean", granule%QA_unfiltered_space_count%mean, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_dev", granule%QA_unfiltered_space_count%dev, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_num", granule%QA_unfiltered_space_count%num, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_num_bad", granule%QA_unfiltered_space_count%num_bad, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_max_track", granule%QA_unfiltered_space_count%max_track, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_max_xtrack", granule%QA_unfiltered_space_count%max_xtrack, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_min_track", granule%QA_unfiltered_space_count%min_track, context)
call nc_get_variable(swid, "QA_unfiltered_space_count_min_xtrack", granule%QA_unfiltered_space_count%min_xtrack, context)
call nc_get_variable(swid, "QA_cal_coef_a0_min", granule%QA_cal_coef_a0%min, context)
call nc_get_variable(swid, "QA_cal_coef_a0_max", granule%QA_cal_coef_a0%max, context)
call nc_get_variable(swid, "QA_cal_coef_a0_mean", granule%QA_cal_coef_a0%mean, context)
call nc_get_variable(swid, "QA_cal_coef_a0_dev", granule%QA_cal_coef_a0%dev, context)
call nc_get_variable(swid, "QA_cal_coef_a0_num", granule%QA_cal_coef_a0%num, context)
call nc_get_variable(swid, "QA_cal_coef_a0_num_bad", granule%QA_cal_coef_a0%num_bad, context)
call nc_get_variable(swid, "QA_cal_coef_a0_max_track", granule%QA_cal_coef_a0%max_track, context)
call nc_get_variable(swid, "QA_cal_coef_a0_max_xtrack", granule%QA_cal_coef_a0%max_xtrack, context)
call nc_get_variable(swid, "QA_cal_coef_a0_min_track", granule%QA_cal_coef_a0%min_track, context)
call nc_get_variable(swid, "QA_cal_coef_a0_min_xtrack", granule%QA_cal_coef_a0%min_xtrack, context)
call nc_get_variable(swid, "QA_cal_coef_a1_min", granule%QA_cal_coef_a1%min, context)
call nc_get_variable(swid, "QA_cal_coef_a1_max", granule%QA_cal_coef_a1%max, context)
call nc_get_variable(swid, "QA_cal_coef_a1_mean", granule%QA_cal_coef_a1%mean, context)
call nc_get_variable(swid, "QA_cal_coef_a1_dev", granule%QA_cal_coef_a1%dev, context)
call nc_get_variable(swid, "QA_cal_coef_a1_num", granule%QA_cal_coef_a1%num, context)
call nc_get_variable(swid, "QA_cal_coef_a1_num_bad", granule%QA_cal_coef_a1%num_bad, context)
call nc_get_variable(swid, "QA_cal_coef_a1_max_track", granule%QA_cal_coef_a1%max_track, context)
call nc_get_variable(swid, "QA_cal_coef_a1_max_xtrack", granule%QA_cal_coef_a1%max_xtrack, context)
call nc_get_variable(swid, "QA_cal_coef_a1_min_track", granule%QA_cal_coef_a1%min_track, context)
call nc_get_variable(swid, "QA_cal_coef_a1_min_xtrack", granule%QA_cal_coef_a1%min_xtrack, context)
call nc_get_variable(swid, "QA_cal_coef_a2_min", granule%QA_cal_coef_a2%min, context)
call nc_get_variable(swid, "QA_cal_coef_a2_max", granule%QA_cal_coef_a2%max, context)
call nc_get_variable(swid, "QA_cal_coef_a2_mean", granule%QA_cal_coef_a2%mean, context)
call nc_get_variable(swid, "QA_cal_coef_a2_dev", granule%QA_cal_coef_a2%dev, context)
call nc_get_variable(swid, "QA_cal_coef_a2_num", granule%QA_cal_coef_a2%num, context)
call nc_get_variable(swid, "QA_cal_coef_a2_num_bad", granule%QA_cal_coef_a2%num_bad, context)
call nc_get_variable(swid, "QA_cal_coef_a2_max_track", granule%QA_cal_coef_a2%max_track, context)
call nc_get_variable(swid, "QA_cal_coef_a2_max_xtrack", granule%QA_cal_coef_a2%max_xtrack, context)
call nc_get_variable(swid, "QA_cal_coef_a2_min_track", granule%QA_cal_coef_a2%min_track, context)
call nc_get_variable(swid, "QA_cal_coef_a2_min_xtrack", granule%QA_cal_coef_a2%min_xtrack, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_min", granule%QA_bb_raw_noise_counts%min, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_max", granule%QA_bb_raw_noise_counts%max, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_mean", granule%QA_bb_raw_noise_counts%mean, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_dev", granule%QA_bb_raw_noise_counts%dev, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_num", granule%QA_bb_raw_noise_counts%num, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_num_bad", granule%QA_bb_raw_noise_counts%num_bad, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_max_track", granule%QA_bb_raw_noise_counts%max_track, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_max_xtrack", granule%QA_bb_raw_noise_counts%max_xtrack, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_min_track", granule%QA_bb_raw_noise_counts%min_track, context)
call nc_get_variable(swid, "QA_bb_raw_noise_counts_min_xtrack", granule%QA_bb_raw_noise_counts%min_xtrack, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_min", granule%QA_sv_raw_noise_counts%min, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_max", granule%QA_sv_raw_noise_counts%max, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_mean", granule%QA_sv_raw_noise_counts%mean, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_dev", granule%QA_sv_raw_noise_counts%dev, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_num", granule%QA_sv_raw_noise_counts%num, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_num_bad", granule%QA_sv_raw_noise_counts%num_bad, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_max_track", granule%QA_sv_raw_noise_counts%max_track, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_max_xtrack", granule%QA_sv_raw_noise_counts%max_xtrack, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_min_track", granule%QA_sv_raw_noise_counts%min_track, context)
call nc_get_variable(swid, "QA_sv_raw_noise_counts_min_xtrack", granule%QA_sv_raw_noise_counts%min_xtrack, context)
call nc_get_variable(swid, "state1", granule%state1, context)
call nc_get_variable(swid, "state2", granule%state2, context)
call nc_get_variable(swid, "cal_coef_a0", granule%cal_coef_a0, context)
call nc_get_variable(swid, "cal_coef_a0_err", granule%cal_coef_a0_err, context)
call nc_get_variable(swid, "cal_coef_a1", granule%cal_coef_a1, context)
call nc_get_variable(swid, "cal_coef_a1_err", granule%cal_coef_a1_err, context)
call nc_get_variable(swid, "cal_coef_a2", granule%cal_coef_a2, context)
call nc_get_variable(swid, "cal_coef_a2_err", granule%cal_coef_a2_err, context)
! call nc_get_variable(swid, "a1_ColdCalPstion", granule%a1_ColdCalPstion, context)
! call nc_get_variable(swid, "a2_ColdCalPstion", granule%a2_ColdCalPstion, context)
! call nc_get_variable(swid, "a1_PLO_Redundncy", granule%a1_PLO_Redundncy, context)
! call nc_get_variable(swid, "a11_mux_temp_used", granule%a11_mux_temp_used, context)
call nc_get_variable(swid, "a11_receiver_temp", granule%a11_receiver_temp, context)
call nc_get_variable(swid, "a11_target_temp", granule%a11_target_temp, context)
! call nc_get_variable(swid, "a12_mux_temp_used", granule%a12_mux_temp_used, context)
call nc_get_variable(swid, "a12_receiver_temp", granule%a12_receiver_temp, context)
call nc_get_variable(swid, "a12_target_temp", granule%a12_target_temp, context)
! call nc_get_variable(swid, "a2_diplexer_temp_used", granule%a2_diplexer_temp_used, context)
call nc_get_variable(swid, "a2_receiver_temp", granule%a2_receiver_temp, context)
call nc_get_variable(swid, "a2_target_temp", granule%a2_target_temp, context)
! call nc_get_variable(swid, "qa_scanline", granule%qa_scanline, context)
! call nc_get_variable(swid, "qa_receiver_a11", granule%qa_receiver_a11, context)
! call nc_get_variable(swid, "qa_receiver_a12", granule%qa_receiver_a12, context)
! call nc_get_variable(swid, "qa_receiver_a2", granule%qa_receiver_a2, context)
! call nc_get_variable(swid, "qa_channel", granule%qa_channel, context)

! Final clean-up
call nc_close_file(swid,context)

return

end subroutine read_amsua_bt_netCDF_granule


!-----------------------------------------------------------------------------
!>  extract the requested AMSUA channel observations from a granule
!>  and convert to DART observation format.  allow caller to specify
!>  a bounding box and only extract data within that region.

subroutine make_obs_sequence (seq, granule, lon1, lon2, lat1, lat2, &
                             use_channels, scan_thin, xtrack_thin)

type(obs_sequence_type), intent(inout) :: seq
type(amsua_bt_granule),  intent(in)    :: granule
real(r8),                intent(in)    :: lon1, lon2, lat1, lat2
logical,                 intent(in)    :: use_channels(:)
integer,                 intent(in)    :: scan_thin, xtrack_thin

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
    call error_handler(E_ERR,routine,string1,source)
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
                            scan_thin, xtrack_thin, obs, prev_obs,  &
                            is_first_obs, pre_time, granule)

end subroutine make_obs_sequence


!-----------------------------------------------------------------------------
!> routine to extract required parts of an AMSUA granule needed for a
!> brightness temperature radiance observation.

subroutine add_granule_observations(seq, lon1, lon2, lat1, lat2, use_channels,  &
                            scan_thin, xtrack_thin, obs, prev_obs, is_first_obs, &
                            pre_time, swath)

type(obs_sequence_type),    intent(inout) :: seq
real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
logical,                    intent(in)    :: use_channels(:)
integer,                    intent(in)    :: scan_thin, xtrack_thin
type(obs_type),             intent(inout) :: obs, prev_obs
logical,                    intent(inout) :: is_first_obs
type(time_type),            intent(inout) :: pre_time
type(amsua_bt_granule),     intent(in)    :: swath

integer :: ipix, iscan, ichan
integer :: days, seconds
integer :: obs_num, key
integer :: which_vert 

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_err
real(r8) :: rqc

real(r8) :: lam1, lam2, phi1, phi2, x, y

real(r8) :: sat_az, sat_ze
integer :: platform_id, sat_id, sensor_id
real(r8) :: mag_field, cosbk
real(r8) :: fastem_p1
real(r8) :: fastem_p2
real(r8) :: fastem_p3
real(r8) :: fastem_p4
real(r8) :: fastem_p5

type(time_type) :: obs_time, time_base
character(len=*), parameter :: routine = 'add_granule_observations:'

integer :: robstype

real(digits12) :: time_offset(size(swath%Time,1),size(swath%Time,2))

! things known (and constant) about the input data and rttov
! consistent with values in rttov_sensor_db.csv

platform_id = 9  ! EOS
sat_id      = 2  ! satellite
sensor_id   = 3  ! AMSUA 

!>@todo FIXME not happy with the time - seems like the swath%start_Time
!       and the time from the data are off by 10 seconds. leaving this
!       code commented out for more exploration later.

!TJH print*,'swath%start_Time ',swath%start_Time
!TJH print*,'min from data    ',minval(swath%Time)
!TJH print*,'max from data    ',maxval(swath%Time)
!TJH print*,'swath%end_Time   ',swath%end_Time
!TJH 
!TJH obs_time = set_date(1993, 1, 1) + set_time(int(swath%start_Time))
!TJH call print_date(obs_time,'start date from 1993')
!TJH call print_time(obs_time,'start time from 1993')
!TJH obs_time = set_date(1993, 1, 1) + set_time(int(swath%end_Time))
!TJH call print_date(obs_time,'end   date from 1993')
!TJH call print_time(obs_time,'end   time from 1993')

! The global metadata has start year.month.day... and start_Time
! The data (Time) have the same base as start_Time.
! So I can convert the data time to an offset to get the number of
! seconds since the start_Time

time_base = set_date(swath%start_year, swath%start_month,      swath%start_day, &
                     swath%start_hour, swath%start_minute, int(swath%start_sec))

time_offset = swath%Time - swath%start_Time

!TJH call print_date(time_base,'start date from global metadata')
!TJH call print_time(time_base,'start time from global metadata')
!TJH obs_time = time_base + set_time(int(minval(time_offset)), 0)
!TJH call print_date(obs_time,'minimum date ')
!TJH call print_time(obs_time,'minimum time ')
!TJH obs_time = time_base + set_time(int(maxval(time_offset)), 0)
!TJH call print_date(obs_time,'maximum date ')
!TJH call print_time(obs_time,'maxiumu time ')

!------------------------------------------------------------------------------
!  loop over all observations within the file

obs_num = 0
which_vert = VERTISUNDEF

! assign each observation the correct observation type
robstype = get_index_for_type_of_obs('EOS_2_AMSUA_TB')
if (robstype < 1) then
   string1 = 'unknown observation type EOS_2_AMSUA_TB'
   call error_handler(E_ERR,routine,string1,source)
endif

!        Latitude(                  AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
!       Longitude(                  AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
! brightness_temp(AMSUA_BT_CHANNEL, AMSUA_BT_GEOXTRACK, AMSUA_BT_GEOTRACK)
! rows are along-track, stepping in the direction the satellite is moving

scanloop:  do iscan=1,size(swath%Latitude,2) ! AKA GEOTRACK

   ! if we're going to subset rows, we will cycle here
   if (scan_thin > 0) then
      if (modulo(iscan, scan_thin) /= 0) cycle scanloop
   endif

   ! Reject bad scans here, if desired.
   ! "Bit field for each scanline (bit 0 set if sun glint in scanline; 
   ! bit 1 set if costal crossing in scanline, bit 2 set if some 
   ! channels had excessive NeDT estimated), dimension (45)"

   rqc = swath%qa_scanline(iscan)
   ! if (rqc /= 0) cycle scanloop ! REJECT IF DESIRED

   ! columns are across-track, varying faster than rows.
   pixloop:  do ipix=1,size(swath%Latitude,1) ! AKA GEOXTRACK

      ! if we're going to subset columns, ditto
      if (xtrack_thin > 0) then
         if (modulo(ipix, xtrack_thin) /= 0) cycle pixloop
      endif

      ! observation lat, lon:
      olat  = swath%Latitude (ipix,iscan) ! valid range [ -90.00,  90.00]
      olon  = swath%Longitude(ipix,iscan) ! valid range [-180.00, 180.00]

      ! verify the location is not outside valid limits.  AIRS  uses -180/180
      if((olon > 180.0_r8) .or. (olon < -180.0_r8) .or.  &
         (olat >  90.0_r8) .or. (olat <  -90.0_r8)) then
         write(*,*)'WARNING : invalid location.  col,row,lon,lat = ', ipix,iscan,olon,olat
         cycle pixloop
      endif

      ! make sure lon is between 0 and 360
      if (olon < 0.0_r8) olon = olon + 360.0_r8

      ! reject observations outside the bounding box (allowing wrapping)
      if(( olat < lat1) .or. ( olat > lat2 ) .or. &
         (.not. is_longitude_between(olon, lon1, lon2))) cycle pixloop

      ! set the zenith angle (aka earth incidence angle)
      sat_ze = swath%satzen(ipix,iscan)

      lam1 = deg2rad*swath%longitude(ipix,iscan)
      lam2 = deg2rad*swath%sat_lon(iscan)
      phi1 = deg2rad*swath%latitude(ipix,iscan)
      phi2 = deg2rad*swath%sat_lat(iscan)

      ! calculate the bearing between the obs lat/lon and the SClat/lon
      y = sin(lam2-lam1)*cos(phi2)
      x = cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lam2-lam1)
      sat_az = rad2deg*atan2(y,x)

      ! convert observation time to DART format 
      obs_time = time_base + set_time(int(time_offset(ipix,iscan)), 0)
      call get_time(obs_time, seconds, days)

      channel_loop: do ichan=1, size(use_channels)

         if (.not. use_channels(ichan)) cycle channel_loop

         ! create the radiance obs for this observation, add to sequence

         ! apparently -9999 is missing data, outside of qc mechanism
         ! FIXME ... there is a '_FillValue' 
         obs_value = swath%brightness_temp(ichan,ipix,iscan)
         if (obs_value < 0.0_r8) cycle channel_loop

         obs_err = swath%brightness_temp_err(ichan,ipix,iscan)

         ! column integrated value, so no vertical location
         vloc = 0.0_r8

         if (get_rttov_option_logical('use_zeeman')) then

            ! From RTTOV v12 Users Guide Section 8.12: "For microwave sensors 
            ! that have high peaking weighting functions in the mesosphere ..., 
            ! channels close to lines of molecular oxygen may be significantly 
            ! affected by the redistribution of line intensity through Zeeman 
            ! splitting as described in the RTTOV v10 Science and Validation Report. 
            ! The absorption for the affected channels will depend on the strength 
            ! and orientation of the magnetic field. You must specify two input 
            ! variables for the geomagnetic field in the rttov_profiles structure, 
            ! these being the magnitude, Be, of the field and the cosine, cosbk, 
            ! of the angle between the field vector and the viewing path considered. 
            ! ...
            ! For AMSU-A, only channel 14 is affected. This channel ... impact is 
            ! therefore much smaller (~0.5K). If a Zeeman coefficient file is used,
            ! then a small set of additional predictors will be included. These will 
            ! contribute for channel 14 but will be nullified for the other channels 
            ! by zero coefficients."
            !
            ! From RTTOV v12 Users Guide Table 12 (p 46):
            ! profiles(i)%Be    is the Earth magnetic field strength in Gauss
            ! profiles(i)%cosbk is the cosine of the angle between the Earth magnetic
            !                   field and wave propagation direction

            write(string1,*) 'Compensating for the Zeeman effect is not supported for AMSU-A observations.'
            write(string2,*) 'See https://github.com/NCAR/DART/issues/99#'
            write(string3,*) 'Also Section 8.12 of the RTTOV v12 user guide.'
            call error_handler(E_ERR, routine, string1, &
                       source, text2=string2, text3=string3)
         else
            mag_field = MISSING_R8
            cosbk = MISSING_R8
         end if

         ! FIXME: consider adding an atlas for looking up FASTEM parameters
         ! From the RTTOV User guide:
         ! Surface type
         ! Typical RTTOV default for land: 3.0, 5.0, 15.0, 0.1, 0.3
         ! Summer land surface:
         !        Forest: 1.7, 1.0, 163.0, 0.0, 0.5
         !        Open grass: 2.2, 1.3, 138.0, 0.0, 0.42
         !        Bare soil: 2.3, 1.9, 21.8, 0.0, 0.5
         ! Winter surface type:
         !        Forest and snow: 2.9, 3.4, 27.0, 0.0, 0.0
         !        Deep dry snow: 3.0, 24.0, 60.0, 0.1, 0.15
         !        Frozen soil: 117.8, 2.0, 0.19, 0.2 ,0.35
         ! Sea ice
         !        Grease ice: 23.7, 7.7, 17.3, 0.0, 0.15
         !        Baltic nilas: 1.6, 3.3, 2.2, 0.0, 0.0
         !        New ice (no snow): 2.9, 3.4, 27.0, 0.0, 0.0
         !        New ice (snow): 2.2, 3.7, 122.0, 0.0, 0.15
         !        Brash ice: 3.0, 5.5, 183.0, 0.0, 0.0
         !        Compact pack ice: 2.0, 1700000.0, 49000000.0, 0.0, 0.0
         !        Fast ice: 1.5, 77.8, 703.0, 0.1, 0.35
         !        Lake ice + snow: 1.8, 67.1, 534.0, 0.1, 0.15
         !        Multi-year ice: 1.5, 85000.0, 4700000.0, 0.0, 0.0
         fastem_p1 = 3.0d0
         fastem_p2 = 5.0d0
         fastem_p3 = 15.0d0
         fastem_p4 = 0.1d0
         fastem_p5 = 0.3d0

         call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, & 
                              ichan, mag_field, cosbk, &
                              fastem_p1, fastem_p2, fastem_p3, fastem_p4, fastem_p5)

         call create_3d_obs(olat, olon, vloc, which_vert, obs_value, robstype, &
                            obs_err, days, seconds, rqc, obs, key)

         call add_obs_to_seq(seq, obs, obs_time, prev_obs, pre_time, is_first_obs)

         obs_num = obs_num + 1

      enddo channel_loop
   enddo pixloop
enddo scanloop

! Print a little summary
if (verbosity > 2) call print_obs_seq_summary(seq)

write(string1,*) 'Converted ',obs_num,' obs for ',trim(swath%instrument), &
                 '; total obs = ',key
call error_handler(E_MSG, routine, string1, source)

end subroutine add_granule_observations


!-----------------------------------------------------------------------------
!> routine to confirm granule metadata

subroutine dump_attributes(granule,filename)
type (amsua_bt_granule), intent(in) :: granule
character(len=*),        intent(in) :: filename

print *,' Summary for file "'//trim(filename)//'" :'
print *,' AutomaticQAFlag               "'//trim(granule%AutomaticQAFlag)//'"'
print *,' DayNightFlag                  "'//trim(granule%DayNightFlag)//'"'
print *,' LatGranuleCen                  ', granule%LatGranuleCen
print *,' LocTimeGranuleCen              ', granule%LocTimeGranuleCen
print *,' LonGranuleCen                  ', granule%LonGranuleCen
print *,' MoonInViewMWCount              ', granule%MoonInViewMWCount
print *,' NumBadData                     ', granule%NumBadData
print *,' NumLandSurface                 ', granule%NumLandSurface
print *,' NumMissingData                 ', granule%NumMissingData
print *,' NumOceanSurface                ', granule%NumOceanSurface
print *,' NumProcessData                 ', granule%NumProcessData
print *,' NumSpecialData                 ', granule%NumSpecialData
print *,' NumTotalData                   ', granule%NumTotalData
print *,' eq_x_longitude                 ', granule%eq_x_longitude
print *,' eq_x_tai                       ', granule%eq_x_tai
print *,' granule_number                 ', granule%granule_number
print *,' instrument                    "'//trim(granule%instrument)//'"'
print *,' node_type                     "'//trim(granule%node_type)//'"'
print *,' num_data_gaps_a1               ', granule%num_data_gaps_a1
print *,' num_data_gaps_a2               ', granule%num_data_gaps_a2
print *,' num_demgeoqa                   ', granule%num_demgeoqa
print *,' num_fpe                        ', granule%num_fpe
print *,' num_ftptgeoqa                  ', granule%num_ftptgeoqa
print *,' num_glintgeoqa                 ', granule%num_glintgeoqa
print *,' num_instr_mode_changes_a1      ', granule%num_instr_mode_changes_a1
print *,' num_instr_mode_changes_a2      ', granule%num_instr_mode_changes_a2
print *,' num_missing_scanlines_a1       ', granule%num_missing_scanlines_a1
print *,' num_missing_scanlines_a2       ', granule%num_missing_scanlines_a2
print *,' num_moongeoqa                  ', granule%num_moongeoqa
print *,' num_satgeoqa                   ', granule%num_satgeoqa
print *,' num_scanlines                  ', granule%num_scanlines
print *,' num_scanlines_not_norm_mode_a1 ', granule%num_scanlines_not_norm_mode_a1
print *,' num_scanlines_not_norm_mode_a2 ', granule%num_scanlines_not_norm_mode_a2
print *,' num_scanlines_rec_cal_prob_a11 ', granule%num_scanlines_rec_cal_prob_a11
print *,' num_scanlines_rec_cal_prob_a12 ', granule%num_scanlines_rec_cal_prob_a12
print *,' num_scanlines_rec_cal_prob_a2  ', granule%num_scanlines_rec_cal_prob_a2
print *,' num_scanlines_sig_coast_xing   ', granule%num_scanlines_sig_coast_xing
print *,' num_scanlines_sig_sun_glint    ', granule%num_scanlines_sig_sun_glint
print *,' num_scansets                   ', granule%num_scansets
print *,' num_zengeoqa                   ', granule%num_zengeoqa
print *,' orbit_path                     ', granule%orbit_path
print *,' orbitgeoqa                     ', granule%orbitgeoqa
print *,' processing_level              "'//trim(granule%processing_level)//'"'
print *,' start_Latitude                 ', granule%start_Latitude
print *,' end_Latitude                   ', granule%end_Latitude
print *,' start_Longitude                ', granule%start_Longitude
print *,' end_Longitude                  ', granule%end_Longitude
print *,' start_Time                     ', granule%start_Time
print *,' end_Time                       ', granule%end_Time
print *,' start_year                     ', granule%start_year
print *,' start_month                    ', granule%start_month
print *,' start_day                      ', granule%start_day
print *,' start_hour                     ', granule%start_hour
print *,' start_minute                   ', granule%start_minute
print *,' start_sec                      ', granule%start_sec
print *,' start_orbit                    ', granule%start_orbit
print *,' end_orbit                      ', granule%end_orbit
print *,' start_orbit_row                ', granule%start_orbit_row
print *,' end_orbit_row                  ', granule%end_orbit_row

end subroutine dump_attributes


!-------------------------------------------------------------------------------
!> 

subroutine define_amsua_dimensions(ncid, context)

integer,                    intent(in) :: ncid
character(len=*), optional, intent(in) :: context

call error_handler(E_ERR,'define_amsua_dimensions','routine not verified',source)

! These names are chosen to match the output of h4tonccf_nc4
call nc_define_dimension(ncid, "GeoXTrack", AMSUA_BT_GEOXTRACK, context)
call nc_define_dimension(ncid, "GeoTrack",  AMSUA_BT_GEOTRACK, context)
call nc_define_dimension(ncid, "Channel",   AMSUA_BT_CHANNEL, context)
call nc_define_unlimited_dimension(ncid, "granule_number", context)

end subroutine define_amsua_dimensions


!------------------------------------------------------------------------------
!> The attributes I am using come from :
!> https://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/AIRS/\
!> 3.3_ScienceDataProductDocumentation/3.3.4_ProductGenerationAlgorithms/\
!> README.AIRABRAD.pdf

subroutine define_amsua_variables(granule, ncid, context)

type(amsua_bt_granule),     intent(in) :: granule
integer,                    intent(in) :: ncid
character(len=*), optional, intent(in) :: context

character(len=NF90_MAX_NAME) :: dimnames(NF90_MAX_VAR_DIMS)
character(len=512) :: string1, string2

call error_handler(E_ERR,'define_amsua_variables','routine not finished',source)

! There are 240 granules per day

call nc_define_integer_variable(ncid,'granule_number','granule_number',context)

! The declarations happen in reverse order from the Fortran.
! The old 'column major' vs. 'row major' argument.
! Fortran slowest is on left, C (and the netCDF libs), slowest on right.

dimnames(1) = 'GeoXTrack' 
dimnames(2) = 'GeoTrack'
dimnames(3) = 'granule_number'

call nc_define_double_variable(ncid,'Latitude' ,dimnames(1:3),context)
call nc_define_double_variable(ncid,'Longitude',dimnames(1:3),context)
call nc_define_double_variable(ncid,'Time'     ,dimnames(1:3),context)

write(string1,*)'AIRS spot boresight geodetic latitude'
write(string2,*)'degrees North'
call nc_add_attribute_to_variable(ncid,'Latitude','long_name',string1,context)
call nc_add_attribute_to_variable(ncid,'Latitude','units'    ,string2,context)
call nc_add_attribute_to_variable(ncid,'Latitude','_FillValue',-9999.0_digits12,context)
call nc_add_attribute_to_variable(ncid,'Latitude','valid_range', &
          (/ -90.0_digits12, 90.0_digits12 /),context)

write(string1,*)'AIRS spot boresight geodetic longitude'
write(string2,*)'degrees East'
call nc_add_attribute_to_variable(ncid,'Longitude','long_name',string1,context)
call nc_add_attribute_to_variable(ncid,'Longitude','units'    ,string2,context)
call nc_add_attribute_to_variable(ncid,'Longitude','_FillValue',-9999.0_digits12,context)
call nc_add_attribute_to_variable(ncid,'Longitude' ,'valid_range', &
          (/ -180.0_digits12, 180.0_digits12 /),context)

write(string1,*)'Footprint (shutter) TAI Time'
write(string2,*)'seconds since 1993-01-01 00:00:00'
call nc_add_attribute_to_variable(ncid,'Time','long_name',string1,context)
call nc_add_attribute_to_variable(ncid,'Time','units'    ,string2,context)
call nc_add_attribute_to_variable(ncid,'Time','_FillValue',-9999.0_digits12,context)

dimnames(1) = 'AMSUA_BT_CHANNEL'
dimnames(2) = 'AMSUA_BT_GEOXTRACK' 
dimnames(3) = 'AMSUA_BT_GEOTRACK'
dimnames(4) = 'granule_number'

call nc_define_real_variable(     ncid,'brightness_temp', dimnames(1:4), context)
call nc_add_attribute_to_variable(ncid,'brightness_temp','units','degrees Kelvin',context)
call nc_add_attribute_to_variable(ncid,'brightness_temp','_FillValue',-9999.0_r4,context)

end subroutine define_amsua_variables


!-------------------------------------------------------------------------------
!> 

subroutine fill_amsua_variables(granule, ncid, context)
type(amsua_bt_granule),     intent(in) :: granule
integer,                    intent(in) :: ncid
character(len=*), optional, intent(in) :: context

call error_handler(E_ERR,'fill_amsua_variables','routine not finished',source)

call nc_put_variable(ncid, 'Latitude',        granule%Latitude,        context)
call nc_put_variable(ncid, 'Longitude',       granule%Longitude,       context)
call nc_put_variable(ncid, 'Time',            granule%Time,            context)
call nc_put_variable(ncid, 'brightness_temp', granule%brightness_temp, context)

end subroutine fill_amsua_variables


!------------------------------------------------------------------------------
!> Determine max possible obs from all files. 
!> Allocating way too much results in excessive run-time and memory use.

function max_possible_obs(num_files, scan_thin, xtrack_thin, channels)

integer, intent(in) :: num_files
integer, intent(in) :: scan_thin
integer, intent(in) :: xtrack_thin
logical, intent(in) :: channels(:)
integer :: max_possible_obs

integer :: thin_factor, num_channels, ichan

num_channels = 0
do ichan = 1, size(channels)
   if (channels(ichan)) num_channels = num_channels + 1
enddo

thin_factor    = 1  ! no thinning

if (  scan_thin > 0) thin_factor = thin_factor *   scan_thin 
if (xtrack_thin > 0) thin_factor = thin_factor * xtrack_thin

max_possible_obs = AMSUA_BT_GEOXTRACK * AMSUA_BT_GEOTRACK * num_channels * &
                   num_files / thin_factor

end function max_possible_obs


!-------------------------------------------------------------------------------
!> Given two observation sequences, insert the observations from the first
!> sequence into the second sequence.

function combine_sequences(seq_in, seq_out, filename)

! The logic of this routine is loosely taken from the obs_sequence_tool

type(obs_sequence_type), intent(in)    :: seq_in
type(obs_sequence_type), intent(inout) :: seq_out
character(len=*),        intent(in)    :: filename
integer                                :: combine_sequences

integer :: num_inserted
logical :: is_there_one, is_this_last

type(obs_type) :: obs_in,  next_obs_in
type(obs_type) :: obs_out, prev_obs_out

! could check to make sure metadata is compatible, skipping for now
integer :: num_copies_in  = 1
integer :: num_qc_in      = 1
integer :: num_copies_out = 1
integer :: num_qc_out     = 1

call init_obs(     obs_in,  num_copies_in,  num_qc_in )
call init_obs(next_obs_in,  num_copies_in,  num_qc_in )
call init_obs(     obs_out, num_copies_out, num_qc_out)
call init_obs(prev_obs_out, num_copies_out, num_qc_out)

! Insert the first observation (which is slow) and then use that
! as the starting point to append the rest of the observations.

num_inserted = 0
is_there_one = get_first_obs(seq_in, obs_in)

if (.not. is_there_one )  then
   write(string1,*)'no first observation in ',trim(filename)
   call error_handler(E_MSG,'combine_sequences', string1)
endif

! insert_obs_in_seq() CHANGES the obs passed in.
! Must pass a copy of incoming obs to insert_obs_in_seq.

obs_out = obs_in
call insert_obs_in_seq(seq_out,obs_out)

prev_obs_out = obs_out            ! records new position in seq_out
num_inserted = num_inserted + 1

call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

ObsLoop : do while (.not. is_this_last)

   obs_in = next_obs_in   ! essentially records position in seq_out

   ! Since the stride through the observation sequence file is always
   ! guaranteed to be in temporally-ascending order, we can use the
   ! 'previous' observation as the starting point to search for the
   ! correct insertion point.  This speeds up the insert code a lot.

   obs_out = obs_in
   call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)

   prev_obs_out = obs_out    ! update position in seq_in for next insert
   num_inserted = num_inserted + 1

   call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

enddo ObsLoop

combine_sequences = num_inserted

end function combine_sequences


!-------------------------------------------------------------------------------
!> either read existing obs_seq or create a new one

subroutine append_or_create(append_output, outputfile, max_num, seq)

logical,                 intent(in)  :: append_output
character(len=*),        intent(in)  :: outputfile
integer,                 intent(in)  :: max_num
type(obs_sequence_type), intent(out) :: seq

logical :: file_exist
integer :: i

! one observation data value and one quality control value
integer, parameter :: NUM_COPIES = 1
integer, parameter :: NUM_QC     = 1

inquire(file=outputfile, exist=file_exist)

if ( file_exist .and. append_output ) then
  call read_obs_seq(outputfile, 0, 0, max_num, seq)
  write(string1,*)'Appending to "'//trim(outputfile)//'"'
  write(string2,*)'Initially has ',get_num_obs(seq),' observations.'
else
  call init_obs_sequence(seq, NUM_COPIES, NUM_QC, max_num)
  do i = 1, NUM_COPIES
    call set_copy_meta_data(seq, i, 'observation')
  end do
  do i = 1, NUM_QC
    call set_qc_meta_data(seq, i, 'QC')
  end do
  write(string1,*)'Creating "'//trim(outputfile)//'" from scratch.'
  write(string2,*)'Initially has ',get_num_obs(seq),' observations.'
endif

if (verbosity > 0) call error_handler(E_MSG,source,string1,text2=string2)

end subroutine append_or_create

!-------------------------------------------------------------------------------

end module amsua_netCDF_support_mod

