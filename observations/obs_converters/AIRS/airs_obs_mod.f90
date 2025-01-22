! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module airs_obs_mod

use types_mod,        only : r4, r8, digits12, deg2rad, rad2deg

use obs_def_mod,      only : obs_def_type, read_obs_def, write_obs_def,        &
                             set_obs_def_time, set_obs_def_type_of_obs,        &
                             set_obs_def_error_variance, set_obs_def_location, &
                             get_obs_def_location

use time_manager_mod, only : time_type, get_date, set_date,                    &
                             get_time, set_time, set_calendar_type,            &
                             GREGORIAN, print_date, print_time,                &
                             operator(+), operator(>=)

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_MSG, is_longitude_between

use     location_mod, only : location_type, set_location, VERTISPRESSURE,      &
                             get_location

use     obs_kind_mod, only : get_index_for_type_of_obs,                        &
                             AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY

use obs_sequence_mod, only : init_obs_sequence, init_obs, obs_sequence_type,   &
                             obs_type, set_copy_meta_data,                     &
                             set_qc_meta_data, print_obs_seq_summary

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use      obs_err_mod, only : rawin_temp_error

use     airs_JPL_mod, only : AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK,            &
                             AIRS_RET_STDPRESSURELEV, AIRS_RET_STDPRESSURELAY, &                 
                             AIRS_RET_AIRSXTRACK, AIRS_RET_AIRSTRACK,          &
                             AIRS_RET_CLOUD, AIRS_RET_H2OPRESSURELEV,          &
                             AIRS_RET_H2OPRESSURELAY, airs_granule_type


implicit none
private

public :: make_obs_sequence, initialize_obs_sequence, compute_thin_factor

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'airs_obs_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.
character(len=512) :: msgstring

logical :: DEBUG = .false.

real(r8), parameter :: mb_to_Pa = 100.0_r8  ! millibars to pascals

! the sizes of the Temperature arrays are:
!   (AIRS_RET_STDPRESSURELAY, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
! the sizes of the MMR arrays are:
!   (AIRS_RET_H2OPRESSURELAY, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
! the rest of the arrays are:
!   (AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)

! MMR is converted into Q:  Q = MMR / (1 + MMR)
real ::   Q    (AIRS_RET_H2OPRESSURELAY, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::   Q_err(AIRS_RET_H2OPRESSURELAY, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)


contains

!-------------------------------------------------------------------------------
!>

subroutine initialize_module

call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module

!-------------------------------------------------------------------------------
!> extract the temperature and humidity observations from a granule
!> and convert to DART observation format.  allow caller to specify
!> a bounding box and only extract data within that region.

subroutine make_obs_sequence ( seq, granule, lon1, lon2, lat1, lat2, &
                             min_MMR_threshold, top_pressure_level, &
                             row_thin, col_thin, use_NCEP_errs, version)

type(obs_sequence_type), intent(inout) :: seq
type(airs_granule_type), intent(in)    :: granule
real(r8), intent(in) :: lon1, lon2, lat1, lat2
real(r8), intent(in) :: min_MMR_threshold, top_pressure_level
integer,  intent(in) :: row_thin, col_thin
logical,  intent(in) :: use_NCEP_errs
integer,  intent(in) :: version

type(obs_type)          :: obs, prev_obs

integer :: i, irow, icol, ivert, istart
integer :: days, seconds
integer :: obs_num, temperature_top_index, humidity_top_index
integer :: which_vert, tobstype, qobstype

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_var
real(r8) :: tqc, qqc
real(r8) :: midpres, log_lower, log_upper

logical :: is_first_obs

type(time_type) :: obs_time, base_time, pre_time, time

if ( .not. module_initialized ) call initialize_module

! Initialize the obs variables
call init_obs(     obs, 1, 1)
call init_obs(prev_obs, 1, 1)

! assign each observation the correct observation type
tobstype = get_index_for_type_of_obs('AIRS_TEMPERATURE')
qobstype = get_index_for_type_of_obs('AIRS_SPECIFIC_HUMIDITY')
if ((tobstype < 1) .or. (qobstype < 1)) then
   msgstring = 'unknown observation type [AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY]'
   call error_handler(E_ERR,'make_obs_sequence',msgstring,source,revision,revdate)
endif

base_time = set_date(1993, 1, 1, 0, 0, 0)   ! reference date: jan 1st, 1993

if (DEBUG) call debug_print_size_check(granule)

! The file contains 'Retrieved Water Vapor Mass Mixing Ratio'.  Convert
! to specific humidity here.  Original units: gm/kg - scale to kg/kg first.

where (granule%H2OMMRStd > min_MMR_threshold) 
   Q = (granule%H2OMMRStd / 1000.0_r8) / ( 1.0_r8 + (granule%H2OMMRStd/1000.0_r8))
elsewhere
   Q = 0.0_r8
endwhere

where (granule%H2OMMRStdErr > 0.0_r8) 
   Q_err = (granule%H2OMMRStdErr / 1000.0_r8) / ( 1.0_r8 + (granule%H2OMMRStdErr/1000.0_r8))
elsewhere
   Q_err = 0.0_r8   ! is this really ok?
endwhere

! see what the top pressure allowed will be, and compute the index number
! that is still below it in each profile.  that will be become the loop end.
! default to doing the whole column
temperature_top_index = size(granule%pressStd)
tloop: do i = 1, size(granule%pressStd)
   if (granule%pressStd(i) <= top_pressure_level) then
      temperature_top_index = i
      exit tloop
   endif
enddo tloop
if (DEBUG) print *, 'temp_top_index = ', temperature_top_index

! temperature obs are on pressure levels.  moisture obs are the mean
! of the layer bounded by the given pressure level below and the next 
! higher layer above, so make sure we are always one level below the
! top.  the loop further down will use i and i+1 as the bounding
! levels for the moisture obs.
humidity_top_index = size(granule%pressH2O) - 1
mloop: do i = 1, size(granule%pressH2O)
   if (granule%pressH2O(i) <= top_pressure_level) then
      humidity_top_index = i - 1
      exit mloop
   endif
enddo mloop
if (DEBUG) print *, 'humd_top_index = ', humidity_top_index

!  loop over all observations within the file
!------------------------------------------------------------------------------

is_first_obs = .true.
obs_num = 1
which_vert = VERTISPRESSURE

! rows are along-track, stepping in the direction the satellite is moving
rowloop:  do irow=1,AIRS_RET_GEOTRACK

   ! if we're going to subset rows, we will cycle here
   if (row_thin > 0) then
      if (modulo(irow, row_thin) /= 1) cycle rowloop
   endif

   ! columns are across-track, varying faster than rows.
   colloop:  do icol=1,AIRS_RET_GEOXTRACK

      ! if we're going to subset columns, ditto
      if (col_thin > 0) then
         if (modulo(icol, col_thin) /= 1) cycle colloop
      endif

      ! observation lat, lon:
      olat  = granule%Latitude (icol,irow) ! valid range [ -90.00,  90.00]
      olon  = granule%Longitude(icol,irow) ! valid range [-180.00, 180.00]

      ! verify the location is not outside valid limits.  AIRS  uses -180/180
      if((olon > 180.0_r8) .or. (olon < -180.0_r8) .or.  &
         (olat >  90.0_r8) .or. (olat <  -90.0_r8)) then
         write(*,*)'WARNING : invalid location.  col,row,lon,lat = ', icol,irow,olon,olat
         cycle colloop
      endif

      ! reject observations outside the bounding box (allowing wrapping)
      if(( olat < lat1) .or. ( olat > lat2 ) .or. &
         (.not. is_longitude_between(olon, lon1, lon2))) cycle colloop

      ! make sure lon is between 0 and 360
      if (olon < 0) olon = olon + 360.0_r8

      obs_time = base_time + set_time(int(granule%Time(icol, irow)))
      call get_time(obs_time, seconds, days)
     
      ! avoid -9999 in the starting vertical level
      istart = granule%nBestStd(icol, irow)
      if (istart < 1) goto 100   ! skip vert_T_loop

      vert_T_loop: do ivert=istart, temperature_top_index

         tqc = 0   ! if we get here, the quality control is 'best' == 0

         ! create the temperature obs for this observation, add to sequence

         ! apparently -9999 is missing data, outside of qc mechanism
         obs_value = granule%TAirStd(ivert, icol, irow)
         if (obs_value == -9999.0_r8) cycle vert_T_loop

         ! temperature values are located directly at the pressure levels
         vloc = granule%pressStd(ivert) * mb_to_Pa

         if (use_NCEP_errs) then
            obs_var = max(rawin_temp_error(vloc / mb_to_Pa), &  ! back to mb for this call
                          granule%TAirStdErr(ivert, icol, irow))
         else
            obs_var = granule%TAirStdErr(ivert, icol, irow)
         endif
         obs_var = obs_var * obs_var

         if (DEBUG) print *, 'obs value, var = ', obs_value, obs_var

         call create_3d_obs(olat, olon, vloc, which_vert, obs_value, AIRS_TEMPERATURE, obs_var, &
                            days, seconds, tqc, obs)
      
         call add_obs_to_seq(seq, obs, time, prev_obs, pre_time, is_first_obs)
   
         if ( DEBUG ) then
            write(*,*)irow, icol, ivert, olat, olon
         endif

         obs_num = obs_num + 1
 
      enddo vert_T_loop

100   continue

      ! avoid -9999 in the starting vertical layer
      istart = granule%nSurfStd(icol, irow)
      if (istart < 1) goto 200  !skip vert_Q_loop

      vert_Q_loop:  do ivert=istart,humidity_top_index

         if ((version.eq.6 .or. version.eq.7) .and. (granule%H2OMMRStd_QC(ivert, icol, irow) > 0)) cycle vert_Q_loop

         qqc = 0   ! if we get here, the quality control is 'Best' == 0

         ! create the moisture obs for this observation, add to sequence
   
         ! if original MMR data was -9999, that is apparently a missing val
         if (granule%H2OMMRStd(ivert, icol, irow) == -9999.0_r8) cycle vert_Q_loop

         obs_value = Q(ivert, icol, irow)

         ! moisture obs are the mean of the layer with the bottom at
         ! the given pressure.  compute the midpoint (in log space)
         ! between this level and the level above it, and use that for
         ! the moisture obs location.  see AIRS.html for more info on
         ! layers vs levels.
         log_lower = log(granule%pressH2O(ivert))
         log_upper = log(granule%pressH2O(ivert+1))
         midpres = exp((log_lower + log_upper) / 2.0_r8)
         vloc = midpres * mb_to_Pa

         if (use_NCEP_errs) then
            obs_var = max(Q_err(ivert, icol, irow), 0.2_r8)  ! ncep routine needs temp & rh but always returns 0.2
         else
            obs_var = Q_err(ivert, icol, irow)
         endif
         obs_var = obs_var * obs_var


         call create_3d_obs(olat, olon, vloc, which_vert, obs_value, AIRS_SPECIFIC_HUMIDITY, obs_var, &
                            days, seconds, qqc, obs)
      
         call add_obs_to_seq(seq, obs, time, prev_obs, pre_time, is_first_obs)
   
         if ( DEBUG ) then
            write(*,*)irow, icol, ivert, olat, olon
         endif

         obs_num = obs_num + 1
 
      enddo vert_Q_loop

200   continue

   enddo colloop
enddo rowloop

! Print a little summary
write(msgstring,*) 'obs used = ', obs_num - 1
call error_handler(E_MSG, ' ', msgstring)

call print_obs_seq_summary(seq)

end subroutine make_obs_sequence


!-------------------------------------------------------------------------------
!>

function compute_thin_factor(along, across)

integer, intent(in) :: along
integer, intent(in) :: across
integer :: compute_thin_factor

if (along > 0 .and. across > 0) then
   compute_thin_factor = along * across
else
   compute_thin_factor = 1
endif

end function compute_thin_factor


!-------------------------------------------------------------------------------
!>

subroutine check_size(size1, size2, varlabel, subrlabel)

integer,          intent(in) :: size1, size2
character(len=*), intent(in) :: varlabel, subrlabel

if (size1 /= size2) then
   write(msgstring, '(A,I6,A,I6)') 'unexpected size '//trim(varlabel)//': ', &
         size1, ' /= ', size2
   call error_handler(E_ERR, trim(subrlabel), msgstring,source,revision,revdate)
endif

end subroutine check_size


!-------------------------------------------------------------------------------
!>

function initialize_obs_sequence(filecount, sample_factor)

integer,          intent(in) :: filecount
integer,          intent(in) :: sample_factor
type(obs_sequence_type) :: initialize_obs_sequence

type(obs_sequence_type) :: seq

! max possible obs from this one granule.
integer :: max_num, num_copies, num_qc

max_num = (AIRS_RET_STDPRESSURELAY * AIRS_RET_GEOXTRACK * AIRS_RET_GEOTRACK + &
           AIRS_RET_H2OPRESSURELAY * AIRS_RET_GEOXTRACK * AIRS_RET_GEOTRACK ) / sample_factor

num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence 
call init_obs_sequence(seq, num_copies, num_qc, max_num * filecount)

! set meta data of obs_seq
call set_copy_meta_data(seq,  1, 'observation')
call set_qc_meta_data(seq, 1, 'AIRS QC')

initialize_obs_sequence = seq

end function initialize_obs_sequence


!-------------------------------------------------------------------------------
!> Check - before we loop over these arrays, make sure
!> they are the size we expect them to be.

subroutine debug_print_size_check(granule)

type(airs_granule_type), intent(in) :: granule

type(time_type) :: obs_time, base_time

base_time = set_date(1993, 1, 1, 0, 0, 0)   ! reference date: jan 1st, 1993

call check_size(size(granule%TAirStd),                           &
   AIRS_RET_STDPRESSURELAY*AIRS_RET_GEOXTRACK*AIRS_RET_GEOTRACK, &
   'TAirStd (T)','make_obs_sequence')

if (DEBUG) print *, 'first  row T', granule%TAirStd(:, 1, 1)
if (DEBUG) print *, 'second row T', granule%TAirStd(:, 2, 1)
if (DEBUG) print *, 'second col T', granule%TAirStd(:, 1, 2)

call check_size(size(granule%TAirStdErr), size(granule%TAirStd), &
   'TAirStdErr (T_err)','make_obs_sequence')

if (DEBUG) print *, 'first  row T err', granule%TAirStdErr(:, 1, 1)

! First (lowest) good pressure level number
call check_size(size(granule%nBestStd), AIRS_RET_GEOXTRACK*AIRS_RET_GEOTRACK, &
   'nBestStd (T QC)','make_obs_sequence')

if (DEBUG) print *, 'first  row nBestStd', granule%nBestStd(:, 1)

call check_size(size(granule%H2OMMRStd),                         &
   AIRS_RET_H2OPRESSURELAY*AIRS_RET_GEOXTRACK*AIRS_RET_GEOTRACK, &
   'H2OMMRStd (MMR)','make_obs_sequence')

if (DEBUG) print *, 'AIRS_RET_H2OPRESSURELEV = ', AIRS_RET_H2OPRESSURELEV
if (DEBUG) print *, 'AIRS_RET_H2OPRESSURELAY = ', AIRS_RET_H2OPRESSURELAY
if (DEBUG) print *, 'first  row MMR', granule%H2OMMRStd(:, 1, 1)
if (DEBUG) print *, 'second col MMR', granule%H2OMMRStd(:, 2, 1)
if (DEBUG) print *, 'second col MMR', granule%H2OMMRStd(:, 1, 2)

call check_size(size(granule%H2OMMRStdErr), size(granule%H2OMMRStd), &
   'H2OMMRStdErr (MMR_err)','make_obs_sequence')

if (DEBUG) print *, 'first  row MMR err', granule%H2OMMRStdErr(:, 1, 1)

call check_size(size(granule%Qual_H2O), AIRS_RET_GEOXTRACK*AIRS_RET_GEOTRACK, &
   'Qual_H2O (Q QC)','make_obs_sequence')

if (DEBUG) print *, 'first  row qual_h2o', granule%Qual_H2O(:, 1)

if (DEBUG) print *, 'first  row Q', Q(:, 1, 1)
if (DEBUG) print *, 'second col Q', Q(:, 2, 1)
if (DEBUG) print *, 'second col Q', Q(:, 1, 2)

call check_size(size(granule%H2OMMRStd_QC), &
   AIRS_RET_H2OPRESSURELAY*AIRS_RET_GEOXTRACK*AIRS_RET_GEOTRACK, &
   'H2OMMRStd_QC (Q QC)','make_obs_sequence')

if (DEBUG) print *, 'first  row H2OMMRStd_QC', granule%H2OMMRStd_QC(:, 1, 1)

if (DEBUG) print *, 'first  row Q', Q(:, 1, 1)
if (DEBUG) print *, 'second col Q', Q(:, 2, 1)
if (DEBUG) print *, 'second col Q', Q(:, 1, 2)

if (DEBUG) print *, 'first  row Q_err', Q_err(:, 1, 1)


if (DEBUG) print *, '1st lat, lon, time', &
           granule%Latitude(1,1), granule%Longitude(1,1), granule%Time(1,1)
if (DEBUG) print *, 'r 2 lat, lon, time', &
           granule%Latitude(2,1), granule%Longitude(2,1), granule%Time(2,1)
if (DEBUG) print *, 'c 2 lat, lon, time', &
           granule%Latitude(1,2), granule%Longitude(1,2), granule%Time(1,2)

! debugging for first/last time of obs in this file
if (DEBUG) print *, 'int start time ', int(granule%Time(1, 1))
if (DEBUG) obs_time = base_time + set_time(int(granule%Time(1, 1)))
if (DEBUG) call print_date(obs_time)

if (DEBUG) print *, 'int end   time ', int(granule%Time(AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK))
if (DEBUG) obs_time = base_time + set_time(int(granule%Time(AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)))
if (DEBUG) call print_date(obs_time)

if (DEBUG) print *, 'pressStd', granule%pressStd
if (DEBUG) print *, 'pressH2O', granule%pressH2O

end subroutine debug_print_size_check



end module airs_obs_mod

