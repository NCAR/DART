! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module airs_obs_mod

use types_mod,        only : r4, r8, digits12, deg2rad, rad2deg

use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def,     &
                             write_obs_def, destroy_obs_def,                   &
                             interactive_obs_def, copy_obs_def,                &
                             set_obs_def_time, set_obs_def_type_of_obs,               &
                             set_obs_def_error_variance, set_obs_def_location, &
                             get_obs_def_location

use time_manager_mod, only : time_type, get_date, set_date,            &
                             get_time, set_time, set_calendar_type,    &
                             GREGORIAN, print_date, print_time,        &
                             operator(+), operator(>=)

use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler,              &
                             E_ERR, E_MSG, is_longitude_between

use     location_mod, only : location_type, set_location, VERTISPRESSURE, &
                             get_location

use     obs_kind_mod, only : get_index_for_type_of_obs, &
                             QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY

use obs_kind_mod,     only : AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type,      &
                             obs_type, copy_obs, set_copy_meta_data,         &
                             set_qc_meta_data, set_obs_def, get_first_obs,   &
                             get_last_obs, get_obs_def

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use obs_seq_utilities_mod, only : print_obs_seq

use airs_JPL_mod, only : AIRS_RET_H2OPRESSURELAY !need ', only' list here

implicit none
private

public :: make_obs_sequence, create_output_filename

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.
character(len=129) :: msgstring

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

!-------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module

!------------------------------------------------------------------------------
!  extract the temperature and humidity observations from a granule
!  and convert to DART observation format.  allow caller to specify
!  a bounding box and only extract data within that region.

function make_obs_sequence ( granule, lon1, lon2, lat1, lat2, &
                             min_MMR_threshold, top_pressure_level, &
                             row_thin, col_thin)
type(airs_granule_type), intent(in) :: granule
real(r8), intent(in) :: lon1, lon2, lat1, lat2
real(r8), intent(in) :: min_MMR_threshold, top_pressure_level
integer,  intent(in) :: row_thin, col_thin
type(obs_sequence_type) :: make_obs_sequence

! max possible obs from this one granule.
integer :: max_num =  &
   AIRS_RET_STDPRESSURELAY * AIRS_RET_GEOXTRACK * AIRS_RET_GEOTRACK + &
   AIRS_RET_H2OPRESSURELAY * AIRS_RET_GEOXTRACK * AIRS_RET_GEOTRACK 

type(obs_def_type)      :: obs_def
type(obs_type)          :: obs, prev_obs
type(location_type)     :: obs_loc

integer :: i, irow, icol, ivert, num_copies, num_qc, istart
integer :: days, seconds
integer :: obs_num, temperature_top_index, humidity_top_index
integer :: which_vert, tobstype, qobstype

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_err
real(r8) :: tqc, qqc, latlon(3)
real(r8) :: midpres, log_lower, log_upper

logical :: is_first_obs
type(time_type) :: obs_time, base_time, pre_time, time

if ( .not. module_initialized ) call initialize_module

num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence 
call init_obs_sequence(make_obs_sequence, num_copies, num_qc, max_num)

! set meta data of obs_seq
do i = 1, num_copies
   call set_copy_meta_data(make_obs_sequence, i, 'observation')
end do

do i = 1, num_qc
   call set_qc_meta_data(make_obs_sequence, i, 'AIRS QC')
end do

! Initialize the obs variables
call init_obs(     obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

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

         obs_err = granule%TAirStdErr(ivert, icol, irow) 

print *, 'obs value, err, var = ', obs_value, obs_err, obs_err*obs_err

         ! temperature values are located directly at the pressure levels
         vloc = granule%pressStd(ivert) * mb_to_Pa

         call create_3d_obs(olat, olon, vloc, which_vert, obs_value, AIRS_TEMPERATURE, obs_err, &
                            days, seconds, tqc, obs)
      
         call add_obs_to_seq(make_obs_sequence, obs, time, prev_obs, pre_time, is_first_obs)
   
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

         if (granule%Qual_H2O(icol, irow) > 0) exit vert_Q_loop

         qqc = 0   ! if we get here, the quality control is 'Best' == 0

         ! create the moisture obs for this observation, add to sequence
   
         ! if original MMR data was -9999, that is apparently a missing val
         if (granule%H2OMMRStd(ivert, icol, irow) == -9999.0_r8) cycle vert_Q_loop

         obs_value = Q(ivert, icol, irow)
         obs_err = Q_err(ivert, icol, irow)

         ! moisture obs are the mean of the layer with the bottom at
         ! the given pressure.  compute the midpoint (in log space)
         ! between this level and the level above it, and use that for
         ! the moisture obs location.  see AIRS.html for more info on
         ! layers vs levels.
         log_lower = log(granule%pressH2O(ivert))
         log_upper = log(granule%pressH2O(ivert+1))
         midpres = exp((log_lower + log_upper) / 2.0_r8)
         vloc = midpres * mb_to_Pa

         call create_3d_obs(olat, olon, vloc, which_vert, obs_value, AIRS_SPECIFIC_HUMIDITY, obs_err, &
                            days, seconds, qqc, obs)
      
         call add_obs_to_seq(make_obs_sequence, obs, time, prev_obs, pre_time, is_first_obs)
   
         if ( DEBUG ) then
            write(*,*)irow, icol, ivert, olat, olon
         endif

         obs_num = obs_num + 1
 
      enddo vert_Q_loop

200   continue

   enddo colloop
enddo rowloop

! Print a little summary
write(msgstring,*) 'obs used = ', obs_num, ' obs skipped = ', max_num - obs_num
call error_handler(E_MSG, ' ', msgstring)

call print_obs_seq(make_obs_sequence, '')

end function make_obs_sequence

!------------------------------------------------------------------------------

subroutine check_size(size1, size2, varlabel, subrlabel)
integer,          intent(in) :: size1, size2
character(len=*), intent(in) :: varlabel, subrlabel

if (size1 /= size2) then
   write(msgstring, '(A,I6,A,I6)') 'unexpected size '//trim(varlabel)//': ', &
         size1, ' /= ', size2
   call error_handler(E_ERR, trim(subrlabel), msgstring,source,revision,revdate)
endif

end subroutine check_size

!-------------------------------------------------
! The L2 filenames have a very long extension that
! records when the data was published - not very interesting
! for our purposes. replace with something DART-y.

subroutine create_output_filename(l2name, ofname)
character(len=*), intent(IN)  :: l2name
character(len=*), intent(OUT) :: ofname

integer :: i, basestart, extstart, strlen

! hardcoded and brittle, but for now...  the first 19 chars
! of the input filename have the date & granule number, which
! seems like the bulk of the useful info.  find the last / and
! copy from there to +19 chars.

strlen = len_trim(l2name)

basestart = 1
slashloop : do i = strlen-1,1,-1
   if (l2name(i:i) == '/' ) then
      basestart = i+1
      exit slashloop
   endif
enddo slashloop

extstart = basestart+19-1

ofname = l2name(basestart:extstart)//'.out'
if (DEBUG) print *, 'output filename = ', ofname

end subroutine create_output_filename

!-------------------------------------------------
! bit of sanity checking - before we loop over these arrays, make sure
! they are the size we expect them to be.

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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
