! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module airs_obs_rad_mod

use types_mod,        only : r4, r8, digits12, deg2rad, rad2deg

use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def,     &
                             write_obs_def, destroy_obs_def,                   &
                             interactive_obs_def, copy_obs_def,                &
                             set_obs_def_time, set_obs_def_type_of_obs,        &
                             set_obs_def_error_variance, set_obs_def_location, &
                             get_obs_def_location

use obs_def_rttov_mod, only : set_rttov_metadata

use time_manager_mod, only : time_type, get_date, set_date,            &
                             get_time, set_time, set_calendar_type,    &
                             GREGORIAN, print_date, print_time,        &
                             operator(+), operator(>=)

use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler,              &
                             E_ERR, E_MSG, is_longitude_between

use     location_mod, only : location_type, set_location, VERTISPRESSURE, &
                             get_location, VERTISUNDEF

use     obs_kind_mod, only : get_index_for_type_of_obs, AQUA_AIRS_AMSU_RADIANCE

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type,      &
                             obs_type, copy_obs, set_copy_meta_data,         &
                             set_qc_meta_data, set_obs_def, get_first_obs,   &
                             get_last_obs, get_obs_def

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use obs_seq_utilities_mod, only : print_obs_seq

use airs_rad_L2_mod   ! need ', only' list here

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

! the sizes of the radiance arrays are:
!   (AIRS_CC_RAD_CHANNEL,AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
! the sizes of the az/ze arrays are:
!   (AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)


contains

!-------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module

!------------------------------------------------------------------------------
!  extract the requested radiance channel observations from a granule
!  and convert to DART observation format.  allow caller to specify
!  a bounding box and only extract data within that region.

function make_obs_sequence ( granule, lon1, lon2, lat1, lat2, &
                             nchans, chanlist, row_thin, col_thin)
type(airs_granule_type), intent(in) :: granule
real(r8), intent(in) :: lon1, lon2, lat1, lat2
integer,  intent(in) :: nchans, chanlist(:)
integer,  intent(in) :: row_thin, col_thin
type(obs_sequence_type) :: make_obs_sequence

! max possible obs from this one granule. in practice if the
! real number of processed channels is very much smaller, make
! another parameter so we don't allocate all these unused obs
! (takes time & space) and then delete them at the end.
integer :: max_num =  &
   AIRS_CC_RAD_CHANNEL * AIRS_CC_RAD_GEOXTRACK * AIRS_CC_RAD_GEOTRACK

type(obs_def_type)      :: obs_def
type(obs_type)          :: obs, prev_obs
type(location_type)     :: obs_loc

integer :: i, irow, icol, ichan, this_chan, num_copies, num_qc, istart
integer :: days, seconds
integer :: obs_num, key
integer :: which_vert, robstype

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_err
real(r8) :: rqc

real(r8) :: sat_az, sat_ze, sun_az, sun_ze
integer  :: platform, sat_id, sensor

logical :: is_first_obs
type(time_type) :: obs_time, base_time, pre_time, time
character(len=*), parameter :: routine = 'make_obs_sequence'

if ( .not. module_initialized ) call initialize_module

! one observation data value and one quality control value 
! per obs.  if you change these you have to set additional
! metadata for them below.
num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence 
call init_obs_sequence(make_obs_sequence, num_copies, num_qc, max_num)

! set meta data of obs_seq
call set_copy_meta_data(make_obs_sequence, 1, 'observation')
call set_qc_meta_data(make_obs_sequence, 1, 'AIRS QC')

! Initialize the obs variables
call init_obs(     obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! assign each observation the correct observation type
robstype = get_index_for_type_of_obs('AQUA_AIRS_AMSU_RADIANCE')
if (robstype < 1) then
   msgstring = 'unknown observation type AQUA_AIRS_AMSU_RADIANCE'
   call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
endif

! things known (and constant) about the input data and rttov

base_time = set_date(1993, 1, 1, 0, 0, 0)   ! Data reference date: jan 1st, 1993

platform = 9   ! EOS, before renamed to AQUA
sat_id   = 2   ! verify this
!sensor   = 11  ! AIRS   (amsu-a is 3)
sensor   = 3  ! AIRS   (amsu-a is 3)


!------------------------------------------------------------------------------
!  loop over all observations within the file

is_first_obs = .true.
obs_num = 1
which_vert = VERTISUNDEF

! rows are along-track, stepping in the direction the satellite is moving
rowloop:  do irow=1,AIRS_CC_RAD_GEOTRACK

   ! if we're going to subset rows, we will cycle here
   if (row_thin > 0) then
      if (modulo(irow, row_thin) /= 1) cycle rowloop
   endif

   ! columns are across-track, varying faster than rows.
   colloop:  do icol=1,AIRS_CC_RAD_GEOXTRACK

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

      ! get the geometry info the radiative transfer model will need
      sat_az = granule%satazi(icol,irow)
      sat_ze = granule%satzen(icol,irow)
      sun_az = granule%solazi(icol,irow)
      sun_ze = granule%solzen(icol,irow)

      obs_time = base_time + set_time(int(granule%Time(icol, irow)))
      call get_time(obs_time, seconds, days)
     
      channel_loop: do ichan=1, nchans

         this_chan = chanlist(ichan)

         ! check channel quality control
         rqc = granule%radiances_QC(this_chan, icol, irow)
         print *, 'rqc = ', rqc
         ! FIXME: reject bad scans here - cycle channel_loop

         ! create the radiance obs for this observation, add to sequence

         ! apparently -9999 is missing data, outside of qc mechanism
         obs_value = granule%radiances(this_chan, icol, irow)
print *, 'obs_value = ', obs_value
         !if (obs_value == -9999.0_r8) cycle channel_loop

         obs_err = granule%radiance_err(this_chan, icol, irow) 

print *, 'obs value, err, var = ', obs_value, obs_err, obs_err*obs_err

         ! column integrated value, so no vertical location
         vloc = 0.0_r8

         ! add additional metadata for this obs type.  returns key to use in create call
         call set_rttov_metadata(key, sat_az, sat_ze, sun_az, sun_ze, platform, sat_id, sensor, this_chan)

         call create_3d_obs(olat, olon, vloc, which_vert, obs_value, robstype, &
                            obs_err, days, seconds, rqc, obs, key)
      
         call add_obs_to_seq(make_obs_sequence, obs, time, prev_obs, pre_time, is_first_obs)
   
         obs_num = obs_num + 1
 
      enddo channel_loop
   enddo colloop
enddo rowloop

! Print a little summary
write(msgstring,*) 'obs used = ', obs_num, ' obs skipped = ', max_num - obs_num
call error_handler(E_MSG, ' ', msgstring)

call print_obs_seq(make_obs_sequence, '')

end function make_obs_sequence

!------------------------------------------------------------------------------
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

end module airs_obs_rad_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
