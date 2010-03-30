! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   text_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISHEIGHT, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use      obs_kind_mod, only : EVAL_U_WIND_COMPONENT, EVAL_V_WIND_COMPONENT, &
                              EVAL_TEMPERATURE

implicit none

character(len=64), parameter :: text_input_file = 'textdata.input'
character(len=64), parameter :: obs_out_file    = 'obs_seq.out'

logical, parameter :: debug = .false.  ! set to .true. to print info

character (len=129) :: input_line

integer :: n, i, oday, osec, rcio, iunit, otype
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs
           
logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc, wdir, wspeed, werr
real(r8) :: lat, lon, vert, uwnd, uerr, vwnd, verr

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

! start of executable code

call initialize_utilities('text_to_obs')

! time setup
call set_calendar_type(GREGORIAN)

!! some times are supplied as number of seconds since some reference
!! date.  This is an example of how to support that.
!! put the reference date into DART format
!comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

! open input text file

iunit = open_file(text_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(text_input_file)


! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 100000
num_copies = 1
num_qc     = 1

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! if you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files 
! once they are in DART obs_seq format.

!  ! existing file found, append to it
!  inquire(file=obs_out_file, exist=file_exist)
!  if ( file_exist ) then
!     call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
!  endif

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

obsloop: do    ! no end limit - have the loop break when input ends

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus model error

   ! assume here a line is a type (1/2), location, time, value, obs error

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=rcio) input_line
   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code from input file, rcio = ', rcio
      exit obsloop
   endif

   ! pull off the first 2 columns as an integer, to decode the type
   read(input_line, "(I2)", iostat=rcio) otype
   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code trying to get obs type, rcio = ', rcio
      exit obsloop
   endif
   
   if (debug) print *, 'next observation type = ', otype

   ! for this example, assume there is an obs type, where otype=1 is
   ! a temperature measured in height, and if otype=2, there's a wind
   ! speed and direction.

   if (otype == 1) then
      read(input_line(3:129), *, iostat=rcio) lat, lon, vert, &
                                 year, month, day, hour, minute, second, &
                                 temp, terr
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code getting rest of temp obs, rcio = ', rcio
         exit obsloop
      endif
   else
      read(input_line(3:129), *, iostat=rcio) lat, lon, vert, &
                                  year, month, day, hour, minute, second, &
                                  wspeed, wdir, werr
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code getting rest of wind obs, rcio = ', rcio
         exit obsloop
      endif
   endif
   
   if (debug) print *, 'next observation located at lat, lon = ', lat, lon

   ! check the lat/lon values to see if they are ok
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle obsloop


   ! if lon comes in between -180 and 180, use these lines instead:
   !if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
   !if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   !! if time is given in seconds since 1/1/1970, here's how to add it.
   !time_obs = comp_day0 + time_obs

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! this example assumes there is an obs type, where otype=1 is
   ! a temperature measured in height, and if otype=2, there's a wind
   ! speed and direction and height is pressure.  any kind of observation
   ! can use any of the vertical types; this is just an example.

   if (otype == 1) then

      ! height is in meters (gph)

      ! make an obs derived type, and then add it to the sequence
      call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
                         EVAL_TEMPERATURE, terr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      if (debug) print *, 'added temperature obs to output seq'
   else

      ! DART usually assimilates wind as 2 separate U and V components
      ! instead of trying to assimilate a vector of speed and direction.
      ! so convert a wind speed & direction into the U and V components
      ! and create 2 obs for it.  assume vert is in mb or hectopascals,
      ! convert to pascals.  DART does assume all pressures are in pascals.
      uwnd = sin(wdir * DEG2RAD) * wspeed
      vwnd = cos(wdir * DEG2RAD) * wspeed
      uerr = sin(wdir * DEG2RAD) * werr
      verr = cos(wdir * DEG2RAD) * werr

      ! convert hectopascals to pascals.
      vert = vert * 100.0_r8

      call create_3d_obs(lat, lon, vert, VERTISPRESSURE, uwnd, &
                         EVAL_U_WIND_COMPONENT, uerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
      call create_3d_obs(lat, lon, vert, VERTISPRESSURE, vwnd, &
                         EVAL_V_WIND_COMPONENT, verr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      if (debug) print *, 'added 2 wind obs to output seq'
    endif

end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.  
!
!       NOTE: assumes the code is using the threed_sphere locations module, 
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, obs)
use        types_mod, only : r8
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 integer,        intent(in)    :: okind, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_3d_obs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_obs_to_seq -- adds an observation to a sequence.  inserts if first
!           obs, inserts with a prev obs to save searching if that's possible.
!
!     seq - observation sequence to add obs to
!     obs - observation, already filled in, ready to add
!     obs_time - time of this observation, in dart time_type format
!     prev_obs - the previous observation that was added to this sequence
!                (will be updated by this routine)
!     prev_time - the time of the previously added observation (will also
!                be updated by this routine)
!     first_obs - should be initialized to be .true., and then will be
!                updated by this routine to be .false. after the first obs
!                has been added to this sequence.
!
!     created Mar 8, 2010   nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_obs_to_seq(seq, obs, obs_time, prev_obs, prev_time, first_obs)
 use        types_mod, only : r8
 use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
 use time_manager_mod, only : time_type, operator(>=)

  type(obs_sequence_type), intent(inout) :: seq
  type(obs_type),          intent(inout) :: obs, prev_obs
  type(time_type),         intent(in)    :: obs_time
  type(time_type),         intent(inout) :: prev_time
  logical,                 intent(inout) :: first_obs

! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.

if(first_obs) then    ! for the first observation, no prev_obs
   call insert_obs_in_seq(seq, obs)
   first_obs = .false.
else               
   if(obs_time >= prev_time) then  ! same time or later than previous obs
      call insert_obs_in_seq(seq, obs, prev_obs)
   else                            ! earlier, search from start of seq
      call insert_obs_in_seq(seq, obs)
   endif
endif

! update for next time
prev_obs = obs
prev_time = obs_time

end subroutine add_obs_to_seq

end program text_to_obs
