! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! An example program to create DART format observations from a text input file
!
! An example program that reads text (ascii) input lines and creates
! DART observations from them.  This is more of a template program.
! You will have to modify it for your own use.  In particular
! you will have to adapt the read routine to match your white-space 
! separated values or fixed-width column data.
!
! This program now includes a small namelist.  The items in the namelist
! can be changed at runtime without having to recompile the program.
! Having the input and output filenames in the namelist should make it
! simpler to script the conversion of a series of files.
!
! history:
!  created 29 Mar 2010   nancy collins NCAR/IMAGe
!  updated  8 Feb 2016   minor changes to wind conversion, and updated
!                        comments to try to be more helpful.
!  updated 29 Aug 2019   added namelist with example values
!

program text_to_obs

use         types_mod, only : r8, PI, DEG2RAD
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG, nmlfileunit,   &
                              do_nml_file, do_nml_term
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISHEIGHT, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data

! Change these to the actual types of observations you have.
use      obs_kind_mod, only : EVAL_U_WIND_COMPONENT, EVAL_V_WIND_COMPONENT, &
                              EVAL_TEMPERATURE

implicit none

! variables that can be changed at runtime from a namelist
! add any other options here that would be useful for your application.
character(len=256) :: text_input_file = '../data/text_input_file'
character(len=256) :: obs_out_file    = 'obs_seq.out'
logical            :: debug           = .false.  ! set to .true. to print info

namelist /text_to_obs_nml/ text_input_file, obs_out_file, debug


! local variables
character (len=129) :: input_line

integer :: oday, osec, rcio, iunit, otype
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs
           
logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc, wdir, wspeed, werr
real(r8) :: lat, lon, vert, uwnd, uerr, vwnd, verr

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: ref_day0, time_obs, prev_time


! start of executable code

call initialize_utilities('text_to_obs')

! Read the namelist entry
call find_namelist_in_file("input.nml", "text_to_obs_nml", iunit)
read(iunit, nml = text_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "text_to_obs_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=text_to_obs_nml)
if (do_nml_term()) write(     *     , nml=text_to_obs_nml)


! time setup
call set_calendar_type(GREGORIAN)

!! Some times are supplied as number of seconds since some reference
!! date.  To support that, set a base/reference time and then add
!! the number of seconds to it.  (time_types support adding two time
!! types together, or adding a scalar to a time_type.)
!! Here is an example of setting a reference date, giving the
!! call: year, month, day, hours, mins, secs
!ref_day0 = set_date(1970, 1, 1, 0, 0, 0)

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
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! this code here assumes an input line has: a type (1/2), location, time, value, obs error
   ! **this is where you will need to adapt this code for your input data values**

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
   ! temperature measured at a location and a vertical height, 
   ! and if otype=2, a wind speed and direction in degrees (0-360)
   ! measured at a location and a vertical pressure.

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

   ! check the lat/lon values to see if they are ok.  we require
   ! longitudes to be 0 to 360.
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle obsloop


   ! if lon comes in between -180 and 180, use these lines instead:
   !if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
   !if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   !! if time is given in seconds since some date, here's how to add it.
   !time_obs = ref_day0 + time_obs

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! this example assumes there is an obs type, where otype=1 is
   ! a temperature measured in height, and if otype=2, a wind speed
   ! and direction with a vertical value in pressure.  any observation
   ! can use any of the vertical types; this is just an example.

   if (otype == 1) then

      ! height is in meters (gph)

      ! make an obs derived type, and then add it to the sequence
      call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
                         EVAL_TEMPERATURE, terr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      if (debug) print *, 'added temperature obs to output seq'
   else

      ! DART assimilates wind as 2 separate U and V components.  assimilating
      ! "direction and speed" is difficult because direction is measured in 
      ! cyclic coordinates so you can't do simple statistics to get a mean value.

      ! convert a wind speed & direction into the U and V components
      ! and create 2 obs for it.  assume vert is in mb or hectopascals,
      ! convert to pascals.  DART assumes all pressures are in pascals.
      ! check your data source; usually wind direction is specified as the 
      ! direction the wind is coming from, increasing degrees going in a 
      ! clockwise circle.  U and V wind components have the opposite sign.
      uwnd = sin(wdir * DEG2RAD) * wspeed * (-1.0_r8)
      vwnd = cos(wdir * DEG2RAD) * wspeed * (-1.0_r8)
      uerr = werr
      verr = werr

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
! inputs:
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    otype - observation type
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
! outputs:
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, otype, oerr, day, sec, qc, obs)
use        types_mod, only : r8
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 integer,        intent(in)    :: otype, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_type_of_obs(obs_def, otype)
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
!     prev_time - the time of the previously added observation 
!                (will also be updated by this routine)
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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
