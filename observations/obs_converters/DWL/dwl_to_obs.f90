! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dwl_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   dwl_to_obs - read in ascii lines with doppler wind lidar data.
!      see below for exact input format.  create 2 observations from
!      each input line - U and V components of the wind.
!
!   3 Sep 2013, nancy collins NCAR, Matic Savli FMF/Univ Ljubljana
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD
use      utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                               check_namelist_read, nmlfileunit, do_nml_file,  &
                               get_next_filename, error_handler, E_ERR, E_MSG, &
                               do_nml_term, finalize_utilities,                &
                               open_file, close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISHEIGHT, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use      obs_kind_mod, only : DWL_U_WIND_COMPONENT, DWL_V_WIND_COMPONENT

implicit none

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

! the max possible number of obs needs to be specified but it will only 
! write out the actual number created.

character(len=128)  :: text_input_file = 'dwldata.input'  ! default input name
character(len=128)  :: obs_out_file    = 'obs_seq.out'    ! default output name
integer             :: max_obs         = 100000   !  max number of obs in one file
logical             :: add_obs_data    = .true.   ! .false. makes empty observations
logical             :: debug           = .false.  ! .true. prints more info

namelist /dwl_to_obs_nml/ max_obs, text_input_file, obs_out_file, add_obs_data, debug


! local variables

character (len=129) :: input_line
character (len=20)  :: date_string

integer :: oday, osec, rcio, iunit, otype
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, lcount
           
logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc, wdir, wspeed, werr
real(r8) :: lat, lon, vert, uwnd, uerr, vwnd, verr

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

! start of executable code

call initialize_utilities('dwl_to_obs')

! namelist handling
call find_namelist_in_file("input.nml", "dwl_to_obs_nml", iunit)
read(iunit, nml = dwl_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "dwl_to_obs_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=dwl_to_obs_nml)
if (do_nml_term()) write(     *     , nml=dwl_to_obs_nml)


! time setup
call set_calendar_type(GREGORIAN)

! open input text file

iunit = open_file(text_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(text_input_file)


if (add_obs_data) then
   print *, 'adding observation values from text file to output'

   ! each observation in this series will have a single observation value 
   ! and a quality control flag. 
   num_copies = 1
   num_qc     = 1
else
   print *, 'creating no-data observations at the requested locations'

   num_copies = 0
   num_qc     = 0
endif

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

if (add_obs_data) then
   ! the first one needs to contain the string 'observation' and the
   ! second needs the string 'QC'.
   call set_copy_meta_data(obs_seq, 1, 'observation')
   call set_qc_meta_data(obs_seq, 1, 'Data QC')
endif

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

! this is just for debugging or errors - count the number of input lines
lcount = 0

obsloop: do    ! no end limit - have the loop break when input ends

   lcount = lcount + 1

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! assume here a line is: a location, time, 2 obs values, obs error

   ! read in entire text line into a buffer.  exit when no more lines.
   read(iunit, "(A)", iostat=rcio) input_line
   if (rcio /= 0) then 
      if (debug) print *, 'line number ', lcount
      if (debug) print *, 'got bad read code from input file, rcio = ', rcio
      if (debug) print *, '(bad read expected if end of file)'
      exit obsloop
   endif

   ! extract the different values from the input line
   read(input_line, *, iostat=rcio) lon, lat, vert, date_string, &
                                    uwnd, vwnd, werr
   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code getting next wind obs, rcio = ', rcio
      if (debug) print *, 'line number ', lcount, ' input line was:'
      if (debug) print *, trim(input_line)
      exit obsloop
   endif
   
   if (debug) print *, 'next observation located at lon, lat = ', lon, lat
   
   ! date format is: dd-mm-ccyy_hh:nn:ss
   read(date_string, "(2(I2,1X),I4,3(1X,I2))", iostat=rcio) day, month, year, hour, minute, second
   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code getting next time value, rcio = ', rcio
      if (debug) print *, 'line number ', lcount, ' input line was:'
      if (debug) print *, trim(input_line)
      exit obsloop
   endif
   
   if (debug) print *, 'next observation is at time ', year, month, day, hour, minute, second
   if (debug) print *, 'next observation values/err ', uwnd, vwnd, werr


   ! check the lat/lon values to see if they are ok
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle obsloop


   ! if lon comes in between -180 and 180, use these lines instead:
   !if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
   !if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! vertical in height in meters

   ! if you want to create different obs types for clear air vs cloud,
   ! change the DWL_x_WIND_COMPONENT below to what you need.  it must be
   ! defined in obs_def/obs_def_dwl_mod.f90, and then quickbuild.sh must
   ! be run again.

   call create_3d_obs(add_obs_data, lat, lon, vert, VERTISHEIGHT, uwnd, &
                         DWL_U_WIND_COMPONENT, werr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
   call create_3d_obs(add_obs_data, lat, lon, vert, VERTISHEIGHT, vwnd, &
                         DWL_V_WIND_COMPONENT, werr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   if (debug) print *, 'added 2 wind obs to output seq'

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
!    add_data - if .false. create the location, time, err, kind only - no data
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
subroutine create_3d_obs(add_data, lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, obs)
use        types_mod, only : r8
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 logical,        intent(in)    :: add_data
 integer,        intent(in)    :: okind, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_type_of_obs(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

if (add_data) then
   obs_val(1) = obsv
   call set_obs_values(obs, obs_val)
   qc_val(1)  = qc
   call set_qc(obs, qc_val)
endif

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

end program dwl_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
