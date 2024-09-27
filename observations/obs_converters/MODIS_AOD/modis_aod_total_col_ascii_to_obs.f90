! Copyright 2019 University Corporation for Atmospheric Research and 
! Colorado Department of Public Health and Environment.
!
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use 
! this file except in compliance with the License. You may obtain a copy of the 
! License at      http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
! CONDITIONS OF ANY KIND, either express or implied. See the License for the 
! specific language governing permissions and limitations under the License.
!
! Development of this code utilized the RMACC Summit supercomputer, which is 
! supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236),
! the University of Colorado Boulder, and Colorado State University.
! The Summit supercomputer is a joint effort of the University of Colorado Boulder
! and Colorado State University.

program modis_aod_total_col_ascii_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   modis_aod_total_col_ascii_to_obs - a program that only needs minor customization to read
!      in a modis_ascii-based dataset - either white-space separated values or
!      fixed-width column data.
!      
!     this is work in progress. MODIS dataset are in HDF format. I do not
!     have HDF libraries for now, so i read the hdf file in matlab and
!     did some processing before i dumped the data in ascii. what you are
!     reading here is a 'processed' dataset of MODIS AOD
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!     modified 23 Nov 2011  ave arellano UArizona
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use types_mod, only : r8, PI, DEG2RAD
   use utilities_mod, only : initialize_utilities, finalize_utilities, &
       open_file, close_file, find_namelist_in_file, &
       check_namelist_read

   use time_manager_mod, only : time_type, set_calendar_type, set_date, &
       operator(>=), increment_time, get_time, &
       operator(-), GREGORIAN, operator(+), print_date
   use location_mod, only : VERTISUNDEF
   use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
       static_init_obs_sequence, init_obs, write_obs_seq, & 
       init_obs_sequence, get_num_obs, & 
       set_copy_meta_data, set_qc_meta_data
   use obs_kind_mod, only : MODIS_AOD_TOTAL_COL

   implicit none

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'modis_aod_total_col_ascii_to_obs.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''

   character(len=64), parameter :: modis_ascii_input_file = 'modis_asciidata.input'
   character(len=64), parameter :: obs_out_file    = 'modis_aod_total_col_obs_seq.out'

   logical, parameter :: debug = .false.  ! set to .true. to print info

   character (len=84) :: input_line

   integer :: n, i, oday, osec, rcio, iunit, otype
   integer :: year, month, day, hour, minute, second
   integer :: num_copies, num_qc, max_obs
           
   logical  :: file_exist, first_obs

   real(r8) :: aod, aoderr, qc
   real(r8) :: lat, lon, vert

   type(obs_sequence_type) :: obs_seq
   type(obs_type)          :: obs, prev_obs
   type(time_type)         :: comp_day0, time_obs, prev_time
!
   integer                 :: beg_year,beg_mon,beg_day,beg_hour,beg_min,beg_sec 
   integer                 :: end_year,end_mon,end_day,end_hour,end_min,end_sec
   integer                 :: io
   real                    :: fac,fac_err,fac_obs_error
   real                    :: lat_mn,lat_mx,lon_mn,lon_mx
   character(len=180)      :: file_in
   logical                 :: use_log_aod
!
   namelist /create_modis_obs_nml/beg_year,beg_mon,beg_day, &
   beg_hour,beg_min,beg_sec,end_year,end_mon,end_day,end_hour,end_min,end_sec, &
   fac_obs_error,file_in,lat_mn,lat_mx,lon_mn,lon_mx,use_log_aod
   
! start of executable code
   call initialize_utilities('modis_aod_total_col_ascii_to_obs')

! initialize some namelist variables
   fac=1.0
   fac_err=1.0
   fac_obs_error=1.0
   
! time setup
   call set_calendar_type(GREGORIAN)

!! some times are supplied as number of seconds since some reference
!! date.  This is an example of how to support that.
!! put the reference date into DART format
!comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

! open and read namelist
! Read the namelist entry
   call find_namelist_in_file("create_modis_obs_nml.nl", "create_modis_obs_nml", iunit)
   read(iunit, nml = create_modis_obs_nml, iostat = io)
   call check_namelist_read(iunit, io, "create_modis_obs_nml")

!   print *, 'beg_year        ',beg_year
!   print *, 'beg_mon         ',beg_mon
!   print *, 'beg_day         ',beg_day
!   print *, 'beg_hour        ',beg_hour
!   print *, 'beg_min         ',beg_min
!   print *, 'beg_sec         ',beg_sec
!   print *, 'end_year        ',end_year
!   print *, 'end_mon         ',end_mon
!   print *, 'end_day         ',end_day
!   print *, 'end_min         ',end_min
!   print *, 'end_sec         ',end_sec
!   print *, 'file_in         ',trim(file_in)
!   print *, 'lat_mn          ',lat_mn
!   print *, 'lat_mx          ',lat_mx
!   print *, 'lon_mn          ',lon_mn
!   print *, 'lon_mx          ',lon_mx
!   print *, ' '
   
! open input modis_ascii file

   iunit = open_file(modis_ascii_input_file, 'formatted', 'read')
   if (debug) print *, 'opened input file ' // trim(modis_ascii_input_file)

! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
   max_obs    = 10000000
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
   call set_copy_meta_data(obs_seq, 1, 'MODIS observation')
   call set_qc_meta_data(obs_seq, 1, 'MODIS QC index')

! if you want to append to existing files (e.g. you have a lot of
! small modis_ascii files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files 
! once they are in DART obs_seq format.

! existing file found, append to it
!   inquire(file=obs_out_file, exist=file_exist)
!   if ( file_exist ) then
!      call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
!   endif

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
   qc = 0.0_r8

   obsloop: do    ! no end limit - have the loop break when input ends

! read in a line from the modis_ascii file.   What you need to create an obs:
!  location: lat, lon, and height in pressure or meters
!  time: when the observation was taken
!  type: from the DART list of obs types
!  error: very important - the instrument error plus representativeness error
!        (see html file for more info)

! assume here a line is a type (1/2), location, time, value, obs error

! read in entire modis_ascii line into a buffer
      read(iunit, "(A)", iostat=rcio) input_line
!      print *, input_line
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
   
! here, otype is fixed to 1 for MODIS AOD
      if (debug) print *, 'next observation type = ', otype

! lon = [-180.,180.]
! lon = [-90.,900.]
      read(input_line(3:84), *, iostat=rcio) lat, lon, vert, &
                              year, month, day, hour, minute, second, &
                              aod, aoderr
      if(lon.lt.0.) lon=lon+360.
!
! APM: adjust obs error for tuning purposes
      aoderr = aoderr * fac_obs_error
!
! temp correction (July 1 comes in as July 0)
      if(month.eq.7 .and. day.eq.0) then
         day = day + 1
      endif
! play with the error for now
      if (debug) print *, 'next observation located at lat, lon = ', lat, lon
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code getting rest of aod obs, rcio = ', rcio
         exit obsloop
      endif
   
      if (debug) print *, 'next observation located at lat, lon = ', lat, lon

! check the lat/lon values to see if they are ok
      if ( lon < 0.0_r8 .or. lon > 360.0_r8 ) cycle obsloop
      if ( lat < -90.0_r8 .or. lat > 90.0_r8 ) cycle obsloop
      if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360

! put date into a dart time format
      time_obs = set_date(year, month, day, hour, minute, second)

      if (debug) call print_date(time_obs, 'next obs time is')

!! if time is given in seconds since 1/1/1970, here's how to add it.
!      time_obs = comp_day0 + time_obs

! extract time of observation into gregorian day, sec.
      call get_time(time_obs, osec, oday)

! this example assumes there is an obs type, where otype=1 is
! a temperature measured in height, and if otype=2, there's a wind
! speed and direction and height is pressure.  any kind of observation
! can use any of the vertical types; this is just an example.

! fixed to otype=1 for MODIS AOD
! no height since it is a column integrated quantity

! make an obs derived type, and then add it to the sequence
      call create_3d_obs(lat, lon, vert, VERTISUNDEF, aod, &
                      MODIS_AOD_TOTAL_COL, aoderr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      if (debug) print *, 'added aod obs to output seq'
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
   use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
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
  call set_obs_def_type_of_obs(obs_def, okind)
  call set_obs_def_time(obs_def, set_time(sec, day))
  call set_obs_def_error_variance(obs_def, oerr * oerr)
  call set_obs_def(obs, obs_def)

  obs_val(1) = obsv
  call set_obs_values(obs, obs_val)
  qc_val(1)  = qc
  call set_qc(obs, qc_val)

end subroutine create_3d_obs
!
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
!
end program modis_aod_total_col_ascii_to_obs
