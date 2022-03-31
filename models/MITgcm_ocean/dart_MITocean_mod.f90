! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module dart_MITocean_mod

use types_mod,        only : r8, rad2deg, PI

use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             set_obs_def_key

use time_manager_mod, only : time_type, get_date, set_time, GREGORIAN, &
                             set_date, set_calendar_type, get_time, &
                             print_date, print_time, operator(>=), &
                             operator(<), operator(>), operator(+)  

use    utilities_mod, only : open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             E_ERR, E_MSG, timestamp, is_longitude_between

use     location_mod, only : location_type, set_location, VERTISHEIGHT, VERTISSURFACE

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, &
                             copy_obs, set_copy_meta_data, set_qc_meta_data, &
                             get_first_obs, get_last_obs, get_obs_def, set_obs_def

use     obs_kind_mod, only : get_index_for_type_of_obs

use obs_def_ocean_mod, only : set_hf_radial_vel

implicit none
private

public :: real_obs_sequence

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'dart_MITocean_mod'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

contains

!-------------------------------------------------
!>  this function is to prepare data to DART sequence format

function real_obs_sequence (obsfile, time1, timeN, max_num, &
          lon1, lon2, lat1, lat2, hfradar)

character(len=*), intent(in) :: obsfile
type(time_type),  intent(in) :: time1, timeN
integer,          intent(in) :: max_num
real(r8),         intent(in) :: lon1, lon2, lat1, lat2
logical,          intent(in) :: hfradar

type(obs_sequence_type) :: real_obs_sequence

type(obs_def_type) :: obs_def
type(obs_type) :: obs, prev_obs
integer :: i, num_copies, num_qc
integer :: days, seconds
integer :: yy, mn, dd, hh, mm, ss
integer :: startdate1, startdate2
integer :: obs_num, calender_type, iskip
integer :: obs_unit, ios, id, defkey
integer :: which_vert, obstype

real (r8) :: lon, lat, vloc, obs_value
real (r8) :: aqc, var2, angle
type(time_type) :: time, pre_time

character(len=32) :: obs_kind_name
character(len=80) :: label
character(len=129) :: copy_meta_data, qc_meta_data

if ( .not. module_initialized ) call initialize_module

num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence 

call init_obs_sequence(real_obs_sequence, num_copies, num_qc, max_num)

! set meta data of obs_seq

do i = 1, num_copies
   copy_meta_data = 'observation'  
   call set_copy_meta_data(real_obs_sequence, i, copy_meta_data)
end do

do i = 1, num_qc
   qc_meta_data = 'QC index'
   call set_qc_meta_data(real_obs_sequence, i, qc_meta_data)
end do

! Initialize the obs variable

call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! set observation time type
calender_type = GREGORIAN
call set_calendar_type(calender_type)

! open observation data file

obs_unit = open_file(obsfile, form='formatted', action='read')
print*, 'input file opened= "'//trim(obsfile)//'"'
rewind (obs_unit)

obs_num = 0
iskip   = 0
defkey  = 1

!  loop over all observations within the file
!------------------------------------------------------------------------------

obsloop:  do

   if (hfradar) then
      read(obs_unit,*,iostat=ios,end=200) lon, lat, vloc, angle, obs_value, &
              which_vert, var2, aqc, obs_kind_name, startdate1, startdate2, id
   else
      read(obs_unit,*,iostat=ios,end=200) lon, lat, vloc, obs_value, &
              which_vert, var2, aqc, obs_kind_name, startdate1, startdate2
   endif

   if (ios /= 0) then
      print*,''
      print*,' Observation ', obs_num+1
      print*,' lon lat vloc obs_value ',lon, lat, vloc, obs_value
      print*,' which_vert var2 aqc ',which_vert, var2, aqc
      print*,' obs_kind_name ',obs_kind_name
      print*,' date1 date2 ',startdate1, startdate2
      if (hfradar) print*, 'angle id ',angle, id
   endif

   ! Calculate the DART time from the observation time 
   yy =     startdate1/10000
   mn = mod(startdate1/100,100)
   dd = mod(startdate1    ,100)
   hh =     startdate2/10000
   mm = mod(startdate2/100,100)
   ss = mod(startdate2    ,100)
   time = set_date(yy,mn,dd,hh,mm,ss)
   call get_time(time,seconds,days)

   ! reject observations outside the time window
   if ( (time < time1) .or. (time > timeN) ) then
      call print_date(time, ' skipping time outside desired window ... ')
      iskip = iskip + 1
      cycle obsloop
   endif

   ! verify the latitude is not outside valid limits
   if ((lat >  90.0_r8) .or. (lat < -90.0_r8)) then
      write(*,*) 'invalid location.  lon,lat = ', lon, lat
      iskip = iskip + 1
      cycle obsloop
   endif

   ! reject observations outside the bounding box
   if(lat < lat1 .or. lat > lat2 .or. &
      .not. is_longitude_between(lon, lon1, lon2)) then
      iskip = iskip + 1
      cycle obsloop
   endif

   ! and now make sure lon is between 0 and 360, so the
   ! dart locations module is happy.  using modulo is important 
   ! here; plain mod() does not handle negative numbers the same.
   lon = modulo(lon, 360.0_r8)


   ! assign each observation the correct observation type
   obstype = get_index_for_type_of_obs(obs_kind_name)
   if(obstype < 1) then
      print*, 'unknown observation type [',trim(obs_kind_name),'] ... skipping ...'
      cycle obsloop
  !else
  !   print*,trim(obs_kind_name),' is ',obstype
   endif

   obs_num = obs_num + 1

   ! print a reassuring message after every Nth processed obs.
   ! if requested, print in the form of a timestamp.  
   ! the default is just a plain string with the current obs count.
   if(mod(obs_num, print_every_Nth) == 0) then
       write(label, *) 'obs count = ', obs_num
       if (print_timestamps) then
          call timestamp(string1=label, pos='')
       else
          write(*,*) trim(label)
       endif
   endif
   if(obs_num == max_num) then
      print*, 'Max limit for observation count reached.  Increase value in namelist'
      stop
   endif
   
!   create the obs_def for this observation, add to sequence
!------------------------------------------------------------------------------
   
   call real_obs(num_copies, num_qc, obs, lon, lat, vloc, obs_value, &
                 var2, aqc, obstype, which_vert, seconds, days,      &
                 hfradar, defkey, id, angle)
   
   if(obs_num == 1) then ! for the first observation 

      call insert_obs_in_seq(real_obs_sequence, obs)

   else  !  not the first observation 

      if(time >= pre_time) then  ! same time or later than prev

         call insert_obs_in_seq(real_obs_sequence, obs, prev_obs)

      else  ! earlier than prev obs, must search from start

         call insert_obs_in_seq(real_obs_sequence, obs)

      endif

   endif

   call copy_obs(prev_obs, obs)
   pre_time = time

end do obsloop

200 continue

close(obs_unit)

! Print a little summary
print*, 'obs used = ', obs_num, ' obs skipped = ', iskip

if ( get_first_obs(real_obs_sequence, obs) ) then
   call get_obs_def(obs, obs_def)
   pre_time = get_obs_def_time(obs_def)
   call print_time(pre_time,' first time in sequence is ')
   call print_date(pre_time,' first date in sequence is ')
endif
if( get_last_obs(real_obs_sequence, obs)) then
   call get_obs_def(obs, obs_def)
   time = get_obs_def_time(obs_def)
   call print_time(time,' last  time in sequence is ')
   call print_date(time,' last  date in sequence is ')
endif
print*, ''

end function real_obs_sequence



subroutine real_obs(num_copies, num_qc, obs, lon, lat, vloc, obs_value, &
                      var2, aqc, obs_kind, which_vert, seconds, days,   &
                      hfradar, defkey, id, angle)
!------------------------------------------------------------------------------
integer,        intent(in)    :: num_copies, num_qc
type(obs_type), intent(inout) :: obs
real(r8),       intent(in)    :: lon, lat, vloc, obs_value, var2, aqc, angle
integer,        intent(in)    :: obs_kind, which_vert, seconds, days, id
logical,        intent(in)    :: hfradar
integer,        intent(inout) :: defkey

integer            :: i
real(r8)           :: aqc01(1), obs_value01(1)
type(obs_def_type) :: obsdef0

if ( .not. module_initialized ) call initialize_module

! Does real initialization of an observation type

call real_obs_def(obsdef0, lon, lat, vloc, &
                    var2, obs_kind, which_vert, seconds, days)
if (hfradar) then
   ! this routine increments defkey and stores the private info in
   ! the ocean module until time to write.
   call set_hf_radial_vel(defkey, id, angle)
   call set_obs_def_key(obsdef0, defkey)
endif

call set_obs_def(obs, obsdef0)

do i = 1, num_copies
   obs_value01(1) = obs_value
   call set_obs_values(obs, obs_value01(1:1) )
end do

do i = 1, num_qc
   aqc01(1) = aqc
   call set_qc(obs, aqc01(1:1))
end do

end subroutine real_obs



subroutine real_obs_def(obs_def, lon, lat, vloc, &
                        var2, obs_kind, which_vert, seconds, days)
!----------------------------------------------------------------------
type(obs_def_type), intent(inout) :: obs_def
real(r8),intent(in) :: lon, lat, vloc, var2
integer, intent(in) :: obs_kind, which_vert, seconds, days

type(location_type) :: loc0

if ( .not. module_initialized ) call initialize_module

! set obs location
loc0 = set_location(lon, lat, vloc, which_vert )
call set_obs_def_location(obs_def, loc0)

! set obs kind
call set_obs_def_type_of_obs(obs_def, obs_kind)

call set_obs_def_time(obs_def, set_time(seconds, days) )
call set_obs_def_error_variance(obs_def, var2)

end subroutine real_obs_def



subroutine initialize_module
!-------------------------------------------------
call register_module(source, revision, revdate)
module_initialized = .true.
end subroutine initialize_module


end module dart_MITocean_mod

