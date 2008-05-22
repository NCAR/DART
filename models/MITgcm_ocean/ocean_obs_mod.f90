! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module ocean_obs_mod

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/ncep_obs/ocean_obs_mod.f90 $
! $Id: ocean_obs_mod.f90 3349 2008-05-16 19:35:31Z thoar $
! $Revision: 3349 $
! $Date: 2008-05-16 13:35:31 -0600 (Fri, 16 May 2008) $

use types_mod,        only : r8, rad2deg, PI
use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use time_manager_mod, only : time_type, get_date, set_time, GREGORIAN, &
                             set_date, set_calendar_type, get_time, &
                             operator(==)
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             E_ERR, E_MSG, timestamp
use     location_mod, only : location_type, set_location, VERTISHEIGHT, VERTISSURFACE
use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, &
                             copy_obs, set_copy_meta_data, set_qc_meta_data, set_obs_def

use     obs_kind_mod, only : KIND_SALINITY, &
                             KIND_TEMPERATURE, &
                             KIND_U_WIND_COMPONENT, &
                             KIND_V_WIND_COMPONENT, &
                             KIND_SEA_SURFACE_HEIGHT

use obs_kind_mod, only : SALINITY
use obs_kind_mod, only : TEMPERATURE
use obs_kind_mod, only : U_CURRENT_COMPONENT
use obs_kind_mod, only : V_CURRENT_COMPONENT
use obs_kind_mod, only : SEA_SURFACE_HEIGHT

implicit none
private

public :: real_obs_sequence

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: http://subversion.ucar.edu/DAReS/DART/trunk/ncep_obs/ocean_obs_mod.f90 $", &
   revision = "$Revision: 3349 $", &
   revdate  = "$Date: 2008-05-16 13:35:31 -0600 (Fri, 16 May 2008) $"

logical, save :: module_initialized = .false.

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

contains

!-------------------------------------------------

function real_obs_sequence (obsfile, year, month, day, max_num, &
          ADCP, DRIFTERS, GLIDERS, TMI, SSH, &
          lon1, lon2, lat1, lat2)
!------------------------------------------------------------------------------
!  this function is to prepare data to DART sequence format
!
character(len=129), intent(in) :: obsfile
integer,            intent(in) :: year, month, day, max_num
logical,            intent(in) :: ADCP, DRIFTERS, GLIDERS, TMI, SSH
real(r8),           intent(in) :: lon1, lon2, lat1, lat2

type(obs_sequence_type) :: real_obs_sequence


type(obs_type) :: obs, prev_obs
integer :: i, num_copies, num_qc
integer :: days, seconds
integer :: yy, mn, dd, hh, mm, ss
integer :: startdate1, startdate2
integer :: obs_num, calender_type, iskip
integer :: obs_unit
integer :: obs_kind, obs_kind_gen, which_vert, obstype

real (r8) :: lon, lat, vloc, obs_value
real (r8) :: aqc, var2, lonc
type(time_type) :: time, pre_time

character(len = 8 ) :: obsdate
character(len = 80) :: label
character(len = 129) :: copy_meta_data, qc_meta_data

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

obs_unit  = get_unit()
open(unit = obs_unit, file = obsfile, form='formatted', status='old')
print*, 'input file opened= ', trim(obsfile)
rewind (obs_unit)

!print*, 'ncep obsdates = ', obsdate

obs_num = 0
iskip   = 0

!  loop over all observations within the file
!------------------------------------------------------------------------------

obsloop:  do

   read(obs_unit,*,end=200) lon, lat, vloc, obs_value, which_vert, var2, aqc, &
                              obstype, startdate1, startdate2

   ! Calculate the DART time from the observation time 
   yy =     startdate1/10000
   mn = mod(startdate1/100,100)
   dd = mod(startdate1    ,100)
   hh =     startdate2/10000
   mm = mod(startdate2/100,100)
   ss = mod(startdate2    ,100)
   time = set_date(yy,mn,dd,hh,mm,ss)
   call get_time(time,seconds,days)

   ! verify the location is not outside valid limits
   if((lon > 360.0_r8) .or. (lon <   0.0_r8) .or.  &
      (lat >  90.0_r8) .or. (lat < -90.0_r8)) then
      write(*,*) 'invalid location.  lon,lat = ', lon, lat
      iskip = iskip + 1
      cycle obsloop
   endif

   lonc = lon
   if (lon2 > 360.0_r8 .and. lon < 180.0_r8) lonc = lon + 360.0_r8

   ! reject observations outside the bounding box
   if(lat < lat1 .or. lat > lat2 .or. lonc < lon1 .or. lonc > lon2) then
     iskip = iskip + 1
     cycle obsloop
   endif

!   assign each observation the correct observation type
   if(obstype == 1) then
     obs_kind_gen = KIND_TEMPERATURE
     obs_kind     =      TEMPERATURE
   elseif(obstype == 2) then
    obs_kind_gen = KIND_SALINITY
    obs_kind     =      SALINITY
   elseif(obstype == 3) then
    obs_kind_gen = KIND_U_WIND_COMPONENT
    obs_kind     =      U_CURRENT_COMPONENT
   elseif(obstype == 4) then
    obs_kind_gen = KIND_V_WIND_COMPONENT
    obs_kind     =      V_CURRENT_COMPONENT
   elseif(obstype == 5) then
    obs_kind_gen = KIND_SEA_SURFACE_HEIGHT
    obs_kind     =      SEA_SURFACE_HEIGHT
   else
      print*, 'unknown observation type ... skipping ...'
      cycle obsloop
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
                 var2, aqc, obs_kind, which_vert, seconds, days)
   
   if(obs_num == 1) then ! for the first observation 

     call insert_obs_in_seq(real_obs_sequence, obs)
     call copy_obs(prev_obs, obs)
     pre_time = time

   else  !  not the first observation 

     if(time == pre_time) then  ! same time as previous observation

       call insert_obs_in_seq(real_obs_sequence, obs, prev_obs)
       call copy_obs(prev_obs, obs)
       pre_time = time

    else  ! not the same time

       call insert_obs_in_seq(real_obs_sequence, obs)
       call copy_obs(prev_obs, obs)
       pre_time = time

    endif

  endif

end do obsloop

200 continue

close(obs_unit)

print*, 'date ', obsdate, ' obs used = ', obs_num, ' obs skipped = ', iskip
!print*, 'obs_num= ',obs_num,' skipped= ',obsdate,iskip

end function real_obs_sequence



subroutine real_obs(num_copies, num_qc, obs, lon, lat, vloc, obs_value, &
                      var2, aqc, obs_kind, which_vert, seconds, days)
!------------------------------------------------------------------------------
integer,        intent(in)    :: num_copies, num_qc
type(obs_type), intent(inout) :: obs
real(r8),       intent(in)    :: lon, lat, vloc, obs_value, var2, aqc
integer,        intent(in)    :: obs_kind, which_vert, seconds, days

integer            :: i
real(r8)           :: aqc01(1), obs_value01(1)
type(obs_def_type) :: obsdef0

if ( .not. module_initialized ) call initialize_module

! Does real initialization of an observation type

call real_obs_def(obsdef0, lon, lat, vloc, &
                    var2, obs_kind, which_vert, seconds, days)
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
call set_obs_def_kind(obs_def, obs_kind)

call set_obs_def_time(obs_def, set_time(seconds, days) )
call set_obs_def_error_variance(obs_def, var2)

end subroutine real_obs_def



subroutine initialize_module
!-------------------------------------------------
call register_module(source, revision, revdate)
module_initialized = .true.
end subroutine initialize_module


end module ocean_obs_mod
