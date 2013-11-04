! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module real_obs_mod

use types_mod,        only : r8, rad2deg, PI
use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use time_manager_mod, only : time_type, operator(>), operator(<), operator(>=), &
                             operator(/=), set_date, set_calendar_type, get_time, &
                             get_date, set_time, GREGORIAN, increment_time
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             E_ERR, E_MSG, timestamp, is_longitude_between
use     location_mod, only : location_type, set_location, &
                             VERTISPRESSURE, VERTISSURFACE
use obs_sequence_mod, only : init_obs_sequence, init_obs, obs_sequence_type, obs_type, &
                             set_copy_meta_data, set_qc_meta_data

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, KIND_SURFACE_PRESSURE, &
                             KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, KIND_RELATIVE_HUMIDITY, &
                             KIND_DEWPOINT, KIND_PRESSURE, KIND_VERTICAL_VELOCITY, &
                             KIND_RAINWATER_MIXING_RATIO, KIND_DENSITY, KIND_VELOCITY, &
                             KIND_1D_INTEGRAL, KIND_RADAR_REFLECTIVITY 

use  obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use  obs_def_altimeter_mod, only: compute_altimeter

use obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT
use obs_kind_mod, only : RADIOSONDE_V_WIND_COMPONENT
use obs_kind_mod, only : RADIOSONDE_SURFACE_PRESSURE
use obs_kind_mod, only : RADIOSONDE_TEMPERATURE
use obs_kind_mod, only : RADIOSONDE_SPECIFIC_HUMIDITY
use obs_kind_mod, only : RADIOSONDE_RELATIVE_HUMIDITY
use obs_kind_mod, only : RADIOSONDE_DEWPOINT
use obs_kind_mod, only : AIRCRAFT_U_WIND_COMPONENT
use obs_kind_mod, only : AIRCRAFT_V_WIND_COMPONENT
use obs_kind_mod, only : AIRCRAFT_TEMPERATURE
use obs_kind_mod, only : AIRCRAFT_SPECIFIC_HUMIDITY
use obs_kind_mod, only : AIRCRAFT_RELATIVE_HUMIDITY
use obs_kind_mod, only : AIRCRAFT_DEWPOINT
use obs_kind_mod, only : ACARS_U_WIND_COMPONENT
use obs_kind_mod, only : ACARS_V_WIND_COMPONENT
use obs_kind_mod, only : ACARS_TEMPERATURE
use obs_kind_mod, only : ACARS_SPECIFIC_HUMIDITY
use obs_kind_mod, only : ACARS_RELATIVE_HUMIDITY
use obs_kind_mod, only : ACARS_DEWPOINT
use obs_kind_mod, only : MARINE_SFC_U_WIND_COMPONENT
use obs_kind_mod, only : MARINE_SFC_V_WIND_COMPONENT
use obs_kind_mod, only : MARINE_SFC_TEMPERATURE
use obs_kind_mod, only : MARINE_SFC_SPECIFIC_HUMIDITY
use obs_kind_mod, only : MARINE_SFC_RELATIVE_HUMIDITY
use obs_kind_mod, only : MARINE_SFC_DEWPOINT
use obs_kind_mod, only : MARINE_SFC_ALTIMETER
use obs_kind_mod, only : LAND_SFC_U_WIND_COMPONENT
use obs_kind_mod, only : LAND_SFC_V_WIND_COMPONENT
use obs_kind_mod, only : LAND_SFC_TEMPERATURE
use obs_kind_mod, only : LAND_SFC_SPECIFIC_HUMIDITY
use obs_kind_mod, only : LAND_SFC_RELATIVE_HUMIDITY
use obs_kind_mod, only : LAND_SFC_DEWPOINT
use obs_kind_mod, only : LAND_SFC_ALTIMETER
use obs_kind_mod, only : RADIOSONDE_SURFACE_ALTIMETER
use obs_kind_mod, only : SAT_U_WIND_COMPONENT
use obs_kind_mod, only : SAT_V_WIND_COMPONENT
use obs_kind_mod, only : ATOV_TEMPERATURE


implicit none
private

public :: real_obs_sequence

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


logical, save :: module_initialized = .false.
! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000
!-------------------------------------------------

contains


function real_obs_sequence (year, month, day, hourt, max_num, select_obs, &
          ObsBase, ADDUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND, &
          obs_U, obs_V, obs_T, obs_PS, obs_QV, inc_specific_humidity,      &
          inc_relative_humidity, inc_dewpoint, bin_beg, bin_end,           &
          lon1, lon2, lat1, lat2, obs_time)
!------------------------------------------------------------------------------
!  this function is to prepare NCEP decoded BUFR data to DART sequence format
!
integer,            intent(in) :: year, month, day, max_num, select_obs
character(len = *), intent(in) :: ObsBase, hourt
logical,            intent(in) :: ADDUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND
logical,            intent(in) :: obs_U, obs_V, obs_T, obs_PS, obs_QV, obs_time
logical,            intent(in) :: inc_specific_humidity, inc_relative_humidity, &
                                  inc_dewpoint
real(r8),           intent(in) :: lon1, lon2, lat1, lat2

type(obs_sequence_type) :: real_obs_sequence


type(obs_type) :: obs, prev_obs
integer :: i, num_copies, num_qc
integer :: days, seconds
integer :: day0, sec0
integer :: hour, imin, sec
integer :: obs_num, calender_type
type(time_type) :: current_day, time_obs, prev_time

integer, parameter :: num_fail_kinds = 6
integer :: iskip(num_fail_kinds)
integer, parameter :: fail_3Z        = 1
integer, parameter :: fail_timerange = 2
integer, parameter :: fail_badloc    = 3
integer, parameter :: fail_areabox   = 4
integer, parameter :: fail_badkind   = 5
integer, parameter :: fail_notwanted = 6
character(len=32) :: skip_reasons(num_fail_kinds) = (/ &
                     'time too early (exactly 03Z)    ', &
                     'time outside bin range          ', &
                     'bad observation location        ', &
                     'observation outside lat/lon box ', &
                     'unrecognized observation kind   ', &
                     'obs type not on select list     ' /)


integer :: obs_unit
integer :: obs_prof, obs_kind, obs_kind_gen, which_vert, iqc, obstype, pc
real (r8) :: obs_err, lon, lat, lev, zob, time, rcount, zob2
real (r8) :: vloc, obs_value, aqc

real (r8) :: bin_beg, bin_end

character(len = 8 ) :: obsdate
character(len = 80) :: obsfile, label
character(len = 6 ) :: subset
logical :: pass, first_obs

if ( .not. module_initialized ) call initialize_module

num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence 

call init_obs_sequence(real_obs_sequence, num_copies, num_qc, max_num)

! set meta data of obs_seq

do i = 1, num_copies
   call set_copy_meta_data(real_obs_sequence, i, 'NCEP BUFR observation')
end do

do i = 1, num_qc
   call set_qc_meta_data(real_obs_sequence, i, 'NCEP QC index')
end do

! Initialize the obs variable

call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! set observation time type
calender_type = GREGORIAN
call set_calendar_type(calender_type)

! get the time (seconds & days) for 00Z of current day (year,month,day)
hour = 0
imin = 0
sec  = 0
current_day = set_date(year, month, day, hour, imin, sec)
call get_time(current_day, sec0, day0)

!   output the day and sec.
print*, 'processing data for day, sec= ', day0, sec0

! open NCEP observation data file

write(obsdate, '(i4.4,i2.2,i2.2)') year, month, day
obsfile   = trim(adjustl(ObsBase))//obsdate//hourt
obs_unit  = get_unit()
open(unit = obs_unit, file = obsfile, form='formatted', status='old')
print*, 'input file opened= ', trim(obsfile)
rewind (obs_unit)

!print*, 'ncep obsdates = ', obsdate

obs_num = 0
iskip(:) = 0
first_obs = .true.

!  loop over all observations within the file
!------------------------------------------------------------------------------

obsloop:  do

   read(obs_unit,880,end=200) obs_err, lon, lat, lev, zob, zob2, rcount, time, &
                              obstype, iqc, subset, pc

 880 format(f4.2,2f9.4,e12.5,f7.2,f7.2,f9.0,f7.3,i4,i2,1x,a6,i2)

!   A 'day' is from 03:01Z of one day through 03Z of the next.
!   skip the observations at exact 03Z of the beginning of the day
!   (obs at 03Z the next day have a time of 27.)
!------------------------------------------------------------------------------
   if(time == 3.0_r8) then
      iskip(fail_3Z) = iskip(fail_3Z) + 1
      cycle obsloop 
   endif 

   !  select the obs for the time window
   if(time < bin_beg .or. time > bin_end) then
      iskip(fail_timerange) = iskip(fail_timerange) + 1
      cycle obsloop
   endif

   ! verify the location is not outside valid limits
   if((lon > 360.0_r8) .or. (lon <   0.0_r8) .or.  &
      (lat >  90.0_r8) .or. (lat < -90.0_r8)) then
      write(*,*) 'invalid location.  lon,lat = ', lon, lat
      iskip(fail_badloc) = iskip(fail_badloc) + 1
      cycle obsloop
   endif

   ! reject observations outside the bounding box
   if(lat < lat1 .or. lat > lat2 .or. & 
     .not. is_longitude_between(lon, lon1, lon2)) then
     iskip(fail_areabox) = iskip(fail_areabox) + 1
     cycle obsloop
   endif

   ! it seems when a machine writes out a number and then reads it back in
   ! you should at least get consistent round-off errors.  perhaps it is more
   ! complex than it seems, but on the ibm xlf compiler, values of exactly
   ! 90.0 and -90.0 are written out as (binary) doubles and when they are
   ! read back in and compared to 90 and -90, they are not equal anymore.
   ! the difference is in a very small but apparently significant digit,
   ! so we error out in our 3d sphere locations module.  this hack is a try
   ! at adjusting the numbers by a small enough value that the location is
   ! still the pole but isn't going to round out of range.
   if      (lat >=  89.9999_r8) then
     lat = lat - 1.0e-12_r8
     !print *, 'lat adjusted down, now ', lat
   else if (lat <= -89.9999_r8) then
     lat = lat + 1.0e-12_r8
     !print *, 'lat adjusted   up, now ', lat
   endif

   obs_prof = rcount/1000000


!   assign each observation the correct observation type
!------------------------------------------------------------------------------

   ! make sure we do not fall through the code below without setting
   ! a valid obs kind (e.g. the obstype is one not listed)
   obs_kind = -1

   if(obs_prof == 1) then
     obs_kind_gen = KIND_TEMPERATURE
     if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_TEMPERATURE
     if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_TEMPERATURE
     if(obstype == 133                    ) obs_kind = ACARS_TEMPERATURE
     if(obstype == 161 .or. obstype == 163) obs_kind = ATOV_TEMPERATURE
     if(obstype == 171 .or. obstype == 173) obs_kind = ATOV_TEMPERATURE
     if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_TEMPERATURE
     if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_TEMPERATURE
   endif

   if(obs_prof == 5) then
     if ( zob2 == 0.0_r8 .and. inc_specific_humidity ) then
       obs_kind_gen = KIND_SPECIFIC_HUMIDITY
       if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_SPECIFIC_HUMIDITY
       if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_SPECIFIC_HUMIDITY
       if(obstype == 133                    ) obs_kind = ACARS_SPECIFIC_HUMIDITY
       if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_SPECIFIC_HUMIDITY
       if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_SPECIFIC_HUMIDITY
     else if ( zob2 == 1.0_r8 .and. inc_relative_humidity ) then
       obs_kind_gen = KIND_RELATIVE_HUMIDITY
       if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_RELATIVE_HUMIDITY
       if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_RELATIVE_HUMIDITY
       if(obstype == 133                    ) obs_kind = ACARS_RELATIVE_HUMIDITY
       if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_RELATIVE_HUMIDITY
       if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_RELATIVE_HUMIDITY
     else if ( zob2 == 2.0_r8 .and. inc_dewpoint ) then
       obs_kind_gen = KIND_DEWPOINT
       if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_DEWPOINT
       if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_DEWPOINT
       if(obstype == 133                    ) obs_kind = ACARS_DEWPOINT
       if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_DEWPOINT
       if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_DEWPOINT
     endif
   endif

   if(obs_prof == 3) then
     obs_kind_gen = KIND_SURFACE_PRESSURE
     if(obstype == 120                    ) obs_kind = RADIOSONDE_SURFACE_ALTIMETER
     if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_ALTIMETER
     if(obstype == 181                    ) obs_kind = LAND_SFC_ALTIMETER
   endif

   if(obs_prof == 2) then
     obs_kind_gen = KIND_U_WIND_COMPONENT
     if(obstype == 220 .or. obstype == 232) obs_kind = RADIOSONDE_U_WIND_COMPONENT
     if(obstype == 221                    ) obs_kind = RADIOSONDE_U_WIND_COMPONENT
     if(obstype == 230 .or. obstype == 231) obs_kind = AIRCRAFT_U_WIND_COMPONENT
     if(obstype == 233                    ) obs_kind = ACARS_U_WIND_COMPONENT
     if(obstype == 242 .or. obstype == 243) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 245 .or. obstype == 246) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 252 .or. obstype == 253) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 255                    ) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 280 .or. obstype == 282) obs_kind = MARINE_SFC_U_WIND_COMPONENT
     if(obstype == 281 .or. obstype == 284) obs_kind = LAND_SFC_U_WIND_COMPONENT
   endif

   if(obs_prof == 9) then
     obs_kind_gen = KIND_V_WIND_COMPONENT
     if(obstype == 220 .or. obstype == 232) obs_kind = RADIOSONDE_V_WIND_COMPONENT
     if(obstype == 221                    ) obs_kind = RADIOSONDE_V_WIND_COMPONENT
     if(obstype == 230 .or. obstype == 231) obs_kind = AIRCRAFT_V_WIND_COMPONENT
     if(obstype == 233                    ) obs_kind = ACARS_V_WIND_COMPONENT
     if(obstype == 242 .or. obstype == 243) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 245 .or. obstype == 246) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 252 .or. obstype == 253) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 255                    ) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 280 .or. obstype == 282) obs_kind = MARINE_SFC_V_WIND_COMPONENT
     if(obstype == 281 .or. obstype == 284) obs_kind = LAND_SFC_V_WIND_COMPONENT
   endif

   if (obs_kind < 0) then
      ! the "real" fix if the record type is not found might actually be to
      ! accept all record types within valid ranges, and depend on the first
      ! preprocessing steps (in the prepbufr converter) to remove obs record
      ! types which are not desired.  for now, avoid giving them the wrong type
      ! and quietly loop.
      !print *, 'unrecognized obs_prof or obstype, skipping', obs_prof, obstype
      iskip(fail_badkind) = iskip(fail_badkind) + 1
      cycle obsloop 
   endif

!   check to see if this observation is desired
!------------------------------------------------------------------------------

   ! if select_obs is 0, we are going to include all observations
   ! and we skip the selection code below.
   if(select_obs /= 0) then

      ! assume we are going to ignore this obs, unless it is
      ! specifically included by one of the selections below.
      pass = .true.

      ! select the specific NCEP obs types
      if( (ADDUPA .and. (subset =='ADPUPA')) .or. &
          (AIRCAR .and. (subset =='AIRCAR')) .or. &
          (AIRCFT .and. (subset =='AIRCFT')) .or. &
          (SATEMP .and. (subset =='SATEMP')) .or. &
          (SFCSHP .and. (subset =='SFCSHP')) .or. &
          (ADPSFC .and. (subset =='ADPSFC')) .or. &
          (SATWND .and. (subset =='SATWND'))       ) then

         ! then select the obs kind requested
         if( (obs_T                 .and. (obs_kind_gen == KIND_TEMPERATURE ))      .or. &
             (obs_U                 .and. (obs_kind_gen == KIND_U_WIND_COMPONENT )) .or. &
             (obs_V                 .and. (obs_kind_gen == KIND_V_WIND_COMPONENT )) .or. &
             (obs_PS                .and. (obs_kind_gen == KIND_SURFACE_PRESSURE))  .or. &
             (obs_QV                .and. (obs_kind_gen == KIND_SPECIFIC_HUMIDITY)) .or. &
             (inc_relative_humidity .and. (obs_kind_gen == KIND_RELATIVE_HUMIDITY)) .or. &
             (inc_dewpoint          .and. (obs_kind_gen == KIND_DEWPOINT)) ) then
             pass = .false.
         endif

      endif

      ! if pass is still true, we want to ignore this obs.
      if(pass) then
         iskip(fail_notwanted) = iskip(fail_notwanted) + 1
         cycle obsloop 
      endif

   endif

!   process this observation
!------------------------------------------------------------------------------

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

   ! set vertical coordinate for upper-air observations
   if (subset == 'AIRCAR' .or. subset == 'AIRCFT' .or. &
       subset == 'SATEMP' .or.                         &
       subset == 'SATWND' .or. subset == 'ADPUPA' ) then
       vloc = lev*100.0_r8          ! convert from mb to Pascal
       which_vert = VERTISPRESSURE
   endif

   ! set vertical coordinate for surface observations
   if (subset == 'ADPSFC' .or. subset == 'SFCSHP') then
     vloc = lev
     ! some obs have elevation of 1e12 - toss those.
     if ( subset == 'SFCSHP' .and. vloc > 4000.0_r8) then
        iskip = iskip + 1
        cycle obsloop
     endif
     which_vert = VERTISSURFACE
   endif

   ! set obs value and error if necessary
   if ( obs_kind == LAND_SFC_ALTIMETER .or. obs_kind == MARINE_SFC_ALTIMETER &
        .or. obs_kind == RADIOSONDE_SURFACE_ALTIMETER ) then
      vloc = lev                  ! station height, not used now for Ps obs
      which_vert = VERTISSURFACE
      obs_value  = compute_altimeter(zob, vloc)  !  altimeter is hPa
   elseif(obs_kind_gen == KIND_SURFACE_PRESSURE) then
      obs_value = zob * 100.0_r8  !  for Ps variable only in Pascal
      vloc = lev                  ! station height, not used now for Ps obs
      which_vert = VERTISSURFACE
   else if(obs_kind_gen == KIND_SPECIFIC_HUMIDITY) then
      obs_err = obs_err*1.0e-3_r8
      obs_value = zob*1.0e-3_r8     !  for Q variable to kg/kg
   else
      obs_value = zob               !  for T, U, V, RH, Tdew
   endif

   aqc = iqc
   days = day0

   if ( obs_time ) then

      seconds = time * 3600

   else

      if ( hourt == '' ) then
         call get_time(increment_time(set_date(year,month,day,0,0,0), &
                        nint(time / 6.0_r8) * 21600, 0), seconds, days)
      else
         read(hourt, fmt='(i2)') hour
         call get_time(increment_time(set_date(year,month,day,0,0,0), & 
                         hour*3600, 0), seconds, days)
      endif

   end if

   ! fill the time_obs variable to pass into add routine
   time_obs = set_time(seconds, days)

   !   create the obs_def for this observation, add to sequence
   !------------------------------------------------------------------------------

   call create_3d_obs(lat, lon, vloc, which_vert, obs_value, &
                      obs_kind, obs_err, days, seconds, aqc, obs)
   call add_obs_to_seq(real_obs_sequence, obs, time_obs, prev_obs, prev_time, first_obs)

end do obsloop

200 continue

close(obs_unit)

print*, 'date ', obsdate
print*, 'num obs used = ', obs_num, ' total obs skipped = ', sum(iskip)
do i=1, num_fail_kinds
  if (iskip(i) >  0) print *, iskip(i), 'skipped because ', skip_reasons(i)
enddo

end function real_obs_sequence


subroutine initialize_module
!-------------------------------------------------
call register_module(source, revision, revdate)
module_initialized = .true.
end subroutine initialize_module


end module real_obs_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
