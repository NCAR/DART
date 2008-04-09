! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module real_obs_mod

! <next five lines automatically updated by CVS, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use types_mod,        only : r8, rad2deg, PI
use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use time_manager_mod, only : time_type, operator(>), operator(<), operator(>=), &
                             operator(/=), set_date, set_calendar_type, get_time, &
                             get_date, set_time
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             E_ERR, E_MSG, timestamp
use     location_mod, only : location_type, set_location, VERTISPRESSURE, VERTISSURFACE
use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, &
                             copy_obs, set_copy_meta_data, set_qc_meta_data, set_obs_def

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, &
                             KIND_V_WIND_COMPONENT, KIND_SURFACE_PRESSURE, &
                             KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, KIND_PRESSURE, &
                             KIND_VERTICAL_VELOCITY, KIND_RAINWATER_MIXING_RATIO, &
                             KIND_DEW_POINT_TEMPERATURE, KIND_DENSITY, KIND_VELOCITY, &
                             KIND_1D_INTEGRAL, KIND_RADAR_REFLECTIVITY 

use  obs_def_altimeter_mod, only: compute_altimeter

use obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT
use obs_kind_mod, only : RADIOSONDE_V_WIND_COMPONENT
use obs_kind_mod, only : RADIOSONDE_SURFACE_PRESSURE
use obs_kind_mod, only : RADIOSONDE_TEMPERATURE
use obs_kind_mod, only : RADIOSONDE_SPECIFIC_HUMIDITY
use obs_kind_mod, only : AIRCRAFT_U_WIND_COMPONENT
use obs_kind_mod, only : AIRCRAFT_V_WIND_COMPONENT
use obs_kind_mod, only : AIRCRAFT_TEMPERATURE
use obs_kind_mod, only : ACARS_U_WIND_COMPONENT
use obs_kind_mod, only : ACARS_V_WIND_COMPONENT
use obs_kind_mod, only : ACARS_TEMPERATURE
use obs_kind_mod, only : MARINE_SFC_U_WIND_COMPONENT
use obs_kind_mod, only : MARINE_SFC_V_WIND_COMPONENT
use obs_kind_mod, only : MARINE_SFC_TEMPERATURE
use obs_kind_mod, only : MARINE_SFC_SPECIFIC_HUMIDITY
use obs_kind_mod, only : MARINE_SFC_ALTIMETER
use obs_kind_mod, only : LAND_SFC_U_WIND_COMPONENT
use obs_kind_mod, only : LAND_SFC_V_WIND_COMPONENT
use obs_kind_mod, only : LAND_SFC_TEMPERATURE
use obs_kind_mod, only : LAND_SFC_SPECIFIC_HUMIDITY
use obs_kind_mod, only : LAND_SFC_ALTIMETER
use obs_kind_mod, only : RADIOSONDE_SURFACE_ALTIMETER
use obs_kind_mod, only : SAT_U_WIND_COMPONENT
use obs_kind_mod, only : SAT_V_WIND_COMPONENT
use obs_kind_mod, only : ATOV_TEMPERATURE


implicit none
private

public :: real_obs_sequence

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, save :: module_initialized = .false.
! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000
!-------------------------------------------------

contains


function real_obs_sequence (year, month, day, hourt, max_num, select_obs, &
          ObsBase, ADDUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND, &
          obs_U, obs_V, obs_T, obs_PS, obs_QV, bin_beg, bin_end,           &
          lon1, lon2, lat1, lat2)
!------------------------------------------------------------------------------
!  this function is to prepare NCEP decoded BUFR data to DART sequence format
!
integer,            intent(in) :: year, month, day, max_num, select_obs
character(len = *), intent(in) :: ObsBase, hourt
logical,            intent(in) :: ADDUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND
logical,            intent(in) :: obs_U, obs_V, obs_T, obs_PS, obs_QV
real(r8),           intent(in) :: lon1, lon2, lat1, lat2

type(obs_sequence_type) :: real_obs_sequence


type(obs_type) :: obs, prev_obs
integer :: i, num_copies, num_qc
integer :: days, seconds
integer :: day0, sec0
integer :: hour, imin, sec
integer :: obs_num, calender_type, iskip
type(time_type) :: current_day, next_day

integer :: obs_unit
integer :: obs_prof, obs_kind, obs_kind_gen, which_vert, iqc, obstype, pc
real (r8) :: obs_err, lon, lat, lev, zob, time, pre_time, rcount, zob2
real (r8) :: vloc, obs_value, aqc, var2, lon2c, lonc

real (r8) :: bin_beg, bin_end

character(len = 8 ) :: obsdate
character(len = 80) :: obsfile, label
character(len = 129) :: copy_meta_data, qc_meta_data
character(len = 6 ) :: subset
logical :: pass

if ( .not. module_initialized ) call initialize_module

num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence 

call init_obs_sequence(real_obs_sequence, num_copies, num_qc, max_num)

! set meta data of obs_seq

do i = 1, num_copies
   copy_meta_data = 'NCEP BUFR observation'  
   call set_copy_meta_data(real_obs_sequence, i, copy_meta_data)
end do

do i = 1, num_qc
   qc_meta_data = 'NCEP QC index'
   call set_qc_meta_data(real_obs_sequence, i, qc_meta_data)
end do

! Initialize the obs variable

call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! set observation time type
calender_type = 3
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
iskip   = 0

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
      iskip = iskip + 1
      cycle obsloop 
   endif 

   !  select the obs for the time window
   if(time < bin_beg .or. time > bin_end) then
      iskip = iskip + 1
      cycle obsloop
   endif

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

   obs_prof = rcount/1000000

!   assign each observation the correct observation type
!------------------------------------------------------------------------------

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
    obs_kind_gen = KIND_SPECIFIC_HUMIDITY
    if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_SPECIFIC_HUMIDITY
    if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_SPECIFIC_HUMIDITY
    if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_SPECIFIC_HUMIDITY
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
         if( (obs_T  .and. (obs_kind_gen == KIND_TEMPERATURE ))      .or. &
             (obs_U  .and. (obs_kind_gen == KIND_U_WIND_COMPONENT )) .or. &
             (obs_V  .and. (obs_kind_gen == KIND_V_WIND_COMPONENT )) .or. &
             (obs_PS .and. (obs_kind_gen == KIND_SURFACE_PRESSURE))  .or. &
             (obs_QV .and. (obs_kind_gen == KIND_SPECIFIC_HUMIDITY))      ) then
             pass = .false.
         endif

      endif
     
      ! if pass is still true, we want to ignore this obs.
      if(pass) then
         iskip = iskip + 1
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
       subset == 'SATWND' .or. subset == 'ADPUPA' ) then
       vloc = lev*100.0_r8          ! convert from mb to Pascal
       which_vert = VERTISPRESSURE
   endif

   ! set vertical coordinate for surface observations
   if (subset == 'ADPSFC' .or. subset == 'SFCSHP') then
     vloc = lev
     if ( subset == 'SFCSHP' ) vloc = 0.0_r8
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
      obs_value = zob               !  for T, U, V
   endif
   var2 = obs_err**2            ! error variance

   aqc = iqc
   seconds = time * 3600
   days = day0

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


end module real_obs_mod
