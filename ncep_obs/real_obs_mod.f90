! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module real_obs_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use types_mod,   only : r8, rad2deg
use obs_def_mod, only : obs_def_type, get_obs_def_time, read_obs_def, &
         write_obs_def, destroy_obs_def, interactive_obs_def, copy_obs_def, &
         set_obs_def_time, set_obs_def_error_variance, &
         set_obs_def_kind, set_obs_def_location
use time_manager_mod, only : time_type, operator(>), operator(<), operator(>=), &
                             operator(/=), set_date, set_calendar_type, get_time
use utilities_mod, only : get_unit, open_file, close_file, file_exist, check_nml_error, &
                          register_module, error_handler, E_ERR, E_MSG
use obs_kind_mod, only : obs_kind_type, set_obs_kind 
use location_mod, only : location_type, set_location

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, copy_obs, &
                             set_copy_meta_data, set_qc_meta_data, copy_obs, set_obs_def
use time_manager_mod, only : time_type, read_time, set_time

implicit none

private

! Public interfaces for obs sequences
public  real_obs_sequence, obs_model_type

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.
!-------------------------------------------------------------      

contains

!----------------------------------------------------------------------
  subroutine initialize_module
!----------------------------------------------------------------------
! subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module


!=================================================
  function real_obs_sequence()
!=================================================
!  this function is to prepare NCEP BUFR data to DART sequence format
!
  type(obs_sequence_type) :: real_obs_sequence
  type(obs_type) :: obs, prev_obs
  integer :: max_num_obs, num_copies, num_qc, i
  integer :: days, seconds, obs_num,  calender_type, obs_day01
  integer :: obs_year, obs_month, obs_day, obs_hour, obs_min, obs_sec
  type(time_type) :: obs_time

  integer :: obs_unit, obs_prof, obs_kind, model_type, which_vert, iqc
  real (r8) :: obs_err, lon, lat, lev, zob, time, type, count, zob2
  real (r8) :: vloc, obs_value, lon01, lat01, aqc, var2

  character(len = 8 ) :: obsdate
  character(len = 80) :: obsfile
  character(len = 129) :: copy_meta_data, qc_meta_data

!-------------------------------------------------------------

!*****************************************************************************
   max_num_obs = 800000    !! for current NCEP daily RA+ACARS+SATWND data only
!*****************************************************************************

   num_copies  = 1
   num_qc      = 1         !! the NCEP data passed NCEP QC procedure by readpb.f
!
! Initialize an obs_sequence structure
!
  call init_obs_sequence(real_obs_sequence, num_copies, num_qc, max_num_obs)

  do i = 1, num_copies
   copy_meta_data = 'NCEP BUFR observation'  
   call set_copy_meta_data(real_obs_sequence, i, copy_meta_data)
  end do

  do i = 1, num_qc
   qc_meta_data = 'NCEP QC index'
   call set_qc_meta_data(real_obs_sequence, i, qc_meta_data)
  end do

!  Initialize the obs variable
!
    call init_obs(obs, num_copies, num_qc)
    call init_obs(prev_obs, num_copies, num_qc)

!-------------------------------------
!   for NCEP real data preparation
!-------------------------------------

!   set observation time type
    calender_type = 3
    call set_calendar_type(calender_type)

!   read in namelist
    obs_unit = get_unit()
    open(obs_unit, file='/home/hliu/newDART/ncep_obs/ncepobs.input', form='formatted')
    read(obs_unit, *) obs_year, obs_month, obs_day
    close(obs_unit)

     if(obs_month .lt. 10) then
       if(obs_day .lt. 10) then
        write(obsdate, '(i4,A1,i1,A1,i1)') obs_year, '0',obs_month, '0',obs_day
       else
        write(obsdate, '(i4,A1,i1,i2)') obs_year, '0',obs_month, obs_day
       endif
     else
       if(obs_day .lt. 10) then
       write(obsdate, '(i4,i2,A1,i1)') obs_year, obs_month, '0',obs_day
       else
       write(obsdate, '(i4,i2,i2)') obs_year, obs_month, obs_day
       endif
     endif
     print*, 'ncep obsdate = ', obsdate

!  open unit for NCEP observation file
    obs_unit = get_unit()
    obsfile  = '/home/hliu/ncepobs/test/temp_obs.'//obsdate
    open(unit = obs_unit, file = obsfile, form='formatted', status='old')

     print*, 'file opened= ', obsfile
    rewind (obs_unit)

    obs_num = 0

! Loop to set up each observation
!
obsloop:  do i = 1, max_num_obs
   
     if(mod(i, 1000) ==0) print*, 'doing obs = ', i
! read in each obs from the daily observation file (include 06, 12, 18, 24Z)

     read(obs_unit, 880, end=200) obs_err, lon, lat, lev, zob, zob2, count,time,type, iqc
     obs_num = obs_num + 1

 880 format(f4.2, 2f7.3, e12.5, f7.2, f7.2, f9.0, f7.3, f5.0, i3)

!  set up observation location

    lon01 = lon*rad2deg            ! in degree
    lat01 = lat*rad2deg            ! in degree
    if(lon01 >= 360.0_r8) lon01 = 360.0_r8
    if(lon01 <=   0.0_r8) lon01 =   0.0_r8

    if(lat01 >=  90.0_r8) lat01 =  90.0_r8
    if(lat01 <= -90.0_r8) lat01 = -90.0_r8

    obs_prof = count/1000000
!   obs_kind = obs_prof*10000 + type 
!   call obs_model_type(obs_kind, model_type)
                                                                                       
    if(obs_prof == 2) obs_kind = 1           ! U
    if(obs_prof == 9) obs_kind = 2           ! V
    if(obs_prof == 3) obs_kind = 3           ! Ps
    if(obs_prof == 1) obs_kind = 4           ! T
    if(obs_prof == 5) obs_kind = 5           ! q

    model_type = obs_kind

    if (model_type == 1 .or. model_type ==2 .or. model_type ==4 .or. model_type==5) then
     vloc = lev*100.0            ! (transfer Pressure coordinate from mb to Pascal) 
     which_vert = 2
    endif

    if (model_type == 3) then    ! for Ps
     vloc = lev            ! station height, not used now for Ps obs
     which_vert = -1
     obs_err = obs_err*100.0     ! convert obs_err to Pa
    endif

    if (model_type == 5) then    ! for Q
     obs_err = obs_err*1.0e-3     ! convert obs_err to kg/kg
    endif

!   set obs value and error covariance

     if(model_type == 3) then
      obs_value = zob*100.0      !  for Ps variable only in Pascal
     else if(model_type == 5) then
      obs_value = zob*1.0e-3     !  for Q variable to kg/kg
     else
      obs_value = zob            !  for T, U, V
     endif

    aqc = iqc
    var2 = obs_err**2                ! error_covariance


!   set obs time

      obs_min  = 60*( time - ifix(time) )
      if(time >= 24) then
       obs_hour   = time -24
       obs_day01  = obs_day + 1
      else
       obs_hour   = time
       obs_day01  = obs_day 
      endif

      obs_sec  = 0

      obs_time = set_date(obs_year, obs_month, obs_day01, obs_hour, obs_min, obs_sec)
   
      call get_time(obs_time, seconds, days)
!     write(*,*)  'obs time', seconds, days

!-----------------------------------------------------------------------------------
   call real_obs(num_copies, num_qc, obs, &
       lon01, lat01, vloc, obs_value, var2, aqc, obs_kind, which_vert, seconds, days)
!-----------------------------------------------------------------------------------

   if(i == 1) then
     call insert_obs_in_seq(real_obs_sequence, obs)
     call copy_obs(prev_obs, obs)
   else
!    for time ordered observation only.
!    call insert_obs_in_seq(real_obs_sequence, obs, prev_obs)
!    call copy_obs(prev_obs, obs)
!
!    for random time order observation, have to do time sort, CPU costly.
     call insert_obs_in_seq(real_obs_sequence, obs)
   endif


end do obsloop

 200 continue
   close(obs_unit)

   print*, 'total obs num= ', obs_num, obsdate

end function real_obs_sequence


!-----------------------------------------------------------------------------------
  subroutine real_obs(num_copies, num_qc, obs, &
       lon01, lat01, vloc, obs_value, var2, aqc, obs_kind, which_vert, seconds, days)
!-----------------------------------------------------------------------------------

  integer, intent(in) :: num_copies, num_qc
  type(obs_type), intent(inout) :: obs
  type(obs_def_type) :: obsdef0

  integer :: i

  integer :: obs_kind, which_vert, days, seconds
  real (r8) :: vloc, obs_value, lon01, lat01, aqc, var2
  real (r8) :: aqc01(1), obs_value01(1)

! Does real initialization of an observation type

  call real_obs_def(obsdef0, &
       lon01, lat01, vloc, obs_value, var2, aqc, obs_kind, which_vert, seconds, days)
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

!------------------------------------------------------------------------------------
subroutine real_obs_def(obs_def, &
       lon01, lat01, vloc, obs_value, var2, aqc, obs_kind, which_vert, seconds, days)
!------------------------------------------------------------------------------------
  type(obs_def_type), intent(inout) :: obs_def

  integer :: ind
  type(obs_kind_type) :: kind0
  type(location_type) :: loc0

  integer :: obs_kind, which_vert, days, seconds
  real (r8) :: vloc, obs_value, lon01, lat01, aqc, var2

  if ( .not. module_initialized ) call initialize_module

! set the location
    loc0 = set_location(lon01, lat01, vloc, which_vert )
  call set_obs_def_location(obs_def, loc0)

! set the obsrvation kind
  kind0 = set_obs_kind(obs_kind)
  call set_obs_def_kind(obs_def, kind0)

  call set_obs_def_time(obs_def, set_time(seconds, days) )
  call set_obs_def_error_variance(obs_def, var2)

end subroutine real_obs_def

!-----------------------------------------------------
subroutine obs_model_type(kind, model_type, ncep_type)
!-----------------------------------------------------
  implicit none

  integer, intent(in)  :: kind
  integer, intent(out) :: model_type
  integer, intent(out), optional :: ncep_type
  integer :: obs_prof

  if ( .not. module_initialized ) call initialize_module

    obs_prof = kind/10000
    if(obs_prof == 2) model_type = 1         !! the model obs types
    if(obs_prof == 9) model_type = 2
    if(obs_prof == 3) model_type = 3
    if(obs_prof == 1) model_type = 4
    if(obs_prof == 5) model_type = 5

    ncep_type = kind - obs_prof * 10000

end subroutine obs_model_type


end module real_obs_mod
