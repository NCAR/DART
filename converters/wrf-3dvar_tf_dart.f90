! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
PROGRAM wrf_3dvar_tf_dart

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use        types_mod, only : r8, missing_r8, missing_data
use    utilities_mod, only : open_file, check_nml_error, close_file, file_exist, initialize_utilities, &
                             finalize_utilities, register_module, logfileunit
use obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, insert_obs_in_seq, &
                             set_copy_meta_data, set_qc_meta_data, write_obs_seq, assignment(=), &
                             init_obs, static_init_obs_sequence, set_obs_def, set_obs_values, set_qc
use      obs_def_mod, only : set_obs_def_location, set_obs_def_error_variance, set_obs_def_kind, &
                             set_obs_def_time, obs_def_type, set_obs_def_platform
use     obs_kind_mod, only : obs_kind_type, set_obs_kind
use     location_mod, only : location_type, set_location
use time_manager_mod, only : time_type, set_date, set_calendar_type, GREGORIAN
use     platform_mod, only : set_platform_location, platform_type

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(obs_kind_type)     :: kind
type(time_type)         :: time
type(platform_type)     :: platform

INTEGER           :: iunit, ierr, iost, io

character(len=19) :: radar_name
real(r8)          :: radar_long, radar_lat, radar_elev
character(len=19) :: start_scan_date
integer           :: num_prof, max_levels

character(len=80) :: dummy

character(len=12) :: platform_name
integer           :: year, month, day, hours, minutes, seconds
real(r8)          :: lat,lon,elv
integer           :: levels

integer           :: ii, num_Radar, rv_qc, rf_qc
integer           :: num_obs, num_copies, num_qc, max_num_obs

real(r8)          :: height, rv_inv, rf_inv, rv_error, rf_error
real(r8)          :: obs_value(1), rstatus(1,1)

!-----------------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: wrf_3dvar_file        = 'qc_radr_3dvar_2002083100.dat', &
                        obs_seq_out_file_name = 'obs_seq.out'
integer              :: calendar_type         = GREGORIAN

namelist /wrf_3dvar_tf_dart_nml/ wrf_3dvar_file, obs_seq_out_file_name, calendar_type

!------------------------------------------------------------------------------

call initialize_utilities
call register_module(source, revision, revdate)
write(logfileunit,*)'STARTING wrf_3dvar_tf_dart ...'
call error_handler(E_MSG,'wrf_3dvar_tf_dart','STARTING ...',source,revision,revdate)

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = wrf_3dvar_tf_dart_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'wrf_3dvar_tf_dart_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! Record the namelist values used for the run ...
call error_handler(E_MSG,'wrf_3dvar_tf_dart_nml','wrf_3dvar_tf_dart_nml values are',' ',' ',' ')
write(logfileunit, nml=wrf_3dvar_tf_dart_nml)
write(     *     , nml=wrf_3dvar_tf_dart_nml)

call set_calendar_type(calendar_type)

iunit = open_file(wrf_3dvar_file, action = 'read')

! -------------------------------------------------------------------
! Initialize the counters:

num_Radar = 0
num_obs = 0

!-----------------------------------------------------------------------------!
! Read the header of a MM5 3D-VAR 2.0 Radar observation file
!-----------------------------------------------------------------------------!

READ (UNIT = iunit, IOSTAT = iost, &
     FMT = '(A19,F8.3,2X,F8.3,F10.1,2X,A19,2I6)' ) &
     radar_name, radar_long, radar_lat, radar_elev, &
     start_scan_date, &
     num_prof, max_levels

READ (UNIT = iunit, IOSTAT = iost, FMT = '(A80)' ) dummy
READ (UNIT = iunit, IOSTAT = iost, FMT = '(A80)' ) dummy

call static_init_obs_sequence()

max_num_obs = num_prof*max_levels*2.0_r8
num_copies = 1
num_qc = 1

! Initialize an obs_sequence structure
call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)

call set_copy_meta_data(seq, 1, 'observations')
call set_qc_meta_data(seq, 1, 'missing')

call init_obs(obs, num_copies, num_qc)

!  READ FORMATS
!  ------------

!  LOOP OVER RECORDS
!  -----------------

reports: &
     DO
!     READ STATION GENERAL INFO
!     =============================

READ (UNIT = iunit, IOSTAT = iost, &
     FMT = '(A12,3X,I4,5(A1,I2),2X,2(F12.3,2X),F8.1,2X,I6)' ) &
     platform_name,  &
     year, dummy, month, dummy, day, dummy, hours, dummy, minutes, dummy, seconds, &
     lat,       &
     lon,       &
     elv,       &
     levels

time = set_date(year, month, day, hours, minutes, seconds)

IF (iost /= 0) THEN
   WRITE (0,'(/,A,I3,/)') ' END OF UNIT: ',iunit
   WRITE (0,'(A,I3)')     ' IOSTAT == ',iost
   EXIT reports
ENDIF

!     READ EACH LEVELS
!     ----------------

loop_level: DO ii = 1, levels

   READ (UNIT = iunit, FMT = '( 3X, F12.1, 2(F12.3,I4,F12.3,2X) )' ) &
        height,         &
        rv_inv,         &
        rv_qc,          &
        rv_error,       &
        rf_inv,         &
        rf_qc,          &
        rf_error

   if (rv_inv /= missing_r8) then

      num_obs = num_obs + 1

      location = set_location(lon, lat, height, 3)
      call set_obs_def_location(obs_def, location)

      location = set_location(radar_long, radar_lat, radar_elev, 3)
      call set_platform_location(platform, location)
      call set_obs_def_platform(obs_def, platform)

      call set_obs_def_time(obs_def, time)

      kind = set_obs_kind(100)
      call set_obs_def_kind(obs_def, kind)

      if (rv_error /= missing_r8) then
         call set_obs_def_error_variance(obs_def, rv_error*rv_error)
      else
         call set_obs_def_error_variance(obs_def, missing_r8)
      endif

      call set_obs_def(obs, obs_def)

      obs_value(1) = rv_inv
      call set_obs_values(obs, obs_value, 1)

      if (rv_inv == missing_r8 .or. &
           rv_error == missing_r8 ) then

         rv_qc = missing_data

      end if

      rstatus(1,1) = rv_qc
      call set_qc(obs, rstatus(1,:), 1)

      if(num_obs == 1) then
         call insert_obs_in_seq(seq, obs)
      else
         call insert_obs_in_seq(seq, obs, prev_obs)
      endif

      prev_obs = obs

   endif

   if (rf_inv /= missing_r8) then

      num_obs = num_obs + 1

      kind = set_obs_kind(101)
      call set_obs_def_kind(obs_def, kind)

      if (rv_error /= missing_r8) then
         call set_obs_def_error_variance(obs_def, rf_error*rf_error)
      else
         call set_obs_def_error_variance(obs_def, missing_r8)
      endif

      call set_obs_def(obs, obs_def)

      obs_value(1) = rf_inv
      call set_obs_values(obs, obs_value, 1)

      if (rf_inv == missing_r8 .or. &
           rf_error == missing_r8 ) then

         rf_qc = missing_data

      end if
 
      rstatus(1,1) = rf_qc
      call set_qc(obs, rstatus(1,:), 1)

      if(num_obs == 1) then
         call insert_obs_in_seq(seq, obs)
      else
         call insert_obs_in_seq(seq, obs, prev_obs)
      endif

      prev_obs = obs

   endif

ENDDO loop_level

num_Radar = num_Radar + 1
   
ENDDO reports

call close_file(iunit)                                                        

!  PRINT OUT
!  =============
 
write(unit=*, fmt='(5x,a,i6,a)') &
     'Read:  ', num_Radar, ' Radar reports,'

! Write out the sequence
call write_obs_seq(seq, obs_seq_out_file_name)

write(logfileunit,*)'FINISHED wrf_3dvar_tf_dart.'
write(logfileunit,*)

call finalize_utilities ! closes the log file.
 
END PROGRAM wrf_3dvar_tf_dart
