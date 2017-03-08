! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
PROGRAM rad_3dvar_to_dart

use         types_mod, only : r8, missing_r8, missing_data, DEG2RAD, earth_radius
use     utilities_mod, only : open_file, close_file, initialize_utilities, &
                              register_module, logfileunit, nmlfileunit, E_MSG, timestamp, &
                              error_handler, find_namelist_in_file, check_namelist_read, &
                              finalize_utilities
use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, insert_obs_in_seq, &
                              set_copy_meta_data, set_qc_meta_data, write_obs_seq, assignment(=), &
                              init_obs, static_init_obs_sequence, set_obs_def, set_obs_values, set_qc
use       obs_def_mod, only : set_obs_def_location, set_obs_def_error_variance, &
                              set_obs_def_type_of_obs, set_obs_def_time, set_obs_def_key, &
                              obs_def_type
use      obs_kind_mod, only : DOPPLER_RADIAL_VELOCITY, RADAR_REFLECTIVITY
use      location_mod, only : location_type, set_location
use  time_manager_mod, only : time_type, set_date, set_calendar_type, GREGORIAN
use obs_def_radar_mod, only : set_radial_vel

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time

INTEGER           :: iunit, iost, io

character(len=19) :: radar_name
real(r8)          :: radar_long, radar_lat, radar_elev
real(r8)          :: orientation(3), nyquist_velocity
character(len=19) :: start_scan_date
integer           :: num_prof, max_levels

character(len=80) :: dummy

character(len=12) :: platform_name
integer           :: year, month, day, hours, minutes, seconds
real(r8)          :: lat,lon,elv
integer           :: levels

integer           :: ii, num_Radar, rv_qc, rf_qc, key, it
integer           :: num_obs, num_copies, num_qc, max_num_obs

real(r8)          :: height, rv_inv, rf_inv, rv_error, rf_error
real(r8)          :: obs_value(1), rstatus(1,1)
real(r8)          :: h, spath, x, y, rad_lon, rad_lat, obs_lon, obs_lat
real(r8)          :: rgate, raz, elev_rad, elev_obs, ae

!-----------------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: var_file              = 'qc_radr_3dvar_2002083100.dat', &
                        obs_seq_out_file_name = 'obs_seq.out'
integer              :: calendar_type         = GREGORIAN

namelist /rad_3dvar_to_dart_nml/ var_file, obs_seq_out_file_name, calendar_type

!------------------------------------------------------------------------------

call initialize_utilities('rad_3dvar_to_dart')
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "rad_3dvar_to_dart_nml", iunit)
read(iunit, nml = rad_3dvar_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "rad_3dvar_to_dart_nml")

! Record the namelist values used for the run ...
write(nmlfileunit, nml=rad_3dvar_to_dart_nml)
write(     *     , nml=rad_3dvar_to_dart_nml)

call set_calendar_type(calendar_type)

iunit = open_file(var_file, action = 'read')

! -------------------------------------------------------------------
! Initialize the counters:

num_Radar = 0
num_obs = 0
key = 0

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

call set_copy_meta_data(seq, 1, 'MM5 3D-VAR 2.0 Radar observation')
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
   WRITE (*,'(/,A,I3,/)') ' END OF UNIT: ',iunit
   WRITE (*,'(A,I3)')     ' IOSTAT == ',iost
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
      key = key + 1

      location = set_location(lon, lat, height, 3)
      call set_obs_def_location(obs_def, location)

      location = set_location(radar_long, radar_lat, radar_elev, 3)
      call set_obs_def_key(obs_def, key)

      obs_lat = lat*DEG2RAD
      rad_lat = radar_lat*DEG2RAD
      obs_lon = lon*DEG2RAD
      rad_lon = radar_long*DEG2RAD

      ae = 1000.0_r8 * earth_radius
      x = ae * cos((obs_lat + rad_lat)/2.0_r8) * (obs_lon - rad_lon)
      y = ae * (obs_lat - rad_lat)
      raz = atan(x/y)
      spath = sqrt(x*x + y*y)
      h = height - radar_elev

      ae = 4000.0_r8 * earth_radius / 3.0_r8
      rgate = spath

      do it=1,10
         elev_rad = asin((h*h + 2.0_r8*ae*h - rgate*rgate)/(2.0_r8*ae*rgate))
         rgate = (ae+h)*sin(spath/ae)/cos(elev_rad)
      enddo

      elev_obs = sqrt(1.5_r8 * h / (1000.0_r8 * earth_radius) + elev_rad*elev_rad)

      orientation(1) = sin(raz)*cos(elev_obs)
      orientation(2) = cos(raz)*cos(elev_obs)
      orientation(3) = sin(elev_obs)

      nyquist_velocity = missing_r8

      call set_radial_vel(key, location, orientation, nyquist_velocity)

      call set_obs_def_time(obs_def, time)

      call set_obs_def_type_of_obs(obs_def, DOPPLER_RADIAL_VELOCITY)

      if (rv_error /= missing_r8) then
         call set_obs_def_error_variance(obs_def, rv_error*rv_error)
      else
         call set_obs_def_error_variance(obs_def, 2.0_r8)
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

      call set_obs_def_type_of_obs(obs_def, RADAR_REFLECTIVITY)

      if (rv_error /= missing_r8) then
         call set_obs_def_error_variance(obs_def, rf_error*rf_error)
      else
         call set_obs_def_error_variance(obs_def, 2.0_r8)
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

write(unit=*, fmt='(5x,a,i6,a)') 'Read:  ', num_Radar, ' Radar reports,'

call write_obs_seq(seq, obs_seq_out_file_name)

call error_handler(E_MSG, 'rad_3dvar_to_dart', 'FINISHED rad_3dvar_to_dart.')
call error_handler(E_MSG, 'rad_3dvar_to_dart', 'Finished successfully.', &
                   source,revision,revdate)
call finalize_utilities()

 
END PROGRAM rad_3dvar_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
