! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module gts_dart_mod

use         types_mod, only : r8, missing_r8, missing_data, DEG2RAD
use     utilities_mod, only : open_file, close_file, initialize_utilities, &
                              register_module, logfileunit, E_MSG, &
                              error_handler, find_namelist_in_file, check_namelist_read
use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, insert_obs_in_seq, &
                              write_obs_seq, init_obs, assignment(=), &
                              static_init_obs_sequence, set_obs_def, set_obs_values, set_qc
use       obs_def_mod, only : set_obs_def_location, set_obs_def_error_variance, &
                              set_obs_def_type_of_obs, set_obs_def_time, set_obs_def_key, &
                              obs_def_type
use   obs_def_gps_mod, only : set_gpsro_ref
!use      obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT, &
!                              RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, &
!                              RADIOSONDE_TEMPERATURE, RADIOSONDE_SPECIFIC_HUMIDITY, &
!                              DEW_POINT_TEMPERATURE, &
!                              METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, METAR_TEMPERATURE_2_METER, &
!                              METAR_SPECIFIC_HUMIDITY_2_METER, METAR_SURFACE_PRESSURE, METAR_POT_TEMP_2_METER
use      location_mod, only : location_type, set_location, get_location, &
                              VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use  time_manager_mod, only : time_type, set_date, set_calendar_type, GREGORIAN

use DA_Constants
use DA_Define_Structures
use module_obs

Implicit none

Interface insert_gts_obs

  module procedure insert_gts_obs_single, insert_gts_obs_multi, insert_gts_obs_gpsref

End interface

contains

!-------------------------------------------------------------------------------
SUBROUTINE insert_gts_obs_single(gtsobs, qc_threshold, &
                          pressure_type, height_type, u_wind_type,v_wind_type, &
                          temperature_type, dew_point_type, rh_type, &
                          which_vert, num_obs, obs, seq)
!-------------------------------------------------------------------------------

   implicit none

   type(single_level_type), intent(in)    :: gtsobs
   integer, optional, intent(in)          :: pressure_type, height_type, &
                                             u_wind_type, v_wind_type, &
                                             temperature_type, dew_point_type, &
                                             rh_type, qc_threshold
   integer, intent(in)                    :: which_vert
   integer, intent(inout)                 :: num_obs
   type(obs_type), intent(inout)          :: obs
   type(obs_sequence_type), intent(inout) :: seq

   type(obs_def_type)      :: obs_def
   type(location_type)     :: location
   type(time_type)         :: time

   character(len=80)       :: dummy
   integer                 :: year, month, day, hours, minutes, seconds
   real(r8)                :: lat,lon,pressure, height

   integer                 :: i, levels
   real(r8)                :: windspd, winddir, uu, vv, tk, td, rh, qv
   real(r8)                :: obs_value(1), obs_err, obs_qc(1)

   logical                 :: valid


   ! time info
   read(gtsobs%info%date_char, '(I4,5(A1,I2))')  &
        year, dummy, month, dummy, day, dummy, hours, dummy, minutes, dummy, seconds
   time = set_date(year, month, day, hours, minutes, seconds)
   call set_obs_def_time(obs_def, time)

   ! number of vertical levels
   levels = gtsobs%info%levels

   ! location info
   lat = gtsobs%info%lat
   lon = gtsobs%info%lon
   ! Dart longitude from 0 to 360
   if(lon < 0.0_r8) lon = lon + 360.0_r8
   height = gtsobs%info%elv
   pressure = gtsobs%each%pressure%data
   if (which_vert.eq.VERTISSURFACE .or. which_vert.eq.VERTISHEIGHT) then
      location = set_location(lon, lat, height, which_vert)
   else if (which_vert.eq.VERTISPRESSURE) then
      location = set_location(lon, lat, pressure, which_vert)
   end if
   call set_obs_def_location(obs_def, location)

   if (present(u_wind_type) .or. present(v_wind_type)) then
      valid = .true. 
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0) &
            valid =       (gtsobs%each%speed%qc .ge. 0) &
                    .and. (gtsobs%each%speed%qc .le. qc_threshold) &
                    .and. (gtsobs%each%direction%qc .ge.0) &
                    .and. (gtsobs%each%direction%qc .le. qc_threshold)
      endif
      !convert wind speed and direction to zonal and meridional components
      windspd = gtsobs%each%speed%data
      winddir = gtsobs%each%direction%data
      call spddir2uv(windspd,winddir,uu,vv,valid)
      if (valid) then
         ! Insert obs into DART obs sequence
         obs_value(1) = uu
         obs_err = gtsobs%each%speed%error
         obs_qc(1)  = real(max(gtsobs%each%speed%qc, gtsobs%each%direction%qc),kind=r8)
         call set_obs_def_type_of_obs(obs_def, u_wind_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
         obs_value(1) = vv
         num_obs = num_obs + 1
         call set_obs_def_type_of_obs(obs_def, v_wind_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(temperature_type)) then
      !temperature in K
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0) &
             valid =      (gtsobs%each%temperature%qc .ge. 0) &
                    .and. (gtsobs%each%temperature%qc .le. qc_threshold)
      endif
      tk = gtsobs%each%temperature%data
      if ( valid .and. tk > 0.0_r8 .and. tk < 1.0e3_r8 ) then
         obs_value(1) = tk
         obs_err = gtsobs%each%temperature%error
         obs_qc(1)  = real(gtsobs%each%temperature%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, temperature_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(dew_point_type)) then
      !dew point in K
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each%dew_point%qc .ge. 0) &
                   .and. (gtsobs%each%dew_point%qc .le. qc_threshold)
      endif
      td = gtsobs%each%dew_point%data
      if ( valid .and. td > 0.0_r8 .and. td < 1.0e3_r8 ) then
         obs_value(1) = td
         obs_err = gtsobs%each%dew_point%error
         obs_qc(1) = real(gtsobs%each%dew_point%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, dew_point_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(rh_type)) then
      !relative humidity rh
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each%rh%qc .ge. 0) &
                   .and. (gtsobs%each%rh%qc .le. qc_threshold)
      endif
      rh = gtsobs%each%rh%data
      if ( valid .and. rh >= 0.0_r8 .and. rh <= 1.2e2_r8 ) then
         obs_value(1) = rh
         obs_err = gtsobs%each%rh%error
         obs_qc(1)  = real(gtsobs%each%rh%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, rh_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(pressure_type)) then
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each%pressure%qc .ge. 0) &
                   .and. (gtsobs%each%pressure%qc .le. qc_threshold)
      endif
      pressure = gtsobs%each%pressure%data
      if ( valid .and. pressure >= -100.0e3_r8 .and. pressure < 100.0e3_r8 ) then
         obs_value(1) = pressure
         obs_err = gtsobs%each%pressure%error
         obs_qc(1)  = real(gtsobs%each%pressure%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, pressure_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(height_type)) then
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each%height%qc .ge. 0) &
                   .and. (gtsobs%each%height%qc .le. qc_threshold)
      endif
      height = gtsobs%each%height%data
      if ( valid .and. height >= -100.0e3_r8 .and. height < 100.0e3_r8 ) then
         obs_value(1) = height
         obs_err = gtsobs%each%height%error
         obs_qc(1)  = real(gtsobs%each%height%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, height_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

END subroutine insert_gts_obs_single

!-------------------------------------------------------------------------------
SUBROUTINE insert_gts_obs_multi(gtsobs, qc_threshold, &
                          pressure_type, height_type, u_wind_type,v_wind_type, &
                          temperature_type, dew_point_type, rh_type, & 
                          which_vert, num_obs, obs, seq)
!-------------------------------------------------------------------------------

   implicit none

   type(multi_level_type), intent(in)     :: gtsobs
   integer, optional, intent(in)          :: pressure_type, height_type, &
                                             u_wind_type, v_wind_type, &
                                             temperature_type, dew_point_type, &
                                             rh_type, qc_threshold
   integer, intent(in)                    :: which_vert
   integer, intent(inout)                 :: num_obs
   type(obs_type), intent(inout)          :: obs
   type(obs_sequence_type), intent(inout) :: seq

   type(obs_def_type)      :: obs_def
   type(location_type)     :: location
   type(time_type)         :: time

   character(len=80)       :: dummy
   integer                 :: year, month, day, hours, minutes, seconds
   real(r8)                :: lat,lon,pressure, height

   integer                 :: i, levels
   real(r8)                :: windspd, winddir, uu, vv, tk, td, rh, qv
   real(r8)                :: obs_value(1), obs_err, obs_qc(1)

   logical                 :: valid


   ! time info
   read(gtsobs%info%date_char, '(I4,5(A1,I2))')  &
        year, dummy, month, dummy, day, dummy, hours, dummy, minutes, dummy, seconds
   time = set_date(year, month, day, hours, minutes, seconds)
   call set_obs_def_time(obs_def, time)

   ! number of vertical levels
   levels = gtsobs%info%levels

!------------------
   do i = 1, levels
!------------------

   ! location info
   lat = gtsobs%info%lat
   lon = gtsobs%info%lon
   ! Dart longitude from 0 to 360
   if(lon < 0.0_r8) lon = lon + 360.0_r8
   height = gtsobs%info%elv
   pressure = gtsobs%each(i)%pressure%data
   if (which_vert.eq.VERTISSURFACE .or. which_vert.eq.VERTISHEIGHT) then
      location = set_location(lon, lat, height, which_vert)
   else if (which_vert.eq.VERTISPRESSURE) then
      location = set_location(lon, lat, pressure, which_vert)
   end if
   call set_obs_def_location(obs_def, location)

   if (present(u_wind_type) .or. present(v_wind_type)) then
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each(i)%speed%qc .ge. 0) &
                   .and. (gtsobs%each(i)%speed%qc .le. qc_threshold) &
                   .and. (gtsobs%each(i)%direction%qc .ge.0) &
                   .and. (gtsobs%each(i)%direction%qc .le. qc_threshold)
      endif
      !convert wind speed and direction to zonal and meridional components
      windspd = gtsobs%each(i)%speed%data
      winddir = gtsobs%each(i)%direction%data
      call spddir2uv(windspd,winddir,uu,vv,valid)
      if (valid) then
         ! Insert obs into DART obs sequence
         obs_value(1) = uu
         obs_err = gtsobs%each(i)%speed%error
         obs_qc(1)  = real(max(gtsobs%each(i)%speed%qc, gtsobs%each(i)%direction%qc),kind=r8)
         call set_obs_def_type_of_obs(obs_def, u_wind_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
         obs_value(1) = vv
         num_obs = num_obs + 1
         call set_obs_def_type_of_obs(obs_def, v_wind_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(temperature_type)) then
      !temperature in K
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each(i)%temperature%qc .ge. 0) &
                   .and. (gtsobs%each(i)%temperature%qc .le. qc_threshold)
      endif
      tk = gtsobs%each(i)%temperature%data
      if ( valid .and. tk > 0.0_r8 .and. tk < 1.0e3_r8 ) then
         obs_value(1) = tk
         obs_err = gtsobs%each(i)%temperature%error
         obs_qc(1)  = real(gtsobs%each(i)%temperature%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, temperature_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(dew_point_type)) then
      !dew point in K
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each(i)%dew_point%qc .ge. 0) &
                   .and. (gtsobs%each(i)%dew_point%qc .le. qc_threshold)
      endif
      td = gtsobs%each(i)%dew_point%data
      if ( valid .and. td > 0.0_r8 .and. td < 1.0e3_r8 ) then
         obs_value(1) = td
         obs_err = gtsobs%each(i)%dew_point%error
         obs_qc(1)  = real(gtsobs%each(i)%dew_point%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, dew_point_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(rh_type)) then
      !relative humidity rh
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each(i)%rh%qc .ge. 0) &
                   .and. (gtsobs%each(i)%rh%qc .le. qc_threshold)
      endif
      rh = gtsobs%each(i)%rh%data
      if ( valid .and. rh >= 0.0_r8 .and. rh <= 1.2e2_r8 ) then
         obs_value(1) = rh
         obs_err = gtsobs%each(i)%rh%error
         obs_qc(1)  = real(gtsobs%each(i)%rh%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, rh_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(pressure_type)) then
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each(i)%pressure%qc .ge. 0) &
                   .and. (gtsobs%each(i)%pressure%qc .le. qc_threshold)
      endif
      pressure = gtsobs%each(i)%pressure%data
      if ( valid .and. pressure >= 0.0_r8 .and. pressure < 105.0e3_r8 ) then
         obs_value(1) = pressure
         obs_err = gtsobs%each(i)%pressure%error
         obs_qc(1)  = real(gtsobs%each(i)%pressure%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, pressure_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

   if (present(height_type)) then
      valid = .true.
      if (present(qc_threshold)) then
         if (qc_threshold .ge. 0 ) &
            valid =      (gtsobs%each(i)%height%qc .ge. 0) &
                   .and. (gtsobs%each(i)%height%qc .le. qc_threshold)
      endif
      height = gtsobs%each(i)%height%data
      if ( valid .and. height >= -100.0e3_r8 .and. height < 100.0e3_r8 ) then
         obs_value(1) = height
         obs_err = gtsobs%each(i)%height%error
         obs_qc(1)  = real(gtsobs%each(i)%height%qc,kind=r8)
         call set_obs_def_type_of_obs(obs_def, height_type)
         call set_obs_def_error_variance(obs_def, obs_err*obs_err)
         call set_obs_def(obs, obs_def)
         call set_obs_values(obs, obs_value)
         call set_qc(obs,obs_qc)
         num_obs = num_obs + 1
         call insert_obs_in_seq(seq, obs)
      endif
   endif

!---------------
   end do 
!---------------

END subroutine insert_gts_obs_multi

!-------------------------------------------------------------------------------
SUBROUTINE insert_gts_obs_gpsref(gtsobs, qc_threshold, &
                                 gpsref_type, gpsref_key, gpsref_form, isgps, &
                                 which_vert, num_obs, obs, seq)
!-------------------------------------------------------------------------------

   implicit none

   type(multi_level_type), intent(in)     :: gtsobs
   integer, optional, intent(in)          :: qc_threshold
   integer, intent(in)                    :: gpsref_type 
   integer, intent(in)                    :: which_vert, gpsref_form
   logical, intent(in)                    :: isgps
   integer, intent(inout)                 :: num_obs, gpsref_key
   type(obs_type), intent(inout)          :: obs
   type(obs_sequence_type), intent(inout) :: seq

   type(obs_def_type)      :: obs_def
   type(location_type)     :: location
   type(time_type)         :: time

   character(len=80)       :: dummy
   integer                 :: year, month, day, hours, minutes, seconds
   real(r8)                :: lat,lon,pressure, height

   character(len=6)        :: subset
   integer                 :: i, levels
   real(r8)                :: td, nx, ny, nz, ds, htop, rfict
   real(r8)                :: obs_value(1), obs_err, obs_qc(1)

   logical                 :: valid


   ! time info
   read(gtsobs%info%date_char, '(I4,5(A1,I2))')  &
        year, dummy, month, dummy, day, dummy, hours, dummy, minutes, dummy, seconds
   time = set_date(year, month, day, hours, minutes, seconds)
   call set_obs_def_time(obs_def, time)

   ! number of vertical levels
   levels = gtsobs%info%levels

!------------------
   do i = 1, levels
!------------------

   gpsref_key = gpsref_key+1

   !gpsref data stored in dew_point

   valid = .true.
   if (present(qc_threshold)) then
      if (qc_threshold .ge. 0 ) &
         valid =      (gtsobs%each(i)%dew_point%qc .ge. 0) &
                .and. (gtsobs%each(i)%dew_point%qc .le. qc_threshold)
   endif
   if (valid) then
      ! location info
      lat = gtsobs%each(i)%speed%data
      lon = gtsobs%each(i)%direction%data
      ! Dart longitude from 0 to 360
      if(lon < 0.0_r8) lon = lon + 360.0_r8
      pressure = gtsobs%each(i)%pressure%data
      height = gtsobs%each(i)%height%data
      if (which_vert.eq.VERTISSURFACE .or. which_vert.eq.VERTISHEIGHT) then
         location = set_location(lon, lat, height, which_vert)
      else if (which_vert.eq.VERTISPRESSURE) then
         location = set_location(lon, lat, pressure, which_vert)
      end if
      call set_obs_def_location(obs_def, location)

      if ( gpsref_form .eq. 1 ) then
         nx = 0.0
         ny = 0.0
         nz = 0.0
         ds = 0.0
         htop = 0.0
         rfict = 0.0
         subset = 'GPSREF'
      else
         subset = 'GPSEXC'
         print*,'Not implemented yet! Return!'
      endif

      call set_gpsro_ref(gpsref_key, nx, ny, nz, rfict, ds, htop, subset)

      obs_value(1) = gtsobs%each(i)%dew_point%data
      obs_err = gtsobs%each(i)%dew_point%error
      obs_qc(1)  = real(gtsobs%each(i)%dew_point%qc,kind=r8)
      call set_obs_def_type_of_obs(obs_def, gpsref_type)
      call set_obs_def_key(obs_def, gpsref_key)
      call set_obs_def_error_variance(obs_def, obs_err*obs_err)
      call set_obs_def(obs, obs_def)
      call set_obs_values(obs, obs_value)
      call set_qc(obs,obs_qc)
      num_obs = num_obs + 1
      call insert_obs_in_seq(seq, obs)

   endif
!---------------
   end do 
!---------------

END subroutine insert_gts_obs_gpsref

!--------------------------------------
subroutine spddir2uv(spd,dir,u,v,valid)
!--------------------------------------
   real(r8)  :: spd, dir, u, v, dir_rad
   logical,intent(inout)   :: valid
!  valid = .true.
   u = missing_r8
   v = missing_r8
   if ( spd /= 0.0_r8 .and. abs(dir) < 1.e5_r8 ) then
      dir_rad = dir*DEG2RAD
      u = -spd*sin(dir_rad)
      v = -spd*cos(dir_rad)
   else
      valid = .false.
   endif
   if ( abs(u) > 1.e3_r8 .or. abs(v) > 1.e3_r8 ) valid=.false.
   return
end subroutine spddir2uv

!--------------------------------------------------------------------------------
! copied from WRFSI/src/mod/module_diagnostic_vars.F
FUNCTION compute_vapor_mixing_ratio(t_k, p_pa, rh_pct, rh_liq_flag) RESULT(wvmr)
!--------------------------------------------------------------------------------

    ! Computes the water vapor mixing ratio (dimensionless) given
    ! temperature (K), pressure (Pa), relative humidity (%), and a flag
    ! that tells the function whether or not the RH being provided is
    ! with respect to liquid everywhere or with respect to liquid only
    ! at temperature > Tref_frz.  This formulation comes from a paper by
    ! Buck (1991)

    IMPLICIT NONE
    REAL(r8)    :: t_k
    REAL(r8)    :: p_pa
    REAL(r8)    :: rh_pct
    LOGICAL     :: rh_liq_flag
    REAL(r8)    :: wvmr

    REAL(r8)    :: e,es
    REAL(r8), parameter :: Tref_frz = 273.15

    !  Compute saturation vapor pressure using an approx.
    !  developed by Buck (1981)

    IF ((rh_liq_flag).OR.((.NOT.rh_liq_flag).AND.(t_k.GE.Tref_frz)))THEN

      ! Compute saturation vapor pressure wrt liquid
    
      es = 100.0 * 6.1121 * EXP( (17.502*t_k-4780.6713)/(t_k-32.18) )

    ELSE 
 
      ! Compute es wrt ice
      ! Factor of 100 converts mb to Pa
      es = 100.0 * 6.1115 * EXP ( (22.452*t_k-6132.7638)/(t_k-0.60) )

    ENDIF

    ! From McKee (1993 CSU AT650 class notes) the values above
    ! need to be multiplied by 1.004 at pressures greater than 800 mb

    IF (p_pa.GE.80000.) es = es * 1.004

    ! Compute e from RH

    e = rh_pct*.01 * es

    ! Compute wvmr from e and p_pa

    wvmr = 0.622 * e / (p_pa - e )
END FUNCTION compute_vapor_mixing_ratio

!--------------------------------------
SUBROUTINE print_gts_obs_single(gtsobs, &
                          pressure_type, height_type, u_wind_type,v_wind_type, &
                          temperature_type, dew_point_type, rh_type, &
                          which_vert, num_obs, obs, seq)
!--------------------------------------
   implicit none
   type(single_level_type), intent(in)    :: gtsobs
   integer, optional, intent(in)          :: pressure_type, height_type, &
                                             u_wind_type, v_wind_type, &
                                             temperature_type, dew_point_type, &
                                             rh_type
   integer, intent(in)                    :: which_vert
   integer, intent(inout)                 :: num_obs
   type(obs_type), intent(inout)          :: obs
   type(obs_sequence_type), intent(inout) :: seq

   print*, gtsobs%info%date_char
   print*, gtsobs%info%levels

   if (present(pressure_type)) print*,'pressure_type=',pressure_type
   if (present(height_type)) print*,'height_type=',height_type
   if (present(u_wind_type)) print*,'u_wind_type=',u_wind_type
   if (present(v_wind_type)) print*,'v_wind_type=',v_wind_type
   if (present(temperature_type)) print*,'temperature_type=',temperature_type
   if (present(dew_point_type)) print*,'dew_point_type=',dew_point_type
   if (present(rh_type)) print*,'rh_type=',rh_type

end subroutine print_gts_obs_single

END module gts_dart_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
