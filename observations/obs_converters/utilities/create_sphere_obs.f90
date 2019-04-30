! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! create synthetic observations evenly spaced on a sphere.
! this program creates observations without data - they have
! a location, type, error variance, a time. run the output file
! through create_fixed_network_seq if you want a time series of
! all these obs. run the output through perfect_model_obs to add data. 
! 
! alter as you wish to change obs types, density, etc.

program create_sphere_obs

use         types_mod, only : r8, missing_r8, pi, rad2deg
use      location_mod, only : VERTISPRESSURE, set_location
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read,         &
                              do_nml_file, do_nml_term, logfileunit, nmlfileunit
use  time_manager_mod, only : time_type, set_calendar_type, set_date, GREGORIAN, &
                              get_time
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, set_obs_def, &
                              set_copy_meta_data, set_qc_meta_data
use      obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, &
                              RADIOSONDE_TEMPERATURE
use      obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             get_obs_def_time, get_obs_def_location,           &
                             get_obs_def_type_of_obs, get_obs_def_error_variance,         &
                             set_obs_def_key
use obs_utilities_mod, only : add_obs_to_seq



! Input is in hectopascals, obs are in pascals
real(r8) :: obs_levels(15) = (/ 1000., 950., 900., 850., 800., 750., 700., 650., &
                                 600., 550., 500., 400., 300., 200., 150. /)

! Date information 
integer :: year = 2008
integer :: month = 10
integer :: day = 1
integer :: hour = 0

! Number of roughly evenly distributed points in horizontal
integer :: number_of_locations, iunit, io
real(r8) :: x, y, z, lat, lon, inc, off, r, phi, dlon, dlat
integer :: num_obs, ob, lvl, nlvls

type(time_type) :: time_obs, prev_time

type(obs_sequence_type) :: seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
logical :: first_obs

namelist /create_sphere_obs_nml/ number_of_locations


! start of executable code

call initialize_utilities('create_sphere_obs')

call find_namelist_in_file("input.nml", "create_sphere_obs_nml", iunit)
read(iunit, nml = create_sphere_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "create_sphere_obs_nml")

if (do_nml_file()) write(nmlfileunit, nml=create_sphere_obs_nml)
if (do_nml_term()) write(     *     , nml=create_sphere_obs_nml)

nlvls = size(obs_levels)

! Total number of observations at single time is levels*n*3
num_obs = nlvls * number_of_locations * 3

! create a new (empty) sequence
! we define only the location, type, time, and error here.
! run perfect_model_obs to fill in the obs values if needed
call init_obs_sequence(seq, 0, 0, num_obs)

call init_obs(obs, 0, 0)
call init_obs(prev_obs, 0, 0)

first_obs = .true.

call set_calendar_type(GREGORIAN)
time_obs = set_date(year, month, day, hour, 0, 0)
call set_obs_def_time(obs_def, time_obs)


inc = pi * (3.0_r8 - sqrt(5.0_r8))
off = 2.0_r8 / number_of_locations

do ob = 1, number_of_locations

   y = (ob-1) * off - 1 + (off / 2)
   r = sqrt(1 - y*y)
   phi = (ob-1) * inc

   x = cos(phi) * r
   z = sin(phi) * r

   ! Now convert to latitude and longitude, radians
   ! first, then degrees

   lon = atan2(y, x) + pi
   lat = asin(z)
   
   dlon = rad2deg * lon
   dlat = rad2deg * lat

   do lvl = 1, nlvls
      call set_obs_def_location(obs_def, set_location(dlon, dlat, obs_levels(lvl), VERTISPRESSURE))
   
      call set_obs_def_error_variance(obs_def, 1.0_r8)
      call set_obs_def_type_of_obs(obs_def, RADIOSONDE_TEMPERATURE)
      call set_obs_def(obs, obs_def)  
      call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
      call set_obs_def_error_variance(obs_def, 4.0_r8)
      call set_obs_def_type_of_obs(obs_def, RADIOSONDE_U_WIND_COMPONENT)
      call set_obs_def(obs, obs_def) 
      call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
      call set_obs_def_type_of_obs(obs_def, RADIOSONDE_V_WIND_COMPONENT)
      call set_obs_def(obs, obs_def) 
      call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
   enddo

enddo

call write_obs_seq(seq, 'even_sphere.in')

call finalize_utilities('create_sphere_obs')

end program create_sphere_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

