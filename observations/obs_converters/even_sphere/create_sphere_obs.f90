! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! Create synthetic observations evenly spaced on a sphere.
!
! This program creates observations without data.  They have
! a location, type, error variance, and a time. Run the output file
! through create_fixed_network_seq if you want to create a time series 
! of these obs. Run the output through perfect_model_obs to add data. 
! 
! Alter as you wish to change obs types, number of types at each loc, etc.
! If you change the obs types either set the obs error explicitly
! or change the call to the obs_error module to match the correct obs type.
!
! There is a related matlab script in the obs_converters/even_sphere directory
! which uses the same algorithm to distribute N points around a sphere.  It has
! additional plotting functions which may be useful.


program create_sphere_obs

use         types_mod, only : r8, pi, rad2deg
use      location_mod, only : VERTISPRESSURE, set_location
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read,         &
                              do_nml_file, do_nml_term, logfileunit, nmlfileunit
use  time_manager_mod, only : time_type, set_calendar_type, set_date, GREGORIAN
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs, write_obs_seq, &
                              init_obs_sequence, set_obs_def
use      obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, &
                              RADIOSONDE_TEMPERATURE
use       obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                              set_obs_def_error_variance, set_obs_def_location
use obs_utilities_mod, only : add_obs_to_seq
use       obs_err_mod, only : rawin_temp_error, rawin_wind_error



integer  :: iunit, io, ntypes
real(r8) :: x, y, z, lat, lon, inc, off, r, phi, deglon, deglat, t_err, w_err
integer  :: num_obs, ob, lvl, nlvls

integer, parameter :: MAX_LEVELS = 40
integer, parameter :: DEF_LEVELS = 21

type(time_type) :: time_obs, prev_time

type(obs_sequence_type) :: seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
logical                 :: first_obs

! These default levels are standard radiosonde reporting pressures in hectopascals.
real(r8), parameter :: default_levels(DEF_LEVELS) = (/ 1000.0_r8, 925.0_r8, 850.0_r8, 700.0_r8, 500.0_r8, &
                                                        400.0_r8, 300.0_r8, 250.0_r8, 200.0_r8, 150.0_r8, &
                                                        100.0_r8,  70.0_r8,  50.0_r8,  30.0_r8,  20.0_r8, &
                                                         10.0_r8,  7.0_r8, 5.0_r8, 3.0_r8, 2.0_r8, 1.0_r8 /)

! namelist variables

! Number of roughly evenly distributed points in horizontal.
! Final number of obs will be N * L times this, where N is the
! number of different obs types at each location and L is the
! number of pressure levels specified.

integer :: number_of_locations

! The level array here is specified in hectopascals (mb).  
! Observations are created and stored in pascals.  
! If this item is not altered by the namelist the default levels are used.

real(r8) :: obs_levels(MAX_LEVELS) =  -1.0_r8

! Date of obs
integer :: year  = 2008
integer :: month = 10
integer :: day   = 1
integer :: hour  = 0


namelist /create_sphere_obs_nml/ number_of_locations, &
                                 obs_levels, &
                                 year, month, day, hour



! start of executable code

call initialize_utilities('create_sphere_obs')

call find_namelist_in_file("input.nml", "create_sphere_obs_nml", iunit)
read(iunit, nml = create_sphere_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "create_sphere_obs_nml")

if (do_nml_file()) write(nmlfileunit, nml=create_sphere_obs_nml)
if (do_nml_term()) write(     *     , nml=create_sphere_obs_nml)

! if the user doesn't set any levels, use the defaults.
! otherwise, stop at the first -1 value or end of the array.

if (all(obs_levels <= 0.0_r8)) then
   obs_levels(1:DEF_LEVELS) = default_levels
   nlvls = DEF_LEVELS
else
   nlvls = findloc(obs_levels, -1.0_r8, 1) - 1
   if (nlvls <= 0) nlvls = size(obs_levels)
endif

! update this if you change the code to have more or fewer 
! obs types created at each location.
ntypes = 3   

! Total number of observations at single time is: nlevels * nlocs * ntypes
num_obs = nlvls * number_of_locations * ntypes

! create a new (empty) sequence.
! this program defines only the location, type, time, and error.
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
   
   deglon = rad2deg * lon
   deglat = rad2deg * lat

   do lvl = 1, nlvls

      ! in dart vertical locations using pressure are specified in units of pascals
      call set_obs_def_location(obs_def, set_location(deglon, deglat, obs_levels(lvl)*100.0_r8, VERTISPRESSURE))
   
      ! the standard obs error routines take vertical locations in hectopascals/mb
      t_err = rawin_temp_error(obs_levels(lvl))
      call set_obs_def_error_variance(obs_def, t_err*t_err)
      call set_obs_def_type_of_obs(obs_def, RADIOSONDE_TEMPERATURE)
      call set_obs_def(obs, obs_def)  
      call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
      w_err = rawin_wind_error(obs_levels(lvl))
      call set_obs_def_error_variance(obs_def, w_err*w_err)
      call set_obs_def_type_of_obs(obs_def, RADIOSONDE_U_WIND_COMPONENT)
      call set_obs_def(obs, obs_def) 
      call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
      call set_obs_def_error_variance(obs_def, w_err*w_err)
      call set_obs_def_type_of_obs(obs_def, RADIOSONDE_V_WIND_COMPONENT)
      call set_obs_def(obs, obs_def) 
      call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
   enddo

enddo

call write_obs_seq(seq, 'even_sphere_seq.in')

call finalize_utilities('create_sphere_obs')

end program create_sphere_obs

