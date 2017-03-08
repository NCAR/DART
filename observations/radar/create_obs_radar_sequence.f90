! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program create_obs_radar_sequence

! This program creates a group of observations as would be returned from a
! radar observation system. The radar location, sweep angle, elevations,
! etc are user inputs.  The output are observation locations and types
! without time information or data.  Either radial velocity, reflectivity,
! or both types of observations can be generated.
!
! To add time, run 'create_fixed_network_seq'.
! To add data, run 'perfect_model_obs'.
!
! This program creates an observation at each possible location, whether there
! is reflectivity there or not.  To cull clear-air observations, some other
! tool would have to be run after 'perfect_model_obs' has added observation
! data values.
!
! If more complex timestepping than 'create_fixed_network_seq' can add is
! wanted, look in the code below for the string 'TIME:'.  There are comments 
! there about how to add this function.  It will require adding new code
! into this program and recompiling it.
 
 
use         types_mod, only : r8, deg2rad, rad2deg, earth_radius, missing_r8
use     utilities_mod, only : timestamp, register_module, finalize_utilities, &
                              open_file, close_file, error_handler, E_ERR, &
                              initialize_utilities, register_module, E_MSG
use  time_manager_mod, only : set_time
use      location_mod, only : location_type, interactive_location, &
                              query_location, set_location, VERTISHEIGHT
use       obs_def_mod, only : obs_def_type, set_obs_def_type_of_obs,           &
                              set_obs_def_location, set_obs_def_time,   &
                              set_obs_def_error_variance, set_obs_def_key
use      obs_kind_mod, only : QTY_RADAR_REFLECTIVITY, QTY_VELOCITY,   &
                              DOPPLER_RADIAL_VELOCITY, RADAR_REFLECTIVITY
use obs_def_radar_mod, only : set_radial_vel
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs,         &
                              write_obs_seq, init_obs_sequence, set_obs_def, &
                              append_obs_to_seq, static_init_obs_sequence,   &
                              destroy_obs, destroy_obs_sequence


implicit none

! version controlled file description for error handling, do not edit.
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq
character(len = 129)    :: file_name
  
! Fixed parameters for WSR-88D radar elevations.
integer,  parameter :: n_elev_clear = 9
real(r8), parameter :: elev_clear(n_elev_clear)  = &
                          (/ 0.5_r8,  1.5_r8,  2.4_r8,  3.4_r8,  4.3_r8, &
                             6.0_r8,  9.9_r8, 14.6_r8, 19.5_r8           /)
integer,  parameter :: n_elev_storm = 14
real(r8), parameter :: elev_storm(n_elev_storm) =  &
                          (/ 0.5_r8,  1.5_r8,  2.4_r8,  3.4_r8,  4.3_r8, &
                             5.3_r8,  6.2_r8,  7.5_r8,  8.7_r8, 10.0_r8, &
                            12.0_r8, 14.0_r8, 16.7_r8, 19.5_r8           /)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! start of executable program code
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

! Record the current time, date, etc. to the logfile
call initialize_utilities('create_obs_radar_sequence')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Create the observation sequence
seq = simulate_radar_sequence()

! Write the sequence to a file
write(*, *) 'Input filename for sequence (  set_def.out   usually works well)'
read(*, *) file_name
call write_obs_seq(seq, file_name)

call destroy_obs_sequence(seq)

call error_handler(E_MSG,'create_obs_radar_sequence','Finished successfully.',&
                   source,revision,revdate)
call finalize_utilities()


contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! start of supporting subroutines
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------

function simulate_radar_sequence()

! Main routine that loops over radar locations, generating a series of
! observation locations based on user input.

type(obs_sequence_type) :: simulate_radar_sequence
type(obs_type)          :: obs
type(obs_def_type)      :: obs_def
type(location_type)     :: radar_location, obs_location

real(r8)                :: first_az, last_az, delta_az, beam_az
real(r8)                :: first_gate, last_gate, delta_gate, beam_gate
real(r8)                :: beam_elev, ceiling
real(r8), allocatable   :: elev(:), azim(:), gate(:)
real(r8)                :: radial_vel_variance, reflectivity_variance

integer                 :: i_elev, num_elev, i_az, num_az
integer                 :: i_gate, num_gate, i_radar, num_radar
integer                 :: i, max_num_obs, select_obstype
integer                 :: num_copies, num_qc, defkey, elev_type, days, secs
integer                 :: which_vert, rc
logical                 :: qc


! This program only generates the observation locations, sets the variance,
! and sets the obs types. It adds no quality control or data values, by design,
! since separate programs already exist to add them. (See top of this file for
! more info.)
num_copies = 0
num_qc = 0

! TIME: This program does not attach any time information to the observations.
! It sets all times to (0,0) and depends on a subsequent program (generally
! 'create_fixed_network_seq') to set the calendar date/times on the obs.
! However if more complex time stepping (e.g. by pairs of observations if
! generating both types, or by scan elevation, or all obs for one cone at the
! same time) then set the days and seconds as desired someplace in this code
! inside the appropriate loops but before the two calls to fillin_obsdef().
! The code will also probably need to prompt the user for initial and final
! times, or delta times, and probably for the actual calendar day of the
! observation, since at the time perfect_model_obs is run, real date
! information is required.  This will require setting a default calendar type.
days = 0
secs = 0

! Max number of obs to allocate space for, total of all obs.
write(*, *) 'Input upper bound on number of observations in sequence (e.g. 1000000)'
read(*, *) max_num_obs

! Initialize an obs_sequence structure
call init_obs_sequence(simulate_radar_sequence, num_copies, num_qc, max_num_obs)

! Initialize an obs variable
call init_obs(obs, num_copies, num_qc)

! The following data applies to all radar observations generated.
! The actual location of the radar will be the only difference between
! each separate group of observations.
write(*, *) 'The following information will apply to all radar sweeps'
write(*, *) 'generated.  Multiple radar locations can be specified but they'
write(*, *) 'all share the other characteristics.'

! Ask user which of the two (or both) observation types to generate at
! each location as the radar sweeps around.
select_obstype = 0
do while (select_obstype < 1 .or. select_obstype > 3)
   write(*,*) 'What type of observations do you want to generate?'
   write(*,*) ' 1=Doppler radial velocity, 2=reflectivity, 3=both'
   read(*,*) select_obstype
end do

! beam elevations
elev_type = -1
do while (elev_type < 1 .or. elev_type > 3)
   write(*, *) 'Input 1 for  9 pre-defined WSR-88D elevations (clear mode)'
   write(*, *) 'Input 2 for 14 pre-defined WSR-88D elevations (storm mode)'
   write(*, *) 'Input 3 to choose number and elevation angles by hand'
   read(*, *) elev_type
end do

! figure out how many
selectcase(elev_type)
   case(1)
      num_elev = n_elev_clear
   case(2)
      num_elev = n_elev_storm
   case(3)
      num_elev = 0
      do while(num_elev < 1)
         write(*, *) 'Input number of elevations'
         read(*, *) num_elev
      end do
end select

! make space
allocate(elev(num_elev), stat=rc)
if (rc /= 0) then
   call error_handler(E_ERR,'simulate_radar_sequence',    &
                            'Allocation error on elev() array', &
                             source, revision, revdate)
endif

! fill in values
selectcase(elev_type)
   case(1)
      elev(:) = elev_clear(:)
   case(2)
      elev(:) = elev_storm(:)
   case(3)
      do i = 1, num_elev
         write(*, "(a,i3)") 'Input elevation angle (degrees) # ',i
         read(*, *) elev(i)
      end do
end select

! convert to radians
elev = deg2rad*elev

! beam azimuth limits
write(*, *) 'Input first azimuth angle (degrees): example value = 0.'
read(*, *) first_az
first_az = deg2rad*first_az

write(*, *) 'Input last azimuth angle (degrees): example value = 359.'
read(*, *) last_az
last_az = deg2rad*last_az

write(*, *) 'Input azimuth angle increment (degrees): example value = 1.'
read(*, *) delta_az
delta_az = deg2rad*delta_az

num_az = int((last_az - first_az)/delta_az) + 1

allocate(azim(num_az), stat=rc)
if (rc /= 0) call error_handler(E_ERR,'simulate_radar_sequence',   &
                                'Allocation error on azim() array', &
                                 source, revision, revdate)
do i_az = 1, num_az
   azim(i_az) = first_az + (i_az-1)*delta_az
enddo

! beam gate limits
write(*, *) 'Input closest gate (m): example value = 3000.'
read(*, *) first_gate

write(*, *) 'Input farthest gate (m): example value = 100000.'
read(*, *) last_gate

write(*, *) 'Input gate length (m): example value = 1000.'
read(*, *) delta_gate
   
num_gate = int((last_gate - first_gate)/delta_gate) + 1

allocate(gate(num_gate), stat=rc)
if (rc /= 0) call error_handler(E_ERR,'simulate_radar_sequence',   &
                                'Allocation error on gate() array', &
                                 source, revision, revdate)
do i_gate = 1, num_gate
   gate(i_gate) = first_gate + (i_gate-1)*delta_gate
enddo

! do not generate observations above this value even if beam sweeps
! up higher than this value.
write(*, *) 'Input highest allowed observation height (m): example value 8500.'
read(*, *) ceiling
   
num_gate = int((last_gate - first_gate)/delta_gate) + 1

! error variances, per obs type
if (select_obstype == 1 .or. select_obstype == 3) then
   write(*, *) 'Input error variance for Doppler velocity (m/s)^2: example value = 4.'
   read(*, *) radial_vel_variance
endif

if (select_obstype == 2 .or. select_obstype == 3) then
   write(*, *) 'Input error variance for reflectivity (dBZ)^2: example value = 4.'
   read(*, *) reflectivity_variance
endif


! All radar locations share the above values
write(*,*) 'How many radar locations do you want to specify?'
read(*,*) num_radar

write(*, *) 'When setting the location of the radar, the vertical type'
write(*, *) '_must_ be specified in height (meters).'

RADARS : do i_radar = 1, num_radar
   write(*, *) ''
   write(*, "(a,i3)") 'Input radar location for radar number ', i_radar

   call interactive_location(radar_location)

   which_vert = nint(query_location(radar_location,'which_vert'))
   if(which_vert /= VERTISHEIGHT) then
      call error_handler(E_ERR,'simulate_radar_sequence',              &
           'Vertical coordinate of the radar location must be height', &
           source, revision, revdate)
   endif


   ELEVS : do i_elev = 1, num_elev
      beam_elev = elev(i_elev)
      AZIMUTHS : do i_az = 1, num_az
         beam_az = azim(i_az)
         GATES : do i_gate = 1, num_gate
            beam_gate = gate(i_gate)

            ! Compute the 3d location of this observation
            call compute_location(radar_location,                   &
                                  beam_gate, beam_az, beam_elev,    &
                                  ceiling, select_obstype, defkey,  &
                                  obs_location, qc)

            ! Bad quality control; most likely above ceiling setting.
            if(.not. qc) cycle GATES
          
            ! Generate a radial velocity obs
            if (select_obstype == 1 .or. select_obstype == 3) then
               call fillin_obsdef(DOPPLER_RADIAL_VELOCITY, obs_location, &
                                  days, secs, radial_vel_variance,       &
                                  defkey, obs_def)

               call set_obs_def(obs, obs_def)

               call append_obs_to_seq(simulate_radar_sequence, obs)

            endif
            ! Generate a reflectivity obs
            if (select_obstype == 2 .or. select_obstype == 3) then
                call fillin_obsdef(RADAR_REFLECTIVITY, obs_location,  &
                                   days, secs, reflectivity_variance, &
                                   0, obs_def)

                call set_obs_def(obs, obs_def)

                call append_obs_to_seq(simulate_radar_sequence, obs)

            endif

         end do GATES
      end do AZIMUTHS
   end do ELEVS
end do RADARS

deallocate(gate, azim, elev)
call destroy_obs(obs)

end function simulate_radar_sequence



!----------------------------------------------------------------------

subroutine compute_location(radar_location, beam_gate, beam_az, beam_elev, &
                            ceiling, select_obstype, defkey, obs_location, qc)

! Compute an actual observation location, given all the input arguments.

type(location_type),    intent(in)    :: radar_location
real(r8),               intent(in)    :: beam_gate, beam_az, beam_elev, ceiling
integer,                intent(in)    :: select_obstype
integer,                intent(inout) :: defkey
type(location_type),    intent(out)   :: obs_location 
logical,                intent(out)   :: qc           ! quality control: 0=good

real(r8) :: h, spath, x, y, radar_lon, radar_lat, obs_lon, obs_lat, vloc
real(r8) :: nyquist_velocity, beam_direction(3), elev_obs 
real(r8) :: effect_radius_m, actual_radius_m
integer  :: which_vert


! Doviak & Zrnic, 1993: Doppler radar and weather observations, eq. 2.28b-c

! Effective earth radius (4/3 times actual radius), converted to meters.
effect_radius_m = 4.0_r8 * earth_radius / 3.0_r8 * 1000.0_r8 

h = sqrt(beam_gate*beam_gate + effect_radius_m*effect_radius_m  &
           + 2.0_r8*beam_gate*effect_radius_m*sin(beam_elev)) - effect_radius_m
spath = effect_radius_m * asin(beam_gate*cos(beam_elev)/(effect_radius_m + h))

x = spath*sin(beam_az)
y = spath*cos(beam_az)

! The query routines return locations in radians.  
radar_lon = query_location(radar_location, 'lon')
radar_lat = query_location(radar_location, 'lat')

! Actual earth radius in meters
actual_radius_m = earth_radius * 1000.0_r8

obs_lat = y/actual_radius_m + radar_lat
obs_lon = x/(actual_radius_m*cos(radar_lat + y/(2.0_r8*actual_radius_m))) + &
          radar_lon

! Already verified to be height in meters.
vloc = query_location(radar_location, 'vloc')

! Add the height increment and check for 'too high'.
vloc = vloc + h

if(vloc.gt.ceiling) then       ! max. height for radar obs = 8.5km
   qc = .false.        ! bad data
   return
endif 

! Convert back to degrees (required by set_location() call).
obs_lon = obs_lon*rad2deg
obs_lat = obs_lat*rad2deg
which_vert = nint(query_location(radar_location,'which_vert'))

obs_location = set_location(obs_lon, obs_lat, vloc, which_vert)


! The elevation angle at the observation location is based on p. 23 of 
! Battan's radar meteorology book (equation 3.18), assuming an effective 
! radius of the earth of 4a/3.

elev_obs = sqrt(1.5_r8 * h / actual_radius_m + beam_elev*beam_elev)

! The radial velocity obs have additional metadata associated with them
! that must be set here.
if(select_obstype .eq. 1 .or. select_obstype .eq. 3) then
   beam_direction(1) = sin(beam_az)*cos(elev_obs)
   beam_direction(2) = cos(beam_az)*cos(elev_obs)
   beam_direction(3) = sin(elev_obs)

   nyquist_velocity = missing_r8

   call set_radial_vel(defkey, radar_location, beam_direction, nyquist_velocity)
endif

qc = .true.    ! all good

end subroutine compute_location

!----------------------------------------------------------------------

subroutine fillin_obsdef(obskind, obsloc, days, secs, variance, defkey, obs_def)
 
! Set the values in a DART obs_def_type derived type.

integer,             intent(in)    :: obskind
type(location_type), intent(in)    :: obsloc
integer,             intent(in)    :: days
integer,             intent(in)    :: secs
real(r8),            intent(in)    :: variance
integer,             intent(in)    :: defkey
type(obs_def_type),  intent(inout) :: obs_def

! set each part of the obs def
call set_obs_def_type_of_obs(obs_def, obskind)

call set_obs_def_location(obs_def, obsloc)

! Note that set_time wants seconds first, days second.
call set_obs_def_time(obs_def, set_time(secs, days) )

call set_obs_def_error_variance(obs_def, variance)

! The key is required in the radial velocity obs; it is unused in the
! reflectivity obs.  (The key is needed only for obs types which have
! auxilliary data associated with each obs.)
call set_obs_def_key(obs_def, defkey)

end subroutine fillin_obsdef

!----------------------------------------------------------------------
!----------------------------------------------------------------------

end program create_obs_radar_sequence

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
