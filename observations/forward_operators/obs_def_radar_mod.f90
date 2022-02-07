! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!-----------------------------------------------------------------------------
! DART radar observation module, including the observation operators for the
! two primary radar-observation types -- Doppler velocity and reflectivity --
! plus other utility subroutines and functions.  A number of simplifications
! are employed for the observation operators.  Most notably, the model state
! is mapped to a "point" observation, whereas a real radar observation is a
! volumetric sample.  The implications of this approximation have not been
! investigated fully, so in the future it might be worth developing and testing
! more sophisticated observation operators that produce volumetric power-
! weighted samples.
!
! This module is able to compute reflectivity and precipitation fall speed
! (needed for computing Doppler radial velocity) from the prognostic model
! fields only for simple single-moment microphysics schemes such as the Kessler
! and Lin schemes.  If a more complicated microphysics scheme is used, then
! reflectivity and fall speed must be accessible instead as diagnostic fields
! in the model state.
!
! Author and Contact information:
!   Radar Science:  David Dowell, ddowell at ucar.edu
!   DART Code:  Nancy Collins, nancy at ucar.edu
!   Original DART/Radar work:  Alain Caya
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! DOPPLER_RADIAL_VELOCITY, QTY_VELOCITY
! RADAR_REFLECTIVITY, QTY_RADAR_REFLECTIVITY
! RADAR_CLEARAIR_REFLECTIVITY, QTY_RADAR_REFLECTIVITY
! PRECIPITATION_FALL_SPEED, QTY_POWER_WEIGHTED_FALL_SPEED
! END DART PREPROCESS TYPE DEFINITIONS
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_radar_mod, only : write_radial_vel, read_radial_vel,           &
!                            interactive_radial_vel, get_expected_radial_vel, &
!                            read_radar_ref, get_expected_radar_ref,          &
!                            get_expected_fall_velocity
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(DOPPLER_RADIAL_VELOCITY)
!     call get_expected_radial_vel(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!  case(RADAR_REFLECTIVITY)
!     call get_expected_radar_ref(state_handle, ens_size, location, expected_obs, istatus)
!  case(RADAR_CLEARAIR_REFLECTIVITY)
!     call get_expected_radar_ref(state_handle, ens_size, location, expected_obs, istatus)
!  case(PRECIPITATION_FALL_SPEED)
!     call get_expected_fall_velocity(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(DOPPLER_RADIAL_VELOCITY)
!      call read_radial_vel(obs_def%key, ifile, fform)
!   case(RADAR_REFLECTIVITY)
!      call read_radar_ref(obs_val, obs_def%key)
!  case(RADAR_CLEARAIR_REFLECTIVITY)
!      call read_radar_ref(obs_val, obs_def%key)
!   case(PRECIPITATION_FALL_SPEED)
!      continue
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(DOPPLER_RADIAL_VELOCITY)
!      call write_radial_vel(obs_def%key, ifile, fform)
!   case(RADAR_REFLECTIVITY)
!      continue
!   case(RADAR_CLEARAIR_REFLECTIVITY)
!      continue
!   case(PRECIPITATION_FALL_SPEED)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(DOPPLER_RADIAL_VELOCITY)
!      call interactive_radial_vel(obs_def%key)
!   case(RADAR_REFLECTIVITY)
!      continue
!   case(PRECIPITATION_FALL_SPEED)
!      continue
!   case(RADAR_CLEARAIR_REFLECTIVITY)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_radar_mod

use        types_mod, only : r8, missing_r8, PI, deg2rad
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location, get_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, &
                             QTY_TEMPERATURE, QTY_VERTICAL_VELOCITY,     &
                             QTY_RAINWATER_MIXING_RATIO, QTY_DENSITY,    &
                             QTY_GRAUPEL_MIXING_RATIO,                    &
                             QTY_SNOW_MIXING_RATIO,                       &
                             QTY_POWER_WEIGHTED_FALL_SPEED,               &
                             QTY_RADAR_REFLECTIVITY

use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: read_radar_ref, get_expected_radar_ref,                          &
          read_radial_vel, write_radial_vel, interactive_radial_vel,       &
          get_expected_radial_vel, get_obs_def_radial_vel, set_radial_vel, &
          get_expected_fall_velocity

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical :: module_initialized = .false.

! Derived type for radial velocity.  Contains auxiliary information stored
! with each obs of this type; used to compute the forward operator.
! See more extensive comments in the interactive_radial_vel() routine for
! expected units, etc.  Technically, the radar location is unused in the
! forward operators currently in this code, but it may be useful for post
! processing or diagnostics, especially if multiple radar locations are
! in the same file.  
type radial_vel_type
   private
   type(location_type) :: radar_location      ! location of radar 
   real(r8)            :: beam_direction(3)   ! direction of beam
   real(r8)            :: nyquist_velocity    ! nyquist velocity
end type radial_vel_type

! Cumulative index into radial velocity metadata array
integer :: velkeycount = 0 

! For error message content
character(len=256) :: msgstring

! Values which are initialized at run time so some can be changed by
! namelist.  After initialization, treated as parameters (values not changed).

! Note that the value of gravity is hardcoded here.  The value for gravity 
! used in the model should match this value.   See additional comments 
! below in the initialize_constants() subroutine.

real(r8) :: param_gravity     ! gravitational acceleration (m s^-2)

                              ! empirical constants:
real(r8) :: param_a           ! 10^^18 * a (const in rain fall speed eqn)
real(r8) :: param_b           !          b (const in rain fall speed eqn)
real(r8) :: param_c           ! 10^^18 * c (const in snow fall speed eqn)
real(r8) :: param_d           !          d (const in snow fall speed eqn)

real(r8) :: param_CD          ! drag coefficient for graupel/hail
real(r8) :: param_rhos0       ! reference density (kg m^-3) in emperical
                              !  dropsize-fall speed eqn
real(r8) :: param_e           ! parameter in graupel/hail fall speed eqn

                              ! results of gamma function applied to:
real(r8) :: param_gam7b       !   (7+b)
real(r8) :: param_gam7d       !   (7+d)
real(r8) :: param_gam7f       !   (7+0.5)

                              ! exponents in equations for:
real(r8) :: param_powr        !   rain fall speed
real(r8) :: param_pows        !   snow fall speed 
real(r8) :: param_powg_dry    !   dry graupel/hail fall speed
real(r8) :: param_powg_wet    !   wet graupel/hail fall speed

                              ! parameters in the equations for:
real(r8) :: param_fs_r        !   rain fall speed
real(r8) :: param_fs_wet_s    !   wet snow fall speed
real(r8) :: param_fs_dry_s    !   dry snow fall speed
real(r8) :: param_fs_wet_g    !   wet graupel/hail fall speed
real(r8) :: param_fs_dry_g    !   dry graupel/hail fall speed

                              ! parameters in the equations for:
real(r8) :: param_refl_r      !   reflectivity from rain
real(r8) :: param_refl_wet_s  !   reflectivity from wet snow
real(r8) :: param_refl_dry_s  !   reflectivity from dry snow
real(r8) :: param_refl_wet_g  !   reflectivity from wet graupel/hail
real(r8) :: param_refl_dry_g  !   reflectivity from dry graupel/hail

!-------------------------------------------------------------
! Namelist with default values
! 
! Obsolete: convert_to_dbz and dbz_threshold
!  convert_to_dbz and dbz_threshold have both been removed from the namelist.
!  Values will always be converted to dBZ, and threshold was only used to 
!  ensure the log() call never saw a real 0.0_r8.  Please remove these 
!  values from your namelist to avoid a run-time error. 
!
! There are two replicated sets of 3 namelist values below.
! In each case, there is 1 logical value and 2 numeric values.
! If the logical is false, the numeric values are ignored.
! If true, then the 2 numeric values are:
!  1) a threshold value.  If the observation or forward operator value is
!     less than this threshold, it will be set to a different value.
!  2) the value it should be set to.
! The value is separate from the threshold to allow, for example, the option
! of setting all values below -20 dBZ to -40 dBZ.
!
! If the observation value or the forward operator already has a value of
! missing_r8 it is assumed either the istatus is marked as failed (for the
! forward operator) or that the QC (quality control) flag is set to not
! assimilate this observation and that value is left unchanged regardless of
! the setting on the apply_ref_limit flag.  Note however that it is not a
! good idea to reset a good but small observation value to missing_r8 -- do
! not use it as the lowest_reflectivity setting.
!
! The next 3 namelist items apply to the incoming observation values.  
! They are in the namelist so they can be changed at runtime, instead 
! of set only when the observation file is originally generated. 
!
! apply_ref_limit_to_obs:    
!  Logical.  If .TRUE. replace all reflectivity values less than 
!  "reflectivity_limit_obs" with "lowest_reflectivity_obs" value. 
!  If .FALSE. leave all values as-is.  Defaults to .FALSE.
!
! reflectivity_limit_obs:    
!  The threshold value.  Observed reflectivity values less than this
!  threshold will be set to the "lowest_reflectivity_obs" value.
!  Units are dBZ.  Defaults to -10.0.
!
! lowest_reflectivity_obs:   
!  The 'set-to' value.  Observed reflectivity values less than the
!  threshold will be set to this value.  Units are dBZ.  Defaults to -10.0.
!
! The next 3 namelist items apply to the forward operator values (the returned
! value from each ensemble member predicting what the observation value should
! be given the current state in this particular member).
!
! apply_ref_limit_to_fwd_op:  
!  Same as "apply_ref_limit_to_obs", but for the forward operator.
!
! reflectivity_limit_fwd_op:  
!  Same as "reflectivity_limit_obs", but for the forward operator values.
!
! lowest_reflectivity_fwd_op: 
!  Same as "lowest_reflectivity_obs", but for the forward operator values.
!
! max_radial_vel_obs:
!  Integer value.  Maximum number of observations of this type to support at
!  run time.  This is combined total of all obs_seq files, for example the
!  observation diagnostic program potentially opens multiple obs_seq.final
!  files, or the obs merge program can also open multiple obs files.
!  Default is 1,000,000 observations.
!
! dielectric_factor:
!  According to Smith (1984), there are two choices for the dielectric
!  factor, depending on how the snow particle sizes are specified.
!  If melted raindrop diameters are used, then the factor is 0.224.  If
!  equivalent ice sphere diameters are used, then the factor is 0.189.
!  Since the common convention is to use melted raindrop diameters, the
!  default here is 0.224.
!
! n0_rain, n0_graupel, n0_snow:
!  Intercept parameters (m^-4) for size distributions of each hydrometeor.
!  The defaults below are for the Lin et al. microphysics scheme
!  with the Hobbs settings for graupel/hail.  (The Hobbs graupel settings
!  are also the default for the Lin scheme in WRF 2.2 and 3.0.)
!
! rho_rain, rho_graupel, rho_snow:
!  Density (kg m^-3) of each hydrometeor type.  The defaults below are for the
!  Lin et al. microphysics scheme with the Hobbs setting for graupel/hail.
!
! allow_wet_graupel:
!  Logical.  It is difficult to predict/diagnose whether graupel/hail has a
!  wet or dry surface.  Even when the temperature is above freezing,
!  evaporation and/or absorption can still result in a dry surface.  This
!  issue is important because the reflectivity from graupel with a wet surface
!  is significantly greater than that from graupel with a dry surface.
!  Currently, the user has two options for how to compute graupel
!  reflectivity.  If allow_wet_graupel is .false. (the default), then graupel
!  is always assumed to be dry.  If allow_wet_graupel is .true., then graupel
!  is assumed to be wet (dry) when the temperature is above (below) freezing.
!  A consequence is that a sharp gradient in reflectivity will be produced at
!  the freezing level.  In the future, it might be better to provide the
!  option of having a transition layer.
!
! microphysics_type:
!  Integer. Tells obs_def_radar_mod what microphysical scheme is employed
!  to enable some smarter decisions in handling fall velocity and radar
!  reflectivity.
!   1  - Kessler scheme
!   2  - Lin et al. microphysics  (default)
!   3  - User selected scheme where 10 cm reflectivity and power weighted fall
!        velocity are expected in the state vector (failure if not found)
!   4  - User selected scheme where only power weighted fall velocity is
!        expected (failure if not found)
!   5  - User selected scheme where only reflectivity is expected (failure
!        if not found)
!  -1  - ASSUME FALL VELOCITY IS ZERO, allows over-riding the failure modes
!        above if reflectivity and/or fall velocity are not available but a
!        result is desired for testing purposes only.
!
! allow_dbztowt_conv:
!  Logical. Flag to enable use of the dbztowt routine where reflectivity is
!  available, but not the power-weighted fall velocity. This scheme uses
!  emperical relations between reflectivity and fall velocity, with poor
!  accuracy for highly reflective, low density particles (such as water coated
!  snow aggregates). Expect questionable accuracy in radial velocity from
!  the forward operator with high elevation angles where ice is present in
!  the model state.
 
logical  :: apply_ref_limit_to_obs     = .false.
real(r8) :: reflectivity_limit_obs     = -10.0_r8
real(r8) :: lowest_reflectivity_obs    = -10.0_r8
logical  :: apply_ref_limit_to_fwd_op  = .false.
real(r8) :: reflectivity_limit_fwd_op  = -10.0_r8
real(r8) :: lowest_reflectivity_fwd_op = -10.0_r8
integer  :: max_radial_vel_obs         = 1000000
logical  :: allow_wet_graupel          = .false.
integer  :: microphysics_type          = 2
logical  :: allow_dbztowt_conv         = .false.


! Constants which depend on the microphysics scheme used.  Should be set
! in the namelist to match the case being run.

real(r8) :: dielectric_factor =  0.224_r8
real(r8) :: n0_rain           =  8.0e6_r8
real(r8) :: n0_graupel        =  4.0e6_r8
real(r8) :: n0_snow           =  3.0e6_r8
real(r8) :: rho_rain          = 1000.0_r8
real(r8) :: rho_graupel       =  400.0_r8
real(r8) :: rho_snow          =  100.0_r8


namelist /obs_def_radar_mod_nml/ apply_ref_limit_to_obs,     &
                                 reflectivity_limit_obs,     &
                                 lowest_reflectivity_obs,    &
                                 apply_ref_limit_to_fwd_op,  &
                                 reflectivity_limit_fwd_op,  &
                                 lowest_reflectivity_fwd_op, &
                                 max_radial_vel_obs,         &
                                 allow_wet_graupel,          &
                                 dielectric_factor,          &
                                 microphysics_type,          &
                                 allow_dbztowt_conv,         &
                                 n0_rain,                    &
                                 n0_graupel,                 &
                                 n0_snow,                    &
                                 rho_rain,                   &
                                 rho_graupel,                &
                                 rho_snow


! Module global storage for auxiliary obs data, allocated in init routine
type(radial_vel_type), allocatable :: radial_vel_data(:)

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine initialize_module

! Called once to set values and allocate space

integer :: iunit, io, rc

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_radar_mod_nml", iunit)
read(iunit, nml = obs_def_radar_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_radar_mod_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_radar_mod_nml)
if (do_nml_term()) write(     *     , nml=obs_def_radar_mod_nml)

! Consistency warning; print a message if the thresholds and lower values
! are going to be used and are different.
call check_namelist_limits(apply_ref_limit_to_obs,     &
                           reflectivity_limit_obs,     &
                           lowest_reflectivity_obs,    &
                           apply_ref_limit_to_fwd_op,  &
                           reflectivity_limit_fwd_op,  &
                           lowest_reflectivity_fwd_op)

! Allocate space for the auxiliary information associated with each obs
! This code must be placed after reading the namelist, so the user can
! increase or decrease the number of obs supported and use more or less
! memory at run time.
allocate(radial_vel_data(max_radial_vel_obs), stat = rc)
if (rc /= 0) then            
   write(msgstring, *) 'initial allocation failed for radial vel obs data,', &
                       'itemcount = ', max_radial_vel_obs
   call error_handler(E_ERR,'initialize_module', msgstring, &
                      source, revision, revdate)
endif                        

! Set the module global values that do not change during the run.
! This code uses some values which are set in the namelist, so this call
! *must* happen after the namelist read above.
call initialize_constants()

! Log the values used for the constants.
call print_constants()

end subroutine initialize_module

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Radial velocity section
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine read_radial_vel(velkey, ifile, fform)

! Main read subroutine for the radial velocity observation auxiliary data.

integer,          intent(out)          :: velkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=8)    :: header
logical             :: is_asciifile
type(location_type) :: radar_location
real(r8)            :: beam_direction(3)
real(r8)            :: nyquist_velocity
integer             :: oldkey

! radar_location: Refers to the lat/lon/height/vertical coordinate option
!   for the radar. Its type is defined through location_type in location_mod. 
!   Uses the read_location function in location_mod.
!
! beam_direction: This is a 3-element array specific to the radar module:
!     beam_direction(1) = sin(azimuth)*cos(elevation)
!     beam_direction(2) = cos(azimuth)*cos(elevation)
!     beam_direction(3) = sin(elevation)

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
      ! Read the character identifier for verbose formatted output
      read(ifile, FMT="(a8)") header
      if(header /= 'platform') then
         call error_handler(E_ERR,'read_radial_vel', &
              "Expected location header 'platform' in input file", &
              source, revision, revdate)
      endif
endif

! read_location is a DART library routine that expects an optional string
! arg, but the other two read routines are local to this module and we can
! tell them exactly what format to be reading because we already know it.
radar_location   = read_location        (ifile, fform)
beam_direction   = read_beam_direction  (ifile, is_asciifile)
nyquist_velocity = read_nyquist_velocity(ifile, is_asciifile)

! Read in the velkey for this particular observation, however, it will
! be discarded and a new, unique key will be generated in the set routine.
if (is_asciifile) then
   read(ifile, *) oldkey
else
   read(ifile) oldkey
endif

! Generate new unique radial velocity observation key, and set the contents
! of the private defined type.
call set_radial_vel(velkey, radar_location, beam_direction, nyquist_velocity)

end subroutine read_radial_vel

!----------------------------------------------------------------------

subroutine write_radial_vel(velkey, ifile, fform)

! Write radial velocity auxiliary information to the obs_seq file.

integer,          intent(in)           :: velkey, ifile
character(len=*), intent(in), optional :: fform

logical             :: is_asciifile
type(location_type) :: radar_location
real(r8)            :: beam_direction(3)
real(r8)            :: nyquist_velocity

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Write the 5 character identifier for verbose formatted output
   write(ifile, "('platform')")
endif

! Extract the values for this key and call the appropriate write routines.
radar_location    = radial_vel_data(velkey)%radar_location
beam_direction(:) = radial_vel_data(velkey)%beam_direction(:)
nyquist_velocity  = radial_vel_data(velkey)%nyquist_velocity

! write_location routine is part of the DART library and wants the optional
! format string argument.  The other two routines are local to this module, 
! and we have already figured out if it is a unformatted/binary file or 
! formatted/ascii, so go ahead and pass that info directly down to the routines.
call         write_location(ifile, radar_location,    fform) 
call   write_beam_direction(ifile, beam_direction(:), is_asciifile) 
call write_nyquist_velocity(ifile, nyquist_velocity,  is_asciifile) 

! Write out the velkey used for this run, however this will be discarded
! when this observation is read in and a new key will be generated.
! (It may be useful for correlating error messages or identifying particular
! observations so we are leaving it as part of the aux data.)
if (is_asciifile) then
   write(ifile, *) velkey
else
   write(ifile) velkey
endif

end subroutine write_radial_vel

!----------------------------------------------------------------------

function read_beam_direction(ifile, is_asciiformat)

! Reads beam_direction from file that was written by write_beam_direction.
! See read_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
real(r8)            :: read_beam_direction(3)

character(len=5)   :: header
real(r8)           :: beam_direction(3)

if ( .not. module_initialized ) call initialize_module


if (is_asciiformat) then
   read(ifile, "(a5)" ) header

   if(header /= 'dir3d') then
      write(msgstring,*)"Expected beam_direction header 'dir3d' in input file, got ", header
      call error_handler(E_ERR, 'read_beam_direction', msgstring, source, revision, revdate)
   endif
   ! Now read the beam_direction data value into temporaries
   read(ifile, *) beam_direction(1), beam_direction(2), beam_direction(3)
else
   ! No header label, just the binary direction values.
   read(ifile)    beam_direction(1), beam_direction(2), beam_direction(3)
endif

! Check for illegal values
if (minval(beam_direction) < -1.0_r8 .or. maxval(beam_direction) > 1.0_r8) then
   write(msgstring,*) "beam_direction value must be between -1 and 1, got: ", &
                       beam_direction(1), beam_direction(2), beam_direction(3)
   call error_handler(E_ERR, 'read_beam_direction', msgstring, &
                      source, revision, revdate)
endif

! set function return value
read_beam_direction(:) = beam_direction(:)

end function read_beam_direction

!----------------------------------------------------------------------

subroutine write_beam_direction(ifile, beam_direction, is_asciiformat)

! Writes beam_direction to obs file.

integer,  intent(in) :: ifile
real(r8), intent(in) :: beam_direction(3)
logical,  intent(in) :: is_asciiformat

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   write(ifile, "('dir3d')" ) 
   write(ifile, *) beam_direction(1), beam_direction(2), beam_direction(3)
else
   write(ifile)    beam_direction(1), beam_direction(2), beam_direction(3)
endif

end subroutine write_beam_direction

!----------------------------------------------------------------------

function read_nyquist_velocity(ifile, is_asciiformat)

! Reads Nyquist velocity from file that was written by write_nyquist_velocity.
! See read_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
real(r8)            :: read_nyquist_velocity

logical, save :: first_time = .true.

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   read(ifile, *) read_nyquist_velocity
else
   read(ifile)    read_nyquist_velocity
endif

! idealized obs leave the nyquist velocity as missing_r8.
! warn the user (once) but allow it.  they will not be unfolded
! at assimilation time, which is usually the desired behavior
! for an idealized ob.

if (read_nyquist_velocity == missing_r8 .and. first_time) then
   write(msgstring, *) "no doppler velocity unfolding can be done on these observations"
   call error_handler(E_MSG, 'read_nyquist_velocity:', &
                      "radar observation(s) with missing nyquist velocities encountered on read", &
                      source, revision, revdate, text2=msgstring)
 
   first_time = .false.
endif

! Check for illegal values; must be non-negative (missing is ok).
if (read_nyquist_velocity < 0.0_r8 .and. read_nyquist_velocity /= missing_r8) then
   write(msgstring,*) "bad value for nyquist velocity: ", read_nyquist_velocity
   call error_handler(E_ERR, 'read_nyquist_velocity:', msgstring, &
                      source, revision, revdate, text2="nyquist velocity must be non-negative")
endif

end function read_nyquist_velocity

!----------------------------------------------------------------------

subroutine write_nyquist_velocity(ifile, nyquist_velocity, is_asciiformat)

! Writes Nyquist velocity to obs file.

integer,  intent(in) :: ifile
real(r8), intent(in) :: nyquist_velocity
logical,  intent(in) :: is_asciiformat

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   write(ifile, *) nyquist_velocity
else
   write(ifile)    nyquist_velocity
endif

end subroutine write_nyquist_velocity

!----------------------------------------------------------------------

subroutine get_obs_def_radial_vel(velkey, radar_location, beam_direction, &
                                  nyquist_velocity)

! Return the auxiliary contents of a given radial velocity observation

integer,             intent(in)  :: velkey
type(location_type), intent(out) :: radar_location
real(r8),            intent(out) :: beam_direction(3)
real(r8),            intent(out) :: nyquist_velocity

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

radar_location    = radial_vel_data(velkey)%radar_location
beam_direction    = radial_vel_data(velkey)%beam_direction(:)
nyquist_velocity  = radial_vel_data(velkey)%nyquist_velocity

end subroutine get_obs_def_radial_vel

!----------------------------------------------------------------------

subroutine set_radial_vel(velkey, radar_location, beam_direction, nyquist_velocity)

! Common code to increment the current key count, and set the private
! contents of this observation's auxiliary data.

integer,             intent(out) :: velkey
type(location_type), intent(in)  :: radar_location
real(r8),            intent(in)  :: beam_direction(3)
real(r8),            intent(in)  :: nyquist_velocity

if ( .not. module_initialized ) call initialize_module

! The total velocity metadata key count from all sequences
velkeycount = velkeycount + 1    
velkey = velkeycount             ! set the return value

! Simple error check on key number before accessing the array
! This errors out if too key value now too large.
call velkey_out_of_range(velkey)

radial_vel_data(velkey)%radar_location    = radar_location
radial_vel_data(velkey)%beam_direction(:) = beam_direction(:)
radial_vel_data(velkey)%nyquist_velocity  = nyquist_velocity

end subroutine set_radial_vel

!----------------------------------------------------------------------

subroutine interactive_radial_vel(velkey)

! Interactively reads in auxiliary information for a radial velocity obs.

integer, intent(out) :: velkey

! Uses interactive_location of DART location_mod, plus the local subroutines
! interactive_beam_direction and set_radial_vel.
! See read_radial_vel for more information.

! velkey is internally incremented in the set routine, and only counts 
! the index for this specialized observation kind. 

type(location_type)  :: location
real(r8)             :: beam_direction(3)
real(r8)             :: nyquist_velocity

if ( .not. module_initialized ) call initialize_module

!Note: Obs location will subsequently be read in by the standard DART module
! "interactive_obs_def". No check is performed here whether radar location,
! obs location, and beam direction are self-consistent. It is the user's
! responsibility to make sure that they are.  This set of information
! does overspecify the problem slightly, but it is expensive to compute
! the true beam direction because of bending of the beam and earth curvature.


write(*, *)
write(*, *) 'Beginning to inquire for information on radar location.'
write(*, *)
write(*, *) 'WARNING!! Make sure that you select 3 (height) for the'
write(*, *) 'vertical coordinate option and enter height in meters.'
write(*, *) 'This location is where the radar source is located.  The later'
write(*, *) 'location question will be asking about where the observation'
write(*, *) 'itself is located.'
write(*, *)

call interactive_location(location)

write(*, *)
write(*, *) 'Beginning to inquire for information on radar beam direction.'
write(*, *)

call interactive_beam_direction(beam_direction)

write(*, *)
write(*, *) 'Beginning to inquire for information on radar Nyquist velocity.'
write(*, *)

call interactive_nyquist_velocity(nyquist_velocity)


call set_radial_vel(velkey, location, beam_direction, nyquist_velocity)

write(*, *)
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_radial_vel

!----------------------------------------------------------------------

subroutine interactive_beam_direction(beam_direction)

! Prompt for beam direction information in azimuth/elevation degrees.

real(r8), intent(out) :: beam_direction(3)

real(r8) :: az, el

az = -1.0
do while (az < 0.0 .or. az > 360.0) 
   write(*, *) 'Input the beam direction azimuth in degrees (0 <= az <= 360):'
   read(*, *) az
end do

el = -1.0
do while (el < 0.0 .or. el > 90.0)
   write(*, *) 'Input the beam direction elevation in degrees (0 <= el <= 90):'
   read(*, *) el
end do

! Convert to radians and compute the actual values stored with the observation.
az = az * deg2rad
el = el * deg2rad

beam_direction(1) = sin(az) * cos(el)
beam_direction(2) = cos(az) * cos(el)
beam_direction(3) = sin(el)

end subroutine interactive_beam_direction

!----------------------------------------------------------------------

subroutine interactive_nyquist_velocity(nyquist_velocity)

! Prompt for Nyquist velocity

real(r8), intent(out) :: nyquist_velocity

nyquist_velocity = -1.0

do while (nyquist_velocity < 0.0)
   write(*, *) 'Input Nyquist velocity for this obs point in m/sec'
   write(*, *) '(Typical values are 10-100 m/s, must be >= 0):'
   read(*, *)  nyquist_velocity
end do

end subroutine interactive_nyquist_velocity

!----------------------------------------------------------------------

subroutine get_expected_radial_vel(state_handle, ens_size, location, velkey, &
                                   radial_vel, istatus)

! This is the main forward operator routine for radar Doppler velocity.

type(ensemble_type),    intent(in) :: state_handle
integer,                intent(in) :: ens_size
type(location_type),    intent(in) :: location
integer,                intent(in) :: velkey
real(r8),              intent(out) :: radial_vel(ens_size)
integer,               intent(out) :: istatus(ens_size)

! Given a location and the state vector from one of the ensemble members,
! compute the model-predicted radial velocity that would be observed
! at that location.  The value is returned in 'radial_vel'.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.
!
! The along-beam component of the 3-d air velocity is computed from the
! u, v, and w fields plus the beam_direction.  The along-beam component
! of power-weighted precipitation fall velocity is added to the result.
!
! Reference: Lin et al., 1983 (J. Climate Appl.Meteor., 1065-1092)

real(r8), dimension(ens_size) :: u, v, w, precip_fall_speed
integer,  dimension(ens_size) :: u_istatus, v_istatus, w_istatus, p_istatus
logical  :: return_now
real(r8) :: debug_location(3)
logical  :: debug = .false.   ! set to .true. to enable debug printout

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)
istatus(:)= 0

call interpolate(state_handle, ens_size, location, QTY_U_WIND_COMPONENT, u, u_istatus)
call track_status(ens_size, u_istatus, radial_vel, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_V_WIND_COMPONENT, v, v_istatus)
call track_status(ens_size, v_istatus, radial_vel, istatus, return_now)
if (return_now) return


call interpolate(state_handle, ens_size, location, QTY_VERTICAL_VELOCITY, w, w_istatus)
call track_status(ens_size, w_istatus, radial_vel, istatus, return_now)
if (return_now) return


call get_expected_fall_velocity(state_handle, ens_size, location, precip_fall_speed, p_istatus)
call track_status(ens_size, p_istatus, radial_vel, istatus, return_now)
if (return_now) return

where (istatus == 0)

   radial_vel = radial_vel_data(velkey)%beam_direction(1) * u +    &
                radial_vel_data(velkey)%beam_direction(2) * v +    &
                radial_vel_data(velkey)%beam_direction(3) * (w - precip_fall_speed)

end where

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'radial velocity key: ', velkey
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated u: ', u
   print *, 'interpolated v: ', v
   print *, 'interpolated w: ', w
   print *, 'interp or derived fall speed: ', precip_fall_speed
   print *, 'final radial_vel: ', radial_vel
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_radial_vel
 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Expected fall velocity section
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine get_expected_fall_velocity(state_handle, &
                ens_size, location, precip_fall_speed, istatus)

! This is the main forward operator routine for the expected
! fall velocity, and it also used as part of computing expected
! radial velocity.

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: precip_fall_speed(ens_size)
integer,             intent(out) :: istatus(ens_size)

! Local vars
logical, save :: first_time = .true.
real(r8), dimension(ens_size) :: qr, qg, qs, rho, temp, refl
integer,  dimension(ens_size) :: qr_istatus, qg_istatus, qs_istatus, rho_istatus, refl_istatus, temp_istatus
integer               :: imem
logical               :: return_now
! If the model can return the precipitation fall speed directly, let it do
! so. Otherwise, attempt to compute if Kessler or Lin type microphysics,
! or see if the dbztowt routine is desired. Note that the computation for
! Lin and Kessler is expected to be reasonably accurate for 10cm radar.

! Should we check microphysics_type var or just go ahead and try to get a value?

istatus(:) = 0
precip_fall_speed(:) = 0.0_r8

call interpolate(state_handle, ens_size, location, QTY_POWER_WEIGHTED_FALL_SPEED, &
                 precip_fall_speed, istatus)

! If able to get value, QTY_POWER_WEIGHT_FALL_SPEED is the
! the state so you can return here.
if (any(istatus == 0) ) return

! If the user explicitly wanted to interpolate in the field, try to complain
! if it could not.  Note that the interp could fail for other reasons.
if (microphysics_type == 3 .or. microphysics_type == 4) then
   ! Could return a specific istatus code here to indicate this condition.
   if (first_time) then
      call error_handler(E_MSG,'get_expected_fall_velocity', &
                         'interpolate failed. Fall Speed may NOT be in state vector', '', '', '')
      first_time = .false.
   endif
   return
endif

! QTY_POWER_WEIGHT_FALL_SPEED is not in the state, try to interpolate
! in a different way
istatus(:) = 0

! If not in the state vector, try to calculate it here based on the
! setting of the microphysics_type namelist.

! if Kessler or Lin we can compute the fall velocity
if (microphysics_type == 1 .or. microphysics_type == 2) then
   call interpolate(state_handle, ens_size, location, QTY_RAINWATER_MIXING_RATIO, qr, qr_istatus)
   call track_status(ens_size, qr_istatus, precip_fall_speed, istatus, return_now)
   if (return_now) return

   if (microphysics_type == 2) then
      call interpolate(state_handle, ens_size, location, QTY_GRAUPEL_MIXING_RATIO, qg, qg_istatus)
      call track_status(ens_size, qg_istatus, precip_fall_speed, istatus, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_SNOW_MIXING_RATIO, qs, qs_istatus)
      call track_status(ens_size, qs_istatus, precip_fall_speed, istatus, return_now)
      if (return_now) return

   endif

   ! Done with Lin et al specific calls
   call interpolate(state_handle, ens_size, location, QTY_DENSITY, rho, rho_istatus)
   call track_status(ens_size, rho_istatus, precip_fall_speed, istatus, return_now)
   if (return_now) return

   call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, temp, temp_istatus)
   call track_status(ens_size, temp_istatus, precip_fall_speed, istatus, return_now)
   if (return_now) return



   do imem = 1, ens_size
      if (istatus(imem) == 0) then
         call get_LK_precip_fall_speed(qr(imem), qg(imem), qs(imem), rho(imem), temp(imem), precip_fall_speed(imem))
      endif
   enddo

   ! Done with Lin et al or Kessler -

else if (microphysics_type == 5 .and. allow_dbztowt_conv) then
   ! Provided reflectivity field - will estimate fall velocity using empirical relations
   call get_expected_radar_ref(state_handle, ens_size, location, refl, refl_istatus)
   call track_status(ens_size, refl_istatus, precip_fall_speed, istatus, return_now)
   if (return_now) return

   call interpolate(state_handle, ens_size, location, QTY_DENSITY, rho, rho_istatus)
   call track_status(ens_size, rho_istatus, precip_fall_speed, istatus, return_now)
   if (return_now) return

   do imem = 1, ens_size ! vectorize this routine?
      precip_fall_speed(imem) = dbztowt(refl(imem), rho(imem), missing_r8)
   enddo

else if (microphysics_type < 0) then
   ! User requested setting fall velocity to zero - use with caution
   precip_fall_speed = 0.0_r8
   istatus(:) = 0
else
   ! Couldn't manage to compute a fall velocity
   precip_fall_speed = missing_r8
   istatus = 2
endif


end subroutine get_expected_fall_velocity

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Radar reflectivity
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine read_radar_ref(obsvalue, refkey)

! Main read subroutine for radar reflectivity observations.
! Reflectivity obs have no auxiliary data to read, but this routine
! may optionally alter the observation value at run time.

real(r8),         intent(inout)        :: obsvalue
integer,          intent(out)          :: refkey


if ( .not. module_initialized ) call initialize_module

! Unused, but set to a known value.
refkey = 0

! Alter observation value, if namelist settings request it.
if ((apply_ref_limit_to_obs) .and. &
    (obsvalue < reflectivity_limit_obs) .and. (obsvalue /= missing_r8)) then
   obsvalue = lowest_reflectivity_obs
endif

end subroutine read_radar_ref

!----------------------------------------------------------------------

subroutine get_expected_radar_ref(state_handle, ens_size, location, ref, istatus)

! The main forward operator routine for radar reflectivity observations.

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(in)    :: location
real(r8),            intent(out)   :: ref(ens_size)
integer,             intent(out)   :: istatus(ens_size)

! Given a location and the state vector from one of the ensemble members,
! compute the model-predicted radar reflectivity that would be observed
! at that location.  The returned value is in dBZ.

! If apply_ref_limit_to_fwd_op, reflectivity values which are below
! reflectivity_limit_fwd_op will be set to lowest_reflectivity_fwd_op.

! "interpolate()" ultimately calls model_mod routine "model_interpolate()"
! to get model values of qr, qg, qs, rho, and temp at the ob location. Then
! the routine "get_LK_reflectivity()" is called to compute the radar reflectivity
! factor value, Z, that corresponds to the hydrometeor and thermodynamic values.

real(r8), dimension(ens_size) :: qr, qg, qs, rho, temp
integer,  dimension(ens_size) :: qr_istatus, qg_istatus, qs_istatus
integer,  dimension(ens_size) :: rho_istatus, temp_istatus

logical, save    :: first_time = .true.
integer          :: imem
logical          :: return_now

if ( .not. module_initialized ) call initialize_module

istatus(:) = 0
! Start with known values before calling interpolate routines.
qr(:)   = 0.0_r8
qg(:)   = 0.0_r8
qs(:)   = 0.0_r8
rho(:)  = 0.0_r8
temp(:) = 0.0_r8

! If the model can return radar reflectivity data directly, give it a chance
! to do so.  Otherwise, compute the various fields individually and then do
! the computation here.  Note that the computation here is accurate only
! for the simple single-moment microphysics schemes (e.g., Kessler or Lin).

! Try to draw from state vector first
call interpolate(state_handle, ens_size, location, QTY_RADAR_REFLECTIVITY, ref, istatus)

! If able to get value, QTY_RADAR_REFLECTIVITY is the
! the state so you can return here.
! Make sure to do apply_ref_limit_to_fwd_op check before "return"ing
if (any(istatus == 0) ) then
   if (apply_ref_limit_to_fwd_op) then
      where ((ref < reflectivity_limit_fwd_op) .and. (istatus == 0))
         ref = lowest_reflectivity_fwd_op
      end where
   endif
   return
endif

! If the user explicitly wanted to interpolate in the field, try to complain
! if it could not.  Note that the interp could fail for other reasons.
if (microphysics_type == 3 .or. microphysics_type == 5) then
   ! Could return a specific istatus code here to indicate this condition.
   if (first_time) then
      call error_handler(E_MSG,'get_expected_radar_ref', &
                            'interpolate failed. Reflectivity may NOT be in state vector', '', '', '')
      first_time = .false.
   endif
   return
endif

istatus(:) = 0

if (microphysics_type == 1 .or. microphysics_type == 2) then

   call interpolate(state_handle, ens_size, location, QTY_RAINWATER_MIXING_RATIO, qr, qr_istatus)
   call track_status(ens_size, qr_istatus, ref, istatus, return_now)
   if (return_now) return

   if (microphysics_type == 2) then
      ! Also need some ice vars
      call interpolate(state_handle, ens_size, location, QTY_GRAUPEL_MIXING_RATIO, qg, qg_istatus)
      call track_status(ens_size, qg_istatus, ref, istatus, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_SNOW_MIXING_RATIO, qs, qs_istatus)
      call track_status(ens_size, qs_istatus, ref, istatus, return_now)
      if (return_now) return

   endif

   call interpolate(state_handle, ens_size, location, QTY_DENSITY, rho, rho_istatus)
   call track_status(ens_size, rho_istatus, ref, istatus, return_now)
   if (return_now) return


   call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, temp, temp_istatus)
   call track_status(ens_size, temp_istatus, ref, istatus, return_now)
   if (return_now) return

   do imem = 1, ens_size
      if (istatus(imem) == 0 ) then
         call get_LK_reflectivity(qr(imem), qg(imem), qs(imem), rho(imem), temp(imem), ref(imem))
         ! Always convert to dbz.  Make sure the value, before taking the logarithm,
         ! is always slightly positive.
         ! tiny() is a fortran intrinsic function that is > 0 by a very small amount.
         ref(imem) = 10.0_r8 * log10(max(tiny(ref(imem)), ref(imem)))
      endif
   enddo

else
   ! not in state vector and not Lin et al or Kessler so can't do reflectivity
   ref(:) = missing_r8
   istatus(:) = 1
endif

if (apply_ref_limit_to_fwd_op) then
   where ((ref < reflectivity_limit_fwd_op) .and. (istatus == 0))
      ref = lowest_reflectivity_fwd_op
   end where
endif


end subroutine get_expected_radar_ref

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Helper routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine get_LK_reflectivity(qr, qg, qs, rho, temp, ref)
 
! Computes the equivalent radar reflectivity factor in mm^6 m^-3 for
! simple single-moment microphysics schemes such as the Kessler and 
! Lin et al. schemes.  (Note:  the equivalent reflectivity factor is usually
! called simply "reflectivity" throughout this code.)  For more details
! about the computations, see the initialize_constants subroutine and 
! the references below.

real(r8), intent(in)  :: qr, qg, qs ! rain,graupel,snow mixing ratios (kg kg^-1)
real(r8), intent(in)  :: rho        ! air density (kg m^-3)
real(r8), intent(in)  :: temp       ! air temperature (K)
real(r8), intent(out) :: ref        ! equivalent reflectivity factor (mm^6 m^-3)

! References: Ferrier, 1994 (JAS, 249-280)
!             Gilmore et al., 2004 (MWR, 2610-2627)
!             Lin et al., 1983 (JCAM, 1065-1092)
!             Smith 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)
!             Caya, "Radar observations in DART" (DART subversion repository)

real(r8) :: precip

if ( .not. module_initialized ) call initialize_module

ref = 0.0_r8

! RAIN
if ( qr >= 1.0e-6_r8 ) then
   precip = rho * qr
   ref = ref + param_refl_r * (precip**1.75_r8)
endif

! HAIL / GRAUPEL
! The exponent 1.6625 is 1.75*0.95.  The 0.95 factor is included as an
! approximation for Mie scattering (Smith et al. 1975).  This approximation 
! is appropriate for a 10-cm wavelength radar.

if ( qg >= 1.0e-6_r8 ) then
    precip = rho * qg
    if (.not. allow_wet_graupel .or. temp < 273.15_r8) then
       ref = ref + param_refl_dry_g * (precip**1.75_r8)
    else
       ref = ref + param_refl_wet_g * (precip**1.6625_r8)
    endif
endif

! SNOW
if ( qs >= 1.0e-6_r8 ) then
   precip = rho * qs
   if ( temp < 273.15_r8 ) then
      ref = ref + param_refl_dry_s * (precip**1.75_r8)
   else
      ref = ref + param_refl_wet_s * (precip**1.75_r8)
   endif
endif

end subroutine get_LK_reflectivity

!----------------------------------------------------------------------

subroutine get_LK_precip_fall_speed(qr, qg, qs, rho, temp, precip_fall_speed)

! Computes power-weighted precipitation fall speed in m s^-1 for simple single-
! moment microphysics schemes such as the Kessler and Lin et al. schemes.
! For more details, see the initialize_constants subroutine and the
! references below.

real(r8), intent(in)  :: qr, qg, qs        ! rain,graupel,snow mixing ratios (kg kg^-1)
real(r8), intent(in)  :: rho               ! air density (kg m^-3)
real(r8), intent(in)  :: temp              ! air temperature (K)
real(r8), intent(out) :: precip_fall_speed ! power-weighted precip. fall speed (m s^-1)

! References: Ferrier, 1994 (JAS, 249-280)
!             Gilmore et al., 2004 (MWR, 2610-2627)
!             Lin et al., 1983 (JCAM, 1065-1092)
!             Smith 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)
!             Caya, "Radar observations in DART" (DART subversion repository)

! alpha is the adjustment for air density in the emperical
! dropsize-fall speed equation
real(r8) :: precip_r, precip_s, precip_g, alpha, ref

if ( .not. module_initialized ) call initialize_module

precip_r = rho * qr
precip_s = rho * qs
precip_g = rho * qg
alpha    = sqrt(param_rhos0/rho)

precip_fall_speed = 0.0_r8

! RAIN
if (qr >= 1.0e-6_r8) then
   precip_fall_speed = precip_fall_speed + alpha * param_fs_r * (precip_r**param_powr)
endif

! HAIL/GRAUPEL
if (qg >= 1.0e-6_r8) then
   if (.not. allow_wet_graupel .or. temp < 273.15_r8) then
      precip_fall_speed = precip_fall_speed + sqrt(param_e/rho) * &
                          param_fs_dry_g * (precip_g**param_powg_dry)
   else
      precip_fall_speed = precip_fall_speed + sqrt(param_e/rho) * &
                          param_fs_wet_g * ((precip_g/(PI*rho_graupel))**param_powg_wet)
   endif
endif

! SNOW
if (qs >= 1.0e-6_r8) then
   if ( temp < 273.15_r8 ) then
      precip_fall_speed = precip_fall_speed + alpha * param_fs_dry_s * &
                          (precip_s**param_pows)
   else
      precip_fall_speed = precip_fall_speed + alpha * param_fs_wet_s * &
                          (precip_s**param_pows)
   endif
endif

if (precip_fall_speed > 0.0_r8) then
   call get_LK_reflectivity(qr, qg, qs, rho, temp, ref)
   if (ref > 0.0_r8) then
      precip_fall_speed = precip_fall_speed/ref
   else
      precip_fall_speed = 0.0_r8
   endif
endif

end subroutine get_LK_precip_fall_speed

!----------------------------------------------------------------------

function dbztowt(rf, rho, spval)

! Convert reflectivity (in DBZ, not Z) to terminal fall speed.
! (Code from the pyncommas system - author D. Dowell)

real(r8), intent(in)  :: rf           ! reflectivity (dBZ)
real(r8), intent(in)  :: rho          ! density (km/m**3)
real(r8), intent(in)  :: spval        ! bad/missing data flag
real(r8)              :: dbztowt

! Local vars
real(r8) :: refl         ! reflectivity (Z)

if ( (rf == spval) .or. (rho == spval) ) then
   dbztowt = spval
else
   ! Convert back to Z for this calculation.
   refl = 10.0**(0.1*rf)
   ! Original code used opposite sign - be careful if updating.
   dbztowt = 2.6 * refl**0.107 * (1.2/rho)**0.4
endif

end function dbztowt

!----------------------------------------------------------------------

subroutine initialize_constants()

! Initialize module global constants. 

! IMPORTANT: Uses namelist values, so this routine cannot be called until
! after the namelist has been read.

! NOTE: gravity is currently hardcoded here.  We should find a way to let the
! model have input if it uses a slightly different value for G, or if it is
! working on a different planet with an entirely different set of constants.
! Question:  how much impact on the results does changing G have?
! 
! Most of the constants below are used for the computation of reflectivity and
! power-weighted precipitation fall speed from the model state.  Only simple
! single-moment microphysics schemes, such as the Kessler and Lin et al.
! schemes, are currently supported.  For these schemes, rain, snow, and
! graupel/hail are assumed to have inverse exponential size distributions:
! 
!    n(D)=n0*exp(-lambda*D),
! 
! where D is the particle diameter, n is the number of particles per unit
! volume and per particle size interval, n0 is the intercept parameter, and
! lambda is the slope parameter.  Lambda is a function of the model-predicted
! hydrometeor mixing ratio and air density.
! 
! For Rayleigh scattering from a spherical raindrop, the reflectivity is
! proportional to the 6th moment of the raindrop's diameter.  The basic form of
! the complicated expressions below comes from integrating this relationship
! between reflectivity and diameter over the entire distribution of hydrometeor
! sizes.  Complications arise because not all hydrometeors are spherical, and
! returned power depends significantly on whether hydrometeors have a wet or
! dry surface.  In addition, for large hailstones, Mie scattering rather than
! Rayleigh scattering must be considered.
!
! References: Ferrier, 1994 (JAS, 249-280)
!             Gilmore et al., 2004 (MWR, 2610-2627)
!             Lin et al., 1983 (JCAM, 1065-1092)
!             Smith 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)
!             Caya, "Radar observations in DART" (DART subversion repository)


! In general, the value used for gravity here should match the value used
! for gravity in the model.
param_gravity  = 9.81_r8        ! Acceleration of gravity (m s^-2)
param_a        = 8.42e20_r8
param_b        = 0.8_r8
param_c        = 4.84e18_r8
param_d        = 0.25_r8
param_CD       = 0.6_r8
param_rhos0    = 1.0_r8
param_e        = 4.0_r8*param_gravity*rho_graupel/(3.0_r8*param_CD)

param_gam7b    = 3376.92_r8
param_gam7d    = 1155.38_r8
param_gam7f    = 1871.25_r8
param_powr     = (7.0_r8 + param_b)/4.0_r8
param_pows     = (7.0_r8 + param_d)/4.0_r8
param_powg_dry = (7.0_r8 + 0.5_r8) /4.0_r8
param_powg_wet = 1.7875_r8


! In the expressions for power-weighted precipitation fall speed for each
! hydrometeor category, the following parameters are computed from the
! constants that do not vary in time and space.  These equations are rather
! complicated, so refer to "Radar observations in DART" by A. Caya for details.
! (A pdf of this paper is in the DART subversion repository.)
! The 0.95 exponenent in the equation for param_fs_wet_g is an approximation
! for the effects of Mie scattering (Smith et al. 1975).  This approximation 
! is appropriate for a 10-cm wavelength radar.

param_fs_r     = n0_rain * param_a * param_gam7b / &
                 (PI * rho_rain * n0_rain)**param_powr

param_fs_wet_s = n0_snow * param_c * param_gam7d / &
                 (PI * rho_snow * n0_snow)**param_pows

param_fs_dry_s = dielectric_factor * ((rho_snow/rho_rain)**2.0_r8) * &
                 param_fs_wet_s

param_fs_wet_g = ((7.2e20_r8)**0.95_r8) * param_gam7f / &
                 (720.0_r8 * (n0_graupel**0.8375_r8))

param_fs_dry_g = 1.0e18_r8 * dielectric_factor *                               &
                 ((rho_graupel/rho_rain)**2.0_r8) * n0_graupel * param_gam7f / &
                 (PI * rho_graupel * n0_graupel)**param_powg_dry

! In the expressions for reflectivity for each hydrometeor category, the
! following parameters are computed from the constants that do not vary in time
! and space.  Computing these parameters now means that the equations in the
! get_LK_reflectivity subroutine will have the following simple form:
! 
!   ref = param_refl_r * (rho*qr)**1.75.
!
! The 0.95 exponenent in the equation for param_refl_wet_g is an approximation
! for the effects of Mie scattering (Smith et al. 1975).  This approximation is
! appropriate for a 10-cm wavelength radar.

param_refl_r     = 7.2e20_r8 / (((PI*rho_rain)**1.75_r8)*(n0_rain**0.75_r8))

param_refl_wet_s = 7.2e20_r8 / (((PI*rho_snow)**1.75_r8)*(n0_snow**0.75_r8))

param_refl_dry_s = dielectric_factor * ((rho_snow/rho_rain)**2.0_r8) * & 
                 param_refl_wet_s

param_refl_wet_g = (7.2e20_r8/(((PI*rho_graupel)**1.75_r8) * &
                 (n0_graupel**0.75_r8)))**0.95_r8

param_refl_dry_g = dielectric_factor * ((rho_graupel/rho_rain)**2.0_r8) * &
                 7.2e20_r8 / (((PI*rho_graupel)**1.75_r8)*(n0_graupel**0.75_r8))



end subroutine initialize_constants

!----------------------------------------------------------------------

subroutine print_constants()

! Log the constants set in the code.
! Prints to both the log file and standard output.

! The values in this list which are also in the namelist will have their
! values written by the write(nml=) code, but this routine includes all
! the fixed constants so they are written in one place, both to standard 
! output and the log file. Using the correct values is critical to doing 
! the appropriate computation, so some duplication is probably a good thing.

write(msgstring, *) 'Constants used in the obs_def_radar module:'
call error_handler(E_MSG,'', msgstring, '', '', '')

call pr_con(dielectric_factor , "dielectric_factor" )
call pr_con(n0_rain           , "n0_rain"           )
call pr_con(n0_graupel        , "n0_graupel"        )
call pr_con(n0_snow           , "n0_snow"           )
call pr_con(rho_rain          , "rho_rain"          )
call pr_con(rho_graupel       , "rho_graupel"       )
call pr_con(rho_snow          , "rho_snow"          )
call pr_con(param_gravity     , "param_gravity"     )
call pr_con(param_a           , "param_a"           )
call pr_con(param_b           , "param_b"           )
call pr_con(param_c           , "param_c"           )
call pr_con(param_d           , "param_d"           )
call pr_con(param_CD          , "param_CD"          )
call pr_con(param_rhos0       , "param_rhos0"       )
call pr_con(param_e           , "param_e"           )
call pr_con(param_gam7b       , "param_gam7b"       )
call pr_con(param_gam7d       , "param_gam7d"       )
call pr_con(param_gam7f       , "param_gam7f"       )
call pr_con(param_powr        , "param_powr"        )
call pr_con(param_pows        , "param_pows"        )
call pr_con(param_powg_dry    , "param_powg_dry"    )
call pr_con(param_powg_wet    , "param_powg_wet"    )
call pr_con(param_fs_r        , "param_fs_r"        )
call pr_con(param_fs_wet_s    , "param_fs_wet_s"    )
call pr_con(param_fs_dry_s    , "param_fs_dry_s"    )
call pr_con(param_fs_wet_g    , "param_fs_wet_g"    )
call pr_con(param_fs_dry_g    , "param_fs_dry_g"    )
call pr_con(param_refl_r      , "param_refl_r"        )
call pr_con(param_refl_wet_s  , "param_refl_wet_s"    )
call pr_con(param_refl_dry_s  , "param_refl_dry_s"    )
call pr_con(param_refl_wet_g  , "param_refl_wet_g"    )
call pr_con(param_refl_dry_g  , "param_refl_dry_g"    )

end subroutine print_constants

!----------------------------------------------------------------------

subroutine pr_con(c_val, c_str)

! Utility routine to print a string and value

real(r8),         intent(in) :: c_val
character(len=*), intent(in) :: c_str

write(msgstring, "(A30,A,ES28.8)") c_str, " = ", c_val
call error_handler(E_MSG,'', msgstring, '', '', '')

end subroutine pr_con

!----------------------------------------------------------------------

subroutine velkey_out_of_range(velkey)

! Range check velkey and trigger a fatal error if larger than allocated array.

integer, intent(in) :: velkey

! fine -- no problem.
if (velkey <= max_radial_vel_obs) return

! Bad news.  Tell the user.
write(msgstring, *) 'velkey (',velkey,') exceeds max_radial_vel_obs (', &
                     max_radial_vel_obs,')'
call error_handler(E_MSG,'set_radial_vel', msgstring, '', '', '')
call error_handler(E_ERR,'set_radial_vel', &
                   'Increase max_radial_vel_obs in namelist', &
                   source, revision, revdate)

end subroutine velkey_out_of_range

!----------------------------------------------------------------------

subroutine check_namelist_limits(apply_ref_limit_to_obs, &
   reflectivity_limit_obs, lowest_reflectivity_obs, apply_ref_limit_to_fwd_op,& 
   reflectivity_limit_fwd_op, lowest_reflectivity_fwd_op)

! Consistency warning; print a message if the thresholds and lower values
! are going to be used and are different.

logical,  intent(in) :: apply_ref_limit_to_obs
real(r8), intent(in) :: reflectivity_limit_obs
real(r8), intent(in) :: lowest_reflectivity_obs
logical,  intent(in) :: apply_ref_limit_to_fwd_op
real(r8), intent(in) :: reflectivity_limit_fwd_op
real(r8), intent(in) :: lowest_reflectivity_fwd_op

! The point here is to gently remind the user if they are setting different
! limits on the actual observation values and the forward operator, since
! that may be what they intend, but it isn't something they should be doing
! by mistake.  But we don't want to be annoying, so this code only prints
! informative messages to the log file and does not stop with an error.

! If neither limit is being enforced, return silently.
if (.not. apply_ref_limit_to_obs .and. .not. apply_ref_limit_to_fwd_op) return

! If both are on, and the limits and set-to values are the same, fine also.
if (apply_ref_limit_to_obs .and. apply_ref_limit_to_fwd_op) then
   if ((reflectivity_limit_obs  == reflectivity_limit_fwd_op) .and. &
       (lowest_reflectivity_obs == lowest_reflectivity_fwd_op)) return
endif

! Either only one of the limits is on, and/or the limits or set-to 
! values do not match.  Print something to the log file to note that they
! are not the same.
if (apply_ref_limit_to_obs) then
   write(msgstring, *) 'reflectivity obs values less than ', &
      reflectivity_limit_obs, 'will be set to', lowest_reflectivity_obs
else
   write(msgstring, *) 'reflectivity obs values will be processed unchanged'
endif
call error_handler(E_MSG,'check_namelist_limits', msgstring, '', '', '')

if (apply_ref_limit_to_fwd_op) then
   write(msgstring, *) 'reflectivity forward operator values less than ', &
      reflectivity_limit_fwd_op, 'will be set to', lowest_reflectivity_fwd_op
else
   write(msgstring, *) 'reflectivity forward operator values will be processed unchanged'
endif
call error_handler(E_MSG,'check_namelist_limits', msgstring, '', '', '')

end subroutine check_namelist_limits

!----------------------------------------------------------------------

end module obs_def_radar_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

