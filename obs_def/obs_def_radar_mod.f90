! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_radar_mod

! BEGIN DART PREPROCESS KIND LIST
! DOPPLER_RADIAL_VELOCITY, KIND_VELOCITY
! RADAR_REFLECTIVITY, KIND_RADAR_REFLECTIVITY
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_radar_mod, only : write_rad_vel, read_rad_vel, &
!                                 interactive_rad_vel, get_expected_rad_vel, &
!                                 write_rad_ref, read_rad_ref, &
!                                 interactive_rad_ref, get_expected_rad_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(DOPPLER_RADIAL_VELOCITY)
!            call get_expected_rad_vel(state, location, obs_def%key, obs_val, istatus)
!         case(RADAR_REFLECTIVITY)
!            call get_expected_rad_ref(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call read_rad_vel(obs_def%key, ifile, fileformat)
!      case(RADAR_REFLECTIVITY)
!         call read_rad_ref(obs_val, obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call write_rad_vel(obs_def%key, ifile, fileformat)
!      case(RADAR_REFLECTIVITY)
!         call write_rad_ref(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call interactive_rad_vel(obs_def%key)
!      case(RADAR_REFLECTIVITY)
!         call interactive_rad_ref(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! <next five lines automatically updated by CVS, do not edit>
! $Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_radar_mod.f90,v $
! $Revision$
! $Date$
! $Author$
! $Name:  $

use        types_mod, only : r8, missing_r8, ps0, PI, gravity, DEG2RAD
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file, &
                             logfileunit, do_output
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_TEMPERATURE, KIND_VERTICAL_VELOCITY, &
                             KIND_RAINWATER_MIXING_RATIO, KIND_DENSITY, &
                             KIND_GRAUPEL_MIXING_RATIO, KIND_SNOW_MIXING_RATIO

implicit none
private

public :: write_rad_vel, read_rad_vel, set_rad_vel, interactive_rad_vel, &
          write_rad_ref, read_rad_ref, interactive_rad_ref, &
          get_expected_rad_vel, get_expected_rad_ref, &
          get_obs_def_rad_ref, get_obs_def_rad_vel

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_radar_mod.f90,v $", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.

! Storage for the special information required for velocity observations of this type
integer, parameter  :: max_rad_vel_obs = 1000000
type(location_type) :: rad_loc(max_rad_vel_obs)
real(r8)            :: direction(3,max_rad_vel_obs)
real(r8)            :: nyquistvel(max_rad_vel_obs)

! Storage for the special information required for reflectivity observations of this type
integer, parameter  :: max_rad_ref_obs = 1000000
!define here        :: globally_defined_reflectivity_metadata_variable(max_rad_ref_obs)

integer             :: velkeycount = 0 ! cumulative index into rad. velocity metadata
integer             :: refkeycount = 0 ! cumulative index into rad. reflectivity metadata

real(r8), parameter :: dief = 0.224_r8

real(r8), parameter :: n0r = 8.0e6_r8, n0g = 4.0e4_r8, n0s = 3.0e6_r8
real(r8), parameter :: rho_r = 1000.0_r8, rho_g = 917.0_r8, rho_s = 100.0_r8

!-------------------------------------------------------------
! Namelist with default values
! 
! convert_to_dbz:           Convert final reflectivity to dB? 
!                           If .TRUE., then make sure dbz_threshold is specified correctly.
!                           (Note that for now, this is reset to .TRUE. in initialize_module
!                           as it is not clear how to assimilate Z.)
!
! dbz_threshold:            Only applied if convert_to_dbz is .TRUE.
!                           Reflectivity below dbz_threshold is set to dbz_threshold, which is chosen
!                           such to convert to the smallest possible dBZ.
!                           The purpose is mainly computational: We want to avoid the occurence of log(0).
!                           Some useful values: Z = 3.163  ->   5.0 dBZ
!                                               Z = 2.512  ->   4.0 dBZ
!                                               Z = 1.259  ->   1.0 dBZ
!                                               Z = 1.0233 ->   0.1 dBZ
!                                               Z = 0.1    -> -10.0 dBZ
!                                               Z = 0.01   -> -20.0 dBZ
!                                               Z = 0.001  -> -30.0 dBZ
!
! apply_ref_limit_to_obs:    If .TRUE., replace reflectivities below "reflectivity_limit" 
!                            with "lowest_reflectivity". If FALSE, ignore "reflectivity_limit"
!                            and "lowest_reflectivity".  
!
! reflectivity_limit_obs:    Reflectivity below this value is set to "lowest_reflectivity". 
!
! lowest_reflectivity_obs:   If reflectivity < reflectivity_limit, reflectivity is set to this value.
!                            The default value of 'missing' is useful for perfect model experiments.
!                            Suggested options: lowest_reflectivity_obs = missing_r8 or 
!                            lowest_reflectivity_obs = any real value.
!
! apply_ref_limit_to_state:  Similar to "apply_ref_limit_to_obs" but applied to the state
!                            (forward operator).
!
! reflectivity_limit_state:  Similar to "reflectivity_limit_obs" for state.
!
! lowest_reflectivity_state: Similar to "lowest_reflectivity_obs" for state.


logical  :: convert_to_dbz            = .true.
real(r8) :: dbz_threshold             = 0.001_r8
logical  :: apply_ref_limit_to_obs    = .false.
real(r8) :: reflectivity_limit_obs    = 0.0_r8
real(r8) :: lowest_reflectivity_obs   = missing_r8
logical  :: apply_ref_limit_to_state  = .false.
real(r8) :: reflectivity_limit_state  = 0.0_r8
real(r8) :: lowest_reflectivity_state = missing_r8

namelist /obs_def_radar_mod_nml/ convert_to_dbz, dbz_threshold, &
                                 apply_ref_limit_to_obs, reflectivity_limit_obs, &
                                 lowest_reflectivity_obs, &
                                 apply_ref_limit_to_state, reflectivity_limit_state, &
                                 lowest_reflectivity_state


contains

!----------------------------------------------------------------------

subroutine initialize_module

implicit none

integer :: iunit, io

call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_radar_mod_nml", iunit)
read(iunit, nml = obs_def_radar_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_radar_mod_nml")

! Record the namelist values used for the run ... 
if (do_output()) call error_handler(E_MSG,'obs_def_radar_mod:initialize_module',&
                         'obs_def_radar_mod_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=obs_def_radar_mod_nml)
if (do_output()) write(     *     , nml=obs_def_radar_mod_nml)

! This is temporary as it is not clear how to assimilate Z
convert_to_dbz = .true.

module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine write_rad_vel(velkey, ifile, fform)

integer,          intent(in)           :: velkey, ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      call write_location(ifile,  rad_loc(velkey), fileformat)
      call write_orientation(ifile, direction(:,velkey), fileformat)
      call write_nyquistvel(ifile, nyquistvel(velkey), fileformat)
      ! Write out the obs_def velkey for this observation
      write(ifile) velkey
   CASE DEFAULT
      ! Write the 5 character identifier for verbose formatted output
      write(ifile, 11)
11    format('platform')
      call write_location(ifile,  rad_loc(velkey), fileformat)
      call write_orientation(ifile, direction(:,velkey), fileformat)
      call write_nyquistvel(ifile, nyquistvel(velkey), fileformat)
      ! Write out the obs_def velkey for this observation
      write(ifile, *) velkey
END SELECT

end subroutine write_rad_vel

!----------------------------------------------------------------------

subroutine get_obs_def_rad_vel(velkey, radarloc, beamdir, velnyquist)

integer, intent(in)              :: velkey
type(location_type), intent(out) :: radarloc
real(r8), intent(out)            :: beamdir(3)
real(r8), intent(out)            :: velnyquist

radarloc   = rad_loc(velkey)
beamdir    = direction(:,velkey)
velnyquist = nyquistvel(velkey)

end subroutine get_obs_def_rad_vel

!----------------------------------------------------------------------

subroutine read_rad_vel(velkey, ifile, fform)

! Main read subroutine for the obs kind radial velocity metadata.

! location: Refers to the lat/lon/height/vertical coordinate option
!           for the radar. Its type is defined through location_type in
!           location_mod. Uses the read_location function in
!           location_mod.
!
! orientation: This is a 3-element array specific to the radar module:
!              orientation(1) = sin(azimuth)*cos(elevation)
!              orientation(2) = cos(azimuth)*cos(elevation)
!              orientation(3) = sin(elevation)

integer,          intent(out)          :: velkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=8)    :: header
character(len=32)   :: fileformat
real(r8)            :: orientation(3)
type(location_type) :: location
real(r8)            :: nyquistv

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      location = read_location(ifile, fileformat)
      orientation = read_orientation(ifile, fileformat)
      nyquistv = read_nyquistvel(ifile, fileformat)
      ! Read in the velkey for this particular observation
      read(ifile) velkey
   CASE DEFAULT
      ! Read the character identifier for verbose formatted output
      read(ifile, FMT='(a8)') header
      if(header /= 'platform') then
         call error_handler(E_ERR,'obs_def_radar_mod:read_rad_vel', &
              'Expected location header "platform" in input file', &
              source, revision, revdate)
      endif
      location = read_location(ifile, fileformat)
      orientation = read_orientation(ifile, fileformat)
      nyquistv = read_nyquistvel(ifile, fileformat)
      ! Read in the velkey for this particular observation
      read(ifile, *) velkey
END SELECT

velkeycount = velkeycount + 1    ! the total velocity metadata key count from all sequences
velkey = velkeycount             ! copied to the output variable

call set_rad_vel(velkey, location, orientation, nyquistv)

end subroutine read_rad_vel

!----------------------------------------------------------------------

subroutine set_rad_vel(velkey, rad_location, rad_orientation, rad_nyquistv)

integer,             intent(in)  :: velkey
real(r8),            intent(in)  :: rad_orientation(3)
type(location_type), intent(in)  :: rad_location
real(r8),            intent(in)  :: rad_nyquistv

if ( .not. module_initialized ) call initialize_module

!Make sure enough space is allocated
if(velkey > max_rad_vel_obs) then
   write(*, *) 'velkey (',velkey,') exceeds max_rad_vel_obs (',max_rad_vel_obs,')'
   call error_handler(E_ERR,'read_rad_vel:set_rad_vel', &
              'Increase max_rad_vel_obs.', source, revision, revdate)
endif

rad_loc(velkey)      = rad_location
direction(:, velkey) = rad_orientation
nyquistvel(velkey)   = rad_nyquistv

end subroutine set_rad_vel

!----------------------------------------------------------------------

subroutine write_rad_ref(refkey, ifile, fform)

integer,          intent(in)           :: refkey, ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! At this point, this is empty as there is no metadata for radar reflectivity
! in the current obs sequence format. In the future, if more metadata pieces
! are added, this part should be modified.

end subroutine write_rad_ref

!----------------------------------------------------------------------

subroutine get_obs_def_rad_ref(refkey, is_this_dbz)

integer, intent(in)  :: refkey
logical, intent(out) :: is_this_dbz
 
is_this_dbz = convert_to_dbz
 
end subroutine get_obs_def_rad_ref

!----------------------------------------------------------------------

subroutine read_rad_ref(obsvalue, refkey, ifile, fform)

! Main read subroutine for the obs kind radar reflectivity metadata.

! reftype: Denotes whether reflectivity obs units are in Z (reftype = .false.)
!                                              or are in dBZ (reftype = .true.)

real(r8),         intent(inout)        :: obsvalue
integer,          intent(out)          :: refkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=8)    :: header
character(len=32)   :: fileformat
logical             :: reftype

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! At this point, no metadata is available for radar reflectivity. For this
! reason, internal counter "refkey" has also no consequence. In the
! future, if more metadata pieces are added, this part should be modified.
refkeycount = refkeycount + 1    ! the total reflectivity metadata key count from all sequences
refkey = refkeycount             ! copied to the output variable

! Convert to Z and appy reflectivity limit if requested through the namelist
! (Note that conversion to Z is suppressed for now by globally resetting
! convert_to_dbz to .TRUE. after reading namelist)
if (obsvalue /= missing_r8) then
   if (.not. convert_to_dbz) then
      obsvalue = 10.0_r8 ** (obsvalue/10.0_r8)
   endif

   if ( (apply_ref_limit_to_obs) .and. (obsvalue < reflectivity_limit_obs) ) then
      obsvalue = lowest_reflectivity_obs
   endif
endif

! Again, since there is no metadata for reflectivity to track, we don't
! need the following call at the moment; thus the
! following commented line and the subroutine "set_rad_ref". In the future,
! if more metadata pieces are added, this part should be modified and
! "set_rad_ref" should be uncommented.

!call set_rad_ref(refkey,some_reflectivity_metadata)

end subroutine read_rad_ref

!----------------------------------------------------------------------
!
!subroutine set_rad_ref(refkey, some_reflectivity_metadata)
!
!integer, intent(in)     :: refkey
!define here, intent(in) :: some_reflectivity_metadata
!
!if ( .not. module_initialized ) call initialize_module
!
!!Make sure enough space is allocated
!if(refkey > max_rad_ref_obs) then
!   write(*, *) 'refkey (',refkey,') exceeds max_rad_ref_obs (',max_rad_ref_obs,')'
!   call error_handler(E_ERR,'obs_def_radar_mod:read_rad_ref:set_rad_ref', &
!              'Increase max_rad_ref_obs.', source, revision, revdate)
!endif
!
!globally_defined_reflectivity_metadata_variable(refkey) = some_reflectivity_metadata
!
!end subroutine set_rad_ref

!----------------------------------------------------------------------

subroutine interactive_rad_vel(velkey)

! Interactively reads in location and orientation components of the obs
! kind radial velocity. See read_rad_vel for more information on these
! two components.

! Uses interactive_location of location_mod and the local subroutines
! interactive_orientation and set_rad_vel.

! velkey is internally incremented and only counts the index for this specialized
! observation kind. It is written in the obs file after the radar location
! and beam direction.

implicit none

integer, intent(out) :: velkey
real(r8)             :: orientation(3)
type(location_type)  :: location
real(r8)             :: nyquistv

if ( .not. module_initialized ) call initialize_module

!Make sure enough space is allocated
if(velkey > max_rad_vel_obs) then
   write(*, *) 'velkey (',velkey,') exceeds max_rad_vel_obs (',max_rad_vel_obs,')'
   call error_handler(E_ERR,'obs_def_radar_mod:interactive_rad_vel', &
              'Increase max_rad_vel_obs.', source, revision, revdate)
endif

!Now interactively obtain radar location and beam direction information
!Note: Obs location is additionally inquired by the standard DART module
!      "interactive_obs_def". No check is performed here whether radar
!      location, obs location, and beam direction match. It is the user's
!      responsibility to make sure that they do.


write(*, *)
write(*, *) 'Beginning to inquire information on radar location.'
write(*, *) 'Make sure that you select 3 (height) for vertical coordinate'
write(*, *) 'option and enter height in units gpm.'
write(*, *)

call interactive_location(location)

write(*, *)
write(*, *) 'Beginning to inquire information on radar beam direction.'
write(*, *)

call interactive_orientation(orientation)

write(*, *)
write(*, *) 'Beginning to inquire information on radar Nyquist velocity.'
write(*, *)

call interactive_nyquistvel(nyquistv)

call set_rad_vel(velkey+1, location, orientation, nyquistv)

! Finally increment velkey to keep track of specialized obs
velkey = velkey + 1

write(*, *)
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_rad_vel

!----------------------------------------------------------------------

subroutine interactive_rad_ref(refkey)

! Interactively (sort of) reads in reftype component of the obs
! kind radar reflectivity. See read_rad_ref for more information on these
! two components.

! Uses the local subroutine interactive_reftype.

! refkey is internally incremented and only counts the index for this specialized
! observation kind. It is written in the obs file after the reflectivity type.

implicit none

integer, intent(out) :: refkey

if ( .not. module_initialized ) call initialize_module

!Make sure enough space is allocated
if(refkey > max_rad_ref_obs) then
   write(*, *) 'refkey (',refkey,') exceeds max_rad_ref_obs (',max_rad_ref_obs,')'
   call error_handler(E_ERR,'obs_def_radar_mod:interactive_rad_ref', &
              'Increase max_rad_ref_obs.', source, revision, revdate)
endif

!Now interactively obtain reflectivity type information
!Note: Reflectivity type will be actually read from the namelist variable
!      convert_to_dbz. Thus, for reftype, the code is not really interactive.
!      Nevertheless, this was chosen to be this way to be consistent with the
!      rest of the code and make future modifications more straightforward.

write(*, *)
write(*, *) 'Beginning to inquire information on reflectivity type.'
write(*, *) 'This will be read from the namelist variable convert_to_dbz.'
write(*, *)

! If reflectivity metadata is added in the future, insert here a call to a
! subroutine to read that metadata interactively (and write the subroutine
! that goes with it!)

! Finally increment refkey to keep track of specialized obs
refkey = refkey + 1

write(*, *)
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_rad_ref

!----------------------------------------------------------------------

subroutine interactive_orientation(orientation)

implicit none

real(r8), intent(out) :: orientation(3)

write(*, *) 'Input first component: sin(azimuth)*cos(elevation)'
read(*, *)  orientation(1)

do while((orientation(1) > 1.0_r8) .or. (orientation(1) < -1.0_r8))
   write(*, *) 'Input value is not allowed, please enter a value -1.0 to 1.0'
   read(*, *) orientation(1)
end do

write(*, *) 'Input second component: cos(azimuth)*cos(elevation)'
read(*, *)  orientation(2)

do while((orientation(2) > 1.0_r8) .or. (orientation(2) < -1.0_r8))
   write(*, *) 'Input value is not allowed, please enter a value -1.0 to 1.0'
   read(*, *) orientation(2)
end do

write(*, *) 'Input third component: sin(elevation)'
read(*, *)  orientation(3)

do while((orientation(3) > 1.0_r8) .or. (orientation(3) < 0.0_r8))
   write(*, *) 'Input value is not allowed, please enter a value 0.0 to 1.0'
   read(*, *) orientation(3)
end do

end subroutine interactive_orientation

!----------------------------------------------------------------------

subroutine interactive_nyquistvel(nyquistv)

implicit none

real(r8), intent(out) :: nyquistv

write(*, *) 'Input Nyquist velocity for this obs point'
read(*, *)  nyquistv

end subroutine interactive_nyquistvel

!----------------------------------------------------------------------

subroutine get_expected_rad_vel(state_vector, location, velkey, vr, istatus)

! This is the main forward operator routine for radar Doppler velocity.
! Given an ob DART location, computes model-predicted radial velocity at
! the same location.

! In addition to u,v,w, from which and the direction of the radar beam is
! calculated the along-beam component of the 3-d wind vector, we also need
! the reflectivity value to compute the terminal velocity of the hydrometeor,
! so that w can be adjusted for it. See also get_expected_rad_ref.

! Reference: Lin et al., 1983 (J. Climate Appl.Meteor., 1065-1092)
! Note that the reflectivity-weighted mean terminal velocities are used here.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: velkey
real(r8),            intent(out) :: vr
integer,             intent(out) :: istatus

real(r8), parameter :: a        = 8.42e20_r8
real(r8), parameter :: b        = 0.8_r8
real(r8), parameter :: c        = 4.84e18_r8
real(r8), parameter :: d        = 0.25_r8
real(r8), parameter :: CD       = 0.6_r8
real(r8), parameter :: rhos0    = 1.0_r8
real(r8), parameter :: e        = 4.0_r8*gravity*rho_g/(3.0_r8*CD)
real(r8), parameter :: f        = 0.5_r8
real(r8), parameter :: gam7b    = 3376.92_r8
real(r8), parameter :: gam7d    = 1155.38_r8
real(r8), parameter :: gam7f    = 1871.25_r8
real(r8), parameter :: powr     = (7.0_r8 + b)/4.0_r8
real(r8), parameter :: pows     = (7.0_r8 + d)/4.0_r8
real(r8), parameter :: powg_dry = (7.0_r8 + f)/4.0_r8
real(r8), parameter :: powg_wet = 1.7875_r8

real(r8) :: ar
real(r8) :: as_wet
real(r8) :: as_dry
real(r8) :: ag_dry

real(r8) :: u, v, w, qr, qg, qs, alpha, wt, rho, temp, precip, ref
real(r8) :: precip_r, precip_s, precip_g

ar       = n0r*a*gam7b / (PI*rho_r*n0r)**powr
as_wet   = n0s*c*gam7d / (PI*rho_s*n0s)**pows
as_dry   = dief*((rho_s/rho_r)**2.0_r8)*as_wet
ag_dry   = 1.0e18_r8*dief*((rho_g/rho_r)**2.0_r8)* &
           n0g*gam7f/(PI*rho_g*n0g)**powg_dry

if ( .not. module_initialized ) call initialize_module

call interpolate(state_vector, location, KIND_U_WIND_COMPONENT, u, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_V_WIND_COMPONENT, v, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_VERTICAL_VELOCITY, w, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO, qr, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, qg, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, qs, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_DENSITY, rho, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_TEMPERATURE, temp, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif

precip_r = rho * qr
precip_s = rho * qs
precip_g = rho * qg
alpha    = sqrt(rhos0/rho)

wt = 0.0_r8

! Computing reflectivity-weighted terminal velocity wt
! RAIN
if (qr >= 1.0e-6_r8) then
   wt = alpha * ar * (precip_r**powr)
endif
! HAIL/GRAUPEL
if (qg >= 1.0e-6_r8) then
   wt = wt + sqrt(e/rho)*ag_dry*(precip_g**powg_dry)
   ! Note that we are assuming dry graupel/hail surface
   ! Thus, no temperature check for hail/graupel (see snow for comparison)
   ! ag_wet = ((7.2e20_r8)**0.95_r8)*gam7f/(720.0_r8*(n0g**0.8375_r8))
   ! wt     = wt + sqrt(e/rho)*ag_wet*(((rho * qg)/(PI*rho_g))**powg_wet)
endif
! SNOW
if (qs >= 1.0e-6_r8) then
   if ( temp < 273.15_r8 ) then
      wt = wt + alpha*as_dry*(precip_s**pows)
   else
      wt = wt + alpha*as_wet*(precip_s**pows)
   endif
endif

if (wt > 0.0_r8) then
   call get_reflectivity(qr, qg, qs, rho, temp, ref)
   wt = wt/ref
endif

vr = direction(1, velkey)*u + direction(2, velkey)*v + direction(3, velkey)*(w-wt)

end subroutine get_expected_rad_vel

!----------------------------------------------------------------------

subroutine get_expected_rad_ref(state_vector, location, ref, istatus)
!
! This is the main forward operator routine for radar reflectivity.
! Given an ob DART location, computes model-predicted radar reflectivity at
! the same location.

! Radar reflectivity = 10 * log_10( radar reflectivity factor) in dBZ.

! If apply_ref_limit_to_state, reflectivities below
! reflectivity_limit_state are set to lowest_reflectivity_state.

! Subroutine "interpolate" ultimately calls model_mod routine "model_interpolate"
! to get model values of qr, qg, qs, rho, and temp at the ob location. Then
! the routine "get_reflectivity" is called to compute the radar reflectivity
! factor value, Z, that corresponds to the hydrometeor and thermodynamic values.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: ref
integer,             intent(out) :: istatus

real(r8) :: qr, qg, qs, rho, temp

if ( .not. module_initialized ) call initialize_module

call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO , qr , istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, qg, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, qs, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_DENSITY, rho, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_TEMPERATURE, temp, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call get_reflectivity(qr, qg, qs, rho, temp, ref)

if (convert_to_dbz) then
   ref = 10.0_r8 * log10(max(dbz_threshold, ref))
endif

if ( (apply_ref_limit_to_state) .and. (ref < reflectivity_limit_state) .and. &
     (ref /= missing_r8) ) then
   ref = lowest_reflectivity_state
endif

if (ref == missing_r8) then
   istatus = 1
endif

end subroutine get_expected_rad_ref

!----------------------------------------------------------------------

subroutine get_reflectivity(qr, qg, qs, rho, temp, ref)
!
! Computes "radar reflectivity factor" (Z) in mm^6 m^-3
!
! References: Ferrier, 1994 (JAS, 249-280)
!             Smith et al., 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)

! According to Smith (1984), there are two choices for the dielectric
! factor (dief), depending on how the snow particle sizes are specified.
! If melted raindrop diameters are used, then the factor is 0.224.  If
! equivalent ice sphere diameters are used, then the factor is 0.189.

real(r8), intent(in)  :: qr, qg, qs, rho, temp
real(r8), intent(out) :: ref

real(r8) :: precip
real(r8) :: ar 
real(r8) :: ag_dry
real(r8) :: as_wet
real(r8) :: as_dry


ar     = 7.2e20_r8  /  (((PI*rho_r)**1.75_r8)*(n0r**0.75_r8))
ag_dry = dief*((rho_g/rho_r)**2.0_r8)*7.2e20_r8/ &
           (((PI*rho_g)**1.75_r8)*(n0g**0.75_r8))
as_wet = 7.2e20_r8/(((PI*rho_s)**1.75_r8)*(n0s**0.75_r8))
as_dry = dief*((rho_s/rho_r)**2.0_r8)*as_wet 

if ( .not. module_initialized ) call initialize_module

ref = 0.0_r8

! RAIN
if ( qr >= 1.0e-6_r8 ) then
   precip = rho * qr
   ref = ref + ar * (precip**1.75_r8)
endif

! HAIL / GRAUPEL
if ( qg >= 1.0e-6_r8 ) then
    precip = rho * qg
    ref = ref + ag_dry * (precip**1.75_r8)
    ! Note that we assume dry surface for hail/graupel
    ! Thus, no temperature check for hail/graupel (see snow for comparison)
    ! ag_wet = (7.2e20_r8/(((PI*rho_g)**1.75_r8)*(n0g**0.75_r8)))**0.95_r8
    ! ref    = ref + ag_wet * (precip**1.6625_r8)
endif

! SNOW
if ( qs >= 1.0e-6_r8 ) then
   precip = rho * qs
   if ( temp < 273.15_r8 ) then
      ref = ref + as_dry * (precip**1.75_r8)
   else
      ref = ref + as_wet * (precip**1.75_r8)
   endif
endif

end subroutine get_reflectivity

!----------------------------------------------------------------------

subroutine write_orientation(ifile, orientation, fform)

! Writes orientation to obs file.
! This variable is only associated with the obs kind radial velocity.
! See read_rad_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: orientation(3)
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) orientation(1), orientation(2), orientation(3)
   CASE DEFAULT
      write(ifile, '(''dir3d'')' ) 
      write(ifile, *) orientation(1), orientation(2), orientation(3)
END SELECT

end subroutine write_orientation

!----------------------------------------------------------------------

subroutine write_nyquistvel(ifile, nyquistv, fform)

! Writes Nyquist velocity to obs file.
! This variable is only associated with the obs kind radial velocity.
! See read_rad_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: nyquistv
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) nyquistv
   CASE DEFAULT
      write(ifile, *) nyquistv
END SELECT

end subroutine write_nyquistvel

!----------------------------------------------------------------------

subroutine write_reftype(ifile, fform)

! Writes reflectivity type to obs file.
! This variable is only associated with the obs kind radar reflectivity.
! See read_rad_ref for additional discussion.

! Note that when writing to the obs file, we already know that this value
! should be equal to the namelist variable convert_to_dbz.

integer,                    intent(in) :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) convert_to_dbz
   CASE DEFAULT
      write(ifile, *) convert_to_dbz
END SELECT

end subroutine write_reftype

!----------------------------------------------------------------------

function read_orientation(ifile, fform)

! Reads orientation from file that was written by write_orientation.
! This variable is only associated with the obs kind radial velocity.
! See read_rad_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8)                               :: read_orientation(3)
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_orientation(1), read_orientation(2), read_orientation(3)
   CASE DEFAULT
      read(ifile, '(a5)' ) header

      if(header /= 'dir3d') then
         write(errstring,*)'Expected orientation header "dir3d" in input file, got ', header
         call error_handler(E_ERR, 'read_orientation', errstring, source, revision, revdate)
      endif
! Now read the orientation data value
      read(ifile, *) read_orientation(1), read_orientation(2), read_orientation(3)
END SELECT

end function read_orientation

!----------------------------------------------------------------------

function read_nyquistvel(ifile, fform)

! Reads Nyquist velocity from file that was written by write_nyquistvel.
! This variable is only associated with the obs kind radial velocity.
! See read_rad_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8)                               :: read_nyquistvel
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_nyquistvel
   CASE DEFAULT
      read(ifile, *) read_nyquistvel
END SELECT

end function read_nyquistvel

!----------------------------------------------------------------------

function read_reftype(ifile, fform)

! Reads reflectivity type from file that was written by write_reftype.
! This variable is only associated with the obs kind radar reflectivity.
! See read_rad_ref for additional discussion.

integer,                    intent(in) :: ifile
logical                                :: read_reftype
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_reftype
   CASE DEFAULT
      read(ifile, *) read_reftype
END SELECT

end function read_reftype



end module obs_def_radar_mod
