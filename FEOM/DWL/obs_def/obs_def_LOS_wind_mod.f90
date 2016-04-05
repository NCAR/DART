! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!-----------------------------------------------------------------------------
! DART Line Of Sight (LOS) or azimuthal wind velocity observations.   These 
! differ from wind speed/angle because ONLY the velocity along the LOS is known.
! There is no information about the velocity of the normal component.  
! These observations require 1 additional piece of metadata - the angle of
! the LOS/azimuth.   The forward operator retrieves the U/V winds and constructs
! the wind speed as measured along the LOS angle.
!
! This converter assumes wind speed is in cm/s, azimuth angle is 0 at due
! north and increases clockwise.  The types of observations include
! Raleigh clear/cloudy and Mie clear/cloudy.
!
! Author and Contact information:
!   LOS Wind Science:  Nedjeljka Zagar, nedjeljka.zagar@fmf.uni-lj.si
!                      Matic Savli, matic.savli@fmf.uni-lj.si
!   DART Code:  Nancy Collins, nancy at ucar.edu
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS KIND LIST
! DWL_RAY_CLEAR_LOS_VELOCITY, KIND_VELOCITY
! DWL_MIE_CLEAR_LOS_VELOCITY, KIND_VELOCITY
! DWL_RAY_CLOUD_LOS_VELOCITY, KIND_VELOCITY
! DWL_MIE_CLOUD_LOS_VELOCITY, KIND_VELOCITY
! END DART PREPROCESS KIND LIST
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_los_vel_mod, only : write_los_vel, read_los_vel,           &
!                           interactive_los_vel, get_expected_los_vel,    &
!                           get_los_angle, set_los_angle
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(DWL_RAY_CLEAR_LOS_VELOCITY, DWL_MIE_CLEAR_LOS_VELOCITY, &
!       DWL_RAY_CLOUD_LOS_VELOCITY, DWL_MIE_CLOUD_LOS_VELOCITY)
!     call get_expected_los_vel(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!  case(DWL_RAY_CLEAR_LOS_VELOCITY, DWL_MIE_CLEAR_LOS_VELOCITY, &
!       DWL_RAY_CLOUD_LOS_VELOCITY, DWL_MIE_CLOUD_LOS_VELOCITY)
!      call read_los_vel(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!  case(DWL_RAY_CLEAR_LOS_VELOCITY, DWL_MIE_CLEAR_LOS_VELOCITY, &
!       DWL_RAY_CLOUD_LOS_VELOCITY, DWL_MIE_CLOUD_LOS_VELOCITY)
!      call write_los_vel(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!  case(DWL_RAY_CLEAR_LOS_VELOCITY, DWL_MIE_CLEAR_LOS_VELOCITY, &
!       DWL_RAY_CLOUD_LOS_VELOCITY, DWL_MIE_CLOUD_LOS_VELOCITY)
!      call interactive_los_vel(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_los_vel_mod

use        types_mod, only : r8, missing_r8, PI, deg2rad
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location, get_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_VELOCITY

implicit none
private

public :: write_los_vel, read_los_vel,                  &
          interactive_los_vel, get_expected_los_vel,    &
          get_los_angle, set_los_angle

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical :: module_initialized = .false.

! Derived type for Line Of Sight (LOS) wind.  Additional metadata required
! for each obs of this type is the angle of the observation.  Azimuth is
! specified in degrees, 0-360.  0 is due north and increases clockwise
! (due east is 90, south 180, west 270).
type los_vel_type
   private
   real(r8)            :: azimuth   ! direction of observation
end type los_vel_type

! Cumulative index into los azimuth angle metadata array
integer :: loskeycount = 0 

! For error message content
character(len=129) :: msgstring

!-------------------------------------------------------------
! Set the max possible obs here.  Increase at run time if needed.
! 
 
integer :: max_los_vel_obs = 100000


namelist /obs_def_los_vel_mod_nml/ max_los_vel_obs


! Module global storage for auxiliary obs data, allocated in init routine
type(los_vel_type), allocatable :: los_vel_data(:)

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
call find_namelist_in_file("input.nml", "obs_def_los_vel_mod_nml", iunit)
read(iunit, nml = obs_def_los_vel_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_los_vel_mod_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_los_vel_mod_nml)
if (do_nml_term()) write(     *     , nml=obs_def_los_vel_mod_nml)

! Allocate space for the auxiliary information associated with each obs
! This code must be placed after reading the namelist, so the user can
! increase or decrease the number of obs supported and use more or less
! memory at run time.
allocate(los_vel_data(max_los_vel_obs), stat = rc)
if (rc /= 0) then            
   write(msgstring, *) 'initial allocation failed for los vel obs data,', &
                       'itemcount = ', max_los_vel_obs
   call error_handler(E_ERR,'initialize_module', msgstring, &
                      source, revision, revdate)
endif                        

end subroutine initialize_module

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LOS velocity section
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine read_los_vel(loskey, ifile, fform)

! Main read subroutine for the los velocity observation auxiliary data.

integer,          intent(out)          :: loskey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=4)    :: header
logical             :: is_asciifile
real(r8)            :: los_angle
integer             :: oldkey

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
      ! Read the character identifier for verbose formatted output
      read(ifile, FMT="(a4)") header
      if(header /= 'LOS') then
         call error_handler(E_ERR,'read_los_vel', &
              "Expected location header 'LOS' in input file", &
              source, revision, revdate)
      endif
endif

los_angle = read_los_angle(ifile, is_asciifile)

! Read in the loskey for this particular observation, however, it will
! be discarded and a new, unique key will be generated in the set routine.
if (is_asciifile) then
   read(ifile, *) oldkey
else
   read(ifile) oldkey
endif

! Generate new unique los velocity observation key, and set the contents
! of the private defined type.
call set_los_angle(loskey, los_angle)

end subroutine read_los_vel

!----------------------------------------------------------------------

subroutine write_los_vel(loskey, ifile, fform)

! Write los velocity auxiliary information to the obs_seq file.

integer,          intent(in)           :: loskey, ifile
character(len=*), intent(in), optional :: fform

logical             :: is_asciifile
real(r8)            :: los_angle

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call loskey_out_of_range(loskey)

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Write the 3 character identifier for verbose formatted output
   write(ifile, "('LOS')")
endif

! Extract the values for this key and call the appropriate write routines.
los_angle  = los_vel_data(loskey)%azimuth

! we have already figured out if it is a unformatted/binary file or 
! formatted/ascii, so go ahead and pass that info directly down to the routine.
call write_los_angle(ifile, los_angle,  is_asciifile) 

! Write out the loskey used for this run, however this will be discarded
! when this observation is read in and a new key will be generated.
! (It may be useful for correlating error messages or identifying particular
! observations so we are leaving it as part of the aux data.)
if (is_asciifile) then
   write(ifile, *) loskey
else
   write(ifile) loskey
endif

end subroutine write_los_vel

!----------------------------------------------------------------------

function read_los_angle(ifile, is_asciiformat)

! Reads LOS angle from file that was written by write_los_angle.
! See read_los_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
real(r8)            :: read_los_angle

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   read(ifile, *) read_los_angle
else
   read(ifile)    read_los_angle
endif

end function read_los_angle

!----------------------------------------------------------------------

subroutine write_los_angle(ifile, los_angle, is_asciiformat)

! Writes LOS angle to obs file.

integer,  intent(in) :: ifile
real(r8), intent(in) :: los_angle
logical,  intent(in) :: is_asciiformat

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   write(ifile, *) los_angle
else
   write(ifile)    los_angle
endif

end subroutine write_los_angle

!----------------------------------------------------------------------

subroutine get_los_angle(loskey, angle)

! Return the auxiliary contents of a given LOS velocity observation
! the input is the private key number for the LOS obs.

integer,  intent(in)  :: loskey
real(r8), intent(out) :: angle

! Simple error check on key number before accessing the array
call loskey_out_of_range(loskey)

angle = los_vel_data(loskey)%azimuth

end subroutine get_los_angle

!----------------------------------------------------------------------

subroutine set_los_angle(loskey, angle)

! Common code to increment the current key count, and set the private
! contents of this observation's auxiliary data.

integer,  intent(out) :: loskey
real(r8), intent(in)  :: angle

if ( .not. module_initialized ) call initialize_module

! The total velocity metadata key count from all sequences
loskeycount = loskeycount + 1    
loskey = loskeycount             ! set the return value

! Simple error check on key number before accessing the array
! This errors out if too key value now too large.
call loskey_out_of_range(loskey)

los_vel_data(loskey)%azimuth = angle

end subroutine set_los_angle

!----------------------------------------------------------------------

subroutine interactive_los_vel(loskey)

! Interactively reads in auxiliary information for a los velocity obs.

integer, intent(out) :: loskey

! Uses the local subroutine set_los_angle.

! loskey is internally incremented in the set routine, and only counts 
! the index for this specialized observation kind. 

real(r8)             :: angle

if ( .not. module_initialized ) call initialize_module

write(*, *)
write(*, *) 'Beginning to inquire for information on LOS wind velocity.'
write(*, *)
write(*, *)

call interactive_los_angle(angle)

call set_los_angle(loskey, angle)

write(*, *)
write(*, *) 'End of specialized section for LOS wind velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_los_vel

!----------------------------------------------------------------------

subroutine interactive_los_angle(angle)

! Prompt for Line of Sight angle.  True north is 0, then increases
! clockwise.  In degrees.

real(r8), intent(out) :: angle

angle = -1.0
do while (angle < 0.0 .or. angle > 360.0) 
   write(*, *) 'Input the Line Of Sight angle in degrees, 0-360. (N=0,E=90,S=180,W=270):'
   read(*, *) angle
end do

end subroutine interactive_los_angle

!----------------------------------------------------------------------

subroutine get_expected_los_vel(state_vector, location, loskey, &
                                los_vel, istatus)

! This is the main forward operator routine for LineOfSight winds

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: loskey
real(r8),            intent(out) :: los_vel
integer,             intent(out) :: istatus


! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted LineOfSight velocity that would be observed
! at that location along the given angle.  The value is returned in 'los_vel'. 
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.
!

real(r8) :: u, v, w, qr, qg, qs, rho, temp, precip_fall_speed
real(r8) :: debug_location(3)
logical  :: debug = .false.   ! set to .true. to enable debug printout

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call loskey_out_of_range(loskey)

call interpolate(state_vector, location, KIND_U_WIND_COMPONENT, u, istatus)
if (istatus /= 0) then
   los_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_V_WIND_COMPONENT, v, istatus)
if (istatus /= 0) then
   los_vel = missing_r8
   return
endif

! FIXME: compute here.  you have u, v (wind components), los_vel_data(loskey)%azimuth
! in degrees, N=0, increasing clockwise.  replace los_vel line below with correct
! computation.  right now it's ignoring the angle and assuming it is always 90.
!  should be something like u * sin(angle * PI / 180.0_r8) +
!                           v * cos(angle * PI / 180.0_r8) +
! and remember angle is 0=N, increasing clockwise.  it is NOT 0=x axis like the
! mathematical definition of polar coordinate angle.

los_vel = u

! Good return code.
istatus = 0

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'los velocity key: ', loskey
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated u: ', u
   print *, 'interpolated v: ', v
   print *, 'final los_vel: ', los_vel
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_los_vel


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Helper routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine loskey_out_of_range(loskey)

! Range check loskey and trigger a fatal error if larger than allocated array.

integer, intent(in) :: loskey

! fine -- no problem.
if (loskey <= max_los_vel_obs) return

! Bad news.  Tell the user.
write(msgstring, *) 'loskey (',loskey,') exceeds max_los_vel_obs (', &
                     max_los_vel_obs,')'
call error_handler(E_MSG,'set_los_angle', msgstring, '', '', '')
call error_handler(E_ERR,'set_los_angle', &
                   'Increase max_los_vel_obs in namelist', &
                   source, revision, revdate)

end subroutine loskey_out_of_range

!----------------------------------------------------------------------

end module obs_def_los_vel_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
