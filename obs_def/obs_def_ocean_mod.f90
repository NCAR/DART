! Data Assimilation Research Testbed -- DART source
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! BEGIN DART PREPROCESS KIND LIST
!SALINITY,                      KIND_SALINITY,              COMMON_CODE
!TEMPERATURE,                   KIND_TEMPERATURE,           COMMON_CODE
!U_CURRENT_COMPONENT,           KIND_U_CURRENT_COMPONENT,   COMMON_CODE
!V_CURRENT_COMPONENT,           KIND_V_CURRENT_COMPONENT,   COMMON_CODE
!SEA_SURFACE_HEIGHT,            KIND_SEA_SURFACE_HEIGHT,    COMMON_CODE
!ARGO_U_CURRENT_COMPONENT,      KIND_U_CURRENT_COMPONENT,   COMMON_CODE
!ARGO_V_CURRENT_COMPONENT,      KIND_V_CURRENT_COMPONENT,   COMMON_CODE
!ARGO_SALINITY,                 KIND_SALINITY,              COMMON_CODE
!ARGO_TEMPERATURE,              KIND_TEMPERATURE,           COMMON_CODE
!ADCP_U_CURRENT_COMPONENT,      KIND_U_CURRENT_COMPONENT,   COMMON_CODE
!ADCP_V_CURRENT_COMPONENT,      KIND_V_CURRENT_COMPONENT,   COMMON_CODE
!ADCP_SALINITY,                 KIND_SALINITY,              COMMON_CODE
!ADCP_TEMPERATURE,              KIND_TEMPERATURE,           COMMON_CODE
!FLOAT_SALINITY,                KIND_SALINITY,              COMMON_CODE
!FLOAT_TEMPERATURE,             KIND_TEMPERATURE,           COMMON_CODE
!DRIFTER_U_CURRENT_COMPONENT,   KIND_U_CURRENT_COMPONENT,   COMMON_CODE
!DRIFTER_V_CURRENT_COMPONENT,   KIND_V_CURRENT_COMPONENT,   COMMON_CODE
!DRIFTER_SALINITY,              KIND_SALINITY,              COMMON_CODE
!DRIFTER_TEMPERATURE,           KIND_TEMPERATURE,           COMMON_CODE
!GLIDER_U_CURRENT_COMPONENT,    KIND_U_CURRENT_COMPONENT,   COMMON_CODE
!GLIDER_V_CURRENT_COMPONENT,    KIND_V_CURRENT_COMPONENT,   COMMON_CODE
!GLIDER_SALINITY,               KIND_SALINITY,              COMMON_CODE
!GLIDER_TEMPERATURE,            KIND_TEMPERATURE,           COMMON_CODE
!MOORING_U_CURRENT_COMPONENT,   KIND_U_CURRENT_COMPONENT,   COMMON_CODE
!MOORING_V_CURRENT_COMPONENT,   KIND_V_CURRENT_COMPONENT,   COMMON_CODE
!MOORING_SALINITY,              KIND_SALINITY,              COMMON_CODE
!MOORING_TEMPERATURE,           KIND_TEMPERATURE,           COMMON_CODE
!MOORING_PRESSURE,              KIND_PRESSURE,              COMMON_CODE
!DOPPLER_U_CURRENT_COMPONENT,   KIND_U_CURRENT_COMPONENT,   COMMON_CODE
!DOPPLER_V_CURRENT_COMPONENT,   KIND_V_CURRENT_COMPONENT,   COMMON_CODE
!DOPPLER_W_CURRENT_COMPONENT,   KIND_W_CURRENT_COMPONENT,   COMMON_CODE
!SATELLITE_MICROWAVE_SST,       KIND_TEMPERATURE,           COMMON_CODE
!SATELLITE_INFRARED_SST,        KIND_TEMPERATURE,           COMMON_CODE
!SATELLITE_SSH,                 KIND_SEA_SURFACE_HEIGHT,    COMMON_CODE
!SATELLITE_SSS,                 KIND_SALINITY,              COMMON_CODE
!CODAR_RADIAL_VELOCITY,         KIND_VELOCITY
! END DART PREPROCESS KIND LIST

! From Ibrahim - 19 May 2009
! "Concerning the radar observations: we can use a format very similar to what we 
! have been using for all the other types of observations:
!
!LAT LON DEPTH ANGLE OBS_value Z_LOC ERR_obs QC_flag TYPE_obs STARTDATE1 STARTDATE2
!
! where here we only add "ANGLE": the angle of the radial observation.
!
! Now what we need is to compute the radial velocity R_model from the state vector 
! (which already contains the horizontal velocities U and V) using 
!
! 'OBS_value' = U*cos(angle) + V*sin(angle)
!
! The angles are in degrees (0 from the east) and the observations in m/s.
! so to compute 'OBS_value' we need to read the observation file to get the ANGLE, 
! and that's basically the only difference with the assimilation of the other 
! types of observations. We also need to introduce the new obs type 
! (CODAR_RADIAL_VELOCITY)."

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_ocean_mod, only : write_radial_vel, read_radial_vel,           &
!                            interactive_radial_vel, get_expected_radial_vel
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(CODAR_RADIAL_VELOCITY)
!     call get_expected_radial_vel(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(CODAR_RADIAL_VELOCITY)
!      call read_radial_vel(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(CODAR_RADIAL_VELOCITY)
!      call write_radial_vel(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(CODAR_RADIAL_VELOCITY)
!      call interactive_radial_vel(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_ocean_mod

use        types_mod, only : r8, missing_r8, PI, deg2rad
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location, get_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_U_CURRENT_COMPONENT, &
                             KIND_V_CURRENT_COMPONENT

implicit none
private

public :: read_radial_vel, write_radial_vel, interactive_radial_vel,       &
          get_expected_radial_vel, get_obs_def_radial_vel, set_radial_vel

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical :: module_initialized = .false.

character(len=129) :: msgstring ! For error message content

! Derived type for radial velocity.  Contains auxiliary information stored
! with each obs of this type; used to compute the forward operator.
! See more extensive comments in the interactive_radial_vel() routine for
! expected units, etc.  The instrument ID is currently unused, but may be
! useful for post processing or diagnostics.

type radial_vel_type
   private
   integer   :: instrument_id       ! ID number for data source
   real(r8)  :: beam_angle          ! angle of beam
end type radial_vel_type

! Cumulative index into radial velocity metadata array
integer :: velkeycount = 0 
                             
! namelist items
integer  :: max_radial_vel_obs         = 1000000

namelist /obs_def_ocean_mod_nml/ max_radial_vel_obs

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
call find_namelist_in_file("input.nml", "obs_def_ocean_mod_nml", iunit)
read(iunit, nml = obs_def_ocean_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_ocean_mod_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_ocean_mod_nml)
if (do_nml_term()) write(     *     , nml=obs_def_ocean_mod_nml)

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
integer             :: instrument_id
real(r8)            :: beam_angle
integer             :: oldkey

! instrument id: Arbitrary number to distinguish different sources of
!  the data.  Not used in any of this code, but carried along in case
!  it is useful in postprocessing or diagnostics.

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
      ! Read the character identifier for verbose formatted output
      read(ifile, FMT="(a5)") header
      if(header /= 'CODAR') then
         call error_handler(E_ERR,'read_radial_vel', &
              "Expected header 'CODAR' in input file", &
              source, revision, revdate)
      endif
endif

! read_location is a DART library routine that expects an optional string
! arg, but the other two read routines are local to this module and we can
! tell them exactly what format to be reading because we already know it.
instrument_id    = read_instrument_id (ifile, is_asciifile)
beam_angle       = read_beam_angle    (ifile, is_asciifile)

! Read in the velkey for this particular observation, however, it will
! be discarded and a new, unique key will be generated in the set routine.
if (is_asciifile) then
   read(ifile, *) oldkey
else
   read(ifile) oldkey
endif

! Generate new unique radial velocity observation key, and set the contents
! of the private defined type.
call set_radial_vel(velkey, instrument_id, beam_angle )

end subroutine read_radial_vel


!----------------------------------------------------------------------


subroutine write_radial_vel(velkey, ifile, fform)

! Write radial velocity auxiliary information to the obs_seq file.

integer,          intent(in)           :: velkey, ifile
character(len=*), intent(in), optional :: fform

logical  :: is_asciifile
integer  :: instrument_id
real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Write the 5 character identifier for verbose formatted output
   write(ifile, "('CODAR')")
endif

! Extract the values for this key and call the appropriate write routines.
instrument_id = radial_vel_data(velkey)%instrument_id
beam_angle    = radial_vel_data(velkey)%beam_angle

! These two write routines are local to this module, and we have already 
! figured out if it is a unformatted/binary file or formatted/ascii, so 
! go ahead and pass that info directly down to the routines.
call write_instrument_id(ifile, instrument_id, is_asciifile)
call write_beam_angle(ifile, beam_angle, is_asciifile)

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


function read_instrument_id(ifile, is_asciiformat)

! Reads instrument id from file that was written by write_instrument_id.
! See read_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
integer             :: read_instrument_id

character(len=5)   :: header
integer            :: instrument_id

if ( .not. module_initialized ) call initialize_module


if (is_asciiformat) then
   read(ifile, "(a5)" ) header

   if(header /= 'InsID') then
      write(msgstring,*)"Expected Instrument ID header 'InsID' in input file, got ", header
      call error_handler(E_ERR, 'read_instrument_id', msgstring, source, revision, revdate)  
   endif
   ! Read the instrument id data value into temporary
   read(ifile, *) instrument_id
else
   ! No header label, just the binary id value. 
   read(ifile)    instrument_id
endif

! set function return value
read_instrument_id = instrument_id

end function read_instrument_id


!----------------------------------------------------------------------


subroutine write_instrument_id(ifile, instrument_id, is_asciiformat)

! Writes instrument_id to obs file.

integer, intent(in) :: ifile
integer, intent(in) :: instrument_id
logical, intent(in) :: is_asciiformat

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   write(ifile, "('InsID')" ) 
   write(ifile, *) instrument_id
else
   write(ifile)    instrument_id
endif

end subroutine write_instrument_id


!----------------------------------------------------------------------


function read_beam_angle(ifile, is_asciiformat)

! Reads beam_angle from file that was written by write_beam_angle.
! See read_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
real(r8)            :: read_beam_angle

character(len=5)   :: header
real(r8)           :: beam_angle

if ( .not. module_initialized ) call initialize_module


if (is_asciiformat) then
   read(ifile, "(a5)" ) header

   if(header /= 'angle') then
      write(msgstring,*)"Expected beam_angle header 'angle' in input file, got ", header
      call error_handler(E_ERR, 'read_beam_angle', msgstring, source, revision, revdate)  
   endif
   ! Now read the beam_angle data value into temporaries
   read(ifile, *) beam_angle
else
   ! No header label, just the binary angle values. 
   read(ifile)    beam_angle
endif

! Check for illegal values 
if (beam_angle < 0.0_r8 .or. beam_angle > 360.0_r8 ) then
   write(msgstring,*) "beam_angle value must be between 0 and 360, got: ", &
                       beam_angle
   call error_handler(E_ERR, 'read_beam_angle', msgstring, &
                      source, revision, revdate)
endif

! set function return value
read_beam_angle = beam_angle

end function read_beam_angle


!----------------------------------------------------------------------


subroutine write_beam_angle(ifile, beam_angle, is_asciiformat)

! Writes beam_angle to obs file.

integer,  intent(in) :: ifile
real(r8), intent(in) :: beam_angle
logical,  intent(in) :: is_asciiformat

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   write(ifile, "('angle')" ) 
   write(ifile, *) beam_angle
else
   write(ifile)    beam_angle
endif

end subroutine write_beam_angle


!----------------------------------------------------------------------


subroutine get_obs_def_radial_vel(velkey, instrument_id, beam_angle)

! Return the auxiliary contents of a given radial velocity observation

integer,   intent(in)  :: velkey
integer,   intent(out) :: instrument_id
real(r8),  intent(out) :: beam_angle

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

instrument_id = radial_vel_data(velkey)%instrument_id
beam_angle    = radial_vel_data(velkey)%beam_angle

end subroutine get_obs_def_radial_vel


!----------------------------------------------------------------------


subroutine set_radial_vel(velkey, instrument_id, beam_angle )

! Common code to increment the current key count, and set the private
! contents of this observation's auxiliary data.

integer,  intent(out) :: velkey
integer,  intent(in)  :: instrument_id
real(r8), intent(in)  :: beam_angle

if ( .not. module_initialized ) call initialize_module

! The total velocity metadata key count from all sequences
velkeycount = velkeycount + 1
velkey = velkeycount             ! set the return value
   
! Simple error check on key number before accessing the array
! This errors out if too key value now too large.
call velkey_out_of_range(velkey)

radial_vel_data(velkey)%instrument_id = instrument_id
radial_vel_data(velkey)%beam_angle    = beam_angle
   
end subroutine set_radial_vel 


!----------------------------------------------------------------------


subroutine interactive_radial_vel(velkey)

! Interactively reads in auxiliary information for a radial velocity obs.

integer, intent(out) :: velkey

! Uses the local subroutines interactive_instrument_id and set_radial_vel.
! See read_radial_vel for more information.

! velkey is internally incremented in the set routine, and only counts 
! the index for this specialized observation kind. 

integer  :: instrument_id
real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

! Get information from user about id, angle.


write(*, *) 
write(*, *) 'Beginning to inquire for information on instrument id.'
write(*, *)

call interactive_instrument_id(instrument_id)


write(*, *) 
write(*, *) 'Beginning to inquire for information on radar beam angle.'
write(*, *)

call interactive_beam_angle(beam_angle)


call set_radial_vel(velkey, instrument_id, beam_angle )

write(*, *) 
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_radial_vel


!----------------------------------------------------------------------


subroutine interactive_instrument_id(instrument_id)

! Prompt for instrument_id.  Not used in this code, but carried along.

integer, intent(out) :: instrument_id

write(*, *) 'Input an integer instrument_id:'
read(*, *) instrument_id

end subroutine interactive_instrument_id


!----------------------------------------------------------------------


subroutine interactive_beam_angle(beam_angle)

! Prompt for beam angle information in azimuth/elevation degrees.

real(r8), intent(out) :: beam_angle

real(r8) :: az

az = -1.0
do while (az < 0.0_r8 .or. az > 360.0_r8) 
   write(*, *) 'Input the beam angle in degrees (0 <= az <= 360):'
   read(*, *) az
end do

! Convert to radians and compute the actual values stored with the observation.
az = az * deg2rad

beam_angle = az

end subroutine interactive_beam_angle


!----------------------------------------------------------------------


subroutine get_expected_radial_vel(state_vector, location, velkey, &
                                   radial_vel, istatus)

! This is the main forward operator routine for radar Doppler velocity.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: velkey
real(r8),            intent(out) :: radial_vel
integer,             intent(out) :: istatus


! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted radial velocity that would be observed
! at that location.  The value is returned in 'radial_vel'. 
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.
!
! The along-beam component of the 3-d ocean surface velocity is computed 
! from the u, v fields plus the beam_angle:
!  return_value = U*cos(angle) + V*sin(angle)
!


real(r8) :: u, v
real(r8) :: debug_location(3) 
logical  :: debug = .false.   ! set to .true. to enable debug printout

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

call interpolate(state_vector, location, KIND_U_CURRENT_COMPONENT, u, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_V_CURRENT_COMPONENT, v, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif

radial_vel = u * cos(radial_vel_data(velkey)%beam_angle) + &
             v * sin(radial_vel_data(velkey)%beam_angle)

! Good return code. 
istatus = 0

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'radial velocity key: ', velkey
   print *, 'instrument id: ', radial_vel_data(velkey)%instrument_id
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated u: ', u
   print *, 'interpolated v: ', v
   print *, 'angle: ', radial_vel_data(velkey)%beam_angle
   print *, 'final radial_vel: ', radial_vel
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_radial_vel

!----------------------------------------------------------------------

function ascii_file_format(fform)

! Common routine for determining input file format. 

character(len=*), intent(in), optional :: fform     
logical                                :: ascii_file_format

! Returns .true. if file is formatted/ascii, .false. if unformatted/binary
! Defaults (if fform not specified) to formatted/ascii.

if ( .not. module_initialized ) call initialize_module

! Default to formatted/ascii. 
if ( .not. present(fform)) then 
   ascii_file_format = .true. 
   return 
endif

SELECT CASE (trim(adjustl(fform)))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED") 
      ascii_file_format = .false.
   CASE DEFAULT
      ascii_file_format = .true.
END SELECT


end function ascii_file_format

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



end module obs_def_ocean_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------
