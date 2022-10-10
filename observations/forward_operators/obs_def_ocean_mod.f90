! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! SALINITY,                      QTY_SALINITY,              COMMON_CODE
! TEMPERATURE,                   QTY_TEMPERATURE,           COMMON_CODE
! U_CURRENT_COMPONENT,           QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! V_CURRENT_COMPONENT,           QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! SEA_SURFACE_HEIGHT,            QTY_SEA_SURFACE_HEIGHT,    COMMON_CODE
! SEA_SURFACE_PRESSURE,          QTY_SEA_SURFACE_PRESSURE,  COMMON_CODE
! ARGO_U_CURRENT_COMPONENT,      QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! ARGO_V_CURRENT_COMPONENT,      QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! ARGO_SALINITY,                 QTY_SALINITY,              COMMON_CODE
! ARGO_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
! ADCP_U_CURRENT_COMPONENT,      QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! ADCP_V_CURRENT_COMPONENT,      QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! ADCP_SALINITY,                 QTY_SALINITY,              COMMON_CODE
! ADCP_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
! FLOAT_SALINITY,                QTY_SALINITY,              COMMON_CODE
! FLOAT_TEMPERATURE,             QTY_TEMPERATURE,           COMMON_CODE
! DRIFTER_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! DRIFTER_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! DRIFTER_SALINITY,              QTY_SALINITY,              COMMON_CODE
! DRIFTER_TEMPERATURE,           QTY_TEMPERATURE,           COMMON_CODE
! GLIDER_U_CURRENT_COMPONENT,    QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! GLIDER_V_CURRENT_COMPONENT,    QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! GLIDER_SALINITY,               QTY_SALINITY,              COMMON_CODE
! GLIDER_TEMPERATURE,            QTY_TEMPERATURE,           COMMON_CODE
! MOORING_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! MOORING_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! MOORING_SALINITY,              QTY_SALINITY,              COMMON_CODE
! MOORING_TEMPERATURE,           QTY_TEMPERATURE,           COMMON_CODE
! MOORING_PRESSURE,              QTY_PRESSURE,              COMMON_CODE
! BOTTLE_SALINITY,               QTY_SALINITY,              COMMON_CODE
! BOTTLE_TEMPERATURE,            QTY_TEMPERATURE,           COMMON_CODE
! CTD_SALINITY,                  QTY_SALINITY,              COMMON_CODE
! CTD_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
! TCTD_SALINITY,                 QTY_SALINITY,              COMMON_CODE
! TCTD_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
! STD_SALINITY,                  QTY_SALINITY,              COMMON_CODE
! STD_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
! XCTD_SALINITY,                 QTY_SALINITY,              COMMON_CODE
! XCTD_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
! MBT_SALINITY,                  QTY_SALINITY,              COMMON_CODE
! MBT_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
! XBT_SALINITY,                  QTY_SALINITY,              COMMON_CODE
! XBT_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
! DBT_SALINITY,                  QTY_SALINITY,              COMMON_CODE
! DBT_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
! APB_SALINITY,                  QTY_SALINITY,              COMMON_CODE
! APB_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
! DOPPLER_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! DOPPLER_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! DOPPLER_W_CURRENT_COMPONENT,   QTY_W_CURRENT_COMPONENT,   COMMON_CODE
! SATELLITE_MICROWAVE_SST,       QTY_TEMPERATURE,           COMMON_CODE
! SATELLITE_INFRARED_SST,        QTY_TEMPERATURE,           COMMON_CODE
! SATELLITE_BLENDED_SST,         QTY_TEMPERATURE,           COMMON_CODE
! SATELLITE_SSH,                 QTY_SEA_SURFACE_HEIGHT,    COMMON_CODE
! SATELLITE_SSS,                 QTY_SALINITY,              COMMON_CODE
! J1_SEA_SURFACE_ANOMALY,        QTY_SEA_SURFACE_ANOMALY,   COMMON_CODE
! EN_SEA_SURFACE_ANOMALY,        QTY_SEA_SURFACE_ANOMALY,   COMMON_CODE
! GFO_SEA_SURFACE_ANOMALY,       QTY_SEA_SURFACE_ANOMALY,   COMMON_CODE
! DRY_LAND,                      QTY_DRY_LAND,              COMMON_CODE
! OI_SEA_SURFACE_TEMPERATURE,    QTY_TEMPERATURE,           COMMON_CODE
! HFRADAR_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
! HFRADAR_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
! HFRADAR_RADIAL_VELOCITY,       QTY_VELOCITY
! FERRYBOX_SALINITY,             QTY_SALINITY,              COMMON_CODE
! FERRYBOX_TEMPERATURE,          QTY_TEMPERATURE,           COMMON_CODE
! OCEAN_COLOR,                   QTY_SURFACE_CHLOROPHYLL,   COMMON_CODE
! END DART PREPROCESS TYPE DEFINITIONS

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
! (HFRADAR_RADIAL_VELOCITY)."
 
! changed all occurences of xxx_radial_vel to xxx_hf_radial_vel (for high 
! frequency) since otherwise the external names here collide with those
! external names in the radar module (for weather radars, e.g. nexrad).  nsc.

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_ocean_mod, only : write_hf_radial_vel, read_hf_radial_vel,           &
!                            interactive_hf_radial_vel, get_expected_hf_radial_vel
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(HFRADAR_RADIAL_VELOCITY)
!     call get_expected_hf_radial_vel(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(HFRADAR_RADIAL_VELOCITY)
!      call read_hf_radial_vel(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(HFRADAR_RADIAL_VELOCITY)
!      call write_hf_radial_vel(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(HFRADAR_RADIAL_VELOCITY)
!      call interactive_hf_radial_vel(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_ocean_mod

use        types_mod, only : r8, missing_r8, PI, deg2rad
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location, get_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_U_CURRENT_COMPONENT, &
                             QTY_V_CURRENT_COMPONENT
use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: read_hf_radial_vel, write_hf_radial_vel, interactive_hf_radial_vel,       &
          get_expected_hf_radial_vel, get_obs_def_hf_radial_vel, set_hf_radial_vel

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical :: module_initialized = .false.

character(len=129) :: msgstring ! For error message content
character(len= 32) :: header, str1

! Derived type to contain the metadata for radial velocity observations.
! Contains auxiliary information used to compute the forward operator.
! See more extensive comments in the interactive_hf_radial_vel() routine for
! expected units, etc.  The instrument ID is currently unused, but may be
! useful for post processing or diagnostics.

type radial_vel_type
   private
   integer   :: instrument_id       ! ID number for data source
   real(r8)  :: beam_angle          ! angle of beam (degrees)
end type radial_vel_type

! Module storage for radial velocity metadata, allocated in init routine.
! Cumulative index into radar metadata array
integer :: metadata_index = 0 
type(radial_vel_type), allocatable :: radial_vel_metadata(:)

! namelist items
integer :: max_radial_vel_obs = 1000000
logical :: debug = .false.   ! set to .true. to enable debug printout

namelist /obs_def_ocean_nml/ max_radial_vel_obs, debug

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
call find_namelist_in_file("input.nml", "obs_def_ocean_nml", iunit)
read(iunit, nml = obs_def_ocean_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_ocean_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_ocean_nml)
if (do_nml_term()) write(     *     , nml=obs_def_ocean_nml)

! Allocate space for the auxiliary information associated with each obs
! This code must be placed after reading the namelist, so the user can
! increase or decrease the number of obs supported and use more or less
! memory at run time.
allocate(radial_vel_metadata(max_radial_vel_obs), stat = rc)
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


subroutine read_hf_radial_vel(velkey, ifile, fform)

! Main read subroutine for the radial velocity observation auxiliary data.
! The beam_angle is expected to have units of degrees.

integer,          intent(out)          :: velkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

logical  :: is_asciifile
integer  :: instrument_id
real(r8) :: beam_angle
integer  :: oldkey

! instrument id: Arbitrary number to distinguish different sources of
!  the data.  Not used in any of this code, but carried along in case
!  it is useful in postprocessing or diagnostics.

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Read the character identifier for verbose formatted output
   read(ifile, FMT="(a)") header
   str1 = adjustl(header)

   if(str1(1:7) /= 'HFRADAR') then
      write(msgstring,*)"Expected header 'HFRADAR' in input file, got ", str1(1:7)
      call error_handler(E_ERR,'read_hf_radial_vel', &
           msgstring, source, revision, revdate)
   endif
endif

! read_location is a DART library routine that expects an optional string
! arg, but the other two read routines are local to this module and we can
! tell them exactly what format to be reading because we already know it.
instrument_id = read_instrument_id (ifile, is_asciifile)
beam_angle    = read_beam_angle    (ifile, is_asciifile)

! Read in the velkey for this particular observation, however, it will
! be discarded and a new, unique key will be generated in the set routine.
if (is_asciifile) then
   read(ifile, *) oldkey
else
   read(ifile) oldkey
endif

! Generate new unique radial velocity observation key, and set the contents
! of the private defined type.
call set_hf_radial_vel(velkey, instrument_id, beam_angle )

end subroutine read_hf_radial_vel


!----------------------------------------------------------------------


subroutine write_hf_radial_vel(velkey, ifile, fform)

! Write radial velocity auxiliary information to the obs_seq file.

integer,          intent(in)           :: velkey, ifile
character(len=*), intent(in), optional :: fform

logical  :: is_asciifile
integer  :: instrument_id
real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey,'write_hf_radial_vel')

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Write the 7 character identifier for verbose formatted output
   write(ifile, "('HFRADAR')")
endif

! Extract the values for this key and call the appropriate write routines.
instrument_id = radial_vel_metadata(velkey)%instrument_id
beam_angle    = radial_vel_metadata(velkey)%beam_angle

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

end subroutine write_hf_radial_vel


!----------------------------------------------------------------------


function read_instrument_id(ifile, is_asciiformat)

! Reads instrument id from file that was written by write_instrument_id.
! See read_hf_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
integer             :: read_instrument_id

integer            :: instrument_id

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   read(ifile, "(a)" ) header
   str1 = adjustl(header)

   if(str1(1:5) /= 'InsID') then
      write(msgstring,*)"Expected Instrument ID header 'InsID' in input file, got ", str1(1:5)
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
! See read_hf_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
real(r8)            :: read_beam_angle

real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

beam_angle = missing_r8

if (is_asciiformat) then
   read(ifile, "(a)" ) header
   str1 = adjustl(header)

   if(str1(1:5) /= 'angle') then
      write(msgstring,*)"Expected beam_angle header 'angle' in input file, got ", str1(1:5)
      call error_handler(E_ERR, 'read_beam_angle', msgstring, source, revision, revdate)  
   endif
   ! Now read the beam_angle data value into temporaries
   read(ifile, *) beam_angle
else
   ! No header label, just the binary angle values. 
   read(ifile)    beam_angle
endif

! Check for illegal values 
if (beam_angle <= 0.0_r8 .or. beam_angle >= 360.0_r8 ) then
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

! Writes beam_angle to obs file. Internally, the angles are in
! radians. When they are written, they are converted to degrees.

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


subroutine get_obs_def_hf_radial_vel(velkey, instrument_id, beam_angle)

! Return the auxiliary contents of a given radial velocity observation

integer,   intent(in)  :: velkey
integer,   intent(out) :: instrument_id
real(r8),  intent(out) :: beam_angle

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey,'get_obs_def_hf_radial_vel')

instrument_id = radial_vel_metadata(velkey)%instrument_id
beam_angle    = radial_vel_metadata(velkey)%beam_angle

end subroutine get_obs_def_hf_radial_vel


!----------------------------------------------------------------------


subroutine set_hf_radial_vel(velkey, instrument_id, beam_angle )

! Common code to increment the current key count, and set the private
! contents of this observation's auxiliary data.

integer,  intent(out) :: velkey
integer,  intent(in)  :: instrument_id
real(r8), intent(in)  :: beam_angle

if ( .not. module_initialized ) call initialize_module

! The total velocity metadata key count from all sequences
metadata_index = metadata_index + 1
velkey         = metadata_index
   
! Simple error check on key number before accessing the array
! This errors out if too key value now too large.
call velkey_out_of_range(velkey,'set_hf_radial_vel')

radial_vel_metadata(velkey)%instrument_id = instrument_id
radial_vel_metadata(velkey)%beam_angle    = beam_angle
   
end subroutine set_hf_radial_vel 


!----------------------------------------------------------------------


subroutine interactive_hf_radial_vel(velkey)

! Interactively reads in auxiliary information for a radial velocity obs.

integer, intent(out) :: velkey

! Uses the local subroutines interactive_instrument_id and set_hf_radial_vel.
! See read_hf_radial_vel for more information.

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

call set_hf_radial_vel(velkey, instrument_id, beam_angle )

write(*, *) 
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_hf_radial_vel


!----------------------------------------------------------------------


subroutine interactive_instrument_id(instrument_id)

! Prompt for instrument_id.  Not used in this code, but carried along.

integer, intent(out) :: instrument_id

write(*, *) 'Input an integer instrument_id:'
read(*, *) instrument_id

end subroutine interactive_instrument_id


!----------------------------------------------------------------------


subroutine interactive_beam_angle(beam_angle)

! Prompt for beam angle information in degrees.

real(r8), intent(out) :: beam_angle

beam_angle = missing_r8
do while (beam_angle <= 0.0_r8 .or. beam_angle >= 360.0_r8) 
   write(*, *) 'Input the beam angle in degrees (0 <= beam_angle <= 360):'
   read(*, *) beam_angle
end do

end subroutine interactive_beam_angle


!----------------------------------------------------------------------


subroutine get_expected_hf_radial_vel(state_handle, ens_size, location, velkey, &
                                      radial_vel, istatus)

! This is the main forward operator routine for radar Doppler velocity.

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: velkey
real(r8),            intent(out) :: radial_vel(ens_size)
integer,             intent(out) :: istatus(ens_size)


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


real(r8) :: u(ens_size), v(ens_size)
real(r8) :: debug_location(3)
integer  :: u_istatus(ens_size), v_istatus(ens_size)
logical  :: return_now

! Start out assuming good istatus
istatus = 0

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey,'get_expected_hf_radial_vel')

call interpolate(state_handle, ens_size, location, QTY_U_CURRENT_COMPONENT, u, u_istatus)
call track_status(ens_size, u_istatus, radial_vel, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_V_CURRENT_COMPONENT, v, v_istatus)
call track_status(ens_size, v_istatus, radial_vel, istatus, return_now)
if (return_now) return

where (istatus == 0)  &
   radial_vel = u * cos(radial_vel_metadata(velkey)%beam_angle*deg2rad) + &
                v * sin(radial_vel_metadata(velkey)%beam_angle*deg2rad)

!> @todo this needs to use the error_handler() and set do_output(.true.)
!> save the original output status and reset at the end
if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'radial velocity key: ', velkey
   print *, 'instrument id: ', radial_vel_metadata(velkey)%instrument_id
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated u: ', u
   print *, 'interpolated v: ', v
   print *, 'angle (deg): ', radial_vel_metadata(velkey)%beam_angle
   print *, 'angle (rad): ', radial_vel_metadata(velkey)%beam_angle*deg2rad
   print *, 'final radial_vel: ', radial_vel
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_hf_radial_vel

!----------------------------------------------------------------------

subroutine velkey_out_of_range(velkey, callingroutine)

! Range check velkey and trigger a fatal error if larger than allocated array.

integer,          intent(in) :: velkey
character(len=*), intent(in) :: callingroutine

! fine -- no problem.
if (velkey <= max_radial_vel_obs) return

! Bad news.  Tell the user.
write(msgstring, *) 'velkey (',velkey,') exceeds max_radial_vel_obs (', &
                     max_radial_vel_obs,')'
call error_handler(E_ERR,'set_hf_radial_vel', &
          'Increase max_radial_vel_obs in namelist', &
          source, revision, revdate, text2=msgstring, text3=callingroutine)

end subroutine velkey_out_of_range

!----------------------------------------------------------------------

end module obs_def_ocean_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

