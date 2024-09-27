! Copyright 2019 University Corporation for Atmospheric Research and 
! Colorado Department of Public Health and Environment.
!
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use 
! this file except in compliance with the License. You may obtain a copy of the 
! License at      http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
! CONDITIONS OF ANY KIND, either express or implied. See the License for the 
! specific language governing permissions and limitations under the License.
!
! Development of this code utilized the RMACC Summit supercomputer, which is 
! supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236),
! the University of Colorado Boulder, and Colorado State University.
! The Summit supercomputer is a joint effort of the University of Colorado Boulder
! and Colorado State University.

! BEGIN DART PREPROCESS KIND LIST
!INSTRUMENT_TEC_NADIR,               QTY_TOTAL_ELECTRON_CONTENT 
!INSTRUMENT_TEC_STATION,             QTY_TOTAL_ELECTRON_CONTENT 
!INSTRUMENT_TEC_LIMB,                QTY_TOTAL_ELECTRON_CONTENT 
! END DART PREPROCESS KIND LIST

! Forward operator for computing Total Electron Content.
!  initial implementation does only nadir integrations (direct downward
!   line from instrument to ground at a single lat/lon), using model levels.
!  additional options which should be implemented are integrations along
!   sight lines from a lat/lon on the ground to an instrument location, and
!   possibly limb measurements (instrument to instrument along a sight line).
!
! calls the model_interpolate() routine for each model level until an error
! (presumably because we ran out of levels) and computes an integral.

! the current code doesn't need additional metadata, but as soon as we
! do anything that's not a nadar observation we'll need the location of
! the satellite, so i'm leaving in the code that would allocate, read/write,
! and manage metadata.

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_tec_mod, only : write_station, read_station, interactive_station, &
!                              get_expected_station,                             &
!                              write_limb, read_limb, interactive_limb,          &
!                              get_expected_limb, get_expected_nadir
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(INSTRUMENT_TEC_NADIR)
!     call get_expected_nadir(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!  case(INSTRUMENT_TEC_STATION)
!     call get_expected_station(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!  case(INSTRUMENT_TEC_LIMB)
!     call get_expected_limb(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!  case(INSTRUMENT_TEC_NADIR)
!      continue
!  case(INSTRUMENT_TEC_STATION)
!      call read_station(obs_def%key, ifile, fform)
!  case(INSTRUMENT_TEC_LIMB)
!      call read_limb(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!  case(INSTRUMENT_TEC_NADIR)
!      continue
!  case(INSTRUMENT_TEC_STATION)
!      call write_station(obs_def%key, ifile, fform)
!  case(INSTRUMENT_TEC_LIMB)
!      call write_limb(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!  case(INSTRUMENT_TEC_NADIR)
!      continue
!  case(INSTRUMENT_TEC_STATION)
!      call interactive_station(obs_def%key)
!  case(INSTRUMENT_TEC_LIMB)
!      call interactive_limb(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_tec_mod


use        types_mod, only : r8, missing_r8, PI, deg2rad
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location, get_location, VERTISLEVEL, &
                             set_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_TOTAL_ELECTRON_CONTENT
use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: read_station, write_station, interactive_station,       &
          get_expected_station, get_obs_def_station, set_station, &
          read_limb, write_limb, interactive_limb, set_limb,      &
          get_expected_limb, get_obs_def_limb, get_expected_nadir

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_tec_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical :: module_initialized = .false.

character(len=512) :: string1, string2, string3 ! For error message content
character(len=32 ) :: header, str1

! Derived type to contain the metadata for more complicated geometries.
! Contains auxiliary information used to compute the forward operator.
! FIXME:  !! CURRENTLY UNUSED  until station or limb obs implemented !!

type station_data_type
   private
! FIXME: probably the location of the instrument, or the angle
! from the ground location (lat/lon) to the instrument
   integer   :: instrument_id       ! ID number for data source
   real(r8)  :: beam_angle          ! angle of beam (degrees)
   type(location_type) :: sat_location  ! where sensor is located
end type station_data_type

type limb_data_type
   private
! FIXME: probably the location of the instrument, or the angle
! from the ground location (lat/lon) to the instrument
   integer   :: instrument_id       ! ID number for data source
   real(r8)  :: beam_angle          ! angle of beam (degrees)
   type(location_type) :: sat_location1  ! where beam start is located
   type(location_type) :: sat_location2  ! where beam end is located
end type limb_data_type

! Module storage for radial velocity metadata, allocated in init routine.
! Cumulative index into radar metadata array
integer :: station_metadata_index = 0 
integer ::    limb_metadata_index = 0 
type(station_data_type), allocatable :: station_metadata(:)
type(limb_data_type),    allocatable ::    limb_metadata(:)

! namelist items
integer :: max_station_obs = 0   ! set to > 0 for testing
integer ::    max_limb_obs = 0   ! set to > 0 for testing
logical :: debug = .false.   ! set to .true. to enable debug printout

namelist /obs_def_tec_nml/ max_station_obs, max_limb_obs, debug

contains


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine initialize_module

! Called once to set values and allocate space

integer :: iunit, io, rc

write(string1,*) 'obs_def_tec_mod is a totally untested routine.'
write(string2,*) 'It is provided as an example,'
write(string3,*) 'but needs to be thoroughly tested for your application.'

call error_handler(E_ERR,'initialize_module',string1, &
           source, revision, revdate, text2=string2, text3=string3)

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_tec_nml", iunit)
read(iunit, nml = obs_def_tec_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_tec_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_tec_nml)
if (do_nml_term()) write(     *     , nml=obs_def_tec_nml)

! Allocate space for the auxiliary information associated with each obs
! This code must be placed after reading the namelist, so the user can
! increase or decrease the number of obs supported and use more or less
! memory at run time.
if (max_station_obs > 0) then
   allocate(station_metadata(max_station_obs), stat = rc)
   if (rc /= 0) then
      write(string1, *) 'initial allocation failed for station TEC obs data,', &
                          'itemcount = ', max_station_obs
      call error_handler(E_ERR,'initialize_module', string1, &
                         source, revision, revdate)
   endif
endif
if (max_limb_obs > 0) then
   allocate(limb_metadata(max_limb_obs), stat = rc)
   if (rc /= 0) then
      write(string1, *) 'initial allocation failed for limb TEC obs data,', &
                          'itemcount = ', max_limb_obs
      call error_handler(E_ERR,'initialize_module', string1, &
                         source, revision, revdate)
   endif
endif

end subroutine initialize_module


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  Nadir-at-tangent-point section
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine get_expected_nadir(state_handle, ens_size, location, teckey, &
                              tec_val, istatus)

! This is the main forward operator routine for nadir total electron content

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: teckey
real(r8),            intent(out) :: tec_val(ens_size)
integer,             intent(out) :: istatus(ens_size)


! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted total electron content that would be in the
! integrated column from an instrument looking straight down at the tangent point.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.
!

integer  :: l, nlevels
real(r8) :: tec(ens_size, 2000)
real(r8) :: loc_vals(3), debug_location(3)
type(location_type) :: probe
logical  :: return_now
integer  :: this_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

! loop over model levels, expect to fail when we reach the top of the model
! give it a (ridiculous) max number of times to loop so in case interpolate
! doesn't return an error, we don't stay here forever.
loc_vals = get_location(location)
nlevels = 0
levelloop: do l = 1, 2000 
   probe = set_location(loc_vals(1), loc_vals(2), real(l,R8), VERTISLEVEL)
   call interpolate(state_handle, ens_size, probe, QTY_TOTAL_ELECTRON_CONTENT, tec(:,l), this_istatus)
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (any(istatus /= 0)) exit levelloop
   nlevels = nlevels + 1
enddo levelloop

! if we didn't get through the loop even once, istatus is /= 0 and we
! can just return
if (nlevels == 0) return

! if we get to the end of the loop, something is probably seriously wrong.
! set a bad istatus and return.
if (nlevels == 2000) then
   tec_val = MISSING_R8
   istatus = 99
   return
endif

! FIXME: do a real integral here.
tec_val = sum(tec(:,1:nlevels)) / nlevels
 
! Good return code. 
istatus = 0

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated tec: ', tec_val
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_nadir


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Station velocity section
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine get_expected_station(state_handle, ens_size, location, teckey, &
                              tec_val, istatus)

! This is the main forward operator routine for station total electron content

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: teckey
real(r8),            intent(out) :: tec_val(ens_size)
integer,             intent(out) :: istatus(ens_size)


! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted total electron content that would be integrated
! along the line of sight from a station on the ground to the instrument.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.
!

integer  :: l, nlevels
real(r8) :: tec(ens_size,2000)
real(r8) :: loc_vals(3), debug_location(3)
type(location_type) :: probe
logical  :: return_now
integer  :: this_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

! loop over model levels, expect to fail when we reach the top of the model
! give it a (ridiculous) max number of times to loop so in case interpolate
! doesn't return an error, we don't stay here forever.
! FIXME: THIS NEEDS TO GET THE SATELLITE LOCATION AND COMPUTE A PATH
! FROM THE LOCATION AND INTEGRATE ALONG THAT LINE.  SO THE LOCATION WILL
! CHANGE WITH EACH HEIGHT.
loc_vals = get_location(location)
nlevels = 0
levelloop: do l = 1, 2000 
   probe = set_location(loc_vals(1), loc_vals(2), real(l,R8), VERTISLEVEL)
   call interpolate(state_handle, ens_size, location, QTY_TOTAL_ELECTRON_CONTENT, tec(:,l), this_istatus)
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (any(istatus /= 0)) exit levelloop
   nlevels = nlevels + 1
enddo levelloop

! if we didn't get through the loop even once, istatus is /= 0 and we
! can just return
if (nlevels == 0) return

! if we get to the end of the loop, something is probably seriously wrong.
! set a bad istatus and return.
if (nlevels == 2000) then
   tec_val = MISSING_R8
   istatus = 99
   return
endif

! FIXME: do a real integral here.
tec_val = sum(tec(:,1:nlevels)) / nlevels
 
! Good return code. 
istatus = 0

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated tec: ', tec_val
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_station

!----------------------------------------------------------------------

subroutine read_station(teckey, ifile, fform)

! Main read subroutine for station TEC data observation

integer,          intent(out)          :: teckey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

logical  :: is_asciifile
integer  :: instrument_id, oldkey
real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Read the character identifier for verbose formatted output
   read(ifile, FMT="(a)") header
   str1 = adjustl(header)

   if(str1(1:7) /= 'STATION') then
      write(string1,*)"Expected header 'STATION' in input file, got ", str1(1:7)
      call error_handler(E_ERR,'read_station', &
           string1, source, revision, revdate)
   endif
endif

! read_location is a DART library routine that expects an optional string
! arg, but the other two read routines are local to this module and we can
! tell them exactly what format to be reading because we already know it.
!instrument_id = read_instrument_id (ifile, is_asciifile)
!beam_angle    = read_beam_angle    (ifile, is_asciifile)

! Read in the teckey for this particular observation, however, it will
! be discarded and a new, unique key will be generated in the set routine.
if (is_asciifile) then
   read(ifile, *) oldkey
else
   read(ifile) oldkey
endif

! Generate new unique station observation key, and set the contents
! of the private defined type.
call set_station(teckey, instrument_id, beam_angle)

end subroutine read_station


!----------------------------------------------------------------------


subroutine write_station(teckey, ifile, fform)

! Write station auxiliary information to the obs_seq file.

integer,          intent(in)           :: teckey, ifile
character(len=*), intent(in), optional :: fform

logical  :: is_asciifile
integer  :: instrument_id
real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call teckey_out_of_range(teckey,'write_station')

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Write the 7 character identifier for verbose formatted output
   write(ifile, "('STATION')")
endif

! Extract the values for this key and call the appropriate write routines.
!instrument_id = station_metadata(teckey)%instrument_id
!beam_angle    = station_metadata(teckey)%beam_angle

! These two write routines are local to this module, and we have already 
! figured out if it is a unformatted/binary file or formatted/ascii, so 
! go ahead and pass that info directly down to the routines.
call write_instrument_id(ifile, instrument_id, is_asciifile)
call write_beam_angle(ifile, beam_angle, is_asciifile)

! Write out the teckey used for this run, however this will be discarded
! when this observation is read in and a new key will be generated.
! (It may be useful for correlating error messages or identifying particular
! observations so we are leaving it as part of the aux data.)
if (is_asciifile) then
   write(ifile, *) teckey
else
   write(ifile) teckey
endif

end subroutine write_station


!----------------------------------------------------------------------


subroutine interactive_station(teckey)

! Interactively reads in auxiliary information for a radial velocity obs.

integer, intent(out) :: teckey

! Uses the local subroutines interactive_instrument_id and set_hf_radial_vel.
! See read_hf_radial_vel for more information.

! teckey is internally incremented in the set routine, and only counts 
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
write(*, *) 'Beginning to inquire for information on beam angle.'
write(*, *)

call interactive_beam_angle(beam_angle)

call set_station(teckey, instrument_id, beam_angle )

write(*, *) 
write(*, *) 'End of specialized section for station TEC.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_station


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Limb velocity section
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine get_expected_limb(state_handle, ens_size, location, teckey, &
                              tec_val, istatus)

! This is the main forward operator routine for limb total electron content

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: teckey
real(r8),            intent(out) :: tec_val(ens_size)
integer,             intent(out) :: istatus(ens_size)


! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted total electron content that would be integrated
! along the line of sight from a limb on the ground to the instrument.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.
!

integer  :: l, nlevels
real(r8) :: tec(ens_size,2000)
real(r8) :: loc_vals(3), debug_location(3)
type(location_type) :: probe
logical  :: return_now
integer  :: this_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

! loop over model levels, expect to fail when we reach the top of the model
! give it a (ridiculous) max number of times to loop so in case interpolate
! doesn't return an error, we don't stay here forever.
! FIXME: THIS NEEDS TO GET THE SATELLITE LOCATION AND COMPUTE A PATH
! FROM THE LOCATION AND INTEGRATE ALONG THAT LINE.  SO THE LOCATION WILL
! CHANGE WITH EACH HEIGHT.
loc_vals = get_location(location)
nlevels = 0
levelloop: do l = 1, 2000 
   probe = set_location(loc_vals(1), loc_vals(2), real(l,R8), VERTISLEVEL)
   call interpolate(state_handle, ens_size, location, QTY_TOTAL_ELECTRON_CONTENT, tec(:,l), this_istatus)
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (any(istatus /= 0)) exit levelloop
   nlevels = nlevels + 1
enddo levelloop

! if we didn't get through the loop even once, istatus is /= 0 and we
! can just return
if (nlevels == 0) return

! if we get to the end of the loop, something is probably seriously wrong.
! set a bad istatus and return.
if (nlevels == 2000) then
   tec_val = MISSING_R8
   istatus = 99
   return
endif

! FIXME: do a real integral here.
tec_val = sum(tec(:,1:nlevels)) / nlevels
 
! Good return code. 
istatus = 0

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated tec: ', tec_val
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_limb

!----------------------------------------------------------------------

subroutine read_limb(teckey, ifile, fform)

! Main read subroutine for limb TEC data observation

integer,          intent(out)          :: teckey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

logical  :: is_asciifile
integer  :: instrument_id, oldkey
real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Read the character identifier for verbose formatted output
   read(ifile, FMT="(a)") header
   str1 = adjustl(header)

   if(str1(1:4) /= 'LIMB') then
      write(string1,*)"Expected header 'LIMB' in input file, got ", str1(1:4)
      call error_handler(E_ERR,'read_limb', &
           string1, source, revision, revdate)
   endif
endif

! read_location is a DART library routine that expects an optional string
! arg, but the other two read routines are local to this module and we can
! tell them exactly what format to be reading because we already know it.
!instrument_id = read_instrument_id (ifile, is_asciifile)
!beam_angle    = read_beam_angle    (ifile, is_asciifile)

! Read in the teckey for this particular observation, however, it will
! be discarded and a new, unique key will be generated in the set routine.
if (is_asciifile) then
   read(ifile, *) oldkey
else
   read(ifile) oldkey
endif

! Generate new unique limb observation key, and set the contents
! of the private defined type.
call set_limb(teckey, instrument_id, beam_angle)

end subroutine read_limb


!----------------------------------------------------------------------


subroutine write_limb(teckey, ifile, fform)

! Write limb auxiliary information to the obs_seq file.

integer,          intent(in)           :: teckey, ifile
character(len=*), intent(in), optional :: fform

logical  :: is_asciifile
integer  :: instrument_id
real(r8) :: beam_angle

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call teckey_out_of_range(teckey,'write_limb')

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Write the 4 character identifier for verbose formatted output
   write(ifile, "('LIMB')")
endif

! Extract the values for this key and call the appropriate write routines.
!instrument_id = limb_metadata(teckey)%instrument_id
!beam_angle    = limb_metadata(teckey)%beam_angle

! These two write routines are local to this module, and we have already 
! figured out if it is a unformatted/binary file or formatted/ascii, so 
! go ahead and pass that info directly down to the routines.
call write_instrument_id(ifile, instrument_id, is_asciifile)
call write_beam_angle(ifile, beam_angle, is_asciifile)

! Write out the teckey used for this run, however this will be discarded
! when this observation is read in and a new key will be generated.
! (It may be useful for correlating error messages or identifying particular
! observations so we are leaving it as part of the aux data.)
if (is_asciifile) then
   write(ifile, *) teckey
else
   write(ifile) teckey
endif

end subroutine write_limb


!----------------------------------------------------------------------


subroutine interactive_limb(teckey)

! Interactively reads in auxiliary information for a radial velocity obs.

integer, intent(out) :: teckey

! Uses the local subroutines interactive_instrument_id and set_hf_radial_vel.
! See read_hf_radial_vel for more information.

! teckey is internally incremented in the set routine, and only counts 
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
write(*, *) 'Beginning to inquire for information on beam angle.'
write(*, *)

call interactive_beam_angle(beam_angle)

call set_limb(teckey, instrument_id, beam_angle )

write(*, *) 
write(*, *) 'End of specialized section for limb TEC.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_limb


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! other supporting routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------


function read_instrument_id(ifile, is_asciiformat)

! Reads instrument id from file that was written by write_instrument_id.
! See read_station for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
integer             :: read_instrument_id

integer            :: instrument_id

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   read(ifile, "(a)" ) header
   str1 = adjustl(header)

   if(str1(1:5) /= 'InsID') then
      write(string1,*)"Expected Instrument ID header 'InsID' in input file, got ", str1(1:5)
      call error_handler(E_ERR, 'read_instrument_id', string1, source, revision, revdate)  
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
! See read_station for additional discussion.

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
      write(string1,*)"Expected beam_angle header 'angle' in input file, got ", str1(1:5)
      call error_handler(E_ERR, 'read_beam_angle', string1, source, revision, revdate)  
   endif
   ! Now read the beam_angle data value into temporaries
   read(ifile, *) beam_angle
else
   ! No header label, just the binary angle values. 
   read(ifile)    beam_angle
endif

! Check for illegal values 
if (beam_angle <= 0.0_r8 .or. beam_angle >= 360.0_r8 ) then
   write(string1,*) "beam_angle value must be between 0 and 360, got: ", &
                       beam_angle
   call error_handler(E_ERR, 'read_beam_angle', string1, &
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


subroutine get_obs_def_station(teckey, instrument_id, beam_angle)

! Return the auxiliary contents of a given radial velocity observation

integer,   intent(in)  :: teckey
integer,   intent(out) :: instrument_id
real(r8),  intent(out) :: beam_angle

! Simple error check on key number before accessing the array
call teckey_out_of_range(teckey,'get_obs_def_station')

!instrument_id = station_metadata(teckey)%instrument_id
!beam_angle    = station_metadata(teckey)%beam_angle

end subroutine get_obs_def_station


!----------------------------------------------------------------------


subroutine set_station(teckey, instrument_id, beam_angle )

! Common code to increment the current key count, and set the private
! contents of this observation's auxiliary data.

integer,  intent(out) :: teckey
integer,  intent(in)  :: instrument_id
real(r8), intent(in)  :: beam_angle

if ( .not. module_initialized ) call initialize_module

! The total velocity metadata key count from all sequences
station_metadata_index = station_metadata_index + 1
teckey                 = station_metadata_index
   
! Simple error check on key number before accessing the array
! This errors out if too key value now too large.
call teckey_out_of_range(teckey,'set_station')

!station_metadata(teckey)%instrument_id = instrument_id
!station_metadata(teckey)%beam_angle    = beam_angle
   
end subroutine set_station 



!----------------------------------------------------------------------


subroutine get_obs_def_limb(teckey, instrument_id, beam_angle)

! Return the auxiliary contents of a given radial velocity observation

integer,   intent(in)  :: teckey
integer,   intent(out) :: instrument_id
real(r8),  intent(out) :: beam_angle

! Simple error check on key number before accessing the array
call teckey_out_of_range(teckey,'get_obs_def_limb')

!instrument_id = limb_metadata(teckey)%instrument_id
!beam_angle    = limb_metadata(teckey)%beam_angle

end subroutine get_obs_def_limb


!----------------------------------------------------------------------


subroutine set_limb(teckey, instrument_id, beam_angle )

! Common code to increment the current key count, and set the private
! contents of this observation's auxiliary data.

integer,  intent(out) :: teckey
integer,  intent(in)  :: instrument_id
real(r8), intent(in)  :: beam_angle

if ( .not. module_initialized ) call initialize_module

! The total velocity metadata key count from all sequences
limb_metadata_index = limb_metadata_index + 1
teckey              = limb_metadata_index
   
! Simple error check on key number before accessing the array
! This errors out if too key value now too large.
call teckey_out_of_range(teckey,'set_limb')

!limb_metadata(teckey)%instrument_id = instrument_id
!limb_metadata(teckey)%beam_angle    = beam_angle
   
end subroutine set_limb 



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

subroutine teckey_out_of_range(teckey, callingroutine)

! Range check teckey and trigger a fatal error if larger than allocated array.

integer,          intent(in) :: teckey
character(len=*), intent(in) :: callingroutine

! fine -- no problem.
if (teckey <= max_station_obs) return

! Bad news.  Tell the user.
write(string1, *) 'teckey (',teckey,') exceeds max_station_obs (', &
                     max_station_obs,')'
call error_handler(E_ERR,'set_station', &
          'Increase max_station_obs in namelist', &
          source, revision, revdate, text2=string1, text3=callingroutine)

end subroutine teckey_out_of_range

!----------------------------------------------------------------------

end module obs_def_tec_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

