! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! RAW_STATE_VARIABLE,    QTY_STATE_VARIABLE,     COMMON_CODE
! RAW_TRACER_CONCENTRATION,       QTY_TRACER_CONCENTRATION,     COMMON_CODE
! RAW_TRACER_SOURCE,     QTY_TRACER_SOURCE, COMMON_CODE
! RAW_STATE_1D_INTEGRAL, QTY_1D_INTEGRAL
! RAW_STATE_VAR_POWER,   QTY_STATE_VAR_POWER
! LARGE_SCALE_STATE,     QTY_LARGE_SCALE_STATE,  COMMON_CODE
! SMALL_SCALE_STATE,     QTY_SMALL_SCALE_STATE,  COMMON_CODE
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &
!                                     interactive_1d_integral, get_expected_1d_integral, &
!                                     write_power, read_power, interactive_power, get_expected_power
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RAW_STATE_1D_INTEGRAL)                                                         
!            call get_expected_1d_integral(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
!         case(RAW_STATE_VAR_POWER)
!            call get_expected_power(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)          
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call read_1d_integral(obs_def%key, ifile, fform)
!      case(RAW_STATE_VAR_POWER)
!         call read_power(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call write_1d_integral(obs_def%key, ifile, fform)
!      case(RAW_STATE_VAR_POWER)
!         call write_power(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call interactive_1d_integral(obs_def%key)
!      case(RAW_STATE_VAR_POWER)
!         call interactive_power(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_1d_state_mod

! This code currently does not require a namelist, but to add one search for
! the string 'NML' and comment in the code lines in the 4 marked sections.

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module, error_handler,               &
                             E_ERR, E_MSG, ascii_file_format !,            & 
                             ! find_namelist_in_file, check_namelist_read, &
                             ! nmlfileunit, do_nml_file, do_nml_term
                             !! Routines only needed for namelist support. NML1
use     location_mod, only : location_type, set_location, get_location 
use     obs_kind_mod, only : QTY_STATE_VARIABLE
use  assim_model_mod, only : interpolate
use   cov_cutoff_mod, only : comp_cov_factor
use ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none

! These are the required interfaces for an obs_def module.  
public :: write_1d_integral, read_1d_integral, &
          interactive_1d_integral, get_expected_1d_integral, &
          set_1d_integral, &
          write_power, read_power, interactive_power, get_expected_power, set_power

! Storage for the special information required for observations of this type
integer               :: num_1d_integral_obs = 0       ! current count of obs
integer               :: max_1d_integral_obs = 100000  ! allocation size limit
real(r8), allocatable :: half_width(:)                 ! metadata storage
integer,  allocatable :: num_points(:)                 ! ditto
integer,  allocatable :: localization_type(:)          ! ditto

! Storage for the power forward operator
integer               :: num_power_obs = 0           ! current count of obs
integer               :: max_power_obs = 220000      ! allocation size limit
real(r8), allocatable :: power(:)                    ! metadata storage


! Set to .true. to get debugging output
logical :: debug = .false.


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

!! To enable the namelist, comment this in.  NML2
!namelist /one_d_integral_nml/  debug, max_1d_integral_obs, max_power_obs

contains

!----------------------------------------------------------------------

subroutine initialize_module

! Anything in this routine is executed exactly once.  It could be used
! to read in a namelist, precompute values needed later, etc. 
! Currently this routine logs the version information and allocates
! space for the metadata arrays.

!! Only needed if you've commented in the namelist code.  NML3
!integer :: iunit, io

! Execute the code here just once
if (module_initialized) return

! Logs version information to file.
call register_module(source, revision, revdate)

! To enable a namelist (e.g to set max number of obs, or to turn
! on the debug messages at run time, or whatever) search for 'NML'
! and comment in the 4 code sections marked with that string.  NML4
!
!! Read the namelist input
!call find_namelist_in_file("input.nml", "one_d_integral_nml", iunit)
!read(iunit, nml = one_d_integral_nml, iostat = io)
!call check_namelist_read(iunit, io, "one_d_integral_nml")
!
!! Record the namelist values used for the run 
!if (do_nml_file()) write(nmlfileunit, nml=one_d_integral_nml)
!if (do_nml_term()) write(     *     , nml=one_d_integral_nml)

! Allocate space for the metadata for the integral observation
allocate(half_width(max_1d_integral_obs),  &
         num_points(max_1d_integral_obs),  &
         localization_type(max_1d_integral_obs))
! Allocate space for metadata for the power observation
allocate(power(max_power_obs))

module_initialized = .true.
if(debug) print*, 'module initialized'

end subroutine initialize_module

!----------------------------------------------------------------------------

subroutine write_1d_integral(igrkey, ifile, fform)
 integer,          intent(in)           :: igrkey, ifile
 character(len=*), intent(in), optional :: fform

! Write out the additional data associated with this observation.
! The obs is identified by the incoming 'key' argument.

logical :: is_ascii

if ( .not. module_initialized ) call initialize_module

! Make sure key value is within valid range -- it will be used as an index below.
call check_valid_key(igrkey, 'GIVEN', 'write_1d_integral')

is_ascii = ascii_file_format(fform)
if(debug) print*, 'write_1d_integral: ascii format = ', is_ascii

! Write out the half_width, num_points, and localization_type for each  
! observation embedded in the observation.  The old key is written out
! for tracking/debug use if needed.

if (is_ascii) then
   write(ifile, *) half_width(igrkey), num_points(igrkey), localization_type(igrkey)
   write(ifile, *) igrkey
else
   write(ifile)    half_width(igrkey), num_points(igrkey), localization_type(igrkey)
   write(ifile)    igrkey
endif
if(debug) print*, 'writing out metadata for 1D integral obs ', igrkey
if(debug) print*, 'metadata values are: ', half_width(igrkey), num_points(igrkey), localization_type(igrkey)

end subroutine write_1d_integral
!----------------------------------------------------------------------------

subroutine write_power(powkey, ifile, fform)
 integer,          intent(in)           :: powkey, ifile
 character(len=*), intent(in), optional :: fform

! Write out the additional data associated with this observation.
! The obs is identified by the incoming 'key' argument.

logical :: is_ascii

if ( .not. module_initialized ) call initialize_module

! Make sure key value is within valid range -- it will be used as an index below.
call check_valid_key_power(powkey, 'GIVEN', 'write_power')

is_ascii = ascii_file_format(fform)
if(debug) print*, 'write_power: ascii format = ', is_ascii

! Write out the power for each observation embedded in the observation.  The old key is written out
! for tracking/debug use if needed.

if (is_ascii) then
   write(ifile, *) power(powkey)
   write(ifile, *) powkey
else
   write(ifile)    power(powkey)
   write(ifile)    powkey
endif
if(debug) print*, 'writing out metadata for power obs ', powkey
if(debug) print*, 'metadata value is: ', power(powkey)

end subroutine write_power

!----------------------------------------------------------------------

subroutine read_1d_integral(igrkey, ifile, fform)
 integer,          intent(out)          :: igrkey
 integer,          intent(in)           :: ifile
 character(len=*), intent(in), optional :: fform

! Read in the additional data associated with this observation.
! The key value in the file will be read and then discarded, and a new key
! will be generated based on the next available index in the metadata arrays.
! Notice that key is intent(out) here, not (in) as in some other routines.

logical            :: is_ascii
integer            :: ignored_igrkey

if ( .not. module_initialized ) call initialize_module

! Increment the counter so all key values are unique
num_1d_integral_obs = num_1d_integral_obs + 1

! Set the return value for the key, and use it as the index below
igrkey = num_1d_integral_obs

! Make sure key is within valid range
call check_valid_key(igrkey, 'GENERATED', 'read_1d_integral')

is_ascii = ascii_file_format(fform)
if(debug) print*, 'read_1d_integral: ascii format = ', is_ascii

! Read in the additional metadata for this observation, and discard the old key.
if (is_ascii) then
   read(ifile, *) half_width(igrkey), num_points(igrkey), localization_type(igrkey) 
   read(ifile, *) ignored_igrkey
else
   read(ifile)    half_width(igrkey), num_points(igrkey), localization_type(igrkey) 
   read(ifile)    ignored_igrkey
endif
if(debug) print*, 'read in metadata for 1D integral obs ', igrkey
if(debug) print*, 'metadata values are: ', half_width(igrkey), num_points(igrkey), localization_type(igrkey)

if(debug) print *, 'ignoring old key', ignored_igrkey
if(debug) print *, 'return key set to ', igrkey

end subroutine read_1d_integral
!----------------------------------------------------------------------

subroutine read_power(powkey, ifile, fform)
 integer,          intent(out)          :: powkey
 integer,          intent(in)           :: ifile
 character(len=*), intent(in), optional :: fform

! Read in the additional data associated with this observation.
! The key value in the file will be read and then discarded, and a new key
! will be generated based on the next available index in the metadata arrays.
! Notice that key is intent(out) here, not (in) as in some other routines.

logical            :: is_ascii
integer            :: ignored_powkey

if ( .not. module_initialized ) call initialize_module

! Increment the counter so all key values are unique
num_power_obs = num_power_obs + 1

! Set the return value for the key, and use it as the index below
powkey = num_power_obs

! Make sure key is within valid range
call check_valid_key_power(powkey, 'GENERATED', 'read_power')

is_ascii = ascii_file_format(fform)
if(debug) print*, 'read_power: ascii format = ', is_ascii

! Read in the additional metadata for this observation, and discard the old key.
if (is_ascii) then
   read(ifile, *) power(powkey)
   read(ifile, *) ignored_powkey
else
   read(ifile)    power(powkey)
   read(ifile)    ignored_powkey
endif
if(debug) print*, 'read in metadata for power obs ', powkey
if(debug) print*, 'metadata value is: ', power(powkey)

if(debug) print *, 'ignoring old key', ignored_powkey
if(debug) print *, 'return key set to ', powkey

end subroutine read_power

!----------------------------------------------------------------------

subroutine interactive_1d_integral(igrkey)
 integer, intent(out) :: igrkey

! Initializes the specialized part of a 1d_integral observation.
! A new key will be generated based on the next available index
! in the metadata arrays.


if ( .not. module_initialized ) call initialize_module

! Increment the counter so all key values are unique
num_1d_integral_obs = num_1d_integral_obs + 1

! Set the return value for the key, and use it as the index below
igrkey = num_1d_integral_obs

! Make sure key is within valid range
call check_valid_key(igrkey, 'GENERATED', 'interactive_1d_integral')

! Prompt for input for the three required metadata items
write(*, *) 'Creating an interactive_1d_integral observation'

write(*, *) 'Input half width of integral '
 read(*, *) half_width(igrkey)

write(*, *) 'Input the number of evaluation points (5-20 recommended) '
 read(*, *) num_points(igrkey)

write(*, *) 'Input localization type: 1=Gaspari-Cohn; 2=Boxcar; 3=Ramped Boxcar'
 read(*, *) localization_type(igrkey)

if(debug) print *, 'return key set to ', igrkey

end subroutine interactive_1d_integral
!----------------------------------------------------------------------

subroutine interactive_power(powkey)
 integer, intent(out) :: powkey

! Initializes the specialized part of a power observation.
! A new key will be generated based on the next available index
! in the metadata arrays.


if ( .not. module_initialized ) call initialize_module

! Increment the counter so all key values are unique
num_power_obs = num_power_obs + 1

! Set the return value for the key, and use it as the index below
powkey = num_power_obs

! Make sure key is within valid range
call check_valid_key_power(powkey, 'GENERATED', 'interactive_power')

! Prompt for input for the metadata item
write(*, *) 'Creating an interactive_power observation'

write(*, *) 'Input power '
 read(*, *) power(powkey)

if(debug) print *, 'return key set to ', powkey

end subroutine interactive_power

!----------------------------------------------------------------------

subroutine get_expected_1d_integral(state_handle, ens_size, location, igrkey, val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: igrkey
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

! The forward operator interface for this type of observation.  It is
! called with a state vector, a location, and a key to identify which
! observation is being processed.  The return 'val' is the expected
! observation value, and istatus is the return code.  0 is ok, 
! > 0 signals an error, and < 0 values are reserved for system use.
! The call to 'interpolate()' below calls the forward operator in
! whatever model this code has been compiled with.

type(location_type)   :: location2
integer               :: i
real(r8)              :: range, loc, bottom, dx, x
integer               :: j !< loop variable
integer               :: point_istatus(ens_size)
logical               :: return_now ! used to return early if an interpoaltion fails
real(r8), dimension(ens_size) :: sum, dist, weight, weight_sum

if ( .not. module_initialized ) call initialize_module

! Make sure key is within valid range
call check_valid_key(igrkey, 'GIVEN', 'get_expected_1d_integral')

! Figure out the total range of the integrated funtion (1 is max)
range = 4.0_r8 * half_width(igrkey)
if(range > 1.0_r8) range = 1.0_r8
if(debug) print*, 'range is ', range

! Get the location value
loc = get_location(location)
if(debug) print*, 'loc base is ', loc

! Compute the bottom and top of the range
bottom = loc - range / 2.0_r8
if(bottom < 0.0_r8) bottom = bottom + 1.0_r8
if(debug) print*, 'bottom is ', bottom

! Next figure out where to put all the points
dx = range / (num_points(igrkey) - 1)
if(debug) print*, 'dx is ', dx

! Loop to compute the value at each point, then multiply by localization
! to get weighted integral
sum = 0.0_r8
weight_sum = 0.0_r8
istatus(:) = 0 ! so you can track istatus from 1 to num_points

do i = 1, num_points(igrkey)
   x = bottom + (i - 1) * dx
   if(x > 1.0_r8) x = x - 1.0_r8
   location2 = set_location(x)
   call interpolate(state_handle, ens_size, location2, QTY_STATE_VARIABLE, val, point_istatus)
   call track_status(ens_size, point_istatus, val, istatus, return_now)
   do j = 1, ens_size
      if (istatus(j) == 0) then
         dist(j) = abs(loc - x)
         if(dist(j) > 0.5_r8) dist(j) = 1.0_r8 - dist(j)
         weight(j) = comp_cov_factor(dist(j), half_width(igrkey), &
            localization_override = localization_type(igrkey))
         sum(j) = sum(j) + weight(j) * val(j)
         weight_sum(j) = weight_sum(j) + weight(j)
      endif
   enddo

enddo

where (istatus == 0) val = sum / weight_sum

if(debug) print*, 'get_expected_1d_integral key is ', igrkey
if(debug) print*, 'metadata values are: ', half_width(igrkey), num_points(igrkey), localization_type(igrkey)
if(debug) print*, 'return value for forward operator is ', val
if(debug) print*, 'return status (0 good; >0 error; <0 reserved for system use) is ', istatus

end subroutine get_expected_1d_integral
!----------------------------------------------------------------------

subroutine get_expected_power(state_handle, ens_size, location, powkey, val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: powkey
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer :: i
real(r8) :: pval(ens_size)

! The forward operator interface for this type of observation.  It is
! called with a state vector, a location, and a key to identify which
! observation is being processed.  The return 'val' is the expected
! observation value, and istatus is the return code.  0 is ok, 
! > 0 signals an error, and < 0 values are reserved for system use.
! The call to 'interpolate()' below calls the forward operator in
! whatever model this code has been compiled with.

if ( .not. module_initialized ) call initialize_module

! Make sure key is within valid range
call check_valid_key_power(powkey, 'GIVEN', 'get_expected_power')

! Interpolate the raw state to the location for each ensemble member
call interpolate(state_handle, ens_size, location, QTY_STATE_VARIABLE, pval, istatus)

if(power(powkey) == int(power(powkey))) then
   ! Integer power, just use standard definition
   val = pval ** power(powkey)
else
   ! For non-integer powers, fix up values for negative bases
   do i = 1, ens_size
      if(pval(i) >= 0.0_r8) then
         val(i) = pval(i) ** power(powkey)
      else
         val(i) = -1.0_r8 * (-1.0_r8 * pval(i)) ** power(powkey)
      endif
   end do
endif

if(debug) print*, 'get_expected_power key is ', powkey
if(debug) print*, 'metadata value is: ', power(powkey)
if(debug) print*, 'return value for forward operator is ', val
if(debug) print*, 'return status (0 good; >0 error; <0 reserved for system use) is ', istatus

end subroutine get_expected_power

!----------------------------------------------------------------------

subroutine set_1d_integral(integral_half_width, num_eval_pts, localize_type, igrkey, istatus)

! inputs are: half width of integral 
!             the number of evaluation points (5-20 recommended) 
!             localization type: 1=Gaspari-Cohn; 2=Boxcar; 3=Ramped Boxcar

 real(r8), intent(in)  :: integral_half_width
 integer,  intent(in)  :: num_eval_pts
 integer,  intent(in)  :: localize_type
 integer,  intent(out) :: igrkey
 integer,  intent(out) :: istatus

! Available to be called by a program creating these types of observations.
! Notice that igrkey is intent(out) here, not (in) as in some other routines.
! Sets the additional metadata for this obs, increments the key, and returns
! the new value.  This key value should be set in the obs_def derived type by
! calling set_obs_def_key().   Notice that this is different from the main
! observation key, which all observation entries have.  This key is specific
! to this observation type and is used to index into the metadata for only
! this type of obs.

if ( .not. module_initialized ) call initialize_module

! Increment the counter so all key values are unique
num_1d_integral_obs = num_1d_integral_obs + 1

! Set the return value for the key, and use it as the index below
igrkey = num_1d_integral_obs

! Make sure key is within valid range
call check_valid_key(igrkey, 'GENERATED', 'set_1d_integral')

! Set the corresponding values in the module global arrays
half_width(igrkey) = integral_half_width
num_points(igrkey) = num_eval_pts
localization_type(igrkey) = localize_type

istatus = 0

if(debug) print*, 'setting metadata for 1D integral obs ', igrkey
if(debug) print*, 'metadata values are: ', half_width(igrkey), num_points(igrkey), localization_type(igrkey)

if(debug) print*, 'return key set to ', igrkey
if(debug) print*, 'return status (0 good; >0 error; <0 reserved for system use) is ', istatus

end subroutine set_1d_integral
!----------------------------------------------------------------------

subroutine set_power(power_in, powkey, istatus)

! input is:  the power to which the state variable interpolation should be raised

 real(r8), intent(in)  :: power_in
 integer,  intent(out) :: powkey
 integer,  intent(out) :: istatus

! Available to be called by a program creating these types of observations.
! Notice that powkey is intent(out) here, not (in) as in some other routines.
! Sets the additional metadata for this obs, increments the key, and returns
! the new value.  This key value should be set in the obs_def derived type by
! calling set_obs_def_key().   Notice that this is different from the main
! observation key, which all observation entries have.  This key is specific
! to this observation type and is used to index into the metadata for only
! this type of obs.

if ( .not. module_initialized ) call initialize_module

! Increment the counter so all key values are unique
num_power_obs = num_power_obs + 1

! Set the return value for the key, and use it as the index below
powkey = num_power_obs

! Make sure key is within valid range
call check_valid_key_power(powkey, 'GENERATED', 'set_power')

! Set the corresponding values in the module global arrays
power(powkey) = power_in

istatus = 0

if(debug) print*, 'setting metadata for power obs ', powkey
if(debug) print*, 'metadata value is: ', power(powkey)

if(debug) print*, 'return key set to ', powkey
if(debug) print*, 'return status (0 good; >0 error; <0 reserved for system use) is ', istatus

end subroutine set_power

!----------------------------------------------------------------------

subroutine check_valid_key(igrkey, what, fromwhere)
 integer, intent(in)          :: igrkey
 character(len=*), intent(in) :: what, fromwhere

! Internal subroutine that verifies that we haven't incremented the key value
! past the size of the allocated space, or that a routine hasn't been called
! with a out-of-range key (which would indicate an internal error of some kind).
! If an error is found, a fatal message is printed and this routine doesn't return.
! The 'what' argument is either 'GIVEN' for a key value that's passed in from
! another routine; or 'GENERATED' for one we have just made and are planning to
! return to the caller.  The 'fromwhere' argument is the name of the calling 
! subroutine so the error message can report where it was called from.

character(len=128) :: msgstring

if (igrkey <= 0 .or. igrkey > max_1d_integral_obs) then
   if (what == 'GENERATED' .and. igrkey > max_1d_integral_obs) then
      ! generating a new key and ran out of space
      write(msgstring, *)'Out of space, max_1d_integral_obs limit ',max_1d_integral_obs
      call error_handler(E_ERR,trim(fromwhere),msgstring,source,revision,revdate, &
                         text2='Increase value of max_1d_integral_obs in obs_def_1d_state_mod')
   else
      ! called with a bad key or a negative key generated somehow. "shouldn't happen".
      write(msgstring, *)'Key is ',igrkey,' must be between 1 and ',max_1d_integral_obs
      call error_handler(E_ERR,trim(fromwhere),msgstring,source,revision,revdate, &
                         text2='Internal error: Invalid key value in RAW_STATE_1D_INTEGRAL obs')
   endif
endif

end subroutine check_valid_key
!----------------------------------------------------------------------

subroutine check_valid_key_power(powkey, what, fromwhere)
 integer, intent(in)          :: powkey
 character(len=*), intent(in) :: what, fromwhere

! Internal subroutine that verifies that we haven't incremented the key value
! past the size of the allocated space, or that a routine hasn't been called
! with a out-of-range key (which would indicate an internal error of some kind).
! If an error is found, a fatal message is printed and this routine doesn't return.
! The 'what' argument is either 'GIVEN' for a key value that's passed in from
! another routine; or 'GENERATED' for one we have just made and are planning to
! return to the caller.  The 'fromwhere' argument is the name of the calling 
! subroutine so the error message can report where it was called from.

character(len=128) :: msgstring

if (powkey <= 0 .or. powkey > max_power_obs) then
   if (what == 'GENERATED' .and. powkey > max_power_obs) then
      ! generating a new key and ran out of space
      write(msgstring, *)'Out of space, max_power_obs limit ',max_power_obs
      call error_handler(E_ERR,trim(fromwhere),msgstring,source,revision,revdate, &
                         text2='Increase value of max_power_obs in obs_def_1d_state_mod')
   else
      ! called with a bad key or a negative key generated somehow. "shouldn't happen".
      write(msgstring, *)'Key is ',powkey,' must be between 1 and ',max_power_obs
      call error_handler(E_ERR,trim(fromwhere),msgstring,source,revision,revdate, &
                         text2='Internal error: Invalid key value in RAW_STATE_VAR_POWER obs')
   endif
endif

end subroutine check_valid_key_power

!----------------------------------------------------------------------

end module obs_def_1d_state_mod

! END DART PREPROCESS MODULE CODE

