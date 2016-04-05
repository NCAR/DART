! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! RAW_STATE_VARIABLE,    KIND_RAW_STATE_VARIABLE, COMMON_CODE
! RAW_STATE_1D_INTEGRAL, KIND_1D_INTEGRAL
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &
!                                     interactive_1d_integral, get_expected_1d_integral
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RAW_STATE_1D_INTEGRAL)                                                         
!            call get_expected_1d_integral(state_ens_handle, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call read_1d_integral(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call write_1d_integral(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call interactive_1d_integral(obs_def%key)
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
use  assim_model_mod, only : interpolate_distrib
use   cov_cutoff_mod, only : comp_cov_factor
use data_structure_mod, only : ensemble_type

implicit none

! These are the required interfaces for an obs_def module.  
public :: write_1d_integral, read_1d_integral, &
          interactive_1d_integral, get_expected_1d_integral, &
          set_1d_integral

! Storage for the special information required for observations of this type
integer               :: num_1d_integral_obs = 0     ! current count of obs
integer               :: max_1d_integral_obs = 1000  ! allocation size limit
real(r8), allocatable :: half_width(:)         ! metadata storage
integer,  allocatable :: num_points(:)         ! ditto
integer,  allocatable :: localization_type(:)  ! ditto


! Set to .true. to get debugging output
logical :: debug = .false.


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

!! To enable the namelist, comment this in.  NML2
!namelist /one_d_integral_nml/  debug, max_1d_integral_obs

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

! Allocate space for the metadata
allocate(half_width(max_1d_integral_obs),  &
         num_points(max_1d_integral_obs),  &
         localization_type(max_1d_integral_obs))

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

subroutine get_expected_1d_integral(state_ens_handle, location, igrkey, val, istatus)
 type(ensemble_type), intent(in)  :: state_ens_handle
 type(location_type), intent(in)  :: location
 integer,             intent(in)  :: igrkey
 real(r8),            intent(out) :: val(:)
 integer,             intent(out) :: istatus(:)

! The forward operator interface for this type of observation.  It is
! called with a state vector, a location, and a key to identify which
! observation is being processed.  The return 'val' is the expected
! observation value, and istatus is the return code.  0 is ok, 
! > 0 signals an error, and < 0 values are reserved for system use.
! The call to 'interpolate_distrib()' below calls the forward operator in 
! whatever model this code has been compiled with.

type(location_type)   :: location2
integer               :: i
real(r8)              :: range, loc, bottom, dx, x
integer               :: j !< loop variable
real(r8), allocatable :: sum(:), dist(:), weight(:), weight_sum(:)
integer               :: n !< number of forward operators to do.

if ( .not. module_initialized ) call initialize_module

! Make sure key is within valid range
call check_valid_key(igrkey, 'GIVEN', 'get_expected_1d_integral')

n = size(val)

allocate(sum(n), dist(n), weight_sum(n), weight(n))

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
do i = 1, num_points(igrkey)
   x = bottom + (i - 1) * dx
   if(x > 1.0_r8) x = x - 1.0_r8
   if(debug) print*, 'location for int ', i, 'is ', x
   location2 = set_location(x)
   call interpolate_distrib(location2, 1, istatus, val, state_ens_handle)
   if(debug) print*, 'model forward operator for ', i, ' returns ', val

   do j = 1, n
      if (istatus(j) /= 0) then
         if(debug) print*, 'forward operator', i, 'returned error'
         val(j) = missing_r8
         !return !not returning for a single failure
      endif
      dist(j) = abs(loc - x)
      if(dist(j) > 0.5_r8) dist(j) = 1.0_r8 - dist(j)
   enddo

   do j = 1, size(val) ! does this vectorize?
      if(debug) print*, 'dist ', i, dist
         weight(j) = comp_cov_factor(dist(j), half_width(igrkey), &
            localization_override = localization_type(igrkey))
      if(debug) print*, 'weight ', i, weight
         sum(j) = sum(j) + weight(j) * val(j)
         weight_sum(j) = weight_sum(j) + weight(j)
   enddo
enddo

val = sum / weight_sum

deallocate(weight, weight_sum, dist, sum)

if(debug) print*, 'get_expected_1d_integral key is ', igrkey
if(debug) print*, 'metadata values are: ', half_width(igrkey), num_points(igrkey), localization_type(igrkey)
if(debug) print*, 'return value for forward operator is ', val
if(debug) print*, 'return status (0 good; >0 error; <0 reserved for system use) is ', istatus

end subroutine get_expected_1d_integral

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

subroutine check_valid_key(igrkey, what, fromwhere)
 integer, intent(in)          :: igrkey
 character(len=*), intent(in) :: what, fromwhere

! Internal subroutine that verifies that we haven't incremented the key value
! past the size of the allocated space, or that a routine hasn't been called
! with a out-of-range key (which would indicate an internal error of some kind).
! If an error is found, an fatal message is printed and this routine doesn't return.
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

end subroutine

!----------------------------------------------------------------------

end module obs_def_1d_state_mod

! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
