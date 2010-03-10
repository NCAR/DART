! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

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
!            call get_expected_1d_integral(state, location, obs_def%key, obs_val, istatus)  
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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location 
use  assim_model_mod, only : interpolate
use   cov_cutoff_mod, only : comp_cov_factor

implicit none

public :: write_1d_integral, read_1d_integral, interactive_1d_integral, &
          get_expected_1d_integral

! Storage for the special information required for observations of this type
integer, parameter                      :: max_1d_integral_obs = 100
integer                                 :: num_1d_integral_obs = 0
real(r8)                                :: half_width(max_1d_integral_obs)
integer, dimension(max_1d_integral_obs) :: num_points, localization_type

! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------

  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module



 subroutine write_1d_integral(key, ifile, fileformat)
!----------------------------------------------------------------------------
!subroutine write_1d_integral(key, ifile, fileformat)

integer, intent(in)             :: key, ifile
character(len=32), intent(in)   :: fileformat

integer :: i
logical :: is_ascii

if ( .not. module_initialized ) call initialize_module

is_ascii = ascii_file_format(fileformat)

! Philosophy, dump ALL information about this special obs_type at once???
! For now, this means you can only write ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this writing
if(.not. already_written) then
   already_written = .true.   
   ! Write out the number of 1d_integral obs descriptions
   if (is_ascii) then
      write(ifile, *) num_1d_integral_obs
   else
      write(ifile)    num_1d_integral_obs
   endif

   ! Write out the half_width, num_points, and localization_type for each  
   do i = 1, num_1d_integral_obs
      if (is_ascii) then
         write(ifile, *) half_width(i), num_points(i), localization_type(i) 
      else
         write(ifile)    half_width(i), num_points(i), localization_type(i) 
      endif
   end do
endif

! Write out the obs_def key for this observation
if (is_ascii) then
   write(ifile, *) key
else
   write(ifile)    key
endif

end subroutine write_1d_integral



 subroutine read_1d_integral(key, ifile, fileformat)
!----------------------------------------------------------------------
!subroutine read_1d_integral(key, ifile, fileformat)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=32), intent(in)   :: fileformat

integer :: i
logical :: is_ascii

if ( .not. module_initialized ) call initialize_module

is_ascii = ascii_file_format(fileformat)

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
if(.not. already_read) then
   already_read = .true.   
   ! Read the number of 1d_integral obs descriptions
   if (is_ascii) then
      read(ifile, *) num_1d_integral_obs
   else
      read(ifile)    num_1d_integral_obs
   endif
   
   ! Read the half_width, num_points, and localization_type for each  
   do i = 1, num_1d_integral_obs
      if (is_ascii) then
         read(ifile, *) half_width(i), num_points(i), localization_type(i) 
      else
         read(ifile)    half_width(i), num_points(i), localization_type(i) 
      endif
   end do
endif

! Read in the key for this particular observation
if (is_ascii) then
   read(ifile, *) key
else
   read(ifile)    key
endif

end subroutine read_1d_integral



 subroutine interactive_1d_integral(key)
!----------------------------------------------------------------------
!subroutine interactive_1d_integral(key)
!
! Initializes the specialized part of a 1d_integral observation
! Passes back up the key for this one

integer, intent(out) :: key

character(len=129) :: msgstring

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_1d_integral_obs >= max_1d_integral_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a 1d_integral_obs.'
   call error_handler(E_MSG,'interactive_1d_integral',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_1d_integral_obs (currently ',max_1d_integral_obs,')'
   call error_handler(E_ERR,'interactive_1d_integral',msgstring,source,revision,revdate)
endif

! Increment the index
num_1d_integral_obs = num_1d_integral_obs + 1
key = num_1d_integral_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_1d_integral observation'
write(*, *) 'Input half width of integral '
read(*, *) half_width(num_1d_integral_obs)
write(*, *) 'Input the number of evaluation points (??? recommended) '
read(*, *) num_points(num_1d_integral_obs)
write(*, *) 'Input localization type: 1=Gaspari-Cohn; 2=Boxcar; 3=Ramped Boxcar'
read(*, *) localization_type(num_1d_integral_obs)

end subroutine interactive_1d_integral



 subroutine get_expected_1d_integral(state, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_1d_integral(state, location, key, val, istatus)

real(r8), intent(in)            :: state(:)
type(location_type), intent(in) :: location
integer, intent(in)             :: key
real(r8), intent(out)           :: val
integer, intent(out)            :: istatus

integer :: i
real(r8) :: range, loc, bottom, dx, x, sum, dist, weight, weight_sum
type(location_type) :: location2

if ( .not. module_initialized ) call initialize_module

! Figure out the total range of the integrated funtion (1 is max)
range = 4.0_r8 * half_width(key)
if(range > 1.0_r8) range = 1.0_r8
!write(*, *) 'range is ', range

! Get the location value
loc = get_location(location)
!write(*, *) 'loc base is ', loc

! Compute the bottom and top of the range
bottom = loc - range / 2.0_r8
if(bottom < 0.0_r8) bottom = bottom + 1.0_r8
!write(*, *) 'bottom is ', bottom

! Next figure out where to put all the points
dx = range / (num_points(key) - 1)
!write(*, *) 'dx is ', dx

! Loop to compute the value at each point, then multiply by localization
! to get weighted integral
sum = 0.0_r8
weight_sum = 0.0_r8
do i = 1, num_points(key)
   x = bottom + (i - 1) * dx
   if(x > 1.0_r8) x = x - 1.0_r8
!write(*, *) 'location for int ', i, 'is ', x
   location2 = set_location(x)
   call interpolate(state, location2, 1, val, istatus)
   dist = abs(loc - x)
   if(dist > 0.5_r8) dist = 1.0_r8 - dist
!write(*, *) 'dist ', i, dist
   weight = comp_cov_factor(dist, half_width(key), &
      localization_override = localization_type(key))
!write(*, *) 'weight ', i, weight
   sum = sum + weight * val
   weight_sum = weight_sum + weight
end do

val = sum / weight_sum

!write(*, *) 'get_expected_1d_integral key is ', key
!write(*, *) half_width(key), num_points(key), localization_type(key)

end subroutine get_expected_1d_integral

!----------------------------------------------------------------------

end module obs_def_1d_state_mod
! END DART PREPROCESS MODULE CODE
