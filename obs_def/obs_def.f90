module obs_def_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use types_mod
use obs_kind_mod, only : obs_kind_type, read_kind, write_kind
use location_mod, only : location_type, read_location, write_location
use obs_model_mod, only : take_obs, interactive_def
use assim_model_mod, only : am_get_close_states=>get_close_states, &
   am_get_num_close_states=>get_num_close_states

private

public init_obs_def, get_expected_obs, get_error_variance, get_obs_location, get_obs_def_kind, &
   get_num_close_states, get_close_states, set_obs_def_location, &
   set_err_var, set_obs_def_kind, read_obs_def, write_obs_def, obs_def_type, &
   interactive_obs_def

type obs_def_type
   private
   type(location_type) :: location
   type(obs_kind_type) :: kind
   real(r8) :: error_variance
end type obs_def_type

contains

!----------------------------------------------------------------------------

function init_obs_def(location, kind, error_variance)

! Constructor for an obs_def.

implicit none

type(obs_def_type) :: init_obs_def
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: kind
real(r8) :: error_variance

init_obs_def%location = location
init_obs_def%kind = kind
init_obs_def%error_variance = error_variance

end function init_obs_def

!----------------------------------------------------------------------------

function get_expected_obs(obs_def, state_vector)

! Given an obs_def and a state vector from a model returns the expected value
! of this observation (forward operator, h).

implicit none

real(r8) :: get_expected_obs
type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: state_vector(:)

! Need to figure out exactly how to expand this
get_expected_obs = take_obs(state_vector, obs_def%location, obs_def%kind)

end function get_expected_obs

!----------------------------------------------------------------------------

function get_error_variance(obs_def)

implicit none

real(r8) :: get_error_variance
type(obs_def_type), intent(in) :: obs_def

get_error_variance = obs_def%error_variance

end function get_error_variance

!----------------------------------------------------------------------------

function get_obs_location(obs_def)

! Returns observation location.

implicit none

type(location_type) :: get_obs_location
type(obs_def_type), intent(in) :: obs_def

get_obs_location = obs_def%location

end function get_obs_location

!----------------------------------------------------------------------------

function get_obs_def_kind(obs_def)

! Returns observation kind

implicit none

type(obs_kind_type) :: get_obs_def_kind
type(obs_def_type), intent(in) :: obs_def

get_obs_def_kind = obs_def%kind

end function get_obs_def_kind

!----------------------------------------------------------------------------

function get_num_close_states(obs_def, radius)

! Returns the number of state variables that are within distance radius of this
! obs_def location. This is a function of the class data for the state, not a 
! particular state so no state vector argument is needed. This limits things to
! one model per executable which might need to be generalized far in the future.
! F90 limitations make this difficult.

implicit none

integer :: get_num_close_states
type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: radius

! Call to assim_model level which knows how to work with locations
get_num_close_states = am_get_num_close_states(obs_def%location, radius)

end function get_num_close_states

!----------------------------------------------------------------------------

subroutine get_close_states(obs_def, radius, number, close_state_list)

! Returns the indices of those state variables that are within distance radius
! of the location of the obs_def along with the number of these states. In the 
! initial quick implementation, close_state_list is a fixed size real array
! and an error is returned by setting number to -number if the number of close states
! is larger than the array. Eventually may want to clean this up and make it 
! more efficient by allowing a dynamic storage allocation return.

implicit none

type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: radius
integer, intent(out) :: number
integer, intent(inout) :: close_state_list(:)

! For now, do this in inefficient redundant way; need to make more efficient soon
! NOTE: Could do the error checking on storage in assim_model if desired, probably
! have to do it there anyway.

number = get_num_close_states(obs_def, radius)

! Check for insufficient storage
if(number > size(close_state_list)) then
   number = -1 * number
   close_state_list = -1
else
   call am_get_close_states(obs_def%location, radius, number, close_state_list)
endif

end subroutine get_close_states

!----------------------------------------------------------------------------

function set_obs_location(obs_def, location)

! Sets the location of an obs_def

implicit none

type(obs_def_type) :: set_obs_location
type(obs_def_type), intent(in) :: obs_def
type(location_type), intent(in) :: location

set_obs_location = obs_def
set_obs_location%location = location

end function set_obs_location

!----------------------------------------------------------------------------

function set_error_variance(obs_def, error_variance)

! Sets the error variance of an obs_def

implicit none

type(obs_def_type) :: set_error_variance
type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: error_variance

set_error_variance = obs_def
set_error_variance%error_variance = error_variance

end function set_error_variance

!----------------------------------------------------------------------------


function set_obs_def_kind(obs_def, kind)

! Sets the kind of an obs_def

implicit none

type(obs_def_type) :: set_obs_def_kind
type(obs_def_type), intent(in) :: obs_def
type(obs_kind_type), intent(in) :: kind

set_obs_def_kind = obs_def
set_obs_def_kind%kind = kind

end function set_obs_def_kind

!----------------------------------------------------------------------------

function read_obs_def(file)

! Reads an obs_def from file which is just an integer unit number in the 
! current preliminary implementation.

implicit none

type(obs_def_type) :: read_obs_def
integer, intent(in) :: file

character*5 :: header

! Begin by reading five character ascii header, then location, kind, error variance

! Need to add additional error checks on read
read(file, 11) header
11 format(a5)
if(header /= 'obdef') then
   write(*, *) 'Error: Expected location header "obdef" in input file'
   stop
endif

! Read the location, kind and error variance
read_obs_def%location = read_location(file)
read_obs_def%kind = read_kind(file)
read(file, *) read_obs_def%error_variance

end function read_obs_def

!----------------------------------------------------------------------------

subroutine write_obs_def(file, obs_def)

! Writes an obs_def to file.

implicit none

integer, intent(in) :: file
type(obs_def_type), intent(in) :: obs_def

! Write the 5 character identifier
write(file, 11)
11 format('obdef')

! Write out the location, kind and error variance
call write_location(file, obs_def%location)
call write_kind(file, obs_def%kind)
write(file, *) obs_def%error_variance

end subroutine write_obs_def



function interactive_obs_def()
!---------------------------------------------------------------------------
!
! Allows interactive creation of an observation

implicit none

type(obs_def_type) :: interactive_obs_def

write(*, *) 'Input error variance for this observation definition '
read(*, *) interactive_obs_def%error_variance

! Call obs_model interfaces to set up kind and location
call interactive_def(interactive_obs_def%location, interactive_obs_def%kind)

end function interactive_obs_def



!----------------------------------------------------------------------------

end module obs_def_mod
!
