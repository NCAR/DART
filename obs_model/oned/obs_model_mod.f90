module obs_model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use types_mod
use location_mod, only : location_type, interactive_location
use assim_model_mod, only : interpolate
use obs_kind_mod, only : obs_kind_type, interactive_kind, get_obs_kind

private

public take_obs, interactive_def

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

contains

!======================================================================

function take_obs(state_vector, location, obs_kind)
!--------------------------------------------------------------------
!

implicit none

real(r8) :: take_obs
real(r8), intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: obs_kind

! Initially, have only raw state observations implemented so obs_kind is
! irrelevant. Just do interpolation and return.

take_obs = interpolate(state_vector, location, get_obs_kind(obs_kind))

end function take_obs



subroutine interactive_def(location, obs_kind)
!---------------------------------------------------------------------
!
! Used for interactive creation of an obs_def. For now, an observation
! is uniquely defined by a location and a kind. The obs_model level
! knows how kinds and locations need to be combined if there are any
! intricacies (although for models to date there are none, 17 April, 2002).

implicit none

type(location_type), intent(out) :: location
type(obs_kind_type), intent(out) :: obs_kind

call interactive_location(location)
call interactive_kind(obs_kind)

end subroutine interactive_def


end module obs_model_mod
