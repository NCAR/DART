module obs_model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use types_mod
use location_mod
use assim_model_mod
use obs_kind_mod

private

public take_obs

contains

!======================================================================

function take_obs(state_vector, location, obs_kind)

implicit none

real(r8) :: take_obs
real(r8), intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: obs_kind

! Initially, have only raw state observations implemented so obs_kind is
! irrelevant. Just do interpolation and return.

take_obs = interpolate(state_vector, location)

end function take_obs


end module obs_model_mod
