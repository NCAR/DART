module obs_set_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

use types_mod
use obs_set_def_mod
use time_manager_mod

private

public obs_set_type

type obs_set_type
   real(r8), pointer :: obs(:)
   integer :: num_copies
   type(time_type) :: time
   type(obs_set_def_type) :: def
end type obs_set_type


contains

!=======================================================

function get_obs_set_time(set)
!-------------------------------------------------------
!
! Returns the time associated with this observation set

implicit none

type(time_type) :: get_obs_set_time
type(obs_set_type), intent(in) :: set

get_obs_set_time = set%time

end function get_obs_set_time


end module obs_set_mod
