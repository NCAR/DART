! Aim: to provide interfaces to state_vector_ip
! for models that don't clamp or fail out of bounds
! variables
module null_clamp

use types_mod, only : r8

implicit none


contains

!-------------------------------------------------------
!> Null version
!> Check whether you need to error out, clamp, or
!> do nothing depending on the variable bounds
function do_clamp_or_fail(var, dom)

integer, intent(in) :: var ! variable index
integer, intent(in) :: dom ! domain index
logical             :: do_clamp_or_fail

do_clamp_or_fail = .false.

end function do_clamp_or_fail

!-------------------------------------------------------
!> Null version
!> Check a variable for out of bounds and clamp or fail if
!> needed
subroutine clamp_or_fail_it(var_index, dom, variable)

integer,     intent(in) :: var_index ! variable index
integer,     intent(in) :: dom ! domain index
real(r8), intent(inout) :: variable(:) ! variable


end subroutine clamp_or_fail_it


end module null_clamp