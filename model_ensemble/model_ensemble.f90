module model_ensemble_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

! Add use only clauses
use types_mod
use assim_model_mod, only : assim_model_type, get_state_vector_ptr, &
   init_assim_model
use time_manager_mod, only : time_type

private

public model_ensemble_type

! For efficiency, we risk the danger of direct pointers into the assim model
! storage for state. 
type model_ensemble_type
   integer :: ens_size
   type(assim_model_type), pointer :: member(:)
   type(state_pointer_type), pointer :: state_ptr(:)
end type model_ensemble_type

type state_pointer_type
   real(r8), pointer :: state(:)
end type state_pointer_type

contains

!================================================================================

subroutine init_model_ensemble(ens, ens_size)
!--------------------------------------------------------------------------------
!
! Initializes storage and direct pointers for an ens of ens_size

implicit none

type(model_ensemble_type), intent(out) :: ens
integer, intent(in) :: ens_size

integer :: i

! Set size
ens%ens_size = ens_size

! Allocate storage
allocate(ens%member(ens_size), ens%state_ptr(ens_size))

! Begin by initializing an assim_model instance for each ensemble member
! CAREFUL about having assim_model initialized, should check?
do i = 1, ens_size
   call init_assim_model(ens%member(i))
! Get direct pointer to model state for this ensemble
   ens%state_ptr(i)%state => get_state_vector_ptr(ens%member(i))
end do

end subroutine init_model_ensemble



subroutine get_ens(ens, i,  x)
!---------------------------------------------------------------------------------
!
! Copies a complete assim_model state into an ensemble member.

implicit none

type(model_ensemble_type), intent(in) :: ens
integer, intent(in) :: i
type(assim_model_type), intent(out) :: x

call copy_assim_model(x, ens%member(i))

end subroutine get_ens



subroutine set_ens(ens, i, x)
!------------------------------------------------------------------------------------
!
! Copies an assim_model complete state into an ensemble.

implicit none

type(model_ensemble_type), intent(inout) :: ens
integer, intent(in) :: i
type(assim_model_stype), intent(in) :: x

call copy_assim_model(ens%member(i), x)

end subroutine set_ens




function ensemble_mean(ens)
!---------------------------------------------------------------------------------
!
! Returns the ensemble mean state

implicit none

type(model_ensemble_type), intent(in) :: ens
! WARNING ON PG COMPILER: IT REQUIRES NEXT LINE TO BE AFTER ABOVE WHICH I THINK
! VIOLATES STANDARD. IS FOLLOWING LINE TECHNICALLY ILLEGAL?
real(r8) :: ensemble_mean(ens%ens_size)

integer :: i

ensemble_mean = ens%state_ptr(1)%state

do i = 2, ens%ens_size
   ensemble_mean = ensemble_mean + ens%state_ptr(i)%state
end do

ensemble_mean = ensemble_mean / ens%ens_size

end function ensemble_mean



function get_ens_state(ens, i, model_size)
!-------------------------------------------------------------------------------------
!
! Set the value of the state for ensemble member i. Might prefer not to have
! the model_size argument at some point, but it made initial implementation
! quick.

implicit none

real(r8) :: get_ens_state(model_size)
type(model_ensemble_type), intent(in) :: ens
integer, intent(in) :: i, model_size

! Do error checks on value of i needed
get_ens_state = ens%state_ptr(i)%state

end function get_ens_state



subroutine set_ens_state(ens, i, x)
!--------------------------------------------------------------------------------------
!
! Get the value of the state for ensemble member i

implicit none

type(model_ensemble_type), intent(inout) :: ens
integer, intent(in) :: i
real(r8), intent(in) :: x(:)

! Do error checks on value of i needed, also need to make sure x is right size?
ens%state_ptr(i)%state = x

end subroutine set_ens_state



function get_ens_members(ens, i)
!-------------------------------------------------------------------------------------
!
! Gets the value of all ensemble members for state variable index i

implicit none

type(model_ensemble_type), intent(in) :: ens
integer, intent(in) :: i
! I'd like to put this as first declaration, but PG won't let me. Is this an illegal
! statement?
real(r8) :: get_ens_members(ens%ens_size)

integer :: j

do j = 1, ens%ens_size
   get_ens_members(j) = ens%state_ptr(j)%state(i)
end do

end function get_ens_members



subroutine set_ens_members(ens, i, x)
!-------------------------------------------------------------------------------------
!
! Gets the value of all ensemble members for state variable index i

implicit none

type(model_ensemble_type), intent(inout) :: ens
integer, intent(in) :: i
real(r8), intent(in) :: x(ens%ens_size)

integer :: j

do j = 1, ens%ens_size
   ens%state_ptr(j)%state(i) = x(j)
end do

end subroutine set_ens_members



subroutine adv_ensemble(ens, time)
!---------------------------------------------------------------------------------------
!
! Advance all members of the ensemble to this time

implicit none

type(model_ensemble_type), intent(in) :: ens
type(time_type), intent(in) :: time

do i = 1, ens%ens_size
   call advance_state(ens%member(i), time)
end do

end subroutine adv_ensemble



function get_ens_time(ens, i)
!---------------------------------------------------------------------------------------
!
! Returns the time associated with the ith member of the ensemble. This raises some
! serious questions about synchronization of time in ensembles, etc.

implicit none

type(time_type) :: get_ens_time
type(model_ensemble_type), intent(in) :: ens
integer, intent(in) :: i

get_ens_time = get_model_time(ens%member(i))

end function get_ens_time



end module model_ensemble_mod

