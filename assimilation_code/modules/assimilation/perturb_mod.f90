! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! Create an ensemble by perturbing a single state instance.
!
! This module has two layers:
!
!   perturb_ensemble  - the handle level entry point.  It owns the policy:
!                       copying member 1 to the other members, saving and
!                       restoring missing_r8, choosing the perturbation
!                       method and falling back when the model provides no
!                       pert_model_copies.  Callers pass their own namelist
!                       values in as arguments; this module reads no namelist
!                       of its own because its callers read different ones.
!
!   perturb_uniform,  - the array level routines.  These take the copies array
!   perturb_bitwise     and the ensemble distribution as plain data so they
!                       can be tested without an ensemble handle.
!
! Note that ens_size is always an explicit argument.  It is NOT size(copies,1),
! because filter carries extra copies (mean, sd, inflation) beyond the ensemble
! members.

module perturb_mod

use types_mod,            only : r8, i8, missing_r8

use utilities_mod,        only : error_handler, E_ERR

use options_mod,          only : get_missing_ok_status

use assim_model_mod,      only : pert_model_copies

use ensemble_manager_mod, only : ensemble_type

use random_seq_mod,       only : random_seq_type, init_random_seq, random_gaussian

implicit none
private

public :: perturb_ensemble, &
          perturb_uniform,  &
          perturb_bitwise,  &
          PERT_MODEL_MOD,   &
          PERT_UNIFORM

character(len=*), parameter :: source = 'perturb_mod.f90'

! Perturbation methods
integer, parameter :: PERT_MODEL_MOD = 1 ! model_mod does it, bitwise fallback
integer, parameter :: PERT_UNIFORM   = 2 ! uniform ladder, for localization tests

contains

!------------------------------------------------------------------
!> Create an ensemble from the single instance held in copy 1.
!>
!> Copy 1 is copied to the other members and then perturbed.  For
!> PERT_MODEL_MOD the model_mod is given the chance to do the perturbing
!> via pert_model_copies; if it does not provide that interface,
!> perturb_bitwise is used instead.
!>   if no model perturb is provided, perturb_copies_task_bitwise is called.
!> Note: Not enforcing a model_mod to produce a
!> pert_model_copies that is bitwise across any number of
!> tasks, although there is enough information in the
!> ens_handle to do this.
!>
!> Some models allow missing_r8 in the state vector.  If missing_r8 is
!> allowed the locations of missing_r8s are stored before the perturb,
!> then the missing_r8s are put back in after the perturb.

subroutine perturb_ensemble(ens_handle, ens_size, amplitude, method)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: amplitude
integer,             intent(in)    :: method

integer               :: i
logical               :: interf_provided ! model does the perturbing
logical, allocatable  :: miss_me(:)
integer               :: partial_state_on_my_task ! the number of elements ON THIS TASK

! Copy from ensemble member 1 to the other copies
do i = 1, ens_handle%my_num_vars
   ens_handle%copies(2:ens_size, i) = ens_handle%copies(1, i)  ! How slow is this?
enddo

! If the state allows missing values, we have to record their locations
! and restore them in all the new perturbed copies.

if (get_missing_ok_status()) then
   partial_state_on_my_task = size(ens_handle%copies,2)
   allocate(miss_me(partial_state_on_my_task))
   miss_me = .false.
   where(ens_handle%copies(1, :) == missing_r8) miss_me = .true.
endif

select case (method)

   case (PERT_UNIFORM)
      call perturb_uniform(ens_handle%copies, ens_size, amplitude)

   case (PERT_MODEL_MOD)
      ! Let model do perturbations if it is prepared to do so
      call pert_model_copies(ens_handle, ens_size, amplitude, interf_provided)
      if (.not. interf_provided) then
         call perturb_bitwise(ens_handle%copies, ens_size, ens_handle%my_vars, &
                              ens_handle%num_vars, amplitude)
      endif

   case default
      call error_handler(E_ERR,'perturb_ensemble', &
         'unknown perturbation method', source)

end select

! Restore the missing_r8
if (get_missing_ok_status()) then
   do i = 1, ens_size
      where(miss_me) ens_handle%copies(i, :) = missing_r8
   enddo
   deallocate(miss_me)
endif

end subroutine perturb_ensemble

!------------------------------------------------------------------
!> Perturb the copies array uniformly for localization tests.
!>
!> Every state variable is given the same set of member offsets, so the
!> ensemble deviations do not depend on the state index.  Every state
!> variable then has the same spread, and the covariance between an
!> observation and any state variable is identical.  That leaves the
!> localization as the only source of variability in the increments.
!>
!> The perturbed values are x_j = x_1 + (j-1)*amplitude, so amplitude is a
!> spacing rather than a standard deviation.  The resulting spread is
!> amplitude*sqrt(N(N+1)/12) and grows with the ensemble size.

subroutine perturb_uniform(copies, ens_size, amplitude)

real(r8), intent(inout) :: copies(:,:) ! (num_copies, my_num_vars)
integer,  intent(in)    :: ens_size
real(r8), intent(in)    :: amplitude

integer :: i, j ! loop variables

do i = 1, size(copies, 2)
   do j = 2, ens_size
      copies(j, i) = copies(1, i) + (j - 1.0_r8) * amplitude
   end do
end do

end subroutine perturb_uniform

!------------------------------------------------------------------
!> Perturb the copies array in a way that is bitwise reproducible
!> no matter how many tasks you run on.
!>
!> Every task walks the whole global state and draws the same random
!> numbers, but only writes the elements it owns. Bitwise across task
!> counts, but the cost is proportional to the global state size on 
!> every task -> so does not scale at all
!>

subroutine perturb_bitwise(copies, ens_size, my_vars, num_vars, amplitude)

real(r8),    intent(inout) :: copies(:,:) ! (num_copies, my_num_vars)
integer,     intent(in)    :: ens_size
integer(i8), intent(in)    :: my_vars(:)  ! global index of each local element
integer(i8), intent(in)    :: num_vars    ! global number of state elements
real(r8),    intent(in)    :: amplitude

integer(i8)           :: i ! global state index
integer               :: j ! loop variable
type(random_seq_type) :: r(ens_size)
real(r8)              :: random_array(ens_size) ! array of random numbers
integer               :: local_index
integer               :: my_num_vars

my_num_vars = size(copies, 2)

! If a task owns no part of the state, so it has nothing to write.  Return
! before local_index is used to index copies, which would be out of bounds.
! I think DART currently does not allow num_procs > num_vars, but this is
! not an expensive safety check for perturbs (called once).
if (my_num_vars < 1) return

! Need ens_size random number sequences.
do j = 1, ens_size
   call init_random_seq(r(j), j)
enddo

local_index = 1 ! same across the ensemble

! Only one task is going to update per i.  This will not scale at all.
do i = 1, num_vars

   do j = 1, ens_size
     ! Can use copies here because the random number
     ! is only relevant to the task than owns element i.
     random_array(j)  =  random_gaussian(r(j), copies(j, local_index), amplitude)
   enddo

   if (my_vars(local_index) == i) then
      copies(1:ens_size, local_index) = random_array(:)
      local_index = local_index + 1 ! task is ready for the next random number
      local_index = min(local_index, my_num_vars)
   endif

enddo

end subroutine perturb_bitwise

!------------------------------------------------------------------

end module perturb_mod
