! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
      
module filter_io_diag_mod
     
!------------------------------------------------------------------------------

use types_mod,             only : r8, missing_r8

use ensemble_manager_mod,  only : ensemble_type, compute_copy_mean_sd

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use io_filenames_mod,      only : file_info_type, COPY_NOT_PRESENT

use assim_model_mod,       only : pert_model_copies

use state_vector_io_mod,   only : get_stage_to_write, write_state

!------------------------------------------------------------------------------

! Identifier for different copies for diagnostic files
integer, parameter :: ENS_START     = 1
integer, parameter :: ENS_END       = 2
integer, parameter :: ENS_MEAN      = 3
integer, parameter :: ENS_SD        = 4
integer, parameter :: PRIOR_INF     = 5
integer, parameter :: PRIOR_INF_SD  = 6
integer, parameter :: POST_INF      = 7
integer, parameter :: POST_INF_SD   = 8

! Data structure for different diag and output stages
integer, parameter :: NUM_SCOPIES    = 8
integer :: DIAG_FILE_COPIES( NUM_SCOPIES )     = COPY_NOT_PRESENT
                                  
!------------------------------------------------------------------------------

public ::  create_ensemble_from_single_file, do_stage_output, &
   ENS_START, ENS_END, ENS_MEAN, ENS_SD, PRIOR_INF, PRIOR_INF_SD, &
   POST_INF, POST_INF_SD, NUM_SCOPIES, DIAG_FILE_COPIES

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!> Produces an ensemble by copying my_vars of the 1st ensemble member
!> and then perturbing the copies array.
!> Mimicks the behaviour of pert_model_state:
!> pert_model_copies is called:
!>   if no model perturb is provided, perturb_copies_task_bitwise is called.
!> Note: Not enforcing a model_mod to produce a
!> pert_model_copies that is bitwise across any number of
!> tasks, although there is enough information in the
!> ens_handle to do this.
!>
!> Some models allow missing_r8 in the state vector.  If missing_r8 is
!> allowed the locations of missing_r8s are stored before the perturb,
!> then the missing_r8s are put back in after the perturb.

subroutine create_ensemble_from_single_file(ens_handle, ens_size, &
   perturbation_amplitude, missing_ok)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: perturbation_amplitude
logical,             intent(in)    :: missing_ok

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

if (missing_ok) then
   partial_state_on_my_task = size(ens_handle%copies,2)
   allocate(miss_me(partial_state_on_my_task))
   miss_me = .false.
   where(ens_handle%copies(1, :) == missing_r8) miss_me = .true.
endif

call pert_model_copies(ens_handle, ens_size, perturbation_amplitude, interf_provided)
if (.not. interf_provided) then
   call perturb_copies_task_bitwise(ens_handle, ens_size, perturbation_amplitude)
endif

! Restore the missing_r8
if (missing_ok) then
   do i = 1, ens_size
      where(miss_me) ens_handle%copies(i, :) = missing_r8
   enddo
   deallocate(miss_me)
endif

end subroutine create_ensemble_from_single_file

!------------------------------------------------------------------------------
! Perturb the copies array in a way that is bitwise reproducible
! no matter how many task you run on.

subroutine perturb_copies_task_bitwise(ens_handle, ens_size, perturbation_amplitude)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: perturbation_amplitude

integer               :: i, j ! loop variables
type(random_seq_type) :: r(ens_size)
real(r8)              :: random_array(ens_size) ! array of random numbers
integer               :: local_index

! Need ens_size random number sequences.
do i = 1, ens_size
   call init_random_seq(r(i), i)
enddo

local_index = 1 ! same across the ensemble

! Only one task is going to update per i.  This will not scale at all.
do i = 1, ens_handle%num_vars

   do j = 1, ens_size
     ! Can use %copies here because the random number
     ! is only relevant to the task than owns element i.
     random_array(j)  =  random_gaussian(r(j), ens_handle%copies(j, local_index), perturbation_amplitude)
   enddo

   if (ens_handle%my_vars(local_index) == i) then
      ens_handle%copies(1:ens_size, local_index) = random_array(:)
      local_index = local_index + 1 ! task is ready for the next random number
      local_index = min(local_index, ens_handle%my_num_vars)
   endif

enddo

end subroutine perturb_copies_task_bitwise

!------------------------------------------------------------------------------

subroutine do_stage_output(stage_name, output_interval, time_step_number, &
   state_ens_handle, file_info, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

character(len = *),   intent(in)    :: stage_name
integer,              intent(in)    :: output_interval, time_step_number
type(ensemble_type),  intent(inout) :: state_ens_handle
type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: ens_size, ENS_MEAN_COPY, ENS_SD_COPY

if(get_stage_to_write(stage_name)) then
   if((output_interval > 0) .and. &
      (time_step_number / output_interval * output_interval == time_step_number)) then

      ! Compute the ensemble mean and standard deviation
      ! For efficiency could make this optional in call
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      call write_state(state_ens_handle, file_info)
   endif
endif

end subroutine do_stage_output

!------------------------------------------------------------------------------

end module filter_io_diag_mod
