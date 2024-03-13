! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!------------------------------------------------------------------------------
!> forward_operator_mod.f90
!>
!> This module contains routines related to the forward operator.
!>
!------------------------------------------------------------------------------
module forward_operator_mod


use types_mod,             only : r8, i8, missing_r8

use time_manager_mod,      only : time_type

use utilities_mod,         only : error_handler, E_ERR

use mpi_utilities_mod,     only : my_task_id

use obs_sequence_mod,      only : init_obs, destroy_obs, obs_sequence_type, &
                                  obs_type, get_obs_values, get_qc,         &
                                  get_obs_def, get_obs_from_key

use obs_def_mod,           only : obs_def_type, get_obs_def_error_variance, &
                                  get_expected_obs_from_def_distrib_state,  &
                                  get_obs_def_type_of_obs 

use obs_kind_mod,          only : assimilate_this_type_of_obs, evaluate_this_type_of_obs

use ensemble_manager_mod,  only : ensemble_type, compute_copy_mean_var, &
                                  get_my_num_copies, copies_in_window,  &
                                  get_allow_transpose, all_vars_to_all_copies, &
                                  all_copies_to_all_vars, allocate_single_copy, &
                                  get_single_copy, put_single_copy, deallocate_single_copy

use distributed_state_mod, only : create_state_window, free_state_window,   &
                                  get_state

use ensemble_manager_mod,    only : copies_in_window

use quality_control_mod,   only : check_outlier_threshold, get_dart_qc, &
                                  input_qc_ok, good_dart_qc, DARTQC_BAD_INCOMING_QC, &
                                  DARTQC_ASSIM_GOOD_FOP, DARTQC_EVAL_GOOD_FOP, &
                                  DARTQC_NOT_IN_NAMELIST

!------------------------------------------------------------------------------

implicit none

private

public :: get_obs_ens_distrib_state, get_expected_obs_distrib_state

character(len=*), parameter :: source = 'forward_operator_mod.f90'

! Module storage for writing error messages
character(len=512) :: string1, string2

contains

!------------------------------------------------------------------------------
!> Subroutine get_obs_ens_distrib_state
!> 
!> Computes the forward observation operators and related quality 
!> control indicators.
!> 
!------------------------------------------------------------------------------
subroutine get_obs_ens_distrib_state(ens_handle, obs_fwd_op_ens_handle, &
   qc_ens_handle, seq, keys, obs_val_index, input_qc_index, &
   OBS_ERR_VAR_COPY,   OBS_VAL_COPY,   OBS_KEY_COPY, &
   OBS_GLOBAL_QC_COPY, OBS_EXTRA_QC_COPY, OBS_MEAN_COPY, &
   OBS_VAR_COPY, isprior, prior_qc_copy)

type(ensemble_type),     intent(inout) :: ens_handle  !! state ensemble handle
type(ensemble_type),     intent(inout) :: obs_fwd_op_ens_handle  !! observation forward operator handle
type(ensemble_type),     intent(inout) :: qc_ens_handle  !! quality control handle
type(obs_sequence_type), intent(in)    :: seq  !! the observation sequence
integer,                 intent(in)    :: keys(:)  !! observation key numbers
integer,                 intent(in)    :: obs_val_index  !! copy index for observation value 
integer,                 intent(in)    :: input_qc_index  !! qc index for incoming QC value
integer,                 intent(in)    :: OBS_ERR_VAR_COPY  !! ensemble copy number for observation error variances
integer,                 intent(in)    :: OBS_VAL_COPY  !! ensemble copy number for observation values
integer,                 intent(in)    :: OBS_KEY_COPY  !! ensemble copy number for observation keys
integer,                 intent(in)    :: OBS_GLOBAL_QC_COPY  !! ensemble copy number for ??
integer,                 intent(in)    :: OBS_EXTRA_QC_COPY  !! ensemble copy number for ??
integer,                 intent(in)    :: OBS_MEAN_COPY  !! ensemble copy number for obs mean values
integer,                 intent(in)    :: OBS_VAR_COPY  !! ensemble copy number for obs variance 
logical,                 intent(in)    :: isprior  !! true for prior eval; false for posterior
real(r8),                intent(inout) :: prior_qc_copy(:)  !! array instead of ensemble copy ??

real(r8) :: input_qc(1), obs_value(1), obs_err_var

real(r8), allocatable :: expected_obs(:) !Also regular obs now?

integer :: i, j, k !! index variables
integer :: thiskey(1)
integer :: global_obs_num, global_qc_value
integer :: num_copies_to_calc
integer :: global_ens_index
integer :: ens_size

integer, allocatable  :: istatus(:)
integer, allocatable  :: var_istatus(:)
integer, allocatable  :: my_copy_indices(:) ! The global ens index for each copy a task has

logical :: evaluate_this_ob, assimilate_this_ob

type(time_type) :: dummy_time

type(obs_def_type) :: obs_def
type(obs_type)     :: observation

! IMPORTANT, IT IS ASSUMED THAT ACTUAL ENSEMBLES COME FIRST
! It is also assumed that the ensemble members are in the same
! order in each of the handles

num_copies_to_calc = copies_in_window(ens_handle)

allocate(istatus(num_copies_to_calc))
allocate(expected_obs(num_copies_to_calc))
allocate(my_copy_indices(num_copies_to_calc))

istatus = 999123
expected_obs = MISSING_R8

! Set up access to the state
call create_state_window(ens_handle, obs_fwd_op_ens_handle, qc_ens_handle)

ens_size = ens_handle%num_copies - ens_handle%num_extras

if(get_allow_transpose(ens_handle)) then ! giant if for transpose or distributed forward op

   my_copy_indices(:) = ens_handle%my_copies(1:num_copies_to_calc) ! var-complete forward operators

   ! Loop through all observations in the set
   ALL_OBSERVATIONS: do j = 1, obs_fwd_op_ens_handle%num_vars
      ! Get the information on this observation by placing it in temporary
      call get_obs_from_key(seq, keys(j), observation)
      call get_obs_def(observation, obs_def)
      ! Check to see if this observation fails input qc test
      call get_qc(observation, input_qc, input_qc_index)

      ! Get the observation value and error variance
      call get_obs_values(observation, obs_value(1:1), obs_val_index)
      obs_err_var = get_obs_def_error_variance(obs_def)

      ! Loop through all copies stored by this process and set values as needed
      do k = 1, obs_fwd_op_ens_handle%my_num_copies
         ! See if this is the copy for error variance or observed value
         global_ens_index = obs_fwd_op_ens_handle%my_copies(k)
         if(global_ens_index == OBS_ERR_VAR_COPY) then
            ! This copy is the instrument observation error variance; read and
            ! store
            obs_fwd_op_ens_handle%vars(j, k) = obs_err_var

         else if(global_ens_index == OBS_VAL_COPY) then
            ! This copy is the observation from the instrument; read and store
            obs_fwd_op_ens_handle%vars(j, k) = obs_value(1)
         else if(global_ens_index == OBS_KEY_COPY) then
            obs_fwd_op_ens_handle%vars(j, k) = keys(j)
         else if(global_ens_index == OBS_GLOBAL_QC_COPY .and. isprior) then
            obs_fwd_op_ens_handle%vars(j, k) = 0.0_r8
         endif
      enddo

      ! If it is bad, set forward operator status value to -99 and return missing_r8 for obs_value
      ! PAR THIS SUBROUTINE SHOULD EVENTUALLY GO IN THE QUALITY CONTROL MODULE
      if(.not. input_qc_ok(input_qc(1), global_qc_value)) then

         qc_ens_handle%vars(j, :) = 0

         do k=1, obs_fwd_op_ens_handle%my_num_copies
            global_ens_index = obs_fwd_op_ens_handle%my_copies(k)
            ! Update prior/post obs values, mean, etc - but leave the key copy
            ! and the QC copy alone.
            if ((global_ens_index /= OBS_KEY_COPY) .and. &
               (global_ens_index /= OBS_EXTRA_QC_COPY) .and. &
               (global_ens_index /= OBS_GLOBAL_QC_COPY)) then
               ! JH test if the global_ens_index is our extra bad qc copy
               obs_fwd_op_ens_handle%vars(j, k) = missing_r8
            endif
            if (global_ens_index == OBS_GLOBAL_QC_COPY) then
               ! JH set extra fwd_op_ens to bad incoming qc 
               obs_fwd_op_ens_handle%vars(j, k) = global_qc_value
            endif
         enddo

         ! No need to do anything else for a failed observation
         cycle ALL_OBSERVATIONS
      endif

      thiskey(1) = keys(j)

      if(qc_ens_handle%my_num_copies > 0) then
         call get_expected_obs_distrib_state(seq, thiskey, &
            dummy_time, isprior, istatus, &
            assimilate_this_ob, evaluate_this_ob, ens_handle, num_copies_to_calc, my_copy_indices, expected_obs)
            obs_fwd_op_ens_handle%vars(j, 1:num_copies_to_calc) = expected_obs
      else ! need to know whether it was assimilate or evaluate this ob.

         ! TJH istatus is not set in here, yet it is in the other half of the if statement.
         ! TJH Consequently, initializing it after it gets allocated.
         call assim_or_eval(seq, thiskey(1), ens_handle%num_vars, assimilate_this_ob, evaluate_this_ob)

      endif

      do k=1, obs_fwd_op_ens_handle%my_num_copies
         global_ens_index = obs_fwd_op_ens_handle%my_copies(k)

         ! Storing assimilate_this_obs, evaluate_this_ob, or ob not in 
         ! namelist in OBS_EXTRA_QC_COPY
         if (global_ens_index == OBS_EXTRA_QC_COPY) then
            if (assimilate_this_ob) then
               obs_fwd_op_ens_handle%vars(j, k) = DARTQC_ASSIM_GOOD_FOP
            else if (evaluate_this_ob) then
               obs_fwd_op_ens_handle%vars(j, k) = DARTQC_EVAL_GOOD_FOP
            else
               obs_fwd_op_ens_handle%vars(j, k) = DARTQC_NOT_IN_NAMELIST
            endif
         endif
      enddo

      ! TJH if qc_ens_handle%my_num_copies <= 0) istatus was not defined
      qc_ens_handle%vars(j, 1:num_copies_to_calc) = istatus(:)

      call check_forward_operator_istatus(num_copies_to_calc, assimilate_this_ob, &
                                 evaluate_this_ob, istatus, expected_obs, keys(j))

   end do ALL_OBSERVATIONS

else ! distributed state

   do i = 1, num_copies_to_calc
      my_copy_indices(i) = i  ! copy-complete fwd operator so indices are 1 to ens_size
   enddo

   ! Loop through all my observations in the set
   MY_OBSERVATIONS: do j = 1,  obs_fwd_op_ens_handle%my_num_vars

   ! convert the local obs number to global obs number
   global_obs_num = obs_fwd_op_ens_handle%my_vars(j) 
   thiskey(1) = keys(global_obs_num)
  
   call get_obs_from_key(seq, keys(global_obs_num), observation)
   call get_obs_def(observation, obs_def)

   ! Get the information on this observation by placing it in temporary
   ! storage. Check to see if this observation fails input qc test
   call get_qc(observation, input_qc, input_qc_index)

   if (isprior) then
      ! this should only need to be done once, in the prior, right?
      obs_fwd_op_ens_handle%copies(OBS_KEY_COPY, j) = thiskey(1)

      ! Check to see if the incoming data qc is bad
      if(.not. input_qc_ok(input_qc(1), global_qc_value)) then
         obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = global_qc_value
         obs_fwd_op_ens_handle%copies(1:ens_size,j) = missing_r8
         ! No need to do anything else for bad incoming qc
         cycle MY_OBSERVATIONS
      endif
   else ! posterior
      ! Check to see if the incoming data qc is bad
      if(.not. input_qc_ok(input_qc(1), global_qc_value)) then
         cycle MY_OBSERVATIONS
      endif
   endif

   ! Get the observation value and error variance
   call get_obs_values(observation, obs_value(1:1), obs_val_index)

   obs_err_var = get_obs_def_error_variance(obs_def)

   call get_expected_obs_distrib_state(seq, thiskey, &
      dummy_time, isprior, istatus, &
      assimilate_this_ob, evaluate_this_ob, ens_handle, num_copies_to_calc, my_copy_indices, expected_obs)

      obs_fwd_op_ens_handle%copies(1:num_copies_to_calc, j) = expected_obs

   ! collect dart qc
   global_qc_value = nint(obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j))

   call get_dart_qc(istatus, num_copies_to_calc, assimilate_this_ob, evaluate_this_ob, &
                  isprior, global_qc_value)

   ! update the dart qc, error variance and for observed value
   obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = global_qc_value
   obs_fwd_op_ens_handle%copies(OBS_ERR_VAR_COPY  , j) = obs_err_var
   obs_fwd_op_ens_handle%copies(OBS_VAL_COPY      , j) = obs_value(1)

   qc_ens_handle%copies(:, j) = istatus

   call check_forward_operator_istatus(num_copies_to_calc, assimilate_this_ob, evaluate_this_ob, &
                               istatus, expected_obs, thiskey(1))

   end do MY_OBSERVATIONS

endif

! End access to the state
call free_state_window(ens_handle, obs_fwd_op_ens_handle, qc_ens_handle)


! QC Section:
! * Consolidate QC for non-distributed forward operator
! * Check outlier threshold (prior only)
! * Set mean and var to missing_r8 if any forward operator failed. 
!   The failure test is any qc_ens_handle%copies = 0.
!   Previously we were testing for any obs_fwd_op_ens_handle%copies = missing_r8 but
!   some models (e.g. CAM) return a non-zero QC, but a forward operator value (not missing_r8) 

if (get_allow_transpose(ens_handle)) then
   ! Extra step for non-distributed to consolidate qc values for
   ! each ensemble member into global_qc_value
   allocate(var_istatus(qc_ens_handle%num_copies))

   MY_OBS: do j = 1,  obs_fwd_op_ens_handle%my_num_vars
      ! collect dart qc
      var_istatus = qc_ens_handle%copies(:,j) 
   
      assimilate_this_ob = (obs_fwd_op_ens_handle%copies(OBS_EXTRA_QC_COPY, j) == DARTQC_ASSIM_GOOD_FOP ) 
      evaluate_this_ob   = (obs_fwd_op_ens_handle%copies(OBS_EXTRA_QC_COPY, j) == DARTQC_EVAL_GOOD_FOP)
      global_qc_value    = nint(obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j))
   
      ! update the status
      call get_dart_qc(var_istatus, qc_ens_handle%num_copies, assimilate_this_ob, evaluate_this_ob, &
                      isprior, global_qc_value)
      obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = global_qc_value
   
   enddo MY_OBS

   deallocate(var_istatus)

endif

call compute_copy_mean_var(obs_fwd_op_ens_handle, 1, qc_ens_handle%num_copies, &
                           OBS_MEAN_COPY, OBS_VAR_COPY)

!>@todo Science quesion: Should groups be considered separately for outlier threshold test?
QC_LOOP: do j = 1,  obs_fwd_op_ens_handle%my_num_vars

   if (isprior) then

      global_qc_value = nint(obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j))
      ! check_outlier_threshold could change the global_qc_value
      call check_outlier_threshold(obs_fwd_op_ens_handle%copies(OBS_MEAN_COPY,    j),       &
                                   obs_fwd_op_ens_handle%copies(OBS_VAR_COPY,     j),       &
                                   obs_fwd_op_ens_handle%copies(OBS_VAL_COPY,     j),       &
                                   obs_fwd_op_ens_handle%copies(OBS_ERR_VAR_COPY, j),  seq, &
                              nint(obs_fwd_op_ens_handle%copies(OBS_KEY_COPY,     j)), global_qc_value)

      obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = global_qc_value

   endif

   ! for either prior or posterior, if any forward operator failed,
   ! reset the mean/var to missing_r8, regardless of the DART QC status
   if (any(qc_ens_handle%copies(:, j) /= 0)) then
      obs_fwd_op_ens_handle%copies(OBS_MEAN_COPY, j) = missing_r8
      obs_fwd_op_ens_handle%copies(OBS_VAR_COPY,  j) = missing_r8
   endif

end do QC_LOOP

if (isprior) call get_single_copy(obs_fwd_op_ens_handle, OBS_GLOBAL_QC_COPY, prior_qc_copy)

deallocate(expected_obs, istatus, my_copy_indices)

end subroutine get_obs_ens_distrib_state

!------------------------------------------------------------------------------
!> Subroutine get_expected_obs_distrib_state
!>               
!> Compute forward operator for set of obs in sequence for distributed state vector. 
!>@todo does this need to be for a set of obs?
!> 
!------------------------------------------------------------------------------
subroutine get_expected_obs_distrib_state(seq, keys, state_time, isprior, &
   istatus, assimilate_this_ob, evaluate_this_ob, state_ens_handle, num_ens, copy_indices, expected_obs)

type(obs_sequence_type), intent(in)    :: seq  !! the observation sequence
integer,                 intent(in)    :: keys(:)  !! list of obs numbers
type(time_type),         intent(in)    :: state_time  !! time state vector data is valid for
logical,                 intent(in)    :: isprior  !! true if prior; false if posterior
integer,                 intent(in)    :: num_ens  !! number of ensemble members (decl must preceed use)
integer,                 intent(out)   :: istatus(num_ens)  !! FO return codes; 0=ok, >0 is error, <0 reserved for system use
logical,                 intent(out)   :: assimilate_this_ob  !! true if assimilated; false otherwise
logical,                 intent(out)   :: evaluate_this_ob  !! true if evaluated; false otherwise
type(ensemble_type),     intent(in)    :: state_ens_handle  !! state ensemble handle
integer,                 intent(in)    :: copy_indices(num_ens)  !! ??
real(r8),                intent(inout) :: expected_obs(num_ens)  !! the array of computed forward operator values

type(obs_type)     :: obs
type(obs_def_type) :: obs_def
character(len=32)  :: state_size_string, obs_key_string, identity_obs_string

integer :: obs_kind_ind
integer :: num_obs, i

num_obs = size(keys)

! NEED to initialize istatus to okay value
istatus = 0

! Initialize the observation type
!!! Can actually init with the correct size here if wanted
call init_obs(obs, 0, 0)

!>@todo do you ever use this with more than one obs?
do i = 1, num_obs
   call get_obs_from_key(seq, keys(i), obs)
   call get_obs_def(obs, obs_def)

   obs_kind_ind = get_obs_def_type_of_obs(obs_def)

   !location = get_obs_def_location(obs_def)
   
   ! Check in kind for negative for identity obs
   if(obs_kind_ind < 0) then
      if ( -obs_kind_ind > state_ens_handle%num_vars ) then
         write(state_size_string, *) state_ens_handle%num_vars
         write(obs_key_string, *) keys(i)
         write(identity_obs_string, *) -obs_kind_ind
         write(string1,  *) 'unable to compute forward operator for obs number '//trim(adjustl(obs_key_string))
         write(string2, *) 'identity index '//trim(adjustl(identity_obs_string))//&
             ' must be between 1 and the state size of '//trim(adjustl(state_size_string))
         call error_handler(E_ERR, 'get_expected_obs', string1, source, text2=string2)
      endif

      expected_obs =  get_state(-1*int(obs_kind_ind,i8), state_ens_handle)

      ! FIXME : we currently have no option to eval only identity obs,
      ! or select to skip their assimilation via namelist.
      assimilate_this_ob = .true.; evaluate_this_ob = .false.
   
   else ! do forward operator for this kind

      call get_expected_obs_from_def_distrib_state(state_ens_handle, num_ens, copy_indices, keys(i), &
         obs_def, obs_kind_ind, state_time, isprior, &
         assimilate_this_ob, evaluate_this_ob, expected_obs, istatus)

   endif
end do

! need to free any observation specific storage that
! might have been allocated.
call destroy_obs(obs)

end subroutine get_expected_obs_distrib_state

!---------------------------------------------------------------------------
!> For a given obs (seq and key) returns the logicals assimilate_this_obs
!> and evaluate_this_obs.
!> This is so a task that does not do a forward operator can still do QC
!> by storing this info in OBS_EXTRA_QC_COPY.
subroutine assim_or_eval(seq, thiskey, numvars, assimilate_this_ob, evaluate_this_ob)

type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: thiskey
integer(i8),             intent(in)  :: numvars
logical,                 intent(out) :: assimilate_this_ob
logical,                 intent(out) :: evaluate_this_ob

type(obs_type)      :: obs
integer             :: obs_kind_ind
type(obs_def_type)  :: obs_def

call init_obs(obs, 0, 0)
call get_obs_from_key(seq, thiskey, obs)
call get_obs_def(obs, obs_def)
obs_kind_ind = get_obs_def_type_of_obs(obs_def)

if (obs_kind_ind < 0) then
   if ( -obs_kind_ind > numvars ) then
      call error_handler(E_ERR,  'get_expected_obs', &
              'identity obs is outside of state vector ', source)
   endif

   ! FIXME : we currently have no option to eval only identity obs,
   ! or select to skip their assimilation via namelist.
   assimilate_this_ob = .true.; evaluate_this_ob = .false.
else

   assimilate_this_ob = assimilate_this_type_of_obs(obs_kind_ind)
   evaluate_this_ob   = evaluate_this_type_of_obs(obs_kind_ind)

endif

end subroutine assim_or_eval

!------------------------------------------------------------------------------
!> Checks for errors from the model_mod forward operator
!> Possible errors:
!>   * Successful istatus but missing_r8 for forward operator
!>   * Negative istatus
!> This routine calls the error handler (E_ERR) if either of these happen.

subroutine check_forward_operator_istatus(num_fwd_ops, assimilate_ob, evaluate_ob, istatus, expected_obs, thiskey)

integer,  intent(in) :: num_fwd_ops
logical,  intent(in) :: assimilate_ob
logical,  intent(in) :: evaluate_ob
integer,  intent(in) :: istatus(num_fwd_ops)
real(r8), intent(in) :: expected_obs(num_fwd_ops)
integer,  intent(in) :: thiskey

integer :: copy

! Check for errors from model_mod forward operator
do copy = 1, num_fwd_ops

   ! Successful istatus but missing_r8 for forward operator
   if(istatus(copy) == 0) then
      if ((assimilate_ob .or. evaluate_ob) .and. (expected_obs(copy) == missing_r8)) then
         write(string1, *) 'istatus was 0 (OK) but forward operator returned missing value.'
         write(string2, *) 'observation number ', thiskey
         call error_handler(E_ERR,'check_forward_operator_istatus', string1, &
                    source, text2=string2)
      endif
   ! Negative istatus
   else if (istatus(copy) < 0) then
      write(string1, *) 'istatus must not be <0 from forward operator. 0=OK, >0 for error'
      write(string2, *) 'observation number ', thiskey
      call error_handler(E_ERR,'check_forward_operator_istatus', string1, &
                    source, text2=string2)
   endif

enddo

end subroutine check_forward_operator_istatus

!------------------------------------------------------------------------------
end module forward_operator_mod

