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


use types_mod,             only : r8, i8, missing_r8, metadatalength

use time_manager_mod,      only : time_type

use utilities_mod,         only : error_handler, E_ERR, open_file

use mpi_utilities_mod,     only : my_task_id

use obs_sequence_mod,      only : init_obs, destroy_obs, obs_sequence_type, &
                                  obs_type, get_obs_values, get_qc,         &
                                  get_obs_def, get_obs_from_key,            &
                                  get_time_range_keys, replace_obs_values,  &
                                  replace_qc, read_obs_seq_header,          &
                                  get_qc_meta_data, add_qc, get_num_qc,     &
                                  get_copy_meta_data, get_num_copies,       &
                                  read_obs_seq, set_qc_meta_data,           &
                                  set_copy_meta_data

use obs_def_mod,           only : obs_def_type, get_obs_def_error_variance, &
                                  get_expected_obs_from_def_distrib_state,  &
                                  get_obs_def_type_of_obs 

use obs_kind_mod,          only : assimilate_this_type_of_obs, evaluate_this_type_of_obs

use ensemble_manager_mod,  only : ensemble_type, compute_copy_mean_var,         &
                                  init_ensemble_manager, map_pe_to_task,        & 
                                  copies_in_window, get_allow_transpose,        &
                                  all_copies_to_all_vars, allocate_single_copy, &
                                  get_single_copy, get_copy
                                  

use distributed_state_mod, only : create_state_window, free_state_window,   &
                                  get_state

use quality_control_mod,   only : check_outlier_threshold, get_dart_qc, input_qc_ok, &
                                  DARTQC_ASSIM_GOOD_FOP, DARTQC_EVAL_GOOD_FOP,       &
                                  DARTQC_NOT_IN_NAMELIST

!------------------------------------------------------------------------------

implicit none

private

! This type keeps track of meta data associated with the ensemble of forward operators
! Explicitly deals with possible group filter application
type forward_op_info_type
   integer               :: in_obs_copy
   integer               :: obs_val_index
   integer               :: input_qc_index
   integer               :: DART_qc_index
   integer               :: prior_obs_mean_index
   integer               :: posterior_obs_mean_index
   integer               :: prior_obs_spread_index
   integer               :: posterior_obs_spread_index
   real(r8), allocatable :: prior_qc_copy(:)
   integer               :: ERR_VAR_COPY
   integer               :: VAL_COPY
   integer               :: KEY_COPY
   integer               :: GLOBAL_QC_COPY
   integer               :: EXTRA_QC_COPY
   integer               :: MEAN_START
   integer               :: MEAN_END
   integer               :: VAR_START
   integer               :: VAR_END
   integer               :: TOTAL_COPIES
end type forward_op_info_type

public :: get_expected_obs_distrib_state, forward_operators, forward_op_info_type, &
          obs_space_sync_QCs, filter_setup_obs_sequence

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
subroutine get_obs_ens_distrib_state(f, ens_handle, obs_fwd_op_ens_handle, qc_ens_handle, &
   seq, keys, isprior)

type(forward_op_info_type), intent(inout) :: f
type(ensemble_type),        intent(inout) :: ens_handle  !! state ensemble handle
type(ensemble_type),        intent(inout) :: obs_fwd_op_ens_handle  !! observation forward operator handle
type(ensemble_type),        intent(inout) :: qc_ens_handle  !! quality control handle
type(obs_sequence_type),    intent(in)    :: seq  !! the observation sequence
integer,                    intent(in)    :: keys(:)  !! observation key numbers
logical,                    intent(in)    :: isprior  !! true for prior eval; false for posterior

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
      call get_qc(observation, input_qc, f%input_qc_index)

      ! Get the observation value and error variance
      call get_obs_values(observation, obs_value(1:1), f%obs_val_index)
      obs_err_var = get_obs_def_error_variance(obs_def)

      ! Loop through all copies stored by this process and set values as needed
      do k = 1, obs_fwd_op_ens_handle%my_num_copies
         ! See if this is the copy for error variance or observed value
         global_ens_index = obs_fwd_op_ens_handle%my_copies(k)
         if(global_ens_index == f%ERR_VAR_COPY) then
            ! This copy is the instrument observation error variance; read and
            ! store
            obs_fwd_op_ens_handle%vars(j, k) = obs_err_var

         else if(global_ens_index == f%VAL_COPY) then
            ! This copy is the observation from the instrument; read and store
            obs_fwd_op_ens_handle%vars(j, k) = obs_value(1)
         else if(global_ens_index == f%KEY_COPY) then
            obs_fwd_op_ens_handle%vars(j, k) = keys(j)
         else if(global_ens_index == f%GLOBAL_QC_COPY .and. isprior) then
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
            if ((global_ens_index /= f%KEY_COPY) .and. &
               (global_ens_index /= f%EXTRA_QC_COPY) .and. &
               (global_ens_index /= f%GLOBAL_QC_COPY)) then
               ! JH test if the global_ens_index is our extra bad qc copy
               obs_fwd_op_ens_handle%vars(j, k) = missing_r8
            endif
            if (global_ens_index == f%GLOBAL_QC_COPY) then
               ! JH set extra fwd_op_ens to bad incoming qc 
               obs_fwd_op_ens_handle%vars(j, k) = global_qc_value
            endif
         enddo

         ! No need to do anything else for a failed observation
         cycle ALL_OBSERVATIONS
      endif

      thiskey(1) = keys(j)

      if(qc_ens_handle%my_num_copies > 0) then
         call get_expected_obs_distrib_state(seq, thiskey, dummy_time, isprior, istatus, &
            assimilate_this_ob, evaluate_this_ob, ens_handle, num_copies_to_calc,        &
            my_copy_indices, expected_obs)
            obs_fwd_op_ens_handle%vars(j, 1:num_copies_to_calc) = expected_obs
      else ! need to know whether it was assimilate or evaluate this ob.

         ! TJH istatus is not set in here, yet it is in the other half of the if statement.
         ! TJH Consequently, initializing it after it gets allocated.
         call assim_or_eval(seq, thiskey(1), ens_handle%num_vars, &
            assimilate_this_ob, evaluate_this_ob)
      endif

      do k=1, obs_fwd_op_ens_handle%my_num_copies
         global_ens_index = obs_fwd_op_ens_handle%my_copies(k)

         ! Storing assimilate_this_obs, evaluate_this_ob, or ob not in 
         ! namelist in OBS_EXTRA_QC_COPY
         if (global_ens_index == f%EXTRA_QC_COPY) then
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
   call get_qc(observation, input_qc, f%input_qc_index)

   if (isprior) then
      ! this should only need to be done once, in the prior, right?
      obs_fwd_op_ens_handle%copies(f%KEY_COPY, j) = thiskey(1)

      ! Check to see if the incoming data qc is bad
      if(.not. input_qc_ok(input_qc(1), global_qc_value)) then
         obs_fwd_op_ens_handle%copies(f%GLOBAL_QC_COPY, j) = global_qc_value
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
   call get_obs_values(observation, obs_value(1:1), f%obs_val_index)

   obs_err_var = get_obs_def_error_variance(obs_def)

   call get_expected_obs_distrib_state(seq, thiskey, dummy_time, isprior, istatus, &
      assimilate_this_ob, evaluate_this_ob, ens_handle, num_copies_to_calc, &
      my_copy_indices, expected_obs)

      obs_fwd_op_ens_handle%copies(1:num_copies_to_calc, j) = expected_obs

   ! collect dart qc
   global_qc_value = nint(obs_fwd_op_ens_handle%copies(f%GLOBAL_QC_COPY, j))

   call get_dart_qc(istatus, num_copies_to_calc, assimilate_this_ob, evaluate_this_ob, &
                  isprior, global_qc_value)

   ! update the dart qc, error variance and for observed value
   obs_fwd_op_ens_handle%copies(f%GLOBAL_QC_COPY, j) = global_qc_value
   obs_fwd_op_ens_handle%copies(f%ERR_VAR_COPY  , j) = obs_err_var
   obs_fwd_op_ens_handle%copies(f%VAL_COPY      , j) = obs_value(1)

   qc_ens_handle%copies(:, j) = istatus

   call check_forward_operator_istatus(num_copies_to_calc, assimilate_this_ob, &
      evaluate_this_ob, istatus, expected_obs, thiskey(1))

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
   
      assimilate_this_ob = &
         (obs_fwd_op_ens_handle%copies(f%EXTRA_QC_COPY, j) == DARTQC_ASSIM_GOOD_FOP ) 
      evaluate_this_ob   = &
         (obs_fwd_op_ens_handle%copies(f%EXTRA_QC_COPY, j) == DARTQC_EVAL_GOOD_FOP)
      global_qc_value    = nint(obs_fwd_op_ens_handle%copies(f%GLOBAL_QC_COPY, j))
   
      ! update the status
      call get_dart_qc(var_istatus, qc_ens_handle%num_copies, assimilate_this_ob, &
         evaluate_this_ob, isprior, global_qc_value)
      obs_fwd_op_ens_handle%copies(f%GLOBAL_QC_COPY, j) = global_qc_value
   
   enddo MY_OBS

   deallocate(var_istatus)

endif

call compute_copy_mean_var(obs_fwd_op_ens_handle, 1, qc_ens_handle%num_copies, &
                           f%MEAN_START, f%VAR_START)

!>@todo Science quesion: Should groups be considered separately for outlier threshold test?
QC_LOOP: do j = 1,  obs_fwd_op_ens_handle%my_num_vars

   if (isprior) then

      global_qc_value = nint(obs_fwd_op_ens_handle%copies(f%GLOBAL_QC_COPY, j))
      ! check_outlier_threshold could change the global_qc_value
      call check_outlier_threshold(obs_fwd_op_ens_handle%copies(f%MEAN_START, j), &
         obs_fwd_op_ens_handle%copies(f%VAR_START, j),                            &
         obs_fwd_op_ens_handle%copies(f%VAL_COPY, j),                             &
         obs_fwd_op_ens_handle%copies(f%ERR_VAR_COPY, j),  seq,                   &
         nint(obs_fwd_op_ens_handle%copies(f%KEY_COPY, j)), global_qc_value)

      obs_fwd_op_ens_handle%copies(f%GLOBAL_QC_COPY, j) = global_qc_value

   endif

   ! for either prior or posterior, if any forward operator failed,
   ! reset the mean/var to missing_r8, regardless of the DART QC status
   if (any(qc_ens_handle%copies(:, j) /= 0)) then
      obs_fwd_op_ens_handle%copies(f%MEAN_START, j) = missing_r8
      obs_fwd_op_ens_handle%copies(f%VAR_START,  j) = missing_r8
   endif

end do QC_LOOP

if (isprior) call get_single_copy(obs_fwd_op_ens_handle, f%GLOBAL_QC_COPY, f%prior_qc_copy)

deallocate(expected_obs, istatus, my_copy_indices)

end subroutine get_obs_ens_distrib_state

!------------------------------------------------------------------------------
!> Subroutine get_expected_obs_distrib_state
!>               
!> Compute forward operator for set of obs in sequence for distributed state vector. 
!>@todo does this need to be for a set of obs?
!> 
!------------------------------------------------------------------------------
subroutine get_expected_obs_distrib_state(seq, keys, state_time, isprior,    &
   istatus, assimilate_this_ob, evaluate_this_ob, state_ens_handle,          &
   num_ens, copy_indices, expected_obs)

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

subroutine check_forward_operator_istatus(num_fwd_ops, assimilate_ob, evaluate_ob, &
   istatus, expected_obs, thiskey)

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

subroutine forward_operators(f, state_ens_handle, obs_fwd_op_ens_handle, qc_ens_handle, &
   seq, ens_size, num_groups, num_obs_in_set, keys, key_bounds, num_output_obs_members, &
   compute_posterior, isprior)

type(forward_op_info_type), intent(inout) :: f
type(ensemble_type),        intent(inout) :: state_ens_handle, obs_fwd_op_ens_handle, qc_ens_handle
type(obs_sequence_type),    intent(inout) :: seq
integer,                    intent(in)    :: ens_size, num_groups, num_obs_in_set
integer, allocatable,       intent(inout) :: keys(:)
integer,                    intent(in)    :: key_bounds(2)
integer,                    intent(in)    :: num_output_obs_members
logical,                    intent(in)    :: compute_posterior
logical,                    intent(in)    :: isprior

! Compute forward operators and handle data structures for the observations and the qc
integer :: obs_mean_index, obs_spread_index, obs_copy_offset

! Initialization required only for prior observation operator computation
if(isprior) then

   f%ERR_VAR_COPY   = ens_size + 1
   f%VAL_COPY       = ens_size + 2
   f%KEY_COPY       = ens_size + 3
   f%GLOBAL_QC_COPY = ens_size + 4
   f%EXTRA_QC_COPY  = ens_size + 5
   f%MEAN_START     = ens_size + 6
   f%MEAN_END       = f%MEAN_START + num_groups -1 
   f%VAR_START      = f%MEAN_END + 1
   f%VAR_END        = f%VAR_START + num_groups - 1
   f%TOTAL_COPIES   = ens_size + 5 + 2*num_groups

   call init_ensemble_manager(obs_fwd_op_ens_handle, f%TOTAL_COPIES, &
                              int(num_obs_in_set,i8), 1, transpose_type_in = 2)
   
   ! Also need a qc field for copy of each observation
   call init_ensemble_manager(qc_ens_handle, ens_size, &
                              int(num_obs_in_set,i8), 1, transpose_type_in = 2)

   ! Allocate storage for the keys for this number of observations
   allocate(keys(num_obs_in_set)) ! This is still var size for writing out the observation sequence

   ! Get all the keys associated with this set of observations
   ! Is there a way to distribute this?
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   ! allocate() space for the prior qc copy
   call allocate_single_copy(obs_fwd_op_ens_handle, f%prior_qc_copy)

endif 

call get_obs_ens_distrib_state(f, state_ens_handle, obs_fwd_op_ens_handle, qc_ens_handle, &
   seq, keys, isprior)
      
! This is where the mean obs ! copy ( + others ) is moved to task 0 so task 0 can update seq.
! There is a transpose (all_copies_to_all_vars(obs_fwd_op_ens_handle)) in obs_space_diagnostics
if(isprior) then
   obs_mean_index = f%prior_obs_mean_index
   obs_spread_index = f%prior_obs_spread_index
   obs_copy_offset = f%in_obs_copy + 1
else
   obs_mean_index = f%posterior_obs_mean_index
   obs_spread_index = f%posterior_obs_spread_index
   obs_copy_offset = f%in_obs_copy + 2
endif

! Do prior observation space diagnostics and associated quality control
call obs_space_diagnostics(f, obs_fwd_op_ens_handle, qc_ens_handle, ens_size, seq, keys, &
   isprior, num_output_obs_members, obs_copy_offset, &
   obs_mean_index, obs_spread_index, num_obs_in_set, compute_posterior)

! Free up the prior_qc_copy if this is being called fopr posterior (all done with it)
if(.not. isprior) deallocate(f%prior_qc_copy)

end subroutine forward_operators

!------------------------------------------------------------------------------
subroutine obs_space_diagnostics(f, obs_fwd_op_ens_handle, qc_ens_handle, ens_size, &
   seq, keys, isprior, num_output_members, members_index, &
   ens_mean_index, ens_spread_index, num_obs_in_set, do_post)

! Do observation space diagnostics on the set of obs corresponding to keys

type(forward_op_info_type), intent(inout) :: f
type(ensemble_type),        intent(inout) :: obs_fwd_op_ens_handle, qc_ens_handle
integer,                    intent(in)    :: ens_size
integer,                    intent(in)    :: num_obs_in_set
logical,                    intent(in)    :: isprior
integer,                    intent(in)    :: keys(num_obs_in_set)
integer,                    intent(in)    :: num_output_members, members_index
integer,                    intent(in)    :: ens_mean_index, ens_spread_index
type(obs_sequence_type),    intent(inout) :: seq
logical,                    intent(in)    :: do_post

integer               :: j, k, ens_offset, copy_factor
integer               :: ivalue, io_task, my_task
real(r8), allocatable :: obs_temp(:)
real(r8)              :: rvalue(1)

! Do verbose forward operator output if requested
!!!JLAif(output_forward_op_errors) call verbose_forward_op_output(qc_ens_handle, isprior, ens_size, keys)

! this is a query routine to return which task has 
! logical processing element 0 in this ensemble.
io_task = map_pe_to_task(obs_fwd_op_ens_handle, 0)
my_task = my_task_id()

! single value per member if no posterior, else 2
if (do_post) then
   copy_factor = 2
else
   copy_factor = 1
endif

! Make var complete for get_copy() calls below.
! Optimize: Could we use a gather instead of a transpose and get copy?
call all_copies_to_all_vars(obs_fwd_op_ens_handle)

! allocate temp space for sending data only on the task that will
! write the obs_seq.final file
if (my_task == io_task) then
   allocate(obs_temp(num_obs_in_set))
else ! TJH: this change became necessary when using Intel 19.0.5 ...
   allocate(obs_temp(1))
endif


! Update the ensemble mean
call get_copy(io_task, obs_fwd_op_ens_handle, f%MEAN_START, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_obs_values(seq, keys(j), rvalue, ens_mean_index)
     end do
  endif

! Update the ensemble spread
call get_copy(io_task, obs_fwd_op_ens_handle, f%VAR_START, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      if (obs_temp(j) /= missing_r8) then
         rvalue(1) = sqrt(obs_temp(j))
      else
         rvalue(1) = obs_temp(j)
      endif
      call replace_obs_values(seq, keys(j), rvalue, ens_spread_index)
   end do
endif

! Update any requested ensemble members
ens_offset = members_index + 2*copy_factor
do k = 1, num_output_members
   call get_copy(io_task, obs_fwd_op_ens_handle, k, obs_temp)
   if(my_task == io_task) then
      ivalue = ens_offset + copy_factor * (k - 1)
      do j = 1, obs_fwd_op_ens_handle%num_vars
         rvalue(1) = obs_temp(j)
         call replace_obs_values(seq, keys(j), rvalue, ivalue)
      end do
   endif
end do

! Update the qc global value
call get_copy(io_task, obs_fwd_op_ens_handle, f%GLOBAL_QC_COPY, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_qc(seq, keys(j), rvalue, f%DART_qc_index)
   end do
endif

deallocate(obs_temp)

end subroutine obs_space_diagnostics

   
!-------------------------------------------------------------------------
!> write out failed forward operators
!> This was part of obs_space_diagnostics
                                  
subroutine verbose_forward_op_output(qc_ens_handle, isprior, ens_size, keys)

type(ensemble_type), intent(inout) :: qc_ens_handle
logical,             intent(in)    :: isprior
integer,             intent(in)    :: ens_size
integer,             intent(in)    :: keys(:) ! I think this is still var size
   
character(len=12) :: task  
integer :: j, i                   
integer :: forward_unit

write(task, '(i6.6)') my_task_id()
   
! all tasks open file?
if(isprior) then 
   forward_unit = open_file('prior_forward_ope_errors' // task, 'formatted', 'append')        
else                              
   forward_unit = open_file('post_forward_ope_errors' // task, 'formatted', 'append')         
endif                             
                                  
! qc_ens_handle is a real representing an integer; values /= 0 get written out                
do i = 1, ens_size                
   do j = 1, qc_ens_handle%my_num_vars
      if(nint(qc_ens_handle%copies(i, j)) /= 0) write(forward_unit, *) i, keys(j), nint(qc_ens_handle%copies(i, j))
   end do
end do                            
                                  
call close_file(forward_unit)     
                                  
end subroutine verbose_forward_op_output
   
!-------------------------------------------------------------------------
      
subroutine obs_space_sync_QCs(f, obs_fwd_op_ens_handle, seq, keys, num_obs_in_set)
   
type(forward_op_info_type), intent(inout) :: f
type(ensemble_type),        intent(inout) :: obs_fwd_op_ens_handle
integer,                    intent(in)    :: num_obs_in_set
integer,                    intent(in)    :: keys(num_obs_in_set)
type(obs_sequence_type),    intent(inout) :: seq
      
integer               :: j
integer               :: io_task, my_task 
real(r8), allocatable :: obs_temp(:)
real(r8)              :: rvalue(1)
   
! this is a query routine to return which task has 
! logical processing element 0 in this ensemble.
io_task = map_pe_to_task(obs_fwd_op_ens_handle, 0)
my_task = my_task_id()            
   
! create temp space for QC values
if (my_task == io_task) then
   allocate(obs_temp(num_obs_in_set))
else 
   allocate(obs_temp(1))
endif
                                  
! Optimize: Could we use a gather instead of a transpose and get copy?
call all_copies_to_all_vars(obs_fwd_op_ens_handle)

! Update the qc global value
call get_copy(io_task, obs_fwd_op_ens_handle, f%GLOBAL_QC_COPY, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_qc(seq, keys(j), rvalue, f%DART_qc_index)
   end do
endif

deallocate(obs_temp)

end subroutine obs_space_sync_QCs

!-------------------------------------------------------------------------

subroutine filter_setup_obs_sequence(f, seq, num_output_obs_members, &
   obs_sequence_in_name, do_post)

type(forward_op_info_type), intent(inout) :: f
type(obs_sequence_type),    intent(inout) :: seq
integer,                    intent(in)    :: num_output_obs_members
character(len = *),         intent(in)    :: obs_sequence_in_name
logical,                    intent(in)    :: do_post
   
character(len=metadatalength) :: no_qc_meta_data = 'No incoming data QC'
character(len=metadatalength) :: dqc_meta_data   = 'DART quality control'
character(len=129) :: obs_seq_read_format
integer :: obs_seq_file_id, copies_num_inc, qc_num_inc
integer :: tnum_copies, tnum_qc, tnum_obs, tmax_num_obs
integer :: my_task, io_task
logical :: pre_I_format
   
! Input file can have one qc field, none, or more.  note that read_obs_seq_header
! does NOT return the actual metadata values, which would be helpful in trying
! to decide if we need to add copies or qcs.
call read_obs_seq_header(obs_sequence_in_name, tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, &
   obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)
   
! return the original number of copies in the obs_seq file
! before we add any copies for diagnostics.
f%in_obs_copy = tnum_copies

! FIXME: this should be called from inside obs_space_diagnostics the first
! time that routine is called, so it has an ensemble handle to query for 
! exactly which task is pe0 (or use a different pe number).  here we 
! have to assume task 0 == pe 0 which is currently true but someday
! we would like to be able to change.
io_task = 0
my_task = my_task_id()

! only the task writing the obs_seq.final file needs space for the
! additional copies/qcs.  for large numbers of individual members
! in the final file this takes quite a bit of memory. 

if (my_task == io_task) then
   ! Determine the number of output obs space fields
   if (do_post) then
      ! 4 is for prior/posterior mean and spread, plus
      ! prior/posterior values for all requested members
      copies_num_inc = 4 + (2 * num_output_obs_members)
   else
      ! 2 is for prior mean and spread, plus
      ! prior values for all requested members
      copies_num_inc = 2 + (1 * num_output_obs_members)
   endif
else
   copies_num_inc = 0
endif

! if there are less than 2 incoming qc fields, we will need
! to make at least 2 (one for the dummy data qc and one for
! the dart qc) on task 0.  other tasks just need 1 for incoming qc.
if (tnum_qc < 2) then
   if (my_task == io_task) then
      qc_num_inc = 2 - tnum_qc
   else
      qc_num_inc = 1 - tnum_qc
   endif
else
   qc_num_inc = 0
endif

! Read in with enough space for diagnostic output values and add'l qc field(s)
! ONLY ADD SPACE ON TASK 0.  everyone else just read in the original obs_seq file.
call read_obs_seq(obs_sequence_in_name, copies_num_inc, qc_num_inc, 0, seq)

! check to be sure that we have an incoming qc field.  if not, look for
! a blank qc field
f%input_qc_index = get_obs_qc_index(seq)
if (f%input_qc_index < 0) then
   f%input_qc_index = get_blank_qc_index(seq)
   if (f%input_qc_index < 0) then
      ! Need 1 new qc field for dummy incoming qc
      call add_qc(seq, 1)
      f%input_qc_index = get_blank_qc_index(seq)
      if (f%input_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen', source)
      endif
   endif
   ! Since we are constructing a dummy QC, label it as such
   call set_qc_meta_data(seq, f%input_qc_index, no_qc_meta_data)
endif

! check to be sure we either find an existing dart qc field and
! reuse it, or we add a new one. only on task 0.
f%DART_qc_index = get_obs_dartqc_index(seq)
if (f%DART_qc_index < 0 .and. my_task == io_task) then
   f%DART_qc_index = get_blank_qc_index(seq)
   if (f%DART_qc_index < 0) then
      ! Need 1 new qc field for the DART quality control
      call add_qc(seq, 1)
      f%DART_qc_index = get_blank_qc_index(seq)
      if (f%DART_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen', source)
      endif
   endif
   call set_qc_meta_data(seq, f%DART_qc_index, dqc_meta_data)
endif

! Determine which copy has actual obs value and return it.
f%obs_val_index = get_obs_copy_index(seq)

! Initialize the output sequences and state files and set their meta data
call filter_generate_copy_meta_data(f, seq, num_output_obs_members, do_post) 

end subroutine filter_setup_obs_sequence

!------------------------------------------------------------------------------
                                  
function get_obs_qc_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_qc_index
                                  
integer :: i                      

! Determine which qc, if any, has the incoming obs qc
! this is tricky because we have never specified what string
! the metadata has to have.  look for 'qc' or 'QC' and the
! first metadata that matches (much like 'observation' above)
! is the winner.
   
do i = 1, get_num_qc(seq)  
   get_obs_qc_index = i           
                                  
   ! Need to avoid 'QC metadata not initialized'
   if(index(get_qc_meta_data(seq, i), 'QC metadata not initialized') > 0) cycle               
                                  
   ! Need to look for 'QC' or 'qc'
   if(index(get_qc_meta_data(seq, i), 'QC') > 0) return                                       
   if(index(get_qc_meta_data(seq, i), 'qc') > 0) return 
   if(index(get_qc_meta_data(seq, i), 'Quality Control') > 0) return
   if(index(get_qc_meta_data(seq, i), 'QUALITY CONTROL') > 0) return
end do
! Falling off end means 'QC' string not found; not fatal!
                                  
get_obs_qc_index = -1             
                                  
end function get_obs_qc_index     

!-------------------------------------------------------------------------
   
function get_obs_dartqc_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_dartqc_index

integer :: i
                                  
! Determine which qc, if any, has the DART qc
                                  
do i = 1, get_num_qc(seq)         
   get_obs_dartqc_index = i       
   ! Need to look for 'DART quality control'
   if(index(get_qc_meta_data(seq, i), 'DART quality control') > 0) return                     
end do                            
! Falling off end means 'DART quality control' not found; not fatal!
                                 
get_obs_dartqc_index = -1  

end function get_obs_dartqc_index 

!-------------------------------------------------------------------------
                                  
function get_blank_qc_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_blank_qc_index
                                  
integer :: i                      

! Determine which qc, if any, is blank
                                  
do i = 1, get_num_qc(seq)
   get_blank_qc_index = i
   ! Need to look for 'QC metadata not initialized'
   if(index(get_qc_meta_data(seq, i), 'QC metadata not initialized') > 0) return
end do
! Falling off end means unused slot not found; not fatal!

get_blank_qc_index = -1

end function get_blank_qc_index

!-------------------------------------------------------------------------

function get_obs_copy_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_copy_index

integer :: i 

! Determine which copy in sequence has actual obs

do i = 1, get_num_copies(seq)
   get_obs_copy_index = i
   ! Need to look for 'observation'
   if(index(get_copy_meta_data(seq, i), 'observation') > 0) return
end do
! Falling of end means 'observations' not found; die
call error_handler(E_ERR,'get_obs_copy_index', &
   'Did not find observation copy with metadata "observation"', source)

end function get_obs_copy_index

!-------------------------------------------------------------------------

subroutine filter_generate_copy_meta_data(f, seq, num_output_obs_members, &
   do_post)

type(forward_op_info_type),  intent(inout) :: f
type(obs_sequence_type),     intent(inout) :: seq
integer,                     intent(in)    :: num_output_obs_members
logical,                     intent(in)    :: do_post

! Figures out the strings describing the output copies for the observation sequence.

character(len=metadatalength) :: prior_meta_data, posterior_meta_data
integer :: i, num_obs_copies

! only PE0 (here task 0) will allocate space for the obs_seq.final
!
! all other tasks should NOT allocate all this space.
! instead, set the copy numbers to an illegal value
! so we'll trap if they're used, and return early.
if (my_task_id() /= 0) then
   f%prior_obs_mean_index  = -1
   f%posterior_obs_mean_index = -1
   f%prior_obs_spread_index  = -1
   f%posterior_obs_spread_index = -1
   return
endif

! Set the metadata for the observations.

! Set up obs ensemble mean
num_obs_copies = f%in_obs_copy

num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
f%prior_obs_mean_index = num_obs_copies

if (do_post) then
   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble mean'
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   f%posterior_obs_mean_index = num_obs_copies
endif


! Set up obs ensemble spread
num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble spread'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
f%prior_obs_spread_index = num_obs_copies

if (do_post) then
   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble spread'
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   f%posterior_obs_spread_index = num_obs_copies
endif

! Make sure there are not too many copies requested - 
! proposed: make this magic number set in 1 place with an accessor
! routine so all parts of the code agree on max values.
if(num_output_obs_members > 10000) then
   write(string1, *)'output metadata in filter needs obs ensemble size < 10000, not ',&
                      num_output_obs_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',string1,source)
endif

! Set up obs ensemble members as requested
do i = 1, num_output_obs_members
   num_obs_copies = num_obs_copies + 1
   write(prior_meta_data, '(a21, 1x, i6)') 'prior ensemble member', i
   call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   if (do_post) then
      num_obs_copies = num_obs_copies + 1
      write(posterior_meta_data, '(a25, 1x, i6)') 'posterior ensemble member', i
      call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   endif
end do


end subroutine filter_generate_copy_meta_data

!-------------------------------------------------------------------------

end module forward_operator_mod
