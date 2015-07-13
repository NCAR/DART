! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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
                                  get_obs_kind 

use ensemble_manager_mod,  only : ensemble_type, compute_copy_mean_var, &
                                  prepare_to_read_from_vars,            &
                                  prepare_to_write_to_vars,             &
                                  get_my_num_copies

use distributed_state_mod, only : create_state_window, free_state_window,   &
                                  get_state

use data_structure_mod,    only : copies_in_window

use quality_control_mod,   only : check_outlier_threshold, get_dart_qc, &
                                  input_qc_ok, good_dart_qc

!------------------------------------------------------------------------------

implicit none

private

public :: get_obs_ens_distrib_state, get_expected_obs_distrib_state

!------------------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!> Subroutine get_obs_ens_distrib_state
!> 
!> Computes the forward observation operators and related quality 
!> control indicators.
!> 
!> @param[in]    ens_handle - state ensemble handle
!> @param[inout] obs_fwd_op_ens_handle - observation forward operator handle
!> @param[inout] qc_ens_handle - quality control handle
!> @param[in]    seq - the observation sequence
!> @param[in]    keys - unique integer key when is an obs_sequence
!> @param[in]    obs_val_index - observation value index
!> @param[in]    input_qc_index - input QC index
!> @param[in]    OBS_ERR_VAR_COPY - observation error variance copy number
!> @param[in]    OBS_VAL_COPY - observation value copy number
!> @param[in]    OBS_KEY_COPY - observation key copy number
!> @param[in]    OBS_GLOBAL_QC_COPY - observation global QC copy number
!> @param[in]    OBS_MEAN_START - start of group observation mean copy number
!> @param[in]    OBS_VAR_START - start of group observation variance copy number
!> @param[in]    isprior - true for prior eval; false for posterior
!------------------------------------------------------------------------------
subroutine get_obs_ens_distrib_state(ens_handle, obs_fwd_op_ens_handle, &
   qc_ens_handle, seq, keys, obs_val_index, input_qc_index, &
   OBS_ERR_VAR_COPY,   OBS_VAL_COPY,   OBS_KEY_COPY, &
   OBS_GLOBAL_QC_COPY, OBS_MEAN_START, OBS_VAR_START, isprior)

type(ensemble_type),     intent(in)    :: ens_handle 
type(ensemble_type),     intent(inout) :: obs_fwd_op_ens_handle
type(ensemble_type),     intent(inout) :: qc_ens_handle 
type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: keys(:)
integer,                 intent(in)    :: obs_val_index
integer,                 intent(in)    :: input_qc_index
integer,                 intent(in)    :: OBS_ERR_VAR_COPY
integer,                 intent(in)    :: OBS_VAL_COPY
integer,                 intent(in)    :: OBS_KEY_COPY
integer,                 intent(in)    :: OBS_GLOBAL_QC_COPY
integer,                 intent(in)    :: OBS_MEAN_START
integer,                 intent(in)    :: OBS_VAR_START 
logical,                 intent(in)    :: isprior

real(r8) :: input_qc(1), obs_value(1), obs_err_var, thisvar(1)
real(r8) :: error, diff_sd, ratio
real(r8) :: obs_prior_mean, obs_prior_var, obs_val

real(r8), allocatable :: expected_obs(:) !Also regular obs now?

integer :: j, k !< index variables 
integer :: thiskey(1)
integer :: global_obs_num, global_qc_value
integer :: forward_min, forward_max !< for global qc
integer :: ens_size, my_num_obs

integer, allocatable  :: istatus(:)

logical :: evaluate_this_ob, assimilate_this_ob
logical :: failed

type(time_type) :: dummy_time

type(obs_def_type) :: obs_def
type(obs_type)     :: observation

! IMPORTANT, IT IS ASSUMED THAT ACTUAL ENSEMBLES COME FIRST
! HK: I think it is also assumed that the ensemble members are in the same 
! order in each of the handles

ens_size = copies_in_window(ens_handle)
my_num_obs = obs_fwd_op_ens_handle%my_num_vars

! make some room for state vectors
allocate(istatus(ens_size)) 
allocate(expected_obs(ens_size))

! FIXME: these no longer do anything?
! call prepare_to_write_to_vars(obs_fwd_op_ens_handle)
! call prepare_to_write_to_vars(qc_ens_handle)
! call prepare_to_read_from_vars(ens_handle)

! create the mpi window for the distributed state
call create_state_window(ens_handle)

! Loop through all observations in the set
ALL_OBSERVATIONS: do j = 1, my_num_obs

   ! convert the local obs number to global obs number
   global_obs_num = obs_fwd_op_ens_handle%my_vars(j) 
   thiskey(1) = keys(global_obs_num)
  
   call get_obs_from_key(seq, keys(global_obs_num), observation)
   call get_obs_def(observation, obs_def)

   if (isprior) then
      ! this should only need to be done once, in the prior, right?
      obs_fwd_op_ens_handle%copies(OBS_KEY_COPY, j) = thiskey(1)

      ! Get the information on this observation by placing it in temporary
      ! storage. Check to see if this observation fails input qc test
      call get_qc(observation, input_qc, input_qc_index)

      ! Check to see if the forward operater failed
      if(.not. input_qc_ok(input_qc(1), global_qc_value)) then
         obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = global_qc_value
         ! No need to do anything else for a failed observation
         cycle ALL_OBSERVATIONS
      endif
   else ! posterior
      global_qc_value = nint(obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) )
      if (.not. good_dart_qc(global_qc_value)) then
         cycle ALL_OBSERVATIONS ! prior forward op failed
      endif
   endif

   ! Get the observation value and error variance
   call get_obs_values(observation, obs_value(1:1), obs_val_index)

   obs_err_var = get_obs_def_error_variance(obs_def)

   ! temporaries to avoid passing array sections which was slow on PGI compiler
   call get_expected_obs_distrib_state(seq, thiskey, &
      dummy_time, isprior, istatus, &
      assimilate_this_ob, evaluate_this_ob, ens_handle, expected_obs)

   obs_fwd_op_ens_handle%copies(1:ens_size, j) = expected_obs

   ! collect dart qc
   global_qc_value = nint(obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j))

   call get_dart_qc(istatus, ens_size, assimilate_this_ob, evaluate_this_ob, &
                     isprior, global_qc_value)

   ! update the dart qc, error variance and for observed value
   obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = global_qc_value
   obs_fwd_op_ens_handle%copies(OBS_ERR_VAR_COPY  , j) = obs_err_var
   obs_fwd_op_ens_handle%copies(OBS_VAL_COPY      , j) = obs_value(1)

   ! do we need this?  nsc: i think no 
   ! qc_ens_handle%copies(:, j) = istatus

   ! redundant, right?
   !obs_fwd_op_ens_handle%copies(OBS_KEY_COPY, j) = thiskey(1)

end do ALL_OBSERVATIONS

!> @todo - don't you have the mean already?
call compute_copy_mean_var(obs_fwd_op_ens_handle, 1, ens_size, &
        OBS_MEAN_START, OBS_VAR_START)

!>@todo Science quesion: Should groups be considered separately for outlier
! we cannot get rid of two loops because we need obs mean and var to do the
! outlier check.  so we'll still need the obs_fwd_op handle to store the istatus,
! or we need a single column compute_mean_var to put inside the first loop.
QC_LOOP: do j = 1, my_num_obs

   global_qc_value = nint(obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j))

   ! If istatus is 0 (successful) then put 0 for assimilate, -1 for evaluate only
   ! and -2 for neither evaluate or assimilate. Otherwise set the istatus
   ! in the forward operator evaluation field
   if (isprior) then
      call check_outlier_threshold(obs_fwd_op_ens_handle%copies(OBS_MEAN_START,j), &
                                   obs_fwd_op_ens_handle%copies(OBS_VAR_START    , j),       &
                                   obs_fwd_op_ens_handle%copies(OBS_VAL_COPY     , j),       &
                                   obs_fwd_op_ens_handle%copies(OBS_ERR_VAR_COPY , j),  seq, &
                              nint(obs_fwd_op_ens_handle%copies(OBS_KEY_COPY     , j)), global_qc_value)
   endif

   obs_fwd_op_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = global_qc_value

   ! for either prior or posterior, if the forward operator failed,
   ! reset the mean/var to missing_r8, regardless of the DART QC status
   ! HK does this fail if you have groups?
   !>@todo Do we want to set all the groups to missing_r8? Not just the start?
   if (.not. good_dart_qc(global_qc_value)) then
      obs_fwd_op_ens_handle%copies(OBS_MEAN_START, j) = missing_r8
      obs_fwd_op_ens_handle%copies(OBS_VAR_START,  j) = missing_r8
   endif

end do QC_LOOP

call free_state_window

deallocate(expected_obs)
deallocate(istatus)
  
end subroutine get_obs_ens_distrib_state

!------------------------------------------------------------------------------
!> Subroutine get_expected_obs_distrib_state
!>               
!> Compute forward operator for set of obs in sequence for distributed state vector. 
!> @todo does this need to be for a set of obs?
!> 
!> @param[in]    seq - the observation sequence 
!> @param[in]    keys - unique integer key into an obs_sequence
!> @param[in]    state_time -  the time of the state vector data
!> @param[in]    isprior - true for prior eval; false for posterior
!> @param[out]   istatus - return code: 0=ok, >0 is error, <0 reserved for system use
!> @param[out]   assimilate_this_ob - true if obs was assimilated, false otherwise 
!> @param[out]   evaluate_this_ob - true if obs was evaluated, false otherwise
!> @param[in]    state_ens_handle - state ensemble handle
!> @param[inout] expected_obs - the computed forward operator value
!------------------------------------------------------------------------------
subroutine get_expected_obs_distrib_state(seq, keys, state_time, isprior, &
   istatus, assimilate_this_ob, evaluate_this_ob, state_ens_handle, expected_obs)

type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: keys(:)
! integer,                 intent(in)    :: ens_index
type(time_type),         intent(in)    :: state_time
logical,                 intent(in)    :: isprior
integer,                 intent(out)   :: istatus(:)
logical,                 intent(out)   :: assimilate_this_ob
logical,                 intent(out)   :: evaluate_this_ob
!HK
type(ensemble_type),     intent(in)    :: state_ens_handle
real(r8), dimension(:),  intent(inout) :: expected_obs !> @todo needs to be 2d for a set of obs

type(obs_type)     :: obs
type(obs_def_type) :: obs_def

integer :: obs_kind_ind
integer :: num_obs, i
integer :: length_of_expected_obs ! HK should this be passed in?


num_obs = size(keys)
length_of_expected_obs = copies_in_window(state_ens_handle)

! NEED to initialize istatus to okay value
istatus = 0

! Initialize the observation type
!!! Can actually init with the correct size here if wanted
call init_obs(obs, 0, 0)
do i = 1, num_obs !> @todo do you ever use this with more than one obs?
   call get_obs_from_key(seq, keys(i), obs)
   call get_obs_def(obs, obs_def)

   obs_kind_ind = get_obs_kind(obs_def)

   !location = get_obs_def_location(obs_def)
   
   ! Check in kind for negative for identity obs
   if(obs_kind_ind < 0) then
      if ( -obs_kind_ind > state_ens_handle%num_vars ) call error_handler(E_ERR, &
         'get_expected_obs', &
         'identity obs is outside of state vector ', &
         source, revision, revdate)

      call get_state(expected_obs, -1*int(obs_kind_ind,i8), state_ens_handle)

      ! FIXME : we currently have no option to eval only identity obs,
      ! or select to skip their assimilation via namelist.
      assimilate_this_ob = .true.; evaluate_this_ob = .false.
   
   else ! do forward operator for this kind
      !> Q. Do we loop around copies here? This would mean ens_size*times the comumication
      !> The alternative is to alter the code in model_mod.f90 to work on arrays of ensemble size.
      !> Currently looping in model_mod.f90 for lorenz_96 

      call get_expected_obs_from_def_distrib_state(keys(i), obs_def, obs_kind_ind, &
         state_time, isprior, istatus, &
         assimilate_this_ob, evaluate_this_ob, expected_obs, state_ens_handle)

   endif
end do

! need to free any observation specific storage that
! might have been allocated.
call destroy_obs(obs)

end subroutine get_expected_obs_distrib_state

!------------------------------------------------------------------------------
end module forward_operator_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$


