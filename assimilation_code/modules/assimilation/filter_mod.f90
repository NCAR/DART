! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module filter_mod

!------------------------------------------------------------------------------
use types_mod,             only : r8, i8, missing_r8, metadatalength, MAX_NUM_DOMS, MAX_FILES

use options_mod,           only : get_missing_ok_status, set_missing_ok_status

use obs_sequence_mod,      only : obs_type, obs_sequence_type, write_obs_seq, & 
                                  static_init_obs_sequence, delete_seq_head,  &        
                                  delete_seq_tail, destroy_obs_sequence
                                 
use obs_def_utilities_mod, only : set_debug_fwd_op

use time_manager_mod,      only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                  operator(-), print_time

use utilities_mod,         only : error_handler, E_ERR, E_MSG,                                &
                                  logfileunit, nmlfileunit,                                   &
                                  find_namelist_in_file, check_namelist_read,      &
                                  do_nml_file, do_nml_term

use assim_model_mod,       only : static_init_assim_model, get_model_size,                    &
                                  end_assim_model, get_state_meta_data

use assim_tools_mod,       only : filter_assim, set_assim_tools_trace
use obs_model_mod,         only : move_ahead, advance_state, set_obs_model_trace

use ensemble_manager_mod,  only : init_ensemble_manager, end_ensemble_manager,                &
                                  ensemble_type, get_copy, get_my_num_copies,                 &
                                  all_vars_to_all_copies, all_copies_to_all_vars,             &
                                  compute_copy_mean, compute_copy_mean_sd,                    &
                                  get_copy_owner_index,                                       &
                                  get_ensemble_time, set_ensemble_time,                       &
                                  map_pe_to_task,                                             &
                                  set_num_extra_copies,                                       &
                                  allocate_vars
                                  
use adaptive_inflate_mod,  only : do_ss_inflate, &
                                  inflate_ens, adaptive_inflate_init,                 &
                                  adaptive_inflate_type, log_inflation_info,          &
                                  do_rtps_inflate, set_inflate_flavor,                &
                                  NO_INFLATION

use mpi_utilities_mod,     only : my_task_id, task_sync, broadcast_send, broadcast_recv,      &
                                  task_count, iam_task0

use state_vector_io_mod,   only : state_vector_io_init, write_state

use io_filenames_mod,      only : file_info_type

use direct_netcdf_mod,     only : finalize_single_file_io

use state_structure_mod,   only : get_num_domains

use forward_operator_mod,  only : forward_operators, forward_op_info_type, obs_space_sync_QCs, &
                                  filter_setup_obs_sequence

use quality_control_mod,   only : initialize_qc

use location_mod,          only : location_type

use probit_transform_mod,  only : transform_to_probit, transform_from_probit

use algorithm_info_mod,    only : probit_dist_info, init_algorithm_info_mod, end_algorithm_info_mod

use distribution_params_mod, only : distribution_params_type

use filter_io_diag_mod,    only : create_ensemble_from_single_file, &
                                  count_state_ens_copies, read_state_and_inflation,      &
                                  output_diagnostics, init_output_file_info,        &
                                  ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY,        &
                                  PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY,&
                                  RTPS_PRIOR_SPREAD_COPY

!------------------------------------------------------------------------------

implicit none
private

public :: filter_sync_keys_time, &
          filter_set_initial_time, &
          filter_main

character(len=*), parameter :: source = 'filter_mod.f90'

! Some convenient global storage items
character(len=512)      :: msgstring

integer :: trace_level, timestamp_level

! Determine if inflation is turned on or off for reading and writing
! inflation restart files
logical :: output_inflation = .false.

! Module Global Variables for inflation
type(adaptive_inflate_type) :: prior_inflate, post_inflate

logical :: has_cycling          = .false. ! filter will advance the model

!----------------------------------------------------------------
! Namelist input with default values
!
! Set of values to control the application of prior and posterior inflation
logical :: do_prior_inflate     = .false.
logical :: do_posterior_inflate = .false.
logical :: prior_inflate_from_restart = .false.
logical :: posterior_inflate_from_restart = .false.

integer  :: async = 0, ens_size = 20
integer  :: tasks_per_model_advance = 1
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days    = 0
integer  :: init_time_seconds = 0
! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days      = -1
integer  :: first_obs_seconds   = -1
integer  :: last_obs_days       = -1
integer  :: last_obs_seconds    = -1
! Assimilation window; defaults to model timestep size.
integer  :: obs_window_days     = -1
integer  :: obs_window_seconds  = -1
! Control diagnostic output for state variables
integer  :: num_output_state_members = 0
integer  :: num_output_obs_members   = 0
integer  :: output_interval     = 1
integer  :: num_groups          = 1
logical  :: output_forward_op_errors = .false.
logical  :: output_timestamps        = .false.
logical  :: trace_execution          = .false.
logical  :: silence                  = .false.
logical  :: distributed_state = .true. ! Default to do distributed forward operators.

! IO options
!>@todo FIXME - how does this work for multiple domains?  ens1d1, ens2d1, ... ens1d2 or
!> ens1d1 ens1d2, ens1d1 ens2d2, etc   i like the latter better.
character(len=256) ::  input_state_files(MAX_FILES) = '' 
character(len=256) :: output_state_files(MAX_FILES) = '' 

! Name of files containing a list of {input,output} restart files, 1 file per domain
character(len=256) ::  input_state_file_list(MAX_NUM_DOMS) = '' 
character(len=256) :: output_state_file_list(MAX_NUM_DOMS) = ''

! Read in a single file and perturb this to create an ensemble
logical  :: perturb_from_single_instance = .false.
real(r8) :: perturbation_amplitude       = 0.2_r8

! File options.  Single vs. Multiple.  really 'unified' or 'combination' vs 'individual'
logical  :: single_file_in  = .false. ! all copies read  from 1 file
logical  :: single_file_out = .false. ! all copies written to 1 file

! optimization option:
logical :: compute_posterior   = .true. ! set to false to not compute posterior values

! Specify which diagnostic stages should be output
logical :: output_forecast_diags  = .true.
logical :: output_preassim_diags  = .true.
logical :: output_postassim_diags = .true.
logical :: output_analysis_diags  = .true.

! What quantities go in the diagnostic files
logical :: output_members   = .true.
logical :: output_mean      = .true.
logical :: output_sd        = .true.

character(len=256) :: obs_sequence_in_name  = "obs_seq.out",    &
                      obs_sequence_out_name = "obs_seq.final",  &
                      adv_ens_command       = './advance_model.csh'

! Some models are allowed to have MISSING_R8 values in the DART state vector.
! If they are encountered, it is not necessarily a FATAL error.
! Most of the time, if a MISSING_R8 is encountered, DART should die.
! CLM should have allow_missing_clm = .true.
logical  :: allow_missing_clm = .false.


namelist /filter_nml/            &
   do_prior_inflate,             &
   do_posterior_inflate,         &
   prior_inflate_from_restart,   &
   posterior_inflate_from_restart, &
   async,                        &
   adv_ens_command,              &
   ens_size,                     &
   tasks_per_model_advance,      &
   output_members,               &
   obs_sequence_in_name,         &
   obs_sequence_out_name,        &
   init_time_days,               &
   init_time_seconds,            &
   first_obs_days,               &
   first_obs_seconds,            &
   last_obs_days,                &
   last_obs_seconds,             &
   obs_window_days,              &
   obs_window_seconds,           &
   num_output_state_members,     &
   num_output_obs_members,       &
   output_interval,              &
   num_groups,                   &
   trace_execution,              &
   output_forward_op_errors,     &
   output_timestamps,            &
   silence,                      &
   distributed_state,            &
   single_file_in,               &
   single_file_out,              &
   perturb_from_single_instance, &
   perturbation_amplitude,       &
   compute_posterior,            &
   input_state_files,            &
   output_state_files,           &
   output_state_file_list,       &
   input_state_file_list,        &
   output_forecast_diags,        & 
   output_preassim_diags,        &
   output_postassim_diags,       &
   output_analysis_diags,        &
   output_mean,                  &
   output_sd,                    &
   allow_missing_clm

!----------------------------------------------------------------

contains

!----------------------------------------------------------------
!> The code does not use %vars arrays except:
!> * Task 0 still writes the obs_sequence file, so there is a transpose (copies to vars)
!>   and sending the obs_fwd_op_ens_handle%vars to task 0. Keys is also size obs%vars.
!> * If you read dart restarts state_ens_handle%vars is allocated.
!> * If you write dart diagnostics state_ens_handle%vars is allocated.
!> * If you are not doing distributed forward operators state_ens_handle%vars is allocated
subroutine filter_main()

type(ensemble_type)         :: state_ens_handle, obs_fwd_op_ens_handle, qc_ens_handle
type(obs_sequence_type)     :: seq
type(forward_op_info_type)  :: forward_op_ens_info
type(time_type)             :: time1, first_obs_time, last_obs_time
type(time_type)             :: curr_ens_time, next_ens_time, window_time

integer,    allocatable :: keys(:)
integer(i8)             :: model_size
integer                 :: iunit, io, time_step_number, num_obs_in_set, ntimes
integer                 :: last_key_used, key_bounds(2)
integer                 :: num_state_ens_copies
logical                 :: read_time_from_file

integer :: num_extras ! the extra ensemble copies

type(file_info_type) :: file_info_read
type(file_info_type) :: file_info_forecast
type(file_info_type) :: file_info_preassim
type(file_info_type) :: file_info_postassim
type(file_info_type) :: file_info_analysis
type(file_info_type) :: file_info_output

logical :: all_gone, allow_missing

call filter_initialize_modules_used() ! static_init_model called in here

! Read the namelist entry
call find_namelist_in_file("input.nml", "filter_nml", iunit)
read(iunit, nml = filter_nml, iostat = io)
call check_namelist_read(iunit, io, "filter_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=filter_nml)
if (do_nml_term()) write(     *     , nml=filter_nml)

if (task_count() == 1) distributed_state = .true.

call set_debug_fwd_op(output_forward_op_errors)
call set_trace(trace_execution, output_timestamps, silence)

! Make sure ensemble size is at least 2 (NEED MANY OTHER CHECKS)
if(ens_size < 2) then
   write(msgstring, *) 'ens_size in namelist is ', ens_size, ': Must be > 1'
   call error_handler(E_ERR,'filter_main', msgstring, source)
endif

! informational message to log
write(msgstring, '(A,I5)') 'running with an ensemble size of ', ens_size
call error_handler(E_MSG,'filter_main:', msgstring, source)

call set_missing_ok_status(allow_missing_clm)
allow_missing = get_missing_ok_status()

! Initialize the adaptive inflation module
call adaptive_inflate_init(prior_inflate)
! Turn it off if not requested in namelist
if(.not. do_prior_inflate) call set_inflate_flavor(prior_inflate, NO_INFLATION)

! Avoid illegal use of RTPS with prior inflation; Should throw an error?
if(do_rtps_inflate(prior_inflate)) call set_inflate_flavor(prior_inflate, NO_INFLATION)

call adaptive_inflate_init(post_inflate)
! Turn it off if not requested in namelist
if(.not. do_posterior_inflate) call set_inflate_flavor(post_inflate, NO_INFLATION)

! Cannot select state space posterior inflation options if not computing posterior
if(.not. compute_posterior .and. do_ss_inflate(post_inflate)) then
   write(msgstring, *) 'cannot use posterior state space inflation if compute_posterior is false'
   call error_handler(E_ERR,'filter_main', msgstring, source, &
           text2='"compute_posterior" is false; cannot have posterior state_space inflation')
endif

! for now, set 'has_cycling' to match 'single_file_out' since we're only supporting
! multi-file output for a single pass through filter, and allowing cycling if we're
! writing to a single file.

has_cycling = single_file_out

! Setup the indices into the ensemble storage:

! Can't output more ensemble members than exist
if(num_output_state_members > ens_size) num_output_state_members = ens_size
if(num_output_obs_members   > ens_size) num_output_obs_members   = ens_size

! Count and set up State copy numbers
call count_state_ens_copies(ens_size, output_mean, output_sd, do_prior_inflate, &
   do_posterior_inflate, prior_inflate, post_inflate, num_state_ens_copies, num_extras)

! Allocate model size storage and ens_size storage for metadata for outputting ensembles
model_size = get_model_size()

if(distributed_state) then
   call init_ensemble_manager(state_ens_handle, num_state_ens_copies, model_size)
   msgstring = 'running with distributed state; model states stay distributed across all tasks for the entire run'
else
   call init_ensemble_manager(state_ens_handle, num_state_ens_copies, model_size, transpose_type_in = 2)
   msgstring = 'running without distributed state; model states are gathered by ensemble for forward operators'
endif
! don't print if running single task.  transposes don't matter in this case.
if (task_count() > 1) &
   call error_handler(E_MSG,'filter_main:', msgstring, source)

call set_num_extra_copies(state_ens_handle, num_extras)

! Don't currently support number of processes > model_size
if(task_count() > model_size) then 
   write(msgstring, *) 'number of MPI processes = ', task_count(), &
                       ' while model size = ', model_size
   call error_handler(E_ERR,'filter_main', &
      'Cannot have number of processes > model size' ,source, text2=msgstring)
endif

if(.not. compute_posterior) then
   msgstring = 'skipping computation of posterior forward operators'
   call error_handler(E_MSG,'filter_main:', msgstring, source)
endif

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time(init_time_days, init_time_seconds, time1, read_time_from_file)

! for now, assume that we only allow cycling if single_file_out is true.
! code in this call needs to know how to initialize the output files.

call init_output_file_info('output', state_ens_handle, num_state_ens_copies, ens_size, &
   file_info_output, output_state_files, output_state_file_list, single_file_out,      &
   has_cycling, do_prior_inflate, do_posterior_inflate)

! Initialize and then read the input file
call read_state_and_inflation(num_state_ens_copies, state_ens_handle, ens_size, single_file_in, &
   perturb_from_single_instance, has_cycling, input_state_files, input_state_file_list, &
   prior_inflate, do_prior_inflate, prior_inflate_from_restart,                          &
   post_inflate, do_posterior_inflate, posterior_inflate_from_restart,              &
   file_info_read, read_time_from_file, time1)

if(iam_task0()) then
   ! Print out info for posterior inflation handle because it can include rtps
   call log_inflation_info(post_inflate)
endif

if (perturb_from_single_instance) then
   ! Only zero has the time, so broadcast the time to all other copy owners
   call broadcast_time_across_copy_owners(state_ens_handle, time1)
   call create_ensemble_from_single_file(state_ens_handle, ens_size, &
      perturbation_amplitude, get_missing_ok_status())
endif

! see what our stance is on missing values in the state vector
allow_missing = get_missing_ok_status()

! Initialize the obs_sequence and metadata; every pe gets a copy for now
call filter_setup_obs_sequence(forward_op_ens_info, seq, num_output_obs_members, &
   obs_sequence_in_name, compute_posterior)

! Need to find first obs with appropriate time, delete all earlier ones
if(first_obs_seconds > 0 .or. first_obs_days > 0) then
   first_obs_time = set_time(first_obs_seconds, first_obs_days)
   call delete_seq_head(first_obs_time, seq, all_gone)
   if(all_gone) then
      msgstring = 'All obs in sequence are before first_obs_days:first_obs_seconds'
      call error_handler(E_ERR,'filter_main',msgstring,source)
   endif
endif

! Start assimilating at beginning of modified sequence
last_key_used = -99

! Also get rid of observations past the last_obs_time if requested
if(last_obs_seconds >= 0 .or. last_obs_days >= 0) then
   last_obs_time = set_time(last_obs_seconds, last_obs_days)
   call delete_seq_tail(last_obs_time, seq, all_gone)
   if(all_gone) then
      msgstring = 'All obs in sequence are after last_obs_days:last_obs_seconds'
      call error_handler(E_ERR,'filter_main',msgstring,source)
   endif
endif

! Time step number is used to do periodic diagnostic output
time_step_number = -1
curr_ens_time = set_time(0, 0)
next_ens_time = set_time(0, 0)
call filter_set_window_time(window_time)

AdvanceTime : do
   time_step_number = time_step_number + 1

   ! Determine how far to advance model to make the window include the next available observation.
   call move_ahead(state_ens_handle, ens_size, seq, last_key_used, window_time, &
      key_bounds, num_obs_in_set, curr_ens_time, next_ens_time)

   ! The last key used is updated to move forward in the observation sequence
   last_key_used = key_bounds(2)

   ! Process 0 broadcast its value of key_bounds so exit condition of no keys can be checked
   call filter_sync_keys_time(state_ens_handle, key_bounds, num_obs_in_set, &
                              curr_ens_time, next_ens_time)

   if(key_bounds(1) < 0) exit AdvanceTime
   
   ! if model state data not at required time, advance model
   if (curr_ens_time /= next_ens_time) then

      ! we are going to advance the model - make sure we're doing single file output
      if (.not. has_cycling) then
         call error_handler(E_ERR,'filter:', &
             'advancing the model inside filter and multiple file output not currently supported', &
             source, text2='support will be added in subsequent releases', &
             text3='set "single_file_out=.true" for filter to advance the model, or advance the model outside filter')
      endif
 
      ! Is this sync still needed given the filter_sync call above
      call task_sync()

      ! make sure storage is allocated in ensemble manager for vars.
      call allocate_vars(state_ens_handle)

      ! Transpose to var complete so models can run on a single process
      call all_copies_to_all_vars(state_ens_handle)

      call advance_state(state_ens_handle, ens_size, next_ens_time, async, &
              adv_ens_command, tasks_per_model_advance, file_info_output, file_info_read)

      call all_vars_to_all_copies(state_ens_handle)

      ! update so curr time is accurate.
      curr_ens_time = next_ens_time
      state_ens_handle%current_time = curr_ens_time
      call set_time_on_extra_copies(state_ens_handle)

      ! only need to sync here since we want to wait for the
      ! slowest task to finish before outputting the time.
      call task_sync()
   endif

   ! Write out forecast diagnostic file(s). 
   if(output_forecast_diags .and. output_diag_now(output_interval, time_step_number)) &
      call output_diagnostics('forecast', state_ens_handle, num_state_ens_copies, ens_size, &
         num_output_state_members, &
         file_info_forecast, single_file_out, has_cycling, output_mean, output_sd,           &
         output_members, do_prior_inflate, do_posterior_inflate, ENS_MEAN_COPY, ENS_SD_COPY)
   
   ! Apply prior inflation 
   if(do_ss_inflate(prior_inflate)) &
      call filter_ensemble_inflate(state_ens_handle, PRIOR_INF_COPY, prior_inflate, ENS_MEAN_COPY)

   ! if relaxation-to-prior-spread inflation, save prior spread in copy RTPS_PRIOR_SPREAD_COPY
   if ( do_rtps_inflate(post_inflate) ) &
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, RTPS_PRIOR_SPREAD_COPY)

   ! Write out preassim diagnostic files if requested.
   if(output_preassim_diags .and. output_diag_now(output_interval, time_step_number)) &
      call output_diagnostics('preassim', state_ens_handle, num_state_ens_copies, ens_size, &
         num_output_state_members, &
         file_info_preassim, single_file_out, has_cycling, output_mean, output_sd,           &
         output_members, do_prior_inflate, do_posterior_inflate, ENS_MEAN_COPY, ENS_SD_COPY)

   ! Compute the forward operators and fill data structures
   call forward_operators(forward_op_ens_info, state_ens_handle, obs_fwd_op_ens_handle, &
      qc_ens_handle, seq, ens_size, num_groups, num_obs_in_set, keys, key_bounds,       &
      num_output_obs_members, compute_posterior, isprior = .true.) 

   call filter_assim(state_ens_handle, obs_fwd_op_ens_handle, forward_op_ens_info, seq, keys, &
      ens_size, num_groups, prior_inflate, ENS_MEAN_COPY, ENS_SD_COPY,         &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, inflate_only = .false.)

   ! Write out postassim diagnostic files if requested.  This contains the assimilated ensemble 
   ! JLA DEVELOPMENT: This used to output the damped inflation. NO LONGER.
   if(output_postassim_diags .and. output_diag_now(output_interval, time_step_number)) &
      call output_diagnostics('postassim', state_ens_handle, num_state_ens_copies, ens_size, &
         num_output_state_members, &
         file_info_postassim, single_file_out, has_cycling, output_mean, output_sd,           &
         output_members, do_prior_inflate, do_posterior_inflate, ENS_MEAN_COPY, ENS_SD_COPY)

   ! This block applies posterior inflation including RTPS if selected
   if(do_ss_inflate(post_inflate) .or. do_rtps_inflate(post_inflate))             &
      call filter_ensemble_inflate(state_ens_handle, POST_INF_COPY, post_inflate, &
                       ENS_MEAN_COPY, RTPS_PRIOR_SPREAD_COPY, ENS_SD_COPY)

   ! this block recomputes the expected obs values for the obs_seq.final file
   if (compute_posterior) then
      ! Compute the forward operators and fill data structures
      call forward_operators(forward_op_ens_info, state_ens_handle, obs_fwd_op_ens_handle, &
         qc_ens_handle, seq, ens_size, num_groups, num_obs_in_set, keys, key_bounds,       & 
         num_output_obs_members, compute_posterior, isprior = .false.) 
   else
      ! Collect any updated QC values that may ! have been set in the assimilation loop
      call obs_space_sync_QCs(forward_op_ens_info, obs_fwd_op_ens_handle, &
         seq, keys, num_obs_in_set)
   endif

   ! Compute the adaptive state space posterior inflation
   if(do_ss_inflate(post_inflate) .and. ( .not. do_rtps_inflate(post_inflate)) )                 &
      call filter_assim(state_ens_handle, obs_fwd_op_ens_handle, forward_op_ens_info, seq, keys, &
         ens_size, num_groups, post_inflate, ENS_MEAN_COPY, ENS_SD_COPY,                         &
         POST_INF_COPY, POST_INF_SD_COPY, inflate_only = .true.)

   ! Free up all the allocated space associated with obs ensemble
   call end_ensemble_manager(obs_fwd_op_ens_handle)
   call end_ensemble_manager(qc_ens_handle)
   deallocate(keys)

   ! Write out analysis diagnostic files if requested. 
   if(output_analysis_diags .and. output_diag_now(output_interval, time_step_number)) &
      call output_diagnostics('analysis', state_ens_handle, num_state_ens_copies, ens_size, &
         num_output_state_members, &
         file_info_analysis, single_file_out, has_cycling, output_mean, output_sd,           &
         output_members, do_prior_inflate, do_posterior_inflate, ENS_MEAN_COPY, ENS_SD_COPY)

end do AdvanceTime

! Output the adjusted ensemble. If cycling only the last timestep is writen out
call write_state(state_ens_handle, file_info_output)

! Only pe 0 outputs the observation space diagnostic file
if(my_task_id() == 0) call write_obs_seq(seq, obs_sequence_out_name)

! close the diagnostic/restart netcdf files
if (single_file_out) then
   call finalize_single_file_io(file_info_output)
   if(output_forecast_diags) call finalize_single_file_io(file_info_forecast)
   if(output_preassim_diags) call finalize_single_file_io(file_info_preassim)
   if(output_postassim_diags) call finalize_single_file_io(file_info_postassim)
   if(output_analysis_diags) call finalize_single_file_io(file_info_analysis)
endif

! Give the model_mod code a chance to clean up.
call end_assim_model()

! deallocate qceff_table_data structures
call end_algorithm_info_mod()

call end_ensemble_manager(state_ens_handle)

! Free up the obs sequence
call destroy_obs_sequence(seq)
if(my_task_id() == 0) then
   write(logfileunit,*)'FINISHED filter.'
   write(logfileunit,*)
endif

end subroutine filter_main

!-----------------------------------------------------------

subroutine filter_initialize_modules_used()

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
call state_vector_io_init()
call initialize_qc()

! Initialize algorothm_info_mod and read in QCF table data
call init_algorithm_info_mod()

end subroutine filter_initialize_modules_used

!-------------------------------------------------------------------------

subroutine filter_set_initial_time(days, seconds, dart_time, read_time_from_file)

integer,         intent(in)  :: days, seconds
type(time_type), intent(out) :: dart_time
logical,         intent(out) :: read_time_from_file

if(days >= 0) then
   dart_time = set_time(seconds, days)
   read_time_from_file = .false.
else
   dart_time = set_time(0, 0)
   read_time_from_file = .true.
endif

end subroutine filter_set_initial_time

!-------------------------------------------------------------------------

subroutine filter_set_window_time(dart_time)

type(time_type), intent(out) :: dart_time


if(obs_window_days >= 0) then
   dart_time = set_time(obs_window_seconds, obs_window_days)
else
   dart_time = set_time(0, 0)
endif

end subroutine filter_set_window_time

!-------------------------------------------------------------------------

subroutine filter_ensemble_inflate(ens_handle, inflate_copy, inflate_handle, ENS_MEAN_COPY, &
                                   RTPS_PRIOR_SPREAD_COPY, ENS_SD_COPY)

type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: inflate_copy, ENS_MEAN_COPY
type(adaptive_inflate_type), intent(inout) :: inflate_handle
integer, optional,           intent(in)    :: RTPS_PRIOR_SPREAD_COPY, ENS_SD_COPY

integer :: j, group, grp_bot, grp_top, grp_size
type(location_type) :: my_state_loc
integer :: my_state_kind
type(distribution_params_type) :: dist_params
real(r8) :: probit_ens(ens_size), probit_ens_mean
logical  :: bounded_below, bounded_above
real(r8) :: lower_bound,   upper_bound
integer  :: dist_type

! Inflate each group separately;  Divide ensemble into num_groups groups
grp_size = ens_size / num_groups

do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1
   ! Compute the mean for this group
   call compute_copy_mean(ens_handle, grp_bot, grp_top, ENS_MEAN_COPY)

   if ( do_rtps_inflate(inflate_handle)) then 
      do j = 1, ens_handle%my_num_vars
         call inflate_ens(inflate_handle, ens_handle%copies(grp_bot:grp_top, j), &
            ens_handle%copies(ENS_MEAN_COPY, j), ens_handle%copies(inflate_copy, j), 0.0_r8, &
            ens_handle%copies(RTPS_PRIOR_SPREAD_COPY, j), ens_handle%copies(ENS_SD_COPY, j)) 
      end do 
   else 

      ! Doing inflation in probit space; do probit probability integral transform
      do j = 1, ens_handle%my_num_vars
         call get_state_meta_data(ens_handle%my_vars(j), my_state_loc, my_state_kind)    

         ! Need to specify what kind of prior to use for each
         call probit_dist_info(my_state_kind, .true., .true., dist_type, &
            bounded_below, bounded_above, lower_bound, upper_bound)

         call transform_to_probit(grp_size, ens_handle%copies(grp_bot:grp_top, j), &
            dist_type, dist_params, probit_ens(1:grp_size), .false., &
               bounded_below, bounded_above, lower_bound, upper_bound)

         ! Compute the ensemble mean in transformed space
         probit_ens_mean = sum(probit_ens(1:grp_size)) / grp_size
         ! Inflate in probit space
         call inflate_ens(inflate_handle, probit_ens(1:grp_size), probit_ens_mean, &
            ens_handle%copies(inflate_copy, j))
         ! Transform back from probit space
         call transform_from_probit(grp_size, probit_ens(1:grp_size), &
            dist_params, ens_handle%copies(grp_bot:grp_top, j))
      end do
   endif
end do

end subroutine filter_ensemble_inflate

!-------------------------------------------------------------------------

subroutine filter_sync_keys_time(ens_handle, key_bounds, num_obs_in_set, time1, time2)

integer,             intent(inout)  :: key_bounds(2), num_obs_in_set
type(time_type),     intent(inout)  :: time1, time2
type(ensemble_type), intent(inout)     :: ens_handle

! Have owner of copy 1 broadcast these values to all other tasks.
! Only tasks which contain copies have this info; doing it this way
! allows ntasks > nens to work.

real(r8) :: rkey_bounds(2), rnum_obs_in_set(1)
real(r8) :: rtime(4)
integer  :: days, secs
integer  :: copy1_owner, owner_index

call get_copy_owner_index(ens_handle, 1, copy1_owner, owner_index)

if( ens_handle%my_pe == copy1_owner) then
   rkey_bounds = key_bounds
   rnum_obs_in_set(1) = num_obs_in_set
   call get_time(time1, secs, days)
   rtime(1) = secs
   rtime(2) = days
   call get_time(time2, secs, days)
   rtime(3) = secs
   rtime(4) = days
   call broadcast_send(map_pe_to_task(ens_handle, copy1_owner), rkey_bounds, rnum_obs_in_set, rtime)
else
   call broadcast_recv(map_pe_to_task(ens_handle, copy1_owner), rkey_bounds, rnum_obs_in_set, rtime)
   key_bounds =     nint(rkey_bounds)
   num_obs_in_set = nint(rnum_obs_in_set(1))
   time1 = set_time(nint(rtime(1)), nint(rtime(2)))
   time2 = set_time(nint(rtime(3)), nint(rtime(4)))
endif

! Every task gets the current time (necessary for the forward operator)
ens_handle%current_time = time1

end subroutine filter_sync_keys_time

!-------------------------------------------------------------------------
! Only copy 1 on task zero has the correct time after reading
! when you read one instance using filter_read_restart.
! perturb_from_single_instance = .true.
! This routine makes the times consistent across the ensemble. 
! Any task that owns one or more state vectors needs the time for
! the move ahead call.
!> @todo This is broadcasting the time to all tasks, not
!> just the tasks that own copies.

subroutine broadcast_time_across_copy_owners(ens_handle, ens_time)

type(ensemble_type), intent(inout) :: ens_handle
type(time_type),     intent(in)    :: ens_time

real(r8) :: rtime(2)
integer  :: days, secs
integer  :: copy1_owner, owner_index
type(time_type) :: time_from_copy1

call get_copy_owner_index(ens_handle, 1, copy1_owner, owner_index)

if( ens_handle%my_pe == copy1_owner) then
   call get_time(ens_time, secs, days)
   rtime(1) = secs
   rtime(2) = days
   call broadcast_send(map_pe_to_task(ens_handle, copy1_owner), rtime)
   ens_handle%time(1:ens_handle%my_num_copies) = ens_time
else
   call broadcast_recv(map_pe_to_task(ens_handle, copy1_owner), rtime)
   time_from_copy1 = set_time(nint(rtime(1)), nint(rtime(2)))
   if (ens_handle%my_num_copies > 0) ens_handle%time(1:ens_handle%my_num_copies) = time_from_copy1
endif

end subroutine broadcast_time_across_copy_owners

!-------------------------------------------------------------------------

subroutine set_trace(trace_execution, output_timestamps, silence)

logical, intent(in) :: trace_execution
logical, intent(in) :: output_timestamps
logical, intent(in) :: silence

! Set whether other modules trace execution with messages
! and whether they output timestamps to trace overall performance

! defaults
trace_level     = 0
timestamp_level = 0

! selectively turn stuff back on
if (trace_execution)   trace_level     = 1
if (output_timestamps) timestamp_level = 1

! turn as much off as possible
if (silence) then
   trace_level     = -1
   timestamp_level = -1
endif

call set_obs_model_trace(trace_level, timestamp_level)
call set_assim_tools_trace(trace_level, timestamp_level)

end subroutine set_trace

!-------------------------------------------------------------------------
!> Set the time on any extra copies that a pe owns
!> Could we just set the time on all copies?

subroutine set_time_on_extra_copies(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

integer :: copy_num, owner, owners_index
integer :: ens_size

ens_size = ens_handle%num_copies - ens_handle%num_extras

do copy_num = ens_size + 1, ens_handle%num_copies
   ! Set time for a given copy of an ensemble
   call get_copy_owner_index(ens_handle, copy_num, owner, owners_index)
   if(ens_handle%my_pe == owner) then
      call set_ensemble_time(ens_handle, owners_index, ens_handle%current_time)
   endif
enddo

end subroutine  set_time_on_extra_copies

!------------------------------------------------------------------

function output_diag_now(output_interval, time_step_number)

logical             :: output_diag_now
integer, intent(in) :: output_interval, time_step_number

output_diag_now = ((output_interval > 0) .and. &  
   (time_step_number / output_interval * output_interval == time_step_number)) 

end function output_diag_now

!------------------------------------------------------------------

end module filter_mod

