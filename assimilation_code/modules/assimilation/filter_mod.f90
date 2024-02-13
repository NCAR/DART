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
                                  do_output, find_namelist_in_file, check_namelist_read,      &
                                  open_file, close_file, do_nml_file, do_nml_term, to_upper,  &
                                  set_multiple_filename_lists

use assim_model_mod,       only : static_init_assim_model, get_model_size,                    &
                                  end_assim_model,  pert_model_copies, get_state_meta_data

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
                                  allocate_single_copy, allocate_vars,                        &
                                  deallocate_single_copy
                                  

use adaptive_inflate_mod,  only : do_ss_inflate, &
                                  inflate_ens, adaptive_inflate_init,                 &
                                  adaptive_inflate_type, log_inflation_info,          &
                                  do_rtps_inflate, set_inflate_flavor,                &
                                  NO_INFLATION
                                  

use mpi_utilities_mod,     only : my_task_id, task_sync, broadcast_send, broadcast_recv,      &
                                  task_count, iam_task0

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use state_vector_io_mod,   only : state_vector_io_init, read_state, write_state, &
                                  set_stage_to_write, get_stage_to_write

use io_filenames_mod,      only : io_filenames_init, file_info_type, &
                                  combine_file_info, set_file_metadata,  &
                                  set_member_file_metadata,  set_io_copy_flag, &
                                  check_file_info_variable_shape, &
                                  query_copy_present, COPY_NOT_PRESENT, &
                                  READ_COPY, WRITE_COPY

use direct_netcdf_mod,     only : finalize_single_file_io

use state_structure_mod,   only : get_num_domains

use forward_operator_mod,  only : forward_operators, forward_op_info_type, obs_space_sync_QCs, &
                                  filter_setup_obs_sequence

use quality_control_mod,   only : initialize_qc

use location_mod,          only : location_type

use probit_transform_mod,  only : transform_to_probit, transform_from_probit

use algorithm_info_mod, only : probit_dist_info, init_algorithm_info_mod, end_algorithm_info_mod

use distribution_params_mod, only : distribution_params_type

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

! Identifier for different copies for diagnostic files
integer, parameter :: MEM_START     = 1
integer, parameter :: MEM_END       = 2
integer, parameter :: ENS_MEAN      = 3
integer, parameter :: ENS_SD        = 4
integer, parameter :: PRIORINF_MEAN = 5
integer, parameter :: PRIORINF_SD   = 6
integer, parameter :: POSTINF_MEAN  = 7
integer, parameter :: POSTINF_SD    = 8

! Number of Stage Copies
integer, parameter :: NUM_SCOPIES    = 8

! Ensemble copy numbers
integer :: ENS_MEM_START                   = COPY_NOT_PRESENT
integer :: ENS_MEM_END                     = COPY_NOT_PRESENT
integer :: ENS_MEAN_COPY                   = COPY_NOT_PRESENT
integer :: ENS_SD_COPY                     = COPY_NOT_PRESENT
integer :: PRIOR_INF_COPY                  = COPY_NOT_PRESENT
integer :: PRIOR_INF_SD_COPY               = COPY_NOT_PRESENT
integer :: POST_INF_COPY                   = COPY_NOT_PRESENT
integer :: POST_INF_SD_COPY                = COPY_NOT_PRESENT

integer :: INPUT_COPIES( NUM_SCOPIES )     = COPY_NOT_PRESENT

integer ::   CURRENT_COPIES( NUM_SCOPIES ) = COPY_NOT_PRESENT
integer ::  FORECAST_COPIES( NUM_SCOPIES ) = COPY_NOT_PRESENT
integer ::  PREASSIM_COPIES( NUM_SCOPIES ) = COPY_NOT_PRESENT
integer :: POSTASSIM_COPIES( NUM_SCOPIES ) = COPY_NOT_PRESENT
integer ::  ANALYSIS_COPIES( NUM_SCOPIES ) = COPY_NOT_PRESENT

integer :: RTPS_PRIOR_SPREAD              = COPY_NOT_PRESENT

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
logical  :: write_obs_every_cycle    = .false.  ! debug only
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

! Stages to write.  Valid values are:
! multi-file:    input, forecast, preassim, postassim, analysis, output
! single-file:          forecast, preassim, postassim, analysis, output
character(len=10)  :: stages_to_write(6) = (/"output    ", "null      ", "null      ", &
                                             "null      ", "null      ", "null      " /)

!>@todo FIXME
!> for preassim and postassim output it might be we should
!> be controlling the writing of individual ensemble members
!> by looking at the num_output_state_member value.  0 means
!> don't write any members, otherwise it's a count.  and for
!> completeness, there could be a count for pre and a count for post.

logical :: output_members   = .true.
logical :: output_mean      = .true.
logical :: output_sd        = .true.
logical :: write_all_stages_at_end = .false.

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
   stages_to_write,              &
   input_state_files,            &
   output_state_files,           &
   output_state_file_list,       &
   input_state_file_list,        &
   output_mean,                  &
   output_sd,                    &
   write_all_stages_at_end,      &
   write_obs_every_cycle,        & 
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

type(file_info_type) :: file_info_input
type(file_info_type) :: file_info_mean_sd
type(file_info_type) :: file_info_forecast
type(file_info_type) :: file_info_preassim
type(file_info_type) :: file_info_postassim
type(file_info_type) :: file_info_analysis
type(file_info_type) :: file_info_output
type(file_info_type) :: file_info_all

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

! don't allow cycling and write all at end - might never be supported
if (has_cycling .and. write_all_stages_at_end) then
   call error_handler(E_ERR,'filter:', &
         'advancing the model inside filter and writing all state data at end not supported', &
          source, text2='delaying write until end only supported when advancing model outside filter', &
          text3='set "write_all_stages_at_end=.false." to cycle and write data as it is computed')
endif

! Setup the indices into the ensemble storage:

! Can't output more ensemble members than exist
if(num_output_state_members > ens_size) num_output_state_members = ens_size
if(num_output_obs_members   > ens_size) num_output_obs_members   = ens_size

! Set up stages to write : input, preassim, postassim, output
call parse_stages_to_write(stages_to_write)

! Count and set up State copy numbers
num_state_ens_copies = count_state_ens_copies(ens_size)
num_extras           = num_state_ens_copies - ens_size

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
call initialize_file_information(num_state_ens_copies ,                     &
                                 file_info_input      , file_info_mean_sd,  &
                                 file_info_forecast   , file_info_preassim, &
                                 file_info_postassim  , file_info_analysis, &
                                 file_info_output)

call check_file_info_variable_shape(file_info_output, state_ens_handle)

call read_state(state_ens_handle, file_info_input, read_time_from_file, time1,      &
                PRIOR_INF_COPY, PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY, &
                prior_inflate, post_inflate,                                        &
                prior_inflate_from_restart, posterior_inflate_from_restart,         &
                perturb_from_single_instance)

if(iam_task0()) then
   ! Print out info for posterior inflation handle because it can include rtps
   call log_inflation_info(post_inflate)
endif

if (perturb_from_single_instance) then
   ! Only zero has the time, so broadcast the time to all other copy owners
   call broadcast_time_across_copy_owners(state_ens_handle, time1)
   call create_ensemble_from_single_file(state_ens_handle)
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
              adv_ens_command, tasks_per_model_advance, file_info_output, file_info_input)

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
   call do_stage_output('forecast', output_interval, time_step_number, write_all_stages_at_end, &
      state_ens_handle, FORECAST_COPIES, file_info_forecast, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)
   
   ! Apply prior inflation 
   if(do_ss_inflate(prior_inflate)) &
      call filter_ensemble_inflate(state_ens_handle, PRIOR_INF_COPY, prior_inflate, ENS_MEAN_COPY)

   ! if relaxation-to-prior-spread inflation, save prior spread in copy RTPS_PRIOR_SPREAD
   if ( do_rtps_inflate(post_inflate) ) &
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, RTPS_PRIOR_SPREAD)

   ! Write out preassim diagnostic files if requested.
   call do_stage_output('preassim', output_interval, time_step_number, write_all_stages_at_end, &
      state_ens_handle, PREASSIM_COPIES, file_info_preassim, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

   ! Compute the forward operators and fill data structures
   call forward_operators(forward_op_ens_info, state_ens_handle, obs_fwd_op_ens_handle, &
      qc_ens_handle, seq, ens_size, num_groups, num_obs_in_set, keys, key_bounds,       &
      num_output_obs_members, compute_posterior, isprior = .true.) 

   call filter_assim(state_ens_handle, obs_fwd_op_ens_handle, forward_op_ens_info, seq, keys, &
      ens_size, num_groups, prior_inflate, ENS_MEAN_COPY, ENS_SD_COPY,         &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, inflate_only = .false.)

   ! Write out postassim diagnostic files if requested.  This contains the assimilated ensemble 
   ! JLA DEVELOPMENT: This used to output the damped inflation. NO LONGER.
   call do_stage_output('postassim', output_interval, time_step_number, write_all_stages_at_end, &
      state_ens_handle, POSTASSIM_COPIES, file_info_postassim, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

   ! This block applies posterior inflation including RTPS if selected
   if(do_ss_inflate(post_inflate) .or. do_rtps_inflate(post_inflate))             &
      call filter_ensemble_inflate(state_ens_handle, POST_INF_COPY, post_inflate, &
                       ENS_MEAN_COPY, RTPS_PRIOR_SPREAD, ENS_SD_COPY)

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

   ! this block computes the adaptive state space posterior inflation
   if(do_ss_inflate(post_inflate) .and. ( .not. do_rtps_inflate(post_inflate)) )                 &
      call filter_assim(state_ens_handle, obs_fwd_op_ens_handle, forward_op_ens_info, seq, keys, &
         ens_size, num_groups, post_inflate, ENS_MEAN_COPY, ENS_SD_COPY,                         &
         POST_INF_COPY, POST_INF_SD_COPY, inflate_only = .true.)

   ! Free up all the allocated space associated with obs ensemble
   call end_ensemble_manager(obs_fwd_op_ens_handle)
   call end_ensemble_manager(qc_ens_handle)
   deallocate(keys)

   ! Write out analysis diagnostic files if requested. 
   call do_stage_output('analysis', output_interval, time_step_number, write_all_stages_at_end, &
      state_ens_handle, ANALYSIS_COPIES, file_info_analysis, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

end do AdvanceTime

! Output the adjusted ensemble. If cycling only the last timestep is writen out
if (get_stage_to_write('output')) then
      ! will write outside loop
      if (.not. write_all_stages_at_end) &
         call write_state(state_ens_handle, file_info_output)
endif

! Only pe 0 outputs the observation space diagnostic file
if(my_task_id() == 0) call write_obs_seq(seq, obs_sequence_out_name)

! Output all restart files if requested
if (write_all_stages_at_end) then
   file_info_all = combine_file_info( &
                    (/file_info_input, file_info_mean_sd, file_info_forecast, &
                      file_info_preassim, file_info_postassim, file_info_analysis, &
                      file_info_output/) )

   call write_state(state_ens_handle, file_info_all)
endif

! close the diagnostic/restart netcdf files
if (single_file_out) then

   if (get_stage_to_write('forecast')) &
      call finalize_single_file_io(file_info_forecast)

   if (get_stage_to_write('preassim')) &
      call finalize_single_file_io(file_info_preassim)

   if (get_stage_to_write('postassim')) &
      call finalize_single_file_io(file_info_postassim)

   if (get_stage_to_write('analysis')) &
      call finalize_single_file_io(file_info_analysis)

   if (get_stage_to_write('output')) &
      call finalize_single_file_io(file_info_output)
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
                                   RTPS_PRIOR_SPREAD, ENS_SD_COPY)

type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: inflate_copy, ENS_MEAN_COPY
type(adaptive_inflate_type), intent(inout) :: inflate_handle
integer, optional,           intent(in)    :: RTPS_PRIOR_SPREAD, ENS_SD_COPY

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
            ens_handle%copies(RTPS_PRIOR_SPREAD, j), ens_handle%copies(ENS_SD_COPY, j)) 
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

subroutine print_ens_time(ens_handle, msg)

type(ensemble_type), intent(in) :: ens_handle
character(len=*), intent(in) :: msg

! Write message to stdout and log file.
type(time_type) :: mtime

if (trace_level <= 0) return

if (do_output()) then
   if (get_my_num_copies(ens_handle) < 1) return
   call get_ensemble_time(ens_handle, 1, mtime)
   call print_time(mtime, ' filter trace: '//msg, logfileunit)
   call print_time(mtime, ' filter trace: '//msg)
endif

end subroutine print_ens_time

!-------------------------------------------------------------------------

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

subroutine create_ensemble_from_single_file(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

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

call pert_model_copies(ens_handle, ens_size, perturbation_amplitude, interf_provided)
if (.not. interf_provided) then
   call perturb_copies_task_bitwise(ens_handle)
endif

! Restore the missing_r8
if (get_missing_ok_status()) then
   do i = 1, ens_size
      where(miss_me) ens_handle%copies(i, :) = missing_r8
   enddo
   deallocate(miss_me)
endif

end subroutine create_ensemble_from_single_file


!------------------------------------------------------------------
! Perturb the copies array in a way that is bitwise reproducible
! no matter how many task you run on.

subroutine perturb_copies_task_bitwise(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

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

!------------------------------------------------------------------
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
!> Copy the current mean, sd, inf_mean, inf_sd to spare copies
!> Assuming that if the spare copy is there you should fill it

subroutine store_copies(ens_handle, STAGE_COPIES)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(inout) :: STAGE_COPIES(NUM_SCOPIES)

integer :: i, offset

if (query_copy_present( STAGE_COPIES(ENS_MEAN)) ) &
   ens_handle%copies(   STAGE_COPIES(ENS_MEAN), :)      = ens_handle%copies(ENS_MEAN_COPY, :)

if (query_copy_present( STAGE_COPIES(ENS_SD)) ) &
   ens_handle%copies(   STAGE_COPIES(ENS_SD), :)        = ens_handle%copies(ENS_SD_COPY, :)

if (query_copy_present( STAGE_COPIES(PRIORINF_MEAN)) ) &
   ens_handle%copies(   STAGE_COPIES(PRIORINF_MEAN), :) = ens_handle%copies(PRIOR_INF_COPY, :)

if (query_copy_present( STAGE_COPIES(PRIORINF_SD)) ) &
   ens_handle%copies(   STAGE_COPIES(PRIORINF_SD), :)   = ens_handle%copies(PRIOR_INF_SD_COPY, :)

if (query_copy_present( STAGE_COPIES(POSTINF_MEAN)) ) &
   ens_handle%copies(   STAGE_COPIES(POSTINF_MEAN), :)  = ens_handle%copies(POST_INF_COPY, :)

if (query_copy_present( STAGE_COPIES(POSTINF_SD)) ) &
   ens_handle%copies(   STAGE_COPIES(POSTINF_SD), :)    = ens_handle%copies(POST_INF_SD_COPY, :)
     
do i = 1, num_output_state_members
   offset = STAGE_COPIES(MEM_START) + i - 1
   if ( query_copy_present(offset) ) ens_handle%copies(offset, :) = ens_handle%copies(i, :)
enddo

end subroutine store_copies


!------------------------------------------------------------------
!> Count the number of copies to be allocated for the ensemble manager

function count_state_ens_copies(ens_size) result(num_copies)

integer,                     intent(in) :: ens_size
integer :: num_copies

integer :: cnum = 0

! Filter Ensemble Members
!   ENS_MEM_XXXX
ENS_MEM_START = next_copy_number(cnum)
ENS_MEM_END   = next_copy_number(cnum, ens_size)

! Filter Extra Copies For Assimilation
!    ENS_MEAN_COPY    
!    ENS_SD_COPY      
!    PRIOR_INF_COPY   
!    PRIOR_INF_SD_COPY
!    POST_INF_COPY    
!    POST_INF_SD_COPY 

ENS_MEAN_COPY     = next_copy_number(cnum)
ENS_SD_COPY       = next_copy_number(cnum)
PRIOR_INF_COPY    = next_copy_number(cnum)
PRIOR_INF_SD_COPY = next_copy_number(cnum)
POST_INF_COPY     = next_copy_number(cnum)
POST_INF_SD_COPY  = next_copy_number(cnum)

! If there are no diagnostic files, we will need to store the
! copies that would have gone in Prior_Diag.nc and Posterior_Diag.nc
! in spare copies in the ensemble.

if (write_all_stages_at_end) then
   if (get_stage_to_write('input')) then
      ! Option to Output Input Mean and SD
      !   INPUT_MEAN
      !   INPUT_SD  
      if (output_mean) then
         INPUT_COPIES(ENS_MEAN) = next_copy_number(cnum)
         if ( do_prior_inflate .and. .not. prior_inflate_from_restart) then
            INPUT_COPIES(PRIORINF_MEAN) = next_copy_number(cnum)
         endif
         if ( do_posterior_inflate .and. .not. posterior_inflate_from_restart) then
            INPUT_COPIES(POSTINF_MEAN)  = next_copy_number(cnum)
         endif
      endif

      if (output_sd) then
         INPUT_COPIES(ENS_SD) = next_copy_number(cnum)
         if ( do_prior_inflate .and. .not. prior_inflate_from_restart) then
            INPUT_COPIES(PRIORINF_SD) = next_copy_number(cnum)
         endif
         if ( do_posterior_inflate .and. .not. posterior_inflate_from_restart) then
            INPUT_COPIES(POSTINF_SD)  = next_copy_number(cnum)
         endif
      endif
   endif

   if (get_stage_to_write('forecast')) &
      call set_copies( cnum, FORECAST_COPIES)

   if (get_stage_to_write('preassim')) &
      call set_copies( cnum, PREASSIM_COPIES)

   if (get_stage_to_write('postassim')) &
      call set_copies( cnum, POSTASSIM_COPIES)

   if (get_stage_to_write('analysis')) &
      call set_copies( cnum, ANALYSIS_COPIES)

else

   ! Write everything in stages
   ! Option to Output Input Mean and SD
   if (output_mean) then
      INPUT_COPIES(ENS_MEAN) = ENS_MEAN_COPY
      if ( do_prior_inflate .and. .not. prior_inflate_from_restart) then
         INPUT_COPIES(PRIORINF_MEAN) = PRIOR_INF_COPY
      endif
      if ( do_posterior_inflate .and. .not. posterior_inflate_from_restart) then
         INPUT_COPIES(POSTINF_MEAN)  = POST_INF_COPY
      endif
   endif

   if (output_sd) then
      INPUT_COPIES(ENS_SD) = ENS_SD_COPY
      if ( do_prior_inflate .and. .not. prior_inflate_from_restart) then
         INPUT_COPIES(PRIORINF_SD) = PRIOR_INF_SD_COPY
      endif
      if ( do_posterior_inflate .and. .not. posterior_inflate_from_restart) then
         INPUT_COPIES(POSTINF_SD)  = POST_INF_SD_COPY
      endif
   endif
  
   FORECAST_COPIES   = (/ ENS_MEM_START, ENS_MEM_END, ENS_MEAN_COPY, ENS_SD_COPY, &
                          PRIOR_INF_COPY, PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY /)

   PREASSIM_COPIES   = (/ ENS_MEM_START, ENS_MEM_END, ENS_MEAN_COPY, ENS_SD_COPY, &
                          PRIOR_INF_COPY, PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY /)

   POSTASSIM_COPIES  = (/ ENS_MEM_START, ENS_MEM_END, ENS_MEAN_COPY, ENS_SD_COPY, &
                         PRIOR_INF_COPY, PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY /)

   ANALYSIS_COPIES   = (/ ENS_MEM_START, ENS_MEM_END, ENS_MEAN_COPY, ENS_SD_COPY, &
                          PRIOR_INF_COPY, PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY /)
  
endif

CURRENT_COPIES    = (/ ENS_MEM_START, ENS_MEM_END, ENS_MEAN_COPY, ENS_SD_COPY, &
                       PRIOR_INF_COPY, PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY /)

! If Whitaker/Hamill (2012) relaxation-to-prior-spread (rpts) inflation
! then we need an extra copy to hold (save) the prior ensemble spread
! ENS_SD_COPY will be overwritten with the posterior spread before
! applying the inflation algorithm; must save the prior ensemble spread in a different copy
if ( do_rtps_inflate(post_inflate) ) then
   RTPS_PRIOR_SPREAD = next_copy_number(cnum)
endif

num_copies = cnum

end function count_state_ens_copies


!------------------------------------------------------------------
!> Set file name information.  For members restarts can be read from
!> an input_state_file_list or constructed using a stage name and
!> num_ens.  The file_info handle knows whether or not there is an
!> associated input_state_file_list. If no list is provided member
!> filenames are written as :
!>    stage_member_####.nc (ex. preassim_member_0001.nc)
!> extra copies are stored as :
!>    stage_basename.nc (ex. preassim_mean.nc)

subroutine set_filename_info(file_info, stage, num_ens, STAGE_COPIES) 

type(file_info_type), intent(inout) :: file_info
character(len=*),     intent(in)    :: stage
integer,              intent(in)    :: num_ens
integer,              intent(inout) :: STAGE_COPIES(NUM_SCOPIES)

call set_member_file_metadata(file_info, num_ens, STAGE_COPIES(MEM_START))


STAGE_COPIES(MEM_END) = STAGE_COPIES(MEM_START) + num_ens - 1

call set_file_metadata(file_info, STAGE_COPIES(ENS_MEAN),      stage, 'mean',          'ensemble mean')
call set_file_metadata(file_info, STAGE_COPIES(ENS_SD),        stage, 'sd',            'ensemble sd')
call set_file_metadata(file_info, STAGE_COPIES(PRIORINF_MEAN), stage, 'priorinf_mean', 'prior inflation mean')
call set_file_metadata(file_info, STAGE_COPIES(PRIORINF_SD),   stage, 'priorinf_sd',   'prior inflation sd')
call set_file_metadata(file_info, STAGE_COPIES(POSTINF_MEAN),  stage, 'postinf_mean',  'posterior inflation mean')
call set_file_metadata(file_info, STAGE_COPIES(POSTINF_SD),    stage, 'postinf_sd',    'posterior inflation sd')

end subroutine set_filename_info

!------------------------------------------------------------------

subroutine set_input_file_info( file_info, num_ens, STAGE_COPIES )

type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: num_ens
integer,              intent(in)    :: STAGE_COPIES(NUM_SCOPIES)

if ( perturb_from_single_instance ) then
   call set_io_copy_flag(file_info, STAGE_COPIES(MEM_START), READ_COPY)
   !>@todo know whether we are perturbing or not
   !#! call set_perturb_members(file_info, MEM_START, num_ens)
else
   call set_io_copy_flag(file_info, STAGE_COPIES(MEM_START), STAGE_COPIES(MEM_START)+num_ens-1, READ_COPY)
endif

if ( do_prior_inflate .and. prior_inflate_from_restart) then
   call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_MEAN), READ_COPY, inherit_units=.false.)
   call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_SD),   READ_COPY, inherit_units=.false.)
endif

if ( do_posterior_inflate .and. posterior_inflate_from_restart) then
   call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_MEAN),  READ_COPY, inherit_units=.false.)
   call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_SD),    READ_COPY, inherit_units=.false.)
endif

! This is for single file augmented state mean and sd if requested
if(single_file_in) then
   if (output_mean) then
     call set_io_copy_flag(file_info,    STAGE_COPIES(ENS_MEAN),  WRITE_COPY, inherit_units=.true.)

      if ( do_prior_inflate .and. .not. prior_inflate_from_restart ) &
        call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_MEAN), WRITE_COPY, inherit_units=.false.)

      if ( do_posterior_inflate .and. .not. posterior_inflate_from_restart ) &
        call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_MEAN),  WRITE_COPY, inherit_units=.false.)
   endif
   
   if (output_sd) then
     call set_io_copy_flag(file_info, STAGE_COPIES(ENS_SD),    WRITE_COPY, inherit_units=.true.)

      if ( do_prior_inflate .and. .not. prior_inflate_from_restart ) &
        call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_SD), WRITE_COPY, inherit_units=.false.)

      if ( do_posterior_inflate .and. .not. posterior_inflate_from_restart ) &
        call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_SD),  WRITE_COPY, inherit_units=.false.)
   endif
endif

end subroutine set_input_file_info

!------------------------------------------------------------------

subroutine set_output_file_info( file_info, num_ens, STAGE_COPIES, do_clamping, force_copy)

type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: num_ens
integer,              intent(in)    :: STAGE_COPIES(NUM_SCOPIES)
logical,              intent(in)    :: do_clamping
logical,              intent(in)    :: force_copy

!>@todo revisit if we should be clamping mean copy for file_info_output
if ( num_ens > 0 .and. output_members ) then
   call set_io_copy_flag(file_info, STAGE_COPIES(MEM_START), STAGE_COPIES(MEM_START)+num_ens-1, WRITE_COPY, &
                         num_output_ens=num_ens, clamp_vars=do_clamping, &
                         force_copy_back=force_copy)
endif

if ( output_mean )          &
   call set_io_copy_flag(file_info, STAGE_COPIES(ENS_MEAN),       WRITE_COPY, &
                         inherit_units=.true., clamp_vars=do_clamping, force_copy_back=force_copy)
if ( output_sd )            &
   call set_io_copy_flag(file_info, STAGE_COPIES(ENS_SD),         WRITE_COPY, &
                         inherit_units=.true., force_copy_back=force_copy)
if ( do_prior_inflate )     &
   call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_MEAN),  WRITE_COPY, &
                         inherit_units=.false., force_copy_back=force_copy)
if ( do_prior_inflate )     &
   call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_SD),    WRITE_COPY, &
                         inherit_units=.false., force_copy_back=force_copy)
if ( do_posterior_inflate ) &
   call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_MEAN),   WRITE_COPY, &
                         inherit_units=.false., force_copy_back=force_copy)
if ( do_posterior_inflate ) &
   call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_SD),     WRITE_COPY, &
                         inherit_units=.false., force_copy_back=force_copy)

end subroutine set_output_file_info

!-----------------------------------------------------------
!> checks the user input and informs the IO modules which files to write.


subroutine parse_stages_to_write(stages)

character(len=*), intent(in) :: stages(:)

integer :: nstages, i
character (len=32) :: my_stage

nstages = size(stages,1)

do i = 1, nstages
   my_stage = stages(i)
   call to_upper(my_stage)
   if (trim(my_stage) /= trim('NULL')) then
   SELECT CASE (my_stage)
      CASE ('INPUT', 'FORECAST', 'PREASSIM', 'POSTASSIM', 'ANALYSIS', 'OUTPUT')
         call set_stage_to_write(stages(i),.true.)
         write(msgstring,*)"filter will write stage : "//trim(stages(i))
         call error_handler(E_MSG,'parse_stages_to_write:',msgstring,source)
      CASE DEFAULT
         write(msgstring,*)"unknown stage : "//trim(stages(i))
         call error_handler(E_ERR,'parse_stages_to_write:',msgstring,source, &
                           text2="currently supported stages include :",&
                           text3="input, forecast, preassim, postassim, analysis, output")
   END SELECT

   endif
enddo

end subroutine parse_stages_to_write

!-----------------------------------------------------------
!> checks the user input and informs the IO modules which files to write.


function next_copy_number(cnum, ncopies)
integer, intent(inout)        :: cnum
integer, intent(in), optional :: ncopies
integer :: next_copy_number

if (present(ncopies)) then
   next_copy_number = cnum + ncopies - 1
else
   next_copy_number = cnum + 1
endif

cnum = next_copy_number

end function next_copy_number

!-----------------------------------------------------------
!> initialize file names and which copies should be read and or written


subroutine initialize_file_information(ncopies, &
                                       file_info_input,     file_info_mean_sd, &
                                       file_info_forecast,  file_info_preassim, &
                                       file_info_postassim, file_info_analysis, &
                                       file_info_output)

integer,              intent(in)  :: ncopies
type(file_info_type), intent(out) :: file_info_input
type(file_info_type), intent(out) :: file_info_mean_sd
type(file_info_type), intent(out) :: file_info_forecast
type(file_info_type), intent(out) :: file_info_preassim
type(file_info_type), intent(out) :: file_info_postassim
type(file_info_type), intent(out) :: file_info_analysis
type(file_info_type), intent(out) :: file_info_output

integer :: noutput_members, ninput_files, noutput_files, ndomains
character(len=256), allocatable :: file_array_input(:,:), file_array_output(:,:)

! local variable to shorten the name for function input
noutput_members = num_output_state_members 
ndomains        = get_num_domains()
noutput_files   = ens_size ! number of incomming ensemble members
ninput_files    = ens_size ! number of incomming ensemble members

! Assign the correct number of input and output files.
if (single_file_in .or. perturb_from_single_instance)  ninput_files = 1
if (single_file_out)                                  noutput_files = 1

! Given either a vector of in/output_state_files or a text file containing
! a list of files, return a vector of files containing the filenames.
call set_multiple_filename_lists(input_state_files(:), &
                                 input_state_file_list(:), &
                                 ndomains, &
                                 ninput_files,     &
                                 'filter','input_state_files','input_state_file_list')
call set_multiple_filename_lists(output_state_files(:), &
                                 output_state_file_list(:), &
                                 ndomains, &
                                 noutput_files, &
                                 'filter','output_state_files','output_state_file_list')

! Allocate space for file arrays.  contains a matrix of files (num_ens x num_domains)
! If perturbing from a single instance the number of input files does not have to
! be ens_size but rather a single file (or multiple files if more than one domain)
allocate(file_array_input(ninput_files, ndomains), file_array_output(noutput_files, ndomains))

file_array_input  = RESHAPE(input_state_files,  (/ninput_files,  ndomains/))
file_array_output = RESHAPE(output_state_files, (/noutput_files, ndomains/))


! Allocate space for the filename handles
call io_filenames_init(file_info_input,                       & 
                       ncopies       = ncopies,               &
                       cycling       = has_cycling,           &
                       single_file   = single_file_in,        &
                       restart_files = file_array_input,      &
                       root_name     = 'input')

! Output Files (we construct the filenames)
call io_filenames_init(file_info_mean_sd,   ncopies, has_cycling, single_file_out, root_name='input')
call io_filenames_init(file_info_forecast,  ncopies, has_cycling, single_file_out, root_name='forecast')
call io_filenames_init(file_info_preassim,  ncopies, has_cycling, single_file_out, root_name='preassim')
call io_filenames_init(file_info_postassim, ncopies, has_cycling, single_file_out, root_name='postassim')
call io_filenames_init(file_info_analysis,  ncopies, has_cycling, single_file_out, root_name='analysis')

! Write restart from output_state_file_list if provided
call io_filenames_init(file_info_output,                       &
                       ncopies       = ncopies,                &
                       cycling       = has_cycling,            &
                       single_file   = single_file_out,        &
                       restart_files = file_array_output,      &
                       root_name     = 'output',               &
                       check_output_compatibility = .true.)


! Set filename metadata information
!   Input Files
call set_filename_info(file_info_input,       'input',     ens_size,          CURRENT_COPIES )

!   Output Files
if (get_stage_to_write('input')) &
   call set_filename_info(file_info_mean_sd,  'input',     0,                   INPUT_COPIES )
if (get_stage_to_write('forecast')) &
   call set_filename_info(file_info_forecast, 'forecast',  noutput_members,  FORECAST_COPIES )
if (get_stage_to_write('preassim')) &
   call set_filename_info(file_info_preassim, 'preassim',  noutput_members,  PREASSIM_COPIES )
if (get_stage_to_write('postassim')) &
   call set_filename_info(file_info_postassim,'postassim', noutput_members, POSTASSIM_COPIES ) 
if (get_stage_to_write('analysis')) &
   call set_filename_info(file_info_analysis, 'analysis',  noutput_members,  ANALYSIS_COPIES )

call set_filename_info(file_info_output,      'output',    ens_size,          CURRENT_COPIES )

! Set file IO information
!   Input Files
call set_input_file_info( file_info_input, ens_size, CURRENT_COPIES ) 

!   Output Files
call set_output_file_info( file_info_mean_sd,           & 
                           num_ens      = 0,            &
                           STAGE_COPIES = INPUT_COPIES, &
                           do_clamping  = .false.,      &
                           force_copy   = .true. )

call set_output_file_info( file_info_forecast,             &
                           num_ens      = noutput_members, &
                           STAGE_COPIES = FORECAST_COPIES, &
                           do_clamping  = .false.,         &
                           force_copy   = .true. )

call set_output_file_info( file_info_preassim,             &
                           num_ens      = noutput_members, &
                           STAGE_COPIES = PREASSIM_COPIES, & 
                           do_clamping  = .false.,         &
                           force_copy   = .true. )

call set_output_file_info( file_info_postassim,             &
                           num_ens      = noutput_members,  &
                           STAGE_COPIES = POSTASSIM_COPIES, &
                           do_clamping  = .false.,          &
                           force_copy   = .true. )

call set_output_file_info( file_info_analysis,             &
                           num_ens      = noutput_members, &
                           STAGE_COPIES = ANALYSIS_COPIES, &
                           do_clamping  = .false.,         &
                           force_copy   = .true. )

call set_output_file_info( file_info_output,              &
                           num_ens      = ens_size,       &
                           STAGE_COPIES = CURRENT_COPIES, &
                           do_clamping  = .true.,         &
                           force_copy   = .false. )

end subroutine initialize_file_information


!-----------------------------------------------------------
!> set copy numbers. this is for when writing all stages at end

subroutine set_copies(cnum, STAGE_COPIES)
integer, intent(inout) :: cnum
integer, intent(inout) :: STAGE_COPIES(NUM_SCOPIES)

! Option to Output Postassim Ensemble Members Before Posterior Inflation
!   MEM_START
!   MEM_END = MEM_START + num_output_state_members - 1
STAGE_COPIES(MEM_START) = next_copy_number(cnum)
STAGE_COPIES(MEM_END)   = next_copy_number(cnum, num_output_state_members)

! Option to Output Input Mean and SD
!   MEAN
!   SD
if (output_mean) then
   STAGE_COPIES(ENS_MEAN) = next_copy_number(cnum)
endif
if (output_sd) then
   STAGE_COPIES(ENS_SD)   = next_copy_number(cnum)
endif

if (output_inflation) then
   ! Option to Output Infation with Damping
   !    PRIORINF_MEAN
   !    PRIORINF_SD
   !    POSTINF_MEAN
   !    POSTINF_SD
   if (do_prior_inflate) then
      STAGE_COPIES(PRIORINF_MEAN) = next_copy_number(cnum)
      STAGE_COPIES(PRIORINF_SD)   = next_copy_number(cnum)
   endif
   if (do_posterior_inflate) then
      STAGE_COPIES(POSTINF_MEAN)  = next_copy_number(cnum)
      STAGE_COPIES(POSTINF_SD)    = next_copy_number(cnum)
   endif
endif

end subroutine set_copies

!------------------------------------------------------------------

subroutine do_stage_output(stage_name, output_interval, time_step_number, &
   write_all_stages_at_end, state_ens_handle, COPIES, file_info, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

character(len = *),   intent(in)    :: stage_name
integer,              intent(in)    :: output_interval, time_step_number
logical,              intent(in)    :: write_all_stages_at_end
type(ensemble_type),  intent(inout) :: state_ens_handle
integer,              intent(inout) :: COPIES(:)
type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: ens_size, ENS_MEAN_COPY, ENS_SD_COPY

if(get_stage_to_write(stage_name)) then
   if((output_interval > 0) .and. &
      (time_step_number / output_interval * output_interval == time_step_number)) then

      ! Compute the ensemble mean and standard deviation
      ! For efficiency could make this optional in call
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      ! Save or output the data
      if (write_all_stages_at_end) then
         call store_copies(state_ens_handle, COPIES)
      else
         call write_state(state_ens_handle, file_info)
      endif
   endif
endif

end subroutine do_stage_output

!-------------------------------------------------------------------
end module filter_mod

