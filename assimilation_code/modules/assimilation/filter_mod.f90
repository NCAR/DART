! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module filter_mod

!------------------------------------------------------------------------------
use types_mod,             only : r8, i8, missing_r8, metadatalength, MAX_NUM_DOMS, MAX_FILES

use options_mod,           only : get_missing_ok_status, set_missing_ok_status

use obs_sequence_mod,      only : read_obs_seq, obs_type, obs_sequence_type,                  &
                                  get_obs_from_key, set_copy_meta_data, get_copy_meta_data,   &
                                  get_obs_def, get_time_range_keys, set_obs_values, set_obs,  &
                                  write_obs_seq, get_num_obs, get_obs_values, init_obs,       &
                                  assignment(=), get_num_copies, get_qc, get_num_qc, set_qc,  &
                                  static_init_obs_sequence, destroy_obs, read_obs_seq_header, &
                                  set_qc_meta_data, get_first_obs, get_obs_time_range,        &
                                  delete_obs_from_seq, delete_seq_head,                       &
                                  delete_seq_tail, replace_obs_values, replace_qc,            &
                                  destroy_obs_sequence, get_qc_meta_data, add_qc
                                 
use obs_def_mod,           only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                                  get_obs_def_type_of_obs

use obs_def_utilities_mod, only : set_debug_fwd_op

use time_manager_mod,      only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                  operator(-), print_time

use utilities_mod,         only : error_handler, E_ERR, E_MSG, E_DBG,                         &
                                  logfileunit, nmlfileunit, timestamp,                        &
                                  do_output, find_namelist_in_file, check_namelist_read,      &
                                  open_file, close_file, do_nml_file, do_nml_term, to_upper,  &
                                  set_multiple_filename_lists, find_textfile_dims

use assim_model_mod,       only : static_init_assim_model, get_model_size,                    &
                                  end_assim_model,  pert_model_copies, get_state_meta_data

use assim_tools_mod,       only : filter_assim, set_assim_tools_trace, test_state_copies
use obs_model_mod,         only : move_ahead, advance_state, set_obs_model_trace

use ensemble_manager_mod,  only : init_ensemble_manager, end_ensemble_manager,                &
                                  ensemble_type, get_copy, get_my_num_copies, put_copy,       &
                                  all_vars_to_all_copies, all_copies_to_all_vars,             &
                                  compute_copy_mean, compute_copy_mean_sd,                    &
                                  compute_copy_mean_var, duplicate_ens, get_copy_owner_index, &
                                  get_ensemble_time, set_ensemble_time, broadcast_copy,       &
                                  map_pe_to_task, copies_in_window, set_num_extra_copies,     &
                                  get_allow_transpose, all_copies_to_all_vars,                &
                                  allocate_single_copy, allocate_vars, get_single_copy,       &
                                  put_single_copy, deallocate_single_copy, print_ens_handle

use adaptive_inflate_mod,  only : do_ss_inflate, mean_from_restart, sd_from_restart,  &
                                  inflate_ens, adaptive_inflate_init,                 &
                                  adaptive_inflate_type, set_inflation_mean_copy ,    &
                                  log_inflation_info, set_inflation_sd_copy,          &
                                  get_minmax_task_zero, do_rtps_inflate,              &
                                  validate_inflate_options, PRIOR_INF, POSTERIOR_INF, &
                                  NO_INFLATION, OBS_INFLATION, VARYING_SS_INFLATION,  &
                                  SINGLE_SS_INFLATION, RELAXATION_TO_PRIOR_SPREAD,    &
                                  ENHANCED_SS_INFLATION

use mpi_utilities_mod,     only : my_task_id, task_sync, broadcast_send, broadcast_recv,      &
                                  task_count

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use state_vector_io_mod,   only : state_vector_io_init, read_state, write_state, &
                                  set_stage_to_write, get_stage_to_write

use io_filenames_mod,      only : io_filenames_init, file_info_type, &
                                  combine_file_info, set_file_metadata,  &
                                  set_member_file_metadata,  set_io_copy_flag, &
                                  check_file_info_variable_shape, &
                                  query_copy_present, COPY_NOT_PRESENT, &
                                  READ_COPY, WRITE_COPY, READ_WRITE_COPY

use direct_netcdf_mod,     only : finalize_single_file_io, write_augmented_state, &
                                  nc_get_num_times

use state_structure_mod,   only : get_num_domains

use forward_operator_mod,  only : get_obs_ens_distrib_state

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

! Defining whether diagnostics are for prior or posterior
integer, parameter :: PRIOR_DIAG = 0, POSTERIOR_DIAG = 2

! Determine if inflation it turned on or off for reading and writing
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

integer :: SPARE_PRIOR_SPREAD              = COPY_NOT_PRESENT

! Module Global Variables for inflation
logical :: do_prior_inflate     = .false.
logical :: do_posterior_inflate = .false.
type(adaptive_inflate_type) :: prior_inflate, post_inflate

logical :: has_cycling          = .false. ! filter will advance the model

! parms for trace/timing messages
integer, parameter :: T_BEFORE  = 1
integer, parameter :: T_AFTER   = 2
integer, parameter :: T_NEITHER = 3
logical, parameter :: P_TIME    = .true.

!----------------------------------------------------------------
! Namelist input with default values
!
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

! The inflation algorithm variables are defined in adaptive_inflate_mod.
! We use the integer parameters for PRIOR_INF and POSTERIOR_INF from 
! adaptive_inflate_mod to index these 'length 2' arrays.

integer  :: inf_flavor(2)                  = 0
logical  :: inf_initial_from_restart(2)    = .false.
logical  :: inf_sd_initial_from_restart(2) = .false.
logical  :: inf_deterministic(2)           = .true.
real(r8) :: inf_initial(2)                 = 1.0_r8
real(r8) :: inf_sd_initial(2)              = 0.0_r8
real(r8) :: inf_sd_max_change(2)           = 1.05_r8
real(r8) :: inf_damping(2)                 = 1.0_r8
real(r8) :: inf_lower_bound(2)             = 1.0_r8
real(r8) :: inf_upper_bound(2)             = 1000000.0_r8
real(r8) :: inf_sd_lower_bound(2)          = 0.0_r8

! Some models are allowed to have MISSING_R8 values in the DART state vector.
! If they are encountered, it is not necessarily a FATAL error.
! Most of the time, if a MISSING_R8 is encountered, DART should die.
! CLM should have allow_missing_clm = .true.
logical  :: allow_missing_clm = .false.


namelist /filter_nml/ async,     &
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
   inf_flavor,                   &
   inf_initial_from_restart,     &
   inf_sd_initial_from_restart,  &
   inf_sd_max_change,            & 
   inf_deterministic,            &
   inf_damping,                  &
   inf_initial,                  &
   inf_sd_initial,               &
   inf_lower_bound,              &
   inf_upper_bound,              &
   inf_sd_lower_bound,           &
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
type(time_type)             :: time1, first_obs_time, last_obs_time
type(time_type)             :: curr_ens_time, next_ens_time, window_time

integer,    allocatable :: keys(:)
integer(i8)             :: model_size
integer                 :: iunit, io, time_step_number, num_obs_in_set, ntimes
integer                 :: last_key_used, key_bounds(2)
integer                 :: in_obs_copy, obs_val_index
integer                 :: prior_obs_mean_index, posterior_obs_mean_index
integer                 :: prior_obs_spread_index, posterior_obs_spread_index
! Global indices into ensemble storage - observations
integer                 :: OBS_VAL_COPY, OBS_ERR_VAR_COPY, OBS_KEY_COPY
integer                 :: OBS_GLOBAL_QC_COPY,OBS_EXTRA_QC_COPY
integer                 :: OBS_MEAN_START, OBS_MEAN_END
integer                 :: OBS_VAR_START, OBS_VAR_END, TOTAL_OBS_COPIES
integer                 :: input_qc_index, DART_qc_index
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

real(r8), allocatable   :: prior_qc_copy(:)

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

call     trace_message('Filter start')
call timestamp_message('Filter start')

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

call trace_message('Before initializing inflation')

call validate_inflate_options(inf_flavor, inf_damping, inf_initial_from_restart, &
   inf_sd_initial_from_restart, inf_deterministic, inf_sd_max_change,            &
   do_prior_inflate, do_posterior_inflate, output_inflation, compute_posterior)

! Initialize the adaptive inflation module
call adaptive_inflate_init(prior_inflate, &
                           inf_flavor(PRIOR_INF), &
                           inf_initial_from_restart(PRIOR_INF), & 
                           inf_sd_initial_from_restart(PRIOR_INF), &
                           output_inflation, &
                           inf_deterministic(PRIOR_INF), &
                           inf_initial(PRIOR_INF), &
                           inf_sd_initial(PRIOR_INF), &
                           inf_lower_bound(PRIOR_INF), &
                           inf_upper_bound(PRIOR_INF), &
                           inf_sd_lower_bound(PRIOR_INF), &
                           inf_sd_max_change(PRIOR_INF), &
                           allow_missing, 'Prior')

call adaptive_inflate_init(post_inflate, &
                           inf_flavor(POSTERIOR_INF), &
                           inf_initial_from_restart(POSTERIOR_INF), &
                           inf_sd_initial_from_restart(POSTERIOR_INF), &
                           output_inflation, &
                           inf_deterministic(POSTERIOR_INF), &
                           inf_initial(POSTERIOR_INF), &
                           inf_sd_initial(POSTERIOR_INF), &
                           inf_lower_bound(POSTERIOR_INF), &
                           inf_upper_bound(POSTERIOR_INF), &
                           inf_sd_lower_bound(POSTERIOR_INF), &
                           inf_sd_max_change(POSTERIOR_INF), &
                           allow_missing, 'Posterior')

if (do_output()) then
   if (inf_flavor(PRIOR_INF) > NO_INFLATION .and. &
            inf_damping(PRIOR_INF) < 1.0_r8) then
      write(msgstring, '(A,F12.6,A)') 'Prior inflation damping of ', &
                                      inf_damping(PRIOR_INF), ' will be used'
      call error_handler(E_MSG,'filter_main:', msgstring)
   endif
   if (inf_flavor(POSTERIOR_INF) > NO_INFLATION .and. &
            inf_damping(POSTERIOR_INF) < 1.0_r8) then
      write(msgstring, '(A,F12.6,A)') 'Posterior inflation damping of ', &
                                      inf_damping(POSTERIOR_INF), ' will be used'
      call error_handler(E_MSG,'filter_main:', msgstring)
   endif
endif

call trace_message('After  initializing inflation')

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
num_state_ens_copies = count_state_ens_copies(ens_size, prior_inflate, post_inflate)
num_extras           = num_state_ens_copies - ens_size

! Observation
OBS_ERR_VAR_COPY     = ens_size + 1
OBS_VAL_COPY         = ens_size + 2
OBS_KEY_COPY         = ens_size + 3
OBS_GLOBAL_QC_COPY   = ens_size + 4
OBS_EXTRA_QC_COPY    = ens_size + 5
OBS_MEAN_START       = ens_size + 6
OBS_MEAN_END         = OBS_MEAN_START + num_groups - 1
OBS_VAR_START        = OBS_MEAN_START + num_groups
OBS_VAR_END          = OBS_VAR_START + num_groups - 1

TOTAL_OBS_COPIES = ens_size + 5 + 2*num_groups

!>@todo FIXME turn trace/timestamp calls into:  
!>
!> integer, parameter :: T_BEFORE = 1
!> integer, parameter :: T_AFTER  = 2
!> integer, parameter :: P_TIME   = 1
!>
!>  call progress(string, T_BEFORE)  ! simple trace msg
!>  call progress(string, T_AFTER)
!>
!>  call progress(string, T_BEFORE, P_TIME)  ! trace plus timestamp
!>  call progress(string, T_AFTER,  P_TIME)

!> DO NOT timestamp every trace message because some are
!> so quick that the timestamps don't impart any info.  
!> we should be careful to timestamp logical *sections* instead.

call     trace_message('Before setting up space for observations')
call timestamp_message('Before setting up space for observations')

! Initialize the obs_sequence; every pe gets a copy for now
call filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, input_qc_index, DART_qc_index, compute_posterior)

call timestamp_message('After  setting up space for observations')
call     trace_message('After  setting up space for observations')

call trace_message('Before setting up space for ensembles')

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

call trace_message('After  setting up space for ensembles')

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

call     trace_message('Before reading in ensemble restart files')
call timestamp_message('Before reading in ensemble restart files')

! for now, assume that we only allow cycling if single_file_out is true.
! code in this call needs to know how to initialize the output files.
call initialize_file_information(num_state_ens_copies ,                     &
                                 file_info_input      , file_info_mean_sd,  &
                                 file_info_forecast   , file_info_preassim, &
                                 file_info_postassim  , file_info_analysis, &
                                 file_info_output)

call check_file_info_variable_shape(file_info_output, state_ens_handle)

call set_inflation_mean_copy( prior_inflate, PRIOR_INF_COPY )
call set_inflation_sd_copy(   prior_inflate, PRIOR_INF_SD_COPY )
call set_inflation_mean_copy( post_inflate,  POST_INF_COPY )
call set_inflation_sd_copy(   post_inflate,  POST_INF_SD_COPY )

call read_state(state_ens_handle, file_info_input, read_time_from_file, time1, &
                prior_inflate, post_inflate, perturb_from_single_instance)

! This must be after read_state
call get_minmax_task_zero(prior_inflate, state_ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
call log_inflation_info(prior_inflate, state_ens_handle%my_pe, 'Prior', single_file_in)
call get_minmax_task_zero(post_inflate, state_ens_handle, POST_INF_COPY, POST_INF_SD_COPY)
call log_inflation_info(post_inflate, state_ens_handle%my_pe, 'Posterior', single_file_in)


if (perturb_from_single_instance) then
   call error_handler(E_MSG,'filter_main:', &
      'Reading in a single member and perturbing data for the other ensemble members')

   ! Only zero has the time, so broadcast the time to all other copy owners
   call broadcast_time_across_copy_owners(state_ens_handle, time1)
   call create_ensemble_from_single_file(state_ens_handle)
else
   call error_handler(E_MSG,'filter_main:', &
      'Reading in initial condition/restart data for all ensemble members from file(s)')
endif

call timestamp_message('After  reading in ensemble restart files')
call     trace_message('After  reading in ensemble restart files')

! see what our stance is on missing values in the state vector
allow_missing = get_missing_ok_status()

call     trace_message('Before initializing output files')
call timestamp_message('Before initializing output files')

! Initialize the output sequences and state files and set their meta data
call filter_generate_copy_meta_data(seq, in_obs_copy, &
      prior_obs_mean_index, posterior_obs_mean_index, &
      prior_obs_spread_index, posterior_obs_spread_index, &
      compute_posterior)

call timestamp_message('After  initializing output files')
call     trace_message('After  initializing output files')

call trace_message('Before trimming obs seq if start/stop time specified')

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

call trace_message('After  trimming obs seq if start/stop time specified')

! Time step number is used to do periodic diagnostic output
time_step_number = -1
curr_ens_time = set_time(0, 0)
next_ens_time = set_time(0, 0)
call filter_set_window_time(window_time)

! Compute mean and spread
call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

! Write out the mean and sd for the input files if requested
if (get_stage_to_write('input')) then

   call     trace_message('Before input state space output')
   call timestamp_message('Before input state space output')

   if (write_all_stages_at_end) then
      call store_input(state_ens_handle, prior_inflate, post_inflate)
   else
      ! if there is only one timestep in your input file insert the mean and sd if requested
      ntimes = nc_get_num_times(file_info_input%stage_metadata%filenames(1,1))
      if (single_file_out) then
        if ( ntimes == 1 ) then
           call write_augmented_state(state_ens_handle, file_info_input)
        else
           call error_handler(E_ERR,'filter_main', &
                   'can not insert mean or spread into input files that have multiple time steps',  &
                   source, text2='please remove "input" from stages_to_write')
        endif
     else ! muti file case
        ! write out input_mean.nc and input_sd.nc if requested
        call write_state(state_ens_handle, file_info_mean_sd)
     endif
   endif

   call timestamp_message('After  input state space output')
   call     trace_message('After  input state space output')

endif


AdvanceTime : do
   call trace_message('Top of main advance time loop')

   time_step_number = time_step_number + 1
   write(msgstring , '(A,I5)') &
      'Main assimilation loop, starting iteration', time_step_number
   call trace_message(' ', ' ', -1)
   call trace_message(msgstring, 'filter: ', -1)

   ! Check the time before doing the first model advance.  Not all tasks
   ! might have a time, so only check on PE0 if running multitask.
   ! This will get broadcast (along with the post-advance time) to all
   ! tasks so everyone has the same times, whether they have copies or not.
   ! If smoothing, we need to know whether the move_ahead actually advanced
   ! the model or not -- the first time through this loop the data timestamp
   ! may already include the first observation, and the model will not need
   ! to be run.  Also, last time through this loop, the move_ahead call
   ! will determine if there are no more obs, not call the model, and return
   ! with no keys in the list, which is how we know to exit.  In both of
   ! these cases, we must not advance the times on the lags.

   ! Figure out how far model needs to move data to make the window
   ! include the next available observation.  recent change is
   ! curr_ens_time in move_ahead() is intent(inout) and doesn't get changed
   ! even if there are no more obs.
   call trace_message('Before move_ahead checks time of data and next obs')

   call move_ahead(state_ens_handle, ens_size, seq, last_key_used, window_time, &
      key_bounds, num_obs_in_set, curr_ens_time, next_ens_time)

   call trace_message('After  move_ahead checks time of data and next obs')

   ! Only processes with an ensemble copy know to exit;
   ! For now, let process 0 broadcast its value of key_bounds
   ! This will synch the loop here and allow everybody to exit
   ! Need to clean up and have a broadcast that just sends a single integer???
   ! PAR For now, can only broadcast real arrays
   call filter_sync_keys_time(state_ens_handle, key_bounds, num_obs_in_set, &
                              curr_ens_time, next_ens_time)

   if(key_bounds(1) < 0) then
      call trace_message('No more obs to assimilate, exiting main loop', 'filter:', -1)
      exit AdvanceTime
   endif

   ! if model state data not at required time, advance model
   if (curr_ens_time /= next_ens_time) then

      ! we are going to advance the model - make sure we're doing single file output
      if (.not. has_cycling) then
         call error_handler(E_ERR,'filter:', &
             'advancing the model inside filter and multiple file output not currently supported', &
             source, text2='support will be added in subsequent releases', &
             text3='set "single_file_out=.true" for filter to advance the model, or advance the model outside filter')
      endif

      call trace_message('Ready to run model to advance data ahead in time', 'filter:', -1)
      call print_ens_time(state_ens_handle, 'Ensemble data time before advance')
      call     trace_message('Before running model')
      call timestamp_message('Before running model', sync=.true.)

      ! make sure storage is allocated in ensemble manager for vars.
      call allocate_vars(state_ens_handle)

      call all_copies_to_all_vars(state_ens_handle)

      call advance_state(state_ens_handle, ens_size, next_ens_time, async, &
              adv_ens_command, tasks_per_model_advance, file_info_output, file_info_input)

      call all_vars_to_all_copies(state_ens_handle)

      ! updated mean and spread after the model advance
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      ! update so curr time is accurate.
      curr_ens_time = next_ens_time
      state_ens_handle%current_time = curr_ens_time
      call set_time_on_extra_copies(state_ens_handle)

      ! only need to sync here since we want to wait for the
      ! slowest task to finish before outputting the time.
      call timestamp_message('After  running model', sync=.true.)
      call     trace_message('After  running model')
      call print_ens_time(state_ens_handle, 'Ensemble data time after  advance')
   else
      call trace_message('Model does not need to run; data already at required time', 'filter:', -1)
   endif

   call trace_message('Before setup for next group of observations')
   write(msgstring, '(A,I7)') 'Number of observations to be assimilated', &
      num_obs_in_set
   call trace_message(msgstring)
   call print_obs_time(seq, key_bounds(1), 'Time of first observation in window')
   call print_obs_time(seq, key_bounds(2), 'Time of last  observation in window')

   ! Create an ensemble for the observations from this time plus
   ! obs_error_variance, observed value, key from sequence, global qc,
   ! then mean for each group, then variance for each group
   call init_ensemble_manager(obs_fwd_op_ens_handle, TOTAL_OBS_COPIES, &
                              int(num_obs_in_set,i8), 1, transpose_type_in = 2)

   ! Also need a qc field for copy of each observation
   call init_ensemble_manager(qc_ens_handle, ens_size, &
                              int(num_obs_in_set,i8), 1, transpose_type_in = 2)

   ! Allocate storage for the keys for this number of observations
   allocate(keys(num_obs_in_set)) ! This is still var size for writing out the observation sequence

   ! Get all the keys associated with this set of observations
   ! Is there a way to distribute this?
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   call trace_message('After  setup for next group of observations')

   ! Write out forecast file(s). This contains the incoming ensemble members and potentially
   ! mean, sd, inflation values if requested.
   if (get_stage_to_write('forecast')) then
      if ((output_interval > 0) .and. &
          (time_step_number / output_interval * output_interval == time_step_number)) then

         call     trace_message('Before forecast state space output')
         call timestamp_message('Before forecast state space output')

         ! save or output the data
         if (write_all_stages_at_end) then
            call store_copies(state_ens_handle, FORECAST_COPIES)
         else
            call write_state(state_ens_handle, file_info_forecast)
         endif

         call timestamp_message('After  forecast state space output')
         call     trace_message('After  forecast state space output')

      endif
   endif

   if(do_ss_inflate(prior_inflate)) then
      call trace_message('Before prior inflation damping and prep')

      if (inf_damping(PRIOR_INF) /= 1.0_r8) then
         state_ens_handle%copies(PRIOR_INF_COPY, :) = 1.0_r8 + &
            inf_damping(PRIOR_INF) * (state_ens_handle%copies(PRIOR_INF_COPY, :) - 1.0_r8)
      endif

      call filter_ensemble_inflate(state_ens_handle, PRIOR_INF_COPY, prior_inflate, &
                                   ENS_MEAN_COPY)

      ! Recompute the the mean and spread as required for diagnostics
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      call trace_message('After  prior inflation damping and prep')
   endif

   ! if relaxation-to-prior-spread inflation, save the prior spread in SPARE_PRIOR_SPREAD
   if ( do_rtps_inflate(post_inflate) ) &
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, &
                                   SPARE_PRIOR_SPREAD)

   call     trace_message('Before computing prior observation values')
   call timestamp_message('Before computing prior observation values')

   ! Compute the ensemble of prior observations, load up the obs_err_var
   ! and obs_values. ens_size is the number of regular ensemble members,
   ! not the number of copies

   ! allocate() space for the prior qc copy
   call allocate_single_copy(obs_fwd_op_ens_handle, prior_qc_copy)

   call get_obs_ens_distrib_state(state_ens_handle, obs_fwd_op_ens_handle, &
           qc_ens_handle, seq, keys, obs_val_index, input_qc_index, &
           OBS_ERR_VAR_COPY, OBS_VAL_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
           OBS_EXTRA_QC_COPY, OBS_MEAN_START, OBS_VAR_START, &
           isprior=.true., prior_qc_copy=prior_qc_copy)

   call timestamp_message('After  computing prior observation values')
   call     trace_message('After  computing prior observation values')

   ! Write out preassim diagnostic files if requested.  This contains potentially 
   ! damped prior inflation values and the inflated ensemble.
   if (get_stage_to_write('preassim')) then
      if ((output_interval > 0) .and. &
          (time_step_number / output_interval * output_interval == time_step_number)) then

         call     trace_message('Before preassim state space output')
         call timestamp_message('Before preassim state space output')

         ! save or output the data
         if (write_all_stages_at_end) then
            call store_copies(state_ens_handle, PREASSIM_COPIES)
         else
            call write_state(state_ens_handle, file_info_preassim)
         endif

         call timestamp_message('After  preassim state space output')
         call     trace_message('After  preassim state space output')

      endif
   endif

   call trace_message('Before observation space diagnostics')

   ! This is where the mean obs
   ! copy ( + others ) is moved to task 0 so task 0 can update seq.
   ! There is a transpose (all_copies_to_all_vars(obs_fwd_op_ens_handle)) in obs_space_diagnostics
   ! Do prior observation space diagnostics and associated quality control
   call obs_space_diagnostics(obs_fwd_op_ens_handle, qc_ens_handle, ens_size, &
           seq, keys, PRIOR_DIAG, num_output_obs_members, in_obs_copy+1, &
           obs_val_index, OBS_KEY_COPY, &
           prior_obs_mean_index, prior_obs_spread_index, num_obs_in_set, &
           OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, &
           OBS_VAL_COPY, OBS_ERR_VAR_COPY, DART_qc_index, compute_posterior)
   call trace_message('After  observation space diagnostics')


   write(msgstring, '(A,I8,A)') 'Ready to assimilate up to', size(keys), ' observations'
   call trace_message(msgstring, 'filter:', -1)

   call     trace_message('Before observation assimilation')
   call timestamp_message('Before observation assimilation')

   call filter_assim(state_ens_handle, obs_fwd_op_ens_handle, seq, keys, &
      ens_size, num_groups, obs_val_index, prior_inflate, &
      ENS_MEAN_COPY, ENS_SD_COPY, &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
      OBS_MEAN_START, OBS_MEAN_END, OBS_VAR_START, &
      OBS_VAR_END, inflate_only = .false.)

   call timestamp_message('After  observation assimilation')
   call     trace_message('After  observation assimilation')

   ! Already transformed, so compute mean and spread for state diag as needed
   call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

   ! This block applies posterior inflation

   if(do_ss_inflate(post_inflate)) then

      call trace_message('Before posterior inflation damping')

      if (inf_damping(POSTERIOR_INF) /= 1.0_r8) then
         state_ens_handle%copies(POST_INF_COPY, :) = 1.0_r8 + &
            inf_damping(POSTERIOR_INF) * (state_ens_handle%copies(POST_INF_COPY, :) - 1.0_r8)
      endif

      call trace_message('After  posterior inflation damping')

   endif


   ! Write out postassim diagnostic files if requested.  This contains the assimilated ensemble 
   ! and potentially damped posterior inflation and updated prior inflation.
   if (get_stage_to_write('postassim')) then
      if ((output_interval > 0) .and. &
          (time_step_number / output_interval * output_interval == time_step_number)) then

         call     trace_message('Before postassim state space output')
         call timestamp_message('Before postassim state space output')

         ! save or output the data
         if (write_all_stages_at_end) then
            call store_copies(state_ens_handle, POSTASSIM_COPIES)
         else
            call write_state(state_ens_handle, file_info_postassim)
         endif

         call timestamp_message('After  postassim state space output')
         call     trace_message('After  postassim state space output')

      endif
   endif

   ! This block applies posterior inflation

   if(do_ss_inflate(post_inflate)) then

      call trace_message('Before posterior inflation applied to state')

      if (do_rtps_inflate(post_inflate)) then   
         call filter_ensemble_inflate(state_ens_handle, POST_INF_COPY, post_inflate, &
                       ENS_MEAN_COPY, SPARE_PRIOR_SPREAD, ENS_SD_COPY)
      else
         call filter_ensemble_inflate(state_ens_handle, POST_INF_COPY, post_inflate, &
                       ENS_MEAN_COPY)
      endif

      ! Recompute the mean or the mean and spread as required for diagnostics
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      call trace_message('After  posterior inflation applied to state')

   endif

   ! this block recomputes the expected obs values for the obs_seq.final file

   if (compute_posterior) then
      call     trace_message('Before computing posterior observation values')
      call timestamp_message('Before computing posterior observation values')
   
      ! Compute the ensemble of posterior observations, load up the obs_err_var
      ! and obs_values.  ens_size is the number of regular ensemble members,
      ! not the number of copies
   
       call get_obs_ens_distrib_state(state_ens_handle, obs_fwd_op_ens_handle, &
                qc_ens_handle, seq, keys, obs_val_index, input_qc_index, &
                OBS_ERR_VAR_COPY, OBS_VAL_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
                OBS_EXTRA_QC_COPY, OBS_MEAN_START, OBS_VAR_START, &
                isprior=.false., prior_qc_copy=prior_qc_copy)
   
      call deallocate_single_copy(obs_fwd_op_ens_handle, prior_qc_copy)
   
      call timestamp_message('After  computing posterior observation values')
      call     trace_message('After  computing posterior observation values')
   
   
      call trace_message('Before posterior obs space diagnostics')
   
      ! Write posterior observation space diagnostics
      ! There is a transpose (all_copies_to_all_vars(obs_fwd_op_ens_handle)) in obs_space_diagnostics
      call obs_space_diagnostics(obs_fwd_op_ens_handle, qc_ens_handle, ens_size, &
              seq, keys, POSTERIOR_DIAG, num_output_obs_members, in_obs_copy+2, &
              obs_val_index, OBS_KEY_COPY, &
              posterior_obs_mean_index, posterior_obs_spread_index, num_obs_in_set, &
              OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, &
              OBS_VAL_COPY, OBS_ERR_VAR_COPY, DART_qc_index, compute_posterior)
   
      call trace_message('After  posterior obs space diagnostics')
   else
      ! call this alternate routine to collect any updated QC values that may
      ! have been set in the assimilation loop and copy them to the outgoing obs seq
      call obs_space_sync_QCs(obs_fwd_op_ens_handle, seq, keys, num_obs_in_set, &
                              OBS_GLOBAL_QC_COPY, DART_qc_index)
      call deallocate_single_copy(obs_fwd_op_ens_handle, prior_qc_copy)
   endif

   ! this block computes the adaptive state space posterior inflation
   ! (it was applied earlier, this is computing the updated values for
   ! the next cycle.)

   ! CSS added condition: Don't update posterior inflation if relaxing to prior spread
   if(do_ss_inflate(post_inflate) .and. ( .not. do_rtps_inflate(post_inflate)) ) then

      ! If not reading the sd values from a restart file and the namelist initial
      !  sd < 0, then bypass this entire code block altogether for speed.
      if ((inf_sd_initial(POSTERIOR_INF) >= 0.0_r8) .or. &
           inf_sd_initial_from_restart(POSTERIOR_INF)) then

         call     trace_message('Before computing posterior state space inflation')
         call timestamp_message('Before computing posterior state space inflation')

         call filter_assim(state_ens_handle, obs_fwd_op_ens_handle, seq, keys, &
                 ens_size, num_groups, obs_val_index, post_inflate, &
                 ENS_MEAN_COPY, ENS_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY, &
                 OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, OBS_MEAN_START, OBS_MEAN_END, &
                 OBS_VAR_START, OBS_VAR_END, inflate_only = .true.)

         call timestamp_message('After  computing posterior state space inflation')
         call     trace_message('After  computing posterior state space inflation')

         ! recalculate standard deviation since this was overwritten in filter_assim
         call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)


      endif  ! sd >= 0 or sd from restart file
   endif  ! if doing state space posterior inflate

   ! Write out analysis diagnostic files if requested.  This contains the 
   ! posterior inflated ensemble and updated {prior,posterior} inflation values
   if (get_stage_to_write('analysis')) then
      if ((output_interval > 0) .and. &
          (time_step_number / output_interval * output_interval == time_step_number)) then

         call     trace_message('Before analysis state space output')
         call timestamp_message('Before analysis state space output')

         ! save or output the data
         if (write_all_stages_at_end) then
            call store_copies(state_ens_handle, ANALYSIS_COPIES)
         else
            call write_state(state_ens_handle, file_info_analysis)
         endif

         call timestamp_message('After  analysis state space output')
         call     trace_message('After  analysis state space output')

      endif
   endif

   ! only intended for debugging when cycling inside filter.
   ! writing the obs_seq file here will be slow - but if filter crashes
   ! you can get partial results by enabling this flag.
   if (write_obs_every_cycle) then
      call     trace_message('Before writing in-progress output sequence file')
      call timestamp_message('Before writing in-progress output sequence file')
      ! Only pe 0 outputs the observation space diagnostic file
      if(my_task_id() == 0) call write_obs_seq(seq, obs_sequence_out_name)
      call timestamp_message('After  writing in-progress output sequence file')
      call     trace_message('After  writing in-progress output sequence file')
   endif

   call trace_message('Near bottom of main loop, cleaning up obs space')
   ! Deallocate storage used for keys for each set
   deallocate(keys)

   ! The last key used is updated to move forward in the observation sequence
   last_key_used = key_bounds(2)

   ! Free up the obs ensemble space; LATER, can just keep it if obs are same size next time
   call end_ensemble_manager(obs_fwd_op_ens_handle)
   call end_ensemble_manager(qc_ens_handle)

   call trace_message('Bottom of main advance time loop')

end do AdvanceTime

call trace_message('End of main filter assimilation loop, starting cleanup', 'filter:', -1)

! Output the adjusted ensemble. If cycling only the last timestep is writen out
if (get_stage_to_write('output')) then
      call     trace_message('Before state space output')
      call timestamp_message('Before state space output')

      ! will write outside loop
      if (.not. write_all_stages_at_end) &
         call write_state(state_ens_handle, file_info_output)
   
      call timestamp_message('After  state space output')
      call     trace_message('After  state space output')

endif

call     trace_message('Before writing output sequence file')
call timestamp_message('Before writing output sequence file')
! Only pe 0 outputs the observation space diagnostic file
if(my_task_id() == 0) call write_obs_seq(seq, obs_sequence_out_name)
call timestamp_message('After  writing output sequence file')
call     trace_message('After  writing output sequence file')

! Output all restart files if requested
if (write_all_stages_at_end) then
   call     trace_message('Before writing all state restart files at end')
   call timestamp_message('Before writing all state restart files at end')

   file_info_all = combine_file_info( &
                    (/file_info_input, file_info_mean_sd, file_info_forecast, &
                      file_info_preassim, file_info_postassim, file_info_analysis, &
                      file_info_output/) )

   call write_state(state_ens_handle, file_info_all)

   call timestamp_message('After  writing all state restart files at end')
   call     trace_message('After  writing all state restart files at end')
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
call trace_message('Before end_model call')
call end_assim_model()
call trace_message('After  end_model call')

! deallocate qceff_table_data structures
call end_algorithm_info_mod()

call trace_message('Before ensemble and obs memory cleanup')
call end_ensemble_manager(state_ens_handle)

! Free up the obs sequence
call destroy_obs_sequence(seq)
call trace_message('After  ensemble and obs memory cleanup')

call     trace_message('Filter done')
call timestamp_message('Filter done')
if(my_task_id() == 0) then
   write(logfileunit,*)'FINISHED filter.'
   write(logfileunit,*)
endif

end subroutine filter_main

!-----------------------------------------------------------
!> This generates the copy meta data for the diagnostic files.
!> And also creates the state space diagnostic file.
!> Note for the state space diagnostic files the order of copies
!> in the diagnostic file is different from the order of copies
!> in the ensemble handle.
subroutine filter_generate_copy_meta_data(seq, in_obs_copy, &
   prior_obs_mean_index, posterior_obs_mean_index, &
   prior_obs_spread_index, posterior_obs_spread_index, &
   do_post)

type(obs_sequence_type),     intent(inout) :: seq
integer,                     intent(in)    :: in_obs_copy
integer,                     intent(out)   :: prior_obs_mean_index
integer,                     intent(out)   :: posterior_obs_mean_index
integer,                     intent(out)   :: prior_obs_spread_index
integer,                     intent(out)   :: posterior_obs_spread_index
logical,                     intent(in)    :: do_post

! Figures out the strings describing the output copies for the three output files.
! THese are the prior and posterior state output files and the observation sequence
! output file which contains both prior and posterior data.

character(len=metadatalength) :: prior_meta_data, posterior_meta_data
integer :: i, num_obs_copies

! only PE0 (here task 0) will allocate space for the obs_seq.final
!
! all other tasks should NOT allocate all this space.
! instead, set the copy numbers to an illegal value
! so we'll trap if they're used, and return early.
if (my_task_id() /= 0) then
   prior_obs_mean_index  = -1
   posterior_obs_mean_index = -1
   prior_obs_spread_index  = -1
   posterior_obs_spread_index = -1
   return
endif
  
! Set the metadata for the observations.

! Set up obs ensemble mean
num_obs_copies = in_obs_copy

num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
prior_obs_mean_index = num_obs_copies

if (do_post) then
   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble mean'
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   posterior_obs_mean_index = num_obs_copies
endif

! Set up obs ensemble spread
num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble spread'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
prior_obs_spread_index = num_obs_copies

if (do_post) then
   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble spread'
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   posterior_obs_spread_index = num_obs_copies
endif

! Make sure there are not too many copies requested - 
! proposed: make this magic number set in 1 place with an accessor
! routine so all parts of the code agree on max values.
if(num_output_obs_members > 10000) then
   write(msgstring, *)'output metadata in filter needs obs ensemble size < 10000, not ',&
                      num_output_obs_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source)
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

subroutine filter_initialize_modules_used()

call trace_message('Before filter_initialize_module_used call')

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
call state_vector_io_init()
call initialize_qc()

! Initialize algorothm_info_mod and read in QCF table data
call init_algorithm_info_mod()

call trace_message('After filter_initialize_module_used call')

end subroutine filter_initialize_modules_used

!-------------------------------------------------------------------------

subroutine filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, &
   input_qc_index, DART_qc_index, do_post)

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(out)   :: in_obs_copy, obs_val_index
integer,                 intent(out)   :: input_qc_index, DART_qc_index
logical,                 intent(in)    :: do_post

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
in_obs_copy = tnum_copies

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
input_qc_index = get_obs_qc_index(seq)
if (input_qc_index < 0) then
   input_qc_index = get_blank_qc_index(seq)
   if (input_qc_index < 0) then
      ! Need 1 new qc field for dummy incoming qc
      call add_qc(seq, 1)
      input_qc_index = get_blank_qc_index(seq)
      if (input_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen', source)
      endif
   endif
   ! Since we are constructing a dummy QC, label it as such
   call set_qc_meta_data(seq, input_qc_index, no_qc_meta_data)
endif

! check to be sure we either find an existing dart qc field and
! reuse it, or we add a new one. only on task 0.
DART_qc_index = get_obs_dartqc_index(seq)
if (DART_qc_index < 0 .and. my_task == io_task) then
   DART_qc_index = get_blank_qc_index(seq)
   if (DART_qc_index < 0) then
      ! Need 1 new qc field for the DART quality control
      call add_qc(seq, 1)
      DART_qc_index = get_blank_qc_index(seq)
      if (DART_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen', source)
      endif
   endif
   call set_qc_meta_data(seq, DART_qc_index, dqc_meta_data)
endif

! Determine which copy has actual obs value and return it.
obs_val_index = get_obs_copy_index(seq)

end subroutine filter_setup_obs_sequence

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

function get_obs_prior_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_prior_index

integer :: i

! Determine which copy in sequence has prior mean, if any.

do i = 1, get_num_copies(seq)
   get_obs_prior_index = i
   ! Need to look for 'prior mean'
   if(index(get_copy_meta_data(seq, i), 'prior ensemble mean') > 0) return
end do
! Falling of end means 'prior mean' not found; not fatal!

get_obs_prior_index = -1

end function get_obs_prior_index

!-------------------------------------------------------------------------

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

subroutine filter_ensemble_inflate(ens_handle, inflate_copy, inflate, ENS_MEAN_COPY, &
                                   SPARE_PRIOR_SPREAD, ENS_SD_COPY)

type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: inflate_copy, ENS_MEAN_COPY
type(adaptive_inflate_type), intent(inout) :: inflate
integer, optional,           intent(in)    :: SPARE_PRIOR_SPREAD, ENS_SD_COPY

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

   if ( do_rtps_inflate(inflate)) then 
      if ( present(SPARE_PRIOR_SPREAD) .and. present(ENS_SD_COPY)) then 
         write(msgstring, *) ' doing RTPS inflation'
         call error_handler(E_MSG,'filter_ensemble_inflate:',msgstring,source)

         !Reset the RTPS factor to the given input.nml value
         ens_handle%copies(inflate_copy, 1:ens_handle%my_num_vars) = inf_initial(POSTERIOR_INF)

         do j = 1, ens_handle%my_num_vars
            call inflate_ens(inflate, ens_handle%copies(grp_bot:grp_top, j), &
               ens_handle%copies(ENS_MEAN_COPY, j), ens_handle%copies(inflate_copy, j), 0.0_r8, &
               ens_handle%copies(SPARE_PRIOR_SPREAD, j), ens_handle%copies(ENS_SD_COPY, j)) 
         end do 
      else 
         write(msgstring, *) 'internal error: missing arguments for RTPS inflation, should not happen'
         call error_handler(E_ERR,'filter_ensemble_inflate',msgstring,source)
      endif 
   else 

      ! This is an initial test of doing inflation in probit space
      ! Note that this appears to work with adaptive inflation, but more research would be good
      ! Probably also shouldn't be used with groups for now although it is coded to do so
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
         call inflate_ens(inflate, probit_ens(1:grp_size), probit_ens_mean, &
            ens_handle%copies(inflate_copy, j))
         ! Transform back from probit space
         call transform_from_probit(grp_size, probit_ens(1:grp_size), &
            dist_params, ens_handle%copies(grp_bot:grp_top, j))
      end do
   endif
end do

end subroutine filter_ensemble_inflate

!-------------------------------------------------------------------------

subroutine obs_space_diagnostics(obs_fwd_op_ens_handle, qc_ens_handle, ens_size, &
   seq, keys, prior_post, num_output_members, members_index, &
   obs_val_index, OBS_KEY_COPY, &
   ens_mean_index, ens_spread_index, num_obs_in_set, &
   OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, OBS_VAL_COPY, &
   OBS_ERR_VAR_COPY, DART_qc_index, do_post)

! Do observation space diagnostics on the set of obs corresponding to keys

type(ensemble_type),     intent(inout) :: obs_fwd_op_ens_handle, qc_ens_handle
integer,                 intent(in)    :: ens_size
integer,                 intent(in)    :: num_obs_in_set
integer,                 intent(in)    :: keys(num_obs_in_set), prior_post
integer,                 intent(in)    :: num_output_members, members_index
integer,                 intent(in)    :: obs_val_index
integer,                 intent(in)    :: OBS_KEY_COPY
integer,                 intent(in)    :: ens_mean_index, ens_spread_index
type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: OBS_MEAN_START, OBS_VAR_START
integer,                 intent(in)    :: OBS_GLOBAL_QC_COPY, OBS_VAL_COPY
integer,                 intent(in)    :: OBS_ERR_VAR_COPY, DART_qc_index
logical,                 intent(in)    :: do_post

integer               :: j, k, ens_offset, copy_factor
integer               :: ivalue, io_task, my_task
real(r8), allocatable :: obs_temp(:)
real(r8)              :: rvalue(1)

! Do verbose forward operator output if requested
if(output_forward_op_errors) call verbose_forward_op_output(qc_ens_handle, prior_post, ens_size, keys)

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
call get_copy(io_task, obs_fwd_op_ens_handle, OBS_MEAN_START, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_obs_values(seq, keys(j), rvalue, ens_mean_index)
     end do
  endif

! Update the ensemble spread
call get_copy(io_task, obs_fwd_op_ens_handle, OBS_VAR_START, obs_temp)
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
call get_copy(io_task, obs_fwd_op_ens_handle, OBS_GLOBAL_QC_COPY, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_qc(seq, keys(j), rvalue, DART_qc_index)
   end do
endif

deallocate(obs_temp)

end subroutine obs_space_diagnostics

!-------------------------------------------------------------------------

subroutine obs_space_sync_QCs(obs_fwd_op_ens_handle,  &
   seq, keys, num_obs_in_set, OBS_GLOBAL_QC_COPY, DART_qc_index)


type(ensemble_type),     intent(inout) :: obs_fwd_op_ens_handle
integer,                 intent(in)    :: num_obs_in_set
integer,                 intent(in)    :: keys(num_obs_in_set)
type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: OBS_GLOBAL_QC_COPY
integer,                 intent(in)    :: DART_qc_index

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
call get_copy(io_task, obs_fwd_op_ens_handle, OBS_GLOBAL_QC_COPY, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_qc(seq, keys(j), rvalue, DART_qc_index)
   end do
endif

deallocate(obs_temp)

end subroutine obs_space_sync_QCs

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

subroutine trace_message(msg, label, threshold)

character(len=*), intent(in)           :: msg
character(len=*), intent(in), optional :: label
integer,          intent(in), optional :: threshold

! Write message to stdout and log file.
integer :: t

t = 0
if (present(threshold)) t = threshold

if (trace_level <= t) return

if (.not. do_output()) return

if (present(label)) then
   call error_handler(E_MSG,trim(label),trim(msg))
else
   call error_handler(E_MSG,' filter trace:',trim(msg))
endif

end subroutine trace_message

!-------------------------------------------------------------------------

subroutine timestamp_message(msg, sync)

character(len=*), intent(in) :: msg
logical, intent(in), optional :: sync

! Write current time and message to stdout and log file.
! if sync is present and true, sync mpi jobs before printing time.

if (timestamp_level <= 0) return

if (present(sync)) then
  if (sync) call task_sync()
endif

if (do_output()) call timestamp(' '//trim(msg), pos='brief')

end subroutine timestamp_message

!-------------------------------------------------------------------------
!>  call progress(string, T_BEFORE, P_TIME, label, threshold, sync)  ! trace plus timestamp
!-------------------------------------------------------------------------

subroutine progress(msg, when, dotime, label, threshold, sync)  ! trace plus timestamp

character(len=*), intent(in)           :: msg
integer,          intent(in)           :: when
logical,          intent(in)           :: dotime
character(len=*), intent(in), optional :: label
integer,          intent(in), optional :: threshold
logical,          intent(in), optional :: sync

! Write message to stdout and log file.
! optionally write timestamp.
integer :: t, lastchar
character(len=40) :: label_to_use

t = 0
if (present(threshold)) t = threshold

if (trace_level <= t) return

if (.not. do_output()) return

if (present(label)) then
   lastchar = min(len_trim(label), len(label_to_use))
   label_to_use = label(1:lastchar)
else
   label_to_use = ' filter_trace: '
endif

select case (when)
  case (T_BEFORE)
    call error_handler(E_MSG, trim(label_to_use)//' Before ', trim(msg))
  case (T_AFTER)
    call error_handler(E_MSG, trim(label_to_use)//' After  ', trim(msg))
  case default
    call error_handler(E_MSG, trim(label_to_use), trim(msg))
end select

if (timestamp_level <= 0) return

! if sync is present and true, sync mpi jobs before printing time.
if (present(sync)) then
  if (sync) call task_sync()
endif

if (do_output()) then
   select case (when)
     case (T_BEFORE)
      call timestamp(' Before '//trim(msg), pos='brief')
     case (T_AFTER)
      call timestamp(' After  '//trim(msg), pos='brief')
     case default
      call timestamp(' '//trim(msg), pos='brief')
   end select
endif

end subroutine progress

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

subroutine print_obs_time(seq, key, msg)

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key
character(len=*), intent(in), optional :: msg

! Write time of an observation to stdout and log file.
type(obs_type) :: obs
type(obs_def_type) :: obs_def
type(time_type) :: mtime

if (trace_level <= 0) return

if (do_output()) then
   call init_obs(obs, 0, 0)
   call get_obs_from_key(seq, key, obs)
   call get_obs_def(obs, obs_def)
   mtime = get_obs_def_time(obs_def)
   call print_time(mtime, ' filter trace: '//msg, logfileunit)
   call print_time(mtime, ' filter trace: '//msg)
   call destroy_obs(obs)
endif

end subroutine print_obs_time

!-------------------------------------------------------------------------
!> write out failed forward operators
!> This was part of obs_space_diagnostics

subroutine verbose_forward_op_output(qc_ens_handle, prior_post, ens_size, keys)

type(ensemble_type), intent(inout) :: qc_ens_handle
integer,             intent(in)    :: prior_post
integer,             intent(in)    :: ens_size
integer,             intent(in)    :: keys(:) ! I think this is still var size

character(len=12) :: task
integer :: j, i
integer :: forward_unit

write(task, '(i6.6)') my_task_id()

! all tasks open file?
if(prior_post == PRIOR_DIAG) then
   forward_unit = open_file('prior_forward_ope_errors' // task, 'formatted', 'append')
else
   forward_unit = open_file('post_forward_ope_errors' // task, 'formatted', 'append')
endif

! qc_ens_handle is a real representing an integer; values /= 0 get written out
do i = 1, ens_size
   do j = 1, qc_ens_handle%my_num_vars
      if(nint(qc_ens_handle%copies(i, j)) /= 0) write(forward_unit, *) i, j, keys(qc_ens_handle%my_vars(j)), nint(qc_ens_handle%copies(i, j))
   end do
end do

call close_file(forward_unit)

end subroutine verbose_forward_op_output

!------------------------------------------------------------------
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

subroutine store_input(ens_handle, prior_inflate, post_inflate)

type(ensemble_type),         intent(inout) :: ens_handle
type(adaptive_inflate_type), intent(in)    :: prior_inflate
type(adaptive_inflate_type), intent(in)    :: post_inflate

if( output_mean ) then
   if (query_copy_present( INPUT_COPIES(ENS_MEAN)) ) &
      ens_handle%copies(   INPUT_COPIES(ENS_MEAN), :) = ens_handle%copies(ENS_MEAN_COPY, :)

   if ( do_prior_inflate .and. .not. mean_from_restart( prior_inflate) ) then
     if (query_copy_present( INPUT_COPIES(PRIORINF_MEAN)) ) &
        ens_handle%copies(   INPUT_COPIES(PRIORINF_MEAN), :) = ens_handle%copies(PRIOR_INF_COPY, :)
   endif 

   if ( do_posterior_inflate .and. .not. mean_from_restart(post_inflate) ) then
      if (query_copy_present( INPUT_COPIES(POSTINF_MEAN)) ) &
         ens_handle%copies(   INPUT_COPIES(POSTINF_MEAN), :) = ens_handle%copies(POST_INF_COPY, :)
   endif
   
endif

if( output_sd ) then

   if (query_copy_present( INPUT_COPIES(ENS_SD)) ) then
      ens_handle%copies(   INPUT_COPIES(ENS_SD), :) = ens_handle%copies(ENS_SD_COPY, :)
   endif

   if ( do_prior_inflate .and. .not. sd_from_restart(prior_inflate) ) then
      if (query_copy_present( INPUT_COPIES(PRIORINF_SD)) ) then
         ens_handle%copies(   INPUT_COPIES(PRIORINF_SD), :) = ens_handle%copies(PRIOR_INF_SD_COPY, :)
      endif
   endif

   if ( do_posterior_inflate .and. .not. sd_from_restart(post_inflate) ) then
      if (query_copy_present( INPUT_COPIES(POSTINF_SD)) ) then
         ens_handle%copies(   INPUT_COPIES(POSTINF_SD), :)  = ens_handle%copies(POST_INF_SD_COPY, :)
      endif
   endif

endif

end subroutine store_input


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

function count_state_ens_copies(ens_size, post_inflate, prior_inflate) result(num_copies)

integer,                     intent(in) :: ens_size
type(adaptive_inflate_type), intent(in) :: prior_inflate
type(adaptive_inflate_type), intent(in) :: post_inflate
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
         if ( do_prior_inflate .and. .not. mean_from_restart(prior_inflate) ) then
            INPUT_COPIES(PRIORINF_MEAN) = next_copy_number(cnum)
         endif
         if ( do_posterior_inflate .and. .not. mean_from_restart(post_inflate) ) then
            INPUT_COPIES(POSTINF_MEAN)  = next_copy_number(cnum)
         endif
      endif

      if (output_sd) then
         INPUT_COPIES(ENS_SD) = next_copy_number(cnum)
         if ( do_prior_inflate .and. .not. sd_from_restart(prior_inflate) ) then
            INPUT_COPIES(PRIORINF_SD) = next_copy_number(cnum)
         endif
         if ( do_posterior_inflate .and. .not. sd_from_restart(post_inflate) ) then
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
   !   INPUT_MEAN
   !   INPUT_SD  
   if (output_mean) then
      INPUT_COPIES(ENS_MEAN) = ENS_MEAN_COPY
      if ( do_prior_inflate     .and. .not. mean_from_restart(prior_inflate) ) then
         INPUT_COPIES(PRIORINF_MEAN) = PRIOR_INF_COPY
      endif
      if ( do_posterior_inflate .and. .not. mean_from_restart(post_inflate) ) then
         INPUT_COPIES(POSTINF_MEAN)  = POST_INF_COPY
      endif
   endif

   if (output_sd) then
      INPUT_COPIES(ENS_SD) = ENS_SD_COPY
      if ( do_prior_inflate     .and. .not. sd_from_restart(prior_inflate) ) then
         INPUT_COPIES(PRIORINF_SD) = PRIOR_INF_SD_COPY
      endif
      if ( do_posterior_inflate .and. .not. sd_from_restart(post_inflate) ) then
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
if ( inf_flavor(POSTERIOR_INF) == RELAXATION_TO_PRIOR_SPREAD ) then
   SPARE_PRIOR_SPREAD = next_copy_number(cnum)
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

if ( do_prior_inflate ) then
   if ( inf_initial_from_restart(PRIOR_INF)    ) &
      call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_MEAN), READ_COPY, inherit_units=.false.)
   if ( inf_sd_initial_from_restart(PRIOR_INF) ) &
      call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_SD),   READ_COPY, inherit_units=.false.)
endif

if ( do_posterior_inflate ) then
   if ( inf_initial_from_restart(POSTERIOR_INF)    ) &
      call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_MEAN),  READ_COPY, inherit_units=.false.)
   if ( inf_sd_initial_from_restart(POSTERIOR_INF) ) &
      call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_SD),    READ_COPY, inherit_units=.false.)
endif

! This is for single file augmented state mean and sd if requested
if(single_file_in) then
   if (output_mean) then
     call set_io_copy_flag(file_info,    STAGE_COPIES(ENS_MEAN),  WRITE_COPY, inherit_units=.true.)

      if ( do_prior_inflate .and. .not. mean_from_restart(prior_inflate) ) &
        call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_MEAN), WRITE_COPY, inherit_units=.false.)

      if ( do_posterior_inflate .and. .not. mean_from_restart(post_inflate) ) &
        call set_io_copy_flag(file_info, STAGE_COPIES(POSTINF_MEAN),  WRITE_COPY, inherit_units=.false.)
   endif
   
   if (output_sd) then
     call set_io_copy_flag(file_info, STAGE_COPIES(ENS_SD),    WRITE_COPY, inherit_units=.true.)

      if ( do_prior_inflate .and. .not. sd_from_restart(prior_inflate) ) &
        call set_io_copy_flag(file_info, STAGE_COPIES(PRIORINF_SD), WRITE_COPY, inherit_units=.false.)

      if ( do_posterior_inflate .and. .not. sd_from_restart(post_inflate) ) &
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

!==================================================================
! TEST FUNCTIONS BELOW THIS POINT
!------------------------------------------------------------------
!> dump out obs_copies to file
subroutine test_obs_copies(obs_fwd_op_ens_handle, information)

type(ensemble_type), intent(in) :: obs_fwd_op_ens_handle
character(len=*),    intent(in) :: information

character(len=20)  :: task_str !! string to hold the task number
character(len=256) :: file_obscopies !! output file name
integer :: i, iunit

write(task_str, '(i10)') obs_fwd_op_ens_handle%my_pe
file_obscopies = TRIM('obscopies_' // TRIM(ADJUSTL(information)) // TRIM(ADJUSTL(task_str)))

iunit = open_file(file_obscopies, 'formatted', 'append')

do i = 1, obs_fwd_op_ens_handle%num_copies - 4
   write(iunit, *) obs_fwd_op_ens_handle%copies(i,:)
enddo

close(iunit)

end subroutine test_obs_copies

!-------------------------------------------------------------------
end module filter_mod

