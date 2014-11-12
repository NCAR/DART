! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> \dir filter  Main program contained here
!> \file filter.f90 Main program

program filter

!------------------------------------------------------------------------------
use types_mod,            only : r8, missing_r8, metadatalength
use obs_sequence_mod,     only : read_obs_seq, obs_type, obs_sequence_type,                  &
                                 get_obs_from_key, set_copy_meta_data, get_copy_meta_data,   &
                                 get_obs_def, get_time_range_keys, set_obs_values, set_obs,  &
                                 write_obs_seq, get_num_obs, get_obs_values, init_obs,       &
                                 assignment(=), get_num_copies, get_qc, get_num_qc, set_qc,  &
                                 static_init_obs_sequence, destroy_obs, read_obs_seq_header, &
                                 set_qc_meta_data, get_first_obs,          &
                                 get_obs_time_range, delete_obs_from_seq, delete_seq_head,   &
                                 delete_seq_tail, replace_obs_values, replace_qc,            &
                                 destroy_obs_sequence, get_qc_meta_data, add_qc,             &
                                 get_expected_obs_distrib_state !HK
use obs_def_mod,          only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                                 get_obs_kind
use time_manager_mod,     only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                 operator(-), print_time
use utilities_mod,        only : register_module,  error_handler, E_ERR, E_MSG, E_DBG,       &
                                 initialize_utilities, logfileunit, nmlfileunit, timestamp,  &
                                 do_output, find_namelist_in_file, check_namelist_read,      &
                                 open_file, close_file, do_nml_file, do_nml_term
use assim_model_mod,      only : static_init_assim_model, get_model_size,                    &
                                 netcdf_file_type, init_diag_output, finalize_diag_output,   & 
                                 aoutput_diagnostics, ens_mean_for_model, end_assim_model
use assim_tools_mod,      only : filter_assim, set_assim_tools_trace, get_missing_ok_status, &
                                 test_state_copies
use obs_model_mod,        only : move_ahead, advance_state, set_obs_model_trace
use ensemble_manager_mod, only : init_ensemble_manager, end_ensemble_manager,                &
                                 ensemble_type, get_copy, get_my_num_copies, put_copy,       &
                                 all_vars_to_all_copies, all_copies_to_all_vars,             &
                                 read_ensemble_restart, write_ensemble_restart,              &
                                 compute_copy_mean, compute_copy_mean_sd,                    &
                                 compute_copy_mean_var, duplicate_ens, get_copy_owner_index, &
                                 get_ensemble_time, set_ensemble_time, broadcast_copy,       &
                                 prepare_to_read_from_vars, prepare_to_write_to_vars, prepare_to_read_from_copies,    &
                                 prepare_to_write_to_copies, get_ensemble_time, set_ensemble_time,    &
                                 map_task_to_pe,  map_pe_to_task, prepare_to_update_copies,  &
                                 get_my_num_vars, allow_complete_state, &
                                 single_restart_file_in

use adaptive_inflate_mod, only : adaptive_inflate_end, do_varying_ss_inflate,                &
                                 do_single_ss_inflate, inflate_ens, adaptive_inflate_init,   &
                                 do_obs_inflate, adaptive_inflate_type,                      &
                                 output_inflate_diagnostics
use mpi_utilities_mod,    only : initialize_mpi_utilities, finalize_mpi_utilities,           &
                                 my_task_id, task_sync, broadcast_send, broadcast_recv,      &
                                 task_count
use smoother_mod,         only : smoother_read_restart, advance_smoother,                    &
                                 smoother_gen_copy_meta_data, smoother_write_restart,        &
                                 init_smoother, do_smoothing, smoother_mean_spread,          &
                                 smoother_assim, filter_state_space_diagnostics,             &
                                 smoother_ss_diagnostics, smoother_end, set_smoother_trace
use distributed_state_mod

use data_structure_mod, only : copies_in_window ! should this be through ensemble_manager?

use state_vector_io_mod,   only : read_transpose, transpose_write, get_state_variable_info,  &
                                  initialize_arrays_for_read, netcdf_filename, state_vector_io_init, &
                                  setup_read_write, turn_read_copy_on, turn_write_copy_on, turn_write_copies_off

use model_mod,            only : variables_domains, fill_variable_list, info_file_name, get_model_time

use io_filenames_mod,         only : io_filenames_init, restart_files_in, &
                                     query_diag_mean, query_diag_spread, query_diag_inf_mean, query_diag_inf_spread

use mpi


!------------------------------------------------------------------------------

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Some convenient global storage items
character(len=129)      :: msgstring
type(obs_type)          :: observation

integer                 :: trace_level, timestamp_level

! Defining whether diagnostics are for prior or posterior
integer, parameter :: PRIOR_DIAG = 0, POSTERIOR_DIAG = 2

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async = 0, ens_size = 20
logical  :: start_from_restart  = .false.
logical  :: output_restart      = .false.
logical  :: output_restart_mean = .false.
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
real(r8) :: outlier_threshold   = -1.0_r8
logical  :: enable_special_outlier_code = .false.
real(r8) :: input_qc_threshold  = 3.0_r8
logical  :: output_forward_op_errors = .false.
logical  :: output_timestamps        = .false.
logical  :: trace_execution          = .false.
logical  :: silence                  = .false.
logical  :: parallel_state_diag      = .true. ! default to write diagnostics in parallel - SKIPING at the moment, now stitching output together post processing
logical  :: direct_netcdf_read = .true. ! default to read from netcdf file

character(len = 129) :: obs_sequence_in_name  = "obs_seq.out",    &
                        obs_sequence_out_name = "obs_seq.final",  &
                        restart_in_file_name  = 'filter_ics',     &
                        restart_out_file_name = 'filter_restart', &
                        adv_ens_command       = './advance_model.csh'

!                  == './advance_model.csh'    -> advance ensemble using a script

! Inflation namelist entries follow, first entry for prior, second for posterior
! inf_flavor is 0:none, 1:obs space, 2: varying state space, 3: fixed state_space
integer              :: inf_flavor(2)             = 0
logical              :: inf_initial_from_restart(2)    = .false.
logical              :: inf_sd_initial_from_restart(2) = .false.
logical              :: inf_output_restart(2)     = .false.
logical              :: inf_deterministic(2)      = .true.
character(len = 129) :: inf_in_file_name(2)       = 'not_initialized',    &
                        inf_out_file_name(2)      = 'not_initialized',    &
                        inf_diag_file_name(2)     = 'not_initialized'
real(r8)             :: inf_initial(2)            = 1.0_r8
real(r8)             :: inf_sd_initial(2)         = 0.0_r8
real(r8)             :: inf_damping(2)            = 1.0_r8
real(r8)             :: inf_lower_bound(2)        = 1.0_r8
real(r8)             :: inf_upper_bound(2)        = 1000000.0_r8
real(r8)             :: inf_sd_lower_bound(2)     = 0.0_r8
logical              :: output_inflation          = .true.

namelist /filter_nml/ async, adv_ens_command, ens_size, tasks_per_model_advance,    &
   start_from_restart, output_restart, obs_sequence_in_name, obs_sequence_out_name, &
   restart_in_file_name, restart_out_file_name, init_time_days, init_time_seconds,  &
   first_obs_days, first_obs_seconds, last_obs_days, last_obs_seconds,              &
   obs_window_days, obs_window_seconds, enable_special_outlier_code,                &
   num_output_state_members, num_output_obs_members, output_restart_mean,           &
   output_interval, num_groups, outlier_threshold, trace_execution,                 &
   input_qc_threshold, output_forward_op_errors, output_timestamps,                 &
   inf_flavor, inf_initial_from_restart, inf_sd_initial_from_restart,               &
   inf_output_restart, inf_deterministic, inf_in_file_name, inf_damping,            &
   inf_out_file_name, inf_diag_file_name, inf_initial, inf_sd_initial,              &
   inf_lower_bound, inf_upper_bound, inf_sd_lower_bound, output_inflation,          &
   silence, parallel_state_diag, direct_netcdf_read


!----------------------------------------------------------------

! Doing this allows independent scoping for subroutines in main program file
call filter_main()

!----------------------------------------------------------------

contains 

!> allow_complete_state() queries the ensemble manager namelist value for
!> no_complete_state.
!> the code is distributed except:
!> task 0 still writes the obs_sequence file, so there is a transpose (copies to vars) and 
!> sending the obs_ens_handle%vars to task 0. Keys is also size obs%vars.
!>

subroutine filter_main()

type(ensemble_type)         :: ens_handle, obs_ens_handle, forward_op_ens_handle
type(obs_sequence_type)     :: seq
type(netcdf_file_type)      :: PriorStateUnit, PosteriorStateUnit
type(time_type)             :: time1, first_obs_time, last_obs_time
type(time_type)             :: curr_ens_time, next_ens_time, window_time
type(adaptive_inflate_type) :: prior_inflate, post_inflate

integer,    allocatable :: keys(:)
integer                 :: i, iunit, io, time_step_number, num_obs_in_set
integer                 :: ierr, last_key_used, model_size, key_bounds(2)
integer                 :: in_obs_copy, obs_val_index
integer                 :: output_state_mean_index, output_state_spread_index
integer                 :: prior_obs_mean_index, posterior_obs_mean_index
integer                 :: prior_obs_spread_index, posterior_obs_spread_index
! Global indices into ensemble storage - why are these in filter?
integer                 :: ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY, PRIOR_INF_SD_COPY
integer                 :: POST_INF_COPY, POST_INF_SD_COPY
! to avoid writing the prior diag
integer                 :: SPARE_COPY_MEAN, SPARE_COPY_SPREAD
integer                 :: SPARE_COPY_INF_MEAN, SPARE_COPY_INF_SPREAD
integer                 :: OBS_VAL_COPY, OBS_ERR_VAR_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer                 :: OBS_MEAN_START, OBS_MEAN_END
integer                 :: OBS_VAR_START, OBS_VAR_END, TOTAL_OBS_COPIES
integer                 :: input_qc_index, DART_qc_index
integer                 :: mean_owner, mean_owners_index
!HK
integer :: owner, owners_index
integer :: num_extras ! the extra ensemble copies

!HK 
doubleprecision start, finish ! for timing with MPI_WTIME

logical                 :: ds, all_gone, allow_missing

! HK
real(r8), allocatable   :: results(:,:)
integer                 :: ii, reps
real(r8), allocatable   :: temp_ens(:)
character*20 task_str, file_obscopies, file_results

!HK debug
logical :: write_flag

call filter_initialize_modules_used() ! static_init_model called in here

! Read the namelist entry
call find_namelist_in_file("input.nml", "filter_nml", iunit)
read(iunit, nml = filter_nml, iostat = io)
call check_namelist_read(iunit, io, "filter_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=filter_nml)
if (do_nml_term()) write(     *     , nml=filter_nml)

call set_trace(trace_execution, output_timestamps, silence)

call     trace_message('Filter start')
call timestamp_message('Filter start')

! Make sure ensemble size is at least 2 (NEED MANY OTHER CHECKS)
if(ens_size < 2) then
   write(msgstring, *) 'ens_size in namelist is ', ens_size, ': Must be > 1'
   call error_handler(E_ERR,'filter_main', msgstring, source, revision, revdate)
endif

! informational message to log
write(msgstring, '(A,I5)') 'running with an ensemble size of ', ens_size
call error_handler(E_MSG,'filter:', msgstring, source, revision, revdate)

! See if smoothing is turned on
ds = do_smoothing()

! Make sure inflation options are legal
do i = 1, 2
   if(inf_flavor(i) < 0 .or. inf_flavor(i) > 3) then
      write(msgstring, *) 'inf_flavor=', inf_flavor(i), ' Must be 0, 1, 2, 3 '
      call error_handler(E_ERR,'filter_main', msgstring, source, revision, revdate)
   endif
   if(inf_damping(i) < 0.0_r8 .or. inf_damping(i) > 1.0_r8) then
      write(msgstring, *) 'inf_damping=', inf_damping(i), ' Must be 0.0 <= d <= 1.0'
      call error_handler(E_ERR,'filter_main', msgstring, source, revision, revdate)
   endif
end do

! if doing something special with outlier threshold, say so
if (enable_special_outlier_code) then
   call error_handler(E_MSG,'filter:', 'special outlier threshold handling enabled', &
      source, revision, revdate)
endif

! Observation space inflation for posterior not currently supported
if(inf_flavor(2) == 1) call error_handler(E_ERR, 'filter_main', &
   'Posterior observation space inflation (type 1) not supported', source, revision, revdate)

! Setup the indices into the ensemble storage
num_extras = 10  ! six plus spare copies

! state
ENS_MEAN_COPY        = ens_size + 1
ENS_SD_COPY          = ens_size + 2
PRIOR_INF_COPY       = ens_size + 3
PRIOR_INF_SD_COPY    = ens_size + 4
POST_INF_COPY        = ens_size + 5
POST_INF_SD_COPY     = ens_size + 6
 ! Aim: to hang on to the prior_inf_copy which would have been written to the Prior_Diag.nc - and others if we need them
SPARE_COPY_MEAN       = ens_size + 7
SPARE_COPY_SPREAD     = ens_size + 8
SPARE_COPY_INF_MEAN   = ens_size + 9
SPARE_COPY_INF_SPREAD = ens_size + 10

! observation
OBS_ERR_VAR_COPY     = ens_size + 1
OBS_VAL_COPY         = ens_size + 2
OBS_KEY_COPY         = ens_size + 3
OBS_GLOBAL_QC_COPY   = ens_size + 4
OBS_MEAN_START       = ens_size + 5
OBS_MEAN_END         = OBS_MEAN_START + num_groups - 1
OBS_VAR_START        = OBS_MEAN_START + num_groups
OBS_VAR_END          = OBS_VAR_START + num_groups - 1

! Can't output more ensemble members than exist
if(num_output_state_members > ens_size) num_output_state_members = ens_size
if(num_output_obs_members   > ens_size) num_output_obs_members   = ens_size

call     trace_message('Before setting up space for observations')
call timestamp_message('Before setting up space for observations')

! Initialize the obs_sequence; every pe gets a copy for now
call filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, input_qc_index, DART_qc_index)

call timestamp_message('After  setting up space for observations')
call     trace_message('After  setting up space for observations')

call trace_message('Before setting up space for ensembles')

! Allocate model size storage and ens_size storage for metadata for outputting ensembles
model_size = get_model_size()

! set up ensemble HK WATCH OUT putting this here.
call init_ensemble_manager(ens_handle, ens_size + num_extras, model_size)

call trace_message('After  setting up space for ensembles')

! Don't currently support number of processes > model_size
if(task_count() > model_size) call error_handler(E_ERR,'filter_main', &
   'Number of processes > model size' ,source,revision,revdate)

call     trace_message('Before reading in ensemble restart files')
call timestamp_message('Before reading in ensemble restart files')

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time(time1)

! set up arrays for which copies to read/write
call setup_read_write(ens_size + num_extras)

! HK Moved initializing inflation to before read of restarts so you can read the restarts
! and inflation files in one step.

! Initialize the adaptive inflation module
! This activates turn_read_copy_on
call adaptive_inflate_init(ens_handle, prior_inflate, inf_flavor(1), inf_initial_from_restart(1), &
   inf_sd_initial_from_restart(1), inf_output_restart(1), inf_deterministic(1),       &
   inf_in_file_name(1), inf_out_file_name(1), inf_diag_file_name(1), inf_initial(1),  &
   inf_sd_initial(1), inf_lower_bound(1), inf_upper_bound(1), inf_sd_lower_bound(1),  &
   ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY, allow_missing, 'Prior',             &
   direct_netcdf_read)

call adaptive_inflate_init(ens_handle, post_inflate, inf_flavor(2), inf_initial_from_restart(2),  &
   inf_sd_initial_from_restart(2), inf_output_restart(2), inf_deterministic(2),       &
   inf_in_file_name(2), inf_out_file_name(2), inf_diag_file_name(2), inf_initial(2),  &
   inf_sd_initial(2), inf_lower_bound(2), inf_upper_bound(2), inf_sd_lower_bound(2),  &
   ens_handle, POST_INF_COPY, POST_INF_SD_COPY, allow_missing, 'Posterior',           &
   direct_netcdf_read)


if (do_output()) then
   if (inf_flavor(1) > 0 .and. inf_damping(1) < 1.0_r8) then
      write(msgstring, '(A,F12.6,A)') 'Prior inflation damping of ', inf_damping(1), ' will be used'
      call error_handler(E_MSG,'filter:', msgstring)
   endif
   if (inf_flavor(2) > 0 .and. inf_damping(2) < 1.0_r8) then
      write(msgstring, '(A,F12.6,A)') 'Posterior inflation damping of ', inf_damping(2), ' will be used'
      call error_handler(E_MSG,'filter:', msgstring)
   endif
endif

call trace_message('After  initializing inflation')

! Read in restart files and initialize the ensemble storage
call turn_read_copy_on(1, ens_size) ! need to read all restart copies

if (direct_netcdf_read) then
   call filter_read_restart_direct(ens_handle, time1, ens_size) 
else ! expecting DART restart files
   call filter_read_restart(ens_handle, time1, model_size)
endif

call test_state_copies(ens_handle, 'after_read')

! Read in or initialize smoother restarts as needed
if(ds) then
   call init_smoother(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)
   call smoother_read_restart(ens_handle, ens_size, model_size, time1, init_time_days)
endif

call timestamp_message('After  reading in ensemble restart files')
call     trace_message('After  reading in ensemble restart files')

! see what our stance is on missing values in the state vector
allow_missing = get_missing_ok_status()

call trace_message('Before initializing inflation')


call     trace_message('Before initializing output files')
call timestamp_message('Before initializing output files')

! Initialize the output sequences and state files and set their meta data
! HK this is where a parallel write of diagnostics differs from the traditional way 
! of having task 0 write all the diagnostics.

! Is there a problem if every task creates the meta data?
call filter_generate_copy_meta_data(seq, prior_inflate, &
      PriorStateUnit, PosteriorStateUnit, in_obs_copy, output_state_mean_index, &
      output_state_spread_index, prior_obs_mean_index, posterior_obs_mean_index, &
      prior_obs_spread_index, posterior_obs_spread_index)

if(ds) call error_handler(E_ERR, 'filter', 'smoother broken by Helen')
if(ds) call smoother_gen_copy_meta_data(num_output_state_members, output_inflation)

call timestamp_message('After  initializing output files')
call     trace_message('After  initializing output files')

! Let user know what settings the outlier threshold has, and data qc.
if (do_output()) then
   if (outlier_threshold <= 0.0_r8) then
      write(msgstring, '(A)') 'No observation outlier threshold rejection will be done'
   else
      write(msgstring, '(A,F12.6,A)') 'Will reject obs values more than', outlier_threshold, ' sigma from mean'
   endif
   call error_handler(E_MSG,'filter:', msgstring)

   write(msgstring, '(A,I4)') 'Will reject obs with Data QC larger than ', nint(input_qc_threshold)
   call error_handler(E_MSG,'filter:', msgstring)
endif

call trace_message('Before trimming obs seq if start/stop time specified')

! Need to find first obs with appropriate time, delete all earlier ones
if(first_obs_seconds > 0 .or. first_obs_days > 0) then
   first_obs_time = set_time(first_obs_seconds, first_obs_days)
   call delete_seq_head(first_obs_time, seq, all_gone)
   if(all_gone) then
      msgstring = 'All obs in sequence are before first_obs_days:first_obs_seconds'
      call error_handler(E_ERR,'filter_main',msgstring,source,revision,revdate)
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
      call error_handler(E_ERR,'filter_main',msgstring,source,revision,revdate)
   endif
endif

call trace_message('After  trimming obs seq if start/stop time specified')

! Time step number is used to do periodic diagnostic output
time_step_number = -1
curr_ens_time = set_time(0, 0)
next_ens_time = set_time(0, 0)
call filter_set_window_time(window_time)

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

   call move_ahead(ens_handle, ens_size, seq, last_key_used, window_time, &
      key_bounds, num_obs_in_set, curr_ens_time, next_ens_time)

   call trace_message('After  move_ahead checks time of data and next obs')

   ! Only processes with an ensemble copy know to exit;
   ! For now, let process 0 broadcast its value of key_bounds
   ! This will synch the loop here and allow everybody to exit
   ! Need to clean up and have a broadcast that just sends a single integer???
   ! PAR For now, can only broadcast real arrays
   call filter_sync_keys_time(ens_handle, key_bounds, num_obs_in_set, curr_ens_time, next_ens_time)
   if(key_bounds(1) < 0) then 
      call trace_message('No more obs to assimilate, exiting main loop', 'filter:', -1)
      exit AdvanceTime
   endif


   ! if model state data not at required time, advance model
   if (curr_ens_time /= next_ens_time) then
      ! Advance the lagged distribution, if needed.
      ! Must be done before the model runs and updates the data.
      if(ds) then
         call     trace_message('Before advancing smoother')
         call timestamp_message('Before advancing smoother')
         call advance_smoother(ens_handle)
         call timestamp_message('After  advancing smoother')
         call     trace_message('After  advancing smoother')
      endif

      call trace_message('Ready to run model to advance data ahead in time', 'filter:', -1)
      call print_ens_time(ens_handle, 'Ensemble data time before advance')
      call     trace_message('Before running model')
      call timestamp_message('Before running model', sync=.true.)
   
      call advance_state(ens_handle, ens_size, next_ens_time, async, &
                         adv_ens_command, tasks_per_model_advance)
   
      ! update so curr time is accurate.
      curr_ens_time = next_ens_time

      ! only need to sync here since we want to wait for the
      ! slowest task to finish before outputting the time.
      call timestamp_message('After  running model', sync=.true.)
      call     trace_message('After  running model')
      call print_ens_time(ens_handle, 'Ensemble data time after  advance')
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
   TOTAL_OBS_COPIES = ens_size + 4 + 2*num_groups
   call init_ensemble_manager(obs_ens_handle, TOTAL_OBS_COPIES, num_obs_in_set, 1)
   ! Also need a qc field for copy of each observation
   call init_ensemble_manager(forward_op_ens_handle, ens_size, num_obs_in_set, 1)

   ! Allocate storage for the keys for this number of observations
   allocate(keys(num_obs_in_set))

   ! Get all the keys associated with this set of observations
   ! Is there a way to distribute this?
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   call trace_message('After  setup for next group of observations')

   ! Compute mean and spread for inflation and state diagnostics
   call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

   if(do_single_ss_inflate(prior_inflate) .or. do_varying_ss_inflate(prior_inflate)) then
      call trace_message('Before prior inflation damping and prep')
      !call test_state_copies(ens_handle, 'before_prior_inflation')

      if (inf_damping(1) /= 1.0_r8) then
         call prepare_to_update_copies(ens_handle)
         ens_handle%copies(PRIOR_INF_COPY, :) = 1.0_r8 + &
            inf_damping(1) * (ens_handle%copies(PRIOR_INF_COPY, :) - 1.0_r8) 
      endif

      call filter_ensemble_inflate(ens_handle, PRIOR_INF_COPY, prior_inflate, ENS_MEAN_COPY)

      !call test_state_copies(ens_handle, 'after_prior_inflation')

      ! Recompute the the mean and spread as required for diagnostics
      call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      call trace_message('After  prior inflation damping and prep')
   endif

   call     trace_message('Before computing prior observation values')
   call timestamp_message('Before computing prior observation values')

   ! Compute the ensemble of prior observations, load up the obs_err_var
   ! and obs_values. ens_size is the number of regular ensemble members,
   ! not the number of copies
   start = MPI_WTIME()

   call get_obs_ens_distrib_state(ens_handle, obs_ens_handle, forward_op_ens_handle, &
     seq, keys, obs_val_index, input_qc_index, &
     OBS_ERR_VAR_COPY, OBS_VAL_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
     OBS_MEAN_START, OBS_VAR_START, isprior=.true.)

   finish = MPI_WTIME()

   if (my_task_id() == 0) print*, 'distributed average ', (finish-start)
   !call test_obs_copies(obs_ens_handle, 'prior')

!   goto 10011 !HK bail out after forward operators

   ! While we're here, make sure the timestamp on the actual ensemble copy
   ! for the mean has the current time.  If the user requests it be written
   ! out, it needs a valid timestamp.
   if (my_task_id() == 0 ) print*, '************ MEAN TIME *****************'
   call get_copy_owner_index(ENS_MEAN_COPY, mean_owner, mean_owners_index)
   if(ens_handle%my_pe == mean_owner) then
      ! Make sure the timestamp for the mean is the current time.
      call set_ensemble_time(ens_handle, mean_owners_index, curr_ens_time)
   endif

   call timestamp_message('After  computing prior observation values')
   call     trace_message('After  computing prior observation values')

   ! Do prior state space diagnostic output as required

!!*********************
! Diagnostic files.

   call trace_message('Before prior state space diagnostics')
   call timestamp_message('Before prior state space diagnostics')

   ! Store inflation mean copy in the spare copy. 
   ! The spare copy is left alone until the end
   ! shoving in four spare copies for now
   ens_handle%copies(SPARE_COPY_MEAN, :)       = ens_handle%copies(ENS_MEAN_COPY, :)
   ens_handle%copies(SPARE_COPY_SPREAD, :)     = ens_handle%copies(ENS_SD_COPY, :)
   ens_handle%copies(SPARE_COPY_INF_MEAN, :)   = ens_handle%copies(PRIOR_INF_COPY, :)
   ens_handle%copies(SPARE_COPY_INF_SPREAD, :) = ens_handle%copies(PRIOR_INF_SD_COPY, :)

   if ((output_interval > 0) .and. &
       (time_step_number / output_interval * output_interval == time_step_number)) then

       ! skeleton just to write the time to diagnostic files.
       call filter_state_space_diagnostics(curr_ens_time, PriorStateUnit, ens_handle, &
           model_size, num_output_state_members, &
           output_state_mean_index, output_state_spread_index, &
           output_inflation, ENS_MEAN_COPY, ENS_SD_COPY, &
           prior_inflate, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)

   endif

   call timestamp_message('After  prior state space diagnostics')
   call trace_message('After  prior state space diagnostics')

   call trace_message('Before observation space diagnostics')

   if (.not. allow_complete_state()) then ! task 0 still updating the sequence.
      allocate(obs_ens_handle%vars(obs_ens_handle%num_vars, obs_ens_handle%my_num_copies))
   endif

   ! This is where the mean obs
   ! copy ( + others ) is moved to task 0 so task 0 can update seq.
   ! There is a transpose (all_copies_to_all_vars(obs_ens_handle)) in obs_space_diagnostics
   ! Do prior observation space diagnostics and associated quality control
   call obs_space_diagnostics(obs_ens_handle, forward_op_ens_handle, ens_size, &
      seq, keys, PRIOR_DIAG, num_output_obs_members, in_obs_copy+1, &
      obs_val_index, OBS_KEY_COPY, &                                 ! new
      prior_obs_mean_index, prior_obs_spread_index, num_obs_in_set, &
      OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, &
      OBS_VAL_COPY, OBS_ERR_VAR_COPY, DART_qc_index)
   call trace_message('After  observation space diagnostics')

   if (.not. allow_complete_state()) then ! task 0 still updating the sequence.
      deallocate(obs_ens_handle%vars)
   endif

!*********************

   ! FIXME:  i believe both copies and vars are equal at the end
   ! of the obs_space diags, so we can skip this. 
   !call all_vars_to_all_copies(obs_ens_handle)

   write(msgstring, '(A,I8,A)') 'Ready to assimilate up to', size(keys), ' observations'
   call trace_message(msgstring, 'filter:', -1)

   call     trace_message('Before observation assimilation')
   call timestamp_message('Before observation assimilation')

   !call test_state_copies(ens_handle, 'before_filter_assim')

   call filter_assim(ens_handle, obs_ens_handle, seq, keys, &
      ens_size, num_groups, obs_val_index, prior_inflate, &
      ENS_MEAN_COPY, ENS_SD_COPY, &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
      OBS_MEAN_START, OBS_MEAN_END, OBS_VAR_START, &
      OBS_VAR_END, inflate_only = .false.)

   !call test_state_copies(ens_handle, 'after_filter_assim')

   call timestamp_message('After  observation assimilation')
   call     trace_message('After  observation assimilation')

   ! Do the update for the smoother lagged fields, too.
   ! Would be more efficient to do these all at once inside filter_assim 
   ! in the future
   if(ds) then
      write(msgstring, '(A,I8,A)') 'Ready to reassimilate up to', size(keys), ' observations in the smoother'
      call trace_message(msgstring, 'filter:', -1)

      call     trace_message('Before smoother assimilation')
      call timestamp_message('Before smoother assimilation')
      call smoother_assim(obs_ens_handle, seq, keys, ens_size, num_groups, &
         obs_val_index, ENS_MEAN_COPY, ENS_SD_COPY, &
         PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
         OBS_MEAN_START, OBS_MEAN_END, OBS_VAR_START, &
         OBS_VAR_END)
      call timestamp_message('After  smoother assimilation')
      call     trace_message('After  smoother assimilation')
   endif

   ! Already transformed, so compute mean and spread for state diag as needed
   call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

!-------- Test of posterior inflate ----------------

   if(do_single_ss_inflate(post_inflate) .or. do_varying_ss_inflate(post_inflate)) then

      call trace_message('Before posterior inflation damping and prep')
      !call test_state_copies(ens_handle, 'before_test_of_inflation')

      if (inf_damping(2) /= 1.0_r8) then
         call prepare_to_update_copies(ens_handle)
         ens_handle%copies(POST_INF_COPY, :) = 1.0_r8 + &
            inf_damping(2) * (ens_handle%copies(POST_INF_COPY, :) - 1.0_r8) 
      endif

    call filter_ensemble_inflate(ens_handle, POST_INF_COPY, post_inflate, ENS_MEAN_COPY)

    !call test_state_copies(ens_handle, 'after_test_of_inflation')

      ! Recompute the mean or the mean and spread as required for diagnostics
      call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      call trace_message('After  posterior inflation damping and prep')
   endif

!-------- End of posterior  inflate ----------------


   call     trace_message('Before computing posterior observation values')
   call timestamp_message('Before computing posterior observation values')

   ! Compute the ensemble of posterior observations, load up the obs_err_var 
   ! and obs_values.  ens_size is the number of regular ensemble members, 
   ! not the number of copies

    call get_obs_ens_distrib_state(ens_handle, obs_ens_handle, forward_op_ens_handle, &
     seq, keys, obs_val_index, input_qc_index, &
     OBS_ERR_VAR_COPY, OBS_VAL_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
     OBS_MEAN_START, OBS_VAR_START, isprior=.false.)

   !call test_obs_copies(obs_ens_handle, 'post')

   call timestamp_message('After  computing posterior observation values')
   call     trace_message('After  computing posterior observation values')

   if(ds) then
      call trace_message('Before computing smoother means/spread')
      call smoother_mean_spread(ens_size, ENS_MEAN_COPY, ENS_SD_COPY)
      call trace_message('After  computing smoother means/spread')
   endif

!***********************
!! Diagnostic files.

   call trace_message('Before posterior state space diagnostics')
   call timestamp_message('Before posterior state space diagnostics')

   ! Do posterior state space diagnostic output as required
   if ((output_interval > 0) .and. &
         (time_step_number / output_interval * output_interval == time_step_number)) then

      ! skeleton just to put time in the diagnostic file
      call filter_state_space_diagnostics(curr_ens_time, PosteriorStateUnit, ens_handle, &
         model_size, num_output_state_members, output_state_mean_index, &
         output_state_spread_index, &
         output_inflation, ENS_MEAN_COPY, ENS_SD_COPY, &
         post_inflate, POST_INF_COPY, POST_INF_SD_COPY)
      ! Cyclic storage for lags with most recent pointed to by smoother_head
      ! ens_mean is passed to avoid extra temp storage in diagnostics

      !> @todo What to do here?
      !call smoother_ss_diagnostics(model_size, num_output_state_members, &
       !  output_inflation, temp_ens, ENS_MEAN_COPY, ENS_SD_COPY, &
        ! POST_INF_COPY, POST_INF_SD_COPY)
   endif

   call timestamp_message('After  posterior state space diagnostics')
   call trace_message('After  posterior state space diagnostics')


   call trace_message('Before posterior obs space diagnostics')

   if (.not. allow_complete_state()) then ! task 0 still updating the sequence.
      allocate(obs_ens_handle%vars(obs_ens_handle%num_vars, obs_ens_handle%my_num_copies))
   endif

   ! Do posterior observation space diagnostics
   ! There is a transpose (all_copies_to_all_vars(obs_ens_handle)) in obs_space_diagnostics
   call obs_space_diagnostics(obs_ens_handle, forward_op_ens_handle, ens_size, &
      seq, keys, POSTERIOR_DIAG, num_output_obs_members, in_obs_copy+2, &
      obs_val_index, OBS_KEY_COPY, &                             ! new
      posterior_obs_mean_index, posterior_obs_spread_index, num_obs_in_set, &
      OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, &
      OBS_VAL_COPY, OBS_ERR_VAR_COPY, DART_qc_index)

   if (.not. allow_complete_state()) then ! task 0 still updating the sequence.
      deallocate(obs_ens_handle%vars)
   endif

!***********************

   call trace_message('After  posterior obs space diagnostics')

!-------- Test of posterior inflate ----------------
 
   if(do_single_ss_inflate(post_inflate) .or. do_varying_ss_inflate(post_inflate)) then

      ! If not reading the sd values from a restart file and the namelist initial
      !  sd < 0, then bypass this entire code block altogether for speed.
      if ((inf_sd_initial(2) >= 0.0_r8) .or. inf_sd_initial_from_restart(2)) then

         call     trace_message('Before computing posterior state space inflation')
         call timestamp_message('Before computing posterior state space inflation')

         call filter_assim(ens_handle, obs_ens_handle, seq, keys, ens_size, num_groups, &
            obs_val_index, post_inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
            POST_INF_COPY, POST_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
            OBS_MEAN_START, OBS_MEAN_END, OBS_VAR_START, &
            OBS_VAR_END, inflate_only = .true.)

         call timestamp_message('After  computing posterior state space inflation')
         call     trace_message('After  computing posterior state space inflation')

      endif  ! sd >= 0 or sd from restart file
   endif  ! if doing state space posterior inflate


!-------- End of posterior  inflate ----------------

   ! If observation space inflation, output the diagnostics
   if(do_obs_inflate(prior_inflate) .and. my_task_id() == 0) &
      call output_inflate_diagnostics(prior_inflate, curr_ens_time)

   call trace_message('Near bottom of main loop, cleaning up obs space')
   ! Deallocate storage used for keys for each set
   deallocate(keys)

   ! The last key used is updated to move forward in the observation sequence
   last_key_used = key_bounds(2)

   ! Free up the obs ensemble space; LATER, can just keep it if obs are same size next time
   call end_ensemble_manager(obs_ens_handle)
   call end_ensemble_manager(forward_op_ens_handle)

   call trace_message('Bottom of main advance time loop')
end do AdvanceTime

!call test_state_copies(ens_handle, 'last')

10011 continue

call trace_message('End of main filter assimilation loop, starting cleanup', 'filter:', -1)

call trace_message('Before finalizing diagnostics files')
! properly dispose of the diagnostics files
if(my_task_id() == 0) then
   ierr = finalize_diag_output(PriorStateUnit)
   ierr = finalize_diag_output(PosteriorStateUnit)
endif
call trace_message('After  finalizing diagnostics files')

call trace_message('Before writing output sequence file')
! Only pe 0 outputs the observation space diagnostic file
if(my_task_id() == 0) call write_obs_seq(seq, obs_sequence_out_name)
call trace_message('After  writing output sequence file')

call trace_message('Before writing inflation restart files if required')
call turn_write_copies_off(1, ens_size + num_extras) ! clean slate

! Output the restart for the adaptive inflation parameters
call adaptive_inflate_end(ens_handle, prior_inflate, ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY, direct_netcdf_read)
call adaptive_inflate_end(ens_handle, post_inflate, ens_handle, POST_INF_COPY, POST_INF_SD_COPY, direct_netcdf_read)
call trace_message('After  writing inflation restart files if required')

! Output a restart file if requested
call trace_message('Before writing state restart files if requested')
call timestamp_message('Before writing state restart files if requested')

call turn_write_copy_on(1,ens_size) ! restarts
! Prior_Diag copies - write spare copies
if (query_diag_mean()) call turn_write_copy_on(SPARE_COPY_MEAN)
if (query_diag_spread()) call turn_write_copy_on(SPARE_COPY_SPREAD)
if (query_diag_inf_mean()) call turn_write_copy_on(SPARE_COPY_INF_MEAN)
if (query_diag_inf_spread()) call turn_write_copy_on(SPARE_COPY_INF_SPREAD)

! Posterior Diag 
call turn_write_copy_on(ENS_MEAN_COPY) ! mean
call turn_write_copy_on(ENS_SD_COPY) ! sd
call turn_write_copy_on(POST_INF_COPY) ! posterior inf mean
call turn_write_copy_on(POST_INF_SD_COPY) ! posterior inf sd

if(direct_netcdf_read) then
   call filter_write_restart_direct(ens_handle)
else ! write
   call filter_write_restart(ens_handle)
endif

if(ds) call smoother_write_restart(1, ens_size)
call trace_message('After  writing state restart files if requested')
call timestamp_message('After  writing state restart files if requested')

! Give the model_mod code a chance to clean up. 
call trace_message('Before end_model call')
call end_assim_model()
call trace_message('After  end_model call')

call trace_message('Before ensemble and obs memory cleanup')
call end_ensemble_manager(ens_handle)

! Free up the observation kind and obs sequence
call destroy_obs(observation)
call destroy_obs_sequence(seq)
call trace_message('After  ensemble and obs memory cleanup')

if(ds) then 
   call trace_message('Before smoother memory cleanup')
   call smoother_end()
   call trace_message('After  smoother memory cleanup')
endif

call     trace_message('Filter done')
call timestamp_message('Filter done')
if(my_task_id() == 0) then 
   write(logfileunit,*)'FINISHED filter.'
   write(logfileunit,*)
endif

! YOU CAN NO LONGER WRITE TO THE LOG FILE BELOW THIS!
! After the call to finalize below, you cannot write to
! any fortran unit number.

! Make this the very last thing done, especially for SGI systems.
! It shuts down MPI and if you try to write after that, some libraries
! choose to discard output that is written after mpi is finalized, or 
! worse, the processes can hang.
call finalize_mpi_utilities(async=async)

end subroutine filter_main

!-----------------------------------------------------------

subroutine filter_generate_copy_meta_data(seq, prior_inflate, PriorStateUnit, &
   PosteriorStateUnit, in_obs_copy, output_state_mean_index, &
   output_state_spread_index, prior_obs_mean_index, posterior_obs_mean_index, &
   prior_obs_spread_index, posterior_obs_spread_index)

type(obs_sequence_type),     intent(inout) :: seq
type(adaptive_inflate_type), intent(in)    :: prior_inflate
type(netcdf_file_type),      intent(inout) :: PriorStateUnit, PosteriorStateUnit
integer,                     intent(out)   :: output_state_mean_index, output_state_spread_index
integer,                     intent(in)    :: in_obs_copy
integer,                     intent(out)   :: prior_obs_mean_index, posterior_obs_mean_index
integer,                     intent(out)   :: prior_obs_spread_index, posterior_obs_spread_index

! Figures out the strings describing the output copies for the three output files.
! THese are the prior and posterior state output files and the observation sequence
! output file which contains both prior and posterior data.

character(len=metadatalength) :: prior_meta_data, posterior_meta_data
! The 4 is for ensemble mean and spread plus inflation mean and spread
! The Prior file contains the prior inflation mean and spread only
! Posterior file contains the posterior inflation mean and spread only
character(len=metadatalength) :: state_meta(num_output_state_members + 4)
integer :: i, ensemble_offset, num_state_copies, num_obs_copies
integer :: ierr ! init_diag return code

! Section for state variables + other generated data stored with them.

! Ensemble mean goes first 
num_state_copies = num_output_state_members + 2
output_state_mean_index = 1
state_meta(output_state_mean_index) = 'ensemble mean'

! Ensemble spread goes second
output_state_spread_index = 2
state_meta(output_state_spread_index) = 'ensemble spread'

! Check for too many output ensemble members
if(num_output_state_members > 10000) then
   write(msgstring, *)'output metadata in filter needs state ensemble size < 10000, not ', &
                      num_output_state_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
endif

! Compute starting point for ensemble member output
ensemble_offset = 2

! Set up the metadata for the output state diagnostic files
do i = 1, num_output_state_members
   write(state_meta(i + ensemble_offset), '(a15, 1x, i6)') 'ensemble member', i
end do

! Next two slots are for inflation mean and sd metadata
! To avoid writing out inflation values to the Prior and Posterior netcdf files,
! set output_inflation to false in the filter section of input.nml 
if(output_inflation) then
   num_state_copies = num_state_copies + 2
   state_meta(num_state_copies-1) = 'inflation mean'
   state_meta(num_state_copies)   = 'inflation sd'
endif

! Have task 0 set up diagnostic output for model state, if output is desired
! I am not using a collective call here, just getting task 0 to set up the files
! - nc_write_model_atts.
if (my_task_id() == 0) then
   PriorStateUnit     = init_diag_output('Prior_Diag', &
                        'prior ensemble state', num_state_copies, state_meta)
   PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                        'posterior ensemble state', num_state_copies, state_meta)
endif

! Set the metadata for the observations.

! Set up obs ensemble mean
num_obs_copies = in_obs_copy
num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
prior_obs_mean_index = num_obs_copies
num_obs_copies = num_obs_copies + 1
posterior_meta_data = 'posterior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
posterior_obs_mean_index = num_obs_copies 

! Set up obs ensemble spread 
num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble spread'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
prior_obs_spread_index = num_obs_copies
num_obs_copies = num_obs_copies + 1
posterior_meta_data = 'posterior ensemble spread'
call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
posterior_obs_spread_index = num_obs_copies

! Make sure there are not too many copies requested
if(num_output_obs_members > 10000) then
   write(msgstring, *)'output metadata in filter needs obs ensemble size < 10000, not ',&
                      num_output_obs_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
endif

! Set up obs ensemble members as requested
do i = 1, num_output_obs_members
   write(prior_meta_data, '(a21, 1x, i6)') 'prior ensemble member', i
   write(posterior_meta_data, '(a25, 1x, i6)') 'posterior ensemble member', i
   num_obs_copies = num_obs_copies + 1
   call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   num_obs_copies = num_obs_copies + 1
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
end do


end subroutine filter_generate_copy_meta_data

!-------------------------------------------------------------------------

subroutine filter_initialize_modules_used()

! Initialize modules used that require it
call initialize_mpi_utilities('Filter')

call register_module(source,revision,revdate)

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call trace_message('Before init_model call')
call static_init_assim_model()
call trace_message('After  init_model call')
call state_vector_io_init()
call trace_message('After  init_stat_vector_io call')


end subroutine filter_initialize_modules_used

!-------------------------------------------------------------------------

subroutine filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, &
   input_qc_index, DART_qc_index)

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(out)   :: in_obs_copy, obs_val_index
integer,                 intent(out)   :: input_qc_index, DART_qc_index

character(len = metadatalength) :: no_qc_meta_data = 'No incoming data QC'
character(len = metadatalength) :: dqc_meta_data   = 'DART quality control'
character(len = 129) :: obs_seq_read_format
integer              :: obs_seq_file_id, num_obs_copies
integer              :: tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, qc_num_inc, num_qc
logical              :: pre_I_format

! Determine the number of output obs space fields
! 4 is for prior/posterior mean and spread, 
! Prior and posterior values for all selected fields (so times 2)
num_obs_copies = 2 * num_output_obs_members + 4

! Input file can have one qc field, none, or more.  note that read_obs_seq_header
! does NOT return the actual metadata values, which would be helpful in trying
! to decide if we need to add copies or qcs.
call read_obs_seq_header(obs_sequence_in_name, tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, &
   obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)


! if there are less than 2 incoming qc fields, we will need
! to make at least 2 (one for the dummy data qc and one for
! the dart qc).
if (tnum_qc < 2) then
   qc_num_inc = 2 - tnum_qc
else
   qc_num_inc = 0
endif

! Read in with enough space for diagnostic output values and add'l qc field(s)
call read_obs_seq(obs_sequence_in_name, num_obs_copies, qc_num_inc, 0, seq)

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
           'error adding blank qc field to sequence; should not happen', &
            source, revision, revdate)
      endif
   endif
   ! Since we are constructing a dummy QC, label it as such
   call set_qc_meta_data(seq, input_qc_index, no_qc_meta_data)
endif

! check to be sure we either find an existing dart qc field and
! reuse it, or we add a new one.
DART_qc_index = get_obs_dartqc_index(seq)
if (DART_qc_index < 0) then
   DART_qc_index = get_blank_qc_index(seq)
   if (DART_qc_index < 0) then
      ! Need 1 new qc field for the DART quality control
      call add_qc(seq, 1)
      DART_qc_index = get_blank_qc_index(seq)
      if (DART_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen', &
            source, revision, revdate)
      endif
   endif
   call set_qc_meta_data(seq, DART_qc_index, dqc_meta_data)
endif

! Get num of obs copies and num_qc
num_qc = get_num_qc(seq)
in_obs_copy = get_num_copies(seq) - num_obs_copies

! Create an observation type temporary for use in filter
call init_obs(observation, get_num_copies(seq), num_qc)

! Set initial DART quality control to 0 for all observations?
! Or leave them uninitialized, since
! obs_space_diagnostics should set them all without reading them

! Determine which copy has actual obs
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
   'Did not find observation copy with metadata "observation"', &
      source, revision, revdate)

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

subroutine filter_set_initial_time(time)

type(time_type), intent(out) :: time


if(init_time_days >= 0) then
   time = set_time(init_time_seconds, init_time_days)
else
   time = set_time(0, 0)
endif

end subroutine filter_set_initial_time

!-------------------------------------------------------------------------

subroutine filter_set_window_time(time)

type(time_type), intent(out) :: time


if(obs_window_days >= 0) then
   time = set_time(obs_window_seconds, obs_window_days)
else
   time = set_time(0, 0)
endif

end subroutine filter_set_window_time

!------------------------------------------------------------------
!> Read the restart information directly from the model output
!> netcdf file
!> Which routine should find model size?
!> 
subroutine filter_read_restart_direct(state_ens_handle, time, ens_size)!, model_size)

type(ensemble_type), intent(inout) :: state_ens_handle
type(time_type),     intent(inout) :: time
integer,             intent(in)    :: ens_size
!integer,             intent(out)   :: model_size

integer                         :: model_size
character(len=256), allocatable :: variable_list(:) !< does this need to be module storage
integer                         :: dart_index !< where to start in state_ens_handle%copies
integer                         :: num_domains !< number of input files to read
integer                         :: domain !< loop index
integer                         :: domain_size !< number of state elements in a domain
integer                         :: num_variables_in_state

! to start with, assume same variables in each domain - this will not always be the case
! If they are not in a domain, just set lengths to zero?
call variables_domains(num_variables_in_state, num_domains)
call io_filenames_init(ens_size, num_domains, inf_in_file_name, inf_out_file_name)

allocate(variable_list(num_variables_in_state))

variable_list = fill_variable_list(num_variables_in_state)

! need to know number of domains
call initialize_arrays_for_read(num_variables_in_state, num_domains) 

if (single_restart_file_in) then ! is this the correct flag to check?
   call get_state_variable_info(num_variables_in_state)
else
   model_size = 0
   do domain = 1, num_domains
      netcdf_filename = info_file_name(domain)
      call get_state_variable_info(num_variables_in_state, variable_list, domain, domain_size)
      model_size = model_size + domain_size
   enddo
endif

! read time from input file if time not set in namelist
if(init_time_days < 0) then
   time = get_model_time(restart_files_in(1,1)) ! Any of the restarts?
endif

state_ens_handle%time = time

! read in the data and transpose
dart_index = 1 ! where to start in state_ens_handle%copies - this is modified by read_transpose
do domain = 1, num_domains
   call read_transpose(state_ens_handle, restart_in_file_name, domain, dart_index)
enddo

deallocate(variable_list)

! Need Temporary print of initial model time?

end subroutine filter_read_restart_direct

!-------------------------------------------------------------------------

subroutine filter_read_restart(state_ens_handle, time, model_size)

type(ensemble_type), intent(inout) :: state_ens_handle
type(time_type),     intent(inout) :: time
integer,             intent(in)    :: model_size

integer :: days, secs

! need to allocate ensemble storage

if(my_task_id()==0) print*, 'regular dart restart files'

if (do_output()) then
   if (start_from_restart) then
      call error_handler(E_MSG,'filter_read_restart:', &
         'Reading in initial condition/restart data for all ensemble members from file(s)')
   else
      call error_handler(E_MSG,'filter_read_restart:', &
         'Reading in a single ensemble and perturbing data for the other ensemble members')
   endif
endif

! allocating storage space in ensemble manager
allocate(state_ens_handle%vars(state_ens_handle%num_vars, state_ens_handle%my_num_copies))

! Only read in initial conditions for actual ensemble members
if(init_time_days >= 0) then
   print*, 'reading ensemble restart task', my_task_id()
   call read_ensemble_restart(state_ens_handle, 1, ens_size, &
      start_from_restart, restart_in_file_name, time)
   if (do_output()) then
      call get_time(time, secs, days)
      write(msgstring, '(A)') 'By namelist control, ignoring time found in restart file.'
      call error_handler(E_MSG,'filter_read_restart:',msgstring,source,revision,revdate)
      write(msgstring, '(A,I6,1X,I5)') 'Setting initial days, seconds to ',days,secs
      call error_handler(E_MSG,'filter_read_restart:',msgstring,source,revision,revdate)
   endif
else
   call read_ensemble_restart(state_ens_handle, 1, ens_size, &
      start_from_restart, restart_in_file_name)
   if (state_ens_handle%my_num_copies > 0) time = state_ens_handle%time(1)
endif

! Temporary print of initial model time
if(state_ens_handle%my_pe == 0) then
   ! FIXME for the future: if pe 0 is not task 0, pe 0 can not print debug messages
   call get_time(time, secs, days)
   write(msgstring, *) 'initial model time of 1st ensemble member (days,seconds) ',days,secs
   call error_handler(E_DBG,'filter_read_restart',msgstring,source,revision,revdate)
endif

call all_vars_to_all_copies(state_ens_handle)

! deallocate whole state storage
deallocate(state_ens_handle%vars)

end subroutine filter_read_restart

!-------------------------------------------------------------------------
!> write the restart information into a DART restart file.
subroutine filter_write_restart(state_ens_handle)

type(ensemble_type) :: state_ens_handle

! allocating storage space in ensemble manager
allocate(state_ens_handle%vars(state_ens_handle%num_vars, state_ens_handle%my_num_copies))

call all_copies_to_all_vars(state_ens_handle)

call write_ensemble_restart(state_ens_handle, restart_out_file_name, 1, ens_size)

! deallocate whole state storage
deallocate(state_ens_handle%vars)

end subroutine filter_write_restart

!-------------------------------------------------------------------------
!> write the restart information directly into the model netcdf file.
subroutine filter_write_restart_direct(state_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle

integer :: dart_index !< where to start in state_ens_handle%copies
integer :: num_variables_in_state
integer :: num_domains
integer :: domain !< loop index

call variables_domains(num_variables_in_state, num_domains)

! transpose and write out the data
dart_index = 1
do domain = 1, num_domains
   call transpose_write(state_ens_handle, restart_out_file_name, domain, dart_index)
enddo

end subroutine filter_write_restart_direct

!-------------------------------------------------------------------------

subroutine filter_ensemble_inflate(ens_handle, inflate_copy, inflate, ENS_MEAN_COPY)

type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: inflate_copy, ENS_MEAN_COPY
type(adaptive_inflate_type), intent(inout) :: inflate

integer :: j, group, grp_bot, grp_top, grp_size

! Assumes that the ensemble is copy complete
call prepare_to_update_copies(ens_handle)

! Inflate each group separately;  Divide ensemble into num_groups groups
grp_size = ens_size / num_groups

do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1
   ! Compute the mean for this group
   call compute_copy_mean(ens_handle, grp_bot, grp_top, ENS_MEAN_COPY)

   do j = 1, ens_handle%my_num_vars
      call inflate_ens(inflate, ens_handle%copies(grp_bot:grp_top, j), &
         ens_handle%copies(ENS_MEAN_COPY, j), ens_handle%copies(inflate_copy, j))
   end do
end do

end subroutine filter_ensemble_inflate

!-------------------------------------------------------------------------


!> Computes forward observation operators and related quality control indicators.
!> @brief
!> Helen is working on this to use a distributed forward operator using MPI remote memeory
!> access
subroutine get_obs_ens_distrib_state(ens_handle, obs_ens_handle, forward_op_ens_handle, seq, keys, &
   obs_val_index, input_qc_index, &
   OBS_ERR_VAR_COPY, OBS_VAL_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, OBS_MEAN_START, OBS_VAR_START, isprior)

type(ensemble_type),     intent(in)    :: ens_handle
type(ensemble_type),     intent(inout) :: obs_ens_handle, forward_op_ens_handle 
type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: keys(:)
integer,                 intent(in)    :: obs_val_index, input_qc_index
integer,                 intent(in)    :: OBS_ERR_VAR_COPY, OBS_VAL_COPY
integer,                 intent(in)    :: OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer,                 intent(in)    :: OBS_MEAN_START, OBS_VAR_START !> @todo groups
logical,                 intent(in)    :: isprior

real(r8)             :: input_qc(1), obs_value(1), obs_err_var, thisvar(1)
integer              :: j, k, my_num_copies, global_ens_index, thiskey(1)
logical              :: evaluate_this_ob, assimilate_this_ob
type(obs_def_type)   :: obs_def
integer, allocatable :: istatus(:)

real(r8), allocatable                   :: expected_obs(:) !Also regular obs now?
integer global_obs_num
type(time_type)                         :: dummy_time
integer :: e
integer :: forward_min, forward_max !< for global qc
real(r8)              :: error, diff_sd, ratio
real(r8)              :: obs_prior_mean, obs_prior_var, obs_val
logical               :: do_outlier, good_forward_op, failed

! IMPORTANT, IT IS ASSUMED THAT ACTUAL ENSEMBLES COME FIRST
! HK: I think it is also assumed that the ensemble members are in the same order in
! each of the handles

allocate(istatus(copies_in_window(ens_handle))) 

! Loop through my copies and compute expected value
my_num_copies = get_my_num_copies(obs_ens_handle)

call prepare_to_write_to_vars(obs_ens_handle)
call prepare_to_write_to_vars(forward_op_ens_handle)
call prepare_to_read_from_vars(ens_handle)

! create the mpi window for the distributed state
call create_state_window(ens_handle)

! make some room for state vectors
allocate(expected_obs(copies_in_window(ens_handle)))

! Loop through all observations in the set
ALL_OBSERVATIONS: do j = 1, obs_ens_handle%my_num_vars

   global_obs_num = obs_ens_handle%my_vars(j) ! convert the local obs number to global obs number
   thiskey(1) = keys(global_obs_num)

   ! Get the information on this observation by placing it in temporary
   call get_obs_from_key(seq, keys(global_obs_num), observation)
   call get_obs_def(observation, obs_def)
   ! Check to see if this observation fails input qc test
   call get_qc(observation, input_qc, input_qc_index)
   ! If it is bad, set forward operator status value to -99 and return missing_r8 for obs_value

   ! PAR THIS SUBROUTINE SHOULD EVENTUALLY GO IN THE QUALITY CONTROL MODULE
   if(.not. input_qc_ok(input_qc(1), input_qc_threshold)) then
      ! The forward operator value is set to -99 if prior qc was failed
      forward_op_ens_handle%copies(:, j) = -99 
      obs_ens_handle%copies(OBS_KEY_COPY, j) = thiskey(1)

      ! No need to do anything else for a failed observation
      cycle ALL_OBSERVATIONS
   endif

   ! Get the observation value and error variance
   call get_obs_values(observation, obs_value(1:1), obs_val_index)
   obs_err_var = get_obs_def_error_variance(obs_def)

   global_ens_index = 1 ! HK where is this used?

   ! temporaries to avoid passing array sections which was slow on PGI compiler
   call get_expected_obs_distrib_state(seq, thiskey, &
     global_ens_index, dummy_time, isprior, &
     istatus, assimilate_this_ob, evaluate_this_ob, ens_handle, expected_obs)

    obs_ens_handle%copies(1:copies_in_window(ens_handle), j) = expected_obs

   ! If istatus is 0 (successful) then put 0 for assimilate, -1 for evaluate only
   ! and -2 for neither evaluate or assimilate. Otherwise pass through the istatus
   ! in the forward operator evaluation field

   do e = 1, copies_in_window(ens_handle) !>@todo this won't always be 6 (groups)

      if(istatus(e) == 0) then
         if ((assimilate_this_ob .or. evaluate_this_ob) .and. (expected_obs(e) == missing_r8))  then
               print*, ' observation ', global_obs_num
               write(msgstring, *) 'istatus was 0 (OK) but forward operator returned missing value.'
            call error_handler(E_ERR,'filter_main', msgstring, source, revision, revdate)
         endif
         if(assimilate_this_ob) then
            forward_op_ens_handle%copies(e, j) = 0
         else if(evaluate_this_ob) then
            forward_op_ens_handle%copies(e, j) = -1
         else
            forward_op_ens_handle%copies(e, j) = -2
         endif
      else if (istatus(e) < 0) then
         write(msgstring, *) 'istatus must not be <0 from forward operator. 0=OK, >0 for error'
         call error_handler(E_ERR,'filter_main', msgstring, source, revision, revdate)
      else
         forward_op_ens_handle%copies(e, j) = istatus(e)
      endif

   enddo

   ! update copy for error variance and for oberved value
   obs_ens_handle%copies(OBS_ERR_VAR_COPY, j) = obs_err_var
   obs_ens_handle%copies(OBS_VAL_COPY, j) = obs_value(1)
   obs_ens_handle%copies(OBS_KEY_COPY, j) = thiskey(1)

end do ALL_OBSERVATIONS

!> @todo - don't you have the mean already?
call compute_copy_mean_var(obs_ens_handle, 1, copies_in_window(ens_handle), OBS_MEAN_START, OBS_VAR_START)

!do_outlier = (prior_post == PRIOR_DIAG .and. outlier_threshold > 0.0_r8)
do_outlier = (isprior .and. outlier_threshold > 0.0_r8)

QC_LOOP: do j = 1, obs_ens_handle%my_num_vars

   good_forward_op = .false.

   ! compute outlier test and consolidate forward operator qc
   ! find the min and max istatus values across all ensemble members.  these are
   ! either set by dart code, or returned by the model-specific model_interpolate() 
   ! routine, or by forward operator code in obs_def_xxx_mod files.
   forward_max = nint(maxval(forward_op_ens_handle%copies(1:copies_in_window(ens_handle), j)))
   forward_min = nint(minval(forward_op_ens_handle%copies(1:copies_in_window(ens_handle), j)))
!--- copied from obs_state_diagnostics ---
   ! Now do a case statement to figure out what the qc result should be
   ! For prior, have to test for a bunch of stuff
   ! FIXME: note that this case statement doesn't cover every possibility;
   ! if there's an error in the code and a minus value gets into the forward
   ! operator istatus without being caught, it will fail all cases below.
   ! add another line for 'internal inconsistency' to be safe.
   if(isprior) then !HK changed from
      if(forward_min == -99) then              ! Failed prior qc in get_obs_ens
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 6
      else if(forward_min == -2) then          ! Observation not used via namelist
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 5
      else if(forward_max > 0) then            ! At least one forward operator failed
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 4
      else if(forward_min == -1) then          ! Observation to be evaluated only
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 1
         good_forward_op = .true.
      else if(forward_min == 0) then           ! All clear, assimilate this ob
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 0
         good_forward_op = .true.
      ! FIXME: proposed enhancement - catchall for cases that we have not caught
      !else   ! 'should not happen'
      !   obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 9  ! inconsistent istatus codes
      endif
        
      ! PAR: THIS SHOULD BE IN QC MODULE 
      ! Check on the outlier threshold quality control: move to QC module?
      ! bug fix: this was using the incoming qc threshold before.
      ! it should only be doing the outlier test on dart qc values of 0 or 1.
      ! if there is already a different qc code set, leave it alone.
      ! only if it is still successful (assim or eval, 0 or 1), then check
      ! for failing outlier test.
      if(do_outlier .and. (obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) < 2)) then
         obs_prior_mean = obs_ens_handle%copies(OBS_MEAN_START, j)
         obs_prior_var = obs_ens_handle%copies(OBS_VAR_START, j)
         obs_val = obs_ens_handle%copies(OBS_VAL_COPY, j)
         obs_err_var = obs_ens_handle%copies(OBS_ERR_VAR_COPY, j)
         error = obs_prior_mean - obs_val
         diff_sd = sqrt(obs_prior_var + obs_err_var)
         if (diff_sd /= 0.0_r8) then
            ratio = abs(error / diff_sd)
         else
            ratio = 0.0_r8
         endif

         ! if special handling requested, pass in the outlier ratio for this obs,
         ! the default outlier threshold value, and enough info to extract the specific 
         ! obs type for this obs.
         ! the function should return .true. if this is an outlier, .false. if it is ok.
         if (enable_special_outlier_code) then
            failed = failed_outlier(ratio, outlier_threshold, obs_ens_handle, &
                                    OBS_KEY_COPY, j, seq)
         else 
            failed = (ratio > outlier_threshold)
         endif

         if (failed) then
            obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 7
         endif
      endif

   else
      ! For failed posterior, only update qc if prior successful
      if(forward_max > 0) then
         ! Both the following 2 tests and assignments were on single executable lines,
         ! but one compiler (gfortran) was confused by this, so they were put in
         ! if/endif blocks.
         if(obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) == 0) then
            obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 2
         endif
         if(obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) == 1) then
            obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 3
         endif
      endif
      ! and for consistency, go through all the same tests as the prior, but
      ! only use the results to set the good/bad forward op flag.
      if ((forward_min == -99) .or. &       ! Failed prior qc in get_obs_ens
          (forward_min == -2)  .or. &       ! Observation not used via namelist
          (forward_max > 0)) then           ! At least one forward operator failed
            continue; 
      else if(forward_min == -1) then       ! Observation to be evaluated only
         good_forward_op = .true.
      else if(forward_min == 0) then        ! All clear, assimilate this ob
         good_forward_op = .true.
      endif
   endif

   ! for either prior or posterior, if the forward operator failed,
   ! reset the mean/var to missing_r8, regardless of the DART QC status
   ! HK does this fail if you have groups?
   if (.not. good_forward_op) then
      obs_ens_handle%copies(OBS_MEAN_START, j) = missing_r8
      obs_ens_handle%copies(OBS_VAR_START,  j) = missing_r8
   endif

! PAR: NEED TO BE ABLE TO OUTPUT DETAILS OF FAILED FORWARD OBSEVATION OPERATORS
! OR FAILED OUTLIER, ETC. NEEDS TO BE DONE BY PE 0 ONLY. PUT IT HERE FOR FIRST
! BUT IN QC MOD FOR THE SECOND???  

!--- end copied from obs_state_diagnostics ---

end do QC_LOOP

call free_state_window

deallocate(expected_obs)
deallocate(istatus)

end subroutine get_obs_ens_distrib_state

!-------------------------------------------------------------------------

subroutine obs_space_diagnostics(obs_ens_handle, forward_op_ens_handle, ens_size, &
   seq, keys, prior_post, num_output_members, members_index, &
   obs_val_index, OBS_KEY_COPY, &
   ens_mean_index, ens_spread_index, num_obs_in_set, &
   OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, OBS_VAL_COPY, &
   OBS_ERR_VAR_COPY, DART_qc_index)

! Do prior observation space diagnostics on the set of obs corresponding to keys

type(ensemble_type),     intent(inout) :: obs_ens_handle, forward_op_ens_handle
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

integer               :: j, k, ens_offset, forward_min, forward_max
integer               :: forward_unit, ivalue
real(r8)              :: error, diff_sd, ratio
real(r8), allocatable :: obs_temp(:)
real(r8)              :: obs_prior_mean, obs_prior_var, obs_val, obs_err_var
real(r8)              :: rvalue(1)

! Do verbose forward operator output if requested
if(output_forward_op_errors) call verbose_forward_op_output(forward_op_ens_handle, prior_post, ens_size, keys)

! Make var complete for get_copy() calls below.
call all_copies_to_all_vars(obs_ens_handle)

! allocate temp space for sending data - surely only task 0 needs to allocate this?
allocate(obs_temp(num_obs_in_set))

! Update the ensemble mean
! Get this copy to process 0
call get_copy(map_task_to_pe(obs_ens_handle, 0), obs_ens_handle, OBS_MEAN_START, obs_temp) 
! Only pe 0 gets to write the sequence
if(my_task_id() == 0) then
     ! Loop through the observations for this time
     do j = 1, obs_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_obs_values(seq, keys(j), rvalue, ens_mean_index)
     end do
  endif

! Update the ensemble spread
! Get this copy to process 0
call get_copy(map_task_to_pe(obs_ens_handle, 0), obs_ens_handle, OBS_VAR_START, obs_temp)
! Only pe 0 gets to write the sequence
if(my_task_id() == 0) then
   ! Loop through the observations for this time
   do j = 1, obs_ens_handle%num_vars
      ! update the spread in each obs
      if (obs_temp(j) /= missing_r8) then
         rvalue(1) = sqrt(obs_temp(j))
      else
         rvalue(1) = obs_temp(j)
      endif
      call replace_obs_values(seq, keys(j), rvalue, ens_spread_index)
   end do
endif

! May be possible to only do this after the posterior call...
! Update any requested ensemble members
ens_offset = members_index + 4
! Update all of these ensembles that are required to sequence file
do k = 1, num_output_members
   ! Get this copy on pe 0
   call get_copy(map_task_to_pe(obs_ens_handle, 0), obs_ens_handle, k, obs_temp)
   ! Only task 0 gets to write the sequence
   if(my_task_id() == 0) then
      ! Loop through the observations for this time
      do j = 1, obs_ens_handle%num_vars
         ! update the obs values 
         rvalue(1) = obs_temp(j)
         ivalue = ens_offset + 2 * (k - 1)
         call replace_obs_values(seq, keys(j), rvalue, ivalue)
      end do
   endif
end do

! Update the qc global value
call get_copy(map_task_to_pe(obs_ens_handle, 0), obs_ens_handle, OBS_GLOBAL_QC_COPY, obs_temp)
! Only task 0 gets to write the observations for this time
if(my_task_id() == 0) then
   ! Loop through the observations for this time
   do j = 1, obs_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_qc(seq, keys(j), rvalue, DART_qc_index)
   end do
endif

! clean up.
deallocate(obs_temp)

end subroutine obs_space_diagnostics

!-------------------------------------------------------------------------

function input_qc_ok(input_qc, threshold)

logical              :: input_qc_ok
real(r8), intent(in) :: input_qc, threshold

! Do checks on input_qc value with namelist control
! Should eventually go in qc module

! NOTE: this code has changed since the original version.
! qc values equal to the threshold are kept now; only qc
! values LARGER THAN the threshold are rejected.  
! e.g. to keep only obs with a data qc of 0, set the 
! threshold to 0.

if(nint(input_qc) <= nint(threshold)) then
   input_qc_ok = .true.
else
   input_qc_ok = .false.
endif

end function input_qc_ok

!-------------------------------------------------------------------------

subroutine filter_sync_keys_time(ens_handle, key_bounds, num_obs_in_set, time1, time2)

integer,             intent(inout)  :: key_bounds(2), num_obs_in_set
type(time_type),     intent(inout)  :: time1, time2
type(ensemble_type), intent(in)     :: ens_handle

! Have owner of copy 1 broadcast these values to all other tasks.
! Only tasks which contain copies have this info; doing it this way
! allows ntasks > nens to work.

real(r8) :: rkey_bounds(2), rnum_obs_in_set(1)
real(r8) :: rtime(4)
integer  :: days, secs
integer  :: copy1_owner, owner_index

call get_copy_owner_index(1, copy1_owner, owner_index)

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

end subroutine filter_sync_keys_time

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

call set_smoother_trace(trace_level, timestamp_level)
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
   call error_handler(E_MSG,'filter trace:',trim(msg))
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

if (do_output()) call timestamp(' '//trim(msg), pos='brief')  ! was debug

end subroutine timestamp_message

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

function failed_outlier(ratio, outlier_threshold, obs_ens_handle, OBS_KEY_COPY, j, seq)

! return true if the observation value is too far away from the ensemble mean
! and should be rejected and not assimilated.

use obs_def_mod, only : get_obs_kind
use obs_kind_mod         ! this allows you to use all the types available


real(r8),                intent(in) :: ratio
real(r8),                intent(in) :: outlier_threshold
type(ensemble_type),     intent(in) :: obs_ens_handle
integer,                 intent(in) :: OBS_KEY_COPY
integer,                 intent(in) :: j
type(obs_sequence_type), intent(in) :: seq
logical                             :: failed_outlier

! the default test is:  if (ratio > outlier_threshold) failed_outlier = .true.
! but you can add code here to do different tests for different observation 
! types.  this function is only called if this namelist item is set to true:
!  enable_special_outlier_code = .true.
! in the &filter_nml namelist.  it is intended to be customized by the user.

integer :: this_obs_key, this_obs_type
type(obs_def_type) :: obs_def
type(obs_type) :: observation
logical :: first_time = .true.

! make sure there's space to hold the observation.   this is a memory
! leak in that we never release this space, but we only allocate it
! one time so it doesn't grow.  if you are going to access the values
! of the observation or the qc values, you must change the 0, 0 below
! to match what's in your obs_seq file.  the example code below uses
! only things in the obs_def derived type and so doesn't need space
! allocated for either copies or qcs.
if (first_time) then
   call init_obs(observation, 0, 0)
   first_time = .false.
endif

call prepare_to_read_from_copies(obs_ens_handle)

! if you want to do something different based on the observation specific type:

this_obs_key = obs_ens_handle%copies(OBS_KEY_COPY, j)
call get_obs_from_key(seq, this_obs_key, observation)
call get_obs_def(observation, obs_def)
this_obs_type = get_obs_kind(obs_def)

! note that this example uses the specific type (e.g. RADIOSONDE_TEMPERATURE)
! to make decisions.  you have the observation so any other part (e.g. the
! time, the value, the error) is available to you as well.

select case (this_obs_type)

! example of specifying a different threshold value for one obs type:
!   case (RADIOSONDE_TEMPERATURE)
!      if (ratio > some_other_value) then
!         failed_outlier = .true.
!      else
!         failed_outlier = .false.
!      endif

! accept all values of this observation type no matter how far
! from the ensemble mean:
!   case (AIRCRAFT_U_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT)
!      failed_outlier = .false.

   case default
      if (ratio > outlier_threshold) then
         failed_outlier = .true.
      else
         failed_outlier = .false.
      endif

end select

end function failed_outlier

!-------------------------------------------------------------------------
!> write out failed forward operators
!> This was part of obs_space_diagnostics
subroutine verbose_forward_op_output(forward_op_ens_handle, prior_post, ens_size, keys)

type(ensemble_type), intent(inout) :: forward_op_ens_handle
integer,             intent(in)    :: prior_post
integer,             intent(in)    :: ens_size
integer,             intent(in)    :: keys(:) ! I think this is still var size

character*12 :: task
integer :: j, i
integer :: forward_unit

write(task, '(i6.6)') my_task_id()

! all tasks open file?
if(prior_post == PRIOR_DIAG) then
   forward_unit = open_file('prior_forward_ope_errors' // task, 'formatted', 'append')
else
   forward_unit = open_file('post_forward_ope_errors' // task, 'formatted', 'append')
endif

! Forward_op_ens_handle is a real representing an integer; values /= 0 get written out
do i = 1, ens_size
   do j = 1, forward_op_ens_handle%my_num_vars
      if(nint(forward_op_ens_handle%copies(i, j)) /= 0) write(forward_unit, *) keys(j), nint(forward_op_ens_handle%copies(i, j))
   end do
end do

call close_file(forward_unit)

end subroutine verbose_forward_op_output

!==================================================================
! TEST FUNCTIONS BELOW THIS POINT
!------------------------------------------------------------------
!> dump out obs_copies to file
subroutine test_obs_copies(obs_ens_handle, information)

type(ensemble_type), intent(in) :: obs_ens_handle
character(len=*),    intent(in) :: information

character*20  :: task_str !< string to hold the task number
character*129 :: file_obscopies !< output file name
integer :: i

write(task_str, '(i10)') obs_ens_handle%my_pe
file_obscopies = TRIM('obscopies_' // TRIM(ADJUSTL(information)) // TRIM(ADJUSTL(task_str)))
open(15, file=file_obscopies, status ='unknown')

do i = 1, obs_ens_handle%num_copies - 4
   write(15, *) obs_ens_handle%copies(i,:)
enddo

close(15)

end subroutine test_obs_copies

!-------------------------------------------------------------------
end program filter

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
