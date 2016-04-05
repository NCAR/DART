! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program inflate_prior_ensemble

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!------------------------------------------------------------------------------
use types_mod,            only : r8, missing_r8, metadatalength
use obs_sequence_mod,     only : read_obs_seq, obs_type, obs_sequence_type,                  &
                                 get_obs_from_key, set_copy_meta_data, get_copy_meta_data,   &
                                 get_obs_def, get_time_range_keys, set_obs_values, set_obs,  &
                                 write_obs_seq, get_num_obs, get_obs_values, init_obs,       &
                                 assignment(=), get_num_copies, get_qc, get_num_qc, set_qc,  &
                                 static_init_obs_sequence, destroy_obs, read_obs_seq_header, &
                                 set_qc_meta_data, get_expected_obs, get_first_obs,          &
                                 get_obs_time_range, delete_obs_from_seq, delete_seq_head,   &
                                 delete_seq_tail, replace_obs_values, replace_qc,            &
                                 destroy_obs_sequence, get_qc_meta_data, add_qc
use obs_def_mod,          only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                                 get_obs_kind
use obs_kind_mod,         only : assimilate_this_obs_kind, evaluate_this_obs_kind, &
                                 max_obs_kinds, get_obs_kind_index
use time_manager_mod,     only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                 operator(-), print_time
use utilities_mod,        only : register_module,  error_handler, E_ERR, E_MSG, E_DBG,       &
                                 initialize_utilities, logfileunit, nmlfileunit, timestamp,  &
                                 do_output, find_namelist_in_file, check_namelist_read,      &
                                 open_file, close_file, do_nml_file, do_nml_term
use assim_model_mod,      only : static_init_assim_model, get_model_size,                    &
                                 netcdf_file_type, init_diag_output, finalize_diag_output,   &
                                 ens_mean_for_model, end_assim_model
use assim_tools_mod,      only : filter_assim, set_assim_tools_trace
use obs_model_mod,        only : move_ahead, advance_state, set_obs_model_trace
use ensemble_manager_mod, only : init_ensemble_manager, end_ensemble_manager,                &
                                 ensemble_type, get_copy, get_my_num_copies, put_copy,       &
                                 all_vars_to_all_copies, all_copies_to_all_vars,             &
                                 read_ensemble_restart, write_ensemble_restart,              &
                                 compute_copy_mean, compute_copy_mean_sd,                    &
                                 compute_copy_mean_var, duplicate_ens, get_copy_owner_index, &
                                 get_ensemble_time, set_ensemble_time, broadcast_copy,       &
                                 prepare_to_read_from_vars, prepare_to_write_to_vars,        &
                                 prepare_to_read_from_copies, prepare_to_write_to_copies,    &
                                 prepare_to_update_copies, map_task_to_pe, map_pe_to_task,   &
                                 get_ensemble_time, set_ensemble_time
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


!------------------------------------------------------------------------------

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Some convenient global storage items
character(len=129)      :: msgstring
type(obs_type)          :: observation

integer                 :: trace_level, timestamp_level

! the public in obs_kind_mod is misnamed; it really defines the number
! of specific types (not generic kinds).  rename it for code clarity.
integer, parameter      :: max_obs_types = max_obs_kinds

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

character(len = 32 ) :: types_use_precomputed_priors(max_obs_types) = 'null'

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

namelist /inflate_prior_ensemble_nml/ async, adv_ens_command, ens_size, tasks_per_model_advance,    &
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
   silence, types_use_precomputed_priors



!----------------------------------------------------------------

! Doing this allows independent scoping for subroutines in main program file
call filter_main()

!----------------------------------------------------------------

contains

subroutine filter_main()

type(ensemble_type)         :: ens_handle, obs_ens_handle, forward_op_ens_handle
type(obs_sequence_type)     :: seq
type(netcdf_file_type)      :: PriorStateUnit, PosteriorStateUnit
type(time_type)             :: time1, first_obs_time, last_obs_time
type(time_type)             :: curr_ens_time, next_ens_time, window_time
type(adaptive_inflate_type) :: prior_inflate, post_inflate

integer, allocatable :: keys(:)

integer :: i, iunit, io, time_step_number, num_obs_in_set
integer :: ierr, last_key_used, model_size, key_bounds(2)
integer :: in_obs_copy, obs_val_index
integer :: output_state_mean_index, output_state_spread_index
integer :: prior_obs_mean_index, posterior_obs_mean_index
integer :: prior_obs_spread_index, posterior_obs_spread_index

! Global indices into ensemble storage
integer :: ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY, PRIOR_INF_SD_COPY
integer :: POST_INF_COPY, POST_INF_SD_COPY
integer :: OBS_VAL_COPY, OBS_ERR_VAR_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer :: OBS_MEAN_START, OBS_MEAN_END
integer :: OBS_VAR_START, OBS_VAR_END, TOTAL_OBS_COPIES
integer :: input_qc_index, DART_qc_index
integer :: mean_owner, mean_owners_index

integer :: num_types_precomputed
logical :: precomputed_priors(max_obs_types)
integer :: precomputed_copy_index

! For now, have model_size real storage for the ensemble mean, don't really want this
! in the long run
real(r8), allocatable :: ens_mean(:)

logical :: ds, all_gone



call filter_initialize_modules_used()

! Read the namelist entry
call find_namelist_in_file("input.nml", "inflate_prior_ensemble_nml", iunit)
read(iunit, nml = inflate_prior_ensemble_nml, iostat = io)
call check_namelist_read(iunit, io, "inflate_prior_ensemble_nml")


! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=inflate_prior_ensemble_nml)
if (do_nml_term()) write(     *     , nml=inflate_prior_ensemble_nml)

call set_trace(trace_execution, output_timestamps, silence)

call     trace_message('Inflate Prior Ensemble start')
call timestamp_message('Inflate Prior Ensemble start')

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

! Observation space inflation for posterior not currently supported
if(inf_flavor(2) == 1) call error_handler(E_ERR, 'filter_main', &
   'Posterior observation space inflation (type 1) not supported', source, revision, revdate)

! Setup the indices into the ensemble storage
ENS_MEAN_COPY        = ens_size + 1
ENS_SD_COPY          = ens_size + 2
PRIOR_INF_COPY       = ens_size + 3
PRIOR_INF_SD_COPY    = ens_size + 4
POST_INF_COPY        = ens_size + 5
POST_INF_SD_COPY     = ens_size + 6

! Can't output more ensemble members than exist
if(num_output_state_members > ens_size) num_output_state_members = ens_size
if(num_output_obs_members   > ens_size) num_output_obs_members   = ens_size

call trace_message('Before setting up space for ensembles')
! Allocate model size storage and ens_size storage for metadata for outputting ensembles
model_size = get_model_size()

! Have ens_mean on all processors for distance computations, really don't want to do this
allocate(ens_mean(model_size))

call trace_message('After  setting up space for ensembles')

! Don't currently support number of processes > model_size
if(task_count() > model_size) call error_handler(E_ERR,'filter_main', &
   'Number of processes > model size' ,source,revision,revdate)

call     trace_message('Before reading in ensemble restart files')
call timestamp_message('Before reading in ensemble restart files')

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time(time1)

! Read in restart files and initialize the ensemble storage
call filter_read_restart(ens_handle, time1, model_size)

! Read in or initialize smoother restarts as needed
if(ds) then
   call init_smoother(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)
   call smoother_read_restart(ens_handle, ens_size, model_size, time1, init_time_days)
endif

call timestamp_message('After  reading in ensemble restart files')
call     trace_message('After  reading in ensemble restart files')

call trace_message('Before initializing inflation')

! Initialize the adaptive inflation module
call adaptive_inflate_init(prior_inflate, inf_flavor(1), inf_initial_from_restart(1), &
   inf_sd_initial_from_restart(1), inf_output_restart(1), inf_deterministic(1),       &
   inf_in_file_name(1), inf_out_file_name(1), inf_diag_file_name(1), inf_initial(1),  &
   inf_sd_initial(1), inf_lower_bound(1), inf_upper_bound(1), inf_sd_lower_bound(1),  &
   ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY, 'Prior')
!call adaptive_inflate_init(post_inflate, inf_flavor(2), inf_initial_from_restart(2),  &
!   inf_sd_initial_from_restart(2), inf_output_restart(2), inf_deterministic(2),       &
!   inf_in_file_name(2), inf_out_file_name(2), inf_diag_file_name(2), inf_initial(2),  &
!   inf_sd_initial(2), inf_lower_bound(2), inf_upper_bound(2), inf_sd_lower_bound(2),  &
!   ens_handle, POST_INF_COPY, POST_INF_SD_COPY, 'Posterior')

if (do_output()) then
   if (inf_flavor(1) > 0 .and. inf_damping(1) < 1.0_r8) then
      write(msgstring, '(A,F12.6,A)') 'Prior inflation damping of ', inf_damping(1), ' will be used'
      call error_handler(E_MSG,'filter:', msgstring)
   endif
endif

call trace_message('After  initializing inflation')

call     trace_message('Before initializing output files')
call timestamp_message('Before initializing output files')

! Initialize the output sequences and state files and set their meta data
if(my_task_id() == 0) then
   call filter_generate_copy_meta_data(seq, prior_inflate, &
      PriorStateUnit, PosteriorStateUnit, in_obs_copy, output_state_mean_index, &
      output_state_spread_index, prior_obs_mean_index, posterior_obs_mean_index, &
      prior_obs_spread_index, posterior_obs_spread_index)
   if(ds) call smoother_gen_copy_meta_data(num_output_state_members, output_inflation)
else
   output_state_mean_index = 0
   output_state_spread_index = 0
   prior_obs_mean_index = 0
   posterior_obs_mean_index = 0
   prior_obs_spread_index = 0
   posterior_obs_spread_index = 0
endif

call timestamp_message('After  initializing output files')
call     trace_message('After  initializing output files')

! Compute mean and spread for inflation and state diagnostics
call all_vars_to_all_copies(ens_handle)

call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

! Output the state before inflation CSS added this here
call all_copies_to_all_vars(ens_handle) 

call trace_message('Before outputting state before inflation')
call filter_state_space_diagnostics(time1, PriorStateUnit, ens_handle, &
   model_size, num_output_state_members, &
   output_state_mean_index, output_state_spread_index, &
   output_inflation, ens_mean, ENS_MEAN_COPY, ENS_SD_COPY, &
   prior_inflate, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
call trace_message('After outputting state before inflation')
! end CSS

if(do_single_ss_inflate(prior_inflate) .or. do_varying_ss_inflate(prior_inflate)) then
      call trace_message('Before prior inflation damping and prep')
      if (inf_damping(1) /= 1.0_r8) then
         ens_handle%copies(PRIOR_INF_COPY, :) = 1.0_r8 + &
            inf_damping(1) * (ens_handle%copies(PRIOR_INF_COPY, :) - 1.0_r8)
      endif

      call filter_ensemble_inflate(ens_handle, PRIOR_INF_COPY, prior_inflate, ENS_MEAN_COPY)

      ! Recompute the the mean and spread as required for diagnostics
      call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      call trace_message('After  prior inflation damping and prep')
endif

! Back to state space for diagnostics if required
call all_copies_to_all_vars(ens_handle) 


! Ship the ensemble mean to the model; some models need this for 
! computing distances
! Who stores the ensemble mean copy
call get_copy_owner_index(ENS_MEAN_COPY, mean_owner, mean_owners_index)
! Broadcast it to everybody else
if(my_task_id() == mean_owner) then
   ens_mean = ens_handle%vars(:, mean_owners_index)
   call broadcast_send(mean_owner, ens_mean)
else
   call broadcast_recv(mean_owner, ens_mean)
endif

! Now send the mean to the model in case it's needed
call ens_mean_for_model(ens_mean)

if ((output_interval > 0) .and. &
       (time_step_number / output_interval * output_interval == time_step_number)) then
      call trace_message('Before prior state space diagnostics')
      call filter_state_space_diagnostics(time1, PriorStateUnit, ens_handle, &
         model_size, num_output_state_members, &
         output_state_mean_index, output_state_spread_index, &
         output_inflation, ens_mean, ENS_MEAN_COPY, ENS_SD_COPY, &
         prior_inflate, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
      call trace_message('After  prior state space diagnostics')
endif


! If observation space inflation, output the diagnostics
if(do_obs_inflate(prior_inflate) .and. my_task_id() == 0) & 
   call output_inflate_diagnostics(prior_inflate, ens_handle%time(1))

call trace_message('End of main filter assimilation loop, starting cleanup', 'filter:', -1)

call trace_message('Before finalizing diagnostics files')
! properly dispose of the diagnostics files
if(my_task_id() == 0) then
   ierr = finalize_diag_output(PriorStateUnit)
   ierr = finalize_diag_output(PosteriorStateUnit)
endif
call trace_message('After  finalizing diagnostics files')

call trace_message('Before writing inflation restart files if required')
! Output the restart for the adaptive inflation parameters
call adaptive_inflate_end(prior_inflate, ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
call adaptive_inflate_end(post_inflate, ens_handle, POST_INF_COPY, POST_INF_SD_COPY)
call trace_message('After  writing inflation restart files if required')

! Output a restart file if requested
call trace_message('Before writing state restart files if requested')
if(output_restart) &
   call write_ensemble_restart(ens_handle, restart_out_file_name, 1, ens_size)
!if(ds) call smoother_write_restart(1, ens_size)
call trace_message('After  writing state restart files if requested')

call trace_message('Before ensemble and obs memory cleanup')
call end_ensemble_manager(ens_handle)

call trace_message('After  ensemble and obs memory cleanup')

!if(ds) then 
!   call trace_message('Before smoother memory cleanup')
!   call smoother_end()
!   call trace_message('After  smoother memory cleanup')
!endif

call     trace_message('Inflate Prior Ensemble done')
call timestamp_message('Inflate Prior Ensemble done')
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


! Set up diagnostic output for model state, if output is desired
!PriorStateUnit     = init_diag_output('Prior_Diag', &
PriorStateUnit     = init_diag_output('before_inflation', &
                        'prior ensemble state', num_state_copies, state_meta)
!PosteriorStateUnit = init_diag_output('Posterior_Diag', &
PosteriorStateUnit = init_diag_output('after_inflation', &
                        'posterior ensemble state', num_state_copies, state_meta)


! Set the metadata for the observations.

! Set up obs ensemble mean
!num_obs_copies = in_obs_copy
!num_obs_copies = num_obs_copies + 1
!prior_meta_data = 'prior ensemble mean'
!call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
!prior_obs_mean_index = num_obs_copies
!num_obs_copies = num_obs_copies + 1
!posterior_meta_data = 'posterior ensemble mean'
!call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
!posterior_obs_mean_index = num_obs_copies 

! Set up obs ensemble spread
!num_obs_copies = num_obs_copies + 1
!prior_meta_data = 'prior ensemble spread'
!call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
!prior_obs_spread_index = num_obs_copies
!num_obs_copies = num_obs_copies + 1
!posterior_meta_data = 'posterior ensemble spread'
!call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
!posterior_obs_spread_index = num_obs_copies

! Make sure there are not too many copies requested
!if(num_output_obs_members > 10000) then
!   write(msgstring, *)'output metadata in filter needs obs ensemble size < 10000, not ',&
!                      num_output_obs_members
!   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
!endif

! Set up obs ensemble members as requested
!do i = 1, num_output_obs_members
!   write(prior_meta_data, '(a21, 1x, i6)') 'prior ensemble member', i
!   write(posterior_meta_data, '(a25, 1x, i6)') 'posterior ensemble member', i
!   num_obs_copies = num_obs_copies + 1
!   call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
!   num_obs_copies = num_obs_copies + 1
!   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
!end do


end subroutine filter_generate_copy_meta_data

!-------------------------------------------------------------------------

subroutine filter_initialize_modules_used()

! Initialize modules used that require it
call initialize_mpi_utilities('Inflate Prior Ensemble')

call register_module(source,revision,revdate)

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()

end subroutine filter_initialize_modules_used

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

subroutine filter_read_restart(ens_handle, time, model_size)

type(ensemble_type), intent(inout) :: ens_handle
type(time_type),     intent(inout) :: time
integer,             intent(in)    :: model_size

integer :: days, secs

! First initialize the ensemble manager storage
! Need enough copies for ensemble plus mean, spread, and any inflates
! Copies are ensemble, ensemble mean, variance inflation and inflation s.d.
! AVOID COPIES FOR INFLATION IF STATE SPACE IS NOT IN USE; NEEDS WORK
!!!if(prior_inflate%flavor >= 2) then
   call init_ensemble_manager(ens_handle, ens_size + 6, model_size, 1)
!!!else
!!!   call init_ensemble_manager(ens_handle, ens_size + 2, model_size, 1)
!!!endif

if (do_output()) then
   if (start_from_restart) then
      call error_handler(E_MSG,'filter_read_restart:', &
         'Reading in initial condition/restart data for all ensemble members from file(s)')
   else
      call error_handler(E_MSG,'filter_read_restart:', &
         'Reading in a single ensemble and perturbing data for the other ensemble members')
   endif
endif

! Only read in initial conditions for actual ensemble members
if(init_time_days >= 0) then
   call read_ensemble_restart(ens_handle, 1, ens_size, &
      start_from_restart, restart_in_file_name, time)
   if (do_output()) then
      call get_time(time, secs, days)
      write(msgstring, '(A)') 'By namelist control, ignoring time found in restart file.'
      call error_handler(E_MSG,'filter_read_restart:',msgstring,source,revision,revdate)
      write(msgstring, '(A,I6,1X,I5)') 'Setting initial days, seconds to ',days,secs
      call error_handler(E_MSG,'filter_read_restart:',msgstring,source,revision,revdate)
   endif
else
   call read_ensemble_restart(ens_handle, 1, ens_size, &
      start_from_restart, restart_in_file_name)
   if (ens_handle%my_num_copies > 0) time = ens_handle%time(1)
endif

! Temporary print of initial model time
if(ens_handle%my_pe == 0) then
   ! FIXME for the future: if pe 0 is not task 0, pe 0 can not print debug messages
   call get_time(time, secs, days)
   write(msgstring, *) 'initial model time of 1st ensemble member (days,seconds) ',days,secs
   call error_handler(E_DBG,'filter_read_restart',msgstring,source,revision,revdate)
endif

end subroutine filter_read_restart

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

end program inflate_prior_ensemble
