! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program filter

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

!-----------------------------------------------------------------------------------------
use types_mod,            only : r8, missing_r8
use obs_sequence_mod,     only : read_obs_seq, obs_type, obs_sequence_type,                  &
                                 get_obs_from_key, set_copy_meta_data, get_copy_meta_data,   &
                                 get_obs_def, get_time_range_keys, set_obs_values, set_obs,  &
                                 write_obs_seq, get_num_obs, get_obs_values, init_obs,       &
                                 assignment(=), get_num_copies, get_qc, get_num_qc, set_qc,  &
                                 static_init_obs_sequence, destroy_obs, read_obs_seq_header, &
                                 set_qc_meta_data, get_expected_obs, get_first_obs,          &
                                 get_obs_time_range, delete_obs_from_seq, delete_seq_head,   &
                                 delete_seq_tail
use obs_def_mod,          only : obs_def_type, get_obs_def_error_variance, get_obs_def_time
use time_manager_mod,     only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                 operator(-)
use utilities_mod,        only : register_module,  error_handler, E_ERR, E_MSG, E_DBG,       &
                                 initialize_utilities, logfileunit, timestamp, do_output,    &
                                 find_namelist_in_file, check_namelist_read, open_file,      &
                                 close_file
use assim_model_mod,      only : static_init_assim_model, get_model_size,                    &
                                 netcdf_file_type, init_diag_output, finalize_diag_output,   & 
                                 aoutput_diagnostics, ens_mean_for_model
use assim_tools_mod,      only : filter_assim
use obs_model_mod,        only : move_ahead
use ensemble_manager_mod, only : init_ensemble_manager, end_ensemble_manager,                &
                                 ensemble_type, get_copy, get_my_num_copies, put_copy,       &
                                 all_vars_to_all_copies, all_copies_to_all_vars,             &
                                 read_ensemble_restart, write_ensemble_restart,              &
                                 compute_copy_mean, compute_copy_mean_sd,                    &
                                 compute_copy_mean_var, duplicate_ens, get_copy_owner_index
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
                                 smoother_ss_diagnostics, smoother_inc_lags,                 &
                                 smoother_end


!-----------------------------------------------------------------------------------------

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! Some convenient global storage items
character(len=129)      :: msgstring
type(obs_type)          :: observation

! Defining whether diagnostics are for prior or posterior
integer, parameter :: PRIOR_DIAG = 0, POSTERIOR_DIAG = 2

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async = 0, ens_size = 20
logical  :: start_from_restart = .false., output_restart = .false.
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days    = 0
integer  :: init_time_seconds = 0
! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days    = -1
integer  :: first_obs_seconds = -1
integer  :: last_obs_days     = -1
integer  :: last_obs_seconds  = -1
! Control diagnostic output for state variables
integer  :: num_output_state_members = 0
integer  :: num_output_obs_members   = 0
integer  :: output_interval = 1
integer  :: num_groups = 1
real(r8) :: outlier_threshold = -1.0_r8
real(r8) :: ncep_qc_threshold = 4.0_r8
logical  :: output_forward_op_errors = .false.

character(len = 129) :: obs_sequence_in_name  = "obs_seq.out",    &
                        obs_sequence_out_name = "obs_seq.final",  &
                        restart_in_file_name  = 'filter_ics',     &
                        restart_out_file_name = 'filter_restart', &
                        adv_ens_command       = './advance_model.csh'

!                  == './advance_model.csh'    -> advance ensemble using a script

! Inflation namelist entries follow, first entry for prior, second for posterior
! inf_flavor is 0:none, 1:obs space, 2: varying state space, 3: fixed state_space
integer              :: inf_flavor(2)             = 0
logical              :: inf_start_from_restart(2) = .false.
logical              :: inf_output_restart(2)     = .false.
logical              :: inf_deterministic(2)      = .true.
character(len = 129) :: inf_in_file_name(2)       = 'not_initialized',    &
                        inf_out_file_name(2)      = 'not_initialized',    &
                        inf_diag_file_name(2)     = 'not_initialized'
real(r8)             :: inf_initial(2)            = 1.0_r8
real(r8)             :: inf_sd_initial(2)         = 0.0_r8
real(r8)             :: inf_lower_bound(2)        = 1.0_r8
real(r8)             :: inf_upper_bound(2)        = 1000000.0_r8
real(r8)             :: inf_sd_lower_bound(2)     = 0.0_r8

namelist /filter_nml/ async, adv_ens_command, ens_size, start_from_restart, &
                      output_restart, obs_sequence_in_name, obs_sequence_out_name, &
                      restart_in_file_name, restart_out_file_name, init_time_days, &
                      init_time_seconds, first_obs_days, first_obs_seconds, last_obs_days, &
                      last_obs_seconds, num_output_state_members, &
                      num_output_obs_members, output_interval, num_groups, outlier_threshold, &
                      ncep_qc_threshold, output_forward_op_errors, inf_flavor, &
                      inf_start_from_restart, &
                      inf_output_restart, inf_deterministic, inf_in_file_name, &
                      inf_out_file_name, inf_diag_file_name, inf_initial, inf_sd_initial, &
                      inf_lower_bound, inf_upper_bound, inf_sd_lower_bound

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
type(adaptive_inflate_type) :: prior_inflate, post_inflate

integer,    allocatable :: keys(:)
integer                 :: i, iunit, io, time_step_number, num_obs_in_set
integer                 :: ierr, last_key_used, model_size, key_bounds(2)
integer                 :: in_obs_copy, obs_val_index
integer                 :: output_state_mean_index, output_state_spread_index
integer                 :: prior_obs_mean_index, posterior_obs_mean_index
integer                 :: prior_obs_spread_index, posterior_obs_spread_index
! Global indices into ensemble storage
integer                 :: ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY, PRIOR_INF_SD_COPY
integer                 :: POST_INF_COPY, POST_INF_SD_COPY
integer                 :: OBS_VAL_COPY, OBS_ERR_VAR_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer                 :: OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END
integer                 :: OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END, TOTAL_OBS_COPIES
integer                 :: input_qc_index, DART_qc_index
integer                 :: mean_owner, mean_owners_index

real(r8)                :: rkey_bounds(2), rnum_obs_in_set(1), small_temp(1)

! For now, have model_size real storage for the ensemble mean, don't really want this
! in the long run
real(r8), allocatable   :: ens_mean(:)

logical                 :: ds, all_gone

call filter_initialize_modules_used(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)

! Read the namelist entry
call find_namelist_in_file("input.nml", "filter_nml", iunit)
read(iunit, nml = filter_nml, iostat = io)
call check_namelist_read(iunit, io, "filter_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'filter_main','filter_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=filter_nml)
if (do_output()) write(     *     , nml=filter_nml)

! Make sure ensemble size is at least 2 (NEED MANY OTHER CHECKS)
if(ens_size < 2) then
   write(msgstring, *) 'ens_size in namelist is ', ens_size, ': Must be > 1'
   call error_handler(E_ERR,'filter_main', msgstring, source, revision, revdate)
endif

! See if smoothing is turned on
ds = do_smoothing()

! Make sure inflation options are legal
do i = 1, 2
   if(inf_flavor(i) < 0 .or. inf_flavor(i) > 3) then
      write(msgstring, *) 'inf_flavor=', inf_flavor(i), ' Must be 0, 1, 2, 3 '
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
OBS_ERR_VAR_COPY     = ens_size + 1
OBS_VAL_COPY         = ens_size + 2
OBS_KEY_COPY         = ens_size + 3
OBS_GLOBAL_QC_COPY   = ens_size + 4
OBS_PRIOR_MEAN_START = ens_size + 5
OBS_PRIOR_MEAN_END   = OBS_PRIOR_MEAN_START + num_groups - 1
OBS_PRIOR_VAR_START  = OBS_PRIOR_MEAN_START + num_groups
OBS_PRIOR_VAR_END    = OBS_PRIOR_VAR_START + num_groups - 1

! Can't output more ensemble members than exist
if(num_output_state_members > ens_size) num_output_state_members = ens_size
if(num_output_obs_members   > ens_size) num_output_obs_members   = ens_size

! Initialize the obs_sequence; every pe gets a copy for now
call filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, input_qc_index, DART_qc_index)

! Allocate model size storage and ens_size storage for metadata for outputting ensembles
model_size = get_model_size()

! Have ens_mean on all processors for distance computations, really don't want to do this
allocate(ens_mean(model_size))

! Don't currently support number of processes > model_size
if(task_count() > model_size) call error_handler(E_ERR,'filter_main', &
   'Number of processes > model size' ,source,revision,revdate)

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time(time1)

! Read in restart files and initialize the ensemble storage
call filter_read_restart(ens_handle, time1, model_size)

! Read in or initialize smoother restarts as needed
if(ds) call smoother_read_restart(ens_handle, ens_size, model_size, time1, init_time_days)

! Initialize the adaptive inflation module
call adaptive_inflate_init(prior_inflate, inf_flavor(1), inf_start_from_restart(1), &
   inf_output_restart(1), inf_deterministic(1), inf_in_file_name(1), inf_out_file_name(1), &
   inf_diag_file_name(1), inf_initial(1), inf_sd_initial(1), inf_lower_bound(1), &
   inf_upper_bound(1), inf_sd_lower_bound(1), ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
call adaptive_inflate_init(post_inflate, inf_flavor(2), inf_start_from_restart(2), &
   inf_output_restart(2), inf_deterministic(2), inf_in_file_name(2), inf_out_file_name(2), &
   inf_diag_file_name(2), inf_initial(2), inf_sd_initial(2), inf_lower_bound(2), &
   inf_upper_bound(2), inf_sd_lower_bound(2), ens_handle, POST_INF_COPY, POST_INF_SD_COPY)

! Initialize the output sequences and state files and set their meta data
if(my_task_id() == 0) then
   call filter_generate_copy_meta_data(seq, prior_inflate, &
   PriorStateUnit, PosteriorStateUnit, in_obs_copy, output_state_mean_index, &
   output_state_spread_index, prior_obs_mean_index, posterior_obs_mean_index, &
   prior_obs_spread_index, posterior_obs_spread_index)
   if(ds) call smoother_gen_copy_meta_data(num_output_state_members)
endif

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

! Time step number is used to do periodic diagnostic output
time_step_number = -1

AdvanceTime : do

   if(my_task_id() == 0) write(*, *) 'starting advance time loop;'
   time_step_number = time_step_number + 1

   ! Advance the lagged distribution
   if(ds) call advance_smoother(ens_handle)

   ! Get the model to a good time to use a next set of observations
   call move_ahead(ens_handle, ens_size, seq, last_key_used, &
      key_bounds, num_obs_in_set, async, adv_ens_command)


   ! Only processes with an ensemble copy know to exit; 
   ! For now, let process 0 broadcast it's value of key_bounds
   ! This will synch the loop here and allow everybody to exit
   ! Need to clean up and have a broadcast that just sends a single integer???
   ! PAR For now, can only broadcast real pair of arrays
   if(my_task_id() == 0) then
      rkey_bounds = key_bounds
      rnum_obs_in_set(1) = num_obs_in_set
      call broadcast_send(0, rkey_bounds, rnum_obs_in_set)
   else
      call broadcast_recv(0, rkey_bounds, rnum_obs_in_set)
      key_bounds = rkey_bounds
      num_obs_in_set = rnum_obs_in_set(1)
   endif
   if(key_bounds(1) < 0) exit AdvanceTime

   ! Create an ensemble for the observations from this time plus &
   ! obs_error_variance, observed value, key from sequence, then mean for each
   ! group, then variance for each group
   TOTAL_OBS_COPIES = ens_size + 4 + 2*num_groups
   call init_ensemble_manager(obs_ens_handle, TOTAL_OBS_COPIES, num_obs_in_set, 1)
   ! Also need a qc field for copy of each observation
   call init_ensemble_manager(forward_op_ens_handle, ens_size, num_obs_in_set, 1)

   ! Allocate storage for the keys for this number of observations
   allocate(keys(num_obs_in_set))

   ! Get all the keys associated with this set of observations
   ! Is there a way to distribute this?
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   ! Compute mean and spread for inflation and state diagnostics
   call all_vars_to_all_copies(ens_handle)
   call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

   if(do_single_ss_inflate(prior_inflate) .or. do_varying_ss_inflate(prior_inflate)) then
      call filter_ensemble_inflate(ens_handle, PRIOR_INF_COPY, prior_inflate, ENS_MEAN_COPY)
      ! Recompute the mean or the mean and spread as required for diagnostics
      call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)
   endif

   ! Back to state space for diagnostics if required
   call all_copies_to_all_vars(ens_handle) 

   ! Compute the ensemble of prior observations, load up the obs_err_var and obs_values
   ! ens_size is the number of regular ensemble members, not the number of copies
   call get_obs_ens(ens_handle, obs_ens_handle, forward_op_ens_handle, seq, keys, &
      obs_val_index, num_obs_in_set, OBS_ERR_VAR_COPY, OBS_VAL_COPY)

   ! Although they are integer, keys are one 'copy' of obs ensemble (the last one?)
   call put_copy(0, obs_ens_handle, OBS_KEY_COPY, keys * 1.0_r8)

   ! Ship the ensemble mean to the model; some models need this for computing distances
   ! Who stores the ensemble mean copy
   call get_copy_owner_index(ENS_MEAN_COPY, mean_owner, mean_owners_index)
   ! Broadcast it to everybody else
   if(my_task_id() == mean_owner) then
      ens_mean = ens_handle%vars(:, mean_owners_index)
      call broadcast_send(mean_owner, ens_mean, small_temp)
   else
      call broadcast_recv(mean_owner, ens_mean, small_temp)
   endif

   ! Now send the mean to the model in case it's needed
   call ens_mean_for_model(ens_mean)

   ! Do prior state space diagnostic output as required
   ! Use ens_mean which is declared model_size for temp storage in diagnostics
   if(time_step_number / output_interval * output_interval == time_step_number) then
      call filter_state_space_diagnostics(PriorStateUnit, ens_handle, model_size, &
         num_output_state_members, output_state_mean_index, output_state_spread_index, &
         ens_mean, ENS_MEAN_COPY, ENS_SD_COPY, &
         prior_inflate, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
   endif
  
   ! Do prior observation space diagnostics and associated quality control
   call obs_space_diagnostics(obs_ens_handle, forward_op_ens_handle, ens_size, seq, keys, &
      PRIOR_DIAG, num_output_obs_members, in_obs_copy + 1, &
      prior_obs_mean_index, prior_obs_spread_index, num_obs_in_set, &
      OBS_PRIOR_MEAN_START, OBS_PRIOR_VAR_START, OBS_GLOBAL_QC_COPY, OBS_VAL_COPY, &
      OBS_ERR_VAR_COPY, DART_qc_index)
  
   ! Need obs to be copy complete for assimilation
   call all_vars_to_all_copies(obs_ens_handle)
   call filter_assim(ens_handle, obs_ens_handle, seq, keys, ens_size, num_groups, &
      obs_val_index, prior_inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
      OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START, &
      OBS_PRIOR_VAR_END, inflate_only = .false.)

   ! Do the update for the smoother lagged fields, too.
   ! Would be more efficient to do these all at once inside filter_assim in the future
   call smoother_assim(obs_ens_handle, seq, keys, ens_size, num_groups, obs_val_index, &
      ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, &
      OBS_GLOBAL_QC_COPY, OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START, &
      OBS_PRIOR_VAR_END)

   ! Already transformed, so compute mean and spread for state diag as needed
   call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

!-------- Test of posterior inflate ----------------

   if(do_single_ss_inflate(post_inflate) .or. do_varying_ss_inflate(post_inflate)) then
      call filter_ensemble_inflate(ens_handle, POST_INF_COPY, post_inflate, ENS_MEAN_COPY) 
      ! Recompute the mean or the mean and spread as required for diagnostics
      call compute_copy_mean_sd(ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)
   endif

!-------- End of posterior  inflate ----------------

   ! Now back to var complete for diagnostics
   call all_copies_to_all_vars(ens_handle)
   
   ! Compute the ensemble of posterior observations, load up the obs_err_var and obs_values
   ! ens_size is the number of regular ensemble members, not the number of copies
   call get_obs_ens(ens_handle, obs_ens_handle, forward_op_ens_handle, seq, keys, &
      obs_val_index, num_obs_in_set, OBS_ERR_VAR_COPY, OBS_VAL_COPY)

   if(ds) call smoother_mean_spread(ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

   ! Do posterior state space diagnostic output as required
   if(time_step_number / output_interval * output_interval == time_step_number) then
      call filter_state_space_diagnostics(PosteriorStateUnit, ens_handle, model_size, &
         num_output_state_members, output_state_mean_index, output_state_spread_index, &
         ens_mean, ENS_MEAN_COPY, ENS_SD_COPY, &
         post_inflate, POST_INF_COPY, POST_INF_SD_COPY)
      ! Cyclic storage for lags with most recent pointed to by smoother_head
      ! ens_mean is passed to avoid extra temp storage in diagnostics
      call smoother_ss_diagnostics(model_size, num_output_state_members, ens_mean, &
         ENS_MEAN_COPY, ENS_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY)
   endif

   ! Increment the number of current lags if smoother set not yet fully spun up
   if(ds) call smoother_inc_lags()
  
   ! Do posterior observation space diagnostics
   call obs_space_diagnostics(obs_ens_handle, forward_op_ens_handle, ens_size, seq, keys, &
      POSTERIOR_DIAG, num_output_obs_members, in_obs_copy + 2, &
      posterior_obs_mean_index, posterior_obs_spread_index, &
      num_obs_in_set, OBS_PRIOR_MEAN_START, OBS_PRIOR_VAR_START, OBS_GLOBAL_QC_COPY, &
       OBS_VAL_COPY, OBS_ERR_VAR_COPY, DART_qc_index)
   

!-------- Test of posterior inflate ----------------
 
   if(do_single_ss_inflate(post_inflate) .or. do_varying_ss_inflate(post_inflate)) then

      ! Ship the ensemble mean to the model; some models need this for computing distances
      ! Who stores the ensemble mean copy
      call get_copy_owner_index(ENS_MEAN_COPY, mean_owner, mean_owners_index)
      ! Broadcast it to everybody else
      if(my_task_id() == mean_owner) then
         ens_mean = ens_handle%vars(:, mean_owners_index)
         call broadcast_send(mean_owner, ens_mean, small_temp)
      else
         call broadcast_recv(mean_owner, ens_mean, small_temp)
      endif
   
      ! Now send the mean to the model in case it's needed
      call ens_mean_for_model(ens_mean)

      ! Need obs to be copy complete for assimilation: IS NEXT LINE REQUIRED???
      call all_vars_to_all_copies(obs_ens_handle)
      call filter_assim(ens_handle, obs_ens_handle, seq, keys, ens_size, num_groups, &
         obs_val_index, post_inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
         POST_INF_COPY, POST_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
         OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START, &
         OBS_PRIOR_VAR_END, inflate_only = .true.)
      call all_copies_to_all_vars(ens_handle)
   endif

!-------- End of posterior  inflate ----------------

   ! If observation space inflation, output the diagnostics
   if(do_obs_inflate(prior_inflate) .and. my_task_id() == 0) & 
      call output_inflate_diagnostics(prior_inflate, ens_handle%time(1))

   ! Deallocate storage used for keys for each set
   deallocate(keys)

   ! The last key used is updated to move forward in the observation sequence
   last_key_used = key_bounds(2)

   ! Free up the obs ensemble space; LATER, can just keep it if obs are same size next time
   call end_ensemble_manager(obs_ens_handle)
   call end_ensemble_manager(forward_op_ens_handle)

end do AdvanceTime

! properly dispose of the diagnostics files
if(my_task_id() == 0) then
   ierr = finalize_diag_output(PriorStateUnit)
   ierr = finalize_diag_output(PosteriorStateUnit)
endif

! Only pe 0 outputs the observation space diagnostic file
if(my_task_id() == 0) call write_obs_seq(seq, obs_sequence_out_name)

! Output the restart for the adaptive inflation parameters
call adaptive_inflate_end(prior_inflate, ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
call adaptive_inflate_end(post_inflate, ens_handle, POST_INF_COPY, POST_INF_SD_COPY)

! Output a restart file if requested
if(output_restart) &
   call write_ensemble_restart(ens_handle, restart_out_file_name, 1, ens_size)
if(ds) call smoother_write_restart(1, ens_size)
call end_ensemble_manager(ens_handle)

if(my_task_id() == 0) then 
   write(logfileunit,*)'FINISHED filter.'
   write(logfileunit,*)
endif

call finalize_mpi_utilities()

! Master task must close the log file
if(my_task_id() == 0) call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

! Free up the observation kind
call destroy_obs(observation)

if(ds) call smoother_end()

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

character(len=129) :: prior_meta_data, posterior_meta_data
! The 4 is for mean, spread plus mean and spread of inflation
character(len=129) :: state_meta(num_output_state_members + 4)
integer :: i, ensemble_offset, num_state_copies, num_obs_copies

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

! Netcdf output diagnostics for inflation; inefficient for single spatial
!!! nsc - temporary fix.  for now always allow space in the netcdf file for
!!! a single copy of the inflation values.  this should be changed to handle
!!! prior, posterior, both, or neither inflation.  but it must match what is
!!! written out in filter_state_space_diagnostics() below.
!!!if(do_single_ss_inflate(prior_inflate) .or. do_varying_ss_inflate(prior_inflate)) then
   num_state_copies = num_state_copies + 2
   state_meta(num_state_copies -1) = 'inflation mean'
   state_meta(num_state_copies) = 'inflation sd'
!!!endif

! Set up diagnostic output for model state, if output is desired
PriorStateUnit     = init_diag_output('Prior_Diag', &
                        'prior ensemble state', num_state_copies, state_meta)
PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                        'posterior ensemble state', num_state_copies, state_meta)


! Set up the metadata for the output ensemble observations space file
! Set up ensemble mean if requested
num_obs_copies = in_obs_copy
num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
prior_obs_mean_index = num_obs_copies
num_obs_copies = num_obs_copies + 1
posterior_meta_data = 'posterior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
posterior_obs_mean_index = num_obs_copies 

! Set up ensemble spread if requested
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

! Set up ensemble members as requested
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

subroutine filter_initialize_modules_used(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in) :: POST_INF_COPY, POST_INF_SD_COPY

! Initialize modules used that require it
call initialize_mpi_utilities()

call initialize_utilities('Filter')
call register_module(source,revision,revdate)

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()

! Initialize the smoother module
call init_smoother(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)

end subroutine filter_initialize_modules_used

!-------------------------------------------------------------------------

subroutine filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, &
   input_qc_index, DART_qc_index)

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(out)   :: in_obs_copy, obs_val_index
integer,                 intent(out)   :: input_qc_index, DART_qc_index

!!!real(r8)             :: qc(1)
character(len = 129) :: qc_meta_data = 'DART quality control'
character(len = 129) :: obs_seq_read_format
integer              :: obs_seq_file_id, num_obs_copies
integer              :: tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, qc_num_inc, num_qc
logical              :: pre_I_format

! Determine the number of output obs space fields
! Prior and posterior values for all selected fields (so times 2)
num_obs_copies = 2 * num_output_obs_members + 4

! Input file can have one qc field or none, not prepared to have more
! The one that exists would be the NCEP and perfect_model_obs generated values in general
call read_obs_seq_header(obs_sequence_in_name, tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, &
   obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)
if(tnum_qc == 0) then
   input_qc_index = 0
   DART_qc_index = 1
else if(tnum_qc == 1) then
   input_qc_index = 1
   DART_qc_index = 2
else
   write(*, msgstring) 'input obs_seq file has ', tnum_qc, ' qc fields; must be < 2'
   call error_handler(E_ERR,'filter_setup_obs_sequence', msgstring, source, revision, revdate)
endif

! Need a new qc field for the DART quality control
qc_num_inc = 1
! Read in with enough space for diagnostic output values and add'l qc field
call read_obs_seq(obs_sequence_in_name, num_obs_copies, qc_num_inc, 0, seq)

! Get num of obs copies and num_qc
num_qc = get_num_qc(seq)
in_obs_copy = get_num_copies(seq) - num_obs_copies

! Create an observation type temporary for use in filter
call init_obs(observation, get_num_copies(seq), num_qc)

! Set initial DART quality control to 0 for all observations
! Leaving them uninitialized, obs_space_diagnostics should set them all without reading them
call set_qc_meta_data(seq, DART_qc_index, qc_meta_data)
!!!qc(1) = 0.0_r8
!!!do i = 1, get_num_obs(seq)
   !!!call get_obs_from_key(seq, i, observation)
   !!!call set_qc(observation, qc, DART_qc_index)
   !!!call set_obs(seq, observation, i)
!!!end do

! Determine which copy has actual obs
obs_val_index = get_obs_copy_index(seq)

end subroutine filter_setup_obs_sequence

!-------------------------------------------------------------------------

function get_obs_copy_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_copy_index

integer :: i

! Determine which copy in sequence has actual obs
!--------
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

! Only read in initial conditions for actual ensemble members
if(init_time_days >= 0) then
   call read_ensemble_restart(ens_handle, 1, ens_size, &
      start_from_restart, restart_in_file_name, time)
else
   call read_ensemble_restart(ens_handle, 1, ens_size, &
      start_from_restart, restart_in_file_name)
endif

! Temporary print of initial model time
if(my_task_id() == 0) then
   call get_time(time, secs, days)
   write(msgstring, *) 'initial model time of 1st ensemble member (days,seconds) ',days,secs
endif
call error_handler(E_DBG,'filter_read_restart',msgstring,source,revision,revdate)

end subroutine filter_read_restart

!-------------------------------------------------------------------------

subroutine filter_ensemble_inflate(ens_handle, inflate_copy, inflate, ENS_MEAN_COPY)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: inflate_copy, ENS_MEAN_COPY
type(adaptive_inflate_type), intent(inout) :: inflate

integer :: j, group, grp_bot, grp_top, grp_size

! Assumes that the ensemble is copy complete

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

subroutine get_obs_ens(ens_handle, obs_ens_handle, forward_op_ens_handle, seq, keys, &
   obs_val_index, num_obs_in_set, OBS_ERR_VAR_COPY, OBS_VAL_COPY)

! Computes forward observation operators and related quality control indicators.

type(ensemble_type),     intent(in)    :: ens_handle
type(ensemble_type),     intent(out)   :: obs_ens_handle, forward_op_ens_handle 
type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: keys(:)
integer,                 intent(in)    ::  obs_val_index, num_obs_in_set, OBS_ERR_VAR_COPY
integer,                 intent(in)    :: OBS_VAL_COPY

real(r8)           :: input_qc(1), obs_value(1), obs_err_var
integer            :: j, k, my_num_copies, istatus , global_ens_index
logical            :: evaluate_this_ob, assimilate_this_ob
type(obs_def_type) :: obs_def

! Assumed that both ensembles are var complete
! Each PE must loop to compute its copies of the forward operators
! IMPORTANT, IT IS ASSUMED THAT ACTUAL ENSEMBLES COME FIRST

! Loop through my copies and compute expected value
my_num_copies = get_my_num_copies(obs_ens_handle)

! Loop through all observations in the set
ALL_OBSERVATIONS: do j = 1, num_obs_in_set
   ! Get the information on this observation by placing it in temporary
   call get_obs_from_key(seq, keys(j), observation)
   call get_obs_def(observation, obs_def)
   ! Check to see if this observation fails input qc test
! PAR WARNING: WHAT IF THERE IS NO INPUT QC FIELD??? NEED TO BE PREPARED FOR THIS
   call get_qc(observation, input_qc, qc_indx = 1)
   ! If it is bad, set forward operator status value to -99 and return missing_r8 for obs_value

   ! PAR THIS SUBROUTINE SHOULD EVENTUALLY GO IN THE QUALITY CONTROL MODULE
   if(.not. input_qc_ok(input_qc(1))) then
      ! The forward operator value is set to -99 if prior qc was failed
      forward_op_ens_handle%vars(j, :) = -99
      obs_ens_handle%vars(j, 1:my_num_copies) = missing_r8 
      ! No need to do anything else for a failed observation
      cycle ALL_OBSERVATIONS
   endif

   ! Get the observation value and error variance
   call get_obs_values(observation, obs_value(1:1), obs_val_index)
   obs_err_var = get_obs_def_error_variance(obs_def)

   ! Loop through all copies stored by this process and set values as needed
   do k = 1, my_num_copies
      global_ens_index = obs_ens_handle%my_copies(k)
      ! If I have a copy that is a standard ensemble member, compute expected value
      if(global_ens_index <= ens_size) then
         ! Compute the expected observation value; direct access to storage
         call get_expected_obs(seq, keys(j:j), ens_handle%vars(:, k), &
            obs_ens_handle%vars(j:j, k), istatus, assimilate_this_ob, evaluate_this_ob)
         ! If istatus is 0 (successful) then put 0 for assimilate, -1 for evaluate only
         ! and -2 for neither evaluate or assimilate. Otherwise pass through the istatus
         ! in the forward operator evaluation field
         if(istatus == 0) then
            if(assimilate_this_ob) then
               forward_op_ens_handle%vars(j, k) = 0
            else if(evaluate_this_ob) then
               forward_op_ens_handle%vars(j, k) = -1
            else
               forward_op_ens_handle%vars(j, k) = -2
            endif
         else
            forward_op_ens_handle%vars(j, k) = istatus
         endif
      ! Otherwise, see if this is the copy for error variance or observed value         
      else if(global_ens_index == OBS_ERR_VAR_COPY) then
         ! This copy is the instrument observation error variance; read and store
         obs_ens_handle%vars(j, k) = obs_err_var
      else if(global_ens_index == OBS_VAL_COPY) then
         ! This copy is the observation from the instrument; read and store
         obs_ens_handle%vars(j, k) = obs_value(1)
      endif
   end do
end do ALL_OBSERVATIONS

end subroutine get_obs_ens

!-------------------------------------------------------------------------

subroutine obs_space_diagnostics(obs_ens_handle, forward_op_ens_handle, ens_size, &
   seq, keys, prior_post, num_output_members, members_index, &
   ens_mean_index, ens_spread_index, num_obs_in_set, &
   OBS_PRIOR_MEAN_START, OBS_PRIOR_VAR_START, OBS_GLOBAL_QC_COPY, OBS_VAL_COPY, &
   OBS_ERR_VAR_COPY, DART_qc_index)

! Do prior observation space diagnostics on the set of obs corresponding to keys

type(ensemble_type),     intent(inout) :: obs_ens_handle, forward_op_ens_handle
integer,                 intent(in)    :: ens_size
integer,                 intent(in)    :: num_obs_in_set
integer,                 intent(in)    :: keys(num_obs_in_set), prior_post
integer,                 intent(in)    :: num_output_members, members_index
integer,                 intent(in)    :: ens_mean_index, ens_spread_index
type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: OBS_PRIOR_MEAN_START, OBS_PRIOR_VAR_START
integer,                 intent(in)    :: OBS_GLOBAL_QC_COPY, OBS_VAL_COPY
integer,                 intent(in)    :: OBS_ERR_VAR_COPY, DART_qc_index

integer               :: j, k, ens_offset, forward_min, forward_max, forward_unit
real(r8)              :: error, diff_sd, ratio, obs_temp(num_obs_in_set)
real(r8)              :: obs_prior_mean, obs_prior_var, obs_val, obs_err_var
real(r8), allocatable :: forward_temp(:)
logical               :: do_outlier

! Assume that mean and spread have been computed if needed???
! Assume that things are copy complete???

! Getting the copies requires communication, so need to only do it
! once per copy. This puts the observation loop on the inside
! which may itself be expensive. Check out cost later.

!PAR: REMEMBER THAT SOME THINGS NEED NOT BE DONE IF QC IS BAD ALREADY!!!

! Compute the ensemble total mean and sd if required for output
! Also compute the mean and spread if qc and outlier_threshold check requested
do_outlier = (prior_post == PRIOR_DIAG .and. outlier_threshold > 0.0_r8)

! Do verbose forward operator output if requested
if(output_forward_op_errors) then
   ! Need to open a file for prior and posterior output
   if(my_task_id() == 0) then
      if(prior_post == PRIOR_DIAG) then
         forward_unit = open_file('prior_forward_op_errors', 'formatted', 'append')
      else
         forward_unit = open_file('post_forward_op_errors', 'formatted', 'append')
      endif
   endif

   allocate(forward_temp(num_obs_in_set))
   do k = 1, ens_size
      ! Get this copy to PE 0
      call get_copy(0, forward_op_ens_handle, k, forward_temp)

      ! Loop through each observation in set for this copy
      ! Forward temp is a real representing an integer; values of 1 or greater get out
      if(my_task_id() == 0) then
         do j = 1, num_obs_in_set
            if(nint(forward_temp(j)) > -1) write(forward_unit, *) keys(j), k, nint(forward_temp(j))
         end do
      endif
   end do
   ! PE 0 does the output for each copy in turn
   deallocate(forward_temp)
   if(my_task_id() == 0) then
      call close_file(forward_unit)
   endif
endif

!PAR DO THESE ALWAYS NEED TO BE DONE? SEE REVERSE ONES AT END, TOO
call all_vars_to_all_copies(obs_ens_handle)
call all_vars_to_all_copies(forward_op_ens_handle)

! Compute mean and spread
call compute_copy_mean_var(obs_ens_handle, &
      1, ens_size, OBS_PRIOR_MEAN_START, OBS_PRIOR_VAR_START)

! At this point can compute outlier test and consolidate forward operator qc
do j = 1, obs_ens_handle%my_num_vars
   forward_max = nint(maxval(forward_op_ens_handle%copies(1:ens_size, j)))
   forward_min = nint(minval(forward_op_ens_handle%copies(1:ens_size, j)))
   ! Now do a case statement to figure out what the qc result should be
   ! For prior, have to test for a bunch of stuff
   if(prior_post == PRIOR_DIAG) then
      if(forward_min == -99) then              ! Failed prior qc in get_obs_ens
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 6
      else if(forward_min == -2) then          ! Observation not used via namelist
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 5
      else if(forward_max > 0) then            ! At least one forward operator failed
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 4
      else if(forward_min == -1) then          ! Observation to be evaluated only
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 1
      else if(forward_min == 0) then           ! All clear, assimilate this ob
         obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 0
      endif
        
      ! PAR: THIS SHOULD BE IN QC MODULE 
      ! Check on the outlier threshold quality control: move to QC module
      if(do_outlier .and. (obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) < ncep_qc_threshold)) then
         obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START, j)
         obs_prior_var = obs_ens_handle%copies(OBS_PRIOR_VAR_START, j)
         obs_val = obs_ens_handle%copies(OBS_VAL_COPY, j)
         obs_err_var = obs_ens_handle%copies(OBS_ERR_VAR_COPY, j)
         error = obs_prior_mean - obs_val
         diff_sd = sqrt(obs_prior_var + obs_err_var)
         ratio = abs(error / diff_sd)
         if(ratio > outlier_threshold) obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 7
      endif

   else
      ! For posterior, only check for failed forward operator if prior successful
      if(forward_min == 0 .and. forward_max > 0) then
         if(obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) == 0) &
            obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 2
         if(obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) == 1) &
            obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, j) = 3
      endif
   endif
enddo

! PAR: NEED TO BE ABLE TO OUTPUT DETAILS OF FAILED FORWARD OBSEVATION OPERATORS
! OR FAILED OUTLIER, ETC. NEEDS TO BE DONE BY PE 0 ONLY. PUT IT HERE FOR FIRST
! BUT IN QC MOD FOR THE SECOND???  

! Is this next call really needed, or does it already exist on input???
call all_copies_to_all_vars(obs_ens_handle)

! Output the ensemble mean
! Get this copy to process 0
call get_copy(0, obs_ens_handle, OBS_PRIOR_MEAN_START, obs_temp)
! Only pe 0 gets to write the sequence
if(my_task_id() == 0) then
     ! Loop through the observations for this time
     do j = 1, obs_ens_handle%num_vars
        call get_obs_from_key(seq, keys(j), observation)
        call set_obs_values(observation, obs_temp(j:j), ens_mean_index)
        ! Store the observation into the sequence
        call set_obs(seq, observation, keys(j))
     end do
  endif

! If requested, output the ensemble spread
! Get this copy to process 0
call get_copy(0, obs_ens_handle, OBS_PRIOR_VAR_START, obs_temp)
! Only pe 0 gets to write the sequence
if(my_task_id() == 0) then
   ! Loop through the observations for this time
   do j = 1, obs_ens_handle%num_vars
      call get_obs_from_key(seq, keys(j), observation)
      call set_obs_values(observation, sqrt(obs_temp(j:j)), ens_spread_index)
      ! Store the observation into the sequence
      call set_obs(seq, observation, keys(j))
   end do
endif

! May be possible to only do this after the posterior call...
! Output any requested ensemble members
ens_offset = members_index + 4
! Output all of these ensembles that are required to sequence file
do k = 1, num_output_members
   ! Get this copy on pe 0
   call get_copy(0, obs_ens_handle, k, obs_temp)
   ! Only pe 0 gets to write the sequence
   if(my_task_id() == 0) then
      ! Loop through the observations for this time
      do j = 1, obs_ens_handle%num_vars
         call get_obs_from_key(seq, keys(j), observation)
         call set_obs_values(observation, obs_temp(j:j), ens_offset + 2 * (k - 1))
         ! Store the observation into the sequence
         call set_obs(seq, observation, keys(j))
      end do
   endif
end do

! Output the qc global value
! First get this copy on pe 0
call get_copy(0, obs_ens_handle, OBS_GLOBAL_QC_COPY, obs_temp)
! Only pe 0 gets to write the observations for this time
if(my_task_id() == 0) then
   ! Loop through the observations for this time
   do j = 1, obs_ens_handle%num_vars
      call get_obs_from_key(seq, keys(j), observation)
      call set_qc(observation, obs_temp(j:j), DART_qc_index)
      ! Store the observation into the sequence
      call set_obs(seq, observation, keys(j))
   end do
endif

end subroutine obs_space_diagnostics

!-------------------------------------------------------------------------

function input_qc_ok(input_qc)

logical              :: input_qc_ok
real(r8), intent(in) :: input_qc

! Should eventually go in qc module

! Do checks on input_qc value with namelist control
! For test always return true for now
if(input_qc < ncep_qc_threshold) then
   input_qc_ok = .true.
else
   input_qc_ok = .false.
endif

end function input_qc_ok

!-------------------------------------------------------------------------

end program filter
