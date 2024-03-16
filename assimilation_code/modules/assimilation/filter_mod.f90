! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module filter_mod

!------------------------------------------------------------------------------
use types_mod,             only : r8, i8, missing_r8, metadatalength, MAX_NUM_DOMS, MAX_FILES

use options_mod,           only : set_missing_ok_status

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
use obs_model_mod,         only : set_obs_model_trace, advance_model

use ensemble_manager_mod,  only : end_ensemble_manager, ensemble_type, &
                                  compute_copy_mean, compute_copy_mean_sd, &
                                  duplicate_state_copies, init_ensemble_manager
                                  
use adaptive_inflate_mod,  only : do_ss_inflate, &
                                  inflate_ens, adaptive_inflate_init,                 &
                                  adaptive_inflate_type, log_inflation_info,          &
                                  do_rtps_inflate, set_inflate_flavor,                &
                                  NO_INFLATION

use mpi_utilities_mod,     only : my_task_id, task_count, iam_task0

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

use filter_io_diag_mod,    only : ens_copies_type, create_ensemble_from_single_file, &
                                  init_state_ens, read_state_and_inflation,      &
                                  output_diagnostics, init_output_file_info

!------------------------------------------------------------------------------

implicit none
private

public :: filter_set_initial_time, filter_main, trim_obs_sequence

character(len=*), parameter :: source = 'filter_mod.f90'

! Some convenient global storage items
character(len=512)      :: msgstring

integer :: trace_level, timestamp_level

! Set to true if forecasts will be made
logical                 :: has_cycling = .false. ! filter will advance the model

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

! Eventually namelist controlled?
integer, parameter :: num_lags = 2

! Ensembles of state, lagged state, and forward operators
type(ensemble_type)         :: state_ens_handle
type(ensemble_type)         :: lag_ens_handle(num_lags)
type(ensemble_type)         :: obs_fwd_op_ens_handle
! Metadata for ensemble copies
type(ens_copies_type) :: ens_copies
type(ens_copies_type) :: lag_copies

! Prior and posterior inflation metadata
type(adaptive_inflate_type) :: prior_inflate, post_inflate
! No inflation for smoothers, but still need a structure that indicates that for now
type(adaptive_inflate_type) :: lag_inflate

! Observation sequence and metadata for copies
type(obs_sequence_type)     :: seq
type(forward_op_info_type)  :: forward_op_ens_info

! Structures to control output for each possible point in algorithms
type(file_info_type) :: file_info_read
type(file_info_type) :: file_info_forecast
type(file_info_type) :: file_info_preassim
type(file_info_type) :: file_info_postassim
type(file_info_type) :: file_info_analysis
type(file_info_type) :: file_info_output
type(file_info_type) :: lag_info_preassim(num_lags)
type(file_info_type) :: lag_info_postassim(num_lags)

type(time_type)         :: time1, curr_ens_time, next_ens_time
integer                 :: time_step_number, num_obs_in_set, i
integer                 :: key_bounds(2)
integer,    allocatable :: keys(:)
logical                 :: read_time_from_file

call filter_initialize_modules_used() ! static_init_model called in here

! Read the nameslist and adjust control variables
call process_namelist_entries()

! Initialize the adaptive inflation module
call init_inflation_options(prior_inflate, post_inflate)

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time(init_time_days, init_time_seconds, time1, read_time_from_file)

! Initialize the ensemble manager and the ens_copies index information
call init_state_ens(state_ens_handle, ens_copies, ens_size, output_mean, output_sd, &
   do_prior_inflate, do_posterior_inflate, post_inflate, num_output_state_members, distributed_state)

! Initialize all lags for fixed lag smoother
do i = 1, num_lags
   call init_state_ens(lag_ens_handle(i), lag_copies, ens_size, output_mean, output_sd, &
      do_prior_inflate, do_posterior_inflate, post_inflate, num_output_state_members, distributed_state)
end do

! for now, assume that we only allow cycling if single_file_out is true.
! code in this call needs to know how to initialize the output files.
call init_output_file_info('output', state_ens_handle, ens_copies, &
   file_info_output, output_state_files, output_state_file_list, single_file_out,      &
   has_cycling, do_prior_inflate, do_posterior_inflate)

! Initialize and then read the input file
call read_state_and_inflation(ens_copies, state_ens_handle, single_file_in, &
   perturb_from_single_instance, has_cycling, input_state_files, input_state_file_list, &
   prior_inflate, do_prior_inflate, prior_inflate_from_restart,                          &
   post_inflate, do_posterior_inflate, posterior_inflate_from_restart,              &
   file_info_read, read_time_from_file, time1)

if (perturb_from_single_instance) &
   call create_ensemble_from_single_file(state_ens_handle, ens_size, &
      perturbation_amplitude, time1, allow_missing_clm)

! Initialize the obs_sequence and metadata; every pe gets a copy for now
call filter_setup_obs_sequence(forward_op_ens_info, seq, num_output_obs_members, &
   obs_sequence_in_name, compute_posterior)

! Remove observations that are too early or too late to be assimilated 
call trim_obs_sequence(seq, key_bounds)

! Initiazlize lagged states
do i = 1, num_lags
   call duplicate_state_copies(state_ens_handle, lag_ens_handle(i), .true.)
end do

! Loop through timesteps until observations are exhausted
AdvanceTime : do time_step_number = 0, huge(time_step_number)

   ! Advance the model to make the window include the next available observation.
   call advance_model(state_ens_handle, ens_size, seq, key_bounds, num_obs_in_set, &
      curr_ens_time, next_ens_time, async, adv_ens_command, tasks_per_model_advance, &
      file_info_output, file_info_read)

   ! No more observations available so exit the time loop
   if(key_bounds(1) < 0) exit AdvanceTime

   call filter_update(state_ens_handle, ens_copies, prior_inflate, post_inflate, seq, &
      forward_op_ens_info, num_obs_in_set, keys, key_bounds, compute_posterior, has_cycling, &
      output_forecast_diags,  'forecast',  file_info_forecast,  &
      output_preassim_diags,  'preassim',  file_info_preassim,  &
      output_postassim_diags, 'postassim', file_info_postassim, &
      output_analysis_diags,  'analysis',  file_info_analysis,  &
      obs_fwd_op_ens_handle)

   ! Propoagate the lags
   call duplicate_state_copies(state_ens_handle, lag_ens_handle(1), .true.)
   do i = 2, num_lags
      call duplicate_state_copies(lag_ens_handle(i-1), lag_ens_handle(i), .true.)
   end do

   ! MOVE THIS LINE UP OUT OF ADVANCE LOOP 
   call set_inflate_flavor(lag_inflate, NO_INFLATION) 

   ! Let current observations impact each lag
   do i = 1, num_lags
      call smoother_update(lag_ens_handle(i), lag_copies, lag_inflate, seq, forward_op_ens_info, &
         obs_fwd_op_ens_handle, num_obs_in_set, keys, key_bounds, i, has_cycling,         &
         output_preassim_diags,  'lag_preassim',  lag_info_preassim(i),             &
         output_postassim_diags, 'lag_postassim', lag_info_postassim(i))
   end do

   deallocate(keys)

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

function output_diag_now(output_interval, time_step_number)

logical             :: output_diag_now
integer, intent(in) :: output_interval, time_step_number

output_diag_now = ((output_interval > 0) .and. &  
   (time_step_number / output_interval * output_interval == time_step_number)) 

end function output_diag_now

!------------------------------------------------------------------

subroutine init_inflation_options(prior_inflate, post_inflate)

type(adaptive_inflate_type), intent(inout) :: prior_inflate, post_inflate

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

end subroutine init_inflation_options

!------------------------------------------------------------------

subroutine process_namelist_entries()

integer :: io, iunit

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
   
! 'has_cycling' set to 'single_file_out'; only allowing cycling if writing to a single file.
has_cycling = single_file_out

! Can't output more ensemble members than exist
if(num_output_obs_members   > ens_size) num_output_obs_members   = ens_size

! Don't currently support number of processes > model_size
if(task_count() > get_model_size()) then
   write(msgstring, *) 'number of MPI processes = ', task_count(), &
                       ' while model size = ', get_model_size()
   call error_handler(E_ERR,'filter_main', &
      'Cannot have number of processes > model size' ,source, text2=msgstring)
endif

end subroutine process_namelist_entries

!------------------------------------------------------------------

subroutine trim_obs_sequence(seq, key_bounds)

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(out)   :: key_bounds(2)

type(time_type) :: first_obs_time, last_obs_time
logical :: all_gone

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
key_bounds(1:2) = -99

! Also get rid of observations past the last_obs_time if requested
if(last_obs_seconds >= 0 .or. last_obs_days >= 0) then
   last_obs_time = set_time(last_obs_seconds, last_obs_days)
   call delete_seq_tail(last_obs_time, seq, all_gone)
   if(all_gone) then
      msgstring = 'All obs in sequence are after last_obs_days:last_obs_seconds'
      call error_handler(E_ERR,'filter_main',msgstring,source)
   endif
endif

end subroutine trim_obs_sequence
!------------------------------------------------------------------

! Does prior inflation, forward operators, assimilation, posterior inflation 

subroutine filter_update(state_ens_handle, ens_copies, prior_inflate, post_inflate, &
   seq, forward_op_ens_info, num_obs_in_set, keys, key_bounds, compute_posterior, has_cycling, &
   output_forecast_diags, filename_forecast, file_info_forecast, &
   output_preassim_diags, filename_preassim, file_info_preassim, &
   output_postassim_diags, filename_postassim, file_info_postassim, &
   output_analysis_diags, filename_analysis, file_info_analysis, &
   save_fwd_op_ens_handle)

type(ensemble_type),         intent(inout) :: state_ens_handle
type(ens_copies_type),       intent(inout) :: ens_copies
type(adaptive_inflate_type), intent(inout) :: prior_inflate, post_inflate
type(forward_op_info_type),  intent(inout) :: forward_op_ens_info
type(obs_sequence_type),     intent(inout) :: seq
integer,                     intent(in)    :: num_obs_in_set
integer, allocatable,        intent(inout) :: keys(:)
integer,                     intent(in)    :: key_bounds(2)
logical,                     intent(in)    :: compute_posterior
logical,                     intent(in)    :: has_cycling
character(len=*),            intent(in)    :: filename_forecast, filename_preassim
character(len=*),            intent(in)    :: filename_postassim, filename_analysis
type(file_info_type),        intent(inout) :: file_info_forecast, file_info_preassim
type(file_info_type),        intent(inout) :: file_info_postassim, file_info_analysis
logical,                     intent(in)    :: output_forecast_diags, output_preassim_diags
logical,                     intent(in)    :: output_postassim_diags, output_analysis_diags
type(ensemble_type),         intent(inout) :: save_fwd_op_ens_handle
!----------
type(ensemble_type)     :: obs_fwd_op_ens_handle, qc_ens_handle
integer                 :: ens_size

ens_size = ens_copies%ens_size
! Write out forecast diagnostic file(s). 
if(output_forecast_diags) &
   call output_diagnostics(filename_forecast, state_ens_handle, ens_copies,            &
      file_info_forecast, single_file_out, has_cycling, output_mean, output_sd, &
      output_members, do_prior_inflate, do_posterior_inflate)
   
! Apply prior inflation 
if(do_ss_inflate(prior_inflate)) &
   call filter_ensemble_inflate(state_ens_handle, ens_copies%PRIOR_INF_COPY, prior_inflate, &
      ens_copies%ENS_MEAN_COPY)

! if relaxation-to-prior-spread inflation, save prior spread in copy RTPS_PRIOR_SPREAD_COPY
if ( do_rtps_inflate(post_inflate) ) &
   call compute_copy_mean_sd(state_ens_handle, 1, ens_copies%ens_size, &
      ens_copies%ENS_MEAN_COPY, ens_copies%RTPS_PRIOR_SPREAD_COPY)

! Write out preassim diagnostic files if requested.
if(output_preassim_diags) &
   call output_diagnostics(filename_preassim, state_ens_handle, ens_copies, &
      file_info_preassim, single_file_out, has_cycling, output_mean, output_sd,           &
      output_members, do_prior_inflate, do_posterior_inflate)

! Compute the forward operators and fill data structures
call forward_operators(forward_op_ens_info, state_ens_handle, obs_fwd_op_ens_handle, &
   qc_ens_handle, seq, ens_size, num_groups, num_obs_in_set, keys, key_bounds,       &
   num_output_obs_members, compute_posterior, isprior = .true.)

! Need to save the forward operator ensemble for smoothers
call init_ensemble_manager(save_fwd_op_ens_handle, forward_op_ens_info%TOTAL_COPIES, &
   int(num_obs_in_set, i8), 1, transpose_type_in = 2)
call duplicate_state_copies(obs_fwd_op_ens_handle, save_fwd_op_ens_handle, .true.)

call filter_assim(state_ens_handle, ens_copies, obs_fwd_op_ens_handle, forward_op_ens_info, &
   seq, keys, num_groups, prior_inflate, ens_copies%PRIOR_INF_COPY, ens_copies%PRIOR_INF_SD_COPY, &
   inflate_only = .false.)

! Write out postassim diagnostic files if requested.  This contains the assimilated ensemble 
! JLA DEVELOPMENT: This used to output the damped inflation. NO LONGER.
if(output_postassim_diags) &
   call output_diagnostics(filename_postassim, state_ens_handle, ens_copies, &
      file_info_postassim, single_file_out, has_cycling, output_mean, output_sd,           &
      output_members, do_prior_inflate, do_posterior_inflate)

! This block applies posterior inflation including RTPS if selected
if(do_ss_inflate(post_inflate) .or. do_rtps_inflate(post_inflate))             &
   call filter_ensemble_inflate(state_ens_handle, ens_copies%POST_INF_COPY, post_inflate, &
      ens_copies%ENS_MEAN_COPY, ens_copies%RTPS_PRIOR_SPREAD_COPY, ens_copies%ENS_SD_COPY)

! this block recomputes the expected obs values for the obs_seq.final file
if (compute_posterior) then
   ! Compute the forward operators and fill data structures
   call forward_operators(forward_op_ens_info, state_ens_handle, obs_fwd_op_ens_handle, &
      qc_ens_handle, seq, ens_size, num_groups, num_obs_in_set, keys, key_bounds,       &
      num_output_obs_members, compute_posterior, isprior = .false.)
else
   ! Collect any updated QC values that may have been set in the assimilation loop
   call obs_space_sync_QCs(forward_op_ens_info, obs_fwd_op_ens_handle, &
      seq, keys, num_obs_in_set)
endif
      
! Compute the adaptive state space posterior inflation
if(do_ss_inflate(post_inflate) .and. ( .not. do_rtps_inflate(post_inflate)) )                  &
   call filter_assim(state_ens_handle, ens_copies, obs_fwd_op_ens_handle, forward_op_ens_info, &
      seq, keys, num_groups, post_inflate, &
      ens_copies%POST_INF_COPY, ens_copies%POST_INF_SD_COPY, inflate_only = .true.)
      
! Free up all the allocated space associated with obs ensemble
call end_ensemble_manager(obs_fwd_op_ens_handle)
call end_ensemble_manager(qc_ens_handle)
      
! Write out analysis diagnostic files if requested. 
if(output_analysis_diags) &
   call output_diagnostics(filename_analysis, state_ens_handle, ens_copies, &
      file_info_analysis, single_file_out, has_cycling, output_mean, output_sd,           &
      output_members, do_prior_inflate, do_posterior_inflate)

end subroutine filter_update

!------------------------------------------------------------------

! Does prior inflation, forward operators, assimilation, posterior inflation 

subroutine smoother_update(state_ens_handle, ens_copies, inflate, seq, forward_op_ens_info, &
   obs_fwd_op_ens_handle, num_obs_in_set, keys, key_bounds, lag, has_cycling,               &
   output_preassim_diags,  filename_preassim,  file_info_preassim,                          &
   output_postassim_diags, filename_postassim, file_info_postassim)

type(ensemble_type),         intent(inout) :: state_ens_handle
type(ens_copies_type),       intent(inout) :: ens_copies
type(adaptive_inflate_type), intent(inout) :: inflate
type(obs_sequence_type),     intent(inout) :: seq
type(forward_op_info_type),  intent(inout) :: forward_op_ens_info
type(ensemble_type),         intent(inout) :: obs_fwd_op_ens_handle
integer,                     intent(in)    :: num_obs_in_set
integer,                     intent(in)    :: keys(:)
integer,                     intent(in)    :: key_bounds(2)
integer,                     intent(in)    :: lag
logical,                     intent(in)    :: has_cycling
character(len=*),            intent(in)    :: filename_preassim, filename_postassim
type(file_info_type),        intent(inout) :: file_info_preassim, file_info_postassim
logical,                     intent(in)    :: output_preassim_diags, output_postassim_diags
!----------

! Storage for filenames with lag extension
character(len=512)      :: filename

! Generate the preassim file name for this lag
write(filename, '(A, I1)') trim(filename_preassim), lag
  
! Write out preassim diagnostic files if requested.
if(output_preassim_diags) &
   call output_diagnostics(filename, state_ens_handle, ens_copies, &
      file_info_preassim, single_file_out, has_cycling, output_mean, output_sd,           &
      output_members, .false., .false.)

call filter_assim(state_ens_handle, ens_copies, obs_fwd_op_ens_handle, forward_op_ens_info, &
   seq, keys, num_groups, inflate, ens_copies%PRIOR_INF_COPY, ens_copies%PRIOR_INF_SD_COPY, &
   inflate_only = .false.)

! Generate the preassim file name for this lag
write(filename, '(A, I1)') trim(filename_postassim), lag
! Write out postassim diagnostic files if requested.  This contains the assimilated ensemble 
! JLA DEVELOPMENT: This used to output the damped inflation. NO LONGER.
if(output_postassim_diags) &
   call output_diagnostics(filename, state_ens_handle, ens_copies, &
      file_info_postassim, single_file_out, has_cycling, output_mean, output_sd,           &
      output_members, .false., .false.)

end subroutine smoother_update

!------------------------------------------------------------------

end module filter_mod


