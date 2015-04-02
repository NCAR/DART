! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program perfect_model_obs

! Program to build an obs_sequence file from simulated observations.

use        types_mod,     only : r8, metadatalength
use    utilities_mod,     only : initialize_utilities, register_module, error_handler, &
                                 find_namelist_in_file, check_namelist_read,           &
                                 E_ERR, E_MSG, E_DBG, nmlfileunit, timestamp,          &
                                 do_nml_file, do_nml_term, logfileunit, &
                                 open_file, close_file, finalize_utilities
use time_manager_mod,     only : time_type, get_time, set_time, operator(/=), print_time,   &
                                 generate_seed
use obs_sequence_mod,     only : read_obs_seq, obs_type, obs_sequence_type,                 &
                                 get_obs_from_key, set_copy_meta_data, get_obs_def,         &
                                 get_time_range_keys, set_obs_values, set_qc, set_obs,      &
                                 write_obs_seq, get_num_obs, init_obs, assignment(=),       &
                                 static_init_obs_sequence, get_num_qc, read_obs_seq_header, &
                                 set_qc_meta_data, get_expected_obs, delete_seq_head,       &
                                 delete_seq_tail, destroy_obs, destroy_obs_sequence

use      obs_def_mod,     only : obs_def_type, get_obs_def_error_variance, get_obs_def_time
use    obs_model_mod,     only : move_ahead, advance_state, set_obs_model_trace
use  assim_model_mod,     only : static_init_assim_model, get_model_size,                    &
                                 aget_initial_condition, netcdf_file_type, init_diag_output, &
                                 aoutput_diagnostics, finalize_diag_output
   
use mpi_utilities_mod,    only : task_count, task_sync

use   random_seq_mod,     only : random_seq_type, init_random_seq, random_gaussian
use ensemble_manager_mod, only : init_ensemble_manager, write_ensemble_restart,              &
                                 end_ensemble_manager, ensemble_type, read_ensemble_restart, &
                                 get_my_num_copies, get_ensemble_time, prepare_to_write_to_vars,      &
                                 prepare_to_read_from_vars

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Module storage for message output
character(len=129) :: msgstring
integer            :: trace_level, timestamp_level

!-----------------------------------------------------------------------------
! Namelist with default values
!
logical  :: start_from_restart = .false.
logical  :: output_restart     = .false.
integer  :: async              = 0
logical  :: trace_execution    = .false.
logical  :: output_timestamps  = .false.
logical  :: silence            = .false.
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days     = 0
integer  :: init_time_seconds  = 0
! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days     = -1
integer  :: first_obs_seconds  = -1
integer  :: last_obs_days      = -1
integer  :: last_obs_seconds   = -1
integer  :: obs_window_days    = -1
integer  :: obs_window_seconds = -1
logical  :: output_forward_op_errors = .false.
integer  :: tasks_per_model_advance = 1
integer  :: output_interval = 1
integer  :: print_every_nth_obs = 0

character(len = 129) :: restart_in_file_name  = 'perfect_ics',     &
                        restart_out_file_name = 'perfect_restart', &
                        obs_seq_in_file_name  = 'obs_seq.in',      &
                        obs_seq_out_file_name = 'obs_seq.out',     &
                        adv_ens_command       = './advance_model.csh'

namelist /perfect_model_obs_nml/ start_from_restart, output_restart, async,         &
                                 init_time_days, init_time_seconds,                 &
                                 first_obs_days, first_obs_seconds,                 &
                                 last_obs_days,  last_obs_seconds, output_interval, &
                                 restart_in_file_name, restart_out_file_name,       &
                                 obs_seq_in_file_name, obs_seq_out_file_name,       &
                                 adv_ens_command, tasks_per_model_advance,          & 
                                 obs_window_days, obs_window_seconds, silence,      &
                                 trace_execution, output_timestamps,                &
                                 print_every_nth_obs, output_forward_op_errors


!------------------------------------------------------------------------------

! Doing this allows independent scoping for subroutines in main program file
call perfect_main()

!------------------------------------------------------------------------------

contains

subroutine perfect_main()

type(obs_sequence_type) :: seq
type(obs_type)          :: obs
type(obs_def_type)      :: obs_def
type(random_seq_type)   :: random_seq
type(ensemble_type)     :: ens_handle
type(netcdf_file_type)  :: StateUnit
type(time_type)         :: first_obs_time, last_obs_time
type(time_type)         :: window_time, curr_ens_time, next_ens_time

integer, allocatable    :: keys(:)
integer                 :: j, iunit, time_step_number, obs_seq_file_id
integer                 :: cnum_copies, cnum_qc, cnum_obs, cnum_max
integer                 :: additional_qc, additional_copies, forward_unit
integer                 :: ierr, io, istatus, num_obs_in_set, nth_obs
integer                 :: model_size, key_bounds(2), num_qc, last_key_used
integer                 :: seed

real(r8)                :: true_obs(1), obs_value(1), qc(1)

character(len=metadatalength) :: copy_meta_data(2), qc_meta_data, obs_seq_read_format
character(len=metadatalength) :: state_meta(1)

logical                 :: assimilate_this_ob, evaluate_this_ob, pre_I_format
logical                 :: all_gone

! Initialize all modules used that require it
call perfect_initialize_modules_used()

! Read the namelist entry
call find_namelist_in_file("input.nml", "perfect_model_obs_nml", iunit)
read(iunit, nml = perfect_model_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "perfect_model_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=perfect_model_obs_nml)
if (do_nml_term()) write(     *     , nml=perfect_model_obs_nml)

! set the level of output
call set_trace(trace_execution, output_timestamps, silence)

call trace_message('Perfect_model start')
call timestamp_message('Perfect_model start')

! Don't let this run with more than one task; just a waste of resource
if(task_count() > 1) then
   write(msgstring, *) 'Only use one mpi process here: ', task_count(), ' were requested'
   call error_handler(E_ERR, 'perfect_main', msgstring,  &
      source, revision, revdate)
endif

! Default to printing nothing
nth_obs = -1

call trace_message('Before setting up space for observations')
call timestamp_message('Before setting up space for observations')

! Find out how many data copies are in the obs_sequence 
call read_obs_seq_header(obs_seq_in_file_name, cnum_copies, cnum_qc, cnum_obs, cnum_max, &
   obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)

! First two copies of output will be truth and observation;
! Will overwrite first two existing copies in file if there are any
additional_copies = 2 - cnum_copies
if(additional_copies < 0) additional_copies = 0

! Want to have a qc field available in case forward op wont work
if(cnum_qc == 0) then
   additional_qc = 1
else
   additional_qc = 0
endif

! Read in definition part of obs sequence; expand to include observation and truth field
call read_obs_seq(obs_seq_in_file_name, additional_copies, additional_qc, 0, seq)

! Initialize an obs type variable
call init_obs(obs, cnum_copies + additional_copies, cnum_qc + additional_qc)

! Need metadata for added qc field
if(additional_qc == 1) then
   qc_meta_data = 'Quality Control'
   call set_qc_meta_data(seq, 1, qc_meta_data)
endif

! Need space to put in the obs_values in the sequence;
copy_meta_data(1) = 'observations'
copy_meta_data(2) = 'truth'
call set_copy_meta_data(seq, 1, copy_meta_data(1))
call set_copy_meta_data(seq, 2, copy_meta_data(2))

call timestamp_message('After  setting up space for observations')
call trace_message('After  setting up space for observations')

! Initialize the model now that obs_sequence is all set up
model_size = get_model_size()
write(msgstring,*)'Model size = ',model_size
call error_handler(E_MSG,'perfect_main',msgstring)

! Set up the ensemble storage and read in the restart file
call trace_message('Before reading in ensemble restart file')
call perfect_read_restart(ens_handle, model_size)
call trace_message('After  reading in ensemble restart file')

call trace_message('Before initializing output diagnostic file')
state_meta(1) = 'true state'
! Set up output of truth for state
StateUnit = init_diag_output('True_State', 'true state from control', 1, state_meta)
call trace_message('After  initializing output diagnostic file')

! Get the time of the first observation in the sequence
write(msgstring, *) 'total number of obs in sequence is ', get_num_obs(seq)
call error_handler(E_MSG,'perfect_main',msgstring)

num_qc = get_num_qc(seq)
write(msgstring, *) 'number of qc values is ',num_qc
call error_handler(E_MSG,'perfect_main',msgstring)

call trace_message('Before trimming obs seq if start/stop time specified')

! Need to find first obs with appropriate time, delete all earlier ones
if(first_obs_seconds > 0 .or. first_obs_days > 0) then
   first_obs_time = set_time(first_obs_seconds, first_obs_days)
   call delete_seq_head(first_obs_time, seq, all_gone)
   if(all_gone) then
      msgstring = 'All obs in sequence are before first_obs_days:first_obs_seconds'
      call error_handler(E_ERR,'perfect_main',msgstring,source,revision,revdate)
   endif
endif

last_key_used = -99

! Also get rid of observations past the last_obs_time if requested
if(last_obs_seconds >= 0 .or. last_obs_days >= 0) then
   last_obs_time = set_time(last_obs_seconds, last_obs_days)
   call delete_seq_tail(last_obs_time, seq, all_gone)
   if(all_gone) then
      msgstring = 'All obs in sequence are after last_obs_days:last_obs_seconds'
      call error_handler(E_ERR,'perfect_main',msgstring,source,revision,revdate)
   endif
endif

! Do verbose forward operator output if requested
if(output_forward_op_errors) forward_unit = open_file('forward_op_errors', 'formatted', 'append')

call trace_message('After  trimming obs seq if start/stop time specified')

! Time step number is used to do periodic diagnostic output
time_step_number = -1
window_time = set_time(0,0)

! Advance model to the closest time to the next available observations
AdvanceTime: do
   time_step_number = time_step_number + 1

   write(msgstring , '(A,I5)') 'Main evaluation loop, starting iteration', time_step_number
   call trace_message(' ', ' ', -1)
   call trace_message(msgstring, 'perfect_model_obs: ', -1)

   call trace_message('Before move_ahead checks time of data and next obs')

   ! Get the model to a good time to use a next set of observations
   call move_ahead(ens_handle, 1, seq, last_key_used, window_time, &
      key_bounds, num_obs_in_set, curr_ens_time, next_ens_time)
   if(key_bounds(1) < 0) then
      call trace_message('No more obs to evaluate, exiting main loop', 'perfect_model_obs:', -1)
      exit AdvanceTime
   endif

   call trace_message('After  move_ahead checks time of data and next obs')

   if (curr_ens_time /= next_ens_time) then
      call trace_message('Ready to run model to advance data ahead in time', 'perfect_model_obs:', -1)
      call print_ens_time(ens_handle, 'Ensemble data time before advance')
      call     trace_message('Before running model')
      call timestamp_message('Before running model', sync=.true.)

      call advance_state(ens_handle, 1, next_ens_time, async, &
         adv_ens_command, tasks_per_model_advance)

      call timestamp_message('After  running model', sync=.true.)
      call     trace_message('After  running model')
      call print_ens_time(ens_handle, 'Ensemble data time after  advance')
   else
      call trace_message('Model does not need to run; data already at required time', 'perfect_model_obs:', -1)
   endif

   ! Initialize a repeatable random sequence for perturbations
   seed = generate_seed(next_ens_time)
   call init_random_seq(random_seq,seed)

   call trace_message('Before setup for next group of observations')
   write(msgstring, '(A,I7)') 'Number of observations to be evaluated', &
      num_obs_in_set
   call trace_message(msgstring)
   call print_obs_time(seq, key_bounds(1), 'Time of first observation in window')
   call print_obs_time(seq, key_bounds(2), 'Time of last  observation in window')

   ! Allocate storage for observation keys for this part of sequence
   allocate(keys(num_obs_in_set))

   ! Get all the keys associated with this set of observations
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   call trace_message('After  setup for next group of observations')

   call prepare_to_read_from_vars(ens_handle)

   ! Output the true state to the netcdf file
   if((output_interval > 0) .and. &
      (time_step_number / output_interval * output_interval == time_step_number)) then

      call trace_message('Before updating truth diagnostics file')
      call aoutput_diagnostics(StateUnit, ens_handle%time(1), ens_handle%vars(:, 1), 1)
      call trace_message('After  updating truth diagnostics file')
   endif

   write(msgstring, '(A,I8,A)') 'Ready to evaluate up to', size(keys), ' observations'
   call trace_message(msgstring, 'perfect_model_obs:', -1)


   ! Compute the forward observation operator for each observation in set
   do j = 1, num_obs_in_set

      
      ! Some compilers do not like mod by 0, so test first.
      if (print_every_nth_obs > 0) nth_obs = mod(j, print_every_nth_obs)

      ! If requested, print out a message every Nth observation
      ! to indicate progress is being made and to allow estimates
      ! of how long the assim will take.
      if (nth_obs == 0) then
         write(msgstring, '(A,1x,I8,1x,A,I8)') 'Processing observation ', j, &
                                            ' of ', num_obs_in_set
         call trace_message(msgstring, 'perfect_model_obs:', -1)
         ! or if you want timestamps:
         !     call timestamp(msgstring, pos="debug")
      endif

      ! Compute the observations from the state
      call get_expected_obs(seq, keys(j:j), &
         1, ens_handle%vars(:, 1), ens_handle%time(1), .true., &
         true_obs(1:1), istatus, assimilate_this_ob, evaluate_this_ob)

      ! Get the observational error covariance (diagonal at present)
      ! Generate the synthetic observations by adding in error samples
      call get_obs_from_key(seq, keys(j), obs)
      call get_obs_def(obs, obs_def)

      ! If observation is not being evaluated or assimilated, skip it
      ! Ends up setting a 1000 qc field so observation is not used again.
      if(istatus == 0 .and. (assimilate_this_ob .or. evaluate_this_ob)) then
         obs_value(1) = random_gaussian(random_seq, true_obs(1), &
            sqrt(get_obs_def_error_variance(obs_def)))

         ! Set qc to 0 if none existed before
         if(cnum_qc == 0) then
            qc(1) = 0.0_r8
            call set_qc(obs, qc, 1)
         endif
      else
         obs_value(1) = true_obs(1)
         qc(1) = 1000.0_r8
         call set_qc(obs, qc, 1)

         ! FIXME: we could set different qc codes, like 1004, 1005, to
         ! indicate why this obs isn't being processed - separate failed 
         ! forward operators from those types not on the assim or eval lists.
         ! the values 4, 5, etc could match the dart QC values + 1000.

         ! if failed forward op logging requested, make sure we're
         ! only writing out obs with real errors and not those that
         ! end up in this code section because their type isn't in the namelist.
         if(output_forward_op_errors) then
            if ((istatus /= 0) .and. (assimilate_this_ob .or. evaluate_this_ob)) &
               write(forward_unit, *) keys(j), istatus
         endif
      endif

      call set_obs_values(obs, obs_value, 1)
      call set_obs_values(obs, true_obs, 2)

      ! Insert the observations into the sequence first copy
      call set_obs(seq, obs, keys(j))

   end do

   ! Deallocate the keys storage
   deallocate(keys)

   ! The last key used is updated to move forward in the observation sequence
   last_key_used = key_bounds(2)

end do AdvanceTime

! if logging errors, close unit
if(output_forward_op_errors) call close_file(forward_unit)

call trace_message('End of main evaluation loop, starting cleanup', 'perfect_model_obs:', -1)

! properly dispose of the diagnostics files
call trace_message('Before finalizing diagnostics file')
ierr = finalize_diag_output(StateUnit)
call trace_message('After  finalizing diagnostics file')

! Write out the sequence
call trace_message('Before writing output sequence file')
call write_obs_seq(seq, obs_seq_out_file_name)
call trace_message('After  writing output sequence file')

! Output a restart file if requested
if(output_restart) then
   call trace_message('Before writing state restart file')
   call write_ensemble_restart(ens_handle, restart_out_file_name, 1, 1, &
      force_single_file = .true.)
   call trace_message('After  writing state restart file')
endif

call trace_message('Before ensemble and obs memory cleanup')

!  Release storage for ensemble
call end_ensemble_manager(ens_handle)

! Free up the observation kind and obs sequence
call destroy_obs(obs)
call destroy_obs_sequence(seq)
call trace_message('After  ensemble and obs memory cleanup')

call trace_message('Perfect_model done')
call timestamp_message('Perfect_model done')

!call error_handler(E_MSG,'perfect_main','FINISHED',source,revision,revdate)

! closes the log file.
call finalize_utilities('perfect_model_obs')

end subroutine perfect_main

!=====================================================================

subroutine perfect_initialize_modules_used()

! Standard initialization (mpi not needed to use ensemble manager
! since we are enforcing that this run as a single task).
call initialize_utilities('perfect_model_obs')

! Initialize modules used that require it
call register_module(source,revision,revdate)

! Initialize the obs sequence module
call static_init_obs_sequence()
! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()

end subroutine perfect_initialize_modules_used

!---------------------------------------------------------------------

subroutine perfect_read_restart(ens_handle, model_size)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: model_size

type(time_type) :: time1
integer         :: secs, days

! First initialize the ensemble manager storage, only 1 copy for perfect
call init_ensemble_manager(ens_handle, 1, model_size, 1)

call prepare_to_write_to_vars(ens_handle)

! If not start_from_restart, use model to get ics for state and time
if(.not. start_from_restart) then
   call error_handler(E_MSG,'perfect_read_restart:', &
         'Using code in model_mod to initialize ensemble')

   call aget_initial_condition(ens_handle%time(1), ens_handle%vars(:, 1))
else
   call error_handler(E_MSG,'perfect_read_restart:', &
         'Reading in initial condition/restart data from file')

   ! Read in initial conditions from restart file
   if(init_time_days >= 0) then
      time1 = set_time(init_time_seconds, init_time_days)
      call read_ensemble_restart(ens_handle, 1, 1, &
         start_from_restart, restart_in_file_name, time1, force_single_file = .true.)

      write(msgstring, '(A)') 'By namelist control, ignoring time found in restart file.'
      call error_handler(E_MSG,'perfect_read_restart:',msgstring,source,revision,revdate)
      write(msgstring, '(A,I6,1X,I5)') 'Setting initial days, seconds to ', &
         init_time_days, init_time_seconds
      call error_handler(E_MSG,'perfect_read_restart:',msgstring,source,revision,revdate)

   else
      call read_ensemble_restart(ens_handle, 1, 1, &
         start_from_restart, restart_in_file_name, force_single_file = .true.)
   endif
endif

! Temporary print of initial model time
call get_time(ens_handle%time(1),secs,days)
write(msgstring, *) 'initial model time of perfect_model member (days,seconds) ',days,secs
call error_handler(E_DBG,'perfect_read_restart',msgstring,source,revision,revdate)

end subroutine perfect_read_restart

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

if (present(label)) then
   call error_handler(E_MSG,trim(label),trim(msg))
else
   call error_handler(E_MSG,'p_m_o trace:',trim(msg))
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

call timestamp(' '//trim(msg), pos='brief')

end subroutine timestamp_message

!-------------------------------------------------------------------------

subroutine print_ens_time(ens_handle, msg)

type(ensemble_type), intent(in) :: ens_handle
character(len=*), intent(in) :: msg

! Write message to stdout and log file.
type(time_type) :: mtime

if (trace_level <= 0) return

if (get_my_num_copies(ens_handle) < 1) return

call get_ensemble_time(ens_handle, 1, mtime)
call print_time(mtime, ' p_m_o trace: '//msg, logfileunit)
call print_time(mtime, ' p_m_o trace: '//msg)

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

call init_obs(obs, 0, 0)
call get_obs_from_key(seq, key, obs)
call get_obs_def(obs, obs_def)
mtime = get_obs_def_time(obs_def)
call print_time(mtime, ' p_m_o trace: '//msg, logfileunit)
call print_time(mtime, ' p_m_o trace: '//msg)
call destroy_obs(obs)

end subroutine print_obs_time

!-------------------------------------------------------------------------


end program perfect_model_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
