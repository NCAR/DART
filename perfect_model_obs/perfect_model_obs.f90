! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program perfect_model_obs

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

! Program to build a simple obs_sequence file for use in testing filters
! for spatial domains with one periodic dimension.

use types_mod,        only : r8
use utilities_mod,    only : open_file, check_nml_error, file_exist, get_unit, close_file, &
                             initialize_utilities, register_module, error_handler, &
                             E_ERR, E_WARN, E_MSG, E_DBG, logfileunit, timestamp
use time_manager_mod, only : time_type, set_time, get_time, operator(/=), operator(*), operator(+)

use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
   get_obs_from_key, set_copy_meta_data, get_copy_meta_data, get_obs_def, get_obs_time_range, &
   get_time_range_keys, set_obs_values, set_qc, set_obs, write_obs_seq, get_num_obs, &
   get_next_obs, get_num_times, init_obs, assignment(=), static_init_obs_sequence, get_num_qc, &
   get_num_copies

use obs_def_mod,      only : obs_def_type, get_obs_def_time, get_obs_def_error_variance

use obs_model_mod,    only : get_expected_obs, move_ahead

use assim_model_mod, only  : static_init_assim_model, get_model_size, &
   aget_initial_condition, get_model_state_vector, set_model_state_vector, &
   set_model_time, get_model_time, &
   netcdf_file_type, init_diag_output, aoutput_diagnostics, finalize_diag_output, &
   init_assim_model, read_state_restart, awrite_state_restart, &
   open_restart_read, open_restart_write, close_restart
use random_seq_mod,  only  : random_seq_type, init_random_seq, random_gaussian
use ensemble_manager_mod, only : init_ensemble_manager, get_ensemble_member, &
   put_ensemble_member, end_ensemble_manager, ensemble_type

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: obs
type(obs_def_type)      :: obs_def
type(time_type)         :: time1, ens_time
type(random_seq_type)   :: random_seq
type(ensemble_type)     :: ens_handle

integer                 :: i, j, iunit

type(netcdf_file_type)  :: StateUnit
integer                 :: ierr, io, istatus, num_obs_in_set
integer                 :: model_size, key_bounds(2), num_qc, last_key_used
integer, allocatable    :: keys(:)
real(r8)                :: true_obs(1), obs_value(1), qc(1)

real(r8), allocatable   :: ens(:)
character(len=129)      :: copy_meta_data(2), msgstring

!-----------------------------------------------------------------------------
! Namelist with default values
!
logical :: start_from_restart = .false., output_restart = .false.
integer :: async = 0
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer :: init_time_days = 0, init_time_seconds = 0, output_interval = 1
character(len = 129) :: restart_in_file_name  = 'perfect_ics',     &
                        restart_out_file_name = 'perfect_restart', &
                        obs_seq_in_file_name  = 'obs_seq.in',      &
                        obs_seq_out_file_name = 'obs_seq.out',     &
                        adv_ens_command       = './advance_ens.csh'


! adv_ens_command  == 'qsub advance_ens.csh' -> system call advances ensemble by
!                                               qsub submission of a batch job
!                                               -l num_nodes can be inserted after qsub
!                  == './advance_ens.csh'    -> advance ensemble using a script which
!                                               explicitly distributes ensemble among nodes
! advance_ens.csh is currently written to handle both batch submissions (qsub) and
!                 non-batch executions.

namelist /perfect_model_obs_nml/ async, adv_ens_command, obs_seq_in_file_name, &
   obs_seq_out_file_name, start_from_restart, output_restart, &
   restart_in_file_name, restart_out_file_name, init_time_days, init_time_seconds, &
   output_interval

!------------------------------------------------------------------------------

! Delete the semaphore files that are used for parallel version 3
call system('rm -f go_advance_model go_end_filter go_assim_regions')

call perfect_initialize_modules_used()

! call perfect_read_namelist()
! Begin by reading the namelist input
! Intel 8.0 quirk that the subroutine does not compile.
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = perfect_model_obs_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'perfect_model_obs_nml')
   enddo
 11 continue
   call close_file(iunit)
endif
write(logfileunit,nml=perfect_model_obs_nml)

! Initialize the two obs type variables
call init_obs(obs, 0, 0)

! Just read in the definition part of the obs sequence; expand to include observation and truth field
call read_obs_seq(obs_seq_in_file_name, 2, 0, 0, seq)

! Want to have error exit if input file has any obs values in it
if(get_num_copies(seq) /= 2) then
   write(msgstring, *) 'Input obs_sequence file should not have any copies of data associated with it'
   call error_handler(E_ERR, 'perfect_model_obs', msgstring, source, revision, revdate)
endif

! Need space to put in the obs_values in the sequence;
copy_meta_data(1) = 'observations'
copy_meta_data(2) = 'truth'
call set_copy_meta_data(seq, 1, copy_meta_data(1))
call set_copy_meta_data(seq, 2, copy_meta_data(2))

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time()

! Initialize the model now that obs_sequence is all set up
model_size = get_model_size()
! Allocate storage for doing advance with ensemble based tools
allocate(ens(model_size))

write(msgstring,*)'Model size = ',model_size
call error_handler(E_MSG,'perfect_model_obs',msgstring,source,revision,revdate)

call perfect_read_restart()

! Set up output of truth for state
StateUnit = init_diag_output('True_State', 'true state from control', 1, (/'true state'/))

! Initialize a repeatable random sequence for perturbations
call init_random_seq(random_seq)

! Get the time of the first observation in the sequence
write(msgstring, *) 'number of obs in sequence is ', get_num_obs(seq)
call error_handler(E_MSG,'perfect_model_obs',msgstring,source,revision,revdate)

num_qc = get_num_qc(seq)
write(msgstring, *) 'number of qc values is ',num_qc
call error_handler(E_MSG,'perfect_model_obs',msgstring,source,revision,revdate)

! Start out with no previously used observations
last_key_used = -99

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).
AdvanceTime: do

   ! Get the model to a good time to use a next set of observations
   call move_ahead(ens_handle, 1, model_size, seq, last_key_used, &
      key_bounds, num_obs_in_set, async, adv_ens_command)
   if(key_bounds(1) < 0) exit AdvanceTime

   ! Allocate storage for the ensemble priors for this number of observations
   allocate(keys(num_obs_in_set))

   ! Get all the keys associated with this set of observations
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

! Output the true state
   if(i / output_interval * output_interval == i) then
      call get_ensemble_member(ens_handle, 1, ens, time1)
      call aoutput_diagnostics(StateUnit, time1, ens, 1)
   endif

! How many observations in this set
   write(msgstring, *) 'num_obs_in_set is ', num_obs_in_set
   call error_handler(E_DBG,'perfect_model_obs',msgstring,source,revision,revdate)

! Can do this purely sequentially in perfect_model_obs for now if desired
   do j = 1, num_obs_in_set
! Compute the observations from the state
      call get_ensemble_member(ens_handle, 1, ens, time1)
      call get_expected_obs(seq, keys(j:j), ens, true_obs(1:1), istatus)
      if(istatus /= 0) qc(1) = 1000.

! Get the observational error covariance (diagonal at present)
      call get_obs_from_key(seq, keys(j), obs)
      call get_obs_def(obs, obs_def)

! Generate the synthetic observations by adding in error samples

      if(istatus == 0) then
         obs_value(1) = random_gaussian(random_seq, true_obs(1), sqrt(get_obs_def_error_variance(obs_def)))
      else
         obs_value(1) = true_obs(1)
      endif

      if (num_qc > 0) call set_qc(obs, qc, 1)

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

! Send a message to the asynchronous version 3 that all is done
! Must be done always because this also terminates option 3 for assim_tools!
call system('echo a > go_end_filter')

! properly dispose of the diagnostics files

ierr = finalize_diag_output(StateUnit)

! Write out the sequence
call write_obs_seq(seq, obs_seq_out_file_name)

! Output a restart file if requested
if(output_restart) then
   iunit = open_restart_write(restart_out_file_name)
   call get_ensemble_member(ens_handle, 1, ens, time1)
   call awrite_state_restart(time1, ens, iunit)
   call close_restart(iunit)
endif

call error_handler(E_MSG,'perfect_model_obs','FINISHED',source,revision,revdate)

! closes the log file.
call timestamp(string1=source,string2=revision,string3=revdate,pos='end')

deallocate(ens)

contains

!=====================================================================

subroutine perfect_initialize_modules_used()

! Initialize modules used that require it
call initialize_utilities
call register_module(source,revision,revdate)
call error_handler(E_MSG,'perfect_model_obs','STARTING',source,revision,revdate)

! Initialize the obs sequence module
call static_init_obs_sequence()
! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()

end subroutine perfect_initialize_modules_used

!---------------------------------------------------------------------

subroutine perfect_read_namelist()
!
! Intel 8.0 compiler chokes on any I/O in this subroutine.
! Consequently, the code block has been duplicated in the main program.
! There is an error report (28Jun2004) to fix this.
!
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
    do while(ierr /= 0)
!      read(iunit, nml = perfect_model_obs_nml, iostat = io, end = 11)
       ierr = check_nml_error(io, 'perfect_model_obs_nml')
    enddo
 11 continue
   call close_file(iunit)
endif
!write(logfileunit,nml=perfect_model_obs_nml)

end subroutine perfect_read_namelist

!---------------------------------------------------------------------

subroutine filter_set_initial_time()

if(init_time_days >= 0) then
   time1 = set_time(init_time_seconds, init_time_days)
else
   time1 = set_time(0, 0)
endif

end subroutine filter_set_initial_time

!---------------------------------------------------------------------

subroutine perfect_read_restart()

! Read restart if requested
if(start_from_restart) then
   if(init_time_days >= 0) then
      call init_ensemble_manager(ens_handle, 1, model_size, restart_in_file_name, time1)
   else
      call init_ensemble_manager(ens_handle, 1, model_size, restart_in_file_name)
   endif

else

   ! Initialize manager but nothing to read in
   call init_ensemble_manager(ens_handle, 1, model_size)
   ! Block to do cold start initialization
   call aget_initial_condition(time1, ens)

   call put_ensemble_member(ens_handle, 1, ens, time1)
   ! End of cold start ensemble initialization block
endif

end subroutine perfect_read_restart

!---------------------------------------------------------------------
!---------------------------------------------------------------------
 
end program perfect_model_obs
