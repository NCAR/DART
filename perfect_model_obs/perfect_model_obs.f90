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
   get_next_obs, get_num_times, init_obs, assignment(=), static_init_obs_sequence, get_num_qc

use obs_def_mod,      only : obs_def_type, get_obs_def_time, get_obs_def_error_variance

use obs_model_mod,    only : get_expected_obs

use assim_model_mod, only  : assim_model_type, static_init_assim_model, get_model_size, &
   get_initial_condition, get_model_state_vector, set_model_state_vector, &
   get_closest_state_time_to, advance_state, set_model_time, get_model_time, &
   netcdf_file_type, init_diag_output, output_diagnostics, finalize_diag_output, &
   init_assim_model, read_state_restart, write_state_restart, binary_restart_files
use random_seq_mod,  only  : random_seq_type, init_random_seq, random_gaussian

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: obs_def
type(time_type)         :: time1, time2, next_time
type(random_seq_type)   :: random_seq

integer                 :: i, j, iunit
integer                 :: days, secs    ! for printing purposes only
type(netcdf_file_type)  :: StateUnit
integer                 :: ierr, io, istatus
integer                 :: model_size, key_bounds(2), num_obs, num_qc
integer, allocatable    :: keys(:)
logical                 :: is_there_one, out_of_range, is_this_last
real(r8)                :: true_obs(1), obs_value(1), rstatus(1,1)

type(assim_model_type)  :: x(1)
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

call initialize_utilities
call register_module(source, revision, revdate)
call error_handler(E_MSG,'perfect_model_obs','STARTING',source,revision,revdate)

! Begin by reading the namelist input
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
call init_obs(next_obs, 0, 0)

! Just read in the definition part of the obs sequence; expand to include observation and truth field
call static_init_obs_sequence()
call read_obs_seq(obs_seq_in_file_name, 2, 0, 0, seq)

! Set a time type for initial time if namelist inputs are not negative
if(init_time_days >= 0) then
   time1 = set_time(init_time_seconds, init_time_days)
else
   time1 = set_time(0, 0)
endif

! Initialize the model now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

write(msgstring,*)'Model size = ',model_size
call error_handler(E_MSG,'perfect_model_obs',msgstring,source,revision,revdate)

!------------------- Read restart if requested ----------------------

if(start_from_restart) then
   call init_assim_model(x(1))
   iunit = get_unit()

   if ( binary_restart_files ) then
      open(unit = iunit, file = restart_in_file_name, form = "unformatted")
      call read_state_restart(x(1), iunit, "unformatted")
   else
      open(unit = iunit, file = restart_in_file_name)
      call read_state_restart(x(1), iunit)
   endif

! If init_time_days an init_time_seconds are not < 0, set time to them
   if(init_time_days >= 0) call set_model_time(x(1) , time1)
   close(iunit)
!-----------------  Restart read in --------------------------------

else

!-----  Block to do cold start initialization ----------
! Initialize the control run
   call init_assim_model(x(1))
   call get_initial_condition(x(1))

! Set time to 0, 0 if none specified, otherwise to specified
   call set_model_time(x(1), time1)
!-------------------- End of cold start ensemble initialization block ------
endif

! Set up output of truth for state
StateUnit = init_diag_output('True_State', 'true state from control', 1, (/'true state'/))

! Initialize a repeatable random sequence for perturbations
call init_random_seq(random_seq)

! Need space to put in the obs_values in the sequence;
copy_meta_data(1) = 'observations'
copy_meta_data(2) = 'truth'
call set_copy_meta_data(seq, 1, copy_meta_data(1))
call set_copy_meta_data(seq, 2, copy_meta_data(2))

! Get the time of the first observation in the sequence
write(msgstring, *) 'number of obs in sequence is ', get_num_obs(seq)
call error_handler(E_MSG,'perfect_model_obs',msgstring,source,revision,revdate)

num_qc = get_num_qc(seq)
write(msgstring, *) 'number of qc values is ',num_qc
call error_handler(E_MSG,'perfect_model_obs',msgstring,source,revision,revdate)

is_there_one = get_first_obs(seq, obs)
! Test for no data at all here?
call get_obs_def(obs, obs_def)
next_time = get_obs_def_time(obs_def)

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).
Advance: do i = 1, get_num_times(seq)
! Verify that starting search at obs works here
   call get_obs_time_range(seq, next_time, next_time, key_bounds, num_obs, out_of_range, obs)
   allocate(keys(num_obs))
   call get_time_range_keys(seq, key_bounds, num_obs, keys)

   call get_time(next_time,secs,days)
   write(msgstring, *) 'time of obs set ', i,' is (d,s) =',days,secs
   call error_handler(E_MSG,'perfect_model_obs',msgstring,source,revision,revdate)

! Figure out time to advance to
   time2 = get_closest_state_time_to(x(1), next_time)
! Advance the state to this time; zero length advance is problem for B-grid so avoid
   if(time2 /= get_model_time(x(1))) call advance_state(x, 1, time2, async, adv_ens_command)

! Output the true state
   if(i / output_interval * output_interval == i) &
      call output_diagnostics( StateUnit, x(1), 1)

! How many observations in this set
   write(msgstring, *) 'num_obs_in_set is ', num_obs
   call error_handler(E_DBG,'perfect_model_obs',msgstring,source,revision,revdate)

! Can do this purely sequentially in perfect_model_obs for now if desired
   do j = 1, num_obs
! Compute the observations from the state
      call get_expected_obs(seq, keys(j:j), get_model_state_vector(x(1)), true_obs(1:1), istatus, rstatus(1:1,1:1))

! Get the observational error covariance (diagonal at present)
      call get_obs_from_key(seq, keys(j), obs)
      call get_obs_def(obs, obs_def)

! Generate the synthetic observations by adding in error samples

      if(istatus == 0) then
         obs_value(1) = random_gaussian(random_seq, true_obs(1), sqrt(get_obs_def_error_variance(obs_def)))
      else
         obs_value(1) = true_obs(1)
      endif

      if (num_qc > 0) call set_qc(obs, rstatus(1,:), 1)

      call set_obs_values(obs, obs_value, 1)
      call set_obs_values(obs, true_obs, 2)

! Insert the observations into the sequence first copy
      call set_obs(seq, obs, keys(j))

   end do

! Get the next time (if any) in the sequence
   call get_next_obs(seq, obs, next_obs, is_this_last)
   if(is_this_last) exit
   call get_obs_def(next_obs, obs_def)
   next_time = get_obs_def_time(obs_def)

! Deallocate the obs size storage
   deallocate(keys)

end do Advance

! properly dispose of the diagnostics files

ierr = finalize_diag_output(StateUnit)

! Write out the sequence
call write_obs_seq(seq, obs_seq_out_file_name)

! Output a restart file if requested
if(output_restart) then
   iunit = get_unit()
   if ( binary_restart_files ) then
      open(unit = iunit, file = restart_out_file_name, form = "unformatted", status = 'replace')
      call write_state_restart(x(1), iunit, "unformatted")
   else
      open(unit = iunit, file = restart_out_file_name, status = 'replace')
      call write_state_restart(x(1), iunit)
   endif
   close(iunit)
endif

call error_handler(E_MSG,'perfect_model_obs','FINISHED',source,revision,revdate)

! closes the log file.
call timestamp(string1=source,string2=revision,string3=revdate,pos='end')
 
end program perfect_model_obs
