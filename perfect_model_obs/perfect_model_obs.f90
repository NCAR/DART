program perfect_model_obs

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

! Program to build a simple obs_sequence file for use in testing filters
! for spatial domains with one periodic dimension.

use types_mod,        only : r8, missing_r
use utilities_mod,    only : open_file, check_nml_error, file_exist, get_unit, close_file
use time_manager_mod, only : time_type, set_time, print_time, operator(/=)

use obs_sequence_mod, only : obs_sequence_type, init_obs_sequence, &
   add_obs_set, write_obs_sequence, read_obs_sequence, associate_def_list, &
   get_num_obs_sets, get_obs_set, get_obs_sequence_time, &
   get_num_obs_in_set, get_expected_obs, get_diag_obs_err_cov, &
   get_obs_values, get_num_obs_copies, inc_num_obs_copies, set_obs_values, &
   read_obs_sequence_def

use obs_def_mod, only : obs_def_type, init_obs_def
use obs_set_def_mod, only : obs_set_def_type, init_obs_set_def, add_obs
use obs_kind_mod, only : set_obs_kind
use set_def_list_mod, only : set_def_list_type, init_set_def_list, &
   add_to_list, write_set_def_list
use obs_set_mod, only : obs_set_type, init_obs_set, set_obs_set_time, write_obs_set, &
   get_obs_set_time, get_num_obs
use assim_model_mod, only : assim_model_type, static_init_assim_model, get_model_size, &
   get_initial_condition, get_model_state_vector, set_model_state_vector, &
   get_closest_state_time_to, advance_state, set_model_time, &
   get_model_time, init_diag_output, output_diagnostics, init_assim_model, &
   read_state_restart, write_state_restart, binary_restart_files
use random_seq_mod, only : random_seq_type, init_random_seq, &
   random_gaussian

use netcdf, only : NF90_close

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_def_type)      :: obs_def
type(set_def_list_type) :: set_def_list
type(obs_set_def_type)  :: obs_set_def
type(obs_set_type)      :: obs_set
type(time_type)         :: time1, time2
type(random_seq_type)   :: random_seq

integer :: i, j, obs_set_def_index, iunit, unit_out, num_obs_in_set
integer :: ierr, state_unit, StateUnit, io

! Need to set up namelists for controlling all of this mess, too!
integer :: model_size, num_obs_sets

type(assim_model_type) :: x(1)
real(r8), allocatable :: obs_err_cov(:), obs(:), true_obs(:)
character(len=129) :: copy_meta_data(2), file_name

!-----------------------------------------------------------------------------
! Namelist with default values
!
logical :: start_from_restart = .false., output_restart = .false.
integer :: async = 0
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer :: init_time_days = -1, init_time_seconds = -1, output_interval = 1
character(len = 129) :: restart_in_file_name = 'perfect_restart_in', &
                        restart_out_file_name = 'perfect_restart_out', &
                        obs_seq_in_file_name = 'obs_seq.in', &
                        obs_seq_out_file_name = 'obs_seq.out'

namelist /perfect_model_obs_nml/ async, obs_seq_in_file_name, &
   obs_seq_out_file_name, start_from_restart, output_restart, &
   restart_in_file_name, restart_out_file_name, init_time_days, init_time_seconds, &
   output_interval

!------------------------------------------------------------------------------

! Change output to diagnostic output block ...
write(*,*)'perfect_model_obs attributes:'
write(*,*)'   ',trim(adjustl(source))
write(*,*)'   ',trim(adjustl(revision))
write(*,*)'   ',trim(adjustl(revdate))
write(*,*)'   '
write(*,*)'    Reading input from input.nml namelist=perfect_model_obs_nml ...'

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file(file = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = perfect_model_obs_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'perfect_model_obs_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! Read in an observation sequence, only definitions part will be used (no data used)

iunit = get_unit()
open(file = obs_seq_in_file_name, unit = iunit)

! Just read in the definition part of the obs sequence
seq = read_obs_sequence_def(iunit)

! Set a time type for initial time if namelist inputs are not negative
if(init_time_days >= 0) then
   time1 = set_time(init_time_seconds, init_time_days)
else
   time1 = set_time(0, 0)
endif

! Initialize the model now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

write(*,*)'Model size = ',model_size

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
StateUnit  = init_diag_output(   'True_State', 'true state from control', 1, (/'true state'/))

! Initialize a repeatable random sequence for perturbations
call init_random_seq(random_seq)

num_obs_sets = get_num_obs_sets(seq)

! Need space to put in the obs_values in the sequence;
copy_meta_data(1) = 'observations'
copy_meta_data(2) = 'truth'

! Really need a way to read in just the definitions part of an obs_sequence
call inc_num_obs_copies(seq, 2, copy_meta_data)

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).
Advance: do i = 1, num_obs_sets
   call get_obs_sequence_time(seq, i, time1)
   write(*, *) ' '
   write(*, *) 'time of obs set ', i
   call print_time(time1)

! Figure out time to advance to
   time2 = get_closest_state_time_to(x(1), time1)
! Advance the state to this time; zero length advance is problem for B-grid so avoid
   if(time2 /= get_model_time(x(1))) call advance_state(x, 1, time2, async)

! Output the true state
   if(i / output_interval * output_interval == i) &
      call output_diagnostics(    StateUnit, x(1), 1)

! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)
   write(*, *) 'num_obs_in_set is ', num_obs_in_set

! Allocate storage for the observations and covariances
   allocate(obs_err_cov(num_obs_in_set), &
      true_obs(num_obs_in_set), obs(num_obs_in_set))

! Compute the observations from the state
   call get_expected_obs(seq, i, get_model_state_vector(x(1)), true_obs)
!   write(*, *) 'exact obs ', true_obs

! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

! Generate the synthetic observations by adding in error samples
   do j = 1, num_obs_in_set
      if (true_obs(j) /= missing_r) then
         obs(j) = random_gaussian(random_seq, true_obs(j), sqrt(obs_err_cov(j)))
      else
         obs(j) = true_obs(j)
      endif
   end do
!   write(*, *) 'obs with error added are ', obs

! Insert the observations into the sequence first copy
   call set_obs_values(seq, i, true_obs, 2)
   call set_obs_values(seq, i,      obs, 1)

! Deallocate the obs size storage
   deallocate(obs_err_cov, true_obs, obs)

end do Advance

! properly dispose of the diagnostics files

ierr = NF90_close(StateUnit)

! Write out the sequence
unit_out = get_unit()
open(file = obs_seq_out_file_name, unit = unit_out, status = 'replace')
call write_obs_sequence(unit_out, seq)

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

end program perfect_model_obs
