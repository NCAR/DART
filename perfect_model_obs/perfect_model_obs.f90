program perfect_model_obs

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

! Program to build a simple obs_sequence file for use in testing filters
! for spatial domains with one periodic dimension.

use types_mod
use utilities_mod,    only : open_file
use time_manager_mod, only : time_type, set_time, print_time, operator(/=)

use obs_sequence_mod, only : init_obs_sequence, obs_sequence_type, &
   add_obs_set, write_obs_sequence, read_obs_sequence, associate_def_list, &
   get_num_obs_sets, get_obs_set, get_obs_sequence_time, &
   get_num_obs_in_set, get_expected_obs, get_diag_obs_err_cov, &
   get_obs_values, get_num_obs_copies, inc_num_obs_copies, set_obs_values, &
   read_obs_sequence_def

use obs_def_mod, only : obs_def_type, init_obs_def
use obs_set_def_mod, only : obs_set_def_type, init_obs_set_def, add_obs
use obs_kind_mod, only : set_obs_kind
use location_mod, only : set_location
use set_def_list_mod, only : set_def_list_type, init_set_def_list, &
   add_to_list, write_set_def_list
use obs_set_mod, only : init_obs_set, obs_set_type, set_obs_set_time, write_obs_set, &
   get_obs_set_time, get_num_obs
use assim_model_mod, only : assim_model_type, static_init_assim_model, get_model_size, &
   get_initial_condition, get_model_state_vector, set_model_state_vector, &
   get_closest_state_time_to, advance_state, set_model_time, &
   get_model_time, init_diag_output, output_diagnostics, init_assim_model
use random_seq_mod, only : random_seq_type, init_random_seq, &
   random_gaussian


implicit none

type(obs_sequence_type) :: seq
type(obs_def_type) :: obs_def
type(set_def_list_type) :: set_def_list
type(obs_set_def_type) :: obs_set_def
type(obs_set_type) :: obs_set
type(time_type) :: time, time2
type(random_seq_type) :: random_seq

integer :: i, j, obs_set_def_index, unit, unit_out, num_obs_in_set, state_unit

! Need to set up namelists for controlling all of this mess, too!
integer :: model_size, num_obs_sets

type(assim_model_type) :: x
real(r8), allocatable :: obs_err_cov(:), obs(:), true_obs(:)
character(len=129) :: copy_meta_data(2), file_name

! Read in an observation sequence, only definitions part will be used (no data used)
write(*, *) 'input file name for obs sequence definition '
read(*, *) file_name
unit = 10
open(file = file_name, unit = 10)
! Just read in the definition part of the obs sequence
seq = read_obs_sequence_def(unit)

! Initialize the model now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

! Get the initial condition
call init_assim_model(x)
call get_initial_condition(x)

! Set up output of truth for state
state_unit = init_diag_output('true_state', 'true state from control', 1, (/'true state'/))

! Advance for a long time (5 days) to get things started?
!time = set_time(0, 5)
!call advance_state(x, time)

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
do i = 1, num_obs_sets
   time = get_obs_sequence_time(seq, i)
   write(*, *) 'time of obs set ', i
   call print_time(time)

! For now, set model initial time to time of first obs_set, want to add much more
! namelist control for this eventually
   if(i == 1) call set_model_time(x, time)

! Figure out time to advance to
   time2 = get_closest_state_time_to(x, time)
! Advance the state to this time; zero length advance is problem for B-grid so avoid
   if(time2 /= get_model_time(x)) call advance_state(x, time2)

! Output the true state
   call output_diagnostics(state_unit, x, 1)

! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)
   write(*, *) 'num_obs_in_set is ', num_obs_in_set

! Allocate storage for the observations and covariances
   allocate(obs_err_cov(num_obs_in_set), &
      true_obs(num_obs_in_set), obs(num_obs_in_set))

! Compute the observations from the state
   call get_expected_obs(seq, i, get_model_state_vector(x), true_obs)
!   write(*, *) 'exact obs ', obs

! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

! Generate the synthetic observations by adding in error samples
   do j = 1, num_obs_in_set
      obs(j) = random_gaussian(random_seq, true_obs(j), sqrt(obs_err_cov(j)))
   end do
!   write(*, *) 'obs with error added are ', obs

! Insert the observations into the sequence first copy
   call set_obs_values(seq, i, true_obs, 2)
   call set_obs_values(seq, i, obs, 1)

! Deallocate the obs size storage
   deallocate(obs_err_cov, true_obs, obs)

end do

10 continue
! Write out the sequence
write(*, *) 'What is file name for output obs sequence?'
read(*, *) file_name
!unit_out = open_file(file_name, action = 'write')
unit_out = 11
open(file = file_name, unit = 11)
call write_obs_sequence(unit_out, seq)

end program perfect_model_obs
