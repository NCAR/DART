program filter

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

! Program to build a simple obs_sequence file for use in testing filters
! for spatial domains with one periodic dimension.

use types_mod
use obs_sequence_mod, only : obs_sequence_type, &
   write_obs_sequence, read_obs_sequence, &
   get_num_obs_sets, get_obs_sequence_time, &
   get_num_obs_in_set, get_expected_obs, get_diag_obs_err_cov, &
   get_obs_values, get_num_close_states, get_close_states, &
   obs_sequence_def_copy, inc_num_obs_copies, set_obs_values, set_single_obs_value
use time_manager_mod, only : time_type, set_time, print_time
use utilities_mod, only :  get_unit
use assim_model_mod, only : assim_model_type, static_init_assim_model, get_model_size, &
   get_initial_condition, get_model_state_vector, set_model_state_vector, &
   get_closest_state_time_to, advance_state, set_model_time, &
   get_model_time, init_diag_output, output_diagnostics, init_assim_model, &
   copy_assim_model
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use assim_tools_mod, only : obs_increment, update_from_obs_inc
use cov_cutoff_mod, only : comp_cov_factor


implicit none

type(obs_sequence_type) :: seq, prior_seq, posterior_seq
type(time_type) :: time, time2
type(random_seq_type) :: random_seq

real(r8), parameter :: cutoff = 0.2, radius = 2.0 * cutoff
character(len = 129) file_name

integer :: i, j, k, ind, unit, prior_obs_unit, posterior_obs_unit, num_obs_in_set
integer :: prior_state_unit, posterior_state_unit

! Need to set up namelists for controlling all of this mess, too!
integer, parameter :: ens_size = 20
integer :: model_size, num_obs_sets

type(assim_model_type) :: x, ens(ens_size)
real(r8) :: obs_inc(ens_size), ens_inc(ens_size), ens_obs(ens_size)
real(r8), allocatable :: ens_temp(:), obs_err_cov(:)
real(r8), allocatable :: obs(:), ens_mean(:)
! Real concern about redundant storage needs for ensemble, must fix
real(r8), allocatable :: ens_state(:, :), dist(:, :)
integer, allocatable :: num_close_states(:), close(:, :)
real(r8) :: rms, sum_rms = 0.0_r8, cov_factor
character(len = 129) :: ens_copy_meta_data(ens_size)

! Input the obs_sequence
write(*, *) 'input name of obs sequence file'
read(*, *) file_name
unit = get_unit()
open(unit = unit, file = file_name)
seq = read_obs_sequence(unit)
close(unit)

! Copy just the definitions part of the sequence to the two output obs sequences
call obs_sequence_def_copy(prior_seq, seq)
call obs_sequence_def_copy(posterior_seq, seq)
! Set up the metadata for the output ensemble observations
do i = 1, ens_size
   write(ens_copy_meta_data(i), *) 'ensemble ', i
end do

! For now output all ensemble members for prior and posterior; add space
call inc_num_obs_copies(prior_seq, ens_size, ens_copy_meta_data)
call inc_num_obs_copies(posterior_seq, ens_size, ens_copy_meta_data)

! Initialize the model now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

! Initialize the control and ensemble states
call init_assim_model(x)
do i = 1, ens_size
   call init_assim_model(ens(i))
end do

! Get the initial condition
call init_assim_model(x)
call get_initial_condition(x)

! Set up diagnostic output for model state
prior_state_unit = init_diag_output('prior_state_diagnostics', &
   'prior ensemble state', ens_size, ens_copy_meta_data)
posterior_state_unit = init_diag_output('posterior_state_diagnostics', &
   'posterior ensemble state', ens_size, ens_copy_meta_data)

! Advance for a long time (5 days) to get things started?
! This should all be parameterized and controlled
time = set_time(0, 5)
write(*, *) 'calling advance state for x 5 days'
call advance_state(x, time)
write(*, *) 'back from advance state for x 5 days'

! Reset the control model time to 0
time = set_time(0, 0)
call set_model_time(x, time)

! Initialize a repeatable random sequence for perturbations
call init_random_seq(random_seq)
! Perturb for ensembles; 
allocate(ens_temp(model_size), ens_mean(model_size), &
   ens_state(model_size, ens_size))
do i = 1, ens_size
   call copy_assim_model(ens(i), x)
   ens_temp = get_model_state_vector(ens(i))
   do j = 1, model_size
      ens_temp(j) = random_gaussian(random_seq, ens_temp(j), 1.0_r8)
   end do
   call set_model_state_vector(ens(i), ens_temp)
end do
deallocate(ens_temp)

num_obs_sets = get_num_obs_sets(seq)

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).
do i = 1, num_obs_sets
   time = get_obs_sequence_time(seq, i)
!   write(*, *) 'time of obs set ', i
   call print_time(time)
   time2 = get_closest_state_time_to(ens(1), time)
! Advance the ensembles to this time
   do j = 1, ens_size
      call advance_state(ens(j), time2)
! Output the prior ensemble state
      call output_diagnostics(prior_state_unit, ens(j), j)
   end do

! Yucky, need to find a less costly way to do this, work directly with pointers?
   do j = 1, ens_size
      ens_state(:, j) = get_model_state_vector(ens(j))
   end do

! Do a covariance inflation for now?
   ens_mean = sum(ens_state, dim=2) / ens_size
   do j = 1, ens_size
      ens_state(:, j) = ens_mean + (ens_state(:, j) - ens_mean) * sqrt(1.05)
   end do

! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)
!   write(*, *) 'num_obs_in_set is ', num_obs_in_set

! Allocate storage for the ensemble priors for this number of observations
   allocate(obs_err_cov(num_obs_in_set), obs(num_obs_in_set), &
      num_close_states(num_obs_in_set))

! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

! Get the observations; from copy 1 for now
   call get_obs_values(seq, i, obs, 1)

! THIS IS HORRENDOUSLY SLOW, NEED TO CACHE AT SOME LEVEL
! Get the number of close states for each of these obs
   call get_num_close_states(seq, i, radius, num_close_states)

! Get the list of close states for this set (look at storage issues)
   allocate(close(num_obs_in_set, maxval(num_close_states)), &
      dist(num_obs_in_set, maxval(num_close_states)))
   call get_close_states(seq, i, radius, num_close_states, &
      close, dist)

! Loop through each observation in the set
   do j = 1, num_obs_in_set

! Compute the ensemble prior for this ob
      do k = 1, ens_size
         call get_expected_obs(seq, i, ens_state(:, k), &
            ens_obs(k:k), j)
      end do

      call obs_increment(ens_obs, ens_size, obs(j), obs_err_cov(j), &
         obs_inc)

! Output the ensemble prior and posterior to diagnostic files
      do k = 1, ens_size
         call set_single_obs_value(prior_seq, i, j, ens_obs(k), k)
         call set_single_obs_value(posterior_seq, i, j, ens_obs(k) + obs_inc(k), k)
      end do

! Now loop through each close state variable for this observation
      do k = 1, num_close_states(j)
         ind = close(j, k)

! Compute distance dependent envelope
         cov_factor = comp_cov_factor(dist(j, k), cutoff)

         call update_from_obs_inc(ens_obs, obs_inc, &
            ens_state(ind, :), ens_size, ens_inc, cov_factor)
         ens_state(ind, :) = ens_state(ind, :) + ens_inc
      end do

   end do

! Put the ensemble storage back into the ens
   do j = 1, ens_size
      call set_model_state_vector(ens(j), ens_state(:, j))
! Output the posterior ensemble state
      call output_diagnostics(posterior_state_unit, ens(j), j)
   end do

! Deallocate the ens_obs storage for this obs set
   deallocate(close, dist)
   deallocate(obs_err_cov, obs, num_close_states)

end do


! Initialize the model state space diagnostic output files

! Output the observation space diagnostic files
prior_obs_unit = get_unit
open(unit = get_unit, file = 'prior_obs_diagnostics')
call write_obs_sequence(prior_obs_unit, prior_seq)
close(prior_obs_unit)
posterior_obs_unit = get_unit
open(unit = get_unit, file = 'posterior_obs_diagnostics')
call write_obs_sequence(posterior_obs_unit, posterior_seq)
close(posterior_obs_unit)

end program filter
