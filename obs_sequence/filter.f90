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
   get_obs_values, get_num_close_states, get_close_states
use time_manager_mod, only : time_type, set_time, print_time
use utilities_mod, only : open_file
use assim_model_mod, only : assim_model_type, initialize_assim_model, get_model_size, &
   get_initial_condition, get_model_state_vector, set_model_state_vector, &
   get_closest_state_time_to, advance_state, set_model_time, &
   get_model_time
use random_seq_mod, only : random_seq_type, init_random_seq, &
   random_gaussian
use assim_tools_mod, only : obs_increment, update_from_obs_inc


implicit none

type(obs_sequence_type) :: seq
type(time_type) :: time, time2
type(random_seq_type) :: random_seq

real(r8), parameter :: cutoff_radius = 0.2
character(len = 129) file_name

integer :: i, j, k, ind, unit, num_obs_in_set

! Need to set up namelists for controlling all of this mess, too!
integer, parameter :: ens_size = 20
integer :: model_size, num_obs_sets

type(assim_model_type) :: x, ens(ens_size)
real(r8) :: obs_inc(ens_size), ens_inc(ens_size), ens_obs(ens_size)
real(r8), allocatable :: ens_temp(:), obs_err_cov(:)
real(r8), allocatable :: obs(:), ens_mean(:)
! Real concern about reduncant storage needs for ensemble, must fix
real(r8), allocatable :: ens_state(:, :)
integer, allocatable :: num_close_states(:), close(:, :)
real(r8) :: rms, sum_rms = 0.0_r8

! Input the obs_sequence
write(*, *) 'input name of obs sequence file'
read(*, *) file_name
unit = open_file(file_name)
seq = read_obs_sequence(unit)
close(unit)


! Initialize the model now that obs_sequence is all set up
call initialize_assim_model()
model_size = get_model_size()

! Get the initial condition
x = get_initial_condition()

! Advance for a long time (5 days) to get things started?
time = set_time(0, 5)
x = advance_state(x, time)

! Reset the control model time to 0
time = set_time(0, 0)
call set_model_time(x, time)

! Initialize a repeatable random sequence for perturbations
call init_random_seq(random_seq)
! Perturb for ensembles; 
allocate(ens_temp(model_size), ens_mean(model_size), &
   ens_state(model_size, ens_size))
do i = 1, ens_size
   ens(i) = x
   ens_temp = get_model_state_vector(ens(i))
   do j = 1, model_size
      ens_temp(j) = ens_temp(j) + random_gaussian(random_seq, 0.0_r8, 1.0_r8)
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
   time2 = get_closest_state_time_to(x, time)
! Advance the state and the ensembles to this time
   x = advance_state(x, time2)
   do j = 1, ens_size
      ens(j) =  advance_state(ens(j), time2)
   end do

! Yucky, need to find a less costly way to do this
   do j = 1, ens_size
      ens_state(:, j) = get_model_state_vector(ens(j))
   end do

! Compute a cheating RMS for now; WARNING: Control x is not
! neccesarily truth!!!
   ens_mean = sum(ens_state, dim=2) / ens_size
   rms = real(sqrt(sum((get_model_state_vector(x) - ens_mean) &
      * (get_model_state_vector(x) - ens_mean))))
   sum_rms = sum_rms + rms
   write(*, *) 'scaled distance to truth is ', rms

! Do a covariance inflation for now?
   do j = 1, ens_size
      ens_state(:, j) = ens_mean + (ens_state(:, j) - ens_mean) * 1.06
   end do

! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)
!   write(*, *) 'num_obs_in_set is ', num_obs_in_set

! Allocate storage for the ensemble priors for this number of observations
   allocate(obs_err_cov(num_obs_in_set), obs(num_obs_in_set), &
      num_close_states(num_obs_in_set))

! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

! Get the observations
   call get_obs_values(seq, i, obs)

! Get the number of close states for each of these obs
   call get_num_close_states(seq, i, cutoff_radius, num_close_states)

! Get the list of close states for this set (look at storage issues)
   allocate(close(num_obs_in_set, maxval(num_close_states)))
   call get_close_states(seq, i, cutoff_radius, num_close_states, close)

! Loop through each observation in the set
   do j = 1, num_obs_in_set

! Compute the ensemble prior for this ob
      do k = 1, ens_size
         call get_expected_obs(seq, i, ens_state(:, k), &
            ens_obs(k:k), j)
      end do



      call obs_increment(ens_obs, ens_size, obs(j), obs_err_cov(j), &
         obs_inc)

! Now loop through each close state variable for this observation
      do k = 1, num_close_states(j)
         ind = close(j, k)
         call update_from_obs_inc(ens_obs, obs_inc, &
            ens_state(ind, :), ens_size, ens_inc, 1.0_r8)
         ens_state(ind, :) = ens_state(ind, :) + ens_inc
      end do

   end do

! Put the ensemble storage back into the ens
   do j = 1, ens_size
      call set_model_state_vector(ens(j), ens_state(:, j))
   end do


! Deallocate the ens_obs storage
   deallocate(close)
   deallocate(obs_err_cov, obs, num_close_states)

end do

write(*, *) 'mean rms is ', sum_rms / num_obs_sets


end program filter
