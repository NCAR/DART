program filter

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$


! WARNING: The question of where to put in covariance inflation is
! still central. The old seq_obs program inflates the whole ensemble
! up front whenever an observation is present. This is clearly incorrect
! in the general case, since a single low quality local observation should
! NOT impact the whole state distribution. However, in some homogeneous
! cases, it may still produce better results. The implementation in the
! filter at the moment only applies covariance inflate when computing the
! update increments (and the regression ???).

! Program to build a simple obs_sequence file for use in testing filters
! for spatial domains with one periodic dimension.

use types_mod
use obs_sequence_mod, only : obs_sequence_type, write_obs_sequence, &
   read_obs_sequence, get_num_obs_sets, get_obs_sequence_time, &
   get_num_obs_in_set, get_expected_obs, get_diag_obs_err_cov, &
   get_obs_values, &
   obs_sequence_def_copy, inc_num_obs_copies, set_obs_values, &
   set_single_obs_value, get_obs_def_index
use time_manager_mod, only : time_type, set_time, print_time
use utilities_mod, only :  get_unit
use assim_model_mod, only : assim_model_type, static_init_assim_model, &
   get_model_size, get_initial_condition, get_closest_state_time_to, &
   advance_state, set_model_time, get_model_time, init_diag_output, &
   output_diagnostics, init_assim_model, get_state_vector_ptr
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use assim_tools_mod, only : obs_increment, update_from_obs_inc
use cov_cutoff_mod, only : comp_cov_factor
use close_state_cache_mod, only : close_state_cache_type, cache_init, &
   get_close_cache


implicit none

! Define a type for doing direct access to ensemble state vectors
type model_state_ptr_type
   real(r8), pointer :: state(:)
end type model_state_ptr_type


type(obs_sequence_type) :: seq, prior_seq, posterior_seq
type(time_type) :: time, time2
type(random_seq_type) :: random_seq

real(r8), parameter :: cutoff = 1.0, radius = 2.0 * cutoff
real(r8), parameter :: cov_inflate = 1.0_r8
integer, parameter :: cache_size = 10
character(len = 129) file_name

integer :: i, j, k, ind, unit, prior_obs_unit, posterior_obs_unit
integer :: prior_state_unit, posterior_state_unit, num_obs_in_set

! Need to set up namelists for controlling all of this mess, too!
integer, parameter :: ens_size = 100
integer :: model_size, num_obs_sets

! Storage for direct access to ensemble state vectors
type(model_state_ptr_type) :: ens_ptr(ens_size), x_ptr

type(assim_model_type) :: x, ens(ens_size)
real(r8) :: obs_inc(ens_size), ens_inc(ens_size), ens_obs(ens_size)
real(r8) ::  swath(ens_size), ens_mean
real(r8), allocatable :: obs_err_cov(:), obs(:)
real(r8) :: rms, sum_rms = 0.0_r8, cov_factor
character(len = 129) :: ens_copy_meta_data(ens_size)

! Storage for caching close states
type(close_state_cache_type) :: cache
integer, pointer :: num_close_ptr(:), close_ptr(:, :)
real(r8), pointer :: dist_ptr(:, :)


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

! Initialize the control and ensemble states and set up direct pointers
call init_assim_model(x)
x_ptr%state => get_state_vector_ptr(x)
do i = 1, ens_size
   call init_assim_model(ens(i))
   ens_ptr(i)%state => get_state_vector_ptr(ens(i))
end do

! Initialize a cache for close state information
call cache_init(cache, cache_size)

! Get the initial condition
call get_initial_condition(x)

! Set up diagnostic output for model state

prior_state_unit = init_diag_output('prior_diag', &
   'prior ensemble state', ens_size, ens_copy_meta_data)
posterior_state_unit = init_diag_output('posterior_diag', &
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
do i = 1, ens_size
   do j = 1, model_size
      ens_ptr(i)%state(j) = random_gaussian(random_seq, x_ptr%state(j), 1.0_r8)
   end do
end do

num_obs_sets = get_num_obs_sets(seq)

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).
AdvanceTime : do i = 1, num_obs_sets

   time = get_obs_sequence_time(seq, i)
   !   write(*, *) 'time of obs set ', i
   call print_time(time)

   time2 = get_closest_state_time_to(ens(1), time)
   ! Advance the ensembles to this time
   do j = 1, ens_size
      call advance_state(ens(j), time2)
<<<<<<< filter.f90
! Output the prior ensemble state
!      call output_diagnostics(prior_state_unit, ens(j), j)
=======
      ! Output the prior ensemble state
      call output_diagnostics(prior_state_unit, ens(j), j)
>>>>>>> 1.11
   end do

   ! Do a covariance inflation for now? 
   ! Inflate the ensemble state estimates
   do k = 1, model_size
      ens_mean = get_ens_mean(ens_ptr, ens_size, k)
      do j = 1, ens_size
         ens_ptr(j)%state(k) = ens_mean + (ens_ptr(j)%state(k) - &
            ens_mean) * sqrt(cov_inflate)
      end do
   end do

   ! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)
   !   write(*, *) 'num_obs_in_set is ', num_obs_in_set

   ! Allocate storage for the ensemble priors for this number of observations
   allocate(obs_err_cov(num_obs_in_set), obs(num_obs_in_set))

   ! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

   ! Get the observations; from copy 1 for now
   call get_obs_values(seq, i, obs, 1)



   ! Try out the cache
   call get_close_cache(cache, seq, i, radius, num_obs_in_set, &
      num_close_ptr, close_ptr, dist_ptr)


   ! Loop through each observation in the set
   Observations : do j = 1, num_obs_in_set

      ! Compute the ensemble prior for this ob
      do k = 1, ens_size
         call get_expected_obs(seq, i, ens_ptr(k)%state, ens_obs(k:k), j)
      end do

      call obs_increment(ens_obs, ens_size, obs(j), obs_err_cov(j), obs_inc)

      ! Output the ensemble prior and posterior to diagnostic files
      do k = 1, ens_size
         call set_single_obs_value(prior_seq, i, j, ens_obs(k), k)
         call set_single_obs_value(posterior_seq, i, j, ens_obs(k) + obs_inc(k), k)
      end do

      ! Now loop through each close state variable for this observation
      do k = 1, num_close_ptr(j)
         ind = close_ptr(j, k)

         ! Compute distance dependent envelope
         cov_factor = comp_cov_factor(dist_ptr(j, k), cutoff)

         ! Get the ensemble elements for this state variable and do regression
         swath = get_ens_swath(ens_ptr, ens_size, ind)

         call update_from_obs_inc(ens_obs, obs_inc, &
                                  swath, ens_size, ens_inc, cov_factor)
         call inc_ens_swath(ens_ptr, ens_size, ind, ens_inc)
      end do

   end do Observations

   ! Put the ensemble storage back into the ens
   do j = 1, ens_size
<<<<<<< filter.f90
! Output the posterior ensemble state
!      call output_diagnostics(posterior_state_unit, ens(j), j)
=======
      ! Output the posterior ensemble state
      call output_diagnostics(posterior_state_unit, ens(j), j)
>>>>>>> 1.11
   end do

   ! Deallocate the ens_obs storage for this obs set
   deallocate(obs_err_cov, obs)

end do AdvanceTime


! Initialize the model state space diagnostic output files

! Output the observation space diagnostic files

prior_obs_unit = get_unit()
open(unit = prior_obs_unit, file = 'prior_obs_diagnostics')
call write_obs_sequence(prior_obs_unit, prior_seq)
close(prior_obs_unit)

posterior_obs_unit = get_unit()
open(unit = posterior_obs_unit, file = 'posterior_obs_diagnostics')
call write_obs_sequence(posterior_obs_unit, posterior_seq)
close(posterior_obs_unit)



!===========================================================

contains


function get_ens_swath(ens_ptr, ens_size, index)
!----------------------------------------------------------
!
! Returns the ensemble for state variable index

implicit none

integer,                    intent(in) :: ens_size, index
type(model_state_ptr_type), intent(in) :: ens_ptr(ens_size)
real(r8)                               :: get_ens_swath(ens_size)

integer :: i

do i = 1, ens_size
   get_ens_swath(i) = ens_ptr(i)%state(index)
end do

end function get_ens_swath

         

subroutine inc_ens_swath(ens_ptr, ens_size, index, ens_inc)
!-----------------------------------------------------------
!
! Increments all ensemble members for a given state variable
! index.

implicit none

integer,                    intent(in)    :: ens_size, index
type(model_state_ptr_type), intent(inout) :: ens_ptr(ens_size)
real(r8),                   intent(in)    :: ens_inc(ens_size)

integer :: i

do i = 1, ens_size
   ens_ptr(i)%state(index) = ens_ptr(i)%state(index) + ens_inc(i)
end do

end subroutine inc_ens_swath



function get_ens_mean(ens_ptr, ens_size, index)
!----------------------------------------------------------
!
! Computes the ensemble mean of the index element of the
! state vector.

implicit none

integer,                    intent(in) :: ens_size, index
type(model_state_ptr_type), intent(in) :: ens_ptr(ens_size)
real(r8)                               :: get_ens_mean

integer :: i

get_ens_mean = ens_ptr(1)%state(index)
do i = 2, ens_size
   get_ens_mean = get_ens_mean + ens_ptr(i)%state(index)
end do

get_ens_mean = get_ens_mean / ens_size

end function get_ens_mean

end program filter





