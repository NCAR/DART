! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim302

! UPDATE from 22 Sept. 2003 for obs space factor correction.
! WARNING:WARNING:WARNING, when mean difference is much less than
! the covariances, the factor hear are NOT independent of the mean
! difference for the mean. BUT, This is just a secondary sampling issue.

! Work done during last week of January 2002 (28 Jan init) to investigate
! impacts of sampling error from small ensembles on update for a single 
! variable that is exactly observed; applicable to observation variable 
! priors. Small sample estimates of variance have approximately normal
! (or is it exactly normal) error distributions. BUT, when one computes
! the updated variance there is a bias. Of course, one must also account
! for the increased uncertainty in the estimate of the mean resulting from
! errors in the computation of the variance.

! This piece looks at what large sample MC gives for correct updated 
! distribution statistics in preparation for correcting EAKF for the small
! sample problems.


use types_mod, only : r8
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

type (random_seq_type) :: r
real(r8) :: growth, prior_var, prior_mean, obs_var, obs, new_var, new_mean
integer :: n_times, i

! Initialize repeatable random sequence
call init_random_seq(r) 

write(*, *) 'Input growth factor'
read(*, *) growth

write(*, *) 'Input number of times'
read(*, *) n_times

prior_var = 1.0
prior_mean = 0.0
obs_var = 1.0
do i = 1, n_times
   write(*, *) prior_mean, prior_var
! Take an observation
   obs = random_gaussian(r, dble(0.0), dble(sqrt(obs_var)))
   call update(prior_mean, prior_var, obs, obs_var, new_mean, new_var)
! Now do equivalent of time advance
   prior_mean = growth * new_mean
   prior_var = new_var * growth**2
end do

contains 

!---------------------------------------------------
subroutine update(prior_mean, prior_var, obs, obs_var, new_mean, new_var)

implicit none

real(r8), intent(in) :: prior_mean, prior_var, obs, obs_var
real(r8), intent(out) :: new_mean, new_var

real(r8) :: error, diff_sd, ratio

! Base computation
new_var = 1.0 / (1.0 / prior_var + 1.0 / obs_var)
new_mean = new_var * (prior_mean / prior_var + obs / obs_var)

error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_var)
ratio = abs(error / diff_sd)

if(ratio < 1.0) then
   new_var = new_var * 0.9
else
   new_var = new_var * 1.1
endif




end subroutine update

end program sys_sim302

