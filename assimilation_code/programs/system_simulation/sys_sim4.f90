! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim4

! Work done during last week of January 2002 (28 Jan init) to investigate
! impacts of sampling error from small ensembles on update for a single 
! variable that is exactly observed; applicable to observation variable 
! priors. Small sample estimates of variance have approximately normal
! (or is it exactly normal) error distributions. BUT, when one computes
! the updated variance there is a bias. Of course, one must also account
! for the increased uncertainty in the estimate of the mean resulting from
! errors in the computation of the variance.

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

type (random_seq_type) :: r
double precision, allocatable :: rnum(:), var_hist(:)
double precision :: sample_mean, sample_sigma_y_p
double precision :: c(2, 2), sum_err_var, sigma_y_p, sigma_y_o, sigma_x_p
double precision :: y_truth, y_o, sample_correl, reg_coef, sample_reg_coef
double precision :: sigma_y_u, y_u, sigma_x_u, x_u, sample_x_u, x_error
double precision :: x_error_var, total_err_var, mean(2), correl, est_err_var
double precision :: temp_sum, pe_2_sum, xx, s_xx, s_xx_2, xx_mean, xx_var
double precision, allocatable :: x_mean(:), x_var(:)
double precision :: sum_sample_mean, sum_sample_sigma, sum_y_u, sum_sigma_y_u, sum_y_u_2
double precision :: sum_err, sum_err_2, inf_sigma_y_u, var_mean, var_var
integer :: n, i, j, n_samples, index

! Initialize repeatable random sequence
call init_random_seq(r) 

! Initialize error accumulation for averaging
sum_err_var = 0.0
mean = 0.0; temp_sum = 0.0; pe_2_sum = 0.0; 

! For now have a set of free parameters that may be too large
write(*, *) 'Input prior variance of observation variable'
read(*, *) sigma_y_p
write(*, *) 'Input variance of observing instrument'
read(*, *) sigma_y_o
write(*, *) 'Input ensemble size'
read(*, *) n

! Allocate storage for computing sample variance 
allocate(rnum(n))

write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples
allocate(var_hist(n_samples))

sum_sample_mean = 0.0
sum_sample_sigma = 0.0
sum_y_u = 0.0; sum_y_u_2 = 0.0
sum_sigma_y_u = 0.0
sum_err = 0.0; sum_err_2 = 0.0

! Loop through the number of samples
do i = 1, n_samples

! Produce a sample of an observation yo by getting truth from prior
! and adding on sample from observatoinal error distribution
   y_truth = random_gaussian(r, dble(0.0), sqrt(sigma_y_p))
   y_o = y_truth + random_gaussian(r, dble(0.0), sqrt(sigma_y_o))

! Compute the sample variance for this ensemble size
   do j = 1, n
      rnum(j) =  random_gaussian(r, dble(0.0), sqrt(sigma_y_p))
   end do

! Compute the sample variance
   call sample_mean_var(rnum, n, sample_mean, sample_sigma_y_p)
   sum_sample_mean = sum_sample_mean + sample_mean
   sum_sample_sigma = sum_sample_sigma + sample_sigma_y_p
   var_hist(i) = sample_sigma_y_p

!   write(*, *) 'sample mean and variance ', real(sample_mean), real(sample_sigma_y_p)

! Compute the updated covariance and mean for y given this observation
   sigma_y_u = (1.0 / (1.0 / sample_sigma_y_p + 1.0 / sigma_y_o))
   y_u = sigma_y_u * (0.0 + y_o / sigma_y_o)

   sum_y_u = sum_y_u + y_u
   sum_y_u_2 = sum_y_u_2 + y_u**2
   sum_sigma_y_u = sum_sigma_y_u + sigma_y_u
   sum_err = sum_err + y_u - y_truth
   sum_err_2 = sum_err_2 + (y_u - y_truth) ** 2
   

end do

write(*, *)'mean sample mean and var ', real(sum_sample_mean / n_samples), real(sum_sample_sigma / n_samples)
write(*, *) 'mean update mean and var ', real(sum_y_u / n_samples), real(sum_sigma_y_u / n_samples)
write(*, *) 'mean err ', real(sum_err / n_samples), real(sum_err_2 / n_samples)

! Correct updated error variance is computed here
inf_sigma_y_u = 1.0 / (1.0 / sigma_y_p + 1.0 / sigma_y_o)
write(*, *) 'infinite sample updated error variance ', real(inf_sigma_y_u)
write(*, *) '---------------------------------'

! Output some ratios of the spread measure
write(*, *) 'ratio of predicted var to correct ', real((sum_sigma_y_u / n_samples)  / inf_sigma_y_u), &
   real(inf_sigma_y_u / (sum_sigma_y_u / n_samples))
write(*, *) 'ratio of actual var to correct ', real((sum_err_2 / n_samples) / inf_sigma_y_u)

! Compute variance of original sampled variance
call sample_mean_var(var_hist, n_samples, var_mean, var_var)
write(*, *) '-----------------------------------'
write(*, *) 'original sample mean and var of var ', real(var_mean), real(var_var)
write(*, *) 'predicted var of var is ', 2.0 * sigma_y_p**4 / (n - 1.0)


if(1 == 1) stop

! Output the expected value of the mean x
write(*, *) 'Expected updated mean is ', sum(x_mean) / n_samples

! Now compute the expected value of the updated variance by doing sampling
! The distribution is a sum of Gaussians with means and variance given
! by x_mean and x_var; each is equally likely. So, do a sample, randomly
! select a 'kernel', then randomly select a value from it.
s_xx = 0.0
s_xx_2 = 0.0
do i = 1, 10 * n_samples
   index = random_uniform(r) * n_samples + 1
   if(index < 1 .or. index > n_samples) then
      write(*, *) 'bad index ', i, index
      stop
   end if
   if(i / 500000 * 500000 == i) then
      write(*, *) 'mean, var ', i, index, real(x_mean(index)), real(x_var(index))
   endif
   xx = random_gaussian(r, x_mean(index), sqrt(x_var(index)))
   s_xx = s_xx + xx
   s_xx_2 = s_xx_2 + xx**2
end do

xx_mean = s_xx / (10.0 * n_samples)
write(*, *) 'Sampled mean is ', xx_mean
write(*, *) 'Expected mean is ', x_u
xx_var = (s_xx_2 - (10.0 * n_samples) * xx_mean**2) / &
   (10.0 * n_samples - 1.0)
write(*, *) 'Sampled variance is ', xx_var
write(*, *) 'Expected raw var is ', sigma_x_u
write(*, *) 'Mean variance is    ', s_xx_2 / (10.0 * n_samples)



if(1 == 1) stop

! Output the mean total_err_var
write(*, *) '----------'
write(*, *) 'Large N Mean total_err_var is ', sigma_x_u
write(*, *) 'Expected sample x variance is ', sum_err_var / n_samples
write(*, *) '----------'

! Using sample for product of correlation and y_0_2
write(*, *) 'Expected error variance is                  ', &
   (temp_sum / n_samples )*(sigma_x_p / sigma_y_p) * (sigma_y_u / sigma_y_o)**2

! Estimate using sample only for correlation
write(*, *) 'Expected error variance, only sample correl ', &
   sigma_x_p * sigma_y_p * (pe_2_sum / n_samples) / (sigma_y_p + sigma_y_o)

est_err_var = ((1.0 - correl**2)**2 / (n - 1.)) * sigma_x_p * sigma_y_p / &
   (sigma_y_p + sigma_y_o)
write(*, *) 'Estimated correl expected variance is       ', est_err_var
write(*, *) '----------'

write(*, *) 'pe_2_mean                   ', pe_2_sum / n_samples
write(*, *) 'expected value of pe_2_mean ', (1.0 - correl**2)**2 / (n - 1.)

write(*, *) '----------'
write(*, *) 'net improvement ', (sigma_x_p - sum_err_var / n_samples) / &
   (sigma_x_p - sigma_x_u)
write(*, *) 'Computed net improvement ', 1. - (pe_2_sum / n_samples) / correl**2

end program sys_sim4

!-----------------------------------------------------
 
subroutine comp_correl(ens, n, correl)
 
implicit none
 
integer, intent(in) :: n
double precision, intent(in) :: ens(2, n)
double precision, intent(out) :: correl
double precision :: sum_x, sum_y, sum_xy, sum_x2, sum_y2


sum_x = sum(ens(2, :))
sum_y = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))
 
! Computation of correlation 
sum_y2 = sum(ens(1, :) * ens(1, :))

correl = (n * sum_xy - sum_x * sum_y) / &
   sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))
 
end subroutine comp_correl

!----------------------------------------------------------------

subroutine sample_mean_var(x, n, mean, var)

implicit none

integer, intent(in) :: n
double precision, intent(in) :: x(n)
double precision, intent(out) :: mean, var

double precision :: sx, s_x2

sx = sum(x)
s_x2 = sum(x * x)
mean = sx / n
var = (s_x2 - sx**2 / n) / (n - 1)

end subroutine sample_mean_var

