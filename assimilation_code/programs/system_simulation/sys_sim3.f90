! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> In this variant modified slightly from system simulation original, 
!> try to determine what fraction of the recommended small sample update 
!> should be used to adjust the mean to get the smallest error variance. 
!> See notes from ~ 19 Dec. 2001.

!> See notes from first two weeks of December, 2001. This program begins attempts 
!> to analyze the value of particular observations. Here, we begin by trying to 
!> determine the value of observations with a given correlation to a state variable
!> using an N member ensemble to compute the correlations.

program sys_sim3

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, twod_gaussians

implicit none

type (random_seq_type) :: r
double precision, allocatable :: rnum(:, :)
double precision :: c(2, 2), sum_err_var, sigma_y_p, sigma_y_o, sigma_x_p
double precision :: y_truth, y_o, sample_correl, reg_coef, sample_reg_coef
double precision :: sigma_y_u, y_u, sigma_x_u, x_u, sample_x_u, x_error
double precision :: x_error_var, total_err_var, mean(2), correl, est_err_var
double precision :: temp_sum, pe_2_sum, fract
integer :: n, i, j, k, n_samples



! For now have a set of free parameters that may be too large
write(*, *) 'Input prior variance of observation variable'
read(*, *) sigma_y_p
write(*, *) 'Input variance of observing instrument'
read(*, *) sigma_y_o
write(*, *) 'Input prior variance of state variable'
read(*, *) sigma_x_p
write(*, *) 'Input ensemble size'
read(*, *) n
write(*, *) 'input an expected correlation'
read(*, *) correl 

! Allocate storage for doing the correlation
allocate(rnum(2, n))

write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples


do k = 75 , 95
   fract = k / 100.0

! Initialize repeatable random sequence
call init_random_seq(r) 

! Initialize error accumulation for averaging
sum_err_var = 0.0
mean = 0.0; temp_sum = 0.0; pe_2_sum = 0.0; 


!!!write(*, *) 'input a fraction to move'
!!!read(*, *) fract

! Set up covariance matrix for computing the sample correlation
c(1, 1) = 1.0; c(2, 2) = 1.0; c(1, 2) = correl; c(2, 1) = correl

! Loop through the number of samples
do i = 1, n_samples

! Produce a sample of an observation yo by getting truth from prior
! and adding on sample from observatoinal error distribution
   y_truth = random_gaussian(r, dble(0.0), sqrt(sigma_y_p))
   y_o = y_truth + random_gaussian(r, dble(0.0), sqrt(sigma_y_o))

! Use random sample to generate erroneous sample correlation
! Generate n pairs of numbers from this distribution
   do j = 1, n
      call twod_gaussians(r, mean, c, rnum(:, j))
   end do

! Compute the sample correlation
   call comp_correl(rnum, n, sample_correl)

! Compute temporary sum to evaluate partial analytic solution
   temp_sum = temp_sum + (sample_correl - correl)**2 * y_o**2
   pe_2_sum = pe_2_sum + (sample_correl - correl)**2

! Compute regression coefficient, erroneous regression coefficient
   reg_coef = correl * sqrt(sigma_x_p) / sqrt(sigma_y_p)
   sample_reg_coef = sample_correl * sqrt(sigma_x_p) / sqrt(sigma_y_p)

! Compute the updated covariance and mean for y given this observation
   sigma_y_u = 1.0 / (1.0 / sigma_y_p + 1.0 / sigma_y_o)
   y_u = sigma_y_u * (0.0 + y_o / sigma_y_o)

! Compute updated covariance for x
! See notes from 11 Dec. for sigma_x_u computation
   sigma_x_u = sigma_x_p - correl**2 * sigma_x_p * sigma_y_p &
      / (sigma_y_o + sigma_y_p)

! Compute correct x_u and sample x_u
   x_u = reg_coef * y_u
!!!   sample_x_u = sample_reg_coef * y_u
! Only move a fraction of the way and see what happens
   sample_x_u = sample_reg_coef * y_u * fract

! Compute error in x
   x_error = sample_x_u - x_u
   x_error_var = x_error **2

! Total expected variance for this one is x_error_var plus sigma_x_u
   total_err_var = x_error_var + sigma_x_u

! Accumulate this for averaging
   sum_err_var = sum_err_var + total_err_var

end do

write(*, *) real(fract), 'Expected sample x variance is ', real(sum_err_var / n_samples)

! This is end of k loop search for best value of fract
end do

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

end program sys_sim3

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

