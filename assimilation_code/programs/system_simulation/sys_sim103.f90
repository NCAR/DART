! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim103

! Now, continue to look at the square error for now, but consider the 
! additional error that comes from the fact that the value of the 
! mean and sd of the regression coefficient themselves come from a 
! sample of size num_groups and so have error. Do the Monte Carlo
! by resampling from a Gaussian with specified base statistics. This
! isn't really right (need to do the inverse sampling problem) but
! is probably a pretty good first cut?


! NOTE: THIS ACTUALLY MINIMIZES THE SQUARE ERROR! MAY WORK WELL, BUT NOT
! NECESSARILY RIGHT QUESTION. GET 1/(x**2 + 1) for the alpha where x
! is the ratio of s.d. to sample mean. Use a single enseble member to
! get appropriate results.
! A look at sampling errors for regression coefficients for ensemble filters
! designed to help with localization from multiple ensembles. Notes (if they
! exist) from late August 2003 describe the problem in more detail. Here, we
! ask the question, if a distance to move comes from a distibution with 
! a given ratio of standard deviation to value, when a given size sample
! is used to evaluate this, how should the resulting distance be scaled
! to extract the most useful information.
!

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

type (random_seq_type) :: r
double precision :: sd_ratio, mean_move, mean_move_2, alpha
double precision, allocatable :: rnum(:), var_hist(:)
double precision :: sample_sigma_y_p, base_mean
double precision :: c(2, 2), sum_err_var, sigma_y_p, sigma_y_o, sigma_x_p
double precision :: y_truth, y_o, sample_correl, reg_coef, sample_reg_coef
double precision :: sigma_y_u, y_u, sigma_x_u, x_u, sample_x_u, x_error
double precision :: x_error_var, total_err_var, mean(2), correl, est_err_var
double precision :: temp_sum, pe_2_sum, xx, s_xx, s_xx_2, xx_mean, xx_var
double precision, allocatable :: sample_mean(:), x_var(:)
double precision :: sum_sample_mean, sum_sample_sigma, sum_y_u, sum_sigma_y_u, sum_y_u_2
double precision :: sum_err, sum_err_2, inf_sigma_y_u, var_mean, var_var
integer :: n, i, j, n_samples, index

! Initialize repeatable random sequence
call init_random_seq(r) 

! Free parameters are ratio and ensemble size
! Just assume that base value is 1.0 (move is 1.0); can't do case where expected move is 0
! but know the answer to this case (DON"T MOVE!)
write(*, *) 'Input ratio of ensemble standard deviation to mean'
read(*, *) sd_ratio
write(*, *) 'Input num groups'
read(*, *) n
write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples

! Allocate storage for the sample means
allocate(sample_mean(n_samples), rnum(n))

! Loop through the number of samples
do i = 1, n_samples

! This first part approximates sampling error from the num groups!!!
! Need an ensemble sized sample of the recursion value
   do j = 1, n
      rnum(j) = random_gaussian(r, dble(1.0), sd_ratio)
   end do

! Compute the sample mean and variance
!!!   call sample_mean_var(rnum, n, sample_mean(i), sample_sigma_y_p)
   call sample_mean_var(rnum, n, base_mean, sample_sigma_y_p)

write(*, *) i, 'ratio ', sqrt(sample_sigma_y_p) / abs(base_mean)

! Now, need to get a single sample from the num_groups sample distribution
   sample_mean(i) = random_gaussian(r, base_mean, sqrt(sample_sigma_y_p))
end do

! Need the sum of the squared sample_mean and sum of the sample mean
mean_move = sum(sample_mean)
mean_move_2 = sum(sample_mean**2)

! Expected value of alpha to minimize squared difference between sample and actual value
alpha = 1.0 * mean_move / mean_move_2

write(*, *) 'alpha is ', alpha

end program sys_sim103

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

