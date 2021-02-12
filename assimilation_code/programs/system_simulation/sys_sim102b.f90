! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim102b

! MAKE SURE THAT RATIO IS ALL THAT MATTERS

! REVISION OF sys_sim101 TO TRY TO DO ABSOLUTE ERROR, NOT SQUARED ERROR
! MINIMIZATION.

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
integer, parameter :: num_tries = 1001
double precision :: mean_abs_error(num_tries)
double precision :: sd_ratio, mean_move, mean_move_2, alpha, sum_abs_error
double precision, allocatable :: rnum(:), var_hist(:)
double precision :: sample_sigma_y_p
double precision :: c(2, 2), sum_err_var, sigma_y_p, sigma_y_o, sigma_x_p
double precision :: y_truth, y_o, sample_correl, reg_coef, sample_reg_coef
double precision :: sigma_y_u, y_u, sigma_x_u, x_u, sample_x_u, x_error
double precision :: x_error_var, total_err_var, mean(2), correl, est_err_var
double precision :: temp_sum, pe_2_sum, xx, s_xx, s_xx_2, xx_mean, xx_var
double precision, allocatable :: sample_mean(:), x_var(:)
double precision :: sum_sample_mean, sum_sample_sigma, sum_y_u, sum_sigma_y_u, sum_y_u_2
double precision :: sum_err, sum_err_2, inf_sigma_y_u, var_mean, var_var
integer :: n, i, j, k, n_samples, index, location(1)

! Initialize repeatable random sequence
call init_random_seq(r) 

! Free parameters are ratio and ensemble size
! Just assume that base value is 0.1 (move is 0.1); can't do case where expected move is 0
! but know the answer to this case (DON"T MOVE!)
write(*, *) 'Input ratio of ensemble standard deviation to mean'
read(*, *) sd_ratio
write(*, *) 'Input group size, NORMALLY ONE FOR THIS COMPUTATION'
read(*, *) n
write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples

! Allocate storage for the sample means
allocate(sample_mean(n_samples), rnum(n))

! Loop through the number of samples
do i = 1, n_samples

! Need an ensemble sized sample of the recursion value
   do j = 1, n
      rnum(j) = random_gaussian(r, dble(0.1), 0.1 * sd_ratio)
   end do

! Compute the sample mean and variance
   call sample_mean_var(rnum, n, sample_mean(i), sample_sigma_y_p)
end do

! For now, don't know how to really compute this derivative for minimization
! Just try a host of values
do k = 1, num_tries 
   alpha = (k - 1.) / (num_tries - 1.)
! Compute the sum of the absolute errors for this value of alpha
   sum_abs_error = 0.0
   do i = 1, n_samples
      sum_abs_error = sum_abs_error + abs(alpha * sample_mean(i) - 0.1) 
   end do
   mean_abs_error(k) = sum_abs_error / n_samples
end do

location = minloc(mean_abs_error)
alpha = (location(1) - 1.) / (num_tries - 1.)
write(*, *) 'alpha is ', alpha, 'min is ', minval(mean_abs_error)



! Need the sum of the squared sample_mean and sum of the sample mean
mean_move = sum(sample_mean)
mean_move_2 = sum(sample_mean**2)

! Expected value of alpha to minimize squared difference between sample and actual value
alpha = 1.0 * mean_move / mean_move_2

write(*, *) 'alpha for squared error is ', alpha

end program sys_sim102b

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

