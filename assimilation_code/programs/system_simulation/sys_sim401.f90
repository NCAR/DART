! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
 
program sys_sim401

implicit none

integer n, n_samples, i, j
double precision :: sd_ratio, alpha

write(*, *) 'input sample size'
read(*, *) n_samples

! Loop through a range of values for alpha until alpha is 0.0
do i = 1, 600
   sd_ratio = (i - 1.) / 100.0
   call sub_sim401(sd_ratio, n_samples, alpha)
   write(*, *) 'sd_ratio, alpha ', sd_ratio, alpha
!!!   if(alpha <= 0.0) stop
end do

end program sys_sim401



subroutine sub_sim401(sd_ratio, n_samples, alpha)

! sys_sim401 adds in an ad hoc measure of additional uncertainty that comes
! from assuming that the sample ratio of s.d. to mean for the regression 
! value is exact by sampling from this distribution with the group size.
! REVISION OF sys_sim102 TO TRY TO DO ABSOLUTE ERROR, NOT SQUARED ERROR
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

double precision, intent(in) :: sd_ratio
integer, intent(in) :: n_samples
double precision, intent(out) :: alpha

type (random_seq_type) :: r
integer, parameter :: num_tries = 1001
double precision :: mean_abs_error(num_tries), base_mean, ratio
double precision :: mean_move, mean_move_2, sum_abs_error
double precision, allocatable :: rnum(:), var_hist(:)
double precision :: sample_sigma_y_p
double precision :: c(2, 2), sum_err_var, sigma_y_p, sigma_y_o, sigma_x_p
double precision :: y_truth, y_o, sample_correl, reg_coef, sample_reg_coef
double precision :: sigma_y_u, y_u, sigma_x_u, x_u, sample_x_u, x_error
double precision :: x_error_var, total_err_var, mean(2), correl, est_err_var
double precision :: temp_sum, pe_2_sum, xx, s_xx, s_xx_2, xx_mean, xx_var
double precision, allocatable :: sample1(:), sample2(:)
double precision :: sum_sample_mean, sum_sample_sigma, sum_y_u, sum_sigma_y_u, sum_y_u_2
double precision :: sum_err, sum_err_2, inf_sigma_y_u, var_mean, var_var
double precision :: factor
integer :: i, j, k, index, location(1)

! Initialize repeatable random sequence
call init_random_seq(r) 

! Free parameters are ratio and ensemble size
! Just assume that base value is 1.0 (move is 1.0); can't do case where expected move is 0
! but know the answer to this case (DON"T MOVE!)

! Allocate storage for the sample means
allocate(sample1(n_samples), sample2(n_samples))

! Loop through the number of samples
do i = 1, n_samples

! Now, draw a single value from this distribution
! This is the revised clean version from 9 Nov., 2003
   sample1(i) = random_gaussian(r, 1.0_r8, sd_ratio)
   sample2(i) = random_gaussian(r, 1.0_r8, sd_ratio)
end do

! For now, don't know how to really compute this derivative for minimization
! Just try a host of values
do k = 1, num_tries 
   alpha = (k - 1.) / (num_tries - 1.)
! Compute the sum of the absolute errors for this value of alpha
   sum_abs_error = 0.0
   do i = 1, n_samples
      sum_abs_error = sum_abs_error + abs(alpha * sample1(i) - sample2(i)) 
   end do
   mean_abs_error(k) = sum_abs_error / n_samples
end do

location = minloc(mean_abs_error)
alpha = (location(1) - 1.) / (num_tries - 1.)
!write(*, *) 'alpha is ', alpha, 'min is ', minval(mean_abs_error)

! Now do a linear damping to zero at 1000
if(sd_ratio > 3.0) then
   factor = (6.0 - sd_ratio) / 3.0
   alpha = alpha * factor
endif

end subroutine sub_sim401

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

