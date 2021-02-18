! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim202

! UPDATE fro 22 September, 2003 looking at sampling ratios in obs space.
!
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

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

type (random_seq_type) :: r
double precision, allocatable :: rnum(:)
double precision :: sigma_y_p, y_p, sigma_y_o, y_o
double precision :: sum_mean, sum_var, var, sigma_y_u
double precision :: y_u, sample_mean, sample_var, exact_mean
double precision :: exact_var, correct_var, mean_ratio, var_ratio
double precision :: sum_fraction, exact_fraction, fraction, fraction_ratio
integer :: n, n_samples, i, j

! Initialize repeatable random sequence
call init_random_seq(r) 

! For now have a set of free parameters that may be too large
write(*, *) 'Input prior variance of observation variable'
read(*, *) sigma_y_p
! Make prior mean 0.0
!write(*, *) 'Input prior value of observation variable'
!read(*, *) y_p
y_p = 0.0
write(*, *) 'Input variance of observing instrument'
read(*, *) sigma_y_o
write(*, *) 'input value of observation '
read(*, *) y_o
write(*, *) 'Input ensemble size'
read(*, *) n

write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples

! What would answer be without sampling error?
exact_var = (1.0 / (1.0 / sigma_y_p + 1.0 / sigma_y_o))
exact_mean = exact_var * (y_p / sigma_y_p + y_o / sigma_y_o)
exact_fraction = (exact_mean - y_p) / (y_o - y_p)

! Allocate storage for computing sample variance
allocate(rnum(n))

! Initialize summations
sum_mean = 0.0
sum_var = 0.0
sum_fraction = 0.0

! Loop through the number of samples
do i = 1, n_samples

! Compute the sample variance for this ensemble size
   do j = 1, n
      rnum(j) =  random_gaussian(r, y_p, sqrt(sigma_y_p))
   end do

! Compute the sample variance
   call sample_mean_var(rnum, n, sample_mean, var)

! Given this var and mean, compute updated var and mean
   sigma_y_u = (1.0 / (1.0 / var + 1.0 / sigma_y_o))
   y_u = sigma_y_u * (sample_mean / var + y_o / sigma_y_o)
   fraction = (y_u - y_p) / (y_o - y_p)

! Need to keep sum of the means, sum of the mean squared, sum of the var
   sum_mean = sum_mean + y_u
   sum_var = sum_var + sigma_y_u
   sum_fraction = sum_fraction + fraction

end do

!write(*, *) 'Exact mean and var ', exact_mean, exact_var
!write(*, *) 'Mean  mean and var ', sum_mean / n_samples, sum_var / n_samples

! Look at fraction that mean moves as another possible statistic
fraction_ratio = (sum_fraction / n_samples) / exact_fraction
!write(*, *) 'exact and mean fraction', exact_fraction, sum_fraction / n_samples
write(*, *) 'fraction ratio ', fraction_ratio

mean_ratio = (sum_mean / n_samples) / exact_mean
var_ratio = (sum_var / n_samples) / exact_var
write(*, *) 'mean var ratio ', mean_ratio, var_ratio

! Verify that uncertainty is not increased
!correct_var = exact_var / var_ratio
!write(*, *) 'prior, exact, corr ', real(sigma_y_p), real(exact_var), &
!   real(correct_var)

end program sys_sim202

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

