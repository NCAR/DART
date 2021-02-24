! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim5

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
double precision :: sigma_y_p, y_p, sigma_y_o, y_o, var_of_var
double precision :: sum_mean, sum_mean_2, sum_var, var, sigma_y_u
double precision :: y_u, sample_mean, sample_var, inf_var, inf_u
double precision :: temp_mean, pre_var_sum, mean_var
integer :: n, n_samples, i, j

! Initialize repeatable random sequence
call init_random_seq(r) 

! For now have a set of free parameters that may be too large
write(*, *) 'Input prior variance of observation variable'
read(*, *) sigma_y_p
write(*, *) 'Input prior value of observation variable'
read(*, *) y_p
write(*, *) 'Input variance of observing instrument'
read(*, *) sigma_y_o
write(*, *) 'input value of observation '
read(*, *) y_o
write(*, *) 'Input ensemble size'
read(*, *) n

write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples

! Allocate storage for computing sample variance
allocate(rnum(n))

! Compute the variance of the variance distribution
var_of_var = 2.0 * sigma_y_p**4 / (n - 1.0)

write(*, *) 'var and var of var', real(sigma_y_p), real(var_of_var)

! Initialize summations
sum_mean = 0.0
sum_mean_2 = 0.0
sum_var = 0.0
pre_var_sum = 0.0

! Loop through the number of samples
do i = 1, n_samples

! We are given the sample variance, what is a sample of the underlying
! variance distribution
! Could do this from distribution, but it doesn't apply for small samples
!!!   var = random_gaussian(r, sigma_y_p, sqrt(var_of_var))

! Compute the sample variance for this ensemble size
   do j = 1, n
      rnum(j) =  random_gaussian(r, dble(0.0), sqrt(sigma_y_p))
   end do

! Compute the sample variance
   call sample_mean_var(rnum, n, temp_mean, var)

   pre_var_sum = pre_var_sum + var



! Given this var, compute updated var and mean
   sigma_y_u = (1.0 / (1.0 / var + 1.0 / sigma_y_o))
   y_u = sigma_y_u * (y_p / var + y_o / sigma_y_o)


if(var < 0.0) write(*, *) 'var is , y_u, si', real(var), real(y_u), real(sigma_y_u)


! Need to keep sum of the means, sum of the mean squared, sum of the var
   sum_mean = sum_mean + y_u
   sum_mean_2 = sum_mean_2 + (y_u)**2
   sum_var = sum_var + sigma_y_u

end do

! Output the mean and variance corresponding to finite sample
! See working notes, 29 Jan 2002
sample_mean = sum_mean / n_samples
sample_var = (SUM_var + sum_mean_2) / n_samples - sample_mean**2
write(*, *) 'sample mean and var ', real(sample_mean), real(sample_var)

! What would happen without these assumptions
inf_var = (1.0 / (1.0 / sigma_y_p + 1.0 / sigma_y_o))
inf_u = inf_var * (y_p / sigma_y_p + y_o / sigma_y_o)
write(*, *) 'infini mean and var ', real(inf_u), real(INF_var)

mean_var = sum_var / n_samples
write(*, *) 'mean var i s ', real(mean_var)
write(*, *) 'pre mean var ', real(pre_var_sum / n_samples)

write(*, *) '-----------------'
write(*, *) 'normalized sample var - mean var ', real((sample_var - sum_var / n_samples) / (y_o - y_p)**2)
write(*, *) 'ratio of mean_var to prior var ', real(mean_var /sigma_y_p) 
write(*, *) 'fraction of movement ', real((sample_mean - y_o) / (y_p - y_o))
write(*, *) 'expected fraction of movement ', real((inf_u - y_o) / (y_p - y_o))


end program sys_sim5

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

