! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim501

! See work from 13 October, 2003. Test to evaluate how to look_for_bias in 
! prior ensemble distribution and observation in observation space.

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

type (random_seq_type) :: r
double precision :: correl, cov(2, 2), zero_2(2) = 0.0, correl_mean, correl_var
double precision :: q, alpha, alpha2, mean_abs_correl
double precision, allocatable :: pairs(:, :), sample_correl(:)
integer :: ens_size, n_samples, i, j, k

! Initialize repeatable random sequence
call init_random_seq(r) 

! For now have a set of free parameters that may be too large
write(*, *) 'Input ensemble size'
read(*, *) ens_size 
write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples

! Allocate storage for sample from this correlation
allocate(pairs(2, ens_size), sample_correl(n_samples))

! Loop through 101 correlaton values
do k = 1, 101
   correl = (k - 1.0) * 0.01
   ! Generate the covariance matrix for this correlation
   cov(1, 1) = 1.0;    cov(2, 2) = 1.0
   cov(1, 2) = correl; cov(2, 1) = correl

   ! Loop for each sample size
   do j = 1, n_samples
      ! Loop to generate an ensemble size sample from this correl
      ! Generate a random sample of size ens_size from something with this correlation
      do i = 1, ens_size
         call twod_gaussians(r, zero_2, cov, pairs(:, i))
      end do

      ! Compute the sample correlation
      call comp_correl(pairs, ens_size, sample_correl(j))
      !!!   write(*, *) 'sample correl is ', sample_correl(j)
   end do

   ! Compute a mean and standard deviation of the sample correlations
   call sample_mean_var(sample_correl, n_samples, correl_mean, correl_var)
   q = sqrt(correl_var) / abs(correl_mean)
   write(*, *) 'correl mean and sd', correl_mean, sqrt(correl_var)
   write(*, *) 'ratio of sd to mean is ', q

   ! Also compute the mean of the absolute value of correl
   mean_abs_correl = sum(abs(sample_correl)) / n_samples
   write(*, *) 'correl, mean_abs_correl ', correl, mean_abs_correl

   ! Value of alpha (regression confidence factor) is computed as limit of
   ! group size m -> infinity from paper: alpha = 1 / (Q^2 + 1)
   alpha = 1.0 / (q**2 + 1)
   if(alpha < 0) alpha = 0.0
   ! Maybe I should really be doing this in limit of m = 1, single sample?
   alpha2 = 1.0 - q**2
   if(alpha2 < 0.0) alpha2 = 0.0
   write(*, *)  correl, alpha, alpha2
end do

end program sys_sim501

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

