! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim502

! Looks at how to deal with correlation errors given null assumption about 
! a uniform distribution for the underlying true correlation distirubtions.
! Generates a huge sample of [-1, 1) uniform correlations, then samples these,
! and keeps track of the sample correlation / true correlation pairs. Then
! it bins these together over small sample correlation ranges and looks at
! the value of alpha that would minimize errors. NOTE: this gives much less
! info than group filter because group filter draws from a much more 
! meaningful prior (i.e. from what the model and assim have generated).

use types_mod, only : r8
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

integer, parameter :: num_bins = 50
real(r8) :: sum(num_bins) = 0.0, sum2(num_bins)
integer  :: num(num_bins) = 0.0
type (random_seq_type) :: r
real(r8) :: correl, cov(2, 2), zero_2(2) = 0.0, correl_mean, correl_var, sd
real(r8) :: q, alpha, alpha2, mean_abs_correl, sample_correl, real_correl
real(r8), allocatable :: pairs(:, :)
real(r8) :: bottom_cor, top_cor, sum_sample2, sum_sample_real, sum_real
real(r8) :: q_factor
integer :: ens_size, n_samples, i, j, k, sample_size, bin

! Generate a bazillion cases where the background correlation is assumed to 
! be U(-1, 1)
! For each case sampled from this, create a SAMPLE correlation
! and store the sample correlation and the actual correlation in a huge
! array.

! Initialize repeatable random sequence
call init_random_seq(r) 

! For now have a set of free parameters that may be too large
write(*, *) 'Input ensemble size'
read(*, *) ens_size 
write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples

! Allocate storage for sample from this correlation
allocate(pairs(2, ens_size))

! Loop through a set of sample draws from uniform (-1, 1)
do k = 1, n_samples
   111 real_correl = -1.0 + 2.0 * random_uniform(r)
   ! Generate the covariance matrix for this correlation
   cov(1, 1) = 1.0;    cov(2, 2) = 1.0
   cov(1, 2) = real_correl; cov(2, 1) = real_correl

   ! Loop to generate an ensemble size sample from this correl
   ! Generate a random sample of size ens_size from something with this correlation
   do i = 1, ens_size
      call twod_gaussians(r, zero_2, cov, pairs(:, i))
   end do

   ! Compute the sample correlation
   call comp_correl(pairs, ens_size, sample_correl) 

   !!!   if(sample_correl > 0.1 .or. sample_correl < -0.1) goto 111
   
   !!!write(*, *) k, real_correl, sample_correl

   ! Now adjust the information for the appropriate sample_correl bin
   bin = (sample_correl + 1.0) / 2.0 * (num_bins) + 1

   num(bin) = num(bin) + 1
   sum(bin) = sum(bin) + real_correl
   sum2(bin) = sum2(bin) + real_correl*real_correl

end do

do i = 1, num_bins
   ! Output sample and mean real for this bin
   sample_correl = ((i - 0.5) / num_bins) * 2.0 - 1.0
   real_correl = sum(i) / num(i)
   sd = sqrt((sum2(i) - num(i) * real_correl**2) / (num(i) - 1.0))
   q = sd / abs(real_correl)
   q_factor = (num(i) - q**2) / ((num(i) - 1) * q + num(i))
   write(*, 11) i, num(i), sample_correl, real_correl, &
      real_correl / sample_correl, sd, q_factor
   11 format(i4, 1x, i8, 1x, 5(e10.4, 1x))
end do

if(1 == 1) stop

!!!   write(*, *) 'optimal alpha is ', sum_sample_real / sum_sample2, sample_size
!!!   write(*, *) 'mean real_correl is ', sum_real / sample_size

end program sys_sim502

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

