! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program sys_sim502

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!



! Looks at how to deal with correlation errors given null assumption about 
! a uniform distribution for the underlying true correlation distirubtions.
! Generates a huge sample of [-1, 1) uniform correlations, then samples these,
! and keeps track of the sample correlation / true correlation pairs. Then
! it bins these together over small sample correlation ranges and looks at
! the value of alpha that would minimize errors. NOTE: this gives much less
! info than group filter because group filter draws from a much more 
! meaningful prior (i.e. from what the model and assim have generated).


use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type (random_seq_type) :: r
double precision :: correl, cov(2, 2), zero_2(2) = 0.0, correl_mean, correl_var
double precision :: q, alpha, alpha2, mean_abs_correl
double precision, allocatable :: pairs(:, :), sample_correl(:), real_correl(:)
double precision :: bottom_cor, top_cor, sum_sample2, sum_sample_real, sum_real
integer :: ens_size, n_samples, i, j, k, sample_size

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
allocate(pairs(2, ens_size), sample_correl(n_samples), real_correl(n_samples))

! Loop through a set of sample draws from uniform (-1, 1)
do k = 1, n_samples
   111 real_correl(k) = -1.0 + 2.0 * random_uniform(r)
   correl = real_correl(k)
   ! Generate the covariance matrix for this correlation
   cov(1, 1) = 1.0;    cov(2, 2) = 1.0
   cov(1, 2) = correl; cov(2, 1) = correl

   ! Loop to generate an ensemble size sample from this correl
   ! Generate a random sample of size ens_size from something with this correlation
   do i = 1, ens_size
      call twod_gaussians(r, zero_2, cov, pairs(:, i))
   end do

   ! Compute the sample correlation
   call comp_correl(pairs, ens_size, sample_correl(k)) 

   if(sample_correl(k) > 0.1 .or. sample_correl(k) < -0.1) goto 111
   
   !!!write(*, *) k, real_correl(k), sample_correl(k)

end do

! Now, gather all the sample_correls from a given range
do j = 1, 200
   bottom_cor = -1.0 + (j - 1) * 0.01 
   top_cor = bottom_cor + 0.01
   write(*, *) '------', j, bottom_cor, top_cor, '--------'
   ! Gather up all the SAMPLE correls in this range
   sum_sample2 = 0.0
   sum_sample_real = 0.0
   sum_real = 0.0
   sample_size = 0
   do k = 1, n_samples
      if(sample_correl(k) >= bottom_cor .and. sample_correl(k) < top_cor) then
         !!!write(*, *) sample_correl(k), real_correl(k)
         sample_size = sample_size + 1
         sum_sample2 = sum_sample2 + sample_correl(k)**2
         sum_sample_real = sum_sample_real + sample_correl(k) * real_correl(k)
         sum_real = sum_real + real_correl(k)
      endif
   end do
   write(*, *) 'optimal alpha is ', sum_sample_real / sum_sample2, sample_size
   write(*, *) 'mean real_correl is ', sum_real / sample_size
end do

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



