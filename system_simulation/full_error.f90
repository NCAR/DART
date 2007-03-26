! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program full_error

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use types_mod
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &          
   twod_gaussians, random_uniform

implicit none

integer, parameter :: num_times = 1, num_points =  100000000
!!!integer, parameter :: num_times = 1, num_points =  10000000
type (random_seq_type) :: ran_id
real(r8), allocatable :: pairs(:, :)
real(r8) :: zero_2(2) = 0.0, cov(2, 2)
real(r8) :: t_correl, correl_mean, sample_correl, alpha, beta
real(r8) :: s_mean(2), s_var(2), reg_mean, reg_sd, t_sd1, t_sd2, true_correl_mean
real(r8) :: tcorrel_sum(201), reg_sum(201), reg_2_sum(201)
integer i, j, k, ens_size, bin_num, bin_count(201)

call init_random_seq(ran_id)


! Set the ensemble size
ens_size = 40

allocate(pairs(2, ens_size))
write(*, *) 'stats for ensemble size ', ens_size

bin_count = 0
tcorrel_sum = 0.0_r8
reg_sum = 0.0_r8
reg_2_sum = 0.0_r8

! Uniformly distributed sequence of true correlations  
do j = 1, num_points
   ! Assume that true correlation is uniform [-1,1]
   t_correl = -1.0_r8 + 2.0_r8 * ((j - 1) / (num_points - 1.0_r8))
   if(j > num_points) t_correl = 0.0_r8

   do k = 1, num_times

      t_sd1 = 1.0_r8
      t_sd2 = 1.0_r8

      ! Generate the covariance matrix for this correlation
      cov(1, 1) = t_sd1**2;    cov(2, 2) = t_sd2**2
      cov(1, 2) = t_correl * t_sd1 * t_sd2; cov(2, 1) = cov(1, 2)
   
      ! Loop to generate an ensemble size sample from this correl
      ! Generate a random sample of size ens_size from something with this correlation
      do i = 1, ens_size
         call twod_gaussians(ran_id, zero_2, cov, pairs(:, i))
      end do
   
      ! Compute the sample correlation
      call comp_correl(pairs, ens_size, sample_correl)

      ! Bin statistics for each 0.01 sample correlation bin; start with just mean true correlation
      bin_num = (sample_correl + 1.0_r8) / 0.01_r8   + 1
      ! Keep a sum of the sample_correl and squared to compute standard deviation
      tcorrel_sum(bin_num) = tcorrel_sum(bin_num) + t_correl
      !-----------------
      ! Also interested in finding out what the spurious variance reduction factor is
      ! First need to compute sample s.d. for the obs and unobs variable
      do i = 1, 2
         call sample_mean_var(pairs(i, :), ens_size, s_mean(i), s_var(i))
      end do

      !-----------------
      reg_sum(bin_num) = reg_sum(bin_num) +t _correl * sqrt(s_var(2) / s_var(1))
      reg_2_sum(bin_num) = reg_2_sum(bin_num) + (t_correl * sqrt(s_var(2) / s_var(1)))**2

      bin_count(bin_num) = bin_count(bin_num) + 1

   end do
end do

deallocate(pairs)


do i = 1, 201
   if(bin_count(i) > 0) then
      ! Compute the standard deviation of the true correlations
      true_correl_mean = tcorrel_sum(i) / bin_count(i)
      reg_mean = reg_sum(i) / bin_count(i)
      reg_sd = sqrt((reg_2_sum(i) - bin_count(i) * reg_mean**2) / (bin_count(i) - 1))

      if(reg_sd <= 0.0_r8) then
         alpha = 1.0_r8
      else
         !!!beta = reg_mean**2 / reg_sd**2
      ! CRAZY IDEA???
         beta = reg_mean**2 / (reg_sd**2 * (1.0_r8 + 1.0_r8 / ens_size))
         alpha = beta / (1.0_r8 + beta)
      endif

      write(*, *) 'bin, count, mean ', i, bin_count(i), true_correl_mean, alpha
   endif
end do
end program full_error



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


