! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sampling_error

use      types_mod, only : r8, digits12
use random_seq_mod, only : random_seq_type, init_random_seq, twod_gaussians

implicit none

integer, parameter :: sample_size = 1000000
type (random_seq_type) :: ran_id
real(r8) :: zero_2(2) = 0.0, cov(2, 2)
real(r8) :: t_correl, correl_mean, sample_correl, ratio_mean, correl_mean2, correl_sd
real(r8), allocatable :: pairs(:, :), obs_inc(:), unobs_inc(:), new_unobs(:)
real(r8) :: s_mean(2), s_var(2)
real(r8) :: reg_sum, reg_mean, reg_sum_2, reg_sd, t_sd1, t_sd2, alpha, beta
integer i, j, k, lji, ens_size

call init_random_seq(ran_id)

! Loop through a range of ensemble sizes
do lji = 2, 2
   ens_size = 10 * 2**(lji - 1)

   ens_size = 10

   allocate(pairs(2, ens_size), obs_inc(ens_size), unobs_inc(ens_size), new_unobs(ens_size))
   write(*, *) 'stats for ensemble size ', ens_size

   ! Loop through real correlations every 0.02
   do j = 0, 50
   !do j = 0, 0

      ! Loop through a large number of samples to get mean
      correl_mean  = 0.0_r8
      correl_mean2 = 0.0_r8
      ratio_mean   = 0.0_r8
      reg_sum      = 0.0_r8
      reg_sum_2    = 0.0_r8

      do k = 1, sample_size
         t_correl = j * 0.02_r8

         ! NOTE: RESULTS DO NOT DEPEND ON SDs of the two priors, just on correl
         t_sd1 = 1.0_r8
         t_sd2 = 1.0_r8

         ! Generate the covariance matrix for this correlation
         cov(1, 1) = t_sd1**2;    cov(2, 2) = t_sd2**2
         cov(1, 2) = t_correl * t_sd1 * t_sd2; cov(2, 1) = cov(1, 2)

         ! Loop to generate an ensemble size sample from this correl
         ! Generate a random sample of size ens_size from something with this correlation
         do i = 1, ens_size
            call twod_gaussians(ran_id, zero_2, cov, pairs(:, i))
         enddo

         ! Compute the sample correlation
         call comp_correl(pairs, ens_size, sample_correl)

         ! Keep a sum of the sample_correl
         correl_mean = correl_mean + abs(sample_correl)
         correl_mean2 = correl_mean2 + sample_correl**2

         !-----------------
         ! Also interested in finding out what the spurious variance reduction factor is
         ! First need to compute sample s.d. for the obs and unobs variable
         do i = 1, 2
            call sample_mean_var(pairs(i, :), ens_size, s_mean(i), s_var(i))
         enddo

         !-----------------
         ! Also compute statistics for regression coefficient
         reg_sum = reg_sum + sample_correl * sqrt(s_var(2) / s_var(1))
         reg_sum_2 = reg_sum_2 + (sample_correl * sqrt(s_var(2) / s_var(1)))**2

      enddo
      correl_mean = correl_mean / sample_size
      correl_sd = sqrt((correl_mean2 - sample_size * correl_mean**2) / (sample_size - 1))
      ratio_mean = ratio_mean / sample_size
      !write(*, *) 'ratio mean ', ratio_mean

      ! Do regression stats, too
      reg_mean = reg_sum / sample_size
      reg_sd = sqrt((reg_sum_2 - sample_size * reg_mean**2) / (sample_size - 1))
      if(reg_sd <= 0.0_r8) then
         alpha = 1.0_r8
      else
         beta = reg_mean**2 / reg_sd**2
         alpha = beta / (1.0_r8 + beta)
      endif
      write(*, 11) t_correl, correl_mean, t_correl / correl_mean, alpha
      11 format(1x, 4(e10.4, 1x))

   enddo
   deallocate(pairs)

enddo

!-----------------------------------------------------
contains
!-----------------------------------------------------

subroutine comp_correl(ens, n, correl)

integer,        intent(in) :: n
real(digits12), intent(in) :: ens(2, n)
real(digits12), intent(out) :: correl

real(digits12) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2

sum_x  = sum(ens(2, :))
sum_y  = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))

! Computation of correlation
sum_y2 = sum(ens(1, :) * ens(1, :))

correl = (n * sum_xy - sum_x * sum_y) / &
   sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))

end subroutine comp_correl

!----------------------------------------------------------------

subroutine sample_mean_var(x, n, mean, var)

integer,        intent(in) :: n
real(digits12), intent(in) :: x(n)
real(digits12), intent(out) :: mean, var

real(digits12) :: sx, s_x2

sx   = sum(x)
s_x2 = sum(x * x)
mean = sx / n
var  = (s_x2 - sx**2 / n) / (n - 1)

end subroutine sample_mean_var


end program sampling_error

