! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program correl_error

use      types_mod, only : r8
use random_seq_mod, only : random_seq_type, init_random_seq, twod_gaussians 

implicit none

!!!integer, parameter :: sample_size = 1000000
integer, parameter :: sample_size = 100000
type (random_seq_type) :: ran_id
real(r8) :: zero_2(2) = 0.0, cov(2, 2)
real(r8) :: t_correl, correl_mean, sample_correl, ratio_mean, correl_mean2, correl_sd
real(r8) :: ratio_mean2, ratio_sd
real(r8), allocatable :: pairs(:, :), obs_inc(:), unobs_inc(:), new_unobs(:)
real(r8) :: s_mean(2), s_var(2), new_mean, new_var, a, expected_ratio
integer i, j, k, lji, ens_size

call init_random_seq(ran_id)

! Loop through a range of ensemble sizes
do lji = 4, 48, 4
   ens_size = lji
   if(lji > 32) ens_size = 32 + (lji - 32) * 4
   if(lji > 40) ens_size = 64 + (lji - 40) * 8
   !!!ens_size = 2**lji

   !!!ens_size = 20
   !!!ens_size = 128

   allocate(pairs(2, ens_size), obs_inc(ens_size), unobs_inc(ens_size), new_unobs(ens_size))
   write(*, *) 'stats for ensemble size ', ens_size
   ! Loop through real correlations every 0.05
   !!!do j = 0, 50
   do j = 0, 0
      ! Loop through a large number of samples to get mean
      correl_mean  = 0.0_r8
      correl_mean2 = 0.0_r8
      ratio_mean   = 0.0_r8
      ratio_mean2  = 0.0_r8
      do k = 1, sample_size
         t_correl = j * 0.02_r8

         ! Generate the covariance matrix for this correlation
         cov(1, 1) = 1.0;    cov(2, 2) = 1.0_r8
         cov(1, 2) = t_correl; cov(2, 1) = t_correl

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
         ! Next, compute the increments for obs variable given reduction in s.d. by a
         a = 0.0_r8
         do i = 1, ens_size
            obs_inc(i) = (a - 1.0_r8) * (pairs(1, i) - s_mean(1))
         enddo

         ! Regress these onto the unobserved variable
         do i = 1, ens_size
            unobs_inc(i) = sample_correl * sqrt(s_var(2) / s_var(1)) * obs_inc(i)
            new_unobs(i) = pairs(2, i) + unobs_inc(i)
         enddo

         ! Now compute the updated variance
         call sample_mean_var(new_unobs, ens_size, new_mean, new_var)

         !write(*, *) 'old and new var ', s_var(2), new_var
         !write(*, *) 'old and new sd ', sqrt(s_var(2)), sqrt(new_var)
         !write(*, *) 'unobs sd ratio ', sqrt(new_var) / sqrt(s_var(2))
         ratio_mean  = ratio_mean  +  sqrt(new_var) / sqrt(s_var(2))
         ratio_mean2 = ratio_mean2 + (sqrt(new_var) / sqrt(s_var(2)))**2

         !-----------------

      enddo
      correl_mean = correl_mean / sample_size
      correl_sd = sqrt((correl_mean2 - sample_size * correl_mean**2) / (sample_size - 1))
      write(*, *) t_correl, correl_mean, correl_sd
      ratio_mean = ratio_mean / sample_size
      ratio_sd = sqrt((ratio_mean2 - sample_size * ratio_mean**2) / (sample_size - 1))
      write(*, *) 'ratio mean ', ratio_mean, (1 - ratio_mean) / (1.0 - a)
      expected_ratio = 1.0 - ((1.0 - a) * t_correl)
      write(44, 11) ens_size, t_correl, ratio_mean, ratio_sd, correl_mean, correl_sd
      11 format(i8, 6(1x, f10.5))
   enddo
   deallocate(pairs, obs_inc, unobs_inc, new_unobs)
enddo

contains

!-----------------------------------------------------

subroutine comp_correl(ens, n, correl)

implicit none

integer, intent(in) :: n
real(r8), intent(in) :: ens(2, n)
real(r8), intent(out) :: correl
real(r8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2


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
real(r8), intent(in) :: x(:)
real(r8), intent(out) :: mean, var

real(r8) :: sx, s_x2

sx = sum(x)
s_x2 = sum(x * x)
mean = sx / n
var = (s_x2 - sx**2 / n) / (n - 1)

end subroutine sample_mean_var

end program correl_error

