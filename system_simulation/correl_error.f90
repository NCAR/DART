program correl_error

use types_mod
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &          
   twod_gaussians, random_uniform

implicit none

integer, parameter :: sample_size = 100000
type (random_seq_type) :: ran_id
real(r8) :: zero_2(2) = 0.0, cov(2, 2)
real(r8) :: t_correl, correl_mean, sample_correl
real(r8), allocatable :: pairs(:, :)
integer i, j, k, lji, ens_size

call init_random_seq(ran_id)

! Loop through a range of ensemble sizes
do lji = 4, 4
ens_size = 10 * 2**(lji - 1)
allocate(pairs(2, ens_size))
write(*, *) 'stats for ensemble size ', ens_size
! Loop through real correlations every 0.10
do j = 7, 20
   ! Loop through a large number of samples to get mean
   correl_mean = 0.0
   do k = 1, sample_size
      t_correl = j * 0.05

      ! Generate the covariance matrix for this correlation
      cov(1, 1) = 1.0;    cov(2, 2) = 1.0
      cov(1, 2) = t_correl; cov(2, 1) = t_correl
   
      ! Loop to generate an ensemble size sample from this correl
      ! Generate a random sample of size ens_size from something with this correlation
      do i = 1, ens_size
         call twod_gaussians(ran_id, zero_2, cov, pairs(:, i))
      end do
   
      ! Compute the sample correlation
      call comp_correl(pairs, ens_size, sample_correl)
   
      ! Keep a sum of the sample_correl
      correl_mean = correl_mean + abs(sample_correl)
   end do
   correl_mean = correl_mean / sample_size
   write(*, *) t_correl, correl_mean
end do
deallocate(pairs)
end do
end program correl_error



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

