! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Test to evaluate how to look_for_bias in prior ensemble distribution 
!> and observation in observation space.
!>
!> See work from 13 October, 2003.

program sys_sim301

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

type (random_seq_type) :: r
double precision :: y_o, sigma_y_o, y_p, sigma_y_p, dist_sq, sample_mean
double precision :: dist_sq_sum, mean_dist_sq
integer :: i, j, n, num, n_samples, count

! Initialize repeatable random sequence
call init_random_seq(r) 

count = 0

! For now have a set of free parameters that may be too large
write(*, *) 'Input prior variance of observation variable'
read(*, *) sigma_y_p
write(*, *) 'Input variance of observing instrument'
read(*, *) sigma_y_o
write(*, *) 'Input ensemble size'
read(*, *) n

write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples
sample_mean = 0.0

! Loop through the number of samples
do i = 1, n_samples

! Generate the observation
   y_o = random_gaussian(r, dble(0.0), sqrt(sigma_y_o))

   dist_sq_sum = 0.0
   do j = 1, n
      y_p = random_gaussian(r, dble(0.0), sqrt(sigma_y_p))
! Compute the distance between the samples and the obs
      dist_sq = (y_p - y_o) ** 2
!!!      dist_sq = abs(y_p - y_o)
      dist_sq_sum = dist_sq_sum + dist_sq
   end do
   mean_dist_sq = dist_sq_sum / n

   sample_mean = sample_mean + mean_dist_sq

end do

write(*, *) 'mean squared distance is ', sample_mean / n_samples

end program sys_sim301

