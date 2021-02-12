! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Looking at sampling error for the obs. space part of the problem
!> Assume that the ratio of the prior sample to obs standard dev. is known
!> Want to compute the following:
!>   1. An unbiased estimate of the updated variance
!>   2. An unbiased estimate of the updated mean which should be
!>      expressed as the fraction of the distance to move between
!>      the prior sample mean and the observed value
!>   3. The standard deviation or variance of the updated mean estimate
!>      Need to add this uncertainty into the sample to avoid loss of
!>      variance (see also full_error.f90 which does this for regression).

program obs_sampling_err

use      types_mod, only : r8, digits12
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none

!!!integer, parameter :: sample_size = 10000000
integer, parameter :: sample_size = 100000
integer :: ens_size
real(r8), allocatable :: sample(:)
real(r8) :: s_var, s_mean, u_var, ratio
real(r8) :: u_var_sum, u_var_mean
real(r8) :: t_var, log_ratio, beta
real(r8) :: s_obs_weight, s_obs_weight_sum, s_obs_weight2_sum
real(r8) :: s_obs_weight_var, s_obs_weight_mean
real(r8) :: s_var_sum, s_var_mean
real(r8) :: mad, mad_sum, mad_mean, sd_sum, sd_mean
type (random_seq_type) :: ran_id
integer :: i, j, k

write(*, *) 'input ensemble size '
read(*, *) ens_size
allocate(sample(ens_size))

write(*, *) 'stats for ensemble size ', ens_size

! Initialize the random number sequence
call init_random_seq(ran_id)

! Loop through a variety of ratios to make a table
do k = 1, 200
   log_ratio = (k - 100.0_r8) / 50.0_r8
   ratio = 10.0_r8**log_ratio

   ! Compute the true updated mean and variance without sampling error
   ! True prior mean is 0.0, 
   ! True prior variance is ratio and true obs variance is 1.0
   t_var = 1.0_r8 / (1.0_r8 / ratio**2 + 1.0_r8)

   ! Initialize the sum arrays
   u_var_sum = 0.0_r8; s_var_sum = 0.0_r8
   s_obs_weight_sum = 0.0_r8; s_obs_weight2_sum = 0.0_r8
   mad_sum = 0.0_r8; sd_sum = 0.0_r8
   
   ! Loop through the sample
   do i = 1, sample_size
      ! Generate a random sample of size ens_size with ratio for variance
      do j = 1, ens_size
         sample(j) = random_gaussian(ran_id, 0.0_r8, ratio)
      enddo
      
      ! Compute the sample mean and variance
      call sample_mean_var(sample, ens_size, s_mean, s_var)



      ! Compute mean of sample variance
      s_var_sum = s_var_sum + s_var
      ! Compute mean of sd
      sd_sum = sd_sum + sqrt(s_var)
      ! Compute the mean absolute deviation, too
      mad = sum(abs(sample - s_mean)) / ens_size
      mad_sum = mad_sum + mad


   
      ! Compute the updated standard deviation and mean; obs var is 1.0
      u_var = 1.0_r8 / (1.0_r8 / s_var + 1.0_r8)
      u_var_sum = u_var_sum + u_var

      ! What weight is given to obs here?
      s_obs_weight = u_var / 1.0_r8
      s_obs_weight_sum = s_obs_weight_sum + s_obs_weight
      s_obs_weight2_sum = s_obs_weight2_sum + s_obs_weight**2
   
   enddo
   
   ! Output the mean stats
   u_var_mean = u_var_sum / sample_size
   s_obs_weight_mean = s_obs_weight_sum / sample_size
   s_obs_weight_var = (s_obs_weight2_sum - s_obs_weight_sum**2 / sample_size) / &
      (sample_size - 1.0_r8)
   s_var_mean = s_var_sum / sample_size
   sd_mean = sd_sum / sample_size
   mad_mean = mad_sum / sample_size
   
   ! Compute beta 
   beta = s_obs_weight_mean**2 / s_obs_weight_var
  
   write(*, 22) log_ratio, ratio, t_var / u_var_mean, t_var / s_obs_weight_mean, beta
   22 format(1x, 6(e10.4, 1x))
   write(*, *) 's_var_mean ', s_var_mean, ratio**2
   write(*, *) 'mad ', mad_mean, mad_mean / sqrt(s_var_mean)
   write(*, *) 'sd_mean ', sd_mean, sqrt(s_var_mean), sd_mean / sqrt(s_var_mean)

enddo

!----------------------------------------------------------------
contains
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

end program obs_sampling_err

