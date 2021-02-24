! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program sys_sim104

! WORKS ON ONE APPROXIMATION FOR THE THIRD MOMENT STATISTICAL CLOSURE
! FIRST TEST IF FOR SQUARED ERROR MINIMIZATION.
! HAVE TO MAKE SOME ASSUMPTIONS ABOUT THE BACKGROUND DISTRIBUTION OF
! REGRESSION S.D. / MEAN TO GET THE ERROR FOR SMALL GROUP SIZES. FOR
! FIRST TRY, JUST ASSUME THAT ONE CAN SAMPLE FROM THE SAMPLE DISTRIBUTION
! WITH THE GROUP SIZE

! USE THE FORMULA 1 / (ratio**2 + 1) to get factor and compute mean factor

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &
   twod_gaussians, random_uniform

implicit none

type (random_seq_type) :: r
double precision :: sd_ratio, mean_move, mean_move_2, alpha, ratio
double precision, allocatable :: rnum(:), var_hist(:)
double precision :: sample_sigma_y_p, base_mean
double precision :: c(2, 2), sum_err_var, sigma_y_p, sigma_y_o, sigma_x_p
double precision :: y_truth, y_o, sample_correl, reg_coef, sample_reg_coef
double precision :: sigma_y_u, y_u, sigma_x_u, x_u, sample_x_u, x_error
double precision :: x_error_var, total_err_var, mean(2), correl, est_err_var
double precision :: temp_sum, pe_2_sum, xx, s_xx, s_xx_2, xx_mean, xx_var
double precision, allocatable :: sample_mean(:), x_var(:)
double precision :: sum_sample_mean, sum_sample_sigma, sum_y_u, sum_sigma_y_u, sum_y_u_2
double precision :: sum_err, sum_err_2, inf_sigma_y_u, var_mean, var_var
integer :: n, i, j, n_samples, index

! Initialize repeatable random sequence
call init_random_seq(r) 

! Free parameters are ratio and ensemble size
! Just assume that base value is 1.0 (move is 1.0); can't do case where expected move is 0
! but know the answer to this case (DON"T MOVE!)
write(*, *) 'Input sample ratio of ensemble standard deviation to mean'
read(*, *) sd_ratio
write(*, *) 'Input num groups'
read(*, *) n
write(*, *) 'Input number of samples for statistics'
read(*, *) n_samples

! Allocate storage for the sample means
allocate(sample_mean(n_samples), rnum(n))

mean_move = 0.0

! Loop through the number of samples
do i = 1, n_samples

! This first part approximates sampling error from the num groups!!!
! Need an ensemble sized sample of the recursion value
   do j = 1, n
      rnum(j) = random_gaussian(r, dble(1.0), sd_ratio)
   end do

! Compute the sample mean and variance
   call sample_mean_var(rnum, n, base_mean, sample_sigma_y_p)
   ratio =  sqrt(sample_sigma_y_p) / abs(base_mean)

!write(*, *) i, 'ratio ', ratio

! Compute the factor for this case and sum
   mean_move = mean_move + 1.0 / ( ratio **2 + 1.0)

end do

! Expected value of alpha to minimize squared difference between sample and actual value
alpha = mean_move / n_samples

write(*, *) 'alpha is ', alpha

end program sys_sim104

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

