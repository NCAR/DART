module correct_corr_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! A module to do updates taking into account errors in correlation computation 
! resulting from finite sample size. See notes from second week of Dec. 2001.

use types_mod
use random_seq_mod, only : random_seq_type, init_random_seq, twod_gaussians

implicit none

private

logical :: first_call = .true.
!!!integer, parameter :: table_size = 5
integer, parameter :: table_size = 100
real(r8) :: table(0:table_size), inverse_table(0:table_size)
real(r8) :: var_table(0:table_size), inverse_var_table(0:table_size)
type (random_seq_type) :: r
integer, parameter :: n_mc = 10000

public init_table, corr_update_from_obs_inc, corr_obs_increment, get_correct_correlation

contains



  subroutine init_table(n)
!==============================================================================
! subroutine init_table(n)
!
! Subroutine to initialize a table with expected values of correlation given
! actual correlation for a given ensemble size. For small n this must be done
! by Monte Carlo since there are large errors in the standard formula.
! The inverse tables allow interpolation to get the expected value of the real
! correlation given the sample correlation and the expected value of the 
! variance of the sample correlation, given that the real variance is 
! given by the first table lookup. At some point, should think very hard about
! this again to make sure that this is the right quantity.

implicit none

integer, intent(in) :: n

real(r8) :: mean(2), c(2, 2), corr, rnum(2, n_mc), sample_correl
real(r8) :: correl_2
integer  :: i, j, k

! Initialize a repeatable random sequence for computing expected correlation

call init_random_seq(r)

table = 0.0_r8       ! Initialize the table for accumulation
mean  = 0.0           ! Mean is 0 for random sampling

do i = 0, table_size

   corr =  1.0_r8 * i / table_size

  ! Round-off can lead to corr greater than 1.0 which messes stuff up

   if(corr > 1.0_r8) corr = 1.0_r8

   ! Load the correlation matrix

   c(1, 1) = 1.0_r8; c(1, 2) =   corr; 
   c(2, 1) =   corr; c(2, 2) = 1.0_r8; 

   correl_2 = 0.0      ! Prepare to accumulate mean of squared correlation

   do j = 1, n_mc
      do k = 1, n      ! Loop through to do MC to get expected value of correlation
         call twod_gaussians(r, mean, c, rnum(:, k))
      end do
      call comp_correl(rnum, n, sample_correl)
      table(i) = table(i) + sample_correl 
      correl_2 = correl_2 + sample_correl**2
   end do

   table(i) = table(i) / n_mc

   if(table(i) < 0.0_r8) table(i) = 0.0_r8
   if(table(i) > 1.0_r8) table(i) = 1.0_r8

   var_table(i) = (correl_2 - n_mc * table(i)**2) / real(n_mc - 1)

   write(*, *) i, real(corr), real(table(i)), real(var_table(i)), &
               real((1 - table(i)**2)**2 / (n - 1))
end do

! Now, want to have a linear interpolation for going from sample to estimated continuous
! Generate table by doing interpolation into the forward table

do i = 0, table_size
   sample_correl = 1.0_r8 * i / real(table_size)
   if(sample_correl > 1.0_r8) sample_correl = 1.0_r8
   do j = 1, 100
      if(sample_correl <= table(j)) then

         inverse_table(i) = (sample_correl - table(j - 1)) / &
            (table(j) - table(j - 1)) / table_size + (j - 1.0_r8) / real(table_size)

         if(inverse_table(i) < 0.0_r8) inverse_table(i) = 0.0_r8
         if(inverse_table(i) > 1.0_r8) inverse_table(i) = 1.0_r8

         inverse_var_table(i) = ((sample_correl - table(j - 1)) / &
            (table(j) - table(j-1))) * (var_table(j) - var_table(j-1)) + var_table(j-1)

         if(inverse_var_table(i) < 0.0_r8) inverse_var_table(i) = 0.0_r8
         if(inverse_var_table(i) > 1.0_r8) inverse_var_table(i) = 1.0_r8

         write(*, *) i, real(sample_correl), real(inverse_table(i)), &
            real(inverse_var_table(i)), (1.0_r8 - inverse_table(i)**2)**2 / real(n - 1)
         goto 10
      endif
   end do
10 continue
end do

end subroutine init_table



  subroutine get_correct_correlation(sample_correl, ens_size, est_correl, var_correl)
!=============================================================================== 
! subroutine get_correct_correlation(sample_correl, ens_size, est_correl, var_correl)

implicit none

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: sample_correl
real(r8), intent(out) :: est_correl, var_correl

real(r8) :: c, sign, fract, interval
integer  :: upper, lower

! There's an intrinsic to do this sign stuff

if(sample_correl >= 0.0_r8) then
   sign = 1.0_r8
else
   sign = -1.0_r8
endif

c = abs(sample_correl)                      ! Lookup the value in the table
interval = 1.0_r8 / real(table_size - 1)    ! Interval in table

lower = int(c / (interval))
upper = lower + 1
fract = (c - lower * interval) / interval
est_correl = sign * &
   (inverse_table(lower) + fract * (inverse_table(upper) - inverse_table(lower)))

var_correl = inverse_var_table(lower) + &
   fract * (inverse_var_table(upper) - inverse_var_table(lower))

! WARNING: Need to make sure lookups are doing exactly what is wanted
!write(*, *) 's_correl ', real(sample_correl), real(est_correl), real(var_correl)

end subroutine get_correct_correlation



  subroutine corr_obs_increment(ens, ens_size, obs, obs_var, mean_inc, spread_inc)
!===============================================================================
! subroutine corr_obs_increment(ens, ens_size, obs, obs_var, mean_inc, spread_inc)
! 
! EAKF version of obs increment: Need updates for mean and spread separately
! To use corrected versions for sample correlation.
 
implicit none
 
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: mean_inc, spread_inc(ens_size)
 
real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8) :: prior_cov, sx, s_x2
integer  :: i

if(first_call) then
   call init_table(ens_size)
   first_call = .false.
end if
 
! Compute mt_rinv_y (obs error normalized by variance)

obs_var_inv = 1.0_r8 / obs_var
 
! Compute prior mean and covariance

sx            = sum(ens)
s_x2          = sum(ens * ens)
prior_mean    = sx / ens_size
prior_cov     = (s_x2 - sx**2 / ens_size) / real(ens_size - 1.0_r8)
prior_cov_inv = 1.0_r8 / prior_cov
new_cov       = 1.0_r8 / (prior_cov_inv + obs_var_inv)
new_mean      = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)
a             = sqrt(new_cov * prior_cov_inv)

! Need increments for both mean and spread for correction

spread_inc = (a - 1.0_r8) * (ens - prior_mean)
mean_inc   = new_mean - prior_mean
 
end subroutine corr_obs_increment


 
  subroutine corr_update_from_obs_inc(obs, mean_inc, spread_inc, state, &
                     ens_size, delta_state, cov_factor)
!===============================================================================
!  subroutine corr_update_from_obs_inc(obs, mean_inc, spread_inc, state, &
!                     ens_size, delta_state, cov_factor)
 
implicit none
 
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: obs(ens_size), mean_inc, spread_inc(ens_size)
real(r8), intent(in)  :: state(ens_size), cov_factor
real(r8), intent(out) :: delta_state(ens_size)
 
real(r8) :: ens(2, ens_size), cov(2, 2)
real(r8) :: sum_x, sum_y, sum_xy, sum_x2, reg_coef, creg_coef
real(r8) :: sum_y2, correl, est_correl, var_correl, varx, sample_factor
real(r8) :: temp_ens(2, ens_size), temp_correl
 
ens(1, :) = state      ! Create combined matrix for covariance
ens(2, :) = obs
 
! For efficiency, just compute regression coefficient here

sum_x  = sum(ens(2, :))
sum_y  = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))
 
reg_coef = (ens_size * sum_xy - sum_x * sum_y) / (ens_size * sum_x2 - sum_x**2)

! Computation of correlation

sum_y2 = sum(ens(1, :) * ens(1, :))
correl = (ens_size * sum_xy - sum_x * sum_y) / &
          sqrt((ens_size * sum_x2 - sum_x**2) * (ens_size * sum_y2 - sum_y**2))

correl = correl * cov_factor     ! TRY MODULATING FIRST

call get_correct_correlation(correl, ens_size, est_correl, var_correl)

!write(*, *) 'correl ', real(correl)

! First addition; don't do anything if expected spread increases
! See notes from 15 Dec. 2001

if(var_correl > est_correl**2 ) then

   !!!IF(Var_correl > 1000000) then
   delta_state = 0.0_r8

else

   ! First pass, just do exactly the same as before but with split terms for deltas
   !!! delta_state = cov_factor * reg_coef * spread_inc + cov_factor * reg_coef * mean_inc 
   creg_coef = est_correl * sqrt(ens_size * sum_y2 - sum_y**2) / &
               sqrt(ens_size * sum_x2 - sum_x **2)

   ! Need to modify distance by sqrt(1 - pe**2 / p**2) where p is correlation

   sample_factor = sqrt(1.0 - var_correl / est_correl **2)

   ! write(*, *) 'cov, sam, correl', real(cov_factor), real(sample_factor), real(correl)

   delta_state = sample_factor * creg_coef * spread_inc + &
                 sample_factor * creg_coef * mean_inc 

endif

! Ouput comparison

!if(var_correl <= est_correl **2) then
!   write(*, *) 'old new ', real(cov_factor * reg_coef), &
!                          real(sample_factor * creg_coef),  &
!                           real((cov_factor * reg_coef) / (sample_factor * creg_coef)), &
!                           real(correl)
!else
!   write(*, *) 'old new ', real(cov_factor * reg_coef), 0.0
!endif

! Temporary computation of updated correlation
!temp_ens = ens
!temp_ens(1, :) = temp_ens(1, :) + delta_state
!call comp_correl(temp_ens, ens_size, temp_correl)
!write(*, *) 'old new correl', real(correl), real(temp_correl), real(temp_correl / correl)

end subroutine corr_update_from_obs_inc



  subroutine comp_correl(ens, n, correl)
!===============================================================================
! subroutine comp_correl(ens, n, correl)
 
implicit none
 
integer,  intent(in)  :: n
real(r8), intent(in)  :: ens(2, n)
real(r8), intent(out) :: correl

real(r8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2
 
sum_x  = sum(ens(2, :))
sum_y  = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))
 
! Computation of correlation

sum_y2 = sum(ens(1, :) * ens(1, :))
 
correl = (n * sum_xy - sum_x * sum_y) / &
         sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))
 
end subroutine comp_correl

!================================================================================
! End of correct_corr_mod.f90
!================================================================================

end module correct_corr_mod
