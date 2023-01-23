! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! A variety of PDFs, CDFs, quantile functions and other tools for working with distributions
! to implement quantile conserving filters in observation space and regression in quantile space.

module quantile_distributions_mod

use types_mod, only : r8, digits12, PI

use sort_mod,  only : sort, index_sort

use utilities_mod, only : E_ERR, error_handler

use algorithm_info_mod, only : probit_dist_info, NORMAL_PRIOR, BOUNDED_NORMAL_RH_PRIOR,  &
                               GAMMA_PRIOR, BETA_PRIOR, LOG_NORMAL_PRIOR, UNIFORM_PRIOR
                               !!!PARTICLE_PRIOR

use normal_distribution_mod, only : norm_cdf, norm_inv, weighted_norm_inv

use gamma_distribution_mod, only : gamma_cdf, inv_gamma_cdf

use beta_distribution_mod,  only : beta_cdf,  inv_beta_cdf

implicit none
private


public :: convert_to_probit, convert_from_probit, convert_all_to_probit, &
   convert_all_from_probit, dist_param_type

type dist_param_type
   integer               :: prior_distribution_type
   real(r8), allocatable :: params(:)
end type

! Saves the ensemble size used in the previous call of obs_inc_bounded_norm_rh
integer :: bounded_norm_rh_ens_size = -99

character(len=512)     :: errstring
character(len=*), parameter :: source = 'quantile_distributions_mod.f90'

! Logical to fix bounds violations for bounded_normal_rh
logical :: fix_bound_violations = .false.

contains

!------------------------------------------------------------------------

subroutine convert_all_to_probit(ens_size, num_vars, state_ens, prior_distribution_type, &
   p, probit_ens, use_input_p, bounded, bounds)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: state_ens(:, :)
integer, intent(in)                  :: prior_distribution_type(num_vars)
type(dist_param_type), intent(inout) :: p(num_vars)
real(r8), intent(out)                :: probit_ens(:, :)
logical, intent(in)                  :: use_input_p
logical, intent(in)                  :: bounded(2)
real(r8), intent(in)                 :: bounds(2)


! NOTE THAT WILL MAKE HELEN CRAZY: THIS WORKS WITH THE INPUT CALLING ARGUMENTS FOR STATE_ENS AND
! PROBIT_ENS BEING THE SAME. A TEMP IS USED TO AVOID OVERWRITING ISSUES. IS THIS YUCKY?

! Note that the input and output arrays may have extra copies (first subscript). Passing sections of a
! leading index could be inefficient for time and storage, so avoiding that for now.

! Assumes that the bounds are the same for any variables that are BNRH for now
! The bounds variables are not used for the normal case or the case where the input p is used

integer  :: i
real(r8) :: temp_ens(ens_size)

do i = 1, num_vars
   call convert_to_probit(ens_size, state_ens(1:ens_size, i), prior_distribution_type(i), &
      p(i), temp_ens, use_input_p, bounded, bounds)
   probit_ens(1:ens_size, i) = temp_ens
end do

end subroutine convert_all_to_probit

!------------------------------------------------------------------------

subroutine convert_to_probit(ens_size, state_ens, prior_distribution_type, p, &
   probit_ens, use_input_p, bounded, bounds)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
integer, intent(in)                  :: prior_distribution_type
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
logical, intent(in)                  :: bounded(2)
real(r8), intent(in)                 :: bounds(2)

! Set the type of the distribution in the parameters defined type
p%prior_distribution_type = prior_distribution_type

if(p%prior_distribution_type == NORMAL_PRIOR) then 
   call to_probit_normal(ens_size, state_ens, p, probit_ens, use_input_p)
elseif(p%prior_distribution_type == LOG_NORMAL_PRIOR) then 
   call to_probit_log_normal(ens_size, state_ens, p, probit_ens, use_input_p)
elseif(p%prior_distribution_type == UNIFORM_PRIOR) then 
   call to_probit_uniform(ens_size, state_ens, p, probit_ens, use_input_p, bounds)
elseif(p%prior_distribution_type == GAMMA_PRIOR) then 
   call to_probit_gamma(ens_size, state_ens, p, probit_ens, use_input_p, bounded, bounds)
elseif(p%prior_distribution_type == BETA_PRIOR) then 
   call to_probit_beta(ens_size, state_ens, p, probit_ens, use_input_p, bounded, bounds)
elseif(p%prior_distribution_type == BOUNDED_NORMAL_RH_PRIOR) then
   call to_probit_bounded_normal_rh(ens_size, state_ens, p, probit_ens, &
      use_input_p, bounded, bounds)
!!!elseif(p%prior_distribution_type == PARTICLE_PRIOR) then
   !!!call to_probit_particle(ens_size, state_ens, p, probit_ens, use_input_p, bounded, bounds)
else
   write(errstring, *) 'Illegal distribution type', p%prior_distribution_type
   call error_handler(E_ERR, 'convert_to_probit', errstring, source)
endif

end subroutine convert_to_probit

!------------------------------------------------------------------------

subroutine to_probit_normal(ens_size, state_ens, p, probit_ens, use_input_p)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p

! Don't need to do anything for normal
probit_ens = state_ens

end subroutine to_probit_normal

!------------------------------------------------------------------------

subroutine to_probit_log_normal(ens_size, state_ens, p, probit_ens, use_input_p)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p

! Taking the logarithm leads directly to a normal distribution
! This normal may not be standard normal, but needs no further adjustment like 
! the regular normal
probit_ens = log(state_ens)

end subroutine to_probit_log_normal

!------------------------------------------------------------------------

subroutine to_probit_uniform(ens_size, state_ens, p, probit_ens, use_input_p, bounds)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
real(r8), intent(in)                 :: bounds(2)

real(r8) :: lower_bound, upper_bound, range, quantile
integer :: i

if(use_input_p) then
   lower_bound = p%params(1)
   upper_bound = p%params(2)
else
   lower_bound = bounds(1)
   upper_bound = bounds(2)   
   if(.not. allocated(p%params)) allocate(p%params(2))
   p%params(1) = lower_bound
   p%params(2) = upper_bound
endif

range = upper_bound - lower_bound
do i = 1, ens_size
   ! Convert to quantile; U(lower_bound, upper_bound) to U(0, 1)
   quantile = (state_ens(i) - lower_bound) / range
   ! Convert to probit space 
   call norm_inv(quantile, probit_ens(i))
end do

end subroutine to_probit_uniform

!------------------------------------------------------------------------

subroutine to_probit_gamma(ens_size, state_ens, p, probit_ens, use_input_p, &
   bounded, bounds)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
logical, intent(in)                  :: bounded(2)
real(r8), intent(in)                 :: bounds(2)

! Probit transform for gamma.
real(r8) :: mean, sd, variance, shape, scale, quantile
integer  :: i

! In full generality, gamma must be bounded either below or above
if(.not. (bounded(1) .neqv. bounded(2))) then
   errstring = 'Gamma distribution requires either bounded above or below to be true'
   call error_handler(E_ERR, 'to_probit_gamma', errstring, source)
endif

! Get parameters
! Representing gamma in terms of shape and scale. 
if(use_input_p) then
   shape = p%params(1)
   scale = p%params(2)
else
   mean = sum(state_ens) / ens_size
   sd  = sqrt(sum((state_ens - mean)**2) / (ens_size - 1))
   variance = sd**2
   ! Get shape and scale
   shape = mean**2 / variance
   scale = variance / mean
   if(.not. allocated(p%params)) allocate(p%params(2))
   p%params(1) = shape
   p%params(2) = scale
endif

do i = 1, ens_size
   ! First, convert the ensemble member to quantile
   quantile = gamma_cdf(state_ens(i), shape, scale)
   ! Convert to probit space 
   call norm_inv(quantile, probit_ens(i))
end do

end subroutine to_probit_gamma

!------------------------------------------------------------------------

subroutine to_probit_beta(ens_size, state_ens, p, probit_ens, use_input_p, &
   bounded, bounds)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
logical, intent(in)                  :: bounded(2)
real(r8), intent(in)                 :: bounds(2)

! Probit transform for beta.
real(r8) :: mean, sd, variance, alpha, beta, quantile, lower_bound, upper_bound
integer  :: i

! For now, check to make sure that distribution is bounded above and below
if(.not. (bounded(1) .and. bounded(2))) then
   errstring = 'Beta distribution requires bounded below and above to be true'
   call error_handler(E_ERR, 'to_probit_beta', errstring, source)
endif

! Get parameters
! Representing beta in terms of alpha and beta
if(use_input_p) then
   alpha = p%params(1)
   beta  = p%params(2)
   ! Bounds for translation and scaling
   lower_bound = p%params(3)
   upper_bound = p%params(4)
   ! Translate and scale the ensemble so it is on [0 1], use the output probit_ens for temp storage
   probit_ens = (state_ens - lower_bound) / (upper_bound - lower_bound)
else
   if(.not. allocated(p%params)) allocate(p%params(4))
   lower_bound = bounds(1)
   upper_bound = bounds(2)
   ! Translate and scale the ensemble so it is on [0 1], use the output probit_ens for temp storage
   probit_ens = (state_ens - lower_bound) / (upper_bound - lower_bound)
   mean = sum(probit_ens) / ens_size
   sd  = sqrt(sum((probit_ens - mean)**2) / (ens_size - 1))
   variance = sd**2
   ! Get alpha and beta
   alpha = mean**2 * (1.0_r8 - mean) / variance - mean
   beta  = alpha * (1.0_r8 / mean - 1.0_r8)
   p%params(1) = alpha
   p%params(2) = beta
   p%params(3) = lower_bound
   p%params(4) = upper_bound
endif

do i = 1, ens_size
   ! First, convert the ensemble member to quantile
   quantile = beta_cdf(probit_ens(i), alpha, beta)
   ! Convert to probit space 
   call norm_inv(quantile, probit_ens(i))
end do

end subroutine to_probit_beta

!------------------------------------------------------------------------

subroutine to_probit_bounded_normal_rh(ens_size, state_ens, p, probit_ens, &
   use_input_p, bounded, bounds)

! Note that this is just for transforming back and forth, not for doing the RHF observation update
! This means that we know a prior that the quantiles associated with the initial ensemble are
! uniformly spaced which can be used to simplify converting.

! How to handle identical ensemble members is an open question for now. This is also a problem
! for ensemble members that are identical to one of the bounds. 

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
logical, intent(in)                  :: bounded(2)
real(r8), intent(in)                 :: bounds(2)

! Probit transform for bounded normal rh.
integer  :: i, j, indx, low_num, up_num
integer  :: ens_index(ens_size)
real(r8) :: x, quantile, q(ens_size)
logical  :: bounded_below, bounded_above, do_uniform_tail_left, do_uniform_tail_right
real(r8) :: lower_bound, tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: upper_bound, tail_amp_right, tail_mean_right, tail_sd_right

! Parameter to control switch to uniform approximation for normal tail
! This defines how many quantiles the bound is from the outermost ensemble member
! If closer than this, can get into precision error problems with F(F-1(X)) on tails
! Can also get other precision problems with the amplitude for the normal in the bounded region
real(r8), parameter :: uniform_threshold = 0.01_r8

! Save to avoid a modestly expensive computation redundancy
real(r8), save :: dist_for_unit_sd
real(r8) :: mean, sd, base_prob, bound_quantile

if(use_input_p) then
   ! Using an existing ensemble for the RH points
   tail_sd_left = p%params(ens_size + 11)
   
   ! Don't know what to do if sd of original ensemble is 0 (or small, work on this later)
   if(tail_sd_left <= 0.0_r8) then
      ! Just return the original ensemble
      probit_ens = state_ens 
      return
   endif

   ! Get rest of variables out of the parameter storage for clarity
   bounded_below = p%params(ens_size + 1) > 0.5_r8
   bounded_above = p%params(ens_size + 2) > 0.5_r8
   lower_bound = p%params(ens_size + 3)
   upper_bound = p%params(ens_size + 4)
   do_uniform_tail_left = p%params(ens_size + 5) > 0.5_r8
   do_uniform_tail_right = p%params(ens_size + 6) > 0.5_r8
   tail_amp_left = p%params(ens_size + 7)
   tail_amp_right = p%params(ens_size + 8)
   tail_mean_left = p%params(ens_size + 9)
   tail_mean_right = p%params(ens_size + 10)
   tail_sd_right = p%params(ens_size + 12)

   ! Get the quantiles for each of the ensemble members in a RH distribution
   call ens_quantiles(p%params(1:ens_size), ens_size, &
      bounded_below, bounded_above, lower_bound, upper_bound, q)

   ! This can be done vastly more efficiently with either binary searches or by first sorting the
   ! incoming state_ens so that the lower bound for starting the search is updated with each ensemble member
   do i = 1, ens_size
      ! Figure out which bin it is in
      x = state_ens(i)

      ! This block of the code is only applied for observation posteriors
      ! Rare round-off issues with the bounded_normal RHF can lead to increments that produce posteriors
      ! that can violate the bound by a tiny amount. The following statements can fix that. For now, they are 
      ! turned off in the default code so that more egregious possible errors can be flagged below.
      if(fix_bound_violations) then
         if(bounded_below .and. x < lower_bound) x = lower_bound
         if(bounded_above .and. x > upper_bound) x = upper_bound
      endif

      if(x < p%params(1)) then
         ! In the left tail
         ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
         if(bounded_below .and. x < lower_bound) then
            write(errstring, *) 'Ensemble member less than lower bound first check(see code)', x, lower_bound
            call error_handler(E_ERR, 'to_probit_bounded_normal_rh', errstring, source)
            ! This error can occur due to roundoff in increment generation from bounded RHF
            ! It may be safe to just set x to the value of the bound here. See fix_bound_violations above
         endif

         if(do_uniform_tail_left) then
            ! Uniform approximation for left tail
            ! The division here could be a concern. However, if p%params(1) == lower_bound, then
            ! x cannot be < p%params(1).
            quantile = (x - lower_bound) / (p%params(1) - lower_bound) * (1.0_r8 / (ens_size + 1.0_r8))
         else
            ! It's a normal tail
            if(bounded_below) then
               quantile = tail_amp_left * (norm_cdf(x, tail_mean_left, tail_sd_left) - &
                  norm_cdf(lower_bound, tail_mean_left, tail_sd_left))
            else        ! Unbounded, tail normal goes all the way down to quantile 0
               quantile = tail_amp_left * norm_cdf(x, tail_mean_left, tail_sd_left)
            endif
            ! Make sure it doesn't sneak past the first ensemble member due to round-off
            quantile = min(quantile, 1.0_r8 / (ens_size + 1.0_r8))
         endif
      elseif(x == p%params(1)) then
         ! This takes care of cases where there are multiple rh values at the bdry or at first ensemble
         quantile = q(1)
      elseif(x > p%params(ens_size)) then
         ! In the right tail
         ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
         if(bounded_above .and. x > upper_bound) then
            write(errstring, *) 'Ensemble member greater than upper bound first check(see code)', x, upper_bound
            call error_handler(E_ERR, 'to_probit_bounded_normal_rh', errstring, source)
            ! This error can occur due to roundoff in increment generation from bounded RHF
            ! It may be safe to just set x to the value of the bound here. See fix_bound_violations above
         endif

         if(do_uniform_tail_right) then
            ! Uniform approximation for right tail
            ! The division here could be a concern. However, if p%params(ens_size) == upper_bound, then
            ! x cannot be > p%params(ens_size).
            quantile = ens_size / (ens_size + 1.0_r8) + &
               (x - p%params(ens_size)) / (upper_bound - p%params(ens_size)) * (1.0_r8 / (ens_size + 1.0_r8))
         else
            ! It's a normal tail
            ! Want to avoid quantiles exceeding 1 due to numerical issues. Do fraction of the normal part
               quantile = ens_size / (ens_size + 1.0_r8) + &
                  tail_amp_right * (norm_cdf(x, tail_mean_right, tail_sd_right) - &
                  norm_cdf(p%params(ens_size), tail_mean_right, tail_sd_right)) * (1.0_r8 / (ens_size + 1.0_r8))
               quantile = min(quantile, 1.0_r8)
         endif

      else
         ! In an interior bin
         do j = 1, ens_size - 1
            if(x < p%params(j+1)) then
               ! The division here could be a concern. 
               ! However, p%params(j)< x < p%params(j+1) so the two cannot be equal
               quantile = (j * 1.0_r8) / (ens_size + 1.0_r8) + &
                  ((x - p%params(j)) / (p%params(j+1) - p%params(j))) * (1.0_r8 / (ens_size + 1.0_r8))
               exit
            elseif(x == p%params(j+1)) then
               quantile = q(j+1)
               exit
            endif
         enddo
      endif
      ! Convert to probit space 
      call norm_inv(quantile, probit_ens(i))
   end do
else

   ! Take care of space for the transform data structure
   if(allocated(p%params)) deallocate(p%params)
   allocate(p%params(ens_size + 2*6))

   ! No pre-existing distribution, create one
   mean = sum(state_ens) / ens_size
   sd  = sqrt(sum((state_ens - mean)**2) / (ens_size - 1))
   
   ! Don't know what to do if sd is 0 (or small, work on this later)
   if(sd <= 0.0_r8) then
      ! Store this info in the left_tail_sd (parameter 11 in structure) for possible subsequent call use
      p%params(ens_size + 11) = sd
      ! Just return the original ensemble
      probit_ens = state_ens 
      return
   endif

   ! Clarity of use for bounds
   lower_bound = bounds(1)
   upper_bound = bounds(2)
   bounded_below = bounded(1)
   bounded_above = bounded(2)

   ! Need to sort. For now, don't worry about efficiency, but may need to somehow pass previous
   ! sorting indexes and use a sort that is faster for nearly sorted data. Profiling can guide the need
   call index_sort(state_ens, ens_index, ens_size)
   p%params(1:ens_size) = state_ens(ens_index)

   ! Get the quantiles for each of the ensemble members in a RH distribution
   call ens_quantiles(p%params(1:ens_size), ens_size, &
      bounded_below, bounded_above, lower_bound, upper_bound, q)

   ! Convert the quantiles to probit space
   do i = 1, ens_size
      indx = ens_index(i)
      call norm_inv(q(i), probit_ens(indx))
   end do 

   ! For BNRH, the required data for inversion is the original ensemble values
   ! Having them in sorted order is useful for subsequent inversion
   ! It is also useful to store additional information regarding the continuous pdf representation of the tails
   ! This includes whether the bounds are defined, the values of the bounds, whether a uniform is used in the outer
   ! bounded bin, the amplitude of the outer continuous normal pdf, the mean of the outer continous
   ! normal pdf, and the standard deviation of the
   ! outer continous. 

   ! Compute the description of the tail continous pdf; 
   ! First two entries are 'logicals' 0 for false and 1 for true indicating if bounds are in use
   if(bounded_below) then
      p%params(ens_size + 1) = 1.0_r8
   else
      p%params(ens_size + 1) = 0.0_r8
   endif

   if(bounded_above) then
      p%params(ens_size + 2) = 1.0_r8
   else
      p%params(ens_size + 2) = 0.0_r8
   endif

   ! Store the bounds (whether used or not) in the probit conversion metadata
   p%params(ens_size + 3) = lower_bound
   p%params(ens_size + 4) = upper_bound

   ! Compute the characteristics of unbounded tail normals

   ! For unit normal, find distance from mean to where cdf is 1/(ens_size+1).
   ! Saved to avoid redundant computation for repeated calls with same ensemble size
   if(bounded_norm_rh_ens_size /= ens_size) then
      call norm_inv(1.0_r8 / (ens_size + 1.0_r8), dist_for_unit_sd)
      ! This will be negative, want it to be a distance so make it positive
      dist_for_unit_sd = -1.0_r8 * dist_for_unit_sd
      ! Keep a record of the ensemble size used to compute dist_for_unit_sd
      bounded_norm_rh_ens_size = ens_size
   endif
   
   ! Fail if lower bound is larger than smallest ensemble member 
   if(bounded_below) then
      ! Do in two ifs in case the bound is not defined
      if(p%params(1) < lower_bound) then
         errstring = 'Ensemble member less than lower bound'
         call error_handler(E_ERR, 'to_probit_bounded_normal_rh', errstring, source)
      endif
   endif
   
   ! Fail if upper bound is smaller than the largest ensemble member 
   if(bounded_above) then
      if(p%params(ens_size) > upper_bound) then
         errstring = 'Ensemble member greater than upper bound'
         call error_handler(E_ERR, 'to_probit_bounded_normal_rh', errstring, source)
      endif
   endif

   ! Find a mean so that 1 / (ens_size + 1) probability is in outer regions
   tail_mean_left = p%params(1) + dist_for_unit_sd * sd
   tail_mean_right = p%params(ens_size) - dist_for_unit_sd * sd

   ! If the distribution is bounded, still want 1 / (ens_size + 1) in outer regions
   ! Put an amplitude term (greater than 1) in front of the tail normals 
   ! Amplitude is 1 if there are no bounds, so start with that
   tail_amp_left  = 1.0_r8
   tail_amp_right = 1.0_r8

   ! DO SOMETHING TO AVOID CASES WHERE THE BOUND AND THE SMALLEST ENSEMBLE ARE VERY CLOSE/SAME
   ! Default: not close
   do_uniform_tail_left = .false.
   base_prob = 1.0_r8 / (ens_size + 1.0_r8)
   if(bounded_below) then
      ! Compute the CDF at the bounds
      bound_quantile = norm_cdf(lower_bound, tail_mean_left, sd)
      if(abs(base_prob - bound_quantile) < uniform_threshold) then
         ! If bound and ensemble member are too close, do uniform approximation
         do_uniform_tail_left = .true.
      else
         ! Compute the left tail amplitude
         tail_amp_left = base_prob / (base_prob - bound_quantile); 
      endif
   endif
   
   ! Default: not close
   do_uniform_tail_right = .false.
   if(bounded_above) then
      ! Compute the CDF at the bounds
      bound_quantile = norm_cdf(upper_bound, tail_mean_right, sd)
      if(abs(base_prob - (1.0_r8 - bound_quantile)) < uniform_threshold) then
         ! If bound and ensemble member are too close, do uniform approximation
         do_uniform_tail_right = .true.
      else
         ! Compute the right tail amplitude
         tail_amp_right = base_prob / (base_prob - (1.0_r8 - bound_quantile))
      endif
   endif

   ! Store the parameters of the tail in the probit data structure
   if(do_uniform_tail_left) then 
      p%params(ens_size + 5) = 1.0_r8
   else
      p%params(ens_size + 5) = 0.0_r8
   endif
   if(do_uniform_tail_right) then 
      p%params(ens_size + 6) = 1.0_r8
   else
      p%params(ens_size + 6) = 0.0_r8
   endif
   p%params(ens_size + 7) = tail_amp_left
   p%params(ens_size + 8) = tail_amp_right
   p%params(ens_size + 9) = tail_mean_left
   p%params(ens_size + 10) = tail_mean_right
   ! Standard deviation of prior tails is prior ensemble standard deviation
   p%params(ens_size + 11) = sd
   p%params(ens_size + 12) = sd
endif

end subroutine to_probit_bounded_normal_rh

!------------------------------------------------------------------------

subroutine to_probit_particle(ens_size, state_ens, p, probit_ens, &
   use_input_p, bounded, bounds)

! Doing a particle filter. Quantiles are (2i-1) / 2n 

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
logical, intent(in)                  :: bounded(2)
real(r8), intent(in)                 :: bounds(2)

integer  :: i, j, indx
integer  :: ens_index(ens_size)
real(r8) :: quantile

! This should fail if any of the input states are not the same as one of the 
! original ensemble states when use_input_p is false. 
if(use_input_p) then
   ! The particles are available from a previous call
   ! The input member gets the same quantile as the corresponding member from the previous call
   ! This can be done vastly more efficiently with either binary searches or by first sorting the
   ! incoming state_ens so that the lower bound for starting the search is updated with each ensemble member
    
   do i = 1, ens_size
      ! Loop through the previous ensemble members
      quantile = -99_r8
      do j = 1, ens_size
         ! Is exact equivalence a problem here?
         if(state_ens(i) == p%params(j)) then
            quantile = 2*(j-1) / (2*ens_size)
            exit
         endif
         ! Test failed to find a match
         if(quantile < 0.0_r8) then
            write(errstring, *) 'Unable to find prior for use_input_p', state_ens(i)
            call error_handler(E_ERR, 'to_probit_particle', errstring, source)
         endif
         ! Do probit transform
         call norm_inv(quantile, probit_ens(i))
      end do
   end do
   
else
   ! Not using a pre-existing distribution
   ! Take care of space for the transform data structure, just need to know sorted prior members
   if(allocated(p%params)) deallocate(p%params)
   allocate(p%params(ens_size))

   ! For particle filter, the required data for inversion is the original ensemble values
   ! Having them in sorted order is useful for subsequent inversion
   call index_sort(state_ens, ens_index, ens_size)
   p%params(1:ens_size) = state_ens(ens_index)

   ! Get the quantiles for each of the ensemble members
   do i = 1, ens_size
      indx = ens_index(i)
      ! The quantiles for a particle filter are just 2(i-1) / 2n
      quantile = 2*(indx - 1) / (2 * ens_size) 

      ! Convert the quantiles to probit space
      call norm_inv(quantile, probit_ens(indx))
   end do 

endif

end subroutine to_probit_particle

!------------------------------------------------------------------------

subroutine convert_all_from_probit(ens_size, num_vars, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: probit_ens(:, :)
type(dist_param_type), intent(inout) :: p(num_vars)
real(r8), intent(out)                :: state_ens(:, :)

! Convert back to the orig
integer  :: i
real(r8) :: temp_ens(ens_size)

do i = 1, num_vars
   call convert_from_probit(ens_size, probit_ens(1:ens_size, i), p(i), temp_ens)
   state_ens(1:ens_size, i) = temp_ens
end do

end subroutine convert_all_from_probit

!------------------------------------------------------------------------

subroutine convert_from_probit(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! Convert back to the orig
if(p%prior_distribution_type == NORMAL_PRIOR) then
   call from_probit_normal(ens_size, probit_ens, p, state_ens)
elseif(p%prior_distribution_type == LOG_NORMAL_PRIOR) then
   call from_probit_log_normal(ens_size, probit_ens, p, state_ens)
elseif(p%prior_distribution_type == UNIFORM_PRIOR) then
   call from_probit_uniform(ens_size, probit_ens, p, state_ens)
elseif(p%prior_distribution_type == GAMMA_PRIOR) then
   call from_probit_gamma(ens_size, probit_ens, p, state_ens)
elseif(p%prior_distribution_type == BETA_PRIOR) then
   call from_probit_beta(ens_size, probit_ens, p, state_ens)
elseif(p%prior_distribution_type == BOUNDED_NORMAL_RH_PRIOR) then
   call from_probit_bounded_normal_rh(ens_size, probit_ens, p, state_ens)
!!!elseif(p%prior_distribution_type == PARTICLE_PRIOR) then
   !!!call from_probit_particle(ens_size, probit_ens, p, state_ens)
else
   write(errstring, *) 'Illegal distribution type', p%prior_distribution_type
   call error_handler(E_ERR, 'convert_from_probit', errstring, source)
   stop
endif


end subroutine convert_from_probit

!------------------------------------------------------------------------

subroutine from_probit_normal(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! Don't do anything for normal
state_ens = probit_ens

end subroutine from_probit_normal


!------------------------------------------------------------------------

subroutine from_probit_log_normal(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! Take the inverse of the log to get back to original space
state_ens = exp(probit_ens)

end subroutine from_probit_log_normal

!------------------------------------------------------------------------

subroutine from_probit_uniform(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

real(r8) :: lower_bound, upper_bound, quantile
integer :: i

! Bounds are the parameters
lower_bound = p%params(1)
upper_bound = p%params(2)

do i = 1, ens_size
   ! First, invert the probit to get a quantile
   quantile = norm_cdf(probit_ens(i), 0.0_r8, 1.0_r8)
   ! Convert from U(0, 1) to U(lower_bound, upper_bound)
   state_ens(i) = lower_bound + quantile * (upper_bound - lower_bound)
end do

! Probably should do an explicit clearing of this storage
! Free the storage
deallocate(p%params)

end subroutine from_probit_uniform

!------------------------------------------------------------------------

subroutine from_probit_gamma(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! Convert back to the orig
real(r8) :: shape, scale, quantile
integer  :: i

! Shape and scale are the distribution parameters
shape = p%params(1)
scale   = p%params(2)

do i = 1, ens_size
   ! First, invert the probit to get a quantile
   quantile = norm_cdf(probit_ens(i), 0.0_r8, 1.0_r8)
   ! Invert the gamma quantiles to get physical space
   state_ens(i) = inv_gamma_cdf(quantile, shape, scale)
end do

! Probably should do an explicit clearing of this storage
! Free the storage
deallocate(p%params)

end subroutine from_probit_gamma

!------------------------------------------------------------------------

subroutine from_probit_beta(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! Convert back to the orig
real(r8) :: alpha, beta, quantile, lower_bound, upper_bound
integer  :: i

! Shape and scale are the distribution parameters
alpha = p%params(1)
beta  = p%params(2)
lower_bound = p%params(3)
upper_bound = p%params(4)

do i = 1, ens_size
   ! First, invert the probit to get a quantile
   quantile = norm_cdf(probit_ens(i), 0.0_r8, 1.0_r8)
   ! Invert the beta quantiles to get scaled physical space
   state_ens(i) = inv_beta_cdf(quantile, alpha, beta)
end do

! Unscale the physical space
state_ens = state_ens * (upper_bound - lower_bound) + lower_bound

! Probably should do an explicit clearing of this storage
! Free the storage
deallocate(p%params)

end subroutine from_probit_beta

!------------------------------------------------------------------------

subroutine from_probit_bounded_normal_rh(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

integer :: i, region
real(r8) :: quantile, target_mass, mass, lower_state, upper_state, lower_q, upper_q
logical  :: bounded_below, bounded_above, do_uniform_tail_left, do_uniform_tail_right
real(r8) :: lower_bound, tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: upper_bound, tail_amp_right, tail_mean_right, tail_sd_right

real(r8) :: fract, lower_mass, upper_mass

! Don't know what to do if original ensemble had all members the same (or nearly so???)
tail_sd_left = p%params(ens_size + 11)
if(tail_sd_left <= 0.0_r8) then
   state_ens = probit_ens
   ! Free the storage; Should do this explicitly?
   deallocate(p%params)
   return
endif

! Get variables out of the parameter storage for clarity
bounded_below = p%params(ens_size + 1) > 0.5_r8
bounded_above = p%params(ens_size + 2) > 0.5_r8
lower_bound = p%params(ens_size + 3)
upper_bound = p%params(ens_size + 4)
do_uniform_tail_left = p%params(ens_size + 5) > 0.5_r8
do_uniform_tail_right = p%params(ens_size + 6) > 0.5_r8
tail_amp_left = p%params(ens_size + 7)
tail_amp_right = p%params(ens_size + 8)
tail_mean_left = p%params(ens_size + 9)
tail_mean_right = p%params(ens_size + 10)
tail_sd_right = p%params(ens_size + 12)

! Convert each probit ensemble member back to physical space
do i = 1, ens_size
   ! First, invert the probit to get a quantile
   ! NOTE: Since we're doing this a ton, may want to have a call specifically for the probit inverse
   quantile = norm_cdf(probit_ens(i), 0.0_r8, 1.0_r8)

   ! Can assume that the quantiles of the original ensemble for the BNRH are uniform
   ! Note that there are some implicit assumptions here about cases where the original 
   ! ensemble had duplicate state members. 
   ! Finding which region this quantile is in is trivial
   region = floor(quantile * (ens_size + 1.0_r8))
   ! Careful about numerical issues moving outside of region [0 ens_size]
   if(region < 0) region = 0
   if(region > ens_size) region = ens_size

   if(region == 0) then
      ! Lower tail
      if(bounded_below .and. do_uniform_tail_left) then
         ! Lower tail uniform
         upper_state = p%params(1)
! NOTE: NEED TO BE CAREFUL OF THE DENOMINATOR HERE AND ON THE PLUS SIDE
         state_ens(i) = lower_bound + &
            (quantile / (1.0_r8 /  (ens_size + 1.0_r8))) * (upper_state - lower_bound)
      !!!elseif(.not. bounded_below) then
      else
         ! Find the mass at the lower bound (which could be unbounded)
         if(bounded_below) then
            lower_mass = tail_amp_left * norm_cdf(lower_bound, tail_mean_left, tail_sd_left)
         else
            lower_mass = 0.0_r8
         endif
         ! Find the mass at the upper bound (ensemble member 1)
         upper_mass = tail_amp_left * norm_cdf(p%params(1), tail_mean_left, tail_sd_left)
         ! What fraction of this mass difference should we go?
         fract = quantile / (1.0_r8 / (ens_size + 1.0_r8))
         target_mass = lower_mass + fract * (upper_mass - lower_mass)
!!!write(*, *) 'first weighted normal call '
         call weighted_norm_inv(tail_amp_left, tail_mean_left, tail_sd_left, target_mass, state_ens(i))
!!!write(*, *) 'back from first weighted normal call '
      endif

   elseif(region == ens_size) then
      ! Upper tail
      if(bounded_above .and. do_uniform_tail_right) then
         ! Upper tail is uniform
         lower_state = p%params(ens_size)
         upper_state = upper_bound
         state_ens(i) = lower_state + & 
            (quantile - (ens_size / (ens_size + 1.0_r8))) * (upper_state - lower_state)
      else
         ! Upper tail is (bounded) normal
         ! Find the mass at the upper bound (which could be unbounded)
         if(bounded_above) then
            upper_mass = tail_amp_right * norm_cdf(upper_bound, tail_mean_right, tail_sd_right)
         else
            upper_mass = 1.0_r8
         endif
         ! Find the mass at the lower bound (ensemble member n)
         lower_mass = tail_amp_right * norm_cdf(p%params(ens_size), tail_mean_right, tail_sd_right)
         ! What fraction of the last interval do we need to move
         fract = (quantile - ens_size / (ens_size + 1.0_r8)) / (1.0_r8 / (ens_size + 1.0_r8))
         target_mass = lower_mass + fract * (upper_mass - lower_mass)
!!!write(*, *) 'first weighted normal call '
         call weighted_norm_inv(tail_amp_right, tail_mean_right, tail_sd_right, target_mass, state_ens(i))
!!!write(*, *) 'back from first weighted normal call '
      endif
         
   else
      ! Interior region; get the quantiles of the region boundary
      lower_q = region / (ens_size + 1.0_r8)
      upper_q = (region + 1.0_r8) / (ens_size + 1.0_r8)
      state_ens(i) = p%params(region) + &
          ((quantile - lower_q) / (upper_q - lower_q)) * (p%params(region + 1) - p%params(region))
   endif
end do

! Check for posterior violating bounds; This may not be needed after development testing
if(bounded_below) then
   do i = 1, ens_size
      if(state_ens(i) < lower_bound) then
         write(errstring, *) 'state_ens ', i, ' less than lower_bound ', state_ens(i)
         call error_handler(E_ERR, 'from_probit_bounded_normal_rh', errstring, source)
      endif
   end do
endif

if(bounded_above) then
   do i = 1, ens_size
      if(state_ens(i) > upper_bound) then
         write(errstring, *) 'state_ens ', i, ' greater than upper_bound ', state_ens(i)
         call error_handler(E_ERR, 'from_probit_bounded_normal_rh', errstring, source)
      endif
   end do
endif

! Probably do this explicitly 
! Free the storage
deallocate(p%params)

end subroutine from_probit_bounded_normal_rh

!------------------------------------------------------------------------

subroutine from_probit_particle(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

integer :: i, indx
real(r8) :: quantile

do i = 1, ens_size
   ! First invert the probit transform to tg
   quantile = norm_cdf(probit_ens(i), 0.0_r8, 1.0_r8)

   ! Invert the quantile for a particle prior
   ! There is a prior ensemble member associated with each 1/ens_size fraction of the quantile 
   ! range
   indx = floor(quantile * ens_size) + 1
   if(indx <= 0) indx = 1
   state_ens(i) = p%params(indx)
end do

! Probably do this explicitly 
! Free the storage
deallocate(p%params)

end subroutine from_probit_particle

!------------------------------------------------------------------------

subroutine ens_quantiles(ens, ens_size, bounded_below, bounded_above, &
                         lower_bound, upper_bound, q)

! Given an ensemble, return information about duplicate values
! in the ensemble. 

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size)
logical,  intent(in)  :: bounded_below, bounded_above
real(r8), intent(in)  :: lower_bound
real(r8), intent(in)  :: upper_bound
real(r8), intent(out) :: q(ens_size)

integer :: i, j, lower_dups, upper_dups, d_start, d_end, series_num
integer :: series_start(ens_size), series_end(ens_size), series_length(ens_size)

! Get number of ensemble members that are duplicates of the lower bound
lower_dups = 0
if(bounded_below) then
   do i = 1, ens_size
      if(ens(i) == lower_bound) then 
         lower_dups = lower_dups + 1
      else
         exit
      endif
   end do
endif

! Get number of ensemble members that are duplicates of the upper bound
upper_dups = 0
if(bounded_above) then
   do i = ens_size, 1, -1
      if(ens(i) == upper_bound) then 
         upper_dups = upper_dups + 1
      else
         exit
      endif
   end do
endif

! If there are duplicate ensemble members away from the boundaries need to revise quantiles
! Make sure not to count duplicates already handled at the boundaries
! Outer loop determines if a series of duplicates starts at sorted index i
d_start = lower_dups + 1 
d_end   = ens_size - upper_dups

! Get start, length, and end of each series of duplicates away from the bounds
series_num = 1
series_start(series_num) = d_start
series_length(series_num) = 1
do i = d_start + 1, d_end
   if(ens(i) == ens(i - 1)) then
      series_length(series_num) = series_length(series_num) + 1
   else
      series_end(series_num) = i-1
      series_num = series_num + 1
      series_start(series_num) = i
      series_length(series_num) = 1
   endif
end do

! Off the end, finish up the last series
series_end(series_num) = d_end

! Now get the value of the quantile for the exact ensemble members
! Start with the lower bound duplicates
do i = 1, lower_dups
   q(i) = lower_dups / (2.0_r8 * (ens_size + 1.0_r8))
end do

! Top bound duplicates next
do i = ens_size - upper_dups + 1, ens_size
   q(i) = upper_dups / (2.0_r8 * (ens_size + 1.0_r8))
end do

! Do the interior series
do i = 1, series_num
   do j = series_start(i), series_end(i)
      q(j) = series_start(i) / (ens_size + 1.0_r8) + (series_length(i) - 1.0_r8) / (2.0_r8 * (ens_size + 1.0_r8))
   end do
end do

end subroutine ens_quantiles

!------------------------------------------------------------------------

end module quantile_distributions_mod
