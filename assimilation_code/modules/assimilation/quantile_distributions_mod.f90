! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! A variety of PDFs, CDFs, quantile functions and other tools for working with distributions
! to implement quantile conserving filters in observation space and regression in quantile space.

module quantile_distributions_mod

use types_mod, only : r8, digits12, PI

use sort_mod,  only : sort, index_sort

use utilities_mod, only : E_ERR, error_handler

implicit none
private

public :: norm_cdf, norm_inv, weighted_norm_inv, convert_to_probit, convert_from_probit, dist_param_type, &
   convert_all_to_probit, convert_all_from_probit


type dist_param_type
   real(r8), allocatable :: params(:)
end type

! Saves the ensemble size used in the previous call of obs_inc_bounded_norm_rhf
integer :: bounded_norm_rhf_ens_size = -99

character(len=512)     :: msgstring
character(len=*), parameter :: source = 'quantile_distributions_mod.f90'

contains

!------------------------------------------------------------------------

subroutine convert_all_to_probit(ens_size, num_vars, state_ens, var_kind, p, probit_ens, use_input_p)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: state_ens(:, :)
type(dist_param_type), intent(inout) :: p(num_vars)
integer, intent(in)                  :: var_kind(num_vars)
real(r8), intent(out)                :: probit_ens(:, :)
logical, intent(in)                  :: use_input_p

! Note that the input and output arrays may have extra copies (first subscript). Passing sections of a
! leading index could be inefficient for time and storage, so avoiding that for now.

integer  :: i

do i = 1, num_vars
   call convert_to_probit(ens_size, state_ens(1:ens_size, i), var_kind(i), p(i), probit_ens(1:ens_size, i), &
      use_input_p)
end do

end subroutine convert_all_to_probit

!------------------------------------------------------------------------

subroutine convert_to_probit(ens_size, state_ens, var_kind, p, probit_ens, use_input_p)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
integer, intent(in)                  :: var_kind
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p

if(var_kind == 0) then
   call to_probit_normal(ens_size, state_ens, p, probit_ens, use_input_p)
elseif(var_kind == 99) then
   ! Need to pass var_kind because different kinds could have different bounds
   call to_probit_bounded_normal_rhf(ens_size, state_ens, var_kind, p, probit_ens, use_input_p)
else
   write(*, *) 'Illegal var_kind in convert_to_probit', var_kind
endif

end subroutine convert_to_probit

!------------------------------------------------------------------------

subroutine to_probit_normal(ens_size, state_ens, p, probit_ens, use_input_p)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p

! Probit transform for nomal. This is just a test since this can be skipped for normals.
real(r8) :: mean, sd

! Initial test is just a bogus thing for normals which require two parameters, mean and sd
if(use_input_p) then
   mean = p%params(1)
   sd   = p%params(2)
else
   mean = sum(state_ens) / ens_size
   sd  = sqrt(sum((state_ens - mean)**2) / (ens_size - 1))
endif

! Do the probit transform for the normal
probit_ens = (state_ens - mean) / sd

! Store these for the inversion
if(.not. use_input_p) then
   if(.not. allocated(p%params)) allocate(p%params(2))
   p%params(1) = mean
   p%params(2) = sd
endif

end subroutine to_probit_normal

!------------------------------------------------------------------------

subroutine to_probit_bounded_normal_rhf(ens_size, state_ens, var_kind, p, probit_ens, use_input_p)

! Note that this is just for transforming back and forth, not for doing the RHF observation update
! This means that we know a prior that the quantiles associated with the initial ensemble are
! uniformly spaced which can be used to simplify converting

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
integer, intent(in)                  :: var_kind
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p

! Probit transform for bounded normal rhf. Need to know the bounds for a given 
integer  :: i, j, indx
integer  :: ens_index(ens_size)
real(r8) :: x, quantile
logical  :: bounded_below, bounded_above, do_uniform_tail_left, do_uniform_tail_right
real(r8) :: lower_bound, tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: upper_bound, tail_amp_right, tail_mean_right, tail_sd_right

! Parameter to control switch to uniform approximation for normal tail
real(r8), parameter :: uniform_threshold = 1e-5_r8

! Save to avoid a modestly expensive computation redundancy
real(r8), save :: dist_for_unit_sd
real(r8) :: mean, sd, base_prob, bound_quantile

if(use_input_p) then
   ! Using an existing ensemble for the RHF points

   ! Get variables out of the parameter storage for clarity
   bounded_below = p%params(1) > 0.5_r8
   bounded_above = p%params(2) > 0.5_r8
   lower_bound = p%params(3)
   upper_bound = p%params(4)
   do_uniform_tail_left = p%params(5) > 0.5_r8
   do_uniform_tail_right = p%params(6) > 0.5_r8
   tail_amp_left = p%params(7)
   tail_amp_right = p%params(8)
   tail_mean_left = p%params(9)
   tail_mean_right = p%params(10)
   tail_sd_left = p%params(11)
   tail_sd_right = p%params(12)

   ! This can be done vastly more efficiently with either binary searches or by first sorting the
   ! incoming state_ens so that the lower bound for starting the search is updated with each ensemble member
   do i = 1, ens_size
      ! Figure out which bin it is in
      x = state_ens(i)
      if(x < p%params(1)) then
         ! In the left tail
         ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
         if(bounded_below .and. x < lower_bound) then
            msgstring = 'Ensemble member less than lower bound first check'
            call error_handler(E_ERR, 'to_probit_bounded_normal_rhf', msgstring, source)
         endif

         if(do_uniform_tail_left) then
            ! Uniform approximation for left tail
            quantile = (x - lower_bound) / (p%params(1) - lower_bound) * (1.0_r8 / (ens_size + 1.0_r8))
         else
            ! It's a normal tail, bounded or not 
            quantile = tail_amp_left * norm_cdf(x, tail_mean_left, tail_sd_left)
         endif

      elseif(x > p%params(ens_size)) then
         ! In the right tail
         ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
         if(bounded_above .and. x > upper_bound) then
            msgstring = 'Ensemble member greater than upper bound first check'
            call error_handler(E_ERR, 'to_probit_bounded_normal_rhf', msgstring, source)
         endif

         if(do_uniform_tail_right) then
            ! Uniform approximation for right tail
            quantile = (ens_size / ens_size + 1.0_r8) + &
               (x - p%params(ens_size)) / (upper_bound - p%params(ens_size)) * (1.0_r8 / (ens_size + 1.0_r8))
         else
            ! It's a normal tail, bounded or not. 
            quantile = tail_amp_right * norm_cdf(x, tail_mean_right, tail_sd_right)
         endif

      else
         ! In an interior bin
         do j = 1, ens_size - 1
            if(x <= p%params(j+1)) then
               quantile = j / (ens_size + 1.0_r8) + &
                  ((x - p%params(j)) / (p%params(j+1) - p%params(j))) * (1.0_r8 / (ens_size + 1.0_r8))
               exit
            endif
         enddo
      endif
      ! Convert to probit space 
      call norm_inv(quantile, probit_ens(i))
   end do
else
   ! No pre-existing distribution, create one
   ! Bounds need to come from somewhere but hard-code here for developmentA
   lower_bound = 0
   upper_bound = 1
   bounded_below = .false.
   bounded_above = .false.

   ! Need to sort. For now, don't worry about efficiency, but may need to somehow pass previous
   ! sorting indexes and use a sort that is faster for nearly sorted data. Profiling can guide the need
   call index_sort(state_ens, ens_index, ens_size)
   do i = 1, ens_size
      quantile = (i * 1.0_r8) / (ens_size + 1.0_r8)
      ! Probit is just the inverse of the standard normal CDF
      call norm_inv(quantile, probit_ens(i))
   end do 

! For RHF, the required data for inversion is the original ensemble values
! Having them in sorted order is useful for subsequent inversion
! It is also useful to store additional information regarding the continuous pdf representation of the tails
! This includes whether the bounds are defined, the values of the bounds, whether a uniform is used in the outer
! bounded bin, the amplitude of the outer continuous normal pdf, the mean of the outer continous
! normal pdf, and the standard deviation of the
! outer continous. 
   ! Do we really want to allow this? Better to deallocate and reallocate?
   if(.not. allocated(p%params)) allocate(p%params(ens_size + 2*6))
   p%params(1:ens_size) = state_ens(ens_index)

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
   if(bounded_norm_rhf_ens_size /= ens_size) then
      call norm_inv(1.0_r8 / (ens_size + 1.0_r8), dist_for_unit_sd)
      ! This will be negative, want it to be a distance so make it positive
      dist_for_unit_sd = -1.0_r8 * dist_for_unit_sd
      ! Keep a record of the ensemble size used to compute dist_for_unit_sd
      bounded_norm_rhf_ens_size = ens_size
   endif
   
   ! Fail if lower bound is larger than smallest ensemble member 
   if(bounded_below) then
      ! Do in two ifs in case the bound is not defined
      if(p%params(1) < lower_bound) then
         msgstring = 'Ensemble member less than lower bound'
         call error_handler(E_ERR, 'to_probit_bounded_normal_rhf', msgstring, source)
      endif
   endif
   
   ! Fail if upper bound is smaller than the largest ensemble member 
   if(bounded_above) then
      if(p%params(ens_size) > upper_bound) then
         msgstring = 'Ensemble member greater than upper bound'
         call error_handler(E_ERR, 'to_probit_bounded_normal_rhf', msgstring, source)
      endif
   endif

   ! Standard deviation of prior tails is prior ensemble standard deviation
   mean = sum(state_ens) / ens_size
   sd  = sqrt(sum((state_ens - mean)**2) / (ens_size - 1))
   ! Find a mean so that 1 / (ens_size + 1) probability is in outer regions
   tail_mean_left = p%params(1) + dist_for_unit_sd * sd
   tail_mean_right = p%params(ens_size) - dist_for_unit_sd * sd

   ! If the distribution is bounded, still want 1 / (ens_size + 1) in outer regions
   ! Put an amplitude term (greater than 1) in front of the tail normals 
   ! Amplitude is 1 if there are no bounds, so start with that
   tail_amp_left  = 1.0_r8
   tail_amp_right = 1.0_r8

   ! DO SOMETHING TO AVOID CASES WHERE THE BOUND AND THE SMALLEST ENSEMBLE ARE VERY CLOSE/SAME
   base_prob = 1.0_r8 / (ens_size + 1.0_r8)
   if(bounded_below) then
      ! Compute the CDF at the bounds
      bound_quantile = norm_cdf(lower_bound, tail_mean_left, sd)
      if(abs(base_prob - bound_quantile) < uniform_threshold) then
         ! If bound and ensemble member are too close, do uniform approximation
         do_uniform_tail_left = .true.
      else
         do_uniform_tail_left = .false.
      endif
   endif
   
   if(bounded_above) then
      ! Compute the CDF at the bounds
      bound_quantile = norm_cdf(upper_bound, tail_mean_right, sd)
      if(abs(base_prob - (1.0_r8 - bound_quantile)) < uniform_threshold) then
         ! If bound and ensemble member are too close, do uniform approximation
         do_uniform_tail_right = .true.
      else
         do_uniform_tail_right = .false.
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
   p%params(ens_size + 11) = sd
   p%params(ens_size + 12) = sd
endif

end subroutine to_probit_bounded_normal_rhf

!------------------------------------------------------------------------

subroutine convert_all_from_probit(ens_size, num_vars, probit_ens, var_kind, p, state_ens)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: probit_ens(:, :)
type(dist_param_type), intent(inout) :: p(num_vars)
integer, intent(in)                  :: var_kind(num_vars)
real(r8), intent(out)                :: state_ens(:, :)

! Convert back to the orig
integer  :: i

do i = 1, num_vars
   call convert_from_probit(ens_size, probit_ens(1:ens_size, i), var_kind(i), p(i), state_ens(1:ens_size, i))
end do

end subroutine convert_all_from_probit

!------------------------------------------------------------------------

subroutine convert_from_probit(ens_size, probit_ens, var_kind, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
integer, intent(in)                  :: var_kind
real(r8), intent(out)                :: state_ens(ens_size)

! Convert back to the orig

if(var_kind == 0) then
   call from_probit_normal(ens_size, probit_ens, p, state_ens)
elseif(var_kind == 99) then
   call from_probit_bounded_normal_rhf(ens_size, probit_ens, var_kind, p, state_ens)
else
   write(*, *) 'Illegal var_kind in convert_from_probit ', var_kind
endif


end subroutine convert_from_probit

!------------------------------------------------------------------------

subroutine from_probit_normal(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! Convert back to the orig
real(r8) :: mean, sd

mean = p%params(1)
sd   = p%params(2)
state_ens = probit_ens * sd + mean

! Free the storage
deallocate(p%params)

end subroutine from_probit_normal

!------------------------------------------------------------------------

subroutine from_probit_bounded_normal_rhf(ens_size, probit_ens, var_kind, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
integer, intent(in)                  :: var_kind
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

integer :: i, region
real(r8) :: quantile, target_mass, delta_q, mass, lower_state, upper_state, lower_q, upper_q
logical  :: bounded_below, bounded_above, do_uniform_tail_left, do_uniform_tail_right
real(r8) :: lower_bound, tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: upper_bound, tail_amp_right, tail_mean_right, tail_sd_right

! Get variables out of the parameter storage for clarity
bounded_below = p%params(1) > 0.5_r8
bounded_above = p%params(2) > 0.5_r8
lower_bound = p%params(3)
upper_bound = p%params(4)
do_uniform_tail_left = p%params(5) > 0.5_r8
do_uniform_tail_right = p%params(6) > 0.5_r8
tail_amp_left = p%params(7)
tail_amp_right = p%params(8)
tail_mean_left = p%params(9)
tail_mean_right = p%params(10)
tail_sd_left = p%params(11)
tail_sd_right = p%params(12)

! Convert each probit ensemble member back to physical space
do i = 1, ens_size
   ! First, invert the probit to get a quantile
   ! NOTE: Since we're doing this a ton, may want to have a call specifically for the probit inverse
   quantile = norm_cdf(probit_ens(i), 0.0_r8, 1.0_r8)

   ! Can assume that the quantiles of the original ensemble for the RHF are uniform
   ! Finding which region this quantile is in is trivial
   region = floor(quantile * (ens_size + 1.0_r8))
   ! Careful about numerical issues moving outside of region [0 ens_size]
   if(region > ens_size) region = ens_size

   if(region == 0) then
      ! Lower tail
      if(do_uniform_tail_left) then
         ! Lower tail uniform
         lower_state = lower_bound
         upper_state = p%params(1)
         state_ens(i) = lower_state + &
            (quantile / (ens_size + 1.0_r8)) * (upper_state - lower_state)
      else
         ! Lower tail is (bounded) normal
         ! What is the value of the weighted normal at the smallest ensemble member
         mass = tail_amp_left * norm_cdf(p%params(1), tail_mean_left, tail_sd_left)
         delta_q = 1.0 / (ens_size + 1.0_r8) - quantile
         target_mass = mass - delta_q
         call weighted_norm_inv(tail_amp_left, tail_mean_left, tail_sd_left, target_mass, state_ens(i))
      endif

   elseif(region == ens_size) then
      ! Upper tail
      if(do_uniform_tail_right) then
         ! Upper tail is uniform
         lower_state = p%params(ens_size)
         upper_state = upper_bound
         state_ens(i) = lower_state + & 
            (quantile - (ens_size / (ens_size + 1.0_r8))) * (upper_state - lower_state)
      else
         ! Upper tail is (bounded) normal
         ! Value of weighted normal at largest ensemble member
         mass = tail_amp_right * norm_cdf(p%params(ens_size), tail_mean_right, tail_sd_right)
         delta_q = quantile - ens_size / (ens_size + 1.0_r8)
         target_mass = mass + delta_q
         call weighted_norm_inv(tail_amp_right, tail_mean_right, tail_sd_right, target_mass, state_ens(i))
      endif
         
   else
      ! Interior region; get the quantiles of the region boundary
      lower_q = region / (ens_size + 1.0_r8)
      upper_q = (region + 1.0_r8) / (ens_size + 1.0_r8)
      state_ens(i) = p%params(region) + &
          ((quantile - lower_q) / (upper_q - lower_q)) * (p%params(region + 1) - p%params(region))
   endif

end do

! Free the storage
deallocate(p%params)

end subroutine from_probit_bounded_normal_rhf

!------------------------------------------------------------------------

function norm_cdf(x_in, mean, sd)

! Approximate cumulative distribution function for normal
! with mean and sd evaluated at point x_in
! Only works for x>= 0.

real(r8)             :: norm_cdf
real(r8), intent(in) :: x_in, mean, sd

real(digits12) :: x, p, b1, b2, b3, b4, b5, t, density, nx

! Convert to a standard normal
nx = (x_in - mean) / sd

x = abs(nx)


! Use formula from Abramowitz and Stegun to approximate
p = 0.2316419_digits12
b1 = 0.319381530_digits12
b2 = -0.356563782_digits12
b3 = 1.781477937_digits12
b4 = -1.821255978_digits12
b5 = 1.330274429_digits12

t = 1.0_digits12 / (1.0_digits12 + p * x)

density = (1.0_digits12 / sqrt(2.0_digits12 * PI)) * exp(-x*x / 2.0_digits12)

norm_cdf = 1.0_digits12 - density * &
   ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t

if(nx < 0.0_digits12) norm_cdf = 1.0_digits12 - norm_cdf

!write(*, *) 'cdf is ', norm_cdf

end function norm_cdf

!------------------------------------------------------------------------

subroutine weighted_norm_inv(alpha, mean, sd, p, x)

! Find the value of x for which the cdf of a N(mean, sd) multiplied times
! alpha has value p.

real(r8), intent(in)  :: alpha, mean, sd, p
real(r8), intent(out) :: x

real(r8) :: np

! Can search in a standard normal, then multiply by sd at end and add mean
! Divide p by alpha to get the right place for weighted normal
np = p / alpha

! Find spot in standard normal
call norm_inv(np, x)

! Add in the mean and normalize by sd
x = mean + x * sd

end subroutine weighted_norm_inv


!------------------------------------------------------------------------

subroutine norm_inv(p, x)

real(r8), intent(in)  :: p
real(r8), intent(out) :: x

! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero

real(r8) :: p_low,p_high
real(r8) :: a1,a2,a3,a4,a5,a6
real(r8) :: b1,b2,b3,b4,b5
real(r8) :: c1,c2,c3,c4,c5,c6
real(r8) :: d1,d2,d3,d4
real(r8) :: q,r
a1 = -39.69683028665376_digits12
a2 =  220.9460984245205_digits12
a3 = -275.9285104469687_digits12
a4 =  138.357751867269_digits12
a5 = -30.66479806614716_digits12
a6 =  2.506628277459239_digits12
b1 = -54.4760987982241_digits12
b2 =  161.5858368580409_digits12
b3 = -155.6989798598866_digits12
b4 =  66.80131188771972_digits12
b5 = -13.28068155288572_digits12
c1 = -0.007784894002430293_digits12
c2 = -0.3223964580411365_digits12
c3 = -2.400758277161838_digits12
c4 = -2.549732539343734_digits12
c5 =  4.374664141464968_digits12
c6 =  2.938163982698783_digits12
d1 =  0.007784695709041462_digits12
d2 =  0.3224671290700398_digits12
d3 =  2.445134137142996_digits12
d4 =  3.754408661907416_digits12
p_low  = 0.02425_digits12
p_high = 1_digits12 - p_low
! Split into an inner and two outer regions which have separate fits
if(p < p_low) then
   q = sqrt(-2.0_digits12 * log(p))
   x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else if(p > p_high) then
   q = sqrt(-2.0_digits12 * log(1.0_digits12 - p))
   x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else
   q = p - 0.5_digits12
   r = q*q
   x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / &
      (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_digits12)
endif

end subroutine norm_inv

!------------------------------------------------------------------------


end module quantile_distributions_mod
