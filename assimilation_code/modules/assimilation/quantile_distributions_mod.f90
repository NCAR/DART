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

integer, parameter :: NORMAL_PRIOR = 1
integer, parameter :: BOUNDED_NORMAL_RH_PRIOR = 2

public :: norm_cdf, norm_inv, weighted_norm_inv, convert_to_probit, &
   convert_from_probit, dist_param_type, convert_all_to_probit, &
   convert_all_from_probit, probit_dist_info, &
   NORMAL_PRIOR, BOUNDED_NORMAL_RH_PRIOR

type dist_param_type
   integer               :: prior_distribution_type
   real(r8), allocatable :: params(:)
end type

! Saves the ensemble size used in the previous call of obs_inc_bounded_norm_rh
integer :: bounded_norm_rh_ens_size = -99

character(len=512)     :: msgstring
character(len=*), parameter :: source = 'quantile_distributions_mod.f90'

contains

!------------------------------------------------------------------------

subroutine probit_dist_info(kind, is_state, is_inflation, dist_type, &
   bounded, bounds)

! Computes the details of the probit transform for initial experiments
! with Molly 

integer,  intent(in)  :: kind
logical,  intent(in)  :: is_state      ! True for state variable, false for obs
logical,  intent(in)  :: is_inflation  ! True for inflation transform
integer,  intent(out) :: dist_type
logical,  intent(out) :: bounded(2)
real(r8), intent(out) :: bounds(2)

! Have input information about the kind of the state or observation being transformed
! along with additional logical info that indicates whether this is an observation
! or state variable and about whether the transformation is being done for inflation
! or for regress. 
! Need to select the appropriate transform. At present, options are NORMAL_PRIOR
! which does nothing or BOUNDED_NORMAL_RH_PRIOR. 
! If the BNRH is selected then information about the bounds must also be set.
! The two dimensional logical array 'bounded' is set to false for no bounds and true
! for bounded. the first element of the array is for the lower bound, the second for the upper.
! If bounded is chosen, the corresponding bound value(s) must be set in the two dimensional 
! real array 'bounds'.
! For example, if my_state_kind corresponds to a sea ice fraction then an appropriate choice
! would be:
! bounded(1) = .true.;  bounded(2) = .true.
! bounds(1)  = 0.0_r8;  bounds(2)  = 1.0_r8

! This logic is consistent with Quantile Regression paper square experiments
if(is_inflation) then 
   dist_type = BOUNDED_NORMAL_RH_PRIOR
   bounded = .false.
elseif(is_state) then
   dist_type = BOUNDED_NORMAL_RH_PRIOR
   bounded = .false.
else
   dist_type = BOUNDED_NORMAL_RH_PRIOR
   bounded(1) = .true.;    bounded(2) = .false.
   bounds(1)  = 0.0_r8
endif

end subroutine probit_dist_info

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
elseif(p%prior_distribution_type == BOUNDED_NORMAL_RH_PRIOR) then
   call to_probit_bounded_normal_rh(ens_size, state_ens, p, probit_ens, &
      use_input_p, bounded, bounds)
else
   write(*, *) 'Illegal distribution in convert_to_probit', p%prior_distribution_type
   stop
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

! Don't need to do anything for normal, but keep code below to show what it could look like
probit_ens = state_ens
return

! Get parameters
if(use_input_p) then
   mean = p%params(1)
   sd   = p%params(2)
else
   mean = sum(state_ens) / ens_size
   sd  = sqrt(sum((state_ens - mean)**2) / (ens_size - 1))
   if(.not. allocated(p%params)) allocate(p%params(2))
   p%params(1) = mean
   p%params(2) = sd
endif

! Do the probit transform for the normal
probit_ens = (state_ens - mean) / sd

end subroutine to_probit_normal

!------------------------------------------------------------------------

subroutine to_probit_bounded_normal_rh(ens_size, state_ens, p, probit_ens, &
   use_input_p, bounded, bounds)

! Note that this is just for transforming back and forth, not for doing the RHF observation update
! This means that we know a prior that the quantiles associated with the initial ensemble are
! uniformly spaced which can be used to simplify converting

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
logical, intent(in)                  :: bounded(2)
real(r8), intent(in)                 :: bounds(2)

! Probit transform for bounded normal rh.
integer  :: i, j, indx
integer  :: ens_index(ens_size)
real(r8) :: x, quantile
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
   tail_sd_left = p%params(ens_size + 11)
   tail_sd_right = p%params(ens_size + 12)

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
            call error_handler(E_ERR, 'to_probit_bounded_normal_rh', msgstring, source)
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
            call error_handler(E_ERR, 'to_probit_bounded_normal_rh', msgstring, source)
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
               quantile = (j * 1.0_r8) / (ens_size + 1.0_r8) + &
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
   lower_bound = bounds(1)
   upper_bound = bounds(2)
   bounded_below = bounded(1)
   bounded_above = bounded(2)

   ! Need to sort. For now, don't worry about efficiency, but may need to somehow pass previous
   ! sorting indexes and use a sort that is faster for nearly sorted data. Profiling can guide the need
   call index_sort(state_ens, ens_index, ens_size)
   do i = 1, ens_size
      indx = ens_index(i)
      quantile = (i * 1.0_r8) / (ens_size + 1.0_r8)
      ! Probit is just the inverse of the standard normal CDF
      call norm_inv(quantile, probit_ens(indx))
   end do 

   ! For BNRH, the required data for inversion is the original ensemble values
   ! Having them in sorted order is useful for subsequent inversion
   ! It is also useful to store additional information regarding the continuous pdf representation of the tails
   ! This includes whether the bounds are defined, the values of the bounds, whether a uniform is used in the outer
   ! bounded bin, the amplitude of the outer continuous normal pdf, the mean of the outer continous
   ! normal pdf, and the standard deviation of the
   ! outer continous. 

   if(allocated(p%params)) deallocate(p%params)
   allocate(p%params(ens_size + 2*6))
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
         msgstring = 'Ensemble member less than lower bound'
         call error_handler(E_ERR, 'to_probit_bounded_normal_rh', msgstring, source)
      endif
   endif
   
   ! Fail if upper bound is smaller than the largest ensemble member 
   if(bounded_above) then
      if(p%params(ens_size) > upper_bound) then
         msgstring = 'Ensemble member greater than upper bound'
         call error_handler(E_ERR, 'to_probit_bounded_normal_rh', msgstring, source)
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
   p%params(ens_size + 11) = sd
   p%params(ens_size + 12) = sd
endif

end subroutine to_probit_bounded_normal_rh

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
elseif(p%prior_distribution_type == BOUNDED_NORMAL_RH_PRIOR) then
   call from_probit_bounded_normal_rh(ens_size, probit_ens, p, state_ens)
else
   write(*, *) 'Illegal distribution in convert_from_probit ', p%prior_distribution_type
   stop
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

! Don't do anything for normal
state_ens = probit_ens
return

mean = p%params(1)
sd   = p%params(2)
state_ens = probit_ens * sd + mean

! Probably should do an explicit clearing of this storage
! Free the storage
deallocate(p%params)

end subroutine from_probit_normal

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

real(r8) :: bound_inv, correction

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
tail_sd_left = p%params(ens_size + 11)
tail_sd_right = p%params(ens_size + 12)

! Convert each probit ensemble member back to physical space
do i = 1, ens_size
   ! First, invert the probit to get a quantile
   ! NOTE: Since we're doing this a ton, may want to have a call specifically for the probit inverse
   quantile = norm_cdf(probit_ens(i), 0.0_r8, 1.0_r8)

   ! Can assume that the quantiles of the original ensemble for the BNRH are uniform
   ! Finding which region this quantile is in is trivial
   region = floor(quantile * (ens_size + 1.0_r8))
   ! Careful about numerical issues moving outside of region [0 ens_size]
   if(region < 0) region = 0
   ! This behavior has been documented 
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
         ! Lower tail is (bounded) normal, work in from the bottom
         ! Value of weighted normal at smallest member
         mass = tail_amp_left * norm_cdf(p%params(1), tail_mean_left, tail_sd_left)
         target_mass = mass - (1.0_r8 / (ens_size + 1.0_r8) - quantile)
         call weighted_norm_inv(tail_amp_left, tail_mean_left, tail_sd_left, target_mass, state_ens(i))

!------------------- The following block can prevent any risk of getting below bounds
!                    results, but is expensive and may be unneeded with thresholds in general
!                    Code is left here, along with elseif 8 lines above in case this becomes an issue.
!                    A similar block would be needed for the upper bounds. Note that there is also
!                    a risk of destroying the sorted order by doing this and that might require further
!                    subtlety
!      elseif(bounded_below .and. .not. do_uniform_tail_left) then
!         ! Work in from the edge??? Have to watch for sorting problems???
!         ! Find mass at the lower bound
!         mass = tail_amp_left * norm_cdf(lower_bound, tail_mean_left, tail_sd_left)
!         ! If the inverse for the boundary gives something less than the bound have to fix it
!         call weighted_norm_inv(tail_amp_left, tail_mean_left, tail_sd_left, mass, bound_inv)
!         correction = abs(min(0.0_r8, bound_inv))
!         target_mass = mass + quantile
!         call weighted_norm_inv(tail_amp_left, tail_mean_left, tail_sd_left, target_mass, state_ens(i))
!         state_ens(i) = state_ens(i) + correction
!------------------- End unused block -------------------------------
       
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
         ! Value of weighted normal at largest ensemble member
         mass = tail_amp_right * norm_cdf(p%params(ens_size), tail_mean_right, tail_sd_right)
         ! How much mass we need past the last ensemble member which has N/(N+1) quantile
         target_mass = mass + quantile - (ens_size / (ens_size + 1.0_r8))
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

! Probably do this explicitly 
! Free the storage
deallocate(p%params)

end subroutine from_probit_bounded_normal_rh

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

if(nx < 0.0_digits12) then
   norm_cdf = 0.5_digits12 * erfc(-nx / sqrt(2.0_digits12))
else
   norm_cdf = 0.5_digits12 * (1.0_digits12 + erf(nx / sqrt(2.0_digits12)))
endif
return

! Old version left for now
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

subroutine norm_inv(p_in, x)

real(r8), intent(in)  :: p_in
real(r8), intent(out) :: x

! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero

real(r8) :: p
real(r8) :: p_low,p_high
real(r8) :: a1,a2,a3,a4,a5,a6
real(r8) :: b1,b2,b3,b4,b5
real(r8) :: c1,c2,c3,c4,c5,c6
real(r8) :: d1,d2,d3,d4
real(r8) :: q,r

! Truncate out of range quantiles, converts them to smallest positive number or largest number <1
! This solution is stable, but may lead to underflows being thrown. May want to 
! think of a better solution. 
p = p_in
if(p <= 0.0_r8) p = tiny(p_in)
if(p >= 1.0_r8) p = nearest(1.0_r8, -1.0_r8)

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
