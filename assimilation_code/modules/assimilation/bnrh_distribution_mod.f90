! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module bnrh_distribution_mod

use types_mod,      only : r8, missing_r8

use utilities_mod,  only : E_ERR, E_MSG, error_handler

use sort_mod, only : index_sort

use normal_distribution_mod, only : normal_cdf, inv_std_normal_cdf, inv_weighted_normal_cdf, &
                                    normal_mean_sd

use distribution_params_mod, only : distribution_params_type

implicit none
private

public :: bnrh_cdf, bnrh_cdf_params, bnrh_cdf_initialized_vector, &
          inv_bnrh_cdf, inv_bnrh_cdf_params, get_bnrh_sd, inv_bnrh_cdf_like

character(len=512)          :: errstring
character(len=*), parameter :: source = 'bnrh_distribution_mod.f90'

! Saves the ensemble size used in the previous call of bnrh_cdf
integer :: saved_ens_size = -99
! Cached value of dist_for_unit_sd for this saved_ens_size
real(r8), save :: dist_for_unit_sd

! Parameter to control switch to uniform approximation for normal tail
! This defines how many quantiles the bound is from the outermost ensemble member
! If closer than this, can get into precision error problems with F(F-1(X)) on tails
! Can also get other precision problems with the amplitude for the normal in the bounded region
real(r8), parameter :: uniform_threshold = 0.01_r8

contains

!-----------------------------------------------------------------------
subroutine bnrh_cdf_params(x, ens_size, bounded_below, bounded_above, &
   lower_bound, upper_bound, p, quantiles)

integer,                        intent(in)    :: ens_size
real(r8),                       intent(in)    :: x(ens_size)
logical,                        intent(in)    :: bounded_below, bounded_above
real(r8),                       intent(in)    :: lower_bound, upper_bound
type(distribution_params_type), intent(inout) :: p
real(r8),                       intent(out)   :: quantiles(ens_size)

real(r8) :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: tail_amp_right, tail_mean_right, tail_sd_right
real(r8) :: sort_ens(ens_size)
logical  :: do_uniform_tail_left, do_uniform_tail_right

call bnrh_cdf(x, ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
   sort_ens, quantiles, &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right)

! Store the info about this cdf in the distribution_params_type
call pack_bnrh_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
   do_uniform_tail_left, do_uniform_tail_right, tail_amp_left, tail_amp_right, &
   tail_mean_left, tail_mean_right, tail_sd_left, tail_sd_right, sort_ens, p)

end subroutine bnrh_cdf_params

!-----------------------------------------------------------------------

subroutine bnrh_cdf(x, ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
   sort_x, quantiles,  &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right)

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: x(ens_size)
logical,  intent(in)  :: bounded_below, bounded_above
real(r8), intent(in)  :: lower_bound,   upper_bound
real(r8), intent(out) :: sort_x(ens_size)
real(r8), intent(out) :: quantiles(ens_size)
real(r8), intent(out) :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8), intent(out) :: tail_amp_right, tail_mean_right, tail_sd_right
logical,  intent(out) :: do_uniform_tail_left, do_uniform_tail_right

real(r8) :: q(ens_size)
real(r8) :: del_q, mean, bound_quantile
integer :: sort_index(ens_size), indx, i

! Computes all information about a rank histogram cdf given the ensemble and bounds

! Get ensemble mean and sd
call normal_mean_sd(x, ens_size, mean, tail_sd_left)
tail_sd_right = tail_sd_left

! Don't know what to do if sd is 0; tail_sds returned are illegal value to indicate this
if(tail_sd_left <= 0.0_r8) then
   tail_sd_left = -99_r8
   tail_sd_right = -99_r8 
   return
endif

! Sort. For now, don't worry about efficiency, but may need to somehow pass previous
! sorting indexes and use a sort that is faster for nearly sorted data. Profiling can guide the need
call index_sort(x, sort_index, ens_size)
sort_x = x(sort_index)

! Fail if lower bound is larger than smallest ensemble member 
if(bounded_below) then
   ! Do in two ifs in case the bound is not defined
   if(sort_x(1) < lower_bound) then
      write(errstring, *) 'Smallest ensemble member less than lower bound', &
         sort_x(1), lower_bound
      call error_handler(E_ERR, 'bnrh_cdf', errstring, source)
   endif
endif
  
! Fail if upper bound is smaller than the largest ensemble member 
if(bounded_above) then
   if(sort_x(ens_size) > upper_bound) then
      write(errstring, *) 'Largest ensemble member greater than upper bound', &
         sort_x(ens_size), upper_bound
      call error_handler(E_ERR, 'bnrh_cdf', errstring, source)
   endif
endif

! The ensemble size array q contains the sorted quantiles corresponding to the sorted ensemble sort_x
call ens_quantiles(sort_x, ens_size, &
   bounded_below, bounded_above, lower_bound, upper_bound, q)
! The quantiles array has the unsorted quantiles corresponding to the unsorted input ensemble, x
do i = 1, ens_size
   indx = sort_index(i)
   quantiles(indx) = q(i)
end do

! Compute the characteristics of tails

! For unit normal, find distance from mean to where cdf is 1/(ens_size+1) (del_q_.
! Saved to avoid redundant computation for repeated calls with same ensemble size
del_q = 1.0_r8 / (ens_size + 1.0_r8)

if(saved_ens_size /= ens_size) then
   dist_for_unit_sd =  inv_std_normal_cdf(del_q)
   ! This will be negative, want it to be a distance so make it positive
   dist_for_unit_sd = -1.0_r8 * dist_for_unit_sd
   ! Keep a record of the ensemble size used to compute dist_for_unit_sd
   saved_ens_size = ens_size
endif

! Find a mean so that 1 / (ens_size + 1) probability is in outer regions
tail_mean_left =  sort_x(1)        + dist_for_unit_sd * tail_sd_left
tail_mean_right = sort_x(ens_size) - dist_for_unit_sd * tail_sd_right

! If the distribution is bounded, still want 1 / (ens_size + 1) (del_q) in outer regions
! Put an amplitude term (greater than 1) in front of the tail normals 
! Amplitude is 1 if there are no bounds, so start with that
tail_amp_left  = 1.0_r8
tail_amp_right = 1.0_r8

! Switch to uniform for cases where bound and outermost ensemble have close quantiles
! Default: not close
do_uniform_tail_left = .false.
if(bounded_below) then
   ! Compute the CDF at the bounds
   bound_quantile = normal_cdf(lower_bound, tail_mean_left, tail_sd_left)
   ! Note that due to roundoff it is possible for del_q - quantile to be slightly negative
   if((del_q - bound_quantile) / del_q < uniform_threshold) then
      ! If bound and ensemble member are too close, do uniform approximation
      do_uniform_tail_left = .true.
   else
      ! Compute the left tail amplitude
      tail_amp_left = del_q / (del_q - bound_quantile);
   endif
endif

! Default: not close
do_uniform_tail_right = .false.
if(bounded_above) then
   ! Compute the CDF at the bounds
   bound_quantile = normal_cdf(upper_bound, tail_mean_right, tail_sd_right)
   ! Note that due to roundoff it is possible for the numerator to be slightly negative
   if((bound_quantile - (1.0_r8 - del_q)) / del_q < uniform_threshold) then
      ! If bound and ensemble member are too close, do uniform approximation
      do_uniform_tail_right = .true.
   else
      ! Compute the right tail amplitude
      tail_amp_right = del_q / (del_q - (1.0_r8 - bound_quantile))
   endif
endif

end subroutine bnrh_cdf

!-----------------------------------------------------------------------

subroutine bnrh_cdf_initialized_vector(x, num, p, quantiles)

integer, intent(in)                        :: num
real(r8), intent(in)                       :: x(num)
type(distribution_params_type), intent(in) :: p
real(r8), intent(out)                      :: quantiles(num)

real(r8) :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: tail_amp_right, tail_mean_right, tail_sd_right
logical  :: do_uniform_tail_left, do_uniform_tail_right

! Given the sorted ensemble (sort_ens) that defines a bnrh CDF and all the corresponding
! information about that distribution, computes the value of the CDF for a vector of num
! elsements (x) and returns those quantiles.

! In the default filter usage, this is only used for doing the probit transform for the
! posterior observation ensemble. In this case, the size of vector x is the same as the
! ensemble size. For hybrid filter applications, the ensemble size defining the BNRH CDF
! might be different from the ensemble that needs updating, so x could have a different size.

real(r8) :: q(p%ens_size)
integer  :: i

! Extract the required information from the distribution_params_type
call unpack_bnrh_params(p, do_uniform_tail_left, do_uniform_tail_right, &
   tail_amp_left, tail_amp_right, tail_mean_left, tail_mean_right, tail_sd_left, tail_sd_right)

! Compute the quantiles of each of the sorted ensemble members that define the BNRH distribution.
! This was all computed when the distribution was originally set up, could choose to cache that
! in the params structure for efficiency. This is only used for the single observation posterior
! so one could only save in that case removing any storage concerns.
call ens_quantiles(p%ens, p%ens_size, p%bounded_below, p%bounded_above, p%lower_bound, p%upper_bound, q)

! Loop through the values in the x vector to compute the CDF at each one.
! This can be done vastly more efficiently with either binary searches or by first sorting the
! vector of values (x) for which the CDF needs to be computed
do i = 1, p%ens_size
   ! Figure out which bin it is in
   call bnrh_cdf_initialized(x(i), p%ens_size, p%ens, p%bounded_below, p%bounded_above, p%lower_bound, p%upper_bound, &
      tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
      tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, q, quantiles(i))
end do

end subroutine bnrh_cdf_initialized_vector

!-----------------------------------------------------------------------

subroutine bnrh_cdf_initialized(x, ens_size, sort_ens, bounded_below, bounded_above,   &
   lower_bound, upper_bound,                                             &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right,&
   q, quantile)

real(r8), intent(in)  :: x
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: sort_ens(ens_size)
logical,  intent(in)  :: bounded_below, bounded_above
real(r8), intent(in)  :: lower_bound, upper_bound
real(r8), intent(in)  :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8), intent(in)  :: tail_amp_right, tail_mean_right, tail_sd_right
logical,  intent(in)  :: do_uniform_tail_left, do_uniform_tail_right
real(r8), intent(in)  :: q(ens_size)
real(r8), intent(out) :: quantile

real(r8) :: upper_q, fract, del_q, q_at_largest_ens
integer :: j

! Quantile increment between ensemble members for bnrh
del_q = 1.0_r8 / (ens_size + 1.0_r8)

if(x < sort_ens(1)) then
   ! In the left tail
   ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
   if(bounded_below .and. x < lower_bound) then
      write(errstring, *) 'Ensemble member less than lower bound', x, lower_bound
      call error_handler(E_ERR, 'bnrh_cdf_initialized', errstring, source)
      ! This error can occur due to roundoff in increment generation from BNRHF
      ! See discussion in function fix_bounds.
   endif

   if(do_uniform_tail_left) then
      ! Uniform approximation for left tail; Note that denominator cannot be 0 but could be small
      quantile = (x - lower_bound) / (sort_ens(1) - lower_bound) * del_q
   else
      ! It's a normal tail
      if(bounded_below) then
         quantile = tail_amp_left * (normal_cdf(x,           tail_mean_left, tail_sd_left) - &
                                     normal_cdf(lower_bound, tail_mean_left, tail_sd_left))
      else        ! Unbounded, tail normal goes all the way down to quantile 0, amplitude is 1
         quantile = (normal_cdf(x,          tail_mean_left, tail_sd_left) / &
                    normal_cdf(sort_ens(1), tail_mean_left, tail_sd_left)) * del_q
      endif
      ! Make sure it doesn't sneak past the quantile of the smallest ensemble member due to round-off
      quantile = min(quantile, q(1))
   endif
elseif(x == sort_ens(1)) then
   ! This takes care of cases where there are multiple bnrh values at the bdry or at first ensemble
   quantile = q(1)
elseif(x > sort_ens(ens_size)) then
   ! In the right tail
   ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
   if(bounded_above .and. x > upper_bound) then
      write(errstring, *) 'Ensemble member greater than upper bound first check(see code)', x, upper_bound
      call error_handler(E_ERR, 'bnrh_cdf_initialized', errstring, source)
      ! This error can occur due to roundoff in increment generation from bounded BNRHF
      ! See discussion in function fix_bounds
   endif

   if(do_uniform_tail_right) then
      ! Uniform approximation for right tail
      ! The division here could be a concern. However, if sort_ens(ens_size) == upper_bound, then
      ! x cannot be > sort_ens(ens_size).
      quantile = ens_size *del_q + &
         (x - sort_ens(ens_size)) / (upper_bound - sort_ens(ens_size)) * del_q
   else
      ! It's a normal tail
      q_at_largest_ens = normal_cdf(sort_ens(ens_size), tail_mean_right, tail_sd_right)
      ! Want to avoid quantiles exceeding 1 due to numerical issues. Do fraction of the normal part
      if(bounded_above) then
         upper_q = tail_amp_right * normal_cdf(upper_bound, tail_mean_right, tail_sd_right)
         fract = (tail_amp_right * normal_cdf(x, tail_mean_right, tail_sd_right) - &
                  tail_amp_right * q_at_largest_ens) / (upper_q - tail_amp_right * q_at_largest_ens)
      else
         ! Normal goes all the way to infinity, amplitude is 1, q at infinity is 1
         fract = (normal_cdf(x, tail_mean_right, tail_sd_right) - q_at_largest_ens) / &
            (1.0_r8 -  q_at_largest_ens)
      endif

      quantile = ens_size * del_q + fract * del_q
      quantile = min(quantile, 1.0_r8)
   endif

else
   ! In an interior bin
   do j = 1, ens_size - 1
      if(x < sort_ens(j+1)) then
         ! The division here could be a concern. 
         ! However, sort_ens(j)< x < sort_ens(j+1) so the two cannot be equal
         quantile = j * del_q + &
            ((x - sort_ens(j)) / (sort_ens(j+1) - sort_ens(j))) * del_q
         exit
      elseif(x == sort_ens(j+1)) then
         quantile = q(j+1)
         exit
      endif
   enddo
endif

end subroutine bnrh_cdf_initialized

!-----------------------------------------------------------------------
subroutine inv_bnrh_cdf_params(quantiles, ens_size, p, x)

integer,                        intent(in)    :: ens_size
real(r8),                       intent(in)    :: quantiles(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out) :: x(ens_size)

real(r8) :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: tail_amp_right, tail_mean_right, tail_sd_right
logical  :: do_uniform_tail_left, do_uniform_tail_right

call unpack_bnrh_params(p, do_uniform_tail_left, do_uniform_tail_right, &
   tail_amp_left, tail_amp_right, tail_mean_left, tail_mean_right, tail_sd_left, tail_sd_right)

call inv_bnrh_cdf(quantiles, ens_size, p%ens,                            &
   p%bounded_below, p%bounded_above, p%lower_bound, p%upper_bound,       &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, x)

end subroutine inv_bnrh_cdf_params

!-----------------------------------------------------------------------

subroutine inv_bnrh_cdf(quantiles, ens_size, sort_ens, &
   bounded_below, bounded_above, lower_bound, upper_bound, &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left,  &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, x)

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: quantiles(ens_size)
real(r8), intent(in)  :: sort_ens(ens_size)
logical,  intent(in)  :: bounded_below, bounded_above
real(r8), intent(in)  :: lower_bound, upper_bound
real(r8), intent(in)  :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8), intent(in)  :: tail_amp_right, tail_mean_right, tail_sd_right
logical,  intent(in)  :: do_uniform_tail_left, do_uniform_tail_right
real(r8), intent(out) :: x(ens_size)

integer :: region, i, j
real(r8) :: lower_state, upper_state, lower_mass, upper_mass, target_mass
real(r8) :: q(ens_size), curr_q, lower_q, upper_q, del_q, fract

! Quantile increment between ensemble members for bnrh
del_q = 1.0_r8 / (ens_size + 1.0_r8)

do i = 1, ens_size
   q(i) = i * del_q
end do

! Loop through each ensemble member to find posterior state
do i = 1, ens_size
   curr_q = quantiles(i)
   ! Which region is this quantile in?
   ! BNRH quantiles are uniform; finding region for this quantile is trivial
   region = floor(curr_q * (ens_size + 1.0_r8))
   ! Careful about numerical issues moving outside of region [0 ens_size]
   if(region < 0) region = 0
   if(region > ens_size) region = ens_size

   if(region == 0) then
      ! Lower tail
      if(bounded_below .and. do_uniform_tail_left) then
         ! Lower tail uniform
         upper_state = sort_ens(1)
         x(i) = lower_bound + (curr_q / q(1)) * (upper_state - lower_bound)
      else
         ! Find the mass at the lower bound (which could be unbounded)
         if(bounded_below) then
            lower_mass = tail_amp_left * &
               normal_cdf(lower_bound, tail_mean_left, tail_sd_left)
         else
            lower_mass = 0.0_r8
         endif
         ! Find the mass at the upper bound (ensemble member 1)
         upper_mass = tail_amp_left * &
            normal_cdf(sort_ens(1), tail_mean_left, tail_sd_left)
         ! What fraction of this mass difference should we go?
         fract = curr_q / q(1)
         target_mass = lower_mass + fract * (upper_mass - lower_mass)
         x(i) = inv_weighted_normal_cdf(tail_amp_left, tail_mean_left, &
            tail_sd_left, target_mass)
      endif
   
   elseif(region == ens_size) then
      ! Upper tail
      if(bounded_above .and. do_uniform_tail_right) then
         ! Upper tail is uniform
         lower_state = sort_ens(ens_size)
         upper_state = upper_bound
         x(i) = lower_state + (curr_q - q(ens_size)) * &
            (upper_state - lower_state) / (1.0_r8 - q(ens_size))
      else
         ! Upper tail is (bounded) normal
         ! Find the mass at the upper bound (which could be unbounded)
         if(bounded_above) then
            upper_mass = tail_amp_right * &
               normal_cdf(upper_bound, tail_mean_right, tail_sd_right)
         else
            upper_mass = 1.0_r8
         endif
         ! Find the mass at the lower edge of the region (ensemble member n)
         lower_mass = tail_amp_right * &
            normal_cdf(sort_ens(ens_size), tail_mean_right, tail_sd_right)
         ! What fraction of the last interval do we need to move
         fract = (curr_q - q(ens_size)) / (1.0_r8 - q(ens_size))
         target_mass = lower_mass + fract * (upper_mass - lower_mass)
         x(i) = inv_weighted_normal_cdf(tail_amp_right, tail_mean_right, &
            tail_sd_right, target_mass)
      endif
   
   else
      ! Interior region; get the quantiles of the region boundary
      lower_q = q(region)
      upper_q = q(region + 1)
      x(i) = sort_ens(region) + ((curr_q - lower_q) / (upper_q - lower_q)) * &
         (sort_ens(region + 1) - sort_ens(region))
   endif
   
   ! Imprecision can lead to x being slightly out of bounds, fix it to bounds
   call check_bounds(x(i), curr_q, bounded_below, lower_bound, &
                              bounded_above, upper_bound, 'inf_bnrh_cdf')
enddo
   
end subroutine inv_bnrh_cdf

!-----------------------------------------------------------------------


subroutine inv_bnrh_cdf_like(quantiles, ens_size, sort_ens, &
   bounded_below, bounded_above, lower_bound, upper_bound, &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left,  &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, x, &
   like)

integer,  intent(in)    :: ens_size
real(r8), intent(in)    :: quantiles(ens_size)
real(r8), intent(in)    :: sort_ens(ens_size)
logical,  intent(in)    :: bounded_below, bounded_above
real(r8), intent(in)    :: lower_bound, upper_bound
real(r8), intent(in)    :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8), intent(in)    :: tail_amp_right, tail_mean_right, tail_sd_right
logical,  intent(in)    :: do_uniform_tail_left, do_uniform_tail_right
real(r8), intent(out)   :: x(ens_size)
real(r8), intent(inout) :: like(ens_size)

! This inverts the cdf which is optionally multiplied by a likelihood.

integer :: region, i, j
real(r8) :: lower_state, upper_state, lower_mass, upper_mass, target_mass
real(r8) :: q(ens_size), curr_q, amp_adj, lower_q, upper_q, del_q, fract

! Quantile increment between ensemble members for bnrh
del_q = 1.0_r8 / (ens_size + 1.0_r8)

! Normalize the likelihood to have a sum of 1
like = like / (sum(like) + like(1) / 2.0_r8 + like(ens_size) / 2.0_r8)

! Go from left to right adjusting the quantiles through the x's
! Assume that the quantiles of the original ensemble for the BNRH are uniform
q(1) = like(1) 
do i = 2, ens_size
   q(i) = q(i - 1) + (like(i-1) + like(i)) / 2.0_r8
end do

! Temporary test to confirm posterior is a pdf
if(abs(q(ens_size) + like(ens_size) - 1.0_r8) > 1.0e-12) then
   write(*, *) 'final q ', q(ens_size) + like(ens_size)
   stop
endif

! Loop through each ensemble member to find posterior state
do i = 1, ens_size
   curr_q = quantiles(i)
   ! Which region is this quantile in?
   ! Find which region this quantile is in
   ! Need to make this more efficient once it is working
   ! Default is that region is the highest one; quantile(i) >= largest q
   region = ens_size
   do j = 1, ens_size
      if(curr_q < q(j)) then 
         region = j - 1
         exit
      endif
   end do

   if(region == 0) then
      ! Lower tail
      if(bounded_below .and. do_uniform_tail_left) then
         ! Lower tail uniform
         upper_state = sort_ens(1)
         x(i) = lower_bound + (curr_q / q(1)) * (upper_state - lower_bound)
      else
         ! Find the mass at the lower bound (which could be unbounded)
         ! The amplitude is changed if there is a non-uniform likelihood
         amp_adj = q(1) / del_q
         if(bounded_below) then
            lower_mass = amp_adj * tail_amp_left * &
               normal_cdf(lower_bound, tail_mean_left, tail_sd_left)
         else
            lower_mass = 0.0_r8
         endif
         ! Find the mass at the upper bound (ensemble member 1)
         upper_mass = amp_adj * tail_amp_left * &
            normal_cdf(sort_ens(1), tail_mean_left, tail_sd_left)
         ! What fraction of this mass difference should we go?
         fract = curr_q / q(1)
         target_mass = lower_mass + fract * (upper_mass - lower_mass)
         x(i) = inv_weighted_normal_cdf(amp_adj*tail_amp_left, tail_mean_left, &
            tail_sd_left, target_mass)
      endif
   
   elseif(region == ens_size) then
      ! Upper tail
      if(bounded_above .and. do_uniform_tail_right) then
         ! Upper tail is uniform
         lower_state = sort_ens(ens_size)
         upper_state = upper_bound
         x(i) = lower_state + (curr_q - q(ens_size)) * &
            (upper_state - lower_state) / (1.0_r8 - q(ens_size))
      else
         ! Upper tail is (bounded) normal
         ! Find the mass at the upper bound (which could be unbounded)
         amp_adj = (1.0_r8 - q(ens_size)) / del_q
         if(bounded_above) then
            upper_mass = amp_adj * tail_amp_right * &
               normal_cdf(upper_bound, tail_mean_right, tail_sd_right)
         else
            upper_mass = amp_adj * 1.0_r8
         endif
         ! Find the mass at the lower edge of the region (ensemble member n)
         lower_mass = amp_adj * tail_amp_right * &
            normal_cdf(sort_ens(ens_size), tail_mean_right, tail_sd_right)
         ! What fraction of the last interval do we need to move
         fract = (curr_q - q(ens_size)) / (1.0_r8 - q(ens_size))
         target_mass = lower_mass + fract * (upper_mass - lower_mass)
         x(i) = inv_weighted_normal_cdf(amp_adj * tail_amp_right, tail_mean_right, &
            tail_sd_right, target_mass)
      endif
   
   else
      ! Interior region; get the quantiles of the region boundary
      lower_q = q(region)
      upper_q = q(region + 1)
      x(i) = sort_ens(region) + ((curr_q - lower_q) / (upper_q - lower_q)) * &
         (sort_ens(region + 1) - sort_ens(region))
   endif
   
   ! Imprecision can lead to x being slightly out of bounds, fix it to bounds
   call check_bounds(x(i), curr_q, bounded_below, lower_bound, &
                              bounded_above, upper_bound, 'inf_bnrh_cdf_like')
enddo

end subroutine inv_bnrh_cdf_like

!-----------------------------------------------------------------------

subroutine ens_quantiles(sorted_ens, ens_size, bounded_below, bounded_above, &
                         lower_bound, upper_bound, q)

! Given sorted ensemble which may have members identical to the bounds or may contain
! duplicates, compute the quantiles for each member in an bounded normal rh distribution

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: sorted_ens(ens_size)
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
      if(sorted_ens(i) == lower_bound) then
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
      if(sorted_ens(i) == upper_bound) then
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
   if(sorted_ens(i) == sorted_ens(i - 1)) then
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
   q(i) = 1.0_r8 - upper_dups / (2.0_r8 * (ens_size + 1.0_r8))
end do

! Do the interior series
do i = 1, series_num
   do j = series_start(i), series_end(i)
      q(j) = series_start(i) / (ens_size + 1.0_r8) + (series_length(i) - 1.0_r8) / (2.0_r8 * (ens_size + 1.0_r8))
   end do
end do

end subroutine ens_quantiles

!-----------------------------------------------------------------------

subroutine pack_bnrh_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
   do_uniform_tail_left, do_uniform_tail_right, tail_amp_left, tail_amp_right, &
   tail_mean_left, tail_mean_right, tail_sd_left, tail_sd_right, sort_ens, p)

integer,                        intent(in)    :: ens_size
logical,                        intent(in)    :: bounded_below,        bounded_above
real(r8),                       intent(in)    :: lower_bound,          upper_bound
logical,                        intent(in)    :: do_uniform_tail_left, do_uniform_tail_right
real(r8),                       intent(in)    :: tail_amp_left,        tail_amp_right
real(r8),                       intent(in)    :: tail_mean_left,       tail_mean_right
real(r8),                       intent(in)    :: tail_sd_left,         tail_sd_right
real(r8),                       intent(in)    :: sort_ens(ens_size)
type(distribution_params_type), intent(inout) :: p

! Set the fixed storage parameters in the distribution_params_type
p%bounded_below = bounded_below;    p%lower_bound = lower_bound
p%bounded_above = bounded_above;    p%upper_bound = upper_bound
p%ens_size = ens_size

! Allocate space needed for the parameters
allocate(p%ens(ens_size))
allocate(p%more_params(2*4))

! Save the sorted bnrh ensemble values
p%ens = sort_ens

! Store the extra information about the distribution in the more_params array
if(do_uniform_tail_left) then 
   p%more_params(1) = 1.0_r8
else
   p%more_params(1) = 0.0_r8
endif
if(do_uniform_tail_right) then 
   p%more_params(2) = 1.0_r8
else
   p%more_params(2) = 0.0_r8
endif

p%more_params(3) = tail_amp_left;    p%more_params(4) = tail_amp_right
p%more_params(5) = tail_mean_left;   p%more_params(6) = tail_mean_right
p%more_params(7) = tail_sd_left;     p%more_params(8) = tail_sd_right

end subroutine pack_bnrh_params

!-----------------------------------------------------------------------

subroutine unpack_bnrh_params(p, do_uniform_tail_left, do_uniform_tail_right, &
   tail_amp_left, tail_amp_right, tail_mean_left, tail_mean_right, tail_sd_left, tail_sd_right)

! Unpack values describing the bnrh distribution from the distribution_params_type p

type(distribution_params_type), intent(in)  :: p
logical,                        intent(out) :: do_uniform_tail_left, do_uniform_tail_right
real(r8),                       intent(out) :: tail_amp_left,        tail_amp_right
real(r8),                       intent(out) :: tail_mean_left,       tail_mean_right
real(r8),                       intent(out) :: tail_sd_left,         tail_sd_right

! Logicals are stored as 1 for true, 0 for false   
do_uniform_tail_left  = p%more_params(1) > 0.5_r8  
do_uniform_tail_right = p%more_params(2) > 0.5_r8

tail_amp_left  = p%more_params(3);    tail_amp_right  = p%more_params(4)
tail_mean_left = p%more_params(5);    tail_mean_right = p%more_params(6)
tail_sd_left   = p%more_params(7);    tail_sd_right   = p%more_params(8)

end subroutine unpack_bnrh_params

!-----------------------------------------------------------------------

function get_bnrh_sd(p)

real(r8)                                   :: get_bnrh_sd
type(distribution_params_type), intent(in) :: p

! Return the standard deviation of this distribution
get_bnrh_sd = p%more_params(7)

end function get_bnrh_sd

!-----------------------------------------------------------------------


subroutine check_bounds(x, q, bounded_below, lower_bound, &
                              bounded_above, upper_bound, msgstring)

real(r8),         intent(inout) :: x
real(r8),         intent(in)    :: q
logical,          intent(in)    :: bounded_below, bounded_above
real(r8),         intent(in)    :: lower_bound,   upper_bound
character(len=*), intent(in)    :: msgstring

! Imprecision in inv_norm could lead to x(i) being out of bounds: check for now
! lower bound. Correct this and output a message. Could be numerically fixed above.
if(bounded_below) then
   if(x < lower_bound) then
      write(errstring, *) 'x less than lower_bound ', x, q
      call error_handler(E_MSG, msgstring, errstring, source)
      x = lower_bound
   endif
endif

! See comment on lower bound in previous code block 
if(bounded_above) then
   if(x > upper_bound) then
      write(errstring, *) 'x greater than upper_bound ', x, q
      call error_handler(E_MSG, msgstring, errstring, source)
      x = upper_bound
   endif
endif

end subroutine check_bounds   

!-----------------------------------------------------------------------

end module bnrh_distribution_mod
