! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module rh_distribution_mod

use types_mod,      only : r8

use utilities_mod,  only : E_ERR, error_handler

use sort_mod, only : index_sort

use normal_distribution_mod, only : norm_cdf, norm_inv, weighted_norm_inv

implicit none
private

public :: rh_cdf_init, rh_cdf, rh_cdf_ens, inv_rh_cdf

character(len=512)          :: errstring
character(len=*), parameter :: source = 'rh_distribution_mod.f90'

! Saves the ensemble size used in the previous call of obs_inc_bounded_norm_rh
integer :: saved_ens_size = -99

! Parameter to control switch to uniform approximation for normal tail
! This defines how many quantiles the bound is from the outermost ensemble member
! If closer than this, can get into precision error problems with F(F-1(X)) on tails
! Can also get other precision problems with the amplitude for the normal in the bounded region
real(r8), parameter :: uniform_threshold = 0.01_r8

contains

!-----------------------------------------------------------------------

subroutine rh_cdf_init(x, ens_size, bounded, bounds, sort_x, quantiles,  &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right)

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: x(ens_size)
logical,  intent(in)  :: bounded(2)
real(r8), intent(in)  :: bounds(2)
! Do we really want to force the sort to happen here?
real(r8), intent(out) :: sort_x(ens_size)
real(r8), intent(out) :: quantiles(ens_size)
real(r8), intent(out) :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8), intent(out) :: tail_amp_right, tail_mean_right, tail_sd_right
logical,  intent(out) :: do_uniform_tail_left, do_uniform_tail_right

real(r8), save :: dist_for_unit_sd
real(r8) :: q(ens_size)
real(r8) :: del_q, mean, bound_quantile
real(r8) :: lower_bound, upper_bound
logical  :: bounded_below, bounded_above
integer :: sort_index(ens_size), indx, i

! Computes all information about a rank histogram cdf given the ensemble and bounds

! Clarity of use for bounds
lower_bound = bounds(1)
upper_bound = bounds(2)
bounded_below = bounded(1)
bounded_above = bounded(2)

! Get ensemble mean and sd
mean = sum(x) / ens_size
tail_sd_left  = sqrt(sum((x - mean)**2) / (ens_size - 1))
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
      call error_handler(E_ERR, 'rh_cdf_init', errstring, source)
   endif
endif
  
! Fail if upper bound is smaller than the largest ensemble member 
if(bounded_above) then
   if(sort_x(ens_size) > upper_bound) then
      write(errstring, *) 'Largest ensemble member greater than upper bound', &
         sort_x(ens_size), upper_bound
      call error_handler(E_ERR, 'rh_cdf_init', errstring, source)
   endif
endif

! Get the quantiles for each of the ensemble members in a RH distribution
call ens_quantiles(sort_x, ens_size, &
   bounded_below, bounded_above, lower_bound, upper_bound, q)
! Put sorted quantiles back into input ensemble order
do i = 1, ens_size
   indx = sort_index(i)
   quantiles(indx) = q(i)
end do

! Compute the characteristics of tails

! For unit normal, find distance from mean to where cdf is 1/(ens_size+1) (del_q_.
! Saved to avoid redundant computation for repeated calls with same ensemble size
del_q = 1.0_r8 / (ens_size + 1.0_r8)

if(saved_ens_size /= ens_size) then
   call norm_inv(del_q, dist_for_unit_sd)
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
   bound_quantile = norm_cdf(lower_bound, tail_mean_left, tail_sd_left)
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
   bound_quantile = norm_cdf(upper_bound, tail_mean_right, tail_sd_right)
   ! Note that due to roundoff it is possible for the numerator to be slightly negative
   if((bound_quantile - (1.0_r8 - del_q)) / del_q < uniform_threshold) then
      ! If bound and ensemble member are too close, do uniform approximation
      do_uniform_tail_right = .true.
   else
      ! Compute the right tail amplitude
      tail_amp_right = del_q / (del_q - (1.0_r8 - bound_quantile))
   endif
endif

end subroutine rh_cdf_init

!-----------------------------------------------------------------------

subroutine rh_cdf_ens(x, ens_size, sort_ens, bounded_below, bounded_above,   &
   lower_bound, upper_bound,                                                 &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left,     &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right,    &
   quantile)

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: x(ens_size)
real(r8), intent(in)  :: sort_ens(ens_size)
logical,  intent(in)  :: bounded_below, bounded_above
real(r8), intent(in)  :: lower_bound, upper_bound
real(r8), intent(in)  :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8), intent(in)  :: tail_amp_right, tail_mean_right, tail_sd_right
logical,  intent(in)  :: do_uniform_tail_left, do_uniform_tail_right
real(r8), intent(out) :: quantile(ens_size)

real(r8) :: q(ens_size)
integer  :: i

! Get the quantiles for each of the ensemble members in a RH distribution
! This was all computed when the distribution was originally set up, could choose to cache that
! in the params structure for efficiency. This is only used for the single observation posterior
! so one could only save in that case removing any storage concerns.
call ens_quantiles(sort_ens, ens_size, &
   bounded_below, bounded_above, lower_bound, upper_bound, q)

! This can be done vastly more efficiently with either binary searches or by first sorting the
! incoming state_ens so that the lower bound for starting the search is updated with each ensemble member
do i = 1, ens_size
   ! Figure out which bin it is in
   call rh_cdf(x(i), ens_size, sort_ens, bounded_below, bounded_above, lower_bound, upper_bound, &
      tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
      tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, q, quantile(i))
end do

end subroutine rh_cdf_ens

!-----------------------------------------------------------------------

subroutine rh_cdf(x, ens_size, sort_ens, bounded_below, bounded_above,   &
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

real(r8) :: upper_q, fract, del_q
integer :: j

! Quantile increment between ensemble members for rh
del_q = 1.0_r8 / (ens_size + 1.0_r8)

if(x < sort_ens(1)) then
   ! In the left tail
   ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
   if(bounded_below .and. x < lower_bound) then
      write(errstring, *) 'Ensemble member less than lower bound', x, lower_bound
      call error_handler(E_ERR, 'rh_cdf', errstring, source)
      ! This error can occur due to roundoff in increment generation from bounded RHF
      ! See discussion in function fix_bounds.
   endif

   if(do_uniform_tail_left) then
      ! Uniform approximation for left tail
      ! The division here could be a concern. However, if sort_ens(1) == lower_bound, then
      ! x cannot be < sort_ens(1).
      quantile = (x - lower_bound) / (sort_ens(1) - lower_bound) * del_q
   else
      ! It's a normal tail
      if(bounded_below) then
         quantile = tail_amp_left * (norm_cdf(x, tail_mean_left, tail_sd_left) - &
            norm_cdf(lower_bound, tail_mean_left, tail_sd_left))
      else        ! Unbounded, tail normal goes all the way down to quantile 0
         quantile = (tail_amp_left * norm_cdf(x, tail_mean_left, tail_sd_left) / &
                    (tail_amp_left * norm_cdf(sort_ens(1), tail_mean_left, tail_sd_left))) &
                    * del_q
      endif
      ! Make sure it doesn't sneak past the first ensemble member due to round-off
      quantile = min(quantile, del_q)
   endif
elseif(x == sort_ens(1)) then
   ! This takes care of cases where there are multiple rh values at the bdry or at first ensemble
   quantile = q(1)
elseif(x > sort_ens(ens_size)) then
   ! In the right tail
   ! Do an error check to make sure ensemble member isn't outside bounds, may be redundant
   if(bounded_above .and. x > upper_bound) then
      write(errstring, *) 'Ensemble member greater than upper bound first check(see code)', x, upper_bound
      call error_handler(E_ERR, 'rh_cdf', errstring, source)
      ! This error can occur due to roundoff in increment generation from bounded RHF
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
      if(bounded_above) then
         upper_q = tail_amp_right * norm_cdf(upper_bound, tail_mean_right, tail_sd_right)
      else
         upper_q = tail_amp_right
      endif

      ! Want to avoid quantiles exceeding 1 due to numerical issues. Do fraction of the normal part
      fract = (tail_amp_right * norm_cdf(x,                  tail_mean_right, tail_sd_right) - &
               tail_amp_right * norm_cdf(sort_ens(ens_size), tail_mean_right, tail_sd_right)) / &
              (upper_q - tail_amp_right * norm_cdf(sort_ens(ens_size), tail_mean_right, tail_sd_right))
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

end subroutine rh_cdf

!-----------------------------------------------------------------------

subroutine inv_rh_cdf(quantile, ens_size, sort_ens, &
   bounded_below, bounded_above, lower_bound, upper_bound, &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left,  &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, x)

real(r8), intent(in)  :: quantile
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: sort_ens(ens_size)
logical,  intent(in)  :: bounded_below, bounded_above
real(r8), intent(in)  :: lower_bound, upper_bound
real(r8), intent(in)  :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8), intent(in)  :: tail_amp_right, tail_mean_right, tail_sd_right
logical,  intent(in)  :: do_uniform_tail_left, do_uniform_tail_right
real(r8), intent(out) :: x

integer :: region
real(r8) :: lower_state, upper_state, lower_mass, upper_mass, target_mass
real(r8) :: lower_q, upper_q, fract, del_q

! Quantile increment between ensemble members for rh
del_q = 1.0_r8 / (ens_size + 1.0_r8)

! Assume that the quantiles of the original ensemble for the BNRH are uniform
! Finding which region this quantile is in is trivial
region = floor(quantile * (ens_size + 1.0_r8))
! Careful about numerical issues moving outside of region [0 ens_size]
if(region < 0) region = 0
if(region > ens_size) region = ens_size

if(region == 0) then
   ! Lower tail
   if(bounded_below .and. do_uniform_tail_left) then
      ! Lower tail uniform
      upper_state = sort_ens(1)
! NOTE: NEED TO BE CAREFUL OF THE DENOMINATOR HERE AND ON THE PLUS SIDE
      x = lower_bound + &
         (quantile / del_q) * (upper_state - lower_bound)
   else
      ! Find the mass at the lower bound (which could be unbounded)
      if(bounded_below) then
         lower_mass = tail_amp_left * norm_cdf(lower_bound, tail_mean_left, tail_sd_left)
      else
         lower_mass = 0.0_r8
      endif
      ! Find the mass at the upper bound (ensemble member 1)
      upper_mass = tail_amp_left * norm_cdf(sort_ens(1), tail_mean_left, tail_sd_left)
      ! What fraction of this mass difference should we go?
      fract = quantile / del_q
      target_mass = lower_mass + fract * (upper_mass - lower_mass)
      call weighted_norm_inv(tail_amp_left, tail_mean_left, tail_sd_left, target_mass, x)
   endif

elseif(region == ens_size) then
   ! Upper tail
   if(bounded_above .and. do_uniform_tail_right) then
      ! Upper tail is uniform
      lower_state = sort_ens(ens_size)
      upper_state = upper_bound
      x = lower_state + (quantile - ens_size *del_q) * &
         (upper_state - lower_state) / del_q

   else
      ! Upper tail is (bounded) normal
      ! Find the mass at the upper bound (which could be unbounded)
      if(bounded_above) then
         upper_mass = tail_amp_right * norm_cdf(upper_bound, tail_mean_right, tail_sd_right)
      else
         upper_mass = 1.0_r8
      endif
      ! Find the mass at the lower bound (ensemble member n)
      lower_mass = tail_amp_right * norm_cdf(sort_ens(ens_size), tail_mean_right, tail_sd_right)
      ! What fraction of the last interval do we need to move
      fract = (quantile - ens_size * del_q) / del_q
      target_mass = lower_mass + fract * (upper_mass - lower_mass)
      call weighted_norm_inv(tail_amp_right, tail_mean_right, tail_sd_right, target_mass, x)
   endif

else
   ! Interior region; get the quantiles of the region boundary
   lower_q = region * del_q
   upper_q = (region + 1.0_r8) / (ens_size + 1.0_r8)
   x = sort_ens(region) + ((quantile - lower_q) / (upper_q - lower_q)) * &
      (sort_ens(region + 1) - sort_ens(region))
endif

! Check for posterior violating bounds; This may not be needed after development testing
if(bounded_below) then
   if(x < lower_bound) then
      write(errstring, *) 'x less than lower_bound ', x
      call error_handler(E_ERR, 'inv_rh_cdf', errstring, source)
   endif
endif

if(bounded_above) then
   if(x > upper_bound) then
      write(errstring, *) 'x greater than upper_bound ', x
      call error_handler(E_ERR, 'inv_rh_cdf', errstring, source)
   endif
endif

end subroutine inv_rh_cdf

!-----------------------------------------------------------------------

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

end module rh_distribution_mod
