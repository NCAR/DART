! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module rh_distribution_mod

use types_mod,      only : r8

use utilities_mod,  only : E_ERR, error_handler

use sort_mod, only : index_sort

use normal_distribution_mod, only : norm_cdf, norm_inv

implicit none
private

public :: ens_quantiles, rh_cdf_init

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
real(r8) :: base_prob, mean, bound_quantile
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
sort_x(1:ens_size) = x(sort_index)

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

! For unit normal, find distance from mean to where cdf is 1/(ens_size+1).
! Saved to avoid redundant computation for repeated calls with same ensemble size
base_prob = 1.0_r8 / (ens_size + 1.0_r8)
if(saved_ens_size /= ens_size) then
   call norm_inv(base_prob, dist_for_unit_sd)
   ! This will be negative, want it to be a distance so make it positive
   dist_for_unit_sd = -1.0_r8 * dist_for_unit_sd
   ! Keep a record of the ensemble size used to compute dist_for_unit_sd
   saved_ens_size = ens_size
endif

! Find a mean so that 1 / (ens_size + 1) probability is in outer regions
tail_mean_left =  sort_x(1)        + dist_for_unit_sd * tail_sd_left
tail_mean_right = sort_x(ens_size) - dist_for_unit_sd * tail_sd_right

! If the distribution is bounded, still want 1 / (ens_size + 1) in outer regions
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
   ! Note that due to roundoff it is possible for base_prob - quantile to be slightly negative
   if((base_prob - bound_quantile) / base_prob < uniform_threshold) then
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
   bound_quantile = norm_cdf(upper_bound, tail_mean_right, tail_sd_right)
   ! Note that due to roundoff it is possible for the numerator to be slightly negative
   if((bound_quantile - (1.0_r8 - base_prob)) / base_prob < uniform_threshold) then
      ! If bound and ensemble member are too close, do uniform approximation
      do_uniform_tail_right = .true.
   else
      ! Compute the right tail amplitude
      tail_amp_right = base_prob / (base_prob - (1.0_r8 - bound_quantile))
   endif
endif

end subroutine rh_cdf_init

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
