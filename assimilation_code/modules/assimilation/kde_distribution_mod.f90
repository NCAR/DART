module kde_distribution_mod

use types_mod,               only : r8, missing_r8

use utilities_mod,           only : E_ERR, E_MSG, error_handler

use sort_mod,                only : sort

use distribution_params_mod, only : distribution_params_type, deallocate_distribution_params, &
                                    KDE_DISTRIBUTION

use rootfinding_mod,         only : inv_cdf

use normal_distribution_mod, only : normal_cdf

implicit none
private

public :: kde_cdf, kde_cdf_params, inv_kde_cdf, inv_kde_cdf_params,       &
          test_kde, obs_dist_types, likelihood_function, pack_kde_params, &
          separate_ensemble

type available_obs_dist_types
   integer  :: uninformative, normal, binomial, gamma, &
               inv_gamma, lognormal, truncated_normal
end type

type(available_obs_dist_types), parameter :: obs_dist_types = &
available_obs_dist_types(uninformative=999, normal=998, binomial=997, &
   gamma=996, inv_gamma=995, lognormal=994, truncated_normal=993)
character(len=512)          :: errstring
character(len=*), parameter :: source = 'kde_distribution_mod.f90'

contains

!---------------------------------------------------------------------------

function likelihood_function(x, y, obs_param, obs_dist_type, &
        bounded_above, bounded_below, upper_bound, lower_bound) result(l)
   real(r8)             :: l  ! likelihood value
   real(r8), intent(in) :: x  ! state value
   real(r8), intent(in) :: y  ! obs value
   real(r8), intent(in) :: obs_param  ! meaning depends on obs_dist_type
   integer,  intent(in) :: obs_dist_type ! obs distribution type
   logical, optional, intent(in)  :: bounded_above, bounded_below
   real(r8), optional, intent(in) :: upper_bound, lower_bound

   ! Evaluates the likelihood pdf(y | x) for various kinds of observation
   ! distributions. The value returned is not equal to the observation
   ! pdf evaluated at y, because normalization constants that don't depend
   ! on x are omitted.

   real(r8) :: gamma_shape, gamma_scale
   real(r8) :: inv_gamma_shape, inv_gamma_scale
   real(r8) :: cdf(2)
   logical  :: bounded_above_val, bounded_below_val

   l = 1._r8 ! Initialize

   select case (obs_dist_type)
      case (obs_dist_types%uninformative)
         ! Uninformative observations have a likelihood equal to one.
         l = 1._r8
      case (obs_dist_types%normal)
         ! For a normal obs distribution, like_param is the obs error variance
         if (obs_param <= 0._r8) then
            write(errstring, *) 'obs error sd <= 0', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            l = exp( -0.5_r8 * (x - y)**2 / obs_param )
         end if
      case (obs_dist_types%binomial)
         ! For a binomial obs distribution 0<=x<=1 is the probability of
         ! observing y successes of obs_param total trials
         if (y < 0._r8) then
            write(errstring, *) 'y value is negative with a binomial obs model ', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (y > obs_param) then
            write(errstring, *) 'successes greater than total trials with a binomial obs model ', y, obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif ((x < 0._r8) .or. (x > 1._r8)) then
            write(errstring, *) 'x outside [0,1] with a binomial obs model ', x
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            l = (x**y) * (1._r8 - x)**(obs_param - y)
         end if
      case (obs_dist_types%gamma)
         ! For a gamma obs distribution, the mean is x and the variance is obs_param * x^2, i.e.
         ! the obs error sd is sqrt(obs_param) times  the true value. If the error sd is p% of x,
         ! set obs_param = (p/100._r8)**2.
         ! For a gamma obs distribution, the likelihood is inverse gamma.
         if (x < 0._r8) then
            write(errstring, *) 'x value is negative with a gamma obs model ', x
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (y <= 0._r8) then
            write(errstring, *) 'y value is non-positive with a gamma obs model ', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (obs_param <= 0._r8) then
            write(errstring, *) 'obs variance is non-positive with a gamma obs model ', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            gamma_shape = 1._r8 / obs_param ! k
            gamma_scale = x * obs_param     ! theta
            if (x .eq. 0._r8) then
               l = 0._r8 ! Technically x must be > 0, but just setting l=0 is convenient.
            else
               l = (y / gamma_scale)**gamma_shape * exp(-y / gamma_scale)
            end if
         end if
      case (obs_dist_types%inv_gamma)
         ! For an inverse gamma obs distribution, the mean is x and the variance is obs_param * x^2,
         ! i.e. the obs error sd is sqrt(obs_param) times  the true value. If the error sd is p% of x,
         ! set obs_param = (p/100._r8)**2.
         ! For an inverse gamma obs distribution, the likelihood is gamma.
         if (x < 0._r8) then
            write(errstring, *) 'x value is negative with an inverse gamma obs model ', x
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (y <= 0._r8) then
            write(errstring, *) 'y value is non-positive with an inverse gamma obs model ', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (obs_param <= 0._r8) then
            write(errstring, *) 'obs variance is non-positive with an inverse gamma obs model ', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            inv_gamma_shape = (1._r8 / obs_param) + 2._r8   ! alpha
            inv_gamma_scale = x * (inv_gamma_shape - 1._r8) ! beta
            if (x .eq. 0._r8) then
               l = 0._r8 ! Technically x must be > 0, but just setting l=0 is convenient.
            else
               l = (inv_gamma_scale**inv_gamma_shape) * exp( - inv_gamma_scale / y )
            end if
         end if
      case (obs_dist_types%lognormal)
         ! For a lognormal obs distribution, ln(y) is normal with mean x and variance obs_param.
         ! The likelihood is normal.
         if (y <= 0._r8) then
            write(errstring, *) 'y value is non-positive with a lognormal obs model', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (obs_param <= 0._r8) then
            write(errstring, *) 'obs error sd <= 0', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            l = exp( -0.5_r8 * (x - log(y))**2 / obs_param )
         end if
      case (obs_dist_types%truncated_normal)
         ! Code based on the version in assim_tools_mod.f90
         ! This implicitly assumes that the bounds on the y distribution are the same
         ! as the bounds on x, which makes sense in observation space. The factor of
         ! sigma * sqrt(2 * pi) is ignored.
         ! A zero observation error variance is a degenerate case
         if(obs_param <= 0.0_r8) then
            if(x == y) then
               l = 1.0_r8
            else
               l = 0.0_r8
            endif
            return
         endif

         cdf = [0._r8, 1._r8]
         if (present(bounded_above)) then
            bounded_above_val = bounded_above
         else
            bounded_above_val = .false.
         end if

         if (present(bounded_below)) then
            bounded_below_val = bounded_below
         else
            bounded_below_val = .false.
         end if

         if (bounded_below_val) cdf(1) = normal_cdf(lower_bound, x, sqrt(obs_param))
         if (bounded_above_val) cdf(2) = normal_cdf(upper_bound, x, sqrt(obs_param))

         l = exp( -0.5_r8 * (x - y)**2 / obs_param ) / (cdf(2) - cdf(1))
      case DEFAULT
         write(errstring, *) 'likelihood called with unrecognized obs_dist_type ', obs_dist_type
         call error_handler(E_MSG, 'kde_distribution_mod:likelihood_function', errstring, source)
   end select

end function likelihood_function

!---------------------------------------------------------------------------

elemental function epanechnikov_kernel(x) ! The fact that this is elemental is not (yet) used
   real(r8) :: epanechnikov_kernel
   real(r8), intent(in) :: x
   real(r8), parameter  :: norm_const = 3._r8 / 4._r8
   epanechnikov_kernel = norm_const * max(0._r8, 1._r8 - x**2)

end function epanechnikov_kernel

!---------------------------------------------------------------------------

elemental function epanechnikov_cdf(x) ! The fact that this is elemental is not (yet) used
   real(r8) :: epanechnikov_cdf
   real(r8), intent(in) :: x
   real(r8), parameter  :: norm_const = 1._r8 / 4._r8
   if (x <= -1._r8) then
      epanechnikov_cdf = 0._r8
   elseif (x <= 1._r8) then
      epanechnikov_cdf = min(1._r8, max(0._r8, norm_const * (2._r8 + 3._r8 * x - x**3)))
   else
      epanechnikov_cdf = 1._r8
   end if

end function epanechnikov_cdf

!---------------------------------------------------------------------------

elemental subroutine boundary_correction(x, lx, mx)
   real(r8), intent(in)  :: x
   real(r8), intent(out) :: lx, mx

   ! Boundary correction for kde, from Jones (1993) and Jones & Foster (1996)
   ! The fact that this is elemental is not (yet) used

   real(r8) :: denom

   denom = ((1._r8 + x)**4) * (19._r8 + 3._r8 * x * (x - 6._r8))
   lx =-64._r8 * (-2._r8 + x * (4._r8 + 3._r8 * x * (x - 2._r8))) / denom
   mx = 240._r8 * (x - 1._r8)**2 / denom

end subroutine boundary_correction

!---------------------------------------------------------------------------

function kde_pdf(x, p)
   real(r8)                                   :: kde_pdf
   real(r8),                       intent(in) :: x
   type(distribution_params_type), intent(in) :: p

   ! Returns the kernel estimate of the pdf evaluated at x.

   integer  :: i
   real(r8) :: u_lower, lx_lower, mx_lower
   real(r8) :: u_upper, lx_upper, mx_upper

   kde_pdf = 0._r8 ! Initialize
   if (p%bounded_below .and. (x <= p%lower_bound)) return
   if (p%bounded_above .and. (x >= p%upper_bound)) return
   do i=1, p%ens_size ! Reduction loop - parallelizable
      u_lower = 0._r8; lx_lower = 1._r8 ; mx_lower = 0._r8
      u_upper = 0._r8; lx_upper = 1._r8 ; mx_upper = 0._r8
      if (p%bounded_below) then ! Bounded below
         u_lower = min( 1._r8, max( 0._r8, (x - p%lower_bound) / p%more_params(i) ) ) ! p%more_params(i) holds kernel width for ensemble member i
         call boundary_correction(u_lower, lx_lower, mx_lower)
      end if
      if (p%bounded_above) then ! Bounded above
         u_upper = min( 1._r8, max( 0._r8, (p%upper_bound - x) / p%more_params(i) ) )
         call boundary_correction(u_upper, lx_upper, mx_upper)
      end if
      u_lower = (x - p%ens(i)) / p%more_params(i) ! Not u_lower any more, just (x-x_i)/h_i
      kde_pdf = kde_pdf + (1._r8 / p%more_params(i)) * &
                          (lx_lower + mx_lower * u_lower) * &
                          (lx_upper - mx_upper * u_lower) * &
                          epanechnikov_kernel( u_lower )
   end do
   kde_pdf = max(0._r8, kde_pdf) / (p%ens_size * p%more_params(p%ens_size + 1)) ! p%more_params(ens_size + 1) normalizes the pdf

end function kde_pdf

!-----------------------------------------------------------------------

subroutine get_kde_bandwidths(ens_size, ens, bandwidths)
   integer,         intent(in)   :: ens_size
   real(r8),        intent(in)   :: ens(ens_size)
   real(r8),        intent(out)  :: bandwidths(ens_size)

   real(r8)                      :: ens_mean, ens_sd, h0, g, d_max
   real(r8), dimension(ens_size) :: f_tilde, d, lambda
   integer                       :: i, j, k


   d(:) = abs( ens(:) - ens(1) )
   d_max = maxval(d(:))
   if (d_max <= 0._r8) then
      errstring = 'Bandwidth = 0 '
      call error_handler(E_ERR, 'get_kde_bandwidths', errstring, source)
   end if
   ens_mean = sum(ens) / ens_size
   ens_sd   = sqrt( sum( (ens - ens_mean)**2 ) / (ens_size - 1._r8) )
   h0 = 2._r8 * ens_sd / (ens_size**0.2_r8) ! This would be the kernel width if the widths were not adaptive.
                                            ! It would be better to use min(sd, iqr/1.34) but don't want to compute iqr
   k = floor( sqrt( real(ens_size) ) ) ! distance to kth nearest neighbor used to set bandwidth
   do i=1,ens_size
      d(:) = sort( abs( ens(:) - ens(i) ) )
      d_max = maxval(d(:))
      where (d < 1.e-3_r8 * d_max) ! set small distances to zero
         d = 0._r8
      end where
      j = 1
      do while ((d(k+j) <= 0._r8) .and. (k+j < ens_size))
         j = j + 1
      end do
      f_tilde(i) = 0.5_r8 * real(k+j-1, r8) / (real(ens_size, r8) * d(k+j)) ! Initial density estimate
   end do
   f_tilde(:) = f_tilde(:) / maxval(f_tilde(:)) ! Avoids overflow in the next line
   g = product( f_tilde )**( 1._r8 / real(ens_size, r8) )
   lambda(:) = sqrt( g / f_tilde(:) )
   bandwidths(:) = h0 * lambda(:)
   if (maxval(bandwidths(:)) <= 0._r8) then
      write(*,*) 'Bandwidth breakdown'
      write(*,*) h0, g, lambda(:), f_tilde(:), ens(:)
   end if

end subroutine get_kde_bandwidths

!---------------------------------------------------------------------------

function gq5(left, right, p) result(q)
   real(r8)                                   :: q
   real(r8),                       intent(in) :: left, right
   type(distribution_params_type), intent(in) :: p

   ! Uses 5-point Gauss quadrature to approximate \int_left^right l(s; y) p(s) ds
   !  where p(x) is the prior pdf and l(x; y) is the likelihood. This only works
   !  correctly if left >= lower_bound and right <= upper_bound. The result is
   ! **Not Normalized**.

   real(r8) :: y
   real(r8) :: obs_param ! See likelihood function for interpretation
   integer  :: obs_dist_type  ! See likelihood function for interpretation
   real(r8) :: xi ! quadrature point
   real(r8), save :: chi(5) = [-sqrt(5._r8 + 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8, &
                               -sqrt(5._r8 - 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8, &
                                0._r8, &
                                sqrt(5._r8 - 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8, &
                                sqrt(5._r8 + 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8] ! Gauss quadrature points
   real(r8), save :: w(5)   = [(322._r8 - 13._r8 * sqrt(70._r8)) / 900._r8, &
                               (322._r8 + 13._r8 * sqrt(70._r8)) / 900._r8, &
                               128._r8/225._r8, &
                               (322._r8 + 13._r8 * sqrt(70._r8)) / 900._r8, &
                               (322._r8 - 13._r8 * sqrt(70._r8)) / 900._r8] ! GQ weights
   real(r8) :: l ! value of the likelihood function
   integer  :: k

   ! Unpack obs info from param struct
   y         = p%more_params(p%ens_size + 2)
   obs_param = p%more_params(p%ens_size + 3)
   obs_dist_type  = nint(p%more_params(p%ens_size + 4))

   q = 0._r8
   do k=1,5
      xi = 0.5_r8 * ((right - left) * chi(k) + left + right)
      if (obs_dist_type .eq. obs_dist_types%truncated_normal) then
         l = likelihood_function(xi, y, obs_param, obs_dist_type, &
            bounded_above=p%bounded_above, bounded_below=p%bounded_below, &
            upper_bound=p%upper_bound, lower_bound=p%lower_bound)
      else
         l = likelihood_function(xi, y, obs_param, obs_dist_type)
      end if
      q  = q + 0.5_r8 * (right - left) * w(k) * kde_pdf(xi, p) * l
   end do

end function gq5

!---------------------------------------------------------------------------

function integrate_pdf(x, p) result(q)
   real(r8)                                   :: q
   real(r8),                       intent(in) :: x
   type(distribution_params_type), intent(in) :: p

   ! Uses quadrature to approximate \int_{-\infty}^x l(x; y) p(s) ds where
   ! p(s) is the prior pdf and l(x; y) is the likelihood. The interval is
   ! broken up into sub-intervals whose boundaries are either one of the bounds,
   ! or the edge of support of one of the kernels, or the value of x. On each
   ! sub-interval the integral is approximated using Gauss-Legendre quadrature
   ! with 5 points. When the likelihood is flat and the boundaries are far
   ! from the ensemble, the result is exact up to roundoff error.

   real(r8) :: y
   real(r8) :: obs_param ! See likelihood function for interpretation
   integer  :: obs_dist_type  ! See likelihood function for interpretation
   real(r8) :: edges(2*p%ens_size)
   real(r8) :: left, right ! edges of current sub-interval, quadrature point
   integer  :: i

   ! Unpack obs info from param struct
   y         = p%more_params(p%ens_size + 2)
   obs_param = p%more_params(p%ens_size + 3)
   obs_dist_type  = nint(p%more_params(p%ens_size + 4))

   edges(1:p%ens_size)                = p%ens(1:p%ens_size) - p%more_params(1:p%ens_size)
   edges(p%ens_size + 1:2*p%ens_size) = p%ens(1:p%ens_size) + p%more_params(1:p%ens_size)
   edges(:) = sort(edges(:)) ! If bandwidths were constant we would not need to sort

   ! If x is outside the support of the pdf then we can skip the quadrature.
   left = edges(1)
   if (p%bounded_below) left = max(left, p%lower_bound)
   if (x <= left) then
      q = 0._r8
      return
   end if

   ! Important to use x > upper_bound here because I use
   ! x = upper_bound to compute the normalization constant.
   if ((p%bounded_above) .and. (x > p%upper_bound)) then
      q = 1._r8
      return
   end if

   ! If we haven't returned yet, then there is at least one subinterval.
   i = 1
   right = min(x, edges(2)) ! left was computed above
   q = gq5(left, right, p)
   do while ((x > right) .and. (i+1 < 2*p%ens_size))
      i     = i + 1
      left  = right
      right = min(x, edges(i+1))
      q     = q + gq5(left, right, p)
   end do
   ! Note that it is possible to have maxval(edges) < x < upper_bound,
   ! but that last sub-interval from maxval(edges) to x has zero integral,
   ! so it can be safely skipped.

end function integrate_pdf

!-----------------------------------------------------------------------

subroutine pack_kde_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
   ens, y, obs_param, obs_dist_type, p)
   integer,                        intent(in)  :: ens_size
   logical,                        intent(in)  :: bounded_below,        bounded_above
   real(r8),                       intent(in)  :: lower_bound,          upper_bound
   real(r8),                       intent(in)  :: ens(ens_size)
   real(r8),                       intent(in)  :: y
   real(r8),                       intent(in)  :: obs_param
   integer,                        intent(in)  :: obs_dist_type
   type(distribution_params_type), intent(out) :: p

   real(r8) :: bandwidths(ens_size)
   real(r8) :: edge, edges(2*ens_size + 2)
   integer  :: i

   ! Set the fixed storage parameters in the distribution_params_type
   p%ens_size = ens_size
   p%distribution_type = KDE_DISTRIBUTION

   ! Allocate space needed for the parameters
   allocate(p%ens(1:ens_size), source=0._r8)
   allocate(p%more_params(1:5*ens_size + 8), source=0._r8)

   ! p%more_params(1:ens_size) are the kernel bandwidths
   ! p%more_params(ens_size + 1) is the normalization constant for the pdf
   ! p%more_params(ens_size + 2) is the observation value y; not used for prior
   ! p%more_params(ens_size + 3) is the observation parameter (see likelihood for details); not used for prior
   ! p%more_params(ens_size + 4) is the observation error distribution type. For prior set to uninformative
   ! p%more_params(ens_size + 5:3*ens_size + 6) are the edges.
   ! p%more_params(3*ens_size+7:5*ens_size+8) are the cdf evaluated at the edges.

   ! Save the ensemble values, sorted now so that they don't have to be re-sorted for the first guess
   p%ens(:) = sort(ens(:))

   ! Store the kernel bandwidths in the more_params array
   ! Important to make sure that p%ens and p%more_params are sorted in same order by passing p%ens rather than ens
   call get_kde_bandwidths(ens_size, p%ens, bandwidths)
   p%more_params(1:ens_size) = bandwidths(:)

   ! Pack obs information
   p%more_params(ens_size + 2) = y
   p%more_params(ens_size + 3) = obs_param
   p%more_params(ens_size + 4) = real(obs_dist_type, r8) ! This is not ideal because it involves a type conversion from int to real

  ! Get the edges of the subintervals on which the pdf is smooth
   edges(1:ens_size)            = p%ens(:) - bandwidths(:)
   edges(ens_size+1:2*ens_size) = p%ens(:) + bandwidths(:)
   if (bounded_below) then
      edges(2*ens_size+1) = lower_bound
   else
      edges(2*ens_size+1) = minval(edges(1:ens_size))
   end if
   if (bounded_above) then
      edges(2*ens_size+2) = upper_bound
   else
      edges(2*ens_size+2) = maxval(edges(ens_size+1:2*ens_size))
   end if
   edges(:) = sort(edges(:))
   p%more_params(ens_size+5:3*ens_size+6) = edges(:)

   ! If the ensemble is sufficiently far from the boundary, that we can use the unbounded code, which is cheaper.
   edge = minval(p%ens(:) - 2*bandwidths(:))
   if (bounded_below .and. (edge > lower_bound)) then
      p%bounded_below = .false.
      p%lower_bound = minval(p%ens(:) - bandwidths(:))
   else
      p%bounded_below = bounded_below
      p%lower_bound = lower_bound
   end if

   ! If the ensemble is sufficiently far from the boundary, that we can use the unbounded code, which is cheaper.
   edge = maxval(p%ens(:) + 2*bandwidths(:))
   if (bounded_above .and. (edge < upper_bound)) then
      p%bounded_above = .false.
      p%upper_bound = maxval(p%ens(:) + bandwidths(:))
   else
      p%bounded_above = bounded_above
      p%upper_bound = upper_bound
   end if

   ! Integrate across all subintervals
   p%more_params(ens_size+1) = 1._r8 ! Temporary value
   p%more_params(3*ens_size+7) = 0._r8 ! 0 at left edge
   do i=2,2*ens_size+2
      if ((p%bounded_below .and. (edges(i) <= p%lower_bound)) .or. &
          (edges(i) <= edges(i-1)) .or. &
          (p%bounded_above .and. (edges(i) > p%upper_bound))) then
         ! The subinterval [edges(i-1), edges(i)] is empty or outside the bounds
         p%more_params(3*ens_size+6+i) = p%more_params(3*ens_size+6+i-1) + 0._r8
      else
         p%more_params(3*ens_size+6+i) = p%more_params(3*ens_size+6+i-1) + gq5(edges(i-1), edges(i), p)
      end if
   end do
   p%more_params(ens_size+1) = maxval(p%more_params(3*ens_size+7:5*ens_size+8))
   p%more_params(3*ens_size+7:5*ens_size+8) = p%more_params(3*ens_size+7:5*ens_size+8) / p%more_params(ens_size+1)

end subroutine pack_kde_params

!-----------------------------------------------------------------------

function kde_cdf_params(x, p) result(quantile)
   real(r8),                       intent(in) :: x
   type(distribution_params_type), intent(in) :: p
   real(r8)                                   :: quantile

   ! Returns the cumulative distribution of the kernel density estimate
   ! at the value x.

   integer  :: i
   logical  :: use_analytical_cdf
   ! It is a waste of memory to unpack p%more_params, but it aids readability
   real(r8) :: bandwidths(p%ens_size)
   real(r8) :: y, obs_param
   integer  :: obs_dist_type
   real(r8), dimension(2*p%ens_size+2) :: edges, cdfs

   bandwidths(:) = p%more_params(1:p%ens_size)
   y             = p%more_params(p%ens_size + 2)
   obs_param     = p%more_params(p%ens_size + 3)
   obs_dist_type = nint(p%more_params(p%ens_size + 4))
   edges(:)      = p%more_params(  p%ens_size+5:3*p%ens_size+6)
   cdfs(:)       = p%more_params(3*p%ens_size+7:5*p%ens_size+8)

   quantile = 0._r8
   ! If the likelihood is uninformative and the distribution is unbounded, we can evaluate
   ! the cdf analytically (instead of using quadrature) to save computation.
   use_analytical_cdf = ( (.not.p%bounded_below) .and. (.not.p%bounded_above) .and. &
                          (obs_dist_type .eq. obs_dist_types%uninformative) )
   if (use_analytical_cdf) then
      do i=1,p%ens_size
         quantile = quantile + epanechnikov_cdf( (x - p%ens(i)) / bandwidths(i) ) / p%ens_size
      end do
   else ! Compute cdf using quadrature
      !quantile  = min(1._r8, max(0._r8, integrate_pdf(x, p)))
      if (x < edges(1)) then
         return
      end if
      do i=2,2*p%ens_size+2
         if ((edges(i-1) <= x) .and. (x <= edges(i))) then
            quantile = min(1._r8, max(0._r8, cdfs(i-1) + gq5(edges(i-1), x, p)))
            return
         else
            cycle
         end if
      end do
      quantile = 1._r8
   end if

end function kde_cdf_params

!---------------------------------------------------------------------------

function kde_cdf(x, ens, ens_size, bounded_below, bounded_above, &
   lower_bound, upper_bound, y, obs_param, obs_dist_type) result(quantile)
   real(r8),                       intent(in)  :: x
   integer,                        intent(in)  :: ens_size
   real(r8),                       intent(in)  :: ens(ens_size)
   logical,                        intent(in)  :: bounded_below,        bounded_above
   real(r8),                       intent(in)  :: lower_bound,          upper_bound
   real(r8),                       intent(in)  :: y
   real(r8),                       intent(in)  :: obs_param
   integer,                        intent(in)  :: obs_dist_type
   real(r8)                                    :: quantile

   ! Returns the cumulative distribution of the kernel density estimate
   ! at the value x.

   type(distribution_params_type)    :: p

   call pack_kde_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
      ens, y, obs_param, obs_dist_type, p)
   quantile = kde_cdf_params(x,p)

end function kde_cdf

!---------------------------------------------------------------------------

function inv_kde_cdf(quantile, ens, ens_size, bounded_below, bounded_above, &
   lower_bound, upper_bound, y, obs_param, obs_dist_type) result(x)
   real(r8),                       intent(in)  :: quantile
   integer,                        intent(in)  :: ens_size
   real(r8),                       intent(in)  :: ens(ens_size)
   logical,                        intent(in)  :: bounded_below,        bounded_above
   real(r8),                       intent(in)  :: lower_bound,          upper_bound
   real(r8),                       intent(in)  :: y
   real(r8),                       intent(in)  :: obs_param
   integer,                        intent(in)  :: obs_dist_type
   real(r8)                                    :: x

   ! Returns the value x such that cdf(x) = quantile.

   type(distribution_params_type) :: p

   call pack_kde_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
      ens, y, obs_param, obs_dist_type, p)

   x = inv_cdf(quantile, kde_cdf_params, inv_kde_first_guess_params, p)

end function inv_kde_cdf

!---------------------------------------------------------------------------

function inv_kde_cdf_params(quantile, p) result(x)
   real(r8)                                   :: x
   real(r8),                       intent(in) :: quantile
   type(distribution_params_type), intent(in) :: p

   ! Returns the value x such that cdf(x) = quantile.

   x = inv_cdf(quantile, kde_cdf_params, inv_kde_first_guess_params, p)

end function inv_kde_cdf_params

!---------------------------------------------------------------------------

function inv_kde_first_guess_params(quantile, p) result(x)
   real(r8)                                   :: x
   real(r8),                       intent(in) :: quantile
   type(distribution_params_type), intent(in) :: p
   real(r8) :: edge ! edge of support of the cdf

   ! This first-guess subroutine evaluates the cdf at the ensemble members,
   ! then finds a pair of ensemble members whose quantiles bracket the
   ! target quantile, then sets the first guess to a convex combination of
   ! these two ensemble members. If the target quantile is outside the
   ! ensemble, the first guess is a combination of the nearest ensemble
   ! member and the edge of the support of the pdf. If the
   ! target quantile is 0 or 1, the appropriate bound is returned.

   integer  :: i
   real(r8) :: dq
   real(r8), dimension(2*p%ens_size+2) :: edges, cdfs

   edges(:)      = p%more_params(  p%ens_size+5:3*p%ens_size+6)
   cdfs(:)       = p%more_params(3*p%ens_size+7:5*p%ens_size+8)

   if (quantile <= 0._r8) then
      x = minval(edges(:))
      if (p%bounded_below) then
         x = max(x, p%lower_bound)
      end if
      return
   end if

   if (quantile >= 1._r8) then
      x = maxval(edges(:))
      if (p%bounded_above) then
         x = min(x, p%upper_bound)
      end if
      return
   end if

   do i=2,2*p%ens_size+2
      if (quantile <= cdfs(i)) then
         dq = cdfs(i) - cdfs(i-1)
         x = (cdfs(i  ) - quantile) * edges(i-1) / dq &
           + (quantile - cdfs(i-1)) * edges(i  ) / dq
         return
      end if
   end do

end function inv_kde_first_guess_params

!-----------------------------------------------------------------------

subroutine separate_ensemble(ens, ens_size, bounded_below, bounded_above, &
   lower_bound, upper_bound, ens_interior, ens_size_interior, &
   q_lower, q_int, q_upper)

   ! Extracts the ensemble members that are not on the boundaries, and
   ! also returns the fraction of ensemble members on each boundary and in
   ! the interior.

   integer,   intent(in) :: ens_size
   real(r8),  intent(in) :: ens(ens_size)
   logical,   intent(in) :: bounded_below, bounded_above
   real(r8),  intent(in) :: lower_bound,   upper_bound
   real(r8), intent(out) :: ens_interior(ens_size)
   integer,  intent(out) :: ens_size_interior
   real(r8), intent(out) :: q_lower, q_int, q_upper

   ! local variables
   integer  :: i, j
   logical  :: is_interior(ens_size)
   real(r8), parameter :: eps = 0._r8 !epsilon(1._r8)

   ! Initialize
   ens_interior(:) = missing_r8

   q_lower = 0._r8
   q_int   = 0._r8
   q_upper = 0._r8
   is_interior(:) = .true.
   j = 0
   do i=1,ens_size
      if (bounded_below) then
         if (ens(i) <= lower_bound + eps) then
            q_lower        = q_lower + 1._r8
            is_interior(i) = .false.
         end if
      end if
      if (bounded_above) then
         if (ens(i) >= upper_bound - eps) then
            q_upper        = q_upper + 1._r8
            is_interior(i) = .false.
         end if
      end if
      if (is_interior(i)) then
         j = j + 1
         ens_interior(j) = ens(i)
      end if
   end do
   ens_size_interior = ens_size - nint(q_lower + q_upper)
   if (j .ne. ens_size_interior) then
      errstring = 'Total number of interior ensemble members incorrect '
      call error_handler(E_ERR, 'separate_ensemble', errstring, source)
   end if
   q_lower = q_lower / real(ens_size, r8)
   q_upper = q_upper / real(ens_size, r8)
   q_int   = 1._r8 - (q_lower + q_upper)

end subroutine separate_ensemble

!-----------------------------------------------------------------------

subroutine test_kde
   ! This routine provides limited tests of the numerics in this module. It tests
   ! the boundary correction function, the likelihood, the bandwidth selection,
   ! and the cdf inversion. It uses an ensemble [-1, 1]. It tests the bandwidth
   ! selection and the cdf inversion with zero, one, or two bounds at [-2, 2].
   ! Failing these tests suggests a serious problem. Passing them does not indicate
   ! that there are acceptable results for all possible inputs.

   integer,             parameter :: ens_size = 2
   real(r8),  dimension(ens_size) :: ensemble, bandwidths, target_bandwidths
   type(distribution_params_type) :: p
   real(r8)                       :: x, y, inv, max_diff, lx, mx, like, obs_param
   integer                        :: i, obs_dist_type

   ! Test the boundary correction code against a Mathematica implementation
   x = 0.5_r8
   call boundary_correction(x, lx, mx)
   write(*, *) '----------------------------'
   write(*, *) 'test boundary correction'
   write(*, *) 'abs difference in lx is ', abs(1.322997416020672_r8 - lx)
   write(*, *) 'abs difference in mx is ', abs(1.102497846683893_r8 - mx)
   write(*, *) 'abs differences should be less than 1e-15'

   ! Test uninformative likelihood
   max_diff = -1.0_r8
   y = 1._r8
   obs_param = 1._r8
   obs_dist_type = obs_dist_types%uninformative
   do i = 0, 1000
      x = (real(i, r8) - 500.0_r8) / 500.0_r8
      like = likelihood_function(x, y, obs_param, obs_dist_type)
      max_diff = max(abs(1._r8 - like), max_diff)
   end do
   write(*, *) '----------------------------'
   write(*, *) 'Uninformative likelihood test'
   write(*, *) 'max difference in likelihood is ', max_diff
   write(*, *) 'max difference should be less than 1e-15'

   ! Test normal likelihood
   x = 0.5_r8
   y = 0._r8
   obs_param = 1._r8
   obs_dist_type = obs_dist_types%normal
   like = likelihood_function(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Normal likelihood test'
   write(*, *) 'abs difference in likelihood is ', abs(0.8824969025845955_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test binomial obs distribution (beta likelihood)
   x = 0.25_r8
   y = 3._r8
   obs_param = 5._r8
   obs_dist_type = obs_dist_types%binomial
   like = likelihood_function(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Binomial obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(0.0087890625_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test gamma obs distribution (inverse gamma likelihood)
   y = 1._r8
   x = 2._r8
   obs_param = 3._r8
   obs_dist_type = obs_dist_types%gamma
   like = likelihood_function(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Gamma obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(0.4658368455179406_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test inverse gamma obs distribution (gamma likelihood)
   y = 3._r8
   x = 2._r8
   obs_param = 4._r8
   obs_dist_type = obs_dist_types%inv_gamma
   like = likelihood_function(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Inverse gamma obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(3.415489474106968_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test lognormal obs distribution (normal likelihood)
   y = 1._r8
   x = 2._r8
   obs_param = 4._r8
   obs_dist_type = obs_dist_types%lognormal
   like = likelihood_function(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Lognormal obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(0.6065306597126334_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test truncated-normal obs distribution
   y = 1._r8
   x = 2._r8
   obs_param = 4._r8
   obs_dist_type = obs_dist_types%truncated_normal
   like = likelihood_function(x, y, obs_param, obs_dist_type, &
         bounded_above=.true.,bounded_below=.true.,upper_bound=4._r8,lower_bound=0._r8)
   write(*, *) '----------------------------'
   write(*, *) 'Truncated-normal obs distribution test, doubly-bounded'
   write(*, *) 'abs difference in likelihood is ', abs(1.292676850528392_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'
   like = likelihood_function(x, y, obs_param, obs_dist_type, &
         bounded_above=.false.,bounded_below=.true.,upper_bound=4._r8,lower_bound=0._r8)
   write(*, *) '----------------------------'
   write(*, *) 'Truncated-normal obs distribution test, bounded below'
   write(*, *) 'abs difference in likelihood is ', abs(1.048912359301403_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'
   like = likelihood_function(x, y, obs_param, obs_dist_type, &
         bounded_above=.true.,bounded_below=.false.,upper_bound=4._r8,lower_bound=0._r8)
   write(*, *) '----------------------------'
   write(*, *) 'Truncated-normal obs distribution test, bounded above'
   write(*, *) 'abs difference in likelihood is ', abs(1.048912359301403_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test bandwidth selection
   target_bandwidths(:) = 2.462288826689833_r8 * [1._r8, 1._r8]

   ensemble(:) = [-1._r8, 1._r8]

   call get_kde_bandwidths(ens_size, ensemble, bandwidths)

   ! Compare computed bandwidths to exact
   write(*, *) '----------------------------'
   write(*, *) 'kde bandwidths test: Absolute value of differences should be less than 1e-15'
   do i = 1, ens_size
      write(*, *) i, abs(bandwidths(i) - target_bandwidths(i))
   end do

   ! Test the inversion of the cdf over the support of the pdf, unbounded
   write(*, *) '----------------------------'
   write(*, *) 'Unbounded cdf/icdf test'
   call pack_kde_params(ens_size, .false., .false., 0._r8, 0._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = (real(i - 500, r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      y = kde_cdf_params(x, p)
      inv = inv_kde_cdf_params(y, p)
      ! write(*,*) 'Quantile: ', y, 'Exact x: ', x, 'Inverted x: ', inv, 'Error: ', abs(inv - x)
      max_diff = max(abs(x-inv), max_diff)
   end do

   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'Errors should be small, max around 1E-8'
   ! The accuracy degrades near the boundary. The slope of the cdf goes to zero on the boundary
   ! But at least the second derivative is nonzero. So near the boundary it behaves like a
   ! double root, and the accuracy of the rootfinding degrades.

   call deallocate_distribution_params(p)

   ! Test the inversion of the cdf over the entire support of the pdf, bounded below
   write(*, *) '----------------------------'
   write(*, *) 'cdf/icdf test bounded below'
   call pack_kde_params(ens_size, .true., .false., -2._r8, 0._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = (real(i - 500, r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      if (x <= -2._r8) then
         cycle
      else
         y = kde_cdf_params(x, p)
         inv = inv_kde_cdf_params(y, p)
         ! write(*,*) 'Quantile: ', y, 'Exact x: ', x, 'Inverted x: ', inv, 'Error: ', abs(inv - x)
         max_diff = max(abs(x-inv), max_diff)
      end if
   end do

   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'Errors should be small, max below 1E-8'

   call deallocate_distribution_params(p)


   ! Test the inversion of the cdf over the entire support of the pdf, bounded above
   write(*, *) '----------------------------'
   write(*, *) 'cdf/icdf test bounded above'
   call pack_kde_params(ens_size, .false., .true., 0._r8, 2._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = ((real(i, r8) - 500.0_r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      if (x >= 2._r8) then
         cycle
      else
         y = kde_cdf_params(x, p)
         inv = inv_kde_cdf_params(y, p)
         ! write(*,*) 'Quantile: ', y, 'Exact x: ', x, 'Inverted x: ', inv, 'Error: ', abs(inv - x)
         max_diff = max(abs(x-inv), max_diff)
      end if
   end do

   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'Errors should be small, max below 1E-8'

   call deallocate_distribution_params(p)

   ! Test the inversion of the cdf over the entire support of the pdf, doubly bounded
   write(*, *) '----------------------------'
   write(*, *) 'cdf/icdf test doubly-bounded'
   call pack_kde_params(ens_size, .true., .true., -2._r8, 2._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = ((real(i, r8) - 500.0_r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      if ((x <= -2._r8) .or. (x >= 2._r8)) then
         cycle
      else
         y = kde_cdf_params(x, p)
         inv = inv_kde_cdf_params(y, p)
         ! write(*,*) 'Quantile: ', y, 'Exact x: ', x, 'Inverted x: ', inv, 'Error: ', abs(inv - x)
         max_diff = max(abs(x-inv), max_diff)
      end if
   end do

   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'Errors should be small, max below 1E-8'

   call deallocate_distribution_params(p)

   ! Test the quadrature: Construct a case with bounds, but where the bounds are
   ! far enough from the data that they are not used. In this case the kernel
   ! density estimate should integrate exactly to one.
   call pack_kde_params(ens_size, .true., .true., -2._r8, 2._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)
   p%lower_bound = -20._r8
   p%upper_bound =  20._r8
   p%more_params(ens_size + 1) = 1._r8

   y = integrate_pdf(0._r8, p)
   write(*, *) '----------------------------'
   write(*, *) 'test quadrature'
   write(*, *) 'abs difference is ', abs(0.5_r8 - y)
   write(*, *) 'abs difference should be less than 1e-15'

   call deallocate_distribution_params(p)

end subroutine test_kde

!------------------------------------------------------------------------

end module kde_distribution_mod
