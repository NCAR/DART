module rootfinding_mod

use types_mod,               only : r8

use utilities_mod,           only : E_ERR, E_MSG, error_handler

use distribution_params_mod, only : distribution_params_type, deallocate_distribution_params

implicit none
private

public :: inv_cdf

character(len=512)          :: errstring
character(len=*), parameter :: source = 'rootfinding_mod.f90'

contains

!------------------------------------------------------------------------

function inv_cdf(quantile_in, cdf, first_guess, p) result(x)

   interface
      function cdf(x, p)
         use types_mod, only : r8
         use distribution_params_mod, only : distribution_params_type
         real(r8)                                   :: cdf
         real(r8), intent(in)                       :: x
         type(distribution_params_type), intent(in) :: p
      end function
   end interface

   interface
      function first_guess(quantile, p)
         use types_mod, only : r8
         use distribution_params_mod, only : distribution_params_type
         real(r8)                                   :: first_guess
         real(r8), intent(in)                       :: quantile
         type(distribution_params_type), intent(in) :: p
      end function
   end interface

   real(r8)                                   :: x
   real(r8), intent(in)                       :: quantile_in
   type(distribution_params_type), intent(in) :: p

   ! Given a quantile q, finds the value of x for which cdf(x) is approximately q
   ! Starts with Newton's method using an approximate gradient until either convergence
   ! or overshoot. If it overshoots then it switches to the bisection-like ITP method
   ! from Oliveira & Takahashi, ACM Trans Math Soft, 2020. Here f(x) = cdf(x) - quantile

   ! Local variables:
   integer,  parameter :: max_iterations = 50 ! Limit on the total number of iterations.
   real(r8), parameter :: min_quantile = 0.0_r8,  max_quantile = 0.999999999999999_r8
   real(r8), parameter :: tol = epsilon(1._r8)**(0.75_r8) ! Absolute on q, relative on x
   real(r8) :: quantile, q_err
   real(r8) :: delta_x, delta_f
   real(r8) :: x_guess, q_guess, x0, x1, f0, f1
   real(r8) :: a, b, fa, fb
   integer  :: iter

   real(r8) :: lower_bound,   upper_bound
   logical  :: bounded_below, bounded_above

   ! Extract the required information from the p type
   bounded_below = p%bounded_below;   bounded_above = p%bounded_above
   lower_bound   = p%lower_bound;     upper_bound   = p%upper_bound

   quantile = quantile_in

   ! Do a test for illegal values on the quantile
   if(quantile <  0.0_r8 .or. quantile > 1.0_r8) then
      write(errstring, *) 'Illegal Quantile input', quantile
      call error_handler(E_ERR, 'inv_cdf', errstring, source)
   endif

   ! If the distribution is bounded, quantiles at the limits have values at the bounds
   if(bounded_below .and. quantile == 0.0_r8) then
      x = lower_bound
      return
   endif
   if(bounded_above .and. quantile == 1.0_r8) then
      x = upper_bound
      return
   endif

   ! If input quantiles are outside the numerically supported range, move them to the extremes
   quantile = min(quantile, max_quantile)
   ! code tests stably for many distributions with min_quantile of 0.0, could remove this
   quantile = max(quantile, min_quantile)

   ! Get first guess from functional approximation
   x_guess = first_guess(quantile, p)

   ! Evaluate the cdf
   q_guess = cdf(x_guess, p)

   ! Check if first guess is close enough
   q_err = quantile - q_guess
   if (abs(q_err) .le. tol) then
      x = x_guess
      return
   end if

   ! If bounded below and target quantile is between 0 and q_guess, go to ITP
   if (bounded_below .and. (quantile .lt. q_guess)) then
      a  = lower_bound
      b  = x_guess
      fa = 0._r8 - quantile
      fb = q_guess - quantile
      x  = inv_cdf_ITP(cdf, quantile, a, b, fa, fb, max_iterations, p)
      return
   end if

   ! If bounded above and target quantile is between q_guess and 1, go to ITP
   if (bounded_above .and. (q_guess .lt. quantile)) then
      a  = x_guess
      b  = upper_bound
      fa = q_guess
      fb = 1._r8
      x  = inv_cdf_ITP(cdf, quantile, a, b, fa, fb, max_iterations, p)
      return
   end if

   ! If we reach this line, then we can't yet bracket the true value of x, so we start
   ! doing steps of the secant method, either until convergence or until we bracket
   ! the root. If we bracket the root then we finish using the ITP method.
   x0 = x_guess
   f0 = q_guess - quantile
   delta_x = max(1e-2_r8, 1e-2_r8 * abs(x_guess))
   if (f0 .gt. 0._r8) then
      delta_x = -delta_x
   end if
   do iter=1,7
      x1 = x_guess + delta_x
      if(bounded_below .and. (x1 .lt. lower_bound)) x1 = lower_bound
      if(bounded_above .and. (x1 .gt. upper_bound)) x1 = upper_bound
      f1 = cdf(x1, p) - quantile
      delta_f = abs(f1 - f0)
      if (delta_f .le. 0._r8) then
         delta_x = 2._r8 * delta_x
      else
         exit
      end if
   end do
   do iter = 1, max_iterations
      ! If f0 * f1 < 0 then x0 and x1 bracket the root, so go to ITP
      if (f0 * f1 .lt. 0._r8 ) then
         a  = min(x0, x1)
         b  = max(x0, x1)
         fa = min(f0, f1)
         fb = max(f0, f1)
         x = inv_cdf_ITP(cdf, quantile, a, b, fa, fb, max_iterations - iter + 1, p)
         return
      else ! Do a secant step
         delta_f = f1 - f0
         if (abs(delta_f) .eq. 0._r8) then ! Stop. Failed step: Flat CDF
            write(errstring, *)  'Failed to converge; flat cdf for quantile ', quantile
            write(*,*) p%bounded_below, p%bounded_above
            write(*,*) 'More params:', p%more_params(:)
            write(*,*) 'iter, x0, x1, f0, f1:', iter, x0, x1, f0, f1
            write(*,*) 'ens:', p%ens(:)
            call error_handler(E_MSG, 'rootfinding_mod:inv_cdf', errstring, source)
            x = x1
            return
         else
            x_guess = x1 - f1 * delta_x / delta_f
            if (bounded_above .and. (x_guess .gt. upper_bound)) x_guess = upper_bound
            if (bounded_below .and. (x_guess .lt. lower_bound)) x_guess = lower_bound
            x0 = x1
            f0 = f1
            x1 = x_guess
            delta_x = x1 - x0
            f1 = cdf(x1, p) - quantile
            x = x1
         end if
      end if

      ! Finished secant step. Check for convergence.
      if (abs(delta_x) .le. tol * max(abs(x0), abs(x1)) .or. (abs(f1) .le. tol)) then
         return
      endif
   end do

   ! For now, have switched a failed convergence to return the latest guess
   ! This has implications for stability of probit algorithms that require further study
   ! Not currently happening for any of the test cases on gfortran
   write(errstring, *)  'Failed to converge for quantile ', quantile
   call error_handler(E_MSG, 'inv_cdf', errstring, source)
   !!!call error_handler(E_ERR, 'inv_cdf', errstring, source)

end function inv_cdf

!-----------------------------------------------------------------------

function inv_cdf_ITP(cdf, quantile, a, b, fa, fb, max_iterations, p) result(x)
   interface
      function cdf(x, p)
         use types_mod, only : r8
         use distribution_params_mod, only : distribution_params_type
         real(r8)                                   :: cdf
         real(r8), intent(in)                       :: x
         type(distribution_params_type), intent(in) :: p
      end function
   end interface

   real(r8),                       intent(in) :: quantile
   real(r8),                    intent(inout) :: a, b, fa, fb
   integer,                        intent(in) :: max_iterations
   type(distribution_params_type), intent(in) :: p
   real(r8)                                   :: x

   ! Given a quantile q, finds the value of x for which cdf(x) is approximately q
   ! Uses the bisection-like ITP method from Oliveira & Takahashi, ACM Trans Math
   ! Soft, 2020. Here f(x) = cdf(x) - quantile. Assumes f(a) < 0 and f(b) > 0 on
   ! input.

   ! Local variables:
   integer,  parameter :: n0 = 1._r8
   real(r8), parameter :: kappa2 = 2._r8
   real(r8), parameter :: phi = 0.5_r8 * (1._r8 + sqrt(5._r8)) ! Golden ratio
   real(r8), parameter :: tol = epsilon(1._r8)**(0.75_r8) ! abs on f, rel on x
   real(r8) :: kappa1
   real(r8) :: delta
   real(r8) :: x_t, x_f, x_half, f_ITP
   real(r8) :: eps, r
   integer  :: i, n_max, n_half

   kappa1 = 0.2_r8 / (b - a)
   eps = tol * max(abs(a), abs(b))
   ! Max number of iterations is either (i) input max, or (ii) number of bisection
   ! iterations (plus n0) needed to achieve the desired relative tolerance on x.
   n_half = max(0, ceiling(log(0.5_r8 * (b - a) / eps) / log(2._r8)))
   n_max = min(max_iterations, n0 + n_half)

   do i=1,n_max
      x_half = 0.5_r8 * (a + b) ! Bisection guess
      x_f    = (b* fa - a * fb) / (fa - fb) ! Secant/Regula-Falsi guess
      delta  = kappa1 * abs(b - a)**kappa2
      if (delta .le. abs(x_half - x_f)) then
         x_t = x_f + sign(delta, x_half - x_f)
      else
         x_t = x_half
      end if
      r = max(0._r8, eps * 2**(n_half - i) - 0.5_r8 * (b - a))
      if (abs(x_t - x_half) .le. r) then
         x = x_t
      else
         x = x_half - sign(r, x_half - x_f)
      end if
      f_ITP = cdf(x, p) - quantile
      if (abs(f_ITP) .le. tol) then
         return
      elseif (f_ITP .gt. 0._r8) then
         b  = x
         fb = f_ITP
      else
         a  = x
         fa = f_ITP
      end if
   end do

end function inv_cdf_ITP

!-----------------------------------------------------------------------

end module rootfinding_mod
