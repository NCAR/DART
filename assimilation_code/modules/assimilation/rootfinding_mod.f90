! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! This module originally written by I. Grooms and described (briefly)
! in "A Quantile-Conserving Ensemble Filter Based on Kernel-Density Estimation"
! by Grooms & Riedel (Remote Sensing, 2024; DOI:10.3390/rs16132377).

module rootfinding_mod

use types_mod,               only : r8

use utilities_mod,           only : E_ERR, E_MSG, E_ALLMSG, error_handler

use distribution_params_mod, only : distribution_params_type, deallocate_distribution_params

implicit none
private

public :: inv_cdf

character(len=512)          :: errstring
character(len=*), parameter :: source = 'rootfinding_mod.f90'

contains

!------------------------------------------------------------------------

function inv_cdf(cdf_in, cdf, first_guess, p) result(quantile)

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
      function first_guess(u, p)
         use types_mod, only : r8
         use distribution_params_mod, only : distribution_params_type
         real(r8)                                   :: first_guess
         real(r8), intent(in)                       :: u
         type(distribution_params_type), intent(in) :: p
      end function
   end interface

   real(r8)                                   :: quantile
   real(r8), intent(in)                       :: cdf_in
   type(distribution_params_type), intent(in) :: p

   ! Given a probability cdf_in, finds the value of x (i.e. the quantile) for which cdf(x)
   ! is approximately cdf_in = u.
   ! Starts with Newton's method using an approximate derivative until either convergence
   ! or overshoot. If it overshoots then it switches to the bisection-like ITP method
   ! from Oliveira & Takahashi, ACM Trans Math Soft, 2020. Here f(x) = cdf(x) - cdf_in

   ! Local variables:
   integer,  parameter :: MAX_ITERATIONS = 50 ! Limit on the total number of iterations.
   real(r8), parameter :: MIN_PROBABILITY = 0.0_r8,  MAX_PROBABILITY = 0.999999999999999_r8
   real(r8), parameter :: TOL = epsilon(1._r8)**(0.75_r8) ! Absolute on q, relative on x
   real(r8) :: u, u_err
   real(r8) :: delta_x, delta_f
   real(r8) :: x_guess, u_guess, x0, x1, f0, f1
   real(r8) :: a, b, fa, fb
   integer  :: iter

   real(r8) :: lower_bound,   upper_bound
   logical  :: bounded_below, bounded_above

   ! Extract the required information from the p type
   bounded_below = p%bounded_below;   bounded_above = p%bounded_above
   lower_bound   = p%lower_bound;     upper_bound   = p%upper_bound

   u = cdf_in

   ! Do a test for illegal values on the probability u
   if(u <  0.0_r8 .or. u > 1.0_r8) then
      write(errstring, *) 'Illegal cdf input', u
      call error_handler(E_ERR, 'inv_cdf', errstring, source)
   endif

   ! If the distribution is bounded, probabilities of 0 and 1 have values at the bounds
   if(bounded_below .and. u == 0.0_r8) then
      quantile = lower_bound
      return
   endif
   if(bounded_above .and. u == 1.0_r8) then
      quantile = upper_bound
      return
   endif

   ! If input probabilities are outside the numerically supported range, move them to the extremes
   u = min(u, MAX_PROBABILITY)
   ! code tests stably for many distributions with min_probability of 0.0, could remove this
   u = max(u, MIN_PROBABILITY)

   ! Get first guess from functional approximation
   x_guess = first_guess(u, p)

   ! Evaluate the cdf
   u_guess = cdf(x_guess, p)

   ! Check if first guess is close enough
   u_err = u - u_guess
   if (abs(u_err) .le. TOL) then
      quantile = x_guess
      return
   end if

!    ! If bounded below and target probability is between 0 and u_guess, go to ITP
!    ! This is inefficient when the ensemble is far from the lower bound
!    if (bounded_below .and. (u .lt. u_guess)) then
!       a  = lower_bound
!       b  = x_guess
!       fa = 0._r8 - u
!       fb = u_guess - u
!       quantile = inv_cdf_ITP(cdf, u, a, b, fa, fb, MAX_ITERATIONS, p)
!       return
!    end if
!
!    ! If bounded above and target probability is between u_guess and 1, go to ITP
!    ! This is inefficient when the ensemble is far from the upper bound
!    if (bounded_above .and. (u_guess .lt. u)) then
!       a  = x_guess
!       b  = upper_bound
!       fa = u_guess
!       fb = 1._r8
!       quantile = inv_cdf_ITP(cdf, u, a, b, fa, fb, MAX_ITERATIONS, p)
!       return
!    end if

   ! If we reach this line, then we can't yet bracket the true value of x, so we start
   ! doing steps of the secant method, either until convergence or until we bracket
   ! the root. If we bracket the root then we finish using the ITP method.
   x0 = x_guess
   f0 = u_guess - u
   delta_x = max(1e-2_r8, 1e-2_r8 * abs(x_guess))
   if (f0 .gt. 0._r8) then
      delta_x = -delta_x
   end if
   do iter=1,7
      x1 = x_guess + delta_x
      if(bounded_below .and. (x1 .lt. lower_bound)) x1 = lower_bound
      if(bounded_above .and. (x1 .gt. upper_bound)) x1 = upper_bound
      f1 = cdf(x1, p) - u
      delta_f = abs(f1 - f0)
      if (delta_f .le. 0._r8) then
         delta_x = 2._r8 * delta_x
      else
         exit
      end if
   end do
   do iter = 1, MAX_ITERATIONS
      ! If f0 * f1 < 0 then x0 and x1 bracket the root, so go to ITP
      if (f0 * f1 .lt. 0._r8 ) then
         a  = min(x0, x1)
         b  = max(x0, x1)
         fa = min(f0, f1)
         fb = max(f0, f1)
         quantile = inv_cdf_ITP(cdf, u, a, b, fa, fb, MAX_ITERATIONS - iter + 1, p)
         return
      else ! Do a secant step
         delta_f = f1 - f0
         if (abs(delta_f) .eq. 0._r8) then ! Stop. Failed step: Flat CDF
            write(errstring, *)  'Failed to converge; flat cdf'
            write(*,*) p%bounded_below, p%bounded_above
            write(*,*) 'More params:', p%more_params(:)
            write(*,*) 'iter, x0, x1, f0, f1:', iter, x0, x1, f0, f1
            write(*,*) 'ens:', p%ens(:)
            call error_handler(E_MSG, 'inv_cdf:', errstring, source)
            quantile = x1
            return
         else
            x_guess = x1 - f1 * delta_x / delta_f
            if (bounded_above .and. (x_guess .gt. upper_bound)) x_guess = upper_bound
            if (bounded_below .and. (x_guess .lt. lower_bound)) x_guess = lower_bound
            x0 = x1
            f0 = f1
            x1 = x_guess
            delta_x = x1 - x0
            f1 = cdf(x1, p) - u
            quantile = x1
         end if
      end if

      ! Finished secant step. Check for convergence.
      if (abs(delta_x) .le. TOL * max(abs(x0), abs(x1)) .or. (abs(f1) .le. TOL)) then
         return
      endif
   end do

   ! For now, have switched a failed convergence to return the latest guess
   ! This has implications for stability of probit algorithms that require further study
   ! Not currently happening for any of the test cases on gfortran
   write(errstring, *)  'Failed to converge for probability ', u
   call error_handler(E_ALLMSG, 'inv_cdf', errstring, source)
   !!!call error_handler(E_ERR, 'inv_cdf', errstring, source)

end function inv_cdf

!-----------------------------------------------------------------------

function inv_cdf_ITP(cdf, u, a, b, fa, fb, max_iterations, p) result(x)
   interface
      function cdf(x, p)
         use types_mod, only : r8
         use distribution_params_mod, only : distribution_params_type
         real(r8)                                   :: cdf
         real(r8), intent(in)                       :: x
         type(distribution_params_type), intent(in) :: p
      end function
   end interface

   real(r8),                       intent(in) :: u
   real(r8),                    intent(inout) :: a, b, fa, fb
   integer,                        intent(in) :: max_iterations
   type(distribution_params_type), intent(in) :: p
   real(r8)                                   :: x

   ! Given a probability u, finds the value of x for which cdf(x) is approximately u
   ! Uses the bisection-like ITP method from Oliveira & Takahashi, ACM Trans Math
   ! Soft, 2020. Here f(x) = cdf(x) - u. Assumes f(a) < 0 and f(b) > 0 on input.

   ! Local variables:
   integer,  parameter :: N0 = 1._r8
   real(r8), parameter :: KAPPA2 = 2._r8
   real(r8), parameter :: TOL = epsilon(1._r8)**(0.75_r8) ! abs on f, rel on x
   real(r8) :: kappa1
   real(r8) :: delta
   real(r8) :: x_t, x_f, x_half, f_ITP
   real(r8) :: eps, r
   real(r8), save :: ln_2 = log(2._r8)
   integer  :: i, n_max, n_half

   kappa1 = 0.2_r8 / (b - a)
   eps = TOL * max(abs(a), abs(b))
   ! Max number of iterations is either (i) input max, or (ii) number of bisection
   ! iterations (plus N0) needed to achieve the desired relative tolerance on x.
   n_half = max(0, ceiling(log(0.5_r8 * (b - a) / eps) / ln_2))
   n_max = min(max_iterations, N0 + n_half)

   do i=0,n_max-1
      x_half = 0.5_r8 * (a + b) ! Bisection guess
      x_f    = (b* fa - a * fb) / (fa - fb) ! Secant/Regula-Falsi guess
      delta  = kappa1 * abs(b - a)**KAPPA2
      if (delta .le. abs(x_half - x_f)) then
         x_t = x_f + sign(delta, x_half - x_f)
      else
         x_t = x_half
      end if
      r = max(0._r8, eps * 2**(n_max - i) - 0.5_r8 * (b - a))
      if (abs(x_t - x_half) .le. r) then
         x = x_t
      else
         x = x_half - sign(r, x_half - x_f)
      end if
      f_ITP = cdf(x, p) - u
      if (abs(f_ITP) .le. TOL) then
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
