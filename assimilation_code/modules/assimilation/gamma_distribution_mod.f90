! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module gamma_distribution_mod

use types_mod, only : r8, PI

use utilities_mod, only : E_ERR, error_handler

use normal_distribution_mod, only : norm_cdf

use random_seq_mod, only : random_seq_type, random_uniform

implicit none
private

public :: gamma_pdf, gamma_cdf, inv_gamma_cdf, random_gamma, test_gamma

character(len=512)     :: errstring
character(len=*), parameter :: source = 'gamma_distribution_mod.f90'

contains

!-----------------------------------------------------------------------

subroutine test_gamma

real(r8) :: x, y, inv
real(r8) :: mean, variance, sd, shape, scale
integer :: i

! Input a mean and variance
mean = 10.0_r8
sd = 1.0_r8
variance = sd**2

! Get shape and scale
shape = mean**2 / variance
scale = variance / mean

! Confirm by going backwards
write(*, *) 'comp mean ', shape * scale
write(*, *) 'comp= sd ', sqrt(shape * scale**2)

write(*, *) 'shape scale ', shape, scale

do i = 0, 1000
   x = mean + ((i - 500.0_r8) / 500.0_r8) * 5.0_r8 * sd
   y = gamma_cdf(x, shape, scale)
   write(*, *) i, x, y
   inv = inv_gamma_cdf(y, shape, scale)
   write(*, *) i, inv
   write(34, *) x, inv, x - inv, y
end do

end subroutine test_gamma

!-----------------------------------------------------------------------

function inv_gamma_cdf(q, shape, scale)

real(r8)             :: inv_gamma_cdf
real(r8), intent(in) :: q
real(r8), intent(in) :: shape
real(r8), intent(in) :: scale

! Given a quantile q, finds the value of x for which the gamma cdf
! with shape and scale has approximately this quantile

! This version uses a Newton method using the fact that the PDF is the derivative of the CDF

real(r8) :: x_guess, old_x_guess, q_guess, dq_dx, del_x, q_err, q_err_new, q_new, x_new
real(r8) :: mn, sd
integer  :: i, j

! Limit on the total iterations; There is no deep thought behind this choice
integer, parameter :: max_iterations = 100
! Limit on number of times to halve the increment; again, no deep thought
integer, parameter :: max_half_iterations = 10

! Unclear what error tolerance is needed for DA applications; 
! A smaller value seems to be possible but leads to more iterations
real(r8), parameter :: xtol = 1.0e-12_r8

! Do a special test for exactly 0
if(q == 0.0_r8) then
   inv_gamma_cdf = 0.0_r8
   return
endif

! Need some sort of first guess, should be smarter here
! For starters, take the mean for this shape and scale
sd = sqrt(shape * scale**2)
mn = shape * scale
x_guess = mn + (q - 0.5_r8) * 6.0_r8 * sd

! If the guess is below zero, just default back to the mean
if(x_guess < 0.0_r8) x_guess = mn
old_x_guess = x_guess

do i = 1, max_iterations
   old_x_guess = x_guess
   q_guess = gamma_cdf(x_guess, shape, scale)
   dq_dx = gamma_pdf(x_guess, shape, scale)
   q_err = q - q_guess
   del_x = q_err / dq_dx
   x_new = x_guess + del_x

   q_new = gamma_cdf(x_new, shape, scale)
   q_err_new = q_new - q

   do j = 1, max_half_iterations
      if(abs(q_err_new) > abs(q_err)) then
         del_x = del_x / 2.0_r8
         x_new = x_guess + del_x
         q_new = gamma_cdf(x_new, shape, scale)
         q_err_new = q_new - q
      else
         ! Inefficient to be in the loop for this
         exit
      endif
   end do

   x_guess = x_new
  
   ! Check for stopping criterion
   if(abs(old_x_guess - x_guess) <= xtol) then
      inv_gamma_cdf = x_guess
      return
   else
      old_x_guess = x_guess
   endif 

enddo

! Fell off the end, should be an error return eventually?
errstring = 'Failed to converge '
call error_handler(E_ERR, 'inv_gamma_cdf', errstring, source)
stop

end function inv_gamma_cdf

!---------------------------------------------------------------------------

function gamma_pdf(x, shape, scale)

! Returns the probability density of a gamma function with shape and scale
! at the value x

! Returns a large negative value if called with illegal values

real(r8) :: gamma_pdf
real(r8), intent(in) :: x, shape, scale

! All inputs must be nonnegative
if(x < 0.0_r8 .or. shape < 0.0_r8 .or. scale < 0.0_r8) then
   gamma_pdf = -99.9_r8
elseif(x == 0.0_r8) then
   gamma_pdf = 0.0_r8
else
   gamma_pdf = x**(shape - 1.0_r8) * exp(-x / scale) / &
      (gamma(shape) * scale**shape)
endif

end function gamma_pdf

!---------------------------------------------------------------------------

function gamma_cdf(x, shape, scale)

! Returns the cumulative distribution of a gamma function with shape and scale
! at the value x

! Returns a large negative value if called with illegal values

real(r8) :: gamma_cdf
real(r8), intent(in) :: x, shape, scale

! All inputs must be nonnegative
if(x < 0.0_r8 .or. shape < 0.0_r8 .or. scale < 0.0_r8) then
   gamma_cdf = -99.9_r8
else
   ! Use definition as incomplete gamma ratio to gamma
   gamma_cdf = gammad(x / scale, shape)
endif

end function gamma_cdf

!---------------------------------------------------------------------------
function gammad (x, p)

implicit none 

real(r8)             :: gammad
real(r8), intent(in) :: x
real(r8), intent(in) :: p

!*****************************************************************************80
!
!! GAMMAD computes the Incomplete Gamma Integral
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by B Shea.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    B Shea,
!    Algorithm AS 239:
!    Chi-squared and Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete 
!    gamma ratio.  0 <= X, and 0 < P.
!
!
!    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete 
!    Gamma integral.
!

real(r8)            :: a, b, c, an, arg, pn(6), rn
real(r8), parameter :: elimit = - 88.0_r8
real(r8), parameter :: oflo   = 1.0e+37_r8
real(r8), parameter :: plimit = 1000.0_r8
real(r8), parameter :: tol    = 1.0e-14_r8
real(r8), parameter :: xbig   = 1.0e+08_r8

! x zero returns zero
if(x == 0.0_r8) then
   gammad = 0.0_r8
elseif(xbig < x) then
   !  If X is large set GAMMAD = 1.
   gammad = 1.0_r8
elseif(plimit < p) then
!  If P is large, use a normal approximation.
   pn(1) = 3.0_r8 * sqrt(p) * ((x / p)**(1.0_r8 / 3.0_r8) + &
       1.0_r8 / (9.0_r8 * p) - 1.0_r8)

   gammad = norm_cdf(pn(1), 0.0_r8, 1.0_r8)
elseif(x <= 1.0_r8 .or. x < p) then
!  Use Pearson's series expansion.
!  Original note: (Note that P is not large enough to force overflow in logAM).
   arg = p * log(x) - x - log(gamma(p + 1.0_r8))
   c = 1.0_r8
   gammad = 1.0_r8
   a = p

   do
      a = a + 1.0_r8
      c = c * x / a
      gammad = gammad + c
      if(c <= tol) exit
   end do

   arg = arg + log(gammad)

   if(elimit <= arg) then
     gammad = exp(arg)
   else
      gammad = 0.0_r8
    end if
else 
   !  Use a continued fraction expansion.
   arg = p * log(x) - x - log(gamma(p))
   a = 1.0_r8 - p
   b = a + x + 1.0_r8
   c = 0.0_r8
   pn(1) = 1.0_r8
   pn(2) = x
   pn(3) = x + 1.0_r8
   pn(4) = x * b
   gammad = pn(3) / pn(4)

   do
      a = a + 1.0_r8
      b = b + 2.0_r8
      c = c + 1.0_r8
      an = a * c
      pn(5) = b * pn(3) - an * pn(1)
      pn(6) = b * pn(4) - an * pn(2)

      if (pn(6) /= 0.0_r8) then
         rn = pn(5) / pn(6)
         if(abs(gammad - rn) <= min(tol, tol * rn)) exit
         gammad = rn
      end if

      pn(1) = pn(3)
      pn(2) = pn(4)
      pn(3) = pn(5)
      pn(4) = pn(6)

      ! Re-scale terms in continued fraction if terms are large.
      if (oflo <= abs(pn(5))) pn(1:4) = pn(1:4) / oflo

   end do

   arg = arg + log(gammad)

   if (elimit <= arg) then
      gammad = 1.0_r8 - exp(arg)
   else
      gammad = 1.0_r8
    endif
endif

end function gammad

!---------------------------------------------------------------------------

function random_gamma(r, rshape, rscale)

! Note that this provides same qualitative functionality as a similarly named
! routine in the random_seq_mod that uses a rejection algorithm. However, once
! we have an inverse cdf function for a distribution, it is possible to generate
! random numbers by first getting a draw from a U(0, 1) and then inverting these
! quantiles to get an actual value

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: rshape
real(r8),              intent(in)    :: rscale
real(r8)                             :: random_gamma

real(r8) :: quantile
if (rshape <= 0.0_r8) then
   write(errstring, *) 'Shape parameter must be positive, was ', rshape
   call error_handler(E_ERR, 'random_gamma', errstring, source)
endif

if (rscale <= 0.0_r8) then
   write(errstring, *) 'Scale parameter (scale=1/rate) must be positive, was ', rscale
   call error_handler(E_ERR, 'random_gamma', errstring, source)
endif

! Draw from U(0, 1) to get a quantile
quantile = random_uniform(r)
! Invert cdf to get a draw from gamma
random_gamma = inv_gamma_cdf(quantile, rshape, rscale)

end function random_gamma

!---------------------------------------------------------------------------

end module gamma_distribution_mod
