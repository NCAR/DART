! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! Thanks to Chris Riedel who developed the methods in this module.

module beta_distribution_mod

use types_mod, only : r8, PI

use utilities_mod, only : E_ERR, error_handler

use random_seq_mod, only : random_seq_type, random_uniform

implicit none
private

public :: beta_pdf, beta_cdf, inv_beta_cdf, random_beta, test_beta

character(len=512)     :: errstring
character(len=*), parameter :: source = 'beta_distribution_mod.f90'

contains

!-----------------------------------------------------------------------

subroutine test_beta

real(r8) :: x, y, p, inv
real(r8) :: mean, variance, sd, alpha, beta
integer :: i

! Input a mean and variance
!!!mean = 10.0_r8
!!!sd = 1.0_r8
!!!variance = sd**2

! Get alpha and beta
!!!shape = mean**2 / variance
!!!scale = variance / mean

! Confirm by going backwards
!!!write(*, *) 'comp mean ', shape * scale
!!!write(*, *) 'comp= sd ', sqrt(shape * scale**2)

alpha = 5.0_r8
beta = 2.0_r8

do i = 0, 1000
   x = i / 1000.0_r8
   p = beta_pdf(x, alpha, beta)
   y = beta_cdf(x, alpha, beta)
   write(33, *) x, y, p
   inv = inv_beta_cdf(y, alpha, beta)
   write(34, *) x, inv, x - inv, y
end do

end subroutine test_beta

!-----------------------------------------------------------------------

function inv_beta_cdf(quantile, alpha, beta)

real(r8)             :: inv_beta_cdf
real(r8), intent(in) :: quantile
real(r8), intent(in) :: alpha
real(r8), intent(in) :: beta

! Given a quantile q, finds the value of x for which the beta cdf
! with alpha and beta has approximately this quantile

integer, parameter :: max_iter = 100
! For beta tests, this loop almost never happens so 25 seems very larg
integer, parameter :: max_half_iterations = 25

real(r8) :: reltol, dq_dx
real(r8) :: x_guess, q_guess, x_new, q_new, del_x, del_q, del_q_old
integer :: iter, j

if (alpha <= 0.0_r8 .or. beta <= 0.0_r8) then
  errstring = 'Negative input beta parameters'
  call error_handler(E_ERR, 'inv_beta_cdf', errstring, source)
endif

if (quantile < 0.0_r8 .or. quantile > 1.0_r8) then
  errstring = 'Bad input quantile value'
  call error_handler(E_ERR, 'inv_beta_cdf', errstring, source)
endif

if (quantile == 0.0_r8) then
    inv_beta_cdf= 0.0_r8
else if (quantile == 1.0_r8) then
    inv_beta_cdf= 1.0_r8
else  
   !Using Newton's Method to find a root of beta_cdf(x, alpha, beta) = quantile
   ! Start with the mean for this alpha and beta as a first guess
   ! Could use information about quantile to refine this and reduce required iterations
   x_guess = alpha/(alpha + beta)
   ! Make sure that the guess isn't too close to 1 or 0 where things can get ugly
   reltol = (EPSILON(x_guess))**(3./4.)
   ! Use information from quantile to refine first guess
   x_guess = max(reltol, min(1.0_r8-reltol, x_guess))

   ! Evaluate the cd 
   q_guess = beta_cdf(x_guess, alpha, beta)
  
   del_q = q_guess - quantile

   ! Iterations of the Newton method to approximate the root
   do iter= 1, max_iter
      ! The PDF is the derivative of the CDF
      dq_dx = beta_pdf(x_guess, alpha, beta)
      ! Linear approximation for how far to move in x
      del_x = del_q / dq_dx 

      ! Avoid moving too much of the fraction towards the bounds at 0 and 1 
      ! because of potential instability there. The factor of 10.0 here is a magic number
      x_new = max(x_guess/10.0_r8, min(1.0_r8 - (1.0_r8 - x_guess)/10.0_r8, x_guess-del_x))

      ! Look for convergence; If the change in x is smaller than approximate precision 
      if (abs(del_x) <= reltol*x_guess) then
         inv_beta_cdf= x_new
         return
      endif
   
      ! If we've gone too far, the new error will be bigger than the old; 
      ! Repeatedly half the distance until this is rectified 
      del_q_old = del_q
      q_new = beta_cdf(x_new, alpha, beta)
      do j = 1, max_half_iterations
         del_q = q_new - quantile
         if (abs(del_q) < abs(del_q_old)) then
            EXIT
         endif
         x_new = (x_guess + x_new)/2.0_r8
         q_new = beta_cdf(x_new, alpha, beta)
      end do

      x_guess = x_new
   end do
   !!!inv_beta_cdf= x_new

   ! Fell off the end, should be an error return eventually?
   errstring = 'Failed to converge '
   call error_handler(E_ERR, 'inv_beta_cdf', errstring, source)

endif

end function inv_beta_cdf

!---------------------------------------------------------------------------

function beta_pdf(x, alpha, beta)

! Returns the probability density of a beta function with alpha and beta
! at the value x

! Returns a large negative value if called with illegal values

real(r8) :: beta_pdf
real(r8), intent(in) :: x, alpha, beta

real(r8) :: gamma_ratio

! Parameters alpha and beta must be positive
if(alpha <= 0.0_r8 .or. beta <= 0.0_r8) then
   beta_pdf = -99.9_r8
elseif(x < 0.0 .or. x > 1.0_r8) then
   beta_pdf = -99.9_r8
elseif(alpha == 1.0_r8 .and. x == 0.0_r8) then
   ! Tricky stuff for x = 0 or 1; 
    beta_pdf = beta
elseif(beta == 1.0_r8 .and. x == 1.0_r8) then
    beta_pdf = alpha
elseif(alpha < 1.0_r8 .and. x == 0.0_r8) then
    beta_pdf = -99.9_r8
elseif(beta < 1.0_r8 .and. x == 1.0_r8) then
    beta_pdf = -99.9_r8
else
    ! Use definition via gammas since this is a F90 intrinsic
    ! Is this numerically robust?
    gamma_ratio = gamma(alpha) * gamma(beta) / gamma(alpha + beta)
    beta_pdf = x**(alpha - 1.0_r8) * (1.0_r8 - x)**(beta - 1.0_r8) / gamma_ratio
endif

end function beta_pdf

!---------------------------------------------------------------------------

function beta_cdf(x, alpha, beta)

! Returns the cumulative distribution of a beta function with alpha and beta
! at the value x

! Returns a large negative value if called with illegal values

real(r8) :: beta_cdf
real(r8), intent(in) :: x, alpha, beta

! Parameters must be positive
if(alpha <= 0.0_r8 .or. beta <= 0.0_r8) then
   beta_cdf = -99.9_r8
elseif(x < 0.0_r8 .or. x > 1.0_r8) then
   ! x must be in 0 1
   beta_cdf = -99.9_r8
elseif(x == 0.0_r8) then 
   beta_cdf = 0.0_r8
elseif(x == 1.0_r8) then
   beta_cdf = 1.0_r8
elseif (x > (alpha + 1.0_r8)/(alpha + beta + 2.0_r8)) then
    beta_cdf = (1.0_r8 - incbeta(beta, alpha, 1.0_r8 - x))
else
    beta_cdf = incbeta(alpha, beta, x)
endif

end function beta_cdf

!---------------------------------------------------------------------------

function random_beta(r, alpha, beta)

! Note that this provides same qualitative functionality as a similarly named
! routine in the random_seq_mod that uses a rejection algorithm. However, once
! we have an inverse cdf function for a distribution, it is possible to generate
! random numbers by first getting a draw from a U(0, 1) and then inverting these
! quantiles to get an actual value

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: alpha
real(r8),              intent(in)    :: beta
real(r8)                             :: random_beta

real(r8) :: quantile
if (alpha <= 0.0_r8) then
   write(errstring, *) 'Alpha parameter must be positive, was ', alpha
   call error_handler(E_ERR, 'random_beta', errstring, source)
endif

if (beta <= 0.0_r8) then
   write(errstring, *) 'Beta parameter must be positive, was ', beta
   call error_handler(E_ERR, 'random_beta', errstring, source)
endif

! Draw from U(0, 1) to get a quantile
quantile = random_uniform(r)
! Invert cdf to get a draw from beta
random_beta = inv_beta_cdf(quantile, alpha, beta)

end function random_beta

!---------------------------------------------------------------------------

function incbeta(a,b,x)
  real(r8), intent(in) :: a,b,x
  real(r8) :: incbeta
  real(r8), parameter :: TINY = 1.0e-30
  real(r8), parameter :: STOP = 1.0e-8
  real(r8) :: lbeta_ab,front,f,c,d,numerator,cd
  integer :: m,i
  integer, parameter :: bot = 2

  if (x < 0 .or. x > 1) then
    errstring = 'Input value for x is not between 0 - 1'
    call error_handler(E_ERR, 'incbeta', errstring, source)
  endif

  call betaln(a,b,lbeta_ab)    
  front = exp(log(x)*a + log(1.0_r8-x)*b - lbeta_ab) / a
  f = 1.0_r8;c=1.0_r8;d=0.0_r8
  do i=0, 200
    m = floor(i/2.0_r8)
    if (i == 0) then
      numerator = 1.0_r8     
    else if (mod(i,2) == 0) then
      numerator = (m*(b-m)*x)/((a+2.0_r8*m-1.0)*(a+2.0_r8*m))
    else
      numerator = -((a+m)*(a+b+m)*x)/((a+2.0_r8*m)*(a+2.0_r8*m+1))
    end if   

    d = 1.0_r8 + (numerator * d)
    if (abs(d) < TINY) d = TINY
    d = 1.0_r8/d

    c = 1.0_r8 + (numerator/c)
    if (abs(c) < TINY) c = TINY

    cd = c*d
    
    f = cd*f

    if (abs(1.0_r8 - cd) < STOP) then
      incbeta = front * (f-1.0_r8)
      return
    end if
  end do
  errstring = 'Alg. did not converge'
  call error_handler(E_ERR, 'incbeta', errstring, source)
end function incbeta

!---------------------------------------------------------------------------

subroutine betaln(a,b,output)
  real(r8), intent(in) :: a,b
  real(r8), intent(out) :: output

  output = gammal(a) + gammal(b) - gammal(a+b)

end subroutine betaln

!---------------------------------------------------------------------------

function gammal(xx)
  real(r8), intent(in) :: xx
  real(r8) :: gammal
  real(r8), parameter :: cov(6) = (/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, -.5395239384953d-5/)
  real(r8), parameter :: stp = 2.5066282746310005d0
  real(r8) :: x,y,tmp,ser
  integer :: j

  x = xx
  y = x

  tmp = x+5.5d0
  tmp = (x+0.5d0)*log(tmp) - tmp
  ser = 1.000000000190015d0
  do j=1, 6
    y = y + 1.0_r8
    ser = ser + cov(j)/y  
  end do 
  gammal = tmp + log(stp*ser/x)
end function gammal

!---------------------------------------------------------------------------

!---------------------------------------------------------------------------


end module beta_distribution_mod
