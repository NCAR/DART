! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Random number and random sequence routines.  Can generate random draws 
!> from a uniform distribution or random draws from differently shaped
!> distributions (e.g. gaussian, gamma, exponential)

module random_seq_mod

! This module now contains both the original contents and the routines
! which used to be in a separate random_nr module.

use     types_mod, only : digits12, i8, r8
use utilities_mod, only : error_handler, E_ERR

implicit none
private

public :: random_seq_type, &
          init_random_seq, &
          random_uniform, &
          random_gaussian, &
          several_random_gaussians, & 
          twod_gaussians, &
          random_gamma, & 
          random_inverse_gamma, & 
          random_exponential, &
          ran_twist

character(len=*), parameter :: source = 'random_seq_mod.f90'

! Gives ability to generate unique repeatable sequences of random numbers
! using random congruential package. Needed to allow different assim algorithms
! that require random numbers to see identical observational sequences.

! Used to give different sequences a different but repeatable start
! There may be problems with incestuous series here; be cautious of this
! in the future.

integer :: seq_number = -1

! the following routines were transcribed from C to F90, originally
! from the GNU scientific library:  init_ran, ran_unif, ran_gauss,
! ran_gamma, ran_twist

integer, parameter :: N = 624   ! period parameters
integer, parameter :: M = 397

! hexadecimal constants
integer(i8), parameter :: UPPER_MASK  = int(z'0000000080000000', i8)
integer(i8), parameter :: LOWER_MASK  = int(z'000000007FFFFFFF', i8)
integer(i8), parameter :: FULL32_MASK = int(z'00000000FFFFFFFF', i8)
integer(i8), parameter :: magic       = int(z'000000009908B0DF', i8)
integer(i8), parameter :: C1          = int(z'000000009D2C5680', i8)
integer(i8), parameter :: C2          = int(z'00000000EFC60000', i8)

type random_seq_type
   private
   integer     :: mti
   integer(i8) :: mt(0:N-1)
   real(r8)    :: lastg
   logical     :: gset
end type random_seq_type

logical, save :: module_initialized = .false.
character(len=128) :: errstring


contains

!========================================================================
! Public entry points
!========================================================================

!------------------------------------------------------------------------

!> An integer seed can be used to get a particular repeatable sequence.

subroutine init_random_seq(r, seed)

type(random_seq_type), intent(inout) :: r
integer, optional,     intent(in)    :: seed

if ( .not. module_initialized ) call initialize_module

! Initialize the generator; use given seed if present, else sequence
if(present(seed)) then
   call init_ran(r, seed)
else
   call init_ran(r, seq_number)
   seq_number = seq_number - 1
endif

end subroutine init_random_seq

!------------------------------------------------------------------------

!> return a random draw from a uniform distribution
!> between 0.0 and 1.0

function random_uniform(r)

type(random_seq_type), intent(inout) :: r
real(r8)                             :: random_uniform

if ( .not. module_initialized ) call initialize_module

random_uniform = ran_unif(r)

end function random_uniform

!------------------------------------------------------------------------

!> return a random draw from a gaussian distribution
!> with the specified mean and standard deviation.

function random_gaussian(r, mean, standard_deviation) 

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: mean, standard_deviation
real(r8)                             :: random_gaussian

if ( .not. module_initialized ) call initialize_module

random_gaussian = ran_gauss(r) * standard_deviation + mean

end function random_gaussian

!------------------------------------------------------------------------

!> return multiple random draws from a gaussian distribution
!> with the specified mean and standard deviation.

subroutine several_random_gaussians(r, mean, standard_deviation, n, rnum)

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: mean
real(r8),              intent(in)    :: standard_deviation
integer,               intent(in)    :: n
real(r8),              intent(out)   :: rnum(n)

integer :: i

if ( .not. module_initialized ) call initialize_module

do i = 1, n
   rnum(i) = ran_gauss(r) * standard_deviation + mean
end do

end subroutine several_random_gaussians

!------------------------------------------------------------------------

!> return a random draw from a 2D multivariate gaussian distribution
!> with the specified 2 means and covariances.  

subroutine twod_gaussians(r, mean, cov, rnum)

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: mean(2)
real(r8),              intent(in)    :: cov(2, 2)
real(r8),              intent(out)   :: rnum(2)

real(r8) :: a11, a21, a22, x1, x2

if ( .not. module_initialized ) call initialize_module

! Use method from Knuth, exercise 13, section 3.4.1 to generate random
! numbers with this mean and covariance

a11 = sqrt(cov(1, 1))
a21 = cov(1, 2) / a11
a22 = sqrt(cov(2, 2) - a21**2)

! Two base independent gaussian deviates
x1 = ran_gauss(r)
x2 = ran_gauss(r)

! Use these to generate correlated
rnum(1) = mean(1) + a11 * x1
rnum(2) = mean(2) + a21 * x1 + a22 * x2

end subroutine twod_gaussians

!------------------------------------------------------------------------

!> return a random draw from a gamma distribution 
!> with the specified shape and scale parameter, often
!> denoted with the symbols kappa and theta.
!>
!> note that there are multiple common formulations of
!> the second parameter to this function.  if you have
!> 'rate' instead of 'scale' specify:  1.0 / rate  
!> for the scale parameter.  shape and rate are often 
!> denoted with the symbols alpha and beta.
!>
!> The distribution mean should be: shape * rate, 
!> and the variance: shape * (rate^2).


function random_gamma(r, rshape, rscale)

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: rshape 
real(r8),              intent(in)    :: rscale 
real(r8)                             :: random_gamma

if ( .not. module_initialized ) call initialize_module

if (rshape <= 0.0_r8) then
   write(errstring, *) 'Shape parameter must be positive, was ', rshape
   call error_handler(E_ERR, 'random_gamma', errstring, source)
endif

if (rscale <= 0.0_r8) then
   write(errstring, *) 'Scale parameter (scale=1/rate) must be positive, was ', rscale
   call error_handler(E_ERR, 'random_gamma', errstring, source)
endif

! internal routine uses rshape, rscale
random_gamma = ran_gamma(r, rshape, rscale)

end function random_gamma

!------------------------------------------------------------------------

!> return a random draw from an inverse gamma distribution 
!> with the specified shape and scale parameter.
!> NOTE THIS IS DIFFERENT FROM THE RANDOM GAMMA FORUMULATION.
!> gamma uses a rate parameter, inverse gamma uses scale.
!> if you have 'rate' instead of 'scale', specify 1.0 / rate

function random_inverse_gamma(r, rshape, rscale)

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: rshape
real(r8),              intent(in)    :: rscale 
real(r8)                             :: random_inverse_gamma

real(r8) :: g

if ( .not. module_initialized ) call initialize_module

if (rshape <= 0.0_r8) then
   write(errstring, *) 'Shape parameter must be positive, was ', rshape
   call error_handler(E_ERR, 'random_inverse_gamma', errstring, source)
endif

if (rscale <= 0.0_r8) then
   write(errstring, *) 'Scale parameter (scale=1/rate) must be positive, was ', rscale
   call error_handler(E_ERR, 'random_inverse_gamma', errstring, source)
endif

g = ran_gamma(r, rshape, rscale)

! ran_gamma can't return 0 so its safe to divide by it
random_inverse_gamma = 1.0_r8/g

end function random_inverse_gamma

!------------------------------------------------------------------------

!> return a random draw from an exponential distribution with
!> the specified rate parameter lambda.  if you have a scale parameter,
!> also called the mean, standard deviation, or survival parameter,
!> specify 1.0 / scale for rate.

function random_exponential(r, rate)

type(random_seq_type), intent(inout) :: r
real(r8),              intent(in)    :: rate 
real(r8)                             :: random_exponential

if ( .not. module_initialized ) call initialize_module

if (rate <= 0.0_r8) then
   write(errstring, *) 'Rate parameter (rate=1/scale) must be positive, was ', rate
   call error_handler(E_ERR, 'random_exponential', errstring, source)
endif

! internal routine uses shape, scale
random_exponential = ran_gamma(r, 1.0_r8, 1.0_r8/rate)

end function random_exponential

!-------------------------------------------------------------------

!========================================================================
! Private, internal routines below here
!========================================================================

!-------------------------------------------------------------------

subroutine initialize_module

if ( .not. module_initialized ) then
   module_initialized = .true.
endif

end subroutine initialize_module

!-------------------------------------------------------------------

!> A random congruential random number generator of the Mersenne Twister type
!> See GNU Scientific Library code.

subroutine init_ran(s, seed)
 integer, intent(in) :: seed
 type(random_seq_type), intent(out) :: s

 integer :: i, sd

! Initialize the generator for use with
! repeatable sequences

sd = seed
if (sd == 0) sd = 4357   ! do not allow seed to be 0, use default

s%mt(0) = iand(int(sd,i8), FULL32_MASK)

! See Knuth's "Art of Computer Programming" Vol. 2, 3rd Ed. p.106 
! for multiplier.

do i=1, N-1
   s%mt(i) = 1812433253_i8 * ieor(s%mt(i-1), ishft(s%mt(i-1), -30_i8)) + i

   s%mt(i) = iand(s%mt(i), FULL32_MASK)
end do

s%mti = N
s%lastg = 0.0_r8
s%gset = .false.

end subroutine init_ran

!-----------------------------------------------------------------

!> A random congruential random number generator from
!> the GNU Scientific Library (The Mersenne Twister MT19937 varient.)

function ran_unif(s)
 type(random_seq_type), intent(inout) :: s
 real(r8) :: ran_unif

integer :: kk
integer(i8) :: k, y, m1

! original c code:
!  define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

if (s%mti >= N) then
   ! generate N words at a time
   do kk = 0, N-M-1
      y = ior(iand(s%mt(kk), UPPER_MASK), iand(s%mt(kk + 1), LOWER_MASK))
      if (iand(y,1_i8) == 1_i8) then
         m1 = magic
      else
         m1 = 0_i8
      endif
      s%mt(kk) = ieor(s%mt(kk + M), ieor(ishft(y,-1_i8), m1))

! original c code:
!          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
!          mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);

   enddo

   do kk = N-M, N-2
      y = ior(iand(s%mt(kk), UPPER_MASK), iand(s%mt(kk + 1), LOWER_MASK))
      if (iand(y,1_i8) == 1_i8) then
         m1 = magic
      else
         m1 = 0_i8
      endif
      s%mt(kk) = ieor(s%mt(kk + (M-N)), ieor(ishft(y,-1_i8), m1))

! original c code:
!          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
!          mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);

   enddo

   y = ior(iand(s%mt(N-1), UPPER_MASK), iand(s%mt(0), LOWER_MASK))
   if (iand(y,1_i8) == 1_i8) then
      m1 = magic
   else
      m1 = 0_i8
   endif
   s%mt(N-1) = ieor(s%mt(M-1), ieor(ishft(y,-1_i8), m1))

! original c code:
!        unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
!        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);

   s%mti = 0
endif

! Tempering

k = s%mt(s%mti)

k = ieor(k, ishft(k, -11_i8))
k = ieor(k, iand(ishft(k, 7_i8),  C1))
k = ieor(k, iand(ishft(k, 15_i8), C2))
k = ieor(k, ishft(k, -18_i8))

! original c code:
!  k ^= (k >> 11);
!  k ^= (k << 7) & 0x9d2c5680UL;
!  k ^= (k << 15) & 0xefc60000UL;
!  k ^= (k >> 18);

s%mti = s%mti + 1

! at this point we have an integer value for k
! this routine returns 0.0 <= real < 1.0, so do
! the divide here.  return range:  [0,1).

ran_unif = real(real(k, digits12) / 4294967296.0_digits12, r8)

end function ran_unif

!------------------------------------------------------------------------

!> Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122 
!> Returns a N(0, 1) random number draw from a gaussian distribution

function ran_gauss(s)

type(random_seq_type), intent(inout) :: s
real(r8) :: ran_gauss

real(digits12) :: x, y, r2, t
integer :: lc

! each pass through this code generates 2 values.  if we have one
! saved from a previous call, return it without doing any new work.
! otherwise, generate 2 new values; return one and save the other.

if (s%gset) then
   ran_gauss = s%lastg
   s%gset = .false.
else

   lc = 0
10 continue
   lc = lc + 1
   
   ! choose x,y in uniform square (-1,-1) to (+1,+1) 
   x = -1.0_digits12 + 2.0_digits12 * ran_unif(s)
   y = -1.0_digits12 + 2.0_digits12 * ran_unif(s)
   
   ! repeat if it is outside the unit circle or at origin
   r2 = x*x + y*y
   if (r2 >= 1.0_digits12 .or. r2 == 0.0_digits12) then
      if (lc > 100) then
         write(errstring, *) 'x, y = ', x, y
         call error_handler(E_ERR, 'ran_gauss', &
                            'if both x and y are -1, random number generator probably not initialized', &
                            source, text2 = errstring);
      endif
      goto 10
   endif

   t = sqrt(-2.0_digits12 * log(r2) / r2)
   s%lastg   = real(x * t, r8)
   s%gset    = .true.
   ran_gauss = real(y * t, r8)
endif

end function ran_gauss

!-------------------------------------------------------------------

!> Gamma(a,b) generator from the GNU Scientific library.
!>
!> There are three different parameterizations in common use:
!>  1. With a shape parameter k and a scale parameter θ.
!>  2. With a shape parameter α = k and an inverse scale parameter β = 1/θ, called a rate parameter.
!>  3. With a shape parameter k and a mean parameter μ = k/β.
!>
!> warning!
!> this is the private, internal interface and it uses the first
!> parameterization, kappa and theta, like the original GSL code.
!> our public interfaces use the words 'shape', 'scale' and 'rate'
!> to try to be more clear about which parameterization we're using.
!>
!> this internal routine assumes the caller has already verified that
!> the shape and scale are positive.

recursive function ran_gamma (r, rshape, rscale) result(g)

type(random_seq_type), intent(inout) :: r
real(r8), intent(in) :: rshape
real(r8), intent(in) :: rscale
real(r8) :: g

real(r8) :: x, v, u
real(r8) :: c, d


! we chose the public interface to take 'rate' since that 
! seems to be more common in the kalman world.  but the original
! code is based on a scale parameter, so that's what we use here.

if (rshape < 1.0) then
   u = random_uniform (r)
   g = ran_gamma (r, 1.0_r8 + rshape, rscale) * u ** (1.0_r8 / rshape)
   return
endif

d = rshape - 1.0_r8 / 3.0_r8
c = (1.0_r8 / 3.0_r8) / sqrt (d)
v = -1.0_r8

OUTER: do while (.true.)
    v = -1.0_r8
    do while (v <= 0.0_r8)
       x = random_gaussian(r, 0.0_r8, 1.0_r8)
       v = 1.0_r8 + c * x
    enddo

    v = v * v * v;
    u = random_uniform (r);

    if (u < 1 - 0.0331_r8 * x * x * x * x) exit OUTER

    if (log (u) < 0.5_r8 * x * x + d * (1 - v + log (v))) exit OUTER

enddo OUTER

g = rscale * d * v

end function ran_gamma

!------------------------------------------------------------------------

!> A random congruential random number generator from
!> the GNU Scientific Library (The Mersenne Twister MT19937 varient.)
!> This routine returns an Integer.

function ran_twist(s)
 type(random_seq_type), intent(inout) :: s
 
integer :: kk
integer(i8) :: ran_twist, k, y, m1

! original c code:
!  define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

if (s%mti >= N) then
   ! generate N words at a time
   do kk = 0, N-M-1
      y = ior(iand(s%mt(kk), UPPER_MASK), iand(s%mt(kk + 1), LOWER_MASK))
      if (iand(y,1_i8) == 1_i8) then
         m1 = magic
      else
         m1 = 0_i8
      endif
      s%mt(kk) = ieor(s%mt(kk + M), ieor(ishft(y,-1_i8), m1))

! original c code:
!          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
!          mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);

   enddo

   do kk = N-M, N-2
      y = ior(iand(s%mt(kk), UPPER_MASK), iand(s%mt(kk + 1), LOWER_MASK))
      if (iand(y,1_i8) == 1_i8) then
         m1 = magic
      else
         m1 = 0_i8
      endif
      s%mt(kk) = ieor(s%mt(kk + (M-N)), ieor(ishft(y,-1_i8), m1))

! original c code:
!          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
!          mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);

   enddo

   y = ior(iand(s%mt(N-1), UPPER_MASK), iand(s%mt(0), LOWER_MASK))
   if (iand(y,1_i8) == 1_i8) then
      m1 = magic
   else
      m1 = 0_i8
   endif
   s%mt(N-1) = ieor(s%mt(M-1), ieor(ishft(y,-1_i8), m1))

! original c code:
!        unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
!        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);

   s%mti = 0
endif

! Tempering

k = s%mt(s%mti)

k = ieor(k, ishft(k, -11_i8))
k = ieor(k, iand(ishft(k, 7_i8),  C1))
k = ieor(k, iand(ishft(k, 15_i8), C2))
k = ieor(k, ishft(k, -18_i8))

! original c code:
!  k ^= (k >> 11);
!  k ^= (k << 7) & 0x9d2c5680UL;
!  k ^= (k << 15) & 0xefc60000UL;
!  k ^= (k >> 18);

s%mti = s%mti + 1

! at this point we have an integer value for k
! this routine returns 0.0 <= real < 1.0, so do
! the divide here.  return range:  [0,1).

ran_twist = k

end function ran_twist

!------------------------------------------------------------------------


end module random_seq_mod

