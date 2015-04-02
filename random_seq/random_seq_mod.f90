! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module random_seq_mod

! This module now contains both the original contents and the routines
! which used to be in a separate random_nr module.

use     types_mod, only : digits12, i8, r8
use utilities_mod, only : register_module, error_handler, E_ERR

implicit none
private

public :: random_seq_type, init_random_seq, random_gaussian, &
   several_random_gaussians, random_uniform, twod_gaussians

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Gives ability to generate unique repeatable sequences of random numbers
! using random congruential package. Needed to allow different assim algorithms
! that require random numbers to see identical observational sequences.

! Used to give different sequences a different but repeatable start
! There may be problems with incestuous series here; be cautious of this
! in the future.

integer :: seq_number = -1

! the following routines were transcribed from C to F90, originally
! from the GNU scientific library:  init_ran, ran_unif, ran_gauss

integer, parameter :: N = 624   ! period parameters
integer, parameter :: M = 397

! hexadecimal constants
integer(i8), parameter :: UPPER_MASK  = z'0000000080000000'
integer(i8), parameter :: LOWER_MASK  = z'000000007FFFFFFF'
integer(i8), parameter :: FULL32_MASK = z'00000000FFFFFFFF'
integer(i8), parameter :: magic       = z'000000009908B0DF'
integer(i8), parameter :: C1          = z'000000009D2C5680'
integer(i8), parameter :: C2          = z'00000000EFC60000'

type random_seq_type
   private
   integer     :: mti
   integer(i8) :: mt(0:N-1)   ! FIXME: move this to 1:N once its working
   real(r8)    :: lastg
   logical     :: gset
end type random_seq_type

logical, save :: module_initialized = .false.
character(len=128) :: errstring


contains



!========================================================================

subroutine init_random_seq(r, seed)
!----------------------------------------------------------------------
! An integer seed can be used to get a particular repeatable sequence.

!
implicit none

type(random_seq_type), intent(inout) :: r
integer, optional,     intent(in)    :: seed

! Initialize the generator; use given seed if present, else sequence
if(present(seed)) then
   call init_ran(r, seed)
else
   call init_ran(r, seq_number)
   seq_number = seq_number - 1
endif

end subroutine init_random_seq

!========================================================================

function random_uniform(r)

implicit none

type(random_seq_type), intent(inout) :: r
real(r8) :: random_uniform

random_uniform = ran_unif(r)

end function random_uniform

!========================================================================

function random_gaussian(r, mean, standard_deviation) 

implicit none

type(random_seq_type), intent(inout) :: r
real(r8), intent(in) :: mean, standard_deviation
real(r8) :: random_gaussian

random_gaussian = ran_gauss(r) * standard_deviation + mean

end function random_gaussian

!========================================================================

subroutine several_random_gaussians(r, mean, standard_deviation, n, rnum)

implicit none

type(random_seq_type), intent(inout) :: r
real(r8), intent(in) :: mean, standard_deviation
integer, intent(in) :: n
real(r8), intent(out) :: rnum(n)

integer :: i

do i = 1, n
   rnum(i) = ran_gauss(r) * standard_deviation + mean
end do

end subroutine several_random_gaussians

!========================================================================

subroutine twod_gaussians(r, mean, cov, rnum)

implicit none

type(random_seq_type), intent(inout) :: r
real(r8), intent(in) :: mean(2), cov(2, 2)
real(r8), intent(out) :: rnum(2)

real(r8) :: a11, a21, a22, x1, x2

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

!========================================================================


!-------------------------------------------------------------------

subroutine initialize_module

if ( .not. module_initialized ) then
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

end subroutine initialize_module

!-------------------------------------------------------------------

! A random congruential random number generator of the Mersenne Twister type
! See GNU Scientific Library code.

subroutine init_ran(s, seed)
 integer, intent(in) :: seed
 type(random_seq_type), intent(out) :: s

 integer :: i, sd

if ( .not. module_initialized ) call initialize_module

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

!  A random congruential random number generator from
! the GNU Scientific Library (The Mersenne Twister MT19937 varient.)

function ran_unif(s)
 type(random_seq_type), intent(inout) :: s
 real(r8) :: ran_unif

integer :: kk
integer(i8) :: k, y, m1

if ( .not. module_initialized ) call initialize_module


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

function ran_gauss(s)

! Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122 

! Returns a N(-1, 1) random number draw from a gaussian distribution

type(random_seq_type), intent(inout) :: s
real(r8) :: ran_gauss

real(digits12) :: x, y, r2, t
integer :: lc

if ( .not. module_initialized ) call initialize_module

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
                            source, revision, revdate, &
                            text2 = errstring);
      endif
      goto 10
   endif

   t = sqrt(-2.0_digits12 * log(r2) / r2)
   s%lastg   = real(x * t, r8)
   s%gset    = .true.
   ran_gauss = real(y * t, r8)
endif

end function ran_gauss

!------------------------------------------------------------------------


end module random_seq_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
