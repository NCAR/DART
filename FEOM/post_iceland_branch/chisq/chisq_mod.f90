! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module chisq_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 
!
! Computes chi square statistics and significance


use     types_mod, only : r8
use utilities_mod, only : error_handler, E_ERR

implicit none
private

public :: chsone

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

integer, parameter :: itmax = 100
real(r8), parameter :: eps = 3.e-7
real(r8), parameter :: fpmin = 1.e-30

real(r8), parameter :: stp = 2.5066282746310005_r8
real(r8) :: cof(6) = (/76.18009172947146_r8, -86.50532032941677_r8, &
                       24.01409824083091_r8, -1.231739572450155_r8, &
                     0.01208650973866179_r8, -0.000005395239384953_r8/)

contains


!  CHISQUARE TEST AND SUPPORTING SUBROUTINES

  subroutine chsone(bins, ebins, knstrn, df, chsq, prob)
!==============================================================================
! subroutine chsone(bins, ebins, knstrn, df, chsq, prob)

implicit none

integer,  intent(in)  :: knstrn
real(r8), intent(in)  :: bins(:), ebins(:)
real(r8), intent(out) :: chsq, df, prob

!  GIVEN THE ARRAY BINS(1:NBINS) CONTAINING THE OBSERVED NUMBER OF EVENTS,
!  AND AN ARRAY EBINS(1:NBINS) CONTAINING THE EXPECTED NUMBERS OF EVENTS,
!  AND GIVEN THE NUMBER OF CONSTRAINTS KNSTRN (NORMALLY ONE), THIS ROUTINE
!  RETURNS (TRIVIALLY) THE NUMBER OF DEGREES OF FREEDOM DF, AND (NONTRIV)
!  THE CHI-SQUARE CHSQ AND THE SIGNIFICANCE PROB.  A SMALL VALUE OF PROB
!  INDICATES A SIGNIFICANT DIFFERENCE BETWEEN THE DISTRIBUTIONS BINS AND
!  EBINS.  NOTE THAT BINS AND EBINS ARE BOTH REAL ARRAYS, ALTHOUGH BINS WILL
!  NORMALLY CONTAIN INTEGER VALUES.

integer :: j
character(len=129) :: errstring

df = size(bins) - knstrn
chsq = 0.0_r8
do j = 1, size(bins)
   if( ebins(j) < 0.0_r8 ) then
      write(errstring,*)'bad expected number (',ebins(j),') in bin ',j
      call error_handler(E_ERR,'chsone', errstring, source, revision, revdate)
   endif
   chsq = chsq + (bins(j) - ebins(j))**2 / ebins(j)
end do 
prob = gammq(0.5*df, 0.5*chsq)

end subroutine chsone



  function gammq(a, x)
!==========================================================================
! function gammq(a, x)
!
!  RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A, X) = 1 - P(A, X)

implicit none

real(r8), intent(in) :: a, x
real(r8)             :: gammq

real(r8) :: gammcf, gamser, gln

character(len=129) :: errstring

if( x < 0.0_r8 .or. a < 0.0_r8) then
   write(errstring,*)'x ',x,' or a ',a,' < 0.0'
   call error_handler(E_ERR,'gammq', errstring, source, revision, revdate)
endif

if(x < a+1.0_r8) then
   call gser(gamser, a, x, gln)
   gammq = 1.0_r8 - gamser
else
   call gcf(gammcf, a, x, gln)
   gammq = gammcf
endif

end function gammq


  subroutine gser(gamser, a, x, gln)
!==========================================================================
! subroutine gser(gamser, a, x, gln)
!
! RETURNS THE INCOMPLETE GAMMA FUNCTION P(A, X) EVALUATED BY ITS SERIES
! REPRESENTATION AS GAMSER.  ALSO RETURN IN(GAMMA(A)) AS GLN.

implicit none

real(r8), intent(in)  :: a, x
real(r8), intent(out) ::  gamser, gln

integer  :: n
real(r8) :: ap, del, sum
character(len=129) :: errstring

gln = gammln(a)
if(x < 0.0_r8) then
   write(errstring,*)'x ',x,' < 0.0'
   call error_handler(E_ERR,'gser', errstring, source, revision, revdate)
endif

if(x == 0._r8) then
   gamser = 0.0_r8
else
   ap = a
   sum = 1.0_r8 / a
   del = sum
   do n = 1, itmax
      ap = ap + 1.0_r8
      del = del * x / ap
      sum = sum + del
      if(abs(del) < abs(sum) * eps) goto 1
   end do
   write(*, *) 'ERROR: a too large, itmax too small in gser in computing chisq'
1  gamser = sum * exp(-x + a * log(x) - gln)
endif

end subroutine gser



  subroutine gcf(gammcf, a, x, gln)
!===========================================================================
! subroutine gcf(gammcf, a, x, gln)
!
! RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A, X) EVALUATED BY ITS CONTINUED 
! FRACTION REPRESENTATION AS GAMMCF.  ALSO RETURNS IN(GAMMA(A)) AS GLN.
! PARAMETERS: ITMAX IS THE MAXIMUM ALLOWED NUMBER OF ITERATIONS; EPS IS
! THE RELATIVE ACCURACY; FPMIN IS A NUMBER NEAR THE SMALLEST REPRESENTABLE
! FLOATING-POINT NUMBER.

implicit none

real(r8), intent(in)  :: a, x
real(r8), intent(out) :: gammcf, gln

integer  :: i
real(r8) :: an, b, c, d, del, h

gln = gammln(a)
b   = x + 1.0_r8 - a
c   = 1.0_r8 / fpmin
d   = 1.0_r8 / b
h   = d

do i = 1, itmax

   an = -i * (i - a)
   b  = b + 2.0_r8
   d  = an * d + b

   if(abs(d) < fpmin) d = fpmin

   c = b + an / c

   if(abs(c) < fpmin) c = fpmin

   d   = 1.0_r8 / d
   del = d * c
   h   = h * del

   if(abs(del - 1.) < eps) goto 1

end do

   call error_handler(E_ERR,'gcf', 'a too large, itmax too small', source, revision, revdate)

1  gammcf = exp(-x + a * log(x) - gln) * h

end subroutine gcf


  function gammln(xx)
!============================================================================
! function gammln(xx)
!
! returns the value  in(gamma(xx)) for xx > 0
!
! internal arithmetic will be done in double precision, a nicety that you
! can omit if five-figure accuracy is good enough.

implicit none

real(r8), intent(in) :: xx
real(r8)             :: gammln

integer  :: j
real(r8) :: ser, tmp, x, y

x   = xx
y   = x
tmp = x + 5.5_r8
tmp = (x + 0.5_r8) * log(tmp) - tmp
ser = 1.000000000190015_r8

do j = 1, 6
   y   = y + 1.0_r8
   ser = ser + cof(j) / y
end do

gammln = tmp + log(stp*ser/x)

end function gammln

!=============================================================================
! end of chisq_mod.f90
!=============================================================================

end module chisq_mod
