module chisq_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Computes chi square statistics and significance, from numerical recipes

implicit none

private
public chsone

integer, parameter :: itmax = 100
double precision, parameter :: eps = 3.e-7
double precision, parameter :: fpmin = 1.e-30

double precision, parameter :: stp = 2.5066282746310005D0
double precision :: cof(6) = (/76.18009172947146D0, -86.50532032941677D0, &
   24.01409824083091D0, -1.231739572450155D0, &
   .1208650973866179D-2, -.5395239384953D-5/)

contains

!==============================================================================

!  CHISQUARE TEST MATERIAL FROM NUMERICAL RECIPES AND ASSOCIATED SUPPORTING
!  SUBROUTINES

subroutine chsone(bins, ebins, knstrn, df, chsq, prob)

implicit none

integer, intent(in) :: knstrn
double precision, intent(in) :: bins(:), ebins(:)
double precision, intent(out) :: chsq, df, prob

!  GIVEN THE ARRAY BINS(1:NBINS) CONTAINING THE OBSERVED NUMBER OF EVENTS,
!  AND AN ARRAY EBINS(1:NBINS) CONTAINING THE EXPECTED NUMBERS OF EVENTS,
!  AND GIVEN THE NUMBER OF CONSTRAINTS KNSTRN (NORMALLY ONE), THIS ROUTINE
!  RETURNS (TRIVIALLY) THE NUMBER OF DEGREES OF FREEDOM DF, AND (NONTRIV)
!  THE CHI-SQUARE CHSQ AND THE SIGNIFICANCE PROB.  A SMALL VALUE OF PROB
!  INDICATES A SIGNIFICANT DIFFERENCE BETWEEN THE DISTRIBUTIONS BINS AND
!  EBINS.  NOTE THAT BINS AND EBINS ARE BOTH REAL ARRAYS, ALTHOUGH BINS WILL
!  NORMALLY CONTAIN INTEGER VALUES.

INTEGER J

df = size(bins) - knstrn
chsq = 0.
do j = 1, size(bins)
   if(ebins(j) < 0.) then
      write(*, *) 'FATAL ERROR: bad expected number in bins in chsone '
      stop
   endif
   chsq = chsq + (bins(j) - ebins(j))**2 / ebins(j)
end do 
prob = gammq(0.5*df, 0.5*chsq)

end subroutine chsone

!==========================================================================

function gammq(a, x)

implicit none

double precision :: gammq
double precision, intent(in) :: a, x

!  RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A, X) = 1 - P(A, X)

double precision :: gammcf, gamser, gln

if(x < 0. .or. a < 0.) then
   write(*, *) 'FATAL ERROR: Bad arguments to gammq in chisq computation'
   stop
endif

if(x < a+1.) then
   call gser(gamser, a, x, gln)
   gammq = 1. - gamser
else
   call gcf(gammcf, a, x, gln)
   gammq = gammcf
endif

end function gammq

!==========================================================================

subroutine gser(gamser, a, x, gln)

implicit none

double precision, intent(in) :: a, x
double precision, intent(out) ::  gamser, gln

!  RETURNS THE INCOMPLETE GAMMA FUNCTION P(A, X) EVALUATED BY ITS SERIES
!  REPRESENTATION AS GAMSER.  ALSO RETURN IN(GAMMA(A)) AS GLN.

integer :: n
double precision :: ap, del, sum

gln = gammln(a)
if(x < 0.0) then
   write(*, *) 'ERROR: x < 0 in gser in computing chisq'
   stop
endif

if(X == 0.) then
   gamser = 0.
else
   ap = a
   sum = 1. / a
   del = sum
   do n = 1, itmax
      ap = ap + 1.
      del = del * x / ap
      sum = sum + del
      if(abs(del) < abs(sum) * eps) goto 1
   end do
   write(*, *) 'ERROR: a too large, itmax too small in gser in computing chisq'
1  gamser = sum * exp(-x + a * log(x) - gln)
endif

end subroutine gser

!===========================================================================

subroutine gcf(gammcf, a, x, gln)

implicit none

double precision, intent(in) :: a, x
double precision, intent(out) :: gammcf, gln

!  RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A, X) EVALUATED BY ITS CONTINUED 
!  FRACTION REPRESENTATION AS GAMMCF.  ALSO RETURNS IN(GAMMA(A)) AS GLN.
!  PARAMETERS: ITMAX IS THE MAXIMUM ALLOWED NUMBER OF ITERATIONS; EPS IS
!  THE RELATIVE ACCURACY; FPMIN IS A NUMBER NEAR THE SMALLEST REPRESENTABLE
!  FLOATING-POINT NUMBER.

integer :: i
double precision :: an, b, c, d, del, h

gln = gammln(a)
b = x + 1. - a
c = 1. / fpmin
d = 1. / b
h = d
do i = 1, itmax
   an = -i * (i - a)
   b = b + 2.
   d = an * d + b
   if(abs(d) < fpmin) d = fpmin
   c = b + an / c
   if(abs(c) < fpmin) c = fpmin
   d = 1. / d
   del = d * c
   h = h * del
   if(abs(del - 1.) < eps) goto 1
end do
write(*, *) 'ERROR: a too large, itmax too small in gcf'
stop

1  gammcf = exp(-x + a * log(x) - gln) * h

end subroutine gcf

!============================================================================

function gammln(xx)

implicit none

double precision gammln
double precision, intent(in) :: xx

!  returns the value  in(gamma(xx)) for xx > 0

integer j
double precision ser, tmp, x, y

!  INTERNAL ARITHMETIC WILL BE DONE IN DOUBLE PRECISION, A NICETY THAT YOU
!  CAN OMIT IF FIVE-FIGURE ACCURACY IS GOOD ENOUGH.

x = xx
y = x
tmp = x + 5.5d0
tmp = (x + 0.5d0) * log(tmp) - tmp
ser = 1.000000000190015d0
do j = 1, 6
   y = y + 1.d0
   ser = ser + cof(j) / y
end do
gammln = tmp + log(stp*ser/x)

end function gammln

!=============================================================================

end module chisq_mod
