! module random_numerical_recipes_mod
module random_nr_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

implicit none

private
public random_seq_type, init_ran1, ran1, gasdev

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer, parameter :: m1 = 259200, ia1 = 7141, ic1 = 54773
integer, parameter :: m2 = 134456, ia2 = 8121, ic2 = 28411
integer, parameter :: m3 = 243000, ia3 = 4561, ic3 = 51349
double precision, parameter :: rm1 = 1./m1, rm2 = 1./m2

type random_seq_type
   private
   integer :: ix1, ix2, ix3, iset
   double precision :: r(97), gset
end type random_seq_type

contains

!-------------------------------------------------------------------

!  package or random number generators from numerical recipes
subroutine init_ran1(s, idum)

implicit none

integer, intent(in) :: idum
type(random_seq_type), intent(out) :: s
integer iff, j

! Initialize the numerical recipes ran1 generator for use with
! repeatable sequences

s%ix1 = mod(ic1 - idum, m1)
s%ix1 = mod(ia1*s%ix1 + ic1, m1)
s%ix2 = mod(s%ix1, m2)
s%ix1 = mod(ia1*s%ix1 + ic1, m1)
s%ix3 = mod(s%ix1, m3)
do j = 1, 97
   s%ix1 = mod(ia1*s%ix1 + ic1, m1)
   s%ix2 = mod(ia2*s%ix2 + ic2, m2)
   s%r(j) = (dble(s%ix1) + dble(s%ix2)*rm2)*rm1
end do

! Initialize the value needed for Gaussian efficiency
s%iset = 0

end subroutine init_ran1

!-----------------------------------------------------------------

!  package or random number generators from numerical recipes
function ran1(s)

implicit none

type(random_seq_type), intent(inout) :: s
double precision :: ran1

integer :: j

!  returns a uniform random deviate between 0.0 and 1.0.

s%ix1 = mod(ia1*s%ix1 + ic1, m1)
s%ix2 = mod(ia2*s%ix2 + ic2, m2)
s%ix3 = mod(ia3*s%ix3 + ic3, m3)
j = 1 + (97*s%ix3) / m3
if(j > 97 .or. j < 1) then
   write(*, *) 'Fatal error in random number generator ran1'
   stop
endif
ran1 = s%r(j)
s%r(j) = (dble(s%ix1) + dble(s%ix2)*rm2)*rm1
return
end function ran1

!---------------------------------------------------------------------

function gasdev(s)
!  returns a normally distributed deviate with zero mean and unit 
!  variance using ran1 as source of uniform deviates.

implicit none

type(random_seq_type), intent(inout) :: s
double precision :: gasdev

double precision :: v1, v2, r, fac

if(s%iset == 0) then
10 v1 = 2. * ran1(s) - 1.
   v2 = 2. * ran1(s) - 1.
   r = v1**2 + v2**2
   if(r >= 1.) goto 10
   fac = sqrt(-2. * log(r) / r)
   s%gset = v1 * fac
   gasdev = v2 * fac
   s%iset = 1
else
   gasdev = s%gset
   s%iset = 0
endif

end function gasdev

!------------------------------------------------------------------------

end module random_nr_mod
! end module random_numerical_recipes_mod
