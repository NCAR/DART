! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module random_nr_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use     types_mod, only : r8, digits12
use utilities_mod, only : register_module, error_handler, E_ERR

implicit none
private

public :: random_seq_type, init_ran1, ran1, gasdev

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer, parameter :: m1 = 259200, ia1 = 7141, ic1 = 54773
integer, parameter :: m2 = 134456, ia2 = 8121, ic2 = 28411
integer, parameter :: m3 = 243000, ia3 = 4561, ic3 = 51349
real(digits12), parameter :: rm1 = 1.0_digits12/m1, rm2 = 1.0_digits12/m2

type random_seq_type
   private
   integer :: ix1, ix2, ix3, iset
   real(digits12) :: r(97), gset
end type random_seq_type

logical, save :: module_initialized = .false.


contains



subroutine initialize_module

   call register_module(source,revision,revdate)
   module_initialized = .true.

end subroutine initialize_module




!-------------------------------------------------------------------

!  A random congruential random number generator (see Knuth)
subroutine init_ran1(s, temp)

implicit none

integer, intent(in) :: temp
type(random_seq_type), intent(out) :: s
integer j

if ( .not. module_initialized ) call initialize_module

! Initialize the generator for use with
! repeatable sequences

s%ix1 = mod(ic1 - temp, m1)
s%ix1 = mod(ia1*s%ix1 + ic1, m1)
s%ix2 = mod(s%ix1, m2)
s%ix1 = mod(ia1*s%ix1 + ic1, m1)
s%ix3 = mod(s%ix1, m3)
do j = 1, 97
   s%ix1 = mod(ia1*s%ix1 + ic1, m1)
   s%ix2 = mod(ia2*s%ix2 + ic2, m2)
   s%r(j) = (s%ix1 + s%ix2*rm2)*rm1
end do

! Initialize the value needed for Gaussian efficiency
s%iset = 0

end subroutine init_ran1

!-----------------------------------------------------------------

!  A random congruential random number generator (see Knuth)
function ran1(s)

implicit none

type(random_seq_type), intent(inout) :: s
real(r8) :: ran1

integer :: j

if ( .not. module_initialized ) call initialize_module

!  Gives a U(0,1) random number

s%ix1 = mod(ia1*s%ix1 + ic1, m1)
s%ix2 = mod(ia2*s%ix2 + ic2, m2)
s%ix3 = mod(ia3*s%ix3 + ic3, m3)
j = 1 + (97*s%ix3) / m3
if(j > 97 .or. j < 1) then
      call error_handler(E_ERR,' ran1', 'Fatal error in random_nr_mod', source, revision, revdate)
endif
ran1 = s%r(j)
s%r(j) = (s%ix1 + s%ix2*rm2)*rm1
return
end function ran1

!---------------------------------------------------------------------

function gasdev(s)

! Returns a N(0, 1) random number

implicit none

type(random_seq_type), intent(inout) :: s
real(r8) :: gasdev

real(digits12) :: v1, v2, r, fac

if ( .not. module_initialized ) call initialize_module

if(s%iset == 0) then
10 v1 = 2.0_digits12 * ran1(s) - 1.0_digits12
   v2 = 2.0_digits12 * ran1(s) - 1.0_digits12
   r = v1**2 + v2**2
   if(r >= 1.0_digits12) goto 10
   fac = sqrt(-2.0_digits12 * log(r) / r)
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
