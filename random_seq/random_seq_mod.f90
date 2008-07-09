! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module random_seq_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use     types_mod, only : r8
use utilities_mod, only : register_module
use random_nr_mod, only : random_seq_type, init_ran1, ran1, gasdev

implicit none
private

public :: random_seq_type, init_random_seq, random_gaussian, &
   several_random_gaussians, random_uniform, twod_gaussians

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Gives ability to generate unique repeatable sequences of random numbers
! using random congruential package. Needed to allow different assim algorithms 
! that require random numbers to see identical observational sequences.

! Used to give different sequences a different but repeatable start
! There may be problems with incestuous series here; be cautious of this
! in the future.

integer :: seq_number = -1

logical, save :: module_initialized = .false.

contains


!========================================================================

subroutine init_random_seq(r, seed)
!----------------------------------------------------------------------
! You cannot generate any random numbers without calling this,
! so this is a sufficient entry point for initializing the module.
! An integer seed can be used to get a particular repeatable sequence.

!
implicit none

type(random_seq_type), intent(inout) :: r
integer, optional,     intent(in)    :: seed

if ( .not. module_initialized ) then
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

! Initialize the generator; use given seed if present, else sequence
if(present(seed)) then
   call init_ran1(r, seed)
else
   call init_ran1(r, seq_number)
   seq_number = seq_number - 1
endif

end subroutine init_random_seq

!========================================================================

function random_uniform(r)

implicit none

type(random_seq_type), intent(inout) :: r
real(r8) :: random_uniform

random_uniform = ran1(r)

end function random_uniform

!========================================================================

function random_gaussian(r, mean, standard_deviation) 

implicit none

type(random_seq_type), intent(inout) :: r
real(r8), intent(in) :: mean, standard_deviation
real(r8) :: random_gaussian

random_gaussian = gasdev(r) * standard_deviation + mean

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
   rnum(i) = gasdev(r) * standard_deviation + mean
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
x1 = gasdev(r)
x2 = gasdev(r)

! Use these to generate correlated
rnum(1) = mean(1) + a11 * x1
rnum(2) = mean(2) + a21 * x1 + a22 * x2

end subroutine twod_gaussians

!========================================================================

end module random_seq_mod
