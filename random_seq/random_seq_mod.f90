module random_seq_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

! use random_numerical_recipes_mod, only : random_seq_type, init_ran1, ran1, gasdev
use random_nr_mod, only : random_seq_type, init_ran1, ran1, gasdev

private
public random_seq_type, init_random_seq, random_gaussian, &
   several_random_gaussians, random_uniform, twod_gaussians


! Gives ability to generate unique repeatable sequences of random numbers
! using numerical recipes package. Needed to allow different assim algorithms 
! that ! require random numbers to see identical observational sequences.

! Used to give different sequences a different but repeatable start
! There may be problems with incestuous series here; be cautious of this
! in the future.
integer :: seq_number = -1

contains

!========================================================================

subroutine init_random_seq(r)

implicit none

type(random_seq_type), intent(inout) :: r

! Initialize the generator
call init_ran1(r, seq_number)
seq_number = seq_number - 1

end subroutine init_random_seq

!========================================================================

function random_uniform(r)

implicit none

type(random_seq_type), intent(inout) :: r
double precision :: random_uniform

random_uniform = ran1(r)

end function random_uniform

!========================================================================

function random_gaussian(r, mean, standard_deviation) 

implicit none

type(random_seq_type), intent(inout) :: r
double precision, intent(in) :: mean, standard_deviation
double precision :: random_gaussian

random_gaussian = gasdev(r) * standard_deviation + mean

end function random_gaussian

!========================================================================

subroutine several_random_gaussians(r, mean, standard_deviation, n, rnum)

implicit none

type(random_seq_type), intent(inout) :: r
double precision, intent(in) :: mean, standard_deviation
integer, intent(in) :: n
double precision, intent(out) :: rnum(n)

integer :: i

do i = 1, n
   rnum(i) = gasdev(r) * standard_deviation + mean
end do

end subroutine several_random_gaussians

!========================================================================

subroutine twod_gaussians(r, mean, cov, rnum)

implicit none

type(random_seq_type), intent(inout) :: r
double precision, intent(in) :: mean(2), cov(2, 2)
double precision, intent(out) :: rnum(2)

double precision :: a11, a21, a22, x1, x2

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
