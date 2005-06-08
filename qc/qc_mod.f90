! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module qc_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

use nag_wrap_mod, only : g01haf_wrap, g01eaf_wrap, g01ecf_wrap
use assim_tools_mod, only : sample_cov

private

public :: single_ob_qc, pair_ob_qc

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

contains

!==========================================================================

function single_ob_qc(ens, ens_size, ob, obs_var)

implicit none

double precision :: single_ob_qc
integer, intent(in) :: ens_size
double precision, intent(in) :: ens(ens_size), ob, obs_var

double precision :: cov(1, 1), cor(1, 1), mean(1), scaled_ob
integer :: ifail = 0

! First compute the covariance and mean of the ensemble
call sample_cov(ens, 1, ens_size, cov, cor, mean)
!write(*, *) '----------------------------'
!write(*, *) 'ensemble ', ens
!write(*, *) 'single ob qc cov, mean, ob ', cov(1, 1), mean(1), ob

! Shift and scale the ob so that it is from a unit normal distribution
! IS TIS THE RIGHT SCALING??? VERIFY
! Need to add in observational variance (sum of independents variance sums)
scaled_ob = (ob - mean(1)) / sqrt(cov(1, 1) + obs_var)

! Compute probability that ob should be more of an outlier than this
single_ob_qc = g01eaf_wrap('S', scaled_ob, ifail)
!write(*, *) 'single ob qc scaled ', scaled_ob, single_ob_qc

end function single_ob_qc

!==========================================================================

function pair_ob_qc(ens, ens_size, obs, obs_variance)

! Given a prior ensemble for two obs distributions, the real obs, and the
! instrumental obs_variance (assumed diagonal here although this isn't
! necessary, determines the probability that the actual obs would be further
! from the mean of the distribution given the covariance structure.
! This is done by first rotating to a space where the expected observation
! covariance structure is diagonal, normalizing to get a bivariate unit
! normal, and then using a chi-square test to determine the chance that
! the obs would have been further away from the mean.  This could easily
! be extended to even higher dimensions if one could think of a reason, in
! order to check trios of obs, etc. For really hardcore obs error correction
! one can think of cases where this might be useful, but in general it's
! going to be just too expensive.

implicit none

double precision :: pair_ob_qc
integer, intent(in) :: ens_size
double precision, intent(in) :: ens(2, ens_size), obs(2), obs_variance(2)

double precision :: cov(2, 2), cor(2, 2), mean(2), lam(2), sv(2, 2), svt(2, 2)
double precision :: rens(2, ens_size), robs(2), rcov(2, 2), dist2
integer :: i, j, ifail = 0

! First compute the covariance and mean of the ensemble
call sample_cov(ens, 2, ens_size, cov, cor, mean)

! Remove the mean from the ensemble
do i = 1, ens_size
   rens(:, i) = ens(:, i) - mean
end do
! Remove the ensemble mean from the observations
robs = obs - mean

! Covariance of expected distribution is ensemble PLUS obs error
do i = 1, 2
   cov(i, i) = cov(i, i) + obs_variance(i)
end do

! Do SVD to get to a frame where this thing has diag cov
call twod_sym_svd(cov, lam, sv)

! Multiply the ensemble members by sv, svt and check out the result
svt = transpose(sv)
rens = matmul(sv, rens)

! Move the obs to this frame and also bring the covariance
robs = matmul(sv, robs)

rcov = matmul(sv, matmul(cov, svt))

! Now scale the obs by this to get standard normal
do i = 1, 2
   robs(i) = robs(i) / sqrt(rcov(i, i))
end do

! Compute distance from origin for observation
dist2 = robs(1)**2 + robs(2)**2

! Try the chi-square distribution
pair_ob_qc = g01ecf_wrap('U', dist2, dble(2.0), ifail)

end function pair_ob_qc

!==========================================================================

subroutine twod_sym_svd(a, lam, x)

! Relatively efficient computation of singular vectors and singular values
! of a 2d symmetric matrix a.

implicit none

double precision, intent(in) :: a(2, 2)
double precision, intent(out) :: lam(2), x(2, 2)

double precision :: b, c, quad_tail

b = -1.0 * (a(1, 1) + a(2, 2))
c = (a(1, 1) * a(2, 2) - a(2, 1)**2)
quad_tail = sqrt(b**2 - 4. * c) / 2.0
lam(1) = -b/2.0 + quad_tail
lam(2) = -b/2.0 - quad_tail

! Compute singular vectors
! Compute normalizing factor to get unit vectors
x(1, 2) = sqrt(1. / (1. + ((lam(1) - a(2, 2)) / a(1, 2))**2))
x(1, 1) = x(1, 2) * (lam(1) - a(2, 2)) / a(1, 2)
x(2, 2) = x(1, 1)
x(2, 1) = -x(1, 2)

end subroutine twod_sym_svd

!========================================================================

end module qc_mod
