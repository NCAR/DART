! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program qc_test

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

! Think of dividing area around distribution into 9 boxes defined by
! the two observations. When the observations have the same sign, boxes
! in lower left and upper right summed (2* either one since they are
! symmetric) is two sided probability that something further out would
! have been picked. When x and y have different signs, the right answer
! is the sum of probabilities in upper left and lower right boxes (again
! 2 * either one by symmetry).

use nag_wrap_mod, only : g01haf_wrap, g05ddf_wrap
use qc_mod

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: i, j, k, ifail = 0
double precision :: x, y, rho, ens(2, 20), obs(2), obs_variance(2), mean(2)

rho = 0.0

do i = 1, 100
   do j = 1, 2
      obs(j) = g05ddf_wrap(dble(0.0), dble(1.0))
   end do
   do j = 1, 20
      do k = 1, 2
         ens(k, j) = g05ddf_wrap(dble(0.0), dble(1.0))
      end do
   end do
   obs_variance = 1.0

   write(*, *) i, pair_ob_qc(ens, 20, obs, obs_variance)
end do
end program qc_test
