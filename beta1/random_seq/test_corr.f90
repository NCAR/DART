! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_random

use random_seq_mod, only : random_seq_type, init_random_seq, twod_gaussians

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
double precision, allocatable :: rnum(:, :)
double precision :: c(2, 2), correl, sample_correl, mean(2), mean_2
double precision :: mean_correl, var
integer :: n, i, j, n_cycles

write(*, *) 'how many pairs'
read(*, *) n
write(*, *) 'input a correlation'
read(*, *) correl 

write(*, *) 'how many cycles'
read(*, *) n_cycles

mean_correl = 0.0
var = 0.0

! Allocate space for the pairs of data
allocate(rnum(2, n))

call init_random_seq(r)

! Mean for generated
mean = 0.0; mean_2 = 0.0

! Generate the covariance matrix
c(1, 1) = 1.0
c(2, 2) = 1.0
c(1, 2) = correl
c(2, 1) = correl

do j = 1, n_cycles
   
   do i = 1, n
      call twod_gaussians(r, mean, c, rnum(:, i))
   end do

! Compute the correlation of the random numbers
   call comp_correl(rnum, n, sample_correl)

   mean_correl = mean_correl + sample_correl
   mean_2 = mean_2 + sample_correl**2
   var = var + (sample_correl - correl) ** 2
end do

mean_correl = mean_correl / n_cycles
write(*, *) 'mean correl is  ', real(mean_correl)
write(*, *) 'variance is     ', real((mean_2 - n_cycles * mean_correl**2) / (n_cycles - 1))
write(*, *) 'predicted is    ', real((1. - correl **2)**2 / (n - 1.))
! Also did this by mistake earlier doesn't use sample mean
write(*, *) 'incorrect var is (uses correl as mean, not sample)', real(var / (n_cycles - 1))

end program test_random

!-----------------------------------------------------
 
subroutine comp_correl(ens, n, correl)
 
implicit none
 
integer, intent(in) :: n
double precision, intent(in) :: ens(2, n)
double precision, intent(out) :: correl
double precision :: sum_x, sum_y, sum_xy, sum_x2, sum_y2


sum_x = sum(ens(2, :))
sum_y = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))
 
! Computation of correlation 
sum_y2 = sum(ens(1, :) * ens(1, :))

correl = (n * sum_xy - sum_x * sum_y) / &
   sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))
 
end subroutine comp_correl

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
