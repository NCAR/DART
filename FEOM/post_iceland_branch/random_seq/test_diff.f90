! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program test_diff

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

integer :: i, n
double precision :: r1, r2, dist, mean_dist, sd1, sd2, mean
type (random_seq_type) :: r

write(*, *) 'input mean'
read(*, *) mean
write(*, *) 'input sd 1'
read(*, *) sd1
write(*, *) 'input sd 2'
read(*, *) sd2



call init_random_seq(r)
n = 1000000
mean_dist = 0.0

do i = 1, n
   r1 = random_gaussian(r, dble(0.0), sd1)
   r2 = random_gaussian(r, mean, sd2)
   dist = dabs(r1 - r2)**2
   mean_dist = mean_dist + dist
end do

write(*, *)  'sample mean distance ', mean_dist / n
write(*, *) 'predicted distance ', sqrt(mean**2 + sd1**2 + sd2**2)


end program test_diff
