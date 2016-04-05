! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program test_random_nr

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

use random_nr_mod, only : random_seq_type, init_ran1, ran1, gasdev

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, n
double precision :: r1, dist, mean_dist

mean_dist = 0.0
call init_ran1(r, -5)

n = 10000000

do i = 1, n
   r1 = gasdev(r)
   dist = dabs(r1)
   mean_dist = mean_dist + dist
end do

write(*, *) 'sd is ', mean_dist / n

end program test_random_nr
