! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program test_random

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  nag_wrap_mod, only : g05ddf_wrap

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: i, n1, n2
double precision :: rb(10), r2(10)
type (random_seq_type) :: r

write(*, *) 'how many in first slice'
read(*, *) n1
write(*, *) 'how many in second slice'
read(*, *) n2


call init_random_seq(r)

! Have a background sequence underway
do i = 1, n1
   rb(i) = g05ddf_wrap(dble(0.0), dble(1.0))
end do

! Also do a second repeatable sequence
do i = 1, n2
   r2(i) = random_gaussian(r, dble(0.0), dble(1.0))
end do

! Have a background sequence underway
do i = n1 + 1, 10
   rb(i) = g05ddf_wrap(dble(0.0), dble(1.0))
end do

! Also do a second repeatable sequence
do i = n2 + 1, 10
   r2(i) = random_gaussian(r, dble(0.0), dble(1.0))
end do

do i = 1, 10
   write(*, *) i, rb(i), r2(i)
end do

end program test_random
