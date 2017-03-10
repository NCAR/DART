! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_random

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  nag_wrap_mod, only : g05ddf_wrap

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
