! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_random_gsl

use      types_mod, only : r4, r8, digits12
use  utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                           initialize_utilities, finalize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, n

double precision :: dpr1, dpdist, dpmean_dist
real(r8)         :: r8r1, r8dist, r8mean_dist
real(r4)         :: r4r1, r4dist, r4mean_dist
real(digits12)   :: d12r1, d12dist, d12mean_dist

call initialize_utilities('test_random')
call register_module(source,revision,revdate)

write(*, *) 'double precision       is ', kind(dpr1)
write(*, *) 'digits12 is defined to be ', digits12
write(*, *) 'r8       is defined to be ', r8
write(*, *) 'r4       is defined to be ', r4

n = 10000000


call init_random_seq(r, -5)
d12mean_dist = 0.0
do i = 1, n
   d12r1        = random_gaussian(r, 0.0_r8, 1.0_r8)
   d12dist      = abs(d12r1)
   d12mean_dist = d12mean_dist + d12dist
end do
write(*, *) 'digits12         sd is ', d12mean_dist / n


call init_random_seq(r, -5)
dpmean_dist = 0.0
do i = 1, n
   dpr1        = random_gaussian(r, 0.0_r8, 1.0_r8)
   dpdist      = dabs(dpr1)
   dpmean_dist = dpmean_dist + dpdist
end do
write(*, *) 'double precision sd is ', dpmean_dist / n


call init_random_seq(r, -5)
r8mean_dist = 0.0_r8
do i = 1, n
   r8r1        = random_gaussian(r, 0.0_r8, 1.0_r8)
   r8dist      = abs(r8r1)
   r8mean_dist = r8mean_dist + r8dist
end do
write(*, *) 'r8               sd is ', r8mean_dist / n


call init_random_seq(r, -5)
r4mean_dist = 0.0_r4
do i = 1, n
   r4r1        = random_gaussian(r, 0.0_r8, 1.0_r8)
   r4dist      = abs(r4r1)
   r4mean_dist = r4mean_dist + r4dist
end do
write(*, *) 'r4               sd is ', r4mean_dist / n

call error_handler(E_MSG, 'test_random_gsl', 'Finished successfully.',&
                   source,revision,revdate)
call finalize_utilities()

end program test_random_gsl

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
