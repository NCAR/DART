! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_diff

use     types_mod,  only : r4, r8, digits12
use utilities_mod,  only : register_module, error_handler, E_ERR, E_MSG, &
                           initialize_utilities, finalize_utilities, &
                           find_namelist_in_file, check_namelist_read
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter :: n = 1000000
integer :: i, iunit, io
type (random_seq_type) :: r

double precision ::  dpr1,  dpr2,  dpdist,  dpmean_dist
real(digits12)   :: d12r1, d12r2, d12dist, d12mean_dist
real(r8)         ::  r8r1,  r8r2,  r8dist,  r8mean_dist
real(r4)         ::  r4r1,  r4r2,  r4dist,  r4mean_dist

!-----------------------------------------------------------------------------
! Namelist with default values

double precision :: mean = 0.0d0, sd1 = 1.0d0, sd2 = 2.0d0
integer :: iseed = 123

namelist /test_diff_nml/ mean, sd1, sd2, iseed

write(*, *) 'double precision       is ', kind(dpr1)
write(*, *) 'digits12 is defined to be ', digits12
write(*, *) 'r8       is defined to be ', r8
write(*, *) 'r4       is defined to be ', r4

call initialize_utilities('test_diff')
call register_module(source,revision,revdate)

! Read the namelist entry  
call find_namelist_in_file("input.nml", "test_diff_nml", iunit)
read(iunit, nml = test_diff_nml, iostat = io)
call check_namelist_read(iunit, io, "test_diff_nml")


call init_random_seq(r,iseed)
!-----------------------------------------------------------------------------
write(*, *) 'double precision test:'
!-----------------------------------------------------------------------------

dpmean_dist = 0.0d0
do i = 1, n
   dpr1        = random_gaussian(r, 0.0d0, sd1)
   dpr2        = random_gaussian(r,  mean, sd2)
   dpdist      = dabs(dpr1 - dpr2)**2
   dpmean_dist = dpmean_dist + dpdist
end do

write(*, *) 'sample mean distance  ', dpmean_dist / n
write(*, *) 'predicted distance    ', sqrt(mean**2 + sd1**2 + sd2**2)
write(*, *) ''


!-----------------------------------------------------------------------------
write(*, *) 'digits12 test:'
!-----------------------------------------------------------------------------

call init_random_seq(r,iseed)
d12mean_dist = 0.0_digits12
do i = 1, n
   d12r1        = random_gaussian(r, 0.0d0, sd1)
   d12r2        = random_gaussian(r,  mean, sd2)
   d12dist      = abs(d12r1 - d12r2)**2
   d12mean_dist = d12mean_dist + d12dist
end do

write(*, *) 'sample mean distance  ', d12mean_dist / n
write(*, *) 'predicted distance    ', sqrt(mean**2 + sd1**2 + sd2**2)
write(*, *) ''

!-----------------------------------------------------------------------------
write(*, *) 'r8 test:'
!-----------------------------------------------------------------------------

call init_random_seq(r,iseed)
r8mean_dist = 0.0_r8
do i = 1, n
   r8r1        = random_gaussian(r, 0.0d0, sd1)
   r8r2        = random_gaussian(r,  mean, sd2)
   r8dist      = abs(r8r1 - r8r2)**2
   r8mean_dist = r8mean_dist + r8dist
end do

write(*, *) 'sample mean distance  ', r8mean_dist / n
write(*, *) 'predicted distance    ', sqrt(mean**2 + sd1**2 + sd2**2)
write(*, *) ''

!-----------------------------------------------------------------------------
write(*, *) 'r4 test:'
!-----------------------------------------------------------------------------

call init_random_seq(r,iseed)
r4mean_dist = 0.0_r4
do i = 1, n
   r4r1        = random_gaussian(r, 0.0d0, sd1)
   r4r2        = random_gaussian(r,  mean, sd2)
   r4dist      = abs(r4r1 - r4r2)**2
   r4mean_dist = r4mean_dist + r4dist
end do

write(*, *) 'sample mean distance  ', r4mean_dist / n
write(*, *) 'predicted distance    ', sqrt(mean**2 + sd1**2 + sd2**2)
write(*, *) ''

call error_handler(E_MSG, 'test_diff', 'Finished successfully.',&
                   source,revision,revdate)
call finalize_utilities()


end program test_diff

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
