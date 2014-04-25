! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
!
! FIXME:  the origins of this test program are lost to time.
! i can see that it is trying to verify that the gaussian distribution
! generator is working as expected, and it's also testing out the
! computation at various (well, 2) precisions.  but i do not
! believe that the computed and expected values are being computed
! right.  the original code was computing the sqrt(distance**2) but
! these are 1d values, not 2d, so the distance is simply abs(r1-r2)
! without squares.  

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

real(r8) :: r1, r2, dist2, mean_dist2

!-----------------------------------------------------------------------------
! Namelist with default values

real(r8) :: mean = 0.0, sd1 = 1.0, sd2 = 2.0, predicted
integer :: iseed = 123

namelist /test_diff_nml/ mean, sd1, sd2, iseed

call initialize_utilities('test_diff')
call register_module(source,revision,revdate)

! Read the namelist entry  
call find_namelist_in_file("input.nml", "test_diff_nml", iunit)
read(iunit, nml = test_diff_nml, iostat = io)
call check_namelist_read(iunit, io, "test_diff_nml")

write(*, *) 'sample size = ', n
write(*, *) 'mean 1, 2 = ', 0.0, mean
write(*, *) 'stddev 1, 2  = ', sd1, sd2
write(*, *) 'predicted distance = ', sqrt(mean**2 + sd1**2 + sd2**2)
write(*, *) ''


call init_random_seq(r,iseed)
mean_dist2 = 0.0d0
do i = 1, n
   r1 = random_gaussian(r, 0.0d0, sd1)
   r2 = random_gaussian(r,  mean, sd2)
   dist2 = (r1 - r2)**2
   mean_dist2 = mean_dist2 + dist2
end do

write(*, *) 'sample mean distance = ', sqrt(mean_dist2 / n)
write(*, *) ''

call finalize_utilities()


end program test_diff

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
