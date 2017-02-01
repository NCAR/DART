! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
!
! do a simple test that the mean distance between sample
! draws from two gaussian distributions are close to the
! expected value.  indirect test of the random number generator
! which is used in the gaussian random numbers.

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
character(len=16) :: wformatI = '(A,I8)'
character(len=16) :: wformat1 = '(A,1(F25.16))'
character(len=16) :: wformat2 = '(A,2(F25.16))'

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

write(*, *) ''
write(*, wformatI) 'sample size = ', n
write(*, *) ''
write(*, wformat2) 'mean, stddev 1 = ', 0.0_r8, sd1
write(*, wformat2) 'mean, stddev 2 = ', mean, sd2
write(*, *) ''


call init_random_seq(r,iseed)
mean_dist2 = 0.0_r8
do i = 1, n
   r1 = random_gaussian(r, 0.0_r8, sd1)
   r2 = random_gaussian(r,  mean,  sd2)
   dist2 = (r1 - r2)**2
   mean_dist2 = mean_dist2 + dist2
end do

write(*, wformat1) 'predicted distance   = ', sqrt(mean**2 + sd1**2 + sd2**2)
write(*, wformat1) 'sample mean distance = ', sqrt(mean_dist2 / n)
write(*, *) ''

call finalize_utilities()


end program test_diff

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
