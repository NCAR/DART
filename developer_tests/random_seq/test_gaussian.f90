! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_gaussian

! test the gaussian distribution random number generator routine
! with different means, standard deviations, and iteration counts.

use      types_mod, only : r8
use  utilities_mod, only : register_module, &
                           open_file, close_file, &
                           initialize_utilities, finalize_utilities, &
                           squeeze_out_blanks
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, j, n, f
real(r8) :: r1, sqdist, mean_sqdist
real(r8) :: mean, sd, var, compvar, compsd, compdiff
logical  :: write_this_one
character(len=256) :: fname, temp
character(len=50) :: formf = '(I12,2(F14.6),1(F24.16),2(F12.6))'
logical :: write_me = .true.       ! if true, write each distribution into a file for later diagnostics
integer :: write_limit = 1000000   ! but only if rep count is not greater than this limit

type t_inputs
 real(r8) :: t_mean
 real(r8) :: t_stddev
 integer  :: t_nreps
end type


! to add more tests or change the parameters, specify a test count
! and update the sets of (mean, stddev, sample size) here:

integer, parameter :: ntests = 21
type(t_inputs) :: t(ntests) = (/ &
    t_inputs(  0.0_r8,  1.0_r8,       100), &
    t_inputs(  0.0_r8,  1.0_r8,     10000), &
    t_inputs(  0.0_r8,  1.0_r8,   1000000), &
    t_inputs(  0.0_r8,  1.0_r8, 100000000), &
    t_inputs(  1.0_r8,  5.0_r8,       100), &
    t_inputs(  1.0_r8,  5.0_r8,     10000), &
    t_inputs(  1.0_r8,  5.0_r8,   1000000), &
    t_inputs(  1.0_r8,  5.0_r8, 100000000), &
    t_inputs(  6.0_r8,  8.0_r8,   1000000), &
    t_inputs(-86.0_r8,  8.0_r8,   1000000), &
    t_inputs(  0.6_r8,  8.0_r8,   1000000), &
    t_inputs( -6.0_r8,  8.0_r8,   1000000), &
    t_inputs(956.0_r8,  8.0_r8,   1000000), &
    t_inputs( 36.0_r8,  8.0_r8,   1000000), &
    t_inputs(  6.0_r8,  1.0_r8,   1000000), &
    t_inputs(  6.0_r8, 18.0_r8,   1000000), &
    t_inputs(  6.0_r8, 99.0_r8,   1000000), &
    t_inputs(  6.0_r8,  0.1_r8,   1000000), &
    t_inputs(  6.0_r8,  0.8_r8,   1000000), &
    t_inputs(  6.0_r8,  0.5_r8,   1000000), &
    t_inputs(  6.0_r8, 33.5_r8,   1000000)  /)


call initialize_utilities('test_gaussian')
call register_module(source,revision,revdate)

write(*, *) ''
write(*, *) 'sample size       input mean & std dev        computed std dev       diff       % diff '
write(*, *) ''

do j=1, ntests

   call init_random_seq(r, 5)

   mean = t(j)%t_mean
   sd = t(j)%t_stddev
   n = t(j)%t_nreps

   ! save all values in a file for post-plotting?
   write_this_one = (write_me .and. n <= write_limit)

   if (write_this_one) then
      write(temp, "(A,F8.3,A,F8.3,A,I10)") "gauss_", mean, "_", sd, "_", n
      call squeeze_out_blanks(temp, fname)
      f = open_file(fname)
   endif

   ! analytical values:
   mean_sqdist = 0.0_r8

   do i = 1, n

      r1 = random_gaussian(r, mean, sd)

      if (write_this_one) write(f,*) r1

      sqdist = (mean-r1)**2
      mean_sqdist = mean_sqdist + sqdist
   end do

   if (write_this_one) call close_file(f)

   ! computed values:
   compvar = mean_sqdist / n
   compsd = sqrt(compvar)
   compdiff = compsd - sd

   write(*, formf) n, mean, sd, compsd, compdiff, abs(compdiff/sd) * 100._r8
   
enddo

call finalize_utilities()

end program test_gaussian

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
