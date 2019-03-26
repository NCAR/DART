! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_exp

! test the exponential distribution random number generator routine
! with different means, standard deviations, and iteration counts.

use      types_mod, only : r8
use  utilities_mod, only : register_module, &
                           open_file, close_file, &
                           initialize_utilities, finalize_utilities, &
                           squeeze_out_blanks
use random_seq_mod, only : random_seq_type, init_random_seq, random_exponential

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, j, n, f
real(r8) :: r1, sum, sumsq, rate
real(r8) :: mean, var, sd, compmean, compvar, compsd, compdiffm, compdiffsd
logical  :: write_this_one
character(len=256) :: fname, temp

character(len=50) :: formf = '(I12,1(F8.3),4(F12.5),2(F16.10),2(F12.5))'
logical :: write_me = .true.       ! if true, write each distribution into a file for later diagnostics
integer :: write_limit = 1000000   ! but only if rep count is not greater than this limit

type t_inputs
 real(r8) :: t_rate
 integer  :: t_nreps
end type


! to add more tests or change the parameters, specify a test count
! and update the sets of inputs here:

integer, parameter :: ntests = 25
type(t_inputs) :: t(ntests) = (/ &
    t_inputs( 1.0_r8,       100), &
    t_inputs( 1.0_r8,     10000), &
    t_inputs( 1.0_r8,   1000000), &
    t_inputs( 1.0_r8, 100000000), &
    t_inputs( 5.0_r8,       100), &
    t_inputs( 5.0_r8,     10000), &
    t_inputs( 5.0_r8,   1000000), &
    t_inputs( 5.0_r8, 100000000), &
    t_inputs( 5.0_r8,   1000000), &
    t_inputs( 0.5_r8,   1000000), &
    t_inputs( 0.8_r8,   1000000), &
    t_inputs( 0.9_r8,   1000000), &
    t_inputs( 8.0_r8,   1000000), &
    t_inputs( 8.2_r8,   1000000), &
    t_inputs( 8.4_r8,   1000000), &
    t_inputs( 8.6_r8,   1000000), &
    t_inputs( 8.8_r8,   1000000), &
    t_inputs( 9.0_r8,   1000000), &
    t_inputs( 1.0_r8,   1000000), &
    t_inputs(18.0_r8,   1000000), &
    t_inputs(99.0_r8,   1000000), &
    t_inputs( 0.1_r8,   1000000), &
    t_inputs( 0.8_r8,   1000000), &
    t_inputs( 0.5_r8,   1000000), &
    t_inputs(33.5_r8,   1000000)  /)


call initialize_utilities('test_exp')
call register_module(source,revision,revdate)

write(*, *) ''
write(*, *) 'sample size   input rate computed mean, sd         actual mean, sd       diff mean         diff sd        % diff mean,  sd'
write(*, *) ''

do j=1, ntests

   call init_random_seq(r, 5)

   rate = t(j)%t_rate
   n = t(j)%t_nreps

   ! save all values in a file for post-plotting?
   write_this_one = (write_me .and. n <= write_limit)

   if (write_this_one) then
      write(temp, "(A,F8.3,A,I10)") "exp_", rate, "_", n
      call squeeze_out_blanks(temp, fname)
      f = open_file(fname)
   endif

   ! analytical values:
   mean = 1.0_r8/rate
   sd = mean

   sum = 0.0_r8
   sumsq = 0.0_r8

   do i = 1, n
      r1 = random_exponential(r, rate)

      if (write_this_one) write(f,*) r1

      sum = sum + r1
      sumsq = sumsq + (r1 - mean)**2
   end do

   if (write_this_one) call close_file(f)

   ! computed values:
   compmean = sum / n
   compsd = sqrt(sumsq/n)

   ! differences
   compdiffm  = compmean - mean
   compdiffsd = compsd - sd

   write(*, formf) n, rate, mean, sd, compmean, compsd, compdiffm, compdiffsd, &
                      abs(compdiffm/mean) * 100._r8, abs(compdiffsd/sd) * 100._r8
   

enddo

call finalize_utilities()

end program test_exp

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
