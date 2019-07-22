! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_gamma

! test the gamma distribution random number generator routine
! with different means, standard deviations, and iteration counts.

use      types_mod, only : r8
use  utilities_mod, only : register_module, &
                           open_file, close_file, &
                           initialize_utilities, finalize_utilities, &
                           squeeze_out_blanks
use random_seq_mod, only : random_seq_type, init_random_seq, random_gamma

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, j, n, f
real(r8) :: r1, sum, sumsq, k, h
real(r8) :: mean, var, sd, compmean, compvar, compsd, compdiffm, compdiffsd
logical  :: write_this_one
character(len=50) :: formf = '(I12,2(F8.3),4(F12.5),2(F16.10),2(F12.5))'
character(len=256) :: fname, temp

logical :: write_me = .true.       ! if true, write each distribution into a file for later diagnostics
integer :: write_limit = 1000000   ! but only if rep count is not greater than this limit

type t_inputs
 real(r8) :: t_kappa
 real(r8) :: t_theta
 integer  :: t_nreps
end type


! to add more tests or change the parameters, specify a test count
! and update the sets of inputs here:

integer, parameter :: ntests = 27
type(t_inputs) :: t(ntests) = (/ &
    t_inputs(  1.0_r8,  1.0_r8,       100), &
    t_inputs(  1.0_r8,  1.0_r8,     10000), &
    t_inputs(  1.0_r8,  1.0_r8,   1000000), &
    t_inputs(  1.0_r8,  1.0_r8, 100000000), &
    t_inputs(  1.0_r8,  1.5_r8,   5000000), &
    t_inputs(  1.0_r8,  2.5_r8,   5000000), &
    t_inputs(  1.0_r8,  5.0_r8,       100), &
    t_inputs(  1.0_r8,  5.0_r8,     10000), &
    t_inputs(  1.0_r8,  5.0_r8,   1000000), &
    t_inputs(  1.0_r8,  5.0_r8, 100000000), &
    t_inputs(  1.0_r8,  5.0_r8,   1000000), &
    t_inputs(  2.0_r8,  5.0_r8,   1000000), &
    t_inputs(  4.0_r8,  5.0_r8,   1000000), &
    t_inputs(  8.0_r8,  5.0_r8,   1000000), &
    t_inputs(  6.0_r8,  8.0_r8,   1000000), &
    t_inputs( 86.0_r8,  8.0_r8,   1000000), &
    t_inputs(  0.6_r8,  8.0_r8,   1000000), &
    t_inputs(  6.0_r8,  8.0_r8,   1000000), &
    t_inputs(956.0_r8,  8.0_r8,   1000000), &
    t_inputs( 36.0_r8,  8.0_r8,   1000000), &
    t_inputs(  6.0_r8,  1.0_r8,   1000000), &
    t_inputs(  6.0_r8, 18.0_r8,   1000000), &
    t_inputs(  6.0_r8, 99.0_r8,   1000000), &
    t_inputs(  6.0_r8,  0.1_r8,   1000000), &
    t_inputs(  6.0_r8,  0.8_r8,   1000000), &
    t_inputs(  6.0_r8,  0.5_r8,   1000000), &
    t_inputs(  6.0_r8, 33.5_r8,   1000000)  /)


call initialize_utilities('test_gamma')
call register_module(source,revision,revdate)

write(*, *) ''
write(*, *) 'sample size     input k & t       computed mean, sd         actual mean, sd       diff mean         diff sd        % diff mean,  sd'
write(*, *) ''

do j=1, ntests

   call init_random_seq(r, 5)

   k = t(j)%t_kappa
   h = t(j)%t_theta
   n = t(j)%t_nreps

   ! save all values in a file for post-plotting?
   write_this_one = (write_me .and. n <= write_limit)

   if (write_this_one) then
      write(temp, "(A,F8.3,A,F8.3,A,I10)") "gamma_", k, "_", h, "_", n
      call squeeze_out_blanks(temp, fname)
      f = open_file(fname)
   endif

   ! analytical values:
   mean = k * h  
   sd = sqrt(k * (h*h))

   sum = 0.0_r8
   sumsq = 0.0_r8

   do i = 1, n
      r1 = random_gamma(r, k, h)

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

   write(*, formf) n, k, h, mean, sd, compmean, compsd, compdiffm, compdiffsd, &
                      abs(compdiffm/mean) * 100._r8, abs(compdiffsd/sd) * 100._r8
   

enddo

call finalize_utilities()

end program test_gamma

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
