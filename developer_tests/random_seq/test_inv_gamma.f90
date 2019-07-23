! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_inv_gamma

! test the gamma distribution random number generator routine
! with different means, standard deviations, and iteration counts.

use      types_mod, only : r8
use  utilities_mod, only : register_module, &
                           open_file, close_file, &
                           initialize_utilities, finalize_utilities, &
                           squeeze_out_blanks
use random_seq_mod, only : random_seq_type, init_random_seq, random_inverse_gamma

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, j, n, f
real(r8) :: r1, sum, sumsq, a, b, l
real(r8) :: mean, var, sd, compmean, compvar, compsd, compdiffm, compdiffsd
logical  :: write_this_one
character(len=256) :: fname, temp

character(len=50) :: formf = '(I12,2(F8.3),4(F12.5),2(F16.10),2(F12.5))'
logical :: write_me = .true.       ! if true, write each distribution into a file for later diagnostics
integer :: write_limit = 1000000   ! but only if rep count is not greater than this limit

type t_inputs
 real(r8) :: t_alpha
 real(r8) :: t_beta
 integer  :: t_nreps
end type


! to add more tests or change the parameters, specify a test count
! and update the sets of inputs here:

integer, parameter :: ntests = 27
type(t_inputs) :: t(ntests) = (/ &
     t_inputs(  3.0_r8,  0.20_r8,       100), &
     t_inputs(  3.0_r8,  0.20_r8,     10000), &
     t_inputs(  3.0_r8,  0.20_r8,   1000000), &
     t_inputs(  3.0_r8,  0.20_r8, 100000000), &
     t_inputs(  5.0_r8,  0.40_r8,   4000000), &
     t_inputs(  5.0_r8,  2.00_r8,   5000000), &
     t_inputs(  5.0_r8,  0.60_r8,       100), &
     t_inputs(  5.0_r8,  0.60_r8,     10000), &
     t_inputs(  5.0_r8,  0.60_r8,   1000000), &
     t_inputs(  5.0_r8,  0.60_r8, 100000000), &
     t_inputs(  3.0_r8,  0.50_r8,   1000000), &
     t_inputs(  4.0_r8,  0.50_r8,   1000000), &
     t_inputs(  5.0_r8,  0.50_r8,   1000000), &
     t_inputs(  6.0_r8,  0.50_r8,   1000000), &
     t_inputs(  6.0_r8,  0.01_r8,   1000000), &
     t_inputs( 86.0_r8,  0.01_r8,   1000000), &
     t_inputs(  6.6_r8,  0.01_r8,   1000000), &
     t_inputs(  6.0_r8,  0.01_r8,   1000000), &
     t_inputs(956.0_r8,  0.01_r8,   1000000), &
     t_inputs( 36.0_r8,  0.01_r8,   1000000), &
     t_inputs(  6.0_r8,  0.03_r8,   1000000), &
     t_inputs(  6.0_r8,  0.09_r8,   1000000), &
     t_inputs(  6.0_r8,  9.00_r8,   1000000), &
     t_inputs(  6.0_r8,  0.10_r8,   1000000), &
     t_inputs(  6.0_r8,  0.80_r8,   1000000), &
     t_inputs(  6.0_r8,  0.50_r8,   1000000), &
     t_inputs(  6.0_r8,  3.50_r8,   1000000)  /)


call initialize_utilities('test_inv_gamma')
call register_module(source,revision,revdate)

write(*, *) ''
write(*, *) 'sample size     input a & b       computed mean, sd         actual mean, sd       diff mean         diff sd        % diff mean,  sd'
write(*, *) ''

do j=1, ntests

   call init_random_seq(r, 5)

   a = t(j)%t_alpha
   b = t(j)%t_beta
   n = t(j)%t_nreps

   ! save all values in a file for post-plotting?
   write_this_one = (write_me .and. n <= write_limit)

   if (write_this_one) then
      write(temp, "(A,F8.3,A,F8.3,A,I10)") "invgamma_", a, "_", b, "_", n
      call squeeze_out_blanks(temp, fname)
      f = open_file(fname)
   endif

   ! analytical values:
   
   ! lambda = scale, beta = rate
   l = 1.0_r8 / b

   if (a > 1.0_r8) then
      mean  = b / (a - 1.0_r8)
   else
      mean  = b  ! ??
   endif
   if (a > 2.0_r8) then
      sd  = sqrt((b*b) / ((a - 1.0)**2 * (a - 2.0_r8)))
   else
      sd  = b    ! ??
   endif

   sum = 0.0_r8
   sumsq = 0.0_r8

   ! inverse gamma wants shape and scale. 
   do i = 1, n
      r1 = random_inverse_gamma(r, a, l)

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

   write(*, formf) n, a, b, mean, sd, compmean, compsd, compdiffm, compdiffsd, &
                      abs(compdiffm/mean) * 100._r8, abs(compdiffsd/sd) * 100._r8
   

enddo

call finalize_utilities()

end program test_inv_gamma

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
