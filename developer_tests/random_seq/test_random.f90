! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_random

! test the uniform random number generator routine

use      types_mod, only : r8
use  utilities_mod, only : register_module, &
                           open_file, close_file, &
                           initialize_utilities, finalize_utilities, &
                           squeeze_out_blanks
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, j, n, f
real(r8) :: r1, sum
real(r8) :: mean, compmean, compdiffm
logical  :: write_this_one
character(len=50) :: formf = '(I12,4(F16.5))'
character(len=256) :: fname, temp

logical :: write_me = .true.       ! if true, write each distribution into a file for later diagnostics
integer :: write_limit = 1000000   ! but only if rep count is not greater than this limit

type t_inputs
 integer  :: t_nreps
end type


! to add more tests or change the parameters, specify a test count
! and update the sets of inputs here:

integer, parameter :: ntests = 11
type(t_inputs) :: t(ntests) = (/ &
    t_inputs(       100), &
    t_inputs(       500), &
    t_inputs(      1000), &
    t_inputs(      5000), &
    t_inputs(     10000), &
    t_inputs(     50000), &
    t_inputs(    100000), &
    t_inputs(    500000), &
    t_inputs(   1000000), &
    t_inputs(   5000000), &
    t_inputs(  10000000) /)


call initialize_utilities('test_random')
call register_module(source,revision,revdate)

write(*, *) ''
write(*, *) 'sample size   computed mean     actual mean       diff mean     % diff mean'
write(*, *) ''

do j=1, ntests

   call init_random_seq(r, 5)

   n = t(j)%t_nreps

   ! save all values in a file for post-plotting?
   write_this_one = (write_me .and. n <= write_limit)

   if (write_this_one) then
      write(temp, "(A,I10)") "random_", n
      call squeeze_out_blanks(temp, fname)
      f = open_file(fname)
   endif

   ! analytical values:
   mean = 0.5

   sum = 0.0_r8

   do i = 1, n
      r1 = random_uniform(r)

      if (write_this_one) write(f,*) r1

      sum = sum + r1
   end do

   if (write_this_one) call close_file(f)

   ! computed values:
   compmean = sum / n

   ! differences
   compdiffm  = compmean - mean

   write(*, formf) n, mean, compmean, compdiffm, &
                   abs(compdiffm/mean) * 100._r8
   

enddo

call finalize_utilities()

end program test_random

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
