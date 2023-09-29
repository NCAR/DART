! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program test_ran_unif

! test the uniform random number generator routine

use      types_mod, only : r8
use  utilities_mod, only : register_module, &
                           open_file, close_file, &
                           initialize_utilities, finalize_utilities, &
                           squeeze_out_blanks
use random_seq_mod, only : random_seq_type, init_random_seq, ran_twist

implicit none


type (random_seq_type) :: r
integer :: j, i, n, f, seq
logical  :: write_this_one
character(len=50) :: formf = '(I12,4(F16.5))'
character(len=256) :: fname, temp

logical :: write_me = .true.       ! if true, write each distribution into a file for later diagnostics
integer :: write_limit = 1000000   ! but only if rep count is not greater than this limit

! to add more tests or change the parameters, specify a test count
! and update the sets of inputs here:

integer, parameter :: ntests = 1

call initialize_utilities('test_ran_unif')

write(*, *) '' 

do j=1, ntests

   call init_random_seq(r, 5)

   n = j

   ! save all values in a file for post-plotting?
   write_this_one = (write_me .and. n <= write_limit)

   if (write_this_one) then
      write(temp, "(A,I10)") "ran_unif_", n
      call squeeze_out_blanks(temp, fname)
      f = open_file(fname)
   endif
   
   seq = ran_twist(r)
   
   if (seq == 953453411) then
      write(*,*) 'ran_unif test: PASS'
      if (write_this_one) write(f,*) 'PASS'
   else
      write(*,*) 'ran_unif test: FAIL'
      if (write_this_one) write(f,*) 'FAIL'
   endif

   if (write_this_one) call close_file(f)

enddo

call finalize_utilities()

end program test_ran_unif

