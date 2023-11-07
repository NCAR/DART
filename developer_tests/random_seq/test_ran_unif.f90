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
integer :: seq

call initialize_utilities('test_ran_unif')

call init_random_seq(r, 5)

seq = ran_twist(r)
   
if (seq == 953453411) then
   write(*,*) 'ran_unif test: PASS'
else
   write(*,*) 'ran_unif test: FAIL'
endif

call finalize_utilities()

end program test_ran_unif

