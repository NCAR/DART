! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_normal_distribution

use  utilities_mod, only : initialize_utilities, finalize_utilities
use normal_distribution_mod, only : test_normal

implicit none

call initialize_utilities()
call test_normal()
call finalize_utilities()

end program test_normal_distribution
