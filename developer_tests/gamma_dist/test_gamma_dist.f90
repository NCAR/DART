! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_gamma_distribution

use  utilities_mod, only : initialize_utilities, finalize_utilities
use gamma_distribution_mod, only : test_gamma

implicit none

call initialize_utilities()
call test_gamma()
call finalize_utilities()

end program test_gamma_distribution
