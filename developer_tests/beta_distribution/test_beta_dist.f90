! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_beta_dist

use  utilities_mod, only : initialize_utilities, finalize_utilities
use beta_distribution_mod,  only : test_beta

call initialize_utilities()

call test_beta()

end program test_beta_dist
