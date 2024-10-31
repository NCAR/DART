! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_kde_mod
   use  utilities_mod, only : initialize_utilities, finalize_utilities
   use kde_distribution_mod, only : test_kde
   call initialize_utilities()
   call test_kde()
   call finalize_utilities()

end program test_kde_mod
