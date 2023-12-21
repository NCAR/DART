! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_table_read

use algorithm_info_mod, only : init_algorithm_info_mod, end_algorithm_info_mod
use utilities_mod,      only : initialize_utilities, finalize_utilities

implicit none

call initialize_utilities('test_table_read')

call init_algorithm_info_mod()
call end_algorithm_info_mod()

call finalize_utilities()

end program test_table_read
