program rttov_unit_tests

use obs_def_rttov_mod, only : test_set_visir_metadata, test_initializations
use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities

implicit none


! DART initialization
call initialize_utilities('rttov_unit_tests')

call test_set_visir_metadata

call test_initializations


end program rttov_unit_tests
