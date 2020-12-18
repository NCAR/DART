! Runs the unit tests for rttov
!
program rttov_unit_tests

use obs_def_rttov_mod, only : test_unit_setup,         &
                              test_set_metadata, &
                              test_initializations,    &
                              test_key_get_expected
use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities

implicit none


! DART initialization
call initialize_utilities('rttov_unit_tests')

! Unit test initialization
call test_unit_setup

call test_set_metadata

call test_key_get_expected

!call test_initializations


end program rttov_unit_tests
