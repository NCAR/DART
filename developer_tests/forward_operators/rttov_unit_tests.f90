! Runs the unit tests for rttov
!
program rttov_unit_tests

use obs_def_rttov_mod, only : test_unit_setup,         &
                              test_set_metadata, &
                              test_initializations,    &
                              test_key_get_expected
use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities
use     assert_mod,    only : assert_equal 

implicit none

integer :: metadata_size(4)

! DART initialization
call initialize_utilities('rttov_unit_tests')

! Unit test initialization
call test_unit_setup

! metadata start from 1 and grow

  metadata_size = test_set_metadata(10,10)

  call assert_equal(metadata_size(1),10, 'one')
  call assert_equal(metadata_size(2),10, 'two')
  call assert_equal(metadata_size(3),10, 'three') 
  call assert_equal(metadata_size(4),10, 'four') 

!call test_key_get_expected

!call test_initializations

end program rttov_unit_tests
