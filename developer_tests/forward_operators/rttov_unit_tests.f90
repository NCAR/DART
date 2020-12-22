! Runs the unit tests for rttov
!
program rttov_unit_tests

use obs_def_rttov_mod, only : test_unit_setup,         &
                              test_set_metadata, &
                              test_initializations,    &
                              test_key_get_expected,   &
                              test_unit_teardown
use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities
use     assert_mod,    only : assert_equal 

implicit none

integer :: metadata_size(3)

! DART initialization
call initialize_utilities('rttov_unit_tests')

! Unit test initialization
if ( test_unit_setup(1) ) then ! metadata start from 1 and grow

   metadata_size = test_set_metadata(10,10) ! visir, mw
   
   call assert_equal(metadata_size(1),32, 'obstype_metadata')
   call assert_equal(metadata_size(2),16, 'visir_obs_metadata')
   call assert_equal(metadata_size(3),16, 'mw_obs_metadata')
   
   !call test_key_get_expected
   
   !call test_initializations
else  ! module is already initialized, unit tests are not reliable
   print*, 'FAIL'
endif

call test_unit_teardown()

if (test_unit_setup(3) ) then 

   metadata_size = test_set_metadata(2,7) ! visir, mw
   
   call assert_equal(metadata_size(1),12, 'obstype_metadata')
   call assert_equal(metadata_size(2),3, 'visir_obs_metadata')
   call assert_equal(metadata_size(3),12, 'mw_obs_metadata')

endif


end program rttov_unit_tests
