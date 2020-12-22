! Runs the unit tests for rttov
!
program rttov_unit_tests

use obs_def_rttov_mod, only : test_unit_setup,    &
                              test_set_metadata,  &
                              test_unit_teardown, &
                              test_metadata
use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities
use     assert_mod,    only : assert_equal 

implicit none

integer :: metadata_size(3)

! metadata SUBTYPE, SUBKEY test 1
integer :: correct1_subtype(32) = (/1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/) 
integer :: correct1_subkey(32) = (/1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
integer :: m1(2,32)

! metadata SUBTYPE, SUBKEY test 2
integer :: correct2_subtype(12) = (/1,1,2,2,2,2,2,2,2,-1,-1,-1/) 
integer :: correct2_subkey(12) = (/1,2,1,2,3,4,5,6,7,-1,-1,-1/)
integer :: m2(2,12)

! DART initialization
call initialize_utilities('rttov_unit_tests')

! Unit test initialization
if ( test_unit_setup(1) ) then ! metadata start from 1 and grow

   metadata_size = test_set_metadata(10,10) ! visir, mw
   
   call assert_equal(metadata_size(1),32, 'obstype_metadata')
   call assert_equal(metadata_size(2),16, 'visir_obs_metadata')
   call assert_equal(metadata_size(3),16, 'mw_obs_metadata')

   call test_metadata(m1)
   call assert_equal(m1(1,:), correct1_subtype, 'subtypes')
   call assert_equal(m1(2,:), correct1_subkey, 'subtypes')
 
else  ! module is already initialized, unit tests are not reliable
   print*, 'FAIL'
endif

call test_unit_teardown()

if ( test_unit_setup(3) ) then 

   metadata_size = test_set_metadata(2,7) ! visir, mw
   
   call assert_equal(metadata_size(1),12, 'obstype_metadata')
   call assert_equal(metadata_size(2),3, 'visir_obs_metadata')
   call assert_equal(metadata_size(3),12, 'mw_obs_metadata')

   call test_metadata(m2)
   call assert_equal(m2(1,:), correct2_subtype, 'subtypes')
   call assert_equal(m2(2,:), correct2_subkey, 'subtypes')

else  ! module is already initialized, unit tests are not reliable
   print*, 'FAIL'
endif


end program rttov_unit_tests
