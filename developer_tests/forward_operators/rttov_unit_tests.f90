! Unit tests for rttov
!
! These unit tests use warn_only_utilities_mod.f90 which 
! will not die with a DART error.  This is to allow multiple
! calls to module_initilize in obs_def_rttov_mod.  
!
! If any of the unit tests fail, the error code from rttov_unit_tests
! is 102.  This is to give an error for test_dart.csh to detect.  
!
! Tests: 
!   metadata 
!    PASS: metadata arrays grow correctly as observations are added
!          metadata arrays contain correct data
!    FAIL: incorrect metadata length
!          incorrect data
!
!
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

logical :: failme = .true.


integer :: metadata_size(3)

! -----------------------------------------------
! PASS correct answers for tests:
! metadata SUBTYPE, SUBKEY test 1
!integer :: correct1_subtype(32) = (/1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/) 
integer :: correct1_subtype(32) = (/14,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/) 
integer :: correct1_subkey(32) = (/1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
integer :: m1(2,32)

! metadata SUBTYPE, SUBKEY test 2
integer :: correct2_subtype(12) = (/1,1,2,2,2,2,2,2,2,-1,-1,-1/) 
integer :: correct2_subkey(12) = (/1,2,1,2,3,4,5,6,7,-1,-1,-1/)
integer :: m2(2,12)
! -----------------------------------------------

! DART initialization
call initialize_utilities('rttov_unit_tests')

! -----------------------------------------------
! Metadata test 1
! -----------------------------------------------
! Unit test initialization
if ( test_unit_setup(1) ) then ! MAXrttovkey starts from 1 and metadata grows

   metadata_size = test_set_metadata(10,10) ! visir, mw
   
   call assert_equal(metadata_size(1),32, 'obstype_metadata')
   call assert_equal(metadata_size(2),16, 'visir_obs_metadata')
   call assert_equal(metadata_size(3),16, 'mw_obs_metadata')

   call test_metadata(m1)
   call assert_equal(m1(1,:), correct1_subtype, 'subtypes')
   call assert_equal(m1(2,:), correct1_subkey, 'subtypes')
 
else  ! module is already initialized, unit tests are not reliable
   print*, 'FAIL: rrtov module is already initialized, unit tests are not reliable'
   failme = .true.
endif

call test_unit_teardown()

! -----------------------------------------------
! Metadata test 2
! -----------------------------------------------
if ( test_unit_setup(3) ) then 

   metadata_size = test_set_metadata(2,7) ! visir, mw
   
   call assert_equal(metadata_size(1),12, 'obstype_metadata')
   call assert_equal(metadata_size(2),3, 'visir_obs_metadata')
   call assert_equal(metadata_size(3),12, 'mw_obs_metadata')

   call test_metadata(m2)
   call assert_equal(m2(1,:), correct2_subtype, 'subtypes')
   call assert_equal(m2(2,:), correct2_subkey, 'subtypes')

else  ! module is already initialized, unit tests are not reliable
   print*, 'FAIL: rrtov module is already initialized, unit tests are not reliable'
   failme = .true.
endif

if (failme) call exit(102)

end program rttov_unit_tests
