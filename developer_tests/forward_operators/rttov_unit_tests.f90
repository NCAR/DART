! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! Unit tests for rttov
!
! These unit tests are best run with a TERMLEVEL of 3, which
! will not die with a DART error.  This is to allow multiple
! calls to module_initialize in obs_def_rttov_mod.  
!
! If any of the unit tests are unable to start, the error code from 
! rttov_unit_tests is 102.  This is to give an error for 
! test_dart.csh to detect.  
!
! Tests: 
!   metadata 
!    PASS: metadata arrays grow correctly as observations are added
!          metadata arrays contain correct data
!    FAIL: incorrect metadata length
!          incorrect data
!
program rttov_unit_tests

use obs_def_rttov_mod, only : test_unit_setup,       &
                              test_set_metadata,     &
                              test_unit_teardown,    &
                              test_metadata,         &
                              test_key_within_range, &
                              test_subkey_within_range
use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities
use     assert_mod,    only : assert_equal 

implicit none

logical :: failme = .false.
integer :: metadata_size(3)

! -----------------------------------------------
! PASS correct answers for tests:
! metadata SUBTYPE, SUBKEY test 1: start MAXrttovkey=1 then add 10 visir obs and 10 mw obs
integer :: correct1_subtype(32) = (/1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/) 
integer :: correct1_subkey(32) = (/1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
integer :: m1(2,32)

! metadata SUBTYPE, SUBKEY test 2: start MAXrttovkey=3 then add 2 visir obs and add 7 mw obs
integer :: correct2_subtype(12) = (/1,1,2,2,2,2,2,2,2,-1,-1,-1/) 
integer :: correct2_subkey(12) = (/1,2,1,2,3,4,5,6,7,-1,-1,-1/)
integer :: m2(2,12)

! key within range test: 2 visir obs, 7 mw obs, MAXrttovkey=100 (9 obs)
integer :: inkeys1(5) = (/-1,100,101,6,2/)  
logical :: correct_range1(5) = (/.false., .false., .false., .true., .true./)

! subkey wthin range tests: 2 visir obs, 7 mw obs 
integer :: inkeys2(5) = (/-1,2,7,8,21/)  
!    test 1: SUBYTPE=NO_OBS
logical :: correct_range2_no_subkey(5) = (/.false., .false., .false., .false., .false./)
!    test 2: SUBYTPE=VISIR
logical :: correct_range2_visir(5) = (/.false., .true., .false., .false., .false./)
!    test 3: SUBYTPE=MW
logical :: correct_range2_mw(5) = (/.false., .true., .true., .false., .false./)
logical :: in_range(5)
 
! -----------------------------------------------

! DART initialization
call initialize_utilities('rttov_unit_tests')

! -----------------------------------------------
! Metadata test 1
! -----------------------------------------------
! Unit test initialization
! MAXrttovkey starts from 1 and metadata grows
if ( test_unit_setup(1) ) then

   metadata_size = test_set_metadata(10,10) ! number of visir, number of mw

   ! test sizes   
   call assert_equal(metadata_size(1),32, 'obstype_metadata')
   call assert_equal(metadata_size(2),16, 'visir_obs_metadata')
   call assert_equal(metadata_size(3),16, 'mw_obs_metadata')

   ! test contents
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
! MAXrttovkey starts from 3 and metadata grows
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

call test_unit_teardown()

! -----------------------------------------------
! key within range
! -----------------------------------------------
! MAXrttovkey=100
if ( test_unit_setup(100) ) then  

  ! KEY
  metadata_size = test_set_metadata(2,7) ! visir, mw
  call test_key_within_range(inkeys1, in_range) 
  call assert_equal(in_range, correct_range1, 'valid keys')

  ! NO SUBTYPE 
  call test_subkey_within_range(inkeys2, 0, in_range)
  call assert_equal(in_range, correct_range2_no_subkey, 'no subtype')

  ! VISIR SUBKEY
  call test_subkey_within_range(inkeys2, 1, in_range)
  call assert_equal(in_range, correct_range2_visir, 'valid visir keys')

  ! MW SUBKEY
  call test_subkey_within_range(inkeys2, 2, in_range)
  call assert_equal(in_range, correct_range2_mw, 'valid mw keys')

else  ! module is already initialized, unit tests are not reliable
   print*, 'FAIL: rrtov module is already initialized, unit tests are not reliable'
   failme = .true.
endif

if (failme) call exit(102)

end program rttov_unit_tests
