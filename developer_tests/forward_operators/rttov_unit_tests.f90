program rttov_unit_tests

use obs_def_rttov_mod, only : set_visir_metadata
use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities

implicit none


! DART initialization
call initialize_utilities('rttov_unit_tests')

call test_set_visir_metadata



contains

subroutine test_set_visir_metadata

! test of set_visir_metadata

integer  :: key
real(r8) :: sat_az, sat_ze, sun_az, sun_ze
integer  :: platform_id, sat_id, sensor_id, channel
real(r8) :: specularity 

call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
      platform_id, sat_id, sensor_id, channel, specularity)

end subroutine test_set_visir_metadata


end program rttov_unit_tests
