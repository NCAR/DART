program assim_region

! Main program to do assimilation of a region with filters for
! parallel multiple executables

use assim_tools_mod, only : async_assim_region, assim_tools_init
use utilities_mod, only : initialize_utilities
use assim_model_mod, only : static_init_assim_model
use obs_sequence_mod, only : static_init_obs_sequence

implicit none

! Need to initialize modules used as appropriate
call initialize_utilities
!  call register_module(source, revision, revdate)
call assim_tools_init()
call static_init_obs_sequence()
call static_init_assim_model()

write(*, *) 'programm assim region calling async_assim_region'
call async_assim_region('filter_assim_region_in', 'filter_assim_region_out')
write(*, *) 'programm assim region back from async_assim_region'

end program assim_region
