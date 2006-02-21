! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program assim_region

! Main program to do assimilation of a region with filters for
! parallel multiple executables
!
! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

use  assim_tools_mod, only : async_assim_region, assim_tools_init
use    utilities_mod, only : initialize_utilities, register_module, timestamp
use  assim_model_mod, only : static_init_assim_model
use obs_sequence_mod, only : static_init_obs_sequence

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! Need to initialize modules used as appropriate
call initialize_utilities('assim_region')
call register_module(source, revision, revdate)
call assim_tools_init(dont_read_restart = .true.)
call static_init_obs_sequence()
call static_init_assim_model()

write(*, *) 'program assim region calling async_assim_region'
call async_assim_region('filter_assim_region_in', 'filter_assim_region_out')
write(*, *) 'program assim region back from async_assim_region'

call timestamp(source,revision,revdate,'end') ! closes the log file.

end program assim_region
