! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program gen_init

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time, get_time, print_time
use  assim_model_mod, only : static_init_assim_model, get_model_size, &
   aget_initial_condition, get_model_state_vector,awrite_state_restart, &
   open_restart_read, open_restart_write, close_restart
use utilities_mod,    only : open_file, file_exist, get_unit, close_file, &
                             initialize_utilities, register_module, error_handler, &
                             E_ERR, E_WARN, E_MSG, E_DBG, timestamp

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical                 :: output_restart
integer                 :: model_size
integer                 :: iunit
type(time_type)         :: time1
real(r8), allocatable   :: ens(:)
character(len=129)      :: restart_out_file_name = "member1_ics"

logical                         :: allocate_wrf = .true.

! Some static initialization
call gen_init_modules_used

! Initialize the model now that obs_sequence is all set up
model_size = get_model_size()

! Allocate storage for doing advance with ensemble based tools
allocate(ens(model_size))

! Grab random IC
call aget_initial_condition(time1, ens)

! Output a restart file
iunit = open_restart_write(restart_out_file_name)
call awrite_state_restart(time1, ens, iunit)
call close_restart(iunit)

contains

!=====================================================================

subroutine gen_init_modules_used()

! Initialize modules used that require it
call initialize_utilities
call register_module(source,revision,revdate)
call error_handler(E_MSG,'gen_init','STARTING',source,revision,revdate)

! Initialize anything else?
call static_init_assim_model()

end subroutine gen_init_modules_used

end program gen_init
