! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program model_to_dart

!----------------------------------------------------------------------
! purpose: interface between ROSE and DART
!
! method: Read ROSE restart file (netCDF format).
!         Reform fields into a DART state vector.
!         Write out state vector in "proprietary" format for DART
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, initialize_utilities, finalize_utilities, &
                             error_handler, E_MSG
use        model_mod, only : model_type, static_init_model, get_model_size, &
                             init_model_instance, read_ROSE_restart, &
                             prog_var_to_vector 
use  assim_model_mod, only : open_restart_write, awrite_state_restart, close_restart
use time_manager_mod, only : time_type

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character (len = 128) ::  &
   file_name = 'rose_restart.nc', & 
   file_out  = 'temp_ud'

! Temporary allocatable storage to read in a native format for ROSE state
type(model_type)       :: var
type(time_type)        :: model_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size

call initialize_utilities(progname='model_to_dart', output_flag=.true.)

! static_init_model reads input.nml, sets the geometry, model size, etc.
call static_init_model()

! Allocate the local state vector
x_size = get_model_size()
allocate(x_state(x_size))

! Allocate an empty instance of the ROSE model type for storage
! This is needed because read_ROSE_restart() requires it.  
call init_model_instance(var)

! Read the ROSE state variables into var and set the model_time
! to reflect the valid time of the ROSE state.
call read_ROSE_restart(file_name, var, model_time)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

! write out state vector in "proprietary" format
file_unit = open_restart_write(file_out)
call awrite_state_restart(model_time, x_state, file_unit)
call close_restart(file_unit)

call error_handler(E_MSG,'model_to_dart','Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
