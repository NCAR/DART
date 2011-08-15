! DART software - Copyright 2004 - 2011 UCAR. This open source software iss
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between TIEGCM and DART
!
! method: Read TIEGCM restart file (netCDF format).
!         Reform fields into a DART state vector.
!         Write out state vector in "proprietary" format for DART
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, initialize_utilities, timestamp
use        model_mod, only : model_type, static_init_model, get_model_size, &
                             init_model_instance, read_TIEGCM_restart,      &
                             read_TIEGCM_secondary,                         &
                             read_TIEGCM_namelist, prog_var_to_vector 
use  assim_model_mod, only : open_restart_write, awrite_state_restart, close_restart
use time_manager_mod, only : time_type

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

character (len = 128) ::  &
   file_namelist = 'tiegcm.nml', &
   file_name1    = 'tiegcm_restart_p.nc', & 
   file_name2    = 'tiegcm_s.nc', & 
   file_out      = 'temp_ud'

! Temporary allocatable storage to read in a native format for TIEGCM state
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

! Allocate an empty instance of the TIEGCM model type for storage
! This is needed because read_TIEGCM_restart() requires it.  
call init_model_instance(var)

! Read the TIEGCM state variables into var and set the model_time
! to reflect the valid time of the TIEGCM state.
call read_TIEGCM_restart(file_name1, var, model_time)
call read_TIEGCM_secondary(file_name2, var)

! Read the TIEGCM input variables into var
call read_TIEGCM_namelist(file_namelist, var)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

! write out state vector in "proprietary" format
file_unit = open_restart_write(file_out)
call awrite_state_restart(model_time, x_state, file_unit)
call close_restart(file_unit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will also call finalize_utilities()
!----------------------------------------------------------------------
call timestamp(string1=source, pos='end')

end program model_to_dart
