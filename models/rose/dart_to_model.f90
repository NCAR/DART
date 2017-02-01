! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_model

!----------------------------------------------------------------------
! purpose: interface between ROSE and DART
!
! method: Read DART state vector ("proprietary" format)
!         Reform state vector back into ROSE fields.
!         Replace those fields on the ROSE restart file with the new values,
!         Replace the 'mtime' variable in the ROSE restart file with
!         the 'valid time' of the DART state vector.
!         Write a new 'ROSE_NML' namelist in file 'rose.nml'.
!
!         Compiler note: Rather curiously, the PG compiler reads 
!         a namelist called 'rose_nml' and then, when writing the
!         namelist - uses uppercase 'ROSE_NML'. Then, the next time
!         you need to read the 'rose_nml' ... it fails! So - we 
!         adopted the uppercase convention from the get-go. TJH
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, initialize_utilities, E_ERR, &
                             error_handler, finalize_utilities, E_MSG
use        model_mod, only : model_type, get_model_size, init_model_instance, &
                             vector_to_prog_var, update_ROSE_restart, &
                             update_ROSE_namelist, static_init_model
use  assim_model_mod, only : assim_model_type, aread_state_restart, &
                             open_restart_read, close_restart
use time_manager_mod, only : time_type, read_time

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(model_type)       :: var
type(time_type)        :: model_time, adv_to_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size, ens_member, io
character (len = 128)  :: file_name = 'rose_restart.nc', file_in = 'temp_ic'

!----------------------------------------------------------------------
! This program has one input argument that is read from STDIN ... 
!----------------------------------------------------------------------

read(*, *, iostat = io )  ens_member
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot read ens_member from STDIN', &
         source,revision,revdate)
endif

!----------------------------------------------------------------------
!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_model', output_flag=.true.)

call static_init_model()        ! reads input.nml, etc., sets the table 
x_size = get_model_size()       ! now that we know how big state vector is ...
allocate(x_state(x_size))       ! allocate space for the (empty) state vector

! Open the DART model state ... 
! Read in the time to which ROSE must advance.  
! Read in the valid time for the model state
! Read in state vector from DART

file_unit = open_restart_read(file_in)

call aread_state_restart(model_time, x_state, file_unit, adv_to_time)
call close_restart(file_unit)

! Parse the vector into ROSE fields (prognostic variables)
call init_model_instance(var, model_time)
call vector_to_prog_var(x_state, var)
deallocate(x_state)

! write fields to the binary ROSE restart file
call update_ROSE_restart(file_name, var)
call update_ROSE_namelist('rose.nml', model_time, adv_to_time, ens_member)

call error_handler(E_MSG,'dart_to_model','Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program dart_to_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
