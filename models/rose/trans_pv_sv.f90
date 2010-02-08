! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

 
program trans_pv_sv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between ROSE and DART
!
! method: Read ROSE restart file (binary format).
!         Reform fields into a state vector.
!         Write out state vector in "proprietary" format for DART
!
!         based on trans_pv_sv for CAM
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit
use        model_mod, only : model_type, init_model_instance, read_ROSE_restart, &
                             prog_var_to_vector
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size , set_model_state_vector, write_state_restart, &
   set_model_time, open_restart_read, open_restart_write, close_restart, &
   aread_state_restart
use time_manager_mod, only : time_type

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

character (len = 128) ::  &
                          file_name = 'rose_restart.nc', & 
                          file_out  = 'temp_ud'

! Temporary allocatable storage to read in a native format for ROSE state
type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: model_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size

! Static init assim model calls static_init_model
PRINT*,'static_init_assim_model in trans_pv_sv'
call static_init_assim_model()

! Initialize the assim_model instance
call init_assim_model(x)

! Allocate the local state vector
x_size = get_model_size()
allocate(x_state(x_size))

! Allocate the instance of the ROSE model type for storage
call init_model_instance(var)

! Read the file ROSE state fragments into var
call read_ROSE_restart(file_name, var, model_time)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

call set_model_state_vector(x, x_state)

call set_model_time(x, model_time)

file_unit = open_restart_write(file_out)
PRINT*,'In trans_pv_sv file_out unit = ',file_unit
PRINT*,' '
! write out state vector in "proprietary" format
call write_state_restart(x, file_unit)
call close_restart(file_unit)

end program trans_pv_sv
