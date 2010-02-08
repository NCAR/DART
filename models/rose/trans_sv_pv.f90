! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program trans_sv_pv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between ROSE and DART
!
! method: Read DART state vector ("proprietary" format)
!         Reform state vector back into ROSE fields.
!         Replace those fields on the ROSE restart file with the new values,
!         preserving all other information on the file.
!
!         based on prog_var_to_vector and vector_to_prog_var for CAM
!
!----------------------------------------------------------------------

use       types_mod, only : r8
use   utilities_mod, only : get_unit
use       model_mod, only : model_type, init_model_instance, &
   vector_to_prog_var, update_ROSE_restart 
use assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size, get_model_state_vector, read_state_restart, &
   open_restart_read, close_restart
use time_manager_mod, only : time_type, read_time

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: adv_to_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size
character (len = 128)  :: file_name = 'rose_restart.nc', file_in = 'temp_ic'

! Static init assim model calls static_init_model
PRINT*,'static_init_assim_model in trans_sv_pv'

call static_init_assim_model()
call init_assim_model(x)

! Allocate the instance of the rose model type for storage
call init_model_instance(var)

file_unit = open_restart_read(file_in)
PRINT*,'In trans_sv_pv file_in unit  = ',file_unit
PRINT*,' '

! Read in time to which ROSE must advance.  
! Neither this, nor time in x (x%time) is used in this program
! read in state vector from DART
call read_state_restart(x, file_unit, adv_to_time)
call close_restart(file_unit)

! Get the state part of the assim_model type x
x_size = get_model_size()
allocate(x_state(x_size))
PRINT*,'(trans_sv_pv) getting model state vector of length ',x_size
x_state = get_model_state_vector(x)

! decompose vector back into ROSE fields
PRINT*,'(trans_sv_pv) converting vector to prog_var'
call vector_to_prog_var (x_state, var)
deallocate (x_state)

! write fields to the binary ROSE restart file
PRINT*,'(trans_sv_pv) updating ',trim(file_name)
call update_ROSE_restart(file_name, var)

end program trans_sv_pv
