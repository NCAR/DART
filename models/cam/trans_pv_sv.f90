! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program trans_pv_sv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read CAM 'initial' file for model state, but not time (netCDF format).
!         Get target time from assim_model_state_ic (temp_ic).
!         Reform fields into a state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, file_exist, &
                             initialize_utilities, finalize_utilities
use        model_mod, only : model_type, init_model_instance, end_model_instance, &
                             prog_var_to_vector, read_cam_init
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size , set_model_state_vector, write_state_restart, &
   set_model_time, open_restart_read, open_restart_write, close_restart, &
   aread_state_restart
! Guam; move time stripping from advance_model to here
use time_manager_mod, only : time_type, read_time

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! character (len = 128) :: file_name = 'caminput.nc', file_out = 'temp_ic'
! Guam; move time stripping from script into here
character (len = 128) :: file_name = 'caminput.nc', file_out = 'temp_ud', &
                         file_time = 'temp_ic'

! Temporary allocatable storage to read in a native format for cam state
type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: model_time, adv_to_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size
logical                :: do_output = .false.

call initialize_utilities('Trans_pv_sv')

if(file_exist('element1')) do_output = .true.

! Static init assim model calls static_init_model
! which now (merge/MPI) calls read_cam_init)
call static_init_assim_model()

! Initialize the assim_model instance
call init_assim_model(x)

! Allocate the local state vector
x_size = get_model_size()
allocate(x_state(x_size))

! Allocate the instance of the cam model type for storage
! Nancy; why did we comment this out?  
!        Do it in read_cam_init?
!        What about end_model_instance?
! I'll just point to the space I need, not;    
call init_model_instance(var)

! Read the file cam state fragments into var, but not time
call read_cam_init(file_name, var)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

call end_model_instance(var)

! Put this in the structure
call set_model_state_vector(x, x_state)

! Integration of model was controlled by the restart file,
! so we use the target time of the restart file (from assim_model_state)
! as the current model state time.
file_unit = open_restart_read(file_time)
! We're done with x_state, so it can be uselessly filled in aread_state_restart,
! while getting model_time.
call aread_state_restart(model_time, x_state, file_unit, adv_to_time)
call set_model_time (x, adv_to_time)
call close_restart(file_unit)

! Get channel for output,
! write out state vector in "proprietary" format
file_unit = open_restart_write(file_out)
call write_state_restart(x, file_unit)
call close_restart(file_unit)

call finalize_utilities()

end program trans_pv_sv
