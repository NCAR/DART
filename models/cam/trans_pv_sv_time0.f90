! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program trans_pv_sv_time0

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read CAM 'initial' file (netCDF format).
!         Reform fields into a state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit
use        model_mod, only : model_type, init_model_instance, read_cam_init, &
                             prog_var_to_vector
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size , set_model_state_vector, write_state_restart, &
   set_model_time, open_restart_read, open_restart_write, close_restart, &
   aread_state_restart
! Guam; move time stripping from advance_model to here
use time_manager_mod, only : time_type, read_time

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
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
real(r8), allocatable  :: x_state(:), x_temp(:)
integer                :: file_unit, x_size

! Static init assim model calls static_init_model
PRINT*,'static_init_assim_model in trans_pv_sv'
call static_init_assim_model()

! Initialize the assim_model instance
call init_assim_model(x)

! Allocate the local state vector
x_size = get_model_size()
allocate(x_state(x_size), x_temp(x_size))

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

! Read the file cam state fragments into var
call read_cam_init(file_name, var)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

! Put this in the structure
call set_model_state_vector(x, x_state)

! Guam; move time stripping from advance_model into here
file_unit = open_restart_read(file_time)
call aread_state_restart(model_time, x_temp, file_unit, adv_to_time)
! for these time 0 _ud files we want the model time (0) to be written,
! not the target time
call set_model_time (x, model_time)
call close_restart(file_unit)

! Get channel for output 
! debug file_unit = 13
file_unit = open_restart_write(file_out)
PRINT*,'In trans_pv_sv file_out unit = ',file_unit
PRINT*,' '
! write out state vector in "proprietary" format
call write_state_restart(x, file_unit)
call close_restart(file_unit)

end program trans_pv_sv_time0
