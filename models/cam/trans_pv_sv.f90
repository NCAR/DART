! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program trans_pv_sv

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
   binary_restart_files, set_model_time
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
character (len = 16)  :: file_form

! Temporary allocatable storage to read in a native format for cam state
type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: adv_to_time
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

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

! Read the file cam state fragments into var
call read_cam_init(file_name, var)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

! Put this in the structure
call set_model_state_vector(x, x_state)

! get form of file output from assim_model_mod
if (binary_restart_files ) then
   file_form = 'unformatted'
else
   file_form = 'formatted'
endif
PRINT*,'In trans_pv_sv binary_restart_files, file_form = ',binary_restart_files, file_form

! Guam; move time stripping from advance_model into here
file_unit = get_unit()
open(unit = file_unit, file = file_time, form=file_form)
adv_to_time = read_time(file_unit, file_form)
call set_model_time (x, adv_to_time)
close(file_unit)

! Get channel for output 
! debug file_unit = 13
file_unit = get_unit()
open(unit = file_unit, file = file_out, form=file_form)
PRINT*,'In trans_pv_sv file_out unit = ',file_unit
PRINT*,' '
! write out state vector in "proprietary" format
call write_state_restart(x, file_unit, file_form)
close(file_unit)

end program trans_pv_sv
