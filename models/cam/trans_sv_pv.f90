program trans_sv_pv

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read DART state vector ("proprietary" format)
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use types_mod
use model_mod, only : model_type, init_model_instance, write_cam_init, &
   vector_to_prog_var
use assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size, get_model_state_vector, read_state_restart
use utilities_mod, only : get_unit

type(assim_model_type) :: x
type(model_type) :: var
real(r8), allocatable :: x_state(:)
integer :: file_unit, x_size
character (len = 128) :: file_name = 'CAM_FILE.nc', file_in = 'temp_ic'

! Static init assim model calls static_init_model
call static_init_assim_model()
call init_assim_model(x)

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

! Get file for DART vector input
file_unit = get_unit()
open(unit = file_unit, file = file_in)

! read in state vector from DART
call read_state_restart(x, file_unit)
close(file_unit)

! Get the state part of the assim_model type x
x_size = get_model_size()
allocate(x_state(x_size))
x_state = get_model_state_vector(x)

! decompose vector back into CAM fields
call vector_to_prog_var (x_state, var)
deallocate (x_state)

! write fields to the netCDF initial file
call write_cam_init(file_name, var)

end program trans_sv_pv
