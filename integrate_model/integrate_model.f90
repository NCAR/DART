program integrate_model

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

! Program to integrate assimilation model forward for asynchronous filter
! execution.

use types_mod
use time_manager_mod, only : time_type, set_time, print_time, operator(/=), &
   operator(>), operator(<), read_time
use utilities_mod,    only :  get_unit, open_file, close_file, check_nml_error, &
   file_exist
use assim_model_mod,  only : assim_model_type, static_init_assim_model, &
   get_model_size, get_initial_condition, get_closest_state_time_to, &
   advance_state, set_model_time, get_model_time, init_diag_output, &
   output_diagnostics, init_assim_model, get_state_vector_ptr, &
   write_state_restart, read_state_restart, &
   init_diag_outputORG, output_diagnosticsORG

implicit none

type(time_type)         :: time, target_time
type(assim_model_type) :: x(1)
integer :: unit, ierr, io, model_size

!----------------------------------------------------------------
! Namelist input with default values
!
integer :: target_time_days = -1, target_time_seconds = -1
character(len = 129) :: ic_file_name = " ", &
                        ud_file_name = ' '

namelist /integrate_model_nml/ target_time_days, target_time_seconds, &
   ic_file_name, ud_file_name
!----------------------------------------------------------------

! Begin by reading the namelist input
if(file_exist('integrate_model.nml')) then
   unit = open_file(file = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(unit, nml = integrate_model_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'integrate_model_nml')
   enddo
 11 continue
   call close_file(unit)
endif

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

!------------------- Read restart from file ----------------------
unit = get_unit()
!!!open(unit = unit, file = restart_in_file_name)
open(unit = unit, file = 'temp_ic')
! Read in the target time
target_time = read_time(unit)
call init_assim_model(x(1))
call read_state_restart(x(1), unit)
close(unit)
!-----------------  Restart read in --------------------------------

! Advance this state to the target time (which comes from namelist)
! If the model time is past the obs set time, just need to skip
call print_time(target_time, 'target time is')
call print_time(get_model_time(x(1)), 'model time is')
if(get_model_time(x(1)) < target_time) then
   call advance_state(x, 1, target_time, .false.)
endif

! Output the restart file if requested
unit = get_unit()
!!!open(unit = unit, file = ud_file_name)
open(unit = unit, file = 'temp_ud')
call write_state_restart(x(1), unit)
close(unit)

end program integrate_model
