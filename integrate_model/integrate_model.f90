! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program integrate_model

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

! Program to integrate assimilation model forward for asynchronous filter
! execution.

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time, print_time, operator(/=), &
                             operator(>), operator(<), read_time
use    utilities_mod, only : get_unit, open_file, close_file, check_nml_error, &
                             file_exist
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   get_model_size, get_initial_condition, get_closest_state_time_to, &
   advance_state, set_model_time, get_model_time, init_diag_output, &
   output_diagnostics, init_assim_model, get_state_vector_ptr, &
   write_state_restart, read_state_restart, binary_restart_files

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(time_type)         :: time, target_time
type(assim_model_type) :: x(1)
integer :: iunit, ierr, io, model_size

character (len=129)    :: adv_ens_command = ''
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
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = integrate_model_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'integrate_model_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

!------------------- Read restart from file ----------------------
iunit = get_unit()
! Read in the target time
if ( binary_restart_files ) then
!!!open(unit = iunit, file = restart_in_file_name, form = "unformatted")
   open(unit = iunit, file = 'temp_ic', form = "unformatted")
   target_time = read_time(iunit, 'unformatted')
   call init_assim_model(x(1))
   call read_state_restart(x(1), iunit, 'unformatted')
else
!!!open(unit = iunit, file = restart_in_file_name)
   open(unit = iunit, file = 'temp_ic')
   target_time = read_time(iunit)
   call init_assim_model(x(1))
   call read_state_restart(x(1), iunit)
endif
close(iunit)
!-----------------  Restart read in --------------------------------

! Advance this state to the target time (which comes from namelist)
! If the model time is past the obs set time, just need to skip
call print_time(target_time, 'target time is')
call print_time(get_model_time(x(1)), 'model time is')
if(get_model_time(x(1)) < target_time) then
   call advance_state(x, 1, target_time, 0, adv_ens_command)
endif

! Output the restart file if requested
iunit = get_unit()
if ( binary_restart_files ) then
!!!open(unit = iunit, file = ud_file_name, form = "unformatted")
   open(unit = iunit, file = 'temp_ud', form = "unformatted")
   call write_state_restart(x(1), iunit, 'unformatted')
else
!!!open(unit = iunit, file = ud_file_name)
   open(unit = iunit, file = 'temp_ud')
   call write_state_restart(x(1), iunit)
endif
close(iunit)

end program integrate_model
