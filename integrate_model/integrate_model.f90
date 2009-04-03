! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program integrate_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Program to integrate assimilation model forward for asynchronous filter
! execution.

use types_mod,           only : r8
use time_manager_mod,    only : time_type, operator(<)
use utilities_mod,       only : initialize_utilities, register_module,              &
                                error_handler, E_MSG, E_ERR, timestamp
                                
use assim_model_mod,     only : static_init_assim_model, get_model_size,              &
                                open_restart_read, open_restart_write, close_restart, &
                                awrite_state_restart, aread_state_restart

use obs_model_mod,        only : advance_state

use ensemble_manager_mod, only : init_ensemble_manager, put_copy, ensemble_type, get_copy

use mpi_utilities_mod,    only : initialize_mpi_utilities, finalize_mpi_utilities, &
                                 task_count


implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(time_type)         :: target_time
integer                 :: iunit, model_size
type(ensemble_type)     :: ens_handle
character (len=129)     :: adv_ens_command = ''

!----------------------------------------------------------------
! Input and output filenames are hardcoded at this point.
!
character(len = 7) :: ic_file_name = "temp_ic", ud_file_name = 'temp_ud'

!----------------------------------------------------------------

! This program should only be run with a single process
call initialize_mpi_utilities('integrate_model')
if(task_count() > 1) &
   call error_handler(E_ERR,'integrate_model','Only use single process', &
   source,revision,revdate)

call register_module(source,revision,revdate)


! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

! Initialize an ensemble manager type with a single copy
call init_ensemble_manager(ens_handle, num_copies=1, num_vars=model_size)

!------------------- Read restart from file ----------------------
iunit = open_restart_read(ic_file_name)
! Read in the target time
call aread_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit, target_time)
call close_restart(iunit)
!-----------------  Restart read in --------------------------------

! Advance this state to the target time
! If the model time is past the obs set time, just need to skip
if(ens_handle%time(1) < target_time) &
   call advance_state(ens_handle, ens_size=1, target_time=target_time, &
      async=0, adv_ens_command=adv_ens_command, tasks_per_model_advance=1)

! Output the restart file if requested; Force to binary for bitwise reproducing
! use in filter and perfect_model obs with shell advance options
iunit = open_restart_write(ud_file_name, "unformatted")
call awrite_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)
call close_restart(iunit)

call timestamp(source,revision,revdate,'end') ! closes the log file.

call finalize_mpi_utilities()

end program integrate_model
