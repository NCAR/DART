! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program integrate_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Program to integrate assimilation model forward for asynchronous filter
! execution.

use types_mod,           only : r8
use time_manager_mod,    only : time_type, operator(<), print_time
use utilities_mod,       only : initialize_utilities, register_module,              &
                                error_handler, E_MSG, E_ERR, timestamp
                                
use assim_model_mod,     only : static_init_assim_model, get_model_size,              &
                                open_restart_read, open_restart_write, close_restart, &
                                awrite_state_restart, aread_state_restart

use obs_model_mod,        only : advance_state

use ensemble_manager_mod, only : init_ensemble_manager, put_copy, ensemble_type, get_copy

use mpi_utilities_mod,    only : initialize_mpi_utilities, finalize_mpi_utilities, &
                                 task_count, iam_task0


implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(time_type)         :: target_time
integer                 :: iunit, model_size
type(ensemble_type)     :: ens_handle

! dummy args for the advance_state call.  presumably this 
! executable was invoked by a script that was originally 
! called by advance_state() in filter, so no need to get 
! recursive here.
!character (len=129)     :: adv_ens_command = ''
!integer                 :: async = 0

! for debugging
logical                 :: trace_execution = .false.

! for speed, accuracy write model_advance files in binary
! both ways.  for debugging make this 'formatted' instead
! of 'unformatted' and you can see what's in the ud file.
character(len = 32) :: advance_restart_format = 'unformatted'

! Input and output filenames are hardcoded at this point.
character(len = 7) :: ic_file_name = "temp_ic", ud_file_name = 'temp_ud'

! NO ONE READS THIS IN RIGHT NOW - it's just a suggested start.
!namelist /integrate_model_nml/ ic_file_name, ud_file_name, &
!                               advance_restart_format,     &
!                               async, adv_ens_command,     &
!                               target_time, trace_execution
!----------------------------------------------------------------

! This program should only be run with a single process
call initialize_mpi_utilities('integrate_model')

! This version is ok to build with MPI and run with more than a single
! task, but only task 0 reads the input, advances the model, and writes 
! the output.  All the other tasks just hang out and exit when task 0 is done.
! FIXME: that could be changed to do multiple ens members in parallel.
! the code does NOT do this now.
if(task_count() > 1) &
   call error_handler(E_MSG,'integrate_model','Only one process doing the work', &
   source,revision,revdate)

call register_module(source,revision,revdate)

if (trace_execution) write(*,*) 'inside integrate_model executable'

! Initialize the model class data
call static_init_assim_model()
model_size = get_model_size()

! Initialize an ensemble manager type with a single copy
call init_ensemble_manager(ens_handle, num_copies=1, num_vars=model_size)

if (iam_task0()) then
   !------------------- Read restart from file ----------------------
   if (trace_execution) write(*,*) 'ready to open input restart file ', trim(ic_file_name)

   iunit = open_restart_read(ic_file_name)

   if (trace_execution) write(*,*) 'opened, iunit = ', iunit

   ! Read in the target time - could make a namelist item that overrides this.
   call aread_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit, target_time)

   if (trace_execution) write(*,*) 'time of data, advance-to are:'
   if (trace_execution) call print_time(ens_handle%time(1))
   if (trace_execution) call print_time(target_time)

   call close_restart(iunit)
   !-----------------  Restart read in --------------------------------

   ! Advance this state to the target time
   ! If the model time is past the obs set time, just need to skip
   if (trace_execution) write(*,*) 'calling advance_state if needed'

   if(ens_handle%time(1) < target_time) &
      call advance_state(ens_handle, ens_size=1, target_time=target_time, &
         async=0, adv_ens_command='', tasks_per_model_advance=1)

   ! Output the restart file if requested; Force to binary for bitwise reproducing
   ! use in filter and perfect_model obs with shell advance options
   if (trace_execution) write(*,*) 'ready to open output restart file ', trim(ud_file_name)

   iunit = open_restart_write(ud_file_name, advance_restart_format)

   if (trace_execution) write(*,*) 'opened, iunit = ', iunit

   call awrite_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)

   if (trace_execution) write(*,*) 'time of data after advance:'
   if (trace_execution) call print_time(ens_handle%time(1))

   call close_restart(iunit)
endif

if (trace_execution) write(*,*) 'end of integrate_model executable'
call finalize_mpi_utilities()

end program integrate_model
