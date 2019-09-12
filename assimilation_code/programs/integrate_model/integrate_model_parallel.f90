! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Program to integrate assimilation model forward without assimilation.
!> Can be used for forecasts after an assimilation, spinning a model up
!> before assimilation, or for advancing subroutine-callable models from
!> a separate executable during an assimilation.
!>
!> Stop time can be set either by the 'advance_to' time in the restart file,
!> or by namelist control.  
!>
!> Usually only advances a single ensemble member, but with namelist control 
!> can cycle through multiple members.

program integrate_model

use time_manager_mod,    only : time_type, operator(<), print_time, get_time, &
                                set_time, set_time_missing
use utilities_mod,       only : register_module,  &
                                error_handler, E_MSG, nmlfileunit, &
                                do_nml_file, do_nml_term,          &
                                find_namelist_in_file, check_namelist_read
use assim_model_mod,     only : static_init_assim_model, get_model_size,              &
                                open_restart_read, open_restart_write, close_restart, &
                                awrite_state_restart, aread_state_restart
use obs_model_mod,        only : advance_state
use ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, &
                                 prepare_to_write_to_vars
use mpi_utilities_mod,    only : initialize_mpi_utilities, finalize_mpi_utilities, &
                                 task_count, iam_task0


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=128), parameter :: id  = "$Id$"

type(ensemble_type) :: ens_handle
type(time_type)     :: target_time, init_time
integer             :: iunit, model_size, rc, member
character(len = 32) :: write_format
logical             :: one_by_one, override_init
character(len=256)  :: ifile, ofile


! to overwrite the target time, set these to something >= 0
! via the namelist.  to overwrite the time in the restart file itself,
! set the init_time to >= 0 (less common).
integer :: target_days    = -1
integer :: target_seconds = -1
integer :: init_days      = -1
integer :: init_seconds   = -1

! args for the advance_state call.  if this executable was 
! invoked by a script that was originally called by advance_state() 
! in filter there is no need to get recursive here.
! if you want to use this outside of the dart advance 
! you could comment these in and add them to the namelist below.
character (len=129) :: adv_ens_command = ''
integer             :: async = 0
integer             :: tasks_per_model_advance = 1   ! unsupported for now
integer             :: ens_size = 1

logical :: single_restart_file_in  = .false.
logical :: single_restart_file_out = .false.

! for use outside dart, or to call after dart finishes and
! writes a restart file (without a target time), set this
! to .false. and supply a target time with target_days, secs
logical :: input_is_model_advance_file = .true.

! for debugging, status
logical             :: trace_execution = .false.
character(len=128)  :: errstring

! for speed, accuracy - write model_advance files in binary
! both ways.  for debugging make this 'formatted' instead
! of 'unformatted' and you can see what's in the ud file.
logical :: write_binary_restart_files = .true.

! Input and output filenames.  The defaults match what filter
! expects when using this to advance ensemble members in async 2.
character(len = 256) :: input_filename = "temp_ic"
character(len = 256) :: output_filename = "temp_ud"

! to enable the use of the namelist, change this to .true. and recompile.
logical :: use_namelist = .true.

! only read in if use_namelist is .true. -- ignored otherwise.
namelist /integrate_model_nml/ &
   ens_size, input_filename, output_filename, write_binary_restart_files, &
   trace_execution, input_is_model_advance_file, target_days, target_seconds, &
   single_restart_file_in, single_restart_file_out, init_days, init_seconds


!----------------------------------------------------------------

! This program should only be run with a single process
call initialize_mpi_utilities('integrate_model')

print *, 'in main program'

! This version is ok to build with MPI and run with more than a single
! task, but only task 0 reads the input, advances the model, and writes 
! the output.  All the other tasks just hang out and exit when task 0 is done.
! FIXME: that could be changed to do multiple ens members in parallel.
! the code does NOT do this now.
if(task_count() > 1) &
   call error_handler(E_MSG,'integrate_model','Only one process doing the work', &
   source,revision,revdate)

call register_module(id)

! this must come AFTER the standard utils are initialized.
! Read the integrate_model_nml namelist from input.nml if 'use_namelist' true.
if (use_namelist) then
   call find_namelist_in_file('input.nml', 'integrate_model_nml', iunit)
   read(iunit, nml = integrate_model_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "integrate_model_nml")
else
   errstring = ' !must edit integrate_model/integrate_model.f90 to enable this namelist'
   if (do_nml_file()) write(nmlfileunit, '(A)') trim(errstring)
   if (do_nml_term()) write(     *     , '(A)') trim(errstring)
endif

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=integrate_model_nml)
if (do_nml_term()) write(     *     , nml=integrate_model_nml)



if (trace_execution) write(*,*) 'inside integrate_model executable'

! Set default target time from namelist, if specified.
if (target_seconds >= 0 .and. target_days >= 0) then
   target_time = set_time(target_seconds, target_days)
else
   target_time = set_time_missing()
endif

! Set override for model state current time from namelist, if specified.
if (init_seconds >= 0 .and. init_days >= 0) then
   init_time = set_time(init_seconds, init_days)
   override_init = .true.
else
   init_time = set_time_missing()
   override_init = .false.
endif

! Initialize the model class data
call static_init_assim_model()
model_size = get_model_size()

! Figure out the output format string.  (Input format automatically detected.)
if (write_binary_restart_files) then
   write_format = "unformatted"
else
   write_format = "formatted"
endif

! Initialize the model so we can get the size.
call static_init_assim_model()
model_size = get_model_size()

! If either the restart in or out data is all together in a single file, 
! then there has to be room for the entire dataset in memory at once --
! which is the size of the state vector times the ensemble count.
one_by_one = .not. (single_restart_file_in .or. single_restart_file_out) 

! if both in and out are single files, we can loop here over all ens members
! doing each one piecemeal; no need for space for all members in memory.
if (one_by_one) then

   ! Initialize the ens manager with enough room for a single ensemble member.
   call init_ensemble_manager(ens_handle, num_copies=1, num_vars=model_size)
   call prepare_to_write_to_vars(ens_handle)

   do member=1, ens_size
 
      ! add member number as a suffix: e.g. base.0000
      write(ifile, "(a,a,i4.4)") trim(input_filename), '.', member
      write(ofile, "(a,a,i4.4)") trim(output_filename), '.', member

      !------------------- Read restart from file ----------------------
      if (trace_execution) write(*,*) 'ready to open input restart file ', trim(ifile)

      iunit = open_restart_read(ifile)

      if (trace_execution) write(*,*) 'opened, iunit = ', iunit

      ! Read in the advance time if present
      if (input_is_model_advance_file) then
         call aread_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit, target_time)
         call get_time(target_time, target_seconds, target_days)
      else
         call aread_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)
      endif
      if (override_init) ens_handle%time(1) = init_time

      if (member == 1) then
         if (trace_execution) call print_time(ens_handle%time(1), 'current model state time')
         if (trace_execution) call print_time(target_time, 'advance-to time')
      endif

      call close_restart(iunit)
      !------------------- Read restart from file ----------------------
      
      !-----------------  Model advance call   --------------------------------
      ! Advance this state to the target time
      ! If the model time is past the obs set time, just need to skip
      if (trace_execution) write(*,*) 'calling advance_state if needed'
   
      if(ens_handle%time(1) < target_time) &
         call advance_state(ens_handle, ens_size=1, target_time=target_time, &
            async=0, adv_ens_command='', tasks_per_model_advance=1)

      !-----------------  Model advance call   --------------------------------

      !------------------- Write restart to file -----------------------
      ! Output the restart file
      if (trace_execution) write(*,*) 'ready to open output restart file ', trim(ofile)

      iunit = open_restart_write(ofile, write_format)

      if (trace_execution) write(*,*) 'opened, iunit = ', iunit

      call awrite_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)

      if (trace_execution) call print_time(ens_handle%time(1), 'time after model advance')
      
      call close_restart(iunit)
      !------------------- Write restart to file -----------------------

   enddo

else

   ! Initialize the ens manager with enough room for all ensemble members.
   ! Either read, write, or both will need this.
   call init_ensemble_manager(ens_handle, num_copies=ens_size, num_vars=model_size)
   call prepare_to_write_to_vars(ens_handle)

   ! make the defaults be a single filename, and overwrite them below if
   ! there are individual files.
   ifile = trim(input_filename)
   ofile = trim(output_filename)
   
   !------------------- Read restart from file ----------------------
   ! If only one restart file on input, read it all up front.  For individuals
   ! read each one in the loop below.
   if (single_restart_file_in) then
      iunit = open_restart_read(ifile)
      ! Read in the advance time if not present
      do member=1, ens_size
         if (input_is_model_advance_file) then
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit, target_time)
            call get_time(target_time, target_seconds, target_days)
         else
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         endif
         if (override_init) ens_handle%time(member) = init_time

         if (member == 1) then
            if (trace_execution) call print_time(ens_handle%time(1), 'current model state time')
            if (trace_execution) call print_time(target_time, 'advance-to time')
         endif
      enddo
      call close_restart(iunit)
    else
      ! loop over each ensemble member, reading in each individually.
      do member=1, ens_size
    
         ! add member number as a suffix: e.g. base.0000
         write(ifile, "(a,a,i4.4)") trim(input_filename), '.', member
   
         iunit = open_restart_read(ifile)
         ! Read in the advance time if present
         if (input_is_model_advance_file) then
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit, target_time)
            call get_time(target_time, target_seconds, target_days)
         else
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         endif
         call close_restart(iunit)

         if (override_init) ens_handle%time(member) = init_time
         if (member == 1) then
            if (trace_execution) call print_time(ens_handle%time(1), 'current model state time')
            if (trace_execution) call print_time(target_time, 'advance-to time')
         endif

      enddo
   endif   
   !------------------- Read restart from file ----------------------
         
   !-----------------  Model advance call   --------------------------------
   ! Advance all states to the target time
   ! If the model time is past the obs set time, just need to skip
   if (trace_execution) write(*,*) 'calling advance_state if needed'

   if(ens_handle%time(1) < target_time) &
      call advance_state(ens_handle, ens_size=ens_size, target_time=target_time, &
         async=0, adv_ens_command='', tasks_per_model_advance=1)

   !-----------------  Model advance call   --------------------------------

   !------------------- Write restart to file -----------------------
   ! If only one restart file on output, write it all in one go.
   ! Otherwise, write each ens in a loop below.
   if (single_restart_file_out) then
      ! Output the restart file
      iunit = open_restart_write(ofile, write_format)
      do member=1, ens_size
         call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
      enddo
      call close_restart(iunit)
      if (member == 1) then 
         if (trace_execution) call print_time(ens_handle%time(1), 'time after model advance')
      endif
    else
      ! loop over each ensemble member, reading in each individually.
      do member=1, ens_size
    
         ! add member number as a suffix: e.g. base.0000
         write(ofile, "(a,a,i4.4)") trim(output_filename), '.', member
   
         ! Output the restart files
         iunit = open_restart_write(ofile, write_format)
         call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         call close_restart(iunit)

         if (member == 1) then
            if (trace_execution) call print_time(ens_handle%time(1), 'time after model advance')
         endif

      enddo
   endif
   !------------------- Write restart to file -----------------------

endif


if (trace_execution) write(*,*) 'end of integrate_model executable'
call finalize_mpi_utilities()

end program integrate_model

