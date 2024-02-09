! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Program to integrate assimilation model forward for asynchronous filter
!> execution.

program integrate_model

!> @todo FIXME:  i think this program needs to have 2 restart lists,
!> one for input, one for output, and possibly an array of ensemble
!> numbers (not sure about that), which it then uses the read_state()
!> and write_state() routines to get and put the data after the
!> model advances.

use time_manager_mod,    only : time_type, operator(<), print_time
use utilities_mod,       only : error_handler, E_MSG, nmlfileunit, &
                                do_nml_file, do_nml_term,          &
                                find_namelist_in_file, check_namelist_read
use assim_model_mod,     only : static_init_assim_model, get_model_size

use obs_model_mod,        only : advance_state
use ensemble_manager_mod, only : init_ensemble_manager, ensemble_type

use mpi_utilities_mod,    only : initialize_mpi_utilities, finalize_mpi_utilities, &
                                 task_count, iam_task0

use state_vector_io_mod,  only : read_state, write_state

use io_filenames_mod,     only : file_info_type, io_filenames_init, set_io_copy_flag,  &
                                 set_file_metadata, READ_COPY, WRITE_COPY, &
                                 get_restart_filename, get_stage_metadata

use types_mod,            only : i8

implicit none

character(len=*), parameter :: source = 'integrate_model.f90'

type(ensemble_type) :: ens_handle
type(time_type)     :: target_time
integer             :: iunit, rc
integer(i8)         :: model_size

type(file_info_type) :: input_file_info, output_file_info

! dummy args for the advance_state call.  presumably this 
! executable was invoked by a script that was originally 
! called by advance_state() in filter, so no need to get 
! recursive here.  but if you want to use this outside of
! the dart advance, you could comment these in and add
! them to the namelist below.
!character (len=129) :: adv_ens_command = ''
!integer             :: async = 0
!integer             :: tasks_per_model_advance = 1

!----------------------------------------------------------------
!>@todo either use the namelist or take it out. 

! to enable the use of the namelist, change this to .true. and recompile.
logical :: use_namelist = .false.

! for use outside dart, or to call after dart finishes and
! writes a restart file (without a target time), set this
! to .false. and supply a target time with target_days, secs
logical :: is_model_advance_file = .true.

! for debugging, status
character(len=256) :: ic_file_name = "temp_ic.nc"
character(len=256) :: ud_file_name = "temp_ud.nc"
logical             :: trace_execution = .false.
character(len=512)  :: errstring

! to overwrite the target time, set these to something >= 0
integer :: target_days = -1, target_seconds = -1

! only read in if use_namelist is .true. -- ignored otherwise.
namelist /integrate_model_nml/ &
   ic_file_name, ud_file_name, &
   trace_execution, is_model_advance_file, target_days, target_seconds

!----------------------------------------------------------------

! This program should only be run with a single process
call initialize_mpi_utilities('integrate_model')

! This version is ok to build with MPI and run with more than a single
! task, but only task 0 reads the input, advances the model, and writes 
! the output.  All the other tasks just hang out and exit when task 0 is done.
! FIXME: that could be changed to do multiple ens members in parallel.
! the code does NOT do this now.
if(task_count() > 1) &
   call error_handler(E_MSG,'integrate_model','Only one process doing the work', source)


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

! Initialize the model class data
call static_init_assim_model()
model_size = get_model_size()

! Initialize an ensemble manager type with a single copy
call init_ensemble_manager(ens_handle, num_copies=1, num_vars=model_size, transpose_type_in = 2)

   !------------------- Read restart from file ----------------------
   
   call io_filenames_init(input_file_info, 1, cycling=.true., single_file=.true.)
   call set_file_metadata(input_file_info, 1, (/ic_file_name/), 'temp_ic', 'initial condition')
   call set_io_copy_flag( input_file_info, 1, READ_COPY) 
   
   ic_file_name = get_restart_filename(get_stage_metadata(input_file_info),1,1)
   if (trace_execution) write(*,*) 'ready to open input restart file ', trim(ic_file_name)
   call read_state(ens_handle, input_file_info, .true., target_time)

   if (trace_execution) write(*,*) 'time of data, advance-to are:'
   if (trace_execution) call print_time(ens_handle%time(1))
   if (trace_execution) call print_time(target_time)

   !-----------------  Restart read in --------------------------------

   ! Advance this state to the target time
   ! If the model time is past the obs set time, just need to skip
   if (trace_execution) write(*,*) 'calling advance_state if needed'

   call io_filenames_init(output_file_info, 1, cycling=.true., single_file=.true.)
   call set_file_metadata(output_file_info, 1, (/ud_file_name/), 'temp_ud', 'advanced member')
   call set_io_copy_flag( output_file_info, 1, WRITE_COPY) 

   if(ens_handle%time(1) < target_time) &
      call advance_state(ens_handle, 1, target_time, 0, '', 1, input_file_info, output_file_info)

   ! Output the restart file if requested; Force to binary for bitwise reproducing
   ! use in filter and perfect_model obs with shell advance options
   ud_file_name = get_restart_filename(get_stage_metadata(output_file_info),1,1)
   if (trace_execution) write(*,*) 'ready to open output restart file ', trim(ud_file_name)

   call write_state(ens_handle, output_file_info)

if (trace_execution) write(*,*) 'end of integrate_model executable'
call finalize_mpi_utilities()

end program integrate_model

