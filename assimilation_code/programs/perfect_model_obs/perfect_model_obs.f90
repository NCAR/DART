! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Program to build an obs_sequence file from simulated observations.

program perfect_model_obs

use        types_mod,     only : r8, i8, metadatalength, MAX_NUM_DOMS
use    utilities_mod,     only : error_handler, &
                                 find_namelist_in_file, check_namelist_read, &
                                 E_ERR, E_MSG, nmlfileunit, timestamp, &
                                 do_nml_file, do_nml_term, open_file, close_file
use time_manager_mod,     only : time_type, get_time, set_time, operator(/=), generate_seed
use obs_sequence_mod,     only : read_obs_seq, obs_type, obs_sequence_type,                 &
                                 get_obs_from_key, set_copy_meta_data, get_obs_def,         &
                                 get_time_range_keys, set_obs_values, set_qc, set_obs,      &
                                 write_obs_seq, get_num_obs, init_obs, assignment(=),       &
                                 static_init_obs_sequence, read_obs_seq_header, &
                                 set_qc_meta_data, destroy_obs, destroy_obs_sequence

use      obs_def_mod,     only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                                 get_obs_def_type_of_obs
use    obs_model_mod,     only : set_obs_model_trace, advance_model
use  assim_model_mod,     only : static_init_assim_model, get_model_size,                    &
                                 get_initial_condition
   
use mpi_utilities_mod,    only : initialize_mpi_utilities, finalize_mpi_utilities

use   random_seq_mod,     only : random_seq_type, init_random_seq, random_gaussian
use ensemble_manager_mod, only : init_ensemble_manager, end_ensemble_manager, ensemble_type,  &
                                 allocate_vars,  all_vars_to_all_copies, all_copies_to_all_vars

use           filter_mod, only : filter_set_initial_time, trim_obs_sequence

use state_vector_io_mod,   only : state_vector_io_init, read_state, write_state

use io_filenames_mod,      only : io_filenames_init, file_info_type, set_file_metadata,  &
                                  set_io_copy_flag, READ_COPY, WRITE_COPY
                                  
use direct_netcdf_mod,     only : finalize_single_file_io

use quality_control_mod,   only : set_input_qc, initialize_qc

use ensemble_manager_mod,  only : set_num_extra_copies

use distributed_state_mod, only : create_state_window, free_state_window

use forward_operator_mod, only : get_expected_obs_distrib_state

use mpi_utilities_mod,    only : my_task_id

use algorithm_info_mod,   only : init_algorithm_info_mod, obs_error_info, end_algorithm_info_mod

implicit none

character(len=*), parameter :: source = 'perfect_model_obs.f90'

! Module storage for message output
character(len=512) :: msgstring
integer            :: trace_level, timestamp_level


!-----------------------------------------------------------------------------
! Namelist with default values
!
logical  :: read_input_state_from_file = .false.
logical  :: write_output_state_to_file = .false.
integer  :: async              = 0
logical  :: trace_execution    = .false.
logical  :: output_timestamps  = .false.
logical  :: silence            = .false.

! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days     = 0
integer  :: init_time_seconds  = 0
! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days     = -1
integer  :: first_obs_seconds  = -1
integer  :: last_obs_days      = -1
integer  :: last_obs_seconds   = -1
integer  :: obs_window_days    = -1
integer  :: obs_window_seconds = -1
logical  :: output_forward_op_errors = .false.
integer  :: tasks_per_model_advance = 1
integer  :: output_interval = 1
integer  :: print_every_nth_obs = 0

logical  :: has_cycling     = .false.
logical  :: single_file_in  = .false.
logical  :: single_file_out = .false.

character(len=256) :: input_state_files(MAX_NUM_DOMS)  = '',               &
                      output_state_files(MAX_NUM_DOMS) = '',                &
                      obs_seq_in_file_name            = 'obs_seq.in',      &
                      obs_seq_out_file_name           = 'obs_seq.out',     &
                      adv_ens_command                 = './advance_model.csh'

namelist /perfect_model_obs_nml/ read_input_state_from_file,&
                                 write_output_state_to_file,                        &
                                 init_time_days, init_time_seconds, async,          &
                                 first_obs_days, first_obs_seconds,                 &
                                 last_obs_days,  last_obs_seconds, output_interval, &
                                 obs_seq_in_file_name, obs_seq_out_file_name,       &
                                 adv_ens_command, tasks_per_model_advance,          & 
                                 obs_window_days, obs_window_seconds, silence,      &
                                 trace_execution, output_timestamps,                &
                                 print_every_nth_obs, output_forward_op_errors,     &
                                 input_state_files, output_state_files,             &
                                 single_file_in, single_file_out

!------------------------------------------------------------------------------

! Doing this allows independent scoping for subroutines in main program file
call perfect_main()

!------------------------------------------------------------------------------

contains

subroutine perfect_main()

type(obs_sequence_type) :: seq
type(obs_type)          :: obs
type(obs_def_type)      :: obs_def
type(random_seq_type)   :: random_seq
type(ensemble_type)     :: ens_handle, fwd_op_ens_handle, qc_ens_handle
type(time_type)         :: curr_ens_time, next_ens_time

integer, allocatable    :: keys(:)
integer                 :: i, j, iunit, time_step_number, obs_seq_file_id
integer                 :: cnum_copies, cnum_qc, cnum_obs, cnum_max
integer                 :: additional_qc, additional_copies, forward_unit
integer                 :: io, num_obs_in_set, nth_obs
integer                 :: key_bounds(2)
integer                 :: seed

integer                 :: istatus(1)

real(r8)                :: true_obs(1), obs_value(1), qc(1)

character(len=metadatalength) :: copy_meta_data(2), qc_meta_data, obs_seq_read_format
character(len=metadatalength) :: state_meta(1)

logical                 :: assimilate_this_ob, evaluate_this_ob, pre_I_format
logical                 :: read_time_from_file

integer                 :: global_obs_num

type(time_type)      :: time1
integer              :: secs, days
character(len=20)    :: task_str ! string to hold the task number
integer              :: ens_size = 1 ! This is to avoid magic number 1s
integer              :: copy_indices(1) = 1

type(file_info_type) :: file_info_input
type(file_info_type) :: file_info_output
type(file_info_type) :: file_info_true

character(len=256), allocatable :: input_filelist(:), output_filelist(:), true_state_filelist(:)
integer :: nfilesin, nfilesout

! Storage for bounded error 
logical  :: bounded_below, bounded_above
real(r8) :: lower_bound,   upper_bound
real(r8) :: error_variance

! Initialize all modules used that require it
call perfect_initialize_modules_used()

! Read the namelist entry
call find_namelist_in_file("input.nml", "perfect_model_obs_nml", iunit)
read(iunit, nml = perfect_model_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "perfect_model_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=perfect_model_obs_nml)
if (do_nml_term()) write(     *     , nml=perfect_model_obs_nml)

! set the level of output
call set_trace(trace_execution, output_timestamps, silence)

!>@todo FIXME: copies and qc should be set using the meta data strings not hard
! coded. 

! Find out how many data copies are in the obs_sequence 
call read_obs_seq_header(obs_seq_in_file_name, cnum_copies, cnum_qc, cnum_obs, cnum_max, &
   obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)

! First two copies of output will be truth and observation;
! Will overwrite first two existing copies in file if there are any
additional_copies = 2 - cnum_copies
if(additional_copies < 0) additional_copies = 0

! Want to have a qc field available in case forward op wont work
if(cnum_qc == 0) then
   additional_qc = 1
else
   additional_qc = 0
endif

! Read in definition part of obs sequence; expand to include observation and truth field
call read_obs_seq(obs_seq_in_file_name, additional_copies, additional_qc, 0, seq)

! Initialize an obs type variable
call init_obs(obs, cnum_copies + additional_copies, cnum_qc + additional_qc)

! Need metadata for added qc field
if(additional_qc == 1) then
   qc_meta_data = 'Quality Control'
   call set_qc_meta_data(seq, 1, qc_meta_data)
endif

! Need space to put in the obs_values in the sequence;
copy_meta_data(1) = 'observations'
copy_meta_data(2) = 'truth'
call set_copy_meta_data(seq, 1, copy_meta_data(1))
call set_copy_meta_data(seq, 2, copy_meta_data(2))

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle, ens_size, get_model_size())

call set_num_extra_copies(ens_handle, 0)

! for now, make the default to allow cycling if we're writing all data
! to a single file, and false if not.
has_cycling = single_file_out

! Initialize file names:
call parse_filenames(input_state_files,  input_filelist,  nfilesin)
call parse_filenames(output_state_files, output_filelist, nfilesout)

!> @todo FIXME  if nfilesout == 0 and write_output_state_to_file is .false.
!> that shouldn't be an error.  if nfilesin == 0 and read_input_state_from_file
!> is false, that also shouldn't be an error.  (unless you're writing the mean
!> and sd, and then maybe we should have a different name for output of input values.)
if (nfilesin == 0 .or. nfilesout == 0 ) then
   msgstring = 'must specify both "input_state_files" and "output_state_files" in the namelist'
   call error_handler(E_ERR,'perfect_main',msgstring,source)
endif

allocate(true_state_filelist(nfilesout))

! mutiple domains ( this is very unlikely to be the case, but in order to
! set_file_metadata we need to have a file list that has the same number
! of files as domains.
if (nfilesout > 1) then
   do i = 1, nfilesout
      write(true_state_filelist(i),'(a,i0.2)') 'true_state_d',i
   enddo
else
   true_state_filelist(1) = 'true_state.nc'
endif

call io_filenames_init(file_info_input,  1, cycling=has_cycling, single_file=single_file_in)
call set_file_metadata(file_info_input,  1, input_filelist, 'perfect_input', 'pmo initial condition')
call set_io_copy_flag( file_info_input,  1, READ_COPY) 

! True State
call io_filenames_init(file_info_true, 1, cycling=has_cycling, single_file=single_file_out)
call set_file_metadata(file_info_true, 1, true_state_filelist, 'true_state', 'true state')
call set_io_copy_flag( file_info_true, 1, 1, WRITE_COPY, num_output_ens=1)

! Perfect Restart
call io_filenames_init(file_info_output, 1, cycling=has_cycling, single_file=single_file_out)
call set_file_metadata(file_info_output, 1, output_filelist, 'perfect_restart', 'pmo restart')
call set_io_copy_flag( file_info_output, 1, 1, WRITE_COPY, num_output_ens=1)

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time(init_time_days, init_time_seconds, time1, read_time_from_file)

if (read_input_state_from_file) then
   call read_state(ens_handle, file_info_input, read_time_from_file, time1)
else ! model spin up
   call allocate_vars(ens_handle)
   if(ens_handle%my_pe == 0) call get_initial_condition(ens_handle%time(1), ens_handle%vars(:, 1))
   call all_vars_to_all_copies(ens_handle)
endif

!>@todo FIXME this block must be supported in the single file loop with time dimension
state_meta(1) = 'true state'

! Remove observations that are too early or too late to be assimilated
call trim_obs_sequence(seq, key_bounds)

! Do verbose forward operator output if requested
write(task_str, '(i10)') ens_handle%my_pe
if(output_forward_op_errors) forward_unit = open_file('forward_op_errors' // TRIM(ADJUSTL(task_str)), 'formatted', 'append')

! Time step number is used to do periodic diagnostic output
curr_ens_time = set_time(0, 0)
next_ens_time = set_time(0, 0)

! Advance model to the closest time to the next available observations
AdvanceTime: do time_step_number = 0, huge(time_step_number)

   ! Advance the model to make the window include the next available observation.
   call advance_model(ens_handle, 1, seq, key_bounds, num_obs_in_set, &
      curr_ens_time, next_ens_time, async, adv_ens_command, tasks_per_model_advance, &
      file_info_output, file_info_input)

   if(key_bounds(1) < 0) exit AdvanceTime

   ! Initialize a repeatable random sequence for perturbations
   seed = generate_seed(next_ens_time)
   call init_random_seq(random_seq,seed)

   ! for multi-core runs, each core needs to store the forward operator and the qc value
   call init_ensemble_manager(fwd_op_ens_handle, ens_size, int(num_obs_in_set,i8), 1, transpose_type_in = 2)
   call init_ensemble_manager(    qc_ens_handle, ens_size, int(num_obs_in_set,i8), 1, transpose_type_in = 2)

   ! Allocate storage for observation keys for this part of sequence
   allocate(keys(num_obs_in_set))

   ! Get all the keys associated with this set of observations
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   ! Output the true state to the netcdf file
   if((output_interval > 0) .and. &
      (time_step_number / output_interval * output_interval == time_step_number)) then

      if(has_cycling) then
         call write_state(ens_handle, file_info_true)
      endif

   endif

   ! Set up access to the state
   call create_state_window(ens_handle)

   ! Compute the forward observation operator for each observation in set
   do j = 1, fwd_op_ens_handle%my_num_vars

      ! Compute the observations from the state
      global_obs_num = fwd_op_ens_handle%my_vars(j)
      call get_expected_obs_distrib_state(seq, keys(global_obs_num:global_obs_num), &
         curr_ens_time, .true., &
         istatus, assimilate_this_ob, evaluate_this_ob, &
         ens_handle, ens_size, copy_indices, true_obs)

      fwd_op_ens_handle%copies(1, j) = true_obs(1)

      qc_ens_handle%copies(1, j) = set_input_qc(istatus(1), assimilate_this_ob, evaluate_this_ob)

      ! FIXME: we could set different qc codes, like 1004, 1005, to
      ! indicate why this obs isn't being processed - separate failed
      ! forward operators from those types not on the assim or eval lists.
      ! the values 4, 5, etc could match the dart QC values + 1000.

      ! if failed forward op logging requested, make sure we're
      ! only writing out obs with real errors and not those that
      ! end up in this code section because their type isn't in the namelist.
      if(output_forward_op_errors) then
         if ((istatus(1) /= 0) .and. (assimilate_this_ob .or. evaluate_this_ob)) &
            write(forward_unit, *) keys(global_obs_num), istatus
      endif

   end do

   ! collect on task 0 and load up the obs_sequence
   call all_copies_to_all_vars(fwd_op_ens_handle)
   call all_copies_to_all_vars(    qc_ens_handle)

   ! Task 0 loads up the obs_sequence
   if(my_task_id() == 0) then

      do i = 1, fwd_op_ens_handle%num_vars

         ! true obs is the forward operator
         ! I'm doing the random noise on task 0 so we remain bitwise the trunk
         true_obs(1) = fwd_op_ens_handle%vars(i, 1)

         ! Get the observational error covariance (diagonal at present)
         ! Generate the synthetic observations by adding in error samples
         call get_obs_from_key(seq, keys(i), obs)
         call get_obs_def(obs, obs_def)

         ! If observation is not being evaluated or assimilated, skip it
         ! Ends up setting a 1000 qc field so observation is not used again.
         if( qc_ens_handle%vars(i, 1) == 0 ) then

            ! Get the information for generating error sample for this observation
            call obs_error_info(obs_def, error_variance, &
                 bounded_below, bounded_above, lower_bound, upper_bound)

            ! Capability to do a bounded normal error
            if(bounded_below .and. bounded_above) then
               ! Bounds on both sides
               obs_value(1) = lower_bound - 1.0_r8
               do while(obs_value(1) < lower_bound .or. obs_value(1) > upper_bound)
                  obs_value(1) = random_gaussian(random_seq, true_obs(1), &
                     sqrt(error_variance))
               end do
            elseif(bounded_below .and. .not. bounded_above) then
               ! Bound on lower side
               obs_value(1) = lower_bound - 1.0_r8
               do while(obs_value(1) < lower_bound)
                  obs_value(1) = random_gaussian(random_seq, true_obs(1), &
                     sqrt(error_variance))
               end do
            elseif(.not. bounded_below .and. bounded_above) then
               ! Bound on upper side
               obs_value(1) = upper_bound + 1.0_r8
               do while(obs_value(1) > upper_bound)
                  obs_value(1) = random_gaussian(random_seq, true_obs(1), &
                     sqrt(error_variance))
               end do
            else
            ! No bounds, regular old normal distribution
               obs_value(1) = random_gaussian(random_seq, true_obs(1), &
                  sqrt(error_variance))
            endif

            ! FIX ME SPINT: if the foward operater passed can we directly set the
            ! qc status?

            ! Set qc to 0 if none existed before
            if(cnum_qc == 0) then
               qc(1) = qc_ens_handle%vars(i, 1)
               call set_qc(obs, qc, 1)
            endif
         else

            obs_value(1) = true_obs(1)
            qc(1) = qc_ens_handle%vars(i, 1)
            call set_qc(obs, qc, 1)

         endif

         ! true obs is fowrard op
         ! obs value is forward op plus noise
         call set_obs_values(obs, obs_value, 1)
         call set_obs_values(obs, true_obs, 2)

         ! Insert the observations into the sequence first copy
         call set_obs(seq, obs, keys(i))

      enddo

   endif

   ! End access to the state
   call free_state_window(ens_handle)

   ! Deallocate the keys storage
   deallocate(keys)

end do AdvanceTime

! if logging errors, close unit
if(output_forward_op_errors) call close_file(forward_unit)

! Write out the sequence
if (ens_handle%my_pe == 0) call write_obs_seq(seq, obs_seq_out_file_name)
call end_ensemble_manager(fwd_op_ens_handle)
call end_ensemble_manager(qc_ens_handle)

! Output a restart file if requested
if(write_output_state_to_file) then
   call write_state(ens_handle, file_info_output)
   if (single_file_out) &
      call finalize_single_file_io(file_info_output)
endif

! Write the last time step to the true state file
if (single_file_out) then
   call write_state(ens_handle, file_info_true)
   call finalize_single_file_io(file_info_true)
endif

!  Release storage for ensemble
call end_ensemble_manager(ens_handle)

! Free up the observation kind and obs sequence
call destroy_obs(obs)
call destroy_obs_sequence(seq)

call end_algorithm_info_mod()

! closes the log file.
call finalize_mpi_utilities()

end subroutine perfect_main

!=====================================================================

subroutine perfect_initialize_modules_used()

! Standard initialization
call initialize_mpi_utilities('perfect_model_obs')

! Initialize modules used that require it

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()

call state_vector_io_init()
call initialize_qc()
call init_algorithm_info_mod()

end subroutine perfect_initialize_modules_used

!-------------------------------------------------------------------------

subroutine set_trace(trace_execution, output_timestamps, silence)

logical, intent(in) :: trace_execution
logical, intent(in) :: output_timestamps
logical, intent(in) :: silence

! Set whether other modules trace execution with messages
! and whether they output timestamps to trace overall performance

! defaults
trace_level     = 0
timestamp_level = 0

! selectively turn stuff back on
if (trace_execution)   trace_level     = 1
if (output_timestamps) timestamp_level = 1

! turn as much off as possible
if (silence) then
   trace_level     = -1
   timestamp_level = -1
endif

call set_obs_model_trace(trace_level, timestamp_level)

end subroutine set_trace

!-------------------------------------------------------------------------

subroutine parse_filenames(file_array, files_out, nfiles)
character(len=*), intent(in)  :: file_array(:)
integer,          intent(out) :: nfiles
character(len=*), allocatable, intent(out) :: files_out(:)
integer :: i

! count the number of valid files
nfiles = 0
do i = 1, size(file_array(:),1)
   if ( file_array(i) == '' ) exit
   nfiles = nfiles + 1 
enddo

! allocate and set output file list
allocate(files_out(nfiles))
do i = 1, nfiles
   files_out(i) = file_array(i)
enddo

end subroutine parse_filenames


end program perfect_model_obs

