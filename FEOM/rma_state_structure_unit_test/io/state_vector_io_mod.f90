! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!----------------------------------------------------------------------
module state_vector_io_mod
!> \defgroup state_vector_io_mod state_vector_io_mod
!> @{ \brief Routines for reading and writing the model state.
!>
!> Read_state() and write_state() are the major routines in this module. filter_write_restart_direct()
!> is used in diagnostic files writes (for certain diagnostic options).
!> They read/write state from an ensmeble_handle%copies array to files described in a file_info_type.
!> The inflation handles (prior and posterior) are optional arguements to read/write_state since
!> inflation is used in filter but not in perfect_model_obs.
!>
!> Two types of state files are supported:
!>     * netcdf
!>     * dart format (this may be ascii or binary set by namelist option: write_binary_restart_files )
!> 
!> netcdf files use \ref io_filenames_mod. See model_mod for construct_file_name_in() and
!> io_filenames_mod for construct_file_name_out(). See \ref io_filenames_mod for filenames 
!> for extra copies (mean, etc.)
!> Each domain is written to a separate file.
!>
!> dart format files use the naming convention from Lanai, e.g. filter_ics.0001 where
!>  filter_ics is the restart_in_base passed to io_filenames_init().
!> 
!> Perturbation (creating an ensemble from a single instance) is done separatly from the read.
!> perturb_from_single_copy is passed as an arguement to read_state. When reading netcdf files
!> only the first ensemble member is read and transposed.  Note for dart format restarts the 
!> whole ensemble is transpose so this is over communicating. create_ensemble_from_single_file()
!> is then called in filter to generate the ensemble.

use adaptive_inflate_mod, only : adaptive_inflate_type, mean_from_restart, sd_from_restart, &
                                 do_single_ss_inflate, do_varying_ss_inflate, get_inflate, &
                                 get_sd, get_inflation_in_filename, get_inflation_out_filename, output_inf_restart, &
                                 get_inflate_mean, get_inflate_sd, do_ss_inflate

use direct_netcdf_mod,    only : read_transpose, transpose_write

use types_mod,            only : r8, MISSING_R8, i4, i8

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id, &
                                 broadcast_send, broadcast_recv

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, &
                                 map_task_to_pe, &
                                 get_copy_owner_index, put_copy, get_allow_transpose, &
                                 all_copies_to_all_vars, all_vars_to_all_copies, &
                                 prepare_to_write_to_vars, get_var_owner_index
                                 

use utilities_mod,        only : error_handler, nc_check, check_namelist_read, &
                                 find_namelist_in_file, nmlfileunit,           &
                                 do_nml_file, do_nml_term,          &
                                 E_MSG, E_ERR, E_DBG, get_unit, ascii_file_format, &
                                 close_file, dump_unit_attributes, &
                                 register_module, set_output

use assim_model_mod,      only : assim_model_type

use time_manager_mod,     only : time_type, read_time, write_time, &
                                 get_time

use io_filenames_mod,     only : get_input_file, file_info_type, get_read_from_netcdf, get_write_to_netcdf, &
                                 get_output_restart, get_output_mean, get_restart_out_base, &
                                 get_restart_in_base, get_read_from_single_file, get_write_to_single_file, &
                                 assert_file_info_initiailzed, restart_names_type, &
                                 assert_restart_names_initiailzed

!> @todo  This should go through assim_model_mod
use model_mod,            only : read_model_time

use copies_on_off_mod,    only : setup_read_write, turn_read_copy_on,      &
                                 turn_write_copy_on, turn_read_copies_off, &
                                 turn_write_copy_off, end_read_write, &
                                 ENS_MEAN_COPY, &
                                 PRIOR_INF_COPY, PRIOR_INF_SD_COPY, &
                                 POST_INF_COPY, POST_INF_SD_COPY, &
                                 SPARE_PRIOR_MEAN, SPARE_PRIOR_SPREAD, &
                                 SPARE_PRIOR_INF_MEAN, SPARE_PRIOR_INF_SPREAD, &
                                 SPARE_POST_INF_MEAN, SPARE_POST_INF_SPREAD, &
                                 query_copy_present
                                 

use state_structure_mod,  only : get_num_domains

use netcdf


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

private

! Initialize, read, write routines for filter and perfect_model_obs
public :: state_vector_io_init, &
          read_state, &
          write_state

! For diagnostic file writes.
public :: turn_write_copy_on, &
          turn_write_copy_off, &
          setup_read_write, &
          end_read_write, &
          filter_write_restart_direct

! obs_model_mod - no netcdf calls in obs_model_mod yet
! model converters - don't use netcdf.
public :: aread_state_restart, &
          awrite_state_restart, &
          open_restart_read, &
          open_restart_write, &
          close_restart

! smoother - this is left over code.
public :: read_ensemble_restart, &
          write_ensemble_restart


! Module storage for writing error messages
character(len=512) :: msgstring

! Logical flag for initialization of module
logical :: module_initialized = .false.

! Global storage for default restart formats
character(len=16) :: read_format = "unformatted", write_format = "unformatted"

! namelist variables with default values
! Aim: to have the regular transpose as the default
integer :: limit_mem = HUGE(1_i4)!< This is the number of elements (not bytes) so you don't have times the number by 4 or 8
logical :: single_precision_output = .false. ! Allows you to write r4 netcdf files even if filter is double precision

! write_binary_restart_files  == .true.  -> use unformatted file format.
!                                     Full precision, faster, smaller,
!                                     but not as portable.
logical  :: write_binary_restart_files = .false.

namelist /  state_vector_io_nml / limit_mem, &
   single_precision_output, &
   write_binary_restart_files

contains

!-----------------------------------------------------------------------
!> Initialize module and read the namelist

subroutine state_vector_io_init()

integer :: iunit, io

if ( .not. module_initialized ) then
   ! Initialize the module with utilities 
   call register_module(source, revision, revdate)
   module_initialized = .true.

   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "state_vector_io_nml", iunit)
   read(iunit, nml = state_vector_io_nml, iostat = io)
   call check_namelist_read(iunit, io, "state_vector_io_nml")

   ! Set the write format for restart files
   if(write_binary_restart_files) then
      write_format = "unformatted"
   else
      write_format = "formatted"
   endif

endif

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=state_vector_io_nml)
if (do_nml_term()) write(     *     , nml=state_vector_io_nml)

end subroutine state_vector_io_init


!-----------------------------------------------------------------------
!> Read in the state vectors into an ensemble handle
!> The copies will end up in the state_ens_handle%copies array.
!> If read_time_from_file = .true. then time is overwritten by the time read from 
!> the restart. If read_time_from_file = .false. then the time associated with each ensemble 
!> member is set to time.
!> Inflation handles are optional (used by filter, not by perfect_model_obs)
!>   - Note the user must give both inflation handles or neither.
!> perturb_from_single_copy is optional (used by filter, not by perfect_model_obs)
!>   - only 1 restart file is read in. For netcdf read only one copy is transposed. Still
!>     doing a full all vars to all copies for dart format read though.

subroutine read_state(state_ens_handle, file_info, read_time_from_file, time, prior_inflate_handle, post_inflate_handle, perturb_from_single_copy)

type(ensemble_type),         intent(inout) :: state_ens_handle
type(file_info_type),        intent(in)    :: file_info
logical,                     intent(in)    :: read_time_from_file ! state time
type(time_type),             intent(inout) :: time
type(adaptive_inflate_type), optional, intent(in) :: prior_inflate_handle
type(adaptive_inflate_type), optional, intent(in) :: post_inflate_handle
logical,                     optional, intent(in) :: perturb_from_single_copy

integer :: ens_size
logical :: inflation_handles
logical :: local_pert

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! check whether file_info handle is initialized
call assert_file_info_initiailzed(file_info, 'read_state')

! check that we either have both inflation handles or neither:
if ( present(prior_inflate_handle) .neqv. present(post_inflate_handle) ) then
   call error_handler(E_ERR, 'read_state', 'must have both inflation handles or neither', &
          source,revision,revdate)
endif

if (present(prior_inflate_handle) .and. present(post_inflate_handle)) then
   inflation_handles = .true.
else
   inflation_handles = .false.
endif

if (present(perturb_from_single_copy)) then
   local_pert = perturb_from_single_copy
else
   local_pert = .false.
endif

ens_size = state_ens_handle%num_copies - state_ens_handle%num_extras

! set up arrays for which copies to read/write
call setup_read_write(state_ens_handle%num_copies)
call turn_read_copies_off(1, state_ens_handle%num_copies)

if (local_pert) then
   call error_handler(E_MSG,'read_state:', &
      'Reading in a single ensemble and perturbing data for the other ensemble members')
else
   call error_handler(E_MSG,'read_state:', &
      'Reading in initial condition/restart data for all ensemble members from file(s)')
endif

if (get_read_from_netcdf(file_info)) then ! netcdf

   ! Read in restart files
   if (local_pert) then
      call turn_read_copy_on(1)
   else
      call turn_read_copy_on(1, ens_size) ! need to read all restart copies
   endif

   ! inflation read from netcdf files only for state_space_inflation
   ! Ideally would only do a read transpose for VARYING state_space_inflation
   if (inflation_handles) then
      if (do_ss_inflate(prior_inflate_handle)) then
         if (mean_from_restart(prior_inflate_handle)) call turn_read_copy_on(PRIOR_INF_COPY)
         if (sd_from_restart(prior_inflate_handle))   call turn_read_copy_on(PRIOR_INF_SD_COPY)
      endif

      if (do_ss_inflate(post_inflate_handle)) then
         if (mean_from_restart(post_inflate_handle)) call turn_read_copy_on(POST_INF_COPY)
         if (sd_from_restart(post_inflate_handle))   call turn_read_copy_on(POST_INF_SD_COPY)
      endif
   endif

   call filter_read_restart_direct(state_ens_handle, file_info, read_time_from_file, time)

   ! Filling copies array for single state space inflation values if required.
   if (inflation_handles) then

      ! If inflation is single state space read from a file, the copies array is filled here.
      call fill_single_ss_inflate_from_read(state_ens_handle, prior_inflate_handle, post_inflate_handle)

      ! If inflation is from a namelist value it is set here.
      call fill_ss_from_nameslist_value(state_ens_handle, prior_inflate_handle, post_inflate_handle)

      ! To match Lanai filter_state_space_diagnostics, if not doing inflation set inf_mean = 1, inf_sd = 0
      if (.not. do_ss_inflate(prior_inflate_handle)) then
         state_ens_handle%copies(PRIOR_INF_COPY, :)    = 1.0_r8
         state_ens_handle%copies(PRIOR_INF_SD_COPY, :) = 0.0_r8
      endif
      if (.not. do_ss_inflate(post_inflate_handle)) then
         state_ens_handle%copies(POST_INF_COPY, :)    = 1.0_r8
         state_ens_handle%copies(POST_INF_SD_COPY, :) = 0.0_r8
      endif

   endif

else ! expecting DART restart files

   ! Allocate storage space in ensemble manager.
   ! First check that %vars is not already allocated.
   if(.not. get_allow_transpose(state_ens_handle) ) allocate(state_ens_handle%vars(state_ens_handle%num_vars, state_ens_handle%my_num_copies))

   ! Read in restart files
   call filter_read_restart(state_ens_handle, file_info, read_time_from_file, local_pert, time)

   ! inflation read (for state_space_inflation only)
   if (inflation_handles) then
      if (do_ss_inflate(prior_inflate_handle)) then
         ! This routine checks mean_from_restart, sd_from_restart
         call read_state_space_inflation(state_ens_handle, prior_inflate_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
      endif
      if (do_ss_inflate(post_inflate_handle)) then
         ! This routine checks mean_from_restart, sd_from_restart
         call read_state_space_inflation(state_ens_handle, post_inflate_handle, POST_INF_COPY, POST_INF_SD_COPY)
      endif
   endif

   call all_vars_to_all_copies(state_ens_handle)
   if (.not. get_allow_transpose(state_ens_handle))deallocate(state_ens_handle%vars)

endif

call end_read_write()

end subroutine read_state


!-----------------------------------------------------------------------
!> Write state vectors from an ensemble_handle
!> The inflation handles are options (used by filter, not by perfect_model_obs)
!> Note both inflation handles must be present (or neither)
!> There is logic in write_state about spare_copies. These are extra state vector 
!> copies that may be in the ensemble handle when doing large model, single timestep runs.

subroutine write_state(state_ens_handle, file_info, prior_inflate_handle, post_inflate_handle)

type(ensemble_type),         intent(inout) :: state_ens_handle
type(file_info_type),        intent(in)    :: file_info
type(adaptive_inflate_type), optional, intent(in)    :: prior_inflate_handle
type(adaptive_inflate_type), optional, intent(in)    :: post_inflate_handle

integer :: ens_size

logical :: inflation_handles

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! check whether file_info handle is initialized
call assert_file_info_initiailzed(file_info, 'write_state')

! check that we either have both inflation handles or neither:
if ( present(prior_inflate_handle) .neqv. present(post_inflate_handle) ) then
   call error_handler(E_ERR, 'read_state', 'must have both inflation handles or neither', &
        source,revision,revdate)
endif

if (present(prior_inflate_handle) .and. present(post_inflate_handle)) then
   inflation_handles = .true.
else
   inflation_handles = .false.
endif

ens_size = state_ens_handle%num_copies - state_ens_handle%num_extras

! set up arrays for which copies to read/write
! clean slate of files
call setup_read_write(state_ens_handle%num_copies)
call turn_write_copy_off(1, state_ens_handle%num_copies)

if (get_write_to_netcdf(file_info)) then ! netcdf

   ! ------ turn on copies for: restarts, mean and spread ------
   if (get_output_restart(file_info)) call turn_write_copy_on(1, ens_size) ! restarts
   if (get_output_mean(file_info))    call turn_write_copy_on(ENS_MEAN_COPY)

   ! These two copies are the mean and spread that would have gone in the Prior_Diag.nc file
   if(.true.) then ! spare copies, single time step
      call turn_write_copy_on(SPARE_PRIOR_MEAN)
      call turn_write_copy_on(SPARE_PRIOR_SPREAD)
   endif
   !--------------------------------------------


   ! ------ turn on copies for: inflation ------
   ! inflation written to netcdf files only for state_space_inflation
   if (inflation_handles) then

      ! PRIOR
      if (do_ss_inflate(prior_inflate_handle)) then

         ! These two copies are input for the next run of filter
         if (output_inf_restart(prior_inflate_handle)) then
            call turn_write_copy_on(PRIOR_INF_COPY)
            call turn_write_copy_on(PRIOR_INF_SD_COPY)
         endif

         ! These two copies are the inflation values that would have gone in the Prior_Diag.nc file
         ! Assume if they are there, they need to be written out
         if (query_copy_present(SPARE_PRIOR_INF_MEAN))    call turn_write_copy_on(SPARE_PRIOR_INF_MEAN)
         if (query_copy_present(SPARE_PRIOR_INF_SPREAD))  call turn_write_copy_on(SPARE_PRIOR_INF_SPREAD)

      endif

      ! POSTERIOR
      if (do_ss_inflate(post_inflate_handle)) then

         ! These two copies are input for the next run of filter
         if (output_inf_restart(post_inflate_handle)) then
            call turn_write_copy_on(POST_INF_COPY)
            call turn_write_copy_on(POST_INF_SD_COPY)
         endif

         ! These two copies are the inflation values that would have gone in the Posterior_Diag.nc file
         ! Assume if they are there, they need to be written out
         if (query_copy_present(SPARE_POST_INF_MEAN))   call turn_write_copy_on(SPARE_POST_INF_MEAN)
         if (query_copy_present(SPARE_POST_INF_SPREAD)) call turn_write_copy_on(SPARE_POST_INF_SPREAD)
      endif
      
   endif
   !----------------------------------------------

   call filter_write_restart_direct(state_ens_handle, file_info%restart_files_out)
   
else ! dart format restarts

   ! allocating storage space in ensemble manager
   if (.not. get_allow_transpose(state_ens_handle)) allocate(state_ens_handle%vars(state_ens_handle%num_vars, state_ens_handle%my_num_copies))
   call all_copies_to_all_vars(state_ens_handle)

   if(get_output_restart(file_info)) then
      call write_ensemble_restart(state_ens_handle, file_info, get_restart_out_base(file_info), 1, ens_size)
   endif

   if(get_output_mean(file_info)) then
      call write_ensemble_restart(state_ens_handle, file_info, trim(get_restart_out_base(file_info))//'.mean', &
                               ENS_MEAN_COPY, ENS_MEAN_COPY, .true.)
   endif

   if (inflation_handles) then
      if (do_ss_inflate(prior_inflate_handle)) then
         ! This routine checks output_inf_restart(inflate_handle)
         call write_state_space_inflation(prior_inflate_handle, file_info, state_ens_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
      endif

      if (do_ss_inflate(post_inflate_handle)) then
         ! This routine checks output_inf_restart(inflate_handle)
         call write_state_space_inflation(post_inflate_handle, file_info, state_ens_handle, POST_INF_COPY, POST_INF_SD_COPY)
      endif
   endif

   ! Destroying storage in ensemble manager.
   if (.not. get_allow_transpose(state_ens_handle)) deallocate(state_ens_handle%vars)

endif

call end_read_write()

end subroutine write_state


!-------------------------------------------------
! DART read write (binary/ascii vector)
!-------------------------------------------------

!-----------------------------------------------------------------------
!> Read in restart files from dart format files (time + vector)
!> Restarts may be:
!>   one file per restart
!>   all restarts in one file
!> The file_info type contains the read options (see get_read_from_single_file(file_info))
subroutine filter_read_restart(state_ens_handle, file_info, read_time_from_file, perturb_from_single_copy,  time)

type(ensemble_type),  intent(inout) :: state_ens_handle
type(file_info_type), intent(in)    :: file_info
logical,              intent(in)    :: read_time_from_file
logical,              intent(in)    :: perturb_from_single_copy
type(time_type),      intent(inout) :: time

integer :: days, secs
integer :: start_copy, end_copy ! 1, ens_size or 1,1 if perturbing from a single instance.

! check whether file_info handle is initialized
call assert_file_info_initiailzed(file_info, 'filter_read_restart')

start_copy = 1
if (perturb_from_single_copy) then
   end_copy = 1
else
   end_copy = state_ens_handle%num_copies - state_ens_handle%num_extras ! ens_size
endif

! Only read in initial conditions for actual ensemble members, assuming these
! are copies 1 to ens_size
if(read_time_from_file) then

   call read_ensemble_restart(state_ens_handle, start_copy, end_copy, get_restart_in_base(file_info), force_single_file=get_read_from_single_file(file_info))
      !> @todo time
   if (state_ens_handle%my_num_copies > 0) time = state_ens_handle%time(1)

else

   call read_ensemble_restart(state_ens_handle, start_copy, end_copy, get_restart_in_base(file_info), time, force_single_file=get_read_from_single_file(file_info))
      call get_time(time, secs, days)
      write(msgstring, '(A)') 'By namelist control, ignoring time found in restart file.'
      call error_handler(E_MSG,'filter_read_restart:',msgstring,source,revision,revdate)
      write(msgstring, '(A,I6,1X,I5)') 'Setting initial days, seconds to ',days,secs
      call error_handler(E_MSG,'filter_read_restart:',msgstring,source,revision,revdate)

endif

! Temporary print of initial model time
if(state_ens_handle%my_pe == 0) then
   ! FIXME for the future: if pe 0 is not task 0, pe 0 can not print debug messages
   call get_time(time, secs, days)
   write(msgstring, *) 'initial model time of 1st ensemble member (days,seconds) ',days,secs
   call error_handler(E_DBG,'filter_read_restart',msgstring,source,revision,revdate)
endif

end subroutine filter_read_restart


!-------------------------------------------------
! Netcdf read write
! Uses direct_netcdf_mod.f90
!-------------------------------------------------

!----------------------------------------------------------------------
!> Read the restart information directly from the model output netcdf file

subroutine filter_read_restart_direct(state_ens_handle, file_info, use_time_from_file, time)

type(ensemble_type),  intent(inout) :: state_ens_handle
type(file_info_type), intent(in)    :: file_info
logical,              intent(in)    :: use_time_from_file
type(time_type),      intent(inout) :: time

integer :: dart_index !< where to start in state_ens_handle%copies
integer :: domain !< loop index

! check whether file_info handle is initialized
call assert_file_info_initiailzed(file_info, 'filter_read_restart_direct')

! read time from input file if time not set in namelist
!> @todo Check time constistency across files? This is assuming they are consistent.
if(use_time_from_file) then
   time = read_model_time(get_input_file(file_info%restart_files_in, 1,1)) ! Any of the restarts?
endif

state_ens_handle%time = time

! read in the data and transpose
dart_index = 1 ! where to start in state_ens_handle%copies - this is modified by read_transpose
do domain = 1, get_num_domains()
   call read_transpose(state_ens_handle, file_info%restart_files_in, domain, dart_index, limit_mem)
enddo

! Need Temporary print of initial model time?

end subroutine filter_read_restart_direct


!-----------------------------------------------------------------------
!> write the restart information directly into the model netcdf file.

subroutine filter_write_restart_direct(state_ens_handle, file_name_handle)

type(ensemble_type),   intent(inout) :: state_ens_handle
type(restart_names_type), intent(in) :: file_name_handle

integer :: dart_index !< where to start in state_ens_handle%copies
integer :: domain !< loop index

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! check whether file_info handle is initialized
call assert_restart_names_initiailzed(file_name_handle, 'filter_write_restart_direct')

! transpose and write out the data
dart_index = 1

! Different filenames for prior vs. posterior vs. diagnostic files
do domain = 1, get_num_domains()
   call transpose_write(state_ens_handle, file_name_handle, domain, dart_index, limit_mem, single_precision_output)
enddo

end subroutine filter_write_restart_direct


!-----------------------------------------------------------------------
!> Single state space inflation from file value is only index 1 of 
!> inflation vars array (mean and sd).
!> This routine find the owner of the 1st element in the vars array
!> The owner broadcasts the values of single state space inflation 
!> to all other tasks who then update their copies array.
!> Note filling both mean and sd values if at least one of mean
!> or sd is read from file.  If one is set from a namelist the copies 
!> array is overwritten in fill_ss_from_namelist value

subroutine fill_single_ss_inflate_from_read(ens_handle, prior_inflate_handle, post_inflate_handle)

type(ensemble_type),         intent(inout) :: ens_handle
type(adaptive_inflate_type), intent(in)    :: prior_inflate_handle
type(adaptive_inflate_type), intent(in)    :: post_inflate_handle

integer :: owner, owners_index
integer(i8) :: first_element
real(r8), allocatable :: inf_array(:) ! 2 or 4 values
integer :: inf_count
logical :: return_me ! flag to return if not read any inflation values from files

! Return if not doing single state space inflation
if (.not. do_single_ss_inflate(prior_inflate_handle) .and. .not. do_single_ss_inflate(post_inflate_handle)) return

return_me = .true.
! Return if not reading any state space inflation values from files
if ( do_single_ss_inflate(prior_inflate_handle)) then
   if (mean_from_restart(prior_inflate_handle)) return_me = .false.
   if (sd_from_restart(prior_inflate_handle))   return_me = .false.
endif
if ( do_single_ss_inflate(post_inflate_handle)) then
   if (mean_from_restart(post_inflate_handle)) return_me = .false.
   if (sd_from_restart(post_inflate_handle))   return_me = .false.
endif
if (return_me) return

! At least some single state space inflation values are being read from files.
inf_count = 0
if (do_single_ss_inflate(prior_inflate_handle)) inf_count = 2
if (do_single_ss_inflate(post_inflate_handle))  inf_count = inf_count + 2

allocate(inf_array(inf_count)) ! for sending and recveiving inflation values

! Find out who owns the first element of vars array
first_element = 1
call get_var_owner_index(first_element, owner, owners_index)

if (ens_handle%my_pe == owner) then
   if (do_single_ss_inflate(prior_inflate_handle) .and. do_single_ss_inflate(post_inflate_handle)) then

      inf_array(1) = ens_handle%copies(PRIOR_INF_COPY, owners_index)
      inf_array(2) = ens_handle%copies(PRIOR_INF_SD_COPY, owners_index)
      inf_array(3) = ens_handle%copies(POST_INF_COPY, owners_index)
      inf_array(4) = ens_handle%copies(POST_INF_SD_COPY, owners_index)

   elseif (do_single_ss_inflate(post_inflate_handle) .and. .not. do_single_ss_inflate(post_inflate_handle)) then

      inf_array(1) = ens_handle%copies(PRIOR_INF_COPY, owners_index)
      inf_array(2) = ens_handle%copies(PRIOR_INF_SD_COPY, owners_index)

   elseif(.not. do_single_ss_inflate(post_inflate_handle) .and. do_single_ss_inflate(post_inflate_handle)) then

      inf_array(1) = ens_handle%copies(POST_INF_COPY, owners_index)
      inf_array(2) = ens_handle%copies(POST_INF_SD_COPY, owners_index)

   endif

   call broadcast_send(map_pe_to_task(ens_handle, owner), inf_array)

else

   call broadcast_recv(map_pe_to_task(ens_handle, owner), inf_array)

   if (do_single_ss_inflate(prior_inflate_handle) .and. do_single_ss_inflate(post_inflate_handle)) then

      ens_handle%copies(PRIOR_INF_COPY, owners_index)     = inf_array(1)
      ens_handle%copies(PRIOR_INF_SD_COPY, owners_index)  = inf_array(2)
      ens_handle%copies(POST_INF_COPY, owners_index)      = inf_array(3)
      ens_handle%copies(POST_INF_SD_COPY, owners_index)   = inf_array(4)

   elseif (do_single_ss_inflate(prior_inflate_handle) .and. .not. do_single_ss_inflate(post_inflate_handle)) then

      ens_handle%copies(PRIOR_INF_COPY, owners_index)    = inf_array(1)
      ens_handle%copies(PRIOR_INF_SD_COPY, owners_index) = inf_array(2)

   elseif(.not. do_single_ss_inflate(prior_inflate_handle) .and. do_single_ss_inflate(post_inflate_handle)) then

      ens_handle%copies(POST_INF_COPY, owners_index)    = inf_array(1)
      ens_handle%copies(POST_INF_SD_COPY, owners_index) = inf_array(2)

   endif

endif

end subroutine fill_single_ss_inflate_from_read


!-----------------------------------------------------------------------
!> Check whether inflation values come from namelist and
!> fill copies array with namelist values for inflation if they do.

subroutine fill_ss_from_nameslist_value(ens_handle, prior_inflate_handle, post_inflate_handle)

type(ensemble_type),         intent(inout) :: ens_handle
type(adaptive_inflate_type), intent(in)    :: prior_inflate_handle
type(adaptive_inflate_type), intent(in)    :: post_inflate_handle

if (do_ss_inflate(prior_inflate_handle)) then
   if (.not. mean_from_restart(prior_inflate_handle)) ens_handle%copies(PRIOR_INF_COPY, :)    = get_inflate_mean(prior_inflate_handle)
   if (.not. sd_from_restart(prior_inflate_handle))   ens_handle%copies(PRIOR_INF_SD_COPY, :) = get_inflate_sd(prior_inflate_handle)
endif

if (do_ss_inflate(post_inflate_handle)) then
   if (.not. mean_from_restart(post_inflate_handle)) ens_handle%copies(POST_INF_COPY, :)    = get_inflate_mean(post_inflate_handle)
   if (.not. sd_from_restart(post_inflate_handle))   ens_handle%copies(POST_INF_SD_COPY, :) = get_inflate_sd(post_inflate_handle)
endif

end subroutine fill_ss_from_nameslist_value


!------------------------------------------------------------------
! Reading/writing DART-format inflation files
!------------------------------------------------------------------

!-----------------------------------------------------------------------
!> Read the dart format state space inflation files

subroutine read_state_space_inflation(ens_handle, inflate_handle, ss_inflate_mean_index, ss_inflate_sd_index)

type(ensemble_type),         intent(inout) :: ens_handle
type(adaptive_inflate_type), intent(in) :: inflate_handle
integer,                     intent(in) :: ss_inflate_mean_index, ss_inflate_sd_index

integer :: owner, owners_index

! Fixed or varying state space inflation types
if(do_ss_inflate(inflate_handle)) then
   ! Initialize state space inflation, copies in ensemble are given
   ! by the inflate and inflate_sd indices. These should be contiguous.

   ! Verify that indices are contiguous
   if(ss_inflate_sd_index /= ss_inflate_mean_index + 1) then
      write(msgstring, *) 'ss_inflate_mean_index = ', ss_inflate_mean_index, &
         ' and ss_inflate_sd_index = ', ss_inflate_sd_index, ' must be continguous'
      call error_handler(E_ERR, 'adaptive_inflate_init', &
         msgstring, source, revision, revdate)
   endif

   ! Read in initial values from file OR get from namelist arguments

   ! If either mean, sd, or both are to be read from the restart file, read them in.
   ! There is no option to read only one; to get either you have to read both.
   ! If one is to be set from the namelist, it gets overwritten in the block below
   ! this one.
   if(mean_from_restart(inflate_handle) .or. sd_from_restart(inflate_handle)) then
      ! the .true. below is 'start_from_restart', which tells the read routine to
      ! read in the full number of ensemble members requested (as opposed to reading
      ! in one and perturbing it).
      call read_ensemble_restart(ens_handle, ss_inflate_mean_index, ss_inflate_sd_index, get_inflation_in_filename(inflate_handle), force_single_file = .true.)

   endif
   ! Now, if one or both values come from the namelist (i.e. is a single static
   ! value), write or overwrite the arrays here.
   if (.not. mean_from_restart(inflate_handle) .or. .not. sd_from_restart(inflate_handle)) then
      ! original code required an expensive transpose which is not necessary.
      ! if setting initial values from the namelist, find out which task has the
      ! inflation and inf sd values and set them only on that task.  this saves us
      ! a transpose.
      if (.not. mean_from_restart(inflate_handle)) then
         call get_copy_owner_index(ss_inflate_mean_index, owner, owners_index)
         if (owner == ens_handle%my_pe) then 
            call prepare_to_write_to_vars(ens_handle)
            ens_handle%vars(:, owners_index) = get_inflate_mean(inflate_handle)
         endif

      endif
      if (.not. sd_from_restart(inflate_handle)) then

        call get_copy_owner_index(ss_inflate_sd_index, owner, owners_index)
        if (owner == ens_handle%my_pe)  then
           call prepare_to_write_to_vars(ens_handle)
           ens_handle%vars(:, owners_index) = get_inflate_sd(inflate_handle)

         endif

      endif
   endif

! Misleading comment - the computation uses the whole array.
   ! Inflation type 3 is spatially-constant.  Make sure the entire array is set to that
   ! value. the computation only uses index 1, but the diagnostics write out the entire
   ! array and it will be misleading if not constant.  the inf values were set above.  
   ! if they were set by namelist, this code changes nothing.  but if they were read in
   ! from a file, then it is possible the values vary across the array.  these lines
   ! ensure the entire array contains a single constant value to match what the code uses.
   if(do_single_ss_inflate(inflate_handle)) then
      call get_copy_owner_index(ss_inflate_mean_index, owner, owners_index)
      if (owner == ens_handle%my_pe) then
         ens_handle%vars(:, owners_index) = ens_handle%vars(1, owners_index)
      endif
      call get_copy_owner_index(ss_inflate_sd_index, owner, owners_index)
      if (owner == ens_handle%my_pe) then
         ens_handle%vars(:, owners_index) = ens_handle%vars(1, owners_index)
      endif
   endif

endif

end subroutine read_state_space_inflation


!-----------------------------------------------------------------------

subroutine write_state_space_inflation(inflate_handle, file_info, state_ens_handle, ss_inflate_mean_index, &
   ss_inflate_sd_index)

type(adaptive_inflate_type), intent(in)    :: inflate_handle
type(file_info_type),        intent(in)    :: file_info
type(ensemble_type),         intent(in)    :: state_ens_handle
integer,                     intent(in)    :: ss_inflate_mean_index, ss_inflate_sd_index

! check whether file_info handle is initialized
call assert_file_info_initiailzed(file_info, 'write_state_space_inflation')

if(output_inf_restart(inflate_handle)) then
   ! Use the ensemble manager to output restart for state space (flavors 2 or 3)
   if(do_ss_inflate(inflate_handle)) then

      call write_ensemble_restart(state_ens_handle, file_info, get_inflation_out_filename(inflate_handle), &
           ss_inflate_mean_index, ss_inflate_sd_index, force_single_file = .true.)
   endif
endif

end subroutine write_state_space_inflation


!-----------------------------------------------------------------------

subroutine read_ensemble_restart(ens_handle, start_copy, end_copy, file_name, init_time, force_single_file)

type(ensemble_type),  intent(inout)           :: ens_handle
integer,              intent(in)              :: start_copy, end_copy
character(len=*),     intent(in)              :: file_name
type(time_type),      intent(in),    optional :: init_time
logical,              intent(in),    optional :: force_single_file

! The ensemble being read from restart is stored in copies start_copy:end_copy
! contiguously (by fiat). So if start_copy is 6, the first ensemble restart
! goes in copy 6, the second in copy 7, etc. This lets higher level determine
! where other copies like mean, inflation, etc., are stored.

! Avoid num_vars size storage on stack; make this allocatable from heap
real(r8), allocatable               :: ens(:) 
integer                             :: iunit, i
character(len=LEN(file_name) + 5)   :: this_file_name
character(len=4)                    :: extension
type(time_type)                     :: ens_time
integer                             :: global_copy_index
logical                             :: single_file_override

! the code reads into the vars array
!ens_handle%valid = VALID_VARS

! Some compilers (absoft, but others also) are particularly unhappy about
! both checking present(N) _and_ evaluating N inside a single if() test.
! (It evaluates both at the same time and blows up on the not present value.)
! The standard says they do not have to evaluate right to left.  Common error
! for anyone with a C programming background.   So -- set a separate local
! logical variable which always has a value, whether the arg is present or not.
if (present(force_single_file)) then
  single_file_override = force_single_file
else
  single_file_override = .false.
endif

!-------- Block for single restart file -----
if(single_file_override) then ! Single restart file is read only by task 0

   if (my_task_id() == 0) then
      allocate(ens(ens_handle%num_vars))
      iunit = open_restart_read(file_name)
   else
      allocate(ens(1)) ! dummy for put_copy
   endif

   do i = start_copy, end_copy
      if(my_task_id() == 0) then
         call aread_state_restart(ens_time, ens, iunit)
      endif

      ! overwrite ens_time if needed
      if(present(init_time)) ens_time = init_time
      ! Store this copy in the appropriate place on the appropriate process
      ! Ensemble manager function so give pe not task.
      call put_copy(map_task_to_pe(ens_handle,0), ens_handle, i, ens, ens_time)

   enddo

   deallocate(ens)
   if (my_task_id() == 0) call close_restart(iunit)

else

!----------- Block for multiple restart files -----------
   ! Loop to read in all my ensemble members
   READ_MULTIPLE_RESTARTS: do i = 1, ens_handle%my_num_copies
      ! Get global index for my ith ensemble
      global_copy_index = ens_handle%my_copies(i)
      ! If this global copy is in the range being read in, proceed
      if(global_copy_index >= start_copy .and. global_copy_index <= end_copy) then
         ! File name extension is the global index number for the copy
         write(extension, '(i4.4)') global_copy_index - start_copy + 1
         this_file_name = trim(file_name) // '.' // extension
         iunit = open_restart_read(this_file_name)

         ! Read the file directly into storage
         call aread_state_restart(ens_handle%time(i), ens_handle%vars(:, i), iunit)
         if(present(init_time)) ens_handle%time(i) = init_time
         ! Close the restart file
         call close_restart(iunit)

      endif

   end do READ_MULTIPLE_RESTARTS
endif

end subroutine read_ensemble_restart


!-----------------------------------------------------------------------

subroutine write_ensemble_restart(ens_handle, file_info, file_name, start_copy, end_copy, &
   force_single_file)

type(ensemble_type),  intent(in)    :: ens_handle
type(file_info_type), intent(in)    :: file_info
character(len=*),     intent(in)    :: file_name
integer,              intent(in)    :: start_copy, end_copy
logical, optional,    intent(in)    :: force_single_file

! Large temporary storage to be avoided if possible
real(r8), allocatable               :: ens(:)
type(time_type)                     :: ens_time
integer                             :: iunit, i, global_index
integer                             :: owner, owners_index
character(len=LEN(file_name) + 10)  :: this_file_name
character(len=4)                    :: extension
logical                             :: single_file_forced

! check whether file_info handle is initialized
call assert_file_info_initiailzed(file_info, 'write_ensemble_restart')

if (present(force_single_file) ) then
   single_file_forced = force_single_file
else
   single_file_forced = .FALSE.
endif

! Error checking
!if (ens_handle%valid /= VALID_VARS .and. ens_handle%valid /= VALID_BOTH) then
!   call error_handler(E_ERR, 'write_ensemble_restart', &
!        'last access not var-complete', source, revision, revdate)
!endif

! For single file, need to send restarts to pe0 and it writes them out.
!-------------- Block for single_restart file -------------
! Need to force single restart file for inflation files
if(get_write_to_single_file(file_info) .or. single_file_forced) then

   ! Single restart file is written only by task 0
   if(my_task_id() == 0) then

      iunit = open_restart_write(file_name)

      ! Loop to write each ensemble member
      allocate(ens(ens_handle%num_vars))   ! used to be on stack.
      do i = start_copy, end_copy
         ! Figure out where this ensemble member is being stored
         call get_copy_owner_index(i, owner, owners_index)
         ! If it's on task 0, just write it
         if(map_pe_to_task(ens_handle, owner) == 0) then
            call awrite_state_restart(ens_handle%time(owners_index), &
               ens_handle%vars(:, owners_index), iunit)
         else
            ! Get the ensemble from the owner and write it out
            ! This communication assumes index numbers are monotonically increasing
            ! and that communications is blocking so that there are not multiple
            ! outstanding messages from the same processors (also see send_to below).
            call receive_from(map_pe_to_task(ens_handle, owner), ens, ens_time)
            call awrite_state_restart(ens_time, ens, iunit)
         endif
      end do
      deallocate(ens)
      call close_restart(iunit)
   else ! I must send my copies to task 0 for writing to file
      do i = 1, ens_handle%my_num_copies
         ! Figure out which global index this is
         global_index = ens_handle%my_copies(i)
         ! Ship this ensemble off to the master
         if(global_index >= start_copy .and. global_index <= end_copy) &
            call send_to(0, ens_handle%vars(:, i), ens_handle%time(i))
      end do
   endif

else
!-------------- Block for multiple restart files -------------
   ! Everyone can just write their own files
   do i = 1, ens_handle%my_num_copies
      ! Figure out which global index this is
      global_index = ens_handle%my_copies(i)
      if(global_index >= start_copy .and. global_index <= end_copy) then
         write(extension, '(i4.4)') ens_handle%my_copies(i)
         this_file_name = trim(file_name) // '.' // extension
         iunit = open_restart_write(this_file_name)
         call awrite_state_restart(ens_handle%time(i), ens_handle%vars(:, i), iunit)
         call close_restart(iunit)
      endif
   end do
endif

end subroutine write_ensemble_restart


!-----------------------------------------------------------------------
!> Write a restart file given a model extended state and a unit number 
!> opened to the restart file. (Need to reconsider what is passed to 
!> identify file or if file can even be opened within this routine).

subroutine awrite_state_restart(model_time, model_state, funit, target_time)

type(time_type), intent(in)           :: model_time
real(r8),        intent(in)           :: model_state(:)
integer,         intent(in)           :: funit
type(time_type), optional, intent(in) :: target_time

integer :: i, io, rc
character(len=16)  :: open_format
character(len=256) :: filename
logical :: is_named

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! Figure out whether the file is opened FORMATTED or UNFORMATTED
inquire(funit, FORM=open_format)

! assume success
io = 0

! Write the state vector
if (ascii_file_format(open_format)) then
   if(present(target_time)) call write_time(funit, target_time, ios_out=io)
   if (io /= 0) goto 10
   call write_time(funit, model_time, ios_out=io)
   if (io /= 0) goto 10
   do i = 1, size(model_state)
      write(funit, *, iostat = io) model_state(i)
      if (io /= 0) goto 10
   end do
else
   if(present(target_time)) call write_time(funit, target_time, form="unformatted", ios_out=io)
   if (io /= 0) goto 10
   call write_time(funit, model_time, form="unformatted", ios_out=io)
   if (io /= 0) goto 10
   write(funit, iostat = io) model_state
   if (io /= 0) goto 10
endif

! come directly here on error. 
10 continue

! if error, use inquire function to extract filename associated with
! this fortran unit number and use it to give the error message context.
if (io /= 0) then
   inquire(funit, named=is_named, name=filename, iostat=rc)
   if ((rc /= 0) .or. (.not. is_named)) filename = 'unknown file'
   write(msgstring,*) 'error writing to restart file ', trim(filename)
   call error_handler(E_ERR,'awrite_state_restart',msgstring,source,revision,revdate)
endif


end subroutine awrite_state_restart


!-----------------------------------------------------------------------
!> Closes a restart file

subroutine close_restart(file_unit)

integer, intent(in) :: file_unit

call close_file(file_unit)

end subroutine close_restart


!-----------------------------------------------------------------------
!> Read a restart file given a unit number

subroutine aread_state_restart(model_time, model_state, funit, target_time)

type(time_type), intent(out)            :: model_time
real(r8),        intent(out)            :: model_state(:)
integer,         intent(in)             :: funit
type(time_type), optional, intent(out) :: target_time

character(len=16) :: open_format
integer :: ios, int1, int2

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

ios = 0

! Figure out whether the file is opened FORMATTED or UNFORMATTED
inquire(funit, FORM=open_format)

if (ascii_file_format(open_format)) then
   if(present(target_time)) target_time = read_time(funit)
   model_time = read_time(funit)
   read(funit,*,iostat=ios) model_state
else
   if(present(target_time)) target_time = read_time(funit, form = "unformatted")
   model_time = read_time(funit, form = "unformatted")
   read(funit,iostat=ios) model_state
endif

! If the model_state read fails ... dump diagnostics.
if ( ios /= 0 ) then
   ! messages are being used as error lines below.  in an MPI filter,
   ! all messages are suppressed that aren't from PE0.  if an error
   ! happens in another task, these lines won't be printed unless we
   ! turn on output.
   call set_output(.true.)

   write(msgstring,*)'dimension of model state is ',size(model_state)
   call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)

   if(present(target_time)) then
      call get_time(target_time, int1, int2)       ! time -> secs/days
      write(msgstring,*)'target_time (secs/days) : ',int1,int2
      call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)
   endif

   call get_time(model_time, int1, int2)       ! time -> secs/days
   write(msgstring,*)'model_time (secs/days) : ',int1,int2
   call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)

   write(msgstring,'(''model max/min/first is'',3(1x,E12.6) )') &
            maxval(model_state), minval(model_state), model_state(1)
   call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)

   call dump_unit_attributes(funit)

   write(msgstring,*)'read error is : ',ios
   call error_handler(E_ERR,'aread_state_restart',msgstring,source,revision,revdate)
endif

end subroutine aread_state_restart


!-----------------------------------------------------------------------
!> Opens a restart file for writing

function open_restart_write(file_name, override_write_format)

character(len=*), intent(in) :: file_name
character(len=*), optional, intent(in) :: override_write_format

integer :: open_restart_write, io

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

open_restart_write = get_unit()
if(present(override_write_format)) then
   open(unit = open_restart_write, file = file_name, form = override_write_format, &
        iostat = io)
else
   open(unit = open_restart_write, file = file_name, form = write_format, iostat = io)
endif
if (io /= 0) then
   write(msgstring,*) 'unable to open restart file ', trim(file_name), ' for writing'
   call error_handler(E_ERR,'open_restart_write',msgstring,source,revision,revdate)
endif

end function open_restart_write


!-----------------------------------------------------------------------
!> Opens a restart file for reading

function open_restart_read(file_name)

integer :: open_restart_read
character(len=*), intent(in) :: file_name

integer :: ios, ios_out
!!logical :: old_output_state
type(time_type) :: temp_time
character(len=64) :: string2

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! DEBUG -- if enabled, every task will print out as it opens the
! restart files.  If questions about missing restart files, first start
! by commenting in only the timestamp line.  If still concerns, then
! go ahead and comment in all the lines.
!!old_output_state = do_output()
!!call set_output(.true.)
!call timestamp("open_restart", "opening restart file "//trim(file_name), pos='')
!!call set_output(old_output_state)
!END DEBUG

! if you want to document which file(s) are being opened before
! trying the open (e.g. in case the fortran runtime library intercepts
! the error and does not return to let us print out the name) then
! comment this in and you can see what files are being opened.
!write(msgstring, *) 'Opening restart file ',trim(adjustl(file_name))
!call error_handler(E_MSG,'open_restart_read',msgstring,source,revision,revdate)

! WARNING: Absoft Pro Fortran 9.0, on a power-pc mac, is convinced
! that certain binary files are, in fact, ascii, because the read_time 
! call is returning what seems like a good time even though it should
! be garbage.  This code works fine on all other platforms/compilers
! we've tried, so we're leaving it as-is.  Best solution if you're
! using absoft on a mac is to set all files to be non-binary in the
! namelist.  You may also have to set the format in both obs_model_mod.f90 
! and interpolate_model.f90 to 'formatted' instead of the hardcoded 
! 'unformatted' for async 2/4 model advance temp_ic and temp_ud files.

! Autodetect format of restart file when opening
! Know that the first thing in here has to be a time, so try to read it.
! If it fails with one format, try the other. If it fails with both, punt.
open_restart_read = get_unit()
read_format = 'formatted'
open(unit   = open_restart_read, &
     file   = trim(file_name),   &
     form   = read_format,       &
     action = 'read',            &
     status = 'old',             &
     iostat = ios)
! An opening error means something is wrong with the file, error and stop
if(ios /= 0) goto 11
temp_time = read_time(open_restart_read, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be formatted, proceed
   rewind open_restart_read
   return
endif

! Next, try to see if an unformatted read works instead
close(open_restart_read)

open_restart_read = get_unit()
read_format = 'unformatted'
open(unit   = open_restart_read, &
     file   = trim(file_name),   &
     form   = read_format,       &
     action = 'read',            &
     status = 'old',             &
     iostat = ios)
! An opening error means something is wrong with the file, error and stop
if(ios /= 0) goto 11
rewind open_restart_read
temp_time = read_time(open_restart_read, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be unformatted, proceed
   rewind open_restart_read
   return
endif

! Otherwise, neither format works. Have a fatal error.
11 continue

write(msgstring, *) 'Problem opening file ',trim(file_name)
write( string2 , *) 'OPEN status was ',ios
call error_handler(E_ERR, 'open_restart_read', msgstring, &
     source, revision, revdate, text2=string2)

end function open_restart_read

!> @}
!-------------------------------------------------------
end module state_vector_io_mod
