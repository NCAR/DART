! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module state_vector_io_mod

!> \defgroup state_vector_io_mod state_vector_io_mod
!> @{ \brief Routines for reading and writing the model state.
!>
!> Read_state() and write_state() are the major routines in this module.
!> write_restart_direct()
!> is used in diagnostic files writes (for certain diagnostic options).
!> They read/write state from an ensmeble_handle%copies array to files described in a file_info_type.
!> The inflation handles (prior and posterior) are optional arguments to read/write_state since
!> inflation is used in filter but not in perfect_model_obs.
!>
!> Only supporting netcdf reads and writes
!> 
!> netcdf files use \ref io_filenames_mod. See model_mod for construct_file_name_in() and
!> io_filenames_mod for construct_file_name_out(). See \ref io_filenames_mod for filenames 
!> for extra copies (mean, etc.)
!> Each domain is written to a separate file.
!>
!> Perturbation (creating an ensemble from a single instance) is done separately from the read.
!> perturb_from_single_copy is passed as an argument to read_state. When reading netcdf files
!> only the first ensemble member is read and transposed.  
!> create_ensemble_from_single_file()
!> is then called in filter to generate the ensemble.

use adaptive_inflate_mod, only : adaptive_inflate_type, mean_from_restart, sd_from_restart, &
                                 do_single_ss_inflate, do_varying_ss_inflate, &
                                 get_inflate_mean, get_inflate_sd, do_ss_inflate, &
                                 get_is_prior, get_is_posterior, get_inflation_mean_copy, &
                                 get_inflation_sd_copy

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
                                 register_module, set_output, to_upper

use assim_model_mod,      only : assim_model_type

use time_manager_mod,     only : time_type, read_time, write_time, &
                                 get_time

use io_filenames_mod,     only : get_restart_filename, file_info_type, &
                                 get_single_file, get_stage_metadata, &
                                 assert_file_info_initialized, stage_metadata_type, &
                                 assert_restart_names_initialized

!>@todo FIXME This should go through assim_model_mod
use model_mod,            only : read_model_time

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
public :: set_stage_to_write, &
          get_stage_to_write

! Module storage for writing error messages
character(len=512) :: msgstring

!>@todo FIXME  these should not be hardcoded this way.
!> the more general solution is to have a get_output_stage() function
!> that takes a string, and compares the string against a list of
!> valid strings.  the strings should be stored in an array and
!> the interpretation of what they mean should NOT be enforced here.
!> the calling code knows what they are for; it should be a black box
!> here.

integer, parameter :: MAX_STAGES = 4

type stages_to_write
   logical            :: write_stage(MAX_STAGES) = .false.
   character(len=32)  ::  stage_name(MAX_STAGES) = 'null'
   integer            :: my_num_stages = 0
end type stages_to_write

type(stages_to_write) :: global_stages_to_write

! Logical flag for initialization of module
logical :: module_initialized = .false.

! namelist variables with default values
! Aim: to have the regular transpose as the default
integer :: limit_mem = HUGE(1_i4) !< This is the number of elements (not bytes) 
                                  !< so you don't have times the number by 4 or 8
logical :: single_precision_output = .false. !< Allows you to write r4 netcdf files 
                                             !< even if filter is double precision

!>@todo FIXME there should be a 'variable_by_variable' logical here
!> that says do i/o a single variable at a time instead of the entire
!> state.  the limit_mem can stay, but it's much harder to use for
!> the average user.

namelist /  state_vector_io_nml / limit_mem, single_precision_output

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

endif

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=state_vector_io_nml)
if (do_nml_term()) write(     *     , nml=state_vector_io_nml)

end subroutine state_vector_io_init


!-----------------------------------------------------------------------
!> Read in the state vectors into an ensemble handle
!> The copies will end up in the state_ens_handle%copies array.
!> If read_time_from_file = .true. then time is overwritten by the time read from 
!> the restart. If read_time_from_file = .false. then the time associated with 
!> each ensemble member is set to time.
!> Inflation handles are optional (used by filter, not by perfect_model_obs)
!>   - Note the user must give both inflation handles or neither.
!> perturb_from_single_copy is optional (used by filter, not by perfect_model_obs)
!>   - only 1 restart file is read in. For netcdf read only one copy is transposed.
!>     Still doing a full all vars to all copies for dart format read though.


subroutine read_state(state_ens_handle, file_info, read_time_from_file, time, &
            prior_inflate_handle, post_inflate_handle)

type(ensemble_type),         intent(inout) :: state_ens_handle
type(file_info_type),        intent(in)    :: file_info
logical,                     intent(in)    :: read_time_from_file ! state time
type(time_type),             intent(inout) :: time
type(adaptive_inflate_type), optional, intent(in) :: prior_inflate_handle
type(adaptive_inflate_type), optional, intent(in) :: post_inflate_handle

logical :: inflation_handles = .false.

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! check whether file_info handle is initialized
call assert_file_info_initialized(file_info, 'read_state')

! check that we either have both inflation handles or neither:
if ( present(prior_inflate_handle) .neqv. present(post_inflate_handle) ) then
   call error_handler(E_ERR, 'read_state', 'must have both inflation handles or neither', &
          source,revision,revdate)
endif

!>@todo all-or-nothing WHY?
if (present(prior_inflate_handle) .and. present(post_inflate_handle)) inflation_handles = .true.

! inflation read from netcdf files only for state_space_inflation
! Ideally would only do a read transpose for VARYING state_space_inflation
!>@todo FIXME do we need a separate file_info for inflation

!>@todo FIXME understand this ... what is it doing and why.
! Filling copies array for single state space inflation values if required.
if (inflation_handles) then

   !>@todo just a print at the moment
   call print_inflation_source(file_info, prior_inflate_handle, 'Prior')
   call print_inflation_source(file_info, post_inflate_handle,  'Posterior')

   ! If inflation is single state space read from a file, the copies array is filled here.
   call fill_single_ss_inflate_from_read(state_ens_handle, prior_inflate_handle, post_inflate_handle)

   ! If inflation is from a namelist value it is set here.
   call fill_ss_from_nameslist_value(state_ens_handle, prior_inflate_handle)
   call fill_ss_from_nameslist_value(state_ens_handle, post_inflate_handle)

endif

if (get_single_file(file_info)) then
   call error_handler(E_ERR,'read_state:', &
      'Writing single file restarts is in progress')
else
   call read_restart_direct(state_ens_handle, file_info, read_time_from_file, time)
endif

end subroutine read_state


!-----------------------------------------------------------------------
!> Write state vectors from an ensemble_handle
!> The inflation handles are options (used by filter, not by perfect_model_obs)
!> Note both inflation handles must be present (or neither)
!> There is logic in write_state about spare_copies. These are extra state vector 
!> copies that may be in the ensemble handle when doing large model, single timestep runs.


subroutine write_state(state_ens_handle, file_info)

type(ensemble_type),                   intent(inout) :: state_ens_handle
type(file_info_type),                  intent(in)    :: file_info

type(stage_metadata_type) :: output_files


if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! check whether file_info handle is initialized
call assert_file_info_initialized(file_info, 'write_state')

! do this once
output_files = get_stage_metadata(file_info)

! write ensemble copies
call write_restart_direct(state_ens_handle, output_files)

end subroutine write_state


!----------------------------------------------------------------------
!> Read the restart information directly from the model netcdf file


subroutine read_restart_direct(state_ens_handle, file_info, use_time_from_file, time)

type(ensemble_type),  intent(inout) :: state_ens_handle
type(file_info_type), intent(in)    :: file_info
logical,              intent(in)    :: use_time_from_file
type(time_type),      intent(inout) :: time

integer :: dart_index ! where to start in state_ens_handle%copies
integer :: domain     ! loop index
type(stage_metadata_type) :: restart_files

! check whether file_info handle is initialized
call assert_file_info_initialized(file_info, 'read_restart_direct')

! do this once
restart_files = get_stage_metadata(file_info)

! read time from input file if time not set in namelist
!>@todo Check time constistency across files? This is assuming they are consistent.
if(use_time_from_file) then
   time = read_model_time(get_restart_filename(restart_files, 1, 1)) ! Any of the restarts?
endif

state_ens_handle%time = time

! read in the data and transpose
dart_index = 1 ! where to start in state_ens_handle%copies - this is modified by read_transpose
do domain = 1, get_num_domains()
   call read_transpose(state_ens_handle, restart_files, domain, dart_index, limit_mem)
enddo

! Need Temporary print of initial model time?

end subroutine read_restart_direct


!-----------------------------------------------------------------------
!> write the restart information directly into the model netcdf file.


subroutine write_restart_direct(state_ens_handle, file_name_handle)

type(ensemble_type),    intent(inout) :: state_ens_handle
type(stage_metadata_type), intent(in) :: file_name_handle

integer :: dart_index !< where to start in state_ens_handle%copies
integer :: domain !< loop index

if ( .not. module_initialized ) call state_vector_io_init() ! to read the namelist

! check whether file_info handle is initialized
call assert_restart_names_initialized(file_name_handle, 'write_restart_direct')

! transpose and write out the data
dart_index = 1

! Different filenames for prior vs. posterior vs. diagnostic files
do domain = 1, get_num_domains()
   call transpose_write(state_ens_handle, file_name_handle, domain, &
                  dart_index, limit_mem, single_precision_output)
enddo

end subroutine write_restart_direct


!-----------------------------------------------------------------------
!> Single state space inflation from file value is only index 1 of 
!> inflation vars array (mean and sd).
!> This routine find the owner of the 1st element in the vars array
!> The owner broadcasts the values of single state space inflation 
!> to all other tasks who then update their copies array.
!> Note filling both mean and sd values if at least one of mean
!> or sd is read from file.  If one is set from a namelist the copies 
!> array is overwritten in fill_ss_from_namelist value

!>@todo We need to refactor this, not sure where it is being used
subroutine fill_single_ss_inflate_from_read(ens_handle, &
                prior_inflate_handle, post_inflate_handle)

type(ensemble_type),         intent(inout) :: ens_handle
type(adaptive_inflate_type), intent(in)    :: prior_inflate_handle
type(adaptive_inflate_type), intent(in)    :: post_inflate_handle

integer :: owner, owners_index
integer(i8) :: first_element
real(r8), allocatable :: inf_array(:) ! 2 or 4 values
integer :: inf_count
logical :: return_me ! flag to return if not read any inflation values from files

integer :: PRIOR_INF_MEAN
integer :: PRIOR_INF_SD
integer :: POST_INF_MEAN
integer :: POST_INF_SD

! Return if not doing single state space inflation
if (.not. do_single_ss_inflate(prior_inflate_handle) .and. &
    .not. do_single_ss_inflate(post_inflate_handle)) return

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

PRIOR_INF_MEAN  = get_inflation_mean_copy(prior_inflate_handle)
PRIOR_INF_SD    = get_inflation_sd_copy(  prior_inflate_handle)
POST_INF_MEAN   = get_inflation_mean_copy(post_inflate_handle)
POST_INF_SD     = get_inflation_sd_copy(  post_inflate_handle)

if (ens_handle%my_pe == owner) then
   if (do_single_ss_inflate(prior_inflate_handle) .and. &
       do_single_ss_inflate(post_inflate_handle)) then
      inf_array(1) = ens_handle%copies(PRIOR_INF_MEAN, owners_index)
      inf_array(2) = ens_handle%copies(PRIOR_INF_SD, owners_index)
      inf_array(3) = ens_handle%copies(POST_INF_MEAN, owners_index)
      inf_array(4) = ens_handle%copies(POST_INF_SD, owners_index)

   elseif (do_single_ss_inflate(post_inflate_handle) .and. &
     .not. do_single_ss_inflate(post_inflate_handle)) then

      inf_array(1) = ens_handle%copies(PRIOR_INF_MEAN, owners_index)
      inf_array(2) = ens_handle%copies(PRIOR_INF_SD, owners_index)

   elseif(.not. do_single_ss_inflate(post_inflate_handle) .and. &
                do_single_ss_inflate(post_inflate_handle)) then

      inf_array(1) = ens_handle%copies(POST_INF_MEAN, owners_index)
      inf_array(2) = ens_handle%copies(POST_INF_SD, owners_index)

   endif

   call broadcast_send(map_pe_to_task(ens_handle, owner), inf_array)

else

   call broadcast_recv(map_pe_to_task(ens_handle, owner), inf_array)

   if (do_single_ss_inflate(prior_inflate_handle) .and. &
       do_single_ss_inflate(post_inflate_handle)) then

      ens_handle%copies(PRIOR_INF_MEAN, owners_index)     = inf_array(1)
      ens_handle%copies(PRIOR_INF_SD, owners_index)  = inf_array(2)
      ens_handle%copies(POST_INF_MEAN, owners_index)      = inf_array(3)
      ens_handle%copies(POST_INF_SD, owners_index)   = inf_array(4)

   elseif (do_single_ss_inflate(prior_inflate_handle) .and. &
     .not. do_single_ss_inflate(post_inflate_handle)) then

      ens_handle%copies(PRIOR_INF_MEAN, owners_index)    = inf_array(1)
      ens_handle%copies(PRIOR_INF_SD, owners_index) = inf_array(2)

   elseif(.not. do_single_ss_inflate(prior_inflate_handle) .and. &
                do_single_ss_inflate(post_inflate_handle)) then

      ens_handle%copies(POST_INF_MEAN, owners_index)    = inf_array(1)
      ens_handle%copies(POST_INF_SD, owners_index) = inf_array(2)

   endif

endif

end subroutine fill_single_ss_inflate_from_read


!-----------------------------------------------------------------------
!> Check whether inflation values come from namelist and
!> fill copies array with namelist values for inflation if they do.


subroutine fill_ss_from_nameslist_value(ens_handle, inflate_handle)

type(ensemble_type),         intent(inout) :: ens_handle
type(adaptive_inflate_type), intent(in)    :: inflate_handle

character(len=32)  :: label
character(len=128) :: nmread
integer            :: INF_MEAN_COPY, INF_SD_COPY
real(r8)           :: inf_initial, sd_initial

INF_MEAN_COPY = get_inflation_mean_copy(inflate_handle)
INF_SD_COPY   = get_inflation_sd_copy(  inflate_handle)

! To match Lanai filter_state_space_diagnostics, 
! if not doing inflation set inf_mean = 1, inf_sd = 0
if (.not. do_ss_inflate(inflate_handle)) then
   ens_handle%copies(INF_MEAN_COPY, :) = 1.0_r8
   ens_handle%copies(INF_SD_COPY, :)   = 0.0_r8
   return
endif

if(get_is_prior(inflate_handle)) then
   label = "Prior"
else if (get_is_posterior(inflate_handle)) then
   label = "Posterior"
else
   write(msgstring, *) "state space inflation but neither prior or posterior"
   call error_handler(E_ERR, 'fill_ss_from_nameslist_value', msgstring, &
      source, revision, revdate)
endif

if (.not. mean_from_restart(inflate_handle)) then
   inf_initial = get_inflate_mean(inflate_handle)
   write(nmread, '(A, F12.6)') 'mean read from namelist as ', inf_initial
   call error_handler(E_MSG, trim(label) // ' inflation:', trim(nmread), &
      source, revision, revdate)
   ens_handle%copies(INF_MEAN_COPY, :) = inf_initial
endif

if (.not. sd_from_restart(inflate_handle)) then
   sd_initial = get_inflate_sd(inflate_handle)
   write(nmread, '(A, F12.6)') 'sd read from namelist as ', sd_initial
   call error_handler(E_MSG, trim(label) // ' inflation:', trim(nmread), &
      source, revision, revdate)
   ens_handle%copies(INF_SD_COPY, :) = sd_initial
endif

end subroutine fill_ss_from_nameslist_value

!-----------------------------------------------------------
!> set stage names


subroutine set_stage_to_write(stage_name_in, output_stage)
character(len=*),  intent(in) :: stage_name_in
logical,           intent(in) :: output_stage

integer :: nstages

nstages = global_stages_to_write%my_num_stages + 1

global_stages_to_write%stage_name( nstages) = stage_name_in
global_stages_to_write%write_stage(nstages) = output_stage
global_stages_to_write%my_num_stages        = nstages

end subroutine set_stage_to_write

!-----------------------------------------------------------
!> set stage names


function get_stage_to_write(stage_name_in) result (write_this_stage)
character(len=*), intent(in) :: stage_name_in
logical :: write_this_stage

character(len=32) :: input_stage, global_stage
integer :: i, nstages

write_this_stage = .false.
nstages     = global_stages_to_write%my_num_stages 

input_stage = stage_name_in
call to_upper(input_stage)

do i = 1, nstages
   global_stage = global_stages_to_write%stage_name(i)
   call to_upper(global_stage)
   if (input_stage == global_stage) then
      write_this_stage = global_stages_to_write%write_stage(i)
      return
   endif
enddo

end function get_stage_to_write

!-----------------------------------------------------------------------


subroutine print_inflation_source(file_info, inflate_handle, label)

type(file_info_type),        intent(in) :: file_info
type(adaptive_inflate_type), intent(in) :: inflate_handle
character(len=*),            intent(in) :: label  

integer :: INF_MEAN_COPY
integer :: INF_SD_COPY
integer :: idom
character(len=256) :: fname

type(stage_metadata_type) :: restart_files

! do this once
restart_files = get_stage_metadata(file_info)

do idom = 1, get_num_domains()
   if (do_ss_inflate(inflate_handle)) then
      if (mean_from_restart(inflate_handle)) then
         ! Format for input prior inflation mean
         !    input_priorinf_mean[_d0X].nc'
         INF_MEAN_COPY = get_inflation_mean_copy(inflate_handle)
         fname = get_restart_filename(restart_files, INF_MEAN_COPY, idom)   
         write(msgstring,*) 'mean read from restart file: ' // trim(fname)
         call error_handler(E_MSG, trim(label) // ' inflation:', trim(msgstring), &
            source, revision, revdate)
      endif

      if (sd_from_restart(inflate_handle)) then  
         ! Format for input prior inflation standard deviation
         !    input_priorinf_sd[_d0X].nc'
         INF_SD_COPY = get_inflation_sd_copy(inflate_handle)
         fname = get_restart_filename(restart_files, INF_SD_COPY, idom)   
         write(msgstring,*) 'sd read from restart file: ' // trim(fname)
         call error_handler(E_MSG, trim(label) // ' inflation:', trim(msgstring), &
            source, revision, revdate)
      endif
   endif
enddo

end subroutine print_inflation_source

!-----------------------------------------------------------------------
!> @}
!-------------------------------------------------------
end module state_vector_io_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
