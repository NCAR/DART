! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module direct_netcdf_mod
!> \defgroup direct_netcdf_mod direct_netcdf_mod
!> @{ \brief Routines for the limited transpose.
!>
!> Netcdf IO for a domain.
!> The idea is to be generic rather than having a converter for each module.
!> Also aim to go to the filesystem only once, e.g.
!> \verbatim
!>      wrf => dart => wrf
!> \endverbatim
!> rather than
!> \verbatim
!>      wrf => wrf_to_dart => dart => dart_to_wrf => wrf
!> \endverbatim
!>
!> Every task needs the dimensions of each state variable to calculate
!> what it is going to recieve in the IO transpose (this is calculated in add_domain )
!> \par Aim of the limited transpose:
!>
!>To limit how much of the state vector you can read at once.
!> You can limit the transpose by memory using <code>limit_mem</code>.
!>
!>What you (potentially) gain from this:
!>
!>* Don't have to have the whole state vector.
!>* Don't have to use a parallel IO library.
!>
!>If limit_mem > state vector size you have the regular transpose + IO, except:
!>  1. You are reading directly from a netcdf file, not a dart state vector file.
!>  2. You only transpose the copies that are being written/read.
!>
!> This code has multiple places where round-robin layout of state onto task is assumed.
!> These are in the section labelled: \n
!> \verbatim
!>    !--------------------------------------------------------
!>    ! Routines that are making the assumption that the ensemble
!>    ! distribution is round-robin (distribution type 1)
!>    !--------------------------------------------------------
!> \endverbatim
!> read_transpose() is the read routine.
!> transpose_write() is the write routine.
!> Note dart_index is an inout variable to read_transpose() and transpose_write(). 
!> This is making the assumption that the calling code is using dart_index in the following way:
!>  * dart_index going in to the subroutines is where the domain starts in the state vector.
!>  * dart_index coming out of the subroutines is where the domain ends.
!>  * The domain is contiguous in the state_vector.
!>
!> Note there was a <code>limit_procs</code> in the limited transpose. This was removed in
!> svn commit 9456.

use types_mod,            only : r4, r8, i4, MISSING_R8, MISSING_R4, MISSING_I

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, &
                                 get_copy_owner_index, get_ensemble_time

use time_manager_mod,     only : time_type

use utilities_mod,        only : error_handler, nc_check, &
                                 E_MSG, E_ALLMSG, E_ERR, file_exist

use state_structure_mod,  only : get_num_variables, get_sum_variables,  &
                                 get_sum_variables_below, &
                                 get_variable_name, get_io_clamping_maxval,   &
                                 get_io_clamping_minval, do_io_clamping,         &
                                 get_io_num_dims, get_io_dim_lengths,            &
                                 get_variable_size, get_io_num_unique_dims,   &
                                 get_io_unique_dim_name, get_dim_name,        &
                                 get_io_unique_dim_length, get_num_variables, &
                                 set_var_id, get_domain_size, do_io_update, &
                                 get_units, get_long_name, get_short_name, &
                                 get_has_missing_value, get_FillValue, &
                                 get_missing_value, get_add_offset, get_scale_factor, &
                                 get_xtype, get_num_domains

use io_filenames_mod,     only : get_input_file, get_output_file, restart_names_type, &
                                 get_file_description

use copies_on_off_mod,    only : query_read_copy, query_write_copy, &
                                 is_inflation_copy, is_mean_copy, is_sd_copy, &
                                 is_ensemble_copy

use model_mod,            only : write_model_time

use netcdf

implicit none
private
public :: read_transpose, transpose_write

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Limit_mem is a namelist item in state_vector_io_mod.

integer :: ret !< netcdf return code

character(len=512) :: msgstring

contains

!-------------------------------------------------
! Limited transpose code
! There are two versions of the limited transpose code:
!   * Single processor (no memory limit is applied)
!   * Multi processor (memory limit applied)
!-------------------------------------------------
subroutine read_transpose(state_ens_handle, name_handle, domain, dart_index, limit_mem)

type(ensemble_type),      intent(inout) :: state_ens_handle
type(restart_names_type), intent(in)    :: name_handle
integer,                  intent(in)    :: domain
integer,                  intent(inout) :: dart_index !< This is for mulitple domains
integer,                  intent(in)    :: limit_mem !< How many state elements you can read at once

if (task_count() == 1) then
   call read_transpose_single_task(state_ens_handle, name_handle, domain, dart_index)
else ! multi task
   call read_transpose_multi_task(state_ens_handle, name_handle, domain, dart_index, limit_mem)
endif

end subroutine read_transpose
!-------------------------------------------------

subroutine transpose_write(state_ens_handle, name_handle, domain, dart_index, limit_mem, write_single_precision)

type(ensemble_type),      intent(inout) :: state_ens_handle
type(restart_names_type), intent(in)    :: name_handle
integer,                  intent(in)    :: domain
integer,                  intent(inout) :: dart_index
integer,                  intent(in)    :: limit_mem !< How many state elements you can write at once
logical,                  intent(in)    :: write_single_precision

if (task_count() == 1) then
   call transpose_write_single_task(state_ens_handle, name_handle, domain, dart_index, write_single_precision)
else ! multi task
   call transpose_write_multi_task(state_ens_handle, name_handle, domain, dart_index, limit_mem, write_single_precision)
endif

end subroutine transpose_write

!-------------------------------------------------

!-------------------------------------------------
! Single processor
!-------------------------------------------------
!> Single processor version of read_transpose.  Reads ens_size whole vectors from
!> netcdf files and fills up a row of %copies for each file.
subroutine read_transpose_single_task(state_ens_handle, name_handle, domain, dart_index)

type(ensemble_type),      intent(inout) :: state_ens_handle
type(restart_names_type), intent(in)    :: name_handle
integer,                  intent(in)    :: domain
integer,                  intent(inout) :: dart_index !< This is for mulitple domains

real(r8), allocatable :: vector(:)

integer :: ncfile !< netcdf input file identifier
character(len=256) :: netcdf_filename

integer :: block_size
integer :: starting_point, ending_point
integer :: copy
integer :: start_var

starting_point = dart_index ! position in state_ens_handle%vars
block_size = 0

! need to read into a tempory array, then fill up copies
allocate(vector(get_domain_size(domain)))

COPIES: do copy = 1, state_ens_handle%my_num_copies

   start_var = 1 ! read first variable first

   ! open netcdf file
   if (query_read_copy(copy)) then
      netcdf_filename = get_input_file(name_handle, copy, domain)
      ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
      call nc_check(ret, 'read_transpose_single_task: opening', netcdf_filename)
   endif

   block_size = get_domain_size(domain)

   ending_point = starting_point + block_size -1

   if (query_read_copy(copy)) then
      call read_variables(ncfile, vector, 1, get_num_variables(domain), domain)
      ! close netcdf file
      ret = nf90_close(ncfile)
      call nc_check(ret, 'read_transpose_single_task: closing', netcdf_filename)
      state_ens_handle%copies(copy, starting_point:ending_point) = vector

   endif

enddo COPIES

! update starting point
starting_point = starting_point + block_size

dart_index = starting_point

deallocate(vector)

end subroutine read_transpose_single_task
!-------------------------------------------------

!> Single processor version of transpose write.  Takes copies array one row
!> at a time and writes copy to a netcdf file.
subroutine transpose_write_single_task(state_ens_handle, name_handle, domain, dart_index, write_single_precision)

type(ensemble_type),      intent(inout) :: state_ens_handle
type(restart_names_type), intent(in)    :: name_handle
integer,                  intent(in)    :: domain
integer,                  intent(inout) :: dart_index
logical,                  intent(in)    :: write_single_precision

integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename_out

real(r8), allocatable :: vector(:)

integer :: block_size
integer :: starting_point, ending_point
integer :: copy
integer :: start_var

type(time_type) :: dart_time
integer :: time_owner, time_owner_index

! need to read into a tempory array to fill with one copies
allocate(vector(get_domain_size(domain)))

starting_point = dart_index ! position in state_ens_handle%vars
block_size = 0

! need to read into a tempory array, then fill up copies

COPIES: do copy = 1, state_ens_handle%my_num_copies

   start_var = 1 ! read first variable first

   ! open netcdf file
   if (query_write_copy(copy)) then
      netcdf_filename_out = get_output_file(name_handle, copy, domain)

      if(file_exist(netcdf_filename_out)) then
         ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
         call nc_check(ret, 'transpose_write: opening', trim(netcdf_filename_out))
         call nc_write_global_att_clamping(ncfile_out, copy, domain)
      else ! create and open file

         !>@todo This is grabbing the time assuming the ensemble is var complete.
         !> Should we instead have all copies time in the ensemble handle?
         call get_copy_owner_index(copy, time_owner, time_owner_index)
         call get_ensemble_time(state_ens_handle, time_owner_index, dart_time)
         ncfile_out = create_and_open_state_output(name_handle, domain, copy, dart_time, write_single_precision)

      endif

   endif

   block_size = get_domain_size(domain)

   ending_point = starting_point + block_size -1

   if (query_write_copy(copy)) then

      vector = state_ens_handle%copies(copy, starting_point:ending_point)

      if (copy <= state_ens_handle%num_copies - state_ens_handle%num_extras) then
         ! actual copy, may need clamping
         call write_variables_clamp(ncfile_out, vector, 1, get_num_variables(domain), domain)
      else ! extra copy, don't clamp
         call write_variables(ncfile_out, vector, 1, get_num_variables(domain), domain)
      endif
      ! close netcdf file
      ret = nf90_close(ncfile_out)
      call nc_check(ret, 'transpose_write closing', netcdf_filename_out)
   endif

enddo COPIES

! update starting point
starting_point = starting_point + block_size

dart_index = starting_point

deallocate(vector)

end subroutine transpose_write_single_task


!-------------------------------------------------
! Mulitple tasks
!-------------------------------------------------
!> Read in variables from model restart file and transpose so that every processor
!> has all copies of a subset of state variables (fill state_ens_handle%copies)
!> Read and transpose data according to the memory limit imposed by
!> limit_mem. Note limit_mem cannot be smaller than a variable.
subroutine read_transpose_multi_task(state_ens_handle, name_handle, domain, dart_index, limit_mem)

type(ensemble_type),      intent(inout) :: state_ens_handle
type(restart_names_type), intent(in)    :: name_handle
integer,                  intent(in)    :: domain
integer,                  intent(inout) :: dart_index !< This is for mulitple domains
integer,                  intent(in)    :: limit_mem !< How many state elements you can read at once

integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: block_size !< number of state elements in a block
integer :: elm_count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: start_rank
integer :: recv_start, recv_end
integer :: send_start, send_end
integer :: ensemble_member !< the ensmeble_member you are receiving.
integer :: dummy_loop
integer :: my_copy !< which copy a pe is reading, from 1 to ens_handle%num_copies
integer :: c !< copies_read loop index
integer :: copies_read
integer :: num_state_variables
logical :: is_reader ! pe is a reader or not

integer :: ncfile !< netcdf input file identifier
character(len=256) :: netcdf_filename !< different for each task

ens_size = state_ens_handle%num_copies ! have the extras, incase you need to read inflation restarts

my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
call get_pe_loops(ens_size, recv_start, recv_end, send_start, send_end)

if (my_pe < ens_size) then
   is_reader = .true.
else
   is_reader = .false.
endif

copies_read = 0

COPIES: do c = 1, ens_size
   if (copies_read >= ens_size) exit

   ! what to do if a variable is larger than the memory limit?
   start_var = 1 ! read first variable first
   starting_point = dart_index ! position in state_ens_handle%copies

   my_copy = copies_read + my_pe + 1

   ! open netcdf file
   ! You have already opened this once to read the variable info. Should you just leave it open
   ! on the readers?
   if (is_reader) then

      if (query_read_copy(my_copy)) then
         netcdf_filename = get_input_file(name_handle, my_copy, domain)
         !print*, 'opening netcdf_filename ', trim(netcdf_filename)
         ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
         call nc_check(ret, 'read_transpose opening', netcdf_filename)
      endif

   endif

   ! Reading of the state variables is broken up into
   VARIABLE_LOOP: do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be read
      end_var = calc_end_var(start_var, domain, limit_mem)
      if ((my_task_id() == 0) .and. (c == 1)) print*, 'start_var, end_var', start_var, end_var
      block_size = get_sum_variables(start_var, end_var, domain)

      if (is_reader) then
         if (query_read_copy(my_copy)) then

            allocate(var_block(block_size))
            call read_variables(ncfile, var_block, start_var, end_var, domain)

         endif
      endif

      start_rank = get_start_rank(start_var, domain)

      ! loop through and post recieves
      RECEIVING_PE_LOOP: do recv_pe = recv_start, recv_end

         ! work out elm_count on the receiving pe
         elm_count = num_elements_on_pe(recv_pe, start_rank, block_size)
         ending_point = starting_point + elm_count -1

         ! work out the start in var_block corresponding to the receiving pe
         i = find_start_point(recv_pe, start_rank)

         if (my_pe == recv_pe) then ! get ready to recieve from each reader

            ensemble_member =  copies_read + 1

            RECEIVE_FROM_EACH: do sending_pe = send_start, send_end

               if (query_read_copy(sending_pe + copies_read + 1)) then

                  if(sending_pe == recv_pe) then ! just copy
                     state_ens_handle%copies(ensemble_member, starting_point:ending_point ) = &
                     var_block(i:elm_count*task_count():task_count())
                  else ! post receive
                     call recv_variables_from_read(state_ens_handle, sending_pe, ensemble_member, starting_point, ending_point)
                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo RECEIVE_FROM_EACH

            ! update starting point

            starting_point = starting_point + elm_count

         elseif (is_reader) then ! sending

            if (query_read_copy(my_copy)) then

               call send_variables_from_read(state_ens_handle, recv_pe, i, elm_count, block_size, var_block)

            endif

         endif

      enddo RECEIVING_PE_LOOP

      start_var = end_var + 1

      if (is_reader) then
         if (query_read_copy(my_copy)) deallocate(var_block)
      endif

   enddo VARIABLE_LOOP

   ! keep track of how many copies have been read.
   copies_read = copies_read + task_count()

   ! close netcdf file
   if (is_reader) then
      if (query_read_copy(my_copy)) then
         ret = nf90_close(ncfile)
         call nc_check(ret, 'read_transpose closing', netcdf_filename)
      endif
   endif

enddo COPIES

dart_index = starting_point

end subroutine read_transpose_multi_task

!-------------------------------------------------
!> Transpose from state_ens_handle%copies to the writers according to 
!> the memory limit imposed by limit_mem.
!> 
!> This is assuming round-robin layout of state on procesors (distribution type 1
!> in the ensemble handle).
subroutine transpose_write_multi_task(state_ens_handle, name_handle, domain, dart_index, limit_mem, write_single_precision)

type(ensemble_type),      intent(inout) :: state_ens_handle
type(restart_names_type), intent(in)    :: name_handle
integer,                  intent(in)    :: domain
integer,                  intent(inout) :: dart_index
integer,                  intent(in)    :: limit_mem !< How many state elements you can read at once
logical,                  intent(in)    :: write_single_precision

integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: block_size !< number of variables in a block
integer :: elm_count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: start_rank
integer :: recv_start, recv_end
integer :: send_start, send_end
integer :: ensemble_member
integer :: dummy_loop
integer :: my_copy !< which copy a pe is reading, starting from 1 to num_copies
integer :: c !< copies_read loop index
integer :: copies_written

integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename_out !< different for each task

integer :: num_state_variables

! single file
type(time_type) :: dart_time
integer :: time_owner, time_owner_index

logical :: is_writer

ens_size = state_ens_handle%num_copies ! have the extras incase you want to read inflation restarts
my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group. 
! Flipped send and recv compared to read_transpose.
call get_pe_loops(ens_size, send_start, send_end, recv_start, recv_end)
if (my_pe < ens_size) then  ! I am a writer
   is_writer = .true.
else
   is_writer = .false.
endif

copies_written = 0

COPIES : do c = 1, ens_size
   if (copies_written >= ens_size) exit

   start_var = 1 ! collect first variable first
   starting_point = dart_index ! position in state_ens_handle%copies

   my_copy = copies_written + my_pe + 1

   ! writers open netcdf output file. This is a copy of the input file
   if (is_writer) then
      if ( query_write_copy(my_copy)) then
         netcdf_filename_out = get_output_file(name_handle, my_copy, domain)

         if(file_exist(netcdf_filename_out)) then
            ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
            call nc_check(ret, 'transpose_write: opening', trim(netcdf_filename_out))
            call nc_write_global_att_clamping(ncfile_out, my_copy, domain)
         else ! create and open output file

            !>@todo This is grabbing the time assuming the ensemble is var complete.
            !> Should we instead have all copies time in the ensemble handle?
            call get_copy_owner_index(my_copy, time_owner, time_owner_index)
            call get_ensemble_time(state_ens_handle, time_owner_index, dart_time)

            ncfile_out = create_and_open_state_output(name_handle, domain, my_copy, dart_time, write_single_precision)
         endif
      endif

   endif

   VARIABLE_LOOP: do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be sent to writer
      end_var = calc_end_var(start_var, domain, limit_mem)
      if ((my_task_id() == 0) .and. (c == 1 )) print*, 'start_var, end_var', start_var, end_var
      block_size = get_sum_variables(start_var, end_var, domain)

      if (is_writer) then
         if (query_write_copy(my_copy)) then
            allocate(var_block(block_size))
         endif
      endif

      start_rank =  get_start_rank(start_var, domain)

      SENDING_PE_LOOP: do sending_pe = send_start, send_end

         ! work out elm_count on the sending pe
         elm_count = num_elements_on_pe(sending_pe, start_rank, block_size)
         ending_point = starting_point + elm_count -1

         ! work out the start in var_block corresponding to the sending_pe
         i = find_start_point(sending_pe, start_rank)

         if (my_pe /= sending_pe ) then ! post recieves

            if (query_write_copy(my_copy)) then
               call recv_variables_to_write(state_ens_handle, sending_pe, i, elm_count, block_size, var_block)
            endif

         else ! send to the collector

            ensemble_member =  copies_written + 1

            do recv_pe = recv_start, recv_end ! no if statement because everyone sends

               if (query_write_copy(recv_pe + copies_written + 1)) then

                  if ( recv_pe /= my_pe ) then
                     call send_variables_to_write(state_ens_handle, recv_pe, ensemble_member, starting_point, ending_point)
                  else ! if sender = receiver just copy
                     var_block(i:elm_count*task_count():task_count()) = state_ens_handle%copies(ensemble_member, starting_point:ending_point)
                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo

            ! update starting point
            starting_point = starting_point + elm_count

         endif

      enddo SENDING_PE_LOOP

      if (is_writer) then ! I am a writer
         if ( query_write_copy(my_copy)) then
            !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
            if (my_copy <= state_ens_handle%num_copies - state_ens_handle%num_extras) then ! actual copy, may need clamping
               call write_variables_clamp(ncfile_out, var_block, start_var, end_var, domain)
            else ! extra copy, don't clamp
               call write_variables(ncfile_out, var_block, start_var, end_var, domain)
            endif
            deallocate(var_block)
         endif
      endif

      start_var = end_var + 1

   enddo VARIABLE_LOOP

   ! keep track of how many copies have been written
   copies_written = copies_written + task_count()

   ! close netcdf file
   if (is_writer) then 
      if (query_write_copy(my_copy)) then
         ret = nf90_close(ncfile_out)
         call nc_check(ret, 'transpose_write', 'closing')
      endif
   endif

enddo COPIES

dart_index = starting_point

end subroutine transpose_write_multi_task

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Check a variable for out of bounds and clamp or fail if needed
!-------------------------------------------------------------------------------
function clamp_variable(dom_id, var_index, variable)

integer,     intent(in) :: dom_id ! domain id
integer,     intent(in) :: var_index ! variable index
real(r8), intent(inout) :: variable(:) ! variable

logical :: clamp_variable

real(r8) :: minclamp, maxclamp
character(len=NF90_MAX_NAME) :: varname ! for debugging only

! assume the variable will not be clamped
clamp_variable = .false.

! is lower bound set
minclamp = get_io_clamping_minval(dom_id, var_index)
if ( minclamp /= missing_r8 ) then
   !>@todo Calculating both minval, and max might be expensive for
   !> large variables.  Should we handle this in a different way?
   !> do we want to have a message if the variable is being clamped?
   if ( minval(variable) < minclamp ) then
      
      varname = get_variable_name(dom_id, var_index) 
      write(msgstring, *) 'min val = ', minval(variable), &
                          'min bounds = ', minclamp
      !@>todo FIXME : do we want to reduce and print out the absolute max/min
      call error_handler(E_ALLMSG, 'clamp_variable', &
                  'Clamping '//trim(varname)//', values out of bounds.', &
                   source,revision,revdate, text2=msgstring)

      !>@todo FIXME : if allow_missing_in_state then we need to be careful
      !>              not to clamp missing values
      variable = max(minclamp, variable)
      clamp_variable = .true.
   endif
endif ! min range set

! is upper bound set
maxclamp = get_io_clamping_maxval(dom_id, var_index)
if ( maxclamp /= missing_r8 ) then
   if ( maxval(variable) > maxclamp ) then
      varname = get_variable_name(dom_id, var_index) 
      write(msgstring, *) 'max val = ', maxval(variable), &
                          'max bounds = ', maxclamp
      call error_handler(E_ALLMSG, 'clamp_variable', &
                  'Clamping '//trim(varname)//', values out of bounds.', &
                   source,revision,revdate, text2=msgstring)

      variable = min(maxclamp, variable)
      clamp_variable = .true.
   endif
endif ! max range set

end function clamp_variable

!-------------------------------------------------------------------------------
!> Read in variables from start_var to end_var
!> FIXME: At the moment, this code is assuming that the variables in the state start
!> at (1,1,1) and that the whole variable is read. This is not the case for 
!> Tiegcm and CLM.  
!-------------------------------------------------------------------------------
subroutine read_variables(ncfile_in, var_block, start_var, end_var, domain)

integer,  intent(in)    :: ncfile_in
real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain

integer :: i
integer :: start_in_var_block, end_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: ret, var_id

start_in_var_block = 1

do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   end_in_var_block = start_in_var_block + var_size - 1

   ! number of dimensions and length of each
   allocate(dims(get_io_num_dims(domain, i)))

   dims = get_io_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_in, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'read_variables: nf90_inq_varid',trim(get_variable_name(domain,i)) )

   ret = nf90_get_var(ncfile_in, var_id, var_block(start_in_var_block:end_in_var_block), count=dims)
   call nc_check(ret, 'read_variables: nf90_get_var',trim(get_variable_name(domain,i)) )

   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine read_variables

!-------------------------------------------------------------------------------
!> Write variables from start_var to end_var no clamping
!-------------------------------------------------------------------------------
subroutine write_variables(ncfile_out, var_block, start_var, end_var, domain)

integer,  intent(in)    :: ncfile_out
real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain 

integer :: start_in_var_block, end_in_var_block
integer, allocatable :: dims(:)
integer :: i, ret, var_id, var_size

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   end_in_var_block = start_in_var_block + var_size - 1
   
   if ( do_io_update(domain, i ) ) then
      ! number of dimensions and length of each
      allocate(dims(get_io_num_dims(domain, i)))

      dims = get_io_dim_lengths(domain, i)

      ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
      call nc_check(ret, 'write_variables', 'getting variable id')

      ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:end_in_var_block), count=dims)
      call nc_check(ret, 'write_variables', 'writing')

      deallocate(dims)
   endif

   start_in_var_block = start_in_var_block + var_size

enddo

end subroutine write_variables

!-------------------------------------------------------------------------------
!> Write variables from start_var to end_var for actual ensemble members
!-------------------------------------------------------------------------------
subroutine write_variables_clamp(ncfile_out, var_block, start_var, end_var, domain)

integer,  intent(in) :: ncfile_out
real(r8), intent(inout) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var
integer,  intent(in) :: domain 

integer :: start_in_var_block, end_in_var_block
integer, allocatable :: dims(:)
integer :: i, ret, var_id, var_size
logical :: clamped

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   end_in_var_block = start_in_var_block + var_size - 1

   if ( do_io_update(domain, i ) ) then
      ! check whether you have to do anything to the variable, clamp or fail
      if ( do_io_clamping(domain, i) ) then
         clamped = clamp_variable(domain, i, var_block(start_in_var_block:end_in_var_block))
      endif

      ! number of dimensions and length of each
      allocate(dims(get_io_num_dims(domain, i)))

      dims = get_io_dim_lengths(domain, i)

      ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
      call nc_check(ret, 'write_variables_clamp', 'getting variable id')

      ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:end_in_var_block), count=dims)
      call nc_check(ret, 'write_variables_clamp', 'writing')

      deallocate(dims)
   endif

   start_in_var_block = start_in_var_block + var_size

enddo

end subroutine write_variables_clamp

!-------------------------------------------------------------------------------
!> Create the output files
!>
!> A 'blank' domain is one variable called state, with dimension = model size.
!> It is used when the model has not suppled any netdcf info but 
!>     direct_netcdf_write = .true.
!>@todo The file is not closed here.
!-------------------------------------------------------------------------------
function create_and_open_state_output(ens_name_handle, dom_id, copy_number, dart_time, single_precision_output) result(ncfile_out)

type(restart_names_type), intent(in) :: ens_name_handle
integer,                  intent(in) :: dom_id !< domain
integer,                  intent(in) :: copy_number
type(time_type),          intent(in) :: dart_time
logical,                  intent(in) :: single_precision_output
integer :: ncfile_out

integer :: ret !> netcdf return code
integer :: create_mode
integer :: i, j ! loop variables
integer :: new_dimid
integer :: new_varid
integer :: ndims
integer :: xtype ! precision for netcdf file
integer :: dimids(NF90_MAX_VAR_DIMS)

character(len=NF90_MAX_NAME) :: filename

filename = get_output_file(ens_name_handle, copy_number, dom_id)

write(msgstring,*) 'Creating output file ', trim(filename)
call error_handler(E_ALLMSG,'create_and_open_state_output:', msgstring)

! What file options do you want?
create_mode = ior(NF90_CLOBBER, NF90_64BIT_OFFSET)
ret = nf90_create(filename, create_mode, ncfile_out)
call nc_check(ret, 'create_and_open_state_output: creating', trim(filename))

! filename discription
call nc_write_file_information(ncfile_out, filename, get_file_description(ens_name_handle, copy_number, dom_id))

! revision information
call nc_write_revision_info(ncfile_out)

! clamping information
call nc_write_global_att_clamping(ncfile_out, copy_number, dom_id, from_scratch=.true.)

! define dimensions, loop around unique dimensions
do i = 1, get_io_num_unique_dims(dom_id)
   ret = nf90_def_dim(ncfile_out, get_io_unique_dim_name(dom_id, i), get_io_unique_dim_length(dom_id, i), new_dimid)
   !> @todo if we already have a unique names we can take this test out
   if(ret /= NF90_NOERR .and. ret /= NF90_ENAMEINUSE) then
      call nc_check(ret, 'create_and_open_state_output', 'defining dimensions')
   endif
enddo

! define variables
do i = 1, get_num_variables(dom_id) ! loop around state variables
   ! double or single precision?
   ndims = get_io_num_dims(dom_id, i)

   if (single_precision_output) then
      xtype = NF90_REAL
   else ! write output that is the precision of filter
      xtype = get_xtype(dom_id, i)
      if (r8 == r4 .and. xtype == NF90_DOUBLE) xtype = NF90_REAL
   endif

   ! query the dimension ids
   do j = 1, ndims
      ret = nf90_inq_dimid(ncfile_out, get_dim_name(dom_id, i, j), dimids(j))
      call nc_check(ret, 'create_and_open_state_output', 'querying dimensions')
   enddo

   !>@todo FIXME : need to include all state variables to the diagnostic files
   !>@             this will include adding some extra logic in write_variable

   ! define variable name and attributes
   ret = nf90_def_var(ncfile_out, trim(get_variable_name(dom_id, i)), &
                      xtype=xtype, dimids=dimids(1:ndims), varid=new_varid)
   call nc_check(ret, 'create_and_open_state_output', 'defining variable')
   !variable_ids(i, dom_id) = new_varid
   call set_var_id(dom_id, i, new_varid)

   call nc_write_attributes(ncfile_out, filename, new_varid, dom_id, i, copy_number)

enddo

ret = nf90_enddef(ncfile_out)
call nc_check(ret, 'create_and_open_state_output', 'end define mode')

call write_model_time(ncfile_out, dart_time)

end function create_and_open_state_output

!-------------------------------------------------
!> Write model attributes if they exist
subroutine nc_write_attributes(ncFileID, filename, ncVarID, domid, varid, copy_number)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid
integer,          intent(in) :: copy_number

if ( get_long_name(domid, varid) /= ' ' ) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'long_name',get_long_name(domid, varid)),&
                'nc_write_attributes','long_name in : '//trim(filename))
endif

if ( get_short_name(domid, varid) /= ' ' ) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'short_name',get_short_name(domid, varid)),&
                'nc_write_attributes','short_name in : '//trim(filename))
endif

! attributes that are only restart files
if ( is_ensemble_copy(copy_number) ) then
   call  nc_write_variable_att_clamping(ncFileID, filename, ncVarID, domid, varid)
endif 

! attributes plus an extra note about clamping for mean and sd files 
if ( is_mean_copy(copy_number) ) then
   call  nc_write_variable_att_clamping(ncFileID, filename, ncVarID, domid, varid)
   write(msgstring,'(A)') 'mean computed prior to clamping'
   call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,'DART_note',msgstring),&
                   'nc_write_attributes','note in : '//trim(filename))
endif

if ( is_sd_copy(copy_number) ) then
   call  nc_write_variable_att_clamping(ncFileID, filename, ncVarID, domid, varid)
   write(msgstring,'(A)') 'sd computed prior to clamping'
   call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,'DART_note',msgstring),&
                   'nc_write_attributes','note in :'//trim(filename))
endif

! attributes that are only important for inflation copies
if( is_inflation_copy(copy_number) ) then 
   call nc_check(nf90_put_att(ncFileID,ncVarID,'units','unitless'),&
                'nc_write_attributes','units in :'//trim(filename))
else if ( get_units(domid, varid) /= ' ' ) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'units',get_units(domid, varid)),&
                'nc_write_attributes','units in :'//trim(filename))
endif
     
! check to see if template file has missing value attributes
if ( get_has_missing_value(domid, varid) ) then
   select case ( get_xtype(domid, varid) )
      case ( NF90_INT )
         call nc_write_missing_value_int(ncFileID, filename, ncVarID, domid, varid)
      case ( NF90_FLOAT )
         call nc_write_missing_value_r4 (ncFileID, filename, ncVarID, domid, varid)
      case ( NF90_DOUBLE )
         call nc_write_missing_value_r8 (ncFileID, filename, ncVarID, domid, varid)
   end select
endif

!>@todo FIXME: also need to have different routines for different different types
!>             of numbers for add_offset and scale_factor.  Keeping it simple for
!>             now since they are not being used.
if (get_scale_factor(domid, varid) /= MISSING_R8) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'scale_factor',get_scale_factor(domid, varid)),&
                'nc_write_attributes','scale_factor '//trim(filename))
endif

if (get_add_offset(domid, varid) /= MISSING_R8) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'add_offset',get_add_offset(domid, varid)),&
                'nc_write_attributes','add_offset '//trim(filename))
endif

end subroutine nc_write_attributes

!-------------------------------------------------
!> Write global clamping attributes to files that already exist

!>@todo FIXME : use derived type for copy number and get long name

subroutine nc_write_file_information(ncFileID, filename, description)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: description

call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,'DART_file_information',description),&
             'nc_write_file_information','file_information'//trim(filename))

end subroutine nc_write_file_information

!-------------------------------------------------
!> Write clamping to variable attributes to files created from scratch
subroutine nc_write_variable_att_clamping(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

real(r8) :: clamp_val

if ( do_io_clamping(domid, varid) ) then
   clamp_val = get_io_clamping_maxval(domid, varid)

   ! max clamping attribute
   if ( clamp_val /= MISSING_R8 ) then
      call nc_check(nf90_put_att(ncFileID,ncVarID,'DART_clamp_max',clamp_val),&
                   'nc_write_variable_att_clamping','DART_clamp_max'//trim(filename))
   endif

   ! min clamping attribute
   clamp_val = get_io_clamping_minval(domid, varid)
   if ( clamp_val /= MISSING_R8 ) then
      call nc_check(nf90_put_att(ncFileID,ncVarID,'DART_clamp_min',clamp_val),&
                   'nc_write_variable_att_clamping','DART_clamp_min'//trim(filename))
   endif

endif


end subroutine nc_write_variable_att_clamping

!-------------------------------------------------
!> Write global clamping attributes for variables that have clamping

subroutine nc_write_global_att_clamping(ncFileID, copy, domid, from_scratch)

integer, intent(in)  :: ncFileID
integer, intent(in)  :: copy
integer, intent(in)  :: domid
logical, intent(in), optional :: from_scratch

integer  :: ivar
real(r8) :: clamp_val
character(len=NF90_MAX_NAME) :: clamp_max, clamp_min, att_name
logical :: need_netcdf_def_mode

need_netcdf_def_mode = .true.

if (present(from_scratch)) then
   if (from_scratch) need_netcdf_def_mode = .false.
endif

if (need_netcdf_def_mode) call nc_check(nf90_Redef(ncFileID),'nc_write_global_att_clamping',   'redef ')

   if ( (is_ensemble_copy(copy)) .or. &
        (is_mean_copy(copy)) ) then
      
   do ivar = 1,get_num_variables(domid)
      if ( do_io_clamping(domid, ivar) ) then
   
        write(clamp_min,*)  'NA'
        write(clamp_max,*)  'NA'
        
        clamp_val = get_io_clamping_maxval(domid, ivar)
        if ( clamp_val /= MISSING_R8 ) write(clamp_max,*)  clamp_val
        
        clamp_val = get_io_clamping_minval(domid, ivar)
        if ( clamp_val /= MISSING_R8 ) write(clamp_min,*)  clamp_val
        
        write(msgstring,'(''min_val = '',A15,'' , max val = '',A15)') trim(clamp_min), trim(clamp_max)
        write(att_name,'(2A)')  'DART_clamp_', trim(get_variable_name(domid, ivar))
        call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,att_name, msgstring), &
                   'nc_write_global_att_clamping','DART_clamping_range')
      endif
   enddo
endif

if (need_netcdf_def_mode) call nc_check(nf90_enddef(ncFileID), 'nc_write_global_att_clamping', 'end define mode')

end subroutine nc_write_global_att_clamping

!-------------------------------------------------
!> Write revision information
subroutine nc_write_revision_info(ncFileID)

integer,          intent(in) :: ncFileID

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'DART_creation_date' ,str1    ), &
           'nc_write_revision_info', 'creation put ')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'DART_source'  ,source  ), &
           'nc_write_revision_info', 'source put ')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'DART_revision',revision), &
           'nc_write_revision_info', 'revision put ')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'DART_revdate' ,revdate ), &
           'nc_write_revision_info', 'revdate put ')

end subroutine nc_write_revision_info

!-------------------------------------------------
!> Write model integer missing_value/_FillValue attributes if they exist
subroutine nc_write_missing_value_int(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

integer :: missingValINT, spvalINT

call get_missing_value(domid, varid, missingValINT)
if (missingValINT /= MISSING_I) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'missing_value',missingValINT), &
                 'nc_write_missing_value_int','missing_value '//trim(filename))
endif
call get_fillValue(domid, varid, spvalINT)
if (spvalINT /= MISSING_I) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'_FillValue',spvalINT), &
                 'nc_write_missing_value_int','_FillValue'//trim(filename))
endif

end subroutine nc_write_missing_value_int

!-------------------------------------------------
!> Write model r4 missing_value/_FillValue attributes if they exist
subroutine nc_write_missing_value_r4(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

real(r4) :: missingValR4, spvalR4

call get_missing_value(domid, varid, missingValR4)
if (missingValR4 /= MISSING_R4) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'missing_value',missingValR4), &
                 'nc_write_missing_value_r4','missing_value '//trim(filename))
endif
call get_fillValue(domid, varid, spvalR4)
if (spValR4 /= MISSING_R4) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'_FillValue',spvalR4), &
                 'nc_write_missing_value_r4','_FillValue'//trim(filename))
endif

end subroutine nc_write_missing_value_r4

!-------------------------------------------------
!> Write model r8 missing_value/_FillValue attributes if they exist
subroutine nc_write_missing_value_r8(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

real(r8) :: missingValR8, spvalR8

call get_missing_value(domid, varid, missingValR8)
if (missingValR8 /= MISSING_R8) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'missing_value',missingValR8), &
                 'nc_write_missing_value_r8','missing_value '//trim(filename))
endif
call get_fillValue(domid, varid, spvalR8)
if (spvalR8 /= MISSING_R8) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'_FillValue',spvalR8), &
                 'nc_write_missing_value_r8','_FillValue'//trim(filename))
endif

end subroutine nc_write_missing_value_r8

!-------------------------------------------------
!> Calculate how many variables to read in one go.
function calc_end_var(start_var, domain, limit_mem)

integer              :: calc_end_var !< end variable index
integer, intent(in)  :: start_var !< start variable index
integer, intent(in)  :: domain
integer, intent(in)  :: limit_mem

integer :: i, var_count
integer :: num_state_variables
integer, allocatable :: num_elements(:) !< cummulative size

num_state_variables = get_num_variables(domain)

allocate(num_elements(num_state_variables - start_var + 1))

calc_end_var = num_state_variables ! assume you can fit them all to start with

var_count = 0

do i = 1, num_state_variables - start_var + 1
   num_elements(i) = get_sum_variables(start_var, start_var + var_count, domain)
   var_count = var_count + 1
enddo

var_count = 1
do i = start_var, num_state_variables

   if (start_var == num_state_variables) then
      calc_end_var = num_state_variables
      exit
   endif

   if (var_count >= num_state_variables) then
      calc_end_var = num_state_variables
      exit
   endif

   if(var_count + 1> size(num_elements)) then
      calc_end_var = num_state_variables
      exit
   endif

   if (num_elements(var_count+1) >= limit_mem ) then
      calc_end_var =  i
      exit
   endif
   var_count = var_count + 1
enddo

deallocate(num_elements)

end function calc_end_var

!-------------------------------------------------
!> Find pes for loop indices
subroutine get_pe_loops(ens_size, recv_start, recv_end, send_start, send_end)

integer, intent(in)    :: ens_size
integer, intent(out)   :: recv_start !< for RECIEVING_PE_LOOP
integer, intent(out)   :: recv_end !< for RECIEVING_PE_LOOP
integer, intent(out)   :: send_start !< for RECEIVE_FROM_EACH_LOOP
integer, intent(out)   :: send_end !< for RECEIVE_FROM_EACH_LOOP


recv_start = 0
recv_end = task_count() -1

send_start  = recv_start

if (ens_size > task_count()) then
   send_end = send_start + task_count() -1
else
   send_end = send_start + ens_size -1
endif

end subroutine get_pe_loops

!--------------------------------------------------------
!--------------------------------------------------------
! Routines that are making the assumption that the ensemble
! distribution is round-robin (disrtibution type 1)
!------------------------------------------------------
!--------------------------------------------------------
!> Send elements of variables to a processor. This routine must be called
!> with a corresponding 'recv_variables_from_read'.
!> The data on the sender are non-contiguous with a stride of task_count. The
!> start is different depending on which pe is the recv and which variables 
!> are being sent (these are calculated in the calling routine).
!> The data on the reciever are contiguous (it is the %copies array)
subroutine send_variables_from_read(state_ens_handle, recv_pe, start, elm_count, block_size, variable_block)

type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in) :: recv_pe ! receiving pe
integer,             intent(in) :: start ! start in variable block on sender.
integer,             intent(in) :: elm_count ! how many elemets
integer,             intent(in) :: block_size ! size of info on sender - the receiver only
                                              ! gets part of this.
real(r8),            intent(in) :: variable_block(block_size) ! variable info

real(r8), allocatable :: buffer(:) ! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   ! MPI vector data type or packing should be used here.
   allocate(buffer(elm_count))
   buffer = variable_block(start:elm_count*task_count():task_count())
   call send_to(map_pe_to_task(state_ens_handle, recv_pe), &
                           buffer)
   deallocate(buffer)

else

   call error_handler(E_ERR, 'send_variables_from_read', 'distributions other than 1 not supported')

endif

end subroutine send_variables_from_read

!--------------------------------------------------------
!> Send the data from a pe to a writer.
!> Note this may be 1 variable or many.  Start is the start index in %copies
!> on the sending pe. Finish is the last index in %copies to send. 
!> If all variables are transposed at once, 
!> start = 1, 
!> finish = ens_handle%my_num_vars  (on sending pe)
!> This routine must be called with a corresponding 'recv_variables_to_write.'
!> The data on the sender is non-contiguous since it is a ROW of %copies.
subroutine send_variables_to_write(state_ens_handle, recv_pe, ensemble_member, start, finish)

type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in) :: recv_pe ! receiving pe
integer,             intent(in) :: ensemble_member
integer,             intent(in) :: start  ! start in copies array on sender.
integer,             intent(in) :: finish ! end in copies array on sender

real(r8), allocatable :: buffer(:) ! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   ! MPI vector data type or packing should be used here.
   allocate(buffer(finish - start + 1))
   buffer = state_ens_handle%copies(ensemble_member, start:finish)
   call send_to(map_pe_to_task(state_ens_handle, recv_pe), buffer)
   deallocate(buffer)

else

   call error_handler(E_ERR, 'send_variables_to_write', 'distributions other than 1 not supported')

endif

end subroutine send_variables_to_write

!--------------------------------------------------------
!> Receieve data from a reader. Start and finish are the local indicies
!> in the %copies array for the data being received.
!> If all variables are transposed at once, 
!> start = 1, 
!> finish = ens_handle%my_num_vars  (on receiveing pe)
!> This routine must be called with a corresponding 'send_variables_from_read'.
!> The data on the sender is non-contiguous since it is a ROW of %copies.
subroutine recv_variables_from_read(state_ens_handle, recv_pe, ensemble_member, start, finish)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: recv_pe ! receiving pe
integer,             intent(in)    :: ensemble_member
integer,             intent(in)    :: start  ! start in copies array on sender.
integer,             intent(in)    :: finish ! end in copies array on sender

real(r8), allocatable :: buffer(:) ! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   ! MPI vector data type or packing should be used here.
   allocate(buffer(finish - start + 1))
   call receive_from(map_pe_to_task(state_ens_handle, recv_pe), buffer)
   state_ens_handle%copies(ensemble_member, start:finish) = buffer
   deallocate(buffer)

else

   call error_handler(E_ERR, 'recv_variables_from_read', 'distributions other than 1 not supported')

endif

end subroutine recv_variables_from_read

!--------------------------------------------------------
!> Receive data to write/collect. The data is put non-contiguously into
!> variable_block.  Variable_block is the block of data writen to a 
!> netcdf file. It may be 1 or more variables.
!> start and elm_count depend on the sending pe. These are calculated in the
!> calling code. The stride is task_count() since we are assuming a round robin 
!> distribution of state vector onto processors.
!> This routine must be called with a corresponding 'send_variables_to_write'
subroutine recv_variables_to_write(state_ens_handle, sending_pe, start, elm_count, block_size, variable_block)

type(ensemble_type), intent(in)    :: state_ens_handle
integer,             intent(in)    :: sending_pe ! sending_pe
integer,             intent(in)    :: start ! start in vars array on reciever.
integer,             intent(in)    :: elm_count ! how many elemets
integer,             intent(in)    :: block_size ! size of info on sender - the receiver only
                                              ! gets part of this.
real(r8),            intent(inout) :: variable_block(block_size) ! variable info

real(r8), allocatable :: buffer(:) ! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   ! MPI vector data type or packing should be used here.
   allocate(buffer(elm_count))
   call receive_from(map_pe_to_task(state_ens_handle, sending_pe), &
                           buffer)
   variable_block(start:elm_count*task_count():task_count()) = buffer
   deallocate(buffer)

else

   call error_handler(E_ERR, 'recv_variables_to_write', 'distributions other than 1 not supported')

endif


end subroutine recv_variables_to_write

!--------------------------------------------------------
!> Calculate number of elements going to the receiving pe (read_transpose)
!> Or being sent from the sending pe (transpose_write) for a given 
!> start_rank and block_size.
!> block_size is the number of elements in a block of variables. There may
!> be 1 variable or all variables depending on limit_mem.
!> start_rank is the pe that owns the 1st element of the 1st variable in 
!> the variable_block.
!> This should go in ensemble manager. 
function num_elements_on_pe(pe, start_rank, block_size) result(elm_count)

integer, intent(in) :: pe
integer, intent(in) :: start_rank
integer, intent(in) :: block_size

integer :: elm_count, remainder

elm_count = block_size/task_count()
remainder = mod(block_size, task_count())

! mop up leftovers CHECK THESE.
! How many elements a pe gets depends on its rank relative to start_rank.
if ( (start_rank <= pe) .and. (pe) < (start_rank + remainder)) elm_count = elm_count + 1
if ( pe < (start_rank + remainder - task_count() )) elm_count = elm_count + 1

end function num_elements_on_pe

!--------------------------------------------------------
!> Give the rank of the processor that owns the start of a variable
function get_start_rank(variable, domain)

integer, intent(in) :: variable
integer, intent(in) :: domain

integer :: get_start_rank

get_start_rank = mod(get_sum_variables_below(variable, domain), task_count())

end function get_start_rank

!-------------------------------------------------------
!> Find i, the start point in var_block for a given recv_pe
!> This is assuming round robin. - will have to query the 
!> ensemble manager to find this for different disrtibutions
function find_start_point(recv_pe, start_rank)

integer, intent(in)  :: recv_pe !< the receiver
integer, intent(in)  :: start_rank !< the pe that owns the 1st element of the var_block
integer              :: find_start_point

if (start_rank < recv_pe) then
   find_start_point = recv_pe - start_rank + 1
elseif(start_rank > recv_pe) then
   find_start_point = recv_pe + task_count() - start_rank + 1
else ! recv_pe = start_rank
   find_start_point = 1
endif

end function find_start_point

!--------------------------------------------------------
!--------------------------------------------------------

!> @}
end module direct_netcdf_mod

