! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> IO for the state vector. \n
!> Idea is to be generic. \n
!> Don't want to go to the filesystem twice wrf => wrf_to_dart => dart => dart_to_wrf \n
!> Get a list of variables in the state \n
!> Get info about their dimenions \n
!> Every tasks needs the dimensions of each state variable to calculate
!> what it is going to recieve in the IO transpose
!> \par Aim:
!>
!>To limit the transpose in two ways:
!>
!>1. Limit how much of the state vector you can read at once.
!>2. Limit the number of tasks involved in a transpose.
!>
!>What you (potentially) gain from this:
!>
!>* Don't have to have the whole state vector.
!>* Don't have to use a parallel IO library.
!>
!>If limit 1 > state vector size and limit 2 > number of tasks, you have the regular transpose, except you are reading directly from a netcdf file, not a dart state vector
!> file.
!>
!> * Reading with limited processors is easy, because you just have muliple readers
!> duplicating the read.
!> * Writing with limitied processors is a bit more involved because it is no longer
!> simply duplicate work.  Every processor has something it needs to contribute to 
!> the write. Thus, there is a second stage of data aggregation if limit_procs <
!> task_count.

module state_vector_io_mod

!> \defgroup state_vector_io state_vector_io
!> Contains all the routines and variables to deal with a limited tranpose
!> You can limit the transpose by memory using <code>limit_mem</code> and by
!> processors using <code> limit_procs </code>
!> @{

use types_mod,            only : r8, MISSING_R8

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id, datasize

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, single_restart_file_in, &
                                 single_restart_file_out

use utilities_mod,        only : error_handler, E_ERR, nc_check, check_namelist_read, &
                                 find_namelist_in_file, nmlfileunit, do_nml_file, do_nml_term, file_exist, &
                                 E_MSG

use assim_model_mod,      only : get_model_size, aread_state_restart, awrite_state_restart, &
                                 open_restart_read, open_restart_write, close_restart
! should you go through assim_model_mod?
!use model_mod,            only : read_file_name, write_file_name

use time_manager_mod,     only : time_type

use netcdf

use mpi

use io_filenames_mod,     only : restart_files_in, restart_files_out, io_filenames_init


implicit none

interface get_state_variable_info
   module procedure get_state_variable_info
   module procedure get_state_variable_info_lorenz96
end interface

interface turn_read_copy_on
   module procedure turn_read_copy_on_single
   module procedure turn_read_copy_on_range
end interface

interface turn_write_copy_on
   module procedure turn_write_copy_on_single
   module procedure turn_write_copy_on_range
end interface


private

public :: state_vector_io_init, &
          initialize_arrays_for_read, &
          netcdf_filename, get_state_variable_info, &
          read_transpose, &
          transpose_write, &
          netcdf_filename_out, &
          setup_read_write, &
          turn_read_copy_on, turn_write_copy_on, &
          turn_read_copies_off, turn_write_copies_off

integer :: ret !< netcdf return code
integer :: ncfile !< netcdf input file identifier
integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename !< needs to be different for each task
character(len=256) :: netcdf_filename_out !< needs to be different for each task
integer, parameter :: MAXDIMS = NF90_MAX_VAR_DIMS 

integer,            allocatable :: variable_ids(:, :)
integer,            allocatable :: variable_sizes(:, :)
integer,   allocatable :: dimensions_and_lengths(:, :, :) !< number of dimensions and length of each dimension
integer :: num_state_variables
integer :: num_domains

! record dimensions for output netcdf file,
! (state_variable, dimension, domain)
integer,            allocatable :: dimIds(:, :, :) !< dimension ids
integer,            allocatable :: copy_dimIds(:, :, :) !< dimension ids copy
integer,            allocatable :: length(:, :, :) !< dimension length
character(len=256), allocatable :: dim_names(:, :, :)

! Stores which copies to read and write
logical, allocatable :: read_copies(:), write_copies(:)

! list of variables names in the state
character(len=256), allocatable :: global_variable_names(:)

! namelist variables with default values
! Aim: to have the regular transpose as the default
integer :: limit_mem = 2147483640!< This is the number of elements (not bytes) so you don't have times the number by 4 or 8
integer :: limit_procs = 100000!< how many processors you want involved in each transpose.
logical :: create_restarts = .true. ! what if the restart files exist?
logical :: time_unlimited = .true. ! You need to keep track of the time.

namelist /  state_vector_io_nml / limit_mem, limit_procs, create_restarts, time_unlimited

contains

!-------------------------------------------------
!> Initialize model 
!> so you can read the namelist
subroutine state_vector_io_init()

integer :: iunit, io

!call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "state_vector_io_nml", iunit)
read(iunit, nml = state_vector_io_nml, iostat = io)
call check_namelist_read(iunit, io, "state_vector_io_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=state_vector_io_nml)
if (do_nml_term()) write(     *     , nml=state_vector_io_nml)

end subroutine state_vector_io_init

!-------------------------------------------------
!> Initialize arrays.  Need to know the number of domains
subroutine initialize_arrays_for_read(n, n_domains)

integer, intent(in) :: n !< number of state variables
integer, intent(in) :: n_domains !< number of domains (and therfore netcdf files) to read

if(allocated(variable_ids))   call error_handler(E_ERR, 'initialize_arrays_for_read', 'already called this routine')
if(allocated(variable_sizes)) call error_handler(E_ERR, 'initialize_arrays_for_read', 'already called this routine')

num_state_variables = n
num_domains = n_domains

allocate(variable_ids(n, num_domains), variable_sizes(n, num_domains))
allocate(dimIds(n, MAXDIMS, num_domains), length(n, MAXDIMS, num_domains), dim_names(n, MAXDIMS, num_domains))
allocate(copy_dimIds(n, MAXDIMS, num_domains))
allocate(dimensions_and_lengths(n, MAXDIMS +1, num_domains))
allocate(global_variable_names(n))

dimensions_and_lengths = -1  ! initialize to a nonsense value
dimIds = -1

end subroutine initialize_arrays_for_read

!-------------------------------------------------
!> Need list of variables in the state
!> Each task grabs the variable dimensions from the netcdf file
subroutine get_state_variable_info(n, variable_names, domain, domain_size)

integer,            intent(in)  :: n !< number of state variables
character(len=256), intent(in)  :: variable_names(n)
integer,            intent(in)  :: domain !< which domain info you are grabbing
integer,            intent(out) :: domain_size

integer :: i

! load up module storage with variable names
global_variable_names = variable_names

! open netcdf file - all restart files have the same info?
ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
call nc_check(ret, 'get_state_variable_info', 'opening')

! get all variable ids
call get_variable_ids(variable_names, domain, variable_ids(:, domain))

! get all variable sizes, only readers store dimensions?
variable_sizes(:, domain) = total_size(n, variable_ids(:, domain), domain)
domain_size = sum(variable_sizes(:, domain))
if(my_task_id() == 0) then
   do i = 1, n
      print*, i, 'variable_sizes', variable_sizes(i, domain), trim(variable_names(i))
   enddo
endif

if ( any(variable_sizes(:, domain)>limit_mem) ) then
   print*, 'memory limit = ', limit_mem
   print*, variable_sizes
   call error_handler(E_ERR, 'get_state_variable_info', 'netcdf variables larger than memory limit')
endif

! close netcdf file
ret = nf90_close(ncfile)
call nc_check(ret, 'get_state_variable_info', 'closing')

end subroutine get_state_variable_info

!-------------------------------------------------
!> tempory lorenz_96 variable info
subroutine get_state_variable_info_lorenz96(n)

integer, intent(in) :: n !< number of state variables

variable_sizes(:,:) = 1

end subroutine get_state_variable_info_lorenz96

!-------------------------------------------------
!> Get netcdf variable ids
subroutine get_variable_ids(variable_names, domain, variable_ids)

character(len=*), intent(in)    :: variable_names(:) !< netcdf variable names
integer,          intent(in)    :: domain ! where is this used.
integer,          intent(inout) :: variable_ids(:) !< netcdf variable ids

integer :: n !< number of variables in the state
integer :: i !< loop variable

n = size(variable_ids)

do i = 1, n

   ret = nf90_inq_varid(ncfile, variable_names(i), variable_ids(i))
   call nc_check(ret, 'get_variable_ids', 'inq_var_id')

enddo

end subroutine get_variable_ids

!-------------------------------------------------
!> Get state variable size, i.e. how many elements are in the variable
!> readers need to store the dimension information
!> Every task is storing the dimension information at the moment.
function total_size(n, varId, domain)

integer, intent(in)  :: n !< number of variables in the state vector
integer, intent(in)  :: varId(n) !< variable ids
integer, intent(in)  :: domain
integer              :: total_size(n) ! NOTE RETURN VALUE LAST

integer              :: ndims !< number of dimensions
integer              :: i, j !< loop variable
integer              :: xtype !< variable type: NF90_FLOAT, NF90_DOUBLE, etc. I don't think we need to store this

do i = 1, n

   ret = nf90_inquire_variable(ncfile, varId(i), ndims=ndims, dimids=dimIds(i, :, domain), xtype=xtype)
   call nc_check(ret, 'totalsize', 'inq_var')

   do j = 1, ndims
      ret = nf90_inquire_dimension(ncfile, dimIds(i, j, domain), name=dim_names(i, j, domain), len=length(i, j, domain))
      call nc_check(ret, 'totalsize', 'inq_dimlen')
   enddo

   total_size(i) = product(length(i, 1:ndims, domain))

   dimensions_and_lengths(i, 1, domain) = ndims
   dimensions_and_lengths(i, 2:ndims+1, domain) = length(i, 1:ndims, domain)

enddo

end function total_size


!-------------------------------------------------
!> Read in variables from model restart file and transpose so that every processor
!> has all copies of a subset of state variables (fill state_ens_handle%copies)
!> Read and transpose data according to the memory limit imposed by
!> limit_mem AND the task limit imposed by limit_procs
!> limit_procs is used to devide the pes into groups.  Note that the
!> groups are not created using mpi_group_incl.
!>
!> Trying to put in single file read for small models.
subroutine read_transpose(state_ens_handle, restart_in_file_name, domain, dart_index)

type(ensemble_type), intent(inout) :: state_ens_handle
character(len=129),  intent(in)    :: restart_in_file_name
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains

integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: block_size !< number of state elements in a block
integer :: count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: remainder
integer :: start_rank
integer :: group_size
integer :: recv_start, recv_end
integer :: send_start, send_end
integer :: ensemble_member !< the ensmeble_member you are receiving.
integer :: dummy_loop
integer :: my_copy !< which copy a pe is reading, starting from 0 to match pe
integer :: c !< copies_read loop index
integer :: copies_read

! single file
integer :: iunit
type(time_type) :: ens_time

ens_size = state_ens_handle%num_copies ! have the extras, incase you need to read inflation restarts

netcdf_filename = restart_in_file_name ! lorenz_96

my_pe = state_ens_handle%my_pe

copies_read = 0

COPIES: do c = 1, ens_size
   if (copies_read >= ens_size) exit

   ! what to do if a variable is larger than the memory limit?
   start_var = 1 ! read first variable first
   starting_point = dart_index ! position in state_ens_handle%copies

   if (single_restart_file_in) then
      my_copy = c -1 ! pe 0 is going to read everything
      group_size = task_count() ! should you be able to limit processors? This assumes one big group
      send_start = 0
      send_end = 0
      recv_start = 0
      recv_end = task_count() -1
   else
      ! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
      if ( task_count() >= ens_size ) then
         my_copy = my_pe
         call get_pe_loops(my_pe, ens_size, group_size, recv_start, recv_end, send_start, send_end)
      else
         my_copy = copies_read + my_pe
         call get_pe_loops(my_pe, task_count(), group_size, recv_start, recv_end, send_start, send_end)
      endif
   endif

   if (single_restart_file_in) then ! assuming not netdf at the moment

      if ( c == 1 .and. my_pe == 0) then ! open the file - do you want task or pe?
         iunit = open_restart_read(netcdf_filename)
      endif

   else

      ! open netcdf file
      ! You have already opened this once to read the variable info. Should you just leave it open
      ! on the readers?
      if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader

         if (query_read_copy(my_copy - recv_start+ 1)) then
            netcdf_filename = restart_files_in((my_copy - recv_start +1), domain)
            !print*, 'opening netcdf_filename ', trim(netcdf_filename)
            ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
            call nc_check(ret, 'read_transpose opening', netcdf_filename)
         endif

      endif

   endif

   ! Reading of the state variables is broken up into
   do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be read
      end_var = calc_end_var(start_var, domain)
      if ((my_task_id() == 0) .and. (c == 1)) print*, 'start_var, end_var', start_var, end_var
      block_size = sum(variable_sizes(start_var:end_var, domain))

      if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader
         if (query_read_copy(my_copy - recv_start + 1)) then

            allocate(var_block(block_size))

            if (single_restart_file_in) then
               call aread_state_restart(ens_time, var_block, iunit)
            else
               call read_variables(var_block, start_var, end_var, domain)
            endif

         endif
      endif

      start_rank = mod(sum_variables_below(start_var,domain), task_count())

      ! loop through and post recieves
      RECEIVING_PE_LOOP: do recv_pe = recv_start, recv_end

         ! work out count on the receiving pe
         count = block_size/task_count()
         remainder = mod(block_size, task_count())

         ! mop up leftovers CHECK THESE.
         if ( (start_rank <= recv_pe) .and. (recv_pe) < (start_rank + remainder)) count = count + 1
         if ( recv_pe < (start_rank + remainder - task_count() )) count = count + 1
         ending_point = starting_point + count -1

         ! work out i for the receiving pe
         i = find_start_point(recv_pe, start_rank)

         if (my_pe == recv_pe) then ! get ready to recieve from each reader

            ensemble_member = 1 + copies_read

            RECEIVE_FROM_EACH: do sending_pe = send_start, send_end ! how do we know ens_size?

               if (query_read_copy(sending_pe + copies_read - recv_start + 1)) then

                  if(sending_pe == recv_pe) then ! just copy
                     ! The row is no longer sending_pe + 1 because it is not just
                     ! the first ens_size processors that are sending
                     state_ens_handle%copies(ensemble_member, starting_point:ending_point ) = &
                     var_block(i:count*task_count():task_count())
                  else ! post receive
                     call receive_from(map_pe_to_task(state_ens_handle, sending_pe), &
                                    state_ens_handle%copies(ensemble_member, starting_point:ending_point))
                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo RECEIVE_FROM_EACH

            ! update starting point

            starting_point = starting_point + count

         elseif ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! sending

            if (query_read_copy(my_copy - recv_start + 1)) then
               call send_to(map_pe_to_task(state_ens_handle, recv_pe), &
                           var_block(i:count*task_count():task_count()))
            endif

         endif

      enddo RECEIVING_PE_LOOP

      start_var = end_var + 1

      !if (.not. single_restart_file_in) then ! do you want to deallocate and reallocate for single restart file in? Probably not
         if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! reader
            if (query_read_copy(my_copy - recv_start + 1)) then
               deallocate(var_block)
            endif
         endif
      !endif

   enddo

   ! keep track of how many copies have been read.
   if (single_restart_file_in) then
      copies_read = c
   else
      copies_read = copies_read + task_count()
   endif

   if (.not. single_restart_file_in) then 
      ! close netcdf file
      if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader
         if (query_read_copy(my_copy - recv_start + 1)) then
            ret = nf90_close(ncfile)
            call nc_check(ret, 'read_transpose closing', netcdf_filename)
         endif
      endif
   endif

enddo COPIES

if (single_restart_file_in) then ! assuming not netdf at the moment
   if ( my_pe == 0) then ! close the file - do you want task or pe?
      call close_restart(iunit)
      !deallocate(var_block)   ! just dellocate after the loop, see line 453
   endif
endif


dart_index = starting_point

end subroutine read_transpose

!-------------------------------------------------
!> Transpose from state_ens_handle%copies to the writers according to 
!> the memory limit imposed by limit_mem AND the task limit imposed by limit_procs
!> limit_procs is used to devide the pes into groups.  Note that the
!> groups are not created using mpi_group_incl.
!> 
!> Two stages of collection.
!> See transpose_write.pdf for explanation of a, k, y.
!> 
subroutine transpose_write(state_ens_handle, restart_out_file_name, domain, dart_index)

type(ensemble_type), intent(inout) :: state_ens_handle
character(len=129),  intent(in)    :: restart_out_file_name
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index

integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: num_vars !< number of variables in a block
integer :: count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: remainder
integer :: start_rank
integer :: group_size
integer :: recv_start, recv_end
integer :: send_start, send_end
integer :: ensemble_member
integer :: g !< group loop index
integer :: num_groups !< number of groups the transpose is split into.  Only relevant if limit_procs < task_count()
integer :: assembling_ensemble !< which ensemble the collectors are assembling
integer :: my_group !< which group my_pe is a member of
integer :: num_in_group !< how many processors are in a group
integer :: dummy_loop, j
integer :: a, k, n, owner, y, sub_block
integer :: my_copy !< which copy a pe is reading, starting from 0 to match pe
integer :: c !< copies_read loop index
integer :: copies_written

! mpi_type variables
integer, allocatable :: array_of_blocks(:)
integer, allocatable :: array_of_displacements(:)
integer              :: num_blocks
integer              :: ierr !< mpi error (all errors are fatal so I don't bother checking this
integer              :: collector_type !< mpi derived datatype
integer status(MPI_STATUS_SIZE)

character(len=256)      :: msgstring

! single file
integer :: iunit
type(time_type) :: ens_time

ens_size = state_ens_handle%num_copies ! have the extras incase you want to read inflation restarts
my_pe = state_ens_handle%my_pe

netcdf_filename_out = restart_out_file_name ! lorenz_96

copies_written = 0

COPIES : do c = 1, ens_size
   if (copies_written >= ens_size) exit

   start_var = 1 ! collect first variable first
   starting_point = dart_index ! position in state_ens_handle%copies
   a = 0 ! start at group 1 element 1


   if (single_restart_file_out) then
      my_copy = c -1! pe 0 is going to write everything
      group_size = task_count()
      recv_start = 0
      recv_end = 0
      send_start = 0
      send_end = task_count() -1
   else
      ! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
      if ( task_count() >= ens_size ) then
         my_copy = my_pe
         call get_pe_loops(my_pe, ens_size, group_size, send_start, send_end, recv_start, recv_end) ! think I can just flip send and recv for transpose_write?
      else
         my_copy = copies_written + my_pe
         call get_pe_loops(my_pe, task_count(), group_size, send_start, send_end, recv_start, recv_end)
      endif
   endif

   if (single_restart_file_out) then
      if ( c == 1 .and. my_pe == 0) then
         iunit = open_restart_write(netcdf_filename_out)
      endif
   else

      ! writers open netcdf output file. This is a copy of the input file
      if (my_pe < ens_size) then
         if ( query_write_copy(my_copy - recv_start + 1)) then
            netcdf_filename_out = restart_files_out((my_copy - recv_start +1), domain)
            if (create_restarts) then ! How do you want to do create restarts
               call create_state_output(netcdf_filename_out, domain)
            else
               if(file_exist(netcdf_filename_out)) then
                  ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
                    call nc_check(ret, 'transpose_write opening', trim(netcdf_filename_out))
               else ! create output file if it does not exist
                  write(msgstring, *) 'Creating output file ', trim(netcdf_filename_out)
                  call error_handler(E_MSG,'state_vector_io_mod:', msgstring)
                  call create_state_output(netcdf_filename_out, domain)
               endif
            endif
         endif

      endif

   endif

   do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be sent to writer
      end_var = calc_end_var(start_var, domain)
      if ((my_task_id() == 0) .and. (c == 1 )) print*, 'start_var, end_var', start_var, end_var
      num_vars = sum(variable_sizes(start_var:end_var, domain))

      if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector
         if (query_write_copy(my_copy - send_start + 1)) then
            allocate(var_block(num_vars))
         endif
      endif

      start_rank =  mod(sum_variables_below(start_var,domain), task_count())

      SENDING_PE_LOOP: do sending_pe = send_start, send_end

         ! work out count on the sending pe
         count = num_vars/task_count()
         remainder = mod(num_vars, task_count())

         ! mop up leftovers CHECK THESE.
         if ( (start_rank <= sending_pe) .and. (sending_pe) < (start_rank + remainder)) count = count + 1
         if ( sending_pe < (start_rank + remainder - task_count() )) count = count + 1
         ending_point = starting_point + count -1

         ! work out i for the sending_pe
         i = find_start_point(sending_pe, start_rank)

         if (my_pe /= sending_pe ) then ! post recieves
            if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector
               if (query_write_copy(my_copy - send_start + 1)) then
                  call receive_from(map_pe_to_task(state_ens_handle, sending_pe), var_block(i:count*task_count():task_count()))
               endif
            endif

         else ! send to the collector

            ensemble_member = 1 + copies_written

            do recv_pe = recv_start, recv_end ! no if statement because everyone sends

               if (query_write_copy(recv_pe + copies_written - send_start + 1)) then

                  if ( recv_pe /= my_pe ) then
                     call send_to(map_pe_to_task(state_ens_handle, recv_pe), state_ens_handle%copies(ensemble_member, starting_point:ending_point))
                  else ! if sender = receiver just copy

                     var_block(i:count*task_count():task_count()) = state_ens_handle%copies(ensemble_member, starting_point:ending_point)

                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo

            ! update starting point
            starting_point = starting_point + count

         endif

      enddo SENDING_PE_LOOP

      if (limit_procs >= task_count()) then ! no need to do second stage

         ! I think for now, the single file should enter this.

         if (my_pe < ens_size) then ! I am a writer
            if ( query_write_copy(my_copy + 1)) then
               !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
               if (single_restart_file_out) then
                  if (my_pe == 0) then
                     ens_time = state_ens_handle%time(c)
                     call awrite_state_restart(ens_time, var_block, iunit)
                     deallocate(var_block)
                  endif
               else
                  call write_variables(var_block, start_var, end_var, domain)
                  deallocate(var_block)
               endif
            endif
         endif

      else ! Need to aggregate onto the writers (ens_size writers) Is there a better way to do this?
           ! I don't think you enter this if task_count < ens_size because limit_procs = task_count()

         if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector

            if (query_write_copy(my_pe - send_start + 1)) then

               assembling_ensemble = my_pe - recv_start + 1
               num_groups = task_count() / limit_procs
               if (mod(task_count(), limit_procs) /= 0) then ! have to do somthing else
                  num_groups = num_groups + 1
                  if (mod(task_count(), limit_procs) < ens_size) num_groups = num_groups - 1 ! last group is big
               endif

               my_group = my_pe / limit_procs + 1
               if (my_group > num_groups) my_group = num_groups ! last group is big

               do g = 2, num_groups ! do whole ensemble at once

                  ! only group sending and first group need to be involved
                  if ((my_group == g) .or. (my_group == 1)) then

                     ! create datatype for the data being sent to the writers - same across the ensemble
                     ! need to find size of group g. Only the last group could be a different size
                     if (g < num_groups) then
                        num_in_group = limit_procs
                     else
                        num_in_group = get_group_size(task_count(), ens_size)
                     endif

                     if (a == 0) then ! group 1, element 1 starts the var_block, g cannot be the owner

                        owner = 1
                        y = limit_procs
                        sub_block = num_vars - (g-1)*limit_procs
                        num_blocks = sub_block / task_count()
                        remainder = mod(sub_block, task_count())
                        if (remainder > 0) num_blocks = num_blocks + 1
                        allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))

                        array_of_displacements(1) = (g-1)*limit_procs
                        array_of_displacements(2) = array_of_displacements(1) + task_count()
                        array_of_blocks(1) = num_in_group

                        if (remainder < num_in_group ) then
                           array_of_blocks(num_blocks) = sub_block - task_count()*(num_blocks-1)
                        else
                           array_of_blocks(num_blocks) = num_in_group
                        endif

                     else ! have to do something else

                        ! calculate which group last element of a is in. a starts at group 1 element 1
                        do j = 1, num_groups
                           if ( a <= cumulative_tasks(j, ens_size) ) then
                              owner = j
                              exit
                           endif
                        enddo

                        ! calulate k:
                        if ( owner == 1 ) then
                           k = a
                        else
                           k = a - cumulative_tasks(owner -1, ens_size)
                        endif

                        if (k == 0) then
                           owner = owner + 1
                           if (owner == num_groups + 1) owner = 1
                        endif

                        y = get_group_size(owner*limit_procs -1, ens_size) - k

                        ! find number of tasks between owner and group 1?
                        n = cumulative_tasks(num_groups, ens_size) - cumulative_tasks(owner, ens_size)

                        !if (my_pe == 0) print*, 'g = ', g, 'owner =', owner, 'n = ', n, 'k = ', k, 'y =', y, 'num_blocks', num_blocks, 'num_vars', num_vars

                        ! find number of blocks:
                        sub_block = num_vars - y - n
                        num_blocks = sub_block / task_count()
                        if ( g >= owner ) num_blocks = num_blocks + 1 ! for y and for any blocks in n

                        remainder = mod( sub_block, task_count() )

                        if (remainder >= cumulative_tasks(g, ens_size) ) then

                           num_blocks = num_blocks + 1
                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = num_in_group

                        elseif ( (cumulative_tasks(g-1, ens_size) < remainder) .and. ( remainder < cumulative_tasks(g, ens_size)) ) then

                           num_blocks = num_blocks + 1 ! ragged end
                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = remainder - cumulative_tasks(g-1, ens_size)

                        else

                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = num_in_group

                        endif

                        if ( g == owner ) then
                           array_of_displacements(1) = 0 ! zero offset for mpi_type_indexed
                           array_of_displacements(2) = task_count() - k  ! zero offset
                           array_of_blocks(1) = y
                        elseif ( g > owner) then
                           array_of_displacements(1) = y + cumulative_tasks(g-1, ens_size) - cumulative_tasks(owner, ens_size)
                           array_of_displacements(2) = array_of_displacements(1) + task_count()
                           array_of_blocks(1) = num_in_group
                        else
                           array_of_displacements(1) = y + n + cumulative_tasks(g-1, ens_size) ! y + n + offest from start of group 1
                           array_of_displacements(2) = array_of_displacements(1) + task_count()
                           array_of_blocks(1) = num_in_group
                        endif

                     endif

                     array_of_blocks(2:num_blocks - 1) = num_in_group

                     do i = 3, num_blocks
                        array_of_displacements(i) = array_of_displacements(i-1) + task_count()
                     enddo
   
                     ! check you are not going over num_vars - you can probably pull this, it was just for debugging.
                     if(my_pe == 0) then
                        if( (array_of_displacements(num_blocks) + array_of_blocks(num_blocks))  > num_vars) then
                           print*, '++++ OVER ++++', num_vars - (array_of_displacements(num_blocks) + array_of_blocks(num_blocks)), 'last block', array_of_blocks(num_blocks), 'last disp', array_of_displacements(num_blocks), 'num_vars', num_vars, 'y', y
                           print*, 'remainder', remainder, cumulative_tasks(g-1, ens_size), cumulative_tasks(g,ens_size), num_blocks
                        endif
                     endif

                     if ( datasize == mpi_real4 ) then

                        call mpi_type_indexed(num_blocks, array_of_blocks, array_of_displacements, mpi_real4, collector_type, ierr)

                     else ! double precision

                        call mpi_type_indexed(num_blocks, array_of_blocks, array_of_displacements, mpi_real8, collector_type, ierr)

                     endif

                     call mpi_type_commit(collector_type, ierr)

                     ! collectors -> writers

                     recv_pe = assembling_ensemble - 1
                     sending_pe = recv_pe + (g-1)*limit_procs
                     if (my_pe == recv_pe) then
                        call mpi_recv(var_block, 1, collector_type, map_pe_to_task(state_ens_handle,sending_pe), 0, mpi_comm_world, status, ierr)
                     elseif (my_pe == sending_pe) then
                        call mpi_send(var_block, 1, collector_type, map_pe_to_task(state_ens_handle,recv_pe), 0, mpi_comm_world, ierr)
                     endif

                     call mpi_type_free(collector_type, ierr)
                     deallocate(array_of_blocks, array_of_displacements)

                  endif

               enddo

            endif

            if (my_pe < ens_size) then ! I am a writer
               if(query_write_copy(my_copy + 1)) then
                  !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
                  call write_variables(var_block, start_var, end_var, domain)
               endif
            endif

            if(query_write_copy(my_copy - send_start + 1)) then
               deallocate(var_block) ! all collectors have var_block
            endif

         endif

      endif

      start_var = end_var + 1
      ! calculate a:
      a = mod(num_vars - (task_count() - a), task_count())

   enddo

   ! keep track of how many copies have been written
   if (single_restart_file_out) then
      copies_written = c
   else
      copies_written = copies_written + task_count()
   endif

   if (.not. single_restart_file_out) then
      ! close netcdf file
      if (my_copy < ens_size ) then ! I am a writer
         if (query_write_copy(my_copy + 1)) then
            ret = nf90_close(ncfile_out)
            call nc_check(ret, 'transpose_write', 'closing')
         endif
      endif
   endif

enddo COPIES

if (single_restart_file_out) then
   if ( my_pe == 0 ) then
      call close_restart(iunit)
      !deallocate(var_block)
   endif
endif

dart_index = starting_point

end subroutine transpose_write

!-------------------------------------------------
!> Calculate how many variables to read in one go.
function calc_end_var(start_var, domain)

integer              :: calc_end_var !< end variable index
integer, intent(in)  :: start_var !< start variable index
integer, intent(in)  :: domain

integer :: i, count
integer, allocatable :: num_elements(:) !< cummulative size

allocate(num_elements(num_state_variables - start_var + 1))

calc_end_var = num_state_variables ! assume you can fit them all to start with

count = 0

do i = 1, num_state_variables - start_var + 1
   num_elements(i) = sum(variable_sizes(start_var:start_var + count, domain))
   count = count + 1
enddo

count = 1
do i = start_var, num_state_variables

   if (start_var == num_state_variables) then
      calc_end_var = num_state_variables
      exit
   endif

   if (count >= num_state_variables) then
      calc_end_var = num_state_variables
      exit
   endif

   if(count + 1> size(num_elements)) then
      calc_end_var = num_state_variables
      exit
   endif

   if (num_elements(count+1) >= limit_mem ) then
      calc_end_var =  i
      exit
   endif
   count = count + 1
enddo

deallocate(num_elements)

end function calc_end_var

!-------------------------------------------------
!> Read in variables from start_var to end_var
subroutine read_variables(var_block, start_var, end_var, domain)

real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain

integer :: i
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id

start_in_var_block = 1

do i = start_var, end_var

   var_size = variable_sizes(i, domain)

   ! number of dimensions and length of each
   allocate(dims(dimensions_and_lengths(i,1, domain)))

   dims = dimensions_and_lengths(i, 2:dimensions_and_lengths(i,1, domain) + 1, domain)
   ret = nf90_inq_varid(ncfile, global_variable_names(i), var_id)
   call nc_check(ret, 'read_variables','inquire variable id')

   ret = nf90_get_var(ncfile, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'read_variables','reading')

   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine read_variables


!-------------------------------------------------
!> Create the output files
!> Assume same variables in each domain
!> Overwriting reading arrays with writing info.  Is this what you want to do?
!> ncfile_out is global - is this ok?
subroutine create_state_output(filename, dom)

character(len=256), intent(in) :: filename
integer,            intent(in) :: dom !< domain, not sure whether you need this?

integer :: ret !> netcdf return code
integer :: create_mode
integer :: i, j ! loop variables
integer :: new_dimid
integer :: new_varid
integer :: ndims
integer :: xtype ! precision for netcdf file
logical :: time_dimension_exists

time_dimension_exists = .false.

! What file options do you want
create_mode = NF90_CLOBBER
ret = nf90_create(filename, create_mode, ncfile_out)
call nc_check(ret, 'create_state_output', 'creating')

! make a copy of the dimIds array
copy_dimIds = dimIds

! define dimensions
do i = 1, num_state_variables ! loop around state variables
   ! check if their dimensions exist
   do j = 1, dimensions_and_lengths(i, 1, dom) ! ndims
      if (time_unlimited .and. (dim_names(i, j, dom) == 'Time')) then ! case sensitive
         ret = nf90_def_dim(ncfile_out, dim_names(i, j, dom), NF90_UNLIMITED, new_dimid) ! does this do nothing if the dimension already exists?
         time_dimension_exists = .true.
      else
         ret = nf90_def_dim(ncfile_out, dim_names(i, j, dom), dimensions_and_lengths(i, j+1, dom), new_dimid) ! does this do nothing if the dimension already exists?
      endif
      if(ret == NF90_NOERR) then ! successfully created, store this dimenion id
         where (dimIds == dimIds(i, j, dom))
            copy_dimIds = new_dimid
         end where
      endif

   enddo
enddo

if ((.not. time_dimension_exists) .and. (time_unlimited)) then ! create unlimlited dimension time
   ret = nf90_def_dim(ncfile_out, 'Time', NF90_UNLIMITED, new_dimid)
   call nc_check(ret, 'create_state_output', 'creating time as the unlimited dimension')

   call shift_dimension_arrays(new_dimid)

endif

! overwrite dimIds
dimIds = copy_dimIds

! define variables
do i = 1, num_state_variables ! loop around state variables
   ! double or single precision?
   ndims = dimensions_and_lengths(i, 1, dom)

   if (datasize == mpi_real4) then
      xtype = nf90_real
   else
      xtype = nf90_double
   endif

   ret = nf90_def_var(ncfile_out, trim(global_variable_names(i)), xtype=xtype, dimids=dimIds(i, 1:ndims, dom), varid=new_varid)
   call nc_check(ret, 'create_state_output', 'defining variable')
      variable_ids(i, dom) = new_varid

enddo

ret = nf90_enddef(ncfile_out)
call nc_check(ret, 'create_state_output', 'end define mode')

end subroutine create_state_output

!-------------------------------------------------
!> Write variables from start_var to end_var
subroutine write_variables(var_block, start_var, end_var, domain)

real(r8), intent(in) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var
integer,  intent(in) :: domain 

integer :: i
integer :: count_displacement
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id 

start_in_var_block = 1
do i = start_var, end_var

   var_size = variable_sizes(i, domain)
   
   ! number of dimensions and length of each
   allocate(dims(dimensions_and_lengths(i, 1, domain)))
   dims = dimensions_and_lengths(i, 2:dimensions_and_lengths(i,1, domain) + 1, domain)

   ret = nf90_inq_varid(ncfile_out, global_variable_names(i), var_id)
   call nc_check(ret, 'write_variables', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'write_variables', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables

!-------------------------------------------------
!> Find pes for loop indices
subroutine get_pe_loops(pe, ens_size, group_size, recv_start, recv_end, send_start, send_end)

integer, intent(in)  :: pe
integer, intent(in)  :: ens_size
integer, intent(out) :: group_size !< size of the group I am in
integer, intent(out) :: recv_start !< for RECIEVING_PE_LOOP
integer, intent(out) :: recv_end !< for RECIEVING_PE_LOOP
integer, intent(out) :: send_start !< for RECEIVE_FROM_EACH_LOOP
integer, intent(out) :: send_end !< for RECEIVE_FROM_EACH_LOOP

group_size = get_group_size(pe, ens_size)

if ( limit_procs > task_count() ) limit_procs = task_count()
if (my_task_id() == 0) print*, 'limit_procs', limit_procs, 'limit_mem', limit_mem

! limit_procs needs to be greater than ens_size

if( group_size > limit_procs ) then ! last group is large because of odd numbers
   recv_start = ((task_count() / limit_procs) - 1) * limit_procs
else
   recv_start = (pe / limit_procs ) * limit_procs
endif

recv_end = recv_start + group_size -1

send_start  = recv_start
send_end = send_start + ens_size -1

end subroutine get_pe_loops

!-------------------------------------------------
!> Find group size for processor limited transpose
!> groups are 1:n, n+1:m, m+1:l ...
function get_group_size(pe, ens_size)

integer, intent(in)  :: pe
integer, intent(in)  :: ens_size
integer              :: get_group_size

integer :: num_groups
integer :: remainder

num_groups = task_count() / limit_procs

remainder = mod(task_count(), limit_procs)

if ( remainder > 0 ) then ! the last group has a diffent size

   if ( remainder < ens_size) then ! need to join the last two groups together

      if ( (pe + 1) > limit_procs*(num_groups - 1) ) then
         get_group_size = remainder + limit_procs
      else
         get_group_size = limit_procs
      endif

   else ! last group is smaller than the rest

      if ( (pe + 1 ) > limit_procs*num_groups ) then ! pe is a member of the last group
         get_group_size = remainder
      else
         get_group_size = limit_procs
      endif

   endif

else ! all same size
   get_group_size = limit_procs
endif

end function get_group_size

!-------------------------------------------------------
!> Find i, the start point in var_block for a given recv_pe
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

!------------------------------------------------------
!> finds number of elements in the state already
function sum_variables_below(start_var, domain)

integer, intent(in) :: start_var
integer, intent(in) :: domain
integer             :: sum_variables_below

integer :: i

sum_variables_below = 0

do i = 1, domain -1
   sum_variables_below = sum(variable_sizes(:, i)) ! whole domain below
enddo

sum_variables_below = sum_variables_below + sum(variable_sizes(1:start_var-1, domain))


end function sum_variables_below

!------------------------------------------------------
!> Given a group, finds the total number of tasks from group 1
!> up to and including that group
function cumulative_tasks(group, ens_size)

integer, intent(in) :: group
integer, intent(in) :: ens_size !< for get group size
integer             :: cumulative_tasks

integer :: i

cumulative_tasks = 0

! what if you give it a group > num_groups? Or a negative group?

do i = 1, group - 1
   cumulative_tasks = cumulative_tasks + limit_procs
enddo

! just in case group is the last group
cumulative_tasks = cumulative_tasks + get_group_size(group*limit_procs -1, ens_size)

end function cumulative_tasks

!-------------------------------------------------------
!> returns true/false depending on whether you should read this copy
function query_read_copy(c)

integer, intent(in) :: c !< copy number
logical             :: query_read_copy

if (c > size(read_copies) ) then
   query_read_copy = .false.
else
   query_read_copy = read_copies(c)
endif

end function query_read_copy

!-------------------------------------------------------
!> returns true/false depending on whether you should read this copy
function query_write_copy(c)

integer, intent(in) :: c !< copy number
logical             :: query_write_copy

if ( c > size(write_copies) ) then
   query_write_copy = .false.
else
   query_write_copy = write_copies(c)
endif

end function query_write_copy

!-------------------------------------------------------
!> Make the arrays for which copies to read and write
!> Default to just the actual copies, no extras
subroutine setup_read_write(num_copies)

integer, intent(in) :: num_copies

if( .not. allocated(read_copies) ) allocate(read_copies(num_copies))
if( .not. allocated(write_copies) ) allocate(write_copies(num_copies))

read_copies(:) = .false.
write_copies(:) = .false.

end subroutine setup_read_write

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_read_copy_on_single(c)

integer, intent(in) :: c !< copy to read

read_copies(c) = .true.

end subroutine turn_read_copy_on_single

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_write_copy_on_single(c)

integer, intent(in) :: c !< copy to write

write_copies(c) = .true.

end subroutine turn_write_copy_on_single

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_read_copy_on_range(c1, c2)

integer, intent(in) :: c1 !< start copy to read
integer, intent(in) :: c2 !< end copy to read

read_copies(c1:c2) = .true.

end subroutine turn_read_copy_on_range

!-------------------------------------------------------
!> Turn on copies to read
subroutine turn_write_copy_on_range(c1, c2)

integer, intent(in) :: c1 !< start copy to write
integer, intent(in) :: c2 !< end copy to write

write_copies(c1:c2) = .true.

end subroutine turn_write_copy_on_range

!-------------------------------------------------------
!> Turn off copies to read
subroutine turn_read_copies_off(c1, c2)

integer, intent(in) :: c1 !< start copy to read
integer, intent(in) :: c2 !< end copy to read

read_copies(c1:c2) = .false.

end subroutine turn_read_copies_off

!-------------------------------------------------------
!> Turn off copies to write
subroutine turn_write_copies_off(c1, c2)

integer, intent(in) :: c1 !< start copy to write
integer, intent(in) :: c2 !< end copy to write

write_copies(c1:c2) = .false.

end subroutine turn_write_copies_off


!-------------------------------------------------------
!> Adding space for an unlimited dimension in the dimesion arrays
subroutine shift_dimension_arrays(unlimited_dimId)

integer, intent(in)  :: unlimited_dimId

integer, allocatable :: local_copy_dimIds(:, :, :)
integer, allocatable :: local_dimensions_and_lengths(:, :, :)

! make space for the copies
allocate(local_dimensions_and_lengths(num_state_variables, MAXDIMS + 1, num_domains))
allocate(local_copy_dimIds(num_state_variables, MAXDIMS, num_domains))

! copy the arrays
local_dimensions_and_lengths = dimensions_and_lengths
local_copy_dimIds = copy_dimIds

! deallocate, and reallocate
deallocate(dimensions_and_lengths, copy_dimIds, dimIds)

allocate(dimensions_and_lengths(num_state_variables, MAXDIMS + 2, num_domains))
allocate(dimIds(num_state_variables, MAXDIMS + 1, num_domains))
allocate(copy_dimIds(num_state_variables, MAXDIMS + 1, num_domains))

! fill the arrays back up
dimensions_and_lengths(:, 3:, :) = local_dimensions_and_lengths(:, 2:, :)
dimIds(:, 2:, :) = local_copy_dimIds(:, :, :)

! add a dimension
dimensions_and_lengths(:, 1, :) = local_dimensions_and_lengths(:, 1, :) + 1
dimensions_and_lengths(:, 2, :) = 1  ! unlimited dimension is length 1?
copy_dimIds(:, 1, :) = unlimited_dimId

deallocate(local_dimensions_and_lengths, local_copy_dimIds)

end subroutine shift_dimension_arrays

!-------------------------------------------------------
!> @}
end module state_vector_io_mod
