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
!> This code has multiple places where round-robin layout of state onto task is assumed.
module direct_netcdf_mod

use types_mod,            only : r8, i4, MISSING_R8

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id,&
                                 datasize

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, &
                                 get_copy_owner_index, get_ensemble_time

use time_manager_mod,     only : time_type

use utilities_mod,        only : error_handler, nc_check, &
                                 E_MSG, E_ERR, file_exist

use state_structure_mod,  only : get_num_variables, get_sum_variables,  &
                                 get_sum_variables_below 

use io_filenames_mod,     only : get_input_file, get_output_file

use copies_on_off_mod,    only : query_read_copy, query_write_copy

use direct_netcdf_common_mod,  only : write_variables, write_variables_clamp, &
                                      read_variables, create_state_output

use mpi
use netcdf

implicit none
private
public :: read_transpose, transpose_write

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Limit_mem, limit_procs are namelist items in state_vector_io_mod.

integer :: ret !< netcdf return code
integer :: ncfile !< netcdf input file identifier
integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename !< needs to be different for each task
character(len=256) :: netcdf_filename_out !< needs to be different for each task


contains
!> \defgroup state_vector_io state_vector_io
!> Contains all the routines and variables to deal with a limited tranpose
!> You can limit the transpose by memory using <code>limit_mem</code> and by
!> processors using <code> limit_procs </code>
!> @{

!-------------------------------------------------
!> Read in variables from model restart file and transpose so that every processor
!> has all copies of a subset of state variables (fill state_ens_handle%copies)
!> Read and transpose data according to the memory limit imposed by
!> limit_mem AND the task limit imposed by limit_procs
!> limit_procs is used to devide the pes into groups.  Note that the
!> groups are not created using mpi_group_incl.
!>
!> Trying to put in single file read for small models.

subroutine read_transpose(state_ens_handle, domain, dart_index, limit_mem, read_limit_procs)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains
integer,             intent(in)    :: limit_mem !< How many state elements you can read at once
integer,             intent(in)    :: read_limit_procs !< How many processors are involved in the transpose


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
integer :: remainder
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
integer :: limit_procs ! local variable that can be modified

ens_size = state_ens_handle%num_copies ! have the extras, incase you need to read inflation restarts

my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

! limit_procs
limit_procs = read_limit_procs

! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
! This call will alter limit_procs based on ens_size and task_count() if necessary.
call get_pe_loops(my_pe, ens_size, limit_procs, recv_start, recv_end, send_start, send_end)

if ((my_pe >= send_start) .and. (my_pe <= send_end)) then
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

   if ( task_count() >= ens_size ) then
      my_copy = my_pe - recv_start + 1
   else
      my_copy = copies_read + my_pe - recv_start + 1
   endif

   ! open netcdf file
   ! You have already opened this once to read the variable info. Should you just leave it open
   ! on the readers?
   if (is_reader) then

      if (query_read_copy(my_copy)) then
         netcdf_filename = get_input_file(my_copy, domain)
         !print*, 'opening netcdf_filename ', trim(netcdf_filename)
         ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
         call nc_check(ret, 'read_transpose opening', netcdf_filename)
      endif

   endif

   ! Reading of the state variables is broken up into
   do dummy_loop = 1, num_state_variables
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
         elm_count = num_elements_on_receiving_pe(recv_pe, start_rank, block_size)
         ending_point = starting_point + elm_count -1

         ! work out i for the receiving pe
         i = find_start_point(recv_pe, start_rank)

         if (my_pe == recv_pe) then ! get ready to recieve from each reader

            ensemble_member =  copies_read + 1

            RECEIVE_FROM_EACH: do sending_pe = send_start, send_end

               if (query_read_copy(sending_pe + copies_read - recv_start + 1)) then

                  if(sending_pe == recv_pe) then ! just copy
                     ! The row is no longer sending_pe + 1 because it is not just
                     ! the first ens_size processors that are sending
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
         if (query_read_copy(my_copy)) then
            deallocate(var_block)
         endif
      endif

   enddo

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
!> This is assuming round-robin layout of state on procesors (distribution type 1
!> in the ensemble handle).
subroutine transpose_write(state_ens_handle, num_extras, domain, dart_index, isprior, limit_mem, write_limit_procs, write_single_precision)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: num_extras ! non restart copies
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index
logical,             intent(in)    :: isprior
integer,             intent(in)    :: limit_mem !< How many state elements you can read at once
integer,             intent(in)    :: write_limit_procs !< How many processors are involved in the transpose
logical,             intent(in)    :: write_single_precision


integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: num_vars !< number of variables in a block
integer :: elm_count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: remainder
integer :: start_rank
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
integer :: my_copy !< which copy a pe is reading, starting from 1 to num_copies
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

integer :: num_state_variables

! single file
type(time_type) :: dart_time
integer :: time_owner, time_owner_index

integer :: create_mode
integer :: limit_procs

logical :: is_writer, is_collector

! Limit_procs may be modified depending on ens_size and task_count()
limit_procs = write_limit_procs

ens_size = state_ens_handle%num_copies ! have the extras incase you want to read inflation restarts
my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group. 
! Flipped send and recv compared to read_transpose. Limit procs is modified depending
! on ens_size and task_count
call get_pe_loops(my_pe, ens_size, limit_procs, send_start, send_end, recv_start, recv_end)

! Note that unlike the read where there is a set of readers per group, only
! the first ens_size tasks actually write.
if (my_pe < ens_size) then  ! I am a writer
   is_writer = .true.
else
   is_writer = .false.
endif

if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then
   is_collector = .true.
else
   is_collector = .false.
endif

copies_written = 0

COPIES : do c = 1, ens_size
   if (copies_written >= ens_size) exit

   start_var = 1 ! collect first variable first
   starting_point = dart_index ! position in state_ens_handle%copies
   a = 0 ! start at group 1 element 1

   if ( task_count() >= ens_size ) then
      my_copy = my_pe - recv_start + 1
   else
      my_copy = copies_written + my_pe - recv_start + 1
   endif

   ! writers open netcdf output file. This is a copy of the input file
   if (is_writer) then
      if ( query_write_copy(my_copy)) then
         netcdf_filename_out = get_output_file((my_copy), domain, isprior)

         if(file_exist(netcdf_filename_out)) then
            ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
            call nc_check(ret, 'transpose_write: opening', trim(netcdf_filename_out))
         else ! create output file if it does not exist

            write(msgstring, *) 'Creating output file ', trim(netcdf_filename_out)
            call error_handler(E_MSG,'state_vector_io_mod:', msgstring)

            ! What file options do you want
            create_mode = ior(NF90_CLOBBER, NF90_64BIT_OFFSET)
            ret = nf90_create(netcdf_filename_out, create_mode, ncfile_out)
            call nc_check(ret, 'transpose_write: creating', trim(netcdf_filename_out))

            !>@todo This is grabbing the time assuming the ensemble is var complete.
            !> Should we instead have all copies time in the ensemble handle?
            call get_copy_owner_index(my_copy, time_owner, time_owner_index)
            call get_ensemble_time(state_ens_handle, time_owner_index, dart_time)

            call create_state_output(ncfile_out, domain, dart_time, write_single_precision)
         endif
      endif

   endif

   do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be sent to writer
      end_var = calc_end_var(start_var, domain, limit_mem)
      if ((my_task_id() == 0) .and. (c == 1 )) print*, 'start_var, end_var', start_var, end_var
      num_vars = get_sum_variables(start_var, end_var, domain)

      if (is_collector) then ! I am a collector
         if (query_write_copy(my_copy)) then
            allocate(var_block(num_vars))
         endif
      endif

      start_rank =  get_start_rank(start_var, domain)

      SENDING_PE_LOOP: do sending_pe = send_start, send_end

         ! work out elm_count on the sending pe
         elm_count = num_vars/task_count()
         remainder = mod(num_vars, task_count())

         ! mop up leftovers CHECK THESE.
         if ( (start_rank <= sending_pe) .and. (sending_pe) < (start_rank + remainder)) elm_count = elm_count + 1
         if ( sending_pe < (start_rank + remainder - task_count() )) elm_count = elm_count + 1
         ending_point = starting_point + elm_count -1

         ! work out i for the sending_pe
         i = find_start_point(sending_pe, start_rank)

         if (my_pe /= sending_pe ) then ! post recieves
            if (is_collector) then ! I am a collector
               if (query_write_copy(my_copy)) then
                  call recv_variables_to_write(state_ens_handle, sending_pe, i, elm_count, num_vars, var_block)
               endif
            endif

         else ! send to the collector

            ensemble_member =  copies_written + 1

            do recv_pe = recv_start, recv_end ! no if statement because everyone sends

               if (query_write_copy(recv_pe + copies_written - send_start + 1)) then

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

      if (limit_procs >= task_count()) then ! no need to do second stage

         ! I think for now, the single file should enter this.

         if (is_writer) then ! I am a writer
            if ( query_write_copy(my_copy)) then
               !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
               if (my_copy <= state_ens_handle%num_copies - num_extras) then ! actual copy, may need clamping
                  call write_variables_clamp(ncfile_out, var_block, start_var, end_var, domain)
               else ! extra copy, don't clamp
                  call write_variables(ncfile_out, var_block, start_var, end_var, domain)
               endif
               deallocate(var_block)
            endif
         endif

      else ! Need to aggregate onto the writers (ens_size writers) Is there a better way to do this?
           ! I don't think you enter this if task_count < ens_size because limit_procs = task_count()

         if (is_collector) then ! I am a collector

            if (query_write_copy(my_copy)) then

               assembling_ensemble = my_copy
               num_groups = task_count() / limit_procs
               if (mod(task_count(), limit_procs) /= 0) then ! have to do somthing else
                  num_groups = num_groups + 1
                  if (mod(task_count(), limit_procs) < ens_size) num_groups = num_groups - 1 ! last group is big
               endif

               my_group = my_pe / limit_procs + 1
               if (my_group > num_groups) my_group = num_groups ! last group is big

               ! Note this loop is around num_groups
               ! Group 1 contains the actual writers (they already have group 1 data)
               ! so loop from 2 to num_groups to collect data from remaining groups.
               do g = 2, num_groups ! do whole ensemble at once

                  ! only group sending and first group need to be involved
                  if ((my_group == g) .or. (my_group == 1)) then

                     ! create datatype for the data being sent to the writers - same across the ensemble
                     ! need to find size of group g. Only the last group could be a different size
                     if (g < num_groups) then
                        num_in_group = limit_procs
                     else
                        num_in_group = get_group_size(task_count(), ens_size, limit_procs)
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
                           if ( a <= cumulative_tasks(j, ens_size, limit_procs) ) then
                              owner = j
                              exit
                           endif
                        enddo

                        ! calulate k:
                        if ( owner == 1 ) then
                           k = a
                        else
                           k = a - cumulative_tasks(owner -1, ens_size, limit_procs)
                        endif

                        if (k == 0) then
                           owner = owner + 1
                           if (owner == num_groups + 1) owner = 1
                        endif

                        y = get_group_size(owner*limit_procs -1, ens_size, limit_procs) - k

                        ! find number of tasks between owner and group 1?
                        n = cumulative_tasks(num_groups, ens_size, limit_procs) - cumulative_tasks(owner, ens_size, limit_procs)

                        !if (my_pe == 0) print*, 'g = ', g, 'owner =', owner, 'n = ', n, 'k = ', k, 'y =', y, 'num_blocks', num_blocks, 'num_vars', num_vars

                        ! find number of blocks:
                        sub_block = num_vars - y - n
                        num_blocks = sub_block / task_count()
                        if ( g >= owner ) num_blocks = num_blocks + 1 ! for y and for any blocks in n

                        remainder = mod( sub_block, task_count() )

                        if (remainder >= cumulative_tasks(g, ens_size, limit_procs) ) then

                           num_blocks = num_blocks + 1
                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = num_in_group

                        elseif ( (cumulative_tasks(g-1, ens_size, limit_procs) < remainder) .and. ( remainder < cumulative_tasks(g, ens_size, limit_procs)) ) then

                           num_blocks = num_blocks + 1 ! ragged end
                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = remainder - cumulative_tasks(g-1, ens_size, limit_procs)

                        else

                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = num_in_group

                        endif

                        if ( g == owner ) then
                           array_of_displacements(1) = 0 ! zero offset for mpi_type_indexed
                           array_of_displacements(2) = task_count() - k  ! zero offset
                           array_of_blocks(1) = y
                        elseif ( g > owner) then
                           array_of_displacements(1) = y + cumulative_tasks(g-1, ens_size, limit_procs) - cumulative_tasks(owner, ens_size, limit_procs)
                           array_of_displacements(2) = array_of_displacements(1) + task_count()
                           array_of_blocks(1) = num_in_group
                        else
                           array_of_displacements(1) = y + n + cumulative_tasks(g-1, ens_size, limit_procs) ! y + n + offest from start of group 1
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
                           print*, 'remainder', remainder, cumulative_tasks(g-1, ens_size, limit_procs), cumulative_tasks(g,ens_size, limit_procs), num_blocks
                        endif
                     endif

                     call mpi_type_indexed(num_blocks, array_of_blocks, array_of_displacements, datasize, collector_type, ierr)

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

            if (is_writer) then ! I am a writer
               if(query_write_copy(my_copy)) then
                  !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
                  if (my_copy <= state_ens_handle%num_copies - num_extras) then ! actual copy, may need clamping
                     call write_variables_clamp(ncfile_out, var_block, start_var, end_var, domain)
                  else ! extra copy, don't clamp
                     call write_variables(ncfile_out, var_block, start_var, end_var, domain)
                  endif
               endif
            endif

            if(query_write_copy(my_copy)) then
               deallocate(var_block) ! all collectors have var_block
            endif

         endif

      endif

      start_var = end_var + 1
      ! calculate a:
      a = mod(num_vars - (task_count() - a), task_count())

   enddo

   ! keep track of how many copies have been written
   copies_written = copies_written + task_count()

   ! close netcdf file
   if (is_writer) then ! I am a writer
      if (query_write_copy(my_copy)) then
         ret = nf90_close(ncfile_out)
         call nc_check(ret, 'transpose_write', 'closing')
      endif
   endif

enddo COPIES

dart_index = starting_point

end subroutine transpose_write

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
subroutine get_pe_loops(pe, ens_size, limit_procs, recv_start, recv_end, send_start, send_end)

integer, intent(in)    :: pe
integer, intent(in)    :: ens_size
integer, intent(inout) :: limit_procs
integer, intent(out)   :: recv_start !< for RECIEVING_PE_LOOP
integer, intent(out)   :: recv_end !< for RECIEVING_PE_LOOP
integer, intent(out)   :: send_start !< for RECEIVE_FROM_EACH_LOOP
integer, intent(out)   :: send_end !< for RECEIVE_FROM_EACH_LOOP


integer :: group_size
! This call will adjust limit_procs:
! If greater than task_count(), limit_procs is set to task_count()
! If less than ens_size, limit_procs
! is set to ens_size because only copy owners read.
group_size = get_group_size(pe, ens_size, limit_procs)

if( group_size > limit_procs ) then ! last group is large because of odd numbers
   recv_start = ((task_count() / limit_procs) - 1) * limit_procs
else
   recv_start = (pe / limit_procs ) * limit_procs
endif

recv_end = recv_start + group_size -1

send_start  = recv_start
if (ens_size > task_count()) then
   send_end = send_start + task_count() -1
else
   send_end = send_start + ens_size -1
endif

end subroutine get_pe_loops

!-------------------------------------------------
!> Find group size for processor limited transpose
!> groups are 1:n, n+1:m, m+1:l ...
function get_group_size(pe, ens_size, limit_procs)

integer, intent(in)  :: pe
integer, intent(in)  :: ens_size ! ensemble size
integer, intent(inout)  :: limit_procs
integer              :: get_group_size

integer :: num_groups
integer :: remainder

! If less than ens_size, limit_procs is set to ens_size because only copy owners read.
if (limit_procs < ens_size) limit_procs = ens_size
! If greater than task_count(), limit_procs is set to task_count()
if (limit_procs > task_count()) limit_procs = task_count()

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

!------------------------------------------------------
!> Given a group, finds the total number of tasks from group 1
!> up to and including that group
function cumulative_tasks(group, ens_size, limit_procs)

integer, intent(in) :: group
integer, intent(in) :: ens_size
integer, intent(inout) :: limit_procs
integer             :: cumulative_tasks

integer :: i

cumulative_tasks = 0

! what if you give it a group > num_groups? Or a negative group?

do i = 1, group - 1
   cumulative_tasks = cumulative_tasks + limit_procs
enddo

! just in case group is the last group
cumulative_tasks = cumulative_tasks + get_group_size(group*limit_procs -1, ens_size, limit_procs)

end function cumulative_tasks

!--------------------------------------------------------
!--------------------------------------------------------
! Routines that are making the assumption that the ensemble
! distribution is round-robin (disrtibution type 1)
!------------------------------------------------------
!--------------------------------------------------------
!> Send elements of variables to correct processor
subroutine send_variables_from_read(state_ens_handle, recv_pe, start, elm_count, block_size, variable_block)

type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in) :: recv_pe ! receiving pe
integer,             intent(in) :: start ! start in copies array on reciever.
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

   call error_handler(E_ERR, 'send_variables_read', 'disrtibutions other than 1 not supported')

endif

end subroutine send_variables_from_read

!--------------------------------------------------------
!> Send to a writer/collector the data from a pe.
!> Note this may be 1 variable or many.  Start is the start index in %copies
!> on the sending pe. Finish is the last index in %copies to send. 
!> If all variables are transposed at once, 
!> start = 1, 
!> finish = ens_handle%my_num_vars  (on sending pe)
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

   call error_handler(E_ERR, 'send_variables_write', 'disrtibutions other than 1 not supported')

endif

end subroutine send_variables_to_write

!--------------------------------------------------------
!> Receieve data from a reader. Start and finish are the local indecies
!> in the %copies array for the data being received.
!> If all variables are transposed at once, 
!> start = 1, 
!> finish = ens_handle%my_num_vars  (on receiveing pe)
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
   buffer = state_ens_handle%copies(ensemble_member, start:finish)
   call receive_from(map_pe_to_task(state_ens_handle, recv_pe), buffer)
   state_ens_handle%copies(ensemble_member, start:finish) = buffer
   deallocate(buffer)

else

   call error_handler(E_ERR, 'send_variables_write', 'disrtibutions other than 1 not supported')

endif

end subroutine recv_variables_from_read

!--------------------------------------------------------
!> Receive data to write/collect. The data is put non-contiguously into
!> variable_block.  Variable_block is the block of data writen to a 
!> netcdf file. It may be 1 or more variables.
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

   call error_handler(E_ERR, 'send_variables_read', 'disrtibutions other than 1 not supported')

endif


end subroutine recv_variables_to_write

!--------------------------------------------------------
!> Calculate number of elements going to the receiving pe
!> This should go in ensemble manager. 
function num_elements_on_receiving_pe(recv_pe, start_rank, block_size) result(count)

integer, intent(in) :: recv_pe
integer, intent(in) :: start_rank
integer, intent(in) :: block_size

integer :: count, remainder

count = block_size/task_count()
remainder = mod(block_size, task_count())

! mop up leftovers CHECK THESE.
if ( (start_rank <= recv_pe) .and. (recv_pe) < (start_rank + remainder)) count = count + 1
if ( recv_pe < (start_rank + remainder - task_count() )) count = count + 1

end function num_elements_on_receiving_pe

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

