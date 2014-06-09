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
!>If limit 1 > state vector size and limit 2 > number of tasks, you have the regular transpose.

module state_vector_io_mod

!> \defgroup state_vector_io state_vector_io
!> Contains all the routines and variables to deal with a limited tranpose
!> You can limit the transpose by memory using <code>limit_mem</code> and by
!> processors using <code> limit_procs </code>
!> @{

use types_mod, only : r8, MISSING_R8
use mpi_utilities_mod, only : task_count, send_to, receive_from, my_task_id
use ensemble_manager_mod, only : ensemble_type, map_pe_to_task
use utilities_mod, only : error_handler, E_ERR, nc_check
use netcdf
use mpi

implicit none

private

public :: netcdf_filename, get_state_variable_info, &
          read_transpose, &
          limit_mem, limit_procs, &
          transpose_write, &
          netcdf_filename_out

integer :: ret !< netcdf return code
integer :: ncfile !< netcdf input file identifier
integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename !< needs to be different for each task
character(len=256) :: netcdf_filename_out !< needs to be different for each task
integer, parameter :: MAXDIMS = 10 ! FIXME netcdf-max_dims?
integer :: limit_mem !< This is the number of elements (not bytes) so you don't have times the number by 4 or 8
integer :: limit_procs

integer,            allocatable :: variable_ids(:)
integer,            allocatable :: variable_sizes(:)
integer,            allocatable :: variable_xtypes(:) !< variable type: NF90_FLOAT, NF90_DOUBLE, etc.
character(len=256), allocatable :: variable_names_local(:) !< local copy of variable names, not sure whether this is the way to go.
integer,   allocatable :: dimensions_and_lengths(:,:) !< number of dimensions and length of each dimension
integer :: num_state_variables

! record dimensions for output netcdf file
!integer              :: dimIds(MAXDIMS) !< dimension ids
!integer              :: length(MAXDIMS) !< dimension length
integer,            allocatable :: dimIds(:, :) !< dimension ids
integer,            allocatable :: length(:, :) !< dimension length
character(len=256), allocatable :: dim_names(:, :)

contains

!-------------------------------------------------
!> Need list of variables in the state
!> Each task grabs the variable dimensions from the netcdf file
subroutine get_state_variable_info(n, variable_names, model_size)

integer,            intent(in)  :: n !< number of state variables
character(len=256), intent(in)  :: variable_names(n)
integer,            intent(out) :: model_size

integer :: i

if(allocated(variable_ids))   call error_handler(E_ERR, 'get_state_variable_info', 'already called this routine')
if(allocated(variable_sizes)) call error_handler(E_ERR, 'get_state_variable_info', 'already called this routine')
allocate(variable_ids(n), variable_sizes(n), variable_xtypes(n), variable_names_local(n))
allocate(dimIds(n, MAXDIMS), length(n, MAXDIMS), dim_names(n, MAXDIMS))

num_state_variables = n
variable_names_local = variable_names

! open netcdf file - all restart files have the same info?
ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
call nc_check(ret, 'get_state_variable_info', 'opening')

! get all variable ids
call get_variable_ids(variable_names, variable_ids)

! get all variable sizes, only readers store dimensions?
variable_sizes = total_size(n, variable_ids)
model_size = sum(variable_sizes)
if(my_task_id() == 0) then
   do i = 1, n
      print*, i, 'variable_sizes', variable_sizes(i)
   enddo
endif

if ( any(variable_sizes>limit_mem) ) then
   print*, 'memory limit = ', limit_mem
   print*, variable_sizes
   call error_handler(E_ERR, 'get_state_variable_info', 'netcdf variables larger than memory limit')
endif

! close netcdf file
ret = nf90_close(ncfile) ! FIXME: Use NC ERROR handler
call nc_check(ret, 'get_state_variable_info', 'closing')

end subroutine get_state_variable_info

!-------------------------------------------------
!> Get netcdf variable ids
subroutine get_variable_ids(variable_names, variable_ids)

character(len=*), intent(in)    :: variable_names(:) !< netcdf variable names
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
function total_size(n, varId)

integer, intent(in)  :: n !< number of variables in the state vector
integer, intent(in)  :: varId(n) !< variable ids
integer              :: total_size(n) ! NOTE RETURN VALUE LAST

integer              :: ndims !< number of dimensions
integer              :: i, j !< loop variable

allocate(dimensions_and_lengths(n, MAXDIMS +1))
dimensions_and_lengths = -1  ! initialize to a nonsense value

do i = 1, n

   ret = nf90_inquire_variable(ncfile, varId(i), ndims=ndims, dimids=dimIds(i, :), xtype=variable_xtypes(i))
   call nc_check(ret, 'totalsize', 'inq_var')

   do j = 1, ndims
      ret = nf90_inquire_dimension(ncfile, dimIds(i, j), name=dim_names(i, j), len=length(i, j))
      call nc_check(ret, 'totalsize', 'inq_dimlen')
   enddo

   total_size(i) = product(length(i, 1:ndims))

   dimensions_and_lengths(i, 1) = ndims
   dimensions_and_lengths(i, 2:ndims+1) = length(i, 1:ndims)

enddo

end function total_size


!-------------------------------------------------
!> Read in variables from model restart file and transpose so that every processor
!> has all copies of a subset of state variables (fill state_ens_handle%copies)
!> Read and transpose data according to the memory limit imposed by
!> limit_mem AND the task limit imposed by limit_procs
!> limit_procs is used to devide the pes into groups.  Note that the
!> groups are not created using mpi_group_incl.
subroutine read_transpose(state_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle

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

ens_size = state_ens_handle%num_copies ! this is not true, this is copies + extras.  
my_pe = state_ens_handle%my_pe

! what to do if a variable is larger than the memory limit?
start_var = 1 ! read first variable first
starting_point = 1 ! position in state_ens_handle%copies

! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
call get_pe_loops(my_pe, ens_size, group_size, recv_start, recv_end, send_start, send_end)

! open netcdf file 
! You have already opened this once to read the variable info. Should you just leave it open
! on the readers?
if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader 
   ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
   call nc_check(ret, 'read_transpose', 'opening')
endif

! Reading of the state variables is broken up into
do dummy_loop = 1, num_state_variables
   if (start_var > num_state_variables) exit ! instead of using do while loop

   ! calculate how many variables will be read
   end_var = calc_end_var(start_var)
   if (my_task_id() == 0) print*, 'start_var, end_var', start_var, end_var
   block_size = sum(variable_sizes(start_var:end_var))

   if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader 
      allocate(var_block(block_size))
      call read_variables(var_block, start_var, end_var)
   endif

   start_rank = mod(sum(variable_sizes(1:start_var-1)), task_count())

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

         ensemble_member = 1

         RECEIVE_FROM_EACH: do sending_pe = send_start, send_end ! how do we know ens_size?

            if(sending_pe == recv_pe) then ! just copy
               ! The row is no longer sending_pe + 1 because it is not just 
               ! the first ens_size processors that are sending
               state_ens_handle%copies(ensemble_member, starting_point:ending_point ) = &
                  var_block(i:count*task_count():task_count())
            else ! post receive
               call receive_from(map_pe_to_task(state_ens_handle, sending_pe), &
                                 state_ens_handle%copies(ensemble_member, starting_point:ending_point))
            endif

            ensemble_member = ensemble_member + 1

         enddo RECEIVE_FROM_EACH

         ! update starting point
         starting_point = starting_point + count

      elseif ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! sending
         call send_to(map_pe_to_task(state_ens_handle, recv_pe), &
                      var_block(i:count*task_count():task_count()))
      endif

   enddo RECEIVING_PE_LOOP

   start_var = end_var + 1

   if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! reader
      deallocate(var_block)
   endif

enddo


! close netcdf file
if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader 
   ret = nf90_close(ncfile)
   call nc_check(ret, 'read_transpose','closing')
endif

end subroutine read_transpose

!-------------------------------------------------
!> Transpose from state_ens_handle%copies to the writers according to 
!> the memory limit imposed by limit_mem AND the task limit imposed by limit_procs
!> limit_procs is used to devide the pes into groups.  Note that the
!> groups are not created using mpi_group_incl.
!> 
!> Start with limited memory only ( use all tasks )
!> Using the term collectors = first stage of amalgamating 
subroutine transpose_write(state_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle

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
integer :: copy
integer :: g !< group loop index
integer :: num_groups !< number of groups the transpose is split into.  Only relevant if limit_procs < task_count()

! mpi_type variables
integer, allocatable :: array_of_blocks(:)
integer, allocatable :: array_of_displacements(:)
integer              :: num_blocks
integer              :: ierr !< mpi error (all errors are fatal so I don't bother checking this
integer              :: collector_type !< mpi derived datatype

ens_size = state_ens_handle%num_copies ! this is not true, this is copies + extras
my_pe = state_ens_handle%my_pe

start_var = 1 ! collect first variable first
starting_point = 1 ! position in state_ens_handle%copies

! post recieves
! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
call get_pe_loops(my_pe, ens_size, group_size, send_start, send_end, recv_start, recv_end) ! think I can just flip send and recv for transpose_write?

! writers open netcdf output file. This is a copy of the input file
if (my_pe < ens_size) then
   ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
   call nc_check(ret, 'transpose_write', 'opening')
endif

do while (start_var <= num_state_variables)

   print*, 'starting point', starting_point

   ! calculate how many variables will be sent to writer
   end_var = calc_end_var(start_var)
   if (my_task_id() == 0) print*, 'start_var, end_var', start_var, end_var
   num_vars = sum(variable_sizes(start_var:end_var))

   if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector
      allocate(var_block(num_vars))
   endif

   start_rank = mod(sum(variable_sizes(1:start_var-1)), task_count())

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
            call receive_from(map_pe_to_task(state_ens_handle, sending_pe), var_block(i:count*task_count():task_count()))
         endif

      else ! send to the collector

         do recv_pe = recv_start, recv_end ! no if statement because everyone sends
            copy = recv_pe + 1
            if ( recv_pe /= my_pe ) then 
               call send_to(map_pe_to_task(state_ens_handle, recv_pe), state_ens_handle%copies(copy, starting_point:ending_point))
            else ! if sender = receiver just copy
                  var_block(i:count*task_count():task_count()) = state_ens_handle%copies(copy, starting_point:ending_point)
            endif

         enddo

         ! update starting point
         starting_point = starting_point + count

      endif

   enddo SENDING_PE_LOOP

   if (limit_procs >= task_count()) then ! no need to do second stage

      if (my_pe < ens_size) then ! I am a writer
         !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
         call write_variables(var_block, start_var, end_var)
         deallocate(var_block)
      endif

   else ! Need to aggregate onto the writers (ens_size writers) Is there a better way to do this?

      if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector
      
         do g = 2, num_groups

            ! create datatype for the data being sent to the writers - same across the ensemble
            num_blocks = 2
            allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
            array_of_blocks = 1
            array_of_displacements = 2

            call mpi_type_indexed(num_blocks, array_of_blocks, array_of_displacements, mpi_real, collector_type, ierr) ! real*4 real*8
            call mpi_type_commit(collector_type, ierr)

            ! collectors -> writers loop
            call mpi_send(var_block, 1, collector_type, recv_pe, mpi_any_tag, mpi_comm_world, ierr)

            deallocate(array_of_blocks, array_of_displacements)
            call mpi_type_free(collector_type, ierr)

         call write_variables(var_block, start_var, end_var)
         deallocate(var_block)

         enddo

      endif

   endif

   start_var = end_var + 1

enddo

! close netcdf file
if (my_pe < ens_size ) then ! I am a writer
   ret = nf90_close(ncfile_out)
   call nc_check(ret, 'transpose_write', 'closing')
endif

end subroutine transpose_write

!-------------------------------------------------
!> Calculate how many variables to read in one go.
function calc_end_var(start_var)

integer :: calc_end_var !< end variable index
integer :: start_var !< start variable index

integer :: i, count
integer, allocatable :: num_elements(:) !< cummulative size

allocate(num_elements(num_state_variables - start_var + 1))

calc_end_var = num_state_variables ! assume you can fit them all to start with

count = 0

do i = 1, num_state_variables - start_var + 1
   num_elements(i) = sum(variable_sizes(start_var:start_var + count))
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
subroutine read_variables(var_block, start_var, end_var)

real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var

integer :: i
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)

start_in_var_block = 1
do i = start_var, end_var

   var_size = variable_sizes(i)

   ! number of dimensions and length of each
   allocate(dims(dimensions_and_lengths(i,1)))

   dims = dimensions_and_lengths(i, 2:dimensions_and_lengths(i,1) + 1)

   ret = nf90_get_var(ncfile, variable_ids(i), var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)

   call nc_check(ret, 'read_variables','reading')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine read_variables

!-------------------------------------------------
!> Write variables from start_var to end_var
subroutine write_variables(var_block, start_var, end_var)

real(r8), intent(in) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var

integer :: i
integer :: count_displacement
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)

start_in_var_block = 1
do i = start_var, end_var

   var_size = variable_sizes(i)
   
   ! number of dimensions and length of each
   allocate(dims(dimensions_and_lengths(i, 1))) 
   dims = dimensions_and_lengths(i, 2:dimensions_and_lengths(i,1) + 1)

   ret = nf90_put_var(ncfile_out, variable_ids(i), var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
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

!-------------------------------------------------------
!> @}
end module state_vector_io_mod
