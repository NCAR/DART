! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program ftest_onesided

! MPI Fortran program that tests one-sided MPI communication.
! this is used in DART for tasks to retrieve data from other
! tasks without having to interrupt them or make them synchronize
! with the requesting task.

use mpi

implicit none

! if your MPI installation does not include an mpi module,
! comment out the 'use mpi' line above and comment in this
! include file.  if you have a choice, the module is better.
!include "mpif.h"

! data kinds
integer, parameter :: r8 = SELECTED_REAL_KIND(12) 
integer, parameter :: datasize = MPI_REAL8

! source, destination arrays, byte sizes, item counts, etc
integer :: our_window 
integer :: window_items = 1000
integer(MPI_OFFSET_KIND) :: window_size
integer(MPI_OFFSET_KIND) :: window_offset

real(r8), allocatable :: contiguous_array(:)
real(r8), allocatable :: target_array(:) 
integer :: target_count = 10

integer :: owner_task  
integer :: ierror
integer :: myrank
integer :: total_tasks

! set this to true to get more output.
! beware - it will print lines from each MPI task
! so it could be very verbose for large task counts.
logical :: verbose = .false.


! start of executable code
call setup_mpi()



! each task makes an MPI window and fills it with
! unique numbers. 

call create_window()

call fill_window()


! allocate and fill a local buffer to transfer
! data into.  make sure the buffer has the original
! contents before the transfer.

allocate(target_array(target_count))

call prefill_target(target_count)

call check_data_before(target_count)


! first test - have each task get 10 items from 
! task 0 at offset 100 and put them into the local
! buffer. the transferred data always goes into 
! offset 0 in the local buffer.
! do the transfer and check results

owner_task = 0
window_offset = 100

call get_data(owner_task, window_offset, target_count, target_array)

call check_data_after(owner_task, window_offset, target_count)



! second test - have each task get data from task N-1,
! wrap around for task 0.  also set different offsets
! in the source window for each task.  reset local buffer
! data first.  do the transfer and check results

call prefill_target(target_count)

call check_data_before(target_count)

owner_task = myrank - 1
if (owner_task < 0) owner_task = total_tasks - 1
window_offset = myrank + 1

call get_data(owner_task, window_offset, target_count, target_array)

call check_data_after(owner_task, window_offset, target_count)


! release resources

deallocate(target_array)

call free_window()

call takedown_mpi()

! end of main program


contains

!-------------------------------------------------------------
! allocate an array on each task and make it available to
! other tasks by creating a window with it.

subroutine create_window()

   integer :: bytesize  ! size in bytes of each item in the window

   call MPI_Type_Size(datasize, bytesize, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Type_Size() did not succeed, error code = ", ierror
      stop
   endif
   
   allocate(contiguous_array(window_items))
   contiguous_array = myrank
   
   window_size = window_items * bytesize
   
   ! Expose local memory to RMA operation by other processes in a communicator.
   call MPI_Win_Create(contiguous_array, window_size, bytesize, MPI_INFO_NULL, &
                       MPI_COMM_WORLD, our_window, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Win_Create() did not succeed, error code = ", ierror
      stop
   endif

end subroutine create_window

!-------------------------------------------------------------
! release the window, and deallocate window array

subroutine free_window()

   call MPI_Win_Free(our_window, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Win_Free() did not succeed, error code = ", ierror
      stop
   endif

   deallocate(contiguous_array)

end subroutine free_window

!-------------------------------------------------------------
! each task fills its local window with unique values of the
! form 1,00x,00y where x is the task number and y is the offset
! into the array.  

subroutine fill_window()

   integer :: i, base
   
   base = 1000000 + (myrank * 1000)
   do i=1, window_items
      contiguous_array(i) = base + i
   enddo

end subroutine fill_window

!-------------------------------------------------------------
! this is the MPI one-sided communication call.  each task
! gets wcount items, starting at an offset windex into the
! remote task's window, and puts the data into local array x.

subroutine get_data(owner, windex, wcount, x)

   integer,  intent(in)  :: owner    ! task that owns the window to copy from
   integer(MPI_OFFSET_KIND), intent(in)  :: windex   ! index into the owning task's window
   integer,  intent(in)  :: wcount   ! how many items to transfer
   real(r8), intent(out) :: x(:)     ! result
   
   ! Note to programmer: The data transfer is not guaranteed
   ! to have occured until the call to MPI_Win_Unlock.
   ! => Don't do anything with target inbetween MPI_Get and MPI_Win_Lock
   
   call MPI_Win_Lock(MPI_LOCK_SHARED, owner, 0, our_window, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Win_Lock() did not succeed, error code = ", ierror
      stop
   endif
   
   call MPI_Get(x, wcount, datasize, owner, windex-1, wcount, datasize, our_window, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Get() did not succeed, error code = ", ierror
      stop
   endif
   
   call MPI_Win_Unlock(owner, our_window, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Win_Unlock() did not succeed, error code = ", ierror
      stop
   endif

   if (verbose) print *, 'MPI_Get() returned on task ', myrank

end subroutine get_data

!-------------------------------------------------------------
! fill the target array on each task.  when
! data is transferred it should overwrite these values.
! target_array is a global.

subroutine prefill_target(item_count)

   integer, intent(in) :: item_count

   integer :: i

   do i=1, item_count
      target_array(i) = -i
   enddo

end subroutine prefill_target

!-------------------------------------------------------------
! make sure the local array has the initial values before
! starting the data transfer
! target_array is a global.

subroutine check_data_before(item_count)

   integer, intent(in) :: item_count

   integer :: i

   do i=1, item_count
      if (target_array(i) /= -i) then
         print *, "target array not initialized consistently on task ", myrank
         print *, "index ", i, " is ", target_array(i), ", expected ", -i
         stop
      endif
   enddo

end subroutine check_data_before

!-------------------------------------------------------------
! after the one-sided communication, verify that the data
! in the local array has changed.  a more exhaustive test
! would ensure only the right data changed and the rest 
! remained fixed.  target_array is a global.

subroutine check_data_after(source_task, offset_into_window, transfer_count)

   integer, intent(in) :: source_task
   integer(MPI_OFFSET_KIND), intent(in) :: offset_into_window
   integer, intent(in) :: transfer_count

   integer :: i, base, expected_val

   base = 1000000 + (source_task * 1000)

   do i=1, transfer_count
      expected_val = base + offset_into_window + i - 1
      if (target_array(i) /= expected_val) then
         print *, "one sided transfer results not as expected on task ", myrank
         print *, "index ", i, " is ", target_array(i), ", expected ", expected_val
         stop
      endif
   enddo

   if (verbose) print *, 'data check succeeded on task ', myrank

end subroutine check_data_after

!-------------------------------------------------------------
! standard mpi initialization

subroutine setup_mpi()

   ierror = -999
   call MPI_Init(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Init() did not succeed, error code = ", ierror
      print *, "If error code is -999, the most likely problem is that"
      print *, "the right MPI libraries were not found at compile time."
      stop
   endif
   
   myrank = -1
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Comm_rank() did not succeed, error code = ", ierror
      stop
   endif
   
   total_tasks = -1
   call MPI_Comm_size(MPI_COMM_WORLD, total_tasks, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Comm_size() did not succeed, error code = ", ierror
      stop
   endif

   if (verbose) print *, 'MPI initialized ok on task ', myrank

end subroutine setup_mpi

!-------------------------------------------------------------
! shut down mpi.  wait here for all tasks to finish.

subroutine takedown_mpi()

   if (verbose) print *, "Ready to finalize MPI and end program on task ", myrank
   
   call MPI_Barrier(MPI_COMM_WORLD, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Barrier() did not succeed, error code = ", ierror
      stop
   endif
   
   if (myrank == 0) print *, 'Test succeeded'
   
   call MPI_Finalize(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Finalize() did not succeed, error code = ", ierror
      stop
   endif
   
   ! shouldn't print after finalize called; hangs on some MPI implementations

end subroutine takedown_mpi

!-------------------------------------------------------------
! 
! the MPI_Get() terminology can be confusing.  
! "origin" is the local buffer and offset where the data will be put, so it
! is the destination.
!
! "target" refers to the remote task rank and offset into the source buffer.
!
! args to MPI_Get():
! origin_addr
!  Initial address of origin buffer (choice).
! origin_count
!  Number of entries in origin buffer (nonnegative integer).
! origin_datatype
!  Data type of each entry in origin buffer (handle).
! target_rank
!  Rank of target (nonnegative integer).
! target_disp
!  Displacement from window start to the beginning of the target buffer (nonnegative integer).
! target_count
!  Number of entries in target buffer (nonnegative integer).
! target datatype
!  datatype of each entry in target buffer (handle)
! win
!  window object used for communication (handle)
! ierror
!  Error status (integer).
! 
!-------------------------------------------------------------


end program ftest_onesided

