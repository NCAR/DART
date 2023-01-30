! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Contains the window information for the state.  Two windows:
!> One for all copies, one for the mean.
!> Not sure whether we should just have one window to avoid multiple synchronizations.

module window_mod

!> \defgroup window window_mod
!> @{
use mpi_utilities_mod,    only : datasize, my_task_id
use types_mod,            only : r8
use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index, &
                                 copies_in_window, init_ensemble_manager, &
                                 get_allow_transpose, end_ensemble_manager, &
                                 set_num_extra_copies, all_copies_to_all_vars, &
                                 all_vars_to_all_copies

use mpi

implicit none

private
public :: create_mean_window, create_state_window, free_mean_window, &
          free_state_window, data_count, mean_win, state_win, current_win, &
          mean_ens_handle, NO_WINDOW, MEAN_WINDOW, STATE_WINDOW

! mpi window handles
integer :: state_win   !! window for the forward operator
integer :: mean_win    !! window for the mean
integer :: current_win !! keep track of current window, start out assuming an invalid window
!>@todo the number of copies in the window is sloppy. You need to make this better.

! parameters for keeping track of which window is open
integer, parameter :: NO_WINDOW    = -1
integer, parameter :: MEAN_WINDOW  = 0 
integer, parameter :: STATE_WINDOW = 2 

integer :: data_count !! number of copies in the window
integer(KIND=MPI_ADDRESS_KIND) window_size
logical :: use_distributed_mean = .false. ! initialize to false

real(r8) :: duplicate_state(*)  ! duplicate array for cray pointer fwd
pointer(a, duplicate_state)

real(r8) :: duplicate_mean(*)  ! duplicate array for cray pointer vert convert
pointer(b, duplicate_mean)
type(ensemble_type) :: mean_ens_handle

contains

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_state_window(state_ens_handle, fwd_op_ens_handle, qc_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle !< state ensemble handle
type(ensemble_type), intent(inout), optional :: fwd_op_ens_handle
type(ensemble_type), intent(inout), optional :: qc_ens_handle

integer :: ii, jj, count, ierr
integer :: bytesize !< size in bytes of each element in the window
integer :: my_num_vars

! Find out how many copies to put in the window
! copies_in_window is not necessarily equal to ens_handle%num_copies
data_count = copies_in_window(state_ens_handle)

if (get_allow_transpose(state_ens_handle)) then
   call all_copies_to_all_vars(state_ens_handle)
   if (present(fwd_op_ens_handle)) then
      call all_copies_to_all_vars(fwd_op_ens_handle)
   endif
   if (present(qc_ens_handle)) then
      call all_copies_to_all_vars(qc_ens_handle)
   endif
else
   ! find how many variables I have
   my_num_vars = state_ens_handle%my_num_vars

   ! allocate some RDMA accessible memory
   ! using MPI_ALLOC_MEM because the MPI standard allows vendors to require MPI_ALLOC_MEM for remote memory access
   call mpi_type_size(datasize, bytesize, ierr)
   window_size = my_num_vars*data_count*bytesize
   a = malloc(my_num_vars*data_count)
   call MPI_ALLOC_MEM(window_size, MPI_INFO_NULL, a, ierr)

   count = 1
   ! create a duplicate copies array for remote memory access
   ! Doing this because you cannot use a cray pointer with an allocatable array
   ! Can't do array assignment with a cray pointer, so you need to loop
   do ii = 1, my_num_vars
      do jj = 1, data_count
         duplicate_state(count) = state_ens_handle%copies(jj,ii)
         count = count + 1
      enddo
   enddo

   ! Expose local memory to RMA operation by other processes in a communicator.
   call mpi_win_create(duplicate_state, window_size, bytesize, MPI_INFO_NULL, mpi_comm_world, state_win, ierr)

endif

! Set the current window to the state window
current_win = STATE_WINDOW

data_count = copies_in_window(state_ens_handle)

end subroutine create_state_window

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_mean_window(state_ens_handle, mean_copy, distribute_mean)

type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in)  :: mean_copy
logical,             intent(in)  :: distribute_mean

integer :: ii, ierr
integer :: bytesize
integer :: my_num_vars

! find out how many variables I have
my_num_vars = state_ens_handle%my_num_vars

! create an ensemble handle of just the mean copy.
use_distributed_mean = distribute_mean

if (use_distributed_mean) then

   call init_ensemble_manager(mean_ens_handle, 1, state_ens_handle%num_vars) ! distributed ensemble
   call set_num_extra_copies(mean_ens_handle, 0)
   mean_ens_handle%copies(1,:) = state_ens_handle%copies(mean_copy, :)

   ! find out how many variables I have
   my_num_vars = mean_ens_handle%my_num_vars
   call mpi_type_size(datasize, bytesize, ierr)
   window_size = my_num_vars*bytesize
   ! allocate some RDMA accessible memory
   ! using MPI_ALLOC_MEM because the MPI standard allows vendors to require MPI_ALLOC_MEM for remote memory access
   ! Have a look at MPI-3, I think this removes cray pointers.
   b = malloc(my_num_vars)
   call MPI_ALLOC_MEM(window_size, MPI_INFO_NULL, b, ierr)

   do ii = 1, my_num_vars
      duplicate_mean(ii) = mean_ens_handle%copies(1, ii)
   enddo

   ! Expose local memory to RMA operation by other processes in a communicator.
   call mpi_win_create(duplicate_mean, window_size, bytesize, MPI_INFO_NULL, mpi_comm_world, mean_win, ierr)

else

   call init_ensemble_manager(mean_ens_handle, 1, state_ens_handle%num_vars, transpose_type_in = 3)
   call set_num_extra_copies(mean_ens_handle, 0)
   mean_ens_handle%copies(1,:) = state_ens_handle%copies(mean_copy, :)
   call all_copies_to_all_vars(mean_ens_handle) ! this is a transpose-duplicate

endif

! Set the current window to the state window
current_win = MEAN_WINDOW

data_count = copies_in_window(mean_ens_handle) ! One

end subroutine create_mean_window

!-------------------------------------------------------------
!> End epoch of state access.
!> Need to transpose qc and fwd operator back to copy complete
subroutine free_state_window(state_ens_handle, fwd_op_ens_handle, qc_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle
type(ensemble_type), intent(inout), optional :: fwd_op_ens_handle
type(ensemble_type), intent(inout), optional :: qc_ens_handle

integer :: ierr

if(get_allow_transpose(state_ens_handle)) then ! the forward operators were done var complete
   !transpose back if allowing transposes
   if (present(fwd_op_ens_handle)) &
      call all_vars_to_all_copies(fwd_op_ens_handle)
   if (present(qc_ens_handle)) &
      call all_vars_to_all_copies(qc_ens_handle)
else
   ! close mpi window
   call mpi_win_free(state_win, ierr)
   call MPI_FREE_MEM(duplicate_state, ierr)
endif

current_win = NO_WINDOW

end subroutine free_state_window

!---------------------------------------------------------
!> Free the mpi window
subroutine free_mean_window()

integer :: ierr

if(get_allow_transpose(mean_ens_handle)) then
   call end_ensemble_manager(mean_ens_handle)
else
   call mpi_win_free(mean_win, ierr)
   call MPI_FREE_MEM(duplicate_mean, ierr)
   call end_ensemble_manager(mean_ens_handle)
endif

current_win = NO_WINDOW

end subroutine free_mean_window

!---------------------------------------------------------
!> @}
end module window_mod

