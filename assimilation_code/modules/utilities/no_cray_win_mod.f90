! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Window without cray pointer. Should you point the window at contigous memory?
module window_mod

!> \defgroup window window_mod
!> @{
use mpi_utilities_mod,    only : datasize, my_task_id, get_dart_mpi_comm
use types_mod,            only : r8, i8
use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index, &
                                 copies_in_window, init_ensemble_manager, &
                                 get_allow_transpose, end_ensemble_manager, &
                                 set_num_extra_copies, all_copies_to_all_vars, &
                                 all_vars_to_all_copies

use mpi_f08

implicit none

private
public :: create_mean_window, create_state_window, free_mean_window, &
          free_state_window, data_count, mean_win, state_win, current_win, &
          mean_ens_handle, NO_WINDOW, MEAN_WINDOW, STATE_WINDOW

! mpi window handles
type(MPI_Win) :: state_win   !< window for the forward operator
type(MPI_Win) :: mean_win    !< window for the mean
integer :: current_win !< keep track of current window, start out assuming an invalid window

! parameters for keeping track of which window is open
!>@todo should this be in the window_mod?  you will have to change in both cray 
!> and non cray versions 
integer, parameter :: NO_WINDOW    = -1
integer, parameter :: MEAN_WINDOW  = 0 
integer, parameter :: STATE_WINDOW = 2 

integer :: data_count !! number of copies in the window
integer(KIND=MPI_ADDRESS_KIND) :: window_size
logical :: use_distributed_mean = .false. ! initialize to false

! Global memory to stick the mpi window to.
! Need a simply contiguous piece of memory to pass to mpi_win_create
! Openmpi 1.10.0 will not compile with ifort 16 if
! you create a window with a 2d array.
real(r8), allocatable :: contiguous_fwd(:)
real(r8), allocatable :: mean_1d(:)

type(ensemble_type) :: mean_ens_handle

contains

!-------------------------------------------------------------
!> For the non-distributed case this is simply a transpose
!> For the distributed case memory is allocated in this module
!> then an mpi window is attached to this memory.
subroutine create_state_window(state_ens_handle, fwd_op_ens_handle, qc_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle
type(ensemble_type), intent(inout), optional :: fwd_op_ens_handle
type(ensemble_type), intent(inout), optional :: qc_ens_handle

integer :: ierr
integer :: bytesize !< size in bytes of each element in the window
integer :: my_num_vars !< my number of vars

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

   call mpi_type_size(datasize, bytesize, ierr)
   window_size = my_num_vars*data_count*bytesize

   allocate(contiguous_fwd(data_count*my_num_vars))
   contiguous_fwd = reshape(state_ens_handle%copies(1:data_count, :), (/my_num_vars*data_count/))

   ! Expose local memory to RMA operation by other processes in a communicator.
   call mpi_win_create(contiguous_fwd, window_size, bytesize, MPI_INFO_NULL, get_dart_mpi_comm(), state_win, ierr)
endif

! Set the current window to the state window
current_win = STATE_WINDOW

data_count = copies_in_window(state_ens_handle)

end subroutine create_state_window

!-------------------------------------------------------------
!> Using a mean ensemble handle.
!> 
subroutine create_mean_window(state_ens_handle, mean_copy, distribute_mean)

type(ensemble_type), intent(in)  :: state_ens_handle
integer,             intent(in)  :: mean_copy
logical,             intent(in)  :: distribute_mean

integer               :: ierr
integer               :: bytesize
integer               :: my_num_vars !< number of elements a task owns

! create an ensemble handle of just the mean copy.
use_distributed_mean = distribute_mean

if (use_distributed_mean) then
   call init_ensemble_manager(mean_ens_handle, 1, state_ens_handle%num_vars) ! distributed ensemble
   call set_num_extra_copies(mean_ens_handle, 0)
   mean_ens_handle%copies(1,:) = state_ens_handle%copies(mean_copy, :)
   allocate(mean_1d(state_ens_handle%my_num_vars))
   mean_1d(:) = mean_ens_handle%copies(1,:)

   ! find out how many variables I have
   my_num_vars = mean_ens_handle%my_num_vars
   call mpi_type_size(datasize, bytesize, ierr)
   window_size = my_num_vars*bytesize

   ! Need a simply contiguous piece of memory to pass to mpi_win_create
   ! Expose local memory to RMA operation by other processes in a communicator.
   call mpi_win_create(mean_1d, window_size, bytesize, MPI_INFO_NULL, get_dart_mpi_comm(), mean_win, ierr)
else
   call init_ensemble_manager(mean_ens_handle, 1, state_ens_handle%num_vars, transpose_type_in = 3)
   call set_num_extra_copies(mean_ens_handle, 0)
   mean_ens_handle%copies(1,:) = state_ens_handle%copies(mean_copy, :)
   call all_copies_to_all_vars(mean_ens_handle) ! this is a transpose-duplicate
endif
   
! grabbing mean directly, no windows are being used
current_win = MEAN_WINDOW

data_count = copies_in_window(mean_ens_handle) ! One.

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
   !transpose back if present
   if (present(fwd_op_ens_handle)) &
      call all_vars_to_all_copies(fwd_op_ens_handle)
   if (present(qc_ens_handle)) &
      call all_vars_to_all_copies(qc_ens_handle)
else
   ! close mpi window
   call mpi_win_free(state_win, ierr)
   deallocate(contiguous_fwd)
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
   deallocate(mean_1d)
   call end_ensemble_manager(mean_ens_handle)
endif

current_win = NO_WINDOW

end subroutine free_mean_window

!-------------------------------------------------------------
!> @}
end module window_mod

