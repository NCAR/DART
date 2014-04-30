!> Aim is to have the distributed forward operator 
!> window contained within this module.
!> We can have a dummy window module for the
!> non-distributed version
module fwd_op_win_mod

use mpi_utilities_mod,  only : datasize, my_task_id
use types_mod,          only : r8
use data_structure_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index
use mpi

implicit none

integer state_win !< mpi window for the forward operator
integer mean_win !< mpi window
integer my_num_vars !< my number of vars
integer copies_in_window !< number of copies in the window - ens_size
!< @todo the number of copies in the window is sloppy. You need to make this better.
integer(KIND=MPI_ADDRESS_KIND) window_size_state
integer(KIND=MPI_ADDRESS_KIND) window_size_mean
integer mean_row !< The row in copies that has the mean copy

interface get_state
   module procedure get_fwd
   module procedure get_mean
end interface

contains

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_state_window(state_ens_handle)

type(ensemble_type) :: state_ens_handle !< state ensemble handle
integer :: ii, jj, count, ierr
integer :: bytesize !< size in bytes of each element in the window

! find how many variables I have
my_num_vars = state_ens_handle%my_num_vars

!> @todo number of copies in the window
copies_in_window = state_ens_handle%num_copies -5

call mpi_type_size(datasize, bytesize, ierr)
window_size_state = my_num_vars*copies_in_window*bytesize

! Expose local memory to RMA operation by other processes in a communicator.
call mpi_win_create(state_ens_handle%copies(1:copies_in_window, :), window_size_state, bytesize, MPI_INFO_NULL, mpi_comm_world, state_win, ierr)

end subroutine create_state_window

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_mean_window(state_ens_handle)

type(ensemble_type) :: state_ens_handle
integer :: ii, ierr
integer :: bytesize

! find out how many variables I have
my_num_vars = state_ens_handle%my_num_vars
mean_row = state_ens_handle%num_copies -5!> @todo How to do this nicely

! allocate some RDMA accessible memory
! using MPI_ALLOC_MEM because the MPI standard allows vendors to require MPI_ALLOC_MEM for remote memory access
! Have a look at MPI-3, I think this removes cray pointers.
call mpi_type_size(datasize, bytesize, ierr)
window_size_mean = my_num_vars*bytesize

! Expose local memory to RMA operation by other processes in a communicator.
call mpi_win_create(state_ens_handle%copies(mean_row, :), window_size_mean, bytesize, MPI_INFO_NULL, mpi_comm_world, mean_win, ierr)

end subroutine create_mean_window


!---------------------------------------------------------
!> Free the mpi window
subroutine free_state_window
integer :: ierr

call mpi_win_free(state_win, ierr)

end subroutine free_state_window

!---------------------------------------------------------
!> Free the mpi window
subroutine free_mean_window
integer :: ierr

call mpi_win_free(mean_win, ierr)

end subroutine free_mean_window

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_fwd(x, index, state_ens_handle)

real(r8), intent(out)            :: x(copies_in_window) !> all copies of an element of the state vector
integer, intent(in)              :: index !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                          :: owner_of_state !> task who owns the state
integer                          :: element_index !> local index of element
integer(KIND=MPI_ADDRESS_KIND)   :: target_disp
integer                          :: ierr

call get_var_owner_index(index, owner_of_state, element_index) ! pe

owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

if (my_task_id() == owner_of_state) then
   x = state_ens_handle%copies(1:copies_in_window, element_index)
else
   target_disp = (element_index - 1) * copies_in_window
   call mpi_win_lock(MPI_LOCK_SHARED, owner_of_state, 0, state_win, ierr)
   call mpi_get(x, copies_in_window, datasize, owner_of_state, target_disp, copies_in_window, datasize, state_win, ierr)
   call mpi_win_unlock(owner_of_state, state_win, ierr)
endif

end subroutine get_fwd

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_mean(x, index, state_ens_handle)

real(r8), intent(out)            :: x !> only grabing the mean
integer, intent(in)              :: index !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                          :: owner_of_state !> task who owns the state
integer                          :: element_index !> local index of element
integer(KIND=MPI_ADDRESS_KIND)   :: target_disp
integer                          :: ierr

call get_var_owner_index(index, owner_of_state, element_index) ! pe

owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

if (my_task_id() == owner_of_state) then
   x = state_ens_handle%copies(mean_row, element_index)
else
   target_disp = (element_index - 1)
   call mpi_win_lock(MPI_LOCK_SHARED, owner_of_state, 0, mean_win, ierr)
   call mpi_get(x, 1, datasize, owner_of_state, target_disp, 1, datasize, mean_win, ierr)
   call mpi_win_unlock(owner_of_state, mean_win, ierr)
endif

end subroutine get_mean

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
! These are test functions/subroutines to test the behaviour of 
! the rest of the module



end module fwd_op_win_mod
