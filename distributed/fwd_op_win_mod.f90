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

integer state_win !> mpi window
integer my_num_vars !> my number of vars
integer copies_in_window !> number of copies in the window - ens_size
!> @todo the number of copies in the window is sloppy. You need to make this better.
integer(KIND=MPI_ADDRESS_KIND) window_size

real(r8) :: duplicate(*)  !> duplicate array for cray pointer
pointer(a, duplicate)

contains

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_state_window(state_ens_handle)

type(ensemble_type) :: state_ens_handle !> state ensemble handle
integer :: ii, jj, count, ierr
integer :: bytesize !> size in bytes of each element in the window

! find how many variables I have
my_num_vars = state_ens_handle%my_num_vars

!> @todo number of copies in the window
copies_in_window = state_ens_handle%num_copies -5

! allocate some RDMA accessible memory
! using MPI_ALLOC_MEM because the MPI standard allows vendors to require MPI_ALLOC_MEM for remote memory access
call mpi_type_size(datasize, bytesize, ierr)
window_size = my_num_vars*copies_in_window*bytesize
a = malloc(my_num_vars)
call MPI_ALLOC_MEM(window_size, MPI_INFO_NULL, a, ierr)

count = 1
! create a duplicate copies array for remote memory access
! Doing this because you cannot use a cray pointer with an allocatable array
! Can't do array assignment with a cray pointer, so you need to loop
!> @todo should you only fill the window with actual ensemble members?
do ii = 1, my_num_vars
   do jj = 1, copies_in_window
      duplicate(count) = state_ens_handle%copies(jj,ii)
      count = count + 1
   enddo
enddo

! Expose local memory to RMA operation by other processes in a communicator.
call mpi_win_create(duplicate, window_size, bytesize, MPI_INFO_NULL, mpi_comm_world, state_win, ierr)

end subroutine create_state_window

!---------------------------------------------------------
!> Free the mpi window
subroutine free_state_window
integer :: ierr

call mpi_win_free(state_win, ierr)
call MPI_FREE_MEM(duplicate, ierr) ! not a

end subroutine free_state_window

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_state(x, index, state_ens_handle)

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

end subroutine get_state

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
! These are test functions/subroutines to test the behaviour of 
! the rest of the module



end module fwd_op_win_mod
