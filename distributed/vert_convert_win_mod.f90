!> Aim to have a module which wraps up the distributed mean state
!> for use in the vertical conversion for localization in get_close_obs
module vert_convert_win_mod

use mpi_utilities_mod,  only : datasize, my_task_id
use types_mod,          only : r8
use data_structure_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index
use mpi

implicit none

integer mean_win !> mpi window
integer my_num_vars !> my number of vars
integer(KIND=MPI_ADDRESS_KIND) window_size
integer mean_row !> The row in copies that has the mean copy

real(r8) :: duplicate(*) !> duplicate array for cray pointer
pointer(a, duplicate)

contains

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_mean_window(state_ens_handle)

type(ensemble_type) :: state_ens_handle
integer :: ii, ierr
integer :: bytesize

! find out how many variables I have
my_num_vars = state_ens_handle%my_num_vars
mean_row = state_ens_handle%num_copies !> @todo How to do this nicely

! allocate some RDMA accessible memory
! using MPI_ALLOC_MEM because the MPI standard allows vendors to require MPI_ALLOC_MEM for remote memory access
call mpi_type_size(datasize, bytesize, ierr)
window_size = my_num_vars*bytesize
a = malloc(my_num_vars)
call MPI_ALLOC_MEM(window_size, MPI_INFO_NULL, a, ierr)

do ii = 1, my_num_vars
   duplicate(ii) = state_ens_handle%copies(mean_row, ii)
enddo

end subroutine create_mean_window

!---------------------------------------------------------
!> Free the mpi window
subroutine free_mean_window
integer :: ierr

call mpi_win_free(mean_win, ierr)
call MPI_FREE_MEM(duplicate, ierr) ! not a

end subroutine free_mean_window

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_state(x, index, state_ens_handle)

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

end subroutine get_state

end module vert_convert_win_mod

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
! These are test functions/subroutines to test the behaviour of 
! the rest of the module

