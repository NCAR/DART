! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Window without cray pointer. Should you point the window at contigous memory?
module window_mod

!> \defgroup window window_mod
!> @{
use mpi_utilities_mod,  only : datasize, my_task_id
use types_mod,          only : r8
use data_structure_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index, &
                               copies_in_window, mean_row

use mpi

implicit none

private
public :: create_mean_window, create_state_window, free_mean_window, free_state_window, mean_win, state_win, row, num_rows

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer state_win !< mpi window for the forward operator
integer mean_win !< mpi window
integer :: num_rows !> number of copies in the window
integer :: row !> mean row index
integer(KIND=MPI_ADDRESS_KIND) window_size_state
integer(KIND=MPI_ADDRESS_KIND) window_size_mean

! Global memory to stick the mpi window to.
! Need a simply contiguous piece of memory to pass to mpi_win_create
! For the forward operator copies, putting the whole of %copies array in
! the window. This is to avoid making a dupilcate of 
! state_ens_handle%copies(1:ens_size)
real(r8), allocatable :: contiguous_mean(:)

contains

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_state_window(state_ens_handle)

type(ensemble_type), intent(in) :: state_ens_handle !< state ensemble handle

integer :: ii, jj, count, ierr
integer :: bytesize !< size in bytes of each element in the window
integer :: my_num_vars !< my number of vars

! find how many variables I have
my_num_vars = state_ens_handle%my_num_vars
! The number of rows in the window is now state_ens_handle%num_copies.
! The number of rows that a call to get_state will grab is still num_rows(ens_size)
! FIXME For now, I am keeping the name copies_in_window().  This should be
! changed to a more sensible name
num_rows = copies_in_window(state_ens_handle)

call mpi_type_size(datasize, bytesize, ierr)
window_size_state = my_num_vars*state_ens_handle%num_copies*bytesize

! Expose local memory to RMA operation by other processes in a communicator.
call mpi_win_create(state_ens_handle%copies, window_size_state, bytesize, MPI_INFO_NULL, mpi_comm_world, state_win, ierr)

end subroutine create_state_window

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_mean_window(state_ens_handle)

type(ensemble_type), intent(in) :: state_ens_handle
integer               :: ii, ierr
integer               :: bytesize
integer               :: my_num_vars !< my number of vars

! find out how many variables I have
my_num_vars = state_ens_handle%my_num_vars
row = mean_row(state_ens_handle)

call mpi_type_size(datasize, bytesize, ierr)
window_size_mean = my_num_vars*bytesize

! Need a simply contiguous piece of memory to pass to mpi_win_create
allocate(contiguous_mean(my_num_vars))
contiguous_mean = state_ens_handle%copies(row, :)

! Expose local memory to RMA operation by other processes in a communicator.
call mpi_win_create(contiguous_mean, window_size_mean, bytesize, MPI_INFO_NULL, mpi_comm_world, mean_win, ierr)

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
deallocate(contiguous_mean)

end subroutine free_mean_window

!-------------------------------------------------------------
!> @}
end module window_mod
