! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Contains the window information for the state.  Two windows:
!> One for all copies, one for the mean.
!> Not sure whether we should just have one window to avoid multiple synchronizations.
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
public :: create_mean_window, create_state_window, free_mean_window, free_state_window, &
          mean_win, state_win, row, num_rows

! Use these if %copies is modified during window access.
public :: get_local_state, get_local_mean

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer state_win !< mpi window for the forward operator
integer mean_win !< mpi window
!< @todo the number of copies in the window is sloppy. You need to make this better.
integer :: num_rows !> number of copies in the window
integer :: row !> mean row index

integer(KIND=MPI_ADDRESS_KIND) window_size_state
integer(KIND=MPI_ADDRESS_KIND) window_size_mean

real(r8) :: duplicate_state(*)  !< duplicate array for cray pointer fwd
pointer(a, duplicate_state)

real(r8) :: duplicate_mean(*)  !< duplicate array for cray pointer vert convert
pointer(b, duplicate_mean)

contains

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_state_window(state_ens_handle)

type(ensemble_type), intent(in) :: state_ens_handle !< state ensemble handle

integer :: ii, jj, count, ierr
integer :: bytesize !< size in bytes of each element in the window
integer :: my_num_vars

! find how many variables I have
my_num_vars = state_ens_handle%my_num_vars
! find out how many rows to put in the window
num_rows = copies_in_window(state_ens_handle)

! allocate some RDMA accessible memory
! using MPI_ALLOC_MEM because the MPI standard allows vendors to require MPI_ALLOC_MEM for remote memory access
call mpi_type_size(datasize, bytesize, ierr)
window_size_state = my_num_vars*num_rows*bytesize
a = malloc(my_num_vars*num_rows)
call MPI_ALLOC_MEM(window_size_state, MPI_INFO_NULL, a, ierr)

count = 1
! create a duplicate copies array for remote memory access
! Doing this because you cannot use a cray pointer with an allocatable array
! Can't do array assignment with a cray pointer, so you need to loop
!> @todo should you only fill the window with actual ensemble members?
do ii = 1, my_num_vars
   do jj = 1, num_rows
      duplicate_state(count) = state_ens_handle%copies(jj,ii)
      count = count + 1
   enddo
enddo

! Expose local memory to RMA operation by other processes in a communicator.
call mpi_win_create(duplicate_state, window_size_state, bytesize, MPI_INFO_NULL, mpi_comm_world, state_win, ierr)

end subroutine create_state_window

!-------------------------------------------------------------
!> Create the window for the ensemble complete state vector
!> I think you have to pass it the state ensemble handle
subroutine create_mean_window(state_ens_handle)

type(ensemble_type), intent(in) :: state_ens_handle

integer :: ii, ierr
integer :: bytesize
integer :: my_num_vars
! find out how many variables I have
my_num_vars = state_ens_handle%my_num_vars
row = mean_row(state_ens_handle)

! allocate some RDMA accessible memory
! using MPI_ALLOC_MEM because the MPI standard allows vendors to require MPI_ALLOC_MEM for remote memory access
! Have a look at MPI-3, I think this removes cray pointers.
call mpi_type_size(datasize, bytesize, ierr)
window_size_mean = my_num_vars*bytesize
b = malloc(my_num_vars)
call MPI_ALLOC_MEM(window_size_mean, MPI_INFO_NULL, b, ierr)

do ii = 1, my_num_vars
   duplicate_mean(ii) = state_ens_handle%copies(row, ii)
enddo

! Expose local memory to RMA operation by other processes in a communicator.
call mpi_win_create(duplicate_mean, window_size_mean, bytesize, MPI_INFO_NULL, mpi_comm_world, mean_win, ierr)

end subroutine create_mean_window

!---------------------------------------------------------
!> Free the mpi window
subroutine free_state_window
integer :: ierr

call mpi_win_free(state_win, ierr)
call MPI_FREE_MEM(duplicate_state, ierr) ! not a

end subroutine free_state_window

!---------------------------------------------------------
!> Free the mpi window
subroutine free_mean_window
integer :: ierr

call mpi_win_free(mean_win, ierr)
call MPI_FREE_MEM(duplicate_mean, ierr) ! not a

end subroutine free_mean_window

!-------------------------------------------------------------
! The functions below are if you need to access data from your
! own local window.
!-------------------------------------------------------------
!>
function get_local_state(index)

integer, intent(in) :: index
real(r8) :: get_local_state(num_rows)

integer :: ierr

call mpi_win_lock(MPI_LOCK_SHARED, my_task_id(), 0, state_win, ierr)
get_local_state = duplicate_state((index-1)*num_rows+1:index*num_rows+1)
call mpi_win_unlock(my_task_id(), state_win, ierr)

end function get_local_state

!-------------------------------------------------------------
function get_local_mean(index)

integer, intent(in) :: index
real(r8) :: get_local_mean

integer :: ierr

call mpi_win_lock(MPI_LOCK_SHARED, my_task_id(), 0, mean_win, ierr)
get_local_mean = duplicate_mean(index)
call mpi_win_unlock(my_task_id(), mean_win, ierr)

end function get_local_mean

!---------------------------------------------------------
!> @}
end module window_mod
