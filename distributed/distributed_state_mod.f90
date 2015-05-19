! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Aim is to have the distributed forward operator 
!> window contained within this module.
!> We can have a dummy window module for the
!> non-distributed version
module distributed_state_mod

!> \defgroup distrib_state distributed_state_mod
!> @{
use mpi_utilities_mod,  only : datasize, my_task_id
use types_mod,          only : r8, i8
use data_structure_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index
use window_mod

use mpi

implicit none

private
public get_state, create_state_window, free_state_window, create_mean_window, free_mean_window

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

interface get_state
   module procedure get_fwd
   module procedure get_mean
end interface

contains

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_fwd(x, index, state_ens_handle)

real(r8),            intent(out) :: x(num_rows) !> all copies of an element of the state vector
integer(i8),         intent(in)  :: index !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                        :: owner_of_state !> task who owns the state
integer                        :: element_index !> local index of element
integer(KIND=MPI_ADDRESS_KIND) :: target_disp
integer                        :: ierr

call get_var_owner_index(index, owner_of_state, element_index) ! pe

owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

if (my_task_id() == owner_of_state) then
   x = state_ens_handle%copies(1:num_rows, element_index)
else
   ! Note all of copies array is in the window, not just the real ensemble members
   target_disp = (element_index - 1) * state_ens_handle%num_copies
   call mpi_win_lock(MPI_LOCK_SHARED, owner_of_state, 0, state_win, ierr)
   call mpi_get(x, num_rows, datasize, owner_of_state, target_disp, num_rows, datasize, state_win, ierr)
   call mpi_win_unlock(owner_of_state, state_win, ierr)
endif

end subroutine get_fwd

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_mean(x, index, state_ens_handle)

real(r8),            intent(out) :: x !> only grabing the mean
integer(i8),         intent(in)  :: index !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                        :: owner_of_state !> task who owns the state
integer                        :: element_index !> local index of element
integer(KIND=MPI_ADDRESS_KIND) :: target_disp
integer                        :: ierr

call get_var_owner_index(index, owner_of_state, element_index) ! pe

owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

if (my_task_id() == owner_of_state) then
   x = state_ens_handle%copies(row, element_index)
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

!> @}

end module distributed_state_mod
