! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Allows a distributed state, or a var complete state
!> for the forward operator or the mean state
!> Note that there are two window modules that can be compiled:
!>    no_cray_win_mod.f90 does not use cray pointers
!>    cray_win_mod.f90 uses cray pointers
module distributed_state_mod

!> \defgroup distrib_state distributed_state_mod
!> @{
use mpi_utilities_mod,    only : datasize, my_task_id
use types_mod,            only : r8, i8
use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index, &
                                 get_allow_transpose
use window_mod,           only : create_mean_window, create_state_window, free_mean_window, &
                                 free_state_window, mean_win, state_win, data_count, &
                                 mean_ens_handle

use mpi

implicit none

private
public :: get_state_array, get_state, create_state_window, free_state_window, create_mean_window, free_mean_window

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
!> Assumes ensemble complete. This differes from get_state as it now works on an
!> array of state indices rather than a single index.
subroutine get_state_array(x, index, state_ens_handle)

real(r8),            intent(out) :: x(data_count) !> all copies of an element of the state vector
integer(i8),         intent(in)  :: index(data_count) !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

real(r8) :: x_tmp(data_count) !> all copies of an element of the state vector
logical  :: not_done(data_count)
integer  :: i, e

not_done = .true.

! only grab state once unique state indices
do e = 1, data_count
   if (not_done(e)) then
      call get_state(x_tmp, index(e), state_ens_handle)
      do i = e, data_count
         if(index(i) == index(e)) then
            x(i) = x_tmp(i)
            not_done(i) = .false. 
         endif
     enddo
   endif
enddo

end subroutine get_state_array

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_fwd(x, index, state_ens_handle)

real(r8),            intent(out) :: x(data_count) !> all copies of an element of the state vector
integer(i8),         intent(in)  :: index !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                        :: owner_of_state !> task who owns the state
integer                        :: element_index !> local index of element
integer(KIND=MPI_ADDRESS_KIND) :: target_disp
integer                        :: ierr

if (get_allow_transpose(state_ens_handle)) then
   x = state_ens_handle%vars(index, 1:data_count)
else

   call get_var_owner_index(index, owner_of_state, element_index) ! pe

   owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

   if (my_task_id() == owner_of_state) then
      !For future use if we need to access the local window
      ! e.g. if %copies is modified while we are doing the forward operator.
      ! Currently, %copies is static => same as what is in the window.
      !x = get_local_state(element_index)
      x = state_ens_handle%copies(1:data_count, element_index)
   else
      ! Note all of copies array is in the window, not just the real ensemble members
      target_disp = (element_index - 1) * data_count
      call mpi_win_lock(MPI_LOCK_SHARED, owner_of_state, 0, state_win, ierr)
      call mpi_get(x, data_count, datasize, owner_of_state, target_disp, data_count, datasize, state_win, ierr)
      call mpi_win_unlock(owner_of_state, state_win, ierr)
   endif
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

if (get_allow_transpose(mean_ens_handle)) then
   x = mean_ens_handle%vars(index, 1)
else

   call get_var_owner_index(index, owner_of_state, element_index) ! pe

   owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

   if (my_task_id() == owner_of_state) then
      x = mean_ens_handle%copies(1, element_index)
   else
      target_disp = (element_index - 1)
      call mpi_win_lock(MPI_LOCK_SHARED, owner_of_state, 0, mean_win, ierr)
      call mpi_get(x, 1, datasize, owner_of_state, target_disp, 1, datasize, mean_win, ierr)
      call mpi_win_unlock(owner_of_state, mean_win, ierr)
   endif

endif

end subroutine get_mean

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
! These are test functions/subroutines to test the behaviour of 
! the rest of the module

!> @}

end module distributed_state_mod
