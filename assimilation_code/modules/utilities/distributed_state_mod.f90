! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module distributed_state_mod

!> \defgroup distrib_state distributed_state_mod
!> @{
!>
!> Allows a distributed state, or a var complete state
!> for the forward operator or the mean state
!> Note that there are two window modules that can be compiled:
!>    no_cray_win_mod.f90 does not use cray pointers
!>    cray_win_mod.f90 uses cray pointers

use mpi_utilities_mod,    only : my_task_id, &
                                 get_from_fwd, get_from_mean
use types_mod,            only : r8, i8
use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index, &
                                 get_allow_transpose
use window_mod,           only : create_mean_window, create_state_window, free_mean_window, &
                                 free_state_window, mean_ens_handle, data_count, &
                                 NO_WINDOW, MEAN_WINDOW, STATE_WINDOW, &
                                 mean_win, state_win, current_win
use utilities_mod,        only : error_handler, E_ERR

implicit none

private
public :: get_state_array, get_state, create_state_window, &
          free_state_window, create_mean_window, free_mean_window

contains

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete. This differes from get_state as it now works on an
!> array of state indices rather than a single index.
subroutine get_state_array(x, my_index, state_ens_handle)

real(r8),            intent(out) :: x(data_count) !! all copies of an element of the state vector
integer(i8),         intent(in)  :: my_index(data_count) !! index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

real(r8) :: x_tmp(data_count) !! all copies of an element of the state vector
logical  :: not_done(data_count)
integer  :: i, e

not_done = .true.

! only grab state once unique state indices
do e = 1, data_count
   if (not_done(e)) then
      x_tmp = get_state(my_index(e), state_ens_handle)
      do i = e, data_count
         if(my_index(i) == my_index(e)) then
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
function get_state(my_index, ens_handle) result (x)

real(r8) :: x(data_count) !! all copies of an element of the state vector

integer(i8),         intent(in)  :: my_index !! index into state vector
type(ensemble_type), intent(in)  :: ens_handle

if (current_win == MEAN_WINDOW) then
   call get_mean(x, my_index, ens_handle)
else if (current_win == STATE_WINDOW) then
   call get_fwd(x, my_index, ens_handle)
else
   call error_handler(E_ERR, 'get_state',' No window currently open')
endif

end function get_state

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_fwd(x, my_index, state_ens_handle)

real(r8),            intent(out) :: x(data_count) !! all copies of an element of the state vector
integer(i8),         intent(in)  :: my_index !! index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                        :: owner_of_state !! task who owns the state
integer                        :: element_index !! local index of element

if (get_allow_transpose(state_ens_handle)) then
   x = state_ens_handle%vars(my_index, 1:data_count)
else

   call get_var_owner_index(state_ens_handle, my_index, owner_of_state, element_index) ! pe

   owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

   if (my_task_id() == owner_of_state) then
      !For future use if we need to access the local window
      ! e.g. if %copies is modified while we are doing the forward operator.
      ! Currently, %copies is static => same as what is in the window.
      !x = get_local_state(element_index)
      x = state_ens_handle%copies(1:data_count, element_index)
   else
     call get_from_fwd(owner_of_state, state_win, element_index, data_count, x)
   endif
endif

end subroutine get_fwd

!---------------------------------------------------------
!> Gets all copies of an element of the state vector from the process who owns it
!> Assumes ensemble complete
subroutine get_mean(x, my_index, state_ens_handle)

real(r8),            intent(out) :: x(1) !! only grabing the mean
integer(i8),         intent(in)  :: my_index !! index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                        :: owner_of_state !! task who owns the state
integer                        :: element_index !! local index of element

if (get_allow_transpose(mean_ens_handle)) then
   x(1) = mean_ens_handle%vars(my_index, 1)
else

   call get_var_owner_index(state_ens_handle, my_index, owner_of_state, element_index) ! pe

   owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

   if (my_task_id() == owner_of_state) then
      x(1) = mean_ens_handle%copies(1, element_index)
   else
      call get_from_mean(owner_of_state, mean_win, element_index, x(1))
   endif

endif

end subroutine get_mean

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
! These are test functions/subroutines to test the behaviour of 
! the rest of the module

!> @}

end module distributed_state_mod

