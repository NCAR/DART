! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Window without cray pointer. Should you point the window at contigous memory?
module window_mod

!> \defgroup window window_mod
!> @{
use types_mod,            only : r8, i8
use ensemble_manager_mod, only : ensemble_type, copies_in_window, &
                                 init_ensemble_manager, &
                                 set_num_extra_copies, &
                                 end_ensemble_manager

implicit none

private
public :: create_mean_window, create_state_window, free_mean_window, &
          free_state_window, data_count, mean_win, state_win, current_win, &
          mean_ens_handle, NO_WINDOW, MEAN_WINDOW, STATE_WINDOW

integer :: data_count !! number of copies in the window
type(ensemble_type) :: mean_ens_handle


! mpi window handles
integer :: state_win   !< window for the forward operator
integer :: mean_win    !< window for the mean
integer :: current_win !< keep track of current window, start out assuming an invalid window

! parameters for keeping track of which window is open
integer, parameter :: NO_WINDOW    = -1
integer, parameter :: MEAN_WINDOW  = 0 
integer, parameter :: STATE_WINDOW = 2 

contains

!-------------------------------------------------------------
!> Always using distributed in non-mpi case
subroutine create_state_window(state_ens_handle, fwd_op_ens_handle, qc_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle
type(ensemble_type), intent(inout), optional :: fwd_op_ens_handle
type(ensemble_type), intent(inout), optional :: qc_ens_handle

! Find out how many copies to put in the window
! copies_in_window is not necessarily equal to ens_handle%num_copies
data_count = copies_in_window(state_ens_handle)

! Set the current window to the state window
current_win = STATE_WINDOW

end subroutine create_state_window

!-------------------------------------------------------------
!> Always using distributed in non-mpi case
subroutine create_mean_window(state_ens_handle, mean_copy, distribute_mean)

type(ensemble_type), intent(in)  :: state_ens_handle
integer,             intent(in)  :: mean_copy
logical,             intent(in)  :: distribute_mean

call init_ensemble_manager(mean_ens_handle, 1, state_ens_handle%num_vars) ! distributed ensemble
call set_num_extra_copies(mean_ens_handle, 0)
mean_ens_handle%copies(1,:) = state_ens_handle%copies(mean_copy, :)

! Set the current window to the state window
current_win = MEAN_WINDOW

data_count = 1

end subroutine create_mean_window

!-------------------------------------------------------------
!> End epoch of state access.
!> Need to transpose qc and fwd operator back to copy complete
subroutine free_state_window(state_ens_handle, fwd_op_ens_handle, qc_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle
type(ensemble_type), intent(inout), optional :: fwd_op_ens_handle
type(ensemble_type), intent(inout), optional :: qc_ens_handle

current_win = NO_WINDOW

end subroutine free_state_window

!---------------------------------------------------------
!> Free the mpi window
subroutine free_mean_window()

call end_ensemble_manager(mean_ens_handle)

current_win = NO_WINDOW

end subroutine free_mean_window

!-------------------------------------------------------------
!> @}
end module window_mod

