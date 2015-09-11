! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Null version of distributed state.
!> I am concerned about programs that call get_state_meta_data_distrib to get
!> state kinds for models like WRF and MPAS, because do vertical conversion
!> inside get_state_meta_data.
module distributed_state_mod

!> \defgroup distrib_state distributed_state_mod
!> @{
use types_mod,          only : r8, i8
use ensemble_manager_mod, only : ensemble_type

implicit none

private
public :: get_state_array, get_state, create_state_window, free_state_window, create_mean_window, free_mean_window

interface get_state
   module procedure get_fwd
   module procedure get_mean
end interface

contains

!---------------------------------------------------------
!> Do nothing
subroutine get_state_array(x, index, state_ens_handle)

real(r8),            intent(out) :: x(:) !> all copies of an element of the state vector
integer(i8),         intent(in)  :: index(:) !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

end subroutine get_state_array

!-------------------------------------------------------------
!> Do nothing
subroutine create_state_window(state_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle

end subroutine create_state_window

!-------------------------------------------------------------
!> Do nothing
subroutine create_mean_window(state_ens_handle, mean_copy, distribute_mean)

type(ensemble_type), intent(in)  :: state_ens_handle
integer,             intent(in)  :: mean_copy
logical,             intent(in)  :: distribute_mean


end subroutine create_mean_window

!-------------------------------------------------------------
!> Do nothing
subroutine free_state_window(state_ens_handle, fwd_op_ens_handle, qc_ens_handle)

type(ensemble_type), intent(inout) :: state_ens_handle
type(ensemble_type), intent(inout) :: fwd_op_ens_handle
type(ensemble_type), intent(inout) :: qc_ens_handle



end subroutine free_state_window

!---------------------------------------------------------
!> Do nothing
subroutine free_mean_window

end subroutine free_mean_window


!---------------------------------------------------------
!> Do nothing
subroutine get_fwd(x, index, state_ens_handle)

real(r8),    intent(out)         :: x(:) !> all copies of an element of the state vector
integer(i8), intent(in)          :: index !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle

integer                          :: owner_of_state !> task who owns the state
integer                          :: element_index !> local index of element
integer                          :: target_disp
integer                          :: ierr


end subroutine get_fwd

!---------------------------------------------------------
!> Do nothing
subroutine get_mean(x, index, state_ens_handle)

real(r8),    intent(out)         :: x !> only grabing the mean
integer(i8), intent(in)          :: index !> index into state vector
type(ensemble_type), intent(in)  :: state_ens_handle


end subroutine get_mean

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
! These are test functions/subroutines to test the behaviour of 
! the rest of the module

!> @}

end module distributed_state_mod
