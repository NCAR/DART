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
use types_mod,          only : r8, i8
use ensemble_manager_mod, only : ensemble_type

implicit none

interface get_state
   module procedure get_fwd
   module procedure get_mean
end interface

contains

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
