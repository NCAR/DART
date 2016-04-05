! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

! One dimensional oscillator model (model size is 2)

use    types_mod, only : r8
use location_mod, only : loc_type, get_dist, set_loc

implicit none
private

public :: init_model, get_model_size, init_conditions, adv_1step, advance, &
   adv_true_state, output, diag_output_index, get_close_pts, state_loc, &
   model_output, model_get_close_states

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

integer,  parameter :: model_size = 2
real(r8), parameter ::    delta_t = 0.01_r8
real(r8), parameter ::      alpha = 0.1_r8

logical :: output_init = .FALSE.

! Define output indices for diagnostics
integer :: diag_output_index(2)

! Define the location of the state variables in module storage
type(loc_type) :: state_loc(model_size)

contains

!==================================================================

subroutine init_model()

! Stub for model initialization, not needed for L96

end subroutine init_model

!==================================================================

subroutine comp_dt(x, dt)

implicit none

! Computes the time tendency of a 1D oscillator

real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: dt(:)

dt(1) = x(2)
dt(2) = - x(1) + alpha * x(2)

end subroutine comp_dt

!
!===================================================================
!
!
!  Does single time step advance for lorenz 96 model
!  using four step rk time step

subroutine adv_1step(x)

implicit none

real(r8), intent(inout) :: x(:)

real(r8) :: dx(size(x))

!  Compute the first intermediate step
call comp_dt(x, dx)

x = x + dx * delta_t

end subroutine adv_1step

!
!----------------------------------------------------------------------
!

subroutine adv_true_state(x)

implicit none

real(r8), intent(inout) :: x(:)

!write(*, *) 'amp = ', sqrt(x(1)**2 + x(2)**2), ' phase = ', atan2(x(1), x(2))

call adv_1step(x)

end subroutine adv_true_state

!
!=====================================================================
!

subroutine init_conditions(x)

!  Initial conditions for lorenz 96
! It is assumed that this is called before any other routines in this
! module. Should probably make that more formal and perhaps enforce for
! more comprehensive models.

implicit none

real(r8), intent(out) :: x(:)

integer  :: i
real(r8) :: x_loc

! Define the interesting indexes for variables to do diag output; span lats
do i = 1, 2
   diag_output_index(i) = i
end do

! Set position and velocity to 1.0
x = 0.0

! Define the locations of the state variables;
do i = 1, model_size
   x_loc = 0.0_r8
   call set_loc(state_loc(i), x_loc)
end do

end subroutine init_conditions

!
!====================================================================
!

!  Advances the lorenz-96 model by a given number of steps
!  Current state in x, new state in xnew, num time steps advanced

subroutine advance(x, num, xnew)

implicit none

real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: xnew(:)
integer,  intent(in)  :: num

integer :: i

!  Copy initial conditions to avoid overwrite
xnew = x

!  Advance the appropriate number of steps
do i = 1, num
   call adv_1step(xnew)
end do

end subroutine advance

!
!---------------------------------------------------------------------
!

subroutine get_close_pts(list, num)

implicit none

integer, intent(in) :: num
integer, intent(inout) :: list(model_size, num)

end subroutine get_close_pts

!
!===================================================================
!

function get_model_size()

! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size




subroutine model_get_close_states(o_loc, radius, numb, indices, dist, x)
!--------------------------------------------------------------------
! 
! Stub for computation of get close states

implicit none

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: numb, indices(:)
real(r8), intent(out) :: dist(:)
real(r8), intent(in) :: x(:)

! Because of F90 limits this stub must be here telling assim_model
! to do exhaustive search (numb = -1 return)
numb = -1

end subroutine model_get_close_states


!
!===================================================================
!
end module model_mod
