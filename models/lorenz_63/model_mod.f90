module model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

! Revised assim_model version of Lorenz-63 3-variable model

use types_mod
use location_mod, only : location_type, get_dist, set_location, get_location
use time_manager_mod

private


public static_init_model, init_conditions, get_model_size, adv_1step, state_loc, &
   init_time,  model_interpolate, get_model_time_step, get_state_meta_data, end_model, &
   model_get_close_states

integer, parameter :: model_size = 3

! Define the location of the state variables in module storage
type(location_type) :: state_loc(model_size)
type(time_type) :: time_step


!  define model parameters

real(r8), parameter ::  sigma = 10.0_r8
real(r8), parameter ::      r = 28.0_r8
real(r8), parameter ::      b = 8.0_r8 / 3.0_r8
real(r8), parameter :: deltat = 0.01_r8

contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for L63 model and outputs I.D.

implicit none

character(len=128) :: source,revision,revdate
integer :: i
real(r8) :: x_loc

! let CVS fill strings ... DO NOT EDIT ...
source   = "$Source$"
revision = "$Revision$"
revdate  = "$Date$"

! Ultimately,  change output to diagnostic output block ...

write(*,*)'assim_model attributes:'
write(*,*)'   ',source
write(*,*)'   ',revision
write(*,*)'   ',revdate

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L93
time_step = set_time(3600, 0)

end subroutine static_init_model



subroutine comp_dt(x, dt)
!==================================================================
! subroutine comp_dt(x, dt)
!
! computes time tendency of the lorenz 1963 3-variable model given 
! current state

implicit none

real(r8), intent(in) :: x(:)
real(r8), intent(out) :: dt(:)

! compute the lorenz model dt from standard equations

dt(1) = sigma * (x(2) - x(1))
dt(2) = -1.0_r8*x(1)*x(3) + r*x(1) - x(2)
dt(3) = x(1)*x(2) - b*x(3)

return
end subroutine comp_dt



subroutine advance(x, num, xnew, time)
!===================================================================
! subroutine advance(x, num, xnew, time)
!
! advances the 3 variable lorenz-63 model by a given number of steps
! current state in x, new state in xnew, num time steps advanced

implicit none

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: num
real(r8), intent(out) :: xnew(:)
type(time_type), intent(in) :: time

integer :: i

xnew = x                  !  copy initial conditions to avoid overwrite
   
do i = 1, num             !  advance the appropriate number of steps
   call adv_1step(xnew, time)
end do

return
end subroutine advance



subroutine init_conditions(x)
!===================================================================
! subroutine init_conditions(x)
!
!  off-attractor initial conditions for lorenz 63

implicit none

real(r8), intent(out) :: x(:)

! Initial conditions that move nicely onto attractor
x = 0.10_r8

end subroutine init_conditions



subroutine linear_dt(x, dx, dt)
!====================================================================
! subroutine linear_dt(x, dx, dt)
!
! old version of linearized lorenz 63 model time tendency computation
   
implicit none
   
real(r8), intent(in)  :: x(:), dx(:)
real(r8), intent(out) :: dt(:)

!  compute linear model lorenz time tendency

dt(1) = -1.0_r8 * sigma * dx(1) + sigma*dx(2)
dt(2) = (r - x(3))*dx(1) - dx(2) - x(1)*dx(3)
dt(3) = x(2)*dx(1) + x(1)*dx(2) - b*dx(3)

end subroutine linear_dt




subroutine adv_1step(x, time)
!====================================================================
! subroutine adv_1step(x, time)
!
! does single time step advance for lorenz convective 3 variable model
! using two step rk time step

implicit none

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8) :: fract

fract = 1.0_r8
call adv_single(x, fract)

return
end  subroutine adv_1step



subroutine adv_single(x, fract)
!====================================================================
! subroutine adv_single(x, fract)
!
! does single time step advance for lorenz convective 3 variable model
! using two step rk time step

implicit none

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(3), x2(3), dx(3)

integer i

call comp_dt(x, dx)            !  compute the first intermediate step
x1 = x + fract * deltat * dx

call comp_dt(x1, dx)           !  compute the second intermediate step
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

return
end subroutine adv_single



function get_model_size()
!=====================================================================
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size




subroutine init_time(time)
!----------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

implicit none

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time




function model_interpolate(x, location, type)
!---------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument type is not used here because there is only one type of variable.
! Type is needed to allow swap consistency with more complex models.

implicit none

real(r8) :: model_interpolate
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type

integer :: lower_index, upper_index
real(r8) :: loc, fraction

! Convert location to real
loc = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
loc = model_size * loc

lower_index = int(loc)
upper_index = lower_index + 1
if(upper_index > model_size) upper_index = 1
if(lower_index == 0) lower_index = model_size

fraction = loc - int(loc)
model_interpolate = (1.0_r8 - fraction) * x(lower_index) + fraction * x(upper_index)

end function model_interpolate



function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step




subroutine get_state_meta_data(index, location)
!---------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?

implicit none

integer, intent(in) :: index
type(location_type), intent(out) :: location

location = state_loc(index)

end subroutine get_state_meta_data




subroutine end_model()
!------------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for L96 for now.


end subroutine end_model



subroutine adv_single_rk4(x, fract)
!=====================================================================
! subroutine adv_single_rk4(x, fract)
!
! does single time step advance for lorenz convective 3 variable model
! using four step rk time step

implicit none

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(3), x2(3), x3(3), x4(3), dx(3), inter(3)

integer i

call comp_dt(x, dx)         !  compute the first intermediate step
x1    = fract * deltat * dx
inter = x + x1 / 2.0

call comp_dt(inter, dx)     !  compute the second intermediate step
x2    = fract * deltat * dx
inter = x + x2 / 2.0

call comp_dt(inter, dx)     !  compute the third intermediate step
x3    = fract * deltat * dx
inter = x + x3

call comp_dt(inter, dx)     !  compute fourth intermediate step
x4 = fract * deltat * dx

!  compute new value for x

x = x + x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8

return
end subroutine adv_single_rk4


subroutine inv_linear_dt(x, dx, px)
!=====================================================================
!  subroutine inv_linear_dt(x, dx, px)
!
!  compute inv linear model lorenz time tendency (see notes 13mar94)
!  for now assumes stupid leap frog, will this be sufficient?

implicit none

real(r8), intent(in)  :: x(:), dx(:)
real(r8), intent(out) :: px(3)

real(r8) a(3, 3), fact, tdx(3)
integer i

a(1, 1) = -sigma * deltat + 1.0_r8
a(1, 2) =  sigma * deltat
a(1, 3) = 0.0_r8

a(2, 1) = (r - x(3)) * deltat
a(2, 2) = -1.0_r8 * deltat + 1.0_r8
a(2, 3) = -x(1) * deltat

a(3, 1) =  x(2) * deltat
a(3, 2) =  x(1) * deltat
a(3, 3) =    -b * deltat + 1.0_r8
   
!  initialize copy of dx

write(*, *) 'this routine is not up to date'

if(1 == 1) stop
!      tdx(i) = dx(i)
   tdx = dx

!  get rid of a(2, 3)

fact = a(2, 3) / a(3, 3)
   a(2, :) = a(2, :) - fact * a(3, :)
tdx(2) = tdx(2) - fact * tdx(3)

!  get rid of a(1, 2)

fact = a(1, 2) / a(2, 2)
   a(1, :) = a(1, :) - fact * a(2, :)
tdx(1) = tdx(1) - fact * tdx(2)

!  solve for the previous step linear perturbation

px(1) = tdx(1) / a(1, 1)
px(2) = (tdx(2) - a(2, 1) * px(1)) / a(2, 2)
px(3) = (tdx(3) - a(3, 1) * px(1) - a(3, 2) * px(2)) / a(3, 3)

end subroutine inv_linear_dt



subroutine linearize(nl, l)
!========================================================================
! subroutine linearize(nl, l)
!
! compute linear operator around state nl

implicit none

real(r8), intent(in)  :: nl(3)
real(r8), intent(out) :: l(3, 3)

l(1, 1) = -1.0_r8 * sigma * deltat + 1.0_r8
l(1, 2) =           sigma * deltat
l(1, 3) =          0.0_r8 * deltat

l(2, 1) = (r - nl(3)) * deltat
l(2, 2) = -1.0_r8 * deltat + 1.0_r8
l(2, 3) = -1.0_r8 * nl(1) * deltat

l(3, 1) =       nl(2) * deltat
l(3, 2) =       nl(1) * deltat
l(3, 3) = -1.0_r8 * b * deltat + 1.0_r8

return
end subroutine linearize




subroutine model_get_close_states(o_loc, radius, number, indices, dist)
!--------------------------------------------------------------------
! 
! Stub for computation of get close states

implicit none

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: number, indices(:)
real(r8), intent(out) :: dist(:)

! Because of F90 limits this stub must be here telling assim_model
! to do exhaustive search (number = -1 return)
number = -1

end subroutine model_get_close_states

   
!=========================================================================
! end module model_mod
!=========================================================================

end module model_mod
