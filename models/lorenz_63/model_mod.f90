module model_mod

use types_mod
use loc_and_dist_mod, only : loc_type, get_dist, set_loc

private
public get_model_size, init_conditions, adv_1step, advance, model_output, &
   diag_output_index, init_model, adv_true_state, state_loc, output

integer, parameter :: model_size = 3


! Define output indices for diagnostics

integer :: diag_output_index(3)

! Define the location of the state variables in module storage

type(loc_type) :: state_loc(model_size)

!  define model parameters

real(r8), parameter ::  sigma = 10.0_r8
real(r8), parameter ::      r = 28.0_r8
real(r8), parameter ::      b = 8.0_r8 / 3.0_r8
real(r8), parameter :: deltat = 0.01_r8

contains

!==================================================================



subroutine init_model()
!==================================================================
! subroutine init_model()
! Stub for model initialization, not needed for L63

end subroutine init_model




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



subroutine advance(x, num, xnew)
!===================================================================
! subroutine advance(x, num, xnew)
!
! advances the 3 variable lorenz-63 model by a given number of steps
! current state in x, new state in xnew, num time steps advanced

implicit none

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: num
real(r8), intent(out) :: xnew(:)

integer :: i

xnew = x                  !  copy initial conditions to avoid overwrite
   
do i = 1, num             !  advance the appropriate number of steps
   call adv_1step(xnew)
end do

return
end subroutine advance


subroutine adv_true_state(x)
!===================================================================
! subroutine adv_true_state(x)
!

implicit none

real(r8), intent(inout) :: x(:)

call adv_1step(x)

end subroutine adv_true_state



subroutine init_conditions(x)
!===================================================================
! subroutine init_conditions(x)
!
!  off-attractor initial conditions for lorenz 63

implicit none

real(r8), intent(out) :: x(:)

integer  :: i
real(r8) :: x_loc

! Define the interesting indexes for variables to do diag output

do i = 1, size(diag_output_index)
   diag_output_index(i) = i
end do

! Define the locations of the state variables; all co-located for L63

do i = 1, model_size
   x_loc = 1.0_r8
   call set_loc(state_loc(i), x_loc)
end do

! Initial conditions that move nicely onto attractor

x = 0.10_r8

return
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

return
end subroutine linear_dt




subroutine adv_1step(x)
!====================================================================
! subroutine adv_1step(x)
!
! does single time step advance for lorenz convective 3 variable model
! using two step rk time step

implicit none

real(r8), intent(inout) :: x(:)
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

   
!=========================================================================
! end module model_mod
!=========================================================================

end module model_mod
