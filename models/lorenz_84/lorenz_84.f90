module model_mod

! Implements Lorenz 84 3 variable model with intermediate level
! attractor.


integer :: model_size = 3

!  define model parameters
double precision, parameter :: a = 0.25, b= 4., f = 8., g = 1.25
double precision, parameter :: deltat = 0.01 

contains

!
!------------------------------------------------------------------
!

subroutine comp_dt(x, dt)

!  compute time tendency of lorenz-84 3-variable model given current state

implicit none

double precision, intent(in) :: x(:)
double precision, intent(out) :: dt(:)


!  compute lorenz 84 model time tendency
dt(1) = -x(2)**2 -x(3)**2 - a*x(1) + a*f
dt(2) = x(1)*x(2) - b*x(1)*x(3) - x(2) + g
dt(3) = b*x(1)*x(2) + x(1)*x(3) - x(3)

end subroutine comp_dt
!
!------------------------------------------------------------------
!

!  advances the 3 variable  model by a given number of steps

subroutine advance(x, num, xnew)

implicit none

double precision, intent(in) :: x(3)
double precision, intent(out) :: xnew(3)
integer, intent(in) :: num

integer i

!  copy initial conditions to avoid overwrite
xnew = x
   
!  advance the appropriate number of steps
do i = 1, num
   call adv_1step(xnew)
end do

end subroutine advance
!
!-------------------------------------------------------------------
!
subroutine init_conditions(x)

!  defines off-attractor initial conditions for lorenz-84

implicit none

double precision, intent(out) :: x(:)

x(1) = 0.0
x(2) = 1.0
x(3) = 0.0

end subroutine init_conditions
!
!--------------------------------------------------------------------
!

!  does single time step advance for lorenz 3 variable model using two step rk
subroutine adv_1step(x)

implicit none

double precision, intent(inout) :: x(:)
double precision :: fract = 1.0

call adv_single(x, fract)

end subroutine adv_1step
!
!--------------------------------------------------------------------
!

!  does single time step advance for lorenz 3 variable model using two step rk
subroutine adv_single(x, fract)

implicit none

double precision, intent(inout) :: x(:)
double precision, intent(in) :: fract
double precision :: x1(3), x2(3), dx(3)
integer i

!  compute the first intermediate step
call comp_dt(x, dx)
x1 = x + fract * deltat * dx
!  compute the second intermediate step
call comp_dt(x1, dx)
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate
x = (x + x2) / 2.0

end subroutine adv_single
!
!-------------------------------------------------------------------
!

        subroutine linearize(nl, l)
        implicit none

!  COMPUTE LINEAR OPERATOR AROUND STATE NL
        DOUBLE PRECISION NL(3), L(3, 3)

        L(1, 1) = -1.0 * A * DELTAT + 1.0
        L(1, 2) = -1.0 * NL(2) * DELTAT
        L(1, 3) = -1.0 * NL(3) * DELTAT
        L(2, 1) = (NL(2) - B * NL(3)) * DELTAT
        L(2, 2) = (-1.0 + NL(1)) * DELTAT + 1.0
        L(2, 3) = -1.0 * B * NL(1) * DELTAT
        L(3, 1) = (B * NL(2) + NL(3)) * DELTAT
        L(3, 2) = B * NL(1) * DELTAT
        L(3, 3) = (NL(1) - 1.0) * DELTAT + 1.0
        RETURN
        END subroutine linearize

!--------------------------------------------------------------------

end module model_mod
