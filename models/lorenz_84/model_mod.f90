module model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Implements Lorenz 84 3 variable model with intermediate level
! attractor.


use types_mod

integer :: model_size = 3

!  define model parameters
real(r8), parameter :: a = 0.25_r8, b= 4.0_r8, f = 8.0_r8, g = 1.25_r8
real(r8), parameter :: deltat = 0.01_r8

contains



  subroutine comp_dt(x, dt)
!------------------------------------------------------------------
! subroutine comp_dt(x, dt)
!
! compute time tendency of lorenz-84 3-variable model given current state

implicit none

real(r8), intent(in) :: x(:)
real(r8), intent(out) :: dt(:)

!  compute lorenz 84 model time tendency

dt(1) = -x(2)**2 -x(3)**2 - a*x(1) + a*f
dt(2) = x(1)*x(2) - b*x(1)*x(3) - x(2) + g
dt(3) = b*x(1)*x(2) + x(1)*x(3) - x(3)

end subroutine comp_dt



subroutine advance(x, num, xnew)
!------------------------------------------------------------------
! subroutine advance(x, num, xnew)
!
! advances the 3 variable  model by a given number of steps

implicit none

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(3)
real(r8), intent(out) :: xnew(3)

integer :: i

xnew = x        !  copy initial conditions to avoid overwrite
   
!  advance the appropriate number of steps

do i = 1, num
   call adv_1step(xnew)
end do

end subroutine advance



  subroutine init_conditions(x)
!-------------------------------------------------------------------
! subroutine init_conditions(x)
!
!  defines off-attractor initial conditions for lorenz-84

implicit none

real(r8), intent(out) :: x(:)

x(1) = 0.0_r8
x(2) = 1.0_r8
x(3) = 0.0_r8

end subroutine init_conditions



  subroutine adv_1step(x)
!--------------------------------------------------------------------
! subroutine adv_1step(x)
!
! does single time step advance for lorenz 3 variable model using 
! two step rk

implicit none

real(r8), intent(inout) :: x(:)

real(r8) :: fract = 1.0_r8

call adv_single(x, fract)

end subroutine adv_1step



subroutine adv_single(x, fract)
!--------------------------------------------------------------------
!
!  does single time step advance for lorenz 3 variable model using two step rk

implicit none

real(r8), intent(in)    :: fract
real(r8), intent(inout) :: x(:)

real(r8) :: x1(3), x2(3), dx(3)
integer  :: i

!  compute the first intermediate step

call comp_dt(x, dx)
x1 = x + fract * deltat * dx

!  compute the second intermediate step

call comp_dt(x1, dx)
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

end subroutine adv_single



  subroutine linearize(nl, l)
!-------------------------------------------------------------------
! subroutine linearize(nl, l)
!
!  compute linear operator around state nl

implicit none

real(r8) :: nl(3), l(3, 3)

l(1, 1) = -1.0_r8 *     a     * deltat + 1.0_r8
l(1, 2) = -1.0_r8 * nl(2)     * deltat
l(1, 3) = -1.0_r8 * nl(3)     * deltat
l(2, 1) = (nl(2) - b * nl(3)) * deltat
l(2, 2) = (-1.0_r8 + nl(1))   * deltat + 1.0_r8
l(2, 3) = -1.0_r8 * b * nl(1) * deltat
l(3, 1) = (b * nl(2) + nl(3)) * deltat
l(3, 2) =  b * nl(1)          * deltat
l(3, 3) = (nl(1) - 1.0_r8)    * deltat + 1.0_r8

return
end subroutine linearize

!--------------------------------------------------------------------
! End of module model_mod.f90
!--------------------------------------------------------------------

end module model_mod
