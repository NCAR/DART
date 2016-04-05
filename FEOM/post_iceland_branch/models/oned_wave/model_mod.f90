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

! Does dF/dt = -cDF/dx with forward in time, centered in space

use    types_mod, only : r8
use ncd_file_mod, only : init_ncd_file, def_axis, def_time_axis, &
                         init_diag_field, output_diag, sum_diag

implicit none
private

type location_type
   private
   real(r8) :: x
end type location_type

public :: model_size, init_conditions, adv_1step, advance, output
public :: barot_to_dp, dp_to_barot
public :: delta_t, adv_true_state
public :: location_type, model_state_location, loc_dist, diag_output_index
public :: model_get_close_state

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

integer, parameter :: reps = 1
integer, parameter :: points_per_rep = 20
integer, parameter :: model_size = reps * points_per_rep

! define model parameters
! Domain is 0 to 1, c is 1, deltat is selected to make Courant = 0.5

real(r8), parameter :: delta_x = 1.0_r8 / (points_per_rep * 1.0_r8)
real(r8), parameter :: delta_t = delta_x / 4.0_r8

logical  :: output_init = .FALSE.
real(r8) :: true_state_time = 0.0_r8

integer  :: diag_output_index(9)       ! Define output indices for diagnostics

contains



  function model_state_location()
!-------------------------------------------------------------------------
! function model_state_location()

implicit none

type (location_type) :: model_state_location(model_size)

integer :: i

do i = 1, model_size
   model_state_location(i)%x = i * delta_x
end do

end function model_state_location



  function loc_dist(a, b)
!-------------------------------------------------------------------------
! function loc_dist(a, b)

implicit none

type (location_type), intent(in) :: a, b
real(r8)                         :: loc_dist

real(r8) :: av, bv, width

av = a%x
bv = b%x

loc_dist = abs(av - bv)
width    = delta_x * model_size
if(loc_dist > width / 2.0_r8) loc_dist = width - loc_dist

end function loc_dist



  subroutine adv_1step(x)
!-------------------------------------------------------------------------
! subroutine adv_1step(x)
!
! does single time step advance for wave equation
! using forward in space
!

implicit none

real(r8), intent(inout) :: x(:)

x = x - (delta_t / (2.0_r8 * delta_x)) *(cshift(x, 1) - cshift(x, -1))

return
end  subroutine adv_1step



  subroutine advance(x, num, xnew)
!-------------------------------------------------------------------------
! subroutine advance(x, num, xnew)
! advances the model by a given number of steps
! current state in x, new state in xnew, num time steps advanced
!

implicit none

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: xnew(:)

integer :: i

xnew = x                 !  copy initial conditions to avoid overwrite

do i = 1, num            !  advance the appropriate number of steps
   call adv_1step(xnew)
end do

return
end subroutine advance



  subroutine init_conditions(x)
!-------------------------------------------------------------------------
! subroutine init_conditions(x)
!
! initial conditions are delta function in middle

implicit none

real(r8), intent(out) :: x(:)
integer :: i

! define the interesting indexes for variables to do diag output; span lats

do i = 1, 9
   diag_output_index(i) = i
end do

true_state_time = 0.0_r8

do i = 1, model_size
   x(i) = cont_ic(i / (points_per_rep * 1.0_r8))
!  x(i) = cont_ic(i / (model_size     * 1.0_r8))
end do

end subroutine init_conditions



  function single_true_state(x, t)
!-------------------------------------------------------------------------
! function single_true_state(x, t)

implicit none

real(r8), intent(in) :: x, t
real(r8)             :: single_true_state

real(r8) :: time

! need g(x-ct) = g(x - t) here
time = x - t

!time = time - floor(time)

time              = modulo(time, 1.0_r8)
single_true_state = cont_ic(time)

end function single_true_state



  subroutine adv_true_state(x)
!-------------------------------------------------------------------------
! subroutine adv_true_state(x)
!

implicit none

real(r8), intent(out) :: x(:)

integer :: i

true_state_time = true_state_time + delta_t      ! advance time

do i = 1, model_size
   x(i) = single_true_state(i / (points_per_rep * 1.0_r8), true_state_time)
!  x(i) = single_true_state(i / (model_size     * 1.0_r8), true_state_time)
end do

end subroutine adv_true_state



  function cont_ic(x)
!-------------------------------------------------------------------------
! function cont_ic(x)

implicit none

real(r8), intent(in) :: x
real(r8)             :: cont_ic

real(r8) :: frac_x

frac_x = x - int(x)

if(frac_x <= 0.5_r8) then
   cont_ic = frac_x
else
   cont_ic = 1.0_r8 - frac_x 
endif

end function cont_ic



  subroutine output(x, time)
!-------------------------------------------------------------------------
! subroutine output(x, time)

implicit none

real, intent(in) :: x(model_size)
real, intent(in) :: time

real    :: points(model_size)
integer :: i

if(.NOT. output_init) then
   output_init = .TRUE.

   ! do netcdf setup for this field

   call def_time_axis('time', 'seconds')
   call init_ncd_file('one_d.nc', 'time', 'NetCDF file for simple one_d model')

   !  Need to get points on axis

   do i = 1,  model_size
      points(i) = delta_x * (i - 1.0_r8)
   end do
   call def_axis('x', points, 'meters', 'x')
   call init_diag_field('f', 'one_d.nc', (/'x'/), 'meters')
endif

! do netcdf output for this field

call sum_diag('f', 'one_d.nc', x)
call output_diag('one_d.nc', time)

end subroutine output



  function dp_to_barot(x)
!-------------------------------------------------------------------------
! function dp_to_barot(x)
!
! Converts from real(r8) physical space array to barotropic
! spherical harmonic stream function

implicit none

real(r8), intent(in) :: x(:)
complex              :: dp_to_barot(1, 1)

dp_to_barot = 0.0

end function dp_to_barot



  function barot_to_dp(psi)
!-------------------------------------------------------------------------
! function barot_to_dp(psi)
!
! Converts from spherical harmonic stream function to real(r8)
! physical space array.

implicit none

complex, intent(in) :: psi(:, :)
real(r8)            :: barot_to_dp(1)

barot_to_dp = 0.0_r8

end function barot_to_dp



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


!-------------------------------------------------------------------------
! End of model_mod.f90
!-------------------------------------------------------------------------

end module model_mod
