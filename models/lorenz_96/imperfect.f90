module model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! A version of the L96 model where the assimilating model doesn't know
! what the value of the forcing is in the true state model. Ends up
! being a 10 variable model; useful for testing model parameter
! estimation.

use ncd_file_mod, only : init_ncd_file, def_axis, def_time_axis, &
   init_diag_field, output_diag, sum_diag

use loc_and_dist_mod, only : loc_type, get_dist, set_loc

private

public init_model, get_model_size, init_conditions, adv_1step, advance, &
   adv_true_state, output, diag_output_index, state_loc, &
   model_output

integer, parameter :: model_size = 41
double precision, parameter :: true_forcing = 8.0
double precision, parameter :: delta_t = 0.05

logical :: output_init = .FALSE.

! Define output indices for diagnostics
integer :: diag_output_index(9)

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

! Computes the time tendency of the lorenz 1996 model given current state

double precision, intent(in) :: x(:)
double precision, intent(out) :: dt(:)

integer :: j, jp1, jm1, jm2

do j = 1, model_size - 1
   jp1 = j + 1
   if(jp1 > model_size - 1) jp1 = 1
   jm2 = j - 2
   if(jm2 < 1) jm2 = (model_size - 1) + jm2
   jm1 = j - 1
   if(jm1 < 1) jm1 = (model_size - 1)

! Forcing is in x(model_size)
   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + x(model_size)
end do

! No change in forcing value
dt(model_size) = 0.0

end subroutine comp_dt

!
!===================================================================
!
!
!  Does single time step advance for lorenz 96 model
!  using four step rk time step

subroutine adv_1step(x)

implicit none

double precision, intent(inout) :: x(:)
double precision, dimension(size(x)) :: x1, x2, x3, x4, dx, inter
integer i

!  Compute the first intermediate step
call comp_dt(x, dx)
x1(1:model_size - 1) = delta_t * dx(1:model_size - 1)
inter(1:model_size - 1) = x(1:model_size - 1) + x1(1:model_size - 1) / 2.0
inter(model_size) = x(model_size)
!  Compute the second intermediate step
call comp_dt(inter, dx)
x2(1:model_size - 1) = delta_t * dx(1:model_size - 1)
inter(1:model_size - 1) = x(1:model_size - 1) + x2(1:model_size - 1) / 2.0
inter(model_size) = x(model_size)
!  Compute the third intermediate step
call comp_dt(inter, dx)
x3(1:model_size - 1) = delta_t * dx(1:model_size - 1)
inter(1:model_size - 1) = x(1:model_size - 1) + x3(1:model_size - 1)
inter(model_size) = x(model_size)
!  Compute fourth intermediate step
call comp_dt(inter, dx)
x4(1:model_size - 1) = delta_t * dx(1:model_size - 1)
x4(model_size) = x(model_size)

!  Compute new value for x
x = x + x1/6.0 + x2/3.0 + x3/3.0 + x4/6.0
x(model_size) = x4(model_size)

end subroutine adv_1step

!
!----------------------------------------------------------------------
!

subroutine adv_true_state(x)

implicit none

double precision, intent(inout) :: x(:)

! Set forcing to true forcing value
x(model_size) = true_forcing

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

double precision, intent(out) :: x(:)
integer :: i
double precision :: x_loc

! Define the interesting indexes for variables to do diag output; span lats
do i = 2, 9
   diag_output_index(i) = i
end do
diag_output_index(1) = model_size

! Define the locations of the state variables;
do i = 1, model_size - 1
   x_loc = (i - 1.0) / (model_size - 1.)
!   write(*, *) 'setting location ', i, x_loc
   call set_loc(state_loc(i), x_loc)
end do

! Need a bogus location for the state variable, F, just set to loc 1 for now
x_loc = 0.0
call set_loc(state_loc(model_size), x_loc)

x = true_forcing
x(1) = 1.001 * true_forcing

end subroutine init_conditions

!
!====================================================================
!

!  Advances the lorenz-96 model by a given number of steps
!  Current state in x, new state in xnew, num time steps advanced

subroutine advance(x, num, xnew)

implicit none

double precision, intent(in) :: x(:)
double precision, intent(out) :: xnew(:)
integer, intent(in) :: num

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

subroutine model_output(x, time)

implicit none

real, intent(in) :: x(model_size - 1)
real, intent(in) :: time
real :: points(model_size - 1)
integer :: i

if(.NOT. output_init) then
   output_init = .TRUE.
! do netcdf setup for this field
   call def_time_axis('time', 'seconds')
   call init_ncd_file('one_d.nc', 'time', 'NetCDF file for l96 model')
!  Need to get points on axis
   do i = 1,  model_size - 1
      points(i) = (i - 1.0) / (1.0 * (model_size - 1.0))
   end do
   call def_axis('x', points, 'meters', 'x')
   call init_diag_field('f', 'one_d.nc', (/'x'/), 'meters')
endif

! do netcdf output for this field
call sum_diag('f', 'one_d.nc', x)
call output_diag('one_d.nc', time)

end subroutine model_output

!
!=======================================================================
!

subroutine get_close_pts(list, num)

implicit none

integer, intent(in) :: num
integer, intent(inout) :: list(model_size, num)

integer :: i, j, offset, index, temp

! Change in ordering of close points needed for debug in new basis_prod
! Change made 27 March, '00
do i = 1, model_size - 1
   list(i, 1) = i
   do j = 2, num
      if(j / 2 * 2 == j) then
         index = i - j / 2
         if(index < 1) index = model_size + index
         list(i, j) = index
      else
         index = i + j / 2
         if(index > model_size) index = index - model_size
         list(i, j) = index
      end if
   end do
end do

end subroutine get_close_pts

!!
!===================================================================
!

function get_model_size()

! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size

!
!===================================================================
!
end module model_mod
