module model_mod

!use ncd_file_mod, only : init_ncd_file, def_axis, def_time_axis, &
!   init_diag_field, output_diag, sum_diag

use types_mod
use loc_and_dist_mod, only : loc_type, get_dist, set_loc

private

public init_model, get_model_size, init_conditions, adv_1step, advance, &
   adv_true_state, output, diag_output_index, get_close_pts, state_loc, &
   model_output

integer(i4), parameter :: model_size =   40_i4
real(r8),    parameter ::    forcing = 8.00_r8
real(r8),    parameter ::    delta_t = 0.05_r8

logical :: output_init = .FALSE.

! Define output indices for diagnostics

integer(i4) :: diag_output_index(9)

! Define the location of the state variables in module storage

type(loc_type) :: state_loc(model_size)

contains

!======================================================================



subroutine init_model()
!----------------------------------------------------------------------
! subroutine init_model()
!
! Stub for model initialization, not needed for L96

end subroutine init_model



subroutine comp_dt(x, dt)
!----------------------------------------------------------------------
! subroutine comp_dt(x, dt)
! 
! Computes the time tendency of the lorenz 1996 model given current state

implicit none

real(r8), intent( in) :: x(:)
real(r8), intent(out) :: dt(:)

integer(i4) :: j, jp1, jm1, jm2

do j = 1, model_size
   jp1 = j + 1
   if(jp1 > model_size) jp1 = 1
   jm2 = j - 2
   if(jm2 < 1) jm2 = model_size + jm2
   jm1 = j - 1
   if(jm1 < 1) jm1 = model_size
   
   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + forcing
end do

end subroutine comp_dt



subroutine adv_1step(x)
!----------------------------------------------------------------------
! subroutine adv_1step(x)
!
! Does single time step advance for lorenz 96 model
! using four-step rk time step

implicit none

real(r8), intent(inout) :: x(:)
real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter
integer(i4) :: i

!  Compute the first intermediate step

call comp_dt(x, dx)
x1    = delta_t * dx
inter = x + x1 / 2.0_r8

!  Compute the second intermediate step

call comp_dt(inter, dx)
x2    = delta_t * dx
inter = x + x2 / 2.0_r8

!  Compute the third intermediate step

call comp_dt(inter, dx)
x3    = delta_t * dx
inter = x + x3

!  Compute fourth intermediate step

call comp_dt(inter, dx)
x4 = delta_t * dx

!  Compute new value for x

x = x + x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8

end subroutine adv_1step



subroutine adv_true_state(x)
!----------------------------------------------------------------------
! subroutine adv_true_state(x)
!

implicit none

real(r8), intent(inout) :: x(:)

call adv_1step(x)

end subroutine adv_true_state



subroutine init_conditions(x)
!----------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for lorenz 96
! It is assumed that this is called before any other routines in this
! module. Should probably make that more formal and perhaps enforce for
! more comprehensive models.

implicit none

real(r8), intent(out) :: x(:)

integer(i4) :: i
real(r8)    :: x_loc

! Define the interesting indexes for variables to do diag output; span lats

do i = 1, size(diag_output_index)
   diag_output_index(i) = i
end do

!do i = 1, 9
!   diag_output_index(i) = model_size * (i - 1) / 9.0 + 1.0
!   write(*, *) 'output index ', i, diag_output_index(i)
!   diag_output_index(i) = 4*i
!end do

! Define the locations of the state variables;
do i = 1, model_size
   x_loc = (i - 1.0) / model_size
!   write(*, *) 'setting location ', i, x_loc
   call set_loc(state_loc(i), x_loc)
end do

x    = forcing
x(1) = 1.001_r8 * forcing

end subroutine init_conditions



subroutine advance(x, num, xnew)
!----------------------------------------------------------------------
! subroutine advance(x, num, xnew)
!
! Advances the lorenz-96 model by a given number of steps
! Current state in x, new state in xnew, num time steps advanced

implicit none

real(r8),  intent( in) :: x(:)
real(r8),  intent(out) :: xnew(:)
integer(i4),intent( in) :: num

integer(i4) :: i

xnew = x      ! Copy initial conditions to avoid overwrite

!  Advance the appropriate number of steps
do i = 1, num
   call adv_1step(xnew)
end do

end subroutine advance



!subroutine model_output(x, time)
!---------------------------------------------------------------------
! subroutine model_output(x, time)
!

!implicit none

!real, intent(in) :: x(model_size)
!real, intent(in) :: time
!real :: points(model_size)
!integer(i4) :: i

!if(.NOT. output_init) then
!   output_init = .TRUE.
!! do netcdf setup for this field
!   call def_time_axis('time', 'seconds')
!   call init_ncd_file('one_d.nc', 'time', 'NetCDF file for l96 model')
!!  Need to get points on axis
!   do i = 1,  model_size
!      points(i) = (i - 1.0) / (1.0 * model_size)
!   end do
!   call def_axis('x', points, 'meters', 'x')
!   call init_diag_field('f', 'one_d.nc', (/'x'/), 'meters')
!endif

!! do netcdf output for this field
!call sum_diag('f', 'one_d.nc', x)
!call output_diag('one_d.nc', time)

!end subroutine model_output



subroutine get_close_pts(list, num)
!---------------------------------------------------------------------
! subroutine get_close_pts(list, num)
!

implicit none

integer(i4), intent(   in) :: num
integer(i4), intent(inout) :: list(model_size, num)

integer(i4) :: i, j, offset, index, temp

!do i = 1, model_size
!   do offset = -num/2, -num/2 + num - 1
!      index = i + offset
!      if(index > model_size) index = index - model_size
!      if(index < 1) index = model_size + index
!      list(i, offset + num/2 + 1) = index
!   end do
! Always need the actual point first in list
!   temp = list(i, 1)
!   list(i, 1) =  list(i, num / 2 + 1)
!   list(i, num / 2 + 1) = temp
!end do

! Change in ordering of close points needed for debug in new basis_prod
! Change made 27 March, '00

do i = 1, model_size
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



function get_model_size()
!---------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer(i4) :: get_model_size

get_model_size = model_size

end function get_model_size

!
!===================================================================
! End of model_mod
!===================================================================
!
end module model_mod
