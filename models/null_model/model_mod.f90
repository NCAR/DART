module model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! NOTE: This was created on 30 November, 2000 to work on evaluating 
! Ensemble Kalman Filter. It is almost identical to the version of the
! Lorenz 96 model in existence at this time except that it has extremely
! simple dynamics. The attractor has the value of 0.0 at all points at
! all times. Anything off the attractor grows exponentially.

use types_mod, only : r8

!use ncd_file_mod, only : init_ncd_file, def_axis, def_time_axis, &
!   init_diag_field, output_diag, sum_diag

use loc_and_dist_mod, only : loc_type, get_dist, set_loc

! TEMPORARY ADDITION OF RANDOM NOISE

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none
private

! Following is for repeatable random numbers

logical :: first_ens_seq = .true.
type (random_seq_type) :: ens_seq

public init_model, get_model_size, init_conditions, adv_1step, advance, &
   adv_true_state, output, diag_output_index, get_close_pts, state_loc, &
   model_output

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

integer,  parameter :: model_size = 40
real(r8), parameter ::    forcing = 8.00
real(r8), parameter ::    delta_t = 0.05

logical :: output_init = .FALSE.

! Define output indices for diagnostics
integer :: diag_output_index(9)

! Define the location of the state variables in module storage
type(loc_type) :: state_loc(model_size)

contains

!==================================================================



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


real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: dt(:)

integer :: j


if(first_ens_seq) then
   call init_random_seq(ens_seq)
   first_ens_seq = .false.
end if

do j = 1, model_size
   dt(j) = x(j)

! TEMPORARY ADDITION OF SOME NOISE
!!!   dt(j) = 1.0 * random_gaussian(ens_seq, 0.0_r8, 1.0_r8)

end do

end subroutine comp_dt



subroutine adv_1step(x)
!----------------------------------------------------------------------
! subroutine adv_1step(x)
!
! Does single time step advance for lorenz 96 model
! using four step rk time step

implicit none

real(r8), intent(inout) :: x(:)

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter
integer :: i

call comp_dt(x, dx)        !  Compute the first intermediate step
x1 = delta_t * dx
inter = x + x1 / 2.0_r8

call comp_dt(inter, dx)    !  Compute the second intermediate step
x2 = delta_t * dx
inter = x + x2 / 2.0_r8

call comp_dt(inter, dx)    !  Compute the third intermediate step
x3 = delta_t * dx
inter = x + x3

call comp_dt(inter, dx)    !  Compute fourth intermediate step
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

integer  :: i
real(r8) :: x_loc

! Define the interesting indexes for variables to do diag output; span lats
!do i = 1, 9
!   diag_output_index(i) = i
!end do
do i = 1, size(diag_output_index)
   diag_output_index(i) = model_size * &
                           (i - 1) / dble(size(diag_output_index)) + 1.0_r8
   write(*, *) 'output index ', i, diag_output_index(i)
!   diag_output_index(i) = 4*i
end do

! Define the locations of the state variables;

do i = 1, model_size
   x_loc = (i - 1.0_r8) / dble(model_size)
!   write(*, *) 'setting location ', i, x_loc
   call set_loc(state_loc(i), x_loc)
end do

x = 0.0_r8

end subroutine init_conditions



subroutine advance(x, num, xnew)
!----------------------------------------------------------------------
! subroutine advance(x, num, xnew)
!
! Advances the lorenz-96 model by a given number of steps
! Current state in x, new state in xnew, num time steps advanced

implicit none

real(r8),  intent(in) :: x(:)
integer,   intent(in) :: num
real(r8), intent(out) :: xnew(:)

integer :: i

xnew = x                 !  Copy initial conditions to avoid overwrite

do i = 1, num            !  Advance the appropriate number of steps
   call adv_1step(xnew)
end do

end subroutine advance



subroutine model_output(x, time)
!----------------------------------------------------------------------
! subroutine model_output(x, time)

implicit none

real(r8), intent(in) :: x(model_size)
real(r8), intent(in) :: time

real(r8) :: points(model_size)
integer  :: i

if(.NOT. output_init) then
   output_init = .TRUE.

! do netcdf setup for this field
!!!   call def_time_axis('time', 'seconds')
!!!   call init_ncd_file('one_d.nc', 'time', 'NetCDF file for l96 model')

!  Need to get points on axis

   do i = 1,  model_size
      points(i) = (i - 1.0_r8) / (1.0_r8 * model_size)
   end do

!!!   call def_axis('x', points, 'meters', 'x')
!!!   call init_diag_field('f', 'one_d.nc', (/'x'/), 'meters')

endif

! do netcdf output for this field
!!!call sum_diag('f', 'one_d.nc', x)
!!!call output_diag('one_d.nc', time)

end subroutine model_output



subroutine get_close_pts(list, num)
!----------------------------------------------------------------------
! subroutine get_close_pts(list, num)

implicit none

integer, intent(in)    :: num
integer, intent(inout) :: list(model_size, num)

integer :: i, j, offset, indx, temp

!do i = 1, model_size
!   do offset = -num/2, -num/2 + num - 1
!      indx = i + offset
!      if(indx > model_size) indx = indx - model_size
!      if(indx < 1) indx = model_size + indx
!      list(i, offset + num/2 + 1) = indx
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
         indx = i - j / 2
         if(indx < 1) indx = model_size + indx
         list(i, j) = indx
      else
         indx = i + j / 2
         if(indx > model_size) indx = indx - model_size
         list(i, j) = indx
      end if
   end do
end do

end subroutine get_close_pts



function get_model_size()
!----------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size





subroutine model_get_close_states(o_loc, radius, numb, indices, dist)
!--------------------------------------------------------------------
! 
! Stub for computation of get close states

implicit none

type(location_type), intent(in)  :: o_loc
real(r8),            intent(in)  :: radius
integer,             intent(out) :: numb, indices(:)
real(r8),            intent(out) :: dist(:)

! Because of F90 limits this stub must be here telling assim_model
! to do exhaustive search (numb = -1 return)
numb = -1

end subroutine model_get_close_states


!===================================================================
! end module model_mod
!===================================================================

end module model_mod
