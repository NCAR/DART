module model_mod

use ncd_file_mod, only : init_ncd_file, def_axis, def_time_axis, &
   init_diag_field, output_diag, sum_diag

!--------------------------------------------------------------
!define a type for locations
type location_type
   double precision :: x
end type location_type

!--------------------------------------------------------------

! Does dF/dt = -cDF/dx with forward in time, centered in space

private
public model_size, init_conditions, adv_1step, advance, output
public barot_to_dp, dp_to_barot
public delta_t, adv_true_state
public location_type, model_state_location, loc_dist, diag_output_index

integer, parameter :: reps = 1
integer, parameter :: points_per_rep = 20

integer, parameter :: model_size = reps * points_per_rep

! define model parameters
! Domain is 0 to 1, c is 1, deltat is selected to make Courant = 0.5
double precision, parameter :: delta_x = 1.0 / points_per_rep
double precision, parameter :: delta_t = delta_x / 4.0

logical :: output_init = .FALSE.

double precision :: true_state_time = 0.0

! Define output indices for diagnostics
integer :: diag_output_index(9)

contains

!!==================================================================

function model_state_location()

implicit none

type (location_type) :: model_state_location(model_size)

integer :: i

do i = 1, model_size
   model_state_location(i)%x = i * delta_x
end do

end function model_state_location

!====================================================================

function loc_dist(a, b)

implicit none

double precision :: loc_dist
type (location_type), intent(in) :: a, b

double precision :: av, bv, width

av = a%x
bv = b%x

loc_dist = dabs(av - bv)
width = delta_x * model_size
if(loc_dist > width / 2) loc_dist = width - loc_dist

end function loc_dist

!!===================================================================
!
!  DOES SINGLE TIME STEP ADVANCE FOR WAVE EQUATION
!  USING FORWARD IN SPACE

   SUBROUTINE ADV_1STEP(X)

   IMPLICIT NONE

   DOUBLE PRECISION, intent(inout) :: X(:)

   x = x - (delta_t / (2. * delta_x)) *(cshift(x, 1) - cshift(x, -1))

   RETURN
   END  subroutine adv_1step

!
!====================================================================
!

!  ADVANCES THE MODEL BY A GIVEN NUMBER OF STEPS
!  CURRENT STATE IN X, NEW STATE IN XNEW, NUM TIME STEPS ADVANCED

subroutine advance(x, num, xnew)

   IMPLICIT NONE

   double precision, intent(in) :: x(:)
   double precision, intent(out) :: xnew(:)
   integer, intent(in) :: num

   integer i

!  COPY INITIAL CONDITIONS TO AVOID OVERWRITE
xnew = x

!  ADVANCE THE APPROPRIATE NUMBER OF STEPS
do i = 1, num
   call adv_1step(xnew)
end do

return
end subroutine advance
!
!===================================================================
!

   subroutine init_conditions(x)

!  INITIAL CONDITIONS ARE DELTA FUNCTION IN MIDDLE

   implicit none

   double precision, INTENT(OUT) :: x(:)
   integer i

! Define the interesting indexes for variables to do diag output; span lats
do i = 1, 9
   diag_output_index(i) = i
end do

   true_state_time = 0.0

   do i = 1, model_size
      x(i) = cont_ic(i / (points_per_rep * dble(1.0)))
!      x(i) = cont_ic(i / (model_size * dble(1.0)))
   end do

   end subroutine init_conditions

!
!===================================================================
!

function single_true_state(x, t)

implicit none

double precision, intent(in) :: x, t
double precision :: single_true_state
double precision :: time

! Need G(x-ct) = G(x - t) here
time = x - t
!time = time - floor(time)
time = modulo(time, dble(1.0))
single_true_state = cont_ic(time)

end function single_true_state

!
!====================================================================
!

subroutine adv_true_state(x)

implicit none

double precision, intent(out) :: x(:)
integer :: i

! Advance time by deltat
true_state_time = true_state_time + delta_t

do i = 1, model_size
   x(i) = single_true_state(i / (points_per_rep * dble(1.0)), true_state_time)
!   x(i) = single_true_state(i / (model_size * dble(1.0)), true_state_time)
end do

end subroutine adv_true_state

!
!====================================================================
!

function cont_ic(x)

implicit none

double precision, intent(in) :: x
double precision :: cont_ic

double precision :: frac_x

frac_x = x - int(x)

if(frac_x <= 0.5) then
   cont_ic = frac_x
else
   cont_ic = 1.0 - frac_x 
endif

end function cont_ic

!
!====================================================================
!

subroutine output(x, time)

implicit none

real, intent(in) :: x(model_size)
real, intent(in) :: time
real :: points(model_size)
integer :: i

if(.NOT. output_init) then
   output_init = .TRUE.
! do netcdf setup for this field
   call def_time_axis('time', 'seconds')
   call init_ncd_file('one_d.nc', 'time', 'NetCDF file for simple one_d model')
!  Need to get points on axis
   do i = 1,  model_size
      points(i) = delta_x * (i - 1.0)
   end do
   call def_axis('x', points, 'meters', 'x')
   call init_diag_field('f', 'one_d.nc', (/'x'/), 'meters')
endif

! do netcdf output for this field
call sum_diag('f', 'one_d.nc', x)
call output_diag('one_d.nc', time)

end subroutine output

!------------------------------------------------------------------------


function dp_to_barot(x)

implicit none

! Converts from double precision physical space array to barotropic
! spherical harmonic stream function

double precision, intent(in) :: x(:)
complex :: dp_to_barot(1, 1)

dp_to_barot = 0.0

end function dp_to_barot

!-------------------------------------------------------------------------
function barot_to_dp(psi)

implicit none

! Converts from spherical harmonic stream function to double precision
! physical space array.

complex, intent(in) :: psi(:, :)
double precision :: barot_to_dp(1)

barot_to_dp = 0.0

end function barot_to_dp

!-------------------------------------------------------------------------


end module model_mod
