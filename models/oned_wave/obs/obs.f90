module obs_mod

use nag_wrap_mod, only : g05ddf_wrap
use model_mod, only : model_size, location_type, model_state_location
   

private
public num_obs, obs_var, take_obs, ens_ics
public num_obs_int, obs_var_int, take_obs_int, obs_location

integer, parameter :: num_obs = model_size
integer, parameter :: num_obs_int = model_size

logical, parameter :: init_obs_locations = .true.
type (location_type) :: location(num_obs)

contains

!=======================================================================

subroutine init_obs

implicit none

type (location_type) :: model_loc(model_size)

integer :: i

! First get model grid locations
model_loc = model_state_location()

do i = 1, model_size - 1
   location(i)%x = (model_loc(i)%x + model_loc(i + 1)%x) / 2.0
end do

location(model_size)%x = model_loc(model_size)%x + &
   0.5 * (model_loc(model_size - 1)%x - model_loc(model_size - 2)%x) 

do i = 1, model_size
   write(*, *) 'in init obs ', i, real(model_loc(i)%x)
end do

end subroutine init_obs

!=======================================================================

function obs_location()

implicit none

type (location_type) :: obs_location(num_obs)

integer :: i

if(init_obs_locations) call init_obs

obs_location = location

end function obs_location

!=======================================================================

function obs_var()

implicit none

double precision :: obs_var(num_obs)

obs_var = 0.1**2

end function obs_var

!=======================================================================

function take_obs(x)

implicit none

double precision :: take_obs(num_obs)
double precision, intent(in) :: x(:)
integer :: i

!=====================
!take_obs = x
!======================

do i = 1, num_obs - 1
   take_obs(i) = (x(i) + x(i + 1))
end do
take_obs(num_obs) = (x(1) + x(num_obs))

!=======================
! This is function h that maps state variables to obs
!do i = 1, num_obs
!   take_obs(i) = x(4*i - 1) + x(4*i)
!   take_obs(i) = x(2*i - 1) + x(2*i)
!end do
!============================

!take_obs = x**2

!===========================
!do i = 1, num_obs
!   take_obs(i) = x(2*i - 1) ** 2
!end do
!===========================

end function take_obs

!========================================================================

subroutine ens_ics(x, as)

implicit none

!  Get initial state for ensemble assimilation

double precision, intent(in) :: x(:)
double precision, intent(out) :: as(:, :)
integer :: i, j

do i = 1, size(x)
   do j = 1, size(as, 2)
! May want the constant to always be the same as the obs_var.
      as(i, j) = x(i) + 0.1 * g05ddf_wrap(dble(0.0), dble(1.0))
   end do
end do

end subroutine ens_ics

!=======================================================================

function obs_var_int()

implicit none

double precision :: obs_var_int(num_obs_int)

!obs_var_int = (2.0)**2
obs_var_int = (10.0)**2

end function obs_var_int

!========================================================================

function take_obs_int(x)

implicit none

double precision take_obs_int(num_obs_int)
double precision, intent(in) :: x(:)


! This is function h for intermediate non-convolution time obs
take_obs_int(1) = x(1) * x(2) + x(2) * x(3)
!take_obs_int = x

end function take_obs_int

!========================================================================

end module obs_mod
