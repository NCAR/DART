module obs_mod

use nag_wrap_mod, only : g05ddf_wrap
! Makes use of NAG, must be modified for NCAR

private
public num_obs, obs_var, take_obs, ens_ics
public num_obs_int, obs_var_int, take_obs_int

integer, parameter :: num_obs = 3, num_obs_int = 1

contains

!=======================================================================

function obs_var()

implicit none

double precision :: obs_var(num_obs)

obs_var = (2.0)**2



!obs_var = (10.0)**2



!obs_var(1) = (1000.0)**2
!obs_var(2) = (1000.0)**2
!obs_var(3) = (1.0)**2


! Following should duplicate Evensen with Variance 2.0 for his errors.
!obs_var = (2.0)

end function obs_var

!=======================================================================

function obs_var_int()

implicit none

double precision :: obs_var_int(num_obs_int)

!obs_var_int = (2.0)**2
obs_var_int = (10.0)**2

end function obs_var_int

!=======================================================================

function take_obs(x)

implicit none

double precision :: take_obs(num_obs)
double precision, intent(in) :: x(:)

! This is function h that maps state variables to obs
take_obs = x

!take_obs(1:2) = x(1:2)

!take_obs(1:3) = x
!take_obs(4:6) = x
!take_obs(7:9) = x
!take_obs(10:12) = x
!take_obs(13:15) = x
!take_obs(16:18) = x

!take_obs(1) = x(1) * x(2) / 10.0
!take_obs(2) = x(2) * x(3) / 10.0
!take_obs(3) = x(3) * x(1) / 10.0

!take_obs(1) = x(1) + x(2) + x(3)

end function take_obs

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

subroutine ens_ics(x, as)

implicit none

!  Get initial state for ensemble assimilation

double precision, intent(in) :: x(:)
double precision, intent(out) :: as(:, :)
integer :: i, j

do i = 1, size(x)
   do j = 1, size(as, 2)
! May want the constant to always be the same as the obs_var.
      as(i, j) = x(i) + 0.5 * g05ddf_wrap(dble(0.0), dble(1.0))
   end do
end do

end subroutine ens_ics

!========================================================================

end module obs_mod
