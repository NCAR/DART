module obs_mod

use new_random_mod, only: gasdev, ran1
use sort_mod

integer, parameter :: num_obs = 3


contains

!-----------------------------------------------------------------------

function obs_sd()

implicit none

double precision :: obs_sd(num_obs)

obs_sd = 0.2
!obs_sd(1) = 2.0
!obs_sd(2) = 2.0
!obs_sd(3) = 2.0

end function obs_sd

!-----------------------------------------------------------------------

function take_obs(x)

implicit none

double precision :: take_obs(num_obs)
double precision, intent(in) :: x(:)

! This is function h that maps state variables to obs
take_obs = x

!take_obs(1) = x(3)

!take_obs(1) = x(1) * x(2)
!take_obs(2) = x(2) * x(3)
!take_obs(3) = x(3) * x(1)

!take_obs(1) = x(1) + x(2) + x(3)

end function take_obs

!------------------------------------------------------------------------

subroutine ens_ics(x, as)

implicit none

!  Get initial state for ensemble assimilation

double precision, intent(in) :: x(:)
double precision, intent(out) :: as(:, :)
integer :: i, j

do i = 1, size(x)
   do j = 1, size(as, 2)
      as(i, j) = dble(gasdev()) * 0.02 + x(i)
   end do
end do

end subroutine ens_ics

!------------------------------------------------------------------------

end module obs_mod
