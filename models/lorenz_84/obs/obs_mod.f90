module obs_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use types_mod
use new_random_mod, only: gasdev, ran1
use sort_mod

integer, parameter :: num_obs = 3


contains



  function obs_sd()
!-----------------------------------------------------------------------
! function obs_sd()

implicit none

real(r8) :: obs_sd(num_obs)

obs_sd = 0.2_r8
!obs_sd(1) = 2.0_r8
!obs_sd(2) = 2.0_r8
!obs_sd(3) = 2.0_r8

end function obs_sd



  function take_obs(x)
!-----------------------------------------------------------------------
! function take_obs(x)
!
! This is function h that maps state variables to obs

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: take_obs(num_obs)

take_obs = x

!take_obs(1) = x(3)
!take_obs(1) = x(1) * x(2)
!take_obs(2) = x(2) * x(3)
!take_obs(3) = x(3) * x(1)

!take_obs(1) = x(1) + x(2) + x(3)

end function take_obs



  subroutine ens_ics(x, as)
!------------------------------------------------------------------------
! subroutine ens_ics(x, as)
!
! Get initial state for ensemble assimilation

implicit none

real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: as(:, :)

integer :: i, j

do i = 1, size(x)
   do j = 1, size(as, 2)
      as(i, j) = gasdev() * 0.02_r8 + x(i)
   end do
end do

end subroutine ens_ics

!------------------------------------------------------------------------
! End of obs_mod.f90
!------------------------------------------------------------------------

end module obs_mod
