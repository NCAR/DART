module obs_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use types_mod

!use nag_wrap_mod, only : g05ddf_wrap
! Makes use of NAG, must be modified for NCAR
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use loc_and_dist_mod, only : loc_type, set_loc, get_loc

private
public num_obs, obs_var, take_obs, ens_ics, take_single_obs, init_obs, &
   get_close_state, obs_loc
public num_obs_int, obs_var_int, take_obs_int

! Following is for repeatable random numbers

logical :: first_ens_seq = .true.
type (random_seq_type) :: ens_seq

! Array of structure for observation locations; static with time in this version
type(loc_type) :: obs_loc(3)

integer, parameter :: num_obs = 3, num_obs_int = 1

contains

!=======================================================================



subroutine init_obs()
!-----------------------------------------------------------------------
! subroutine init_obs()
!
! Stub for init_obs

implicit none

integer  :: i
real(r8) :: loc


! Make all observations co-located with all state variables (no ramping)

do i = 1, num_obs
   loc = 1.0_r8
   call set_loc(obs_loc(i), loc)
end do

end subroutine init_obs





function obs_var()
!-----------------------------------------------------------------------
! function obs_var()
!

implicit none

real(r8) :: obs_var(num_obs)

obs_var = (2.0)**2

!obs_var = (10.0)**2

!obs_var(1) = (1000.0)**2
!obs_var(2) = (1000.0)**2
!obs_var(3) = (1.0)**2

! Following should duplicate Evensen with Variance 2.0 for his errors.
!obs_var = (2.0)

end function obs_var



function obs_var_int()
!-----------------------------------------------------------------------
! function obs_var_int()
!

implicit none

real(r8) :: obs_var_int(num_obs_int)

!obs_var_int = (2.0)**2
obs_var_int = (10.0)**2

end function obs_var_int



function take_obs(x)
!-----------------------------------------------------------------------
! function take_obs(x)
!

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: take_obs(num_obs)

take_obs = x      ! This is function h that maps state variables to obs

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



function take_single_obs(x, index)
!-----------------------------------------------------------------------
! function take_single_obs(x, index)
!

implicit none

real(r8), intent(in) :: x(:)
integer,  intent(in) :: index
real(r8)             :: take_single_obs

take_single_obs = x(index)

end function take_single_obs



function take_obs_int(x)
!-----------------------------------------------------------------------
! function take_obs_int(x)
!

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: take_obs_int(num_obs_int)


! This is function h for intermediate non-convolution time obs

take_obs_int(1) = x(1) * x(2) + x(2) * x(3)

! take_obs_int = x

end function take_obs_int



subroutine ens_ics(x, as)
!-----------------------------------------------------------------------
! subroutine ens_ics(x, as)
!

implicit none

!  Get initial state for ensemble assimilation

real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: as(:, :)

integer :: i, j

if (first_ens_seq) then
   call init_random_seq(ens_seq)
   first_ens_seq = .false.
endif

do i = 1, size(x)
   do j = 1, size(as, 2)

      ! May want the constant to always be the same as the obs_var.

      as(i, j) = x(i) + 0.5_r8 * random_gaussian(ens_seq, 0.0_r8, 1.0_r8)

   end do
end do

end subroutine ens_ics


subroutine get_close_state(obs_num, list, max_list, num)
!-----------------------------------------------------------------------
! subroutine get_close_state(obs_num, list, max_list, num)
!

implicit none

integer, intent(in)  :: obs_num, max_list
integer, intent(out) :: list(max_list), num

integer :: i

num = 3               ! Everything is close here

do i = 1, 3
   list(i) = obs_num + i - 1
   if(list(i) > 3) list(i) = list(i) - 3
end do

end subroutine get_close_state



!========================================================================
! end module obs_mod
!========================================================================

end module obs_mod
