module obs_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use types_mod
use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian
use model_mod,        only : get_model_size, balance_init

use loc_and_dist_mod, only : loc_type, get_dist, set_loc

private

public num_obs, obs_var, take_obs, ens_ics, init_obs, take_single_obs, &
       get_close_state, state_to_obs, get_num_obs, get_obs_m_rinv, &
       max_num_pos_obs, num_obs_int, obs_var_int, take_obs_int, obs_loc

integer, parameter :: num_obs = 9

! Following is to allow initialization of obs_def_type

  logical :: obs_init = .false.

! Following is for repeatable random numbers

  logical :: first_ens_seq = .true.
  type (random_seq_type) :: ens_seq

! Define a type for model operator

  type obs_def_type
     integer           :: num
     integer,  pointer :: state_var_index(:)
     real(r8), pointer :: coef(:)
  end type obs_def_type

  type (obs_def_type) :: obs_def(num_obs)

  type state_obs_dep_type
     integer          :: num
     integer, pointer :: obs_index(:)
  end type state_obs_dep_type

! type (state_obs_dep_type) :: state_obs_dep(model_size)
  type (state_obs_dep_type) :: state_obs_dep(9)

! Size of arrays for upper limit on number of obs associated with set of data

  integer, parameter :: max_num_pos_obs = 9

! Array of structure for observation locations; static with time in this version
! Not used in 9var as all variables assumed to be equidistant

  type(loc_type) :: obs_loc(num_obs)

contains



  subroutine init_obs
!-----------------------------------------------------------------------
! subroutine init_obs

implicit none

integer  :: i
real(r8) :: loc

obs_init = .true.       ! Initialization for identity observations

!--------------- block for diagonal observations operator -------

do i = 1, num_obs

   obs_def(i)%num = 1

   allocate(obs_def(i)%state_var_index(1), obs_def(i)%coef(1))

   obs_def(i)%state_var_index(1) = i
   obs_def(i)%coef(1) = 1.0_r8

!!!   obs_def(i)%coef(1) = 1.0_r8 * i

   ! Set the location to 0 for all obs

   loc = 0.0_r8
   call set_loc(obs_loc(i), loc)
end do

call init_state_obs_dep

end subroutine init_obs



  function obs_var()
!-----------------------------------------------------------------------
! function obs_var()

implicit none

real(r8) :: obs_var(num_obs)

obs_var(:) = (0.05 ** 2)      ! Default values of all obs for 9var_paper

!obs_var(:) = 1.0
!obs_var(:) = 0.2

!---
!obs_var(:)   = 10.0 ** 2
!obs_var(4:6) = (0.05 **2)
!---

!-- Following is for divergence only obs

!obs_var(:)   = -100.0 **2
!obs_var(1:3) = 0.0025 ** 2
!--

!---
!obs_var(:) = 10.0 **2
!obs_var(4) = (0.05** 2)
!---

end function obs_var



  subroutine get_close_state(obs_num, list, max_list, num)
!----------------------------------------------------------------------
! subroutine get_close_state(obs_num, list, max_list, num)

implicit none

integer, intent(in)  :: obs_num, max_list
integer, intent(out) :: list(max_list), num
integer :: i

   ! Every ob impacts every state

   num = 9
   do i = 1, 9
      list(i) = i
   end do

! Only impact vars 4 to 9
!   num = 6
!   do i = 1, 6
!      list(i) = i + 3
!   end do

 
! Obs only impact same class of variable (see below for type)
!num = 3
!if(obs_num == 1 .or. obs_num == 4 .or. obs_num == 7) then
!   do i = 1, 3
!      list(i) = 3 * (i - 1) + 1
!   end do
!else if(obs_num == 2 .or. obs_num == 5 .or. obs_num == 8) then
!   do i = 1, 3
!      list(i) = 3 * (i - 1) + 2
!   end do
!else
!   do i = 1, 3
!      list(i) = 3*i
!   end do
!endif


! Every ob only impacts its own state variable
!num = 1
!list(1) = obs_num

! Only obs of same type are used
!num = 3
!if(obs_num < 4) then
!   do i = 1, 3
!      list(i) = i
!   end do
!else if(obs_num >=4 .and. obs_num <= 6) then
!   do i = 1, 3
!      list(i) = 3 + i
!   end do
!else
!   do i = 1, 3
!      list(i) = 6 + i
!   end do
!end if

end subroutine get_close_state



function obs_var_int()
!-----------------------------------------------------------------------
! function obs_var_int()

implicit none

real(r8) :: obs_var_int(num_obs)

!obs_var_int = (0.20)**2
 obs_var_int = (0.01)**2

end function obs_var_int



function take_obs(x)
!-----------------------------------------------------------------------
! function take_obs(x)

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: take_obs(num_obs)

! Initialize observation state

if(.not. obs_init) call init_obs

! This is function h that maps state variables to obs

!take_obs(:) = x(4:9)
take_obs = x

end function take_obs



  function take_single_obs(x, index)
!-----------------------------------------------------------------------
! function take_single_obs(x, index)

implicit none

real(r8), intent(in) :: x(:)
integer,  intent(in) :: index
real(r8)             :: take_single_obs

take_single_obs = x(index)

end function take_single_obs



  function state_to_obs(x)
!-----------------------------------------------------------------------
!????????????????????????????????????????????????????????
! function state_to_obs(x)

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: state_to_obs(num_obs)

integer :: i, j

if(.not. obs_init) call init_obs

state_to_obs = 0.0_r8
do i = 1, num_obs
   do j = 1, obs_def(i)%num
      state_to_obs(i) = state_to_obs(i) + &
         obs_def(i)%coef(j) * x(obs_def(i)%state_var_index(j))
   end do
end do

end function state_to_obs



  function take_obs_int(x)
!-----------------------------------------------------------------------
! function take_obs_int(x)

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: take_obs_int(num_obs)

! This is function h that maps state variables to obs

take_obs_int = x

end function take_obs_int



  subroutine ens_ics(x, as)
!------------------------------------------------------------------------
! subroutine ens_ics(x, as)
!
! Get initial state for ensemble assimilation

implicit none


real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: as(:, :)

integer  :: i, j
real(r8) :: var(size(x)), obs(size(x))

if ( first_ens_seq ) then
   call init_random_seq(ens_seq)
   first_ens_seq = .false.
endif

! Added as part of revised ICS for 13 Feb. as outlined below

var = obs_var()

! For large variances want this to not cause overflow

do i = 1, size(x)
   if(var(i) < 0.0_r8) var(i) = 0.1_r8
end do

do i = 1, size(x)
   obs(i) = random_gaussian(ens_seq, x(i), var(i))
end do

do i = 1, size(x)
   do j = 1, size(as, 2)

      ! Results from 13 Feb. 2001 show that centered ens_ics as above can lead
      ! to VERY persistent features in ensemble spread ratios. Revised ens_ics
      ! to generate a distribution from a single ob and obs_var. Only works if
      ! obs are identity as they are in all cases here.

      as(i, j) = random_gaussian(ens_seq, obs(i), var(i))
   end do
end do

! Need to balance these?
!do j = 1, size(as, 2)
!   call balance_init(as(:, j), as(:, j))
!end do

end subroutine ens_ics



  subroutine get_num_obs(list, num_state, obs_list, num_obs)
!------------------------------------------------------------------------
! subroutine get_num_obs(list, num_state, obs_list, num_obs)
!

implicit none

integer, intent(in)  :: num_state, list(num_state)
integer, intent(out) :: obs_list(max_num_pos_obs), num_obs

integer :: i, j, num_pos_obs, ind, ob, k

! First get a non-redundant list of all the obs that are candidates

num_pos_obs = 0

StateLoop : do i = 1, num_state
   ind = list(i)

   DepLoop : do j = 1, state_obs_dep(ind)%num

      ob = state_obs_dep(ind)%obs_index(j)

      ! See if this is in list of obs

      do k = 1, num_pos_obs
         if(obs_list(k) == ob) cycle DepLoop
      end do

      num_pos_obs = num_pos_obs + 1

      if(num_pos_obs > max_num_pos_obs) then
         write(*, *) 'error in get_num_obs; max_num_pos_obs exceeded'
         stop
      endif

      obs_list(num_pos_obs) = ob
 
   end do DepLoop
end do StateLoop

! Loop through possible obs and see which have all required state vars in list

num_obs = 0

do i = 1, num_pos_obs

   ob = obs_list(i)

   do j = 1, obs_def(ob)%num
      if(.not. any(list == obs_def(ob)%state_var_index(j))) goto 20
   end do

   num_obs = num_obs + 1

   ! Writing over list at bottom; should be safe

   obs_list(num_obs) = ob

20 continue
end do

end subroutine get_num_obs



  subroutine get_obs_m_rinv(state_list, n_state, obs_list, n_obs, r_inv, m)
!------------------------------------------------------------
! subroutine get_obs_m_rinv(state_list, n_state, obs_list, n_obs, r_inv, m)
!

implicit none

integer,  intent(in)  :: n_state, n_obs, state_list(n_state), obs_list(n_obs)
real(r8), intent(out) :: m(n_state, n_obs), r_inv(n_obs, n_obs)

integer  :: i, ob, sv, j, k
real(r8) :: obs_variance(num_obs)

! Get inverse of covariance matrix for this set of obs; for now assues
! that the full observational covariance matrix is diagonal (could be extended)

obs_variance = obs_var()
r_inv        = 0.0_r8

do i = 1, n_obs
   ob          = obs_list(i)
   r_inv(i, i) = 1.0_r8 / obs_variance(ob)
end do

! Get M, must be a more efficient way (sort at start?), table??

m = 0.0_r8

ILoop : do i = 1, n_state
   sv = state_list(i)
   JLoop : do j = 1, n_obs
      ob = obs_list(j)

      KLoop : do k = 1, obs_def(ob)%num
	 if(obs_def(ob)%state_var_index(k) == sv) then
            m(i, j) = obs_def(ob)%coef(k)
       !    goto 10
            cycle JLoop
         endif

      end do KLoop

10    continue

   end do JLoop
end do ILoop

end subroutine get_obs_m_rinv



  subroutine init_state_obs_dep
!------------------------------------------------------------
! subroutine init_state_obs_dep
!

implicit none

integer :: i, j, num, ind, model_size

! Initialize count of obs per state variable to 0

model_size = get_model_size()

do i = 1, model_size
   state_obs_dep(i)%num = 0
end do

! Get count of obs per state variable

do i = 1, num_obs
   do j = 1, obs_def(i)%num
      ind = obs_def(i)%state_var_index(j)
      state_obs_dep(ind)%num = state_obs_dep(ind)%num + 1
   end do
end do

! Initialize space for list of these obs; reset num to be pointer for loading

do i = 1, model_size
   num = state_obs_dep(i)%num
   allocate(state_obs_dep(i)%obs_index(num))
   state_obs_dep(i)%num = 0
end do

do i = 1, num_obs
   do j = 1, obs_def(i)%num
      ind = obs_def(i)%state_var_index(j)
      num = state_obs_dep(ind)%num + 1
      state_obs_dep(ind)%num = num
      state_obs_dep(ind)%obs_index(num) = i
   end do
end do

end subroutine init_state_obs_dep


!------------------------------------------------------------
! End of 9var/obs/obs_mod.f90 Module 
!------------------------------------------------------------

end module obs_mod
