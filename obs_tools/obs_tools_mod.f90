module obs_tools_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Provides a set of tools that are used by observations modules to generate
! data structures and information needed to do ensemble assimilation.
! These should be common to a variety of models and observational operator
! types.

use types_mod

private

public conv_state_to_obs, obs_def_type, obs_m_rinv, state_obs_dep_type, &
   dep_obs, init_state_obs_dep, def_single_obs, get_state_var_indices


! Define a type for linear model operator

type obs_def_type
   private
   integer :: num
   integer, pointer :: state_var_index(:)
   real(r8), pointer :: coef(:)
end type obs_def_type

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Type that keeps list of which observations depend on which state variables

type state_obs_dep_type
   private
   integer :: num
   integer, pointer :: obs_index(:)
end type state_obs_dep_type

contains




  function conv_state_to_obs(x, obs_def, num_obs)
!=======================================================================
! function conv_state_to_obs(x, obs_def, num_obs)
!
! Given a state vector x, returns corresponding observations given the
! observations data definition structure obs_def.

implicit none

real(r8),            intent(in) :: x(:)
integer,             intent(in) :: num_obs
type (obs_def_type), intent(in) :: obs_def(num_obs)
real(r8)                        :: conv_state_to_obs(num_obs)

integer :: i, j

conv_state_to_obs = 0.0_r8
do i = 1, num_obs
   do j = 1, obs_def(i)%num
      conv_state_to_obs(i) = conv_state_to_obs(i) + &
         obs_def(i)%coef(j) * x(obs_def(i)%state_var_index(j))
   end do
end do

end function conv_state_to_obs


  subroutine dep_obs(state_obs_dep, model_size, obs_def, num_obs, list, &
                     num_state, obs_list, max_num_pos_obs, num_dep_obs)
!=======================================================================
! subroutine dep_obs(state_obs_dep, model_size, obs_def, num_obs, list, &
!                    num_state, obs_list, max_num_pos_obs, num_dep_obs)
!
! Given a state_obs_dep strucutre and an obs_def structure along with a
! list of state variables and observations, returns the number of observations
! that depend exclusively on this set of state variables.

implicit none

integer,                  intent(in)  :: num_state, list(num_state), model_size
integer,                  intent(in)  :: max_num_pos_obs, num_obs
type(state_obs_dep_type), intent(in)  :: state_obs_dep(model_size)
type(obs_def_type),       intent(in)  :: obs_def(num_obs)
integer,                  intent(out) :: obs_list(max_num_pos_obs), num_dep_obs

integer :: i, j, num_pos_obs, ind, ob, k

! First get a non-redundant list of all the obs that are candidates

num_pos_obs = 0

StateLoop : do i = 1, num_state

   ind = list(i)

   DepLoop : do j = 1, state_obs_dep(ind)%num
      ob = state_obs_dep(ind)%obs_index(j)

      ! See if this is in list of obs

      do k = 1, num_pos_obs
       ! if(obs_list(k) == ob) goto 10
         if(obs_list(k) == ob) cycle DepLoop
      end do

      num_pos_obs = num_pos_obs + 1
      if(num_pos_obs > max_num_pos_obs) then
         write(*, *) 'error in dep_obs; max_num_pos_obs exceeded'
         stop
      endif
      obs_list(num_pos_obs) = ob

 10   continue
   end do DepLoop
end do StateLoop

! Loop through possible obs and see which have all required state vars in list

num_dep_obs = 0

ILoop : do i = 1, num_pos_obs

   ob = obs_list(i)
   JLoop : do j = 1, obs_def(ob)%num
    ! if(.not. any(list == obs_def(ob)%state_var_index(j))) goto 20
      if(.not. any(list == obs_def(ob)%state_var_index(j))) cycle ILoop
   end do JLoop
   num_dep_obs = num_dep_obs + 1

   !Writing over list at bottom; should be safe
   obs_list(num_dep_obs) = ob

20 continue
end do ILoop

end subroutine dep_obs



  subroutine obs_m_rinv(obs_def, total_num_obs, obs_variance, state_list, &
                        n_state, obs_list, n_obs, r_inv, m)
!==========================================================================
! subroutine obs_m_rinv(obs_def, total_num_obs, obs_variance, state_list, &
!                       n_state, obs_list, n_obs, r_inv, m)
!
! Computes the observations operator, M, to operate on a set of state variables
! to get a set of observation variables, both provided in lists. Also returns
! the inverse of the observational error covariance matrix for the set of obs.

implicit none

integer,             intent(in)  :: total_num_obs, n_state, n_obs
type (obs_def_type), intent(in)  :: obs_def(total_num_obs)
integer,             intent(in)  :: state_list(n_state), obs_list(n_obs)

! Should we be using full covariance for obs yet???

real(r8),            intent(in)  :: obs_variance(total_num_obs)
real(r8),            intent(out) :: m(n_state, n_obs), r_inv(n_obs, n_obs)

integer :: i, ob, sv, j, k

! Get inverse of covariance matrix for this set of obs; for now assues
! that the full observational covariance matrix is diagonal (could be extended)

r_inv = 0.0_r8
do i = 1, n_obs
   ob = obs_list(i)
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
         !  goto 10
            cycle JLoop
         endif
      end do KLoop

10    continue
   end do JLoop
end do ILoop

end subroutine obs_m_rinv



  subroutine init_state_obs_dep(obs_def, num_obs, state_obs_dep, model_size)
!==========================================================================
! subroutine init_state_obs_dep(obs_def, num_obs, state_obs_dep, model_size)
!
! Given a completed obs_def structure, sets up the state_obs_dep structure
! to be used by obs_m_rinv and get_num_obs.

implicit none

integer,                  intent(in)    :: num_obs, model_size
type(obs_def_type),       intent(in)    :: obs_def(num_obs)
type(state_obs_dep_type), intent(inout) :: state_obs_dep(model_size)

integer :: i, j, num, ind

! Initialize count of obs per state variable to 0

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



  subroutine def_single_obs(num_state, state_ind, coef, obs_def)
!=======================================================================
! subroutine def_single_obs(num_state, state_ind, coef, obs_def)
!
! Puts definition of a single observation into an obs_def data structure.

implicit none

integer,            intent(in)    :: num_state, state_ind(num_state)
real(r8),           intent(in)    :: coef(num_state)
type(obs_def_type), intent(inout) :: obs_def

integer :: i

! Set aside storage for defining this ob

obs_def%num = num_state
allocate(obs_def%state_var_index(num_state), obs_def%coef(num_state))

! Load the state variable index and coefficient for each state variable

do i = 1, num_state
   obs_def%state_var_index(i) = state_ind(i)
   obs_def%coef(i) = coef(i)
end do

end subroutine def_single_obs



  subroutine get_state_var_indices(obs_def, state_ind, n)
!=======================================================================
! subroutine get_state_var_indices(obs_def, state_ind, n)
!
! Gets list of state variables upon which a given ob depends

implicit none

type(obs_def_type), intent(in)  :: obs_def
integer,            intent(out) :: n, state_ind(:)

integer :: i

n = obs_def%num

if(size(state_ind, 1) < n) then
   write(*, *) 'fatal: get_state_var_indices: state_ind array too small '
   stop
endif

do i = 1, n
   state_ind(i) = obs_def%state_var_index(i)
end do

end subroutine get_state_var_indices

!=======================================================================
! End of obs_tools/obs_tools_mod.f90
!=======================================================================

end module obs_tools_mod
