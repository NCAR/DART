module obs_mod

! MODIFIED FROM REGULAR L96 OBS, 10 JULY, 2000

use types_mod
!use nag_wrap_mod,    only : g05ddf_wrap
use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian
use model_mod,        only : model_size
use obs_tools_mod,    only : conv_state_to_obs, obs_def_type, obs_m_rinv, &
                             dep_obs, def_single_obs, &
                             state_obs_dep_type, init_state_obs_dep
use loc_and_dist_mod, only : loc_type, set_loc, get_loc, get_dist
use utilities_mod,    only : open_file, close_file

private

public num_obs, obs_var, take_obs, ens_ics, init_obs, take_single_obs, &
       state_to_obs, max_num_pos_obs, get_close_obs

!  Following are not needed for lim_obs_pattern
!  public get_num_obs, get_obs_m_rinv

   integer :: num_obs = 0

!  Following is to allow initialization of obs_def_type

   logical :: obs_init = .false.

!  Following is for repeatable random numbers

   logical :: first_ens_seq = .true.
   type (random_seq_type) :: ens_seq

!  The obs_def structure defines a linear M operator from state vars to obs

   type (obs_def_type), allocatable :: obs_def(:)

!  Array of structure for observation locations; static w/time in this version

   type(loc_type), allocatable :: obs_loc(:)

!  Storage for the observational variance

   real(r8), allocatable :: obs_variance(:)

!  Set a cut-off for distance on which obs can be used; only look 4 points away

   real(r8), parameter :: max_close_distance = 50.0_r8 / dble(model_size)
!  real(r8), parameter :: max_close_distance = 10.0_r8 / dble(model_size)
!  real(r8), parameter :: max_close_distance =  5.0_r8 / dble(model_size)
!  real(r8), parameter :: max_close_distance = 0.25_r8



! The state_obs_dep operator indicates which obs depend on each state variable

type (state_obs_dep_type) :: state_obs_dep(model_size)

! Size of arrays for upper limit on number of obs associated with set of data

integer, parameter :: max_num_pos_obs = 50



contains



  subroutine init_obs
!=======================================================================
! subroutine init_obs
!
! Initializes the description a linear observations operator. For each
! observation, a list of the state variables on which it depends and the
! coefficients for each state variable is passed to def_single_obs
! which establishes appropriate data structures. init_state_obs_dep
! is then called to set up the state observations relation data structure.
 
implicit none

integer  :: i, state_index(2), unit_num, model_size
real(r8) :: coef(2), x, loc, pos

obs_init = .true.        ! Initialization for identity observations

! Read in the lorenz_96_obs_def file

!!!unit_num = Get_Unit()
unit_num = open_file(file = 'lorenz_96_obs_def', action = 'read')
!!!open(unit = unit_num, file = 'lorenz_96_obs_def')
read(unit_num, *) num_obs
allocate(obs_def(num_obs), obs_loc(num_obs), obs_variance(num_obs))

do i = 1, num_obs
   read(unit_num, *) loc, obs_variance(i)
   if(loc < 0.0_r8 .or. loc .gt. 1.0_r8) then
      write(*, *) 'illegal obs location read in init_obs for lorenz_96'
      stop
   end if

   call set_loc(obs_loc(i), loc) 
end do

close(unit = unit_num)

! Now loop through the obs and figure out how to get linear interpolation M

model_size = get_model_size()

do i = 1, num_obs
   call get_loc(obs_loc(i), x)
   pos = x * (model_size - 1)
   state_index(1) = int(pos) + 1
   state_index(2) = state_index(1) + 1
   if(state_index(2) > model_size - 1) state_index(2) = 1
   coef(1) = 1.0_r8 - (pos - int(pos))
   coef(2) = 1.0_r8 - coef(1)
   call def_single_obs(2, state_index(1:2), coef(1:2), obs_def(i))
end do

! Set up data structure of dependencies of obs on state; 
! Not needed for lim_obs_pattern development
!call init_state_obs_dep(obs_def, num_obs, state_obs_dep, model_size)

end subroutine init_obs



  function obs_var()
!=======================================================================
! function obs_var()
!
! Defines the observational error variance. Eventually may need to be 
! generalized to covariance.

implicit none

real(r8) :: obs_var(num_obs)

obs_var = obs_variance

end function obs_var



  function take_obs(x)
!=======================================================================
! function take_obs(x)
!
! Given a model state, x, returns observations for assimilation.
! For perfect model, take_obs is just state_to_obs

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: take_obs(num_obs)

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

take_obs = conv_state_to_obs(x, obs_def, num_obs)

end function take_obs



  function take_single_obs(x, index)
!========================================================================
! function take_single_obs(x, index)
!
! Given a model state, x, returns observations for assimilation.
! For perfect model, take_obs is just state_to_obs

implicit none

real(r8), intent(in) :: x(:)
integer,  intent(in) :: index
real(r8)             :: take_single_obs

real(r8) :: take(1)

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

take = conv_state_to_obs(x, obs_def(index:index), 1)
take_single_obs = take(1)

end function take_single_obs



  function state_to_obs(x)
!========================================================================
! function state_to_obs(x)
!
! Given a model state, returns the associated 'observations' by applying
! the observational operator to the state.

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: state_to_obs(num_obs)

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

state_to_obs = conv_state_to_obs(x, obs_def, num_obs)

end function state_to_obs



  subroutine ens_ics(x, as)
!========================================================================
! subroutine ens_ics(x, as)
!
! Get initial ensemble states for ensemble assimilation. This may not be
! the appropriate module for this in the long run.

implicit none

real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: as(:, :)

integer :: i, j

if(first_ens_seq) then
   call init_random_seq(ens_seq)
   first_ens_seq = .false.
endif

! TJH 15.03.2002 leaving explicit typing to dble for g05ddf_wrap, unable
!                to determine if this is an actual library routine or one
!                under our control.

do i = 1, size(x)
   do j = 1, size(as, 2)
! May want the constant to always be the same as the obs_var.
!     as(i, j) = x(i) + 0.1_r8 * g05ddf_wrap(dble(0.0), dble(1.0))
      as(i, j) = x(i) + 0.1_r8 * random_gaussian(ens_seq, 0.0_r8, 1.0_r8)
   end do
end do

end subroutine ens_ics


TJH bomb note 15.03.2002  herein lies the problem ...
The other lorenz_96 obs modules have a "get_close_state" routine
and use sort_mod, this obs module has several more pieces.



  subroutine get_num_obs(list, num_state, obs_list, num_dep_obs)
!=======================================================================
! subroutine get_num_obs(list, num_state, obs_list, num_dep_obs)
!
! Given a list of state variables, returns a list of observations that
! depend exclusively on these state variables in obs_list and the number
! of the dependent obs in num_dep_obs.

implicit none

integer, intent(in)  :: num_state, list(num_state)
integer, intent(out) :: obs_list(max_num_pos_obs), num_dep_obs

call dep_obs(state_obs_dep, model_size, obs_def, num_obs, list, num_state, &
   obs_list, max_num_pos_obs, num_dep_obs)

end subroutine get_num_obs



  subroutine get_obs_m_rinv(state_list, n_state, obs_list, n_obs, r_inv, m)
!==========================================================================
! subroutine get_obs_m_rinv(state_list, n_state, obs_list, n_obs, r_inv, m)
!
! Given a list of state variables and observations, returns the observations
! operator, M, that operates on these state variables to get these obs and
! the inverse of the appropriate observational covariance matrix.

implicit none

integer, intent(in) :: n_state, n_obs, state_list(n_state), obs_list(n_obs)
real(r8), intent(out) :: m(n_state, n_obs), r_inv(n_obs, n_obs)

call obs_m_rinv(obs_def, num_obs, obs_variance, state_list, n_state, obs_list, &
   n_obs, r_inv, m)

end subroutine get_obs_m_rinv



  subroutine get_close_obs(list, num_close, max_num, state_loc)
!==========================================================================
! subroutine get_close_obs(list, num_close, max_num, state_loc)
!
! Computes a list of observations 'close' to each state variable and
! returns the list sorted by the measure of 'closeness' from closest to
! least close. Num is the number of close obs requested for each state
! variable, list is the returned array of close points, and state_loc is
! an array of structure containing information about the location of the
! state variables. This information may not be used in some cases of 
! computing the close obs. In big models, careful attention to efficiency
! will be required in computing this list to avoid large O(n**2) computations.

implicit none

integer, intent(in) :: max_num
integer, intent(inout) :: list(model_size, max_num), num_close(model_size)

type (loc_type) :: state_loc(model_size)

integer :: i, j, k, index, temp_ind
real(r8) :: temp_dist, dist(max_num), d

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

! This will be far too inefficient for big models: what to do about it?
! Upper bound on dist; Should check for this in computation and die

! Loop through every pair of state/obs and compute distance. Do downward
! bubble sort of close obs for each state. For lorenz-96, know that distance
! is in range 0 to 1.
list = -99
do i = 1, model_size
! Make initial distance for this state variable huge and replace below
   dist = huge(dist)
   do j = 1, num_obs
      d = get_dist(state_loc(i), obs_loc(j))
! If distance exceeds max distance allowed, move on to check next obs
      if(d > max_close_distance) goto 131
! Downward bubble sort for now
      if(d < dist(max_num)) then
         dist(max_num) = d
         list(i, max_num) = j
         do k = max_num - 1, 1, -1
            if(dist(k) > dist(k + 1)) then
               temp_dist = dist(k + 1)
               temp_ind = list(i, k + 1)
               dist(k + 1) = dist(k)
               list(i, k + 1) = list(i, k)
               dist(k) = temp_dist
               list(i, k) = temp_ind 
            end if
         end do
      end if
131 continue
   end do
end do

! Make sure that all entries have been filled. List value of -99 is fatal

do i = 1, model_size
   do j = 1, max_num
      if(list(i, j) == -99) then
         num_close(i) = j -1
         goto 141
      endif
   end do
! Fall off end, number of close obs if max_num
   num_close(i) = max_num
141 continue
end do

end subroutine get_close_obs

!==========================================================================
! end of lorenz_96/obs/imperfect/obs_mod.f90
!==========================================================================

end module obs_mod
