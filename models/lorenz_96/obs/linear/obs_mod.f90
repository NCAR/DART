module obs_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use types_mod
!use nag_wrap_mod,    only : g05ddf_wrap
use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian
use model_mod,        only : get_model_size
use obs_tools_mod,    only : conv_state_to_obs, obs_def_type, obs_m_rinv, &
                             dep_obs, def_single_obs
use loc_and_dist_mod, only : loc_type, set_loc, get_loc
use utilities_mod,    only : open_file, close_file
use sort_mod,         only : index_sort

private

public num_obs, obs_var, take_obs, ens_ics, init_obs, take_single_obs, &
       get_close_state


!  TEMPORARY ADDITION TO ALLOW RAMPED IMPACT

public obs_loc

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

!  Set a cut-off for distance on which obs can be used

!  real(r8), parameter :: max_close_distance =     1.0_r8
   real(r8), parameter :: max_close_distance =   0.175_r8
!  real(r8), parameter :: max_close_distance =    1.00_r8
!  real(r8), parameter :: max_close_distance = 0.00011_r8

contains

!=======================================================================



  subroutine init_obs
!=======================================================================
! subroutine init_obs
!
! Initializes the description a linear observations operator. For each
! observation, a list of the state variables on which it depends and the
! coefficients for each state variable is passed to def_single_obs
! which establishes appropriate data structures. 
 
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

call close_file(unit_num)

! Now loop through the obs and figure out how to get linear interpolation M

model_size = get_model_size()

do i = 1, num_obs
   call get_loc(obs_loc(i), x)
   pos = x * model_size
   state_index(1) = int(pos) + 1
   state_index(2) = state_index(1) + 1
   if(state_index(2) > model_size) state_index(2) = 1
   coef(1) = 1.0_r8 - (pos - int(pos))
   coef(2) = 1.0_r8 - coef(1)
   call def_single_obs(2, state_index(1:2), coef(1:2), obs_def(i))
end do

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

! TJH 14.03.2002 leaving explicit typing to dble for g05ddf_wrap, unable
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



  subroutine get_close_state(obs_num, list, max_list, num)
!=======================================================================
! subroutine get_close_state(obs_num, list, max_list, num)
!

implicit none

integer, intent( in) :: obs_num, max_list
integer, intent(out) :: list(max_list), num

real(r8) :: x, x_min, x_max, dist(max_list)
integer  :: ind_min, ind_max, i, model_size
integer  :: index(max_list), list_temp(max_list)

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

! Get location of this ob

call get_loc(obs_loc(obs_num), x)
model_size = get_model_size()

! Model specific computation of which state variables are close
! THe 1.000001 factor is a quick fix to avoid round-off problems
! for obs that are on evenly spaced grid points.

x_min   = x - 1.000001_r8 * max_close_distance
x_max   = x + 1.000001_r8 * max_close_distance
ind_min = ceiling(x_min * model_size)
ind_max = floor(x_max * model_size)
num     = ind_max - ind_min + 1
if(num > model_size) num = model_size

if(num > max_list) then
   write(*, *) 'max list length is too short in get_close_state'
   stop
endif

ind_min = modulo(ind_min, model_size)

! Generate index list

do i = 1, num
   list(i) = modulo(ind_min + i - 1, model_size) + 1
end do

if(1 == 1) return

write(*, *) '                          SORTING'

! TEST CODE ADDED 10 FEB, 2001, sort list by distance, closest first?
! HAS NO IMPACT SANS 5 DECEMBER RESAMPLING OF OBS
do i = 1, num
   dist(i) = abs(x - (list(i) - 1.0_r8) / model_size)
   if(dist(i) > 0.5_r8) dist(i) = 1.0_r8 - dist(i)
!   write(*, *) 'list/dist ', i, x, list(i), dist(i)
end do

call index_sort(dist(1:num), index(1:num), num)
do i = 1, num
   list_temp(i) = list(index(i))
end do
list(1:num) = list_temp(1:num)
! END OF SORT BLOCK

end subroutine get_close_state

!==========================================================================
! end module lorenz_96/obs/linear/obs_mod.f90
!==========================================================================

end module obs_mod
