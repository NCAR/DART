! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

!use nag_wrap_mod, only : g05ddf_wrap
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use      model_mod, only : get_model_size 
use   location_mod, only : loc_type, get_dist, set_loc

implicit none
private

public :: num_obs, obs_var, take_obs, ens_ics, get_close_state, &
   max_num_pos_obs, init_obs, take_single_obs, obs_loc

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

integer, parameter :: num_obs = 2
!!!integer, parameter :: num_obs = 1

! Following is to allow initialization of obs_def_type
logical :: obs_init = .false.

! Following is for repeatable random numbers
logical :: first_ens_seq = .true.
type (random_seq_type) :: ens_seq

! Define a type for model operator
type obs_def_type
   integer :: num
   integer, pointer :: state_var_index(:)
   double precision, pointer :: coef(:)
end type obs_def_type

type (obs_def_type) :: obs_def(num_obs)

type state_obs_dep_type
   integer :: num
   integer, pointer :: obs_index(:)
end type state_obs_dep_type

!type (state_obs_dep_type) :: state_obs_dep(model_size)
type (state_obs_dep_type) :: state_obs_dep(2)

! Size of arrays for upper limit on number of obs associated with set of data
integer, parameter :: max_num_pos_obs = 2

! Array of structure for observation locations; static with time in this version
! Not used in 9var as all variables assumed to be equidistant
type(loc_type) :: obs_loc(num_obs)


contains

!-----------------------------------------------------------------------

subroutine init_obs

implicit none

integer :: i
double precision :: loc

! Initialization for identity observations
obs_init = .true.

if(1 == 1) return

end subroutine init_obs

!-----------------------------------------------------------------------

function obs_var()

implicit none

double precision :: obs_var(num_obs)

! Default values of all obs for 9var_paper
obs_var(:) = 1.0

end function obs_var

!----------------------------------------------------------------------

subroutine get_close_state(obs_num, list, max_list, num)

implicit none

integer, intent(in) :: obs_num, max_list
integer, intent(out) :: list(max_list), num
integer :: i

! Every ob impacts evry state
   num = 2
   do i = 1, 2
      list(i) = i
   end do

end subroutine get_close_state

!-----------------------------------------------------------------------

function take_obs(x)

implicit none

double precision :: take_obs(num_obs)
double precision, intent(in) :: x(:)

! Initialize observation state
if(.not. obs_init) call init_obs

! This is function h that maps state variables to obs
take_obs = x

! Single ob of position
!take_obs = x(1)


! Obs are amplitude and phase in this experiment
!!!take_obs(1) = sqrt(x(1)**2 + x(2)**2)
!!!take_obs(2) = atan2(x(1), x(2))

end function take_obs

!-----------------------------------------------------------------------

function take_single_obs(x, index)

implicit none

double precision :: take_single_obs
double precision, intent(in) :: x(:)
integer, intent(in) :: index

! Default case is identity obs
take_single_obs = x(index)


! Only single ob of position
!if(index == 1) then
!   take_single_obs = x(1)
!else
!   write(*, *) 'index is not 1 in take_single_obs'
!   stop
!endif

! Experimental case is amplitude and phase
!if(index == 1) then
!   take_single_obs = sqrt(x(1)**2 + x(2)**2)
!else
!   take_single_obs = atan2(x(1), x(2))
!      continue
!endif

end function take_single_obs

!-----------------------------------------------------------------------

subroutine ens_ics(x, as)

implicit none

!  Get initial state for ensemble assimilation

double precision, intent(in) :: x(:)
double precision, intent(out) :: as(:, :)
integer :: i, j
double precision :: var(size(x)), obs(size(x))

if(first_ens_seq) then
   call init_random_seq(ens_seq)
   first_ens_seq = .false.
endif

! Added as part of revised ICS for 13 Feb. as outlined below
var = obs_var()
! For large variances want this to not cause overflow
do i = 1, size(x)
   if(var(i) < 0.0) var(i) = 0.1
end do

do i = 1, size(x)
   obs(i) = random_gaussian(ens_seq, x(i), var(i))
end do

do i = 1, size(x)
   do j = 1, size(as, 2)
      as(i, j) = random_gaussian(ens_seq, obs(i), var(i))
   end do
end do

end subroutine ens_ics

!------------------------------------------------------------------------

end module obs_mod
