program create_obs_sequence

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

! Program to build a simple obs_sequence file for use in testing filters
! for spatial domains with one periodic dimension.

use types_mod
use obs_sequence_mod, only : init_obs_sequence, obs_sequence_type, &
   add_obs_set
use obs_def_mod, only : obs_def_type, init_obs_def
use obs_set_def_mod, only : obs_set_def_type, init_obs_set_def, add_obs
use obs_kind_mod, only : set_obs_kind
use location_mod, only : set_location
use set_def_list_mod, only : set_def_list_type, init_set_def_list, &
   add_to_list
use obs_set_mod, only : init_obs_set, obs_set_type, set_obs_set_time
use time_manager_mod, only : time_type, set_time

implicit none

type(obs_sequence_type) :: seq
type(obs_def_type) :: obs_def
type(set_def_list_type) :: set_def_list
type(obs_set_def_type) :: obs_set_def
type(obs_set_type) :: obs_set

integer, parameter :: num_steps = 100
integer, parameter :: max_sets = 5
integer, parameter :: num_obs = 10
real(r8), parameter :: error_variance = 4.0
integer :: i, obs_set_def_index

! Initialize the obs_set_def
obs_set_def = init_obs_set_def(num_obs)

! Build a set of obs_defs and put them into an obs_set_def
do i = 1, num_obs
! Create an obs_def with specificed location, kind and error variance
   obs_def = init_obs_def(set_location(i * 1.0_r8 / num_obs), set_obs_kind(1), error_variance)
! Insert this obs_def into an obs_set
   call add_obs(obs_set_def, obs_def)
end do

! Insert this obs_set_def into the 
set_def_list = init_set_def_list(max_sets)
obs_set_def_index = add_to_list(set_def_list, obs_set_def)

! Initialize a sequence
seq = init_obs_sequence(num_steps, 0)

! Now generate a long obs_sequence with regular occurences of the obs_set_def
do i = 1, num_steps
   obs_set = init_obs_set(set_def_list, obs_set_def_index, 0)
   call set_obs_set_time(obs_set, set_time(i * 3600, 0))
! Put this obs_set into the sequence
   call add_obs_set(seq, obs_set)
end do

end program create_obs_sequence
