program filter

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
   add_obs_set, write_obs_sequence, read_obs_sequence, associate_def_list, &
   get_num_obs_sets, get_obs_set
use obs_def_mod, only : obs_def_type, init_obs_def
use obs_set_def_mod, only : obs_set_def_type, init_obs_set_def, add_obs
use obs_kind_mod, only : set_obs_kind
use location_mod, only : set_location
use set_def_list_mod, only : set_def_list_type, init_set_def_list, &
   add_to_list, write_set_def_list
use obs_set_mod, only : init_obs_set, obs_set_type, set_obs_set_time, write_obs_set, &
   get_obs_set_time
use time_manager_mod, only : time_type, set_time
use utilities_mod, only : open_file
use assim_model_mod, only : assim_model_type, initialize_assim_model, get_model_size, &
   get_initial_condition, get_model_state_vector, set_model_state_vector, &
   get_closest_state_time_to, advance_state

implicit none

type(obs_sequence_type) :: seq, seq2
type(obs_def_type) :: obs_def
type(set_def_list_type) :: set_def_list
type(obs_set_def_type) :: obs_set_def
type(obs_set_type) :: obs_set
type(time_type) :: time, time2

integer, parameter :: num_steps = 100
integer, parameter :: max_sets = 5
integer, parameter :: num_obs = 10
real(r8), parameter :: error_variance = 4.0
integer :: i, j, obs_set_def_index, unit, test_unit, obs_set_unit

! Need to set up namelists for controlling all of this mess, too!
integer, parameter :: ens_size = 20
integer :: model_size, num_obs_sets

type(assim_model_type) :: x, ens(ens_size)
real(r8), allocatable :: ens_temp(:)

! Initialize the obs_set_def
obs_set_def = init_obs_set_def(num_obs)

! Build a set of obs_defs and put them into an obs_set_def
do i = 1, num_obs
! Create an obs_def with specificed location, kind and error variance
   obs_def = init_obs_def(set_location(i * 1.0_r8 / num_obs), set_obs_kind(1), error_variance)
! Insert this obs_def into an obs_set
   call add_obs(obs_set_def, obs_def)
end do

! Insert this obs_set_def into the list
set_def_list = init_set_def_list(max_sets)
obs_set_def_index = add_to_list(set_def_list, obs_set_def)

! Initialize a sequence
seq = init_obs_sequence(num_steps, 0)

! Put the obs_set_def into the sequence
call associate_def_list(seq, set_def_list)

! Now generate a long obs_sequence with regular occurences of the obs_set_def
do i = 1, num_steps
   obs_set = init_obs_set(set_def_list, obs_set_def_index, 0)
   call set_obs_set_time(obs_set, set_time(i * 3600, 0))
! Put this obs_set into the sequence
   call add_obs_set(seq, obs_set)
end do

! Output the obs_sequence to a file
unit = open_file('seq_out')
call write_obs_sequence(unit, seq)
close(unit)


! Initialize the model now that obs_sequence is all set up
call initialize_assim_model()
model_size = get_model_size

! Get the initial condition
x = get_initial_condition()

! Perturb for ensembles; (done with no random for now)
allocate(ens_temp(model_size))
do i = 1, ens_size
   ens_temp = get_model_state_vector(ens(i))
   ens_temp(i) = ens_temp(i) + 0.1
   call set_model_state_vector(ens(i), ens_temp)
end do
deallocate(ens_temp)

num_obs_sets = get_num_obs_sets(seq)

do i = 1, num_obs_sets
   obs_set = get_obs_set(seq, i)
   time = get_obs_set_time(obs_set)
   write(*, *) 'time of obs set ', i, time
   time2 = get_closest_state_time_to(x, time)
   write(*, *) 'closest model time is ', time2
! Advance the state and the ensembles to this time
   x = advance_state(x, time2)
   do j = 1, ens_size
      ens(j) =  advance_state(ens(j), time2)
   end do
end do


end program filter
