! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_obs_set_def

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

! Main program to interactively create a set_def_list description file for a 
! particular spatial domain and obs_kind set. This is a prototype for more
! user friendly GUI driven methods of creating set_def_lists. For now, there
! is no support for nested observation definition subsets, but that will be
! needed in the long run.

use      obs_def_mod, only : obs_def_type, init_obs_def, interactive_obs_def
use  obs_set_def_mod, only : obs_set_def_type, init_obs_set_def, add_obs
use set_def_list_mod, only : set_def_list_type, init_set_def_list, &
   add_to_list, write_set_def_list, read_set_def_list
use    utilities_mod, only : open_file, register_module, error_handler, E_MSG, E_WARN, E_ERR, &
                             initialize_utilities, finalize_utilities, logfileunit
use  assim_model_mod, only : static_init_assim_model, get_model_size

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_def_type) :: obs_def
type(set_def_list_type) :: set_def_list
type(obs_set_def_type) :: obs_set_def

integer :: max_sets, num_obs
integer :: i, j, obs_set_def_index, unit
character(len = 129) :: file_name

! Write program control information, once to logfile, once to stdout
call initialize_utilities()
write(logfileunit,*)'STARTING create_obs_set_def ...'

call register_module(source, revision, revdate)

! Get output filename
write(*, *) 'Input the filename for output of observation set_def_list? [set_def.out]'
read(*, *) file_name

! With potential for perfect model observations, must do static_init of model
call static_init_assim_model()
write(*, *) 'model size is ', get_model_size()

! Initialize a set_def_list; for now need to specify max_sets
write(*, *) 'Input the number of unique observation sets you might define'
read(*, *) max_sets
set_def_list = init_set_def_list(max_sets)

! Loop through to get definitions for each set
SetDefLoop : do i = 1, max_sets

   write(*, *) 'How many observations in set ', i
   read(*, *) num_obs

   ! Initialize the obs_set_def
   obs_set_def = init_obs_set_def(num_obs)

   ! Build a set of obs_defs and put them into an obs_set_def
   do j = 1, num_obs
      ! Make calls to obs_def to define observations
      write(*, *) 'Defining observation ', j
      obs_def = interactive_obs_def()
      ! Insert this obs_def into an obs_set
      call add_obs(obs_set_def, obs_def)
   end do

   ! Insert this obs_set_def into the list
   obs_set_def_index = add_to_list(set_def_list, obs_set_def)

end do SetDefLoop

! write(*,*)'DEBUG(obs_set_def:create_obs_set_def): finished definition for each set'
! write(*,*)'DEBUG(obs_set_def:create_obs_set_def): opening file ',trim(adjustl(file_name))

! Output the set_def_list
unit = open_file(file_name, action = 'write')
! write(*,*)'DEBUG(obs_set_def:create_obs_set_def): as unit ',unit
call write_set_def_list(unit, set_def_list)
close(unit)

! Read and rewrite as test
!unit = open_file(file_name)
!set_def_list = read_set_def_list(unit)
!close(unit)

!unit = open_file(file_name)
!call write_set_def_list(unit, set_def_list)
!close(unit)

write(*, *)trim(adjustl(file_name)),' successfully created.'
write(*, *)'Terminating normally.'

write(logfileunit, *)trim(adjustl(file_name)),' successfully created.'
write(logfileunit, *)'Terminating normally.'
write(logfileunit, *)'FINISHED create_obs_set_def ...'
write(logfileunit, *)
call finalize_utilities()       ! gracefully closes logfile.

end program create_obs_set_def
