! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program identity_obs_set_def

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
! needed in the long run. This version is used to create sets with many identity
! observations in them.

use        types_mod, only : r8
use      obs_def_mod, only : obs_def_type, init_obs_def, interactive_obs_def
use  obs_set_def_mod, only : obs_set_def_type, init_obs_set_def, add_obs
use set_def_list_mod, only : set_def_list_type, init_set_def_list, &
                             add_to_list, write_set_def_list, read_set_def_list
use    utilities_mod, only : open_file
use  assim_model_mod, only : static_init_assim_model, get_model_size, get_state_meta_data
use     location_mod, only : location_type, get_location, set_location
use        model_mod, only : get_close_states_devel, get_state_meta_data

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_def_type) :: obs_def
type(set_def_list_type) :: set_def_list
type(obs_set_def_type) :: obs_set_def
type(location_type) :: loc

integer :: max_sets, num_obs, model_size
integer :: i, j, obs_set_def_index, iunit, var_type, input_type
character(len = 129) :: file_name
real(r8) :: location(3), err_var

! For test of get_close_state_devel
real(r8) :: o_lon, o_lat, o_lev, radius, lon_lat_lev(3)
type(location_type) :: o_loc
integer :: var_type
integer, allocatable :: indices(:)
real(r8), allocatable :: dist(:)

! Get output filename
write(*, *) 'Input the filename for output of observation set_def_list'
read(*, *) file_name

! With potential for perfect model observations, must do static_init of model
call static_init_assim_model()
model_size = get_model_size()

! Initialize a set_def_list; for now need to specify max_sets
max_sets = 1
set_def_list = init_set_def_list(max_sets)


! For now, assume that no more than one observation for each state variable
num_obs = model_size

! Initialize the obs_set_def
obs_set_def = init_obs_set_def(num_obs)


! Loop for different variable types
do i = 1, 1000
   write(*, *) 'input integer type index, negative for finished'
   read(*, *) input_type
   if(input_type < 0) goto 21
   write(*, *) 'input error variance for var type ', input_type
   read(*, *) err_var
   do j = 1, model_size
      call get_state_meta_data(j, loc, var_type)
      if(var_type == input_type) then
         ! Add this obs_def to the set
         write(*, *) 'adding var index ', j
         obs_def = init_obs_def(j, err_var)
         call add_obs(obs_set_def, obs_def)
      endif
   end do
 
end do


! Insert this obs_set_def into the list
21   obs_set_def_index = add_to_list(set_def_list, obs_set_def)


! Output the set_def_list
iunit = open_file(file_name, action = 'write')
call write_set_def_list(iunit, set_def_list)
close(iunit)

! Read and rewrite as test
!iunit = open_file(file_name)
!set_def_list = read_set_def_list(iunit)
!close(iunit)

!iunit = open_file(file_name)
!call write_set_def_list(iunit, set_def_list)
!close(iunit)

end program identity_obs_set_def
