program new_create_obs_sequence

use types_mod, only : r8
use utilities_mod, only : open_file, close_file
use obs_sequence_mod, only : init_obs_sequence, obs_sequence_type, &
   insert_obs_in_seq, interactive_obs, write_obs_seq, &
   interactive_obs_sequence, read_obs_seq, obs_type, get_first_obs, &
   get_obs_def, get_num_obs, get_next_obs, init_obs, assignment(=), static_init_obs_sequence
use obs_def_mod, only : get_obs_def_location, obs_def_type
use location_mod, only : location_type, get_location

type(obs_sequence_type) :: seq, seq2
type(obs_type) :: obs, next_obs
type(obs_def_type) :: obs_def
type(location_type) :: location
integer :: unit_num
character(len = 129) :: file_name
logical :: is_there_one, is_this_last

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Initialize the obs types
call init_obs(obs, 0, 0)
call init_obs(next_obs, 0, 0)

! Interactive creation of an observation sequence
seq = interactive_obs_sequence()

! Write the sequence to a file
write(*, *) 'Input filename for sequence'
read(*, *) file_name
call write_obs_seq(seq, file_name)

end program new_create_obs_sequence
