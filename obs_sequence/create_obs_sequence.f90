program create_obs_sequence

use utilities_mod, only : open_file, close_file
use obs_sequence_mod, only : obs_sequence_type, interactive_obs, &
   write_obs_seq, interactive_obs_sequence, static_init_obs_sequence

type(obs_sequence_type) :: seq
character(len = 129) :: file_name

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Interactive creation of an observation sequence
seq = interactive_obs_sequence()

! Write the sequence to a file
write(*, *) 'Input filename for sequence'
read(*, *) file_name
call write_obs_seq(seq, file_name)

end program create_obs_sequence
