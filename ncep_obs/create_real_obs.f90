program create_real_obs

use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence

use real_obs_mod, only : real_obs_sequence

type(obs_sequence_type) :: seq
character(len = 129) :: file_name = 'obs_seq.out'

! Initialize the obs_sequence module
!
  call static_init_obs_sequence()

! Creation of an real observation sequence
!
  seq = real_obs_sequence()

! Write the sequence to a file
!
  call write_obs_seq(seq, file_name)

end program create_real_obs
