program create_fixed_network_seq

use types_mod, only : r8
use utilities_mod, only : open_file, close_file
use obs_def_mod, only : obs_def_type, get_obs_def_time, set_obs_def_time
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
   get_num_obs, init_obs_sequence, get_first_obs, write_obs_seq, set_copy_meta_data, &
   get_obs_def, set_obs_def, append_obs_to_seq, get_next_obs, insert_obs_in_seq, init_obs, &
   assignment(=), static_init_obs_sequence
use time_manager_mod, only : time_type, operator(*), operator(+), set_time

implicit none

type(obs_sequence_type) :: seq, seq_in
type(obs_type) :: obs, next_obs, new_obs
type(obs_def_type) :: obs_def
integer :: unit_num
character(len = 129) :: file_name
logical :: is_there_one, is_this_last
type(time_type) :: ob_time, init_time, this_time, period
integer :: seconds, days, i, j, network_size, option, num_times
character(len = 129) :: in_string

! Initialize the obs_sequence module
call static_init_obs_sequence

! Initialize the obs_type variables
call init_obs(obs, 0, 0)
call init_obs(next_obs, 0, 0)
call init_obs(new_obs, 0, 0)

! Write the sequence to a file
write(*, *) 'Input filename for network definition sequence '
read(*, *) file_name
call read_obs_seq(file_name, 0, 0, 0, seq_in)

! Find out how many obs there are
network_size = get_num_obs(seq_in)

! Get the time information 

20 write(*, *) 'To input a regularly repeating time sequence enter 1'
write(*, *) 'To enter an irregular list of times enter 2'
read(*, *) option

! Should also allow both regular and irregular for same set,
! but too much work for now

if(option == 1) then

   write(*, *) 'Input number of observation times in sequence'
   read(*, *) num_times

   write(*, *) 'Input initial time in sequence in days and seconds'
   read(*, *) days, seconds
   init_time = set_time(seconds, days)

   write(*, *) 'Input period of obs in sequence in days and seconds'
   read(*, *) days, seconds
   period = set_time(seconds, days)

! Initialize the output sequence
   call init_obs_sequence(seq, 0, 0, network_size * num_times)

   do j = 1, num_times
      write(*, *) j
      ob_time = init_time + (j - 1) * period

      is_there_one = get_first_obs(seq_in, obs)
      
      do i = 1, network_size
         new_obs = obs
! Set the time
         call get_obs_def(new_obs, obs_def)
         call set_obs_def_time(obs_def, ob_time) 
         call set_obs_def(new_obs, obs_def)
! Append it to the sequence 
         call append_obs_to_seq(seq, new_obs)
! Find the next observation in the input set
         call get_next_obs(seq_in, obs, next_obs, is_there_one)
         obs = next_obs
      end do

   enddo

!-------------------------------------------------------------------------

else if(option == 2) then

   ! The irregular input section does minimal error checking.
   ! non-monotonic input times will cause an error in add_obs_set().
   write(*, *) 'Input an upper bound on the number of observation times'
   read(*, *) num_times

! Initialize the output sequence
   call init_obs_sequence(seq, 0, 0, network_size * num_times)

   IRREGULAR : do

      write(*, *) 'Input time in days and seconds, negative days if finished with this set'
      read(*, *) days, seconds

      if ( days < 0 ) exit IRREGULAR         ! Done with this set.

      this_time = set_time(seconds, days)    ! Set the time

      ! Input this on a long list for later sorting
      ! (fixed storage for now should be changed)

      ob_time = this_time

! Put all the observations in the output sequence with this time      
      is_there_one = get_first_obs(seq_in, obs)
      
      do i = 1, network_size
         new_obs = obs
! Set the time
         call get_obs_def(new_obs, obs_def)
         call set_obs_def_time(obs_def, ob_time) 
         call set_obs_def(new_obs, obs_def)
! Append it to the sequence 
         call append_obs_to_seq(seq, new_obs)
! Find the next observation in the input set
         call get_next_obs(seq_in, obs, next_obs, is_there_one)
         obs = next_obs
      end do

   enddo IRREGULAR

else
   write(*, *) 'option must be 1 or 2, try again'
   goto 20
endif

write(*, *) 'What is output file name for sequence'
read(*, *) file_name
call write_obs_seq(seq, file_name)


end program create_fixed_network_seq
