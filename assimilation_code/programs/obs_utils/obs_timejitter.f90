! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program obs_timejitter

use        types_mod, only : r8
use    utilities_mod, only : open_file, close_file, &
                             initialize_utilities, finalize_utilities
use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian
use      obs_def_mod, only : obs_def_type, get_obs_def_time, set_obs_def_time
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
   get_num_obs, init_obs_sequence, get_first_obs, write_obs_seq, set_copy_meta_data, &
   get_obs_def, set_obs_def, get_next_obs, insert_obs_in_seq, init_obs, &
   assignment(=), static_init_obs_sequence, get_num_copies, get_num_qc, &
   get_copy_meta_data, get_qc_meta_data, set_qc_meta_data
use time_manager_mod, only : time_type, set_time, get_time, interactive_time, &
                             operator(*), operator(+), operator(-), operator(<)
use        model_mod, only : static_init_model

implicit none

type(obs_sequence_type) :: seq, seq_in
type(obs_type)          :: obs, next_obs, new_obs, last_obs
type(obs_def_type)      :: obs_def
character(len = 129)    :: file_name
logical                 :: is_there_one, is_this_last
type(time_type)         :: ob_time, last_time
integer                 :: seconds, days, i, num_obs, num_copies, num_qc
type(random_seq_type)   :: random_seq
real(r8)                :: perturbation_amplitude, delseconds  ! size of perturbations



! Record the current time, date, etc. to the logfile
call initialize_utilities('obs_timejitter')

! Call the underlying model's static initialization for calendar info
call static_init_model()

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Initialize the random number sequence
call init_random_seq(random_seq, 1)

! Write the sequence to a file
write(*, *) 'Input filename for network definition sequence (usually  set_def.out  )'
read(*, *) file_name
call read_obs_seq(file_name, 0, 0, 0, seq_in)

! Find out how many obs there are
num_obs = get_num_obs(seq_in)

! Initialize the obs_type variables
num_copies = get_num_copies(seq_in)
num_qc     = get_num_qc(seq_in)
call init_obs(     obs, num_copies, num_qc)
call init_obs(next_obs, num_copies, num_qc)
call init_obs( new_obs, num_copies, num_qc)
call init_obs(last_obs, num_copies, num_qc)

! Get the time spread information 

write(*, *) 'Input the spread in time in seconds'
read(*, *) perturbation_amplitude

! Initialize the output sequence
call init_obs_sequence(seq, num_copies, num_qc, num_obs)
do i = 1, num_copies
   call set_copy_meta_data(seq, i, get_copy_meta_data(seq_in, i))
end do
do i = 1, num_qc
   call set_qc_meta_data(seq, i, get_qc_meta_data(seq_in, i))
end do

last_time = set_time(0, 0)
is_there_one = get_first_obs(seq_in, obs)

do i = 1, num_obs
   new_obs = obs

   ! Set the time
   call get_obs_def(new_obs, obs_def)
   ob_time = get_obs_def_time(obs_def)

   ! jitter here
   call get_time(ob_time, seconds, days)
   delseconds = int(random_gaussian(random_seq, real(seconds, r8), perturbation_amplitude))
   if (delseconds > 0) then
      ob_time = ob_time + set_time(int(delseconds), 0)
   else
      if (days == 0 .and. delseconds > seconds) then
         ob_time = set_time(0, 0)
      else
         ob_time = ob_time - set_time(-int(delseconds), 0)
      endif
   endif

   call set_obs_def_time(obs_def, ob_time) 
   call set_obs_def(new_obs, obs_def)

   ! Insert it in the new sequence.  If time has not been
   ! moved back before the time of the last obs, use it
   ! as the start of the insert to save search time.
   ! Otherwise, start at the first obs in the sequence.
   if (i == 1 .or. ob_time < last_time) then
      call insert_obs_in_seq(seq, new_obs)
   else
      call insert_obs_in_seq(seq, new_obs, last_obs) 
   endif
   last_obs = new_obs
   last_time = ob_time

   ! Find the next observation in the input set
   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if(.not. is_this_last) obs = next_obs
end do

write(*, *) 'What is output file name for sequence (  obs_seq.in   is recommended )'
read(*, *) file_name
call write_obs_seq(seq, file_name)

! Clean up
call finalize_utilities()

end program obs_timejitter

