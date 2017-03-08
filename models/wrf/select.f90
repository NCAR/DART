! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
PROGRAM select

use        types_mod, only : r8, metadatalength
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             register_module, error_handler, E_MSG
use obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, &
                             insert_obs_in_seq, get_first_obs, get_next_obs, &
                             write_obs_seq, &
                             assignment(=), &
                             init_obs, static_init_obs_sequence, &
                             get_num_obs, get_num_copies, get_num_qc, &
                             get_obs_def, read_obs_seq, &
                             get_copy_meta_data, set_copy_meta_data, &
                             get_qc_meta_data, set_qc_meta_data
use     obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT, &
                             RADIOSONDE_V_WIND_COMPONENT, &
                             RADIOSONDE_SURFACE_PRESSURE, &
                             RADIOSONDE_TEMPERATURE, &
                             RADIOSONDE_SPECIFIC_HUMIDITY
use      obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, &
                             get_obs_def_time, get_obs_def_location
use     location_mod, only : location_type, get_location
use time_manager_mod, only : time_type, operator(/=), get_time, print_time, &
                             set_calendar_type, GREGORIAN

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq, real_seq
type(obs_type)          :: obs, prev_obs, real_obs
type(obs_def_type)      :: real_obs_def
type(location_type)     :: location
type(time_type)         :: time, prev_time

integer                 :: seconds, days, kind, delta
real(r8), dimension(3)  :: loc

integer           :: i
integer           :: num_obs, num_copies, num_qc, real_seq_num_obs

character(len = metadatalength) :: meta_data

character(len = 129) :: out_file_name = 'obs_seq.out', &
                         in_file_name = 'obs_seq.in'
integer              :: calendar_type         = GREGORIAN

!------------------------------------------------------------------------------

logical  :: is_this_last, is_there_one

call initialize_utilities('Select')
call register_module(source, revision, revdate)

call set_calendar_type(calendar_type)

call static_init_obs_sequence()

! -------------------------------------------------------------------
! Initialize the counters:

num_obs = 0

delta = 10800

call read_obs_seq(in_file_name, 0, 0, 0, real_seq)

num_copies = get_num_copies(real_seq)
real_seq_num_obs = get_num_obs(real_seq)
num_qc = get_num_qc(real_seq)

! Initialize an obs_sequence structure
call init_obs_sequence(seq, num_copies, num_qc, real_seq_num_obs)

do i = 1, num_copies
   meta_data = get_copy_meta_data(real_seq, i)
   call set_copy_meta_data(seq, i, meta_data)
enddo

do i = 1, num_qc
   meta_data = get_qc_meta_data(real_seq, i)
   call set_qc_meta_data(seq, i, meta_data)
enddo

call init_obs(obs, num_copies, num_qc)

print*,'Number of observations in original file: ',real_seq_num_obs

do i = 1, real_seq_num_obs

   if(i > 1) then

      call get_next_obs(real_seq, real_obs, real_obs, is_this_last)
      prev_time = time

   else

      is_there_one = get_first_obs(real_seq, real_obs)

   endif

   call get_obs_def(real_obs, real_obs_def)

   time = get_obs_def_time(real_obs_def)

   call get_time(time, seconds, days)

!!$   if(time /= prev_time) call print_time(time)

   location = get_obs_def_location(real_obs_def)
   loc = get_location(location)

   kind = get_obs_def_type_of_obs(real_obs_def)

   if (  &
!!$         seconds >= (86401 - delta) .or. &
!!$         seconds <= delta           .or. &
!!$        (seconds >= (43201 - delta) .and. seconds <= (43200 + delta)) .and. &
        (loc(2) >= 0.0_r8) .and. &  ! NH
        ((kind == RADIOSONDE_U_WIND_COMPONENT) .or. &
         (kind == RADIOSONDE_V_WIND_COMPONENT) .or. &
         (kind == RADIOSONDE_TEMPERATURE)      .or. &
         (kind == RADIOSONDE_SPECIFIC_HUMIDITY))    &
        ) then

      obs = real_obs

      num_obs = num_obs + 1

      if(num_obs == 1) then
         call insert_obs_in_seq(seq, obs)
      else
         call insert_obs_in_seq(seq, obs, prev_obs)
      endif

      prev_obs = obs

   endif

enddo

write(unit=*, fmt='(5x,a,i6,a)') &
     'Total number of observations:  ', num_obs

call write_obs_seq(seq, out_file_name)

call error_handler(E_MSG,'select','FINISHED select.')
call error_handler(E_MSG,'select','Finished successfully.',source,revision,revdate)
call finalize_utilities()
 
END PROGRAM select

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
