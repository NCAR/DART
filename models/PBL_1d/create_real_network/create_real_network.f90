! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_real_network_seq

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! JPH
! This code originated from create_fixed_network.  It uses module_wrf to get
! obs from smos file, with file, date, and interval controlled via the wrf1d
! namelist.  Note that an obs_def is still required to control which
! obs are actually written out.  Normally, this would be created with 
! create_obs_sequence.  This would be run in place of both create_fixed_network
! and perfect_model_obs.

use        types_mod, only : r8, missing_r8, missing_i
use    utilities_mod, only : timestamp, register_module, open_file, &
                             close_file, find_namelist_in_file, &
                             error_handler, check_namelist_read, &
                             initialize_utilities, E_ERR
use     obs_kind_mod, only : assimilate_this_obs_kind, evaluate_this_obs_kind
use      obs_def_mod, only : obs_def_type, get_obs_def_time, set_obs_def_time,&
                             get_obs_kind, get_obs_name
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                             get_num_obs, init_obs_sequence, get_first_obs, &
                             write_obs_seq, set_copy_meta_data, &
                             get_obs_def, set_obs_def, append_obs_to_seq, &
                             get_next_obs, insert_obs_in_seq, init_obs, &
                             assignment(=), static_init_obs_sequence, &
                             get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, &
                             set_qc_meta_data, read_obs_seq_header, &
                             set_obs_values, set_qc, get_qc
use time_manager_mod, only : time_type, operator(*), operator(+), set_time, &
                             set_date, increment_time, get_time, print_time, &
                             operator(==), operator(/), operator(<), operator(-)
use        model_mod, only : static_init_model, real_obs_period, start_real_obs
use        module_wrf, only : static_init_wrf, init_wrf, nt_f_smos, &
                              start_year_f, start_month_f,start_day_f, &
                              start_hour_f, start_minute_f, &
                              start_forecast, interval_smos, &
                              init_f_type, u10_init_f, v10_init_f, &
                              q2_init_f, t2_init_f, forecast_length


implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"


type(obs_sequence_type) :: seq, seq_in, seq_out
type(obs_type)          :: obs, next_obs, new_obs
type(obs_def_type)      :: obs_def
character(len = 129)    :: file_name, obs_seq_in_file_name
logical                 :: is_there_one, is_this_last
type(time_type)         :: ob_time, init_time, this_time
type(time_type)         :: obs_seq_period, obs_list_period
type(time_type),dimension(:), allocatable :: obs_list_time
type(time_type),dimension(:), allocatable :: obs_seq_time
type(time_type)         :: start_seq_time, flen_time, end_time
integer                 :: seconds, days, i, j, network_size,  num_times, num_copies, num_qc
integer                 :: obs_seq_file_id, iunit, io
integer                 :: cnum_copies, cnum_qc, cnum_obs, cnum_max
integer                 :: additional_qc, additional_copies
integer                 :: last_key_used, time_step_number
integer                 :: num_obs, obs_kind_ind
real(r8)                :: this_obs_val, this_qc_val
real(r8), dimension(:), allocatable :: obs_vals, qc_vals, qc_sequence
logical                 :: assimilate_this_ob, evaluate_this_ob, pre_I_format
character(len=129)      :: copy_meta_data(2), qc_meta_data, obs_seq_read_format
integer                 :: wrf_rnd_seed = -1

! Record the current time, date, etc. to the logfile
call initialize_utilities('Create_real_network_seq')
call register_module(source,revision,revdate)

! The only necessary namelist variables come from the model

! Call the underlying model's static initialization for calendar info
call static_init_model()

! Initialize the obs_sequence module
call static_init_obs_sequence

! fail if we are not initializing from OBS (this could be easily modified
! to get values from WRF, and may come in handy later!
if ( init_f_type == 'WRF' ) then
  call error_handler(E_ERR, 'create_real_network', &
     'CANNOT PRODUCE OBS SEQUENCE FROM WRF OUTPUT YET', source, revision, revdate)
endif

! Write the sequence to a file
write(*, *) 'Input filename for network definition sequence (usually  set_def.out  )'
read(*, *) file_name
call read_obs_seq(file_name, 0, 0, 0, seq_in)

! Find out how many obs there are
network_size = get_num_obs(seq_in)

! Initialize the obs_type variables
num_copies = get_num_copies(seq_in)
num_qc = get_num_qc(seq_in)
call init_obs(obs, num_copies, num_qc)
call init_obs(next_obs, num_copies, num_qc)
call init_obs(new_obs, num_copies, num_qc)

! set init time and period, including increment for forecast start and
! increment to the proper time of day.  One might want to start the assimilation
! later then when constrained to start the obs list.
   init_time = set_date(start_year_f, start_month_f, start_day_f, &
                        start_hour_f, start_minute_f, 0)
   init_time = increment_time(init_time,start_forecast,0)
   flen_time = set_time(forecast_length,0) 
   end_time  = init_time + flen_time

   call get_time(init_time,seconds,days)

   start_seq_time = increment_time(init_time,start_real_obs,0)
   obs_seq_period = set_time(real_obs_period, 0)
   obs_list_period = set_time(interval_smos, 0)

   num_times = (end_time - start_seq_time) / set_time(real_obs_period,0) + 1

   ! time information comes from the wrf1d_namelist.input
   ! only supports regularly-repeating obs right now
   allocate(obs_seq_time(num_times))
   allocate(obs_list_time(nt_f_smos))

   ! associate a time with each obs in the input list
   do j = 1, nt_f_smos
      obs_list_time(j) = init_time + (j - 1) * obs_list_period
   enddo

   ! Initialize the output sequence
   call init_obs_sequence(seq, num_copies, &
      num_qc, network_size * num_times)

   ! Get the metadata (might want a call in obs_sequence to do this)
   do i = 1, num_copies
      call set_copy_meta_data(seq, i, get_copy_meta_data(seq_in, i))
   end do
   do i = 1, num_qc
      call set_qc_meta_data(seq, i, get_qc_meta_data(seq_in, i))
   end do

   ! while looping through the times, generate a list of obs times
   ! and qc values
   do j = 1, num_times
      write(*, *) j
      ob_time = start_seq_time + (j - 1) * obs_seq_period
      obs_seq_time(j) = ob_time
      call print_time(obs_seq_time(j))

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
         call get_next_obs(seq_in, obs, next_obs, is_this_last)
         if(.not. is_this_last) obs = next_obs
      end do

   enddo

!-------------------------------------------------------------------------
! write to a temporary file for ingestion into the next block
file_name = 'real_obs_seq.in'

call write_obs_seq(seq, file_name)

! Clean up
call timestamp(string1=source,string2=revision,string3=revdate,pos='end')

!-------------------------------------------------------------------------
! Now the part that replaces perfect_model_obs.  There are some 
! assumptions in here about what type of obs we are ingesting:
! 1.  pressure and vapor pressure are used to get mixing ratio
! 2.  T and winds are in correct units (K and m/s)

call init_wrf(wrf_rnd_seed)

!do i = 1,num_times
!  print*,t2_init_f(i),u10_init_f(i),v10_init_f(i),q2_init_f(i)
!enddo

obs_seq_in_file_name = file_name

call read_obs_seq_header(obs_seq_in_file_name, cnum_copies, cnum_qc, &
                         cnum_obs, cnum_max, obs_seq_file_id, &
                         obs_seq_read_format, pre_I_format, &
                         close_the_file = .true.)

! First two copies of output will be truth and observation;
! Will overwrite first two existing copies in file if there are any
! Note that truth=obs for this case of real obs
additional_copies = 2 - cnum_copies
if(additional_copies < 0) additional_copies = 0

! currently no need for additional qc field
additional_qc = 0

! Just read in the definition part of the obs sequence; expand to include
! observation and truth field
call read_obs_seq(obs_seq_in_file_name, additional_copies, additional_qc, &
                  0, seq)

! Initialize an obs type variable
call init_obs(obs, cnum_copies + additional_copies, cnum_qc + additional_qc)

! Need metadata for added qc field (here in case needed later)
if(additional_qc == 1) then
   qc_meta_data = 'Quality Control'
   call set_qc_meta_data(seq, 1, qc_meta_data)
endif

time_step_number = 0
num_qc = get_num_qc(seq)
num_copies = get_num_copies(seq)
num_obs = get_num_obs(seq)

! init output obs sequence
call init_obs_sequence(seq_out, num_copies, num_qc, num_obs)
call init_obs(obs, num_copies, num_qc)
call init_obs(next_obs, num_copies, num_qc)
call init_obs(new_obs, num_copies, num_qc)

! Need space to put in the obs_values in the sequence;
copy_meta_data(1) = 'observations'
copy_meta_data(2) = 'truth'
call set_copy_meta_data(seq_out, 1, copy_meta_data(1))
call set_copy_meta_data(seq_out, 2, copy_meta_data(2))
do i = 1, num_qc
  call set_qc_meta_data(seq_out, i, get_qc_meta_data(seq, i))
end do

! simply look through obs one-by-one and pull from the proper vector
allocate(obs_vals(num_copies), qc_vals(num_qc))
allocate(qc_sequence(num_obs))

is_there_one = get_first_obs(seq, obs)
if ( is_there_one ) then
   do i = 1, num_obs
     new_obs = obs

     ! Set the time
     call get_obs_def(new_obs, obs_def)
     ob_time = get_obs_def_time(obs_def)
     obs_kind_ind = get_obs_kind(obs_def)
     assimilate_this_ob = assimilate_this_obs_kind(obs_kind_ind)
     evaluate_this_ob = evaluate_this_obs_kind(obs_kind_ind)

     this_obs_val = get_obs_from_input(ob_time,obs_kind_ind,num_times)
     this_qc_val  = get_qc_from_obs(obs_kind_ind,this_obs_val)
     if ( num_qc > 0 ) then
       call get_qc(new_obs,qc_sequence(i:i),1)
       this_qc_val = max(this_qc_val,qc_sequence(i))
     endif

     ! for input, all copies are the same
     obs_vals = this_obs_val
     call set_obs_values(new_obs,obs_vals)
     qc_vals = this_qc_val

     if ( num_qc > 0 ) then
       call set_qc(new_obs,qc_vals)
     endif
       
     call set_obs_def(new_obs, obs_def)

     ! Append it to the sequence
     call append_obs_to_seq(seq_out, new_obs)

     ! Find the next observation in the input set
     call get_next_obs(seq, obs, next_obs, is_this_last)
     if(.not. is_this_last) obs = next_obs
  
   end do ! obs

file_name = 'real_obs_seq.out'

call write_obs_seq(seq_out, file_name)
stop


else 

   print*, "could not find any obs in the input sequence"

endif


!--------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------


real(r8) function get_obs_from_input(ob_time,obs_kind_in,num_times)

implicit none

type(time_type), intent(in)      :: ob_time
integer, intent(in)              :: obs_kind_in, num_times

integer                          :: seconds, days, i
integer                          :: this_time_ind
real(r8)                         :: obs_val

get_obs_from_input = missing_r8

this_time_ind = missing_i
do i = 1, nt_f_smos
  if ( obs_list_time(i) == ob_time ) this_time_ind = i
enddo

if ( this_time_ind == missing_i ) return

select case ( trim(get_obs_name(obs_kind_in)) )
  case ('METAR_U_10_METER_WIND')
    obs_val = u10_init_f(this_time_ind)
  case ('METAR_V_10_METER_WIND')
    obs_val = v10_init_f(this_time_ind)
  case ('METAR_TEMPERATURE_2_METER')
    obs_val = t2_init_f(this_time_ind)
  case ('METAR_SPECIFIC_HUMIDITY_2_METER')
    obs_val = q2_init_f(this_time_ind)
  case default
    return
end select

get_obs_from_input = obs_val

end function get_obs_from_input

!--------------------------------------------------------------

real(r8) function get_qc_from_obs(obs_kind_in,obs_val)
! simple gross error check on qc

implicit none

integer, intent(in)              :: obs_kind_in
real(r8), intent(in)             :: obs_val

get_qc_from_obs = 0.0_r8
if ( obs_val == missing_r8 ) then
  get_qc_from_obs = 9.0_r8
  return
end if

! no real qc yet
select case ( trim(get_obs_name(obs_kind_in)) )
  case default
    return
end select

end function get_qc_from_obs

end program create_real_network_seq
