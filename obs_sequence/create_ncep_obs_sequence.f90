! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_ncep_obs_sequence

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! Program to build a simple ncep_obs_sequence file for use in testing filters. Reads
! in a previously created set_def_list that has spatial information about 
! available obs_set_defs. This then creates a time-ordered set of observations.
!
use        types_mod, only : r8
use obs_sequence_mod, only : obs_sequence_type, init_obs_sequence, &
   add_obs_set, write_obs_sequence, associate_def_list, read_obs_sequence
use set_def_list_mod, only : set_def_list_type, get_num_sets_in_list, read_set_def_list
use      obs_set_mod, only : obs_set_type, init_obs_set, set_obs_set_time
use time_manager_mod, only : time_type, set_time, get_time, operator(*), operator(+)
use    utilities_mod, only : open_file
use         sort_mod, only : index_sort
use time_manager_mod, only : set_date, set_calendar_type

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"


type(obs_sequence_type) :: seq
type(set_def_list_type) :: set_def_list
type(obs_set_type)      :: obs_set

! Fixed storage leads to use of parameter for array, fix this
integer, parameter :: max_obs = 10000
integer            :: set_index(max_obs), sort_ind(max_obs)
type(time_type)    :: time(max_obs), init_time, period, this_time

! Kluge for sorting time types
real(r8) :: time_secs(max_obs)

integer :: in_unit, out_unit, in_unit2, out_unit2
integer :: i, j, obs_set_def_index, index = 0, days, seconds, option
integer :: num_obs_set_defs, num_obs
character(len = 129) :: file_name

integer :: obs_year, obs_month, obs_day, obs_hour, obs_min, obs_sec
integer :: calender_type, calender_type_test
type(time_type)    :: obs_time

! Change output to diagnostic output block ...
write(*,*)'create_obs_sequence attributes:'
write(*,*)'   ',trim(adjustl(source))
write(*,*)'   ',trim(adjustl(revision))
write(*,*)'   ',trim(adjustl(revdate))
write(*,*)'   '

! Get file name for input set_def_list
  file_name = 'set_def.out'
write(*, *) 'name of set_def_list is ', file_name
in_unit = open_file(file_name, action = 'read')

set_def_list = read_set_def_list(in_unit)     ! Read in the set_def_list
seq = init_obs_sequence(max_obs, 0)           ! Initialize a sequence
call associate_def_list(seq, set_def_list)    ! Put the obs_set_def into the sequence

num_obs_set_defs = get_num_sets_in_list(set_def_list)


    calender_type = 3
    call set_calendar_type(calender_type)

    print*, 'input calendar type for B-grid test: 1=Gregrian, 0=Arbitrary'
    read(*,*) calender_type_test

    open(800, file='/home/hliu/DART/ncep_obs/ncepobs.input', form='formatted') 
    read(800, *) obs_year, obs_month, obs_day            
    close(800)          
                                                                      
    print*, 'obs_date in ncepobs = ', obs_year, obs_month, obs_day

      obs_hour = 0     !! begin from 00Z data set
      obs_min  = 0
      obs_sec  = 0

      obs_time = set_date(obs_year,obs_month,obs_day,obs_hour,obs_min,obs_sec)

      call get_time(obs_time, seconds, days)
      write(*,*)  'obs time', seconds, days

! Loop through the individual defs 
 do i = 1, num_obs_set_defs
  
!     Relative time of the current obs set in sequence in days and seconds
      days = 0
      seconds = 3600* i
      init_time = set_time(seconds, days)

!     add obs time to the init time
      if(calender_type_test .ne. 0) then
      init_time = init_time + obs_time
      endif
      
      call get_time(init_time, seconds, days)
      write(*,*)  'obs + init time', seconds, days

!     put the obs in the list
      index = index + 1
      set_index(index) = i
      time(index) = init_time 

 enddo

! Next sort the list;  
 do i = 1, index
   call get_time(time(i), seconds, days)
   time_secs(i) = days * 86400 + seconds
 end do
   call index_sort(time_secs, sort_ind, index) 

! Now generate a long obs_sequence with regular occurences of the obs_set_def

 do i = 1, index
   obs_set = init_obs_set(set_def_list, set_index(sort_ind(i)), 0)
   call get_time(time(sort_ind(i)), seconds, days)
   write(*, *) 'time ', i, ' is ', seconds, days
   call set_obs_set_time(obs_set, time(sort_ind(i)))
   call add_obs_set(seq, obs_set)       ! Put this obs_set into the sequence
 enddo

! Output the obs_sequence to a file

  file_name = 'obs_seq.in'
  out_unit = open_file(file_name, action = 'write')
  call write_obs_sequence(out_unit, seq)


end program create_ncep_obs_sequence
