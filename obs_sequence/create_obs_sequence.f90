program create_obs_sequence

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

! Program to build a simple obs_sequence file for use in testing filters. Reads
! in a previously created set_def_list that has spatial information about 
! available obs_set_defs. This then creates a time-ordered set of observations.
! Easiest way to do this is to simply specify a set of times for each observation
! set in the definition. Initially allow regularly spaced in time or ... Will
! eventually need to embelish to allow arbitrary initial time, time read from
! files, etc.
!
! TJH Sun Jun 30 15:20:04 MDT 2002: the sets must be input in a monotonic ascending
!    fashion. At some point we need to sort these sets.

use types_mod
use obs_sequence_mod, only : obs_sequence_type, init_obs_sequence, &
   add_obs_set, write_obs_sequence, associate_def_list, read_obs_sequence
use set_def_list_mod, only : set_def_list_type, get_num_sets_in_list, read_set_def_list
use obs_set_mod, only : obs_set_type, init_obs_set, set_obs_set_time
use time_manager_mod, only : time_type, set_time, get_time, operator(*), operator(+)
use utilities_mod, only : open_file
use sort_mod, only : index_sort

implicit none

! Everybody needs to know these ... TJH Feb 10, 2003
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
type(time_type)    :: ob_time(max_obs), init_time, period, this_time

! Kluge for sorting time types
real(r8) :: time_secs(max_obs)

integer :: in_unit, out_unit, in_unit2, out_unit2
integer :: i, j, obs_set_def_index, indx = 0, days, seconds, option
integer :: num_obs_set_defs, num_obs
character(len = 129) :: file_name

! Change output to diagnostic output block ...
write(*,*)'create_obs_sequence attributes:'
write(*,*)'   ',trim(adjustl(source))
write(*,*)'   ',trim(adjustl(revision))
write(*,*)'   ',trim(adjustl(revdate))
write(*,*)'   '

! Get file name for input set_def_list
write(*, *) 'What is name of set_def_list? [set_def.out]'
read(*, *) file_name
in_unit = open_file(file_name, action = 'read')

set_def_list = read_set_def_list(in_unit)     ! Read in the set_def_list

seq = init_obs_sequence(max_obs, 0)           ! Initialize a sequence

call associate_def_list(seq, set_def_list)    ! Put the obs_set_def into the sequence

! Loop through the individual defs (easy with no subsets)

num_obs_set_defs = get_num_sets_in_list(set_def_list)

do i = 1, num_obs_set_defs

   write(*, *) 'Setting times for obs_def ', i
   ! Would want to output some information about this obs_def, too, 
   ! should come from obs_def_set level
  
   ! Input of time should be pushed into time manager eventually, but I'm lazy 

20 continue

   write(*, *) 'To input a regularly repeating time sequence enter 1'
   write(*, *) 'To enter an irregular list of times enter 2'
   read(*, *) option

   ! Should also allow both regular and irregular for same set, 
   ! but too much work for now

   if(option == 1) then

      write(*, *) 'Input number of observations in sequence'
      read(*, *) num_obs

      write(*, *) 'Input time of initial ob in sequence in days and seconds'
      read(*, *) days, seconds
      init_time = set_time(seconds, days)

      write(*, *) 'Input period of obs in days and seconds'
      read(*, *) days, seconds
      period = set_time(seconds, days)

      ! Loop to put these obs in the list

      do j = 1, num_obs
         indx = indx + 1
         if(indx > max_obs) then
            write(*, *) 'Number of observations exceeds max_obs'
            stop
         endif
         set_index(indx) = i
         ob_time(indx) = init_time + (j - 1) * period
      enddo

   else if(option == 2) then

      ! The irregular input section does minimal error checking.
      ! non-monotonic input times will cause an error in add_obs_set().

      IRREGULAR : do 

         write(*, *) 'Input time in days and seconds, negative days if finished with this set'
         read(*, *) days, seconds

         if ( days < 0 ) exit IRREGULAR         ! Done with this set.

         this_time = set_time(seconds, days)    ! Set the time

         ! Input this on a long list for later sorting 
         ! (fixed storage for now should be changed)

         indx = indx + 1
         set_index(indx) = i
         if(indx > max_obs) then
            write(*, *) 'Number of observations exceeds max_obs'
            stop
         endif
         ob_time(indx) = this_time 

      enddo IRREGULAR

   else 
      write(*, *) 'option must be 1 or 2, try again'
      goto 20
   endif
enddo

! Next would have to sort the list;  This is a quick fix ; should have
! time sorting routines
! May be that this is okay with most floating point precisions???
do i = 1, indx
   call get_time(ob_time(i), seconds, days)
   time_secs(i) = days * 86400 + seconds
end do
call index_sort(time_secs, sort_ind, indx) 

! Now generate a long obs_sequence with regular occurences of the obs_set_def

do i = 1, indx
   obs_set = init_obs_set(set_def_list, set_index(sort_ind(i)), 0)
   call get_time(ob_time(sort_ind(i)), seconds, days)
   write(*, *) 'time ', i, ' is ', seconds, days
   call set_obs_set_time(obs_set, ob_time(sort_ind(i)))
   call add_obs_set(seq, obs_set)       ! Put this obs_set into the sequence
enddo

! Output the obs_sequence to a file

write(*, *) 'Input file name for output of obs_sequence? [obs_seq.in]'
read( *, *) file_name
out_unit = open_file(file_name, action = 'write')
call write_obs_sequence(out_unit, seq)

! Temporary read and write test
! THIS DOESN'T WORK ON PGF90 COMPILER: WHY???
!close(out_unit)
!write(*, *) 're-reading from ', file_name
!in_unit2 = open_file(file_name)
!seq = read_obs_sequence(in_unit2)

!write(*, *) 'rewriting to  garb_out' 
!out_unit2 = open_file('garb_out')
!call write_obs_sequence(out_unit2, seq)

end program create_obs_sequence
