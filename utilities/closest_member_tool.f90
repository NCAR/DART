! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program closest_member_tool

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Program to overwrite the time on each ensemble in a restart file.

use types_mod,           only : r8
use time_manager_mod,    only : time_type, set_time_missing,               &
                                operator(==), print_time, print_date,      &
                                set_calendar_type, GREGORIAN, NO_CALENDAR, &
                                get_calendar_type

use utilities_mod,       only : register_module, do_output,                &
                                error_handler, nmlfileunit, E_MSG, E_ERR,  &
                                timestamp, find_namelist_in_file,          &
                                check_namelist_read, logfileunit,          &
                                do_nml_file, do_nml_term, open_file, close_file
                                
use       sort_mod,      only : slow_index_sort

use assim_model_mod,     only : static_init_assim_model, get_model_size,   &
                                open_restart_read, open_restart_write,     &
                                awrite_state_restart, aread_state_restart, &
                                close_restart

use mpi_utilities_mod,    only : initialize_mpi_utilities, task_count,     &
                                 finalize_mpi_utilities


implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer               :: iunit, model_size, io, ens
integer, allocatable  :: index_list(:)
character(len = 128)  :: ifile, msgstring
real(r8), allocatable :: mean(:), member(:), diffs(:)
type(time_type)       :: mean_time, member_time, mean_advance_time, advance_time

character(len=64)     :: method_name(4) = (/     &
   "Simple Difference    ", &
   "Normalized Difference", &
   "Simple RMSE          ", &
   "Normalized RMSE      "  /)

!----------------------------------------------------------------
! These variables are namelist-controllable.
!
character(len = 128) :: input_file_name        = "filter_restart"
character(len = 128) :: output_file_name       = "closest_restart"
integer              :: ens_size               = 1
logical              :: single_restart_file_in = .true.
integer              :: difference_method      = 4

!----------------------------------------------------------------
! different methods to compute 'distance' from mean:
!  1 = simple absolute difference
!  2 = normalized absolute difference
!  3 = simple rmse difference
!  4 = normalized rmse difference
!
!  (suggest more...)
!----------------------------------------------------------------

namelist /closest_member_tool_nml/  &
   input_file_name,              &
   output_file_name,             &
   ens_size,                     &
   single_restart_file_in,       &
   difference_method


logical :: input_is_model_advance_file  = .false.
! FIXME: could add this to namelist:
!      input_is_model_advance_file 
! right now we don't output model_advance means, so for now
! it stays out.  but if you did have the mean, you could
! use this tool on it.


!----------------------------------------------------------------
!----------------------------------------------------------------

! This program should only be run with a single process
call initialize_mpi_utilities('closest_member_tool')

if(task_count() > 1) &
   call error_handler(E_ERR,'closest_member_tool','Only use single process', &
                      source,revision,revdate)

call register_module(source,revision,revdate)

! Read the namelist entry and print it
call find_namelist_in_file("input.nml", "closest_member_tool_nml", iunit)
read(iunit, nml = closest_member_tool_nml, iostat = io)
call check_namelist_read(iunit, io, "closest_member_tool_nml")

if (do_nml_file()) write(nmlfileunit, nml=closest_member_tool_nml)
if (do_nml_term()) write(     *     , nml=closest_member_tool_nml)


! Initialize the model so we can get the size.
call static_init_assim_model()
model_size = get_model_size()

write(msgstring, *) 'Model size/restart data length =', model_size
call error_handler(E_MSG,'',msgstring)
write(msgstring, *) 'Ensemble member count = ', ens_size
call error_handler(E_MSG,'',msgstring)
write(msgstring, *) 'Computing difference using method: '//trim(method_name(difference_method))
call error_handler(E_MSG,'',msgstring)


! make space for the mean and a single member, plus place to sort list for output
allocate(mean(model_size), member(model_size), index_list(ens_size), diffs(ens_size))

mean_time         = set_time_missing()
member_time       = set_time_missing()
mean_advance_time = set_time_missing()
mean_advance_time = set_time_missing()

! read in the mean - always in a separate file
iunit = open_restart_read(trim(input_file_name)//'.mean')
! Read in the advance time if present
if (input_is_model_advance_file) then
   call aread_state_restart(mean_time, mean, iunit, mean_advance_time)
else
   call aread_state_restart(mean_time, mean, iunit)
endif
call close_restart(iunit)


! either loop over individual files or open a single file and read a member
! at a time.  same functionality; where the file open/close happens differs.
if (single_restart_file_in) then

   ! One restart file - open once and loop on the read

   iunit = open_restart_read(trim(input_file_name))

   ! Read in the advance time if present
   do ens=1, ens_size
      if (input_is_model_advance_file) then
         call aread_state_restart(member_time, member, iunit, advance_time)
      else
         call aread_state_restart(member_time, member, iunit)
      endif

      ! FIXME: should put in tests for time matching mean to avoid comparing
      ! against the wrong file.

      !------------------- Compute scaled RMSE   -----------------------

      diffs(ens) = compute_diff(mean, member, model_size) 
     
      !------------------- Compute scaled RMSE   -----------------------

   enddo
   call close_restart(iunit)

else
   do ens=1, ens_size
 
      ! add member number as a suffix: e.g. base.0000
      write(ifile, "(a,a,i4.4)") trim(input_file_name), '.', ens

      !------------------- Read restart from file ----------------------
      iunit = open_restart_read(ifile)
      ! Read in the advance time if present
      if (input_is_model_advance_file) then
         call aread_state_restart(member_time, member, iunit, advance_time)
      else
         call aread_state_restart(member_time, member, iunit)
      endif
      call close_restart(iunit)
      !------------------- Read restart from file ----------------------
      
      !------------------- Compute scaled RMSE   -----------------------

      diffs(ens) = compute_diff(mean, member, model_size) 
     
      !------------------- Compute scaled RMSE   -----------------------

   enddo

endif

!------------------- Print out results     -----------------------

call slow_index_sort(diffs, index_list, ens_size)
call error_handler(E_MSG, '', ' ')
write(msgstring, "(A,I5)") 'Member with the minimum difference from the mean is ', index_list(1)
call error_handler(E_MSG, '', msgstring)
call error_handler(E_MSG, '', ' ')

do ens=1, ens_size
   write(msgstring, "(A,I5,A,F18.6)") "Member ", index_list(ens), " difference ", diffs(index_list(ens))
   call error_handler(E_MSG, '', msgstring)
enddo

!------------------- Print out results     -----------------------

!------------------- Write results to file -----------------------

! if the input is a single file, write the ensemble member number to a file.
! if the input is separate files, write the full filename to a file.

iunit = open_file(output_file_name, 'formatted', 'write')

if (single_restart_file_in) then
   write(iunit, "(I4)") index_list(1)
else
   write(iunit, "(A,A,I4.4)") trim(input_file_name), '.', index_list(1)
endif

call close_file(iunit)
  
call error_handler(E_MSG, '', ' ')
write(msgstring, *) 'Writing closest member information to file: ', trim(output_file_name)
call error_handler(E_MSG, '', msgstring)

!------------------- Write results to file -----------------------

deallocate(mean, member, index_list, diffs)

call finalize_mpi_utilities()   ! now closes log file, too

!----------------------------------------------------------------
!----------------------------------------------------------------

contains

function compute_diff(target, candidate, arraysize)
 real(r8), intent(in) :: target(:)
 real(r8), intent(in) :: candidate(:)
 integer,  intent(in) :: arraysize

 real(r8) :: compute_diff

integer  :: i
real(r8) :: val, r, diff, biggest

select case (difference_method)

   ! simple absolute difference
   case (1)

      val = 0.0
      do i = 1, arraysize  
         val = val + abs(target(i) - candidate(i))
      enddo 

      compute_diff = val

   ! normalized absolute difference
   case (2)

      val = 0.0
      do i = 1, arraysize  
         if (target(i) == 0.0) then
            if (candidate(i) == 0.0) then
               diff = 0.0
            else
               diff = candidate(i)
            endif
         else
            diff = abs((target(i) - candidate(i)) / target(i))
         endif
         val = val + diff
      enddo
      
      compute_diff = val

   ! simple rmse difference
   case (3)

      val = 0.0
      do i = 1, arraysize  
         diff = target(i) - candidate(i)
         val = val + ( diff * diff )
      enddo
      
      if (val > 0) then
         compute_diff = sqrt(val)
      else
         compute_diff = val
      endif

   ! normalized rmse difference
   case (4)

      val = 0.0
      do i = 1, arraysize  
         if (target(i) == 0.0) then
            if (candidate(i) == 0.0) then
               diff = 0.0
            else
               diff = candidate(i)
            endif
         else
            diff = (target(i) - candidate(i)) / target(i)
         endif
         val = val + (diff*diff)
      enddo
      
      if (val > 0) then
         compute_diff = sqrt(val)
      else
         compute_diff = val
      endif

   case default
      write(msgstring, *) 'Valid values for difference_method are 1-4, value is', difference_method
      call error_handler(E_ERR,'closest_member_tool','Bad value for difference_method', &
                         source,revision,revdate, text2=msgstring)
end select


end function compute_diff

end program closest_member_tool
