! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program closest_member_tool

! Program to overwrite the time on each ensemble in a restart file.

use types_mod,         only : r8
use time_manager_mod,  only : time_type, set_time_missing,               &
                              operator(/=), print_time
 
use utilities_mod,     only : register_module, find_namelist_in_file,        &
                              error_handler, nmlfileunit, E_MSG, E_ERR,      &
                              check_namelist_read, do_nml_file, do_nml_term, &
                              open_file, close_file

use  location_mod,     only : location_type

use  obs_kind_mod,     only : get_num_raw_obs_kinds, get_raw_obs_kind_index, &
                              paramname_length, get_raw_obs_kind_name

use  sort_mod,         only : slow_index_sort

use assim_model_mod,   only : static_init_assim_model, get_model_size,   &
                              open_restart_read, aread_state_restart,    &
                              close_restart, get_state_meta_data

use mpi_utilities_mod, only : initialize_mpi_utilities, task_count,     &
                              finalize_mpi_utilities


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer               :: iunit, model_size, io, ens, i, j, kindindex
integer, allocatable  :: index_list(:)
character(len = 128)  :: ifile, msgstring, msgstring1
real(r8), allocatable :: mean(:), member(:), diffs(:)
logical               :: allkinds, done
logical, allocatable  :: usekind(:), useindex(:)
type(location_type)   :: loc
integer               :: num_kinds, stype
type(time_type)       :: mean_time, member_time, mean_advance_time, advance_time
integer, parameter    :: max_list_len = 500

character(len=64)     :: method_name(4) = (/     &
   "Simple Difference    ", &
   "Normalized Difference", &
   "Simple RMS Diff      ", &
   "Normalized RMS Diff  "  /)

!----------------------------------------------------------------
! These variables are namelist-controllable.
!
character(len = 128) :: input_file_name        = "filter_restart"
character(len = 128) :: output_file_name       = "closest_restart"
integer              :: ens_size               = 1
logical              :: single_restart_file_in = .true.
integer              :: difference_method      = 4
character(len = paramname_length) :: use_only_kinds(max_list_len) = ''

!----------------------------------------------------------------
! different methods to compute 'distance' from mean:
!  1 = simple absolute difference
!  2 = normalized absolute difference
!  3 = simple rms difference
!  4 = normalized rms difference
!
!  (suggest more...)
!----------------------------------------------------------------

namelist /closest_member_tool_nml/  &
   input_file_name,              &
   output_file_name,             &
   ens_size,                     &
   single_restart_file_in,       &
   difference_method,            &
   use_only_kinds         


logical :: input_is_model_advance_file  = .false.
! FIXME: could add this to namelist:
!      input_is_model_advance_file 
! right now we don't output model_advance means, so for now
! it stays out.  but if you did have the mean, you could
! use this tool on it.


!----------------------------------------------------------------
! program start
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
allocate(mean(model_size), member(model_size))
allocate(index_list(ens_size), diffs(ens_size))

! are we adding up only differences from particular kinds, or the entire vector?
if (use_only_kinds(1) /= '') then
   allkinds = .false.

   num_kinds = get_num_raw_obs_kinds()
   allocate(usekind(num_kinds))
   usekind = .false.

   done = .false.
   KindList:do i=1, max_list_len
      if (use_only_kinds(i) == '') then
         done = .true.
         exit KindList
      endif
      kindindex = get_raw_obs_kind_index(use_only_kinds(i))   
      if (kindindex < 0) then
         write(msgstring, *) 'unrecognized KIND string: '//trim(use_only_kinds(i))
         call error_handler(E_ERR,'closest_member_tool', msgstring, &
                            source,revision,revdate)
      endif
      usekind(kindindex) = .true.

   enddo KindList

   if (.not. done) then
      write(msgstring, *) 'cannot have more than ', max_list_len, ' kinds'
      call error_handler(E_ERR,'closest_member_tool', msgstring, &
                         source,revision,revdate)
   endif

   write(msgstring, *) 'Computing difference based only on items in state vector items of kind:'
   call error_handler(E_MSG,'',msgstring)
   do i=1, num_kinds
      if (usekind(i)) then
         write(msgstring, *) '   ', trim(get_raw_obs_kind_name(i))
         call error_handler(E_MSG,'',msgstring)
      endif
   enddo

else
   allkinds = .true.
endif

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

! if we are not processing all kinds of state vector items, set up a mask
! for the kinds we are.  do this once at the start so we don't replicate
! work.  the usekind(number_of_different_kinds) array says whether we are
! going to add differences for this type.   the useindex(state_vector_len)
! array is being precomputed here so it's fast to loop over the ensemble
! members and only add up differences for the kinds of interest.
allocate(useindex(model_size))
if (.not. allkinds) then
   useindex(:) = .false.

   j = 0
   do i=1, model_size
      call get_state_meta_data(i, loc, stype)
      if (stype < 1 .or. stype > num_kinds) then
         write(msgstring, *) 'bad KIND from get_state_meta_data, ', stype, ' for index ', i 
         write(msgstring1, *) 'must be between 1 and ', num_kinds
         call error_handler(E_ERR,'closest_member_tool', msgstring, &
                            source,revision,revdate, text2=msgstring1)

      endif
      if (usekind(stype)) then 
         useindex(i) = .true.
         j = j + 1
      endif
   enddo

   write(msgstring, *) 'using ', j, ' of ', model_size, ' items in the state vector'
   call error_handler(E_MSG,'closest_member_tool', msgstring)
else
   ! use everything.
   useindex(:) = .true.
endif

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

      if (mean_time /= member_time) then
         call print_time(mean_time, "time of ensemble mean data")
         call print_time(member_time, "time of ensemble member data")
         write(msgstring, *) 'member ', ens, ' has a different timestamp than mean'
         call error_handler(E_ERR,'closest_member_tool', msgstring)
      endif

      !------------------- Compute difference    -----------------------

      diffs(ens) = compute_diff(mean, member, model_size) 
     
      !------------------- Compute difference    -----------------------

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
      
      if (mean_time /= member_time) then
         call print_time(mean_time, "time of ensemble mean data")
         call print_time(member_time, "time of ensemble member data")
         write(msgstring, *) 'member ', ens, ' has a different timestamp than mean'
         call error_handler(E_ERR,'closest_member_tool', msgstring)
      endif

      !------------------- Compute difference    -----------------------

      diffs(ens) = compute_diff(mean, member, model_size) 
     
      !------------------- Compute difference    -----------------------

   enddo

endif

!------------------- Print out results     -----------------------

call slow_index_sort(diffs, index_list, ens_size)
call error_handler(E_MSG, '', ' ')
write(msgstring, "(A,I5)") 'Member with the minimum difference from the mean is ', index_list(1)
call error_handler(E_MSG, '', msgstring)
call error_handler(E_MSG, '', ' ')

do ens=1, ens_size
   write(msgstring, "(A,I5,A,G18.6)") "Member ", index_list(ens), " difference ", diffs(index_list(ens))
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

deallocate(mean, member, index_list, diffs, useindex)
if (.not. allkinds) deallocate(usekind)

call finalize_mpi_utilities()   ! now closes log file, too

!----------------------------------------------------------------
!----------------------------------------------------------------

contains

function compute_diff(target, candidate, arraysize)
 real(r8), intent(in) :: target(:)
 real(r8), intent(in) :: candidate(:)
 integer,  intent(in) :: arraysize

 real(r8) :: compute_diff

real(r8), allocatable :: adiff(:)

! new strategy:  compute an array of differences and sum them at the end.
! try to use array operations when possible.  useindex() is a logical array
! that can be used as a mask if only some kinds are going to be used.
! it is set to all .true. if no kinds were set, so it can be used in 
! either case.

allocate(adiff(arraysize)) 
adiff = 0.0_r8

select case (difference_method)

   ! simple absolute difference
   case (1)
      where (useindex) adiff = abs(target - candidate)
      
      compute_diff = sum(adiff)

   ! normalized absolute difference
   case (2)

      where (useindex) 
         where (target /= 0.0_r8) 
            adiff = abs((target - candidate) / target)
         elsewhere
            adiff = abs(candidate)
         endwhere
      endwhere

      compute_diff = sum(adiff)

   ! simple rms difference
   case (3)

      where (useindex) adiff = target - candidate

      compute_diff = sum(adiff * adiff)
      if (compute_diff > 0.0_r8) compute_diff = sqrt(compute_diff)

   ! normalized rms difference
   case (4)

      where (useindex) 
         where (target /= 0.0_r8) 
            adiff = (target - candidate) / target
         elsewhere
            adiff = candidate
         endwhere
      endwhere

      compute_diff = sum(adiff * adiff)
      if (compute_diff > 0.0_r8) compute_diff = sqrt(compute_diff)

   case default
      write(msgstring, *) 'Valid values for difference_method are 1-4, value is', difference_method
      call error_handler(E_ERR,'closest_member_tool','Bad value for difference_method', &
                         source,revision,revdate, text2=msgstring)
end select


deallocate(adiff)

end function compute_diff


end program closest_member_tool

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
