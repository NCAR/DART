! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


!> This program removes duplicate observations.

!> This file contains 1 module and 1 program.  the module code has
!> to come first, so page down for the main program. 
!> The module is a custom sort routine needed to compare 2 observations, 
!> and then the program uses the sort to match duplicate observations 
!> in an obs_seq file.


!> Module:
!> Routine to sort observations at the same time into a consistent order.  
!> Obs sequence files are guarenteed to be traversed in time order; 
!> running them through the standard obs_sequence_tool will physically 
!> order them in the file in the same way a linked-list traversal would.  
!>
!> But for multiple obs at the same time there is no way to indicate or
!> control how to order them relative to each other.  This tool sorts 
!> same-time obs based on location, then kind, then variance.  If there 
!> are duplicate obs in the same file this helps get them together and
!> pass through only one copy of the observation.

module special_sort

use        types_mod, only : r8
use         sort_mod, only : index_sort
use     location_mod, only : location_type, get_location, set_location,        &
                             LocationName, read_location, operator(/=),        &
                             write_location, LocationDims
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs,     &
                             get_obs_def_location, get_obs_def_error_variance, operator(==)
use time_manager_mod, only : time_type, operator(/=), print_time
use obs_sequence_mod, only : get_obs_def, obs_type

implicit none
private

public :: obssort, obs_this_bin

type(obs_type), allocatable :: obs_this_bin(:)

contains


!---------------------------------------------------------------------

function obssort(i, j)
 integer, intent(in) :: i, j
 integer :: obssort

! this is requesting a compare of obs_this_bin(i) and obs_this_bin(j)
! they should have identical times, so the compare needs to be by 
! location, type, value, etc.   return -1 if i < j ; 0 if == ; 1 if i > j

type(obs_def_type)  :: this_obs_def1, this_obs_def2
integer             :: this_kind1, this_kind2
type(location_type) :: this_loc1, this_loc2
type(time_type)     :: this_time1, this_time2
real(r8)            :: this_var1, this_var2
real(r8)            :: loc1(LocationDims), loc2(LocationDims)   ! try for general?
integer             :: ndim
character(len=129)  :: locstring1, locstring2
logical             :: local_debug

local_debug = .false.

call get_obs_def(obs_this_bin(i), this_obs_def1)
call get_obs_def(obs_this_bin(j), this_obs_def2)

this_time1 = get_obs_def_time(this_obs_def1)
this_time2 = get_obs_def_time(this_obs_def2)

this_loc1 = get_obs_def_location(this_obs_def1)
this_loc2 = get_obs_def_location(this_obs_def2)

loc1 = get_location(this_loc1)
loc2 = get_location(this_loc2)

this_kind1 = get_obs_def_type_of_obs(this_obs_def1)
this_kind2 = get_obs_def_type_of_obs(this_obs_def2)

this_var1 = get_obs_def_error_variance(this_obs_def1)
this_var2 = get_obs_def_error_variance(this_obs_def2)

if (this_time1 /= this_time2) then
  print *, 'error, times not the same'
  print *, 'comparing items ', i, j
  call print_time(this_time1, 'time1')
  call print_time(this_time2, 'time2')
  stop
endif

if (local_debug) then
  print *, 'comparing items ', i, j
  call print_time(this_time1, 'time: ')
  print *, 'kinds ', this_kind1, this_kind2
  print *, 'vars ', this_var1, this_var2
  call write_location(0, this_loc1, charstring=locstring1)
  call write_location(0, this_loc2, charstring=locstring2)
  print *, 'locs: '
  print *, trim(locstring1)
  print *, trim(locstring2)
  print *, ''
endif

! try for a general location solution
do ndim=1, LocationDims
 
 if (loc1(ndim) > loc2(ndim)) then
    obssort = 1
    return
 else if (loc1(ndim) < loc2(ndim)) then
    obssort = -1
    return
 endif

enddo
 
! locations the same, so try types
if (this_kind1 > this_kind2) then
   obssort = 1
   return
else if (this_kind1 < this_kind2) then
   obssort = -1
   return
endif

! same up to now, so try errors (variance)
if (this_var1 > this_var2) then
   obssort = 1
   return
else if (this_var1 < this_var2) then
   obssort = -1
   return
endif

! ok, i give up.  they're the same.
! or enough for us to not try to resort them.
obssort = 0

if (local_debug) then
  print *, 'decided items ', i, j, ' are same'
  call print_time(this_time1, 'time: ')
  print *, 'kinds ', this_kind1, this_kind2
  print *, 'vars ', this_var1, this_var2
  call write_location(0, this_loc1, charstring=locstring1)
  call write_location(0, this_loc2, charstring=locstring2)
  print *, 'locs: '
  print *, trim(locstring1)
  print *, trim(locstring2)
  print *, ''
endif

end function obssort

end module special_sort

!---------------------------------------------------------------------

!> simple program that opens an obs_seq file and loops over the obs

!> program that opens an obs_seq file and loops over the obs
!> and copies them to a new output file.   this is intended to be a
!> template for programs that want to alter existing obs in some simple way.

program obs_remove_dups

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             open_file, close_file, finalize_utilities
use         sort_mod, only : index_sort
use     location_mod, only : location_type, get_location, set_location,        &
                             LocationName, read_location, operator(/=),        &
                             write_location
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs, &
                             operator(==) !, print_obs_def
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs
use time_manager_mod, only : time_type, operator(>), print_time, set_time,     &
                             print_date, set_calendar_type, operator(==),      &
                             operator(/=), get_calendar_type, NO_CALENDAR,     &
                             operator(-), set_time_missing
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=), get_obs_def,             &
                             init_obs_sequence, static_init_obs_sequence,      &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_last_obs, get_next_obs,        &
                             insert_obs_in_seq, get_num_copies, get_num_qc,    &
                             get_copy_meta_data, get_qc_meta_data,             &
                             set_copy_meta_data, set_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence,                &
                             get_num_key_range, get_obs_key, get_qc,           &
                             copy_partial_obs, get_next_obs_from_key,          &
                             get_obs_def, set_obs_def, operator(==), operator(/=)

use special_sort

implicit none

character(len=*), parameter :: source = 'obs_remove_dups.f90'

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs_in, next_obs_in, last_obs
type(obs_type)          :: obs_out, prev_obs_out
type(time_type)         :: this_time, prev_time
logical                 :: is_this_last
integer                 :: size_seq_in, size_seq_out
integer                 :: num_copies_in, num_qc_in
integer                 :: num_inserted, iunit, io, i, j
integer                 :: max_num_obs, file_id, sort_count, last_in
integer, allocatable    :: index(:)
character(len = 129)    :: read_format
logical                 :: pre_I_format, cal
character(len = 256)    :: msgstring, msgstring1, msgstring2
type(obs_def_type)      :: this_obs_def, that_obs_def

character(len = metadatalength) :: meta_data

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 5000

! lazy, pick big number.  make it bigger if too small.
integer, parameter :: max_obs_input_types = 500

!----------------------------------------------------------------
! Namelist input with default values


character(len = 160) :: filename_in = ''
character(len = 160) :: filename_out = ''

! if true, only compare the obs_defs and not the obs values
logical              :: ignore_values    = .false.

logical              :: print_only    = .false.
character(len=32)    :: calendar      = 'Gregorian'

! true for more output
logical              :: debug = .false.

namelist /obs_remove_dups_nml/ &
         filename_in, filename_out, &
         ignore_values, print_only, calendar, debug

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_remove_dups_nml", iunit)
read(iunit, nml = obs_remove_dups_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_remove_dups_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_remove_dups_nml)
if (do_nml_term()) write(     *     , nml=obs_remove_dups_nml)

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! set a logial to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! if you add anything to the namelist, you can process it here.

! end of namelist processing and setup


! single pass algorithm (unlike other obs tools).

call read_obs_seq_header(filename_in, num_copies_in, num_qc_in, &
   size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
   close_the_file = .true.)

if (max_num_obs == 0) then
   write(msgstring,*) 'No obs in input sequence file ', trim(filename_in)
   call error_handler(E_ERR,'obs_remove_dups',msgstring)
endif

write(msgstring, *) 'Starting to process input sequence file: '
write(msgstring1,*)  trim(filename_in)
call error_handler(E_MSG,'obs_remove_dups',msgstring, &
                   text2=msgstring1)

call read_obs_seq(filename_in, 0, 0, 0, seq_in)

! sanity check - ensure the linked list times are in increasing time order
call validate_obs_seq_time(seq_in, filename_in)

! output is same size (or less) than input, generally.
! if this program is going to dup obs, account for it here.
size_seq_out = max_num_obs

! blank line, start of actually creating output file
call error_handler(E_MSG,' ',' ')

! Initialize individual observation variables
call init_obs(     obs_in,  num_copies_in, num_qc_in)
call init_obs(next_obs_in,  num_copies_in, num_qc_in)
call init_obs(     obs_out, num_copies_in, num_qc_in)
call init_obs(prev_obs_out, num_copies_in, num_qc_in)
   
! space for sorting obs with the same timestamp
allocate(obs_this_bin(max_num_obs), index(max_num_obs))
do i=1, max_num_obs
   call init_obs(obs_this_bin(i), num_copies_in, num_qc_in)
enddo

! create the output sequence here 
call init_obs_sequence(seq_out, num_copies_in, num_qc_in, size_seq_out)
do j=1, num_copies_in
   meta_data = get_copy_meta_data(seq_in, j)
   call set_copy_meta_data(seq_out, j, meta_data)
enddo
do j=1, num_qc_in
   meta_data = get_qc_meta_data(seq_in, j)
   call set_qc_meta_data(seq_out, j, meta_data)
enddo

call print_obs_seq(seq_in, filename_in)

!! is this needed?
!if (print_only) call print_obs_seq(seq_in, filename_in)

!-------------------------------------------------------------
! Start to insert obs from sequence_in into sequence_out
!
! NOTE: insert_obs_in_seq CHANGES the obs passed in.
!       Must pass a copy of incoming obs to insert_obs_in_seq.
!--------------------------------------------------------------
num_inserted = 0

if ( get_first_obs(seq_in, obs_in) )  then

   is_this_last = .false.
   next_obs_in = obs_in
   call get_obs_def(obs_in, this_obs_def)
   prev_time = get_obs_def_time(this_obs_def)

   ObsLoop : do while ( .not. is_this_last )

      obs_in = next_obs_in

      ! obs_out will be modified when it is inserted in the output sequence
      ! so we have to make a copy of obs_in before modifiying it.
      obs_out = obs_in

      ! see if this obs is the same time as the prev obs
      ! if not, carry on by putting it into the output.
      ! if it's the same time, we have to sort first.

      call get_obs_def(obs_out, this_obs_def)
      this_time = get_obs_def_time(this_obs_def)
      if (debug) print *, 'next observation: '
      if (debug) call print_time(this_time, 'obs_in this_time')
      if (debug) call print_time(prev_time, 'obs_in prev_time')

      if (prev_time == this_time) then

         if (debug) print *, 'matched prev_time'
         sort_count = 0

         SortObsLoop : do while ( .not. is_this_last )

            obs_in = next_obs_in

            sort_count = sort_count + 1
            obs_this_bin(sort_count) = obs_in

            call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

            call get_obs_def(next_obs_in, this_obs_def)
            this_time = get_obs_def_time(this_obs_def)
            if (debug) call print_time(this_time, 'next_obs_in')
            if (debug) print *, 'sort_count = ', sort_count

            if (prev_time /= this_time) exit SortObsLoop

         enddo SortObsLoop

         if (debug) print *, 'out of loop, sort_count = ', sort_count
         ! sort obs here
         call index_sort(index, sort_count, obssort)
         if (debug) print *, 'sorted index:'
         if (debug) print *, index(1:sort_count)

         ! the first obs in this time bin can't be a dup, so output it
         if (num_inserted > 0) then
            call insert_obs_in_seq(seq_out, obs_this_bin(index(1)), prev_obs_out)
         else
            call insert_obs_in_seq(seq_out, obs_this_bin(index(1)))
         endif

         prev_obs_out = obs_this_bin(index(1))
         last_in = 1
         call get_obs_def(obs_this_bin(index(1)), this_obs_def)
         SAMEOBS: do i=2, sort_count
            if (debug) print *, i, 'comparing obs ', index(last_in), index(i)
            if (ignore_values) then
               call get_obs_def(obs_this_bin(index(i)), that_obs_def)
               if (debug) print *, 'obs ', index(last_in)
               !if (debug) call print_obs_def(this_obs_def)
               if (debug) print *, 'obs ', index(i)
               !if (debug) call print_obs_def(that_obs_def)
               if (this_obs_def == that_obs_def) then
                  if (debug) print *, 'same in obs_def - dup being ignored.'
                  if (debug) print *, ''
                  cycle SAMEOBS
               endif
            else
               if (debug) print *, 'obs ', index(last_in)
               !if (debug) call print_obs(obs_this_bin(index(last_in)))
               if (debug) print *, 'obs ', index(i)
               !if (debug) call print_obs(obs_this_bin(index(i)))
               if (obs_this_bin(index(last_in)) == obs_this_bin(index(i))) then
                  if (debug) print *, 'same in obs vals - dup being ignored.'
                  if (debug) print *, ''
                  cycle SAMEOBS
               endif
            endif 
            if (debug) print *, 'obs differ - next obs being added to output.'
            if (debug) print *, ''
            last_in = i
            call get_obs_def(obs_this_bin(index(i)), this_obs_def)

            call insert_obs_in_seq(seq_out, obs_this_bin(index(i)), prev_obs_out)
            prev_obs_out = obs_this_bin(index(i))
         enddo SAMEOBS
   
         num_inserted = num_inserted + sort_count

         prev_time = this_time

         if (print_every > 0) then
            if ((mod(num_inserted,print_every) == 0) .or. &
                (num_inserted > print_every)) then
               print*, 'inserted number ',num_inserted,' of ',size_seq_out
            endif
         endif
   
         ! no call to get_next_obs() because we've already done it

      else

         ! Since the stride through the observation sequence file is always
         ! guaranteed to be in temporally-ascending order, we can use the
         ! 'previous' observation as the starting point to search for the
         ! correct insertion point.  This speeds up the insert code a lot.
   
         if (num_inserted > 0) then
            call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)
         else
            call insert_obs_in_seq(seq_out, obs_out)
         endif
   
         prev_obs_out = obs_out  ! update position in seq for next insert
         num_inserted = num_inserted + 1

         prev_time = this_time

         if (print_every > 0) then
            if (mod(num_inserted,print_every) == 0) then
               print*, 'inserted number ',num_inserted,' of ',size_seq_out
            endif
         endif
   
         call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)
      endif

   enddo ObsLoop

else
   write(msgstring, *)'no first observation in ',trim(filename_in)
   call error_handler(E_MSG,'obs_remove_dups', msgstring)
endif

if (.not. print_only) then
   print*, '---------  Obs seqs '
   print*, 'Number of obs input sequence    :         ', size_seq_in
   print*, 'Number of obs copied to output  :         ', num_inserted
   print*, '---------------------------------------------------------'
endif



write(msgstring, *) 'Starting to process output sequence file ', &
                            trim(filename_out)
call error_handler(E_MSG,'obs_remove_dups',msgstring)

print*, 'Number of obs in the output seq file :', get_num_key_range(seq_out)

call print_obs_seq(seq_out, filename_out)
if (.not. print_only) then
   call write_obs_seq(seq_out, filename_out)
else
   write(msgstring,*) 'Output sequence file not created; print_only in namelist is .true.'
   call error_handler(E_MSG,'', msgstring)
endif

! clean up

call destroy_obs_sequence(seq_in)
call destroy_obs_sequence(seq_out)
call destroy_obs(     obs_in )
call destroy_obs(next_obs_in )
call destroy_obs(     obs_out)
!call destroy_obs(prev_obs_out)  ! copy of something already deleted

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_remove_dups')
call static_init_obs_sequence()

end subroutine setup


!---------------------------------------------------------------------
subroutine shutdown()

call finalize_utilities('obs_remove_dups')

end subroutine shutdown


!---------------------------------------------------------------------
subroutine print_obs_seq(seq_in, filename)

! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

type(obs_sequence_type), intent(in) :: seq_in
character(len=*), intent(in)        :: filename

type(obs_type)      :: obs, next_obs
type(obs_def_type)  :: this_obs_def
logical             :: is_there_one, is_this_last
integer             :: size_seq_in
integer             :: i
integer             :: this_obs_type
integer             :: type_count(0:max_defined_types_of_obs), identity_count


! Initialize input obs_types
type_count(:) = 0
identity_count = 0

! make sure there are obs left to process before going on.
! num_obs should be ok since we just constructed this seq so it should
! have no unlinked obs.  if it might for some reason, use this instead:
! size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_remove_dups',msgstring)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq_in), get_num_qc(seq_in))
call init_obs(next_obs, get_num_copies(seq_in), get_num_qc(seq_in))

! blank line
call error_handler(E_MSG,'',' ')

write(msgstring,*) 'Processing sequence file ', trim(filename)
call error_handler(E_MSG,'',msgstring)

call print_metadata(seq_in, filename)

!-------------------------------------------------------------
! Start to process obs from seq_in
!--------------------------------------------------------------
is_there_one = get_first_obs(seq_in, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in ',trim(filename)
   call error_handler(E_MSG,'obs_remove_dups', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
! does not work with NO_CALENDAR
if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_type = get_obs_def_type_of_obs(this_obs_def)
   if (this_obs_type < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_type) = type_count(this_obs_type) + 1
   endif
!   print *, 'obs type index = ', this_obs_type
!   if(this_obs_type > 0)print *, 'obs name = ', get_name_for_type_of_obs(this_obs_type)

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def), '  Last timestamp: ')
      if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')
   endif

enddo ObsLoop


write(msgstring, *) 'Number of obs processed  :          ', size_seq_in
call error_handler(E_MSG, '', msgstring)
write(msgstring, *) '---------------------------------------------------------'
call error_handler(E_MSG, '', msgstring)
do i = 0, max_defined_types_of_obs
   if (type_count(i) > 0) then 
      write(msgstring, '(a32,i8,a)') trim(get_name_for_type_of_obs(i)), &
                                     type_count(i), ' obs'
      call error_handler(E_MSG, '', msgstring)
   endif
enddo
if (identity_count > 0) then 
   write(msgstring, '(a32,i8,a)') 'Identity observations', &
                                  identity_count, ' obs'
   call error_handler(E_MSG, '', msgstring)
endif

! another blank line
call error_handler(E_MSG, '', ' ')

! Time to clean up

call destroy_obs(     obs)
call destroy_obs(next_obs)

end subroutine print_obs_seq


!---------------------------------------------------------------------
subroutine validate_obs_seq_time(seq, filename)

! this eventually belongs in the obs_seq_mod code, but for now
! try it out here.  we just fixed a hole in the interactive create
! routine which would silently let you create out-of-time-order
! linked lists, which gave no errors but didn't assimilate the
! right obs at the right time when running filter.   this runs
! through the times in the entire sequence, ensuring they are
! monotonically increasing in time.  this should help catch any
! bad files which were created with older versions of code.

type(obs_sequence_type), intent(in) :: seq
character(len=*),        intent(in) :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq, obs_count
integer                 :: key
type(time_type)         :: last_time, this_time


! make sure there are obs left to process before going on.
size_seq = get_num_obs(seq) 
if (size_seq == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_remove_dups:validate',msgstring)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(next_obs, get_num_copies(seq), get_num_qc(seq))

obs_count = 0

!-------------------------------------------------------------
! Start to process obs from seq
!--------------------------------------------------------------
is_there_one = get_first_obs(seq, obs)

! we already tested for 0 obs above, so there should be a first obs here.
if ( .not. is_there_one )  then
   write(msgstring,*)'no first obs in sequence ' // trim(filename)
   call error_handler(E_ERR,'obs_remove_dups:validate', msgstring, source)
   return
endif

is_this_last = .false.
last_time = set_time(0, 0)
ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_time = get_obs_def_time(this_obs_def)

   if (last_time > this_time) then
      ! bad time order of observations in linked list
      call print_time(last_time, ' previous timestamp: ')
      if (cal) call print_date(last_time, '      calendar date: ')
      call print_time(this_time, '     next timestamp: ')
      if (cal) call print_date(this_time, '      calendar date: ')

      key = get_obs_key(obs)
      write(msgstring1,*)'obs number ', key, ' has earlier time than previous obs'
      write(msgstring2,*)'observations must be in increasing time order, file ' // trim(filename)
      call error_handler(E_ERR,'obs_remove_dups:validate', msgstring2, &
                         source, text2=msgstring1)
   endif

   last_time = this_time
   obs_count = obs_count + 1

   call get_next_obs(seq, obs, next_obs, is_this_last)
   if (.not. is_this_last) obs = next_obs

enddo ObsLoop

! clean up
call destroy_obs(     obs)
call destroy_obs(next_obs)

! technically not a time validation, but easy to check.  obs_count should never
! be larger than size_seq - that's a fatal error.  obs_count < size_seq would 
! suggest there are obs in the file that aren't part of the linked list.  
! this does not necessarily indicate a fatal error but it's not a common 
! situation and might indicate someone should check on the file.
if (obs_count /= size_seq) then
   write(msgstring,*) 'input sequence ', trim(filename)
   call error_handler(E_MSG,'obs_remove_dups:validate', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'obs_remove_dups:validate', msgstring, &
                         source, text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'obs_remove_dups:validate', msgstring, &
                         source, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time


!---------------------------------------------------------------------
subroutine print_metadata(seq, fname)

!
! print out the metadata strings, trimmed
!

type(obs_sequence_type), intent(in) :: seq
character(len=*), optional :: fname

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str
character(len=255) :: msgstring3

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(    seq)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring3,*)' illegal copy or obs count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_remove_dups', msgstring3, source)
endif

MetaDataLoop : do i=1, num_copies
   str = get_copy_meta_data(seq,i)

   write(msgstring,*)'Data Metadata: ',trim(str)
   call error_handler(E_MSG, '', msgstring)

enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str = get_qc_meta_data(seq,i)

   write(msgstring,*)'  QC Metadata: ', trim(str)
   call error_handler(E_MSG, '', msgstring)

enddo QCMetaData

end subroutine print_metadata

end program obs_remove_dups

