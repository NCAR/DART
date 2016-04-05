! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This file contains a module and a program. The module exists only to enclose
!> a comparison routine that's needed by the sort module.  The actual program
!> follows.

!> Utility to sort observations at the same time into a consistent order.  
!> Obs sequence files are guarenteed to be traversed in time order, so 
!> running them through the standard obs_sequence_tool will physically 
!> order them in the file in the same way a linked-list traversal would.  
!> But for multiple obs at the same time there is no way to indicate how to order
!> them relative to each other.  This tool tries to sort same-time
!> obs based on location, then kind, then variance.  If there are
!> duplicate obs in the same file this might be a way to get them
!> together and possibly add code to remove them.  Users have in the past
!> also requested the option to order different kinds of obs at the same
!> time in a particular order (e.g. all the radar obs before the others)
!> and this might allow us to support that.

module special_sort

use        types_mod, only : r8
use     location_mod, only : location_type, get_location, LocationDims,        &
                             write_location 
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind,     &
                             get_obs_def_location, get_obs_def_error_variance
use time_manager_mod, only : time_type, operator(>), operator(==),             &
                             operator(<), operator(/=), print_time
use obs_sequence_mod, only : obs_type, get_obs_def

implicit none
private

! obssort function needed by the sort module.  will be called to do
! a comparison of 2 observations.  returns -1 for a<b, 0 for a==b, 1 for a>b
! obs_this_bin needed for accessing from the program below
public :: obssort, obs_this_bin

type(obs_type), allocatable :: obs_this_bin(:)

contains


!-----------------------------------------------------------------------
!> The sort module code will call this function with i, j, which is a
!> request to compare obs_this_bin(i) and obs_this_bin(j).
!> It will sort first by time, then by location, type, value, etc.
!>
!> Must return -1 if i < j ; 0 if i == j ; 1 if i > j
!>
!> If you want to change what values are used to order the obs do it here.
!> This order was picked to try to minimize the number of comparisons
!> needed, and to keep obs at the same location together.  In filter
!> at assimilation time there is an advantage to having obs at the same 
!> location be together in the input file - distances can be cached and 
!> reused for better performance.
!>
!> @param i index of first obs to compare
!> @param j index of other obs to compare
!>

function obssort(i, j)

integer, intent(in) :: i
integer, intent(in) :: j
integer :: obssort

type(obs_def_type)  :: this_obs_def1, this_obs_def2
integer             :: this_kind1, this_kind2
type(location_type) :: this_loc1, this_loc2
type(time_type)     :: this_time1, this_time2
real(r8)            :: this_var1, this_var2
real(r8)            :: loc1(LocationDims), loc2(LocationDims)   ! try for general
integer             :: ndim
character(len=129)  :: locstring1, locstring2

logical             :: local_debug = .false.

call get_obs_def(obs_this_bin(i), this_obs_def1)
call get_obs_def(obs_this_bin(j), this_obs_def2)


this_time1 = get_obs_def_time(this_obs_def1)
this_time2 = get_obs_def_time(this_obs_def2)

this_loc1 = get_obs_def_location(this_obs_def1)
this_loc2 = get_obs_def_location(this_obs_def2)

loc1 = get_location(this_loc1)
loc2 = get_location(this_loc2)

this_kind1 = get_obs_kind(this_obs_def1)
this_kind2 = get_obs_kind(this_obs_def2)

this_var1 = get_obs_def_error_variance(this_obs_def1)
this_var2 = get_obs_def_error_variance(this_obs_def2)

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

! sort first by time
if (this_time1 > this_time2) then
   obssort = 1
   return
else if (this_time1 < this_time2) then
   obssort = -1
   return
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
! (they could have different metadata and we
! have no way to know about that.)
obssort = 0

if (local_debug) print *, 'decided items ', i, j, ' are same'

end function obssort

end module special_sort

!-----------------------------------------------------------------------
!> Open an obs sequence file and reorder observations that have the same
!> timestamp so they have a reproducible order.  Could be used to make
!> two obs_sequence files which have had different processing have their
!> observations in the same order, or look for duplicate observations.
!>

program obs_sort

use         types_mod, only : r8, metadatalength
use     utilities_mod, only : register_module, initialize_utilities,            &
                              find_namelist_in_file, check_namelist_read,       &
                              error_handler, E_ERR, E_MSG, nmlfileunit,         &
                              do_nml_file, do_nml_term, finalize_utilities
use          sort_mod, only : index_sort
use       obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind
use      obs_kind_mod, only : max_obs_kinds, get_obs_kind_name
use  time_manager_mod, only : time_type, operator(>), print_time, set_time,     &
                              print_date, set_calendar_type, operator(==),      &
                              operator(/=), get_calendar_type, NO_CALENDAR,     &
                              operator(-)
use  obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                              init_obs, assignment(=),                          &
                              init_obs_sequence, static_init_obs_sequence,      &
                              read_obs_seq_header, read_obs_seq, get_num_obs,   &
                              get_first_obs, get_next_obs, get_obs_def,         &
                              insert_obs_in_seq, get_num_copies, get_num_qc,    &
                              get_copy_meta_data, get_qc_meta_data,             &
                              set_copy_meta_data, set_qc_meta_data,             &
                              destroy_obs, destroy_obs_sequence,                &
                              get_num_key_range, get_obs_key
use obs_utilities_mod, only : print_obs_seq, validate_obs_seq_time, print_metadata

use special_sort

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs_in, next_obs_in
type(obs_type)          :: obs_out, prev_obs_out
type(time_type)         :: this_time, prev_time
logical                 :: is_this_last
integer                 :: size_seq_in, size_seq_out
integer                 :: num_copies_in, num_qc_in
integer                 :: num_inserted, iunit, io, i, j
integer                 :: max_num_obs, file_id, sort_count
integer                 :: num_rejected_badqc, num_rejected_diffqc
integer                 :: num_rejected_other
integer, allocatable    :: indices(:)
character(len=129)      :: read_format
logical                 :: pre_I_format, cal
character(len=512)      :: msgstring, msgstring1, msgstring2
type(obs_def_type)      :: this_obs_def

character(len=metadatalength) :: meta_data

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 5000

! lazy, pick big number.  make it bigger if too small.
integer, parameter :: max_obs_input_types = 500

!-----------------------------------------------------------------------
! Namelist input with default values

character(len=160) :: filename_in  = ''
character(len=160) :: filename_out = ''
logical            :: print_only   = .false.
character(len=32)  :: calendar     = 'Gregorian'
logical            :: debug        = .false.  ! true for more output

namelist /obs_sort_nml/ &
         filename_in, filename_out, &
         print_only, calendar, debug

!-----------------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!-----------------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_sort_nml", iunit)
read(iunit, nml = obs_sort_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_sort_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_sort_nml)
if (do_nml_term()) write(     *     , nml=obs_sort_nml)

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
   call error_handler(E_ERR,'obs_sort',msgstring)
endif

write(msgstring, *) 'Starting to process input sequence file: '
write(msgstring1,*)  trim(filename_in)
call error_handler(E_MSG,'obs_sort',msgstring, &
                   text2=msgstring1)

call read_obs_seq(filename_in, 0, 0, 0, seq_in)

! sanity check - ensure the linked list times are in increasing time order
! Some observation sequences that have not been created with the DART tools
! have been known to violate this assumption. 
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
allocate(obs_this_bin(max_num_obs), indices(max_num_obs))
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

! is this needed?
if (print_only) call print_obs_seq(seq_in, filename_in)

!-----------------------------------------------------------------------
! Start to insert obs from sequence_in into sequence_out
!
! NOTE: insert_obs_in_seq CHANGES the obs passed in.
!       Must pass a copy of incoming obs to insert_obs_in_seq.
!-----------------------------------------------------------------------

num_inserted = 0
num_rejected_badqc = 0
num_rejected_diffqc = 0
num_rejected_other = 0

! initialize "obs_in" here by calling get_first_obs() so we 
! can get the first obs_def to prime the loop below

if (.not. get_first_obs(seq_in, obs_in) )  then
   write(msgstring, *)'error getting first observation in ',trim(filename_in)
   call error_handler(E_MSG,'obs_sort', msgstring)
endif

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
      call index_sort(indices, sort_count, obssort)
      if (debug) print *, 'sorted indices:'
      if (debug) print *, indices(1:sort_count)

      if (num_inserted > 0) then
         call insert_obs_in_seq(seq_out, obs_this_bin(indices(1)), prev_obs_out)
      else
         call insert_obs_in_seq(seq_out, obs_this_bin(indices(1)))
      endif

      prev_obs_out = obs_this_bin(indices(1))
      do i=2, sort_count
         call insert_obs_in_seq(seq_out, obs_this_bin(indices(i)), prev_obs_out)
         prev_obs_out = obs_this_bin(indices(i))
      enddo

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

if (.not. print_only) then
   print*, '---------  Obs seqs '
   print*, 'Number of obs input sequence    :         ', size_seq_in
   print*, 'Number of obs copied to output  :         ', num_inserted
   print*, '---------------------------------------------------------'
endif


write(msgstring, *) 'Starting to process output sequence file ', &
                            trim(filename_out)
call error_handler(E_MSG,'obs_sort',msgstring)

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
! do not destroy prev_obs_out; deleting obs_out destroys it as well


deallocate(indices, obs_this_bin)

call shutdown()

!-----------------------------------------------------------------------
! end of main program.
!-----------------------------------------------------------------------


contains


!-----------------------------------------------------------------------
!> Initialize modules used

subroutine setup()

call initialize_utilities('obs_sort')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

end subroutine setup


!-----------------------------------------------------------------------
!> Finalize modules used
subroutine shutdown()

call finalize_utilities('obs_sort')

end subroutine shutdown



end program obs_sort


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
