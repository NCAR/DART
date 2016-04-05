! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program obs_info

! read the input obs_seq file from stdin; count the number and first/last
! timestamps for each obs type and print out.  does not require a namelist;
! cannot do a list of files - one at a time only.

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : register_module, initialize_utilities,            &
                             error_handler, E_ERR, E_MSG, finalize_utilities
use     location_mod, only : location_type, get_location, set_location,        &
                             LocationName, read_location, operator(/=),        &
                             write_location
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind,     &
                             get_obs_def_location, read_obs_def,               &
                             set_obs_def_time, get_obs_kind
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name,                 &
                             get_obs_kind_index, read_obs_kind
use time_manager_mod, only : time_type, operator(>), print_time, set_time,     &
                             print_date, set_calendar_type,                    &
                             operator(==), get_calendar_type, NO_CALENDAR,     &
                             operator(-), set_time_missing, operator(<)
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=), get_obs_def,             &
                             init_obs_sequence, static_init_obs_sequence,      &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_last_obs, get_next_obs,        &
                             insert_obs_in_seq, get_num_copies, get_num_qc,    &
                             get_copy_meta_data, get_qc_meta_data,             &
                             set_copy_meta_data, set_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence,                &
                             delete_seq_head, delete_seq_tail,                 &
                             get_num_key_range, get_obs_key, get_qc,           &
                             copy_partial_obs, get_next_obs_from_key,          &
                             get_obs_def, set_obs_def

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=128), parameter :: id  = "$Id$"

type(obs_sequence_type) :: seq_in
type(obs_type)          :: obs_in, next_obs_in
logical                 :: is_this_last
integer                 :: size_seq_in
integer                 :: num_copies_in, num_qc_in
integer                 :: iunit, io, i
integer                 :: max_num_obs, file_id
character(len=129)      :: read_format
logical                 :: pre_I_format, cal
character(len=512)      :: msgstring, msgstring1, msgstring2, msgstring3
type(obs_def_type)      :: this_obs_def
type(time_type)         :: obs_time

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 5000

! lazy, pick big number.  make it bigger if too small.
integer, parameter :: max_obs_input_types = 600

type obs_info_type
   integer :: count
   type(time_type) :: first_time
   type(time_type) :: last_time
end type

! in spite of the name, this is the number of specific types.
! also one for all obs types.
type(obs_info_type) :: oinfo(max_obs_kinds)
type(obs_info_type) :: all_obs

type(location_type) :: location
integer :: obs_type_ind
character(len=256) :: string

!----------------------------------------------------------------

! calendar is fixed at gregorian at this point, and the filename
! is read from the console/stdin

character(len=512)   :: filename_in = ''
character(len=32)    :: calendar      = 'Gregorian'


!----------------------------------------------------------------
! Start of the program:
!
! Process the input observations in the input sequence file
!----------------------------------------------------------------

call setup()

! NO NAMELIST!  read from the console
print *, 'Enter name of input obs_seq file: '
read *, filename_in

! only support a gregorian calendar.  if you are using a different type
! edit the line above and recompile.  this only affects the output
! format for the first/last obs times.
call set_calendar_type(calendar)

! set a logial to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! initialize the bookkeeping structures
do i=1, max_obs_kinds
   call initialize(oinfo(i))
enddo
call initialize(all_obs)


! single pass algorithm (unlike other obs tools).

call read_obs_seq_header(filename_in, num_copies_in, num_qc_in, &
   size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
   close_the_file = .true.)

if (max_num_obs == 0) then
   write(msgstring,*) 'No obs in input sequence file ', trim(filename_in)
   call error_handler(E_ERR,'obs_info',msgstring)
endif

write(msgstring, *) 'Starting to process input sequence file: '
write(msgstring1,*)  trim(filename_in)
call error_handler(E_MSG,'obs_info',msgstring, &
                   text2=msgstring1)

call read_obs_seq(filename_in, 0, 0, 0, seq_in)

! sanity check - ensure the linked list times are in increasing time order
call validate_obs_seq_time(seq_in, filename_in)

! blank line
call error_handler(E_MSG,' ',' ')

! Initialize individual observation variables
call init_obs(     obs_in,  num_copies_in, num_qc_in)
call init_obs(next_obs_in,  num_copies_in, num_qc_in)
   
!-------------------------------------------------------------
! Start of obs loop
!
!--------------------------------------------------------------

if ( get_first_obs(seq_in, obs_in) )  then

   is_this_last = .false.
   next_obs_in = obs_in

   ObsLoop : do while ( .not. is_this_last )

      obs_in = next_obs_in

      ! get obs_def info
      call get_obs_def(obs_in, this_obs_def)
      location = get_obs_def_location(this_obs_def)
      obs_type_ind = get_obs_kind(this_obs_def)
      obs_time = get_obs_def_time(this_obs_def)

      call update(all_obs, obs_time)
      call update(oinfo(obs_type_ind), obs_time)

      call write_location(0, location, charstring = string)
      !write(*, *) trim(string) // '  ' // trim(get_obs_kind_name(obs_type_ind))


      call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

   enddo ObsLoop

else
   write(msgstring, *)'no first observation in ',trim(filename_in)
   call error_handler(E_MSG,'obs_info', msgstring)
endif

print*, 'Totals for all obs types:'
print*, '  Count: ', all_obs%count
call print_date(all_obs%first_time, '.  First obs:')
call print_date(all_obs%last_time,  '.   Last obs:')
print*, '---------------------------------------------------------'


! print out the results
ALLTYPES: do i=1, max_obs_kinds
   if (oinfo(i)%count == 0) cycle ALLTYPES
   write(msgstring, '(A,I8)') get_obs_kind_name(i), oinfo(i)%count
   call error_handler(E_MSG, '', msgstring)
   call print_date(oinfo(i)%first_time, '.  First obs:')
   call print_date(oinfo(i)%last_time,  '.   Last obs:')
enddo ALLTYPES

call destroy_obs_sequence(seq_in)
call destroy_obs(     obs_in )
call destroy_obs(next_obs_in )

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_info')
call register_module(source, revision, revdate)
call static_init_obs_sequence()

end subroutine setup


!---------------------------------------------------------------------
subroutine shutdown()

call finalize_utilities('obs_info')

end subroutine shutdown


!---------------------------------------------------------------------
subroutine initialize(op)
 
! set everything to 0 or missing

type(obs_info_type), intent(inout) :: op

op%count = 0

op%first_time = set_time_missing()
op%last_time = set_time_missing()

end subroutine initialize

!---------------------------------------------------------------------
subroutine update(op, otime)
 
! add one to the count and set the time if either are outliers

type(obs_info_type), intent(inout) :: op
type(time_type),     intent(in)    :: otime

op%count = op%count + 1

if (op%first_time == set_time_missing()) then
   op%first_time = otime
else
   if (otime < op%first_time) op%first_time = otime
endif

if (op%last_time == set_time_missing()) then
   op%last_time = otime
else
   if (otime > op%last_time) op%last_time = otime
endif

end subroutine update

!---------------------------------------------------------------------
subroutine print_obs_seq(seq_in, filename)

! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

type(obs_sequence_type), intent(in) :: seq_in
character(len=*),        intent(in) :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in
integer                 :: i
integer                 :: this_obs_kind
! max_obs_kinds is a public from obs_kind_mod.f90 and really is
! counting the max number of types, not kinds
integer                 :: type_count(max_obs_kinds), identity_count


! Initialize input obs_types
do i = 1, max_obs_kinds
   type_count(i) = 0
enddo
identity_count = 0

! make sure there are obs left to process before going on.
! num_obs should be ok since we just constructed this seq so it should
! have no unlinked obs.  if it might for some reason, use this instead:
! size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_info',msgstring)
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
   call error_handler(E_MSG,'obs_info', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
! does not work with NO_CALENDAR
if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_kind = get_obs_kind(this_obs_def)
   if (this_obs_kind < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_kind) = type_count(this_obs_kind) + 1
   endif
!   print *, 'obs kind index = ', this_obs_kind
!   if(this_obs_kind > 0)print *, 'obs name = ', get_obs_kind_name(this_obs_kind)

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
do i = 1, max_obs_kinds
   if (type_count(i) > 0) then 
      write(msgstring, '(a32,i8,a)') trim(get_obs_kind_name(i)), &
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
   call error_handler(E_MSG,'obs_info:validate',msgstring)
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
   call error_handler(E_ERR,'obs_info:validate', &
                      msgstring, source, revision, revdate)
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
      call error_handler(E_ERR,'obs_info:validate', msgstring2, &
                         source, revision, revdate, &
                         text2=msgstring1)
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
   call error_handler(E_MSG,'obs_info:validate', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'obs_info:validate', msgstring, &
                         source, revision, revdate, &
                         text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'obs_info:validate', msgstring, &
                         source, revision, revdate, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time


!---------------------------------------------------------------------
subroutine print_metadata(seq, fname)

!
! print out the metadata strings, trimmed
!

type(obs_sequence_type),    intent(in) :: seq
character(len=*), optional, intent(in) :: fname

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(    seq)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring3,*)' illegal copy or obs count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_info', msgstring3, &
                      source, revision, revdate)
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


!---------------------------------------------------------------------
end program obs_info

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
