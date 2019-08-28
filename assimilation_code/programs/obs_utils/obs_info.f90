! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!> print out information about observation sequence file(s).
!> summarizes obs types, times, counts.

program obs_info

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : register_module, initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             open_file, close_file, finalize_utilities
use   parse_args_mod, only : get_args_from_string
use     location_mod, only : location_type, get_location, set_location,        &
                             write_location
use      obs_def_mod, only : obs_def_type, get_obs_def_time,                   &
                             get_obs_def_type_of_obs, get_obs_def_location
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs
use time_manager_mod, only : time_type, operator(>), print_time, set_time,     &
                             print_date, set_calendar_type,                    &
                             operator(==), get_calendar_type, NO_CALENDAR,     &
                             operator(-), set_time_missing, operator(<)
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=), get_obs_def,             &
                             static_init_obs_sequence,                         &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_next_obs,                      &
                             get_num_copies, get_num_qc,                       &
                             get_copy_meta_data, get_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence,                &
                             get_num_key_range, get_obs_key, get_qc,           &
                             get_obs_def

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
integer                 :: iunit, io, i, fnum
integer                 :: num_input_files = 0
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
type(obs_info_type) :: oinfo(0:max_defined_types_of_obs)
type(obs_info_type) :: all_obs

type(location_type) :: location
integer :: obs_type_ind
character(len=256) :: string

!----------------------------------------------------------------
! Namelist input with default values

integer, parameter :: MAX_IN_FILES = 5000

character(len=256)   :: filename_in(MAX_IN_FILES) = ''
character(len=256)   :: filelist_in = ''
character(len=32)    :: calendar = 'Gregorian'
logical              :: filenames_from_terminal = .false.


namelist /obs_info_nml/ filename_in, filelist_in, &
                        calendar, filenames_from_terminal

!----------------------------------------------------------------
! Start of the program:
!
! Process the input observations in the input sequence file
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_info_nml", iunit)
read(iunit, nml = obs_info_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_info_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_info_nml)
if (do_nml_term()) write(     *     , nml=obs_info_nml)

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! set a logial to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! after this call, filelist_in() has the names
if (filenames_from_terminal) then
   call parse_filenames_from_stdin(num_input_files, filename_in)
else
   call handle_filenames(filename_in, filelist_in, num_input_files)
endif

! end of namelist processing and setup

! for each file...

do fnum = 1, num_input_files

   ! initialize the bookkeeping structures
   do i=0, max_defined_types_of_obs
      call initialize(oinfo(i))
   enddo
   call initialize(all_obs)
   
   ! single pass algorithm (unlike other obs tools).
   
   call read_obs_seq_header(filename_in(fnum), num_copies_in, num_qc_in, &
      size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)
   
   if (max_num_obs == 0) then
      write(msgstring,*) 'No obs in input sequence file ', trim(filename_in(fnum))
      call error_handler(E_ERR,'obs_info',msgstring)
   endif
   
   write(msgstring, *) 'Starting to process input sequence file: '
   write(msgstring1,*)  trim(filename_in(fnum))
   call error_handler(E_MSG,'obs_info',msgstring, &
                      text2=msgstring1)
   
   call read_obs_seq(filename_in(fnum), 0, 0, 0, seq_in)
   
   ! sanity check - ensure the linked list times are in increasing time order
   call validate_obs_seq_time(seq_in, filename_in(fnum))
   
   ! blank line
   call error_handler(E_MSG,' ',' ')
   
   ! Initialize individual observation variables
   call init_obs(     obs_in,  num_copies_in, num_qc_in)
   call init_obs(next_obs_in,  num_copies_in, num_qc_in)
      
   !-------------------------------------------------------------
   ! Start of obs loop
   !--------------------------------------------------------------
   
   if ( get_first_obs(seq_in, obs_in) )  then
   
      is_this_last = .false.
      next_obs_in = obs_in
   
      ObsLoop : do while ( .not. is_this_last )
   
         obs_in = next_obs_in
   
         ! get obs_def info
         call get_obs_def(obs_in, this_obs_def)
         location = get_obs_def_location(this_obs_def)
         obs_type_ind = get_obs_def_type_of_obs(this_obs_def)
         obs_time = get_obs_def_time(this_obs_def)
   
         call update(all_obs, obs_time)
         call update(oinfo(obs_type_ind), obs_time)
   
         call write_location(0, location, charstring = string)
         !write(*, *) trim(string) // '  ' // trim(get_name_for_type_of_obs(obs_type_ind))
   
         call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)
   
      enddo ObsLoop
   
   else
      write(msgstring, *)'no first observation in ',trim(filename_in(fnum))
      call error_handler(E_MSG,'obs_info', msgstring)
   endif
   
   print*, 'Totals for all obs types:'
   print*, '  Count: ', all_obs%count
   call print_date(all_obs%first_time, '.  First obs:')
   call print_date(all_obs%last_time,  '.   Last obs:')
   print*, '---------------------------------------------------------'
   
   
   ! print out the results
   ALLTYPES: do i=0, max_defined_types_of_obs
      if (oinfo(i)%count == 0) cycle ALLTYPES
      write(msgstring, '(A,I8)') get_name_for_type_of_obs(i), oinfo(i)%count
      call error_handler(E_MSG, '', msgstring)
      call print_date(oinfo(i)%first_time, '.  First obs:')
      call print_date(oinfo(i)%last_time,  '.   Last obs:')
   enddo ALLTYPES
   
   call destroy_obs_sequence(seq_in)
   call destroy_obs(     obs_in )
   call destroy_obs(next_obs_in )
   
enddo

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_info')
call register_module(source,revision,revdate)
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
integer                 :: this_obs_type
integer                 :: type_count(max_defined_types_of_obs), identity_count


! Initialize input obs_types
do i = 0, max_defined_types_of_obs
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
subroutine parse_filenames_from_stdin(num_in, filenames)
integer,          intent(out) :: num_in
character(len=*), intent(out) :: filenames(:)

character(len=512) :: inline

! let the user know, if they don't type anything, why we are
! waiting for input.  comment this out if it gets annoying.
print *, 'reading input filename(s) from terminal'

read (*, '(A512)'), inline
call get_args_from_string(inline, num_in, filenames)

end subroutine parse_filenames_from_stdin

!---------------------------------------------------------------------
subroutine handle_filenames(filename_in, filelist_in, num_input_files)
! sort out the input lists, set the length if not given by user,
! make sure what's specified is consistent.
character(len=*), intent(inout) :: filename_in(:)
character(len=*), intent(in)    :: filelist_in
integer,          intent(inout) :: num_input_files

integer :: index
logical :: from_file
character(len=32) :: fsource

! ok, here's the new logic:
! if the user specifies neither filename_in nor filelist_in, we
! default to trying 'obs_seq.out' and num files = 1.
! if the user specifies both, it's an error.
! if the user gives a filelist, we make sure the length is not more
!   than maxfiles and read it into the explicit list and continue.
! if num_input_files = 0, we count up the non-null items in the list
!   and set num_input_files to that.
! if the user specifies num_input_files but it doesn't match the list length,
!   we give an error (and maybe suggest 0 if they don't want to keep
!   updating the num_input_files.)

! default case - input file is 'obs_seq.out' and count is 1.
if (filename_in(1) == '' .and. filelist_in == '') then

   if (num_input_files /= 0 .and. num_input_files /= 1) then
      call error_handler(E_ERR,'obs_info', &
          'if no filenames specified, num_input_files must be 0 or 1', &
          source,revision,revdate)
   endif

   num_input_files = 1
   filename_in(1) = 'obs_seq.out'
   return
endif

! make sure the namelist specifies one or the other but not both
if (filename_in(1) /= '' .and. filelist_in /= '') then
   call error_handler(E_ERR,'obs_info', &
       'cannot specify both filename_in and filelist_in', &
       source,revision,revdate)
endif

! if they have specified a file which contains a list, read it into
! the filename_in array and set the count.
if (filelist_in /= '') then
   fsource = 'filelist_in'
   from_file = .true.
else
   fsource = 'filename_in'
   from_file = .false.
endif

do index = 1, MAX_IN_FILES
   if (from_file) &
      filename_in(index) = get_next_filename(filelist_in, index)

   if (filename_in(index) == '') then
      if (index == 1) then
         call error_handler(E_ERR,'obs_info', &
             'namelist item '//trim(fsource)//' contains no filenames', &
             source,revision,revdate)
      endif
      ! leaving num_input_files unspecified (or set to 0) means use
      ! whatever number of files is in the list.
      if (num_input_files == 0) then
         num_input_files = index - 1
         return
      else
         ! if they do give a count, make it match.
         if (num_input_files == (index - 1)) return

         write(msgstring, *) 'num_input_files is ', num_input_files, &
                     ' but namelist item '//trim(fsource)//' has filecount ', index - 1

         write(msgstring2, *) 'if num_input_files is 0, the number of files will be automatically computed'
         write(msgstring3, *) 'if num_input_files is not 0, it must match the number of filenames specified'

         call error_handler(E_ERR,'obs_info', msgstring, &
            source,revision,revdate, text2=msgstring2, text3=msgstring3)

      endif
   endif
enddo

write(msgstring, *) 'cannot specify more than ',MAX_IN_FILES,' files'
call error_handler(E_ERR,'obs_info', msgstring, &
     source,revision,revdate)

end subroutine handle_filenames


!---------------------------------------------------------------------
end program obs_info

