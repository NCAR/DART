! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> print out information about observation sequence file(s).
!> summarizes obs types, times, counts.
!>
!> If 'csv_style_output = .true.' and 'output_file="summary.csv"', 
!> the resulting comma-separated-values file
!> can be trivially read into matlab with:
!>
!> X = readtable('summary.csv')
!>
!> The first line of the output file describes the columns.
!>
!>@todo FIXME This routine should use print_obs_seq_summary. 

program obs_info

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             open_file, close_file, finalize_utilities
use   parse_args_mod, only : get_args_from_string
use     location_mod, only : location_type
use      obs_def_mod, only : obs_def_type, get_obs_def_time,                   &
                             get_obs_def_type_of_obs, get_obs_def_location
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs
use time_manager_mod, only : time_type, operator(>), print_time, set_time,     &
                             print_date, set_calendar_type, operator(+),       &
                             operator(==), get_calendar_type, NO_CALENDAR,     &
                             operator(-), set_time_missing, operator(<),       &
                             operator(/), get_time, month_name, get_date
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

character(len=*), parameter :: source = 'obs_info.f90'

type(obs_sequence_type) :: seq_in
type(obs_type)          :: obs_in, next_obs_in
logical                 :: is_this_last
integer                 :: size_seq_in
integer                 :: num_copies_in, num_qc_in
integer                 :: iunit, io, i, fnum, ounit
integer                 :: num_input_files = 0
integer                 :: max_num_obs, file_id
character(len=129)      :: read_format
logical                 :: pre_I_format, cal
character(len=512)      :: msgstring, msgstring1, msgstring2, msgstring3
type(obs_def_type)      :: this_obs_def
type(time_type)         :: obs_time
!integer                 :: mid_day, mid_sec
character(len=32)       :: mid_string

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 5000

! lazy, pick big number.  make it bigger if too small.
integer, parameter :: max_obs_input_types = 600

type obs_info_type
   integer :: count
   type(time_type) :: first_time
   type(time_type) :: last_time
end type

! an array to hold counts for each obs type.
! also one for all obs types.
type(obs_info_type) :: oinfo(max_defined_types_of_obs)
type(obs_info_type) :: identity_obs
type(obs_info_type) :: all_obs

type(location_type) :: location
integer :: obs_type_ind

!----------------------------------------------------------------
! Namelist input with default values

integer, parameter :: MAX_IN_FILES = 5000

character(len=256)   :: filename_in(MAX_IN_FILES) = ''
character(len=256)   :: filelist_in = ''
character(len=32)    :: calendar = 'Gregorian'
logical              :: filenames_from_terminal = .false.
logical              :: csv_style_output = .false.
character(len=256)   :: output_file = ''


namelist /obs_info_nml/ filename_in, filelist_in, csv_style_output, &
                        calendar, filenames_from_terminal, output_file

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

if (output_file /= '') then
   ounit = open_file(output_file, action='write')
   write(msgstring, *) 'Output counts will be written to text file: '
   write(msgstring1,*)  trim(output_file)
   call error_handler(E_MSG,'obs_info',msgstring, text2=msgstring1)

   if (csv_style_output) then
      ! The first row of the csv file describes the columns.
      write(ounit,'(''Filename, ObsTypeInteger, ObsType, Count, AverageTime'')')
   endif
else
   ounit = 0
endif

! end of namelist processing and setup

! for each file...

do fnum = 1, num_input_files

   ! initialize the bookkeeping structures
   do i=1, max_defined_types_of_obs
      call initialize(oinfo(i))
   enddo
   call initialize(identity_obs)
   call initialize(all_obs)
   
   ! single pass algorithm (unlike other obs tools).
   
   call read_obs_seq_header(filename_in(fnum), num_copies_in, num_qc_in, &
      size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)
   
   if (max_num_obs == 0) then
      write(msgstring,*) 'No obs in input sequence file ', trim(filename_in(fnum))
      call error_handler(E_ERR,'obs_info',msgstring)
   endif
   
   if (.not. csv_style_output) then
      write(msgstring,  *) '--------------------------------------------------'
      write(msgstring1, *) 'Starting to process input sequence file: '
      write(msgstring2, *)  trim(filename_in(fnum))
      call error_handler(E_MSG,'obs_info',msgstring, &
                         text2=msgstring1, text3=msgstring2)
   endif
   
   call read_obs_seq(filename_in(fnum), 0, 0, 0, seq_in)
   
   ! sanity check - ensure the linked list times are in increasing time order
   call validate_obs_seq_time(seq_in, filename_in(fnum))
   
   ! blank line
   if (.not. csv_style_output) call error_handler(E_MSG,' ',' ')
   
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
         if (obs_type_ind < 0) then
            call update(identity_obs, obs_time)
         else
            call update(oinfo(obs_type_ind), obs_time)
         endif
   
         call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)
   
      enddo ObsLoop
   
   else
      write(msgstring, *)'no first observation in ',trim(filename_in(fnum))
      call error_handler(E_MSG,'obs_info', msgstring)
   endif
   
   if (.not. csv_style_output) then
      write(ounit, *) 'Totals for all obs types:'
      write(ounit, *) '  Count: ', all_obs%count
      call print_date(all_obs%first_time, '.  First obs:', ounit)
      call print_date(all_obs%last_time,  '.   Last obs:', ounit)
   endif
   
   ! print out the results
   ALLTYPES: do i=1, max_defined_types_of_obs
      if (oinfo(i)%count == 0) cycle ALLTYPES
      if (csv_style_output) then
         call compute_times(oinfo(i)%first_time, oinfo(i)%last_time, avg_string=mid_string)
         write(ounit, '(A,'','',I8,'','',A36,'','',I8,'','',A)') &
               trim(filename_in(fnum)), i, trim(get_name_for_type_of_obs(i)), &
               oinfo(i)%count, trim(mid_string)
      else
         write(ounit, '(A32,I8)') get_name_for_type_of_obs(i), oinfo(i)%count
         call print_date(oinfo(i)%first_time, '.  First obs:', ounit)
         call print_date(oinfo(i)%last_time,  '.   Last obs:', ounit)
      endif
   enddo ALLTYPES
   if (identity_obs%count > 0) then
      if (csv_style_output) then
         call compute_times(identity_obs%first_time, identity_obs%last_time, avg_string=mid_string)
         write(ounit, '(A,'','',I8,'','',A36,'','',I8,'','',A)') &
               trim(filename_in(fnum)), -1, 'IDENTITY_OBSERVATIONS', &
               identity_obs%count, trim(mid_string)
      else
         write(ounit, '(A32,I8)') "IDENTITY_OBSERVATIONS          ", identity_obs%count
         call print_date(identity_obs%first_time, '.  First obs:', ounit)
         call print_date(identity_obs%last_time,  '.   Last obs:', ounit)
      endif
   endif
   
   call destroy_obs_sequence(seq_in)
   call destroy_obs(     obs_in )
   call destroy_obs(next_obs_in )
   
   ! blank line only if not doing the CSV output
   if (.not. csv_style_output) call error_handler(E_MSG,' ',' ')

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
call static_init_obs_sequence()

end subroutine setup


!---------------------------------------------------------------------
subroutine shutdown()

call close_file(ounit)
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
   call error_handler(E_ERR,'obs_info:validate', msgstring, source)
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
      call error_handler(E_ERR,'obs_info:validate', msgstring2, source, text2=msgstring1)
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
      call error_handler(E_ERR,'obs_info:validate', msgstring, source, text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'obs_info:validate', msgstring, source, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time


!---------------------------------------------------------------------
subroutine parse_filenames_from_stdin(num_in, filenames)
integer,          intent(out) :: num_in
character(len=*), intent(out) :: filenames(:)

character(len=512) :: inline

! let the user know, if they don't type anything, why we are
! waiting for input.  comment this out if it gets annoying.
print *, 'reading input filename(s) from terminal'

read (*, *) inline
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
          'if no filenames specified, num_input_files must be 0 or 1', source)
   endif

   num_input_files = 1
   filename_in(1) = 'obs_seq.out'
   return
endif

! make sure the namelist specifies one or the other but not both
if (filename_in(1) /= '' .and. filelist_in /= '') then
   call error_handler(E_ERR,'obs_info', &
       'cannot specify both filename_in and filelist_in', source)
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
             'namelist item '//trim(fsource)//' contains no filenames', source)
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
            source, text2=msgstring2, text3=msgstring3)

      endif
   endif
enddo

write(msgstring, *) 'cannot specify more than ',MAX_IN_FILES,' files'
call error_handler(E_ERR,'obs_info', msgstring, source)

end subroutine handle_filenames

!---------------------------------------------------------------------

subroutine compute_times(first_time, last_time, avg_day, avg_sec, avg_string)
type(time_type),  intent(in)  :: first_time
type(time_type),  intent(in)  :: last_time
integer,          intent(out), optional :: avg_day
integer,          intent(out), optional :: avg_sec
character(len=*), intent(out), optional :: avg_string

type(time_type)  :: avg_time
integer          :: yr, mo, dy, hr, mn, sc
integer          :: aday, asec
character(len=9) :: mon_name

avg_time = (first_time + last_time) / 2

call get_time(avg_time, asec, aday)

if (present(avg_day) .and. present(avg_sec)) then
   avg_day = aday
   avg_sec = asec
else if (present(avg_sec)) then
   call get_time(avg_time, avg_sec)
endif

if (present(avg_string)) then
   if (cal) then
      call get_date(avg_time, yr,mo,dy,hr,mn,sc)
      mon_name = month_name(mo) 
      write(avg_string, "(I2.2,'-',A3,'-',I4,' ',I2.2,':',I2.2,':',I2.2 )") dy, mon_name, yr, hr, mn, sc 
   else 
      write(avg_string, "('day ', I8, ',', ' sec ', I8)") aday, asec
   endif
endif

end subroutine compute_times

!---------------------------------------------------------------------
end program obs_info

