! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Prints out a quick table of obs types and counts, overall start and
!> stop times, and metadata strings and counts.
!>
!> You can also get more/different info by running the obs_diag program.
!>
!> This program works on any obs_seq file (e.g. input to perfect_model_obs,
!> input to filter, output from filter).  It can print out the first/last
!> obs times in the file, and the count of each type of obs.
!>
!> If the DART QC data is available it can print out the total number
!> of obs for each QC category, and it can print out the QC counts per
!> obs type.  This program is intended to be a quick look at an unknown
!> obs_seq file to see what types and how many obs it contains. if it
!> has been through filter, it prints what happened to the obs (assimilated, 
!> failed forward operator, all rejected by outlier threshold, etc).
!>
!> It has an option to open a separate output file and write single line 
!> summaries of the QC values.  This output is intended to be ingested into
!> matlab or other plotting programs. (pro tip: see the matlab 'cell'
!> documentation for how to read in data that contains strings as well as
!> numbers.)  
!>
!> A missing piece for this purpose is that obs_seq.final files do not 
!> contain the analysis time of the model itself.  This may be what you 
!> want to use for the X axis on a plot.  Tim suggested some use of the 
!> schedule module might make it possible to generate periodic analysis 
!> times although they'd have to be sync'd with the input filenames.  
!> This is still TBD.

program obs_assim_count

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             open_file, close_file, finalize_utilities,        &
                             file_exist
use     location_mod, only : location_type, get_location, set_location,        &
                             LocationName, read_location, operator(/=),        &
                             write_location
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs,&
                             get_obs_def_location, read_obs_def,               &
                             set_obs_def_time
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs, &
                             get_index_for_type_of_obs, read_type_of_obs_table
use time_manager_mod, only : time_type, operator(>), print_time, set_time,     &
                             print_date, set_calendar_type,                    &
                             operator(/=), get_calendar_type, NO_CALENDAR,     &
                             operator(-)
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

character(len=*), parameter :: source   = 'obs_assim_count.f90'

type(obs_sequence_type) :: seq_in
integer                 :: size_seq_in
integer                 :: num_copies_in, num_qc_in
integer                 :: iunit, io, ifile
integer                 :: max_num_obs, file_id
character(len = 32)     :: read_format
character(len = 256)    :: filename_in
logical                 :: pre_I_format, cal
character(len = 512)    :: msgstring, msgstring1, msgstring2

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 5000

! lazy, pick big number.  make it bigger if too small.
integer, parameter :: max_obs_input_types = 800

! 8 is now failed vert convert
integer, parameter :: MAX_DART_QC = 8

integer, parameter :: QC_STR_LEN = 20

character(len=QC_STR_LEN) :: qc_meaning(0:MAX_DART_QC) = &
 (/ "Assimilated_OK      ", &
    "Evaluated_OK        ", &
    "Prior_only_Assim_OK ", &
    "Prior_only_Eval_OK  ", &
    "Forward_Op_failed   ", &
    "Namelist_Excluded   ", &
    "Bad_incoming_QC     ", &
    "Outlier_Failed      ", &
    "Bad_Vert_Convert    " /)

! if the user gives us an output file and asks for counts_only,
! print some stuff to the log & stdout too, but open this file
! and print the counts there as well.
integer :: output_unit = -1

! first time through a loop 
logical :: first_time

!----------------------------------------------------------------
! Namelist input with default values

character(len = 256) :: obs_sequence_name = ''
character(len = 256) :: obs_sequence_list = ''
character(len = 256) :: output_file = ''

logical :: stats_by_obs_type = .false.
logical :: counts_only       = .false.
logical :: print_metadata    = .false.
logical :: print_time_info   = .true.

character(len=32) :: calendar    = 'Gregorian'


namelist /obs_assim_count_nml/ &
         obs_sequence_name, obs_sequence_list, output_file, &
         calendar, stats_by_obs_type, counts_only, print_metadata, &
         print_time_info

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, collecting
! and printing stats from each input file.
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_assim_count_nml", iunit)
read(iunit, nml = obs_assim_count_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_assim_count_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_assim_count_nml)
if (do_nml_term()) write(     *     , nml=obs_assim_count_nml)

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! set a logical to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! Check the user input for sanity
if ((obs_sequence_name /= '') .and. (obs_sequence_list /= '')) then
   write(msgstring1,*)'specify "obs_sequence_name" or "obs_sequence_list"'
   write(msgstring2,*)'set other to an empty string ... i.e. ""'
   call error_handler(E_ERR, 'obs_assim_count', msgstring1, &
                   source, text2=msgstring2)
endif

if (counts_only) &
   call error_handler(E_MSG, 'obs_assim_count', 'printing summary counts only')

if (output_file /= '') &
   output_unit = open_file(output_file, action='write', form='formatted')

! if you add anything to the namelist, you can process it here.

! end of namelist processing and setup

ifile = 0
ObsFileLoop : do      ! until we run out of names

   ifile = ifile + 1

   if (obs_sequence_list == '') then
      if (ifile > 1) exit ObsFileLoop
      filename_in = obs_sequence_name
   else
      filename_in = get_next_filename(obs_sequence_list,ifile)
      if (filename_in == '') exit ObsFileLoop
   endif

   if ( .not. file_exist(filename_in) ) then
      write(msgstring1,*)'cannot open ', trim(filename_in)
      call error_handler(E_ERR,'obs_assim_count:',msgstring1,source)
   endif

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   call read_obs_seq_header(filename_in, num_copies_in, num_qc_in, &
      size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)
   
   if (max_num_obs == 0) then
      write(msgstring,*) 'No obs in input sequence file ', trim(filename_in)
      call error_handler(E_MSG,'obs_assim_count',msgstring)
      cycle ObsFileLoop
   endif
   
   write(msgstring, *) 'Starting to process input sequence file: '
   write(msgstring1,*)  trim(filename_in)
   call error_handler(E_MSG,'obs_assim_count',msgstring, &
                      text2=msgstring1)
   
   call read_obs_seq(filename_in, 0, 0, 0, seq_in)
   
   ! sanity check - ensure the linked list times are in increasing time order
   call validate_obs_seq_time(seq_in, filename_in)
   
   ! the counting up is done here now.
   call compute_and_print_obs_seq_info(seq_in, filename_in, output_unit)
   
   ! clean up
   call destroy_obs_sequence(seq_in)

enddo ObsFileLoop

! see comment where the output_unit variable is declared
if (output_file /= '') &
   call close_file(output_unit)

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_assim_count')
call static_init_obs_sequence()

end subroutine setup


!---------------------------------------------------------------------
subroutine shutdown()

call finalize_utilities('obs_assim_count')

end subroutine shutdown


!---------------------------------------------------------------------
! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

! and this one counts up, if there is a 'posterior ensemble mean' copy,
! how many are missing_r8 and how many are not.   it could also count
! up the 'DART quality control' settings?  start with the latter for now.

subroutine compute_and_print_obs_seq_info(seq_in, filename, ounit)

type(obs_sequence_type), intent(in) :: seq_in
character(len=*),        intent(in) :: filename
integer,                 intent(in) :: ounit

type(obs_type)     :: obs, next_obs
type(obs_def_type) :: this_obs_def
logical            :: is_there_one, is_this_last
integer            :: size_seq_in
integer            :: i, qc_int, obt
integer            :: this_obs_type
integer            :: type_count(0:max_defined_types_of_obs), identity_count, qc_count(0:MAX_DART_QC), qcindex
integer            :: qc_count_by_type(0:MAX_DART_QC, 0:max_defined_types_of_obs)
real(r8)           :: qcval(1)
character(len=32)  :: this_obs_name


! Initialize counters
type_count(:) = 0
identity_count = 0
qc_count(:) = 0
qc_count_by_type(:,:) = 0


size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty, skipping.'
   call error_handler(E_MSG,'obs_assim_count',msgstring)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq_in), get_num_qc(seq_in))
call init_obs(next_obs, get_num_copies(seq_in), get_num_qc(seq_in))

! find the dart qc copy, if there is one
qcindex = get_dartqc_index(seq_in, filename)

! blank line
call error_handler(E_MSG,'',' ')

write(msgstring,*) 'Processing sequence file: ', trim(filename)
call error_handler(E_MSG,'',msgstring)

if (print_metadata) &
   call print_seq_metadata(seq_in, filename)

!-------------------------------------------------------------
! Start to process obs from seq_in
!--------------------------------------------------------------
is_there_one = get_first_obs(seq_in, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in ',trim(filename)
   call error_handler(E_MSG,'obs_assim_count', msgstring)
endif

! process it here
is_this_last = .false.

! blank line
call error_handler(E_MSG, '', ' ')

call get_obs_def(obs, this_obs_def)
if (print_time_info) then
   call print_time(get_obs_def_time(this_obs_def), ' First obs time: ')
   ! should, but does not work with NO_CALENDAR
   if (cal) call print_date(get_obs_def_time(this_obs_def), ' First obs date: ')
endif

! collect data counts
ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_type = get_obs_def_type_of_obs(this_obs_def)
   if (this_obs_type < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_type) = type_count(this_obs_type) + 1
   endif
   if (qcindex > 0) then
      call get_qc(obs, qcval, qcindex)
      qc_int = nint(qcval(1))
      qc_count(qc_int) = qc_count(qc_int) + 1
      if (this_obs_type >= 0) &
         qc_count_by_type(qc_int,this_obs_type) = qc_count_by_type(qc_int,this_obs_type) + 1
   endif

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
      if (print_time_info) then
         call print_time(get_obs_def_time(this_obs_def), '  Last obs time: ')
         if (cal) call print_date(get_obs_def_time(this_obs_def), '  Last obs date: ')
         call error_handler(E_MSG, '', ' ')
      endif
   endif

enddo ObsLoop

! how to output - depends on 'counts_only', 'stats_by_obs_type' settings and
! whether input file has DART QCs or not.

! FIXME: this is very messy.  get it working and then refactor so there 
! isn't such a deep set of nested if/else/thens.  the cases are:
!  - are you asking for counts only or a full summary?
!  - are you asking for stats by obs type or only totals?
!  - does this file contain DART QC info?
! which makes 6 different code blocks below, each with slightly
! different output strings and formats.

if (counts_only) then
   if (qcindex > 0) then
      if (stats_by_obs_type) then
         first_time = .true.
         do obt = 0, max_defined_types_of_obs
            if (sum(qc_count_by_type(:,obt)) <= 0) cycle
            this_obs_name = get_name_for_type_of_obs(obt)
   
            ! print a header, and then try to make the rest of the lines have
            ! all the info on a single line.  (easier to grep, and to read into
            ! matlab for stats, i hope.)
            if (first_time) call print_header('DART QC results:', &
               'Obs_type, QC value, QC meaning, count_of_that_QC_value, total_obs_of_this_type', &
               first_time)

            do i=0, MAX_DART_QC
               if (qc_count_by_type(i,obt) > 0) then
                  write(msgstring, '(a34,a4,i1,a22,2(i12))') trim(this_obs_name), ' QC ', i, qc_meaning(i), &
                                                          qc_count_by_type(i,obt), sum(qc_count_by_type(:,obt))
                  call error_handler(E_MSG, '', msgstring)
                  if (ounit >= 0) write(ounit, '(A)') trim(msgstring)
               endif
            enddo
         enddo
      else
         first_time = .true.

         ! print a header, and then try to make the rest of the lines have
         ! all the info on a single line.  (easier to grep, and to read into
         ! matlab for stats, i hope.)
         if (first_time) call print_header('DART QC results:', &
            'QC value, count_of_that_QC_value, total_obs', &
            first_time)

         do i=0, MAX_DART_QC
            if (qc_count(i) > 0) then
               write(msgstring, '(a22,a4,i1,2(i12))') qc_meaning(i), ' QC ', i, qc_count(i), sum(qc_count(:))
               call error_handler(E_MSG, '', msgstring)
               if (ounit >= 0) write(ounit, '(A)') trim(msgstring)
            endif
         enddo
      endif
   else
      if (stats_by_obs_type) then
         first_time = .true.
         do obt = 0, max_defined_types_of_obs
            if (type_count(obt) <= 0) cycle
            this_obs_name = get_name_for_type_of_obs(obt)
   
            ! print a header, and then try to make the rest of the lines have
            ! all the info on a single line.  (easier to grep, and to read into
            ! matlab for stats, i hope.)
            if (first_time) call print_header('DART Obs results:', &
               'Obs_type, total_obs_of_this_type', &
               first_time)

            write(msgstring, '(a34,1(i12))') trim(this_obs_name), type_count(obt)
            call error_handler(E_MSG, '', msgstring)
            if (ounit >= 0) write(ounit, '(A)') trim(msgstring)
         enddo
      else
         first_time = .true.
         call print_header('DART QC results:', &
            'QC value, count_of_that_QC_value, total_obs', &
            first_time)

         do i=0, MAX_DART_QC
            if (qc_count(i) > 0) then
               write(msgstring, '(3(i12))') i, qc_count(i), sum(qc_count(:))
               call error_handler(E_MSG, '', msgstring)
               if (ounit >= 0) write(ounit, '(A)') trim(msgstring)
            endif
         enddo
      endif
   endif
else
   
   write(msgstring, *) 'Number of obs processed  :          ', size_seq_in
   call error_handler(E_MSG, '', msgstring)
   write(msgstring, *) '---------------------------------------------------------'
   call error_handler(E_MSG, '', msgstring)
   do i = 1, max_defined_types_of_obs
      if (type_count(i) > 0) then 
         write(msgstring, '(a32,i12,a)') trim(get_name_for_type_of_obs(i)), &
                                        type_count(i), ' obs'
         call error_handler(E_MSG, '', msgstring)
      endif
   enddo
   if (identity_count > 0) then 
      write(msgstring, '(a32,i12,a)') 'Identity observations', &
                                     identity_count, ' obs'
      call error_handler(E_MSG, '', msgstring)
   endif
   if (qcindex > 0) then
      call error_handler(E_MSG, '', ' ')
      write(msgstring, *) 'DART QC results for Total obs: ', sum(qc_count(:))
      call error_handler(E_MSG, '', msgstring)
      do i=0, MAX_DART_QC
         if (qc_count(i) > 0) then
            write(msgstring, '(a22,a4,i1,i12)') qc_meaning(i), ' QC ', i, qc_count(i)
            call error_handler(E_MSG, '', msgstring)
         endif
      enddo
   
      if (stats_by_obs_type) then
         do obt = 0, max_defined_types_of_obs
            if (sum(qc_count_by_type(:,obt)) <= 0) cycle
            this_obs_name = get_name_for_type_of_obs(obt)
   
            call error_handler(E_MSG, '', ' ')
            write(msgstring, *) 'Total '//trim(this_obs_name)//' obs: ', sum(qc_count_by_type(:,obt))
            call error_handler(E_MSG, '', msgstring)
            do i=0, MAX_DART_QC
               if (qc_count_by_type(i,obt) > 0) then
                  write(msgstring, '(a22,a4,i1,i12)') qc_meaning(i), ' QC ', i, &
                                                      qc_count_by_type(i,obt)
                  call error_handler(E_MSG, '', msgstring)
               endif
            enddo
         enddo
      endif
   endif
   
   ! another blank line
   call error_handler(E_MSG, '', ' ')
endif

! Time to clean up

call destroy_obs(     obs)
call destroy_obs(next_obs)

end subroutine compute_and_print_obs_seq_info


!---------------------------------------------------------------------

! something like this eventually belongs in the obs_seq_mod code, but 
! for now try it out here.  we just fixed a hole in the interactive create
! routine which would silently let you create out-of-time-order
! linked lists, which gave no errors but didn't assimilate the
! right obs at the right time when running filter.   this runs
! through the times in the entire sequence, ensuring they are
! monotonically increasing in time.  this should help catch any
! bad files which were created with older versions of code, or
! created by users with python...

subroutine validate_obs_seq_time(seq, filename)

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
   call error_handler(E_MSG,'obs_assim_count:validate',msgstring)
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
   call error_handler(E_ERR,'obs_assim_count:validate', msgstring, source)
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
      call error_handler(E_ERR,'obs_assim_count:validate', msgstring2, &
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
   call error_handler(E_MSG,'obs_assim_count:validate', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'obs_assim_count:validate', msgstring, &
                         source, text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'obs_assim_count:validate', msgstring, &
                         source, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time


!---------------------------------------------------------------------
! print out the metadata strings, trimmed

subroutine print_seq_metadata(seq, fname)

type(obs_sequence_type), intent(in) :: seq
character(len=*) :: fname

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str
character(len=255) :: msgstring3

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(    seq)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring3,*)' illegal copy or obs count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_assim_count', msgstring3, source)
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

end subroutine print_seq_metadata


!---------------------------------------------------------------------
! return the index number of the dart qc copy (-1 if none)

function get_dartqc_index(seq, fname)

type(obs_sequence_type), intent(in) :: seq
character(len=*) :: fname
integer :: get_dartqc_index

integer :: num_qc, i
character(len=metadatalength) :: str
character(len=255) :: msgstring3

num_qc = get_num_qc(seq)

if ( num_qc < 0 ) then
   write(msgstring3,*)' illegal qc metadata count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_assim_count', msgstring3, source)
endif

QCMetaData : do i=1, num_qc
   str = get_qc_meta_data(seq,i)

   if (str == 'DART quality control') then
      get_dartqc_index = i
      return 
   endif

enddo QCMetaData

get_dartqc_index = -1

end function get_dartqc_index

!---------------------------------------------------------------------

subroutine print_header(string1, string2, first)

character(len=*), intent(in)    :: string1
character(len=*), intent(in)    :: string2
logical,          intent(inout) :: first
 
if (.not. first) return

call error_handler(E_MSG, '', ' ')
call error_handler(E_MSG, '', string1)
call error_handler(E_MSG, '', string2)
call error_handler(E_MSG, '', ' ')

first = .false.

end subroutine print_header

!---------------------------------------------------------------------

end program obs_assim_count

