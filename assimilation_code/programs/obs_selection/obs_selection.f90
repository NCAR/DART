! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> This specialized tool selects a subset of input observations
!> from an observation sequence file.

program obs_selection

!> nsc 12apr2012 -
!> was too slow for large lists and large obs_seq files.
!> sorted the obs_def list by time and then started the search
!> at the time of the next obs, and quit looping when past that time.
!> really speeded up - may not have to do anything more complicated
!> at this time.
!> if more speed is needed, next likely place to pick up speed is
!> by sorting all obs at the same time by locations (maybe just
!> by the x coord for starters), or binning in spatial bins if
!> still too slow.  but the current fixes should go a long way
!> to making the performance acceptable.

! this latest addition has select by list of obs types.

use        types_mod, only : r8, missing_r8, metadatalength
use    utilities_mod, only : timestamp, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit,   &
                             do_nml_file, do_nml_term, get_next_filename, &
                             open_file, close_file, finalize_utilities
use     location_mod, only : location_type, get_location, set_location, &
                             LocationName, read_location, operator(==), &
                             write_location, VERTISSURFACE, VERTISHEIGHT, &
                             VERTISLEVEL, VERTISPRESSURE, VERTISSCALEHEIGHT, &
                             VERTISUNDEF, query_location
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs, &
                             get_obs_def_location, read_obs_def
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs, get_index_for_type_of_obs, &
                             read_type_of_obs_table
use time_manager_mod, only : time_type, operator(>), print_time, set_time, &
                             print_date, set_calendar_type, GREGORIAN,     &
                             operator(/=), operator(<=), NO_CALENDAR,      &
                             get_calendar_type, time_index_sort, operator(==)
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq, &
                             init_obs, assignment(=), get_obs_def, &
                             init_obs_sequence, static_init_obs_sequence, &
                             read_obs_seq_header, read_obs_seq, get_num_obs, &
                             get_first_obs, get_last_obs, get_next_obs, &
                             insert_obs_in_seq, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, &
                             set_copy_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence, &
                             delete_seq_head, delete_seq_tail, &
                             get_num_key_range, get_obs_key,  &
                             copy_partial_obs, get_next_obs_from_key

! call the get_close routines?

implicit none

character(len=*), parameter :: source = 'obs_selection.f90'

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs_in, next_obs_in
type(obs_type)          :: obs_out, prev_obs_out
type(time_type)         :: t1, t2
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in, num_copies_in, num_qc_in
integer                 :: size_seq_out, num_copies_out, num_qc_out
integer                 :: num_inserted, iunit, io, i, j, ocount
integer                 :: total_num_inserted, base_index, next_base_index
integer                 :: max_num_obs, file_id, this_type
logical, allocatable    :: type_wanted(:)
integer                 :: first_seq, num_good_called, num_good_searched
character(len=129)      :: read_format, meta_data
logical                 :: pre_I_format, cal
character(len=512)      :: msgstring, msgstring1, msgstring2, msgstring3


!----------------------------------------------------------------
! Namelist input with default values

! max_num_input_files : maximum number of input sequence files to be processed
! lazy, pick big numbers.  make them bigger if too small.
integer, parameter               :: max_num_input_files = 1000
integer, parameter               :: max_obs_input_types = 500
integer                          :: num_input_files = 0
type(obs_def_type),  allocatable :: obs_def_list(:)
integer                          :: obs_def_count
integer                          :: print_every_nth_obs = 100


character(len=256) :: filename_seq(max_num_input_files) = ''
character(len=256) :: filename_seq_list = ''
character(len=256) :: filename_out  = 'obs_seq.processed'
logical            :: process_file(max_num_input_files)

character(len=256) :: selections_file = 'obsdef_mask.txt'
logical            :: selections_is_obs_seq = .false.

! max differences allowed when deciding a location is the same
real(r8) :: latlon_tolerance      = 0.000001 ! horizontal, degrees
real(r8) :: surface_tolerance     = 0.0001   ! vertical, meters
real(r8) :: pressure_tolerance    = 0.001    ! vertical, pascals
real(r8) :: height_tolerance      = 0.0001   ! vertical, meters
real(r8) :: scaleheight_tolerance = 0.001    ! vertical, pressure ratio
real(r8) :: level_tolerance       = 0.00001  ! vertical, fractional model levels

logical  :: match_vertical        = .false.
logical  :: print_only            = .false.
logical  :: partial_write         = .false.
logical  :: print_timestamps      = .false.
character(len=32) :: calendar     = 'Gregorian'


namelist /obs_selection_nml/ &
         num_input_files, filename_seq, filename_seq_list, filename_out, &
         selections_file, selections_is_obs_seq, print_only, calendar,   &
         print_timestamps, partial_write, latlon_tolerance,              &
         surface_tolerance, pressure_tolerance, height_tolerance,        &
         scaleheight_tolerance, level_tolerance, match_vertical

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call obs_seq_modules_used()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_selection_nml", iunit)
read(iunit, nml = obs_selection_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_selection_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_selection_nml)
if (do_nml_term()) write(     *     , nml=obs_selection_nml)

! ok, here's the new logic:
! if the user specifies neither filename_seq nor filename_seq_list, we
! default to trying 'obs_seq.out' and num files = 1.
! if the user specifies both, it's an error.
! if the user gives a filelist, we make sure the length is not more
!   than maxfiles and read it into the explicit list and continue.
! if num_input_files = 0, we count up the non-null items in the list
!   and set num_input_files to that.
! if the user specifies num_input_files but it doesn't match the list length,
!   we give an error (and maybe suggest 0 if they don't want to keep
!   updating the num_input_files.)

call handle_filenames(filename_seq, filename_seq_list, num_input_files)

! set the calendar type, and set a flag to say whether it is something
! other than no_calendar.  note that cal is set as a logical.
call set_calendar_type(calendar)
cal = (get_calendar_type() /= NO_CALENDAR)

call read_selection_list(selections_file, selections_is_obs_seq, obs_def_list, obs_def_count, type_wanted)

! end of namelist processing and setup

! some statistics for timing/debugging
num_good_called  = 0
num_good_searched = 0

! Read header information for the sequences to see if we need
! to accomodate additional copies or qc values from subsequent sequences.
! Also, calculate how many observations to be added to the first sequence.
num_copies_out   = 0
num_qc_out       = 0
size_seq_out     = 0

! TWO PASS algorithm; open each file, trim it if requested, and count
! the number of actual observations.  then the output file can be
! created with the correct size, and as observations are put into it
! they'll be sorted, and unused obs will be removed.

! pass 1:

if (print_timestamps) call timestamp(string1='start of pass1', pos='brief')

first_seq = -1
do i = 1, num_input_files

   if ( len(filename_seq(i)) .eq. 0 .or. filename_seq(i) .eq. "" ) then
      call error_handler(E_ERR,'obs_selection', &
         'num_input_files and filename_seq mismatch',source)
   endif

   ! count up the max number of observations possible if every obs in the
   ! input file was copied to the output.

   call read_obs_seq_header(filename_seq(i), num_copies_in, num_qc_in, &
      size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)

   if (max_num_obs == 0) then
      process_file(i) = .false.
      write(msgstring,*) 'No obs in input sequence file ', trim(filename_seq(i))
      call error_handler(E_MSG,'obs_selection',msgstring)
      cycle
   endif

   process_file(i) = .true.
   size_seq_in = max_num_obs

   ! keep counts of first file to compare with others.
   if ( first_seq < 0 ) then
      first_seq = i
      num_copies_out = num_copies_in
      num_qc_out = num_qc_in
      size_seq_out = size_seq_in
   else
      size_seq_out = size_seq_out + size_seq_in
   endif

enddo

if (print_timestamps) call timestamp(string1='end   of pass1', pos='brief')
 
! no valid obs found?  if the index value is still negative, we are
! still waiting to process the first one and never found one.
if (first_seq < 0 .or. size_seq_out == 0) then
   msgstring = 'All input files are empty'
   call error_handler(E_MSG,'obs_by_selection',msgstring,source)
endif


! pass 2:

if (print_timestamps) call timestamp(string1='start of pass2', pos='brief')

! blank line, start of actually creating output file
call error_handler(E_MSG,' ',' ')

! Initialize individual observation variables
call init_obs(     obs_in,  num_copies_in,  num_qc_in )
call init_obs(next_obs_in,  num_copies_in,  num_qc_in )
call init_obs(     obs_out, num_copies_out, num_qc_out)
call init_obs(prev_obs_out, num_copies_out, num_qc_out)

total_num_inserted = 0

! Read obs seq to be added, and insert obs from it to the output seq
first_seq = -1
FILES: do i = 1, num_input_files

   if (.not. process_file(i)) cycle

   write(msgstring, *) 'input file ', i
   if (print_timestamps) call timestamp(string1='start of '//trim(msgstring), pos='brief')

   write(msgstring,*) 'Starting to process input sequence file ', trim(filename_seq(i))
   call error_handler(E_MSG,'obs_selection',msgstring)
   call timestamp(' at ', pos='brief')

   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)

   ! create the output sequence here based on the first input file
   if (first_seq < 0) then
      call init_obs_sequence(seq_out, num_copies_out, num_qc_out, size_seq_out)
      do j=1, num_copies_out
         meta_data = get_copy_meta_data(seq_in, j)
         call set_copy_meta_data(seq_out, j, meta_data)
      enddo
      do j=1, num_qc_out
         meta_data = get_qc_meta_data(seq_in, j)
         call set_qc_meta_data(seq_out, j, meta_data)
      enddo
      first_seq = i
   else
      ! we have an existing output sequence already.  make sure the next one
      ! is completely compatible.

      ! Compare metadata between the observation sequences.
      ! This routine exits if they do not match.
      call compare_metadata(seq_out, seq_in, filename_seq(first_seq), filename_seq(i))
   endif

   size_seq_out = get_num_key_range(seq_out)   !current size of seq_out
   size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

   size_seq_out = get_num_key_range(seq_out)   !current size of seq_out
   size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

   ! ensure the linked list times are in increasing time order
   call validate_obs_seq_time(seq_in, filename_seq(i))

   if (print_only) call print_obs_seq(seq_in, filename_seq(i))

   !-------------------------------------------------------------
   ! Start to insert obs from sequence_in into sequence_out
   !
   ! NOTE: insert_obs_in_seq CHANGES the obs passed in.
   !       Must pass a copy of incoming obs to insert_obs_in_seq.
   !--------------------------------------------------------------
   num_good_called = 0
   num_good_searched = 0

   num_inserted = 0
   is_there_one = get_first_obs(seq_in, obs_in)

   if ( .not. is_there_one )  then
      write(msgstring2,*)  'no valid observations in ',trim(filename_seq(i))
      call error_handler(E_MSG,'obs_selection', 'skipping to next input file', &
              source, text2=msgstring2)

      call destroy_obs_sequence(seq_in)

      write(msgstring, *) 'input file ', i
      if (print_timestamps) call timestamp(string1='end   of '//trim(msgstring), pos='brief')

      cycle FILES
   endif

   if (print_timestamps) call timestamp(string1='start of first obs', pos='brief')
   ocount = 1

   ! figure out the time of the first obs in this file, and set the
   ! offset for the first item in the selection file that's at this time.
   t1 = get_time_from_obs(obs_in)
   base_index = set_base(t1, obs_def_list, obs_def_count)

   if (base_index < 0) then
      write(msgstring2, *) 'skipping all obs in ', trim(filename_seq(i))
      call error_handler(E_MSG, 'obs_selection', &
          'first time in obs_sequence file is after all times in selection file', &
           source, text2=msgstring2)
      call destroy_obs_sequence(seq_in)
      write(msgstring, *) 'input file ', i
      if (print_timestamps) call timestamp(string1='cyc 1 of '//trim(msgstring), pos='brief')
      cycle FILES      
   endif
   next_base_index = base_index

   if (good_selection(obs_in, obs_def_list, obs_def_count, next_base_index)) then
      obs_out = obs_in

      call insert_obs_in_seq(seq_out, obs_out)  ! new_obs linked list info changes

      prev_obs_out = obs_out            ! records new position in seq_out
      num_inserted = num_inserted + 1
   endif

   call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)
   ObsLoop : do while ( .not. is_this_last)

      ocount = ocount + 1

      ! before we fool with checking times and setting offsets into
      ! the selection list, skip until we are handling an obs type 
      ! that we care about.
      this_type = get_type_from_obs(next_obs_in)

      ! support identity obs - cannot exit early for them.
      if (this_type < 0 .or. type_wanted(this_type)) then

         ! if the time of the next obs is different from this one, bump
         ! up the start of the search index offset.
         t2 = get_time_from_obs(next_obs_in)
         if (t1 /= t2) then
            base_index = next_base_index
            next_base_index = set_base(t2, obs_def_list, obs_def_count, base_index)
            t1 = t2
  
            if (next_base_index < 0) then
               ! next obs in selection file is after all the rest of the obs in this
               ! input obs_seq file, so we can move on now.
               write(msgstring, *) 'done with ', trim(filename_seq(i))
               call error_handler(E_MSG, 'obs_selection', msgstring, source)
               call destroy_obs_sequence(seq_in)
               write(msgstring, *) 'input file ', i
               if (print_timestamps) call timestamp(string1='cyc 2 of '//trim(msgstring), pos='brief')
               cycle FILES      
            endif
         endif

         if (good_selection(next_obs_in, obs_def_list, obs_def_count, next_base_index)) then
            obs_out = next_obs_in

            ! Since the stride through the observation sequence file is always
            ! guaranteed to be in temporally-ascending order, we can use the
            ! 'previous' observation as the starting point to search for the
            ! correct insertion point.  This speeds up the insert code a lot.

            if (num_inserted > 0) then
               call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)
            else
               call insert_obs_in_seq(seq_out, obs_out)
            endif

            prev_obs_out = obs_out    ! update position in seq_in for next insert
            num_inserted = num_inserted + 1

            if (print_every_nth_obs > 0) then
               if (mod(num_inserted,print_every_nth_obs) == 0) then
                  write(msgstring,*) 'inserted number ',num_inserted,' of ',size_seq_in, ' possible obs'
                  call error_handler(E_MSG, 'obs_selection', msgstring, source)
               endif
            endif

         endif

      endif

      obs_in     = next_obs_in   ! next obs is the current one now

      call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

      if (print_timestamps) then
         if (mod(ocount,10000) == 0) then
            write(msgstring, *) 'processed obs ', ocount, ' of ', size_seq_in
            if (print_timestamps) call timestamp(string1=trim(msgstring), pos='brief')
         endif
      endif

   enddo ObsLoop

   total_num_inserted = total_num_inserted + num_inserted


   if (.not. print_only) then
      print*, '--------------  Obs seq file # :          ', i
      print*, 'Number of obs in previous seq  :          ', size_seq_out
      print*, 'Number of obs possible to get  :          ', size_seq_in
      print*, 'Number of obs really accepted  :          ', num_inserted
      print*, '---------------------------------------------------------'
      if (partial_write) call write_obs_seq(seq_out, filename_out)
   endif

   call destroy_obs_sequence(seq_in)

   write(msgstring, *) 'input file ', i
   if (print_timestamps) call timestamp(string1='end   of '//trim(msgstring), pos='brief')

   ! DEBUG diagnostics
   !print *, 'size of mask file = ',  obs_def_count
   !print *, 'average search length = ', num_good_searched / num_good_called
   !print *, 'number of time search called = ', num_good_called

enddo FILES

write(msgstring,*) 'Starting to process output sequence file ', trim(filename_out)
call error_handler(E_MSG,'obs_selection',msgstring)

if (.not. print_only) then
   print*, 'Total number of obs  inserted  :          ', total_num_inserted
   print*, 'Actual number of obs in the new seq file :', get_num_key_range(seq_out)
else
   print*, 'Total number of selected obs in all files :', get_num_key_range(seq_out)
endif

call print_obs_seq(seq_out, filename_out)
if (.not. print_only) then
   call write_obs_seq(seq_out, filename_out)
else
   write(msgstring,*) 'Output sequence file not created; print_only in namelist is .true.'
   call error_handler(E_MSG,'', msgstring)
endif

! Time to clean up

call destroy_obs_sequence(seq_out)
call destroy_obs(     obs_in )
call destroy_obs(next_obs_in )
call destroy_obs(     obs_out)
!call destroy_obs(prev_obs_out)

call finalize_utilities()

contains


!---------------------------------------------------------------------
subroutine obs_seq_modules_used()

! Initialize modules used that require it
call initialize_utilities('obs_selection')
call static_init_obs_sequence()

end subroutine obs_seq_modules_used

!---------------------------------------------------------------------
subroutine handle_filenames(filename_seq, filename_seq_list, num_input_files)
! sort out the input lists, set the length if not given by user,
! make sure what's specified is consistent.
character(len=*), intent(inout) :: filename_seq(:)
character(len=*), intent(in)    :: filename_seq_list
integer,          intent(inout) :: num_input_files

integer :: index
logical :: from_file
character(len=32) :: fsource

! ok, here's the new logic:
! if the user specifies neither filename_seq nor filename_seq_list, we
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
if (filename_seq(1) == '' .and. filename_seq_list == '') then

   if (num_input_files /= 0 .and. num_input_files /= 1) then
      call error_handler(E_ERR,'obs_selection', &
          'if no filenames specified, num_input_files must be 0 or 1', source)
   endif
   
   num_input_files = 1
   filename_seq(1) = 'obs_seq.out'
   return
endif

! make sure the namelist specifies one or the other but not both
if (filename_seq(1) /= '' .and. filename_seq_list /= '') then
   call error_handler(E_ERR,'obs_selection', &
       'cannot specify both filename_seq and filename_seq_list', source)
endif

! if they have specified a file which contains a list, read it into
! the filename_seq array and set the count.
if (filename_seq_list /= '') then
   fsource = 'filename_seq_list'
   from_file = .true.
else
   fsource = 'filename_seq'
   from_file = .false.
endif

do index = 1, max_num_input_files
   if (from_file) &
      filename_seq(index) = get_next_filename(filename_seq_list, index)

   if (filename_seq(index) == '') then
      if (index == 1) then
         call error_handler(E_ERR,'obs_selection', &
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

         call error_handler(E_ERR,'obs_selection', msgstring, &
            source, text2=msgstring2, text3=msgstring3)
         
      endif
   endif
enddo

write(msgstring, *) 'cannot specify more than ',max_num_input_files,' files'
call error_handler(E_ERR,'obs_selection', msgstring, source)

end subroutine handle_filenames

!---------------------------------------------------------------------
subroutine compare_metadata(seq1, seq2, fname1, fname2)

!
! This subroutine compares the metadata for two different observation
! sequences and terminates the program if they are not conformable.
! In order to be merged, the two observation sequences must have the same
! number of qc values, the same number of copies ... 
!
! FIXME:  this routine uses several globals now from the namelist that
! should be arguments to this routine.  i'm being lazy (or expedient) here
! but this should be fixed at some point soon.
!
! and now there is an additional restriction -- if editing the copies or
! qc entries, seq1 has already been edited and only 2 needs the editing
! applied.  before they were completely symmetric.

type(obs_sequence_type),    intent(in) :: seq1, seq2
character(len=*), optional, intent(in) :: fname1, fname2

integer :: num_copies1, num_qc1
integer :: num_copies2, num_qc2
integer :: num_copies , num_qc, i
character(len=metadatalength) :: str1, str2

num_copies1 = get_num_copies(seq1)
num_qc1     = get_num_qc(    seq1)

num_copies2 = get_num_copies(seq2)
num_qc2     = get_num_qc(    seq2)

num_copies  = num_copies2
num_qc      = num_qc2

! get this ready in case we have to use it below.  do not overwrite it!
if (present(fname1) .and. present(fname2)) then
   write(msgstring1,*)'Sequence files ', trim(fname1), ' and ', trim(fname2), &
                      ' are not compatible'
else
  msgstring1 = 'Sequence files cannot be merged because they are not compatible'
endif

if ( num_copies1 /= num_copies2 ) then
   write(msgstring2,*)'Different numbers of data copies found: ', &
                      num_copies1, ' vs ', num_copies2 
   call error_handler(E_ERR, 'obs_selection', msgstring1, source, text2=msgstring2)
endif
if ( num_qc1 /= num_qc2 ) then
   write(msgstring2,*)'Different different numbers of QCs found: ', &
                      num_qc1, ' vs ', num_qc2
   call error_handler(E_MSG, 'obs_selection', msgstring1, source, text2=msgstring2)
endif

! watch the code flow in this loop and the one below it.
! if a match is found, cycle.  if there is no match
! a single set of (fatal) error messages is called.
CopyMetaData : do i=1, num_copies
   str1 = get_copy_meta_data(seq1,i)
   str2 = get_copy_meta_data(seq2,i) 

   ! they match.  cycle to next copy.
   if( str1 == str2 ) then
      ! for now, don't print out if things are ok.  this could become
      ! part of a verbose option.
      !write(msgstring2,*)'metadata ',trim(str1), ' in both.'
      !call error_handler(E_MSG, 'obs_selection', msgstring2)
      cycle CopyMetaData
   endif

   ! if you get here, the metadata is not the same and the user has not
   ! given us strings that are ok to match.  fail.
   write(msgstring2,*)'metadata value mismatch. seq1: ', trim(str1)
   write(msgstring3,*)'metadata value mismatch. seq2: ', trim(str2)
   call error_handler(E_ERR, 'obs_selection', msgstring1, &
       source, text2=msgstring2, text3=msgstring3)

enddo CopyMetaData

QCMetaData : do i=1, num_qc
   str1 = get_qc_meta_data(seq1,i)
   str2 = get_qc_meta_data(seq2,i) 

   ! they match.  cycle to next copy.
   if( str1 == str2 ) then
      ! see comment in copy section above about a verbose option.
      !write(msgstring2,*)'metadata ',trim(str1), ' in both.'
      !call error_handler(E_MSG, 'obs_selection', msgstring2)
      cycle QCMetaData
   endif

   ! if you get here, the metadata is not the same and the user has not
   ! given us strings that are ok to match.  fail.
   write(msgstring2,*)'qc metadata value mismatch. seq1: ', trim(str1)
   write(msgstring3,*)'qc metadata value mismatch. seq2: ', trim(str2)
   call error_handler(E_ERR, 'obs_selection', msgstring1, &
       source, text2=msgstring2, text3=msgstring3)

enddo QCMetaData

end subroutine compare_metadata

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
! max_defined_types_of_obs is a public from obs_kind_mod.f90 and really is
! counting the max number of types, not kinds.
integer                 :: type_count(max_defined_types_of_obs), identity_count


! Initialize input obs_types
type_count(:) = 0
identity_count = 0

! num_obs should be ok since we just constructed this seq so it should
! have no unlinked obs.  if it might for some reason, use this instead:
! size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_selection',msgstring)
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
   call error_handler(E_MSG,'obs_selection', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
         call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
if (cal) call print_date(get_obs_def_time(this_obs_def), '   which is date: ')

ObsLoop : do while ( .not. is_this_last)

   this_obs_type = get_type_from_obs(obs)
   if (this_obs_type < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_type) = type_count(this_obs_type) + 1
   endif

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
               call print_time(get_obs_def_time(this_obs_def), '  Last timestamp: ')
      if (cal) call print_date(get_obs_def_time(this_obs_def), '   which is date: ')
   endif

enddo ObsLoop


write(msgstring, *) 'Number of obs processed  :          ', size_seq_in
call error_handler(E_MSG, '', msgstring)
write(msgstring, *) '---------------------------------------------------------'
call error_handler(E_MSG, '', msgstring)
do i = 1, max_defined_types_of_obs
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
logical                 :: is_there_one, is_this_last
integer                 :: size_seq
integer                 :: key
type(time_type)         :: last_time, this_time


! make sure there are obs left to process before going on.
size_seq = get_num_obs(seq) 
if (size_seq == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_selection:validate',msgstring)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(next_obs, get_num_copies(seq), get_num_qc(seq))

!-------------------------------------------------------------
! Start to process obs from seq
!--------------------------------------------------------------
is_there_one = get_first_obs(seq, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in sequence ' // trim(filename)
   call error_handler(E_MSG,'obs_selection:validate', msgstring, source)
endif

last_time = get_time_from_obs(obs)

is_this_last = .false.
ObsLoop : do while ( .not. is_this_last)

   this_time = get_time_from_obs(obs)

   if (last_time > this_time) then
      ! bad time order of observations in linked list
               call print_time(last_time, ' previous timestamp: ')
      if (cal) call print_date(last_time, '      which is date: ')
               call print_time(this_time, '     next timestamp: ')
      if (cal) call print_date(this_time, '      which is date: ')

      key = get_obs_key(obs)
      write(msgstring2,*)'obs number ', key, ' has earlier time than previous obs'
      write(msgstring,*)'observations must be in increasing time order, file ' // trim(filename)
      call error_handler(E_ERR,'obs_selection:validate', msgstring, &
                source, text2=msgstring2)
   endif

   last_time = this_time

   call get_next_obs(seq, obs, next_obs, is_this_last)
   if (.not. is_this_last) obs = next_obs

enddo ObsLoop

! clean up
call destroy_obs(     obs)
call destroy_obs(next_obs)

end subroutine validate_obs_seq_time

!---------------------------------------------------------------------
subroutine print_metadata(seq1, fname1)

!
! print out the metadata strings, trimmed
!

type(obs_sequence_type),    intent(in) :: seq1
character(len=*), optional, intent(in) :: fname1

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str1

num_copies = get_num_copies(seq1)
num_qc     = get_num_qc(    seq1)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring1,*)' illegal copy or obs count in file '//trim(fname1)
   call error_handler(E_ERR, 'obs_selection', msgstring1, source)
endif

MetaDataLoop : do i=1, num_copies
   str1 = get_copy_meta_data(seq1,i)

   write(msgstring1,*)'Data Metadata: ',trim(str1)
   call error_handler(E_MSG, '', msgstring1)

enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str1 = get_qc_meta_data(seq1,i)

   write(msgstring1,*)'  QC Metadata: ', trim(str1)
   call error_handler(E_MSG, '', msgstring1)

enddo QCMetaData

end subroutine print_metadata

!---------------------------------------------------------------------
subroutine read_selection_list(select_file, select_is_seq, &
                               selection_list, selection_count, type_wanted)
 character(len=*),                intent(in)  :: select_file
 logical,                         intent(in)  :: select_is_seq
 type(obs_def_type), allocatable, intent(out) :: selection_list(:)
 integer,                         intent(out) :: selection_count
 logical,            allocatable, intent(out) :: type_wanted(:)

 ! the plan:
 ! open file
 ! read count
 ! call read_obs_def right number of times
 ! close file

 integer :: iunit, count, i, copies, qcs, this_type
 character(len=15) :: label    ! must be 'num_definitions'
 type(obs_type) :: obs, prev_obs
 type(obs_sequence_type) :: seq_in
 logical :: is_this_last
 real(r8) :: dummy
 type(obs_def_type), allocatable :: temp_sel_list(:)
 type(time_type), allocatable :: temp_time(:)
 integer, allocatable :: sort_index(:)

 ! if the list of which obs to select comes from the coverage tool,
 ! it's a list of obs_defs.   if it's a full obs_seq file, then
 ! use the normal tools and just pull out the list of obs_defs.
 if (.not. select_is_seq) then
     iunit = open_file(select_file, form='formatted', action='read')
    
     read(iunit, *) label, count
     if (label /= 'num_definitions') then
         call error_handler(E_ERR,'obs_selection', &
           'bad format, expected "num_definitions" at start of selection file', source)
     endif
    
     ! set up the mapping table for the kinds here
     call read_type_of_obs_table(iunit, .false.)
    
     ! these ones stay around and are returned from this subroutine
     allocate(selection_list(count), type_wanted(max_defined_types_of_obs))
  
     ! these are temporaries for sorting into time order
     allocate(temp_sel_list(count), temp_time(count), sort_index(count))
    
     ! assume no types are wanted
     type_wanted(:) = .false.

     ! read into array in whatever order these are in the file
     ! and bookkeep what types are encountered
     do i = 1, count
         call read_obs_def(iunit, temp_sel_list(i), 0, dummy)
         this_type = get_obs_def_type_of_obs(temp_sel_list(i))
         if (this_type > 0) type_wanted(this_type) = .true.
     enddo
    
     ! extract just the times into an array
     do i = 1, count
         temp_time(i) = get_obs_def_time(temp_sel_list(i))
     enddo
   
     ! sort into time order
     call time_index_sort(temp_time, sort_index, count)

     ! copy into output array using that order
     do i = 1, count
         selection_list(i) = temp_sel_list(sort_index(i))
     enddo


     deallocate(temp_sel_list, temp_time, sort_index)

     call close_file(iunit)
 else
    
     call read_obs_seq(select_file, 0, 0, 0, seq_in)

     count  = get_num_obs(seq_in)
     copies = get_num_copies(seq_in)
     qcs    = get_num_qc(seq_in)

     call init_obs(obs,      copies, qcs)
     call init_obs(prev_obs, copies, qcs)

     allocate(selection_list(count), type_wanted(max_defined_types_of_obs))
    
     if (.not. get_first_obs(seq_in, obs)) then
         call error_handler(E_ERR,'obs_selection', &
                 'empty obs_seq for selection file', source)
     endif

     ! in an obs_seq file, using get_next_obs you will
     ! be guarenteed to get these in increasing time order.

     is_this_last = .false.
     do i = 1, count
         if (is_this_last) exit 

         call get_obs_def(obs, selection_list(i))
         this_type = get_obs_def_type_of_obs(selection_list(i))
         if (this_type > 0) type_wanted(this_type) = .true.

         prev_obs = obs
         call get_next_obs(seq_in, prev_obs, obs, is_this_last)
     enddo
    
    call destroy_obs(obs)
    call destroy_obs_sequence(seq_in)
 endif

 selection_count = count

write(msgstring, *) 'selection file contains ', count, ' entries'
call error_handler(E_MSG, 'obs_selection', msgstring, source)
         call print_time(get_obs_def_time(selection_list(1)),     'time of first selection:')
if (cal) call print_date(get_obs_def_time(selection_list(1)),     '          which is date:')
         call print_time(get_obs_def_time(selection_list(count)), 'time of last  selection:')
if (cal) call print_date(get_obs_def_time(selection_list(count)), '          which is date:')

end subroutine read_selection_list

!---------------------------------------------------------------------
function set_base(obs_time, selection_list, selection_count, startindex)

! find offset of first obs_def entry that has a time >= to the given one
! if optional startindex is given, start looking there.

 type(time_type),    intent(in) :: obs_time
 type(obs_def_type), intent(in) :: selection_list(:)
 integer,            intent(in) :: selection_count
 integer, optional,  intent(in) :: startindex
 integer :: set_base

 integer :: i, s
 type(time_type) :: def_time

 ! if we know an offset to start from, use it.  otherwise, start at 1.
 if (present(startindex)) then
    s = startindex
 else
    s = 1
 endif
 
 ! find the index of the first item in this list which is >= given one
 ! we will start subsequent searches at this offset.
 do i = s, selection_count

    def_time = get_obs_def_time(selection_list(i))

    ! if we have looped through the obs_def list far enough so
    ! the next obs_def entry has a time more than or equal to the
    ! one of the next observation, stop and return this index.
    if (obs_time <= def_time) then
       set_base = i
       return
    endif

 enddo

! we got all the way through the file and the last obs_def entry
! was still before the observation time we are looking for.
set_base = -1

end function set_base



! compare horiz location, time, type - ignores vertical
!---------------------------------------------------------------------
function good_selection(obs_in, selection_list, selection_count, startindex)
 type(obs_type),     intent(in) :: obs_in
 type(obs_def_type), intent(in) :: selection_list(:)
 integer,            intent(in) :: selection_count
 integer,            intent(in) :: startindex
 logical :: good_selection

 ! first pass, iterate list.
 ! if too slow, pull in the get_close code


 integer :: i
 type(obs_def_type)  :: base_obs_def
 integer             :: base_obs_type, test_obs_type
 type(time_type)     :: base_obs_time, test_obs_time
 type(location_type) :: base_obs_loc,  test_obs_loc

 call get_obs_def(obs_in, base_obs_def)
 base_obs_loc  = get_obs_def_location(base_obs_def)
 base_obs_time = get_obs_def_time(base_obs_def)
 base_obs_type = get_obs_def_type_of_obs(base_obs_def)

 ! this program now time-sorts the selection list first, so we
 ! are guarenteed the obs_defs will be encountered in time order.
 ! now optimize in two ways:  first, the caller will pass in an
 ! offset that is less than or equal to the first obs for this time,
 ! and second this routine will return when the obs_def time is
 ! larger than the time to match.  that should save a lot of looping.

 if (startindex < 1 .or. startindex > selection_count) then
    write(msgstring2, *) 'startindex ', startindex, ' is not between 1 and ', selection_count
    call error_handler(E_ERR, 'good_selection', &
       'invalid startindex, internal error should not happen', &
       source, text2=msgstring2) 
 endif
 
 ! statistics for timing/debugging
 num_good_called = num_good_called + 1

 ! the obs_def list is time-sorted now.  if you get to a larger
 ! time without a match you can bail early.
 good_selection = .false.

 ! first select on time - it is an integer comparison
 ! and quicker than location test.  then type, and
 ! finally location.
 do i = startindex, selection_count

    test_obs_time = get_obs_def_time(selection_list(i))

    ! if past the possible times, return now.
    if (base_obs_time > test_obs_time) then
       num_good_searched = i - startindex + 1
       return
    endif

    if (base_obs_time /= test_obs_time) cycle

    test_obs_type = get_obs_def_type_of_obs(selection_list(i))
    if (base_obs_type /= test_obs_type) cycle

    test_obs_loc = get_obs_def_location(selection_list(i))
    if ( .not. same_location(base_obs_loc, test_obs_loc, match_vertical)) cycle

    ! all match - good return.
    num_good_searched = i - startindex + 1
    good_selection = .true.
    return
 enddo

 ! if you get here, no match, and return value is already false.
num_good_searched = selection_count - startindex + 1

end function good_selection

!---------------------------------------------------------------------
subroutine destroy_selections(selection_list, selection_count)
 type(obs_def_type), allocatable, intent(out) :: selection_list(:)
 integer,                         intent(out) :: selection_count

  ! clean up

  deallocate(selection_list, type_wanted)
  selection_count = 0

end subroutine destroy_selections

!---------------------------------------------------------------------------
function same_location(loc1,loc2,vert_flag)

! compares two locations.
! the vert_flag specifies whether to consider only the horizontal values
! of lat/lon, or whether to include the vertical coordinate.
!
! if including the vertical, the two coordinate systems must be
! the same (pressure, height, level, surface, scaleheight, undef).
! the tolerances for saying the two locations are the same is a namelist
! item setting.  horizontal values are in degrees.  if that tolerance 
! is specified as <= 0, allow the lowest bit to differ from roundoff error.
! if vert_flag is true, vert must also match within tolerance to return same.
! if vert_flag is false, just use the horizontal.

type(location_type), intent(in) :: loc1, loc2
logical, intent(in)             :: vert_flag
logical                         :: same_location

real(r8) :: larray1(3), larray2(3)
integer  :: whichv1, whichv2

larray1 = get_location(loc1)
larray2 = get_location(loc2)

! strategy is to assume the locations don't match, and when
! you know this you can immediately return.  if you get all
! the way to the end of this routine, then they are the same.

same_location = .false.

if (latlon_tolerance <= 0.0_r8) then
   if ( abs(larray1(1) - larray2(1)) > epsilon(larray1(1)) ) return
   if ( abs(larray1(2) - larray2(2)) > epsilon(larray1(2)) ) return
else
   if ( abs(larray1(1) - larray2(1)) > latlon_tolerance ) return
   if ( abs(larray1(2) - larray2(2)) > latlon_tolerance ) return
endif

if (vert_flag) then
   whichv1 = nint(query_location(loc1))
   whichv2 = nint(query_location(loc2))
   if (whichv1 /= whichv2) return

   select case(whichv1) 
      case (VERTISUNDEF)
         continue
      case (VERTISSURFACE)
         if ( abs(larray1(3) - larray2(3)) > surface_tolerance ) return
      case (VERTISLEVEL)
         if ( abs(larray1(3) - larray2(3)) > level_tolerance ) return
      case (VERTISPRESSURE)
         if ( abs(larray1(3) - larray2(3)) > pressure_tolerance ) return
      case (VERTISHEIGHT)
         if ( abs(larray1(3) - larray2(3)) > height_tolerance ) return
      case (VERTISSCALEHEIGHT)
         if ( abs(larray1(3) - larray2(3)) > scaleheight_tolerance ) return
      case default
         write(msgstring, *) 'unrecognized key for vertical type: ', whichv1
         call error_handler(E_ERR, 'same_location', msgstring, source)
   end select
 
endif

same_location = .true.

end function same_location

!---------------------------------------------------------------------------
function get_time_from_obs(this_obs)
  type(obs_type), intent(in) :: this_obs
  type(time_type) :: get_time_from_obs

type(obs_def_type) :: this_obs_def

call get_obs_def(this_obs, this_obs_def)
get_time_from_obs = get_obs_def_time(this_obs_def)

end function get_time_from_obs

!---------------------------------------------------------------------------
function get_type_from_obs(this_obs)
  type(obs_type), intent(in) :: this_obs
  integer :: get_type_from_obs

type(obs_def_type) :: this_obs_def

call get_obs_def(this_obs, this_obs_def)
get_type_from_obs = get_obs_def_type_of_obs(this_obs_def)

end function get_type_from_obs

!---------------------------------------------------------------------
end program obs_selection

