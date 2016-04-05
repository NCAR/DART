! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program obs_selection

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! this latest addition has select by list of obs types.

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : timestamp, register_module, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit,   &
                             do_nml_file, do_nml_term, get_next_filename, &
                             open_file, close_file
use     location_mod, only : location_type, get_location, set_location, &
                             LocationName, read_location, operator(==), &
                             write_location
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind, &
                             get_obs_def_location, read_obs_def
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name, get_obs_kind_index, &
                             read_obs_kind
use time_manager_mod, only : time_type, operator(>), print_time, set_time, &
                             print_date, set_calendar_type, GREGORIAN,     &
                             operator(/=), NO_CALENDAR, get_calendar_type
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

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs_in, next_obs_in
type(obs_type)          :: obs_out, prev_obs_out
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in, num_copies_in, num_qc_in
integer                 :: size_seq_out, num_copies_out, num_qc_out
integer                 :: num_inserted, iunit, io, i, j, total_num_inserted
integer                 :: max_num_obs, file_id
integer                 :: first_seq
character(len = 129)    :: read_format, meta_data
logical                 :: pre_I_format, cal
character(len = 129)    :: msgstring

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 20

!----------------------------------------------------------------
! Namelist input with default values

! max_num_input_files : maximum number of input sequence files to be processed
! lazy, pick big numbers.  make them bigger if too small.
integer, parameter               :: max_num_input_files = 1000
integer, parameter               :: max_obs_input_types = 500
integer                          :: num_input_files = 0
type(obs_def_type),  allocatable :: obs_def_list(:)
integer                          :: obs_def_count


character(len = 129) :: filename_seq(max_num_input_files) = ''
character(len = 129) :: filename_seq_list = ''
character(len = 129) :: filename_out  = 'obs_seq.processed'
logical              :: process_file(max_num_input_files)

character(len = 129) :: selections_file = 'obsdef_mask.txt'

logical  :: selections_is_obs_seq = .false.
logical  :: print_only            = .false.
character(len=32) :: calendar     = 'Gregorian'


namelist /obs_selection_nml/ &
         num_input_files, filename_seq, filename_seq_list, filename_out, &
         selections_file, selections_is_obs_seq, print_only, calendar

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

call read_selection_list(selections_file, selections_is_obs_seq, obs_def_list, obs_def_count)

! end of namelist processing and setup

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

first_seq = -1
do i = 1, num_input_files

   if ( len(filename_seq(i)) .eq. 0 .or. filename_seq(i) .eq. "" ) then
      call error_handler(E_ERR,'obs_selection', &
         'num_input_files and filename_seq mismatch',source,revision,revdate)
   endif

   ! count up the number of observations we are going to eventually have.
   ! if all the observations in a file are not part of the linked list, the
   ! output number of observations might be much smaller than the total size in
   ! the header.  it is slower, but go ahead and read in the entire sequence
   ! and count up the real number of obs - trim_seq will do the count even if
   ! it is not trimming in time.  this allows us to create an empty obs_seq
   ! output file of exactly the right size.

   call read_obs_seq_header(filename_seq(i), num_copies_in, num_qc_in, &
      size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)

   call destroy_obs_sequence(seq_in)
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

! no valid obs found?  if the index value is still negative, we are
! still waiting to process the first one and never found one.
if (first_seq < 0 .or. size_seq_out == 0) then
   msgstring = 'All input files are empty'
   call error_handler(E_MSG,'obs_by_selection',msgstring,source,revision,revdate)
endif


! pass 2:

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
do i = 1, num_input_files

   if (.not. process_file(i)) cycle

   write(msgstring,*) 'Starting to process input sequence file ', trim(filename_seq(i))
   call error_handler(E_MSG,'obs_selection',msgstring)

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
   num_inserted = 0
   is_there_one = get_first_obs(seq_in, obs_in)

   if ( is_there_one )  then

      if (good_selection(obs_in, obs_def_list, obs_def_count)) then
         obs_out = obs_in

         call insert_obs_in_seq(seq_out, obs_out)  ! new_obs linked list info changes

         prev_obs_out = obs_out            ! records new position in seq_out
         num_inserted = num_inserted + 1
      endif

      call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)
      ObsLoop : do while ( .not. is_this_last)

         obs_in     = next_obs_in   ! essentially records position in seq_out

         if (good_selection(obs_in, obs_def_list, obs_def_count)) then
            obs_out = obs_in

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

            if (print_every > 0) then
               if (mod(num_inserted,print_every) == 0) then
                  print*, 'inserted number ',num_inserted,' of ',size_seq_in
               endif
            endif

         endif

         call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

      enddo ObsLoop

      total_num_inserted = total_num_inserted + num_inserted

   else
      write(msgstring,*)'no first observation in ',trim(filename_seq(i))
      call error_handler(E_MSG,'obs_selection', msgstring)
   endif

   if (.not. print_only) then
      print*, '--------------  Obs seq file # :          ', i
      print*, 'Number of obs in previous seq  :          ', size_seq_out
      print*, 'Number of obs to be  inserted  :          ', size_seq_in
      print*, 'Number of obs really inserted  :          ', num_inserted
      print*, '---------------------------------------------------------'
   endif

   call destroy_obs_sequence(seq_in)

enddo

write(msgstring,*) 'Starting to process output sequence file ', trim(filename_out)
call error_handler(E_MSG,'obs_sequence_tool',msgstring)

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

call timestamp(source,revision,revdate,'end')

contains


!---------------------------------------------------------------------
subroutine obs_seq_modules_used()

! Initialize modules used that require it
call initialize_utilities('obs_sequence_tool')
call register_module(source,revision,revdate)
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
character(len=32) :: source

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
      call error_handler(E_ERR,'obs_sequence_tool', &
          'if no filenames specified, num_input_files must be 0 or 1', &
          source,revision,revdate)
   endif
   
   num_input_files = 1
   filename_seq(1) = 'obs_seq.out'
   return
endif

! make sure the namelist specifies one or the other but not both
if (filename_seq(1) /= '' .and. filename_seq_list /= '') then
   call error_handler(E_ERR,'obs_sequence_tool', &
       'cannot specify both filename_seq and filename_seq_list', &
       source,revision,revdate)
endif

! if they have specified a file which contains a list, read it into
! the filename_seq array and set the count.
if (filename_seq_list /= '') then
   source = 'filename_seq_list'
   from_file = .true.
else
   source = 'filename_seq'
   from_file = .false.
endif

do index = 1, max_num_input_files
   if (from_file) &
      filename_seq(index) = get_next_filename(filename_seq_list, index)

   if (filename_seq(index) == '') then
      if (index == 1) then
         call error_handler(E_ERR,'obs_sequence_tool', &
             trim(source)//' contains no filenames', &
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

         write(msgstring, *) 'if num_input_files is 0, the number of files will be automatically computed'
         call error_handler(E_MSG,'obs_sequence_tool', msgstring)
         write(msgstring, *) 'if num_input_files is not 0, it must match the number of filenames specified'
         call error_handler(E_MSG,'obs_sequence_tool', msgstring)
         write(msgstring, *) 'num_input_files is ', num_input_files, &
                     ' but '//trim(source)//' has filecount ', index - 1
         call error_handler(E_ERR,'obs_sequence_tool', msgstring, &
            source,revision,revdate)
         
      endif
   endif
enddo

write(msgstring, *) 'cannot specify more than ',max_num_input_files,' files'
call error_handler(E_ERR,'obs_sequence_tool', msgstring, &
     source,revision,revdate)

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

 type(obs_sequence_type), intent(IN) :: seq1, seq2
 character(len=*), optional :: fname1, fname2

integer :: num_copies1, num_qc1
integer :: num_copies2, num_qc2
integer :: num_copies , num_qc, i
character(len=129) :: str1, str2
character(len=255) :: msgstring1, msgstring2

num_copies1 = get_num_copies(seq1)
num_qc1     = get_num_qc(    seq1)

num_copies2 = get_num_copies(seq2)
num_qc2     = get_num_qc(    seq2)

num_copies  = num_copies2
num_qc      = num_qc2

! get this ready in case we have to use it
if (present(fname1) .and. present(fname2)) then
   write(msgstring1,*)'Sequence files ', trim(fname1), ' and ', trim(fname2), &
                      ' are not compatible'
else
  msgstring1 = 'Sequence files cannot be merged because they are not compatible'
endif

if ( num_copies1 /= num_copies2 ) then
   write(msgstring2,*)'Different numbers of data copies found: ', &
                      num_copies1, ' vs ', num_copies2 
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   num_copies = -1
endif
if ( num_qc1 /= num_qc2 ) then
   write(msgstring2,*)'Different different numbers of QCs found: ', &
                      num_qc1, ' vs ', num_qc2
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   num_qc = -1
endif
if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)
endif

! watch the code flow in this loop and the one below it.
! the smoothest code order is to determine what the strings are,
! and then try different things to match them.  if a match is found,
! cycle.  if you get to the bottom of the loop, there is no match
! and a single set of (fatal) error messages is called.
CopyMetaData : do i=1, num_copies
   str1 = get_copy_meta_data(seq1,i)
   str2 = get_copy_meta_data(seq2,i) 

   ! easy case - they match.  cycle to next copy.
   if( str1 == str2 ) then
      write(msgstring2,*)'metadata ',trim(str1), ' in both.'
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
      cycle CopyMetaData
   endif

   ! if you get here, the metadata is not the same and the user has not
   ! given us strings that are ok to match.  fail.
   write(msgstring2,*)'metadata value mismatch. seq1: ', trim(str1)
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   write(msgstring2,*)'metadata value mismatch. seq2: ', trim(str2)
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)

enddo CopyMetaData

QCMetaData : do i=1, num_qc
   str1 = get_qc_meta_data(seq1,i)
   str2 = get_qc_meta_data(seq2,i) 

   ! easy case - they match.  cycle to next copy.
   if( str1 == str2 ) then
      write(msgstring2,*)'metadata ',trim(str1), ' in both.'
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
      cycle QCMetaData
   endif

   ! if you get here, the metadata is not the same and the user has not
   ! given us strings that are ok to match.  fail.
   write(msgstring2,*)'qc metadata value mismatch. seq1: ', trim(str1)
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   write(msgstring2,*)'qc metadata value mismatch. seq2: ', trim(str2)
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)

enddo QCMetaData

end subroutine compare_metadata

!---------------------------------------------------------------------
subroutine print_obs_seq(seq_in, filename)

! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

type(obs_sequence_type), intent(in) :: seq_in
character(len=*), intent(in)        :: filename

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
   call error_handler(E_MSG,'obs_sequence_tool',msgstring)
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
   call error_handler(E_MSG,'obs_sequence_tool', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
         call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
if (cal) call print_date(get_obs_def_time(this_obs_def), '   which is date: ')

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
      if (cal) call print_date(get_obs_def_time(this_obs_def), '   which is date: ')
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
integer                 :: size_seq
integer                 :: key
type(time_type)         :: last_time, this_time


! make sure there are obs left to process before going on.
size_seq = get_num_obs(seq) 
if (size_seq == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_sequence_tool:validate',msgstring)
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
   call error_handler(E_MSG,'obs_sequence_tool:validate', msgstring, source, revision, revdate)
endif

call get_obs_def(obs, this_obs_def)
last_time = get_obs_def_time(this_obs_def)


is_this_last = .false.
ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_time = get_obs_def_time(this_obs_def)

   if (last_time > this_time) then
      ! bad time order of observations in linked list
               call print_time(last_time, ' previous timestamp: ')
      if (cal) call print_date(last_time, '      which is date: ')
               call print_time(this_time, '     next timestamp: ')
      if (cal) call print_date(this_time, '      which is date: ')

      key = get_obs_key(obs)
      write(msgstring,*)'obs number ', key, ' has earlier time than previous obs'
      call error_handler(E_MSG,'obs_sequence_tool:validate', msgstring)
      write(msgstring,*)'observations must be in increasing time order, file ' // trim(filename)
      call error_handler(E_ERR,'obs_sequence_tool:validate', msgstring, source, revision, revdate)
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

type(obs_sequence_type), intent(in) :: seq1
character(len=*), optional :: fname1

integer :: num_copies , num_qc, i
character(len=129) :: str1
character(len=255) :: msgstring1

num_copies = get_num_copies(seq1)
num_qc     = get_num_qc(    seq1)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring1,*)' illegal copy or obs count in file '//trim(fname1)
   call error_handler(E_ERR, 'obs_selection', msgstring1, source, revision, revdate)
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
                               selection_list, selection_count)
 character(len=*),                intent(in)  :: select_file
 logical,                         intent(in)  :: select_is_seq
 type(obs_def_type), allocatable, intent(out) :: selection_list(:)
 integer,                         intent(out) :: selection_count

 ! the plan:
 ! open file
 ! read count
 ! call read_obs_def right number of times
 ! close file

 integer :: iunit, count, i, copies, qcs
 character(len=15) :: label    ! must be 'num_definitions'
 type(obs_type) :: obs, prev_obs
 type(obs_sequence_type) :: seq_in
 logical :: is_this_last
 real(r8) :: dummy

 ! if the list of which obs to select comes from the coverage tool,
 ! it's a list of obs_defs.   if it's a full obs_seq file, then
 ! use the normal tools and just pull out the list of obs_defs.
 if (.not. select_is_seq) then
     iunit = open_file(select_file, form='formatted', action='read')
    
     read(iunit, *) label, count
     if (label /= 'num_definitions') then
         call error_handler(E_ERR,'obs_selection', &
           'bad format, expected "num_definitions" at start of selection file', &
            source,revision,revdate)
     endif
    
     allocate(selection_list(count))
    
     ! set up the mapping table for the kinds here
     call read_obs_kind(iunit, .false.)
    
     do i = 1, count
         call read_obs_def(iunit, selection_list(i), 0, dummy)
     enddo
    
     call close_file(iunit)
 else
    
     call read_obs_seq(select_file, 0, 0, 0, seq_in)

     count  = get_num_obs(seq_in)
     copies = get_num_copies(seq_in)
     qcs    = get_num_qc(seq_in)

     call init_obs(obs,      copies, qcs)
     call init_obs(prev_obs, copies, qcs)

     allocate(selection_list(count))
    
     if (.not. get_first_obs(seq_in, obs)) then
         call error_handler(E_ERR,'obs_selection', &
           'empty obs_seq for selection file', &
            source,revision,revdate)
     endif

     is_this_last = .false.
     do i = 1, count
         if (is_this_last) exit 

         call get_obs_def(obs, selection_list(i))

         prev_obs = obs
         call get_next_obs(seq_in, prev_obs, obs, is_this_last)
     enddo
    
    call destroy_obs(obs)
    call destroy_obs_sequence(seq_in)
 endif

 selection_count = count

end subroutine read_selection_list


! compare horiz location, time, type - ignores vertical
!---------------------------------------------------------------------
function good_selection(obs_in, selection_list, selection_count)
 type(obs_type),     intent(in) :: obs_in
 type(obs_def_type), intent(in) :: selection_list(:)
 integer,            intent(in) :: selection_count
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
 base_obs_type = get_obs_kind(base_obs_def)

 ! first select on time - it is an integer comparison
 ! and quicker than location test.  then type, and
 ! finally location.
 do i = 1, selection_count

    test_obs_type = get_obs_kind(selection_list(i))
    if (base_obs_type /= test_obs_type) cycle

    test_obs_time = get_obs_def_time(selection_list(i))
    if (base_obs_time /= test_obs_time) cycle

    test_obs_loc = get_obs_def_location(selection_list(i))
    if ( .not. horiz_location_equal(base_obs_loc, test_obs_loc)) cycle

    ! all match - good return.
    good_selection = .true.
    return
 enddo

 ! got to end of selection list without a match, cycle.
 good_selection = .false.

end function good_selection

!---------------------------------------------------------------------
subroutine destroy_selections(selection_list, selection_count)
 type(obs_def_type), allocatable, intent(out) :: selection_list(:)
 integer,                         intent(out) :: selection_count

  ! clean up

  deallocate(selection_list)
  selection_count = 0

end subroutine destroy_selections

!---------------------------------------------------------------------------
function horiz_location_equal(loc1,loc2)

! function to compare only the lat & lon and ignore the vert location.

type(location_type), intent(in) :: loc1, loc2
logical                         :: horiz_location_equal

real(r8) :: l1(3), l2(3)

l1 = get_location(loc1)
l2 = get_location(loc2)

horiz_location_equal = .false.

if ( abs(l1(1)  - l2(1) ) > epsilon(l1(1) ) ) return
if ( abs(l1(2)  - l2(2) ) > epsilon(l1(2) ) ) return

horiz_location_equal = .true.

end function horiz_location_equal

!---------------------------------------------------------------------
end program obs_selection

