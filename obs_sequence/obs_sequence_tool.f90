! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! this latest addition has select by list of obs types.

program obs_sequence_tool

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision
! $Date

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : timestamp, register_module, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit,   &
                             do_nml_file, do_nml_term, get_next_filename
use     location_mod, only : location_type, get_location, set_location2, &
                             LocationName !! , vert_is_height 
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind, &
                             get_obs_def_location
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name, get_obs_kind_index
use time_manager_mod, only : time_type, operator(>), print_time, set_time, &
                             print_date, set_calendar_type, GREGORIAN
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
                             get_num_key_range, delete_obs_by_typelist, &
                             get_obs_key, copy_partial_obs, &
                             delete_obs_from_seq, get_next_obs_from_key, &
                             delete_obs_by_qc, delete_obs_by_copy, &
                             select_obs_by_location 
implicit none

! <next few lines under version control, do not edit>
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
integer                 :: max_num_obs, file_id, remaining_obs_count
integer                 :: first_seq
character(len = 129)    :: read_format, meta_data
logical                 :: pre_I_format, all_gone
logical                 :: trim_first, trim_last
character(len = 129)    :: msgstring

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 20000

!----------------------------------------------------------------
! Namelist input with default values

! max_num_input_files : maximum number of input sequence files to be processed
integer, parameter :: max_num_input_files = 500
integer :: num_input_files = 0

! lazy, pick a big number
integer, parameter :: max_obs_input_types = 500
character(len = obstypelength) :: obs_types(max_obs_input_types) = ''
logical :: restrict_by_obs_type
integer :: num_obs_input_types
logical :: restrict_by_location
logical :: restrict_by_qc, restrict_by_copy, restrict_by_height
logical :: editing_obs, matching_copy_metadata, matching_qc_metadata
integer :: matching_copy_limit = 0
integer :: matching_qc_limit   = 0


character(len = 129) :: filename_seq(max_num_input_files) = ''
character(len = 129) :: filename_seq_list = ''
character(len = 129) :: filename_out  = 'obs_seq.processed'
logical              :: process_file(max_num_input_files)

! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days    = -1
integer  :: first_obs_seconds = -1
integer  :: last_obs_days     = -1
integer  :: last_obs_seconds  = -1
! must be location independent...
type(location_type) :: min_loc, max_loc
real(r8) :: min_box(4) = missing_r8
real(r8) :: max_box(4) = missing_r8
real(r8) :: min_lat = -90.0_r8
real(r8) :: max_lat =  90.0_r8
real(r8) :: min_lon =   0.0_r8
real(r8) :: max_lon = 360.0_r8
type(time_type) :: first_obs_time, last_obs_time
real(r8) :: min_qc = missing_r8
real(r8) :: max_qc = missing_r8
real(r8) :: min_copy = missing_r8
real(r8) :: max_copy = missing_r8
character(len = obstypelength) :: copy_type = ''
character(len=metadatalength)  :: copy_metadata = ''
character(len=metadatalength)  :: qc_metadata = ''
! 256 is an arb max of number of copies for data and qc
logical  :: edit_copy_metadata = .false.
character(len=metadatalength) :: new_copy_metadata(256) = ''
logical  :: edit_copies = .false.
integer  :: new_copy_index(256) = -1
integer  :: copy_index_len = 0
logical  :: edit_qc_metadata = .false.
character(len=metadatalength) :: new_qc_metadata(256)   = ''
logical  :: edit_qcs = .false.
integer  :: new_qc_index(256) = -1
integer  :: qc_index_len = 0
character(len=metadatalength) :: synonymous_copy_list(256) = ''
character(len=metadatalength) :: synonymous_qc_list(256)   = ''
logical  :: keep_types = .true.
logical  :: print_only = .false.
logical  :: gregorian_cal = .true.
real(r8) :: min_gps_height = missing_r8


namelist /obs_sequence_tool_nml/ &
         num_input_files, filename_seq, filename_seq_list, filename_out,     &
         first_obs_days, first_obs_seconds, last_obs_days, last_obs_seconds, &
         obs_types, keep_types, min_box, max_box, print_only,                &
         min_lat, max_lat, min_lon, max_lon, min_qc, max_qc, qc_metadata,    &
         min_copy, max_copy, copy_metadata, copy_type, gregorian_cal,        &
         min_gps_height, edit_copy_metadata, new_copy_index,                 &
         edit_qc_metadata, new_qc_index, synonymous_copy_list,               &
         synonymous_qc_list, edit_copies, edit_qcs, new_copy_metadata,       &
         new_qc_metadata

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call obs_seq_modules_used()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_sequence_tool_nml", iunit)
read(iunit, nml = obs_sequence_tool_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_sequence_tool_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_sequence_tool_nml)
if (do_nml_term()) write(     *     , nml=obs_sequence_tool_nml)

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


! if you are not using a gregorian cal, set this to false in the namelist.
! if users need it, we could add a calendar type integer to the namelist,
! if users want to specify a particular calendar which is not gregorian.
! (earlier versions of this file had the test before the namelist read - duh.)
if (gregorian_cal) then
   call set_calendar_type(GREGORIAN)
endif

! See if the user is restricting the obs types to be processed, and set up
! the values if so.
num_obs_input_types = 0
do i = 1, max_obs_input_types
   if ( len(obs_types(i)) .eq. 0 .or. obs_types(i) .eq. "" ) exit
   num_obs_input_types = i
enddo
if (num_obs_input_types == 0) then
   restrict_by_obs_type = .false.
else
   restrict_by_obs_type = .true.
endif

! See if the user is restricting the obs locations to be processed, and set up
! the values if so.  Note that if any of the values are not missing_r8, all of
! them must have good (and consistent) values.  i don't have code in here yet
! to idiotproof this, but i should add it.
if ((minval(min_box).ne.missing_r8) .or. (maxval(min_box).ne.missing_r8) .or. &
    (minval(max_box).ne.missing_r8) .or. (maxval(max_box).ne.missing_r8)) then
   restrict_by_location = .true.
   min_loc = set_location2(min_box)
   max_loc = set_location2(max_box)
else
   restrict_by_location = .false.
endif

! 3d sphere box check - locations module dependent, but an important one.
if ((min_lat /= -90.0_r8) .or. (max_lat /=  90.0_r8) .or. &
    (min_lon /=   0.0_r8) .or. (max_lon /= 360.0_r8)) then
   ! do not allow namelist to set BOTH min/max box and by lat/lon.
   if (restrict_by_location) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'use either lat/lon box or min/max box but not both', &
                         source,revision,revdate)
   endif
   restrict_by_location = .true.
   if (trim(LocationName) /= 'loc3Dsphere') then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'can only use lat/lon box with 3d sphere locations', &
                         source,revision,revdate)
   endif
   ! simple err checks before going on; try to catch radians vs degrees or
   ! just plain junk or typos.
   if (min_lat >= max_lat) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lat must be less than max_lat', &
                         source,revision,revdate)
   endif
   if (min_lat < -90.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lat cannot be less than -90.0 degrees', &
                         source,revision,revdate)
   endif
   if (max_lat >  90.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'max_lat cannot be greater than  90.0 degrees', &
                         source,revision,revdate)
   endif
   ! this is ok - e.g.-180 to 180
   !if (min_lon < 0.0_r8) then
   !   call error_handler(E_ERR,'obs_sequence_tool', &
   !                      'min_lon cannot be less than 0.0 degrees', &
   !                      source,revision,revdate)
   !endif
   if (min_lon > 360.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lon cannot be greater than 360.0 degrees', &
                         source,revision,revdate)
   endif
   if (max_lon > 360.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'max_lon cannot be greater than 360.0 degrees', &
                         source,revision,revdate)
   endif

   ! it is risky to allow == operations on floating point numbers.
   ! if someone wants to select only obs exactly on a particular lon,
   ! force them to do an interval [min-epsilon, max-epsilon].
   ! move this test BEFORE the modulo so that if afterwards they are
   ! the same (e.g. 0 and 360) that's ok -- it will mean all longitudes 
   ! will be accepted.
   if (min_lon == max_lon) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lon cannot exactly equal max_lon', &
                         source,revision,revdate)
   endif
  
   ! do not test for min < max, because lon can wrap around 0.  ensure that
   ! the values after here are [0, 360) and the compare routine knows about wrap.
   ! must use modulo() and not just mod() to handle negative values correctly;
   ! the result of modulo() is positive even if the input is negative.
   min_lon = modulo(min_lon,360.0_r8)
   max_lon = modulo(max_lon,360.0_r8)

   ! create real location_type objects, knowing we are running with the
   ! 3d sphere locations module.
   min_box(1) = min_lon
   min_box(2) = min_lat
   min_box(3) = 0.0_r8
   min_box(4) = -2.0_r8 !! FIXME: VERTISUNDEF, but not all loc mods have it
   min_loc = set_location2(min_box)
   max_box(1) = max_lon
   max_box(2) = max_lat
   max_box(3) = 0.0_r8
   max_box(4) = -2.0_r8 !! FIXME: VERTISUNDEF, but not all loc mods have it
   max_loc = set_location2(max_box)
else
   restrict_by_location = .false.
endif
  
!%! ! SPECIAL: cut off all GPS obs below the given height
!%! if (min_gps_height /= missing_r8) then
!%!    restrict_by_height = .true.
!%! else
   restrict_by_height = .false.
!%! endif


! See if the user is restricting the data values or qc to be processed, 
! and if so, set up
if ((min_qc /= missing_r8) .or. (max_qc /= missing_r8)) then
   if (len(trim(qc_metadata)) == 0) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'must specify the metadata name of a QC field', &
                         source,revision,revdate)
   endif
   restrict_by_qc = .true.
else
   restrict_by_qc = .false.
endif

if ((min_copy /= missing_r8) .or. (max_copy /= missing_r8)) then
   if (len(trim(copy_metadata)) == 0) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'must specify the metadata name of a copy field', &
                         source,revision,revdate)
   endif
   restrict_by_copy = .true.
else
   restrict_by_copy = .false.
endif

! See if the user is requesting edit of either the number of copies
! or qc values, or reordering them.
editing_obs = .false.
if (edit_copies) then
   ! count up how many copies and qc there are going out.
   do i = 1, size(new_copy_index)
      if (new_copy_index(i) < 0) exit
      copy_index_len = copy_index_len + 1
   enddo
   editing_obs = .true.
endif
if (edit_qcs) then
   do i = 1, size(new_qc_index)
      if (new_qc_index(i) < 0) exit
      qc_index_len = qc_index_len + 1
   enddo
   editing_obs = .true.
endif

! check to see if we are going to trim the sequence by time
if(first_obs_seconds >= 0 .or. first_obs_days >= 0) then
   first_obs_time = set_time(first_obs_seconds, first_obs_days)
   trim_first = .true.
else
   trim_first = .false.
endif
if(last_obs_seconds >= 0 .or. last_obs_days >= 0) then
   last_obs_time = set_time(last_obs_seconds, last_obs_days)
   trim_last = .true.
else
   trim_last = .false.
endif
if (trim_first .and. trim_last) then
   if (first_obs_time > last_obs_time) then
      call error_handler(E_ERR,'obs_sequence_tool', 'first time cannot be later than last time', &
                         source,revision,revdate)
   endif
endif

! check to see if user gave us lists of metadata that should be
! considered the same across different obs_seq files.  count the
! length of each list.  used in the compare_metadata() routine.
matching_copy_metadata = .false.
if (synonymous_copy_list(1) /= '') then
   matching_copy_metadata = .true.
   do i = 1, size(synonymous_copy_list)
      if (synonymous_copy_list(i) == '') exit
      matching_copy_limit = i
   enddo
endif
matching_qc_metadata = .false.
if (synonymous_qc_list(1) /= '') then
   matching_qc_metadata = .true.
   do i = 1, size(synonymous_qc_list)
      if (synonymous_qc_list(i) == '') exit
      matching_qc_limit = i
   enddo
endif


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
      call error_handler(E_ERR,'obs_sequence_tool', &
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
   
   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)
   call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time, &
                 filename_seq(i), .true., remaining_obs_count)
   call destroy_obs_sequence(seq_in)
   if (remaining_obs_count == 0) then
      process_file(i) = .false.
      write(msgstring,*) 'No obs used from input sequence file ', trim(filename_seq(i))
      call error_handler(E_MSG,'obs_sequence_tool',msgstring)
      cycle
   else
      process_file(i) = .true.
      size_seq_in = remaining_obs_count
   endif

   ! First sequence with valid data gets to set up the format of the output.
   ! This only matters if we keep the 'non-exact match of metadata' in the
   ! compare_metadata() subroutine below.
   if ( first_seq < 0 ) then
      first_seq = i
      ! set up allowing number of copies or qcs to be changed or reordered.
      if (edit_copies) then
         ! do a tiny bit of sanity checking here.
         do j = 1, copy_index_len
            if ((new_copy_index(j) > num_copies_in) .or. &
                (new_copy_index(j) < 1)) then
               write(msgstring,*)'new_copy_index values must be between 1 and ', num_copies_in
               call error_handler(E_ERR,'obs_sequence_tool', msgstring, &
                                  source,revision,revdate)

            endif
         enddo
         if (copy_index_len > num_copies_in) then
            write(msgstring,*)'WARNING: more output copies than input; some are being replicated'
            call error_handler(E_MSG,'obs_sequence_tool', msgstring)
         endif
         num_copies_out = copy_index_len
      else
         num_copies_out = num_copies_in
         do j = 1, num_copies_in
            new_copy_index(j) = j
         enddo
      endif
      if (edit_qcs) then
         do j = 1, qc_index_len
            if ((new_qc_index(j) > num_qc_in) .or. &
                (new_qc_index(j) < 1)) then
               write(msgstring,*)'new_qc_index values must be between 1 and ', num_qc_in
               call error_handler(E_ERR,'obs_sequence_tool', msgstring, &
                                  source,revision,revdate)

            endif
         enddo
         if (qc_index_len > num_qc_in) then
            write(msgstring,*)'WARNING: more output qcs than input; some are being replicated'
            call error_handler(E_MSG,'obs_sequence_tool', msgstring)
         endif
         num_qc_out = qc_index_len
      else
         num_qc_out = num_qc_in
         do j = 1, num_qc_in
            new_qc_index(j) = j
         enddo
      endif
      size_seq_out = size_seq_in
   else
      size_seq_out = size_seq_out + size_seq_in
   endif

enddo

! no valid obs found?  if the index value is still negative, we are
! still waiting to process the first one and never found one.
if (first_seq < 0 .or. size_seq_out == 0) then
   msgstring = 'All input files are empty or all obs excluded by time/type/location'
   call error_handler(E_ERR,'obs_sequence_tool',msgstring,source,revision,revdate)
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
   call error_handler(E_MSG,'obs_sequence_tool',msgstring)

   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)

   ! If you get here, there better be observations in this file which
   ! are going to be used (the process_file flag wouldn't be set otherwise.)
   call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time, &
                 filename_seq(i), .false., remaining_obs_count)

   ! This would be an error at this point.
   if(remaining_obs_count == 0) then
      call destroy_obs_sequence(seq_in) 
      write(msgstring, *) 'Internal error trying to process file ', trim(filename_seq(i))
      call error_handler(E_ERR,'obs_sequence_tool',msgstring,source,revision,revdate)
   endif

   ! create the output sequence here based on the first input file
   if (first_seq < 0) then
      call init_obs_sequence(seq_out, num_copies_out, num_qc_out, size_seq_out) 
      do j=1, num_copies_out
         if (edit_copy_metadata) then
            meta_data = new_copy_metadata(j)
            write(msgstring, *)'replacing original copy meta_data: ' // &
                                trim(get_copy_meta_data(seq_in, new_copy_index(j)))
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            write(msgstring, *)' with: ' // trim(meta_data)
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            if (new_copy_index(j) /= j) then
               write(msgstring, *)'original index was ', new_copy_index(j), ', now ', j
               call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            endif
         else
            meta_data = get_copy_meta_data(seq_in, new_copy_index(j)) 
            if (new_copy_index(j) /= j) then
               write(msgstring, *)'copy with meta_data: ' // trim(meta_data)
               call error_handler(E_MSG,'obs_sequence_tool',msgstring)
               write(msgstring, *)'moved from index ', new_copy_index(j), ' to ', j
               call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            endif
         endif
         call set_copy_meta_data(seq_out, j, meta_data)
      enddo 
      do j=1, num_qc_out
         if (edit_qc_metadata) then
            meta_data = new_qc_metadata(j)
            write(msgstring, *)'replacing original qc meta_data: ' // &
                                trim(get_qc_meta_data(seq_in, new_qc_index(j)))
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            write(msgstring, *)' with: ' // trim(meta_data)
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            if (new_qc_index(j) /= j) then
               write(msgstring, *)'original index was ', new_qc_index(j), ' now ', j
               call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            endif
         else
            meta_data = get_qc_meta_data(seq_in, new_qc_index(j)) 
            if (new_qc_index(j) /= j) then
               write(msgstring, *)'qc with meta_data: ' // trim(meta_data)
               call error_handler(E_MSG,'obs_sequence_tool',msgstring)
               write(msgstring, *)'moved from index ', new_qc_index(j), ' now ', j
               call error_handler(E_MSG,'obs_sequence_tool',msgstring)
            endif
         endif
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

      if (editing_obs) then
         ! copy subset or reordering
         call copy_partial_obs(obs_out, obs_in,                &
                               num_copies_out, new_copy_index, &
                               num_qc_out, new_qc_index) 
      else
         obs_out = obs_in
      endif

!#!      call change_variance(obs_out)

      call insert_obs_in_seq(seq_out, obs_out)  ! new_obs linked list info changes

      prev_obs_out = obs_out            ! records new position in seq_out
      num_inserted = num_inserted + 1
   
      call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)
      ObsLoop : do while ( .not. is_this_last)

         if (print_every > 0) then
            if (mod(num_inserted,print_every) == 0) then
               print*, 'inserted number ',num_inserted,' of ',size_seq_in
            endif
         endif

         obs_in     = next_obs_in   ! essentially records position in seq_out

         if (editing_obs) then
            ! copy subset or reordering
            call copy_partial_obs(obs_out, obs_in,                &
                                  num_copies_out, new_copy_index, &
                                  num_qc_out, new_qc_index) 
         else
            obs_out = obs_in
         endif

!#!         call change_variance(obs_out)

         ! Since the stride through the observation sequence file is always 
         ! guaranteed to be in temporally-ascending order, we can use the
         ! 'previous' observation as the starting point to search for the
         ! correct insertion point.  This speeds up the insert code a lot.

         call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)

         prev_obs_out = obs_out    ! update position in seq_in for next insert
         num_inserted = num_inserted + 1

         call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

      enddo ObsLoop

      total_num_inserted = total_num_inserted + num_inserted

   else
      write(msgstring,*)'no first observation in ',trim(filename_seq(i))
      call error_handler(E_MSG,'obs_sequence_tool', msgstring)
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
call destroy_obs(prev_obs_out)

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
subroutine compare_metadata(seq1, seq2, fname1, fname2)

!
! This subroutine compares the metadata for two different observation
! sequences and terminates the program if they are not conformable.
! In order to be merged, the two observation sequences must have the same
! number of qc values, the same number of copies ... 
!
! The messages might be a bit misleading 'warning', 'warning', 'error' ...
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
integer :: num_copies , num_qc, i, j
logical :: have_match1, have_match2
character(len=129) :: str1, str2
character(len=255) :: msgstring1, msgstring2

num_copies1 = get_num_copies(seq1)
num_qc1     = get_num_qc(    seq1)

num_copies2 = get_num_copies(seq2)
num_qc2     = get_num_qc(    seq2)

if (edit_copies) num_copies2 = copy_index_len
if (edit_qcs   ) num_qc2     =   qc_index_len

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

   if (edit_copy_metadata) then
      str2 = new_copy_metadata(i)
   else
      str2 = get_copy_meta_data(seq2, new_copy_index(i)) 
   endif

   ! easy case - they match.  cycle to next copy.
   if( str1 == str2 ) then
      write(msgstring2,*)'metadata ',trim(str1), ' in both.'
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
      cycle CopyMetaData
   endif

   ! see if user provided a list of metadata strings that are
   ! the same values and can be considered a match.
   if (matching_copy_metadata) then
      have_match1 = .false.
      do j=1, matching_copy_limit
         if (trim(str1) == trim(synonymous_copy_list(j))) then
            have_match1 = .true.
            exit
         endif
      enddo

      have_match2 = .false.
      do j=1, matching_copy_limit
         if (trim(str2) == trim(synonymous_copy_list(j))) then
            have_match2 = .true.
            exit
         endif
      enddo

      ! if both are true, you found both strings in the list and
      ! it is ok to proceed.
      if (have_match1 .and. have_match2) then
         write(msgstring2,*)'different copy metadata strings ok because both on synonymous list'
         call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
         write(msgstring2,*)'one is: ', trim(str1)
         call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
         write(msgstring2,*)'one is: ', trim(str2)
         call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
         cycle CopyMetaData
      endif

      ! if no match, fall out of the if.
   endif
    
   !! FIXME: this could be dangerous - it allows any metadata string with
   !! the substring 'observation' to match any other.  if there are multiple
   !! strings with 'observation', it will allow them to match.  for now it
   !! allows 'NCEP BUFR observations' to match 'observation', for example,
   !! but it's dangerous.
   !if ((index(str1, 'observation') > 0) .and. &
   !    (index(str2, 'observation') > 0)) then
   !   write(msgstring2,*)'observation metadata in both ',trim(str1), ' and ', trim(str2)
   !   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   !   write(msgstring2,*)'ALLOWING NON-EXACT MATCH'
   !   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   !   cycle CopyMetaData
   !! end FIXME

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

   if (edit_qc_metadata) then
      str2 = new_qc_metadata(i)
   else
      str2 = get_qc_meta_data(seq2, new_qc_index(i)) 
   endif


   ! easy case - they match.  cycle to next copy.
   if( str1 == str2 ) then
      write(msgstring2,*)'metadata ',trim(str1), ' in both.'
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
      cycle QCMetaData
   endif

   ! see if user provided a list of metadata strings that are
   ! the same values and can be considered a match.
   if (matching_copy_metadata) then
      have_match1 = .false.
      do j=1, matching_qc_limit
         if (trim(str1) == trim(synonymous_qc_list(j))) then
            have_match1 = .true.
            exit
         endif
      enddo

      have_match2 = .false.
      do j=1, matching_qc_limit
         if (trim(str2) == trim(synonymous_qc_list(j))) then
            have_match2 = .true.
            exit
         endif
      enddo

      ! if both are true, you found both strings in the list and
      ! it is ok to proceed.
      if (have_match1 .and. have_match2) then
         write(msgstring2,*)'different qc metadata strings ok because both on synonymous list'
         call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
         write(msgstring2,*)'one is: ', trim(str1)
         call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
         write(msgstring2,*)'one is: ', trim(str2)
         call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
         cycle QCMetaData
      endif

      ! if no match, fall out of the if.
   endif
    
   !! FIXME: this is even more dangerous than the obs - better to make the
   !! user give an explicit list of strings that are ok to match.  but here
   !! is the code if you wanted to make it more mindless.
   !if (((index(str1, 'QC') > 0).or.(index(str1, 'quality control') > 0)).and. &
   !    ((index(str2, 'QC') > 0).or.(index(str2, 'quality control') > 0))) then
   !   write(msgstring2,*)'QC metadata in both ',trim(str1), ' and ', trim(str2)
   !   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   !   write(msgstring2,*)'ALLOWING NON-EXACT MATCH'
   !   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
   !   cycle QCMetaData
   !endif
   !! end FIXME

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
! pass in an already opened sequence and a start/end time.  this routine
! really trims the observations out of the sequence, and returns a count
! of how many remain.
subroutine trim_seq(seq, trim_first, first_time, trim_last, last_time, &
                    seqfilename, print_msg, remaining_obs_count)
 type(obs_sequence_type), intent(inout) :: seq
 logical, intent(in)                    :: trim_first, trim_last
 type(time_type), intent(in)            :: first_time, last_time
 character(len = *), intent(in)         :: seqfilename
 logical, intent(in)                    :: print_msg
 integer, intent(out)                   :: remaining_obs_count

 integer :: i
 logical :: found

   ! Need to find first obs with appropriate time, delete all earlier ones
   if(trim_first) then
      call delete_seq_head(first_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are before first_obs time'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
         endif
         remaining_obs_count = 0
         return
      endif
   endif
   
   ! Also get rid of observations past the last_obs_time if requested
   if(trim_last) then
      call delete_seq_tail(last_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are after last_obs time'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
         endif
         remaining_obs_count = 0
         return
      endif
   endif
   
   ! optionally select only a list of obs types
   if (restrict_by_obs_type) then
      call delete_obs_by_typelist(num_obs_input_types, obs_types, &
                                  keep_types, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            if (keep_types) then
               msgstring = 'Skipping: no obs in ' // trim(seqfilename) // &
                           ' are on the keep obs_types list'
            else
               msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                           ' are on the discard obs_types list'
            endif
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   ! optionally select only obs inside a bounding box
   ! most common expected use is a lat/lon box with the 3d sphere locations
   ! mod, but this should work with any locations mod if the namelist uses
   ! the right number of corners appropriate for the dimensions of the location
   ! when setting min_box and max_box.
   if (restrict_by_location) then
      call select_obs_by_location(min_loc, max_loc, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: no obs in ' // trim(seqfilename) // &
                        ' are above the given height'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   ! optionally select only a range of qc and copies
   if (restrict_by_qc) then
      ! validate the metadata string
      found = .false.
      do i=1, num_qc_in
         meta_data = get_qc_meta_data(seq_in, i) 
         if (trim(qc_metadata) == trim(meta_data)) then 
            call delete_obs_by_qc(i, min_qc, max_qc, seq, all_gone)
            found = .true.
            exit
         endif
      enddo 
      if (.not. found) then 
         msgstring = 'QC metadata string: ' // trim(qc_metadata) // &
                     ' not found in obs_seq file'
         call error_handler(E_MSG,'obs_sequence_tool',msgstring)
      endif
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are outside the qc min/max range'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
         endif
         remaining_obs_count = 0
         return
      endif
   endif
   if (restrict_by_copy) then
      ! validate the metadata string
      found = .false.
      do i=1, num_copies_in
         meta_data = get_copy_meta_data(seq_in, i) 
         if (trim(copy_metadata) == trim(meta_data)) then 
            call delete_obs_by_copy(i, min_copy, max_copy, copy_type, &
                                    seq, all_gone)
            found = .true.
            exit
         endif
      enddo 
      if (.not. found) then 
         msgstring = 'Copy metadata string: ' // trim(copy_metadata) // &
                     ' not found in obs_seq file'
         call error_handler(E_MSG,'obs_sequence_tool',msgstring)
      endif
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are outside the copy min/max range'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   ! SPECIAL: optionally restrict GPS obs to above a height
   if (restrict_by_height) then
      call select_gps_by_height(min_gps_height, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: no obs in ' // trim(seqfilename) // &
                        ' are above the GPS height threshold'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   remaining_obs_count = get_num_key_range(seq)

end subroutine trim_seq


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
if (gregorian_cal) then
   call print_date(get_obs_def_time(this_obs_def), '   Gregorian day: ')
endif

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
      if (gregorian_cal) then
         call print_date(get_obs_def_time(this_obs_def), '   Gregorian day: ')
      endif
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
      if (gregorian_cal) call print_date(last_time, '   Gregorian day: ')
      call print_time(this_time, ' next timestamp: ')
      if (gregorian_cal) call print_date(this_time, '   Gregorian day: ')

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
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)
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
subroutine select_gps_by_height(min_height, seq, all_gone)

! CURRENTLY COMMENTED OUT

! Delete all gps observations in the sequence which are below the given ht.
! If there are no obs left afterwards return that the sequence is all_gone.

real(r8),                intent(in)    :: min_height
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

all_gone = .false.
return

!%! type(obs_def_type)   :: obs_def
!%! type(obs_type)       :: obs, prev_obs
!%! integer              :: i, key, gps_type_index, this_obs_type
!%! type(location_type)  :: location
!%! logical              :: out_of_range, is_this_last, above, first_obs
!%! real(r8)             :: ll(3), vloc
!%! 
!%! ! figure out what index number is gps
!%! gps_type_index = get_obs_kind_index('GPSRO_REFRACTIVITY')
!%! if (gps_type_index < 0) then
!%!    write(msgstring,*) 'obs_type GPSRO not found'
!%!    call error_handler(E_ERR,'select_gps_by_height', msgstring, &
!%!                       source, revision, revdate)
!%! endif
!%! 
!%! 
!%! ! Initialize an observation type with appropriate size
!%! call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
!%! call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))
!%! 
!%! ! Iterate entire sequence, deleting obs which are not in the box.
!%! ! First, make sure there are obs to delete, and initialize first obs.
!%! if(.not. get_first_obs(seq, obs)) then
!%!    all_gone = .true.
!%!    call destroy_obs(obs)
!%!    call destroy_obs(prev_obs)
!%!    return
!%! endif
!%! 
!%! first_obs = .true.
!%! prev_obs = obs
!%! 
!%! ! This is going to be O(n), n=num obs in seq
!%! is_this_last = .false.
!%! allobs : do while (.not. is_this_last)
!%! 
!%!    ! verify GPS obs first, then do height
!%! 
!%!    call get_obs_def(obs, obs_def)
!%!    this_obs_type = get_obs_kind(obs_def)
!%!    if (this_obs_type /= gps_type_index) then
!%!       ! we are going to keep this
!%!       above = .true.   
!%!    else
!%!    
!%!       ! must check height.  at this point, all gps obs are be height only. 
!%!       location = get_obs_def_location(obs_def)
!%!    
!%!       ! this makes the tool sphere_3d dependent.  also assumes height as vert.
!%!       ! should verify. 
!%!       if (.not. vert_is_height(location)) then
!%!          write(msgstring,*) 'obs_type GPSRO vertical location not height'
!%!          call error_handler(E_ERR,'select_gps_by_height', msgstring, &
!%!                             source, revision, revdate)
!%!       endif
!%! 
!%!       ll = get_location(location)
!%!       vloc = ll(3)
!%! 
!%!       if (vloc < min_height) then
!%!          ! delete this one
!%!          above = .false.
!%!       else
!%!          ! will be kept
!%!          above = .true.
!%!       endif
!%!    endif 
!%! 
!%!    ! same code as delete/keep by obstype; do any code fixes both places
!%!    if (.not. above) then
!%!       if (first_obs) then
!%!          call delete_obs_from_seq(seq, obs)
!%!          if(.not. get_first_obs(seq, obs)) exit allobs
!%!       else
!%! !print *, 'going to del obs key ', obs%key
!%! !print *, 'prev key is ', prev_obs%key
!%!          call delete_obs_from_seq(seq, obs)
!%!          ! cannot simply use prev_obs; cached copy out of sync with seq one
!%!          key = get_obs_key(prev_obs)
!%!          call get_next_obs_from_key(seq, key, obs, is_this_last)
!%! !print *, 'next obs now is key ', obs%key
!%!       endif
!%!    else
!%! !print *, 'no del, keep this obs key ', obs%key
!%!       first_obs = .false.
!%!       prev_obs = obs
!%! !print *, 'prev obs now is key ', prev_obs%key
!%! !print *, 'obs was key ', obs%key
!%!       call get_next_obs(seq, prev_obs, obs, is_this_last)
!%! !print *, 'obs now is key ', obs%key
!%!    endif
!%!    
!%! end do allobs
!%! 
!%! ! Figure out if there are no more obs left in the sequence.
!%! if(.not. get_first_obs(seq, obs)) then
!%!    all_gone = .true.
!%! else
!%!    all_gone = .false.
!%! endif
!%! 
!%! ! Done.  delete temp storage and return.
!%! call destroy_obs(obs)
!%! call destroy_obs(prev_obs)

end subroutine select_gps_by_height

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
!#! subroutine change_variance(this_obs)
!#! 
!#! 
!#! use obs_kind_mod
!#! use obs_def_mod
!#! use obs_sequence_mod
!#! 
!#! ! change the variance on specific kinds.
!#! 
!#! type(obs_type), intent(inout) :: this_obs
!#! 
!#! type(obs_def_type) :: this_obs_def
!#! integer            :: this_obs_kind
!#! real(r8)           :: oldvar, newvar
!#! 
!#! 
!#! call get_obs_def(obs, this_obs_def)
!#! this_obs_kind = get_obs_kind(this_obs_def)
!#! 
!#! ! ignore identity obs
!#! if (this_obs_kind < 0) return
!#! 
!#! oldvar = get_obs_def_error_variance(this_obs_def)
!#! 
!#! !print *, 'kind, var = ', this_obs_kind, oldvar
!#! ! SOYOUNG: Change the code here for what you want.
!#! 
!#! ! Set the new variance here based on the type
!#! select case (this_obs_kind)
!#!   case (LAND_SFC_TEMPERATURE)
!#!      newvar = oldvar * 0.5
!#!   case (RADIOSONDE_TEMPERATURE)
!#!      newvar = 1.3
!#!   case (GPSRO_REFRACTIVITY)
!#!      newvar = 1.666
!#!   ! etc
!#!   case default
!#!      newvar = oldvar
!#! end select
!#! 
!#! !print *, 'newvar = ', newvar
!#! if (newvar /= oldvar) then
!#!    call set_obs_def_error_variance(this_obs_def, newvar)
!#!    call set_obs_def(this_obs, this_obs_def)
!#! endif
!#! 
!#! end subroutine change_variance

!---------------------------------------------------------------------
end program obs_sequence_tool
