! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Change observation sequence files by adding or removing observations
!> based on type, location, values, time.  Long list of options based on 
!> namelist settings.  See the html documention for extensive examples.

program obs_sequence_tool

use        types_mod, only : r8, missing_r8, metadatalength, obstypelength, &
                             MISSING_I

use    utilities_mod, only : finalize_utilities, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit,   &
                             do_nml_file, do_nml_term, set_filename_list, &
                             to_upper

use     location_mod, only : location_type, get_location, set_location, &
                             LocationName !%! , vert_is_height 
                             ! comment in select_gps_by_height() explains !%!

use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs, &
                             get_obs_def_location, set_obs_def_write_external_FO

use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs, &
                             get_index_for_type_of_obs

use time_manager_mod, only : time_type, operator(>), print_time, set_time, &
                             print_date, set_calendar_type, GREGORIAN

use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq, &
                             init_obs, assignment(=), get_obs_def, set_obs_def, &
                             init_obs_sequence, static_init_obs_sequence, &
                             read_obs_seq_header, read_obs_seq, get_num_obs, &
                             get_first_obs, get_next_obs, &
                             insert_obs_in_seq, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, &
                             set_copy_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence, &
                             delete_seq_head, delete_seq_tail, &
                             get_num_key_range, delete_obs_by_typelist, &
                             get_obs_key, copy_partial_obs, &
                             delete_obs_from_seq, get_next_obs_from_key, &
                             delete_obs_by_qc, delete_obs_by_copy, &
                             select_obs_by_location, set_obs_values, set_qc

implicit none

character(len=*), parameter :: source = 'obs_sequence_tool.f90'

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type) :: obs_in, next_obs_in
type(obs_type) :: obs_out, prev_obs_out
logical :: is_there_one, is_this_last
integer :: size_seq_in, num_copies_in, num_qc_in
integer :: size_seq_out, num_copies_out, num_qc_out
integer :: num_inserted, iunit, io, i, j, total_num_inserted
integer :: max_num_obs, file_id, remaining_obs_count
integer :: first_seq
character(len=metadatalength) :: read_format, meta_data
logical :: pre_I_format, all_gone
logical :: trim_first, trim_last
character(len=512) :: msgstring1, msgstring2, msgstring3

! could go into namelist if you wanted runtime control
integer, parameter :: print_every = 20000

! specify sizes of allocated arrays; increase and recompile if too small
integer, parameter :: copy_qc_listlen = 256           ! max number of data or qc copies
integer, parameter :: max_num_input_files = 1000
integer, parameter :: max_obs_input_types = 500

logical :: process_file(max_num_input_files)
logical :: restrict_by_obs_type
logical :: restrict_by_location
logical :: restrict_by_qc
logical :: restrict_by_copy
!%! logical :: restrict_by_height
integer :: num_obs_input_types
type(location_type) :: min_loc, max_loc
logical :: editing_obs, matching_copy_metadata, matching_qc_metadata
type(time_type) :: first_obs_time, last_obs_time
integer :: matching_copy_limit = 0
integer :: matching_qc_limit   = 0
integer :: copy_index_len      = 0
integer :: qc_index_len        = 0


! variables to manage whether or not the precomputed forward operator values
! are written out

integer :: FO_write_types(max_obs_input_types) = MISSING_I
integer :: num_FO_types_to_suppress

!----------------------------------------------------------------
! Namelist input with default values

integer :: num_input_files = 0    ! DEPRECATED; count set by number of input files specified.
character(len=256) :: filename_seq(max_num_input_files) = ''
character(len=256) :: filename_seq_list = ''
character(len=256) :: filename_out  = 'obs_seq.processed'

! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days    = -1
integer  :: first_obs_seconds = -1
integer  :: last_obs_days     = -1
integer  :: last_obs_seconds  = -1

character(len=obstypelength) :: remove_precomputed_FO_values(max_obs_input_types) = ''
character(len=obstypelength) :: obs_types(max_obs_input_types) = ''
logical  :: keep_types = .true.

! must be location independent...
real(r8) :: min_box(4) = missing_r8
real(r8) :: max_box(4) = missing_r8
real(r8) :: min_lat = -90.0_r8
real(r8) :: max_lat =  90.0_r8
real(r8) :: min_lon =   0.0_r8
real(r8) :: max_lon = 360.0_r8

real(r8) :: min_copy = missing_r8
real(r8) :: max_copy = missing_r8
character(len=obstypelength)  :: copy_type = ''
character(len=metadatalength) :: copy_metadata = ''
character(len=metadatalength) :: new_copy_metadata(copy_qc_listlen) = ''
logical  :: edit_copy_metadata = .false.
logical  :: edit_copies        = .false.
integer  :: new_copy_index(copy_qc_listlen) = -1
real(r8) :: new_copy_data(copy_qc_listlen)  = MISSING_R8
character(len=metadatalength) :: synonymous_copy_list(copy_qc_listlen) = ''

real(r8) :: min_qc   = missing_r8
real(r8) :: max_qc   = missing_r8
character(len=metadatalength) :: qc_metadata = ''
character(len=metadatalength) :: new_qc_metadata(copy_qc_listlen)   = ''
logical  :: edit_qc_metadata = .false.
logical  :: edit_qcs = .false.
integer  :: new_qc_index(copy_qc_listlen) = -1
real(r8) :: new_qc_data(copy_qc_listlen) = MISSING_R8
character(len=metadatalength) :: synonymous_qc_list(copy_qc_listlen)   = ''

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
         new_qc_metadata, new_copy_data, new_qc_data, remove_precomputed_FO_values

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

! this is a deprecated namelist item.  ignored, and eventually will be removed.
! check it first before overwriting it in the call to set_filename_list().
if (num_input_files > 0) then
   write(msgstring1, *) '"num_input_files" is a DEPRECATED namelist item and is ignored.'
   write(msgstring2, *) 'the count of input files is set by the number of filenames specified.'
   call error_handler(E_MSG,'obs_sequence_tool', msgstring1, source, text2=msgstring2)
endif

! when set_filename_list() returns, filename_seq contains all the filenames
! whether they were specified directly in the namelist or in a separate file.
num_input_files = set_filename_list(filename_seq, filename_seq_list, 'obs_sequence_tool')

! this is a logical.  if the calendar is gregorian it prints out times in
! both (day,second) format and date (year/month/day/hour/minute/second) format.
! otherwise it just prints (day,second).  if users want other calendar type
! support this should be renamed to 'calendar_type' and be a character string.
if (gregorian_cal) call set_calendar_type(GREGORIAN)

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

! See if the user wants to suppress the precomputed forward operator values 
! for any/all observation types with precomputed forward operator values

num_FO_types_to_suppress = 0
do i = 1, max_obs_input_types
   if ( len(remove_precomputed_FO_values(i)) == 0  .or.  &
            remove_precomputed_FO_values(i)  == "" ) exit
   num_FO_types_to_suppress = i
enddo
if (num_FO_types_to_suppress > 0) call setup_FO_suppress_list(FO_write_types)

! See if the user is restricting the obs locations to be processed, and set up
! the values if so.  Note that if any of the values are not missing_r8, all of
! them must have good (and consistent) values.  i don't have code in here yet
! to idiotproof this, but i should add it.

restrict_by_location = .false.
if ((minval(min_box).ne.missing_r8) .or. (maxval(min_box).ne.missing_r8) .or. &
    (minval(max_box).ne.missing_r8) .or. (maxval(max_box).ne.missing_r8)) then
   restrict_by_location = .true.
   min_loc = set_location(min_box)
   max_loc = set_location(max_box)
endif

! 3d sphere box check - locations module dependent, but an important one.
if ((min_lat /= -90.0_r8) .or. (max_lat /=  90.0_r8) .or. &
    (min_lon /=   0.0_r8) .or. (max_lon /= 360.0_r8)) then
   ! do not allow namelist to set BOTH min/max box and by lat/lon.
   if (restrict_by_location) then
      call error_handler(E_ERR,'obs_sequence_tool', &
              'use either lat/lon box or min/max box but not both', source)
   endif
   restrict_by_location = .true.
   if ((trim(LocationName) /= 'loc3Dsphere') .and.  &
       (trim(LocationName) /= 'loc2Dsphere')) then  
      call error_handler(E_ERR,'obs_sequence_tool', &
              'can only use lat/lon box with 2d/3d sphere locations', source)
   endif
   ! simple err checks before going on; try to catch radians vs degrees or
   ! just plain junk or typos.
   if (min_lat >= max_lat) then
      call error_handler(E_ERR,'obs_sequence_tool', &
              'min_lat must be less than max_lat', source)
   endif
   if (min_lat < -90.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
              'min_lat cannot be less than -90.0 degrees', source)
   endif
   if (max_lat >  90.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
              'max_lat cannot be greater than  90.0 degrees', source)
   endif
   ! this is ok - e.g.-180 to 180
   !if (min_lon < 0.0_r8) then
   !   call error_handler(E_ERR,'obs_sequence_tool', &
   !           'min_lon cannot be less than 0.0 degrees', source)
   !endif
   if (min_lon > 360.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lon cannot be greater than 360.0 degrees', source)
   endif
   if (max_lon > 360.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'max_lon cannot be greater than 360.0 degrees', source)
   endif

   ! it is risky to allow == operations on floating point numbers.
   ! if someone wants to select only obs exactly on a particular lon,
   ! force them to do an interval [min-epsilon, max-epsilon].
   ! move this test BEFORE the modulo so that if afterwards they are
   ! the same (e.g. 0 and 360) that's ok -- it will mean all longitudes 
   ! will be accepted.
   if (min_lon == max_lon) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lon cannot exactly equal max_lon', source)
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
   min_loc = set_location(min_box)
   max_box(1) = max_lon
   max_box(2) = max_lat
   max_box(3) = 0.0_r8
   max_box(4) = -2.0_r8 !! FIXME: VERTISUNDEF, but not all loc mods have it
   max_loc = set_location(max_box)
endif
  
!%! ! SPECIAL: cut off all GPS obs below the given height
!%! ! see comments in select_gps_by_height() for more info.
!%! if (min_gps_height /= missing_r8) then
!%!    restrict_by_height = .true.
!%! else
!%!    restrict_by_height = .false.
!%! endif


! See if the user is restricting the data values or qc to be processed, 
! and if so, set up
if ((min_qc /= missing_r8) .or. (max_qc /= missing_r8)) then
   if (len(trim(qc_metadata)) == 0) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'must specify the metadata name of a QC field', source)
   endif
   restrict_by_qc = .true.
else
   restrict_by_qc = .false.
endif

if ((min_copy /= missing_r8) .or. (max_copy /= missing_r8)) then
   if (len(trim(copy_metadata)) == 0) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'must specify the metadata name of a copy field', source)
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

if (trim_first) then
   call print_time(first_obs_time,    'Excluding observations before: ')
   if (gregorian_cal) &
      call print_date(first_obs_time, '       which is Gregorian day: ')
endif
if (trim_last) then
   call print_time(last_obs_time,     'Excluding observations  after: ')
   if (gregorian_cal) &
      call print_date(last_obs_time,  '       which is Gregorian day: ')
endif

if (trim_first .and. trim_last) then
   if (first_obs_time > last_obs_time) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'first time cannot be later than last time', source)
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
         'num_input_files and filename_seq mismatch',source)
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
      write(msgstring1,*) 'No obs used from input sequence file ', trim(filename_seq(i))
      call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
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
                (new_copy_index(j) < 0)) then   ! was < 1 here and below
               write(msgstring1,*)'new_copy_index values must be between 0 and ', num_copies_in
               call error_handler(E_ERR,'obs_sequence_tool', msgstring1, source)

            endif
         enddo
         if (copy_index_len > num_copies_in) then
            write(msgstring1,*)'WARNING: more output copies than input; some are being replicated'
            call error_handler(E_MSG,'obs_sequence_tool', msgstring1)
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
                (new_qc_index(j) < 0)) then  ! was < 1 here and below
               write(msgstring1,*)'new_qc_index values must be between 0 and ', num_qc_in
               call error_handler(E_ERR,'obs_sequence_tool', msgstring1, source)

            endif
         enddo
         if (qc_index_len > num_qc_in) then
            write(msgstring1,*)'WARNING: more output qcs than input; some are being replicated'
            call error_handler(E_MSG,'obs_sequence_tool', msgstring1)
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
   msgstring1 = 'All input files are empty or all obs excluded by time/type/location'
   call error_handler(E_ERR,'obs_sequence_tool',msgstring1,source)
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
 
   write(msgstring1,*) 'Starting to process input sequence file ', trim(filename_seq(i))
   call error_handler(E_MSG,'obs_sequence_tool',msgstring1)

   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)

   ! If you get here, there better be observations in this file which
   ! are going to be used (the process_file flag wouldn't be set otherwise.)
   call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time, &
                 filename_seq(i), .false., remaining_obs_count)

   ! This would be an error at this point.
   if(remaining_obs_count == 0) then
      call destroy_obs_sequence(seq_in) 
      write(msgstring1, *) 'Internal error trying to process file ', trim(filename_seq(i))
      call error_handler(E_ERR,'obs_sequence_tool',msgstring1,source)
   endif

   ! create the output sequence here based on the first input file
   if (first_seq < 0) then
      call init_obs_sequence(seq_out, num_copies_out, num_qc_out, size_seq_out) 
      do j=1, num_copies_out
         if (edit_copy_metadata) then
            meta_data = new_copy_metadata(j)
            if (new_copy_index(j) > 0) then
               write(msgstring1, *)'replacing original copy meta_data: ' // &
                                   trim(get_copy_meta_data(seq_in, new_copy_index(j)))
            else
               write(msgstring1, *)'replacing original copy meta_data: ' // &
                                   trim(get_copy_meta_data(seq_out, j))
            endif
            write(msgstring2, *)' with: ' // trim(meta_data)
            if (new_copy_index(j) /= j) then
               write(msgstring3, *)'original index was ', new_copy_index(j), ', now ', j
            else
               msgstring3 = 'copy index has not changed'
            endif
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1,&
                               text2=msgstring2,text3=msgstring3)

         else
            if (new_copy_index(j) > 0) then
               meta_data = get_copy_meta_data(seq_in, new_copy_index(j)) 
            else
               meta_data = get_copy_meta_data(seq_out, j) 
            endif
            if (new_copy_index(j) /= j) then
               write(msgstring1, *)'copy with meta_data: ' // trim(meta_data)
               write(msgstring2, *)'moved from index ', new_copy_index(j), ' to ', j
               call error_handler(E_MSG,'obs_sequence_tool',msgstring1,text2=msgstring2)
            endif
         endif
         call set_copy_meta_data(seq_out, j, meta_data)
      enddo 
      do j=1, num_qc_out
         if (edit_qc_metadata) then
            meta_data = new_qc_metadata(j)
            if (new_qc_index(j) > 0) then
               write(msgstring1, *)'replacing original qc meta_data: ' // &
                                   trim(get_qc_meta_data(seq_in, new_qc_index(j)))
            else
               write(msgstring1, *)'replacing original qc meta_data: ' // &
                                   trim(get_qc_meta_data(seq_out, j))
            endif
            write(msgstring2, *)' with: ' // trim(meta_data)
            if (new_qc_index(j) /= j) then
               write(msgstring3, *)'original index was ', new_qc_index(j), ' now ', j
            else
               msgstring3 = 'QC index has not changed'
            endif
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1,&
                               text2=msgstring2,text3=msgstring3)

         else
            if (new_qc_index(j) > 0) then
               meta_data = get_qc_meta_data(seq_in, new_qc_index(j)) 
            else
               meta_data = get_qc_meta_data(seq_out, j) 
            endif
            if (new_qc_index(j) /= j) then
               write(msgstring1, *)'qc with meta_data: ' // trim(meta_data)
               write(msgstring2, *)'moved from index ', new_qc_index(j), ' now ', j
               call error_handler(E_MSG,'obs_sequence_tool',msgstring1,text2=msgstring2)
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

         ! if new copies or qcs were added, set the initial data values
         call set_new_data(obs_out, edit_copies, new_copy_index, new_copy_data, &
                                    edit_qcs,    new_qc_index,   new_qc_data)
      else
         obs_out = obs_in
      endif

      ! changing metadata to output precomputed FOs based on namelist
      call set_FO_metadata(obs_out)

!#!      ! see comment in subroutine code for an explanation
!#!      call change_obs(obs_out)

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

         obs_in = next_obs_in   ! essentially records position in seq_out

         !setting flag to output precomputed FOs based on namelist
         call set_FO_metadata(obs_in)

         if (editing_obs) then
            ! copy subset or reordering
            call copy_partial_obs(obs_out, obs_in,                &
                                  num_copies_out, new_copy_index, &
                                  num_qc_out, new_qc_index) 

            ! if new copies or qcs were added, set the initial data values
            call set_new_data(obs_out, edit_copies, new_copy_index, new_copy_data, &
                                       edit_qcs,    new_qc_index,   new_qc_data)
         else
            obs_out = obs_in
         endif

!#!         ! see comment in subroutine code for an explanation
!#!         call change_obs(obs_out)

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
      write(msgstring1,*)'no first observation in ',trim(filename_seq(i))
      call error_handler(E_MSG,'obs_sequence_tool', msgstring1)
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


write(msgstring1,*) 'Starting to process output sequence file ', trim(filename_out)
call error_handler(E_MSG,'obs_sequence_tool',msgstring1)

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
   write(msgstring1,*) 'Output sequence file not created; print_only in namelist is .true.'
   call error_handler(E_MSG,'', msgstring1)
endif

! Time to clean up

call destroy_obs_sequence(seq_out)
call destroy_obs(     obs_in )
call destroy_obs(next_obs_in )
call destroy_obs(     obs_out)
call destroy_obs(prev_obs_out)

call error_handler(E_MSG, 'obs_sequence_tool', 'Finished successfully.',source)
call finalize_utilities()


contains


!---------------------------------------------------------------------
subroutine obs_seq_modules_used()

! Initialize modules used that require it
call initialize_utilities('obs_sequence_tool')
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
! change the error to use the 2 additional text lines rather than
! making some warnings and the last one an error.
!
! FIXME:  this routine uses several globals now from the namelist that
! should be arguments to this routine.  i'm being lazy (or expedient) here
! but this should be fixed at some point soon.
!
! and now there is an additional restriction -- if editing the copies or
! qc entries, seq1 has already been edited and only 2 needs the editing
! applied.  before they were completely symmetric.
!
! in this case, seq1 must be the new out sequence, and seq2 is one
! of the input sequences

type(obs_sequence_type), intent(in) :: seq1, seq2
character(len=*), optional :: fname1, fname2

integer :: num_copies1, num_qc1
integer :: num_copies2, num_qc2
integer :: num_copies , num_qc, i, j
logical :: have_match1, have_match2
character(len=metadatalength) :: str1, str2

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
   num_copies = -1
else
   write(msgstring2,*)'Number of data copies match'
endif
if ( num_qc1 /= num_qc2 ) then
   write(msgstring3,*)'Different different numbers of QCs found: ', &
                      num_qc1, ' vs ', num_qc2
   num_qc = -1
else
   write(msgstring3,*)'Number of qc copies match'
endif
if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, &
                      source, text2=msgstring2, text3=msgstring3)
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
      if (new_copy_index(i) > 0) then
         str2 = get_copy_meta_data(seq2, new_copy_index(i)) 
      else
         str2 = get_copy_meta_data(seq1, i) 
      endif
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
         write(msgstring2,*)'one is: ', trim(str1)
         write(msgstring3,*)'one is: ', trim(str2)
         call error_handler(E_MSG, 'obs_sequence_tool', &
                              'different copy metadata strings ok because both on synonymous list', &
                              text2=msgstring2, text3=msgstring3)
         cycle CopyMetaData
      endif

      ! if no match, fall out of the if.
   endif
    
   ! if you get here, the metadata is not the same and the user has not
   ! given us strings that are ok to match.  fail.
   write(msgstring2,*)'copy metadata mismatch, file 1: ', trim(str1)
   write(msgstring3,*)'copy metadata mismatch, file 2: ', trim(str2)
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, &
                      source, &
                      text2=msgstring2, text3=msgstring3)

enddo CopyMetaData

QCMetaData : do i=1, num_qc
   str1 = get_qc_meta_data(seq1,i)

   if (edit_qc_metadata) then
      str2 = new_qc_metadata(i)
   else
      if (new_qc_index(i) > 0) then
         str2 = get_qc_meta_data(seq2, new_qc_index(i)) 
      else
         str2 = get_qc_meta_data(seq1, i) 
      endif
   endif


   ! easy case - they match.  cycle to next copy.
   if( str1 == str2 ) then
      write(msgstring2,*)'metadata ',trim(str1), ' in both.'
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2)
      cycle QCMetaData
   endif

   ! see if user provided a list of metadata strings that are
   ! the same values and can be considered a match.
   if (matching_qc_metadata) then
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
         write(msgstring2,*)'one is: ', trim(str1)
         write(msgstring3,*)'one is: ', trim(str2)
         call error_handler(E_MSG, 'obs_sequence_tool', &
                              'different qc metadata strings ok because both on synonymous list', &
                              text2=msgstring2, text3=msgstring3)
         cycle QCMetaData
      endif

      ! if no match, fall out of the if.
   endif
    
   ! if you get here, the metadata is not the same and the user has not
   ! given us strings that are ok to match.  fail.
   write(msgstring2,*)'qc metadata mismatch, file 1: ', trim(str1)
   write(msgstring3,*)'qc metadata mismatch, file 2: ', trim(str2)
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, &
                      source, text2=msgstring2, text3=msgstring3)

enddo QCMetaData

end subroutine compare_metadata

!---------------------------------------------------------------------
! pass in an already opened sequence and a start/end time.  this routine
! really trims the observations out of the sequence, and returns a count
! of how many remain.
subroutine trim_seq(seq, trim_first, first_time, trim_last, last_time, &
                    seqfilename, print_msg, remaining_obs_count)
 type(obs_sequence_type), intent(inout) :: seq
 logical,                 intent(in)    :: trim_first, trim_last
 type(time_type),         intent(in)    :: first_time, last_time
 character(len=*),        intent(in)    :: seqfilename
 logical,                 intent(in)    :: print_msg
 integer,                 intent(out)   :: remaining_obs_count

 integer :: i
 logical :: found

   ! Need to find first obs with appropriate time, delete all earlier ones
   if(trim_first) then
      call delete_seq_head(first_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring1 = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are before first_obs time'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
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
            msgstring1 = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are after last_obs time'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
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
               msgstring1 = 'Skipping: no obs in ' // trim(seqfilename) // &
                           ' are on the keep obs_types list'
            else
               msgstring1 = 'Skipping: all obs in ' // trim(seqfilename) // &
                           ' are on the discard obs_types list'
            endif
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
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
            msgstring1 = 'Skipping: no obs in ' // trim(seqfilename) // &
                        ' are within the given min/max location bounding box'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
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
         msgstring1 = 'QC metadata string: ' // trim(qc_metadata) // &
                     ' not found in obs_seq file'
         call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
      endif
      if(all_gone) then
         if (print_msg) then
            msgstring1 = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are outside the qc min/max range'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
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
         msgstring1 = 'Copy metadata string: ' // trim(copy_metadata) // &
                     ' not found in obs_seq file'
         call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
      endif
      if(all_gone) then
         if (print_msg) then
            msgstring1 = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are outside the copy min/max range'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

!%!    ! SPECIAL: optionally restrict GPS obs to above a given height
!%!    ! see comments in select_gps_by_height() for more info.
!%!    if (restrict_by_height) then
!%!       call select_gps_by_height(min_gps_height, seq, all_gone)
!%!       if(all_gone) then
!%!          if (print_msg) then
!%!             msgstring1 = 'Skipping: no obs in ' // trim(seqfilename) // &
!%!                         ' are above the GPS height threshold'
!%!             call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
!%!          endif
!%!          remaining_obs_count = 0
!%!          return
!%!       endif
!%!    endif

   remaining_obs_count = get_num_key_range(seq)

end subroutine trim_seq


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
! max_defined_types_of_obs is a public from obs_kind_mod.f90 and really is
! counting the max number of types, not kinds
integer                 :: type_count(max_defined_types_of_obs), identity_count


! Initialize input obs_types
do i = 1, max_defined_types_of_obs
   type_count(i) = 0
enddo
identity_count = 0

! make sure there are obs left to process before going on.
! num_obs should be ok since we just constructed this seq so it should
! have no unlinked obs.  if it might for some reason, use this instead:
! size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring1 = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_sequence_tool',msgstring1)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq_in), get_num_qc(seq_in))
call init_obs(next_obs, get_num_copies(seq_in), get_num_qc(seq_in))

! blank line
call error_handler(E_MSG,'',' ')

write(msgstring1,*) 'Processing sequence file ', trim(filename)
call error_handler(E_MSG,'',msgstring1)

call print_metadata(seq_in, filename)

!-------------------------------------------------------------
! Start to process obs from seq_in
!--------------------------------------------------------------
is_there_one = get_first_obs(seq_in, obs)

if ( .not. is_there_one )  then
   write(msgstring1,*)'no first observation in ',trim(filename)
   call error_handler(E_MSG,'obs_sequence_tool', msgstring1)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def),    '  First obs time: ')
if (gregorian_cal) &
   call print_date(get_obs_def_time(this_obs_def), '   Gregorian day: ')

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_kind = get_obs_def_type_of_obs(this_obs_def)
   if (this_obs_kind < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_kind) = type_count(this_obs_kind) + 1
   endif
!   print *, 'obs kind index = ', this_obs_kind
!   if(this_obs_kind > 0)print *, 'obs name = ', get_name_for_type_of_obs(this_obs_kind)

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def),    '   Last obs time: ')
      if (gregorian_cal) &
         call print_date(get_obs_def_time(this_obs_def), '   Gregorian day: ')
   endif

enddo ObsLoop


write(msgstring1, *) 'Number of obs processed  :          ', size_seq_in
call error_handler(E_MSG, '', msgstring1)
write(msgstring1, *) '---------------------------------------------------------'
call error_handler(E_MSG, '', msgstring1)
do i = 1, max_defined_types_of_obs
   if (type_count(i) > 0) then 
      write(msgstring1, '(a32,i8,a)') trim(get_name_for_type_of_obs(i)), &
                                     type_count(i), ' obs'
      call error_handler(E_MSG, '', msgstring1)
   endif
enddo
if (identity_count > 0) then 
   write(msgstring1, '(a32,i8,a)') 'Identity observations', &
                                  identity_count, ' obs'
   call error_handler(E_MSG, '', msgstring1)
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
   msgstring1 = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_sequence_tool:validate',msgstring1)
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
   write(msgstring1,*)'no first observation in sequence ' // trim(filename)
   call error_handler(E_MSG,'obs_sequence_tool:validate', msgstring1, source)
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
      write(msgstring1,*)'obs number ', key, ' has earlier time than previous obs'
      write(msgstring2,*)'observations must be in increasing time order, file ' // trim(filename)
      call error_handler(E_ERR,'obs_sequence_tool:validate', msgstring1, &
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

type(obs_sequence_type), intent(in) :: seq1
character(len=*), optional :: fname1

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str1

num_copies = get_num_copies(seq1)
num_qc     = get_num_qc(    seq1)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring1,*)' illegal copy or obs count in file '//trim(fname1)
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source)
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
subroutine set_new_data(obs, edit_copies, new_copy_index, new_copy_data, &
                             edit_qcs,    new_qc_index,   new_qc_data)

! if new copies or qcs were added, set the initial data values

type(obs_type), intent(inout) :: obs
logical,        intent(in)    :: edit_copies
integer,        intent(in)    :: new_copy_index(:)
real(r8),       intent(in)    :: new_copy_data(:)
logical,        intent(in)    :: edit_qcs
integer,        intent(in)    :: new_qc_index(:)
real(r8),       intent(in)    :: new_qc_data(:)

integer :: i, j

if (edit_copies) then     
   j = 1
   copy_loop: do i = 1, size(new_copy_index)
      if (new_copy_index(i) == -1) exit copy_loop
      if (new_copy_index(i) /= 0) cycle copy_loop
      call set_obs_values(obs, new_copy_data(j:j), i)
      j = j + 1
   enddo copy_loop
endif

if (edit_qcs) then     
   j = 1
   qc_loop: do i = 1, size(new_qc_index)
      if (new_qc_index(i) == -1) exit qc_loop
      if (new_qc_index(i) /= 0) cycle qc_loop
      call set_qc(obs, new_qc_data(j:j), i)
      j = j + 1
   enddo qc_loop
endif

end subroutine set_new_data


!-------------------------------------------------------------------------------
! Construct list of observation types whose precomputed forward operator (FO) 
! values are to be removed.

subroutine setup_FO_suppress_list(type_list)

integer, intent(out) :: type_list(:)

integer :: i
character(len=obstypelength) :: upperstring

LOOP : do i = 1, num_FO_types_to_suppress
   upperstring = trim(remove_precomputed_FO_values(i))
   call to_upper(upperstring)

   if (upperstring == 'ALL') then
      ! If all possible observations with FO values are selected 
      ! just indicate that by setting num_FO_types_to_suppress
      ! to value that is larger than the largest possible
      type_list = 1
      num_FO_types_to_suppress = max_obs_input_types + 1
      exit LOOP
   else
      type_list(i) = get_index_for_type_of_obs(upperstring)
      if (type_list(i) < 0) then
         write(msgstring1,*) 'obs_type "'//trim(remove_precomputed_FO_values(i))//'"'
         write(msgstring2,*) 'not a valid observation type.'
         write(msgstring3,*) 'check entries for "remove_precomputed_FO_values"'
         call error_handler(E_ERR,'setup_FO_suppress_list', msgstring1, source, &
                            text2=msgstring2, text3=msgstring3)
      endif
   endif
enddo LOOP

end subroutine setup_FO_suppress_list


!-------------------------------------------------------------------------------
! If the observation type matches one specified by the remove_precomputed_FO_values 
! variable, adjust the metadata for that observation so that write_obs_def() 
! does the right thing. Since read_obs_def() sets the flag to suppress the writing
! of the precomputed values, this routine must set the flag to enable the writing
! of the precomputed values.

subroutine set_FO_metadata(obs)

type(obs_type), intent(inout) :: obs

type(obs_def_type) :: obs_def
integer            :: this_obs_type
logical            :: write_FO_values

write_FO_values = .true.

call get_obs_def(obs, obs_def)
this_obs_type = get_obs_def_type_of_obs(obs_def)

! If all types are to be suppressed, just match - and avoid using
! num_FO_types_to_suppress as an index ...

if (num_FO_types_to_suppress > max_obs_input_types) then
   write_FO_values = .false.
elseif ( any(this_obs_type == FO_write_types(1:num_FO_types_to_suppress)) ) then
   write_FO_values = .false.
endif

if (write_FO_values) then
   call set_obs_def_write_external_FO(obs_def,write_FO_values)
   call set_obs_def(obs,obs_def)
endif

end subroutine set_FO_metadata


!%! !---------------------------------------------------------------------
!%! subroutine select_gps_by_height(min_height, seq, all_gone)
!%! 
!%! ! This code is intentionally commented out.  Enabling this code requires 
!%! ! you are compiling with the threed_sphere locations module (which is true for
!%! ! most large models).  To do selections in the vertical for observations
!%! ! remove all !%! characters here and in the code above.
!%! 
!%! ! Adapt the code as needed - select on a different observation type or ignore
!%! ! the types and process all observations.  Test the type of the vertical
!%! ! coordinate before using the value; vertical in the threed_sphere module
!%! ! can be pressure, height, model level, surface, or intentionally undefined.
!%! ! There is a namelist item 'min_gps_height' which can be repurposed for
!%! ! setting run-time values and is passed into this routine.
!%! 
!%! ! The code below is working code.  It processes only GPS Radio Occultation
!%! ! observation types.  They are defined with a vertical of 'height' in meters.
!%! ! The code removes any observations below a given threshold.  If all obs are
!%! ! removed the logical 'all_gone' is returned true.
!%! 
!%! real(r8),                intent(in)    :: min_height
!%! type(obs_sequence_type), intent(inout) :: seq
!%! logical,                 intent(out)   :: all_gone
!%! 
!%! type(obs_def_type)   :: obs_def
!%! type(obs_type)       :: obs, prev_obs
!%! integer              :: i, key, gps_type_index, this_obs_type
!%! type(location_type)  :: location
!%! logical              :: out_of_range, is_this_last, above, first_obs
!%! real(r8)             :: ll(3), vloc
!%! 
!%! ! figure out what index number GPS obs are
!%! gps_type_index = get_index_for_type_of_obs('GPSRO_REFRACTIVITY')
!%! if (gps_type_index < 0) then
!%!    write(msgstring1,*) 'obs_type GPSRO_REFRACTIVITY not found'
!%!    call error_handler(E_ERR,'select_gps_by_height', msgstring1, source)
!%! endif
!%! 
!%! ! Initialize an observation type with appropriate size
!%! call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
!%! call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))
!%! 
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
!%!    this_obs_type = get_obs_def_type_of_obs(obs_def)
!%!    if (this_obs_type /= gps_type_index) then
!%!       ! we are going to keep this
!%!       above = .true.   
!%!    else
!%!    
!%!       ! must check height.  at this point, all gps obs are by height only. 
!%!       location = get_obs_def_location(obs_def)
!%!    
!%!       ! this makes the tool locations/threed_sphere dependent.
!%!       if (.not. vert_is_height(location)) then
!%!          write(msgstring1,*) 'obs_type GPSRO_REFRACTIVITY vertical location not height'
!%!          call error_handler(E_ERR,'select_gps_by_height', msgstring1, source)
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
!%!    ! remove obs if selected for delete, and set the prev_obs to the right value
!%!    if (.not. above) then
!%!       if (first_obs) then
!%!          call delete_obs_from_seq(seq, obs)
!%!          if(.not. get_first_obs(seq, obs)) exit allobs
!%!       else
!%!          call delete_obs_from_seq(seq, obs)
!%!          ! cannot simply use prev_obs; cached copy out of sync with seq one
!%!          key = get_obs_key(prev_obs)
!%!          call get_next_obs_from_key(seq, key, obs, is_this_last)
!%!       endif
!%!    else
!%!       first_obs = .false.
!%!       prev_obs = obs
!%!       call get_next_obs(seq, prev_obs, obs, is_this_last)
!%!    endif
!%!    
!%! end do allobs
!%! 
!%! ! Return indicator of whether all obs were deleted by this routine or not.
!%! if(.not. get_first_obs(seq, obs)) then
!%!    all_gone = .true.
!%! else
!%!    all_gone = .false.
!%! endif
!%! 
!%! ! Done.  delete temp storage and return.
!%! call destroy_obs(obs)
!%! call destroy_obs(prev_obs)
!%! 
!%! end subroutine select_gps_by_height
!%! !---------------------------------------------------------------------

!#! !---------------------------------------------------------------------
!#! subroutine change_obs(this_obs)
!#! 
!#! ! This code is intentionally commented out.  It is here as an example
!#! ! of how to make a systematic change to observations in an obs_seq file.
!#! ! To use it, remove all instances of !#! in the code, here and above.
!#! ! Adapt as needed.  Use with care so you don't run it on obs_sequences
!#! ! you don't want to modify in this way.
!#! 
!#! ! The code below is working code.  It changes the variance of the selected
!#! ! observation types.
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
!#! integer            :: this_obs_type
!#! real(r8)           :: oldvar, newvar
!#! logical, save      :: first_call = .true.
!#! 
!#! 
!#! if (first_call) then
!#!    call error_handler(E_MSG, 'obs_sequence_tool', 'special: using the change_obs routine to alter variances')
!#!    first_call = .false.
!#! endif
!#! 
!#! call get_obs_def(this_obs, this_obs_def)
!#! this_obs_type = get_obs_def_type_of_obs(this_obs_def)
!#! 
!#! ! do not alter identity obs
!#! if (this_obs_type < 0) return
!#! 
!#! oldvar = get_obs_def_error_variance(this_obs_def)
!#! 
!#! ! Set the new variance here based on the type
!#! select case (this_obs_type)
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
!#! !print *, 'type, original variance, new variance = ', this_obs_type, oldvar, newvar
!#!
!#! if (newvar /= oldvar) then
!#!    call set_obs_def_error_variance(this_obs_def, newvar)
!#!    call set_obs_def(this_obs, this_obs_def)
!#! endif
!#! 
!#! end subroutine change_obs
!#! !---------------------------------------------------------------------

!---------------------------------------------------------------------
end program obs_sequence_tool

