! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program merge_obs_seq

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8
use    utilities_mod, only : timestamp, register_module, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, logfileunit
use time_manager_mod, only : time_type, operator(>), print_time, set_time
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq, &
                             init_obs, init_obs_sequence, static_init_obs_sequence, &
                             read_obs_seq_header, read_obs_seq, assignment(=), &
                             get_num_obs, get_first_obs, get_last_obs, get_next_obs, &
                             insert_obs_in_seq, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence, delete_seq_head, &
                             delete_seq_tail, get_num_key_range, set_copy_meta_data

implicit none

! <next few lines under version control, do not edit>
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs, prev_obs, next_obs, new_obs
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

!----------------------------------------------------------------
! Namelist input with default values

! max_num_input_files : maximum number of input sequence files to be merged
integer  :: max_num_input_files, num_input_files
parameter( max_num_input_files = 50) 

character(len = 129) :: filename_seq(max_num_input_files)
character(len = 129) :: filename_out  = 'obs_seq.merged'
logical              :: process_file(max_num_input_files)

! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days    = -1
integer  :: first_obs_seconds = -1
integer  :: last_obs_days     = -1
integer  :: last_obs_seconds  = -1
type(time_type) :: first_obs_time, last_obs_time


namelist /merge_obs_seq_nml/ num_input_files, filename_seq, filename_out, &
         first_obs_days, first_obs_seconds, last_obs_days, last_obs_seconds

!----------------------------------------------------------------
! Start of the routine.
! This routine basically opens the second observation sequence file
! and 'inserts' itself into the first observation sequence file.
! Each observation in the second file is independently inserted 
! or appended into the first file.
! Inserting takes time, appending is much faster when possible.
!----------------------------------------------------------------

call merge_obs_seq_modules_used()

! Initialize input obs_seq filenames
do i = 1, max_num_input_files
   filename_seq(i) = ""
enddo

! Read the namelist entry
call find_namelist_in_file("input.nml", "merge_obs_seq_nml", iunit)
read(iunit, nml = merge_obs_seq_nml, iostat = io)
call check_namelist_read(iunit, io, "merge_obs_seq_nml")

if(num_input_files .gt. max_num_input_files) then
  call error_handler(E_ERR,'merge_obs_seq', &
     'num_input_files > max_num_input_files. change max_num_input_files in source file', &
     source,revision,revdate)
endif

! Record the namelist values used for the run ...
call error_handler(E_MSG,'merge_obs_seq','merge_obs_seq_nml values are',' ',' ',' ')
write(logfileunit, nml=merge_obs_seq_nml)
write(     *     , nml=merge_obs_seq_nml)

! Read header information for the sequences to see if we need
! to accomodate additional copies or qc values from subsequent sequences.
! Also, calculate how many observations to be added to the first sequence.
num_copies_out   = 0
num_qc_out       = 0
size_seq_out     = 0

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
      call error_handler(E_ERR,'merge_obs_seq', 'first time cannot be later than last time', &
                         source,revision,revdate)
   endif
endif

! TWO PASS algorithm; open each file, trim it if requested, and count
! the number of actual observations.  then the output file can be
! created with the correct size, and as observations are put into it
! they'll be sorted, and unused obs will be removed.

! pass 1:

first_seq = -1
do i = 1, num_input_files

   if ( len(filename_seq(i)) .eq. 0 .or. filename_seq(i) .eq. "" ) then
      call error_handler(E_ERR,'merge_obs_seq', &
         'num_input_files and filename_seq mismatch',source,revision,revdate)
   endif

   ! count up the number of observations we are going to eventually have.
   ! if all the observations in a file are not part of the linked list, the
   ! output number of observations might be much smaller than the total size in 
   ! the header.  it is slower, but go ahead and read in the entire sequence
   ! and count up the real number of obs - trim_seq will do the count even if
   ! it is not trimming in time.  this allows us to create an empty obs_seq
   ! output file of exactly the right size.

   call read_obs_seq_header(filename_seq(i), num_copies_in, num_qc_in, size_seq_in, &
      max_num_obs, file_id, read_format, pre_I_format, close_the_file = .true.)
   
   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)
   call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time,   &
                 filename_seq(i), .true., remaining_obs_count)
   call destroy_obs_sequence(seq_in)
   if (remaining_obs_count == 0) then
      process_file(i) = .false.
      cycle
   else
      process_file(i) = .true.
      size_seq_in = remaining_obs_count
   endif

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
   msgstring = 'All input files are empty or all obs outside the first/last times'
   call error_handler(E_ERR,'merge_obs_seq',msgstring,source,revision,revdate)
endif

! pass 2:

! Initialize individual observation variables 
call init_obs(     obs, num_copies_out, num_qc_out)
call init_obs( new_obs, num_copies_out, num_qc_out)
call init_obs(next_obs, num_copies_out, num_qc_out)
call init_obs(prev_obs, num_copies_out, num_qc_out)

total_num_inserted = 0

! Read obs seq to be added, and insert obs from it to the output seq
first_seq = -1
do i = 1, num_input_files

   if (.not. process_file(i)) cycle
 
   write(msgstring,*) 'Starting to process input sequence file ', trim(filename_seq(i))
   call error_handler(E_MSG,'merge_obs_seq',msgstring,source,revision,revdate)

   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)

   ! If you get here, there better be observations in this file which
   ! are going to be used (the process_file flag wouldn't be set otherwise.)
   call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time,   &
                 filename_seq(i), .false., remaining_obs_count)

   ! This would be an error at this point.
   if(remaining_obs_count == 0) then
      call destroy_obs_sequence(seq_in) 
      write(msgstring, *) 'Internal error trying to process file ', trim(filename_seq(i))
      call error_handler(E_ERR,'merge_obs_seq',msgstring,source,revision,revdate)
   endif

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

   !-------------------------------------------------------------
   ! Start to insert obs from sequence_in into sequence_out
   !
   ! NOTE: insert_obs_in_seq CHANGES the obs passed in.
   !       Must pass a copy of incoming obs to insert_obs_in_seq.
   !--------------------------------------------------------------
   num_inserted = 0
   is_there_one = get_first_obs(seq_in, obs)

   if ( is_there_one )  then

      new_obs      = obs           ! obs records position in seq_out

      call insert_obs_in_seq(seq_out, new_obs)  ! new_obs linked list info changes

      prev_obs     = new_obs       ! records new position in seq_in
      num_inserted = num_inserted + 1
   
      call get_next_obs(seq_in, obs, next_obs, is_this_last)
      ObsLoop : do while ( .not. is_this_last)

         if (mod(num_inserted,1000) == 0) then
            print*, 'inserted number ',num_inserted,' of ',size_seq_in
         endif

         obs     = next_obs   ! essentially records position in seq_out
         new_obs = obs        ! will be modified w/ position in seq_in

         ! Since the stride through the observation sequence file is always 
         ! guaranteed to be in temporally-ascending order, we can use the
         ! 'previous' observation as the starting point to search for the
         ! correct insertion point.  This speeds up the insert code a lot.

         call insert_obs_in_seq(seq_out, new_obs, prev_obs)

         prev_obs     = new_obs    ! update position in seq_in for next insert
         num_inserted = num_inserted + 1

         call get_next_obs(seq_in, obs, next_obs, is_this_last)

      enddo ObsLoop

      total_num_inserted = total_num_inserted + num_inserted

   else
      write(msgstring,*)'no first observation in ',trim(adjustl(filename_seq(i)))
      call error_handler(E_MSG,'merge_obs_seq', msgstring, source, revision, revdate)
   endif

   print*, '--------------  Obs seq file # :          ', i
   print*, 'Number of obs in previous seq  :          ', size_seq_out
   print*, 'Number of obs to be  inserted  :          ', size_seq_in
   print*, 'Number of obs really inserted  :          ', num_inserted
   print*, '---------------------------------------------------------'

   call destroy_obs_sequence(seq_in)

enddo

print*, 'Total number of obs  inserted  :          ', total_num_inserted
print*, 'Actual number of obs in the new seq file :', get_num_key_range(seq_out)

call write_obs_seq(seq_out, filename_out)

! Time to clean up

call destroy_obs_sequence(seq_out)
call destroy_obs(     obs)
call destroy_obs( new_obs)
call destroy_obs(next_obs)
call destroy_obs(prev_obs)

call timestamp(source,revision,revdate,'end')

contains


!=====================================================================


subroutine merge_obs_seq_modules_used()

! Initialize modules used that require it
call initialize_utilities('merge_obs_seq')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

end subroutine merge_obs_seq_modules_used


  subroutine compare_metadata(seq1, seq2, fname1, fname2)
!
! This subroutine compares the metadata for two different observation
! sequences and terminates the program if they are not conformable.
! In order to be merged, the two observation sequences must have the same
! number of qc values, the same number of copies ... 
!
! The messages might be a bit misleading 'warning', 'warning', 'error' ...
!

 type(obs_sequence_type), intent(IN) :: seq1, seq2
 character(len=*), optional :: fname1, fname2

integer :: num_copies1, num_qc1
integer :: num_copies2, num_qc2
integer :: num_copies , num_qc, i
character(len=129) :: str1, str2
character(len=255) :: msgstring1, msgstring2

num_qc1     = get_num_qc(    seq1)
num_copies1 = get_num_copies(seq1)

num_qc2     = get_num_qc(    seq2)
num_copies2 = get_num_copies(seq2)

num_copies  = num_copies1
num_qc      = num_qc1

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
   call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
   num_copies = -1
endif
if ( num_qc1 /= num_qc2 ) then
   write(msgstring2,*)'Different different numbers of QCs found: ', &
                      num_qc1, ' vs ', num_qc2
   call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
   num_qc = -1
endif
if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'merge_obs_seq', msgstring1, source, revision, revdate)
endif

MetaDataLoop : do i=1, num_copies
   str1 = trim(adjustl(get_copy_meta_data(seq1,i)))
   str2 = trim(adjustl(get_copy_meta_data(seq2,i)))

   if( str1 == str2 ) then
      write(msgstring2,*)'metadata ',trim(adjustl(str1)), ' in both.'
      call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
   else
      write(msgstring2,*)'metadata value mismatch. seq1: ', trim(adjustl(str1))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
      write(msgstring2,*)'metadata value mismatch. seq2: ', trim(adjustl(str2))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
      call error_handler(E_ERR, 'merge_obs_seq', msgstring1, source, revision, revdate)
   endif
enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str1 = trim(adjustl(get_qc_meta_data(seq1,i)))
   str2 = trim(adjustl(get_qc_meta_data(seq2,i)))

   if( str1 == str2 ) then
      write(msgstring2,*)'qc metadata ', trim(adjustl(str1)), ' in both.'
      call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
   else
      write(msgstring2,*)'qc metadata value mismatch. seq1: ', trim(adjustl(str1))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
      write(msgstring2,*)'qc metadata value mismatch. seq2: ', trim(adjustl(str2))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring2, source, revision, revdate)
      call error_handler(E_ERR, 'merge_obs_seq', msgstring1, source, revision, revdate)
   endif
enddo QCMetaData

end subroutine compare_metadata

! pass in an already opened sequence and a start/end time.  this routine
! really trims the observations out of the sequence, and returns a count
! of how many remain.
subroutine trim_seq(seq, trim_first, first_time, trim_last, last_time, seqfilename, &
                    print_msg, remaining_obs_count)
 type(obs_sequence_type), intent(inout) :: seq
 logical, intent(in)                    :: trim_first, trim_last
 type(time_type), intent(in)            :: first_time, last_time
 character(len = *), intent(in)         :: seqfilename
 logical, intent(in)                    :: print_msg
 integer, intent(out)                   :: remaining_obs_count

   ! Need to find first obs with appropriate time, delete all earlier ones
   if(trim_first) then
      call delete_seq_head(first_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are before first_obs_days:first_obs_seconds'
            call error_handler(E_MSG,'merge_obs_seq',msgstring,source,revision,revdate)
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
                        ' are after last_obs_days:last_obs_seconds'
            call error_handler(E_MSG,'merge_obs_seq',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif
   
   remaining_obs_count = get_num_key_range(seq)

end subroutine trim_seq

end program merge_obs_seq
