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
use time_manager_mod, only : time_type, operator(>), print_time
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq, &
                             init_obs, init_obs_sequence, static_init_obs_sequence, &
                             read_obs_seq_header, read_obs_seq, assignment(=), &
                             get_num_obs, get_first_obs, get_last_obs, get_next_obs, &
                             insert_obs_in_seq, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence

implicit none

! <next few lines under version control, do not edit>
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq1, seq2
type(obs_type)          :: obs, prev_obs, next_obs, new_obs
logical                 :: is_there_one, is_this_last
integer                 :: size_seq1, num_copies1, num_qc1
integer                 :: size_seq2, num_copies2, num_qc2
integer                 :: size_seq,  num_copies,  num_qc
integer                 :: add_size_seq = 0
integer                 :: add_copies = 0, add_qc = 0
integer                 :: num_inserted, iunit, io, i, total_num_inserted
integer                 :: max_num_obs, file_id
character(len = 129)    :: read_format
logical                 :: pre_I_format
character(len = 129)    :: msgstring

!----------------------------------------------------------------
! Namelist input with default values

! max_num_input_files : maximum number of input sequence files to be merged
integer  :: max_num_input_files, num_input_files
parameter( max_num_input_files = 50) 

character(len = 129) :: filename_seq(max_num_input_files)
character(len = 129) :: filename_out  = 'obs_seq.merged'

namelist /merge_obs_seq_nml/ num_input_files, filename_seq, filename_out

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
! Also, calculate how many observation to be added to the first sequence.
add_copies   = 0
add_qc       = 0
add_size_seq = 0

do i = 1, num_input_files

   if ( len(filename_seq(i)) .eq. 0 .or. filename_seq(i) .eq. "" ) then
      call error_handler(E_ERR,'merge_obs_seq', &
         'num_input_files and filename_seq mismatch',source,revision,revdate)
   endif

   if ( i .eq. 1) then
      ! for the header of the first seq
      call read_obs_seq_header(filename_seq(i), num_copies1, num_qc1, size_seq1, &
         max_num_obs, file_id, read_format, pre_I_format, close_the_file = .true.)
   else
      call read_obs_seq_header(filename_seq(i), num_copies2, num_qc2, size_seq2, &
         max_num_obs, file_id, read_format, pre_I_format, close_the_file = .true.)
      add_size_seq = add_size_seq + size_seq2
      if (num_copies2 > num_copies1) then 
         add_copies  = num_copies2 - num_copies1
         num_copies1 = num_copies2
      endif
      if (num_qc2 > num_qc1) then
         add_qc  = num_qc2 - num_qc1
         num_qc1 = num_qc2
      endif
   endif

enddo

! Read 1st obs seq, expand its size for insertion
call read_obs_seq(filename_seq(1), add_copies, add_qc, add_size_seq, seq1)
size_seq    = get_num_obs(seq1)      ! does not include size_seq2
num_copies  = get_num_copies(seq1)   ! already includes add_copies
num_qc      = get_num_qc(seq1)       ! already includes add_qc

! Initialize individual observation variables 
call init_obs(     obs, num_copies, num_qc)
call init_obs( new_obs, num_copies, num_qc)
call init_obs(next_obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

total_num_inserted = 0

! Read obs seq to be added, and insert obs from it to the first seq
do i = 2, num_input_files

   call read_obs_seq(filename_seq(i), 0, 0, 0, seq2)

   size_seq1 = get_num_obs(seq1)     !current size of seq1
   size_seq2 = get_num_obs(seq2)     !current size of seq2

   ! Compare metadata between the observation sequences.
   ! If it is compatible, keep going, if not ... it will terminate here.
   call compare_metadata(seq1, seq2)

   call error_handler(E_MSG,'merge_obs_seq', &
       'Good news - observation sequence files compatible ...', source,revision,revdate)

   ! Getting the time of the last observation in the first sequence to see
   ! if we can append instead of insert. Appending is MUCH faster.
   is_there_one = get_last_obs(seq1, obs)
   if (.not. is_there_one .and. size_seq1 >= 1) then
      call error_handler(E_ERR,'merge_obs_seq', &
       'BAD news - first obs_seq is neverending ...', source,revision,revdate) 
   endif

   !-------------------------------------------------------------
   ! Start to insert obs from sequence 2 into sequence 1
   !
   ! NOTE: insert_obs_in_seq CHANGES the obs passed in.
   !       Must pass a copy of incoming obs to insert_obs_in_seq.
   !--------------------------------------------------------------
   num_inserted = 0
   is_there_one = get_first_obs(seq2, obs)

   if ( is_there_one )  then

      new_obs      = obs           ! obs records position in seq2

      call insert_obs_in_seq(seq1, new_obs)  ! new_obs linked list info changes

      prev_obs     = new_obs       ! records new position in seq1
      num_inserted = num_inserted + 1
   
      call get_next_obs(seq2, obs, next_obs, is_this_last)
      ObsLoop : do while ( .not. is_this_last)

         if (mod(num_inserted,1000) == 0) then
            print*, 'inserted number ',num_inserted,' of ',size_seq2
         endif

         obs     = next_obs   ! essentially records position in seq2
         new_obs = obs        ! will be modified w/ position in seq1

         ! Since the stride through the observation sequence file is always 
         ! guaranteed to be in temporally-ascending order, we can use the
         ! 'previous' observation as the starting point to search for the
         ! correct insertion point. 

         call insert_obs_in_seq(seq1, new_obs, prev_obs)

         prev_obs     = new_obs    ! update position in seq1 for next insert
         num_inserted = num_inserted + 1

         call get_next_obs(seq2, obs, next_obs, is_this_last)

      enddo ObsLoop

      total_num_inserted = total_num_inserted + num_inserted

   else
      write(msgstring,*)'no first observation in ',trim(adjustl(filename_seq(i)))
      call error_handler(E_MSG,'merge_obs_seq', msgstring, source, revision, revdate)
   endif

   print*, '--------------  Obs seq file # :          ', i
   print*, 'Number of obs in previous seq  :          ', size_seq1
   print*, 'Number of obs to be  inserted  :          ', size_seq2
   print*, 'Number of obs really inserted  :          ', num_inserted
   print*, '---------------------------------------------------------'

   call destroy_obs_sequence(seq2)

enddo

print*, 'Number of obs in first seq     :          ', size_seq
print*, 'Total number of obs  inserted  :          ', total_num_inserted
print*, 'Target number of obs in the new seq file :', size_seq + total_num_inserted
print*, 'Actual number of obs in the new seq file :', get_num_obs(seq1)

call write_obs_seq(seq1, filename_out)

! Time to clean up

call destroy_obs_sequence(seq1)
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


  subroutine compare_metadata(seq1, seq2)
! subroutine compare_metadata(seq1, seq2)
!
! This subroutine compares the metadata for two different observation
! sequences and terminates the program if they are not conformable.
! In order to be merged, the two observation sequences must have the same
! number of qc values, the same number of copies ... 
!
! The messages might be a bit misleading 'warning', 'warning', 'error' ...
!

 type(obs_sequence_type), intent(IN) :: seq1, seq2

integer :: num_copies1, num_qc1
integer :: num_copies2, num_qc2
integer :: num_copies , num_qc, i
character(len=129) :: str1, str2
character(len=255) :: msgstring

num_qc1     = get_num_qc(    seq1)
size_seq1   = get_num_obs(   seq1)
num_copies1 = get_num_copies(seq1)

num_qc2     = get_num_qc(    seq2)
size_seq2   = get_num_obs(   seq2)
num_copies2 = get_num_copies(seq2)

num_copies  = num_copies1
num_qc      = num_qc1

if ( num_copies1 /= num_copies2 ) then
   write(msgstring,*)'The obs_sequences have incompatible numbers of copies', &
                      num_copies1, num_copies2 
   call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
   num_copies = -1
endif
if ( num_qc1 /= num_qc2 ) then
   write(msgstring,*)'The obs_sequences have incompatible numbers of QC metadata', &
                      num_qc1, num_qc2
   call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
   num_qc = -1
endif
if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'merge_obs_seq', &
        'obs_sequence files incompatible ... stopping.', source, revision, revdate)
endif

MetaDataLoop : do i=1, num_copies
   str1 = trim(adjustl(get_copy_meta_data(seq1,i)))
   str2 = trim(adjustl(get_copy_meta_data(seq2,i)))

   if( str1 == str2 ) then
      write(msgstring,*)'metadata ',trim(adjustl(str1)), ' in both.'
      call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
   else
      write(msgstring,*)'metadata seq1 ', trim(adjustl(str1))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
      write(msgstring,*)'metadata seq2 ', trim(adjustl(str2))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
      call error_handler(E_ERR, 'merge_obs_seq', &
        'obs_sequence files incompatible ... stopping.', source, revision, revdate)
   endif
enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str1 = trim(adjustl(get_qc_meta_data(seq1,i)))
   str2 = trim(adjustl(get_qc_meta_data(seq2,i)))

   if( str1 == str2 ) then
      write(msgstring,*)'qc metadata ', trim(adjustl(str1)), ' in both.'
      call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
   else
      write(msgstring,*)'qc metadata seq1 ', trim(adjustl(str1))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
      write(msgstring,*)'qc metadata seq2 ', trim(adjustl(str2))
      call error_handler(E_MSG, 'merge_obs_seq', msgstring, source, revision, revdate)
      call error_handler(E_ERR, 'merge_obs_seq', &
        'obs_sequence files incompatible ... stopping.', source, revision, revdate)
   endif
enddo QCMetaData

end subroutine compare_metadata


end program merge_obs_seq
