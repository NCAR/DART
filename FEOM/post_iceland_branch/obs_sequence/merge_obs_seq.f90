! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, 2006,
! Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program merge_obs_seq

use        types_mod, only : r8
use    utilities_mod, only : timestamp, register_module, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, logfileunit
use time_manager_mod, only : time_type, operator(>), print_time
use      obs_def_mod, only : obs_def_type, get_obs_def_time
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq, &
                             init_obs, init_obs_sequence, static_init_obs_sequence, &
                             read_obs_seq_header, read_obs_seq, assignment(=), &
                             get_num_obs, get_first_obs, get_last_obs, get_next_obs, &
                             insert_obs_in_seq, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence, get_obs_def, &
                             append_obs_to_seq

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq1, seq2
type(obs_type)          :: obs, next_obs, new_obs
type(obs_def_type)      :: obs_def
type(time_type)         :: last_time, this_time
logical                 :: is_there_one, is_this_last
integer                 :: size_seq1, num_copies1, num_qc1
integer                 :: size_seq2, num_copies2, num_qc2
integer                 :: size_seq,  num_copies,  num_qc
integer                 :: add_copies = 0, add_qc = 0
integer                 :: num_inserted, iunit, io
integer                 :: max_num_obs, file_id
character(len = 129)    :: read_format
logical                 :: pre_I_format
character(len = 129)    :: msgstring

!----------------------------------------------------------------
! Namelist input with default values

character(len = 129) :: filename_seq1 = "obs_seq.one",  &
                        filename_seq2 = "obs_seq.two",  &
                        filename_out  = 'obs_seq.merged'

namelist /merge_obs_seq_nml/ filename_seq1, filename_seq2, filename_out

!----------------------------------------------------------------
! Start of the routine.
! This routine basically opens the second observation sequence file
! and 'inserts' itself into the first observation sequence file.
! Each observation in the second file is independently inserted 
! or appended into the first file.
! Inserting takes time, appending is much faster when possible.
!----------------------------------------------------------------

call merge_obs_seq_modules_used()

! Read the namelist entry
call find_namelist_in_file("input.nml", "merge_obs_seq_nml", iunit)
read(iunit, nml = merge_obs_seq_nml, iostat = io)
call check_namelist_read(iunit, io, "merge_obs_seq_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'merge_obs_seq','merge_obs_seq_nml values are',' ',' ',' ')
write(logfileunit, nml=merge_obs_seq_nml)
write(     *     , nml=merge_obs_seq_nml)

! Read header information for the sequences to see if we need
! to accomodate additional copies or qc values from sequence two.

call read_obs_seq_header(filename_seq1, num_copies1, num_qc1, size_seq1, &
   max_num_obs, file_id, read_format, pre_I_format, close_the_file = .true.)

call read_obs_seq_header(filename_seq2, num_copies2, num_qc2, size_seq2, &
   max_num_obs, file_id, read_format, pre_I_format, close_the_file = .true.)

add_copies = 0
add_qc     = 0

if (num_copies2 > num_copies1) add_copies = num_copies2 - num_copies1
if (    num_qc2 >     num_qc1) add_qc     =     num_qc2 -     num_qc1

! Read 1st obs seq, expand its size for insertion
call read_obs_seq(filename_seq1, add_copies, add_qc, size_seq2, seq1)
size_seq    = get_num_obs(seq1)      ! does not include size_seq2
num_copies  = get_num_copies(seq1)   ! already includes add_copies
num_qc      = get_num_qc(seq1)       ! already includes add_qc

! Read 2nd obs seq
call read_obs_seq(filename_seq2, 0, 0, 0, seq2)

! Compare metadata between the observation sequences.
! If it is compatible, keep going, if not ... it will terminate here.
call compare_metadata(seq1, seq2)

call error_handler(E_MSG,'merge_obs_seq', &
    'Good news - observation sequence files compatible ...', source,revision,revdate)

! Initialize individual observation variables 
call init_obs(     obs, num_copies, num_qc)
call init_obs( new_obs, num_copies, num_qc)
call init_obs(next_obs, num_copies, num_qc)

! Getting the time of the last observation in the first sequence to see
! if we can append instead of insert. Appending is MUCH faster.
is_there_one = get_last_obs(seq1, obs)
if ( is_there_one ) then
   call get_obs_def(obs, obs_def)
   last_time = get_obs_def_time(obs_def)
else
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

   new_obs = obs    ! using interface procedure

   call        get_obs_def(new_obs, obs_def)
   this_time = get_obs_def_time(obs_def)

   if ( this_time > last_time ) then
      call append_obs_to_seq(seq1, new_obs)
      last_time = this_time
   else
      call insert_obs_in_seq(seq1, new_obs)
   endif

   num_inserted = num_inserted + 1
   
   call get_next_obs(seq2, obs, next_obs, is_this_last)
   ObsLoop : do while ( .not. is_this_last)

      if ( mod(num_inserted, 1000) == 0 ) then
         print*, 'inserted number ',num_inserted,' of ',size_seq2
      endif

      obs     = next_obs
      new_obs = obs

      ! If this obs is later than the last obs, we can append ... much faster
      call        get_obs_def(obs, obs_def)
      this_time = get_obs_def_time(obs_def)

      if ( this_time > last_time ) then
         call append_obs_to_seq(seq1, new_obs)
         last_time = this_time
      else
         call insert_obs_in_seq(seq1, new_obs)
      endif

      num_inserted = num_inserted + 1
      call get_next_obs(seq2, obs, next_obs, is_this_last)
   enddo ObsLoop

else
   write(msgstring,*)'no first observation in ',trim(adjustl(filename_seq2))
   call error_handler(E_MSG,'merge_obs_seq', msgstring, source, revision, revdate)
endif

print*, 'Number of obs in first file   :          ', size_seq
print*, 'Number of obs to be  inserted :          ', size_seq2
print*, 'Number of obs really inserted :          ', num_inserted
print*, 'Target number of obs in the new seq file:', size_seq + num_inserted
print*, 'Actual number of obs in the new seq file:', get_num_obs(seq1)

call write_obs_seq(seq1, filename_out)

! Time to clean up

call destroy_obs_sequence(seq1)
call destroy_obs_sequence(seq2)
call destroy_obs(     obs)
call destroy_obs( new_obs)
call destroy_obs(next_obs)

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

num_copies = num_copies1
num_qc     = num_qc1

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
