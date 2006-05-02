! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, 2006,
! Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program merge_obs_seq

use        types_mod, only : r8
use    utilities_mod, only : timestamp, register_module, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, logfileunit
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                             get_num_obs, init_obs_sequence, get_first_obs, &
                             write_obs_seq, set_copy_meta_data, get_obs_def, &
                             set_obs_def, append_obs_to_seq, get_next_obs, &
                             insert_obs_in_seq, init_obs, assignment(=), &
                             static_init_obs_sequence, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence

! use obs_def_gps_mod

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq1, seq2
type(obs_type)          :: obs, next_obs, new_obs
logical                 :: is_there_one, is_this_last
integer                 :: num_copies, num_qc
integer                 :: size_seq1, size_seq2
integer                 :: num_inserted, iunit, io
character(len = 129)    :: msgstring

!----------------------------------------------------------------
! Namelist input with default values

character(len = 129) :: filename_seq1 = "obs_seq.out",    &
                        filename_seq2 = "obs_seq.final",  &
                        filename_out  = 'filter_ics'

namelist /merge_obs_seq_nml/ filename_seq1, filename_seq2, filename_out

!----------------------------------------------------------------
! Start of the routine.
! This routine basically opens the second observation sequence file
! and 'inserts' itself into the first observation sequence file.
! Each observation in the second file is independently inserted 
! into the first file. This takes time ... Might want to throw
! in some checking to see if the observations in the second file 
! are AFTER the last time of the first file ... appending would
! be a lot faster -- less searching and comparing dates.
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

! Read 2nd obs seq first, and find out its num_obs
call read_obs_seq(filename_seq2, 0, 0, 0, seq2)
size_seq2 = get_num_obs(seq2)

! Read 1st obs seq, expand its size for insertion
call read_obs_seq(filename_seq1, 0, 0, size_seq2, seq1)
size_seq1 = get_num_obs(seq1)

! Compare metadata between the observation sequences.
! If it is compatible, keep going, if not ... it will terminate here.
 call compare_metadata(seq1, seq2)

num_copies  = get_num_copies(seq1)   ! By this time, num_copies, num_qc are
num_qc      = get_num_qc(seq1)       ! the same for both seq1 and seq2

call error_handler(E_MSG,'merge_obs_seq', &
    'Good news - observation sequence files compatible ...', source,revision,revdate)

! Initialize individual observation variables 
call init_obs(     obs, num_copies, num_qc)
call init_obs( new_obs, num_copies, num_qc)
call init_obs(next_obs, num_copies, num_qc)

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
   call insert_obs_in_seq(seq1, new_obs)
   num_inserted = num_inserted + 1
   
   call get_next_obs(seq2, obs, next_obs, is_this_last)
   ObsLoop : do while ( .not. is_this_last)

      if ( mod(num_inserted, 1000) == 0 ) then
         print*, 'inserted number ',num_inserted,' of ',size_seq2
      endif

      obs = next_obs
      new_obs = obs
      call insert_obs_in_seq(seq1, new_obs)
      num_inserted = num_inserted + 1
      call get_next_obs(seq2, obs, next_obs, is_this_last)
   enddo ObsLoop

else
   write(msgstring,*)'no first observation in ',trim(adjustl(filename_seq2))
   call error_handler(E_MSG,'merge_obs_seq', msgstring, source, revision, revdate)
endif

print*, 'Number of obs to be inserted :          ', size_seq2
print*, 'Number of inserted obs:                 ', num_inserted
print*, 'Total number of obs in the new seq file:', num_inserted+size_seq1

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
! 'observation types' header.  
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
