! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Prints out a quick summary number of the total error in the mean
!> and spread.  This should be run on an obs_seq file which has been
!> through both perfect_model_obs and filter, so it has both a copy
!> of the 'truth' as well as 'ensemble mean' and 'ensemble spread'.
!> You can get more info by running the obs_diag program 

program obs_total_error

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
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs, &
                             get_obs_def_location, read_obs_def,               &
                             set_obs_def_time
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs
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
                             delete_seq_head, delete_seq_tail, get_obs_values, &
                             get_num_key_range, get_obs_key, get_qc,           &
                             copy_partial_obs, get_next_obs_from_key,          &
                             get_obs_def, set_obs_def

implicit none

character(len=*), parameter :: source = 'obs_total_error.f90'

type(obs_sequence_type) :: seq_in
integer                 :: size_seq_in
integer                 :: num_copies_in, num_qc_in
integer                 :: iunit, io, ifile
integer                 :: max_num_obs, file_id
character(len = 129)    :: read_format, filename_in
logical                 :: pre_I_format, cal
character(len = 256)    :: msgstring, msgstring1, msgstring2

!----------------------------------------------------------------
! Namelist input with default values


character(len = 160) :: obs_sequence_name = ''
character(len = 160) :: obs_sequence_list = ''

character(len=32)    :: calendar    = 'Gregorian'


namelist /obs_total_error_nml/ &
         obs_sequence_name, obs_sequence_list, calendar

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_total_error_nml", iunit)
read(iunit, nml = obs_total_error_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_total_error_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_total_error_nml)
if (do_nml_term()) write(     *     , nml=obs_total_error_nml)

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! set a logial to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! Check the user input for sanity
if ((obs_sequence_name /= '') .and. (obs_sequence_list /= '')) then
   write(msgstring1,*)'specify "obs_sequence_name" or "obs_sequence_list"'
   write(msgstring2,*)'set other to an empty string ... i.e. ""'
   call error_handler(E_ERR, 'obs_assim_count', msgstring1, source, text2=msgstring2)
endif

! if you add anything to the namelist, you can process it here.

! end of namelist processing and setup

ifile = 0
ObsFileLoop : do      ! until we run out of names
!-----------------------------------------------------------------------

   ifile = ifile + 1

   ! Determine the next input filename ...

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

call read_obs_seq_header(filename_in, num_copies_in, num_qc_in, &
      size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)
   
if (max_num_obs == 0) then
      write(msgstring,*) 'No obs in input sequence file ', trim(filename_in)
   call error_handler(E_ERR,'obs_total_error',msgstring)
endif
   
write(msgstring, *) 'Starting to process input sequence file: '
write(msgstring1,*)  trim(filename_in)
call error_handler(E_MSG,'obs_total_error',msgstring, &
                      text2=msgstring1)
   
call read_obs_seq(filename_in, 0, 0, 0, seq_in)
   
! sanity check - ensure the linked list times are in increasing time order
call validate_obs_seq_time(seq_in, filename_in)
   
! the computation is done here.
call print_obs_seq_info(seq_in, filename_in)
   
   ! clean up
   
   call destroy_obs_sequence(seq_in)

enddo ObsFileLoop

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_total_error')
call static_init_obs_sequence()

end subroutine setup


!---------------------------------------------------------------------
subroutine shutdown()

call finalize_utilities('obs_total_error')

end subroutine shutdown


!---------------------------------------------------------------------
subroutine print_obs_seq_info(seq_in, filename)

! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

! and this one counts up, if there is a 'posterior ensemble mean' copy,
! how many are missing_r8 and how many are not.   it could also count
! up the 'DART quality control' settings?  start with the latter for now.

type(obs_sequence_type), intent(in) :: seq_in
character(len=*), intent(in)        :: filename

type(obs_type)     :: obs, next_obs
type(obs_def_type) :: this_obs_def
logical            :: is_there_one, is_this_last
integer            :: size_seq_in
integer            :: i
integer            :: this_obs_type
integer            :: type_count(max_defined_types_of_obs), identity_count, qc_count(0:7), qcindex
integer            :: trueindex, meanindex, spreadindex, erritems
real(r8)           :: qcval(1), copyval(1), truth, mean, spread
real(r8)           :: tminusm, spreadsq
logical            :: can_do_error


! Initialize counters and sums
type_count(:) = 0
identity_count = 0
qc_count(:) = 0
tminusm = 0.0_r8
spreadsq = 0.0_r8
erritems = 0

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_total_error',msgstring)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq_in), get_num_qc(seq_in))
call init_obs(next_obs, get_num_copies(seq_in), get_num_qc(seq_in))

! find the dart qc copy, if there is one
qcindex = get_dartqc_index(seq_in, filename)

! find the true_state, ensemble_mean, ensemble_spread if they are there
trueindex = get_dartcopy_index(seq_in, filename, 'truth')
meanindex = get_dartcopy_index(seq_in, filename, 'posterior ensemble mean')
spreadindex = get_dartcopy_index(seq_in, filename, 'posterior ensemble spread')

can_do_error = all( (/trueindex, meanindex, spreadindex/) >= 0 )

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
   call error_handler(E_MSG,'obs_total_error', msgstring)
endif

! process it here
is_this_last = .false.

! blank line
call error_handler(E_MSG, '', ' ')

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
! does not work with NO_CALENDAR
if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_type = get_obs_def_type_of_obs(this_obs_def)
   if (this_obs_type < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_type) = type_count(this_obs_type) + 1
   endif
!   print *, 'obs type index = ', this_obs_type
!   if(this_obs_type > 0)print *, 'obs name = ', get_name_for_type_of_obs(this_obs_type)
   if (qcindex > 0) then
      call get_qc(obs, qcval, qcindex)
      qc_count(nint(qcval(1))) = qc_count(nint(qcval(1))) + 1
   endif

   if (can_do_error) then
      call get_obs_values(obs, copyval, trueindex)
      truth = copyval(1)
      call get_obs_values(obs, copyval, meanindex)
      mean = copyval(1)
      call get_obs_values(obs, copyval, spreadindex)
      spread = copyval(1)

      tminusm = tminusm + (truth - mean)**2.0
      spreadsq = spreadsq + spread**2.0
      erritems = erritems + 1
   endif

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def), '  Last timestamp: ')
      if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')
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
if (qcindex > 0) then
   call error_handler(E_MSG, '', ' ')
   write(msgstring, *) 'DART QC results: '
   call error_handler(E_MSG, '', msgstring)
   do i=0, 7
      if (qc_count(i) > 0) then
         write(msgstring, '(a16,2(i8))') 'DART QC value', i, &
                                        qc_count(i)
         call error_handler(E_MSG, '', msgstring)
      endif
   enddo
   write(msgstring, *) 'Total obs: ', sum(qc_count(:))
   call error_handler(E_MSG, '', msgstring)
endif

call error_handler(E_MSG, '', ' ')
if (can_do_error) then
   write(msgstring, *) 'Error results: '
   call error_handler(E_MSG, '', msgstring)
   write(msgstring, *) 'Total  error: ', sqrt(tminusm)
   call error_handler(E_MSG, '', msgstring)
   write(msgstring, *) 'Total spread: ', sqrt(spreadsq)
   call error_handler(E_MSG, '', msgstring)
else
   write(msgstring, *) 'No error results printed because one or more of "truth",'
   call error_handler(E_MSG, '', msgstring)
   write(msgstring, *) '"posterior ensemble mean", and/or "posterior ensemble spread" are missing'
   call error_handler(E_MSG, '', msgstring)
endif

! another blank line
call error_handler(E_MSG, '', ' ')

! Time to clean up

call destroy_obs(     obs)
call destroy_obs(next_obs)

end subroutine print_obs_seq_info


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
   call error_handler(E_MSG,'obs_total_error:validate',msgstring)
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
   call error_handler(E_ERR,'obs_total_error:validate', msgstring, source)
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
      call error_handler(E_ERR,'obs_total_error:validate', msgstring2, &
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
   call error_handler(E_MSG,'obs_total_error:validate', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'obs_total_error:validate', msgstring, &
                         source, text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'obs_total_error:validate', msgstring, &
                         source, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time


!---------------------------------------------------------------------
subroutine print_metadata(seq, fname)

!
! print out the metadata strings, trimmed
!

type(obs_sequence_type), intent(in) :: seq
character(len=*) :: fname

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str
character(len=255) :: msgstring3

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(    seq)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring3,*)' illegal copy or obs count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_total_error', msgstring3, source)
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

end subroutine print_metadata


!---------------------------------------------------------------------
function get_dartqc_index(seq, fname)

!
! return the index number of the dart qc copy (-1 if none)
!

type(obs_sequence_type), intent(in) :: seq
character(len=*) :: fname
integer :: get_dartqc_index

integer :: num_qc, i
character(len=metadatalength) :: str
character(len=255) :: msgstring3

num_qc = get_num_qc(seq)

if ( num_qc < 0 ) then
   write(msgstring3,*)' illegal qc metadata count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_total_error', msgstring3, source)
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

function get_dartcopy_index(seq, fname, copyname)

!
! return the index number of the given dart copy (-1 if none)
!

type(obs_sequence_type), intent(in) :: seq
character(len=*) :: fname, copyname
integer :: get_dartcopy_index

integer :: num_copies, i
character(len=metadatalength) :: str
character(len=255) :: msgstring3

num_copies = get_num_copies(seq)

if ( num_copies < 0 ) then
   write(msgstring3,*)' illegal copy metadata count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_total_error', msgstring3, source)
endif

MetaData : do i=1, num_copies
   str = get_copy_meta_data(seq,i)

   if (str == copyname) then
      get_dartcopy_index = i
      return 
   endif

enddo MetaData

get_dartcopy_index = -1

end function get_dartcopy_index


!---------------------------------------------------------------------
end program obs_total_error

