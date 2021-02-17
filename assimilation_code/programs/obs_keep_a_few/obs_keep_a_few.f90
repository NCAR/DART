! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> open an obs_seq file and copy over the first N of each obs type to
!> the output file. the value of N is namelist-settable.  intended
!> to subset a large obs_seq file for testing or other purposes but
!> still keep examples of each type of obs from the input.

program obs_keep_a_few

use        types_mod, only : r8, missing_r8, metadatalength
use    utilities_mod, only : initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             open_file, close_file, finalize_utilities
use     location_mod, only : location_type, get_location, set_location,        &
                             LocationName, read_location, operator(/=),        &
                             write_location
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs, &
                             get_obs_def_location, read_obs_def,               &
                             set_obs_def_time
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs, &
                             get_index_for_type_of_obs, read_type_of_obs_table,  &
                             get_num_types_of_obs
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

character(len=*), parameter :: source = 'obs_keep_a_few.f90'

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs_in, next_obs_in
type(obs_type)          :: obs_out, prev_obs_out
logical                 :: is_this_last
integer                 :: size_seq_in, size_seq_out
integer                 :: num_copies_in, num_qc_in
integer                 :: num_inserted, iunit, io, j
integer                 :: max_num_obs, file_id
integer                 :: num_rejected_badqc, num_rejected_diffqc
integer                 :: num_rejected_other
character(len=129)      :: read_format
logical                 :: pre_I_format, cal
character(len=512)      :: msgstring, msgstring1, msgstring2, msgstring3
type(obs_def_type)      :: this_obs_def

integer, allocatable    :: n_this_type(:)
integer                 :: this_type

character(len=metadatalength) :: meta_data

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 5000

! lazy, pick big number.  make it bigger if too small.
integer, parameter :: max_obs_input_types = 500

!----------------------------------------------------------------
! Namelist input with default values


character(len=256) :: filename_in = ''
character(len=256) :: filename_out = ''

integer            :: max_count_per_type = 10
integer            :: max_total_count    = -1

logical            :: print_only    = .false.
character(len=32)  :: calendar      = 'Gregorian'


namelist /obs_keep_a_few_nml/ &
         filename_in, filename_out, &
         max_count_per_type, max_total_count, &
         print_only, calendar

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_keep_a_few_nml", iunit)
read(iunit, nml = obs_keep_a_few_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_keep_a_few_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_keep_a_few_nml)
if (do_nml_term()) write(     *     , nml=obs_keep_a_few_nml)

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! set a logial to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! if you add anything to the namelist, you can process it here.

! end of namelist processing and setup

! make space for the counts.  0 is for all identity obs.
allocate(n_this_type(0:get_num_types_of_obs()))
n_this_type(:) = 0

! single pass algorithm (unlike other obs tools).

call read_obs_seq_header(filename_in, num_copies_in, num_qc_in, &
   size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
   close_the_file = .true.)

if (max_num_obs == 0) then
   write(msgstring,*) 'No obs in input sequence file ', trim(filename_in)
   call error_handler(E_ERR,'obs_keep_a_few',msgstring)
endif

write(msgstring, *) 'Starting to process input sequence file: '
write(msgstring1,*)  trim(filename_in)
call error_handler(E_MSG,'obs_keep_a_few',msgstring, &
                   text2=msgstring1)

call read_obs_seq(filename_in, 0, 0, 0, seq_in)

! sanity check - ensure the linked list times are in increasing time order
call validate_obs_seq_time(seq_in, filename_in)

! output is same size (or less) than input, generally.
! if this program is going to dup obs, account for it here.
size_seq_out = max_num_obs

! blank line, start of actually creating output file
call error_handler(E_MSG,' ',' ')

! Initialize individual observation variables
call init_obs(     obs_in,  num_copies_in, num_qc_in)
call init_obs(next_obs_in,  num_copies_in, num_qc_in)
call init_obs(     obs_out, num_copies_in, num_qc_in)
call init_obs(prev_obs_out, num_copies_in, num_qc_in)
   
! create the output sequence here 
call init_obs_sequence(seq_out, num_copies_in, num_qc_in, size_seq_out)
do j=1, num_copies_in
   meta_data = get_copy_meta_data(seq_in, j)
   call set_copy_meta_data(seq_out, j, meta_data)
enddo
do j=1, num_qc_in
   meta_data = get_qc_meta_data(seq_in, j)
   call set_qc_meta_data(seq_out, j, meta_data)
enddo

! is this needed?
if (print_only) call print_obs_seq(seq_in, filename_in)

!-------------------------------------------------------------
! Start to insert obs from sequence_in into sequence_out
!
! NOTE: insert_obs_in_seq CHANGES the obs passed in.
!       Must pass a copy of incoming obs to insert_obs_in_seq.
!--------------------------------------------------------------
num_inserted = 0
num_rejected_badqc = 0
num_rejected_diffqc = 0
num_rejected_other = 0

if ( get_first_obs(seq_in, obs_in) )  then

   is_this_last = .false.
   next_obs_in = obs_in

   ObsLoop : do while ( .not. is_this_last )

      obs_in = next_obs_in

      ! obs_out will be modified when it is inserted in the output sequence
      ! so we have to make a copy of obs_in before modifiying it.
      obs_out = obs_in

      ! count up how many of this type you already have
      ! and skip it if you've got enough.

      call get_obs_def(obs_out, this_obs_def)

      this_type = get_obs_def_type_of_obs(this_obs_def)

      if (this_type < 0) this_type = 0   ! identity obs

      if (n_this_type(this_type) < max_count_per_type .or. max_count_per_type < 0) then
 
         ! copy to output obs_seq and increment the count for this type
         n_this_type(this_type) = n_this_type(this_type) + 1
   
         if (num_inserted > 0) then
            call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)
         else
            call insert_obs_in_seq(seq_out, obs_out)
         endif
   
         prev_obs_out = obs_out  ! update position in seq for next insert
         num_inserted = num_inserted + 1
   
         if (print_every > 0) then
            if (mod(num_inserted,print_every) == 0) then
               print*, 'inserted number ',num_inserted,' of ',size_seq_out
            endif
         endif

      endif

      if (max_total_count > 0 .and. num_inserted >= max_total_count) exit ObsLoop

      call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

   enddo ObsLoop

else
   write(msgstring, *)'no first observation in ',trim(filename_in)
   call error_handler(E_MSG,'obs_keep_a_few', msgstring)
endif

if (.not. print_only) then
   print*, '---------  Obs seqs '
   print*, 'Number of obs input sequence    :         ', size_seq_in
   print*, 'Number of obs copied to output  :         ', num_inserted
   print*, '---------------------------------------------------------'
endif


write(msgstring, *) 'Starting to process output sequence file ', &
                            trim(filename_out)
call error_handler(E_MSG,'obs_keep_a_few',msgstring)

print*, 'Number of obs in the output seq file :', get_num_key_range(seq_out)

call print_obs_seq(seq_out, filename_out)
if (.not. print_only) then
   call write_obs_seq(seq_out, filename_out)
else
   write(msgstring,*) 'Output sequence file not created; print_only in namelist is .true.'
   call error_handler(E_MSG,'', msgstring)
endif

! clean up

call destroy_obs_sequence(seq_in)
call destroy_obs_sequence(seq_out)
call destroy_obs(     obs_in )
call destroy_obs(next_obs_in )
call destroy_obs(     obs_out)
!call destroy_obs(prev_obs_out)  ! copy of something already deleted
deallocate(n_this_type)

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_keep_a_few')
call static_init_obs_sequence()

end subroutine setup


!---------------------------------------------------------------------
subroutine shutdown()

call finalize_utilities('obs_keep_a_few')

end subroutine shutdown


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
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_keep_a_few',msgstring)
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
   call error_handler(E_MSG,'obs_keep_a_few', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
! does not work with NO_CALENDAR
if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')

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
integer                 :: size_seq, obs_count
integer                 :: key
type(time_type)         :: last_time, this_time


! make sure there are obs left to process before going on.
size_seq = get_num_obs(seq) 
if (size_seq == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_keep_a_few:validate',msgstring)
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
   call error_handler(E_ERR,'obs_keep_a_few:validate', msgstring, source)
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
      call error_handler(E_ERR,'obs_keep_a_few:validate', msgstring2, &
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
   call error_handler(E_MSG,'obs_keep_a_few:validate', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'obs_keep_a_few:validate', msgstring, &
                         source, text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'obs_keep_a_few:validate', msgstring, &
                         source, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time


!---------------------------------------------------------------------
subroutine print_metadata(seq, fname)

!
! print out the metadata strings, trimmed
!

type(obs_sequence_type),    intent(in) :: seq
character(len=*), optional, intent(in) :: fname

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(    seq)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring3,*)' illegal copy or obs count in file '//trim(fname)
   call error_handler(E_ERR, 'obs_keep_a_few', msgstring3, source)
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
end program obs_keep_a_few

