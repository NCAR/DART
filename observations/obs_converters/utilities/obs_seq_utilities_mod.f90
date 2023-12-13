! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> utility routines that work on observation sequences
!>

module obs_seq_utilities_mod


use        types_mod, only : r8, MISSING_R8, MISSING_I
use    utilities_mod, only : E_MSG, E_ERR, error_handler
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             get_obs_def_time, get_obs_def_location,           &
                             get_obs_kind, get_obs_def_error_variance,         &
                             set_obs_def_key
use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq, &
                             set_obs_values, set_qc, set_obs_def, get_qc,    &
                             get_obs_values, get_obs_def
use time_manager_mod, only : time_type, operator(>=), set_time, get_time
use     location_mod, only : set_location, location_type, get_location, query_location

use        types_mod, only : r8, metadatalength
use    utilities_mod, only : register_module, initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, finalize_utilities
use         sort_mod, only : index_sort
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name
use time_manager_mod, only : time_type, operator(>), print_time, set_time,     &
                             print_date, set_calendar_type, operator(==),      &
                             operator(/=), get_calendar_type, NO_CALENDAR,     &
                             operator(-)
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=),                          &
                             init_obs_sequence, static_init_obs_sequence,      &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_next_obs, get_obs_def,         &
                             insert_obs_in_seq, get_num_copies, get_num_qc,    &
                             get_copy_meta_data, get_qc_meta_data,             &
                             set_copy_meta_data, set_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence,                &
                             get_num_key_range, get_obs_key

use netcdf

implicit none
private

public :: print_obs_seq,           &
          validate_obs_seq_time,   &
          print_metadata

!>@todo FIXME there is no documentation for this module

! module global storage
character(len=NF90_MAX_NAME) :: missing_name = ''
character(512) :: msgstring, msgstring1, msgstring2

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

contains

!-----------------------------------------------------------------------
!> prints out a quick table of obs types and counts, overall start and
!> stop times, and metadata strings and counts.
!> you can get more info by running the obs_diag program.
!>
!> @param seq observation sequence of interest 
!> @param filename the file that was the origin of the observation sequence.
!>                 Used for messages.

subroutine print_obs_seq(seq, filename)

type(obs_sequence_type), intent(in) :: seq
character(len=*),        intent(in) :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq
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
! size_seq = get_num_key_range(seq)     !current size of seq

size_seq = get_num_obs(seq)
if (size_seq == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'obs_sort',msgstring)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(next_obs, get_num_copies(seq), get_num_qc(seq))

! blank line
call error_handler(E_MSG,'',' ')

write(msgstring,*) 'Processing sequence file ', trim(filename)
call error_handler(E_MSG,'',msgstring)

call print_metadata(seq, filename)

! Start to process obs from seq

is_there_one = get_first_obs(seq, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in ',trim(filename)
   call error_handler(E_MSG,'obs_sort', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')

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

   call get_next_obs(seq, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def), '  Last timestamp: ')
      call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')
   endif

enddo ObsLoop


write(msgstring, *) 'Number of obs processed  :          ', size_seq
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


!-----------------------------------------------------------------------
!> ensures the observation sequence is temporally ascending when 
!> traversed as a linked list.
!>
!> @param seq observation sequence of interest 
!> @param filename the file that was the origin of the observation sequence.
!>                 Used for messages.

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
   call error_handler(E_MSG,'obs_sort:validate',msgstring)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(next_obs, get_num_copies(seq), get_num_qc(seq))

obs_count = 0

! Start to process obs from seq

is_there_one = get_first_obs(seq, obs)

! we already tested for 0 obs above, so there should be a first obs here.
if ( .not. is_there_one )  then
   write(msgstring,*)'no first obs in sequence ' // trim(filename)
   call error_handler(E_ERR,'obs_sort:validate', &
                      msgstring, source, revision, revdate)
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
      call print_date(last_time, '      calendar date: ')
      call print_time(this_time, '     next timestamp: ')
      call print_date(this_time, '      calendar date: ')

      key = get_obs_key(obs)
      write(msgstring1,*)'obs number ', key, ' has earlier time than previous obs'
      write(msgstring2,*)'observations must be in increasing time order, file ' // trim(filename)
      call error_handler(E_ERR,'obs_sort:validate', msgstring2, &
                         source, revision, revdate, &
                         text2=msgstring1)
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
   call error_handler(E_MSG,'obs_sort:validate', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'obs_sort:validate', msgstring, &
                         source, revision, revdate, &
                         text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'obs_sort:validate', msgstring, &
                         source, revision, revdate, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time


!-----------------------------------------------------------------------
!> print the (trimmed) metadata strings
!>
!> @param seq observation sequence of interest 
!> @param filename the file that was the origin of the observation sequence.
!>                 Used for messages.

subroutine print_metadata(seq, filename)

type(obs_sequence_type),    intent(in) :: seq
character(len=*), optional, intent(in) :: filename

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(    seq)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring,*)' illegal copy or obs count in file '//trim(filename)
   call error_handler(E_ERR, 'obs_sort', msgstring, &
                      source, revision, revdate)
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

!-----------------------------------------------------------------------

end module obs_seq_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
