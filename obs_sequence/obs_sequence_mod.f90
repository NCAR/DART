module obs_sequence_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

use types_mod
use utilities_mod,    only : open_file
use time_manager_mod, only : time_type, set_time, operator(<=)

use obs_set_mod, only : obs_set_type, read_obs_set, write_obs_set, get_obs_set_time, &
   obs_set_copy, get_num_obs, read_obs_set_time, obs_set_time_copy, &
   os_get_obs_def_index    => get_obs_def_index, &
   os_get_obs_values       => get_obs_values, &
   os_set_obs_values       => set_obs_values, &
   os_set_single_obs_value => set_single_obs_value, &
   os_inc_num_obs_copies   => inc_num_obs_copies

use set_def_list_mod, only : set_def_list_type, &
   read_set_def_list, write_set_def_list, set_def_list_copy, &
   sd_get_expected_obs     => get_expected_obs, &
   sd_get_diag_obs_err_cov => get_diag_obs_err_cov, &
   sd_get_num_close_states => get_num_close_states, &
   sd_get_close_states     => get_close_states

private

public obs_sequence_type, init_obs_sequence, &
   get_num_obs_copies, get_copy_meta_data, &
   get_obs_set, add_obs_set, associate_def_list, read_obs_sequence, &
   write_obs_sequence, obs_sequence_copy, get_num_obs_sets, &
   get_obs_sequence_time, get_num_obs_in_set, &
   get_expected_obs, get_diag_obs_err_cov, get_obs_values, &
   inc_num_obs_copies, set_obs_values, get_num_close_states, &
   get_close_states, read_obs_sequence_def, obs_sequence_def_copy, &
   set_single_obs_value, get_obs_def_index
	

type obs_sequence_type 
   private
   integer                       :: num_copies
   character(len = 129), pointer :: copy_meta_data(:)
   integer                       :: max_obs_sets, num_obs_sets
   type(set_def_list_type)       :: def_list
   type(obs_set_type),   pointer :: obs_sets(:)
end type obs_sequence_type

contains

!===============================================================================

function init_obs_sequence(max_obs_sets, num_copies_in, copy_meta_data)
!------------------------------------------------------------------------------
!
! Returns an obs_sequence type that is initialized and ready for output. The num_copies 
! is optional with a default of 1 and may have a value of 0 for an obs_def file with
! no associated observations. For now the storage size is fixed with a size
! limit of max_obs_sets. Need to make this dynamic and allow for some sort of
! streaming input and output with windows eventually.

implicit none

type(obs_sequence_type)        :: init_obs_sequence
integer,            intent(in) :: max_obs_sets
integer, optional,  intent(in) :: num_copies_in
character(len = *), intent(in), optional :: copy_meta_data(:)

integer :: file_id, num_copies, i

! Get the actual number of copies
num_copies = 1
if(present(num_copies_in)) num_copies = num_copies_in
init_obs_sequence%num_copies = num_copies

! Should later verify that the appropriate amount of copy_meta_data is available

! Allocate space for copy meta data
allocate(init_obs_sequence%copy_meta_data(num_copies))

! Set the copy meta_data
if(present(copy_meta_data)) then
   do i = 1, num_copies
      init_obs_sequence%copy_meta_data(i) = copy_meta_data(i)
   end do
endif

init_obs_sequence%max_obs_sets = max_obs_sets
init_obs_sequence%num_obs_sets = 0

! Initialize the storage
allocate(init_obs_sequence%obs_sets(max_obs_sets))

end function init_obs_sequence




subroutine obs_sequence_copy(seq_out, seq_in)
!---------------------------------------------------------------------------------
!
! Comprehensive copy for obs_sequence_type

implicit none

type(obs_sequence_type), intent(out) :: seq_out
type(obs_sequence_type), intent(in)  :: seq_in

integer :: i

seq_out%num_copies   = seq_in%num_copies
seq_out%max_obs_sets = seq_in%max_obs_sets
seq_out%num_obs_sets = seq_in%num_obs_sets

! Allocate and copy the meta_data

allocate(seq_out%copy_meta_data(seq_in%num_copies), &
         seq_out%obs_sets(seq_in%max_obs_sets))
seq_out%copy_meta_data = seq_in%copy_meta_data

! Copy the set_def_list
call set_def_list_copy(seq_out%def_list, seq_in%def_list)

! Copy the obs sets one by one
do i = 1, seq_in%num_obs_sets
   call obs_set_copy(seq_out%obs_sets(i), seq_in%obs_sets(i))
end do

end subroutine obs_sequence_copy



subroutine obs_sequence_def_copy(seq_out, seq_in)
!---------------------------------------------------------------------------------
!
! Comprehensive copy for obs_sequence_type

implicit none

type(obs_sequence_type), intent(out) :: seq_out
type(obs_sequence_type), intent(in)  :: seq_in

integer :: i

! Copy just the definition part of the sequence setting copies to 0
seq_out%num_copies   = 0
seq_out%max_obs_sets = seq_in%max_obs_sets
seq_out%num_obs_sets = seq_in%num_obs_sets

! Allocate the copy metadata space to 0 and the obs space
allocate(seq_out%copy_meta_data(0), seq_out%obs_sets(seq_in%max_obs_sets))

! Copy the set_def_list
call set_def_list_copy(seq_out%def_list, seq_in%def_list)

! Copy the obs sets one by one, just the def info
do i = 1, seq_in%num_obs_sets
   call obs_set_time_copy(seq_out%obs_sets(i), seq_in%obs_sets(i))
end do

end subroutine obs_sequence_def_copy




subroutine inc_num_obs_copies(seq, inc, copy_meta_data)
!---------------------------------------------------------
!
! Increments the number of copies of the actual data associated
! with an obs sequence while retaining all the other information.
! This requires an inefficient re-allocation in the current
! fixed storage version.

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: inc
character(len=129),      intent(in)    :: copy_meta_data(inc)

character(len=129) :: temp_meta_data(inc + seq%num_copies)
integer :: old_num, new_num, i

! Check to make sure the increment is positive
if(inc < 0) then
   write(*, *) 'Error: increment must be positive in inc_num_obs_copies'
   stop
endif

! Get current num copies and adjust the metadata in sequnce
old_num = seq%num_copies
new_num = old_num + inc
seq%num_copies = new_num

! Need temporary storage for metadata
temp_meta_data(1:old_num) = seq%copy_meta_data
temp_meta_data(old_num + 1:new_num) = copy_meta_data

! Deallocate and reallocate (really want to avoid this for fragmentation
! possibilities?)
deallocate(seq%copy_meta_data)
allocate(seq%copy_meta_data(new_num))
seq%copy_meta_data = temp_meta_data

! Now need to loop through all the obs_sets and add extra space there
do i = 1, seq%num_obs_sets
   call os_inc_num_obs_copies(seq%obs_sets(i), inc)
end do

end subroutine inc_num_obs_copies



function get_num_obs_copies(sequence)
!--------------------------------------------------------------------------------
!
! Returns global meta_data string for this sequence

implicit none

integer :: get_num_obs_copies
type(obs_sequence_type), intent(in) :: sequence

get_num_obs_copies = sequence%num_copies

end function get_num_obs_copies




function get_num_obs_sets(sequence)
!---------------------------------------------------------------------------------
!
! Returns number of obs_sets in this sequence (want this to be dynamic eventually)

implicit none

integer :: get_num_obs_sets
type(obs_sequence_type) :: sequence

get_num_obs_sets = sequence%num_obs_sets

end function get_num_obs_sets



subroutine get_obs_values(seq, index, obs, copy_in)
!---------------------------------------------------
!
! Returns the observed values associated with the index set
! in the sequence.

implicit none

type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: index
real(r8),                intent(out) :: obs(:)
integer, optional,       intent(in)  :: copy_in

integer :: copy

copy = 1
if(present(copy_in)) copy = copy_in

call os_get_obs_values(seq%obs_sets(index), obs, copy)

end subroutine get_obs_values



function get_obs_def_index(seq, index)
!----------------------------------------------------
!
! Returns the definition index of the index-th observation 
! set in the sequence.

implicit none

integer                             :: get_obs_def_index
type(obs_sequence_type), intent(in) :: seq
integer,                 intent(in) :: index

get_obs_def_index = os_get_obs_def_index(seq%obs_sets(index))

end function get_obs_def_index




subroutine set_obs_values(seq, index, obs, copy_in)
!----------------------------------------------------
!
! Sets the index obs_set in the sequence to the obs input.

implicit none

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: index
real(r8),                intent(in)    :: obs(:)
integer, optional,       intent(in)    :: copy_in

integer :: copy

! Default value for copy is 1
copy = 1
if(present(copy_in)) copy = copy_in

call os_set_obs_values(seq%obs_sets(index), obs, copy)

end subroutine set_obs_values




subroutine set_single_obs_value(seq, index, num_obs, obs, copy_in)
!----------------------------------------------------
!
! Sets the num_obs obs of the index obs_set in the sequence to the obs input.

implicit none

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: index, num_obs
real(r8),                intent(in)    :: obs
integer, optional,       intent(in)    :: copy_in

integer :: copy

! Default value for copy is 1
copy = 1
if(present(copy_in)) copy = copy_in

call os_set_single_obs_value(seq%obs_sets(index), num_obs, obs, copy)

end subroutine set_single_obs_value




function get_copy_meta_data(sequence, index_in)
!---------------------------------------------------------------------------------
! Returns meta data for the index-th copy of obs. Default for index is 1

implicit none

character(len = 129) :: get_copy_meta_data
type(obs_sequence_type) :: sequence
integer, optional, intent(in) :: index_in

integer :: index

! Get the proper index
index = 1
if(present(index_in)) index = index_in
! Should do some error checks on index

get_copy_meta_data = sequence%copy_meta_data(index)

end function get_copy_meta_data



function get_obs_set(sequence, index)
!----------------------------------------------------------------------------------
!
! Returns the index obs_set in this sequence. Many additonal embelishments to this 
! interface will be needed for other assimilation applications.

implicit none

type(obs_set_type)                  :: get_obs_set
type(obs_sequence_type), intent(in) :: sequence
integer,                 intent(in) :: index

call obs_set_copy(get_obs_set, sequence%obs_sets(index))
!!!get_obs_set = sequence%obs_sets(index)

end function get_obs_set



subroutine get_diag_obs_err_cov(seq, index, cov)
!-------------------------------------------------------------
!
! Returns the diagonal part of observational error covariance
! for index observation in the sequence.

implicit none

type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: index
real(r8),                intent(out) :: cov(:)

type(obs_set_type) :: obs_set
integer :: def_index

! Get the set_def_list index for this obs_set
obs_set = get_obs_set(seq, index)
def_index = os_get_obs_def_index(obs_set)

call sd_get_diag_obs_err_cov(seq%def_list, def_index, cov)

end subroutine get_diag_obs_err_cov




subroutine get_num_close_states(seq, index, radius, num)
!------------------------------------------------------
! 
! Gets the number of state variables within distance radius
! of each of the observations in the index set in sequence.

implicit none

type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: index
real(r8),                intent(in)  :: radius
integer,                 intent(out) :: num(:)

type(obs_set_type) :: obs_set
integer :: def_index

! Get the set_def_list index for this obs_set
obs_set = get_obs_set(seq, index)
def_index = os_get_obs_def_index(obs_set)

call sd_get_num_close_states(seq%def_list, def_index, radius, num)

end subroutine get_num_close_states




subroutine get_close_states(seq, index, radius, num, indices, dist)
!------------------------------------------------------
!
! Gets the number and index of state variables within distance radius
! of each of the observations in the index set in sequence.

implicit none

type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: index
real(r8),                intent(in)  :: radius
integer,                 intent(out) :: num(:), indices(:, :)
real(r8),                intent(out) :: dist(:, :)

type(obs_set_type) :: obs_set
integer :: def_index

! Get the set_def_list index for this obs_set
obs_set = get_obs_set(seq, index)
def_index = os_get_obs_def_index(obs_set)

call sd_get_close_states(seq%def_list, def_index, radius, num, &
   indices, dist)

end subroutine get_close_states




subroutine get_expected_obs(sequence, index, state, obs, num)
!------------------------------------------------------
!
! Gets the expected value of observations for the obs_set_def
! corresponding to the index obs_set in the sequence with
! the state given. num allows selection of a single observation
! from the set.

implicit none

type(obs_sequence_type), intent(in)  :: sequence
integer,                 intent(in)  :: index
real(r8),                intent(in)  :: state(:)
real(r8),                intent(out) :: obs(:)
integer, optional,       intent(in)  :: num

type(obs_set_type) :: obs_set
integer :: def_index

! Get the set_def_list index for this obs_set
obs_set = get_obs_set(sequence, index)
def_index = os_get_obs_def_index(obs_set)

if(present(num)) then
   call sd_get_expected_obs(sequence%def_list, def_index, state, obs, num)
else
   call sd_get_expected_obs(sequence%def_list, def_index, state, obs)
endif

end subroutine get_expected_obs



function get_obs_sequence_time(sequence, index)
!-----------------------------------------------------
!
! Returns the time of the index-th obs_set in the sequence.

implicit none

type(time_type)                     :: get_obs_sequence_time
type(obs_sequence_type), intent(in) :: sequence
integer,                 intent(in) :: index

type(obs_set_type) :: obs_set

obs_set = get_obs_set(sequence, index)
get_obs_sequence_time = get_obs_set_time(obs_set)

end function get_obs_sequence_time



function get_num_obs_in_set(sequence, index)
!-------------------------------------------------
!
! Gets the number of obs in the index obs_set in the sequence

implicit none

integer                             :: get_num_obs_in_set
type(obs_sequence_type), intent(in) :: sequence
integer,                 intent(in) :: index

type(obs_set_type) :: obs_set

obs_set = get_obs_set(sequence, index)
get_num_obs_in_set = get_num_obs(obs_set)

end function get_num_obs_in_set



subroutine add_obs_set(sequence, obs_set)
!-----------------------------------------------------------------------------------
!

implicit none

type(obs_sequence_type), intent(inout) :: sequence
type(obs_set_type),      intent(in)    :: obs_set

type(time_type) :: time

! Make sure the times are ascending
time = get_obs_set_time(sequence%obs_sets(sequence%num_obs_sets))
if(get_obs_set_time(obs_set) <= time .and. sequence%num_obs_sets > 0) then
   write(*, *) 'Error: obs_set being added to sequence is not in time order: add_obs_set'
   stop
endif

! Make sure there's room in the fixed storage implementation
if(sequence%num_obs_sets == sequence%max_obs_sets) then
   write(*, *) 'Error: No more room in sequence in add_obs_set'
   stop
endif

! Increment counters and add
! Is Default copy sufficient
sequence%num_obs_sets = sequence%num_obs_sets + 1
call obs_set_copy(sequence%obs_sets(sequence%num_obs_sets), obs_set)
!!!sequence%obs_sets(sequence%num_obs_sets) = obs_set

end subroutine add_obs_set




subroutine associate_def_list(sequence, def_list)
!------------------------------------------------------------------
!
! Set the def_list associated with the sequence. May just want to
! make this a pointer transaction in the future or work directly
! with sequences def_list to avoid copy of huge amounts of data?

implicit none

type(obs_sequence_type), intent(inout) :: sequence
type(set_def_list_type), intent(in)    :: def_list

! Is default copy sufficient?
call set_def_list_copy(sequence%def_list, def_list)
!!!sequence%def_list = def_list

end subroutine associate_def_list



subroutine write_obs_sequence(file_id, sequence)
!------------------------------------------------------------------
!

implicit none

integer,                 intent(in) :: file_id
type(obs_sequence_type), intent(in) :: sequence

integer :: i

! Write number of copies of data 
write(file_id, *) sequence%num_copies

! Write number of obs sets
write(file_id, *) sequence%num_obs_sets

! Write the per copy metadata
do i = 1, sequence%num_copies
   write(file_id, *) sequence%copy_meta_data(i)
end do

! Write out the observation definition list
call write_set_def_list(file_id, sequence%def_list)

! Write out the obs sets 
do i = 1, sequence%num_obs_sets
   call write_obs_set(file_id, sequence%obs_sets(i))
end do

end subroutine write_obs_sequence



function read_obs_sequence(file_id)
!----------------------------------------------------------
!

implicit none

type(obs_sequence_type) :: read_obs_sequence 
integer,     intent(in) :: file_id

integer :: i, num_copies, num_obs_sets

! Read number of copies of data 
read(file_id, *) num_copies

! Read the number of observations in set
read(file_id, *) num_obs_sets

! Initialize the sequence
read_obs_sequence = init_obs_sequence(num_obs_sets, num_copies)
read_obs_sequence%num_obs_sets = num_obs_sets

! Read the per copy metadata
do i = 1, num_copies
   read(file_id, *) read_obs_sequence%copy_meta_data(i)
end do

! Read the observation definition list
call set_def_list_copy(read_obs_sequence%def_list, read_set_def_list(file_id))
!!!read_obs_sequence%def_list =  read_set_def_list(file_id)

! Read the obs sets 
do i = 1, read_obs_sequence%num_obs_sets
   call obs_set_copy(read_obs_sequence%obs_sets(i), read_obs_set(file_id))
!!!   read_obs_sequence%obs_sets(i) = read_obs_set(file_id)
end do

end function read_obs_sequence



function read_obs_sequence_def(file_id)
!----------------------------------------------------------
!

implicit none

type(obs_sequence_type) :: read_obs_sequence_def 
integer,     intent(in) :: file_id

integer :: i, num_copies, num_obs_sets

! Read number of copies of data 
read(file_id, *) num_copies

! Read the number of observations in set
read(file_id, *) num_obs_sets

! Initialize the sequence
read_obs_sequence_def = init_obs_sequence(num_obs_sets, num_copies)
read_obs_sequence_def%num_obs_sets = num_obs_sets

! Read the per copy metadata
do i = 1, num_copies
   read(file_id, *) read_obs_sequence_def%copy_meta_data(i)
end do

! Read the observation definition list
call set_def_list_copy(read_obs_sequence_def%def_list, &
   read_set_def_list(file_id))

! Read the obs sets definition portion but not the values
do i = 1, read_obs_sequence_def%num_obs_sets
   call obs_set_copy(read_obs_sequence_def%obs_sets(i), &
      read_obs_set_time(file_id))
end do

end function read_obs_sequence_def






end module obs_sequence_mod
