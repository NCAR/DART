module obs_sequence_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

use types_mod
use obs_set_mod, only : obs_set_type, read_obs_set, write_obs_set, get_obs_set_time
use set_def_list_mod, only : set_def_list_type, &
   read_set_def_list, write_set_def_list
use utilities_mod, only : open_file
use time_manager_mod, only : time_type, set_time, operator(<=)

private

public obs_sequence_type, init_obs_sequence, &
   get_num_obs_copies, get_copy_meta_data, &
   get_obs_set, add_obs_set, set_def_list, read_obs_sequence, &
   write_obs_sequence
	

type obs_sequence_type 
   private
   integer :: num_copies
   character(len = 129), pointer :: copy_meta_data(:)
   integer :: max_obs_sets, num_obs_sets
   type(set_def_list_type) :: def_list
   type(obs_set_type), pointer :: obs_sets(:)
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

type(obs_sequence_type) :: init_obs_sequence
integer, intent(in) :: max_obs_sets
integer, optional, intent(in) :: num_copies_in
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




function get_num_obs_copies(sequence)
!--------------------------------------------------------------------------------
!
! Returns global meta_data string for this sequence

implicit none

integer :: get_num_obs_copies
type(obs_sequence_type), intent(in) :: sequence

get_num_obs_copies = sequence%num_copies

end function get_num_obs_copies




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

type(obs_set_type) :: get_obs_set
type(obs_sequence_type), intent(in) :: sequence
integer, intent(in) :: index

get_obs_set = sequence%obs_sets(index)

end function get_obs_set



subroutine add_obs_set(sequence, obs_set)
!-----------------------------------------------------------------------------------
!

implicit none

type(obs_sequence_type), intent(inout) :: sequence
type(obs_set_type), intent(in) :: obs_set

type(time_type) :: time

! Make sure the times are ascending
time = get_obs_set_time(sequence%obs_sets(sequence%num_obs_sets))
if(get_obs_set_time(obs_set) <= time) then
   write(*, *) 'Error: obs_set being added to sequence is not in time order: output_obs_set'
   stop
endif

! Make sure there's room in the fixed storage implementation
if(sequence%num_obs_sets == sequence%max_obs_sets) then
   write(*, *) 'Error: No more room in sequence in add_obs_set'
   stop
endif

! Increment counters and add
sequence%num_obs_sets = sequence%num_obs_sets + 1
sequence%obs_sets(sequence%num_obs_sets) = obs_set

end subroutine add_obs_set




subroutine set_def_list(sequence, def_list)
!------------------------------------------------------------------
!
! Set the def_list associated with the sequence. May just want to
! make this a pointer transaction in the future or work directly
! with sequences def_list to avoid copy of huge amounts of data?

implicit none

type(obs_sequence_type), intent(inout) :: sequence
type(set_def_list_type), intent(in) :: def_list

sequence%def_list = def_list

end subroutine set_def_list



subroutine write_obs_sequence(file_id, sequence)
!------------------------------------------------------------------
!

implicit none

integer, intent(in) :: file_id
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
do i = 1, sequence%num_copies
   call write_obs_set(file_id, sequence%obs_sets(i))
end do

end subroutine write_obs_sequence



function read_obs_sequence(file_id)
!----------------------------------------------------------
!

implicit none

type(obs_sequence_type) :: read_obs_sequence 
integer, intent(in) :: file_id

integer :: i, num_copies, num_obs_sets

! Read number of copies of data 
read(file_id, *) num_copies

! Read the number of observations in set
read(file_id, *) num_obs_sets

! Initialize the sequence
read_obs_sequence = init_obs_sequence(num_obs_sets, num_copies)

! Read the per copy metadata
do i = 1, num_copies
   read(file_id, *) read_obs_sequence%copy_meta_data(i)
end do

! Read out the observation definition list
read_obs_sequence%def_list =  read_set_def_list(file_id)

! Read the obs sets 
do i = 1, read_obs_sequence%num_copies
   read_obs_sequence%obs_sets(i) = read_obs_set(file_id)
end do

end function read_obs_sequence



end module obs_sequence_mod
