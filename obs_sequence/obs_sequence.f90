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
use utilities_mod, only : open_file
use time_manager_mod, only : time_type, set_time, operator(<=)

private

public obs_sequence_type, init_output_obs_sequence, init_input_obs_sequence, &
   get_global_meta_data, get_num_obs_copies, get_copy_meta_data, &
   get_next_obs_set, output_obs_set 

type obs_sequence_type
   integer :: file_id
   character(len = 128) :: file_meta_data
   integer :: num_copies
   character(len = 129), pointer :: copy_meta_data(:)
   type(time_type) :: time
end type obs_sequence_type

contains

!===============================================================================

function init_output_obs_sequence(file_name, file_meta_data, num_copies_in, copy_meta_data)
!------------------------------------------------------------------------------
!
! Returns an obs_sequence type that is initialized and ready for output. The num_copies 
! is optional with a default of 1 and may have a value of 0 for an obs_def file with
! no associated observations.

implicit none

type(obs_sequence_type) :: init_output_obs_sequence
character(len = *), intent(in) :: file_name, file_meta_data
integer, optional, intent(in) :: num_copies_in
character(len = *), intent(in), optional :: copy_meta_data(:)

integer :: file_id, num_copies, i

! Get the actual number of copies
num_copies = 1
if(present(num_copies_in)) num_copies = num_copies_in
init_output_obs_sequence%num_copies = num_copies

! Should later verify that the appropriate amount of copy_meta_data is available

! Open the output file
file_id = open_file(file_name)
init_output_obs_sequence%file_id = file_id

! Write the file meta_data
init_output_obs_sequence%file_meta_data = file_meta_data
write(file_id, *) file_meta_data

! Write the number of copies of data
write(file_id, *) num_copies

! Allocate space for copy meta data
allocate(init_output_obs_sequence%copy_meta_data(num_copies))

! Write the copy meta_data
do i = 1, num_copies
   init_output_obs_sequence%copy_meta_data(i) = copy_meta_data(i)
   write(file_id, *) copy_meta_data(i) 
end do

! Initialize the time to 0
init_output_obs_sequence%time = set_time(0, 0)

end function init_output_obs_sequence




function init_input_obs_sequence(file_name)
!------------------------------------------------------------------------------
!
! Returns an obs_sequence type that is initialized and ready for reading input. 

implicit none

type(obs_sequence_type) :: init_input_obs_sequence
character(len = *), intent(in) :: file_name

integer :: unit, i, num_copies

! Open the input file
unit = open_file(file_name)
init_input_obs_sequence%file_id = unit

! Read the file meta_data
read(unit, *) init_input_obs_sequence%file_meta_data

! Read the number of copies
read(unit, *) num_copies
init_input_obs_sequence%num_copies = num_copies

! Allocate space for copy meta data
allocate(init_input_obs_sequence%copy_meta_data(num_copies))

! Read the copy meta_data
do i = 1, num_copies
   read(unit, *) init_input_obs_sequence%copy_meta_data(i)
end do

! Initialize the time to 0
init_input_obs_sequence%time = set_time(0, 0)

end function init_input_obs_sequence



function get_global_meta_data(sequence)
!--------------------------------------------------------------------------------
!
! Returns global meta_data string for this sequence

implicit none

character(len = *) :: get_global_meta_data
type(obs_sequence_type), intent(in) :: sequence

get_global_meta_data = sequence%file_meta_data

end function get_global_meta_data




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

character(len = *) :: get_copy_meta_data
type(obs_sequence_type) :: sequence
integer, optional, intent(in) :: index_in

integer :: index

! Get the proper index
index = 1
if(present(index_in)) index = index_in
! Should do some error checks on index

get_copy_meta_data = sequence%copy_meta_data(index)

end function get_copy_meta_data



function get_next_obs_set(sequence)
!----------------------------------------------------------------------------------
!
! Returns the next obs_set in this sequence. Many additonal embelishments to this 
! interface will be needed for other assimilation applications.

implicit none

type(obs_set_type) :: get_next_obs_set
type(obs_sequence_type), intent(in) :: sequence

get_next_obs_set = read_obs_set(sequence%file_id)

end function get_next_obs_set



subroutine output_obs_set(sequence, obs_set)
!-----------------------------------------------------------------------------------
!
! Outputs the obs_set to the sequence. Makes error checks for time ordering. Need to
! add fields to make time checking straightforward here.

implicit none

type(obs_sequence_type), intent(inout) :: sequence
type(obs_set_type), intent(in) :: obs_set

type(time_type) :: time

! Make sure the times are ascending
time = get_obs_set_time(obs_set)
if(time <= sequence%time) then
   write(*, *) 'Error: obs_set being added to sequence is not in time order: output_obs_set'
   stop
endif

call write_obs_set(sequence%file_id, obs_set)

end subroutine output_obs_set





end module obs_sequence_mod
