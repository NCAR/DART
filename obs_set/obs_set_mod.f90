! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_set_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

use        types_mod, only : r8
use set_def_list_mod, only : set_def_list_type, get_total_num_obs
use time_manager_mod, only : time_type, read_time, write_time, &
                             get_time, set_time

implicit none
private

public obs_set_type, init_obs_set, get_obs_set_time, get_obs_values,&
   set_obs_values, set_single_obs_value, set_obs_set_time, &
   get_single_obs_value, &
   contains_data, obs_value_missing, &
   read_obs_set, write_obs_set, obs_set_copy, get_num_obs, &
   get_obs_def_index, inc_num_obs_copies, read_obs_set_time, &
   obs_set_time_copy

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type obs_set_type
   private
! The two arrays will be allocated as (num_obs, num_copies)
   real(r8), pointer :: obs(:, :)
   logical, pointer :: missing(:, :)
   integer :: num_copies
   integer :: num_obs
   type(time_type) :: time
   integer :: def_index
end type obs_set_type


! NOTE: There are a variety of places where some no change lock mechanism
! needs to be implemented. Once obs_set_defs are in the table they should
! not have internals changed (except links). Once they are used by an 
! obs_set they should not be changed at all. Same thing with obs_sets,
! once they are in a sequence they should not be changed. Will need to
! implement lock semaphore in structures.


contains

!=======================================================


function init_obs_set(set_def_list, ind, num_copies_in)
!--------------------------------------------------------------------
!
! Creates and initializes storage for an obs_set associated
! with the obs_set_def index in obs_set_def_list and with num_copies of the
! observations associated (default is 1). 

type(obs_set_type) :: init_obs_set
type(set_def_list_type), intent(in) :: set_def_list
integer, intent(in) :: ind
integer, optional, intent(in) :: num_copies_in

integer :: num_copies, num_obs

! Begin by getting the number of copies
num_copies = 1
if(present(num_copies_in)) num_copies = num_copies_in
if(num_copies < 0) then
   write(*, *) 'Error: Negative num_copies in init_obs_set'
   stop
endif

! Set the number of copies
init_obs_set%num_copies = num_copies

! Assign the set definition
init_obs_set%def_index = ind

! Get the total number of obs in the set
num_obs = get_total_num_obs(set_def_list, ind)
init_obs_set%num_obs = num_obs

! Allocate storage for obs; need to verify that the 0 size allocation
! is legal F90
allocate(init_obs_set%obs(num_obs, num_copies), &
     init_obs_set%missing(num_obs, num_copies))

! Initialize all obs to present
init_obs_set%missing = .FALSE.

end function init_obs_set



subroutine obs_set_copy(set_out, set_in)
!-----------------------------------------------------------
!
! Comprehensive copy for obs_set
! TJH Jan 25, 2003 -- if we can reuse the same storage, we should

implicit none

! Copies with intent(out) for set_out lead to memory leak; 7 Oct. 2002
!!!type(obs_set_type), intent(out) :: set_out
type(obs_set_type), intent(inout) :: set_out
type(obs_set_type), intent(in) :: set_in

integer :: days, seconds

! Set the sizes
set_out%num_copies = set_in%num_copies
set_out%num_obs = set_in%num_obs

! The intel compiler doesn't do what I believe is default type copy!!!
! NEED TO CHECK FOR THIS THROUGHOUT, YUCK   11 Nov., 2002
call get_time(set_in%time, seconds, days)
set_out%time = set_time(seconds, days)

set_out%def_index = set_in%def_index

! Allocate storage for obs and missing
! INTEL 7.1 compiler with full obs checking dies on this for uninitialized
! HOWEVER, memory leaks here without this. Try leaving for now.
! TJH Jan 25, 2003 ... tried to fix it, have not checked.

! This is needed for the Intel 7.1 compiler.
! I am not sure if it hurts with the PG compiler.
! It does break with the Intel 8.0 compiler.

if( associated(set_out%obs) ) then
   if (size(set_out%obs) > 0 ) deallocate(set_out%obs)
else
   nullify(set_out%obs)
endif
if( associated(set_out%missing) ) then
   if (size(set_out%missing) > 0) deallocate(set_out%missing)
else
   nullify(set_out%missing)
endif

allocate(set_out%obs(set_in%num_obs, set_in%num_copies), &
     set_out%missing(set_in%num_obs, set_in%num_copies))

! Copy the values in obs and missing
set_out%obs     = set_in%obs
set_out%missing = set_in%missing

end subroutine obs_set_copy




subroutine obs_set_time_copy(set_out, set_in)
!-----------------------------------------------------------
!
! Copy just the definition (time and set_def portion) of an 
! obs set. This could really have been set as an intermediate
! object (the obs_set part). NEED TO WATCH TERMINOLOGY.

implicit none

type(obs_set_type), intent(out) :: set_out
type(obs_set_type), intent(in) :: set_in

! Set the sizes, num_copies is reset to 0
set_out%num_copies = 0
set_out%num_obs = set_in%num_obs
set_out%time = set_in%time
set_out%def_index = set_in%def_index

! Allocate storage for obs and missing
allocate(set_out%obs(set_in%num_obs, 0), &
     set_out%missing(set_in%num_obs, 0))

end subroutine obs_set_time_copy



subroutine inc_num_obs_copies(set, inc)
!------------------------------------------------------
!
! Increments the number of obs_copies and creates space.
! Current hard storage version requires ugly release
! and reallocation which may lead to fragmentation issues.

implicit none

type(obs_set_type), intent(inout) :: set
integer, intent(in) :: inc

logical :: temp_missing(set%num_obs, set%num_copies)
real(r8) :: temp_obs(set%num_obs, set%num_copies)
integer :: old_num, new_num

! Increment the number of copies
old_num = set%num_copies
new_num = old_num + inc
set%num_copies = new_num

! Copy to temporary storage for obs and missing
temp_missing = set%missing
temp_obs     = set%obs

! Deallocate and reallocate
allocate(set%obs(set%num_obs, new_num), &
     set%missing(set%num_obs, new_num))

! Copy in the existing data
set%missing(:, 1:old_num) = temp_missing
set%obs(:, 1:old_num) = temp_obs

! Set new missing to false
set%missing(:, old_num + 1 : new_num) = .FALSE.

end subroutine inc_num_obs_copies




function get_obs_def_index(obs_set)
!-------------------------------------------------------
!
! Returns the index of the obs_set_def in the set_def_list
! associated with this obs_set.

implicit none

integer :: get_obs_def_index
type(obs_set_type), intent(in) :: obs_set

get_obs_def_index = obs_set%def_index

end function get_obs_def_index





function get_obs_set_time(set)
!-------------------------------------------------------
!
! Returns the time associated with this observation set

implicit none

type(time_type) :: get_obs_set_time
type(obs_set_type), intent(in) :: set

get_obs_set_time = set%time

end function get_obs_set_time



function get_num_obs(set)
!-------------------------------------------------------
!
! Returns the number of obs in the set

implicit none

integer :: get_num_obs
type(obs_set_type), intent(in) :: set

get_num_obs = set%num_obs

end function get_num_obs




subroutine get_obs_values(set, obs, index_in)
!-----------------------------------------------------------------
!
! Returns the values of observations from the set. If the set has
! multiple observations per definition, the optional argument index
! selects which set to return.

implicit none

type(obs_set_type), intent(in) :: set
real(r8), intent(out) :: obs(:)
integer, optional, intent(in) :: index_in

integer :: ind

! Get the appropriate index
ind = 1
if(present(index_in)) ind = index_in
if((ind < 1) .or. (ind > set%num_copies)) then
   write(*, *) 'Error: Out of range index in get_obs_values'
   stop
endif

! Next make sure there's enough room in obs
if(size(obs) < set%num_obs) then
   write(*, *) 'Error: obs array too small in get_obs_values'
   stop
endif

!Copy the obs
obs = set%obs(:, ind)

end subroutine get_obs_values



subroutine set_obs_values(set, obs, copy_in)
!--------------------------------------------------------------------
!
! Sets the obs values; copy_in is optional with default 1.

implicit none

type(obs_set_type), intent(inout) :: set
real(r8), intent(in) :: obs(:)
integer, optional, intent(in) :: copy_in

integer :: copy

! Get the appropriate copycopy
copy = 1
if(present(copy_in)) copy = copy_in
if(copy < 1 .or. copy > set%num_copies) then
   write(*, *) 'Error: Out of range copy in set_obs_values'
   stop
endif

! Make sure the obs array is the right size
if(size(obs) /= set%num_obs) then
   write(*, *) 'Error: obs array wrong size in set_obs_values'
   stop
endif

!Copy the obs
set%obs(:, copy) = obs

end subroutine set_obs_values




subroutine set_single_obs_value(set, num_obs, obs, copy_in)
!--------------------------------------------------------------------
!
! Sets the obs value number num_obs;  copy_in is optional with default 1.

implicit none

type(obs_set_type), intent(inout) :: set
integer, intent(in) :: num_obs
real(r8), intent(in) :: obs
integer, optional, intent(in) :: copy_in

integer :: copy

! Get the appropriate copy
copy = 1
if(present(copy_in)) copy = copy_in 
if(copy < 1 .or. copy > set%num_copies) then
   write(*, *) 'Error(set_single_obs_value): Out of range "copy"'
   write(*, *) 'range is [1,',set%num_copies,'], "copy" is ',copy
   stop
endif

! Make sure the obs index is legal
if(num_obs > set%num_obs .or. num_obs < 1) then
   write(*, *) 'Error(set_single_obs_value): "num_obs" wrong size'
   write(*, *) 'num_obs is ',num_obs,' must be between [1,',set%num_obs,']'
   stop
endif

! Copy the obs
set%obs(num_obs, copy) = obs

end subroutine set_single_obs_value




subroutine set_obs_set_missing(set, missing, index_in)
!--------------------------------------------------------------------
!
! Sets the obs values; index is optional with default 1.

implicit none

type(obs_set_type), intent(inout) :: set
logical, intent(in) :: missing(:)
integer, optional, intent(in) :: index_in

integer :: ind

! Get the appropriate index
ind = 1
if(present(index_in)) ind = index_in
if((ind < 1) .or. (ind > set%num_copies)) then
   write(*, *) 'Error: Out of range index in init_obs_set'
   stop
endif

! Make sure the missing array is the right size
if(size(missing) /= set%num_obs) then
   write(*, *) 'Error: missing array wrong size in set_obs_set_missing'
   stop
endif

! Set the data missing
set%missing(:, ind) = missing

end subroutine set_obs_set_missing



subroutine set_obs_set_time(set, time)
!----------------------------------------------------
!
! Set the time for the obs_set

implicit none

type(obs_set_type), intent(inout) :: set
type(time_type), intent(in) :: time

set%time = time

end subroutine set_obs_set_time




function contains_data(set)
!--------------------------------------------------------------------
!
! Returns true if the number of associated data copies with the set is
! not 0.

implicit none

logical :: contains_data
type(obs_set_type), intent(in) :: set

contains_data = set%num_copies > 0

end function contains_data




subroutine obs_value_missing(set, missing, index_in)
!--------------------------------------------------------------------
!
! Returns true if the data associated with the copy is missing. Default
! for index_in is 1.

implicit none

type(obs_set_type), intent(in) :: set
logical, intent(out) :: missing(:)
integer, optional, intent(in) :: index_in

integer :: ind

! Get the appropriate index
ind = 1
if(present(index_in)) ind = index_in
if(ind < 1 .or. ind > set%num_copies) then
   write(*, *) 'Error: Out of range index in init_obs_set'
   stop
endif

! Next make sure there's enough room in obs
if(size(missing) /= set%num_obs) then
   write(*, *) 'Error: missing array too small in obs_value_missing'
   stop
endif

! Copy the missing data
missing = set%missing(:, ind)

end subroutine obs_value_missing




function read_obs_set(file_id)
!------------------------------------------------------------------------
!
! Reads an obs_set from a file

implicit none

type(obs_set_type) :: read_obs_set
integer, intent(in) :: file_id

character(len=5) :: header
integer :: num_obs, num_copies, i

! Read the header and verify 
read(file_id, *) header
if(header /= 'obset') then
   write(*, *) 'Error: Expected "obset" in header in read_obs_set' 
   stop
end if

! Read the obs_set def index
read(file_id, *) read_obs_set%def_index

! Read the number of obs and the number of copies
read(file_id, *) num_obs, num_copies
read_obs_set%num_obs = num_obs
read_obs_set%num_copies = num_copies

! Allocate space
allocate(read_obs_set%obs(num_obs, num_copies), &
     read_obs_set%missing(num_obs, num_copies))

! Read the data for each copy in turn
do i = 1, num_copies
   read(file_id, *) read_obs_set%obs(:, i)
end do

! Read the missing fields for each copy in turn
do i = 1, num_copies
   read(file_id, *) read_obs_set%missing(:, i)
end do

! Read the time 
read_obs_set%time = read_time(file_id)

end function read_obs_set



function read_obs_set_time(file_id)
!------------------------------------------------------------------------
!
! Reads an obs_set from a file but keep only the time info, acting as
! if 0 copies of the data were available.

implicit none

type(obs_set_type) :: read_obs_set_time
integer, intent(in) :: file_id

character(len=5) :: header
integer :: num_obs, num_copies, i
real(r8), allocatable :: obs(:)
logical, allocatable :: missing(:)

! Read the header and verify 
read(file_id, *) header
if(header /= 'obset') then
   write(*, *) 'Error: Expected "obset" in header in read_obs_set_time' 
   stop
end if

! Read the obs_set def index
read(file_id, *) read_obs_set_time%def_index

! Read the number of obs and the number of copies
read(file_id, *) num_obs, num_copies
read_obs_set_time%num_obs = num_obs
! Set the num_copies to 0
read_obs_set_time%num_copies = 0

! Allocate space with 0 for num_copies
allocate(read_obs_set_time%obs(num_obs, 0), &
         read_obs_set_time%missing(num_obs, 0))

! Allocate some storage to read stuff that is discarded
allocate(obs(num_obs), missing(num_obs))

! Read the data for each copy in turn
do i = 1, num_copies
   read(file_id, *) obs
end do

! Read the missing fields for each copy in turn
do i = 1, num_copies
   read(file_id, *) missing
end do

! Read the time 
read_obs_set_time%time = read_time(file_id)

deallocate(obs, missing)

end function read_obs_set_time



subroutine write_obs_set(file_id, set)
!-------------------------------------------------------------------------
!
! Writes an obs_set with all copies of its observations to file, currently
! represented by integer unit number.

implicit none

type(obs_set_type), intent(in) :: set
integer, intent(in) :: file_id

integer :: i, j

! First write ascii header saying set is coming
write(file_id, *) 'obset'

! Write the obs_set_def_index
write(file_id, *) set%def_index

! The number of obs followed by the number of copies
write(file_id, *) set%num_obs, set%num_copies

! Write out the data for each copy in turn
do i = 1, set%num_copies
   write(file_id, *) set%obs(:, i)
end do

! Write out the missing fields for each copy in turn
do i = 1, set%num_copies
   write(file_id, *) set%missing(:, i)
end do

! Write out the time associated with this observation
call write_time(file_id, set%time)

end subroutine write_obs_set


subroutine get_single_obs_value(set, num_obs, obs, copy_in)
!--------------------------------------------------------------------
!
! Sets the obs value number num_obs;  copy_in is optional with default 1.

implicit none

type(obs_set_type), intent(in) :: set
integer, intent(in) :: num_obs
real(r8), intent(out) :: obs
integer, optional, intent(in) :: copy_in

integer :: copy

! Get the appropriate copy
copy = 1
if(present(copy_in)) copy = copy_in 
if(copy < 1 .or. copy > set%num_copies) then
   write(*, *) 'Error: Out of range copy in get_single_obs_value'
   stop
endif

! Make sure the obs index is legal
if(num_obs > set%num_obs .or. num_obs < 1) then
   write(*, *) 'Error: num_obs wrong size in get_single_obs_value'
   stop
endif

! Copy the obs
   obs=set%obs(num_obs, copy)

end subroutine get_single_obs_value


end module obs_set_mod
