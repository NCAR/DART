! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module set_def_list_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! Provides a facility to keep permanent storage for a set of 
! obs_set_defs and allows construction of nested subsets if desired.

use       types_mod, only :  r8
use    location_mod, only : location_type
use   utilities_mod, only : register_module, error_handler, E_ERR
use obs_set_def_mod, only : obs_set_def_type, get_num_obs, &
                            read_obs_set_def, write_obs_set_def, obs_set_def_copy, &
                            os_get_expected_obs => get_expected_obs, &
                            os_get_seq_loc => get_seq_loc, &
                            os_get_obs_location3 => get_obs_location3, &
                            os_get_obs_kind3 => get_obs_kind3, &
                            os_get_diag_obs_err_cov => get_diag_obs_err_cov, &
                            os_get_num_close_states => get_num_close_states, &
                            os_get_close_states => get_close_states

implicit none
private

public set_def_list_type, list_element_type, get_number_obs_subsets, &
   init_set_def_list, add_to_list, get_total_num_obs, &
   write_set_def_list, read_set_def_list, write_list_element, &
   read_list_element, list_element_copy, set_def_list_copy, &
   get_expected_obs, get_diag_obs_err_cov, get_num_close_states, &
   get_seq_loc, get_obs_location2, get_obs_kind2, &
   get_close_states, get_num_sets_in_list

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! For now set up with fixed size array storage declared at allocation
! time. Eventually want a linked list or linked arrays .
type set_def_list_type
    private
   integer :: max_sets
   integer :: num_sets
   type(list_element_type), pointer :: sets(:)
end type set_def_list_type

type list_element_type
    private
   type(obs_set_def_type) :: obs_set
   integer :: index
   integer :: total_num_obs
   integer :: num_subsets, max_subsets
   integer, pointer :: subset(:)
end type list_element_type
! How to really do total_num_obs and keep hierarchical subset
! count. Need back pointers, just ignore???

logical, save :: module_initialized = .false.

contains

!==============================================================


subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



function init_set_def_list(max_sets)
!--------------------------------------------------------------
!
! Initialize a set_def_list with room for a total of max_sets.

implicit none

type(set_def_list_type) :: init_set_def_list
integer, intent(in) :: max_sets

if ( .not. module_initialized ) call initialize_module

init_set_def_list%max_sets = max_sets
init_set_def_list%num_sets = 0
allocate(init_set_def_list%sets(max_sets))

end function init_set_def_list




subroutine set_def_list_copy(list_out, list_in)
!--------------------------------------------------------------
! 
! Full copy of set_def_list_type

implicit none

type(set_def_list_type), intent(in) :: list_in
type(set_def_list_type), intent(out) :: list_out

integer :: i

if ( .not. module_initialized ) call initialize_module

list_out%max_sets = list_in%max_sets
list_out%num_sets = list_in%num_sets
allocate(list_out%sets(list_in%max_sets))

do i = 1, list_in%num_sets
   call list_element_copy(list_out%sets(i), list_in%sets(i))
end do

end subroutine set_def_list_copy




subroutine get_expected_obs(list, list_index, state, obs, num)
!---------------------------------------------------------------
!
! Returns the expected value of the observations in this list_element
! hierarchy of obs_sets given the state_vector. For now, sub_sets are
! not implemented but need to do this.  num allows picking a single
! obs from the set for efficiency.

implicit none

type(set_def_list_type), intent(in) :: list
integer, intent(in) :: list_index
real(r8), intent(in) :: state(:)
real(r8), intent(out) :: obs(:)
integer, intent(in), optional :: num

if ( .not. module_initialized ) call initialize_module

! For now, just call os_get_expected_obs for the single set in the list_element
if(present(num)) then
   call os_get_expected_obs(list%sets(list_index)%obs_set, state, obs, num)
else
   call os_get_expected_obs(list%sets(list_index)%obs_set, state, obs)
endif

! Note that list is not currently used, included for later subset
! applications, make sure this is appropriate

end subroutine get_expected_obs
   



subroutine list_element_copy(list_out, list_in)
!---------------------------------------------------------------
!
! Overloaded to = for list_element

implicit none

type(list_element_type) :: list_out
type(list_element_type), intent(in) :: list_in

integer :: i

if ( .not. module_initialized ) call initialize_module

list_out%obs_set = list_in%obs_set
list_out%index = list_in%index
list_out%total_num_obs = list_in%total_num_obs
list_out%num_subsets = list_in%num_subsets
list_out%max_subsets = list_in%max_subsets
allocate(list_out%subset(list_in%max_subsets))
do i = 1, list_in%max_subsets
   list_out%subset(i) = list_in%subset(i)
end do

end subroutine list_element_copy




function add_to_list(list, set, max_subsets_in)
!--------------------------------------------------------------
!
! Add an obs_set_def to the list. The index of the set in the
! list is returned (this is the handle).
! 

implicit none

integer :: add_to_list
type(set_def_list_type), intent(inout) :: list
type(obs_set_def_type), intent(in) :: set
integer, intent(in), optional :: max_subsets_in

integer :: max_subsets, indx

if ( .not. module_initialized ) call initialize_module

! Default for max_subsets is 0
max_subsets = 0
if(present(max_subsets_in)) max_subsets = max_subsets_in

if(list%num_sets == list%max_sets) call error_handler(E_ERR, 'add_to_list', &
      'Too many sets in list in add_to_list', source, revision, revdate)

! Increment the number of sets in the list
list%num_sets = list%num_sets + 1
indx          = list%num_sets

! Create the list_element
list%sets(indx)%index = indx
list%sets(indx)%total_num_obs = get_num_obs(set)
list%sets(indx)%num_subsets = 0
list%sets(indx)%max_subsets = max_subsets
allocate(list%sets(indx)%subset(max_subsets))

! Copy in the obs_set_def
call obs_set_def_copy(list%sets(indx)%obs_set, set)

! Return the index
add_to_list = indx

end function add_to_list



function get_def_from_list(list, indx)
!--------------------------------------------------------------
!
! Returns an obs_set_def that is the index-th element in the list.
! This does not return info on the subsets which will require 
! extra interfaces. Need error checks on indexing

implicit none

 
type(obs_set_def_type) :: get_def_from_list
type(set_def_list_type), intent(in) :: list
integer, intent(in) :: indx

if ( .not. module_initialized ) call initialize_module

if(indx > list%num_sets .or. indx < 1) call error_handler(E_ERR, 'get_def_from_list',&
    'bad index in get_def_from_list',source, revision, revdate)

call obs_set_def_copy(get_def_from_list, list%sets(indx)%obs_set)

end function get_def_from_list



function get_num_sets_in_list(list)
!---------------------------------------------------------------
!
! Returns the number of total elements in this list.

implicit none

integer :: get_num_sets_in_list
type(set_def_list_type), intent(in) :: list

if ( .not. module_initialized ) call initialize_module

get_num_sets_in_list = list%num_sets

end function get_num_sets_in_list




function get_total_num_obs(list, indx)
!------------------------------------------------------------
!
! Returns the total number of observations in the index set
! from this list.

implicit none

integer :: get_total_num_obs
type(set_def_list_type), intent(in) :: list
integer, intent(in) :: indx

if ( .not. module_initialized ) call initialize_module

if(indx > list%num_sets .or. indx < 1) call error_handler(E_ERR,'get_total_num_obs',&
   'bad index in get_total_num_obs', source, revision, revdate)

get_total_num_obs = list%sets(indx)%total_num_obs

end function get_total_num_obs



subroutine get_diag_obs_err_cov(list, indx, cov)
!----------------------------------------------------
!
! Returns the diagonal part of the observational error
! covariance for this list_element. For now, with no subsets
! implemented, this can just call for this individual set.

implicit none

type(set_def_list_type), intent(in) :: list
integer, intent(in) :: indx
real(r8), intent(out) :: cov(:)

if ( .not. module_initialized ) call initialize_module

call os_get_diag_obs_err_cov(list%sets(indx)%obs_set, cov)

end subroutine get_diag_obs_err_cov



subroutine get_num_close_states(list, indx, radius, num)
!-------------------------------------------------------
!
! Returns the number of close states for the index list_element
! in the list. For now, with no subsets implemented, this can
! just call for this individual set.

implicit none

type(set_def_list_type), intent(in) :: list
integer, intent(in) :: indx
real(r8), intent(in) :: radius
integer, intent(out) :: num(:)

if ( .not. module_initialized ) call initialize_module

call os_get_num_close_states(list%sets(indx)%obs_set, radius, num)

end subroutine get_num_close_states



subroutine get_close_states(list, indx, radius, num, indices, dist, obs_num)
!-------------------------------------------------------
!
! Returns the list of indices of the close states for this
! set_def_list.

implicit none

type(set_def_list_type), intent(in) :: list
integer, intent(in) :: indx
real(r8), intent(in) :: radius
integer, intent(out) :: num(:), indices(:, :)
real(r8), intent(out) :: dist(:, :)
integer, optional, intent(in) :: obs_num

if ( .not. module_initialized ) call initialize_module

if(present(obs_num)) then
   call os_get_close_states(list%sets(indx)%obs_set, radius, num, &
   indices, dist, obs_num)
else
   call os_get_close_states(list%sets(indx)%obs_set, radius, num, &
   indices, dist)
endif

end subroutine get_close_states




function get_number_obs_subsets(list, indx)
!-----------------------------------------------------------
!
! Returns the number of obs subsets included in this set.

implicit none

integer :: get_number_obs_subsets
type(set_def_list_type), intent(in) :: list
integer, intent(in) :: indx

if ( .not. module_initialized ) call initialize_module

! Return 0 for limited implementation
! Eventually need to do recursive pass (or forbid adding
! subsets below existing subsets?)
get_number_obs_subsets = 0

end function get_number_obs_subsets





subroutine write_set_def_list(file_id, list)
!----------------------------------------------------------
!

implicit none

integer, intent(in) :: file_id
type(set_def_list_type), intent(in) :: list

integer :: i

if ( .not. module_initialized ) call initialize_module

! write(*,*)'DEBUG(set_def_list_mod:write_set_def_list) file id is ',file_id

! Write a header
write(file_id, *) 'defls'

! Write the number of sets
write(file_id, *) list%num_sets

! Loop through each set and output its info
do i = 1, list%num_sets
   call write_list_element(file_id, list%sets(i))
end do

! write(*,*)'DEBUG(set_def_list_mod:write_set_def_list) completed.',file_id

end subroutine write_set_def_list




function read_set_def_list(file_id)
!---------------------------------------------------------
!

implicit none

type(set_def_list_type):: read_set_def_list
integer, intent(in) :: file_id

character(len=5) :: header
integer :: i, num_sets
character(len=129) :: stringerror

if ( .not. module_initialized ) call initialize_module

! Read header
read(file_id, *) header
write(stringerror,*)'expected "defls" in read_set_def_list, got ',header
if(header /= 'defls') call error_handler(E_ERR,'read_set_def_list',&
     stringerror, source, revision, revdate)

! Read the number of sets
read(file_id, *) num_sets

! Initialize storage
read_set_def_list = init_set_def_list(num_sets)
read_set_def_list%num_sets = num_sets


! Need to verify that storage is okay through here, fear null pointer
! drops??

! Loop through each set and read its info (eventually will have subsets here)
do i = 1, num_sets
! May want to re-examine function syntax for reads and perhaps make then copy
! subroutines ???
   call list_element_copy(read_set_def_list%sets(i), read_list_element(file_id))
end do

end function read_set_def_list



subroutine write_list_element(file_id, element)
!----------------------------------------------------------
!

implicit none

integer, intent(in) :: file_id
type(list_element_type), intent(in) :: element

integer :: i

if ( .not. module_initialized ) call initialize_module

! write(*,*)'DEBUG(set_def_list_mod:write_list_element): file_id is ',file_id

! Write out a sequence of integer arguments
write(file_id, *) element%index
write(file_id, *) element%total_num_obs
write(file_id, *) element%num_subsets
do i = 1, element%num_subsets
   write(file_id, *) element%subset(i)
end do

! write(*,*)'DEBUG(set_def_list_mod:write_list_element): before write_obs_set_def ',file_id
! Write out the obs_set_def for this set
call write_obs_set_def(file_id, element%obs_set)

end subroutine write_list_element



function read_list_element(file_id)
!----------------------------------------------------------
!

implicit none

type(list_element_type) :: read_list_element
integer, intent(in) :: file_id

integer :: i

if ( .not. module_initialized ) call initialize_module

! Read in a sequence of integer arguments
read(file_id, *) read_list_element%index
read(file_id, *) read_list_element%total_num_obs  
read(file_id, *) read_list_element%num_subsets
read_list_element%max_subsets = read_list_element%num_subsets

! Allocate the storage for the subsets
allocate(read_list_element%subset(read_list_element%num_subsets))
do i = 1, read_list_element%num_subsets
   read(file_id, *) read_list_element%subset(i)
end do

! Read the obs_set_def for this set
call obs_set_def_copy(read_list_element%obs_set, read_obs_set_def(file_id))
!!!read_list_element%obs_set = read_obs_set_def(file_id)


end function read_list_element


subroutine get_seq_loc(list, list_index, obsloc0, num)
!---------------------------------------------------------------
!
implicit none

type(set_def_list_type), intent(in) :: list
integer, intent(in) :: list_index
real(r8), intent(out) :: obsloc0(:)
integer, intent(in), optional :: num

if ( .not. module_initialized ) call initialize_module

! For now, just call os_get_seq_loc for the single set in the list_element
if(present(num)) then
   call os_get_seq_loc(list%sets(list_index)%obs_set, obsloc0, num)
endif

end subroutine get_seq_loc


   
subroutine get_obs_location2(list, list_index, obsloc)
!---------------------------------------------------------------
!
implicit none

type(set_def_list_type), intent(in) :: list
integer, intent(in) :: list_index
real(r8), intent(out) :: obsloc(:, :)

if ( .not. module_initialized ) call initialize_module

! For now, just call os_get_seq_loc for the single set in the list_element
   call os_get_obs_location3(list%sets(list_index)%obs_set, obsloc)

end subroutine get_obs_location2
   
   
subroutine get_obs_kind2(list, list_index, obskind)
!---------------------------------------------------------------
!
implicit none

type(set_def_list_type), intent(in) :: list
integer, intent(in) :: list_index
real(r8), intent(out) :: obskind(:)

if ( .not. module_initialized ) call initialize_module

! For now, just call os_get_seq_loc for the single set in the list_element
   call os_get_obs_kind3(list%sets(list_index)%obs_set, obskind)

end subroutine get_obs_kind2



end module set_def_list_mod
