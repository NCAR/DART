! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_set_def_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR, register_module
use   obs_def_mod, only : obs_def_type, get_obs_location, read_obs_def, &
                        write_obs_def, get_error_variance, &
                        od_get_expected_obs => get_expected_obs, &
                        od_get_close_states=>get_close_states, &
                        od_get_num_close_states=>get_num_close_states, &
                        od_get_seq_loc => get_seq_loc, &
                        od_get_obs_location4 => get_obs_location4, &
                        od_get_obs_kind4 => get_obs_kind4
use   location_mod, only: location_type, vert_is_level

implicit none
private

public obs_set_def_type, init_obs_set_def, get_diag_obs_err_cov, &
       get_expected_obs, get_obs_def, get_num_obs, get_obs_locations, &
       get_close_states, get_num_close_states, add_obs, &
       diag_obs_err_cov, read_obs_set_def, write_obs_set_def, obs_set_def_copy, &
       get_seq_loc, get_obs_location3, get_obs_kind3

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type obs_set_def_type
   private
   real(r8), pointer  :: diag_cov(:)
   type(obs_def_type), pointer :: obs_defs(:)
   integer :: max_num_obs_defs
   integer :: num_obs_defs
end type obs_set_def_type


logical, save :: module_initialized = .false.

contains

!===============================================================

subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



function init_obs_set_def(max_num_obs)
!--------------------------------------------------------------
!
! Returns a handle to a structure for an obs_set_def. Specifying
! the max_num_obs up front is annoying but 
! precludes need for dynamic data structures at this point. Might
! want to make them optional at a much later time

implicit none

type(obs_set_def_type) :: init_obs_set_def
integer, intent(in) :: max_num_obs

if ( .not. module_initialized ) call initialize_module

! Allocate storage for the inividual observations
allocate(init_obs_set_def%diag_cov(max_num_obs), &
   init_obs_set_def%obs_defs(max_num_obs))
init_obs_set_def%max_num_obs_defs = max_num_obs
init_obs_set_def%num_obs_defs = 0

end function init_obs_set_def




subroutine obs_set_def_copy(set_out, set_in)
!--------------------------------------------------------------
!
! Needs to be used in place of assignement for obs_set_defs to
! get full copy functionality.

implicit none

type(obs_set_def_type), intent(out) :: set_out
type(obs_set_def_type), intent(in) :: set_in

integer :: i

if ( .not. module_initialized ) call initialize_module

set_out = init_obs_set_def(set_in%max_num_obs_defs)

! Loop to add in the obs_defs
do i = 1, set_in%num_obs_defs
   call add_obs(set_out, set_in%obs_defs(i))
end do

end subroutine obs_set_def_copy




subroutine get_diag_obs_err_cov(set, cov)
!-------------------------------------------------------------
!
! Returns the diagonal of the observation error covariance for the
! observation set. Dies if cov isn't big enough?

implicit none

type(obs_set_def_type), intent(in) :: set
real(r8), intent(out) :: cov(:)

if ( .not. module_initialized ) call initialize_module

if(size(cov) < set%num_obs_defs) then
      call error_handler(E_ERR,'get_diag_obs_err_cov', 'cov array too small', source, revision, revdate)
endif

cov(1:set%num_obs_defs) = set%diag_cov

end subroutine get_diag_obs_err_cov



subroutine get_expected_obs(set, state, obs, num)
!--------------------------------------------------------------
!
! Returns the expected value of the observations in this set given
! the state vector

implicit none

type(obs_set_def_type), intent(in) :: set
real(r8), intent(in) :: state(:)
real(r8), intent(out) :: obs(:)
integer, intent(in), optional :: num

integer :: i

if ( .not. module_initialized ) call initialize_module

if(present(num)) then
   obs(1) = od_get_expected_obs(set%obs_defs(num), state)
else

! Check for sufficient storage
   if(size(obs) < set%num_obs_defs) then
      call error_handler(E_ERR,'get_expected_obs', 'obs array too small', source, revision, revdate)
   endif


   do i = 1, set%num_obs_defs
      obs(i) = od_get_expected_obs(set%obs_defs(i), state)
   end do
endif

end subroutine get_expected_obs





function get_obs_def(set, index)
!------------------------------------------------------------
!
! Returns the index-th obs_def in the set.

implicit none

type(obs_def_type) :: get_obs_def
type(obs_set_def_type), intent(in) :: set
integer, intent(in) :: index

if ( .not. module_initialized ) call initialize_module

! Check for index exceeding total number of obs available
if(index > get_num_obs(set)) then
      call error_handler(E_ERR,'get_obs_def', 'index too large', source, revision, revdate)
endif

get_obs_def = set%obs_defs(index)

end function get_obs_def



function get_num_obs(set)
!-----------------------------------------------------------
!
! Returns the number of obs included in this set.

implicit none

integer :: get_num_obs
type(obs_set_def_type), intent(in) :: set

if ( .not. module_initialized ) call initialize_module

get_num_obs = set%num_obs_defs

end function get_num_obs





subroutine get_obs_locations(set, locations)
!--------------------------------------------------------
!
! Gets locations for all observations in the set. 

implicit none

type(obs_set_def_type), intent(in) :: set
type(location_type), intent(out) :: locations(:)

integer i

if ( .not. module_initialized ) call initialize_module

! Make sure locations array is big enough
if(size(locations) < set%num_obs_defs) then
   call error_handler(E_ERR,'get_obs_locations','locations array too small', source, revision, revdate)
endif

do i = 1, set%num_obs_defs
   locations(i) = get_obs_location(set%obs_defs(i))
end do

end subroutine get_obs_locations



subroutine get_close_states(set, radius, number, indices, dist, obs_num)
!--------------------------------------------------------
!
! Returns the number of close states for each observation in
! the set followed by the indices for those. Storage fixed at
! call for now.

type(obs_set_def_type), intent(in) :: set
real(r8), intent(in) :: radius
integer, intent(out) :: number(:)
integer, intent(out) :: indices(:, :)
real(r8), intent(out) :: dist(:, :)
integer, optional, intent(in) :: obs_num

integer :: i

if ( .not. module_initialized ) call initialize_module

! Make sure that number array is big enough to hold all obs in set
if(present(obs_num)) then
   if(size(number) < 1) then
      call error_handler(E_ERR,'get_close_states', 'number array too small', source, revision, revdate)
   endif
else
   if(size(number) < set%num_obs_defs) then
      call error_handler(E_ERR,'get_close_states', 'number array too small', source, revision, revdate)
   endif
endif

! Loop through the observations and get their close states
if(present(obs_num)) then
   call  od_get_close_states(set%obs_defs(obs_num), radius, number(1), &
         indices(1, :), dist(1, :))
else
   do i = 1, set%num_obs_defs
      call od_get_close_states(set%obs_defs(i), radius, number(i), &
         indices(i, :), dist(i, :))
   end do
endif

end subroutine get_close_states




subroutine get_num_close_states(set, radius, number)
!--------------------------------------------------------
!
! Returns the number of close states for each observation in
! the set.

implicit none

type(obs_set_def_type), intent(in) :: set
real(r8), intent(in) :: radius
integer, intent(out) :: number(:)

integer :: i

if ( .not. module_initialized ) call initialize_module

! Make sure that number array is big enough to hold all obs in set
if(size(number) < set%num_obs_defs) then
      call error_handler(E_ERR,'get_num_close_states', 'number array too small', source, revision, revdate)
endif

! Loop through the observations and get their close states
do i = 1, set%num_obs_defs
   number(i) = od_get_num_close_states(set%obs_defs(i), radius)
end do

end subroutine get_num_close_states



subroutine add_obs(set, obs_def)
!-------------------------------------------------------
!
! Adds this obs_def to the set; if the max number of 
! obs would be exceeded it's an error.

implicit none

type(obs_set_def_type), intent(inout) :: set
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

if(set%num_obs_defs == set%max_num_obs_defs) then
      call error_handler(E_ERR,'add_obs', 'no room for another obs_def', source, revision, revdate)
endif

set%num_obs_defs = set%num_obs_defs + 1
set%obs_defs(set%num_obs_defs) = obs_def

! Update the covariance array (which is redundant storage for diag
set%diag_cov(set%num_obs_defs) = get_error_variance(obs_def)

end subroutine add_obs




function diag_obs_err_cov(set)
!-----------------------------------------------------
!
! Returns true if covariance for this set if diagonal
! For now, all are diagonal so always returns true.

implicit none

logical :: diag_obs_err_cov
type(obs_set_def_type), intent(in) :: set

if ( .not. module_initialized ) call initialize_module

diag_obs_err_cov = .true.

end function diag_obs_err_cov



function read_obs_set_def(file_id)
!------------------------------------------------------
! 
! Reads an obs_set_def from the file and sets up storage, etc.

implicit none

type(obs_set_def_type) :: read_obs_set_def
integer, intent(in) :: file_id

character(len=5) :: header
integer :: num, i

if ( .not. module_initialized ) call initialize_module

! Read the header for a set_def
read(file_id, *) header
if(header /= 'obset') then
    call error_handler(E_ERR,'read_obs_set_deF', 'expected string "obset" in input file', &
                       source, revision, revdate)
endif

! Read the number of observations in the set
read(file_id, *) num

! Initialize the obs_set for this many obs
read_obs_set_def = init_obs_set_def(num)
read_obs_set_def%num_obs_defs = num

! Loop to read each obs in the set and insert it
do i = 1, num
   read_obs_set_def%obs_defs(i) = read_obs_def(file_id)
end do

end function read_obs_set_deF





subroutine write_obs_set_def(file_id, set)
!------------------------------------------------------
!
! Outputs the obs_set_def to the file, currently represented
! by an integer unit number.

implicit none

integer, intent(in) :: file_id
type(obs_set_def_type), intent(in) :: set

integer :: i

if ( .not. module_initialized ) call initialize_module

!write(*,*)'DEBUG(obs_set_def_mod:write_obs_set_def), fid/N is ', &
!          file_id, set%num_obs_defs

! Write a header starting the obs_set_def
write(file_id, *) 'obset'

! Write the number of obs to follow
write(file_id, *) set%num_obs_defs

! For now no need to write out covariance

! Loop to write out each observation
do i = 1, set%num_obs_defs
!   write(*,*)'DEBUG(obs_set_def_mod:write_obs_set_def), i is ',i
   call write_obs_def(file_id, set%obs_defs(i))
end do

! write(*,*)'DEBUG(obs_set_def_mod:write_obs_set_def), leaving ...'

end subroutine write_obs_set_def


subroutine get_seq_loc(set, obsloc0, num)
!--------------------------------------------------------------
!
implicit none

type(obs_set_def_type), intent(in) :: set
real(r8), intent(out) :: obsloc0(:)
integer, intent(in), optional :: num

integer :: i

if ( .not. module_initialized ) call initialize_module

if(present(num)) then
   call od_get_seq_loc(set%obs_defs(num), obsloc0)
endif

end subroutine get_seq_loc


subroutine get_obs_location3(set, obsloc)
!--------------------------------------------------------------
!
implicit none

type(obs_set_def_type), intent(in) :: set
real(r8), intent(out) :: obsloc(:, :)
real(r8) :: obsloc0(3)

integer :: i

   if ( .not. module_initialized ) call initialize_module

   do i = 1, set%num_obs_defs
      call od_get_obs_location4(set%obs_defs(i), obsloc0)
      obsloc(i,1) = obsloc0(1)
      obsloc(i,2) = obsloc0(2)
      obsloc(i,3) = obsloc0(3)
   end do

end subroutine get_obs_location3


subroutine get_obs_kind3(set, obskind)
!--------------------------------------------------------------
!
implicit none

type(obs_set_def_type), intent(in) :: set
real(r8), intent(out) :: obskind(:)
real(r8) :: obskind0

integer :: i

   if ( .not. module_initialized ) call initialize_module

   do i = 1, set%num_obs_defs
      call od_get_obs_kind4(set%obs_defs(i), obskind0)
      obskind(i) = obskind0
   end do

end subroutine get_obs_kind3


end module obs_set_def_mod
