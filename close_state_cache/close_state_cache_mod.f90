! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module close_state_cache_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use obs_sequence_mod, only : obs_sequence_type, get_obs_def_index, &
                             get_num_close_states, get_close_states

implicit none
private

public :: close_state_cache_type, cache_init, get_close_cache

! CVS Generated file description for error handling, do not edit 
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type close_state_cache_type
   private
   integer :: size, num
   integer, pointer :: def_index(:), num_obs(:), order(:)
   real(r8), pointer :: radius(:)
   type(cache_element_type), pointer :: elt(:)
end type close_state_cache_type

type cache_element_type
   private
   integer, pointer :: num_close(:), close(:, :)
   real(r8), pointer :: dist(:, :)
end type cache_element_type

logical, save :: module_initialized = .false.

contains

!===============================================

subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module




subroutine cache_init(cache, size)
!-----------------------------------------------------
!
! Initializes storage for a close state cache with size elements

implicit none

type(close_state_cache_type), intent(out) :: cache
integer, intent(in) :: size

if ( .not. module_initialized ) call initialize_module

cache%size = size
cache%num = 0
allocate(cache%order(size), cache%def_index(size), &
   cache%num_obs(size), cache%radius(size), cache%elt(size))

cache%order = 0
cache%def_index = 0
cache%num_obs = 0
cache%radius = 0.0

end  subroutine cache_init



subroutine get_close_cache(cache, seq, index, radius, &
   num_obs_in_set, num_close, close, dist)
!----------------------------------------------------
!
! Gets the close state info for index-th observation in 
! sequence seq with radius for cutoff distance. The cache
! is used and / or update appropriately.

implicit none

type(close_state_cache_type), intent(inout) :: cache
type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: index, num_obs_in_set
real(r8), intent(in) :: radius
integer, pointer :: num_close(:), close(:, :)
real(r8), pointer :: dist(:, :)

integer :: def_index, i, ind

if ( .not. module_initialized ) call initialize_module

! Get the index for the associated obs_set_def
def_index = get_obs_def_index(seq, index)

! See if this index is in the cache
do i = 1, cache%num
   if(cache%def_index(i) == def_index .and. cache%radius(i) == radius) then
! Already in cache at index i, just return the pointers
      num_close => cache%elt(i)%num_close
      close => cache%elt(i)%close
      dist => cache%elt(i)%dist
      return
   endif
end do

! Need to put this in the cache; need to work on update
! For more general use want a last touched update order
ind = cache%num + 1
cache%num = cache%num + 1
cache%def_index(i) = def_index
cache%radius(i) = radius
cache%num_obs(i) = num_obs_in_set

if(ind > cache%size) then
      call error_handler(E_ERR,'get_close_cache', 'cache is full', source, revision, revdate)
endif  

! Can I combine the allocations to avoid multiple memory stabs
allocate(cache%elt(ind)%num_close(num_obs_in_set))
call get_num_close_states(seq, index, radius, cache%elt(ind)%num_close)
num_close => cache%elt(ind)%num_close

! Loop to allocate and get the close states for each observation
allocate(cache%elt(ind)%close(num_obs_in_set, maxval(num_close)), &
   cache%elt(ind)%dist(num_obs_in_set, maxval(num_close)))

! It would be a savings of a factor of 2 more or less if the 
! get_num_close_states could be eliminated.
call get_close_states(seq, index, radius, cache%elt(ind)%num_close, &
   cache%elt(ind)%close, cache%elt(ind)%dist)


! CAN COMBINE THESEE WITH FIRST CASE
close => cache%elt(ind)%close
dist => cache%elt(ind)%dist

end subroutine get_close_cache




end module close_state_cache_mod
