! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> \file obs_sequence_mod.f90  modifing this to have distributed identity obs
!> \dir obs_sequence modifing obs_sequence to have distributed identity obs

!> @brief For observations sequence stuff
!> get expected obs is in here
module obs_sequence_mod

! WARNING OPERATOR OVERLOAD FOR EQUIVALENCE???
! FURTHER WARNING: Compiler problems exist with the use of assignment(=) in
! use only statement. First, can only use it at the level above if the internals
! of the type are not private. Second, if I inherit assignment(=) from obs_def
! and also define one in obs_sequence, I get an error if I try to make it public
! to a module that uses obs_sequence but not obs_def with the intel compiler. No
! obvious workaround exists. For now, make modules at higher levels use explicit
! copy subroutines. USERS MUST BE VERY CAREFUL TO NOT DO DEFAULT ASSIGNMENT
! FOR THESE TYPES THAT HAVE COPY SUBROUTINES.

use        types_mod, only : r8, DEG2RAD, MISSING_R8, metadatalength
use     location_mod, only : location_type, interactive_location, &
                             is_location_in_region
use      obs_def_mod, only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, copy_obs_def, &
                             interactive_obs_def, get_obs_def_location, &
                             get_expected_obs_from_def, get_obs_kind, &
                             get_obs_def_key, &
                             get_expected_obs_from_def_distrib_state !HK
use     obs_kind_mod, only : write_obs_kind, read_obs_kind, max_obs_kinds, &
                             get_obs_kind_index
use time_manager_mod, only : time_type, operator(>), operator(<), &
                             operator(>=), operator(/=), set_time, &
                             operator(-), operator(+), operator(==)
use    utilities_mod, only : get_unit, close_file,                       &
                             register_module, error_handler,             &
                             find_namelist_in_file, check_namelist_read,   &
                             E_ERR, E_WARN, E_MSG, nmlfileunit, do_output, &
                             do_nml_file, do_nml_term
!HK
use mpi_utilities_mod, only : task_count, my_task_id
use ensemble_manager_mod, only: get_var_owner_index, map_pe_to_task, &
                                ensemble_type

use mpi

implicit none
private

interface assignment(=)
   module procedure copy_obs
end interface

! Public interfaces for obs sequences
public :: obs_sequence_type, init_obs_sequence, interactive_obs_sequence, &
   get_num_copies, get_num_qc, get_num_obs, get_max_num_obs, &
   get_copy_meta_data, get_qc_meta_data, get_next_obs, get_prev_obs, &
   insert_obs_in_seq, delete_obs_from_seq, set_copy_meta_data, &
   set_qc_meta_data, get_first_obs, get_last_obs, add_copies, add_qc, &
   write_obs_seq, read_obs_seq, set_obs, append_obs_to_seq, &
   get_obs_from_key, get_obs_time_range, get_time_range_keys, &
   get_num_times, get_num_key_range, &
   static_init_obs_sequence, destroy_obs_sequence, read_obs_seq_header, &
   get_expected_obs, delete_seq_head, delete_seq_tail, &
   get_next_obs_from_key, get_prev_obs_from_key, delete_obs_by_typelist, &
   select_obs_by_location, delete_obs_by_qc, delete_obs_by_copy,         &
   get_expected_obs_distrib_state !HK

! Public interfaces for obs
public :: obs_type, init_obs, destroy_obs, get_obs_def, set_obs_def, &
   get_obs_values, set_obs_values, replace_obs_values, get_qc, set_qc, &  
   read_obs, write_obs, replace_qc, interactive_obs, copy_obs, assignment(=), &
   get_obs_key, copy_partial_obs

! Public interfaces for obs covariance modeling
public :: obs_cov_type

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type obs_sequence_type
   private
   integer :: num_copies
   integer :: num_qc
   integer :: num_obs
   integer :: max_num_obs
   ! F95 allows pointers to be initialized to a known value.
   ! However, if you get an error on the following lines from your
   ! compiler, remove the => NULL() from the end of the 5 lines below.
   character(len = metadatalength), pointer :: copy_meta_data(:)  => NULL()
   character(len = metadatalength), pointer :: qc_meta_data(:)    => NULL()
   integer :: first_time
   integer :: last_time
!   integer :: first_avail_time, last_avail_time
   type(obs_type), pointer :: obs(:)   => NULL()
! What to do about groups
end type obs_sequence_type

type obs_type
   private
! The key is needed to indicate the element number in the storage for the obs_sequence
! Do I want to enforce the identity of the particular obs_sequence?
   integer :: key
   type(obs_def_type) :: def
   real(r8), pointer :: values(:)  => NULL()
   real(r8), pointer :: qc(:)      => NULL()
   ! Put sort indices directly into the data structure
   integer :: prev_time, next_time
   integer :: cov_group
end type obs_type

type obs_cov_type
   private
   integer :: num_cov_groups
! ??????
end type obs_cov_type

! for errors
character(len=129) :: msgstring, string1

!-------------------------------------------------------------
! Namelist with default values
! write_binary_restart_files  == .true.  -> use unformatted file format.
!                                     Full precision, faster, smaller,
!                                     but not as portable.

logical :: write_binary_obs_sequence = .false.

namelist /obs_sequence_nml/ write_binary_obs_sequence

!--------------------------------------------------------------


contains

!--------------------------------------------------------------

subroutine static_init_obs_sequence

! reads namelist and registers module
! Read the namelist input

integer :: iunit, io

call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_sequence_nml", iunit)
read(iunit, nml = obs_sequence_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_sequence_nml")

if (do_nml_file()) write(nmlfileunit,nml=obs_sequence_nml)
if (do_nml_term()) write(     *     ,nml=obs_sequence_nml)

end subroutine static_init_obs_sequence

!--------------------------------------------------------------

!WHAT ABOUT PASS THROUGHS TO THE OBS_DEF???
! WHAT ABOUT copy_obs_sequence similar to read.
!-------------------------------------------------
subroutine init_obs_sequence(seq, num_copies, num_qc, &
   expected_max_num_obs)

! Constructor for an obs_sequence

type(obs_sequence_type), intent(out) :: seq
integer, intent(in) :: num_copies, num_qc, expected_max_num_obs

integer :: i

seq%num_copies  = num_copies
seq%num_qc      = num_qc
seq%num_obs     = 0
seq%max_num_obs = expected_max_num_obs

allocate(seq%copy_meta_data(seq%num_copies), &
         seq%qc_meta_data(seq%num_qc), &
         seq%obs(seq%max_num_obs) )

do i = 1, seq%num_copies
   seq%copy_meta_data(i) = 'Copy metadata not initialized'
end do

do i = 1, seq%num_qc
   seq%qc_meta_data(i) = 'QC metadata not initialized'
end do

! Initialize the pointers to allocated and initialize to something benign
! (Go ahead and allocated even in the case the counts are 0.)
do i = 1, seq%max_num_obs
   allocate(seq%obs(i)%values(num_copies))
   if (num_copies > 0) seq%obs(i)%values = MISSING_R8
   allocate(seq%obs(i)%qc(num_qc))
   if (num_qc > 0) seq%obs(i)%qc = 0.0_r8
end do
seq%first_time = -1
seq%last_time  = -1
!seq%first_avail_time = -1
!seq%last_avail_time = -1

end subroutine init_obs_sequence


!--------------------------------------------------------------


subroutine destroy_obs_sequence(seq)
! Destructor for an obs_sequence

type(obs_sequence_type), intent(inout) :: seq

integer :: i

if ( seq%max_num_obs > 0 ) then

   if (associated(seq%copy_meta_data)) then
      deallocate(seq%copy_meta_data)
      nullify(seq%copy_meta_data)
   endif
   if (associated(seq%qc_meta_data)) then
      deallocate(seq%qc_meta_data)
      nullify(seq%qc_meta_data)
   endif
          
   do i = 1, seq%max_num_obs
   ! seq%obs is a derived type, not a pointer.
   !    if (associated(seq%obs(i))) call destroy_obs( seq%obs(i) )
      call destroy_obs( seq%obs(i) )
   end do

   ! Also free up the obs storage in the sequence
   if(associated(seq%obs)) then 
      deallocate(seq%obs)
      nullify(seq%obs)
   else
      print *, 'destroy_obs_sequence called but seq%obs not associated'
   endif

   seq%first_time  = -1
   seq%last_time   = -1
   seq%num_copies  = -1                                                       
   seq%num_qc      = -1
   seq%num_obs     = -1
   seq%max_num_obs = -1                                                       

endif


end subroutine destroy_obs_sequence


!--------------------------------------------------------------

function interactive_obs_sequence()

! Interactive creation of an observation sequence
type(obs_sequence_type) :: interactive_obs_sequence

type(obs_type)     :: obs, prev_obs
type(obs_def_type) :: obs_def
type(time_type)    :: obs_time, prev_time
integer            :: max_num_obs, num_copies, num_qc, i, end_it_all

write(*, *) 'Input upper bound on number of observations in sequence'
read(*, *) max_num_obs

write(*, *) 'Input number of copies of data (0 for just a definition)'
read(*, *) num_copies

write(*, *) 'Input number of quality control values per field (0 or greater)'
read(*, *) num_qc

! Initialize an obs_sequence structure
call init_obs_sequence(interactive_obs_sequence, num_copies, num_qc, max_num_obs)

do i = 1, num_copies
   write(*, *) 'input meta data for data copy ', i
   read(*, *) interactive_obs_sequence%copy_meta_data(i)
end do

do i = 1, num_qc
   write(*, *) 'input meta data for qc field ', i
   read(*, *) interactive_obs_sequence%qc_meta_data(i)
end do

! Initialize the obs variable
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! Loop to initialize each observation in turn; terminate by -1
do i = 1, max_num_obs
   write(*, *) 'input a -1 if there are no more obs'
   read(*, *) end_it_all
   if(end_it_all == -1) exit
   ! Need to have key available for specialized observation modules
   call interactive_obs(num_copies, num_qc, obs, i)
   if(i == 1) then
      call insert_obs_in_seq(interactive_obs_sequence, obs)
   else
      ! if this is not the first obs, make sure the time is larger
      ! than the previous observation.  if so, we can start the
      ! linked list search at the location of the previous obs.
      ! otherwise, we have to start at the beginning of the entire
      ! sequence to be sure the obs are ordered correctly in
      ! monotonically increasing times. 
      call get_obs_def(obs, obs_def)
      obs_time = get_obs_def_time(obs_def)
      call get_obs_def(prev_obs, obs_def)
      prev_time = get_obs_def_time(obs_def)
      if(prev_time > obs_time) then
         call insert_obs_in_seq(interactive_obs_sequence, obs)
      else
         call insert_obs_in_seq(interactive_obs_sequence, obs, prev_obs)
      endif
   endif
   prev_obs = obs
end do

call destroy_obs(obs)
call destroy_obs(prev_obs)

end function interactive_obs_sequence


!---------------------------------------------------------

!> @brief Compute forward operator for set of obs in sequence for distributed state vector. 
!> @todo does this need to be for a set of obs?
subroutine get_expected_obs_distrib_state(seq, keys, ens_index, state_time, isprior, &
   istatus, assimilate_this_ob, evaluate_this_ob, state_ens_handle, win, expected_obs)

use mpi_utilities_mod, only : datasize ! This is here rather than at the top because the mpi_get calls will most likely get wraped up inside mpi_utilities

type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: keys(:)
integer,                 intent(in)    :: ens_index
type(time_type),         intent(in)    :: state_time
logical,                 intent(in)    :: isprior
integer,                 intent(out)   :: istatus(:)
logical,                 intent(out)   :: assimilate_this_ob, evaluate_this_ob
!HK
type(ensemble_type),     intent(in)    :: state_ens_handle
real(r8), dimension(:),  intent(inout) :: expected_obs !> @todo needs to be 2d for a set of obs
integer, intent(in)                    :: win !< window for mpi remote memory access

integer              :: num_obs, i
!type(location_type) :: location
type(obs_type)       :: obs
type(obs_def_type)   :: obs_def
integer              :: obs_kind_ind

! HK
integer                        :: ierr
integer(KIND=MPI_ADDRESS_KIND) :: target_disp ! must be mpi_address_kind to avoid seg faults on some systems
integer owner_of_state
integer element_index

num_obs = size(keys)

! NEED to initialize istatus to okay value
istatus = 0

! Initialize the observation type
!!! Can actually init with the correct size here if wanted
call init_obs(obs, 0, 0)

do i = 1, num_obs !> @todo do you ever use this with more than one obs?
   call get_obs_from_key(seq, keys(i), obs)
   call get_obs_def(obs, obs_def)
   !location = get_obs_def_location(obs_def)
   obs_kind_ind = get_obs_kind(obs_def)

   ! Check in kind for negative for identity obs
   if(obs_kind_ind < 0) then
      if ( -obs_kind_ind > state_ens_handle%num_vars ) call error_handler(E_ERR, &
         'get_expected_obs', &
         'identity obs is outside of state vector ', &
         source, revision, revdate)

      ! Find which task has the element of state vector
      call get_var_owner_index(-1*obs_kind_ind, owner_of_state, element_index) ! pe
      owner_of_state = map_pe_to_task(state_ens_handle, owner_of_state)        ! task

      if (my_task_id() == owner_of_state) then

         expected_obs = state_ens_handle%copies(:, element_index)

      else

         target_disp = ( element_index - 1) * state_ens_handle%num_copies

         call mpi_win_lock(MPI_LOCK_SHARED, owner_of_state, 0 , win, ierr)
         call mpi_get(expected_obs, state_ens_handle%num_copies, datasize, owner_of_state, target_disp, state_ens_handle%num_copies, datasize, win, ierr)
         call mpi_win_unlock(owner_of_state, win, ierr)

      endif

      assimilate_this_ob = .true.; evaluate_this_ob = .false.
   
   else ! do forward operator for this kind
      !> Q. Do we loop around copies here? This would mean ens_size*times the comumication
      !> The alternative is to alter the code in model_mod.f90 to work on arrays of ensemble size.
      !> Currently looping in model_mod.f90 for lorenz_96 

      call get_expected_obs_from_def_distrib_state(keys(i), obs_def, obs_kind_ind, &
         ens_index, state_time, isprior, istatus, &
         assimilate_this_ob, evaluate_this_ob, expected_obs, state_ens_handle, win)


   endif
end do

! need to free any observation specific storage that
! might have been allocated.
call destroy_obs(obs)

end subroutine get_expected_obs_distrib_state

!---------------------------------------------------------

subroutine get_expected_obs(seq, keys, ens_index, state, state_time, isprior, &
   obs_vals, istatus, assimilate_this_ob, evaluate_this_ob)

! Compute forward operator for set of obs in sequence

type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: keys(:)
integer,                 intent(in)  :: ens_index
real(r8),                intent(in)  :: state(:)
type(time_type),         intent(in)  :: state_time
logical,                 intent(in)  :: isprior
real(r8),                intent(out) :: obs_vals(:)
integer,                 intent(out) :: istatus
logical,                 intent(out) :: assimilate_this_ob, evaluate_this_ob

integer              :: num_obs, i
!type(location_type) :: location
type(obs_type)       :: obs
type(obs_def_type)   :: obs_def
integer              :: obs_kind_ind

num_obs = size(keys)

! NEED to initialize istatus to okay value
istatus = 0

! Initialize the observation type
!!! Can actually init with the correct size here if wanted
call init_obs(obs, 0, 0)

do i = 1, num_obs
   call get_obs_from_key(seq, keys(i), obs)
   call get_obs_def(obs, obs_def)
   !location = get_obs_def_location(obs_def)
   obs_kind_ind = get_obs_kind(obs_def)
! Check in kind for negative for identity obs
   if(obs_kind_ind < 0) then
      if ( -obs_kind_ind > size(state) ) call error_handler(E_ERR, &
         'get_expected_obs', &
         'identity obs is outside of state vector ', &
         source, revision, revdate)
      obs_vals(i) = state(-1 * obs_kind_ind)
      assimilate_this_ob = .true.; evaluate_this_ob = .false.
! Otherwise do forward operator for this kind
   else
      call get_expected_obs_from_def(keys(i), obs_def, obs_kind_ind, &
         ens_index, state, state_time, isprior, obs_vals(i), istatus, &
         assimilate_this_ob, evaluate_this_ob)
   endif
end do

! need to free any observation specific storage that
! might have been allocated.
call destroy_obs(obs)

end subroutine get_expected_obs

!---------------------------------------------------------

function get_num_copies(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_num_copies

get_num_copies = seq%num_copies

end function get_num_copies

!-------------------------------------------------

function get_num_qc(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_num_qc

get_num_qc= seq%num_qc

end function get_num_qc

!-------------------------------------------------

function get_num_obs(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_num_obs

get_num_obs = seq%num_obs

end function get_num_obs

!-------------------------------------------------

function get_max_num_obs(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_max_num_obs

get_max_num_obs = seq%max_num_obs

end function get_max_num_obs
!-------------------------------------------------

function get_copy_meta_data(seq, copy_num)


type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: copy_num
character(len=metadatalength) :: get_copy_meta_data

! Should have an error check for copy_num range
get_copy_meta_data = seq%copy_meta_data(copy_num)

end function get_copy_meta_data

!-------------------------------------------------
function get_qc_meta_data(seq, qc_num)


type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: qc_num
character(len=metadatalength) :: get_qc_meta_data

! Should have an error check for qc_num range
get_qc_meta_data = seq%qc_meta_data(qc_num)

end function get_qc_meta_data

!-------------------------------------------------

subroutine get_next_obs(seq, obs, next_obs, is_this_last)


type(obs_sequence_type), intent(in) :: seq
type(obs_type), intent(in) :: obs
type(obs_type), intent(out) :: next_obs
logical, intent(out) :: is_this_last

integer :: next_index

! Get index of the next observation
next_index = obs%next_time
if(next_index == -1) then
   is_this_last = .true.
   return
else
   is_this_last = .false.
   next_obs = seq%obs(next_index)
endif
!print *, 'next index = ', next_index

end subroutine get_next_obs

!-------------------------------------------------

subroutine get_prev_obs(seq, obs, prev_obs, is_this_first)


type(obs_sequence_type), intent(in) :: seq
type(obs_type), intent(in) :: obs
type(obs_type), intent(out) :: prev_obs
logical, intent(out) :: is_this_first

integer :: prev_index

! Get index of the next observation
prev_index = obs%prev_time
if(prev_index == -1) then
   is_this_first= .true.
   return
else
   is_this_first= .false.
   prev_obs = seq%obs(prev_index)
endif

end subroutine get_prev_obs

!-------------------------------------------------------------

subroutine get_obs_from_key(seq, key, obs)

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key
type(obs_type) :: obs

obs = seq%obs(key)

end subroutine get_obs_from_key

!-------------------------------------------------

subroutine get_next_obs_from_key(seq, last_key_used, next_obs, is_this_last)


type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: last_key_used
type(obs_type), intent(out) :: next_obs
logical, intent(out) :: is_this_last

integer :: next_index

! Get index of the next observation
next_index = seq%obs(last_key_used)%next_time
if(next_index == -1) then
   is_this_last = .true.
   return
else
   is_this_last = .false.
   next_obs = seq%obs(next_index)
endif

end subroutine get_next_obs_from_key

!-------------------------------------------------

subroutine get_prev_obs_from_key(seq, last_key_used, prev_obs, is_this_first)


type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: last_key_used
type(obs_type), intent(out) :: prev_obs
logical, intent(out) :: is_this_first

integer :: prev_index

! Get index of the next observation
prev_index = seq%obs(last_key_used)%prev_time
if(prev_index == -1) then
   is_this_first= .true.
   return
else
   is_this_first= .false.
   prev_obs = seq%obs(prev_index)
endif

end subroutine get_prev_obs_from_key

!-----------------------------------------------------------------

subroutine set_obs(seq, obs, key_in)

! Copies the obs into the key element of sequence where key is the key field
! in obs. If the integer argument key is present, the obs is copied into
! the key-th element of the sequence.

type(obs_sequence_type), intent(inout) :: seq
type(obs_type), intent(in) :: obs
integer, intent(in), optional :: key_in

integer :: key

! Get the key to copy into
if(present(key_in)) then 
   key = key_in
else
   key = obs%key
endif

seq%obs(key) = obs

! Make sure the key in sequence is set properly
seq%obs(key)%key = key

end subroutine set_obs

!-------------------------------------------------------------------

subroutine get_obs_time_range(seq, time1, time2, key_bounds, num_keys, out_of_range, obs)

! Add other options for getting the first time to minimize search
type(obs_sequence_type),  intent(in) :: seq
type(time_type),          intent(in) :: time1, time2
integer,                 intent(out) :: key_bounds(2)
integer,                 intent(out) :: num_keys
logical,                 intent(out) :: out_of_range
type(obs_type), intent(in), optional :: obs

type(time_type)    :: cur_time
type(obs_def_type) :: obs_def
integer            :: current, last_key

! Returns the first key and last key of sequence of obs between time1 and
! time2 along with the total number.
! A complete list of the keys can be obtained by call to get_time_range_keys
! Logical out_of_range is true if the time range is all past the end of sequence times

num_keys = 0
out_of_range = .false.

! The optional argument obs says the search can be started at this observation

! Figure out where to begin search
if(present(obs)) then
   current = obs%key
else
   current = seq%first_time
endif

! Check for all observations after the last time in the window
call get_obs_def(seq%obs(current), obs_def)
cur_time = get_obs_def_time(obs_def)
if(cur_time > time2) then
   out_of_range = .true.
   return
endif

! Find the first element in the time window
do while(current /= -1)
   call get_obs_def(seq%obs(current), obs_def)
   cur_time = get_obs_def_time(obs_def)
   if(cur_time >= time1) goto 10
   current = seq%obs(current)%next_time
end do
! Falling off the end means there are no times greater than time1
out_of_range = .true.
return

10 continue
! current is pointer to first

! First pass, count the keys for storage requirements
key_bounds(1) = current
last_key = current
do while(current /= -1)
   call get_obs_def(seq%obs(current), obs_def)
   cur_time = get_obs_def_time(obs_def)
   if(cur_time > time2) goto 20
! Found a time in the range
   num_keys = num_keys + 1
   last_key = current
   current = seq%obs(current)%next_time
end do

20 continue
key_bounds(2) = last_key

end subroutine get_obs_time_range

!---------------------------------------------------------------

subroutine get_time_range_keys(seq, key_bounds, num_keys, keys)

! Given bounds from get_obs_time_range and an array keys big enough to hold
! all the keys in the range, returns the keys in the range

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key_bounds(2), num_keys
integer, intent(out) :: keys(num_keys)

integer :: current, i

! Now loop through again to get these keys
current = key_bounds(1)
do i = 1, num_keys
   keys(i) = seq%obs(current)%key
   current = seq%obs(current)%next_time
end do

end subroutine get_time_range_keys


!-------------------------------------------------

subroutine insert_obs_in_seq(seq, obs, prev_obs)

type(obs_sequence_type), intent(inout) :: seq
type(obs_type),          intent(inout) :: obs
type(obs_type),   intent(in), optional :: prev_obs

type(time_type) :: obs_time, current_time
integer :: prev, next, current

! Inserts an observation into a sequence, optional argument
! prev_obs says that this was the predecessor in time.
! This avoids time search in cases where one is building
! a sequence from scratch.

! Make sure there is room, fail for now if not
if(seq%num_obs >= seq%max_num_obs) then
   ! Later do an increase of space and copy
   write(msgstring,*) 'ran out of room, num_obs (',seq%num_obs, &
                               ') > max_num_obs (',seq%max_num_obs,')'
   call error_handler(E_ERR,'insert_obs_in_seq',msgstring,source,revision,revdate)
endif

! Set the key for the observation
obs%key     = seq%num_obs + 1
seq%num_obs = seq%num_obs + 1

! Get the time for the observation
obs_time = get_obs_def_time(obs%def)

! Assume we're starting at the beginning.
! If we make this smarter eventually, here is where
! we'd set the initial key number for a search.

! If given an existing obs, be sure the new obs time is
! consistent - later or equal to the given previous obs. 
if(present(prev_obs)) then
   prev = prev_obs%key
   current = prev
   next = prev_obs%next_time
   
   ! it is an error to try to insert an observation after an
   ! existing obs which has a smaller timestamp.
   if (prev /= -1) then
       current_time = get_obs_def_time(seq%obs(prev)%def)
       if (obs_time < current_time) then
          !! or, do the insert searching from the start
          !prev = -1
          !current = -1
          !next = seq%first_time
          ! error out 
          write(msgstring,*) 'time of prev_obs cannot be > time of new obs'
          call error_handler(E_ERR,'insert_obs_in_seq',msgstring,source,revision,revdate)
       endif
    endif
   
    ! the insert code will search forward starting at the
    ! given obs, so it is not an error to give an obs which
    ! has a larger time than the next obs.
else
   ! Start search at beginning
   prev = -1
   current = -1
   next = seq%first_time
endif

! Have to search through the linked list to find last member
! already in with a time less than or equal to obs time
do while(next /= -1)
   prev = current
   current = next
   next = seq%obs(current)%next_time
   current_time = get_obs_def_time(seq%obs(current)%def)
! If the time of the observation in the sequence is >, stop
   if(current_time > obs_time) then 
! The observation that will follow the one being inserted is current
      next = current
      goto 10 
   endif
end do

! Falling off the end means that next is -1, so current should be previous for insertion
prev = current

! If the time check occured, previous is already pointing to previous
10 continue

! prev now holds the key of the previous observation, next holds the one after

! Link into the foward moving pointer chain
! If prev is -1, new observation goes at the start
if(prev == -1) then
   obs%next_time = seq%first_time
   obs%prev_time = -1
   seq%first_time = obs%key
else
   obs%prev_time = prev
   obs%next_time = next
   seq%obs(prev)%next_time = obs%key
endif

! Link into the backward moving pointer chain
if(next == -1) then
   obs%prev_time = seq%last_time
   obs%next_time = -1
   seq%last_time = obs%key
else
   seq%obs(next)%prev_time = obs%key
endif

! Finally, copy this obs structure into the sequence
seq%obs(obs%key) = obs

end subroutine insert_obs_in_seq

!----------------------------------------------------------------------

subroutine append_obs_to_seq(seq, obs)

! Appends an observation to an existing sequence; Error if new obs is 
! not later than time of last obs already in seq

type(obs_sequence_type), intent(inout) :: seq
type(obs_type), intent(inout) :: obs

type(obs_type) :: last_obs
type(time_type) :: obs_time, last_time

! Initialize obs_type before using
call init_obs(last_obs, 0, 0)

! If this is first, just put it in
if(.not. get_last_obs(seq, last_obs)) then
   call insert_obs_in_seq(seq, obs)
else

! Otherwise, get last obs from sequence and do insert with it as
! the previous after checking times

! Get the time for the observation
   obs_time = get_obs_def_time(obs%def)
   last_time = get_obs_def_time(last_obs%def)
   if(obs_time < last_time) then
      write(msgstring, *) 'time of appended obs cannot be < time of last obs in sequence'
      call error_handler(E_ERR,'append_obs_to_seq',msgstring,source,revision,revdate)
   endif

!!!   call insert_obs_in_seq(seq, obs)
!!!   if(1 == 1) return

! Make sure there is room, fail for now if not
   if(seq%num_obs >= seq%max_num_obs) then
! Later do an increase of space and copy
      write(msgstring,*) 'ran out of room, max_num_obs = ',seq%max_num_obs
      call error_handler(E_ERR,'append_obs_to_seq',msgstring,source,revision,revdate)
   endif

! Set the key for the observation
   obs%key = seq%num_obs + 1
   seq%num_obs = seq%num_obs + 1
! Link into the pointer chains
! Previous last points to this one, this one points back to previous last
   obs%prev_time = seq%last_time
   seq%obs(seq%last_time)%next_time = obs%key
   seq%last_time = obs%key
! Appended is at end, put a -1 for the next
   obs%next_time = -1

! Put this obs into the sequence's last slot
   seq%obs(seq%num_obs) = obs

endif

! free any space allocated at init time.
call destroy_obs(last_obs)

end subroutine append_obs_to_seq

!---------------------------------------------------------------

!subroutine insert_obs_group_in_seq(seq, obs_grp, prev_obs)

! Insert a group of observations from the same time into a sequence
!type(obs_sequence_type), intent(inout) :: seq
!type(obs_type), intent(inout) :: obs
!type(obs_type), intent(in), optional :: prev_obs
!
!end subroutine insert_obs_group_in_seq

!-------------------------------------------------

subroutine delete_obs_from_seq(seq, obs)

! Removes this observation from the sequence, does not free storage in this implementation
type(obs_sequence_type), intent(inout) :: seq
type(obs_type), intent(inout) :: obs

integer :: prev, next

prev = obs%prev_time
next = obs%next_time

!print *, 'del key, initial prev,next=', obs%key, prev, next

! update obs count??  i think this should be done, but other code
! is not prepared to deal with it.
!seq%num_obs = seq%num_obs - 1

! If only one obs, seq first_time and last_time to -1
if(prev == -1 .and. next == -1) then
  seq%first_time = -1
  seq%last_time  = -1
  return
endif

! Previous should now point to next; if deleted was first update sequence first_time
if(prev /= -1) then
   seq%obs(prev)%next_time = next
else
   seq%obs(next)%prev_time = -1
   seq%first_time = next
endif

! Next should point to previous; if deleted is last, set previous next_time to -1
if(next /= -1) then
   seq%obs(next)%prev_time = prev
else
   seq%obs(prev)%next_time = -1
   seq%last_time = prev
endif


!print *, 'prev key, next = ', prev, seq%obs(prev)%next_time
!print *, 'next key, prev = ', next, seq%obs(next)%prev_time
!print *, 'seq entire first/last = ', seq%first_time, seq%last_time

end subroutine delete_obs_from_seq

!-------------------------------------------------
subroutine set_copy_meta_data(seq, copy_num, meta_data)

! Need all sorts of error checking to avoid silly stuff eventually

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: copy_num
character(len = *),      intent(in)    :: meta_data

character(len=len(meta_data)) :: lj_meta_data ! left justified version

lj_meta_data = adjustl(meta_data)

if (len_trim(lj_meta_data) > metadatalength) then
   write(msgstring,*) 'metadata string [', trim(lj_meta_data),']'
   write(string1,*) 'must be shorter than ',metadatalength
   call error_handler(E_ERR, 'set_copy_meta_data', msgstring, &
                      source, revision, revdate, text2=string1)
endif

if (copy_num > seq%num_copies) then
   write(msgstring,*) 'trying to set copy (', copy_num, &
                      ') which is larger than num_copies (', seq%num_copies, ')'
   call error_handler(E_ERR,'set_copy_meta_data',msgstring,source,revision,revdate)
endif

seq%copy_meta_data(copy_num) = trim(lj_meta_data)

end subroutine set_copy_meta_data

!-------------------------------------------------

subroutine set_qc_meta_data(seq, qc_num, meta_data)

! Need error checks
type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: qc_num
character(len = *),      intent(in)    :: meta_data

character(len=len(meta_data)) :: lj_meta_data ! left justified version

lj_meta_data = adjustl(meta_data)

if (len_trim(lj_meta_data) > metadatalength) then
   write(msgstring,*) 'metadata string [', trim(lj_meta_data),']'
   write(string1,*) 'must be shorter than ',metadatalength
   call error_handler(E_ERR, 'set_qc_meta_data', msgstring, &
                      source, revision, revdate, text2=string1)
endif

if (qc_num > seq%num_qc) then
   write(msgstring,*) 'trying to set qc (', qc_num, &
                      ') which is larger than num_qc (', seq%num_qc, ')'
   call error_handler(E_ERR,'set_qc_meta_data',msgstring,source,revision,revdate)
endif

seq%qc_meta_data(qc_num) = trim(lj_meta_data)

end subroutine set_qc_meta_data

!-------------------------------------------------

function get_first_obs(seq, obs)

type(obs_sequence_type), intent(in) :: seq
type(obs_type),         intent(out) :: obs
logical                             :: get_first_obs

if(seq%num_obs == 0 .or. seq%first_time <= 0) then
   get_first_obs = .false.
else
   get_first_obs = .true.
   obs = seq%obs(seq%first_time)
endif

end function get_first_obs

!-------------------------------------------------

function get_last_obs(seq, obs)

type(obs_sequence_type), intent(in) :: seq
type(obs_type), intent(out) :: obs
logical :: get_last_obs

if(seq%num_obs == 0 .or. seq%last_time <=0) then
   get_last_obs = .false.
   return
else
   get_last_obs = .true.
   obs = seq%obs(seq%last_time)
endif

end function get_last_obs

!-------------------------------------------------

subroutine add_copies(seq, num_to_add)

! This requires a complete recreation of the entire obs sequence
! Add additional copies to an observation sequence. This increases
! the space for copy meta_data and goes through the whole string of
! observations deallocating and allocating (yuck), to add space.
! In the long run, may want a smoother way to do this globally.

type(obs_sequence_type), intent(inout) :: seq
integer, intent(in) :: num_to_add

character(len = metadatalength) :: meta_temp(seq%num_copies)
real(r8) :: values_temp(seq%num_copies)
integer :: i, old_num

old_num = seq%num_copies
seq%num_copies = old_num + num_to_add

! Copy the old copy metadata to temp storage, reallocate and copy
if(old_num > 0) then
   meta_temp = seq%copy_meta_data
endif

! Deallocate and reallocate with enhanced length
deallocate(seq%copy_meta_data)
allocate(seq%copy_meta_data(old_num + num_to_add))
seq%copy_meta_data(1:old_num) = meta_temp
seq%copy_meta_data(old_num+1 : old_num + num_to_add) = 'Copy metadata not initialized'

! Loop through all the observations, copy and increase size
do i = 1, seq%max_num_obs

! Copy the existing values
   if(old_num > 0) values_temp = seq%obs(i)%values

! Deallocate, reallocate and copy
   deallocate(seq%obs(i)%values)
   allocate(seq%obs(i)%values(old_num + num_to_add))
   seq%obs(i)%values(1:old_num) = values_temp
   seq%obs(i)%values(old_num+1:old_num+num_to_add) = MISSING_r8

end do

end subroutine add_copies

!-------------------------------------------------

subroutine add_qc(seq, num_to_add)

! This requires a complete recreation of the entire obs sequence
! Add additional copies to an observation sequence. This increases
! the space for copy meta_data and goes through the whole string of
! observations deallocating and allocating (yuck), to add space.
! In the long run, may want a smoother way to do this globally.

type(obs_sequence_type), intent(inout) :: seq
integer,                    intent(in) :: num_to_add

character(len = metadatalength) ::     qc_temp(seq%num_qc)
real(r8)                        :: values_temp(seq%num_qc)
integer                         :: i, old_num

old_num = seq%num_qc
seq%num_qc = old_num + num_to_add

! Copy the old copy metadata to temp storage, reallocate and copy
if(old_num > 0) then
   qc_temp = seq%qc_meta_data
endif

! Deallocate and reallocate with enhanced length
deallocate(seq%qc_meta_data)
allocate(seq%qc_meta_data(old_num + num_to_add))
seq%qc_meta_data(1:old_num) = qc_temp
seq%qc_meta_data(old_num+1 : old_num + num_to_add) = 'QC metadata not initialized'

! Loop through all the observations, copy and increase size
do i = 1, seq%max_num_obs

! Copy the existing values
   if(old_num > 0) values_temp = seq%obs(i)%qc

! Deallocate, reallocate and copy
   deallocate(seq%obs(i)%qc)
   allocate(seq%obs(i)%qc(old_num + num_to_add))
   seq%obs(i)%qc(1:old_num) = values_temp
   seq%obs(i)%qc(old_num+1:old_num+num_to_add) = 0.0_r8

end do

end subroutine add_qc

!------------------------------------------------------------------

subroutine write_obs_seq(seq, file_name)

type(obs_sequence_type), intent(in) :: seq
character(len = *),      intent(in) :: file_name

integer :: i, file_id, rc
integer :: have(max_obs_kinds)
character(len=11) :: useform


if(write_binary_obs_sequence) then
   useform = 'unformatted'
else
   useform = 'formatted'
endif

! Open the file. nsc - why is this not using open_file()?
file_id = get_unit()
write(msgstring, *) 'opening '// trim(useform) // ' file ',trim(file_name)
call error_handler(E_MSG,'write_obs_seq',msgstring)

open(unit = file_id, file = file_name, form = useform, &
     action='write', position='rewind', iostat=rc)
if (rc /= 0) then
   write(msgstring, *) 'unable to create file '//trim(file_name)
   call error_handler(E_ERR,'write_obs_seq',msgstring,source,revision,revdate)
endif

! Write the initial string for help in figuring out binary
if(write_binary_obs_sequence) then
   write(file_id) 'obs_sequence'
else
   write(file_id, *) 'obs_sequence'
endif

! Figure out which of the total possible kinds (really types) exist in this
! sequence, and set the array values to 0 for no, 1 for yes.
call set_used_kinds(seq, have)

! Write the TOC, with only the kinds that exist in this seq.
call write_obs_kind(file_id, useform, have)

! First inefficient ugly pass at writing an obs sequence, need to 
! update for storage size.  CHANGE - use num_obs for the max_num_obs, to
! limit the amount of memory needed when this sequence is read in.
if(write_binary_obs_sequence) then
   write(file_id) seq%num_copies, seq%num_qc, seq%num_obs, seq%num_obs
else
   write(file_id, *) ' num_copies: ',seq%num_copies, ' num_qc: ',     seq%num_qc
   write(file_id, *) ' num_obs: ',   seq%num_obs,    ' max_num_obs: ',seq%num_obs
endif 

do i = 1, seq%num_copies
   if(write_binary_obs_sequence) then
      write(file_id) seq%copy_meta_data(i)
   else
      write(file_id, '(a)') seq%copy_meta_data(i)
   endif
end do

do i = 1, seq%num_qc
   if(write_binary_obs_sequence) then
      write(file_id) seq%qc_meta_data(i)
   else
      write(file_id, '(a)') seq%qc_meta_data(i)
   endif
end do

if(write_binary_obs_sequence) then
   write(file_id) seq%first_time, seq%last_time
else
   write(file_id, *) ' first: ',seq%first_time, ' last: ',seq%last_time
endif

do i = 1, seq%num_obs
   if(.not. write_binary_obs_sequence) write(file_id, *) 'OBS ',seq%obs(i)%key
   call write_obs(seq%obs(i), file_id, seq%num_copies, seq%num_qc)
end do

! Close up the file
call close_file(file_id)

write(msgstring, *) 'closed file '//trim(file_name)
call error_handler(E_MSG,'write_obs_seq',msgstring)

end subroutine write_obs_seq

!------------------------------------------------------------------

subroutine read_obs_seq(file_name, add_copies, add_qc, add_obs, seq)

! Be able to increase size at read in time for efficiency

character(len = *),      intent(in)  :: file_name
integer,                 intent(in)  :: add_copies, add_qc, add_obs
type(obs_sequence_type), intent(out) :: seq

integer :: i, num_copies, num_qc, num_obs, max_num_obs, file_id, io
character(len = 16) :: label(2)
logical :: pre_I_format
character(len = 129) :: read_format

! Use read_obs_seq_header to get file format and header info
call read_obs_seq_header(file_name, num_copies, num_qc, num_obs, &
   max_num_obs, file_id, read_format, pre_I_format)

call init_obs_sequence(seq, num_copies + add_copies, &
   num_qc + add_qc, num_obs + add_obs)

! Set the number of obs available at present
seq%num_obs = num_obs

! Get the available copy_meta_data
do i = 1, num_copies
   if(read_format == 'unformatted') then
      read(file_id, iostat=io) seq%copy_meta_data(i)
   else
      read(file_id, '(a)', iostat=io) seq%copy_meta_data(i)
   endif
   if (io /= 0) then
      ! Read error of some type
      write(msgstring, *) 'Read error in copy metadata ', i, ' rc= ', io
      call error_handler(E_ERR, 'read_obs_seq', msgstring, &
         source, revision, revdate)
   endif
end do

! Get the available qc_meta_data
do i = 1, num_qc
   if(read_format == 'unformatted') then
      read(file_id, iostat=io) seq%qc_meta_data(i)
   else
      read(file_id, '(a)', iostat=io) seq%qc_meta_data(i)
   endif
   if (io /= 0) then
      ! Read error of some type
      write(msgstring, *) 'Read error in qc metadata ', i, ' rc= ', io
      call error_handler(E_ERR, 'read_obs_seq', msgstring, &
         source, revision, revdate)
   endif
end do

! Read the first and last avail_time pointers
if(read_format == 'unformatted') then
   read(file_id, iostat=io) seq%first_time, seq%last_time
else
   read(file_id, *, iostat=io) label(1),seq%first_time,label(2), seq%last_time
endif
if (io /= 0) then
   ! Read error of some type
   write(msgstring, *) 'Read error in first/last times, rc= ', io
   call error_handler(E_ERR, 'read_obs_seq', msgstring, &
      source, revision, revdate)
endif

if (seq%first_time < -1 .or. seq%first_time > max_num_obs) then
   write(msgstring, *) 'Bad value for first', seq%first_time, ', min is -1, max is ', max_num_obs 
   call error_handler(E_ERR, 'read_obs_seq', msgstring, source, revision, revdate)
endif
if (seq%last_time < -1 .or. seq%last_time > max_num_obs) then
   write(msgstring, *) 'Bad value for last', seq%last_time, ', min is -1, max is ', max_num_obs 
   call error_handler(E_ERR, 'read_obs_seq', msgstring, source, revision, revdate)
endif

! Now read in all the previously defined observations
do i = 1, num_obs
   if(.not. read_format == 'unformatted') read(file_id,*, iostat=io) label(1)
   if (io /= 0) then
      ! Read error of some type
      write(msgstring, *) 'Read error in obs label', i, ' rc= ', io
      call error_handler(E_ERR, 'read_obs_seq', msgstring, &
         source, revision, revdate)
   endif
   call read_obs(file_id, num_copies, add_copies, num_qc, add_qc, i, seq%obs(i), &
      read_format, num_obs)
! Also set the key in the obs
   seq%obs(i)%key = i
end do

! Close up the file
call close_file(file_id)

end subroutine read_obs_seq

!------------------------------------------------------------------

subroutine read_obs_seq_header(file_name, num_copies, num_qc, num_obs, &
   max_num_obs, file_id, read_format, pre_I_format, close_the_file)

! Be able to increase size at read in time for efficiency

character(len = *),     intent(in) :: file_name
integer,               intent(out) :: num_copies, num_qc, num_obs, max_num_obs, file_id
character(len = *),    intent(out) :: read_format
logical,               intent(out) :: pre_I_format
logical,      intent(in), optional :: close_the_file

character(len = 16) label(2)
character(len = 12) header
integer :: ios

! Determine the format for an obs_sequence file to be read. Options are:
! 1. Formatted, I-release format
! 2. Formatted, Hawaii-release format
! 3. Unformatted, I-release format
! 4. Unformatted, Hawaii-release format
!
! Also return the num_copies, num_qc, num_obs and max_num_obs along
! with the read_format (formatted or unformatted) and release version
! (pre_I_format).

! Have to be backwards compatible: assume new format for now
pre_I_format = .false.

! Try opening the file as formatted
file_id = get_unit()
read_format = 'formatted'
open(unit = file_id, file = file_name, form = read_format, &
   action = 'read', status = 'old', iostat = ios)
! If opening error, move to unformatted; else try to find lines
if(ios == 0) then
   ! Try to read in the I-format file header
   read(file_id, *, iostat = ios) header

   ! If read succeeds and header is 'obs_sequence' it is I-format, formatted
   if(ios == 0 .and. header == 'obs_sequence') goto 31

   ! Maybe it's formatted old H-format
   rewind file_id
   read(file_id, *, iostat = ios) label(1), num_copies, label(2), num_qc

   ! If read succeeds and label(1) is 'num_copies:' it is pre_I, formatted
   if(ios == 0 .and. label(1) == 'num_copies:') then
      pre_I_format = .true.
      ! Can read next line to extract rest of header information and return
      read(file_id, *, iostat = ios) label(1), num_obs, label(2), max_num_obs
      ! Also call read_obs_kind to initialize default obs_kind mapping
      call read_obs_kind(file_id, pre_I_format, read_format)
      goto 51
   endif

endif

! Next try unformatted open and reads
close(file_id)
read_format = 'unformatted'
pre_I_format = .false.
open(unit = file_id, file = file_name, form = read_format, &
   action = 'read', status = 'old', iostat = ios)
! If no opening error try to detect pre_i or I format, else error
if(ios == 0) then
   ! Try reading the 'obs_sequence' header
   read(file_id, iostat = ios) header
   if(header == 'obs_sequence') goto 41

   ! Maybe it's unformatted but in the old format
   rewind file_id
   read(file_id, iostat = ios) num_copies, num_qc, num_obs, max_num_obs
   ! If it's pre-i unformatted, we've read in the header and we're done
   ! Initialize the default mapping for obs_kind by calling read_obs_kind
   ! If you have a binary file on the wrong endian machine, you end up
   ! here.  the test for num_copies and num_qc will catch bad binary
   ! numbers and kick you out here -- nsc
   if(ios == 0 .and. num_copies < 1000000 .and. num_qc < 1000000) then
      pre_I_format = .true.
      call read_obs_kind(file_id, pre_I_format, read_format)
      goto 51
   endif

else
   ! Unable to figure out what to do with file or it doesn't exist
   write(msgstring, *) 'Unable to open file ', trim(file_name)
   call error_handler(E_ERR, 'read_obs_seq_header', msgstring, &
      source, revision, revdate)
endif

! Falling off the end here means file didn't correspond with any 
! expected format
write(msgstring, *) 'Unable to determine format of file ', trim(file_name)
call error_handler(E_ERR, 'read_obs_seq_header', msgstring, &
   source, revision, revdate)


! Format is I and formatted
31 continue
! Read in the obs_kind mapping table
call read_obs_kind(file_id, pre_I_format, read_format)
! Read in the rest of the header information
read(file_id, *) label(1), num_copies, label(2), num_qc
read(file_id, *) label(1), num_obs, label(2), max_num_obs
goto 51

! Format is I and unformatted
41 continue
! Read in the obs_kind mapping table
call read_obs_kind(file_id, pre_I_format, read_format)
! Read in the rest of the header information
read(file_id) num_copies, num_qc, num_obs, max_num_obs


! Ready to exit, close the file if requested by optional argument
51 continue
if(present(close_the_file)) then
   if(close_the_file) close(file_id)
endif
!write(*, *) 'pre_I is ', pre_I_format
!write(*, *) 'format is ', read_format
!write(*, *) num_copies, num_qc, num_obs, max_num_obs

end subroutine read_obs_seq_header
!-------------------------------------------------


subroutine delete_seq_head(first_time, seq, all_gone)

! Deletes all observations in the sequence with times before first_time. 
! If no observations remain, return all_gone as .true.

type(time_type),         intent(in)    :: first_time
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

type(obs_def_type)   :: obs_def
type(obs_type)       :: obs
type(time_type)      :: pre_first_time, time0, seq_start_time
integer              :: key_bounds(2), num_keys, i
integer, allocatable :: keys(:)
logical              :: out_of_range

! Initialize an observation type with appropriate size
call init_obs(obs, get_num_copies(seq), get_num_qc(seq))

! Set lowest possible time
time0 = set_time(0, 0)

! Get time of first observation in sequence; if there isn't one, return all_gone
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
   call destroy_obs(obs)
   return
else
   call get_obs_def(obs, obs_def)
   seq_start_time = get_obs_def_time(obs_def)
endif

! If first_time is lowest possible time no need to delete
if(first_time == time0) then
   all_gone = .false.
   call destroy_obs(obs)
   return
end if

! Get last possible time for observations that should NOT be used
pre_first_time = first_time - set_time(1, 0)

! Get bounds of keys in sequence that are before the first time
call get_obs_time_range(seq, time0, pre_first_time, key_bounds, num_keys, out_of_range)

! If it is out_of_range could be because all obs are after or all are before
if(out_of_range) then
   if(seq_start_time > pre_first_time) then
      ! Whole sequence is after
      all_gone = .false.
   else
      ! Whole sequence is before; but sequence is not altered?
      all_gone = .true.
   endif
   ! Destroy temp storage and return
   call destroy_obs(obs)
   return
endif

! compare num_keys with all possible keys in file; if equal, you have
! also removed all obs and should return all_gone = .true.  
if (num_keys == get_num_key_range(seq)) then
   all_gone = .true.
   ! Destroy temp storage and return
   call destroy_obs(obs)
   return
endif

! If here, then there are a set of observations that are not being used at beginning
! Delete them from the sequence
all_gone = .false.
allocate(keys(num_keys))
call get_time_range_keys(seq, key_bounds, num_keys, keys)

! Loop through the keys and delete these observations
do i = 1, num_keys
   call get_obs_from_key(seq, keys(i), obs)
   call delete_obs_from_seq(seq, obs)
end do

! Free up storage before returning
deallocate(keys)
call destroy_obs(obs)

end subroutine delete_seq_head


!-------------------------------------------------


subroutine delete_seq_tail(last_time, seq, all_gone)

! Delete all observations in the sequence with times after last_time.
! If there are none before this time return that the sequence is all_gone.

type(time_type),         intent(in)    :: last_time
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

type(obs_def_type)   :: obs_def
type(obs_type)       :: obs
type(time_type)      :: post_last_time, end_of_seq_time
integer              :: key_bounds(2), num_keys, i
integer, allocatable :: keys(:)
logical              :: out_of_range

! Initialize an observation type with appropriate size
call init_obs(obs, get_num_copies(seq), get_num_qc(seq))

! Get earliest time of observations that should be deleted
post_last_time = last_time + set_time(1, 0)

! Get time of last observation in sequence; if there are none, return all_gone
if(.not. get_last_obs(seq, obs)) then
   all_gone = .true.
   call destroy_obs(obs)
   return
endif
call get_obs_def(obs, obs_def)
end_of_seq_time = get_obs_def_time(obs_def)

! Get bounds of keys in sequence that are after the last_time
call get_obs_time_range(seq, post_last_time, end_of_seq_time, &
   key_bounds, num_keys, out_of_range)

! If it is out_of_range could be because all obs are before or all are after (none left)
if(out_of_range) then
   if(end_of_seq_time < post_last_time) then
      ! Whole sequence is after, start at beginning
      all_gone = .false.
   else
      ! Whole sequence is before
      all_gone = .true.
   endif
   ! Free storage and return
   call destroy_obs(obs)
   return
endif

! compare num_keys with all possible keys in file; if equal, you have
! also removed all obs and should return all_gone = .true.  
if (num_keys == get_num_key_range(seq)) then
   all_gone = .true.
   ! Destroy temp storage and return
   call destroy_obs(obs)
   return
endif

! If here, then there are a set of observations that are not being used at the end
! Delete them from the sequence
all_gone = .false.
allocate(keys(num_keys))
call get_time_range_keys(seq, key_bounds, num_keys, keys)

! Loop through the keys and delete these observations
do i = 1, num_keys
   call get_obs_from_key(seq, keys(i), obs)
   call delete_obs_from_seq(seq, obs)
end do

! Free storage before ending
deallocate(keys)
call destroy_obs(obs)

end subroutine delete_seq_tail


!-------------------------------------------------


subroutine delete_obs_by_typelist(num_obs_input_types, obs_input_types, &
                                  keep_list, seq, all_gone)

! Delete all observations in the sequence which either are or are not
! in the given obs types list; the sense depends on the keep flag.
! If there are no obs left afterwards return that the sequence is all_gone.

integer,                 intent(in)    :: num_obs_input_types
character(len=*),        intent(in)    :: obs_input_types(:)
logical,                 intent(in)    :: keep_list
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

type(obs_def_type)   :: obs_def
type(obs_type)       :: obs, prev_obs
integer              :: i
logical              :: is_this_last, remove_me, first_obs
integer              :: obs_type_index(num_obs_input_types), this_obs_type

! Some sanity checking on the input args.
if (num_obs_input_types <= 0) then
   write(msgstring,*) 'num_obs_input_types must be > 0'
   call error_handler(E_ERR,'delete_obs_by_typelist', msgstring, &
                      source, revision, revdate)
endif
! Ok for list to be longer; only first N items will be used.  But list
! cannot be shorter.
if (size(obs_input_types) < num_obs_input_types) then
   write(msgstring,*) 'num_obs_input_types must be >= length of list'
   call error_handler(E_ERR,'delete_obs_by_typelist', msgstring, &
                      source, revision, revdate)
endif


! Get index numbers for each type string
do i=1, num_obs_input_types
   obs_type_index(i) = get_obs_kind_index(obs_input_types(i))
   if (obs_type_index(i) < 0) then
      write(msgstring,*) 'obs_type ', trim(obs_input_types(i)), ' not found'
      call error_handler(E_ERR,'delete_obs_by_typelist', msgstring, &
                         source, revision, revdate)
   endif
enddo

! Initialize an observation type with appropriate size
call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))

! Iterate entire sequence, deleting obs which are (are not) on the list
! First, make sure there are obs to delete, and initialize first obs.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
   call destroy_obs(obs)
   call destroy_obs(prev_obs)
   return
endif

first_obs = .true.
prev_obs = obs

! This is going to be O(n*m), n=num obs in seq, m=typelist length
is_this_last = .false.
allobs : do while (.not. is_this_last)

   call get_obs_def(obs, obs_def)
   this_obs_type = get_obs_kind(obs_def)
!print *, 'this_obs_key, type = ', obs%key, this_obs_type

   ! Do we keep things on the list, or toss them?
   if (keep_list) then
      ! Assume we are going to delete the obs unless we find it on the list
      ! (can exit do loop early this way).
      remove_me = .true.
      do i=1, num_obs_input_types
         if (obs_type_index(i) == this_obs_type) then
            remove_me = .false.
            exit
         endif
      end do
   else
      ! Assume we are going to keep the obs unless we find it on the list
      ! (can exit do loop early this way).
      remove_me = .false.
      do i=1, num_obs_input_types
         if (obs_type_index(i) == this_obs_type) then
            remove_me = .true.
            exit
         endif
      end do
   endif

   ! either remove the obs and update prev, or move to next obs
   ! must be careful here; wrong order == wrong output
   if (remove_me) then
      if (first_obs) then
         call delete_obs_from_seq(seq, obs)
         if(.not. get_first_obs(seq, obs)) exit allobs
      else
         call delete_obs_from_seq(seq, obs)
         ! cannot simply use prev_obs; cached copy out of sync with seq one
         call get_next_obs_from_key(seq, prev_obs%key, obs, is_this_last)
      endif
   else
      first_obs = .false.
      prev_obs = obs
      call get_next_obs(seq, prev_obs, obs, is_this_last)
   endif
   
end do allobs

! Figure out if there are no more obs left in the sequence.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
else
   all_gone = .false.
endif

! Done.  delete temp storage and return.
call destroy_obs(obs)
call destroy_obs(prev_obs)

end subroutine delete_obs_by_typelist


!-------------------------------------------------

subroutine delete_obs_by_qc(qc_index, qc_min, qc_max, seq, all_gone)

! Delete all observations in the sequence which are outside min/max range.
! missing_r8 means infinity in that direction.
! If there are no obs left afterwards return that the sequence is all_gone.

integer,                 intent(in)    :: qc_index
real(r8),                intent(in)    :: qc_min, qc_max
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

type(obs_type)       :: obs, prev_obs
logical              :: is_this_last, remove_me, first_obs
real(r8)             :: qcval(1)

! Some sanity checking on the input args.
if (qc_index > seq%num_qc) then
   write(msgstring,*) 'qc_index must be <', seq%num_qc
   call error_handler(E_ERR,'delete_obs_by_qc', msgstring, &
                      source, revision, revdate)
endif
! Ok for min/max to be missing_r8; if both specified, min must be <= max.
if (qc_min /= missing_r8 .and. qc_max /= missing_r8 .and. qc_min > qc_max) then
   write(msgstring,*) 'qc_min must be less than or equal qc_max'
   call error_handler(E_ERR,'delete_obs_by_qc', msgstring, &
                      source, revision, revdate)
endif

! Initialize an observation type with appropriate size
call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))
 
! Iterate entire sequence, deleting obs which have a qc outside the range
! First, make sure there are obs to delete, and initialize first obs.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
   call destroy_obs(obs)
   call destroy_obs(prev_obs)
   return
endif

first_obs = .true.
prev_obs = obs

! This is going to be O(n), n=num obs in seq
is_this_last = .false.
allobs : do while (.not. is_this_last)

   call get_qc(obs, qcval, qc_index)
!print *, 'this_obs_key, qc = ', obs%key, qcval(1)

   remove_me = .false.
   if (qc_min /= missing_r8 .and. qcval(1) < qc_min) remove_me = .true.
   if (qc_max /= missing_r8 .and. qcval(1) > qc_max) remove_me = .true.

   ! either remove the obs and update prev, or move to next obs
   if (remove_me) then
      if (first_obs) then
         call delete_obs_from_seq(seq, obs)
         if(.not. get_first_obs(seq, obs)) exit allobs
      else
         call delete_obs_from_seq(seq, obs)
         ! cannot simply use prev_obs; cached copy out of sync with seq one
         call get_next_obs_from_key(seq, prev_obs%key, obs, is_this_last)
      endif
   else
      first_obs = .false.
      prev_obs = obs
      call get_next_obs(seq, prev_obs, obs, is_this_last)
   endif
   
end do allobs

! Figure out if there are no more obs left in the sequence.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
else
   all_gone = .false.
endif

! Done.  delete temp storage and return.
call destroy_obs(obs)
call destroy_obs(prev_obs)

end subroutine delete_obs_by_qc


!-------------------------------------------------

subroutine delete_obs_by_copy(copy_index, copy_min, copy_max, obs_type_name, &
                              seq, all_gone)

! Delete all observations in the sequence which are outside min/max range.
! missing_r8 means infinity in that direction.
! If there are no obs left afterwards return that the sequence is all_gone.

integer,                 intent(in)    :: copy_index
real(r8),                intent(in)    :: copy_min, copy_max
character(len=*),        intent(in)    :: obs_type_name
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

type(obs_def_type)   :: obs_def
type(obs_type)       :: obs, prev_obs
integer              :: obs_type_index, this_obs_type
logical              :: is_this_last, remove_me, first_obs
real(r8)             :: copyval(1)

! Some sanity checking on the input args.
if (copy_index > seq%num_copies) then
   write(msgstring,*) 'copy_index must be <', seq%num_copies
   call error_handler(E_ERR,'delete_obs_by_copy', msgstring, &
                      source, revision, revdate)
endif
! Ok for min/max to be missing_r8; if both specified, min must be <= max.
if (copy_min /= missing_r8 .and. copy_max /= missing_r8 .and. &
    copy_min > copy_max) then
   write(msgstring,*) 'copy_min must be less than or equal copy_max'
   call error_handler(E_ERR,'delete_obs_by_copy', msgstring, &
                      source, revision, revdate)
endif

! Get index number for the type
if (len(trim(obs_type_name)) > 0) then
   obs_type_index = get_obs_kind_index(obs_type_name)
   if (obs_type_index < 0) then
      write(msgstring,*) 'obs_type ', trim(obs_type_name), ' not found'
      call error_handler(E_ERR,'delete_obs_by_copy', msgstring, &
                         source, revision, revdate)
   endif
else
   obs_type_index = -1
endif

! Initialize an observation type with appropriate size
call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))
 
! Iterate entire sequence, deleting obs which have a copyval outside the range
! First, make sure there are obs to delete, and initialize first obs.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
   call destroy_obs(obs)
   call destroy_obs(prev_obs)
   return
endif

first_obs = .true.
prev_obs = obs

! This is going to be O(n), n=num obs in seq
is_this_last = .false.
allobs : do while (.not. is_this_last)

   call get_obs_values(obs, copyval, copy_index)
!print *, 'this_obs_key, val = ', obs%key, copyval(1)

   remove_me = .false.

   ! need to check type here, or below?
   if (obs_type_index > 0) then
      call get_obs_def(obs, obs_def)
      this_obs_type = get_obs_kind(obs_def)
      !print *, 'this_obs_key, type = ', obs%key, this_obs_type
      if (this_obs_type /= obs_type_index) remove_me = .true.
   endif

   if (copy_min /= missing_r8 .and. copyval(1) < copy_min) remove_me = .true.
   if (copy_max /= missing_r8 .and. copyval(1) > copy_max) remove_me = .true.

   ! either remove the obs and update prev, or move to next obs
   if (remove_me) then
      if (first_obs) then
         call delete_obs_from_seq(seq, obs)
         if(.not. get_first_obs(seq, obs)) exit allobs
      else
         call delete_obs_from_seq(seq, obs)
         ! cannot simply use prev_obs; cached copy out of sync with seq one
         call get_next_obs_from_key(seq, prev_obs%key, obs, is_this_last)
      endif
   else
      first_obs = .false.
      prev_obs = obs
      call get_next_obs(seq, prev_obs, obs, is_this_last)
   endif
   
end do allobs

! Figure out if there are no more obs left in the sequence.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
else
   all_gone = .false.
endif

! Done.  delete temp storage and return.
call destroy_obs(obs)
call destroy_obs(prev_obs)


end subroutine delete_obs_by_copy

!-------------------------------------------------

! To be portable between different location types (i.e. 1D, 3D sphere)
! this can only refer to the location type and the actual comparison must
! be inside the locations module itself.
subroutine select_obs_by_location(min_loc, max_loc, seq, all_gone)

! Delete all observations in the sequence which are outside the bounding box.
! If there are no obs left afterwards return that the sequence is all_gone.

type(location_type),     intent(in)    :: min_loc, max_loc
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

type(obs_def_type)   :: obs_def
type(obs_type)       :: obs, prev_obs
type(location_type)  :: location
logical              :: is_this_last, inside, first_obs


! Initialize an observation type with appropriate size
call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))

! Iterate entire sequence, deleting obs which are (are not) on the list
! First, make sure there are obs to delete, and initialize first obs.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
   call destroy_obs(obs)
   call destroy_obs(prev_obs)
   return
endif

first_obs = .true.
prev_obs = obs

! This is going to be O(n*m), n=num obs in seq, m=typelist length
is_this_last = .false.
allobs : do while (.not. is_this_last)

   call get_obs_def(obs, obs_def)
   location = get_obs_def_location(obs_def)

   ! each diff locations mod has a different one of these
   inside = is_location_in_region(location, min_loc, max_loc)
   
   ! same code as delete/keep by obstype; do any code fixes both places
   if (.not. inside) then
      if (first_obs) then
         call delete_obs_from_seq(seq, obs)
         if(.not. get_first_obs(seq, obs)) exit allobs
      else
!print *, 'going to del obs key ', obs%key
!print *, 'prev key is ', prev_obs%key
         call delete_obs_from_seq(seq, obs)
         ! cannot simply use prev_obs; cached copy out of sync with seq one
         call get_next_obs_from_key(seq, prev_obs%key, obs, is_this_last)
!print *, 'next obs now is key ', obs%key
      endif
   else
!print *, 'no del, keep this obs key ', obs%key
      first_obs = .false.
     prev_obs = obs
!print *, 'prev obs now is key ', prev_obs%key
!print *, 'obs was key ', obs%key
      call get_next_obs(seq, prev_obs, obs, is_this_last)
!print *, 'obs now is key ', obs%key
   endif
   
end do allobs

! Figure out if there are no more obs left in the sequence.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
else
   all_gone = .false.
endif

! Done.  delete temp storage and return.
call destroy_obs(obs)
call destroy_obs(prev_obs)

end subroutine select_obs_by_location


!------------------------------------------------------------------
! Figure out which of the total possible kinds (really types) exist in this
! sequence, and set the array values to 0 for no, 1 for yes.

subroutine set_used_kinds(seq, have)
type(obs_sequence_type), intent(in) :: seq
integer, intent(out) :: have(:)

integer :: i, num_copies, num_qc
integer :: num_obs
type(obs_type) :: obs
type(obs_def_type) :: obs_def
integer :: obs_kind_ind

! Get existing header info
num_copies  = get_num_copies(seq)
num_qc      = get_num_qc(seq)
num_obs     = get_num_obs(seq)

call init_obs(obs, num_copies, num_qc)

! start with no types
have(:) = 0
do i=1, num_obs
   ! cheating here, i know.  iterate the list in order the obs occur in
   ! the file, not linked list order.  i just want to know about the type
   ! of each obs, nothing about time or anything else.
   call get_obs_from_key(seq, i, obs)
   call get_obs_def(obs, obs_def)
   obs_kind_ind = get_obs_kind(obs_def)
   if (obs_kind_ind < 0) cycle   ! ignore identity obs
   have(obs_kind_ind) = 1
enddo

call destroy_obs(obs)

end subroutine set_used_kinds


!------------------------------------------------------------------
! Follow the linked list entries to copy only the linked observations
! from one sequence to the other.

subroutine copy_obs_seq(oldseq, newseq, time1, time2)

type(obs_sequence_type),   intent(in)  :: oldseq
type(obs_sequence_type),   intent(out) :: newseq
type(time_type), optional, intent(in)  :: time1, time2

integer :: i, num_copies, num_qc, max_num_obs
integer :: num_obs, num_keys, key_bounds(2)
integer, pointer :: keylist(:)
type(obs_type) :: obs
type(time_type) :: first_time, last_time
logical :: out_of_range

! Get existing header info
num_copies  = get_num_copies(oldseq)
num_qc      = get_num_qc(oldseq)
num_obs     = get_num_obs(oldseq)
max_num_obs = get_max_num_obs(oldseq)

call init_obs(obs, num_copies, num_qc)

! Really count how many obs are in the linked list, with
! optional time starts and ends.
if (present(time1)) then
   first_time = time1
else
   call get_obs_from_key(oldseq, oldseq%first_time, obs)
   first_time = get_obs_def_time(obs%def)
endif
if (present(time2)) then
   last_time = time2
else
   call get_obs_from_key(oldseq, oldseq%last_time, obs)
   last_time = get_obs_def_time(obs%def)
endif

call destroy_obs(obs)

call get_obs_time_range(oldseq, first_time, last_time, &
                        key_bounds, num_keys, out_of_range)
if (out_of_range) then
   write(msgstring, *) 'All keys out of range'
   call error_handler(E_ERR, 'copy_obs_seq', msgstring, &
                      source, revision, revdate)
endif


call init_obs_sequence(newseq, num_copies, num_qc, num_keys)

allocate(keylist(num_keys))
call get_time_range_keys(oldseq, key_bounds, num_keys, keylist)

call init_obs(obs, num_copies, num_qc)

do i=1, num_keys
   call get_obs_from_key(oldseq, keylist(i), obs)
   call set_obs(newseq, obs, i)
enddo

! Release the temp storage
deallocate(keylist)
call destroy_obs(obs)

end subroutine copy_obs_seq


!=================================================

! Functions for the obs_type
!-------------------------------------------------
subroutine init_obs(obs, num_copies, num_qc)

! Sort of a constructor for obs_type
! Should this be public or private just for sequence?

integer, intent(in) :: num_copies, num_qc
type(obs_type), intent(out) :: obs

! Intentionally allocate even 0 copies.  This creates an 
! associated pointer with an explicit size of zero.
allocate(obs%values(num_copies))
if (num_copies > 0) obs%values = missing_r8

allocate(obs%qc(num_qc))
if (num_qc > 0) obs%qc = 0.0_r8

obs%key = -1
obs%prev_time = -1
obs%next_time = -1
obs%cov_group = -1

end subroutine init_obs

!-----------------------------------------------------

subroutine copy_obs(obs1, obs2)

! This routine is overloaded with the = operator

!type(obs_type), intent(out) :: obs1
type(obs_type), intent(inout) :: obs1
type(obs_type), intent(in) :: obs2

obs1%key = obs2%key
call copy_obs_def(obs1%def, obs2%def)

if (associated(obs1%values)) then
   if (size(obs1%values) /= size(obs2%values)) then
      deallocate(obs1%values)
      allocate(obs1%values(size(obs2%values)))
   endif
else
   allocate(obs1%values(size(obs2%values)))
endif

if (associated(obs1%qc)) then
   if (size(obs1%qc) /= size(obs2%qc)) then
      deallocate(obs1%qc)
      allocate(obs1%qc(size(obs2%qc)))
   endif
else
   allocate(obs1%qc(size(obs2%qc)))
endif

obs1%values = obs2%values
obs1%qc = obs2%qc
obs1%prev_time = obs2%prev_time
obs1%next_time = obs2%next_time
obs1%cov_group = obs2%cov_group

end subroutine copy_obs

!-------------------------------------------------

subroutine destroy_obs(obs)

! Free up allocated storage in an observation type
type(obs_type), intent(inout) :: obs

if (associated(obs%values)) then
   deallocate(obs%values)
   nullify(obs%values)
endif
if (associated(obs%qc)) then
   deallocate(obs%qc)
   nullify(obs%qc)
endif
!if pointers are nullified() then this is safe (and simpler).
!deallocate(obs%values, obs%qc)
call destroy_obs_def(obs%def)

end subroutine destroy_obs

!-----------------------------------------------------

subroutine copy_partial_obs(obs1, obs2, numcopies, copylist, &
                            numqc, qclist)

! Copy from obs2 to obs1, the entire contents of the
! obs def, but only the copies and qcs as listed (in order)
! Special value (0) means leave space but there is
! no existing value to copy.

type(obs_type), intent(inout) :: obs1
type(obs_type), intent(in) :: obs2
integer, intent(in) :: numcopies, copylist(:), numqc, qclist(:)

integer :: i, ival

! only basic idiotproofing - detect bad indices in the lists
! without too much expense in time.  no checks here that length
! of lists are >= num sizes.

! numcopies and numqc are the new outgoing sizes in obs1.
! check the index lists to be sure they are >= 0 and <= size
! of existing data in obs2.  
ival = min(minval(copylist(1:numcopies)), minval(qclist(1:numqc)))
if (ival < 0) then
   write(msgstring, '(A,I8,A)') 'index list value, ', ival, ' must be >= 0'
   call error_handler(E_ERR, 'copy_partial_obs:', msgstring, &
               source, revision, revdate)
endif
ival = maxval(copylist(1:numcopies))
if (ival > size(obs2%values)) then
   write(msgstring, '(A,I8,A,I8)') 'index list value, ', ival, &
      ' is larger than copies length, ', size(obs2%values)
   call error_handler(E_ERR, 'copy_partial_obs:', msgstring, &
               source, revision, revdate)
endif
ival = maxval(qclist(1:numqc))
if (ival > size(obs2%qc)) then
   write(msgstring, '(A,I8,A,I8)') 'index list value, ', ival, &
      ' is larger than qc length, ', size(obs2%qc)
   call error_handler(E_ERR, 'copy_partial_obs:', msgstring, &
               source, revision, revdate)
endif

obs1%key = obs2%key
call copy_obs_def(obs1%def, obs2%def)

if (associated(obs1%values)) then
   if (size(obs1%values) /= numcopies) then
      deallocate(obs1%values)
      allocate(obs1%values(numcopies))
   endif
else
   allocate(obs1%values(numcopies))
endif

if (associated(obs1%qc)) then
   if (size(obs1%qc) /= numqc) then
      deallocate(obs1%qc)
      allocate(obs1%qc(numqc))
   endif
else
   allocate(obs1%qc(numqc))
endif

do i = 1, numcopies
   if (copylist(i) == 0) then
       obs1%values(i) = MISSING_R8
   else
       obs1%values(i) = obs2%values(copylist(i))
   endif
enddo
do i = 1, numqc
   if (qclist(i) == 0) then
      obs1%qc(i) = 0.0_r8
   else
      obs1%qc(i) = obs2%qc(qclist(i))
   endif
enddo

obs1%prev_time = obs2%prev_time
obs1%next_time = obs2%next_time
obs1%cov_group = obs2%cov_group

end subroutine copy_partial_obs

!-------------------------------------------------
subroutine get_obs_def(obs, obs_def)

type(obs_type), intent(in) :: obs
type(obs_def_type), intent(out) :: obs_def

! WARNING: NEED TO DEFINE A COPY ROUTINE FOR OBS_DEF !!!
call copy_obs_def(obs_def, obs%def)

end subroutine get_obs_def

!-------------------------------------------------
subroutine set_obs_def(obs, obs_def)

type(obs_type), intent(inout) :: obs
type(obs_def_type), intent(in) :: obs_def

call copy_obs_def(obs%def, obs_def)

end subroutine set_obs_def

!-------------------------------------------------

subroutine get_obs_values(obs, values, copy_indx)


type(obs_type),    intent(in)  :: obs
real(r8),          intent(out) :: values(:)
integer, optional, intent(in)  :: copy_indx

if(present(copy_indx)) then
   values(1) = obs%values(copy_indx)
else
   values = obs%values
endif

end subroutine get_obs_values

!-------------------------------------------------

subroutine set_obs_values(obs, values, copy_indx)

type(obs_type), intent(inout) :: obs
real(r8), intent(in) :: values(:)
integer, optional, intent(in) :: copy_indx

if(present(copy_indx)) then
   obs%values(copy_indx) = values(1)
else
   obs%values = values
endif

end subroutine set_obs_values

!-------------------------------------------------

subroutine replace_obs_values(seq, key, values, copy_indx)

type(obs_sequence_type), intent(inout) :: seq
integer, intent(in) :: key
real(r8), intent(in) :: values(:)
integer, optional, intent(in) :: copy_indx

if(present(copy_indx)) then
   seq%obs(key)%values(copy_indx) = values(1)
else
   seq%obs(key)%values = values
endif

end subroutine replace_obs_values

!-------------------------------------------------
subroutine get_qc(obs, qc, qc_indx)


type(obs_type),    intent(in) :: obs
real(r8),         intent(out) :: qc(:)
integer, optional, intent(in) :: qc_indx

if(present(qc_indx)) then
   qc(1) = obs%qc(qc_indx)
else
   qc = obs%qc
endif

end subroutine get_qc

!-------------------------------------------------
subroutine set_qc(obs, qc, qc_indx)

type(obs_type),   intent(inout) :: obs
real(r8),          intent(in) :: qc(:)
integer, optional, intent(in) :: qc_indx

if(present(qc_indx)) then
   obs%qc(qc_indx) = qc(1)
else
   obs%qc = qc
endif

end subroutine set_qc

!-------------------------------------------------

subroutine replace_qc(seq, key, qc, qc_indx)

type(obs_sequence_type), intent(inout) :: seq
integer, intent(in) :: key
real(r8), intent(in) :: qc(:)
integer, optional, intent(in) :: qc_indx

if(present(qc_indx)) then
   seq%obs(key)%qc(qc_indx) = qc(1)
else
   seq%obs(key)%qc = qc
endif

end subroutine replace_qc

!-------------------------------------------------------------

function get_obs_key(obs)

type(obs_type), intent(in) :: obs
integer :: get_obs_key

get_obs_key = obs%key

end function get_obs_key

!-------------------------------------------------

subroutine write_obs(obs, file_id, num_copies, num_qc)

! Write out an observation to file, inefficient

type(obs_type), intent(in) :: obs
integer,        intent(in) :: file_id, num_copies, num_qc

integer :: i

do i = 1, num_copies
   if(write_binary_obs_sequence) then
      write(file_id) obs%values(i)
   else
      write(file_id, *) obs%values(i)
   endif
end do

do i = 1, num_qc
   if(write_binary_obs_sequence) then
      write(file_id) obs%qc(i)
   else
      write(file_id, *) obs%qc(i)
   endif
end do

if(write_binary_obs_sequence) then
   write(file_id) obs%prev_time, obs%next_time, obs%cov_group
   call write_obs_def(file_id, obs%def, obs%key, 'unformatted')
else
   write(file_id, *) obs%prev_time, obs%next_time, obs%cov_group
   call write_obs_def(file_id, obs%def, obs%key)
endif

end subroutine write_obs

!-------------------------------------------------

subroutine read_obs(file_id, num_copies, add_copies, num_qc, add_qc, key, &
                    obs, read_format, max_obs)

! Read in observation from file, watch for allocation of storage
! This RELIES on the fact that obs%values(1) is ALWAYS the observation value
! (as opposed to the prior or mean or ...)
!
! Are the checks for num_copies == 0 or <0 necessary? 
! Yes, they happen in create_fixed_network_sequence

integer,              intent(in)    :: file_id, num_copies, add_copies
integer,              intent(in)    :: num_qc, add_qc, key
character(len = *),   intent(in)    :: read_format
type(obs_type),       intent(inout) :: obs
integer, optional,    intent(in)    :: max_obs

integer  :: i, io
real(r8) :: temp_val

! Read in values and qc
if(num_copies > 0) then
   if(read_format == 'unformatted') then
      do i = 1, num_copies
         read(file_id, iostat=io) obs%values(i)
         if (io /= 0) then
            ! Read error of some type
            write(msgstring, *) 'Read error in obs values, obs ', i, ' rc= ', io
            call error_handler(E_ERR, 'read_obs', msgstring, &
               source, revision, revdate)
         endif
      end do
   else
      read(file_id, *, iostat=io) obs%values(1:num_copies)
      if (io /= 0) then
         ! Read error of some type
         write(msgstring, *) 'Read error in obs values, rc= ', io
         call error_handler(E_ERR, 'read_obs', msgstring, &
            source, revision, revdate)
      endif
   endif
endif

if(num_qc > 0) then
   if(read_format == 'unformatted') then
      do i = 1, num_qc
         read(file_id, iostat=io) obs%qc(i)
         if (io /= 0) then
            ! Read error of some type
            write(msgstring, *) 'Read error in qc values, obs ', i, ' rc= ', io
            call error_handler(E_ERR, 'read_obs', msgstring, &
               source, revision, revdate)
         endif
      end do
   else
      read(file_id, *, iostat=io) obs%qc(1:num_qc)
      if (io /= 0) then
         ! Read error of some type
         write(msgstring, *) 'Read error in qc values, rc= ', io
         call error_handler(E_ERR, 'read_obs', msgstring, &
            source, revision, revdate)
      endif
   endif
endif

! Need to pass the value if available
if(num_copies > 0) then
   temp_val = obs%values(1)
else
   temp_val = missing_r8
endif 

! Read in linked list pointers and error check
if(read_format == 'unformatted') then
   read(file_id, iostat=io) obs%prev_time, obs%next_time, obs%cov_group
else
   read(file_id, *, iostat=io) obs%prev_time, obs%next_time, obs%cov_group
endif
if (io /= 0) then
   ! Read error of some type
   write(msgstring, *) 'Read error in linked list or cov grp, rc= ', io
   call error_handler(E_ERR, 'read_obs', msgstring, &
      source, revision, revdate)
endif

! if max_obs specified, do additional error checking
if (present(max_obs)) then
   ! -1 is ok; used for first and last entries.
   if (obs%prev_time < -1 .or. obs%prev_time > max_obs) then
      write(msgstring, *) 'Bad value for previous obs, ', obs%prev_time, ', in obs ', key 
      call error_handler(E_ERR, 'read_obs', msgstring, &
         source, revision, revdate)
   endif
   if (obs%next_time < -1 .or. obs%next_time > max_obs) then
      write(msgstring, *) 'Bad value for next obs, ', obs%next_time, ', in obs ', key
      call error_handler(E_ERR, 'read_obs', msgstring, &
         source, revision, revdate)
   endif
endif

! Get model-dependent values
if(read_format == 'unformatted') then
   call read_obs_def(file_id, obs%def, key, temp_val, 'unformatted')
else
   call read_obs_def(file_id, obs%def, key, temp_val)
endif

! Copy the temp_val back to obs%values(1) if there are copies of data
if(num_copies > 0) obs%values(1) = temp_val

end subroutine read_obs

!------------------------------------------------------------------------------

subroutine interactive_obs(num_copies, num_qc, obs, key)

integer,           intent(in) :: num_copies, num_qc, key
type(obs_type), intent(inout) :: obs

integer :: i

! Does interactive initialization of an observation type

call interactive_obs_def(obs%def, key)
do i = 1, num_copies
   write(*, *) 'Enter value ', i, 'for this observation'
   read(*, *) obs%values(i)
end do

do i = 1, num_qc
   write(*, *) 'Enter quality control value ', i, 'for this observation'
   read(*, *) obs%qc(i)
end do

! WHAT ABOUT THE COVARIANCE GROUPING???

end subroutine interactive_obs


!---------------------------------------------------------

function get_num_times(seq)

! Returns number of different times for observations in sequence
! Could also be computed as sequence is built?

type(obs_sequence_type), intent(in) :: seq
integer :: get_num_times

integer :: next
type(obs_def_type) :: obs_def
type(time_type) :: this_time, prev_time

! Just loop through the time sorted sequence and look for different times
get_num_times = 0
next = seq%first_time

do while (next /= -1)
   call get_obs_def(seq%obs(next), obs_def)
   this_time = get_obs_def_time(obs_def)
   if(get_num_times == 0) then
      get_num_times = 1
   else if(this_time /= prev_time) then
      get_num_times = get_num_times + 1
   endif
   prev_time = this_time
   next = seq%obs(next)%next_time
end do

end function get_num_times

!---------------------------------------------------------

function get_num_key_range(seq, key1, key2)

! Returns number of observations between the two given keys

type(obs_sequence_type), intent(in) :: seq
integer, optional, intent(in) :: key1, key2
integer :: get_num_key_range

integer :: next, last


if (present(key1)) then
   if (key1 < seq%first_time .or. key1 > seq%last_time) then
      write(msgstring, *) 'Bad value for key1, must be between ', &
                            seq%first_time, ' and ', seq%last_time
      call error_handler(E_ERR, 'get_num_key_range', msgstring, &
         source, revision, revdate)
   endif
   next = key1
else
   next = seq%first_time
endif
if (present(key2)) then
   if (key2 < seq%first_time .or. key2 > seq%last_time) then
      write(msgstring, *) 'Bad value for key2, must be between ', &
                            seq%first_time, ' and ', seq%last_time
      call error_handler(E_ERR, 'get_num_key_range', msgstring, &
         source, revision, revdate)
   endif
   last = key2
else
   last = seq%last_time
endif

! count them up
get_num_key_range = 0
do while (next /= -1)
   get_num_key_range = get_num_key_range + 1
   if (next == last) exit
   next = seq%obs(next)%next_time
end do

end function get_num_key_range


!-------------------------------------------------
!subroutine get_cov_group
!-------------------------------------------------
!subroutine set_cov_group ???

!=================================================


end module obs_sequence_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
