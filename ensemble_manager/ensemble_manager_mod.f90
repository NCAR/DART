! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module ensemble_manager_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Manages a data structure that is composed of copies of a vector of variables.
! The general data structure simply represents this two-dimensional array but
! provides tools for storing it with one dimension (copies) or the other 
! (variables) complete on each process of a multiprocess implementation. Some
! operations that are specifically aimed at supporting ensemble filter applications
! have been placed here for efficiency even though they might be more 
! appropriately abstracted at a higher level of code.

use types_mod,         only : r8, MISSING_R8
use utilities_mod,     only : register_module, do_nml_file, do_nml_term, &
                              error_handler, E_ERR, E_MSG, do_output, &
                              nmlfileunit, find_namelist_in_file,        &
                              check_namelist_read
use assim_model_mod,   only : aread_state_restart, awrite_state_restart, &
                              open_restart_read, open_restart_write,     &
                              close_restart, pert_model_state
use time_manager_mod,  only : time_type, set_time
use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian
use mpi_utilities_mod, only : task_count, my_task_id, send_to, receive_from

implicit none
private

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!PAR other storage option control can be implemented here. In particular, want to find
!PAR some way, either allocating or multiple addressing, to use same chunk of storage
!PAR for both copy and var complete representations.

public :: init_ensemble_manager,      end_ensemble_manager,     get_ensemble_time,      &
          ensemble_type,              duplicate_ens,            get_var_owner_index,    &
          get_my_num_copies,          get_my_copies,            get_my_num_vars,        &
          get_my_vars,                compute_copy_mean,        compute_copy_mean_sd,   &
          get_copy,                   put_copy,                 all_vars_to_all_copies, &
          all_copies_to_all_vars,     read_ensemble_restart,    write_ensemble_restart, &
          compute_copy_mean_var,      get_copy_owner_index

type ensemble_type
   !DIRECT ACCESS INTO STORAGE IS USED TO REDUCE COPYING: BE CAREFUL
   !!!private
   integer                  :: num_copies, num_vars, my_num_copies, my_num_vars
   integer,         pointer :: my_copies(:), my_vars(:)
   ! Storage in next line is to be used when each pe has all copies of subset of vars
   real(r8),        pointer :: copies(:, :)         ! Dimensioned (num_copies, my_num_vars)
   ! Storage on next line is used when each pe has subset of copies of all vars
   real(r8),        pointer :: vars(:, :)           ! Dimensioned (num_vars, my_num_copies)
   ! Time is only related to var complete
   type(time_type), pointer :: time(:)
   integer                  :: distribution_type
end type ensemble_type

! Logical flag for initialization of module
logical              :: module_initialized = .false.

! Module storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! Module storage for writing error messages
character(len = 129) :: errstring

! Module storage for pe information for this process avoids recomputation
integer              :: my_pe, num_pes

!-----------------------------------------------------------------
!
! namelist with default values

!PAR: Might want this to come from above now so multiple ensembles could have
! different restart types.
! If true a single restart file containing all ensemble vars is read/written
! If false, one restart file is used for each ensemble member
logical  :: single_restart_file_in  = .true.
logical  :: single_restart_file_out = .true.
! Size of perturbations for creating ensembles when model won't do it
real(r8) :: perturbation_amplitude  = 0.2_r8

namelist / ensemble_manager_nml / single_restart_file_in,  &
                                  single_restart_file_out, &
                                  perturbation_amplitude

!-----------------------------------------------------------------

contains

!-----------------------------------------------------------------

subroutine init_ensemble_manager(ens_handle, num_copies, &
   num_vars, distribution_type_in)

type(ensemble_type), intent(out)            :: ens_handle
integer,             intent(in)             :: num_copies, num_vars
integer,             intent(in), optional   :: distribution_type_in

integer :: iunit, io

! Distribution type controls pe layout of storage; Default is 1. 1 is only one implemented.
if(.not. present(distribution_type_in)) then
   ens_handle%distribution_type = 1
else
   ! Check for error: only type 1 implemented for now
   if(distribution_type_in /= 1) call error_handler(E_ERR, 'init_ensemble_manager', &
      'only distribution_type 1 is implemented', source, revision, revdate)
   ens_handle%distribution_type = distribution_type_in
endif

! First call to init_ensemble_manager must initialize module and read namelist
if ( .not. module_initialized ) then
   ! Initialize the module with utilities 
   call register_module(source, revision, revdate)
   module_initialized = .true.

   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "ensemble_manager_nml", iunit)
   read(iunit, nml = ensemble_manager_nml, iostat = io)
   call check_namelist_read(iunit, io, "ensemble_manager_nml")

   if (do_nml_file()) write(nmlfileunit, nml=ensemble_manager_nml)
   if (do_nml_term()) write(     *     , nml=ensemble_manager_nml)

   ! Get mpi information for this process; it's stored in module storage
   num_pes = task_count()
   my_pe = my_task_id()
endif

! Set the global storage bounds for the number of copies and variables
ens_handle%num_copies = num_copies
ens_handle%num_vars = num_vars

! Figure out how the ensemble copies are partitioned
call set_up_ens_distribution(ens_handle)

end subroutine init_ensemble_manager

!-------------------------------------------------------------------------------

subroutine read_ensemble_restart(ens_handle, start_copy, end_copy, &
   start_from_restart, file_name, init_time, force_single_file)

type(ensemble_type),  intent(inout)           :: ens_handle
integer,              intent(in)              :: start_copy, end_copy
logical,              intent(in)              :: start_from_restart
character(len = *),   intent(in)              :: file_name
type(time_type),      intent(in),    optional :: init_time
logical,              intent(in),    optional :: force_single_file

! The ensemble being read from restart is stored in copies start_copy:end_copy
! contiguously (by fiat). So if start_copy is 6, the first ensemble restart
! goes in copy 6, the second in copy 7, etc. This lets higher level determine
! where other copies like mean, inflation, etc., are stored.

! Would like to avoid num_vars size storage
! LARGE ARRAY ON STACK -- make this a pointer and allocate it from heap?
real(r8)                            :: ens(ens_handle%num_vars)
integer                             :: iunit, i, j
character(len = LEN(file_name) + 5) :: this_file_name
character(len = 4)                  :: extension
type(time_type)                     :: ens_time
integer                             :: global_copy_index
logical                             :: interf_provided
logical                             :: single_file_override

! Does not make sense to have start_from_restart and single_restart_file_in BOTH false
if(.not. start_from_restart .and. .not. single_restart_file_in) then
   write(errstring, *) 'start_from_restart in filter_nml and single_restart_file_in in &
      &ensemble_manager_nml cannot both be false'
   call error_handler(E_ERR,'read_ensemble_restart', errstring, source, revision, revdate)
endif


! Some compilers (absoft, but others also) are particularly unhappy about
! both checking present(N) _and_ evaluating N inside a single if() test.
! (It evaluates both at the same time and blows up on the not present value.)
! The standard says they do not have to evaluate right to left.  Common error
! for anyone with a C programming background.   So -- set a separate local
! logical variable which always has a value, whether the arg is present or not.
if (present(force_single_file)) then
  single_file_override = force_single_file
else
  single_file_override = .false.
endif

!-------- Block for single restart file or single member  being perturbed -----
if(single_restart_file_in .or. .not. start_from_restart .or. &
   single_file_override) then 
   ! Single restart file is read only by master_pe and then distributed
   if(my_pe == 0) iunit = open_restart_read(file_name)

   ! Loop through the total number of copies
   do i = start_copy, end_copy
      ! Only master_pe does reading. Everybody can do their own perturbing
      if(my_pe == 0) then
         ! Read restarts in sequence; only read once for not start_from_restart
         if(start_from_restart .or. i == start_copy) &
            call aread_state_restart(ens_time, ens, iunit)
         ! Override this time if requested by namelist
         if(present(init_time)) ens_time = init_time
      endif

      ! Store this copy in the appropriate place on the appropriate process
      call put_copy(0, ens_handle, i, ens, ens_time)

   end do

   ! Master pe must close the file it's been reading
   if(my_pe == 0) call close_restart(iunit)

else
!----------- Block that follows is for multiple restart files -----------
   ! Loop to read in all my ensemble members
   READ_MULTIPLE_RESTARTS: do i = 1, ens_handle%my_num_copies
      ! Get global index for my ith ensemble
      global_copy_index = ens_handle%my_copies(i)
      ! If this global copy is in the range being read in, proceed
      if(global_copy_index >= start_copy .and. global_copy_index <= end_copy) then
         ! File name extension is the global index number for the copy
         write(extension, '(i4.4)') global_copy_index - start_copy + 1
         this_file_name = trim(file_name) // '.' // extension
         iunit = open_restart_read(this_file_name)

         ! Read the file directly into storage
         call aread_state_restart(ens_handle%time(i), ens_handle%vars(:, i), iunit)
         if(present(init_time)) ens_handle%time(i) = init_time
      
         ! Close the restart file
         call close_restart(iunit)
      endif
   end do READ_MULTIPLE_RESTARTS
endif

!------------------- Block that follows perturbs single base restart --------
if(.not. start_from_restart) then
   PERTURB_MY_RESTARTS: do i = 1, ens_handle%my_num_copies
      global_copy_index = ens_handle%my_copies(i)
      ! If this is one of the actual ensemble copies, then perturb
      if(global_copy_index >= start_copy .and. global_copy_index <= end_copy) then
         ! See if model has an interface to perturb
         call pert_model_state(ens_handle%vars(:, i), ens_handle%vars(:, i), interf_provided)
         ! If model does not provide a perturbing interface, do it here
         if(.not. interf_provided) then
            ! To reproduce for varying pe count, need  fixed sequence for each copy
            call init_random_seq(random_seq, global_copy_index)
            do j = 1, ens_handle%num_vars
               ens_handle%vars(j, i) = random_gaussian(random_seq, ens_handle%vars(j, i), &
                  perturbation_amplitude)
            end do
         endif
      endif
   end do PERTURB_MY_RESTARTS
endif

end subroutine read_ensemble_restart

!-----------------------------------------------------------------

subroutine write_ensemble_restart(ens_handle, file_name, start_copy, end_copy, &
   force_single_file)

type(ensemble_type),  intent(inout) :: ens_handle
character(len = *),   intent(in)    :: file_name
integer,              intent(in)    :: start_copy, end_copy
logical, optional,    intent(in)    :: force_single_file

! Large temporary storage to be avoided if possible
! LARGE ARRAY ON STACK -- make this a pointer and allocate it from heap?
real(r8)                            :: ens(ens_handle%num_vars)
type(time_type)                     :: ens_time
integer                             :: iunit, i, global_index
integer                             :: owner, owners_index
character(len = LEN(file_name) + 5) :: this_file_name
character(len = 4)                  :: extension
logical                             :: single_file_forced

if (present(force_single_file) ) then
        single_file_forced = force_single_file
else
        single_file_forced = .FALSE.
endif

! For single file, need to send restarts to pe0 and it writes them out.
!-------------- Block for single_restart file -------------
! Need to force single restart file for inflation files
if(single_restart_file_out .or. single_file_forced) then

   ! Single restart file is written only by the master_pe
   if(my_pe == 0) then
      iunit = open_restart_write(file_name)
      
      ! Loop to write each ensemble member
      do i = start_copy, end_copy
         ! Figure out where this ensemble member is being stored
         call get_copy_owner_index(i, owner, owners_index)
         ! If it's on the master pe, just write it
         if(owner == 0) then
            call awrite_state_restart(ens_handle%time(owners_index), &
               ens_handle%vars(:, owners_index), iunit)
         else
            ! Get the ensemble from the owner and write it out
            ! This communication assumes index numbers are monotonically increasing
            ! and that communications is blocking so that there are not multiple 
            ! outstanding messages from the same processors (also see send_to below).
            call receive_from(owner, ens, ens_time)
            call awrite_state_restart(ens_time, ens, iunit)
         endif
      end do
      call close_restart(iunit)
   else
      ! If I'm not the master pe with a single restart, I must send my copies
      ! to master pe for writing to file
      do i = 1, ens_handle%my_num_copies
         ! Figure out which global index this is
         global_index = ens_handle%my_copies(i)
         ! Ship this ensemble off to the master
         if(global_index >= start_copy .and. global_index <= end_copy) &
            call send_to(0, ens_handle%vars(:, i), ens_handle%time(i))
      end do
   endif

else
!-------------- Block for multiple restart files -------------
   ! Everyone can just write their own files
   do i = 1, ens_handle%my_num_copies
      ! Figure out which global index this is
      global_index = ens_handle%my_copies(i)
      if(global_index >= start_copy .and. global_index <= end_copy) then
         write(extension, '(i4.4)') ens_handle%my_copies(i)
         this_file_name = trim(file_name) // '.' // extension
         iunit = open_restart_write(this_file_name)
         call awrite_state_restart(ens_handle%time(i), ens_handle%vars(:, i), iunit)
         call close_restart(iunit)
      endif
   end do
endif

end subroutine write_ensemble_restart

!-----------------------------------------------------------------

subroutine get_copy(receiving_pe, ens_handle, copy, vars, mtime)

! The requested copy of the vars from the ensemble is stored into
! vars on the receiving_pe process. Only the receiving pe and the pe
! that actually stores the copy do anything. The vars storage only need be 
! allocated on the receiving_pe. Time only transferred if mtime present.
! An example application of this routine would be for writing state 
! space diagnostics. The filter only uses a single writer to an open file,
! currently process 0. To write out the ensemble mean state to the 
! prior and posterior netcdf state space diagnostic files, process 0
! would need to request and receive the copy containing the ensemble
! mean from whatever process stores it when processes have a subset
! of copies of all variables.

integer,             intent(in)              :: receiving_pe
type(ensemble_type), intent(in)              :: ens_handle
integer,             intent(in)              :: copy
real(r8),            intent(out)             :: vars(:)
type(time_type),     intent(out),  optional  :: mtime

integer :: owner, owners_index

! Verify that requested copy exists
if(copy < 1 .or. copy > ens_handle%num_copies) then
   write(errstring, *) 'Requested copy: ', copy, ' is > maximum copy: ', ens_handle%num_copies 
   call error_handler(E_ERR,'get_copy', errstring, source, revision, revdate)
endif

! Make sure that vars has enough space to handle the answer
if(size(vars) < ens_handle%num_vars) then
   write(errstring, *) 'Size of vars: ', size(vars), ' Must be at least ', ens_handle%num_vars
   call error_handler(E_ERR,'get_copy', errstring, source, revision, revdate)
endif

! Figure out which PE stores this copy and what its local storage index is
call get_copy_owner_index(copy, owner, owners_index)

!----------- Block of code that must be done by receiving pe -----------------------------
if(my_pe == receiving_pe) then
   ! If PE that stores is the same, just copy and return
   if(my_pe == owner) then
      vars = ens_handle%vars(:, owners_index)
      if(present(mtime)) mtime = ens_handle%time(owners_index)
      ! If I'm the receiving PE and also the owner, I'm all finished; return
      return
   endif
 
   ! Otherwise, must wait to receive vars and time from storing pe
      call receive_from(owner, vars, mtime)
endif

!----- Block of code that must be done by PE that stores the copy IF it is NOT receiver -----
if(my_pe == owner) then
   ! Send copy to receiving pe
   if(present(mtime)) then
      call send_to(receiving_pe, ens_handle%vars(:, owners_index), ens_handle%time(owners_index))
   else
      call send_to(receiving_pe, ens_handle%vars(:, owners_index))
   endif
endif
!------ End of block ---------------------------------------------------------------------

end subroutine get_copy

!-----------------------------------------------------------------

subroutine put_copy(sending_pe, ens_handle, copy, vars, mtime)

! The vars on the sending_pe is stored on the pe
! that stores this copy. Only the sending_pe and storing pe do anything. The
! vars storage only needs to be allocated on the sending_pe. An example
! use is when process 0 in the filter reads the observation keys from the 
! obs_sequence file and then needs to place this as the appropriate copy
! in the obs_ensemble file. This copy is stored by some process that may
! not be 0.

integer,             intent(in)           :: sending_pe
type(ensemble_type), intent(inout)        :: ens_handle
integer,             intent(in)           :: copy
real(r8),            intent(in)           :: vars(:)
type(time_type),     intent(in), optional :: mtime

integer :: owner, owners_index

if(copy < 1 .or. copy > ens_handle%num_copies) then
   write(errstring, *) 'Requested copy: ', copy, ' is > maximum copy: ', ens_handle%num_copies 
   call error_handler(E_ERR,'put_copy', errstring, source, revision, revdate)
endif

! Make sure that vars has enough space to handle the answer
if(size(vars) < ens_handle%num_vars) then
   write(errstring, *) 'Size of vars: ', size(vars), ' Must be at least ', ens_handle%num_vars
   call error_handler(E_ERR,'put_copy', errstring, source, revision, revdate)
endif

! What PE stores this copy and what is its local storage index
call get_copy_owner_index(copy, owner, owners_index)

! Block of code that must be done by PE that is to send the copy
if(my_pe == sending_pe) then
   ! If PE that stores is the same, just copy and return
   if(my_pe == owner) then
      ens_handle%vars(:, owners_index) = vars
      if(present(mtime)) ens_handle%time(owners_index) = mtime
      ! If I'm the sending PE and also the owner, I'm all finished; return
      return
   endif
 
   ! Otherwise, must send vars and possibly time to storing pe
      call send_to(owner, vars, mtime)
endif

! Block of code that must be done by PE that stores the copy IF it is NOT sender
if(my_pe == owner) then
   ! Need to receive copy from sending_pe
   if(present(mtime)) then
      call receive_from(sending_pe, ens_handle%vars(:, owners_index), ens_handle%time(owners_index))
   else
      call receive_from(sending_pe, ens_handle%vars(:, owners_index))
   endif
endif

end subroutine put_copy

!-----------------------------------------------------------------

subroutine get_ensemble_time(ens_handle, indx, mtime)

! Returns the time of an ensemble member indexed by local storage on this pe.

type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: indx
type(time_type),     intent(out) :: mtime

if(indx < 1 .or. indx > ens_handle%my_num_copies) then
   write(errstring, *) 'indx: ', indx, ' cannot exceed ', ens_handle%my_num_copies
   call error_handler(E_ERR,'get_ensemble_time', errstring, source, revision, revdate)
endif

mtime = ens_handle%time(indx)

end subroutine get_ensemble_time

!-----------------------------------------------------------------

subroutine end_ensemble_manager(ens_handle)

! Frees up allocated storage for an ensemble handle

type(ensemble_type), intent(inout) :: ens_handle

! Free up the allocated storage
deallocate(ens_handle%my_copies, ens_handle%time, ens_handle%my_vars, &
           ens_handle%vars,    ens_handle%copies)

end subroutine end_ensemble_manager

!-----------------------------------------------------------------

subroutine duplicate_ens(ens1, ens2, duplicate_time)

type(ensemble_type), intent(in)             :: ens1
type(ensemble_type), intent(inout)          :: ens2
logical,             intent(in)             :: duplicate_time

! Name was changed from copy_ens to avoid confusion in naming.
! Should only be used if the ens%vars storage is current (each pe has subset of
! copies of all vars).
! If duplicate_time is true, also copies the time information from ens1 to ens2. 
! If duplicate_time is false, the times in ens2 are left unchanged.

! Check to make sure that the ensembles are compatible
if(ens1%num_copies /= ens2%num_copies) then
   write(errstring, *) 'num_copies ', ens1%num_copies, ' and ', ens2%num_copies, &
      'must be equal'
   call error_handler(E_ERR,'duplicate_ens', errstring, source, revision, revdate)
endif
if(ens1%num_vars /= ens2%num_vars) then
   write(errstring, *) 'num_vars ', ens1%num_vars, ' and ', ens2%num_vars, &
      'must be equal'
   call error_handler(E_ERR,'duplicate_ens', errstring, source, revision, revdate)
endif
if(ens1%distribution_type /= ens2%distribution_type) then
   write(errstring, *) 'distribution_type ', ens1%distribution_type, ' and ', &
      ens2%distribution_type, 'must be equal'
   call error_handler(E_ERR,'duplicate_ens', errstring, source, revision, revdate)
endif

! Duplicate each copy that is stored locally on this process.
ens2%vars = ens1%vars

! Duplicate time if requested
if(duplicate_time) ens2%time = ens1%time

end subroutine duplicate_ens

!-----------------------------------------------------------------

function get_my_num_copies(ens_handle)

! Returns the number of copies stored on my processor when var complete.
! Same as num_copies if single process is in use.

integer                          :: get_my_num_copies
type (ensemble_type), intent(in) :: ens_handle

get_my_num_copies = ens_handle%my_num_copies

end function get_my_num_copies

!-----------------------------------------------------------------

subroutine get_my_copies(ens_handle, copies)

! Returns a list of all copies stored on this processor when var complete.
! Requires copies to be dimensioned with at least my_num_copies.
! Returns 1...num_copies for single process.

type (ensemble_type), intent(in)  :: ens_handle
integer,              intent(out) :: copies(:)

if(size(copies) < ens_handle%my_num_copies) then
   write(errstring, *) 'Array copies only has size ', size(copies), &
      ' but must be at least ', ens_handle%my_num_copies
   call error_handler(E_ERR, 'get_my_copies', errstring, source, revision, revdate)
endif

copies(1:ens_handle%my_num_copies) = ens_handle%my_copies(1:ens_handle%my_num_copies)

end subroutine get_my_copies

!-----------------------------------------------------------------

function get_my_num_vars(ens_handle)

! Returns the number of vars stored on my processor when copy complete.
! Same as num_vars if single process is in use.

integer                          :: get_my_num_vars
type (ensemble_type), intent(in) :: ens_handle

get_my_num_vars = ens_handle%my_num_vars

end function get_my_num_vars

!-----------------------------------------------------------------

subroutine get_my_vars(ens_handle, vars)

! Returns a list of all vars stored on this processor when copy complete.
! Requires vars to be dimensioned with at least my_num_vars.
! Returns 1...num_vars for single process.

type (ensemble_type), intent(in)  :: ens_handle
integer,              intent(out) :: vars(:)

if(size(vars) < ens_handle%my_num_vars) then
   write(errstring, *) 'Array vars only has size ', size(vars), &
      ' but must be at least ', ens_handle%my_num_vars
   call error_handler(E_ERR,'get_my_vars', errstring, source, revision, revdate)
endif

vars(1:ens_handle%my_num_vars) = ens_handle%my_vars(1:ens_handle%my_num_vars)

end subroutine get_my_vars

!-----------------------------------------------------------------

subroutine set_up_ens_distribution(ens_handle)

! Figures out how to lay-out the copy complete and vars complete
! distributions. The distribution_type identifies 
! different options. Only distribution_type 1 is implemented.
! This puts every nth var or copy on a given processor where n is the
! total number of processes.

type (ensemble_type),  intent(inout)  :: ens_handle

integer :: num_per_pe_below, num_left_over, i

! Option 1: Maximum separation for both vars and copies
! Compute the total number of copies I'll get for var complete
num_per_pe_below = ens_handle%num_copies / num_pes
num_left_over = ens_handle%num_copies - num_per_pe_below * num_pes
if(num_left_over >= (my_pe + 1)) then
   ens_handle%my_num_copies = num_per_pe_below + 1
else
   ens_handle%my_num_copies = num_per_pe_below
endif

! Do the same thing for copy complete: figure out which vars I get
num_per_pe_below = ens_handle%num_vars / num_pes
num_left_over = ens_handle%num_vars - num_per_pe_below * num_pes
if(num_left_over >= (my_pe + 1)) then
   ens_handle%my_num_vars = num_per_pe_below + 1
else
   ens_handle%my_num_vars = num_per_pe_below
endif

!Allocate the storage for copies and vars all at once
allocate(ens_handle%my_copies(ens_handle%my_num_copies),              &
         ens_handle%time     (ens_handle%my_num_copies),              &
         ens_handle%my_vars  (ens_handle%my_num_vars),                &
         ens_handle%vars     (ens_handle%num_vars,   ens_handle%my_num_copies),  &
         ens_handle%copies   (ens_handle%num_copies, ens_handle%my_num_vars))

! Set everything to missing value
ens_handle%vars = MISSING_R8
ens_handle%copies = MISSING_R8

! Fill out the number of my members
call get_copy_list(ens_handle%num_copies, my_pe, ens_handle%my_copies, i)

! Initialize times to missing
do i = 1, ens_handle%my_num_copies
   ens_handle%time(i) = set_time(0, 0)
end do

! Fill out the number of my vars
call get_var_list(ens_handle%num_vars, my_pe, ens_handle%my_vars, i)

end subroutine set_up_ens_distribution

!-----------------------------------------------------------------

subroutine get_copy_owner_index(copy_number, owner, owners_index)
!!!subroutine get_copy_owner_index(copy_number, owner, owners_index, distribution_type)

! Given the copy number, returns which PE stores it when copy complete
! and its index in that pes local storage. Depends on distribution_type
! with only option 1 currently implemented.

integer, intent(in)  :: copy_number
integer, intent(out) :: owner, owners_index
!!!integer, intent(in) :: distribution_type

integer :: div

div = (copy_number - 1) / num_pes
owner = copy_number - div * num_pes - 1
owners_index = div + 1

end subroutine get_copy_owner_index

!-----------------------------------------------------------------

subroutine get_var_owner_index(var_number, owner, owners_index)
!!!subroutine get_var_owner_index(var_number, owner, owners_index, distribution_type)

! Given the var number, returns which PE stores it when var complete
! and its index in that pes local storage. Depends on distribution_type
! with only option 1 currently implemented.

integer, intent(in)  :: var_number
integer, intent(out) :: owner, owners_index
!!!integer, intent(in) :: distribution_type

integer :: div

div = (var_number - 1) / num_pes
owner = var_number - div * num_pes - 1
owners_index = div + 1

end subroutine get_var_owner_index

!-----------------------------------------------------------------

function get_max_num_vars(num_vars)
!!!function get_max_num_vars(num_vars, distribution_type)

! Returns the largest number of vars that are on any pe when copy complete.
! Depends on distribution_type with only option 1 currently implemented.
! Used to get size for creating storage to receive a list of the vars on a pe.

integer             :: get_max_num_vars
integer, intent(in) :: num_vars
!!!integer, intent(in) :: distribution_type

get_max_num_vars = num_vars / num_pes + 1

end function get_max_num_vars

!-----------------------------------------------------------------

function get_max_num_copies(num_copies)
!!!function get_max_num_copies(num_copies, distribution_type)

! Returns the largest number of copies that are on any pe when var complete.
! Depends on distribution_type with only option 1 currently implemented.
! Used to get size for creating storage to receive a list of the copies on a pe.

integer             :: get_max_num_copies
integer, intent(in) :: num_copies
!!!integer, intent(in) :: distribution_type

get_max_num_copies = num_copies / num_pes + 1

end function get_max_num_copies

!-----------------------------------------------------------------

subroutine get_var_list(num_vars, pe, var_list, pes_num_vars)
!!!subroutine get_var_list(num_vars, pe, var_list, pes_num_vars, distribution_type)

! Returns a list of the vars stored by process pe when copy complete
! and the number of these vars.
! var_list must be dimensioned large enough to hold all vars.
! Depends on distribution_type with only option 1 currently implemented.

integer,   intent(in)     :: num_vars, pe
integer,   intent(out)    :: var_list(:), pes_num_vars
!!!integer, intent(in) :: distribution_type

integer :: num_per_pe_below, num_left_over, i

! Figure out number of vars stored by pe
num_per_pe_below = num_vars / num_pes
num_left_over = num_vars - num_per_pe_below * num_pes
if(num_left_over >= (pe + 1)) then
   pes_num_vars = num_per_pe_below + 1
else
   pes_num_vars = num_per_pe_below
endif

! Fill out the pe's vars
do i = 1, pes_num_vars
   var_list(i) = (pe + 1) + (i - 1) * num_pes
end do

end subroutine get_var_list

!-----------------------------------------------------------------

subroutine get_copy_list(num_copies, pe, copy_list, pes_num_copies)
!!!subroutine get_copy_list(num_copies, pe, copy_list, pes_num_copies, distribution_type)

! Returns a list of the copies stored by process pe when var complete.
! copy_list must be dimensioned large enough to hold all copies.
! Depends on distribution_type with only option 1 currently implemented.

integer,   intent(in)    :: num_copies, pe
integer,   intent(out)   :: copy_list(:), pes_num_copies
!!!integer, intent(in) :: distribution_type

integer :: num_per_pe_below, num_left_over, i

! Figure out which copies stored by pe
num_per_pe_below = num_copies / num_pes
num_left_over = num_copies - num_per_pe_below * num_pes
if(num_left_over >= (pe + 1)) then
   pes_num_copies = num_per_pe_below + 1
else
   pes_num_copies = num_per_pe_below
endif

! Fill out the pe's copies
do i = 1, pes_num_copies
   copy_list(i) = (pe + 1) + (i - 1) * num_pes
end do

end subroutine get_copy_list

!-----------------------------------------------------------------

subroutine all_vars_to_all_copies(ens_handle)

! Converts from having subset of copies of all variables to having
! all copies of a subset of variables on a given PE.

type (ensemble_type), intent(inout) :: ens_handle

integer,  allocatable :: var_list(:), copy_list(:)
real(r8), allocatable :: transfer_temp(:)
integer               :: num_copies, num_vars, my_num_vars, my_num_copies
integer               :: max_num_vars, max_num_copies, num_copies_to_receive
integer               :: sending_pe, recv_pe, k, sv, copy, num_vars_to_send
integer               :: global_ens_index

! Accelerated version for single process
if(num_pes == 1) then
   ens_handle%copies = transpose(ens_handle%vars)
   return
end if

! Short var definitions
num_copies    = ens_handle%num_copies
num_vars      = ens_handle%num_vars
my_num_vars   = ens_handle%my_num_vars
my_num_copies = ens_handle%my_num_copies

! What is maximum number of vars stored on a copy complete pe?
max_num_vars = get_max_num_vars(num_vars)

! What is maximum number of copies stored on a var complete pe?
max_num_copies = get_max_num_copies(num_copies)

! Allocate storage to hold lists of vars and copies
allocate(var_list(max_num_vars), transfer_temp(max_num_vars), &
   copy_list(max_num_copies))

! Loop to give each pe a turn to send its vars
SENDING_PE_LOOP: do sending_pe = 0, num_pes - 1
   ! If I'm the sending pe, do this block
   if(my_pe == sending_pe) then
      ! Figure out what piece to send to each other PE and send it
      SEND_TO_EACH: do recv_pe = 0, num_pes - 1
         call get_var_list(num_vars, recv_pe, var_list, num_vars_to_send)

         if (num_vars_to_send > 0) then
            ! Loop to send these vars for each copy stored on my_pe
            ALL_MY_COPIES: do k = 1, my_num_copies
   
               ! Fill up the transfer array
               do sv = 1, num_vars_to_send
                  transfer_temp(sv) = ens_handle%vars(var_list(sv), k)
               end do
   
               ! If sending_pe is receiving_pe, just copy
               if(sending_pe == recv_pe) then
                  global_ens_index = ens_handle%my_copies(k)
                  ens_handle%copies(global_ens_index, :) = transfer_temp(1:num_vars_to_send)
               else
                  ! Otherwise, ship this off
                  call send_to(recv_pe, transfer_temp(1:num_vars_to_send)) 
               endif
            end do ALL_MY_COPIES
         endif
      end do SEND_TO_EACH
   else
      ! I'm not the sending PE, figure out what copies of my vars I'll receive from sending_pe
      call get_copy_list(num_copies, sending_pe, copy_list, num_copies_to_receive)
       
      do copy = 1, num_copies_to_receive
         if (my_num_vars > 0) then
            ! Have to  use temp because %copies section is not contiguous storage
            call receive_from(sending_pe, transfer_temp(1:my_num_vars))
            ! Figure out which global ensemble member this is
            global_ens_index = copy_list(copy)
            ! Store this chunk in my local storage
            ens_handle%copies(global_ens_index, :) = transfer_temp(1:my_num_vars)
         endif
      end do
      
   endif
end do SENDING_PE_LOOP

! Free up the temporary storage
deallocate(var_list, transfer_temp, copy_list)

end subroutine all_vars_to_all_copies

!-----------------------------------------------------------------

subroutine all_copies_to_all_vars(ens_handle)

! Converts from having subset of copies of all variables to having
! all copies of a subset of variables on a given PE.

type (ensemble_type), intent(inout) :: ens_handle

integer,  allocatable :: var_list(:), copy_list(:)
real(r8), allocatable :: transfer_temp(:)
integer               :: num_copies, num_vars, my_num_vars, my_num_copies
integer               :: max_num_vars, max_num_copies, num_vars_to_receive
integer               :: sending_pe, recv_pe, k, sv, copy, num_copies_to_send
integer               :: global_ens_index

! Accelerated version for single process
if(num_pes == 1) then
   ens_handle%vars = transpose(ens_handle%copies)
   return
end if

! Short var definitions
num_copies    = ens_handle%num_copies
num_vars      = ens_handle%num_vars
my_num_vars   = ens_handle%my_num_vars
my_num_copies = ens_handle%my_num_copies

! What is maximum number of vars stored on a copy complete pe?
max_num_vars = get_max_num_vars(num_vars)

! What is maximum number of copies stored on a var complete pe?
max_num_copies = get_max_num_copies(num_copies)

allocate(var_list(max_num_vars), transfer_temp(max_num_vars), &
   copy_list(max_num_copies))

! Loop to give each pe a turn to receive its vars
RECEIVING_PE_LOOP: do recv_pe = 0, num_pes - 1
   ! If I'm the receiving pe, do this block
   if(my_pe == recv_pe) then
      ! Figure out what piece to receive from each other PE and receive it
      RECEIVE_FROM_EACH: do sending_pe = 0, num_pes - 1
         call get_var_list(num_vars, sending_pe, var_list, num_vars_to_receive)

         ! Loop to receive these vars for each copy stored on my_pe
         ALL_MY_COPIES: do k = 1, my_num_copies

            ! If sending_pe is receiving_pe, just copy
            if(sending_pe == recv_pe) then
               global_ens_index = ens_handle%my_copies(k)
               do sv = 1, num_vars_to_receive
                  ens_handle%vars(var_list(sv), k) = ens_handle%copies(global_ens_index, sv)
               end do
            else
               if (num_vars_to_receive > 0) then
                  ! Otherwise, receive this part from the sending pe
                  call receive_from(sending_pe, transfer_temp(1:num_vars_to_receive)) 
   
                  ! Copy the transfer array to my local storage
                  do sv = 1, num_vars_to_receive
                     ens_handle%vars(var_list(sv), k) = transfer_temp(sv)
                  end do
               endif
            endif
         end do ALL_MY_COPIES
      end do RECEIVE_FROM_EACH
   else
      ! I'm the sending PE, figure out what copies of my vars I'll send.
      call get_copy_list(num_copies, recv_pe, copy_list, num_copies_to_send)
       
      do copy = 1, num_copies_to_send
         if (my_num_vars > 0) then
            transfer_temp(1:my_num_vars) = ens_handle%copies(copy_list(copy), :)
            ! Have to  use temp because %copies section is not contiguous storage
            call send_to(recv_pe, transfer_temp(1:my_num_vars))
         endif
      end do
      
   endif
end do RECEIVING_PE_LOOP

! Free up the temporary storage
deallocate(var_list, transfer_temp, copy_list)

end subroutine all_copies_to_all_vars

!-----------------------------------------------------------------

subroutine compute_copy_mean(ens_handle, start_copy, end_copy, mean_copy)

! Assumes that ens_handle%copies is current; each pe has all copies of subset of vars
! Computes the mean of ensemble copies start_copy:end_copy and stores
! the result in mean_copy

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: start_copy, end_copy, mean_copy

integer :: num_copies, i

! Should check to make sure that start, end and mean are all legal

num_copies = end_copy - start_copy + 1

do i = 1, ens_handle%my_num_vars
   ens_handle%copies(mean_copy, i) = sum(ens_handle%copies(start_copy:end_copy, i)) / num_copies
end do

end subroutine compute_copy_mean

!--------------------------------------------------------------------------------

subroutine compute_copy_mean_sd(ens_handle, start_copy, end_copy, mean_copy, sd_copy)

! Assumes that ens_handle%copies is current; each pe has all copies of subset of vars
! Computes the mean and sd of ensemble copies start_copy:end_copy and stores
! mean in copy mean_copy and sd in copy sd_copy.

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: start_copy, end_copy, mean_copy, sd_copy

integer :: num_copies, i

! Should check to make sure that start, end, mean and sd are all legal copies

num_copies = end_copy - start_copy + 1

do i = 1, ens_handle%my_num_vars
   ens_handle%copies(mean_copy, i) = sum(ens_handle%copies(start_copy:end_copy, i)) / num_copies
   ens_handle%copies(sd_copy, i)   = sqrt((sum((ens_handle%copies(start_copy:end_copy, i) - &
      ens_handle%copies(mean_copy, i))**2) / (num_copies - 1)))
end do

end subroutine compute_copy_mean_sd
!--------------------------------------------------------------------------------

subroutine compute_copy_mean_var(ens_handle, start_copy, end_copy, mean_copy, var_copy)

! Assumes ensemble complete. 
! Computes the mean and variance of ensemble copies start_copy:end_copy and stores
! mean in copy mean_copy and variance in copy var_copy.

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: start_copy, end_copy, mean_copy, var_copy

integer :: num_copies, i

! Should check to make sure that start, end, mean and var are all legal copies

num_copies = end_copy - start_copy + 1

do i = 1, ens_handle%my_num_vars
   ens_handle%copies(mean_copy, i) = sum(ens_handle%copies(start_copy:end_copy, i)) / num_copies
   ens_handle%copies(var_copy, i)   = (sum((ens_handle%copies(start_copy:end_copy, i) - &
      ens_handle%copies(mean_copy, i))**2) / (num_copies - 1))
end do

end subroutine compute_copy_mean_var

!--------------------------------------------------------------------------------

end module ensemble_manager_mod
