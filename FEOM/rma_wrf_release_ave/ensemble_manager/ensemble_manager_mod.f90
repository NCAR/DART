! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module ensemble_manager_mod

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
                              check_namelist_read, timestamp, set_output
use assim_model_mod,   only : aread_state_restart, awrite_state_restart, &
                              open_restart_read, open_restart_write,     &
                              close_restart, pert_model_state
use time_manager_mod,  only : time_type, set_time
use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian
use mpi_utilities_mod, only : task_count, my_task_id, send_to, receive_from, &
                              task_sync, broadcast_send, broadcast_recv
use sort_mod,          only : index_sort

use data_structure_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index

use copies_on_off_mod,  only : query_read_copy, query_write_copy

use io_filenames_mod,   only : restart_files_out

implicit none
private

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!PAR other storage option control can be implemented here. In particular, want to find
!PAR some way, either allocating or multiple addressing, to use same chunk of storage
!PAR for both copy and var complete representations.

public :: init_ensemble_manager,      end_ensemble_manager,     get_ensemble_time,      &
          ensemble_type,              duplicate_ens,            get_var_owner_index,    &
          get_my_num_copies,          get_my_copies,            get_my_num_vars,        &
          get_my_vars,                compute_copy_mean,        compute_copy_mean_sd,   &
          get_copy,                   put_copy,                 all_vars_to_all_copies, &
          all_copies_to_all_vars,     read_ensemble_restart,    write_ensemble_restart, &
          compute_copy_mean_var,      get_copy_owner_index,     set_ensemble_time,      &
          broadcast_copy,             prepare_to_write_to_vars, prepare_to_write_to_copies, &
          prepare_to_read_from_vars,  prepare_to_read_from_copies, prepare_to_update_vars,  &
          prepare_to_update_copies,   print_ens_handle,                                 &
          map_task_to_pe,             map_pe_to_task, single_restart_file_in, single_restart_file_out

! track if copies modified last, vars modified last, both are in sync
! (and therefore both valid to be used r/o), or unknown.
integer, parameter :: VALID_UNKNOWN = -1
integer, parameter :: VALID_BOTH    =  0      ! vars & copies have same data
integer, parameter :: VALID_VARS    =  1      ! vars(:,:) modified last
integer, parameter :: VALID_COPIES  =  2      ! copies(:,:) modified last

! unique counter per ensemble handle
integer              :: global_counter = 1

! Logical flag for initialization of module
logical              :: module_initialized = .false.

! Module storage for writing error messages
character(len = 255) :: msgstring

! Module storage for pe information for this process avoids recomputation
integer              :: num_pes

! Control order of communication loops in the transpose routines
logical  :: use_copy2var_send_loop = .true.
logical  :: use_var2copy_rec_loop  = .true.

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
! Complain if unneeded transposes are done
logical  :: flag_unneeded_transposes = .false.
! Communication configuration:
!  1 = usual default, 2 - 4 are valid and depend on the machine, ensemble count, and task count
integer  :: communication_configuration = 1
! task layout options:
integer  :: layout = 1 ! default to my_pe = my_task_id(). Layout2 assumes that the user knows the correct tasks_per_node
integer  :: tasks_per_node = 1 ! default to 1 if the user does not specify a number of tasks per node.
logical  :: debug = .false.

namelist / ensemble_manager_nml / single_restart_file_in,  &
                                  single_restart_file_out, &
                                  perturbation_amplitude,  &
                                  communication_configuration, &
                                  layout, tasks_per_node,  &
                                  debug, flag_unneeded_transposes
                                  
!-----------------------------------------------------------------

contains

!-----------------------------------------------------------------

subroutine init_ensemble_manager(ens_handle, num_copies, &
   num_vars, distribution_type_in, layout_type)

type(ensemble_type), intent(out)            :: ens_handle
integer,             intent(in)             :: num_copies, num_vars
integer,             intent(in), optional   :: distribution_type_in
integer,             intent(in), optional   :: layout_type

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
endif

! Optional layout_type argument to assign how my_pe is related to my_task_id
! layout_type can be set individually for each ensemble handle. It is not advisable to do this
! because get_obs_ens assumes that the layout is the same for each ensemble handle.
if(.not. present(layout_type) ) then
   ens_handle%layout_type = layout ! namelist option
else
   ens_handle%layout_type = layout_type
endif

! Check for error: only layout_types 1 and 2 are implemented
if (ens_handle%layout_type /= 1 .and. ens_handle%layout_type /=2 ) then
   call error_handler(E_ERR, 'init_ensemble_manager', 'only layout values 1 (standard), 2 (round-robin) allowed ', &
                      source, revision, revdate)
endif

allocate(ens_handle%task_to_pe_list(num_pes))
allocate(ens_handle%pe_to_task_list(num_pes))

call assign_tasks_to_pes(ens_handle, num_copies, ens_handle%layout_type)
ens_handle%my_pe = map_task_to_pe(ens_handle, my_task_id())

! Set the global storage bounds for the number of copies and variables
ens_handle%num_copies = num_copies
ens_handle%num_vars = num_vars

! For debugging, error checking
ens_handle%id_num = global_counter
global_counter = global_counter + 1

! This controls which way the two transpose routines order their
! loops:  single sender with multiple receivers, or single receiver
! with multiple senders.  For small problems it doesn't matter;
! for very large ensemble sizes and very large MPI task counts you
! will want to profile each combination and pick the fastest.
if (communication_configuration == 1) then
   use_copy2var_send_loop = .true.
   use_var2copy_rec_loop  = .true.
else if (communication_configuration == 2) then
   use_copy2var_send_loop = .false.
   use_var2copy_rec_loop  = .true.
else if (communication_configuration == 3) then
   use_copy2var_send_loop = .true.
   use_var2copy_rec_loop  = .false.
else if (communication_configuration == 4) then
   use_copy2var_send_loop = .false.
   use_var2copy_rec_loop  = .false.
else
   write(msgstring, *) 'communication_configuration is', communication_configuration
   call error_handler(E_ERR, 'init_ensemble_manager', &
      'communication_configuration must be between 1 and 4', &
       source, revision, revdate, text2=msgstring)
endif

! initially no data
ens_handle%valid = VALID_BOTH

if(debug .and. my_task_id()==0) then
   print*, 'pe_to_task_list', ens_handle%pe_to_task_list
   print*, 'task_to_pe_list', ens_handle%task_to_pe_list
endif

! Figure out how the ensemble copies are partitioned
call set_up_ens_distribution(ens_handle)

if (debug) call print_ens_handle(ens_handle)

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

! Avoid num_vars size storage on stack; make this allocatable from heap
real(r8), allocatable               :: ens(:) 
integer                             :: iunit, i, j
character(len = LEN(file_name) + 5) :: this_file_name
character(len = 4)                  :: extension
type(time_type)                     :: ens_time
integer                             :: global_copy_index
logical                             :: interf_provided
logical                             :: single_file_override
type(random_seq_type)               :: random_seq
integer                             :: my_num_vars

! timing variables
double precision :: time_at_start

! Does not make sense to have start_from_restart and single_restart_file_in BOTH false
if(.not. start_from_restart .and. .not. single_restart_file_in) then
   write(msgstring, *) 'start_from_restart in filter_nml and single_restart_file_in in &
        &ensemble_manager_nml cannot both be false'
   call error_handler(E_ERR,'read_ensemble_restart', msgstring, source, revision, revdate)
endif

! the code reads into the vars array
ens_handle%valid = VALID_VARS

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
   ! Single restart file is read only by task 0 and then distributed
   if(my_task_id() == 0) iunit = open_restart_read(file_name)
   allocate(ens(ens_handle%num_vars))   ! used to be on stack.

   ! Loop through the total number of copies
   do i = start_copy, end_copy
       ! Only task 0 does reading. Everybody can do their own perturbing
       if(my_task_id() == 0) then
       ! Read restarts in sequence; only read once for not start_from_restart
          if(start_from_restart .or. i == start_copy) &
               call aread_state_restart(ens_time, ens, iunit)
         ! Override this time if requested by namelist
         if(present(init_time)) ens_time = init_time
      endif

      ! Store this copy in the appropriate place on the appropriate process
      ! map from my_pe to physical task number all done in send and recieves only
      call put_copy(map_task_to_pe(ens_handle,0), ens_handle, i, ens, ens_time)
   end do
   
   deallocate(ens)
   ! Task 0 must close the file it's been reading
   if(my_task_id() == 0) call close_restart(iunit)

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
               if (ens_handle%vars(j,i) /= MISSING_R8) &
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
real(r8), allocatable               :: ens(:)
type(time_type)                     :: ens_time
integer                             :: iunit, i, global_index
integer                             :: owner, owners_index
character(len = LEN(file_name) + 10) :: this_file_name
character(len = 4)                  :: extension
logical                             :: single_file_forced

! timing variables
double precision :: start_at_time

if (present(force_single_file) ) then
   single_file_forced = force_single_file
else
   single_file_forced = .FALSE.
endif

! Error checking
!if (ens_handle%valid /= VALID_VARS .and. ens_handle%valid /= VALID_BOTH) then
!   call error_handler(E_ERR, 'write_ensemble_restart', &
!        'last access not var-complete', source, revision, revdate)
!endif

! For single file, need to send restarts to pe0 and it writes them out.
!-------------- Block for single_restart file -------------
! Need to force single restart file for inflation files
if(single_restart_file_out .or. single_file_forced) then

   ! Single restart file is written only by task 0
   if(my_task_id() == 0) then

      iunit = open_restart_write(file_name)

      ! Loop to write each ensemble member
      allocate(ens(ens_handle%num_vars))   ! used to be on stack.
      do i = start_copy, end_copy
         ! Figure out where this ensemble member is being stored
         call get_copy_owner_index(i, owner, owners_index)
         ! If it's on task 0, just write it
         if(map_pe_to_task(ens_handle, owner) == 0) then
            call awrite_state_restart(ens_handle%time(owners_index), &
               ens_handle%vars(:, owners_index), iunit)
         else
            ! Get the ensemble from the owner and write it out
            ! This communication assumes index numbers are monotonically increasing
            ! and that communications is blocking so that there are not multiple
            ! outstanding messages from the same processors (also see send_to below).
            call receive_from(map_pe_to_task(ens_handle, owner), ens, ens_time)
            call awrite_state_restart(ens_time, ens, iunit)
         endif
      end do
      deallocate(ens)
      call close_restart(iunit)
   else ! I must send my copies to task 0 for writing to file
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

! Error checking
if (ens_handle%valid /= VALID_VARS .and. ens_handle%valid /= VALID_BOTH) then
   call error_handler(E_ERR, 'get_copy', &
        'last access not var-complete', source, revision, revdate)
endif

! Verify that requested copy exists
if(copy < 1 .or. copy > ens_handle%num_copies) then
   write(msgstring, *) 'Requested copy: ', copy, ' is > maximum copy: ', ens_handle%num_copies 
   call error_handler(E_ERR,'get_copy', msgstring, source, revision, revdate)
endif

! Make sure that vars has enough space to handle the answer
if(ens_handle%my_pe == receiving_pe) then !HK I think only the reciever needs the space
   if(size(vars) < ens_handle%num_vars) then
      write(msgstring, *) 'Size of vars: ', size(vars), ' Must be at least ', ens_handle%num_vars
      call error_handler(E_ERR,'get_copy', msgstring, source, revision, revdate)
   endif
endif

! Figure out which PE stores this copy and what its local storage index is
call get_copy_owner_index(copy, owner, owners_index)

!----------- Block of code that must be done by receiving pe -----------------------------
if(ens_handle%my_pe == receiving_pe) then
   ! If PE that stores is the same, just copy and return
   if(ens_handle%my_pe == owner) then
      vars = ens_handle%vars(:, owners_index)
      if(present(mtime)) mtime = ens_handle%time(owners_index)
      ! If I'm the receiving PE and also the owner, I'm all finished; return
      return
   endif
 
   ! Otherwise, must wait to receive vars and time from storing pe
      call receive_from(map_pe_to_task(ens_handle, owner), vars, mtime)
endif

!----- Block of code that must be done by PE that stores the copy IF it is NOT receiver -----
if(ens_handle%my_pe == owner) then
   ! Send copy to receiving pe

   if(present(mtime)) then
      call send_to(map_pe_to_task(ens_handle, receiving_pe), ens_handle%vars(:, owners_index), ens_handle%time(owners_index))
   else
      call send_to(map_pe_to_task(ens_handle, receiving_pe), ens_handle%vars(:, owners_index))
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

! Error checking
if (ens_handle%valid /= VALID_VARS .and. ens_handle%valid /= VALID_BOTH) then
   call error_handler(E_ERR, 'put_copy', &
        'last access not var-complete', source, revision, revdate)
endif

if(copy < 1 .or. copy > ens_handle%num_copies) then
   write(msgstring, *) 'Requested copy: ', copy, ' is > maximum copy: ', ens_handle%num_copies 
   call error_handler(E_ERR,'put_copy', msgstring, source, revision, revdate)
endif

! Make sure that num_vars has enough space to handle the answer
if(ens_handle%num_vars < size(vars)) then
   write(msgstring, *) 'Size of vars: ', size(vars), ' Cannot be more than ', ens_handle%num_vars
   call error_handler(E_ERR,'put_copy', msgstring, source, revision, revdate)
endif

! What PE stores this copy and what is its local storage index
call get_copy_owner_index(copy, owner, owners_index)

! Block of code that must be done by PE that is to send the copy
if(ens_handle%my_pe == sending_pe) then
   ! If PE that stores is the same, just copy and return
   if(ens_handle%my_pe == owner) then
      ens_handle%vars(:, owners_index) = vars
      if(present(mtime)) ens_handle%time(owners_index) = mtime
      ! If I'm the sending PE and also the owner, I'm all finished; return
      return
   endif
 
   ! Otherwise, must send vars and possibly time to storing pe
   call send_to(map_pe_to_task(ens_handle, owner), vars, mtime)

endif

! Block of code that must be done by PE that stores the copy IF it is NOT sender
if(ens_handle%my_pe == owner) then
   ! Need to receive copy from sending_pe
   if(present(mtime)) then
      call receive_from(map_pe_to_task(ens_handle, sending_pe), ens_handle%vars(:, owners_index), ens_handle%time(owners_index))
   else
      call receive_from(map_pe_to_task(ens_handle, sending_pe), ens_handle%vars(:, owners_index))
   endif
endif

ens_handle%valid = VALID_VARS

end subroutine put_copy

!-----------------------------------------------------------------

subroutine broadcast_copy(ens_handle, copy, arraydata)

! find which PE has the global copy number and have it broadcast 
! that copy to all the other PEs.  arraydata is an output on
! all PEs, even on the PE which is the owner it is separate
! storage from the vars array in the ensemble handle.

type(ensemble_type), intent(in)           :: ens_handle
integer,             intent(in)           :: copy
real(r8),            intent(out)          :: arraydata(:)

integer :: owner, owners_index

! Error checking
if (ens_handle%valid /= VALID_VARS .and. ens_handle%valid /= VALID_BOTH) then
   call error_handler(E_ERR, 'broadcast_copy', &
        'last access not var-complete', source, revision, revdate)
endif

if(copy < 1 .or. copy > ens_handle%num_copies) then
   write(msgstring, *) 'Requested copy: ', copy, ' is > maximum copy: ', ens_handle%num_copies 
   call error_handler(E_ERR,'broadcast_copy', msgstring, source, revision, revdate)
endif

! Make sure that arraydata has enough space to handle the answer
if(size(arraydata) < ens_handle%num_vars) then
   write(msgstring, *) 'Size of arraydata: ', size(arraydata), ' must be at least ', ens_handle%num_vars
   call error_handler(E_ERR,'broadcast_copy', msgstring, source, revision, revdate)
endif

! What PE stores this copy and what is its local storage index
call get_copy_owner_index(copy, owner, owners_index)

! First block of code that must be done by PE that is to send the copy
if(ens_handle%my_pe == owner) then
   arraydata = ens_handle%vars(:, owners_index) 
   call broadcast_send(map_pe_to_task(ens_handle, owner), arraydata)
else 
   call broadcast_recv(map_pe_to_task(ens_handle, owner), arraydata)
endif

end subroutine broadcast_copy

!-----------------------------------------------------------------

subroutine prepare_to_write_to_vars(ens_handle)

! Warn ens manager that we're going to directly update the %vars array

type(ensemble_type), intent(inout) :: ens_handle

!ens_handle%valid = VALID_VARS

end subroutine prepare_to_write_to_vars

!-----------------------------------------------------------------

subroutine prepare_to_write_to_copies(ens_handle)

! Warn ens manager that we're going to directly update the %copies array

type(ensemble_type), intent(inout) :: ens_handle

!ens_handle%valid = VALID_COPIES

end subroutine prepare_to_write_to_copies

!-----------------------------------------------------------------

subroutine prepare_to_read_from_vars(ens_handle)

! Check to be sure that the vars array is current

type(ensemble_type), intent(in) :: ens_handle

!if (ens_handle%valid /= VALID_VARS .and. ens_handle%valid /= VALID_BOTH) then
!   call error_handler(E_ERR, 'prepare_to_read_from_vars', &
 !       'last access not var-complete', source, revision, revdate)
!endif

end subroutine prepare_to_read_from_vars

!-----------------------------------------------------------------

subroutine prepare_to_read_from_copies(ens_handle)

! Check to be sure that the copies array is current

type(ensemble_type), intent(in) :: ens_handle

!if (ens_handle%valid /= VALID_COPIES .and. ens_handle%valid /= VALID_BOTH) then
!   call error_handler(E_ERR, 'prepare_to_read_from_copies', &
!        'last access not copy-complete', source, revision, revdate)
!endif

end subroutine prepare_to_read_from_copies

!-----------------------------------------------------------------

subroutine prepare_to_update_vars(ens_handle)

! We need read/write access, so it has to start valid for vars or both,
! and then is going to be vars only going out.

type(ensemble_type), intent(inout) :: ens_handle

!if (ens_handle%valid /= VALID_VARS .and. ens_handle%valid /= VALID_BOTH) then
!   call error_handler(E_ERR, 'prepare_to_update_vars', &
 !       'last access not var-complete', source, revision, revdate)
!endif
!ens_handle%valid = VALID_VARS

end subroutine prepare_to_update_vars

!-----------------------------------------------------------------

subroutine prepare_to_update_copies(ens_handle)

! We need read/write access, so it has to start valid for copies or both,
! and then is going to be copies only going out.

type(ensemble_type), intent(inout) :: ens_handle

!if (ens_handle%valid /= VALID_COPIES .and. ens_handle%valid /= VALID_BOTH) then
!   call error_handler(E_ERR, 'prepare_to_update_copies', &
!        'last access not copy-complete', source, revision, revdate)
!endif
!ens_handle%valid = VALID_COPIES

end subroutine prepare_to_update_copies

!-----------------------------------------------------------------

subroutine set_ensemble_time(ens_handle, indx, mtime)

! Sets the time of an ensemble member indexed by local storage on this pe.

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: indx
type(time_type),     intent(in)    :: mtime

if(indx < 1 .or. indx > ens_handle%my_num_copies) then
   write(msgstring, *) 'indx: ', indx, ' cannot exceed ', ens_handle%my_num_copies
   call error_handler(E_ERR,'get_ensemble_time', msgstring, source, revision, revdate)
endif

ens_handle%time(indx) = mtime

end subroutine set_ensemble_time

!-----------------------------------------------------------------

subroutine get_ensemble_time(ens_handle, indx, mtime)

! Returns the time of an ensemble member indexed by local storage on this pe.

type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: indx
type(time_type),     intent(out) :: mtime

if(indx < 1 .or. indx > ens_handle%my_num_copies) then
   write(msgstring, *) 'indx: ', indx, ' cannot exceed ', ens_handle%my_num_copies
   call error_handler(E_ERR,'get_ensemble_time', msgstring, source, revision, revdate)
endif

mtime = ens_handle%time(indx)

end subroutine get_ensemble_time

!-----------------------------------------------------------------

subroutine end_ensemble_manager(ens_handle)

! Frees up allocated storage for an ensemble handle

type(ensemble_type), intent(inout) :: ens_handle

! Free up the allocated storage
deallocate(ens_handle%my_copies, ens_handle%time, ens_handle%my_vars, &
           ens_handle%copies, ens_handle%task_to_pe_list, ens_handle%pe_to_task_list)

end subroutine end_ensemble_manager

!-----------------------------------------------------------------

subroutine duplicate_ens(ens1, ens2, duplicate_time)

type(ensemble_type), intent(in)             :: ens1
type(ensemble_type), intent(inout)          :: ens2
logical,             intent(in)             :: duplicate_time

! Name was changed from copy_ens to avoid confusion in naming.
! Should only be used if the ens1%vars storage is current (each pe has subset of
! copies of all vars).
! If duplicate_time is true, also copies the time information from ens1 to ens2. 
! If duplicate_time is false, the times in ens2 are left unchanged.

! Error checking
if (ens1%valid /= VALID_VARS .and. ens1%valid /= VALID_BOTH) then
   call error_handler(E_ERR, 'duplicate_ens', &
        'last access not var-complete for source ensemble', source, revision, revdate)
endif

! Check to make sure that the ensembles are compatible
if(ens1%num_copies /= ens2%num_copies) then
   write(msgstring, *) 'num_copies ', ens1%num_copies, ' and ', ens2%num_copies, &
      'must be equal'
   call error_handler(E_ERR,'duplicate_ens', msgstring, source, revision, revdate)
endif
if(ens1%num_vars /= ens2%num_vars) then
   write(msgstring, *) 'num_vars ', ens1%num_vars, ' and ', ens2%num_vars, &
      'must be equal'
   call error_handler(E_ERR,'duplicate_ens', msgstring, source, revision, revdate)
endif
if(ens1%distribution_type /= ens2%distribution_type) then
   write(msgstring, *) 'distribution_type ', ens1%distribution_type, ' and ', &
      ens2%distribution_type, 'must be equal'
   call error_handler(E_ERR,'duplicate_ens', msgstring, source, revision, revdate)
endif

! Duplicate each copy that is stored locally on this process.
ens2%vars = ens1%vars

ens2%valid = VALID_VARS

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
   write(msgstring, *) 'Array copies only has size ', size(copies), &
      ' but must be at least ', ens_handle%my_num_copies
   call error_handler(E_ERR, 'get_my_copies', msgstring, source, revision, revdate)
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
   write(msgstring, *) 'Array vars only has size ', size(vars), &
      ' but must be at least ', ens_handle%my_num_vars
   call error_handler(E_ERR,'get_my_vars', msgstring, source, revision, revdate)
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
if(num_left_over >= (ens_handle%my_pe + 1)) then
   ens_handle%my_num_copies = num_per_pe_below + 1
else
   ens_handle%my_num_copies = num_per_pe_below
endif

! Do the same thing for copy complete: figure out which vars I get
num_per_pe_below = ens_handle%num_vars / num_pes
num_left_over = ens_handle%num_vars - num_per_pe_below * num_pes
if(num_left_over >= (ens_handle%my_pe + 1)) then
   ens_handle%my_num_vars = num_per_pe_below + 1
else
   ens_handle%my_num_vars = num_per_pe_below
endif

!Allocate the storage for copies and vars all at once
allocate(ens_handle%my_copies(ens_handle%my_num_copies),              &
         ens_handle%time     (ens_handle%num_copies),                 &
         ens_handle%my_vars  (ens_handle%my_num_vars),                &
         ens_handle%copies   (ens_handle%num_copies, ens_handle%my_num_vars))

! Set everything to missing value
ens_handle%copies = MISSING_R8

! Fill out the number of my members
call get_copy_list(ens_handle%num_copies, ens_handle%my_pe, ens_handle%my_copies, i)

! Initialize times to missing
ens_handle%time(:) = set_time(0, 0) ! every task has the time for every copy

! Fill out the number of my vars
call get_var_list(ens_handle%num_vars, ens_handle%my_pe, ens_handle%my_vars, i)

end subroutine set_up_ens_distribution

!-----------------------------------------------------------------

subroutine get_copy_owner_index(copy_number, owner, owners_index)

! Given the copy number, returns which PE stores it when var complete
! and its index in that pes local storage. Depends on distribution_type
! with only option 1 currently implemented.


integer, intent(in)  :: copy_number
integer, intent(out) :: owner, owners_index

integer :: div

! Asummes distribution type 1
div = (copy_number - 1) / num_pes
owner = copy_number - div * num_pes - 1
owners_index = div + 1

end subroutine get_copy_owner_index

!!-----------------------------------------------------------------
!
!subroutine get_var_owner_index(var_number, owner, owners_index)
!
!! Given the var number, returns which PE stores it when copy complete
!! and its index in that pes local storage. Depends on distribution_type
!! with only option 1 currently implemented.
!
!! Assumes that all tasks are used in the ensemble
!
!integer, intent(in)  :: var_number
!integer, intent(out) :: owner, owners_index
!
!integer :: div
!
!! Asummes distribution type 1
!div = (var_number - 1) / num_pes
!owner = var_number - div * num_pes - 1
!owners_index = div + 1
!
!end subroutine get_var_owner_index

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

subroutine all_vars_to_all_copies(ens_handle, label)

! Converts from having subset of copies of all variables to having
! all copies of a subset of variables on a given PE.
!
! updated version of the original routine.  here, all tasks send
! while 1 receives.  original was all tasks receiving while 1 sends.
! apparently the sends can overlap and the execution time is less.
! HK Added namelist option (use_var2copy_rec_loop) to select updated
! or original version of the routine. 
!   Default: use updated version

type (ensemble_type), intent(inout)        :: ens_handle
character (len=*),    intent(in), optional :: label

integer,  allocatable :: var_list(:), copy_list(:)
real(r8), allocatable :: transfer_temp(:)
integer               :: num_copies, num_vars, my_num_vars, my_num_copies, my_pe
integer               :: max_num_vars, max_num_copies, num_copies_to_receive
integer               :: sending_pe, recv_pe, k, sv, num_vars_to_send, copy
integer               :: global_ens_index

! only output if there is a label
if (present(label)) then
   call timestamp_message('vars_to_copies start: '//label, alltasks=.true.)
endif

! Error checking, but can't return early in case only some of the
! MPI tasks need to transpose.  Only if all N tasks say this is an
! unneeded transpose can we skip it.
!if (ens_handle%valid == VALID_BOTH) then
!   if (flag_unneeded_transposes) then
!      write(msgstring, *) 'task ', my_task_id(), ' ens_handle ', ens_handle%id_num
!      call error_handler(E_MSG, 'all_vars_to_all_copies', &
!           'vars & copies both valid, transpose not needed for this task', &
!            source, revision, revdate, text2=msgstring)
!   endif
!else if (ens_handle%valid /= VALID_VARS) then
!   write(msgstring, *) 'ens_handle ', ens_handle%id_num
!   call error_handler(E_ERR, 'all_vars_to_all_copies', &
!        'last access not var-complete', source, revision, revdate, &
!         text2=msgstring)
!endif

ens_handle%valid = VALID_BOTH

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
my_pe         = ens_handle%my_pe

! What is maximum number of vars stored on a copy complete pe?
max_num_vars = get_max_num_vars(num_vars)

! What is maximum number of copies stored on a var complete pe?
max_num_copies = get_max_num_copies(num_copies)

allocate(var_list(max_num_vars), transfer_temp(max_num_vars), &
   copy_list(max_num_copies))

if ( use_var2copy_rec_loop .eqv. .true. ) then ! use updated version

  ! Loop to give each pe a turn to receive its copies
  RECEIVING_PE_LOOP: do recv_pe = 0, num_pes - 1
     ! If I'm the receiving pe, do this block
     if(my_pe == recv_pe) then

        ! Figure out what piece to receive from each other PE and receive it
        RECEIVE_FROM_EACH: do sending_pe = 0, num_pes - 1
           call get_copy_list(num_copies, sending_pe, copy_list, num_copies_to_receive)

           ! Loop to receive for each copy stored on my_pe
           ALL_MY_COPIES: do k = 1, num_copies_to_receive

              global_ens_index = copy_list(k)

              ! If sending_pe is receiving_pe, just copy
              if(sending_pe == recv_pe) then
                 do sv = 1, my_num_vars
                    ens_handle%copies(global_ens_index, sv) = ens_handle%vars(ens_handle%my_vars(sv), k)
                 end do
              else
                 if (num_copies_to_receive > 0) then
                    ! Otherwise, receive this part from the sending pe
                    call receive_from(map_pe_to_task(ens_handle, sending_pe), transfer_temp(1:my_num_vars))
   
                    ! Copy the transfer array to my local storage
                    ens_handle%copies(global_ens_index, :) = transfer_temp(1:my_num_vars)
                 endif
              endif
           end do ALL_MY_COPIES
        end do RECEIVE_FROM_EACH
     else
        ! I'm the sending PE, figure out what vars of my copies I'll send.
        call get_var_list(num_vars, recv_pe, var_list, num_vars_to_send)
       
        do k = 1, my_num_copies
           do sv = 1, num_vars_to_send
              ! Have to use temp because %var section is not contiguous storage
              transfer_temp(sv) = ens_handle%vars(var_list(sv), k)
           enddo
           call send_to(map_pe_to_task(ens_handle, recv_pe), transfer_temp(1:num_vars_to_send))
        end do
      
     endif
  end do RECEIVING_PE_LOOP

else ! use older version

  ! Loop to give each pe a turn to send its vars
  SENDING_PE_LOOP: do sending_pe = 0, num_pes - 1
    ! If I'm the sending pe, do this block
    if(my_pe == sending_pe) then
         ! Figure out what piece to send to each other PE and send it
         SEND_TO_EACH: do recv_pe = 0, num_pes - 1
           call get_var_list(num_vars, recv_pe, var_list, num_vars_to_send)

           if (num_vars_to_send > 0) then
             ! Loop to send these vars for each copy stored on my_pe
             ALL_MY_COPIES_SEND_LOOP: do k = 1, my_num_copies

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
                 call send_to(map_pe_to_task(ens_handle, recv_pe), transfer_temp(1:num_vars_to_send))
               endif
             end do ALL_MY_COPIES_SEND_LOOP
           endif
         end do SEND_TO_EACH

    else
       ! I'm not the sending PE, figure out what copies of my vars I'll receive from sending_pe
        call get_copy_list(num_copies, sending_pe, copy_list, num_copies_to_receive)

        do copy = 1, num_copies_to_receive
          if (my_num_vars > 0) then
            ! Have to  use temp because %copies section is not contiguous storage
            call receive_from(map_pe_to_task(ens_handle, sending_pe), transfer_temp(1:my_num_vars))
            ! Figure out which global ensemble member this is
            global_ens_index = copy_list(copy)
            ! Store this chunk in my local storage
            ens_handle%copies(global_ens_index, :) = transfer_temp(1:my_num_vars)
         endif
        end do

    endif

  end do SENDING_PE_LOOP

endif


! Free up the temporary storage
deallocate(var_list, transfer_temp, copy_list)

! only output if there is a label
if (present(label)) then
   call timestamp_message('vars_to_copies   end: '//label, alltasks=.true.)
endif

end subroutine all_vars_to_all_copies

!-----------------------------------------------------------------

subroutine all_copies_to_all_vars(ens_handle, label)

! Converts from having subset of copies of all variables to having
! all copies of a subset of variables on a given PE.

type (ensemble_type), intent(inout) :: ens_handle
character (len=*),    intent(in), optional :: label

integer,  allocatable :: var_list(:), copy_list(:)
real(r8), allocatable :: transfer_temp(:)
integer               :: num_copies, num_vars, my_num_vars, my_num_copies, my_pe
integer               :: max_num_vars, max_num_copies, num_vars_to_receive
integer               :: sending_pe, recv_pe, k, sv, copy, num_copies_to_send
integer               :: global_ens_index

! only output if there is a label
if (present(label)) then
   call timestamp_message('copies_to_vars start: '//label, alltasks=.true.)
endif

! Error checking, but can't return early in case only some of the
! MPI tasks need to transpose.  Only if all N tasks say this is an
! unneeded transpose can we skip it.
!if (ens_handle%valid == VALID_BOTH) then
!   if (flag_unneeded_transposes) then
!      write(msgstring, *) 'task ', my_task_id(), ' ens_handle ', ens_handle%id_num
!      call error_handler(E_MSG, 'all_copies_to_all_vars', &
!           'vars & copies both valid, transpose not needed for this task', &
!            source, revision, revdate, text2=msgstring)
!   endif
!else if (ens_handle%valid /= VALID_COPIES) then
!   write(msgstring, *) 'ens_handle ', ens_handle%id_num
!   call error_handler(E_ERR, 'all_copies_to_all_vars', &
!        'last access not copy-complete', source, revision, revdate, &
!         text2=msgstring)
!endif

ens_handle%valid = VALID_BOTH

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
my_pe         = ens_handle%my_pe

! What is maximum number of vars stored on a copy complete pe?
max_num_vars = get_max_num_vars(num_vars)

! What is maximum number of copies stored on a var complete pe?
max_num_copies = get_max_num_copies(num_copies)

allocate(var_list(max_num_vars), transfer_temp(max_num_vars), &
   copy_list(max_num_copies))


if (use_copy2var_send_loop .eqv. .true. ) then
! Switched loop index from receiving_pe to sending_pe
! Aim: to make the commication scale better on Yellowstone, as num_pes >> ens_size
! For small numbers of tasks (32 or less) the receiving_pe loop may be faster.
! Namelist option use_copy2var_send_loop can be used to select which
! communication pattern to use
!    Default: use sending_pe loop (use_copy2var_send_loop = .true.)

SENDING_PE_LOOP: do sending_pe = 0, num_pes - 1
 
   if (my_pe /= sending_pe ) then

      ! figure out what piece to recieve from each other PE and recieve it
      call get_var_list(num_vars, sending_pe, var_list, num_vars_to_receive)

      if( num_vars_to_receive > 0 ) then
         ! Loop to receive these vars for each copy stored on my_pe
         ALL_MY_COPIES_RECV_LOOP: do k = 1, my_num_copies

            call receive_from(map_pe_to_task(ens_handle, sending_pe), transfer_temp(1:num_vars_to_receive))
            ! Copy the transfer array to my local storage
            do sv = 1, num_vars_to_receive
               ens_handle%vars(var_list(sv), k) = transfer_temp(sv)
            enddo

         enddo ALL_MY_COPIES_RECV_LOOP
      endif

   else

      do recv_pe = 0, num_pes - 1
      ! I'm the sending PE, figure out what copies of my vars I'll send
      call get_copy_list(num_copies, recv_pe, copy_list, num_copies_to_send)

         SEND_COPIES: do copy = 1, num_copies_to_send
            if (my_pe /= recv_pe ) then
               if (my_num_vars > 0) then
                  transfer_temp(1:my_num_vars) = ens_handle%copies(copy_list(copy), :)
                  ! Have to  use temp because %copies section is not contiguous storage
                  call send_to(map_pe_to_task(ens_handle, recv_pe), transfer_temp(1:my_num_vars))
               endif

            else

               ! figure out what piece to recieve from myself and recieve it
               call get_var_list(num_vars, sending_pe, var_list, num_vars_to_receive)
               do k = 1,  my_num_copies
                  ! sending to yourself so just copy
                  global_ens_index = ens_handle%my_copies(k)
                  do sv = 1, num_vars_to_receive
                     ens_handle%vars(var_list(sv), k) = ens_handle%copies(global_ens_index, sv)
                  end do
               enddo
            endif
         enddo SEND_COPIES
      enddo

   endif

enddo SENDING_PE_LOOP

else ! use old communication pattern

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
                  call receive_from(map_pe_to_task(ens_handle, sending_pe), transfer_temp(1:num_vars_to_receive))
   
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
            call send_to(map_pe_to_task(ens_handle, recv_pe), transfer_temp(1:my_num_vars))
         endif
      end do
      
   endif
end do RECEIVING_PE_LOOP

endif


! Free up the temporary storage
deallocate(var_list, transfer_temp, copy_list)

! only output if there is a label
if (present(label)) then
   call timestamp_message('copies_to_vars   end: '//label, alltasks=.true.)
endif

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

! Error checking
if (ens_handle%valid /= VALID_COPIES .and. ens_handle%valid /= VALID_BOTH) then
   call error_handler(E_ERR, 'compute_copy_mean', &
        'last access not copy-complete', source, revision, revdate)
endif

num_copies = end_copy - start_copy + 1

MYLOOP : do i = 1, ens_handle%my_num_vars
   if (any(ens_handle%copies(start_copy:end_copy, i) == MISSING_R8)) then
      ens_handle%copies(mean_copy, i) = MISSING_R8
   else
      ens_handle%copies(mean_copy, i) = sum(ens_handle%copies(start_copy:end_copy, i)) / num_copies
   endif
end do MYLOOP

ens_handle%valid = VALID_COPIES

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

! Error checking
!if (ens_handle%valid /= VALID_COPIES .and. ens_handle%valid /= VALID_BOTH) then
!   call error_handler(E_ERR, 'compute_copy_mean_sd', &
!        'last access not copy-complete', source, revision, revdate)
!endif

num_copies = end_copy - start_copy + 1

MYLOOP : do i = 1, ens_handle%my_num_vars

   if (any(ens_handle%copies(start_copy:end_copy, i) == MISSING_R8)) then
      ens_handle%copies(mean_copy, i) = MISSING_R8
      ens_handle%copies(  sd_copy, i) = MISSING_R8
   else
      ens_handle%copies(mean_copy, i) = sum(ens_handle%copies(start_copy:end_copy, i)) / num_copies
      ens_handle%copies(  sd_copy, i) = sqrt((sum((ens_handle%copies(start_copy:end_copy, i) - &
                                        ens_handle%copies(mean_copy, i))**2) / (num_copies - 1)))
   endif

end do MYLOOP

ens_handle%valid = VALID_COPIES

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

! Error checking
if (ens_handle%valid /= VALID_COPIES .and. ens_handle%valid /= VALID_BOTH) then
   call error_handler(E_ERR, 'compute_copy_mean_var', &
        'last access not copy-complete', source, revision, revdate)
endif

num_copies = end_copy - start_copy + 1

MYLOOP : do i = 1, ens_handle%my_num_vars
   if (any(ens_handle%copies(start_copy:end_copy, i) == MISSING_R8)) then
      ens_handle%copies(mean_copy, i) = MISSING_R8
      ens_handle%copies( var_copy, i) = MISSING_R8
   else
      ens_handle%copies(mean_copy, i) = sum(ens_handle%copies(start_copy:end_copy, i)) / num_copies
      ens_handle%copies( var_copy, i) = (sum((ens_handle%copies(start_copy:end_copy, i) - &
         ens_handle%copies(mean_copy, i))**2) / (num_copies - 1))
   endif
end do MYLOOP

ens_handle%valid = VALID_COPIES

end subroutine compute_copy_mean_var

!--------------------------------------------------------------------------------

subroutine timestamp_message(msg, sync, alltasks)

character(len=*), intent(in) :: msg
logical, intent(in), optional :: sync
logical, intent(in), optional :: alltasks

integer, save :: timestamp_level = 1
logical :: should_output, old_flag
character (len=129) :: tbuf

! Write current time and message to stdout and log file.
! if sync is present and true, sync mpi jobs before printing time.

if (timestamp_level <= 0) return

if (present(sync)) then
  if (sync) call task_sync()
endif

should_output = do_output()
if (present(alltasks)) then
   if (alltasks) should_output = .true.
endif

if (should_output) then
   old_flag = do_output()
   call set_output(.true.)
   write(tbuf, "(A,I4,A)") 'Task', my_task_id(), ': '//trim(msg)
   call timestamp(trim(tbuf), pos='brief')  ! was debug
   call set_output(old_flag)
endif

end subroutine timestamp_message

!--------------------------------------------------------------------------------

subroutine print_ens_handle(ens_handle, force, label)
 type(ensemble_type),        intent(in) :: ens_handle
 logical,          optional, intent(in) :: force
 character(len=*), optional, intent(in) :: label

logical :: print_anyway
logical :: has_label

print_anyway = .false.
if (present(force)) then
   print_anyway = force
endif

has_label = .false.
if (present(label)) then
   has_label = .true.
endif

! print out contents of an ensemble handle derived type
if (.not. debug .and. .not. print_anyway) return

if (has_label) then
   call error_handler(E_MSG, 'ensemble handle: ', label, source, revision, revdate)
endif
write(msgstring, *) 'handle num: ',          ens_handle%id_num 
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
write(msgstring, *) 'number of    copies: ', ens_handle%num_copies
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
write(msgstring, *) 'number of    vars  : ', ens_handle%num_vars
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
write(msgstring, *) 'number of my_copies: ', ens_handle%my_num_copies
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
write(msgstring, *) 'number of my_vars  : ', ens_handle%my_num_vars
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
write(msgstring, *) 'valid              : ', ens_handle%valid
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
write(msgstring, *) 'distribution_type  : ', ens_handle%distribution_type
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
write(msgstring, *) 'my_pe number       : ', ens_handle%my_pe
call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
if (allocated(ens_handle%pe_to_task_list)) then
   write(msgstring, *) 'task_to_pe_list    : ', ens_handle%task_to_pe_list
   call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
   write(msgstring, *) 'pe_to_task_list    : ', ens_handle%pe_to_task_list
   call error_handler(E_MSG, 'ensemble handle: ', msgstring, source, revision, revdate)
endif


end subroutine print_ens_handle

!--------------------------------------------------------------------------------

subroutine assign_tasks_to_pes(ens_handle, nEns_members, layout_type)

! Calulate the task layout based on the tasks per node and the total number of tasks.
! Allows the user to spread out the ensemble members as much as possible to balance 
! memory usage between nodes.
!
! Possible options:
!   1. Standard task layout - first n tasks have the ensemble members my_pe = my_task_id()
!   2. Round-robin on the nodes

type(ensemble_type), intent(inout)    :: ens_handle
integer,             intent(in)       :: nEns_members
integer,             intent(inout)    :: layout_type

if (layout_type /= 1 .and. layout_type /=2) call error_handler(E_ERR,'assign_tasks_to_pes', &
    'not a valid layout_type, must be 1 (standard) or 2 (round-robin)',source,revision,revdate)

if (tasks_per_node >= num_pes) then ! all tasks are on one node, don't try to spread them out
   call simple_layout(ens_handle, num_pes)
   return
endif

if (layout_type == 1) then 
   call simple_layout(ens_handle, num_pes)
else
   call round_robin(ens_handle)
endif

end subroutine assign_tasks_to_pes

!------------------------------------------------------------------------------

subroutine round_robin(ens_handle)

! Round-robin MPI task layout starting at the first node.  
! Starting on the first node forces pe 0 = task 0. 
! The smoother code assumes task 0 has an ensemble member.
! If you want to break the assumption that pe 0 = task 0, this routine is a good 
! place to start. Test with the smoother. 

type(ensemble_type), intent(inout)   :: ens_handle
integer                              :: last_node_task_number, num_nodes
integer                              :: i, j
integer, allocatable                 :: count(:)


! Find number of nodes and find number of tasks on last node
call calc_tasks_on_each_node(num_nodes, last_node_task_number)

allocate(count(num_nodes))

count(:) = 1  ! keep track of the pes assigned to each node
i = 0         ! keep track of the # of pes assigned

do while (i < num_pes)   ! until you run out of processors
   do j = 1, num_nodes   ! loop around the nodes

      if(j == num_nodes) then  ! special case for the last node - it could have fewer tasks than the other nodes
         if(count(j) <= last_node_task_number) then
            ens_handle%task_to_pe_list(tasks_per_node*(j-1) + count(j)) = i
            count(j) = count(j) + 1
            i = i + 1
         endif
      else
         if(count(j) <= tasks_per_node) then
            ens_handle%task_to_pe_list(tasks_per_node*(j-1) + count(j)) = i
            count(j) = count(j) + 1
            i = i + 1
         endif
      endif

   enddo
enddo

deallocate(count)

call create_pe_to_task_list(ens_handle)

end subroutine round_robin

!-------------------------------------------------------------------------------

subroutine create_pe_to_task_list(ens_handle)

! Creates the ens_handle%pe_to_task_list
! ens_handle%task_to_pe_list must have been assigned first, otherwise this 
! routine will just return nonsense. 

!FIXME set ens_handle%task_to_pe_list to -1 when it is allocated, then test if has been changed

type(ensemble_type), intent(inout)   :: ens_handle
integer                              :: temp_sort(num_pes), idx(num_pes)
integer                              :: ii

temp_sort = ens_handle%task_to_pe_list
call sort_task_list(temp_sort, idx, num_pes)

do ii = 1, num_pes
   ens_handle%pe_to_task_list(ii) = temp_sort(idx(ii))
enddo

end subroutine create_pe_to_task_list

!-------------------------------------------------------------------------------

subroutine calc_tasks_on_each_node(nodes, last_node_task_number)

! Finds the of number nodes and how many tasks are on the last node, given the 
! number of tasks and the tasks_per_node (ptile).
! The total number of tasks is num_pes = task_count()
! The last node may have fewer tasks, for example, if ptile = 16 and the number of
! mpi tasks = 17

integer, intent(out)  :: last_node_task_number, nodes

if ( mod(num_pes, tasks_per_node) == 0) then
   nodes = num_pes / tasks_per_node
   last_node_task_number = tasks_per_node
else
   nodes = num_pes / tasks_per_node + 1
   last_node_task_number = tasks_per_node - (nodes*tasks_per_node - num_pes)
endif

end subroutine calc_tasks_on_each_node

!-----------------------------------------------------------------------------

subroutine simple_layout(ens_handle, n)

! assigns the arrays task_to_pe_list and pe_to_task list for the simple layout
! where my_pe = my_task_id()

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: n
integer                            :: ii

do ii = 0, num_pes - 1
   ens_handle%task_to_pe_list(ii + 1) = ii
enddo

ens_handle%pe_to_task_list = ens_handle%task_to_pe_list

end subroutine simple_layout

!------------------------------------------------------------------------------
!> sorts an array and returns the sorted array, and the index of the original
!> array
subroutine sort_task_list(x, idx, n)

integer, intent(in)    :: n 
integer, intent(inout) :: x(n)   ! array to be sorted
integer, intent(out)   :: idx(n) ! index of sorted array

integer                :: xcopy(n), i

xcopy = x

call index_sort(x, idx, n)

do i = 1, n
   x(i) = xcopy(idx(i))
enddo

end subroutine sort_task_list

!!--------------------------------------------------------------------------------
!
!function map_pe_to_task(ens_handle, p)
!
!! Return the physical task for my_pe
!
!type(ensemble_type), intent(in) :: ens_handle
!integer,             intent(in)    :: p
!integer                            :: map_pe_to_task
!
!map_pe_to_task = ens_handle%pe_to_task_list(p + 1)
!
!end function map_pe_to_task

!--------------------------------------------------------------------------------
!> ! Return my_pe corresponding to the physical task
function map_task_to_pe(ens_handle, t)

type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: t
integer                         :: map_task_to_pe

map_task_to_pe = ens_handle%task_to_pe_list(t + 1)

end function map_task_to_pe

!---------------------------------------------------------------------------------

end module ensemble_manager_mod


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
