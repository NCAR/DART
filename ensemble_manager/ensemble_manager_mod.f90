! Data Assimilation Research Testbed -- DART

! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module ensemble_manager_mod

use types_mod, only : r8
use    utilities_mod, only :  get_unit, open_file, close_file, register_module, &
                              check_nml_error, file_exist, error_handler, &
                              E_ERR, E_WARN, E_MSG, E_DBG, initialize_utilities, &
                              logfileunit, timestamp
use assim_model_mod, only : aread_state_restart, awrite_state_restart, open_restart_read, &
   open_restart_write, close_restart, get_model_time_step, adv_1step
use time_manager_mod, only : time_type, get_time, read_time, write_time, get_calendar_type, &
                             THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP, NO_CALENDAR, &
                             operator(<), operator(>), operator(+), operator(-), &
                             operator(/), operator(*), operator(==), operator(/=)

implicit none
private

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

public init_ensemble_manager, get_ensemble_member, put_ensemble_member, &
   update_ens_mean, update_ens_mean_spread, end_ensemble_manager, &
   get_ensemble_region, put_ensemble_region, get_ensemble_time, Aadvance_state, &
   ensemble_type

! This type gives a handle to an ensemble, not currently playing a role but
! allows later implementations to possibly support multiple ensembles open at once
type ensemble_type
   logical :: null_variable
end type ensemble_type

! Flags for allocating mean and spread storage
logical :: mean_allocated = .false., spread_allocated = .false.
   
! Global in core storage for the ensemble
real(r8), allocatable :: ens(:, :), ens_mean(:), ens_spread(:)
type(time_type), allocatable :: ens_time(:)
integer :: ens_size, model_size

integer :: seq_unit
character(len = 129) :: errstring
! File name for the temporary files; has extensions added
character(len = 21) :: seq_file_name = 'ens_manager_temp_file'

!-----------------------------------------------------------------
!
! namelist with default values

! If true, then ensemble is held in core storage; if false, only
! a single ensemble member worth of storage is used and ensemble
! is kept in a direct access temporary file.
logical :: in_core = .true.
! If true a single restart file containing all ensemble states is read/written
! If false, one restart file is used for each ensemble member
logical :: single_restart_file_in = .true.
logical :: single_restart_file_out = .true.

namelist / ensemble_manager_nml / in_core, single_restart_file_in, single_restart_file_out

!-----------------------------------------------------------------

contains

!-----------------------------------------------------------------

subroutine init_ensemble_manager(ens_handle, ens_size_in, &
   model_size_in, file_name, init_time)

type(ensemble_type), intent(out)            :: ens_handle
integer, intent(in)                         :: ens_size_in, model_size_in
character(len = 129), intent(out), optional :: file_name
type(time_type), intent(in), optional       :: init_time

integer :: iunit, i, ierr, io
character(len = 129) :: msgstring, this_file_name
character(len = 4) :: extension

! Initialize the module with utilities 
call register_module(source, revision, revdate)

! Set the global storage bounds
ens_size = ens_size_in
model_size = model_size_in

! Read namelist for run time control
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1

   READBLOCK: do while(ierr /= 0)
      read(iunit, nml = ensemble_manager_nml, iostat = io)
      if ( io < 0 ) exit READBLOCK          ! end-of-file
      ierr = check_nml_error(io, 'ensemble_manager_nml')
   enddo READBLOCK

   call close_file(iunit)
endif

call error_handler(E_MSG,'init_ensemble_manager','ensemble_manager_nml values are',' ',' ',' ')
write(logfileunit,nml=ensemble_manager_nml)
write(     *     ,nml=ensemble_manager_nml)
! Initialize the storage and read in from restart file if needed
if(in_core) then 
   ! Don't allocate space for mean and spread until needed to conserve space
   allocate(ens(ens_size, model_size), ens_time(-1:ens_size))
else
   ! For out of core will only have storage for one ensemble member
   allocate(ens(1, model_size), ens_time(-1:ens_size))
endif

! If a file_name is present, read in ensemble from (a) restart file(s)
if(present(file_name)) then
   if(single_restart_file_in) iunit = open_restart_read(file_name)

   do i = 1, ens_size
      ! If multiple restart files, open one for each ensemble member
      if(.not. single_restart_file_in) then
         write(extension, 99) i
         99 format(i4.4)
         this_file_name = trim(file_name) // '.' // extension
         iunit = open_restart_read(this_file_name)
      endif

      write(msgstring, *) 'trying to read restart ', i
      call error_handler(E_DBG,'init_ensemble_manager',msgstring,source,revision,revdate)

      ! If in core, can just read the ensemble members straight into storage
      if(in_core) then 
         call aread_state_restart(ens_time(i), ens(i, :), iunit)
      else
         ! If out of core, read in from restart file and write out to temp file
         call aread_state_restart(ens_time(i), ens(1, :), iunit)
         ! Open the file for this ensemble member
         call get_seq_file_name(this_file_name, i)
         seq_unit = get_unit()
         open(unit = seq_unit, file = this_file_name, access = 'sequential', &
            status = 'replace', form = 'unformatted')
         write(seq_unit) ens(1, :)
         close(seq_unit)
      endif

      ! If the inital file time is being overridden, do it
      if(present(init_time)) ens_time(i) = init_time

      if(.not. single_restart_file_in) call close_restart(iunit)
   end do
   
   if(single_restart_file_in) call close_restart(iunit)
endif

end subroutine init_ensemble_manager

!-----------------------------------------------------------------

subroutine get_ensemble_member(ens_handle, index, member, mtime)

type(ensemble_type), intent(in) :: ens_handle
integer, intent(in) :: index
real(r8), intent(out) :: member(:)
type(time_type), intent(out) :: mtime

! Return one member of an ensemble by index
! Index 0 gets ensemble mean, index -1 gets ensemble spread

character(len = 129) :: this_file_name

if(index < -1 .or. index > ens_size) then
   write(errstring, *) 'index out of range'
   call error_handler(E_ERR,'get_ensemble_member', errstring, source, revision, revdate)
endif

! For in memory, just copy and return
if(in_core) then
   if(index > 0) then
      member = ens(index, :)
   else if(index == 0) then
      member = ens_mean
   else if(index == -1) then 
      member = ens_spread
   endif
else
   ! Out of core, need to read in from file
   call get_seq_file_name(this_file_name, index)
   ! Open the file and read in the field      
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   read(seq_unit) member
   close(unit = seq_unit)
endif

mtime = ens_time(index)

end subroutine get_ensemble_member

!-----------------------------------------------------------------

subroutine put_ensemble_member(ens_handle, index, member, mtime)

! Set one member of an ensemble by index
! Index 0 gets ensemble mean, index -1 gets ensemble spread

type(ensemble_type), intent(in) :: ens_handle
integer, intent(in) :: index
real(r8), intent(in) :: member(:)
type(time_type), intent(in) :: mtime

character(len = 129) :: this_file_name

if(index < -1 .or. index > ens_size) then
   write(errstring, *) 'index out of range ', index
   call error_handler(E_ERR,'put_ensemble_member', errstring, source, revision, revdate)
endif

! For in core, just copy
if(in_core) then
   if(index > 0) then
      ens(index, :) = member
   else if(index == 0) then
      ens_mean = member
   else
      ens_spread = member
   endif
else
   ! For out of core, open file and write out new values
   call get_seq_file_name(this_file_name, index)
   ! Open the file and write the field
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   write(seq_unit) member
   close(unit = seq_unit)
end if

ens_time(index) = mtime

end subroutine put_ensemble_member

!-----------------------------------------------------------------

subroutine get_ensemble_time(ens_handle, index, mtime)

! Returns the time of an ensemble member; 0 for ensemble mean, -1 for spread

type(ensemble_type), intent(in) :: ens_handle
integer, intent(in) :: index
type(time_type), intent(out) :: mtime

if(index < -1 .or. index > ens_size) then
   write(errstring, *) 'index out of range ', index
   call error_handler(E_ERR,'get_ensemble_time', errstring, source, revision, revdate)
endif

mtime = ens_time(index)

end subroutine get_ensemble_time

!-----------------------------------------------------------------

subroutine update_ens_mean(ens_handle)

! Update the value of the ensemble mean from the current ensemble members

type(ensemble_type), intent(in) :: ens_handle

integer :: i, j
character(len = 129) :: this_file_name
! Sequential version read efficiency requires having second model_size storage
real(r8) :: ens_member(model_size)

! Update the ensemble mean (copy 0) from ensemble members
if(in_core) then
   ! Make sure that space has been allocated
   if(.not. mean_allocated) then
      allocate(ens_mean(model_size))  
      mean_allocated = .true.
   endif
   do i = 1, model_size
      ens_mean(i) = sum(ens(1:ens_size, i)) / ens_size
   end do
else
   ! Stored on disk, read in each member and sum
   ens = 0.0
   do j = 1, ens_size
      call get_seq_file_name(this_file_name, j)
      ! Open the file and read in the field      
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
      read(seq_unit) ens_member
      close(unit = seq_unit)
      ! Add this member to sum
      ens(1, :) = ens(1, :) + ens_member
   end do
   ens = ens / ens_size
   ! Put the ensemble mean into the file; ensemble mean is 0 member
   call get_seq_file_name(this_file_name, 0)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   write(seq_unit) ens
   close(unit = seq_unit)
endif

! Update the ensemble mean time, assume everything else stays in synch
ens_time(0) = ens_time(1)

end subroutine update_ens_mean

!-----------------------------------------------------------------

subroutine update_ens_mean_spread(ens_handle)

! Update the ensemble spread from the current ensemble members,
! Also does an update of the ensemble mean
type(ensemble_type), intent(in) :: ens_handle

integer :: i, j
character(len = 129) this_file_name
! For efficient disk use now need three copies, if this is too much can alter
real(r8) :: tens_mean(model_size), tens_spread(model_size)

! Update both the mean and spread (copy -1) from ensemble members
call update_ens_mean(ens_handle)

if(in_core) then
   ! Make sure storage has been allocated
   if(.not. spread_allocated) then
      allocate(ens_spread(model_size))
      spread_allocated = .true.
   endif
   do i = 1, model_size
      ens_spread(i) = sqrt(sum((ens(1:, i) - ens_mean(i))**2) / (ens_size - 1))
   end do
else
   ! Get the ensemble mean from disk file; mean is copy 0
   call get_seq_file_name(this_file_name, 0)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   read(seq_unit) tens_mean
   close(unit = seq_unit)

   tens_spread = 0.0
   ! Loop through and add squared difference from mean
      do j = 1, ens_size
      call get_seq_file_name(this_file_name, j)
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
      read(seq_unit) ens
      close(unit = seq_unit)
      tens_spread = tens_spread + (ens(1, :) - tens_mean)**2
   end do
   
   tens_spread = sqrt(tens_spread / (ens_size - 1))
   ! Write to the ensemble_spread file, index is -1
   call get_seq_file_name(this_file_name, -1)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   write(seq_unit) tens_spread
   close(unit = seq_unit)
endif

! Update the spread time
ens_time(-1) = ens_time(1)

end subroutine update_ens_mean_spread

!-----------------------------------------------------------------

subroutine get_ensemble_region(ens_handle, region, rtime, state_vars_in, ens_members_in)

type(ensemble_type), intent(in) :: ens_handle
real(r8), intent(out) :: region(:, :)
type(time_type), intent(out) :: rtime(:)
integer, intent(in), optional :: state_vars_in(:), ens_members_in(:)

integer :: state_vars(model_size), ens_members(ens_size + 2), i, j
integer :: num_state_vars, num_ens_members
character(len = 129) :: this_file_name
! Returns a subset of the state variables and ensemble members
! The ensemble members to be returned are specified by number in ens_members_in
! The ensemble mean is number 0 and ensemble spread is -1
! The state variables to be returned are specified by number in state_var_in

! NOTE: Leaving out ens_members_in results in return of all base
! ensemble members but NOT the ensemble mean or spread

! Determine size and elements of state variable subset
if(present(state_vars_in)) then
   num_state_vars = size(state_vars_in)
   state_vars(1:num_state_vars) = state_vars_in
else
   num_state_vars = model_size
   do i = 1, model_size
      state_vars(i) = i
   end do
endif

! Determine size and elements of ensemble subset
if(present(ens_members_in)) then
   num_ens_members = size(ens_members_in)
   ens_members(1:num_ens_members) = ens_members_in
else
   num_ens_members = ens_size
   do i = 1, ens_size
      ens_members(i) = i
   end do
endif

do i = 1, num_ens_members
   rtime(i) = ens_time(ens_members(i))

   ! If not in core, open the file for this ensemble member
   if(.not. in_core) then
      call get_seq_file_name(this_file_name, ens_members(i))
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
      read(seq_unit) ens
      close(seq_unit)
   endif

   do j = 1, num_state_vars
      if(in_core) then
         if(ens_members(i) == 0) then
            region(i, j) = ens_mean(state_vars(j))
         else if(ens_members(i) == -1) then
            region(i, j) = ens_spread(state_vars(j))
         else
            region(i, j) = ens(ens_members(i), state_vars(j))
         endif
      else
         ! Already have proper ensemble member read in for out of core
         region(i, j) = ens(1, state_vars(j))
      endif
   end do
end do

end subroutine get_ensemble_region

!-----------------------------------------------------------------

subroutine put_ensemble_region(ens_handle, region, rtime, state_vars_in, ens_members_in)

type(ensemble_type), intent(in) :: ens_handle
real(r8), intent(in) :: region(:, :)
type(time_type), intent(in) :: rtime(:)
integer, intent(in), optional :: state_vars_in(:), ens_members_in(:)

integer :: state_vars(model_size), ens_members(ens_size + 2), i, j
integer :: num_state_vars, num_ens_members
character(len = 129) :: this_file_name
! Puts a subset of state variables and ensemble members

! Determine size and elements of state variable subset
if(present(state_vars_in)) then
   num_state_vars = size(state_vars_in)
   state_vars(1:num_state_vars) = state_vars_in
else
   num_state_vars = model_size
   do i = 1, model_size
      state_vars(i) = i
   end do
endif

! Determine size and elements of ensemble subset
if(present(ens_members_in)) then
   num_ens_members = size(ens_members_in)
   ens_members(1:num_ens_members) = ens_members_in
else
   num_ens_members = ens_size
   do i = 1, ens_size
      ens_members(i) = i
   end do
endif

do i = 1, num_ens_members
   ens_time(ens_members(i)) = rtime(i)

   ! If not in core, open the file for this ensemble member
   if(.not. in_core) then
      call get_seq_file_name(this_file_name, ens_members(i))
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
      read(seq_unit) ens
      close(seq_unit)
   endif

   do j = 1, num_state_vars
      if(in_core) then 
         if(ens_members(i) == 0) then
            ens_mean(state_vars(j)) = region(i, j)
         else if(ens_members(i) == -1) then
            ens_spread(state_vars(j)) = region(i, j)
         else
            ens(ens_members(i), state_vars(j)) = region(i, j)
         endif
      else
         ens(1, state_vars(j)) = region(i, j)
         endif
   end do

   ! If not in core, open the file for this ensemble member
   if(.not. in_core) then
      call get_seq_file_name(this_file_name, ens_members(i))
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
      write(seq_unit) ens
      close(seq_unit)
   endif
end do

end subroutine put_ensemble_region

!-----------------------------------------------------------------

subroutine end_ensemble_manager(ens_handle, file_name)

type(ensemble_type), intent(in) :: ens_handle
character(len = 129), intent(in), optional :: file_name

integer :: iunit, i, seq_unit
character(len = 129) :: this_file_name, command_string
character(len = 4) :: extension

if(single_restart_file_out) iunit = open_restart_write(file_name)

do i = 1, ens_size
   
   if(.not. single_restart_file_out) then
      write(extension, 99) i
      99 format(i4.4)
      this_file_name = trim(file_name) // '.' // extension
      iunit = open_restart_write(this_file_name)
   endif

   if(in_core) then
      call awrite_state_restart(ens_time(i), ens(i, :), iunit)
   else
      call get_seq_file_name(this_file_name, i)
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
      read(seq_unit) ens
      close(seq_unit)
      call awrite_state_restart(ens_time(i), ens(1, :), iunit)
   endif

   if(.not. single_restart_file_out) call close_restart(iunit)

end do

if(single_restart_file_out) call close_restart(iunit)

! Free up the allocated storage
deallocate(ens, ens_time)
if(mean_allocated) deallocate(ens_mean)
if(spread_allocated) deallocate(ens_spread)

! LETS GET RID OF ALL THE TEMP FILES, TOO!
command_string = 'rm -f ' // trim(seq_file_name) // '.*'
!write(*, *) 'command string is ', command_string
call system(command_string)

end subroutine end_ensemble_manager

!-----------------------------------------------------------------






subroutine Aadvance_state(ens_handle, target_time, asynch, adv_ens_command)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.

implicit none

type(ensemble_type),intent(in)    :: ens_handle
type(time_type),    intent(in)    :: target_time
integer,            intent(in)    :: asynch
character(len=129), intent(in)    :: adv_ens_command                                                      

type(time_type) :: time_step

integer :: seconds, days, i, control_unit, ic_file_unit, ud_file_unit

character(len = 26), dimension(ens_size) :: ic_file_name, ud_file_name 
character(len = 128) :: input_string
integer :: is1,is2,id1,id2
type(time_type) :: smodel_time
real(r8) :: smodel_state(model_size)

! Loop through each model state and advance
do i = 1, ens_size

   ! Get the state and time for this ensemble member
   call get_ensemble_member(ens_handle, i, smodel_state, smodel_time)

   ! On first one, make sure that there is a need to advance, else return
   if(smodel_time == target_time) return

   ! Check for time error
   if(smodel_time > target_time) then
      call get_time(smodel_time, is1, id1)
      call get_time(target_time, is2, id2)
      write(errstring,*)'target time ',is2,id2,' is before model_time ',is1,id1
      call error_handler(E_ERR,'Aadvance_state', errstring, source, revision, revdate)
   endif

! Two blocks here for now, one for single executable, one for asynch multiple execs

!------------- Block for single executable ----------------------------
! At some point probably need to push the determination of the time back
! into the model itself and out of assim_model

   if(asynch == 0) then

      time_step = get_model_time_step()
      call get_time(time_step, seconds, days)

      do while(smodel_time < target_time)
         call adv_1step(smodel_state, smodel_time)
         smodel_time = smodel_time + time_step
      end do

      ! Put the updated ensemble member back into storage
      call put_ensemble_member(ens_handle, i, smodel_state, smodel_time)

   !-------------- End single executable block ------------------------
   else 
   !-------------- Block for multiple asynch executables --------------

   ! Loop to write out state for each member to separate file, preface with target time
      if(i < 10) then
         write(ic_file_name(i), 11) 'assim_model_state_ic', i
         write(ud_file_name(i), 11) 'assim_model_state_ud', i
      else if(i < 100) then
         write(ic_file_name(i), 21) 'assim_model_state_ic', i
         write(ud_file_name(i), 21) 'assim_model_state_ud', i
      else if(i < 1000) then
         write(ic_file_name(i), 31) 'assim_model_state_ic', i
         write(ud_file_name(i), 31) 'assim_model_state_ud', i
      else if(i < 10000) then
         write(ic_file_name(i), 41) 'assim_model_state_ic', i
         write(ud_file_name(i), 41) 'assim_model_state_ud', i
      else 
         write(errstring,*)'Trying to use ',ens_size,' model states -- too many.'
         call error_handler(E_MSG,'Aadvance_state',errstring,source,revision,revdate)
         call error_handler(E_ERR,'Aadvance_state','Use less than 10000 model states.',source,revision,revdate)
      endif

 11   format(a21, i1)
 21   format(a21, i2)
 31   format(a21, i3)
 41   format(a21, i4)
      write(*, *) 'ic and ud files ', i, ic_file_name(i), ud_file_name(i)
      write(logfileunit, *) 'ic and ud files ', i, ic_file_name(i), ud_file_name(i)

      ! Output the destination time followed by assim_model_state 
      ! Write the time to which to advance
      ! Write the assim model extended state   

      ! Open a restart file and write out with a target time
      ic_file_unit = open_restart_write(ic_file_name(i))
      call awrite_state_restart(smodel_time, smodel_state, ic_file_unit, target_time)
      call close_restart(ic_file_unit)
     

   endif

   !-------------- End of multiple async executables block ------------

end do


! Also need synchronization block at the end for the asynch

if(asynch /= 0) then

   ! Write out the file names to a control file

   control_unit = get_unit()
   open(unit = control_unit, file = 'filter_control')
   write(control_unit, *) ens_size
   do i = 1, ens_size
      write(control_unit, '(a26)' ) ic_file_name(i)
      write(control_unit, '(a26)' ) ud_file_name(i)
   end do
!!!   call write_time(control_unit, target_time)
   close(control_unit)

   if(asynch == 1) then

      ! We sleep 1 second here for the following reason:
      ! The first time async_filter.csh is executed,
      ! it removes the file async_may_go that may be left over.
      ! With small models, perfect_model_obs and filter can come
      ! to this point before the script async_filter.csh has time to
      ! remove the file async_may_go.

      call system('sleep 1')

      ! Create the file async_may_go to allow the async model integrations
      control_unit = get_unit()
      open(unit = control_unit, file = 'async_may_go')
      call write_time(control_unit, target_time)
      close(control_unit)

! Suspend on a read from standard in for integer value
      do
         read(*, '(a80)') input_string
         if(trim(input_string) == 'All_done:Please_proceed') exit
! Following line can allow diagnostic pass through of output
         write(*, *) input_string
      end do

      ! This sleep is necessary to insure that all files coming from
      ! remote systems are completely written.

      call system('sleep 1')

   elseif(asynch == 2) then

!      call system('./sync_filter.csh ; sleep 1')
      call system('echo go > batchflag; '//adv_ens_command//' ; sleep 1')

      do
         if(file_exist('batchflag')) then
            call system ('sleep 10')
         else
!            call system ('echo "All_done:Please_proceed"')
            exit
         endif
      end do
   
   elseif(asynch == 3) then
   ! Script is running all the time, waiting for clearance to proceed
      call system ('echo a > go_advance_model')
      do
         if(file_exist('go_advance_model')) then
            call system('sleep 1')
         else
            exit
         endif
      end do

   else
   ! Unsupported option for async error
      write(errstring,*)'input.nml - async is ',asynch,' must be 0, 1, 2, or 3' 
      call error_handler(E_ERR,'Aadvance_state', errstring, source, revision, revdate)

   endif

   write(*, *) 'got clearance to proceed in Aadvance_state'

   ! All should be done, read in the states and proceed
   do i = 1, ens_size
      ud_file_unit = open_restart_read(ud_file_name(i))
      call aread_state_restart(smodel_time, smodel_state, ud_file_unit)
      call close_restart(ud_file_unit)
      ! Put the updated state in the ensemble storage
      call put_ensemble_member(ens_handle, i, smodel_state, smodel_time)
   end do

end if

end subroutine Aadvance_state


!-----------------------------------------------------------------

subroutine get_seq_file_name(file_name, index)

character(len = *), intent(out) :: file_name
integer, intent(in) :: index

character(len = 4) extension

if(index > 0) then
   ! Open the direct access file for this ensemble member
   write(extension, 99) index
   99 format(i4.4)
   file_name = trim(seq_file_name) // '.' // extension
else if(index == 0) then
   file_name = trim(seq_file_name) // '.mean'
else if(index == -1) then
   file_name = trim(seq_file_name) // '.spread'
endif

end subroutine get_seq_file_name


end module ensemble_manager_mod

