! Data Assimilation Research Testbed -- DART

! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module ensemble_manager_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use types_mod, only : r8
use    utilities_mod, only :  get_unit, open_file, close_file, register_module, &
                              file_exist, error_handler, &
                              E_ERR, E_WARN, E_MSG, E_DBG, initialize_utilities, &
                              logfileunit, timestamp, find_namelist_in_file, check_namelist_read
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

public :: init_ensemble_manager, get_ensemble_member, put_ensemble_member, &
   update_ens_mean, update_ens_mean_spread, end_ensemble_manager, &
   get_ensemble_region, put_ensemble_region, get_ensemble_time, Aadvance_state, &
   ensemble_type, get_region_by_number, put_region_by_number, &
   transpose_ens_to_regions, transpose_regions_to_ens, &
   is_ens_in_core, copy_ens

! This type gives a handle to an ensemble, not currently playing a role but
! allows later implementations to possibly support multiple ensembles open at once
type ensemble_type
   ! Allow direct access into this for efficient in_core use
   !!!private
   integer :: id
   integer :: ens_size, model_size
   real(r8), pointer :: ens(:, :), mean(:), spread(:)
   type(time_type), pointer :: time(:)
   logical :: in_core
end type ensemble_type

! Counter for unique id's for ensembles
integer :: num_ens_init = 0

! Logical flag for initialization of module
logical :: module_initialized = .false.

character(len = 129) :: errstring
! File name for the temporary files; has extensions added
character(len = 20) :: ens_file_name = 'ens_manager_ens_file'
character(len = 20) :: reg_file_name = 'ens_manager_reg_file'

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

subroutine init_ensemble_manager(ens_handle, ens_size, &
   model_size, file_name, init_time)

type(ensemble_type), intent(out)            :: ens_handle
integer, intent(in)                         :: ens_size, model_size
character(len = 129), intent(in), optional  :: file_name
type(time_type), intent(in), optional       :: init_time

real(r8) :: ens(model_size)
integer :: iunit, i, io, seq_unit
character(len = 129) :: msgstring, this_file_name
character(len = 4) :: extension

! integer :: secs, days       ! for print statement below
! type(time_type) :: ensemble_valid_time

if ( .not. module_initialized ) then
   ! Initialize the module with utilities 
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

! Set the global storage bounds
ens_handle%ens_size = ens_size
ens_handle%model_size = model_size

! Set the unique id
num_ens_init = num_ens_init + 1
ens_handle%id = num_ens_init

! Read the namelist entry
call find_namelist_in_file("input.nml", "ensemble_manager_nml", iunit)
read(iunit, nml = ensemble_manager_nml, iostat = io)
call check_namelist_read(iunit, io, "ensemble_manager_nml")

call error_handler(E_MSG,'init_ensemble_manager','ensemble_manager_nml values are',' ',' ',' ')
write(logfileunit,nml=ensemble_manager_nml)
write(     *     ,nml=ensemble_manager_nml)

! Namelist controls the in_core status 
ens_handle%in_core = in_core

! Initialize the storage and read in from restart file if needed
if(ens_handle%in_core) then 
   allocate(ens_handle%ens(ens_size, model_size), ens_handle%mean(model_size), &
      ens_handle%spread(model_size), ens_handle%time(-1:ens_size))
else
   ! For out of core will only have storage for one ensemble member
   allocate(ens_handle%time(-1:ens_size))
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
      if(ens_handle%in_core) then 
         call aread_state_restart(ens_handle%time(i), ens_handle%ens(i, :), iunit)
      else
         ! If out of core, read in from restart file and write out to temp file
         call aread_state_restart(ens_handle%time(i), ens, iunit)
         ! Open the file for this ensemble member
         this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, i)
         seq_unit = get_unit()
         open(unit = seq_unit, file = this_file_name, access = 'sequential', &
            status = 'replace', form = 'unformatted')

         !open(unit = seq_unit, file = this_file_name, access = 'direct', &
         !   form = 'unformatted', status = 'replace', recl = req_rec_length)

         write(seq_unit) ens
         close(seq_unit)
      endif

      ! If the initial file time is being overridden, do it
      if(present(init_time)) ens_handle%time(i) = init_time

      if(.not. single_restart_file_in) call close_restart(iunit)

      ! Print of initial model time
 !    ensemble_valid_time = ens_handle%time(i)
 !    call get_time(ensemble_valid_time,secs,days)
 !    write(msgstring,*)'ensemble file ',trim(adjustl(this_file_name)), &
 !                      ' (days, seconds)',days,secs
 !    call error_handler(E_MSG,'init_ensemble_manager',msgstring,source,revision,revdate)

   end do
   
   if(single_restart_file_in) call close_restart(iunit)
endif

end subroutine init_ensemble_manager

!-----------------------------------------------------------------

subroutine get_ensemble_member(ens_handle, indx, member, mtime)

type(ensemble_type), intent(in) :: ens_handle
integer, intent(in) :: indx
real(r8), intent(out) :: member(:)
type(time_type), intent(out) :: mtime

! Return one member of an ensemble by indx
! indx 0 gets ensemble mean, indx -1 gets ensemble spread

integer :: seq_unit
character(len = 129) :: this_file_name

if(indx < -1 .or. indx > ens_handle%ens_size) then
   write(errstring, *) 'indx out of range'
   call error_handler(E_ERR,'get_ensemble_member', errstring, source, revision, revdate)
endif

! For in memory, just copy and return
if(ens_handle%in_core) then
   if(indx > 0) then
      member = ens_handle%ens(indx, :)
   else if(indx == 0) then
      member = ens_handle%mean
   else if(indx == -1) then 
      member = ens_handle%spread
   endif
else
   ! Out of core, need to read in from file
   this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, indx)
   ! Open the file and read in the field      
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   
    !open(unit = seq_unit, file = this_file_name, access = 'direct', &
    !        form = 'unformatted', recl = req_rec_length)

   read(seq_unit) member
   close(unit = seq_unit)
endif

mtime = ens_handle%time(indx)

end subroutine get_ensemble_member

!-----------------------------------------------------------------

subroutine put_ensemble_member(ens_handle, indx, member, mtime)

! Set one member of an ensemble by index
! indx 0 gets ensemble mean, indx -1 gets ensemble spread

type(ensemble_type), intent(inout) :: ens_handle
integer, intent(in) :: indx
real(r8), intent(in) :: member(:)
type(time_type), intent(in) :: mtime

integer :: seq_unit
character(len = 129) :: this_file_name

if(indx < -1 .or. indx > ens_handle%ens_size) then
   write(errstring, *) 'indx out of range ', indx
   call error_handler(E_ERR,'put_ensemble_member', errstring, source, revision, revdate)
endif

! For in core, just copy
if(ens_handle%in_core) then
   if(indx > 0) then
      ens_handle%ens(indx, :) = member
   else if(indx == 0) then
      ens_handle%mean = member
   else
      ens_handle%spread = member
   endif
else
   ! For out of core, open file and write out new values
   this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, indx)
   ! Open the file and write the field
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
  
   !open(unit = seq_unit, file = this_file_name, access = 'direct', &
   !         form = 'unformatted', recl = req_rec_length)

   write(seq_unit) member
   close(unit = seq_unit)
end if

ens_handle%time(indx) = mtime

end subroutine put_ensemble_member

!-----------------------------------------------------------------

subroutine get_ensemble_time(ens_handle, indx, mtime)

! Returns the time of an ensemble member; 0 for ensemble mean, -1 for spread

type(ensemble_type), intent(in) :: ens_handle
integer, intent(in) :: indx
type(time_type), intent(out) :: mtime

if(indx < -1 .or. indx > ens_handle%ens_size) then
   write(errstring, *) 'indx out of range ', indx
   call error_handler(E_ERR,'get_ensemble_time', errstring, source, revision, revdate)
endif

mtime = ens_handle%time(indx)

end subroutine get_ensemble_time

!-----------------------------------------------------------------

subroutine update_ens_mean(ens_handle)

! Update the value of the ensemble mean from the current ensemble members

type(ensemble_type), intent(inout) :: ens_handle

integer :: i, j, seq_unit
character(len = 129) :: this_file_name
! Sequential version read efficiency requires having second model_size storage
real(r8) :: ens_member(ens_handle%model_size), ens_mean(ens_handle%model_size)

! Update the ensemble mean (copy 0) from ensemble members
if(ens_handle%in_core) then
   do i = 1, ens_handle%model_size
      ens_handle%mean(i) = sum(ens_handle%ens(1:ens_handle%ens_size, i)) / ens_handle%ens_size
   end do
else
   ! Stored on disk, read in each member and sum
   ens_mean = 0.0
   do j = 1, ens_handle%ens_size
      this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, j)
      ! Open the file and read in the field      
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
      !open(unit = seq_unit, file = this_file_name, access = 'direct', &
      !      form = 'unformatted', recl = req_rec_length)

      read(seq_unit) ens_member
      close(unit = seq_unit)
      ! Add this member to sum
      ens_mean = ens_mean + ens_member
   end do
   ens_mean = ens_mean / ens_handle%ens_size
   ! Put the ensemble mean into the file; ensemble mean is 0 member
   this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, 0)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

   !open(unit = seq_unit, file = this_file_name, access = 'direct', &
   !         form = 'unformatted', recl = req_rec_length)

   write(seq_unit) ens_mean
   close(unit = seq_unit)
endif

! Update the ensemble mean time, assume everything else stays in synch
ens_handle%time(0) = ens_handle%time(1)

end subroutine update_ens_mean

!-----------------------------------------------------------------

subroutine update_ens_mean_spread(ens_handle)

! Update the ensemble spread from the current ensemble members,
! Also does an update of the ensemble mean
type(ensemble_type), intent(inout) :: ens_handle

integer :: i, j, seq_unit
character(len = 129) this_file_name
! For efficient disk use now need three copies, if this is too much can alter
real(r8) :: tens_mean(ens_handle%model_size), tens_spread(ens_handle%model_size), &
   ens(ens_handle%model_size)

! Update both the mean and spread (copy -1) from ensemble members
call update_ens_mean(ens_handle)

if(ens_handle%in_core) then
   do i = 1, ens_handle%model_size
      ens_handle%spread(i) = sqrt(sum((ens_handle%ens(1:, i) - ens_handle%mean(i))**2) / &
         (ens_handle%ens_size - 1))
      ! Substitute output of variance IMPORTANT FOR LOW ORDER DISCUSSION
      !!!ens_handle%spread(i) = (sum((ens_handle%ens(1:, i) - ens_handle%mean(i))**2) / &
      !!!   (ens_handle%ens_size - 1))
   end do
else
   ! Get the ensemble mean from disk file; mean is copy 0
   this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, 0)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

   !open(unit = seq_unit, file = this_file_name, access = 'direct', &
   !         form = 'unformatted', recl = req_rec_length)

   read(seq_unit) tens_mean
   close(unit = seq_unit)

   tens_spread = 0.0
   ! Loop through and add squared difference from mean
   do j = 1, ens_handle%ens_size
      this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, j)
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

      !open(unit = seq_unit, file = this_file_name, access = 'direct', &
      !      form = 'unformatted', recl = req_rec_length)

      read(seq_unit) ens
      close(unit = seq_unit)
      tens_spread = tens_spread + (ens - tens_mean)**2
   end do
   
   tens_spread = sqrt(tens_spread / (ens_handle%ens_size - 1))
   ! Write to the ensemble_spread file, index is -1
   this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, -1)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

   !open(unit = seq_unit, file = this_file_name, access = 'direct', &
   !         form = 'unformatted', recl = req_rec_length)

   write(seq_unit) tens_spread
   close(unit = seq_unit)
endif

! Update the spread time
ens_handle%time(-1) = ens_handle%time(1)

end subroutine update_ens_mean_spread

!-----------------------------------------------------------------

subroutine get_ensemble_region(ens_handle, region, rtime, state_vars_in, ens_members_in)

type(ensemble_type), intent(inout) :: ens_handle
real(r8), intent(out) :: region(:, :)
type(time_type), intent(out) :: rtime(:)
integer, intent(in), optional :: state_vars_in(:), ens_members_in(:)

real(r8) :: ens(ens_handle%model_size)
integer :: state_vars(ens_handle%model_size), ens_members(ens_handle%ens_size + 2), i, j
integer :: num_state_vars, num_ens_members, seq_unit
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
   num_state_vars = ens_handle%model_size
   do i = 1, ens_handle%model_size
      state_vars(i) = i
   end do
endif

! Determine size and elements of ensemble subset
if(present(ens_members_in)) then
   num_ens_members = size(ens_members_in)
   ens_members(1:num_ens_members) = ens_members_in
else
   num_ens_members = ens_handle%ens_size
   do i = 1, ens_handle%ens_size
      ens_members(i) = i
   end do
endif

do i = 1, num_ens_members
   rtime(i) = ens_handle%time(ens_members(i))

   ! If not in core, open the file for this ensemble member
   if(.not. ens_handle%in_core) then
      this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, ens_members(i))
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

      !open(unit = seq_unit, file = this_file_name, access = 'direct', &
      !      form = 'unformatted', recl = req_rec_length)

      read(seq_unit) ens
      close(seq_unit)
   endif

   do j = 1, num_state_vars
      if(ens_handle%in_core) then
         if(ens_members(i) == 0) then
            region(i, j) = ens_handle%mean(state_vars(j))
         else if(ens_members(i) == -1) then
            region(i, j) = ens_handle%spread(state_vars(j))
         else
            region(i, j) = ens_handle%ens(ens_members(i), state_vars(j))
         endif
      else
         ! Already have proper ensemble member read in for out of core
         region(i, j) = ens(state_vars(j))
      endif
   end do
end do

end subroutine get_ensemble_region

!-----------------------------------------------------------------

subroutine put_ensemble_region(ens_handle, region, rtime, state_vars_in, ens_members_in)

type(ensemble_type), intent(inout) :: ens_handle
real(r8), intent(in) :: region(:, :)
type(time_type), intent(in) :: rtime(:)
integer, intent(in), optional :: state_vars_in(:), ens_members_in(:)

real(r8) :: ens(ens_handle%model_size)
integer :: state_vars(ens_handle%model_size), ens_members(ens_handle%ens_size + 2), i, j
integer :: num_state_vars, num_ens_members, seq_unit
character(len = 129) :: this_file_name
! Puts a subset of state variables and ensemble members

! Determine size and elements of state variable subset
if(present(state_vars_in)) then
   num_state_vars = size(state_vars_in)
   state_vars(1:num_state_vars) = state_vars_in
else
   num_state_vars = ens_handle%model_size
   do i = 1, ens_handle%model_size
      state_vars(i) = i
   end do
endif

! Determine size and elements of ensemble subset
if(present(ens_members_in)) then
   num_ens_members = size(ens_members_in)
   ens_members(1:num_ens_members) = ens_members_in
else
   num_ens_members = ens_handle%ens_size
   do i = 1, ens_handle%ens_size
      ens_members(i) = i
   end do
endif

do i = 1, num_ens_members
   ens_handle%time(ens_members(i)) = rtime(i)

   ! If not in core, open the file for this ensemble member
   if(.not. ens_handle%in_core) then
      this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, ens_members(i))
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

      !open(unit = seq_unit, file = this_file_name, access = 'direct', &
      !      form = 'unformatted', recl = req_rec_length)

      read(seq_unit) ens
      close(seq_unit)
   endif

   do j = 1, num_state_vars
      if(ens_handle%in_core) then 
         if(ens_members(i) == 0) then
            ens_handle%mean(state_vars(j)) = region(i, j)
         else if(ens_members(i) == -1) then
            ens_handle%spread(state_vars(j)) = region(i, j)
         else
            ens_handle%ens(ens_members(i), state_vars(j)) = region(i, j)
         endif
      else
         ens(state_vars(j)) = region(i, j)
      endif
   end do

   ! If not in core, open the file for this ensemble member
   if(.not. ens_handle%in_core) then
      this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, ens_members(i))
      seq_unit = get_unit()
      open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

      !open(unit = seq_unit, file = this_file_name, access = 'direct', &
      !      form = 'unformatted', recl = req_rec_length)

      write(seq_unit) ens
      close(seq_unit)
   endif
end do

end subroutine put_ensemble_region

!-----------------------------------------------------------------

subroutine end_ensemble_manager(ens_handle, file_name)

type(ensemble_type), intent(inout) :: ens_handle
character(len = 129), intent(in), optional :: file_name

real(r8) :: ens(ens_handle%model_size)
integer :: iunit, i, seq_unit
character(len = 129) :: this_file_name, command_string
character(len = 4) :: extension

if(present(file_name)) then

   if(single_restart_file_out) iunit = open_restart_write(file_name)

   do i = 1, ens_handle%ens_size
   
      if(.not. single_restart_file_out) then
         write(extension, 99) i
99       format(i4.4)
         this_file_name = trim(file_name) // '.' // extension
         iunit = open_restart_write(this_file_name)
      endif

      if(ens_handle%in_core) then
         call awrite_state_restart(ens_handle%time(i), ens_handle%ens(i, :), iunit)
      else
         this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, i)
         seq_unit = get_unit()
         open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')

         read(seq_unit) ens
         close(seq_unit)
         call awrite_state_restart(ens_handle%time(i), ens, iunit)
      endif

      if(.not. single_restart_file_out) call close_restart(iunit)

   end do

   if(single_restart_file_out) call close_restart(iunit)

endif

! Free up the allocated storage
if(ens_handle%in_core) then
   deallocate(ens_handle%ens, ens_handle%mean, ens_handle%spread, ens_handle%time)
else
   deallocate(ens_handle%time)
endif

! First extension is ens_id
write(extension, 99) ens_handle%id

if(.not. ens_handle%in_core) then
   command_string = 'rm -f ' // trim(ens_file_name) // '.' // extension //'.*'
   call system(command_string)
   command_string = 'rm -f ' // trim(reg_file_name) // '.' // extension // '.*'
   call system(command_string)
endif

end subroutine end_ensemble_manager

!-----------------------------------------------------------------






subroutine Aadvance_state(ens_handle, target_time, asynch, adv_ens_command)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.

implicit none

type(ensemble_type),intent(inout) :: ens_handle
type(time_type),    intent(in)    :: target_time
integer,            intent(in)    :: asynch
character(len=129), intent(in)    :: adv_ens_command

type(time_type) :: time_step

integer :: seconds, days, i, control_unit, ic_file_unit, ud_file_unit

character(len = 26), dimension(ens_handle%ens_size) :: ic_file_name, ud_file_name 
character(len = 128) :: input_string
integer :: is1,is2,id1,id2
type(time_type) :: smodel_time
real(r8) :: smodel_state(ens_handle%model_size)

! Loop through each model state and advance
do i = 1, ens_handle%ens_size

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
         write(errstring,*)'Trying to use ',ens_handle%ens_size,' model states -- too many.'
         call error_handler(E_MSG,'Aadvance_state',errstring,source,revision,revdate)
         call error_handler(E_ERR,'Aadvance_state','Use less than 10000 model states.',source,revision,revdate)
      endif

 11   format(a20, i1)
 21   format(a20, i2)
 31   format(a20, i3)
 41   format(a20, i4)
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

! Following block needed for WRF only
! Write out the ensemble mean for the calculation of flow-dependent BC of WRF.
!call update_ens_mean(ens_handle)
!ic_file_unit = open_restart_write('assim_model_state_ic_mean')
!call awrite_state_restart(smodel_time, ens_mean, ic_file_unit, target_time)
!call close_restart(ic_file_unit)

! Also need synchronization block at the end for the asynch

if(asynch /= 0) then

   ! Write out the file names to a control file

   control_unit = get_unit()
   open(unit = control_unit, file = 'filter_control')
   write(control_unit, *) ens_handle%ens_size
   do i = 1, ens_handle%ens_size
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
   do i = 1, ens_handle%ens_size
      ud_file_unit = open_restart_read(ud_file_name(i))
      call aread_state_restart(smodel_time, smodel_state, ud_file_unit)
      call close_restart(ud_file_unit)
      ! Put the updated state in the ensemble storage
      call put_ensemble_member(ens_handle, i, smodel_state, smodel_time)
   end do

end if

end subroutine Aadvance_state

!-----------------------------------------------------------------

subroutine transpose_ens_to_regions(ens_handle, num_regions, region_id, region_size)

type(ensemble_type), intent(inout) :: ens_handle
integer, intent(in) :: num_regions
integer, intent(in) :: region_id(ens_handle%model_size), region_size(num_regions)

! Temporary storage to reorder into regions
real(r8) :: reorder(ens_handle%model_size), ens(ens_handle%model_size)
integer :: start(num_regions), indx(num_regions), i, j, this_region
integer :: reg_unit, seq_unit, req_rec_length
character(len = 129) :: this_file_name


! Don't do anything if things are in core
if(ens_handle%in_core) return

! Compute where in reordered storage each region should start
start(1) = 1
do i = 2, num_regions
   start(i) = start(i - 1) + region_size(i - 1)
end do

! Loop through each of the ensembles and extract the regions
! Then write the regions to region files
do i = 1, ens_handle%ens_size

   ! Initialize where each region should write
   indx = start

   ! Open the ensemble file
   this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, i)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   ! Read in this ensemble member
   read(seq_unit) ens 
   close(seq_unit)

   ! Can now delete the ensemble file for storage efficiency
   ! SYSTEM CALLS ARE VERY EXPENSIVE!!!
   !!!call system('rm -f ' // trim(this_file_name))

   ! Now loop through and extract the regions
   do j = 1, ens_handle%model_size
      this_region = region_id(j)
      reorder(indx(this_region)) = ens(j)
      indx(this_region) = indx(this_region) + 1
   end do

   ! Now the whole beasts are there, write them to the region files
   reg_unit = get_unit()
   do j = 1, num_regions

      ! The regions file must be random access for reverse transpose
      ! Open a file to this regions file
      this_file_name = get_disk_file_name(reg_file_name, ens_handle%id, j)
      ! Get the length of the record (one ensemble worth)
      inquire(iolength = req_rec_length) reorder(start(j):start(j) + region_size(j) - 1)
      open(unit = reg_unit, file = this_file_name, access = 'direct', &
         form = 'unformatted', recl = req_rec_length)
      ! Write state for this region, this ensemble member
      write(reg_unit, rec = i) reorder(start(j):start(j) + region_size(j) - 1)
      close(reg_unit)
   end do
end do

end subroutine transpose_ens_to_regions


!-----------------------------------------------------------------

subroutine transpose_regions_to_ens(ens_handle, num_regions, region_id, region_size)

type(ensemble_type), intent(inout) :: ens_handle
integer, intent(in) :: num_regions
integer, intent(in) :: region_id(ens_handle%model_size), region_size(num_regions)

! Temporary storage to reorder into regions
real(r8) :: reorder(ens_handle%model_size), ens(ens_handle%model_size)
integer :: start(num_regions), indx(num_regions), i, j, this_region
integer :: reg_unit, seq_unit, req_rec_length
character(len = 129) :: this_file_name

! Don't do anything if things are in core
if(ens_handle%in_core) return

! Compute where in reordered storage each region should start
start(1) = 1
do i = 2, num_regions
   start(i) = start(i - 1) + region_size(i - 1)
end do

! Loop through to read the ensembles from each region in turn
do i = 1, ens_handle%ens_size
   ! Loop through each regions file
   do j = 1, num_regions
      this_file_name = get_disk_file_name(reg_file_name, ens_handle%id, j)
      reg_unit = get_unit()
      ! Get the record length for this region (one ensemble member)
      inquire(iolength = req_rec_length) reorder(start(j):start(j) + region_size(j) -1)
      open(unit = reg_unit, file = this_file_name, access = 'direct', &
         form = 'unformatted', recl = req_rec_length)
      ! Read the region into temporary storage
      read(reg_unit, rec = i) reorder(start(j):start(j) + region_size(j) - 1)
      close(reg_unit)
      ! Can now delete the region file if this is last ensemble
      ! SYSTEM CALLS ARE VERY EXPENSIVE!!!
      !!!if(i == ens_size) call system('rm -f ' // this_file_name)
   end do

   ! All regions are in, move them to ensemble order
   ! Initialize where each region should write
   indx = start
   do j = 1, ens_handle%model_size
      this_region = region_id(j)
      ens(j) = reorder(indx(this_region))
      indx(this_region) = indx(this_region) + 1
   end do

   ! Now write this ensemble member to disk
   this_file_name = get_disk_file_name(ens_file_name, ens_handle%id, i)
   seq_unit = get_unit()
   open(unit = seq_unit, file = this_file_name, access = 'sequential', form = 'unformatted')
   write(seq_unit) ens
   close(seq_unit)
end do

end subroutine transpose_regions_to_ens

!-----------------------------------------------------------------

subroutine get_region_by_number(ens_handle, rnum, rsize, region, region_id)

! Returns all ensemble members for region number rnum. When the ensemble is
! being stored in-core, the region_id is used to find all the state variables
! in this region and they are copied into region. When the ensemble is NOT
! being stored in core, it is ASSUMED that the ensemble has been transposed to
! regional files and that the size and ids of the regions for the transpose
! are consistent with the region size and region_id in this call. This could
! all be made safer, but with some overhead.

type(ensemble_type), intent(inout) :: ens_handle
integer, intent(in) :: rnum, rsize
real(r8), intent(out) :: region(ens_handle%ens_size, rsize)
integer, intent(in) :: region_id(rsize)

integer :: i, j, indx, reg_unit, req_rec_length
character(len = 129) :: this_file_name

! For in-core, just do copying
if(ens_handle%in_core) then
   do i = 1, rsize
      indx = region_id(i)
      do j = 1, ens_handle%ens_size
         region(j, i) = ens_handle%ens(j, indx)
      end do
   end do

else
   !For out of core, need to open and read the appropriate file
   this_file_name = get_disk_file_name(reg_file_name, ens_handle%id, rnum)
   reg_unit = get_unit()
   inquire(iolength = req_rec_length) region(1, :)
   open(unit = reg_unit, file = this_file_name, access = 'direct', &
      form = 'unformatted', recl = req_rec_length)
   do j = 1, ens_handle%ens_size
      read(reg_unit, rec = j) region(j, :)
   end do
   close(reg_unit)
   
endif

end subroutine get_region_by_number

!-----------------------------------------------------------------

subroutine put_region_by_number(ens_handle, rnum, rsize, region, region_id)

type(ensemble_type), intent(inout) :: ens_handle
integer, intent(in) :: rnum, rsize
real(r8), intent(in) :: region(ens_handle%ens_size, rsize)
integer, intent(in) :: region_id(rsize)

integer :: i, j, indx, reg_unit, req_rec_length
character(len = 129) :: this_file_name

! For in-core, just do copying
if(ens_handle%in_core) then
   do i = 1, rsize
      indx = region_id(i)
      do j = 1, ens_handle%ens_size
         ens_handle%ens(j, indx) = region(j, i)
      end do
   end do

else
   !For out of core, need to open and read the appropriate file
   this_file_name = get_disk_file_name(reg_file_name, ens_handle%id, rnum)
   reg_unit = get_unit()
   inquire(iolength = req_rec_length) region(1, :)
   open(unit = reg_unit, file = this_file_name, access = 'direct', &
      form = 'unformatted', recl = req_rec_length)
   do j = 1, ens_handle%ens_size
      write(reg_unit, rec = j) region(j, :)
   end do
   close(reg_unit)
   
endif

end subroutine put_region_by_number

!-----------------------------------------------------------------

function get_disk_file_name(base_file, ens_id, indx)

character(len = 129)  :: get_disk_file_name
character(len = *), intent(in) :: base_file
integer, intent(in) :: ens_id, indx

character(len = 4) extension1, extension2

! First extension is ens_id
write(extension1, 99) ens_id
99 format(i4.4)

if(indx > 0) then
   ! Open the direct access file for this ensemble member
   write(extension2, 99) indx
   get_disk_file_name = trim(base_file) // '.' // extension1 // '.' // extension2
else if(indx == 0) then
   get_disk_file_name = trim(base_file) // '.' // extension1 // '.mean'
else if(indx == -1) then
   get_disk_file_name = trim(base_file) // '.' // extension1 // '.spread'
endif

end function get_disk_file_name

!-----------------------------------------------------------------

function is_ens_in_core(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

! Possible concerns about direct access to a namelist variable,
! Use a function to return whether the ensemble is stored in core or not.

logical :: is_ens_in_core

is_ens_in_core = ens_handle%in_core

end function is_ens_in_core


!-----------------------------------------------------------------

subroutine copy_ens(ens1, ens2, copy_time)

type(ensemble_type), intent(in)    :: ens1
type(ensemble_type), intent(inout) :: ens2
logical, intent(in)              :: copy_time

character(len = 129) :: errstring
real(r8)             :: ens(ens1%model_size)
type(time_type)      :: ens_time
integer              :: i

! Copies all data from ens1 to ens2. If copy_time is true, also 
! copies the time information from ens1 to ens2. If copy_time is
! false, the times in ens2 are left unchanged.
! WARNING: DOES NOT COPY THE ENSEMBLE MEAN AND SPREAD.

! Check to make sure that the ensembles are compatible
if(ens1%ens_size /= ens2%ens_size .or. ens1%model_size /= ens2%model_size)then
   write(errstring, *) 'ensemble size and/or model size incompatible for copy', &
      ens1%ens_size, ens2%ens_size, ens1%model_size, ens2%model_size
   call error_handler(E_ERR,'copy_ens', errstring, source, revision, revdate)
endif

! Loop through each ensemble member, get, then put
do i = 1, ens1%ens_size
   call get_ensemble_member(ens1, i, ens, ens_time)

   ! Use the current time from ens2 if copy_time is false
   if(.not. copy_time) call get_ensemble_time(ens2, i, ens_time)

   call put_ensemble_member(ens2, i, ens, ens_time)
end do

end subroutine copy_ens

!-----------------------------------------------------------------


end module ensemble_manager_mod
