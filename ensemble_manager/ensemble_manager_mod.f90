! Data Assimilation Research Testbed -- DART

   ! Open the direct access file for on disk ensemble
   
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module ensemble_manager_mod

use types_mod, only : r8
use    utilities_mod, only :  get_unit, open_file, close_file, register_module, &
                              check_nml_error, file_exist, error_handler, &
                              E_ERR, E_WARN, E_MSG, E_DBG, initialize_utilities, &
                              logfileunit, timestamp
use time_manager_mod, only : time_type
use assim_model_mod, only : aread_state_restart, awrite_state_restart, open_restart_read, &
   open_restart_write, close_restart, assim_model_type, get_model_time_step, adv_1step
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
   get_ensemble_region, put_ensemble_region, get_ensemble_time, Aadvance_state
   
! Global in core storage for the ensemble
real(r8), allocatable :: ens(:, :)
type(time_type), allocatable :: ens_time(:)
integer :: ens_size, model_size, req_rec_length

integer :: direct_unit

!-----------------------------------------------------------------
!
! namelist with default values

logical :: in_core = .true.
! If true, then ensemble is held in core storage; if false, only
! a single ensemble member worth of storage is used and ensemble
! is kept in a direct access temporary file.

namelist / ensemble_manager_nml / in_core

!-----------------------------------------------------------------

contains

!-----------------------------------------------------------------

subroutine init_ensemble_manager(ens_size_in, model_size_in, file_name, &
   init_time)

integer, intent(in) :: ens_size_in, model_size_in
character(len = 129), intent(out), optional :: file_name
type(time_type), intent(in), optional :: init_time

real(r8) :: x
integer :: iunit, i, j, ilength, rlength, rnum, ierr, io
character(len = 129) :: msgstring

call register_module(source, revision, revdate)

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

if(.not. in_core) then
   ! Find out the required record length for real(r8)
   inquire(iolength = req_rec_length) x

   ! Open the direct access file for on disk ensemble
   direct_unit = get_unit()
   open(unit = direct_unit, file = "temp_ens_file", access = 'direct', &
      status = 'replace', form = 'unformatted', recl = req_rec_length)
   
   ! Write out the first two 0 filled ensemble members for ensemble spread and mean
   do i = 1, 2
      do j = 1, model_size
         rnum = (i - 1) * model_size + j
         write(direct_unit, rec = rnum) -99.99_r8
      end do
   end do
endif

! Set the global storage bounds
ens_size = ens_size_in
model_size = model_size_in

! Initialize the storage and read in from restart file if needed
if(in_core) then 
   allocate(ens(-1:ens_size, model_size), ens_time(-1:ens_size))
else
   allocate(ens(1, model_size), ens_time(-1:ens_size))
endif

if(present(file_name)) then
   iunit = open_restart_read(file_name)

   do i = 1, ens_size
      write(msgstring, *) 'trying to read restart ', i
      call error_handler(E_DBG,'init_ensemble_manager',msgstring,source,revision,revdate)
      if(in_core) then 
         call aread_state_restart(ens_time(i), ens(i, :), iunit)
      else
         call aread_state_restart(ens_time(i), ens(1, :), iunit)
         ! Write out the ensemble to direct access file
         do j = 1, model_size
            rnum = (i + 1) * model_size + j
            write(direct_unit, rec = rnum) ens(j, 1) 
         end do
      endif

      ! If the inital file time is being overridden, do it
      if(present(init_time)) ens_time(i) = init_time
   end do
   call close_restart(iunit)
endif


end subroutine init_ensemble_manager

!-----------------------------------------------------------------

subroutine get_ensemble_member(index, member, mtime)

integer, intent(in) :: index
real(r8), intent(out) :: member(:)
type(time_type), intent(out) :: mtime

integer :: i, rnum

if(index < -1 .or. index > ens_size) then
   write(*, *) 'index out of range in get_ensemble_member'
   stop
endif

if(in_core) then
   member = ens(index, :)
else
   do i = 1, model_size
      rnum = (index + 1) * model_size + i
      read(direct_unit, rec = rnum) member(i)
   end do
endif

mtime = ens_time(index)

end subroutine get_ensemble_member

!-----------------------------------------------------------------

subroutine put_ensemble_member(index, member, mtime)

integer, intent(in) :: index
real(r8), intent(in) :: member(:)
type(time_type), intent(in) :: mtime

integer :: i, rnum

if(index < -1 .or. index > ens_size) then
   write(*, *) 'index out of range in put_ensemble_member'
   write(*, *) 'index is ', index
   stop
endif

if(in_core) then
   ens(index, :) = member
else
   do i = 1, model_size
      rnum = (index + 1) * model_size + i
      write(direct_unit, rec = rnum) member(i)
   end do
end if

ens_time(index) = mtime

end subroutine put_ensemble_member

!-----------------------------------------------------------------

subroutine get_ensemble_time(index, mtime)

integer, intent(in) :: index
type(time_type), intent(out) :: mtime

if(index < -1 .or. index > ens_size) then
   write(*, *) 'index out of range in get_ensemble_member'
   stop
endif

mtime = ens_time(index)

end subroutine get_ensemble_time

!-----------------------------------------------------------------

subroutine update_ens_mean()

integer :: i, j, rnum
real(r8) :: val

! Update the ensemble mean (copy 0) from ensemble members
if(in_core) then
    do i = 1, model_size
      ens(0, i) = sum(ens(1:ens_size, i)) / ens_size
   end do
else
   ens = 0.0
   do j = 1, ens_size
      do i = 1, model_size
         rnum = (j + 1) * model_size + i
         read(direct_unit, rec = rnum) val
         ens(1, i) = ens(1, i) + val
      end do
   end do
   ens = ens / ens_size
   ! Put the ensemble mean into the file
   do i = 1, model_size
      rnum = model_size + i
      write(direct_unit, rec = rnum) ens(1, i)
   end do
endif

! Update the ensemble mean time, assume everything else stays in synch
ens_time(0) = ens_time(1)

end subroutine update_ens_mean

!-----------------------------------------------------------------

subroutine update_ens_mean_spread()

integer :: i, rnum, j
real(r8) :: temp, ens_spread, ens_mean, ens_element

! Update both the mean and spread (copy -1) from ensemble members
call update_ens_mean()

if(in_core) then
   do i = 1, model_size
      ens(-1, i) = sqrt(sum((ens(1:, i) - ens(0, i))**2) / (ens_size - 1))
   end do
else
   ! Loop through and add squared difference from mean
   ! This is going to wreak havoc on disks potentially
   do i = 1, model_size
      ens_spread = 0.0
      rnum = model_size + i
      read(direct_unit, rec = rnum) ens_mean
      do j = 1, ens_size
         rnum = (j + 1) * model_size + i
         read(direct_unit, rec = rnum) ens_element
         ens_spread = ens_spread + (ens_element - ens_mean)**2
      end do
      write(direct_unit, rec = i) sqrt(ens_spread / (ens_size - 1))
   end do
endif

! Update the spread time
ens_time(-1) = ens_time(1)

end subroutine update_ens_mean_spread

!-----------------------------------------------------------------

subroutine get_ensemble_region(region, rtime, state_vars_in, ens_members_in)

real(r8), intent(out) :: region(:, :)
type(time_type), intent(out) :: rtime(:)
integer, intent(in), optional :: state_vars_in(:), ens_members_in(:)

integer :: state_vars(model_size), ens_members(-1:ens_size), i, j, rnum
integer :: num_state_vars, num_ens_members
! Returns a subset of the state variables and ensemble members

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
   do j = 1, num_state_vars
      if(in_core) then
         region(i, j) = ens(ens_members(i), state_vars(j))
      else
         rnum = (ens_members(i) + 1) * model_size + state_vars(j)
         read(direct_unit, rec = rnum) region(i, j)
      endif
   end do
end do

end subroutine get_ensemble_region

!-----------------------------------------------------------------

subroutine put_ensemble_region(region, rtime, state_vars_in, ens_members_in)

real(r8), intent(in) :: region(:, :)
type(time_type), intent(in) :: rtime(:)
integer, intent(in), optional :: state_vars_in(:), ens_members_in(:)

integer :: state_vars(model_size), ens_members(-1:ens_size), i, j, rnum
integer :: num_state_vars, num_ens_members
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
   do j = 1, num_state_vars
      if(in_core) then 
         ens(ens_members(i), state_vars(j)) = region(i, j)
      else
         rnum = (ens_members(i) + 1) * model_size + state_vars(j)
         write(direct_unit, rec = rnum) region(i, j)
      endif
      ens_time(ens_members(i)) = rtime(i)
   end do
end do

end subroutine put_ensemble_region

!-----------------------------------------------------------------

subroutine end_ensemble_manager(file_name)

character(len = 129), intent(in), optional :: file_name

integer :: iunit, i, j, rnum

iunit = open_restart_write(file_name)

do i = 1, ens_size
   if(in_core) then
      call awrite_state_restart(ens_time(i), ens(i, :), iunit)
   else
      do j = 1, model_size
         rnum = (i + 1) * ens_size + j
         read(direct_unit, rec = rnum) ens(1, j)
      end do
      call awrite_state_restart(ens_time(i), ens(1, :), iunit)
   endif
end do
call close_restart(iunit)
close(direct_unit)

! Free up the allocated storage
deallocate(ens, ens_time)

end subroutine end_ensemble_manager

!-----------------------------------------------------------------






subroutine Aadvance_state(ens_size, model_size, target_time, asynch, adv_ens_command)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.

implicit none

integer,            intent(in)    :: ens_size, model_size
type(time_type),    intent(in)    :: target_time
integer,            intent(in)    :: asynch
character(len=129), intent(in)    :: adv_ens_command                                                      

type(time_type) :: time_step

integer :: seconds, days, i, control_unit, ic_file_unit, ud_file_unit

character(len = 26), dimension(ens_size) :: ic_file_name, ud_file_name 
character(len = 128) :: input_string
character(len = 129) :: errstring
integer :: is1,is2,id1,id2
type(time_type) :: smodel_time
real(r8) :: smodel_state(model_size)

! Loop through each model state and advance
do i = 1, ens_size

   ! Get the state and time for this ensemble member
   call get_ensemble_member(i, smodel_state, smodel_time)

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
      call put_ensemble_member(i, smodel_state, smodel_time)

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

   else

      write(errstring,*)'input.nml - async is ',asynch,' must be 0, 1, or 2' 
      call error_handler(E_ERR,'Aadvance_state', errstring, source, revision, revdate)

   endif

   write(*, *) 'got clearance to proceed in Aadvance_state'

   ! All should be done, read in the states and proceed
   do i = 1, ens_size
      ud_file_unit = open_restart_read(ud_file_name(i))
      call aread_state_restart(smodel_time, smodel_state, ud_file_unit)
      call close_restart(ud_file_unit)
      ! Put the updated state in the ensemble storage
      call put_ensemble_member(i, smodel_state, smodel_time)
   end do

end if

end subroutine Aadvance_state


end module ensemble_manager_mod

