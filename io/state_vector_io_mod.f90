! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Routines for reading the model state.
!> The limited transpose routines are in direct_netcdf_mod.
module state_vector_io_mod

use direct_netcdf_mod,    only : read_transpose, transpose_write

use types_mod,            only : r8, MISSING_R8, i4

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, &
                                 is_single_restart_file_in, map_task_to_pe, &
                                 get_copy_owner_index, put_copy
                                 

use utilities_mod,        only : error_handler, nc_check, check_namelist_read, &
                                 find_namelist_in_file, nmlfileunit,           &
                                 do_nml_file, do_nml_term,          &
                                 E_MSG, E_ERR, get_unit, ascii_file_format, &
                                 close_file, dump_unit_attributes, &
                                 register_module, set_output

use assim_model_mod,      only : get_model_size, pert_model_state, &
                                 assim_model_type

use time_manager_mod,     only : time_type, read_time, write_time, &
                                 get_time

use io_filenames_mod,     only : get_input_file


use random_seq_mod,       only : random_seq_type, init_random_seq, random_gaussian

!> @todo  This should go through assim_model_mod
use model_mod,            only : read_model_time

use copies_on_off_mod,    only : setup_read_write, turn_read_copy_on,      &
                                 turn_write_copy_on, turn_read_copies_off, &
                                 turn_write_copy_off

use state_structure_mod,  only : get_num_domains

use netcdf


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

private

public :: state_vector_io_init, &
          read_transpose,       &
          transpose_write,      &
          setup_read_write,     &
          turn_read_copy_on,    &
          turn_write_copy_on,   &
          turn_read_copies_off, &
          turn_write_copy_off, &
          aread_state_restart, &
          awrite_state_restart, &
          write_ensemble_restart, &
          read_ensemble_restart, &
          open_restart_read, &
          open_restart_write, &
          close_restart, &
          filter_read_restart_direct, &
          filter_write_restart_direct

integer :: ret !< netcdf return code
integer :: ncfile !< netcdf input file identifier
integer :: ncfile_out !< netcdf output file handle

! Module storage for writing error messages
character(len = 255) :: msgstring

! Logical flag for initialization of module
logical              :: module_initialized = .false.

! Global storage for default restart formats
character(len = 16) :: read_format = "unformatted", write_format = "unformatted"

! namelist variables with default values
! Aim: to have the regular transpose as the default
integer :: limit_mem = HUGE(1_i4)!< This is the number of elements (not bytes) so you don't have times the number by 4 or 8
logical :: time_unlimited = .true. ! You need to keep track of the time.
logical :: single_precision_output = .false. ! Allows you to write r4 netcdf files even if filter is double precision

logical  :: single_restart_file_out = .true.
! Size of perturbations for creating ensembles when model won't do it
real(r8) :: perturbation_amplitude  = 0.2_r8
! write_binary_restart_files  == .true.  -> use unformatted file format. 
!                                     Full precision, faster, smaller,
!                                     but not as portable.
logical  :: write_binary_restart_files = .false.


namelist /  state_vector_io_nml / limit_mem, time_unlimited, &
   single_precision_output, single_restart_file_out, perturbation_amplitude, &
   write_binary_restart_files

contains

!-------------------------------------------------
!> Initialize model 
!> so you can read the namelist
subroutine state_vector_io_init()

integer :: iunit, io

if ( .not. module_initialized ) then
   ! Initialize the module with utilities 
   call register_module(source, revision, revdate)
   module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "state_vector_io_nml", iunit)
read(iunit, nml = state_vector_io_nml, iostat = io)
call check_namelist_read(iunit, io, "state_vector_io_nml")

   ! Set the write format for restart files
   if(write_binary_restart_files) then
      write_format = "unformatted"
   else
      write_format = "formatted"
   endif

endif


! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=state_vector_io_nml)
if (do_nml_term()) write(     *     , nml=state_vector_io_nml)

end subroutine state_vector_io_init

!-------------------------------------------------
! DART read write (binary/ascii vector)
!-------------------------------------------------


subroutine write_state_restart(assim_model, funit, target_time)
!----------------------------------------------------------------------
!
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

implicit none

type (assim_model_type), intent(in)           :: assim_model
integer,                 intent(in)           :: funit
type(time_type),         optional, intent(in) :: target_time

if(present(target_time)) then
!   call awrite_state_restart(assim_model%time, assim_model%state_vector, funit, target_time)
else
!   call awrite_state_restart(assim_model%time, assim_model%state_vector, funit)
endif

end subroutine write_state_restart




subroutine awrite_state_restart(model_time, model_state, funit, target_time)
!----------------------------------------------------------------------
!
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

implicit none

type(time_type), intent(in)           :: model_time
real(r8),        intent(in)           :: model_state(:)
integer,         intent(in)           :: funit
type(time_type), optional, intent(in) :: target_time

integer :: i, io, rc
character(len = 16) :: open_format
character(len=128) :: filename
logical :: is_named

if ( .not. module_initialized ) call state_vector_io_init()

! Figure out whether the file is opened FORMATTED or UNFORMATTED
inquire(funit, FORM=open_format)

! assume success
io = 0

! Write the state vector
if (ascii_file_format(open_format)) then
   if(present(target_time)) call write_time(funit, target_time, ios_out=io)
   if (io /= 0) goto 10
   call write_time(funit, model_time, ios_out=io)
   if (io /= 0) goto 10
   do i = 1, size(model_state)
      write(funit, *, iostat = io) model_state(i)
      if (io /= 0) goto 10
   end do
else
   if(present(target_time)) call write_time(funit, target_time, form="unformatted", ios_out=io)
   if (io /= 0) goto 10
   call write_time(funit, model_time, form="unformatted", ios_out=io)
   if (io /= 0) goto 10
   write(funit, iostat = io) model_state
   if (io /= 0) goto 10
endif

! come directly here on error. 
10 continue

! if error, use inquire function to extract filename associated with
! this fortran unit number and use it to give the error message context.
if (io /= 0) then
   inquire(funit, named=is_named, name=filename, iostat=rc)
   if ((rc /= 0) .or. (.not. is_named)) filename = 'unknown file'
   write(msgstring,*) 'error writing to restart file ', trim(filename)
   call error_handler(E_ERR,'awrite_state_restart',msgstring,source,revision,revdate)
endif


end subroutine awrite_state_restart

!----------------------------------------------------------------------
!
! Closes a restart file
subroutine close_restart(file_unit)

integer, intent(in) :: file_unit

call close_file(file_unit)

end subroutine close_restart

!---------------------------------------------------------------------------------

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
if(.not. start_from_restart .and. .not. is_single_restart_file_in()) then
   write(msgstring, *) 'start_from_restart in filter_nml and single_restart_file_in in &
        &ensemble_manager_nml cannot both be false'
   call error_handler(E_ERR,'read_ensemble_restart', msgstring, source, revision, revdate)
endif

! the code reads into the vars array
!ens_handle%valid = VALID_VARS

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
if(is_single_restart_file_in() .or. .not. start_from_restart .or. &
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
      ! map from my_pe to physical task number all done in send and receives only
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

subroutine read_state_restart(assim_model, funit, target_time)
!----------------------------------------------------------------------
!
! Read a restart file given a unit number (see write_state_restart)

implicit none

type(assim_model_type), intent(inout)         :: assim_model
integer,                intent(in)            :: funit
type(time_type),        optional, intent(out) :: target_time

if ( .not. module_initialized ) call state_vector_io_init()

if(present(target_time)) then
   !call aread_state_restart(assim_model%time, assim_model%state_vector, funit, target_time)
else
   !call aread_state_restart(assim_model%time, assim_model%state_vector, funit)
endif

end subroutine read_state_restart

!----------------------------------------------------------------------

subroutine aread_state_restart(model_time, model_state, funit, target_time)
!----------------------------------------------------------------------
!
! Read a restart file given a unit number (see write_state_restart)

implicit none

type(time_type), intent(out)            :: model_time
real(r8),        intent(out)            :: model_state(:)
integer,         intent(in)             :: funit
type(time_type), optional, intent(out) :: target_time

character(len = 16) :: open_format
integer :: ios, int1, int2

if ( .not. module_initialized ) call state_vector_io_init()

ios = 0

! Figure out whether the file is opened FORMATTED or UNFORMATTED
inquire(funit, FORM=open_format)

if (ascii_file_format(open_format)) then
   if(present(target_time)) target_time = read_time(funit)
   model_time = read_time(funit)
   read(funit,*,iostat=ios) model_state
else
   if(present(target_time)) target_time = read_time(funit, form = "unformatted")
   model_time = read_time(funit, form = "unformatted")
   read(funit,iostat=ios) model_state
endif

! If the model_state read fails ... dump diagnostics.
if ( ios /= 0 ) then
   ! messages are being used as error lines below.  in an MPI filter,
   ! all messages are suppressed that aren't from PE0.  if an error
   ! happens in another task, these lines won't be printed unless we
   ! turn on output.
   call set_output(.true.)

   write(msgstring,*)'dimension of model state is ',size(model_state)
   call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)

   if(present(target_time)) then
      call get_time(target_time, int1, int2)       ! time -> secs/days
      write(msgstring,*)'target_time (secs/days) : ',int1,int2
      call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)
   endif

   call get_time(model_time, int1, int2)       ! time -> secs/days
   write(msgstring,*)'model_time (secs/days) : ',int1,int2
   call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)

   write(msgstring,'(''model max/min/first is'',3(1x,E12.6) )') &
            maxval(model_state), minval(model_state), model_state(1)
   call error_handler(E_MSG,'aread_state_restart',msgstring,source,revision,revdate)

   call dump_unit_attributes(funit)

   write(msgstring,*)'read error is : ',ios
   call error_handler(E_ERR,'aread_state_restart',msgstring,source,revision,revdate)
endif

end subroutine aread_state_restart

!----------------------------------------------------------------------


function open_restart_write(file_name, override_write_format)
!----------------------------------------------------------------------
!
! Opens a restart file for writing

character(len = *), intent(in) :: file_name
character(len = *), optional, intent(in) :: override_write_format

integer :: open_restart_write, io

if ( .not. module_initialized ) call state_vector_io_init()

open_restart_write = get_unit()
if(present(override_write_format)) then
   open(unit = open_restart_write, file = file_name, form = override_write_format, &
        iostat = io)
else
   open(unit = open_restart_write, file = file_name, form = write_format, iostat = io)
endif
if (io /= 0) then
   write(msgstring,*) 'unable to open restart file ', trim(file_name), ' for writing'
   call error_handler(E_ERR,'open_restart_write',msgstring,source,revision,revdate)
endif

end function open_restart_write


function open_restart_read(file_name)
!----------------------------------------------------------------------
!
! Opens a restart file for reading

integer :: open_restart_read
character(len = *), intent(in) :: file_name

integer :: ios, ios_out
!!logical :: old_output_state
type(time_type) :: temp_time
character(len=64) :: string2

if ( .not. module_initialized ) call state_vector_io_init()

! DEBUG -- if enabled, every task will print out as it opens the
! restart files.  If questions about missing restart files, first start
! by commenting in only the timestamp line.  If still concerns, then
! go ahead and comment in all the lines.
!!old_output_state = do_output()
!!call set_output(.true.)
!call timestamp("open_restart", "opening restart file "//trim(file_name), pos='')
!!call set_output(old_output_state)
!END DEBUG

! if you want to document which file(s) are being opened before
! trying the open (e.g. in case the fortran runtime library intercepts
! the error and does not return to let us print out the name) then
! comment this in and you can see what files are being opened.
!write(msgstring, *) 'Opening restart file ',trim(adjustl(file_name))
!call error_handler(E_MSG,'open_restart_read',msgstring,source,revision,revdate)

! WARNING: Absoft Pro Fortran 9.0, on a power-pc mac, is convinced
! that certain binary files are, in fact, ascii, because the read_time 
! call is returning what seems like a good time even though it should
! be garbage.  This code works fine on all other platforms/compilers
! we've tried, so we're leaving it as-is.  Best solution if you're
! using absoft on a mac is to set all files to be non-binary in the
! namelist.  You may also have to set the format in both obs_model_mod.f90 
! and interpolate_model.f90 to 'formatted' instead of the hardcoded 
! 'unformatted' for async 2/4 model advance temp_ic and temp_ud files.

! Autodetect format of restart file when opening
! Know that the first thing in here has to be a time, so try to read it.
! If it fails with one format, try the other. If it fails with both, punt.
open_restart_read = get_unit()
read_format = 'formatted'
open(unit   = open_restart_read, &
     file   = trim(file_name),   &
     form   = read_format,       &
     action = 'read',            &
     status = 'old',             &
     iostat = ios)
! An opening error means something is wrong with the file, error and stop
if(ios /= 0) goto 11
temp_time = read_time(open_restart_read, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be formatted, proceed
   rewind open_restart_read
   return
endif

! Next, try to see if an unformatted read works instead
close(open_restart_read)

open_restart_read = get_unit()
read_format = 'unformatted'
open(unit   = open_restart_read, &
     file   = trim(file_name),   &
     form   = read_format,       &
     action = 'read',            &
     status = 'old',             &
     iostat = ios)
! An opening error means something is wrong with the file, error and stop
if(ios /= 0) goto 11
rewind open_restart_read
temp_time = read_time(open_restart_read, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be unformatted, proceed
   rewind open_restart_read
   return
endif

! Otherwise, neither format works. Have a fatal error.
11 continue

write(msgstring, *) 'Problem opening file ',trim(file_name)
write( string2 , *) 'OPEN status was ',ios
call error_handler(E_ERR, 'open_restart_read', msgstring, &
     source, revision, revdate, text2=string2)

end function open_restart_read


!-----------------------------------------------------------------

!-------------------------------------------------
! Netcdf read write
! Uses direct_netcdf_mpi_mod.f90 or direct_netcdf_no_mpi_mod.f90
!-------------------------------------------------

!------------------------------------------------------------------
!> Read the restart information directly from the model output
!> netcdf file
!> Which routine should find model size?
!> 
subroutine filter_read_restart_direct(state_ens_handle, time, num_extras, use_time_from_file)

type(ensemble_type), intent(inout) :: state_ens_handle
type(time_type),     intent(inout) :: time
integer,             intent(in)    :: num_extras
logical,             intent(in)    :: use_time_from_file

integer                         :: model_size
character(len=256), allocatable :: variable_list(:) !< does this need to be module storage
integer                         :: dart_index !< where to start in state_ens_handle%copies
integer                         :: domain !< loop index

if(get_num_domains()==0) then
   call error_handler(E_ERR, 'filter_read_restart_direct', 'model needs to call add_domain for direct_netcdf_read = .true.')
endif

! read time from input file if time not set in namelist
!> @todo get_model_time should be read_model_time, and should be a namelist to set the filename
!> also need a write_model_time for creating netcdf files
if(use_time_from_file) then
   time = read_model_time(get_input_file(1,1)) ! Any of the restarts?
endif

state_ens_handle%time = time

! read in the data and transpose
dart_index = 1 ! where to start in state_ens_handle%copies - this is modified by read_transpose
do domain = 1, get_num_domains()
   call read_transpose(state_ens_handle, domain, dart_index, limit_mem)
enddo

! Need Temporary print of initial model time?

end subroutine filter_read_restart_direct

!-------------------------------------------------------------------------
!> write the restart information directly into the model netcdf file.
subroutine filter_write_restart_direct(state_ens_handle, num_extras, isprior)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: num_extras ! non restart copies
logical,             intent(in)    :: isprior

integer :: dart_index !< where to start in state_ens_handle%copies
integer :: domain !< loop index
integer :: component_id

!> @todo should we add a blank domain if there is not domains?

! transpose and write out the data
dart_index = 1
do domain = 1, get_num_domains()
   call transpose_write(state_ens_handle, num_extras, domain, dart_index, isprior, limit_mem, single_precision_output)
enddo

end subroutine filter_write_restart_direct

!-------------------------------------------------------
end module state_vector_io_mod
