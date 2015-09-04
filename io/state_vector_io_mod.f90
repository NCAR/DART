! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> IO for the state vector. \n
!> Idea is to be generic. \n
!> Don't want to go to the filesystem twice wrf => wrf_to_dart => dart => dart_to_wrf \n
!> Get a list of variables in the state \n
!> Get info about their dimenions \n
!> Every tasks needs the dimensions of each state variable to calculate
!> what it is going to recieve in the IO transpose
!> \par Aim:
!>
!>To limit the transpose in two ways:
!>
!>1. Limit how much of the state vector you can read at once.
!>2. Limit the number of tasks involved in a transpose.
!>
!>What you (potentially) gain from this:
!>
!>* Don't have to have the whole state vector.
!>* Don't have to use a parallel IO library.
!>
!>If limit 1 > state vector size and limit 2 > number of tasks, you have the regular transpose, except you are reading directly from a netcdf file, not a dart state vector
!> file.
!>
!> * Reading with limited processors is easy, because you just have muliple readers
!> duplicating the read.
!> * Writing with limitied processors is a bit more involved because it is no longer
!> simply duplicate work.  Every processor has something it needs to contribute to 
!> the write. Thus, there is a second stage of data aggregation if limit_procs <
!> task_count.

module state_vector_io_mod

!> \defgroup state_vector_io state_vector_io
!> Contains all the routines and variables to deal with a limited tranpose
!> You can limit the transpose by memory using <code>limit_mem</code> and by
!> processors using <code> limit_procs </code>
!> @{

use types_mod,            only : r8, MISSING_R8, digits12, i4

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id,&
                                 datasize

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, &
                                 is_single_restart_file_in, map_task_to_pe, &
                                 get_copy_owner_index, put_copy
                                 

use utilities_mod,        only : error_handler, nc_check, check_namelist_read, &
                                 find_namelist_in_file, nmlfileunit,           &
                                 do_nml_file, do_nml_term, file_exist,         &
                                 E_MSG, E_ERR, get_unit, ascii_file_format, &
                                 close_file, dump_unit_attributes, &
                                 register_module, set_output

use assim_model_mod,      only : get_model_size, clamp_or_fail_it,             &
                                do_clamp_or_fail, &
                               pert_model_state, assim_model_type

use time_manager_mod,     only : time_type, read_time, write_time, &
                                 get_time

use io_filenames_mod,     only : get_input_file, get_output_file, io_filenames_init

use state_structure_mod,  only : get_domain_size, get_num_variables,           &
                                 get_dim_lengths, get_num_dims, get_dim_ids,   &
                                 get_variable_name, get_variable_size,         &
                                 get_dim_name, get_dim_length,                 &
                                 get_unique_dim_name, get_num_unique_dims,     &
                                 get_unique_dim_length, get_sum_variables,     &
                                 get_sum_variables_below, set_var_id,          &
                                 get_num_domains

use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian

!> @todo  This should go through assim_model_mod
use model_mod, only : write_model_time, read_model_time

use copies_on_off_mod

use netcdf

use mpi


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

private

public :: state_vector_io_init, &
          netcdf_filename,      &
          read_transpose,       &
          transpose_write,      &
          netcdf_filename_out,  &
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
character(len=256) :: netcdf_filename !< needs to be different for each task
character(len=256) :: netcdf_filename_out !< needs to be different for each task

! Module storage for writing error messages
character(len = 255) :: msgstring

! Logical flag for initialization of module
logical              :: module_initialized = .false.

! Global storage for default restart formats
character(len = 16) :: read_format = "unformatted", write_format = "unformatted"


! keep track of whether the initial conditions were netcdf files
logical :: netcdf_input = .false.

! namelist variables with default values
! Aim: to have the regular transpose as the default
integer :: limit_mem = HUGE(1_i4)!< This is the number of elements (not bytes) so you don't have times the number by 4 or 8
integer :: limit_procs = 100000!< how many (~maximum) processors you want involved in each transpose.
logical :: time_unlimited = .true. ! You need to keep track of the time.
logical :: single_precision_output = .false. ! Allows you to write r4 netcdf files even if filter is double precision

logical  :: single_restart_file_out = .true.
! Size of perturbations for creating ensembles when model won't do it
real(r8) :: perturbation_amplitude  = 0.2_r8
! write_binary_restart_files  == .true.  -> use unformatted file format. 
!                                     Full precision, faster, smaller,
!                                     but not as portable.
logical  :: write_binary_restart_files = .false.


namelist /  state_vector_io_nml / limit_mem, limit_procs, time_unlimited, &
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
   call read_transpose(state_ens_handle, domain, dart_index)
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
   call transpose_write(state_ens_handle, num_extras, domain, dart_index, isprior)
enddo

end subroutine filter_write_restart_direct


!-------------------------------------------------
!> Read in variables from model restart file and transpose so that every processor
!> has all copies of a subset of state variables (fill state_ens_handle%copies)
!> Read and transpose data according to the memory limit imposed by
!> limit_mem AND the task limit imposed by limit_procs
!> limit_procs is used to devide the pes into groups.  Note that the
!> groups are not created using mpi_group_incl.
!>
!> Trying to put in single file read for small models.
subroutine read_transpose(state_ens_handle, domain, dart_index)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains

integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: block_size !< number of state elements in a block
integer :: count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: remainder
integer :: start_rank
integer :: group_size
integer :: recv_start, recv_end
integer :: send_start, send_end
integer :: ensemble_member !< the ensmeble_member you are receiving.
integer :: dummy_loop
integer :: my_copy !< which copy a pe is reading, starting from 0 to match pe
integer :: c !< copies_read loop index
integer :: copies_read
integer :: num_state_variables


! single file
integer :: iunit
type(time_type) :: ens_time

netcdf_input = .true.

ens_size = state_ens_handle%num_copies ! have the extras, incase you need to read inflation restarts

my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

copies_read = 0

COPIES: do c = 1, ens_size
   if (copies_read >= ens_size) exit

   ! what to do if a variable is larger than the memory limit?
   start_var = 1 ! read first variable first
   starting_point = dart_index ! position in state_ens_handle%copies

   ! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
   if ( task_count() >= ens_size ) then
      my_copy = my_pe
      call get_pe_loops(my_pe, ens_size, group_size, recv_start, recv_end, send_start, send_end)
   else
      my_copy = copies_read + my_pe
      call get_pe_loops(my_pe, task_count(), group_size, recv_start, recv_end, send_start, send_end)
   endif

   ! open netcdf file
   ! You have already opened this once to read the variable info. Should you just leave it open
   ! on the readers?
   if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader

      if (query_read_copy(my_copy - recv_start+ 1)) then
         netcdf_filename = get_input_file((my_copy - recv_start +1), domain)
         print*, 'opening netcdf_filename ', trim(netcdf_filename)
         ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
         call nc_check(ret, 'read_transpose opening', netcdf_filename)
      endif

   endif

   ! Reading of the state variables is broken up into
   do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be read
      end_var = calc_end_var(start_var, domain)
      if ((my_task_id() == 0) .and. (c == 1)) print*, 'start_var, end_var', start_var, end_var
      block_size = get_sum_variables(start_var, end_var, domain)

      if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader
         if (query_read_copy(my_copy - recv_start + 1)) then

            allocate(var_block(block_size))
            call read_variables(var_block, start_var, end_var, domain)

         endif
      endif

      start_rank = mod(get_sum_variables_below(start_var,domain), task_count())

      ! loop through and post recieves
      RECEIVING_PE_LOOP: do recv_pe = recv_start, recv_end

         ! work out count on the receiving pe
         count = block_size/task_count()
         remainder = mod(block_size, task_count())

         ! mop up leftovers CHECK THESE.
         if ( (start_rank <= recv_pe) .and. (recv_pe) < (start_rank + remainder)) count = count + 1
         if ( recv_pe < (start_rank + remainder - task_count() )) count = count + 1
         ending_point = starting_point + count -1

         ! work out i for the receiving pe
         i = find_start_point(recv_pe, start_rank)

         if (my_pe == recv_pe) then ! get ready to recieve from each reader

            ensemble_member = 1 + copies_read

            RECEIVE_FROM_EACH: do sending_pe = send_start, send_end ! how do we know ens_size?

               if (query_read_copy(sending_pe + copies_read - recv_start + 1)) then

                  if(sending_pe == recv_pe) then ! just copy
                     ! The row is no longer sending_pe + 1 because it is not just
                     ! the first ens_size processors that are sending
                     state_ens_handle%copies(ensemble_member, starting_point:ending_point ) = &
                     var_block(i:count*task_count():task_count())
                  else ! post receive
                     call receive_from(map_pe_to_task(state_ens_handle, sending_pe), &
                                    state_ens_handle%copies(ensemble_member, starting_point:ending_point))
                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo RECEIVE_FROM_EACH

            ! update starting point

            starting_point = starting_point + count

         elseif ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! sending

            if (query_read_copy(my_copy - recv_start + 1)) then
               call send_to(map_pe_to_task(state_ens_handle, recv_pe), &
                           var_block(i:count*task_count():task_count()))
            endif

         endif

      enddo RECEIVING_PE_LOOP

      start_var = end_var + 1

      if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! reader
         if (query_read_copy(my_copy - recv_start + 1)) then
            deallocate(var_block)
         endif
      endif

   enddo

   ! keep track of how many copies have been read.
   copies_read = copies_read + task_count()

   ! close netcdf file
   if ((my_pe >= send_start) .and. (my_pe <= send_end)) then ! I am a reader
      if (query_read_copy(my_copy - recv_start + 1)) then
         ret = nf90_close(ncfile)
         call nc_check(ret, 'read_transpose closing', netcdf_filename)
      endif
   endif

enddo COPIES

dart_index = starting_point

end subroutine read_transpose

!-------------------------------------------------
!> Transpose from state_ens_handle%copies to the writers according to 
!> the memory limit imposed by limit_mem AND the task limit imposed by limit_procs
!> limit_procs is used to devide the pes into groups.  Note that the
!> groups are not created using mpi_group_incl.
!> 
!> Two stages of collection.
!> See transpose_write.pdf for explanation of a, k, y.
!> 
subroutine transpose_write(state_ens_handle, num_extras, domain, dart_index, isprior)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: num_extras ! non restart copies
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index
logical,             intent(in)    :: isprior


integer :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer :: num_vars !< number of variables in a block
integer :: count !< number of elements to send
integer :: starting_point!< position in state_ens_handle%copies
integer :: ending_point
integer :: ens_size !< ensemble size
integer :: remainder
integer :: start_rank
integer :: group_size
integer :: recv_start, recv_end
integer :: send_start, send_end
integer :: ensemble_member
integer :: g !< group loop index
integer :: num_groups !< number of groups the transpose is split into.  Only relevant if limit_procs < task_count()
integer :: assembling_ensemble !< which ensemble the collectors are assembling
integer :: my_group !< which group my_pe is a member of
integer :: num_in_group !< how many processors are in a group
integer :: dummy_loop, j
integer :: a, k, n, owner, y, sub_block
integer :: my_copy !< which copy a pe is reading, starting from 0 to match pe
integer :: c !< copies_read loop index
integer :: copies_written

! mpi_type variables
integer, allocatable :: array_of_blocks(:)
integer, allocatable :: array_of_displacements(:)
integer              :: num_blocks
integer              :: ierr !< mpi error (all errors are fatal so I don't bother checking this
integer              :: collector_type !< mpi derived datatype
integer status(MPI_STATUS_SIZE)

character(len=256)      :: msgstring

integer :: num_state_variables

! single file
integer :: iunit
type(time_type) :: dart_time


ens_size = state_ens_handle%num_copies ! have the extras incase you want to read inflation restarts
my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

copies_written = 0

COPIES : do c = 1, ens_size
   if (copies_written >= ens_size) exit

   start_var = 1 ! collect first variable first
   starting_point = dart_index ! position in state_ens_handle%copies
   a = 0 ! start at group 1 element 1

   ! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
   if ( task_count() >= ens_size ) then
      my_copy = my_pe
      call get_pe_loops(my_pe, ens_size, group_size, send_start, send_end, recv_start, recv_end) ! think I can just flip send and recv for transpose_write?
   else
      my_copy = copies_written + my_pe
      call get_pe_loops(my_pe, task_count(), group_size, send_start, send_end, recv_start, recv_end)
   endif

   ! writers open netcdf output file. This is a copy of the input file
   if (my_pe < ens_size) then
      if ( query_write_copy(my_copy - recv_start + 1)) then
         netcdf_filename_out = get_output_file((my_copy - recv_start +1), domain, isprior)

         if(file_exist(netcdf_filename_out)) then
            ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
            call nc_check(ret, 'transpose_write: opening', trim(netcdf_filename_out))
         else ! create output file if it does not exist
            write(msgstring, *) 'Creating output file ', trim(netcdf_filename_out)
            call error_handler(E_MSG,'state_vector_io_mod:', msgstring)
            dart_time = state_ens_handle%time(my_copy - recv_start + 1)
            call create_state_output(netcdf_filename_out, domain, dart_time)
         endif
      endif

   endif

   do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be sent to writer
      end_var = calc_end_var(start_var, domain)
      if ((my_task_id() == 0) .and. (c == 1 )) print*, 'start_var, end_var', start_var, end_var
      num_vars = get_sum_variables(start_var, end_var, domain)

      if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector
         if (query_write_copy(my_copy - send_start + 1)) then
            allocate(var_block(num_vars))
         endif
      endif

      start_rank =  mod(get_sum_variables_below(start_var,domain), task_count())

      SENDING_PE_LOOP: do sending_pe = send_start, send_end

         ! work out count on the sending pe
         count = num_vars/task_count()
         remainder = mod(num_vars, task_count())

         ! mop up leftovers CHECK THESE.
         if ( (start_rank <= sending_pe) .and. (sending_pe) < (start_rank + remainder)) count = count + 1
         if ( sending_pe < (start_rank + remainder - task_count() )) count = count + 1
         ending_point = starting_point + count -1

         ! work out i for the sending_pe
         i = find_start_point(sending_pe, start_rank)

         if (my_pe /= sending_pe ) then ! post recieves
            if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector
               if (query_write_copy(my_copy - send_start + 1)) then
                  call receive_from(map_pe_to_task(state_ens_handle, sending_pe), var_block(i:count*task_count():task_count()))
               endif
            endif

         else ! send to the collector

            ensemble_member = 1 + copies_written

            do recv_pe = recv_start, recv_end ! no if statement because everyone sends

               if (query_write_copy(recv_pe + copies_written - send_start + 1)) then

                  if ( recv_pe /= my_pe ) then
                     call send_to(map_pe_to_task(state_ens_handle, recv_pe), state_ens_handle%copies(ensemble_member, starting_point:ending_point))
                  else ! if sender = receiver just copy

                     var_block(i:count*task_count():task_count()) = state_ens_handle%copies(ensemble_member, starting_point:ending_point)

                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo

            ! update starting point
            starting_point = starting_point + count

         endif

      enddo SENDING_PE_LOOP

      if (limit_procs >= task_count()) then ! no need to do second stage

         ! I think for now, the single file should enter this.

         if (my_pe < ens_size) then ! I am a writer
            if ( query_write_copy(my_copy + 1)) then
               !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
               if (my_copy <= state_ens_handle%num_copies - num_extras) then ! actual copy, may need clamping
                  call write_variables_clamp(var_block, start_var, end_var, domain)
               else ! extra copy, don't clamp
                  call write_variables(var_block, start_var, end_var, domain)
               endif
               deallocate(var_block)
            endif
         endif

      else ! Need to aggregate onto the writers (ens_size writers) Is there a better way to do this?
           ! I don't think you enter this if task_count < ens_size because limit_procs = task_count()

         if ((my_pe >= recv_start) .and. (my_pe <= recv_end)) then ! I am a collector

            if (query_write_copy(my_pe - send_start + 1)) then

               assembling_ensemble = my_pe - recv_start + 1
               num_groups = task_count() / limit_procs
               if (mod(task_count(), limit_procs) /= 0) then ! have to do somthing else
                  num_groups = num_groups + 1
                  if (mod(task_count(), limit_procs) < ens_size) num_groups = num_groups - 1 ! last group is big
               endif

               my_group = my_pe / limit_procs + 1
               if (my_group > num_groups) my_group = num_groups ! last group is big

               do g = 2, num_groups ! do whole ensemble at once

                  ! only group sending and first group need to be involved
                  if ((my_group == g) .or. (my_group == 1)) then

                     ! create datatype for the data being sent to the writers - same across the ensemble
                     ! need to find size of group g. Only the last group could be a different size
                     if (g < num_groups) then
                        num_in_group = limit_procs
                     else
                        num_in_group = get_group_size(task_count(), ens_size)
                     endif

                     if (a == 0) then ! group 1, element 1 starts the var_block, g cannot be the owner

                        owner = 1
                        y = limit_procs
                        sub_block = num_vars - (g-1)*limit_procs
                        num_blocks = sub_block / task_count()
                        remainder = mod(sub_block, task_count())
                        if (remainder > 0) num_blocks = num_blocks + 1
                        allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))

                        array_of_displacements(1) = (g-1)*limit_procs
                        array_of_displacements(2) = array_of_displacements(1) + task_count()
                        array_of_blocks(1) = num_in_group

                        if (remainder < num_in_group ) then
                           array_of_blocks(num_blocks) = sub_block - task_count()*(num_blocks-1)
                        else
                           array_of_blocks(num_blocks) = num_in_group
                        endif

                     else ! have to do something else

                        ! calculate which group last element of a is in. a starts at group 1 element 1
                        do j = 1, num_groups
                           if ( a <= cumulative_tasks(j, ens_size) ) then
                              owner = j
                              exit
                           endif
                        enddo

                        ! calulate k:
                        if ( owner == 1 ) then
                           k = a
                        else
                           k = a - cumulative_tasks(owner -1, ens_size)
                        endif

                        if (k == 0) then
                           owner = owner + 1
                           if (owner == num_groups + 1) owner = 1
                        endif

                        y = get_group_size(owner*limit_procs -1, ens_size) - k

                        ! find number of tasks between owner and group 1?
                        n = cumulative_tasks(num_groups, ens_size) - cumulative_tasks(owner, ens_size)

                        !if (my_pe == 0) print*, 'g = ', g, 'owner =', owner, 'n = ', n, 'k = ', k, 'y =', y, 'num_blocks', num_blocks, 'num_vars', num_vars

                        ! find number of blocks:
                        sub_block = num_vars - y - n
                        num_blocks = sub_block / task_count()
                        if ( g >= owner ) num_blocks = num_blocks + 1 ! for y and for any blocks in n

                        remainder = mod( sub_block, task_count() )

                        if (remainder >= cumulative_tasks(g, ens_size) ) then

                           num_blocks = num_blocks + 1
                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = num_in_group

                        elseif ( (cumulative_tasks(g-1, ens_size) < remainder) .and. ( remainder < cumulative_tasks(g, ens_size)) ) then

                           num_blocks = num_blocks + 1 ! ragged end
                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = remainder - cumulative_tasks(g-1, ens_size)

                        else

                           allocate(array_of_blocks(num_blocks), array_of_displacements(num_blocks))
                           array_of_blocks(num_blocks) = num_in_group

                        endif

                        if ( g == owner ) then
                           array_of_displacements(1) = 0 ! zero offset for mpi_type_indexed
                           array_of_displacements(2) = task_count() - k  ! zero offset
                           array_of_blocks(1) = y
                        elseif ( g > owner) then
                           array_of_displacements(1) = y + cumulative_tasks(g-1, ens_size) - cumulative_tasks(owner, ens_size)
                           array_of_displacements(2) = array_of_displacements(1) + task_count()
                           array_of_blocks(1) = num_in_group
                        else
                           array_of_displacements(1) = y + n + cumulative_tasks(g-1, ens_size) ! y + n + offest from start of group 1
                           array_of_displacements(2) = array_of_displacements(1) + task_count()
                           array_of_blocks(1) = num_in_group
                        endif

                     endif

                     array_of_blocks(2:num_blocks - 1) = num_in_group

                     do i = 3, num_blocks
                        array_of_displacements(i) = array_of_displacements(i-1) + task_count()
                     enddo
   
                     ! check you are not going over num_vars - you can probably pull this, it was just for debugging.
                     if(my_pe == 0) then
                        if( (array_of_displacements(num_blocks) + array_of_blocks(num_blocks))  > num_vars) then
                           print*, '++++ OVER ++++', num_vars - (array_of_displacements(num_blocks) + array_of_blocks(num_blocks)), 'last block', array_of_blocks(num_blocks), 'last disp', array_of_displacements(num_blocks), 'num_vars', num_vars, 'y', y
                           print*, 'remainder', remainder, cumulative_tasks(g-1, ens_size), cumulative_tasks(g,ens_size), num_blocks
                        endif
                     endif

                     if ( datasize == mpi_real4 ) then

                        call mpi_type_indexed(num_blocks, array_of_blocks, array_of_displacements, mpi_real4, collector_type, ierr)

                     else ! double precision

                        call mpi_type_indexed(num_blocks, array_of_blocks, array_of_displacements, mpi_real8, collector_type, ierr)

                     endif

                     call mpi_type_commit(collector_type, ierr)

                     ! collectors -> writers

                     recv_pe = assembling_ensemble - 1
                     sending_pe = recv_pe + (g-1)*limit_procs
                     if (my_pe == recv_pe) then
                        call mpi_recv(var_block, 1, collector_type, map_pe_to_task(state_ens_handle,sending_pe), 0, mpi_comm_world, status, ierr)
                     elseif (my_pe == sending_pe) then
                        call mpi_send(var_block, 1, collector_type, map_pe_to_task(state_ens_handle,recv_pe), 0, mpi_comm_world, ierr)
                     endif

                     call mpi_type_free(collector_type, ierr)
                     deallocate(array_of_blocks, array_of_displacements)

                  endif

               enddo

            endif

            if (my_pe < ens_size) then ! I am a writer
               if(query_write_copy(my_copy + 1)) then
                  !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
                  if (my_copy <= state_ens_handle%num_copies - num_extras) then ! actual copy, may need clamping
                     call write_variables_clamp(var_block, start_var, end_var, domain)
                  else ! extra copy, don't clamp
                     call write_variables(var_block, start_var, end_var, domain)
                  endif
               endif
            endif

            if(query_write_copy(my_copy - send_start + 1)) then
               deallocate(var_block) ! all collectors have var_block
            endif

         endif

      endif

      start_var = end_var + 1
      ! calculate a:
      a = mod(num_vars - (task_count() - a), task_count())

   enddo

   ! keep track of how many copies have been written
   copies_written = copies_written + task_count()

   ! close netcdf file
   if (my_copy < ens_size ) then ! I am a writer
      if (query_write_copy(my_copy + 1)) then
         ret = nf90_close(ncfile_out)
         call nc_check(ret, 'transpose_write', 'closing')
      endif
   endif

enddo COPIES

dart_index = starting_point

end subroutine transpose_write

!-------------------------------------------------
!> Calculate how many variables to read in one go.
function calc_end_var(start_var, domain)

integer              :: calc_end_var !< end variable index
integer, intent(in)  :: start_var !< start variable index
integer, intent(in)  :: domain

integer :: i, count
integer :: num_state_variables
integer, allocatable :: num_elements(:) !< cummulative size

num_state_variables = get_num_variables(domain)

allocate(num_elements(num_state_variables - start_var + 1))

calc_end_var = num_state_variables ! assume you can fit them all to start with

count = 0

do i = 1, num_state_variables - start_var + 1
   num_elements(i) = get_sum_variables(start_var, start_var + count, domain)
   count = count + 1
enddo

count = 1
do i = start_var, num_state_variables

   if (start_var == num_state_variables) then
      calc_end_var = num_state_variables
      exit
   endif

   if (count >= num_state_variables) then
      calc_end_var = num_state_variables
      exit
   endif

   if(count + 1> size(num_elements)) then
      calc_end_var = num_state_variables
      exit
   endif

   if (num_elements(count+1) >= limit_mem ) then
      calc_end_var =  i
      exit
   endif
   count = count + 1
enddo

deallocate(num_elements)

end function calc_end_var

!-------------------------------------------------
!> Read in variables from start_var to end_var
!> FIXME: At the moment, this code is assuming that the variables in the state start
!> at (1,1,1) and that the whole variable is read. This is not the case for 
!> Tiegcm and CLM.  
subroutine read_variables(var_block, start_var, end_var, domain)

real(r8),           intent(inout) :: var_block(:)
integer,            intent(in)    :: start_var
integer,            intent(in)    :: end_var
integer,            intent(in)    :: domain

integer :: i
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id

start_in_var_block = 1

do i = start_var, end_var

   var_size = get_variable_size(domain, i)

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'read_variables','inquire variable id')

   ret = nf90_get_var(ncfile, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'read_variables','reading')

   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine read_variables


!-------------------------------------------------
!> Create the output files
!> ncfile_out is global - is this ok?
!> I have removed fresh_netcdf_file, since filter_write_restart_direct
!> can add a blank domain.
!> A 'blank' domain is one variable called
!> state, with dimension = model size.
!> It is used when the model has not suppled any netdcf info but direct_netcdf_write = .true.

subroutine create_state_output(filename, dom, dart_time)

character(len=256), intent(in) :: filename
integer,            intent(in) :: dom !< domain, not sure whether you need this?
type(time_type),    intent(in) :: dart_time

integer :: ret !> netcdf return code
integer :: create_mode
integer :: i, j ! loop variables
integer :: new_dimid
integer :: new_varid
integer :: ndims
integer :: xtype ! precision for netcdf file
logical :: time_dimension_exists
integer :: dimids(NF90_MAX_VAR_DIMS)

time_dimension_exists = .false.

! What file options do you want
create_mode = ior(NF90_CLOBBER, NF90_64BIT_OFFSET)
ret = nf90_create(filename, create_mode, ncfile_out)
call nc_check(ret, 'create_state_output', 'creating')

! define dimensions, loop around unique dimensions
do i = 1, get_num_unique_dims(dom)
   ret = nf90_def_dim(ncfile_out, get_unique_dim_name(dom, i), get_unique_dim_length(dom, i), new_dimid)
   !> @todo if we already have a unique names we can take this test out
   if(ret /= NF90_NOERR .and. ret /= NF90_ENAMEINUSE) then
      call nc_check(ret, 'create_state_output', 'defining dimensions')
   endif
enddo

! define variables
do i = 1, get_num_variables(dom) ! loop around state variables
   ! double or single precision?
   ndims = get_num_dims(dom, i)

   if (single_precision_output) then
      xtype = nf90_real
   else ! write output that is the precision of filter
      if (r8 == digits12) then ! datasize = MPI_REAL8  ! What should we be writing?
         xtype = nf90_double
      else
         xtype = nf90_real
      endif
   endif

   ! query the dimension ids
   do j = 1, get_num_dims(dom, i)
      ret = nf90_inq_dimid(ncfile_out, get_dim_name(dom, i, j), dimids(j))
      call nc_check(ret, 'create_state_output', 'querying dimensions')
   enddo

   ret = nf90_def_var(ncfile_out, trim(get_variable_name(dom, i)), xtype=xtype, dimids=dimids(1:get_num_dims(dom, i)), varid=new_varid)
   call nc_check(ret, 'create_state_output', 'defining variable')
      !variable_ids(i, dom) = new_varid
   call set_var_id(dom, i, new_varid)
enddo

ret = nf90_enddef(ncfile_out)
call nc_check(ret, 'create_state_output', 'end define mode')

call write_model_time(ncfile_out, dart_time)

end subroutine create_state_output

!-------------------------------------------------
!> Write variables from start_var to end_var
!> no clamping
subroutine write_variables(var_block, start_var, end_var, domain)

real(r8), intent(inout) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var
integer,  intent(in) :: domain 

integer :: i
integer :: count_displacement
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id 

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'write_variables', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'write_variables', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables

!-------------------------------------------------
!> Write variables from start_var to end_var
!> For actual ensemble members
subroutine write_variables_clamp(var_block, start_var, end_var, domain)

real(r8), intent(inout) :: var_block(:)
integer,  intent(in) :: start_var
integer,  intent(in) :: end_var
integer,  intent(in) :: domain 

integer :: i
integer :: count_displacement
integer :: start_in_var_block
integer :: var_size
integer, allocatable :: dims(:)
integer :: var_id 

start_in_var_block = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)

   ! check whether you have to do anything to the variable, clamp or fail
   if (do_clamp_or_fail(i, domain)) then
      call clamp_or_fail_it(i, domain, var_block(start_in_var_block:start_in_var_block+var_size-1))
   endif

   ! number of dimensions and length of each
   allocate(dims(get_num_dims(domain, i)))

   dims = get_dim_lengths(domain, i)

   ret = nf90_inq_varid(ncfile_out, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'write_variables_clamp', 'getting variable id')

   ret = nf90_put_var(ncfile_out, var_id, var_block(start_in_var_block:start_in_var_block+var_size-1), count=dims)
   call nc_check(ret, 'write_variables_clamp', 'writing')
   start_in_var_block = start_in_var_block + var_size

   deallocate(dims)

enddo

end subroutine write_variables_clamp

!-------------------------------------------------
!> Find pes for loop indices
subroutine get_pe_loops(pe, ens_size, group_size, recv_start, recv_end, send_start, send_end)

integer, intent(in)  :: pe
integer, intent(in)  :: ens_size
integer, intent(out) :: group_size !< size of the group I am in
integer, intent(out) :: recv_start !< for RECIEVING_PE_LOOP
integer, intent(out) :: recv_end !< for RECIEVING_PE_LOOP
integer, intent(out) :: send_start !< for RECEIVE_FROM_EACH_LOOP
integer, intent(out) :: send_end !< for RECEIVE_FROM_EACH_LOOP

group_size = get_group_size(pe, ens_size)

if ( limit_procs > task_count() ) limit_procs = task_count()
!if (my_task_id() == 0) print*, 'limit_procs', limit_procs, 'limit_mem', limit_mem

! limit_procs needs to be greater than ens_size

if( group_size > limit_procs ) then ! last group is large because of odd numbers
   recv_start = ((task_count() / limit_procs) - 1) * limit_procs
else
   recv_start = (pe / limit_procs ) * limit_procs
endif

recv_end = recv_start + group_size -1

send_start  = recv_start
send_end = send_start + ens_size -1

end subroutine get_pe_loops

!-------------------------------------------------
!> Find group size for processor limited transpose
!> groups are 1:n, n+1:m, m+1:l ...
function get_group_size(pe, ens_size)

integer, intent(in)  :: pe
integer, intent(in)  :: ens_size
integer              :: get_group_size

integer :: num_groups
integer :: remainder

num_groups = task_count() / limit_procs

remainder = mod(task_count(), limit_procs)

if ( remainder > 0 ) then ! the last group has a diffent size

   if ( remainder < ens_size) then ! need to join the last two groups together

      if ( (pe + 1) > limit_procs*(num_groups - 1) ) then
         get_group_size = remainder + limit_procs
      else
         get_group_size = limit_procs
      endif

   else ! last group is smaller than the rest

      if ( (pe + 1 ) > limit_procs*num_groups ) then ! pe is a member of the last group
         get_group_size = remainder
      else
         get_group_size = limit_procs
      endif

   endif

else ! all same size
   get_group_size = limit_procs
endif

end function get_group_size

!-------------------------------------------------------
!> Find i, the start point in var_block for a given recv_pe
function find_start_point(recv_pe, start_rank)

integer, intent(in)  :: recv_pe !< the receiver
integer, intent(in)  :: start_rank !< the pe that owns the 1st element of the var_block
integer              :: find_start_point

if (start_rank < recv_pe) then
   find_start_point = recv_pe - start_rank + 1
elseif(start_rank > recv_pe) then
   find_start_point = recv_pe + task_count() - start_rank + 1
else ! recv_pe = start_rank
   find_start_point = 1
endif

end function find_start_point

!------------------------------------------------------
!> Given a group, finds the total number of tasks from group 1
!> up to and including that group
function cumulative_tasks(group, ens_size)

integer, intent(in) :: group
integer, intent(in) :: ens_size !< for get group size
integer             :: cumulative_tasks

integer :: i

cumulative_tasks = 0

! what if you give it a group > num_groups? Or a negative group?

do i = 1, group - 1
   cumulative_tasks = cumulative_tasks + limit_procs
enddo

! just in case group is the last group
cumulative_tasks = cumulative_tasks + get_group_size(group*limit_procs -1, ens_size)

end function cumulative_tasks

!-------------------------------------------------------
!> @}
end module state_vector_io_mod
