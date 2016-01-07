! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module mpi_utilities_mod

!-----------------------------------------------------------------------------
!
!   USED FOR IMPLEMENTATIONS THAT DO NOT WANT TO USE MPI.
!   ALSO REQUIRED FOR INTEGRATE_MODEL PROGRAM TO WORK WITH SHELL
!   DRIVEN MODEL ADVANCES. ACTS AS IF THERE IS ONLY A SINGLE TASK.
!
!   Programs using this module instead of the actual MPI routines do not
!   need to be compiled with the MPI wrapper commands (e.g. mpif90).
!   In most cases it will be better to compile with the fortran compiler
!   directly.  (On some platforms this is required.)
!
!   These routines mimic the actual interfaces to the MPI (Message Passing
!   Interface) library but do not use MPI.
!
!    # initialize_mpi_utilities()  Subroutine that initializes MPI and sets
!                                  local values needed later.  Must be called
!                                  before any other routine here.
!
!    # finalize_mpi_utilities()  Subroutine that shuts down MPI cleanly.
!                                Must be called before program exits, and no
!                                other routines here can be used afterwards.
!
!    # task_count()       Function that returns the total number of MPI tasks.
!
!    # my_task_id()       Function that returns my task number.  Note that
!                         in the MPI world task numbers run from 0 to N-1.
!
!    # send_to()          Subroutine which sends a 1D data array
!                         synchronously to another task (point-to-point).
!
!    # receive_from()     Subroutine which receives a 1D data array
!                         synchronously from another task (point-to-point).
!
!    # task_sync()        Subroutine that only returns after every task has
!                         reached this same location in the code.
!        
!    # array_broadcast()  Subroutine that sends a copy of the entire data 
!                         array to all other tasks. 
!                
!  *** exit_all()         Subroutine that substitutes for the intrinsic exit.
!                         It calls MPI_Abort() to force other MPI tasks to
!                         exit as well in case of error.
!
!    * transpose_array()  Subroutine that transposes a 2D array
!                         from column-major to row-major or back.
!
!   ** array_distribute() Subroutine that distributes a data array across the
!                         other tasks, so each task gets a non-overlapping 
!                         subset of the data.
!                
!   MPI cover routines more specific for DART and hopefully more useful.
!
!    # iam_task0()        Function which returns .TRUE. if task id is 0,
!                         .FALSE. for anything else.
!
!    # broadcast_send()   Subroutine which takes two r8 arrays and broadcasts
!                         them to all other tasks.  If sending ID is not the
!                         same as the local task ID, an error is returned.
!                         Does not return until all other tasks have called
!                         recv to pick up the data.
!
!    # broadcast_recv()   Subroutine which receives two r8 arrays from the 
!                         sending task ID.  If the sending ID is the same as
!                         the local task ID, an error is returned.  All other
!                         tasks must call recv before this routine returns.
!
!    * sum_across_tasks() Subroutine which takes a single integer argument
!                         from each task, and returns the sum of all integers
!                         across all tasks back to all tasks.  All tasks must
!                         call this routine before it can compute and return
!                         the value.
!
!   Lower level utility routines which interact with the utilities_mod.f90 
!   code to open a named pipe per MPI task, read and write from them, and
!   close and/or remove them.
!
!  *** make_pipe()        Function that creates a named pipe (fifo), opens it,
!                         and returns the unit number.  Ok to call if the pipe
!                         already exists or is already open; it will skip
!                         those steps and just return the unit number.  The 
!                         name argument is used as a base and a filename
!                         in the form 'base.NNNN' is generated, where the N's
!                         are the MPI rank number, 0 padded.
!
!  *** destroy_pipe()     The unit number is closed and the pipe file is 
!                         removed.
!
!    * read_pipe()        The character string is read from the pipe.
!                         (Can be overloaded to read ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine blocks until data is available.
!
!    * write_pipe()       The character string is written to the pipe.
!                         (Can be overloaded to write ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine writes and returns immediately.
!
!   Wrappers for system functions.  Covers differences if you run with
!   or without MPI.
!
!  *** shell_execute()    Use the system() command to execute a command string.
!                         Will wait for the command to complete and returns an
!                         error code unless you end the command with & to put
!                         it into background.   Function which returns the rc
!                         of the command, 0 being all is ok.
!
!  *** sleep_seconds()    Wrapper for the sleep command.  Argument is a real
!                         in seconds.  Different systems have different lower
!                         resolutions for the minimum time it will sleep.
!                         Subroutine, no return value.
!
!
!   # code done and tested
! *** both code and interface are done (but untested so far)
!  ** interface with proposed arguments exists but code not complete
!   * interface name only; no arg list devised yet 
!
!-----------------------------------------------------------------------------
! 
! these do not exist - i believe a single transpose will work.  but if not,
! they can be separated into these two, which can either work on a real
! 2D array or a single linearized array which is logically 2D but in reality
! stored in a 1D fortran array:
!
!      transpose_row_major()  Subroutine which transposes a logical 2D array
!                             from column-major to row-major.  The source and
!                             destination arrays must be stored in 1D arrays 
!                             of length (nrows * ncols).
!
!      transpose_col_major()  Subroutine which transposes a logical 2D array
!                             from row-major to column-major.  The source and
!                             destination arrays must be stored in 1D arrays 
!                             of length (nrows * ncols).
!
!-----------------------------------------------------------------------------

use types_mod, only        : r8, digits12
use utilities_mod, only    : register_module, error_handler,             &
                             initialize_utilities, get_unit, close_file, & 
                             E_ERR, E_WARN, E_MSG, E_DBG, finalize_utilities
use time_manager_mod, only : time_type, set_time


implicit none
private


! BUILD TIP 
! Some compilers require an interface block for the system() function;
! some fail if you define one.  If you get an error at link time (something
! like 'undefined symbol _system_') try running the fixsystem script in
! this directory.  It is a sed script that comments in and out the interface
! block below.  Please leave the BLOCK comment lines unchanged.

! !!SYSTEM_BLOCK_EDIT START COMMENTED_OUT
! ! interface block for getting return code back from system() routine
! interface
!  function system(string)    
!   character(len=*) :: string
!   integer :: system         
!  end function system
! end interface
! ! end block                 
! !!SYSTEM_BLOCK_EDIT END COMMENTED_OUT

!   ---- private data for mpi_utilities ----

integer :: myrank          ! my mpi number
integer :: total_tasks     ! total mpi tasks/procs
integer :: my_local_comm   ! duplicate communicator private to this file
integer :: comm_size       ! if ens count < tasks, only the first N participate

public :: initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, reduce_min_max, &
          get_from_fwd, get_from_mean

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len = 129) :: saved_progname = ''

character(len = 129) :: errstring

! Namelist input - placeholder for now; no options yet in this module.
!namelist /mpi_utilities_nml/ x


contains

!-----------------------------------------------------------------------------
! mpi cover routines
!-----------------------------------------------------------------------------

subroutine initialize_mpi_utilities(progname, alternatename)
 character(len=*), intent(in), optional :: progname
 character(len=*), intent(in), optional :: alternatename

! Initialize MPI and query it for global information.  Make a duplicate
! communicator so that any user code which wants to call MPI will not 
! interfere with any outstanding asynchronous requests, accidental tag
! matches, etc.  This routine must be called before any other routine in
! this file, and it should not be called more than once (but it does have
! defensive code in case that happens.)

if ( module_initialized ) then
   ! return without calling the code below multiple times
   write(errstring, *) 'initialize_mpi_utilities has already been called'
   call error_handler(E_WARN,'initialize_mpi_utilities', errstring, source, revision, revdate)
   return
endif

call initialize_utilities(progname, alternatename)

if (present(progname)) then
   if (len_trim(progname) < len(saved_progname)) then
      saved_progname = trim(progname)
   else
      saved_progname = progname(1:len(saved_progname))
   endif
endif

if ( .not. module_initialized ) then
   ! Initialize the module with utilities
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

myrank = 0
total_tasks = 1

! TODO: if there are fewer ensembles than tasks, all the collective routines
! need to take that into account and not participate if they are > comm_size.
comm_size = total_tasks

if (r8 /= digits12) then
   write(errstring, *) "Using real * 4 for datasize of r8"
   call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring,source,revision,revdate)
endif

! non-MPI successfully initialized.
call error_handler(E_MSG,'initialize_mpi_utilities: ','Running single process', &
                   source, revision, revdate)

end subroutine initialize_mpi_utilities

!-----------------------------------------------------------------------------

subroutine finalize_mpi_utilities(callfinalize, async)
 logical, intent(in), optional :: callfinalize
 integer, intent(in), optional :: async

! Shut down cleanly.  Call normal utilities finalize if we have actually
! ever called initialize.  Otherwise there is nothing to do in the null case.


if ( .not. module_initialized ) return

if (saved_progname /= '') then
   call finalize_utilities(saved_progname)
else
   call finalize_utilities()
endif

end subroutine finalize_mpi_utilities


!-----------------------------------------------------------------------------

function task_count()

! Return the total number of MPI tasks.  e.g. if the number of tasks is 4,
! it returns 4.  (The actual task numbers are 0-3.)  For the null mpi utils,
! this always returns 1.

integer :: task_count

if ( .not. module_initialized ) call initialize_mpi_utilities()

task_count = total_tasks

end function task_count


!-----------------------------------------------------------------------------

function my_task_id()

! Return my unique task id.  Values run from 0 to N-1 (where N is the
! total number of MPI tasks.  For the null mpi utils, this is always 0.

integer :: my_task_id

if ( .not. module_initialized ) call initialize_mpi_utilities()

my_task_id = myrank

end function my_task_id


!-----------------------------------------------------------------------------

subroutine task_sync()

! Synchronize all tasks.  This subroutine does not return until all tasks
! execute this line of code.

if ( .not. module_initialized ) call initialize_mpi_utilities()


end subroutine task_sync


!-----------------------------------------------------------------------------

subroutine send_to(dest_id, srcarray, time, label)
 integer, intent(in) :: dest_id
 real(r8), intent(in) :: srcarray(:)
 type(time_type), intent(in), optional :: time
 character(len=*), intent(in), optional :: label

! Send the srcarray to the destination id.
! If time is specified, it is also sent in a separate communications call.  
! This is a synchronous call; it will not return until the destination has 
! called receive to accept the data.  If the send_to/receive_from calls are 
! not paired correctly the code will hang.

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if ((dest_id < 0) .or. (dest_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "destination task id ", dest_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

! this style cannot be easily simulated correctly with one task.
! always throw an error.
write(errstring, '(a)') "cannot call send_to() in the single process case"
call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)

end subroutine send_to


!-----------------------------------------------------------------------------

subroutine receive_from(src_id, destarray, time, label)
 integer, intent(in) :: src_id
 real(r8), intent(out) :: destarray(:)
 type(time_type), intent(out), optional :: time
 character(len=*), intent(in), optional :: label

! Receive data into the destination array from the src task.
! If time is specified, it is received in a separate communications call.  
! This is a synchronous call; it will not return until the source has 
! sent the data.  If the send_to/receive_from calls are not paired correctly 
! the code will hang.

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if ((src_id < 0) .or. (src_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "source task id ", src_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

! this style cannot be easily simulated correctly with one task.
! always throw an error.
write(errstring, '(a)') "cannot call receive_from() in the single process case"
call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)

destarray = 0
if (present(time)) time = set_time(0, 0)

end subroutine receive_from



!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 1D, 2D, 3D inputs?


subroutine transpose_array

! not implemented here yet.  will have arguments -- several of them.

if ( .not. module_initialized ) call initialize_mpi_utilities()

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)

end subroutine transpose_array


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 2D inputs?

subroutine array_broadcast(array, root)
 real(r8), intent(inout) :: array(:)
 integer, intent(in) :: root

! The data array values on the root task will be broadcast to every other
! task.  When this routine returns, all tasks will have the contents of the
! root array in their own arrays.  Thus 'array' is intent(in) on root, and
! intent(out) on all other tasks.

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

! array already has the values, nothing to do.

end subroutine array_broadcast


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 2D inputs?

subroutine array_distribute(srcarray, root, dstarray, dstcount, how, which)
 real(r8), intent(in) :: srcarray(:)
 integer, intent(in) :: root
 real(r8), intent(out) :: dstarray(:)
 integer, intent(out) :: dstcount
 integer, intent(in) :: how
 integer, intent(out) :: which(:)

! 'srcarray' on the root task will be distributed across all the tasks
! into 'dstarray'.  dstarray must be large enough to hold each task's share
! of the data.  The actual number of values returned on each task will be
! passed back in the 'count' argument.  'how' is a flag to select how to
! distribute the data (round-robin, contiguous chunks, etc).  'which' is an
! integer index array which lists which of the original values were selected
! and put into 'dstarray'.

integer :: i

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

dstarray = srcarray
dstcount = size(srcarray)
which = (/ ((i), i=1,size(srcarray))  /)

end subroutine array_distribute

!-----------------------------------------------------------------------------
! DART-specific cover utilities
!-----------------------------------------------------------------------------

function iam_task0()

! Return .TRUE. if my local task id is 0, .FALSE. otherwise.
! (Task numbers in MPI start at 0, contrary to the rules of polite fortran.)

logical :: iam_task0

if ( .not. module_initialized ) call initialize_mpi_utilities()

iam_task0 = (myrank == 0)

end function iam_task0

!-----------------------------------------------------------------------------
subroutine broadcast_send(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)
 integer, intent(in) :: from
 ! really only intent(in) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:)
 real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' matches my local task id.  also, these arrays are
! intent(in) here, but they call a routine which is intent(inout) so they
! must be the same here.

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if (from /= myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "must be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! this does nothing, because the array already has the data.
! the subroutine does validate 'from' to be sure it's a valid task id.
call array_broadcast(array1, from)

end subroutine broadcast_send

!-----------------------------------------------------------------------------
subroutine broadcast_recv(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)
 integer, intent(in) :: from
 ! really only intent(out) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:)
 real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' is not the same as my local task id.  these arrays are
! intent(out) here, but they call a routine which is intent(inout) so they
! must be the same here.

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if (from == myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "cannot be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! this does nothing, because the array already has the data.
! the subroutine does validate 'from' to be sure it's a valid task id.
call array_broadcast(array1, from)

end subroutine broadcast_recv

!-----------------------------------------------------------------------------
subroutine sum_across_tasks(addend, sum)
 integer, intent(in) :: addend
 integer, intent(out) :: sum

! cover routine for MPI all-reduce

if ( .not. module_initialized ) call initialize_mpi_utilities()

sum = addend

end subroutine sum_across_tasks


!-----------------------------------------------------------------------------
! pipe-related utilities
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine block_task()


if ( .not. module_initialized ) call initialize_mpi_utilities()
 
 
end subroutine block_task

!-----------------------------------------------------------------------------
subroutine restart_task()
   
   
if ( .not. module_initialized ) call initialize_mpi_utilities()


end subroutine restart_task


!-----------------------------------------------------------------------------
!    * make_pipe()        Function that creates a named pipe (fifo), opens it,
!                         and returns the unit number.  Ok to call if the pipe
!                         already exists or is already open; it will skip
!                         those steps and just return the unit number.  The 
!                         name argument is used as a base and a filename
!                         in the form 'base.NNNN' is generated, where the N's
!                         are the MPI rank number, 0 padded.
!
!-----------------------------------------------------------------------------

function make_pipe(pipename, exists) result (iunit)
 character(len=*), intent(in) :: pipename
 logical, intent(in), optional :: exists
 integer :: iunit

! Create, open, and return a fortran unit number for a named pipe.
! The local MPI rank number will be appended to the given name to create
! a file of the form 'base.NNNN', where N's are the MPI rank number, 0 padded.
! TODO: based on the total number of tasks get extra style points for
! creating the shortest name necessary; e.g. base.N, base.NN, base.NNN, etc.
!
! If the optional 'exists' flag is not present, then it is not an error
! whether the pipe already exists or not.  It is made if it does not exist, 
! it is opened if not already opened, and the fortran unit number is returned.
! If 'exists' is present then it forces the issue of whether the pipe file
! must exist already or not.  The error handler is called if things aren't 
! as expected.  apologies to tim hoar for the intentional imitation of
! the generic file utilities_mod.f90 code.

logical :: open, there
character(len=128) :: fname
character(len=11) :: format
integer :: rc

if ( .not. module_initialized ) call initialize_mpi_utilities()

write(fname, "(a,i4.4)") trim(pipename)//".", myrank
print *, "fname now = ", trim(fname)

! check to see if the pipe already exists; if so, we've got the unit number
! (directly into the output value) and we're done.  otherwise, make it and
! open it.
inquire (file=fname, exist=there, opened=open, number=iunit, form=format)

if (.not. open) then

   if (.not. there) then
      ! make pipe file; mkfifo should be standard on any unix/linux system.
      rc = system('mkfifo '//trim(fname)//' '//char(0))

      ! and check to be sure it was made
      inquire (file=fname, exist=there)

      if (.not. there) then
        write(errstring, *) 'mkfifo command failed to create '//trim(fname)
        call error_handler(E_ERR,'make_pipe', errstring, source, revision, revdate)
      endif
   endif

   ! open pipe using an available unit number
   iunit = get_unit()
   open(unit=iunit, file=fname)

endif

! iunit contains the function return value.

end function make_pipe


!-----------------------------------------------------------------------------
!    * destroy_pipe()     The unit number is closed and the pipe file is 
!                         removed.
!
subroutine destroy_pipe(iunit)
 integer, intent(in) :: iunit

character(len=128) :: pipename
integer :: ios, rc

if ( .not. module_initialized ) call initialize_mpi_utilities()

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'destroy_pipe', errstring, source, revision, revdate)


! general idea is:

! call inquire to get name
inquire(unit=iunit, name=pipename, iostat=ios)
if (ios /= 0) then
   write(errstring, '(a,i4)') 'failure trying to inquire about unit ', iunit
   call error_handler(E_ERR,'destroy_pipe', errstring, source, revision, revdate)
endif

call close_file(iunit)

! remove echo when we trust this command.
rc = system('echo rm -f '//trim(pipename)//' '//char(0))


end subroutine destroy_pipe

!-----------------------------------------------------------------------------
!    * read_pipe()        The character string is read from the pipe.
!                         (Can be overloaded to read ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine blocks until data is available.
!
subroutine read_pipe(iunit, chardata)
 integer, intent(in) :: iunit
 character(len=*), intent(out) :: chardata

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'read_pipe', errstring, source, revision, revdate)

chardata = " "
 
end subroutine

!-----------------------------------------------------------------------------
!    * write_pipe()       The character string is written to the pipe.
!                         (Can be overloaded to write ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine writes and returns immediately.
!
subroutine write_pipe(iunit, chardata)
 integer, intent(in) :: iunit
 character(len=*), intent(in) :: chardata

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'write_pipe', errstring, source, revision, revdate)
 
end subroutine


!-----------------------------------------------------------------------------
! general system util wrappers.
!-----------------------------------------------------------------------------
function shell_execute(execute_string, serialize)
 character(len=*), intent(in) :: execute_string
 logical, intent(in), optional :: serialize
 integer :: shell_execute

! Use the system() command to execute a command string.
! Will wait for the command to complete and returns an
! error code unless you end the command with & to put
! it into background.   Function which returns the rc
! of the command, 0 being all is ok.

! on some platforms/mpi implementations, the system() call
! does not seem to be reentrant.  if serialize is set and
! is true, do each call serially.

character(len=255) :: doit

   !print *, "in-string is: ", trim(execute_string)

   write(doit, "(a, 1x, a1)") trim(execute_string), char(0)

   !print *, "about to run: ", trim(doit)
   !print *, "input string length = ", len(trim(doit))

   shell_execute = system(doit)
   print *, "execution returns, rc = ", shell_execute

end function shell_execute

!-----------------------------------------------------------------------------
subroutine sleep_seconds(naplength)
 real(r8), intent(in) :: naplength

! Wrapper for the sleep command.  Argument is a real
! in seconds.  Different systems have different lower
! resolutions for the minimum time it will sleep.
! Subroutine, no return value.

 integer :: sleeptime

 sleeptime = floor(naplength)
 if (sleeptime <= 0) sleeptime = 1

 call sleep(sleeptime)

end subroutine sleep_seconds

!-----------------------------------------------------------------------------
! Collect min and max on task. This is for adaptive_inflate_mod
subroutine reduce_min_max(minmax, task, global_val)

real(r8), intent(in)  :: minmax(2) ! min max on each task
integer,  intent(in)  :: task ! task to collect on
real(r8), intent(out) :: global_val(2) ! only concerned with this on task collecting result

integer :: errcode

global_val(:) = minmax(:) ! only one task.

end subroutine reduce_min_max

!-----------------------------------------------------------------------------
! One sided communication

subroutine get_from_mean(owner, window, index, x)

integer,  intent(in)  :: owner  ! task in the window that owns the memory
integer,  intent(in)  :: window ! window object
integer,  intent(in)  :: index  ! index in the tasks memory
real(r8), intent(out) :: x ! result

call error_handler(E_ERR,'get_from_mean', 'should not call this code', source, revision, revdate)

end subroutine get_from_mean

!-----------------------------------------------------------------------------

subroutine get_from_fwd(owner, window, index, num_rows, x)

integer, intent(in)  :: owner    ! task in the window that owns the memory
integer,  intent(in)  :: window   ! window object
integer,  intent(in)  :: index    ! index in the tasks memory
integer,  intent(in)  :: num_rows ! number of rows in the window
real(r8), intent(out) :: x(:)     ! result

call error_handler(E_ERR,'get_from_fwd', 'should not call this code', source, revision, revdate)

end subroutine get_from_fwd


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

end module mpi_utilities_mod

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! NOTE -- non-module code, so this subroutine can be called from the
!  utilities module, which this module uses (and cannot have circular refs)
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

subroutine exit_all(exit_code)
 integer, intent(in) :: exit_code

! Call exit with the specified code.

   call exit(exit_code)

end subroutine exit_all

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
