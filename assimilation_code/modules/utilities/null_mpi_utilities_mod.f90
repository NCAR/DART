! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> A collection of interfaces that bypass calling MPI and
!> allows programs to be compiled in a serial configuration.
!> Uses only a single task.  Does NOT require actual MPI libs.

module mpi_utilities_mod

use types_mod, only        : i8, r8, digits12
use utilities_mod, only    : register_module, error_handler,             &
                             initialize_utilities, get_unit, close_file, & 
                             E_ERR, E_WARN, E_MSG, E_DBG, finalize_utilities
use time_manager_mod, only : time_type, set_time


!#ifdef __NAG__
 !use F90_unix_proc, only : sleep, system, exit
 !! block for NAG compiler
 !  PURE SUBROUTINE SLEEP(SECONDS,SECLEFT)
 !    INTEGER,INTENT(IN) :: SECONDS
 !    INTEGER,OPTIONAL,INTENT(OUT) :: SECLEFT
 !
 !  SUBROUTINE SYSTEM(STRING,STATUS,ERRNO)
 !    CHARACTER*(*),INTENT(IN) :: STRING
 !    INTEGER,OPTIONAL,INTENT(OUT) :: STATUS,ERRNO
 !
 !!also used in exit_all outside this module
 !  SUBROUTINE EXIT(STATUS)
 !    INTEGER,OPTIONAL :: STATUS
 !! end block
!#endif


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


interface sum_across_tasks
   module procedure sum_across_tasks_int4
   module procedure sum_across_tasks_int8
   module procedure sum_across_tasks_real
end interface
!   ---- private data for mpi_utilities ----

integer :: myrank          ! my mpi number
integer :: total_tasks     ! total mpi tasks/procs
integer :: comm_size       ! if ens count < tasks, only the first N participate
integer :: datasize        ! should be an accessor function, not a public

public :: initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, send_minmax_to, datasize,                        &
          get_from_fwd, get_from_mean, broadcast_minmax, broadcast_flag,     &
          start_mpi_timer, read_mpi_timer, send_sum_to,                      &
          all_reduce_min_max   ! deprecated, replace with broadcast_minmax

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
   datasize = 4
   write(errstring, *) "Using real * 4 for datasize of r8"
   call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring,source,revision,revdate)
else
   datasize = 8
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
real(r8), intent(in)  :: srcarray(:)
integer,  intent(in)  :: root
real(r8), intent(out) :: dstarray(:)
integer,  intent(out) :: dstcount
integer,  intent(in)  :: how
integer,  intent(out) :: which(:)

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
subroutine sum_across_tasks_int4(addend, sum)
 integer, intent(in) :: addend
 integer, intent(out) :: sum

sum = addend

end subroutine sum_across_tasks_int4

!-----------------------------------------------------------------------------
subroutine sum_across_tasks_int8(addend, sum)
 integer(i8), intent(in) :: addend
 integer(i8), intent(out) :: sum

sum = addend

end subroutine sum_across_tasks_int8

!-----------------------------------------------------------------------------
subroutine sum_across_tasks_real(addend, sum)
 real(r8), intent(in) :: addend
 real(r8), intent(out) :: sum

sum = addend

end subroutine sum_across_tasks_real


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

!> start a time block.  call with different argument
!> to start multiple or nested timers.  same argument
!> must be supplied to read_timer function to get
!> elapsed time since that timer was set.  contrast this with
!> 'start_timer/read_timer' in the utils module which returns
!> elapsed seconds.  this returns whatever units the mpi wtime()
!> function returns.
!>
!> usage:
!>  real(digits12) :: base, time_elapsed
!>
!>  call start_mpi_timer(base)
!>  time_elapsed = read_mpi_timer(base)

subroutine start_mpi_timer(base)

real(digits12), intent(out) :: base

integer :: temp

call system_clock(temp)
base = real(temp, digits12)

end subroutine start_mpi_timer

!-----------------------------------------------------------------------------

!> return the time since the last call to start_timer().
!> can call multiple times to get running times.
!> call with a different base for nested timers.

function read_mpi_timer(base)

real(digits12), intent(in) :: base
real(digits12) :: read_mpi_timer

real(digits12) :: now
integer :: temp

call system_clock(temp)
now = real(temp, digits12)

read_mpi_timer = now - base

end function read_mpi_timer

!-----------------------------------------------------------------------------

function get_dart_mpi_comm()
 integer :: get_dart_mpi_comm

! return dummy value
get_dart_mpi_comm = 0

end function get_dart_mpi_comm

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Collect sum across tasks for a given array.
subroutine send_sum_to(local_val, task, global_val)

real(r8), intent(in)  :: local_val(:) !> min max on each task
integer,  intent(in)  :: task !> task to collect on
real(r8), intent(out) :: global_val(:) !> only concerned with this on task collecting result

integer :: errcode

! collect values on a single given task 
global_val(:) = local_val(:) ! only one task.

end subroutine send_sum_to

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Collect min and max on task. This is for adaptive_inflate_mod
subroutine send_minmax_to(minmax, task, global_val)

real(r8), intent(in)  :: minmax(2) ! min max on each task
integer,  intent(in)  :: task ! task to collect on
real(r8), intent(out) :: global_val(2) ! only concerned with this on task collecting result

global_val(:) = minmax(:) ! only one task.

end subroutine send_minmax_to

!-----------------------------------------------------------------------------
! cover routine which is deprecated.  when all user code replaces this
! with broadcast_minmax(), remove this.
subroutine all_reduce_min_max(min_var, max_var, num_elements)

integer,  intent(in)    :: num_elements
real(r8), intent(inout) :: min_var(num_elements)
real(r8), intent(inout) :: max_var(num_elements)

call broadcast_minmax(min_var, max_var, num_elements)

end subroutine all_reduce_min_max


!-----------------------------------------------------------------------------
! Find min and max of each element of an array across tasks, put the result on every task.
! For this null_mpi_version min_var and max_var are unchanged because there is
! only 1 task.
subroutine broadcast_minmax(min_var, max_var, num_elements)

integer,  intent(in)    :: num_elements
real(r8), intent(inout) :: min_var(num_elements)
real(r8), intent(inout) :: max_var(num_elements)

end subroutine broadcast_minmax

!-----------------------------------------------------------------------------
! One sided communication

subroutine get_from_mean(owner, window, mindex, x)

integer,  intent(in)  :: owner  ! task in the window that owns the memory
integer,  intent(in)  :: window ! window object
integer,  intent(in)  :: mindex ! index in the tasks memory
real(r8), intent(out) :: x ! result

call error_handler(E_ERR,'get_from_mean', 'cannot be used in serial mode', source, revision, revdate)

! NOT REACHED
x = 0.0_r8

end subroutine get_from_mean

!-----------------------------------------------------------------------------

subroutine get_from_fwd(owner, window, mindex, num_rows, x)

integer,  intent(in)  :: owner    ! task in the window that owns the memory
integer,  intent(in)  :: window   ! window object
integer,  intent(in)  :: mindex   ! index in the tasks memory
integer,  intent(in)  :: num_rows ! number of rows in the window
real(r8), intent(out) :: x(:)     ! result

call error_handler(E_ERR,'get_from_fwd', 'cannot be used in serial mode', source, revision, revdate)

! NOT REACHED
x(:) = 0.0_r8

end subroutine get_from_fwd


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

! Broadcast logical
subroutine broadcast_flag(flag, root)

logical, intent(inout) :: flag
integer, intent(in)    :: root ! relative to get_dart_mpi_comm()

end subroutine broadcast_flag


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
