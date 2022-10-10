! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> See the mpi_utilities_mod.f90 documentation for more information on
!> this file.  If you change either file you must make the corresponding
!> changes in the other if it affects a public interface.  They must stay
!> in sync.
!>
!> Substitute this code for mpi_utilities_mod.f90 if you do not want to
!> have to link in an MPI library, and you only want to run single task.
!> Many of the single task DART utility programs use this file instead of
!> the parallel version.  Note that this file has the same module name
!> and the same external entry points as the real mpi_utilities_mod.f90
!> file, so it will link correctly as a replacement file.
!>
!> Programs using this module instead of the actual MPI routines do not
!> need to be compiled with the MPI wrapper commands (e.g. mpif90).
!> In most cases it will be better to compile with the fortran compiler
!> directly.  (On some platforms this is required.)

module mpi_utilities_mod

use        types_mod, only : i8, r8, digits12
use    utilities_mod, only : error_handler, E_ERR, E_WARN, E_MSG, &
                             initialize_utilities, finalize_utilities
use time_manager_mod, only : time_type, set_time

! We build on case-insensitive systems so we cannot reliably
! count on having the build system run the fortran preprocessor
! since the usual distinction is between bob.F90 and bob.f90
! to decide what needs preprocessing.  instead we utilize a
! script we provide called 'fixsystem' which looks for the
! special XXX_BLOCK_EDIT comment lines and comments the blocks
! in and out depending on the target compiler.

! the NAG compiler needs these special definitions enabled.
! the #ifdef lines are only there in case someday we can use
! the fortran preprocessor.  they need to stay commented out.

! !!NAG_BLOCK_EDIT START COMMENTED_OUT
! !#ifdef __NAG__
!
! use F90_unix_proc, only : sleep, system, exit
!
! !! NAG only needs the use statement above, but
! !! these are the calling sequences if you need
! !! to use these routines additional places in code.
! !  PURE SUBROUTINE SLEEP(SECONDS,SECLEFT)
! !    INTEGER,INTENT(IN) :: SECONDS
! !    INTEGER,OPTIONAL,INTENT(OUT) :: SECLEFT
! !
! !  SUBROUTINE SYSTEM(STRING,STATUS,ERRNO)
! !    CHARACTER*(*),INTENT(IN) :: STRING
! !    INTEGER,OPTIONAL,INTENT(OUT) :: STATUS,ERRNO
! !
! !!also used in exit_all outside this module
! !  SUBROUTINE EXIT(STATUS)
! !    INTEGER,OPTIONAL :: STATUS
! !! end block
!
!  !#endif
! !!NAG_BLOCK_EDIT END COMMENTED_OUT


implicit none
private


! BUILD TIP 
! Some compilers require an interface block for the system() function;
! some fail if you define one.  If you get an error at link time (something
! like 'undefined symbol _system_') try running the fixsystem script in
! this directory.  It is a sed script that comments in and out the interface
! block below.  Please leave the BLOCK comment lines unchanged.

 !!SYSTEM_BLOCK_EDIT START COMMENTED_IN
 !#if .not. defined (__GFORTRAN__) .and. .not. defined(__NAG__)
 ! interface block for getting return code back from system() routine
 interface
  function system(string)
   character(len=*) :: string
   integer :: system
  end function system
 end interface
 ! end block
 !#endif
 !!SYSTEM_BLOCK_EDIT END COMMENTED_IN


! allow global sum to be computed for integers, r4, and r8s
interface sum_across_tasks
   module procedure sum_across_tasks_int4
   module procedure sum_across_tasks_int8
   module procedure sum_across_tasks_real
end interface

!   ---- private data for mpi_utilities ----

integer :: myrank        = 0  ! my mpi number
integer :: total_tasks   = 1  ! total mpi tasks/procs
integer :: my_local_comm = 0  ! duplicate communicator private to this file
integer :: datasize      = 8  ! which MPI type corresponds to our r8 definition



public :: initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, get_dart_mpi_comm, datasize, send_minmax_to,     &
          get_from_fwd, get_from_mean, broadcast_minmax, broadcast_flag,     &
          start_mpi_timer, read_mpi_timer, send_sum_to, get_global_max,      &
          all_reduce_min_max  ! deprecated, replace by broadcast_minmax

character(len=*), parameter :: source = 'null_mpi_utilities_mod.f90'

logical :: module_initialized   = .false.

character(len = 256) :: saved_progname = ''
character(len = 128) :: shell_name = ''   ! if needed, add ksh, tcsh, bash, etc

character(len = 256) :: errstring

! Namelist input - placeholder for now; no options yet in this module.
!namelist /mpi_utilities_nml/ x


contains

!-----------------------------------------------------------------------------
! mpi cover routines
!-----------------------------------------------------------------------------

!> Initialize the utilities module, and print out a message including the 
!> program name.
subroutine initialize_mpi_utilities(progname, alternatename, communicator)

character(len=*), intent(in), optional :: progname
character(len=*), intent(in), optional :: alternatename
integer,          intent(in), optional :: communicator

if ( module_initialized ) then
   ! return without calling the code below multiple times.  Print out a warning each
   ! time this is called again because it may indicate an error in logic.  In a well-
   ! constructed program the initialize routine will only be called once.
   write(errstring, *) 'initialize_mpi_utilities has already been called'
   call error_handler(E_WARN,'initialize_mpi_utilities', errstring, source)
   return
endif

module_initialized = .true.

! Initialize the module with utilities
call initialize_utilities(progname, alternatename)

if (present(progname)) then
   if (len_trim(progname) <= len(saved_progname)) then
      saved_progname = trim(progname)
   else
      saved_progname = progname(1:len(saved_progname))
   endif
endif

! log info if requested

if (r8 /= digits12) then
   datasize = 4
   write(errstring, *) "Using real * 4 for datasize of r8"
   call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring, source)
else
   datasize = 8
endif

! non-MPI successfully initialized.
call error_handler(E_MSG,'initialize_mpi_utilities: ','Running single process', source)

end subroutine initialize_mpi_utilities

!-----------------------------------------------------------------------------

!> Shut down cleanly.  Call normal utilities finalize if we have actually
!> ever called initialize.  Otherwise there is nothing to do in the null case.

subroutine finalize_mpi_utilities(callfinalize, async)
 logical, intent(in), optional :: callfinalize
 integer, intent(in), optional :: async

if ( .not. module_initialized ) return

if (saved_progname /= '') then
   call finalize_utilities(saved_progname)
else
   call finalize_utilities()
endif

end subroutine finalize_mpi_utilities


!-----------------------------------------------------------------------------

!> Return the number of MPI tasks.  For this code this is always 1.

function task_count()

integer :: task_count

if ( .not. module_initialized ) call initialize_mpi_utilities()

task_count = total_tasks

end function task_count


!-----------------------------------------------------------------------------

!> Return my unique task id.  For this code this is always 0.

function my_task_id()

integer :: my_task_id

if ( .not. module_initialized ) call initialize_mpi_utilities()

my_task_id = myrank

end function my_task_id


!-----------------------------------------------------------------------------

!> A no-op for this code.

subroutine task_sync()

end subroutine task_sync


!-----------------------------------------------------------------------------

!> Send the srcarray to the destination task id.
!> This communication style cannot be easily simulated correctly with one task.
!> If called, always throw an error.

subroutine send_to(dest_id, srcarray, time, label)
 integer, intent(in) :: dest_id
 real(r8), intent(in) :: srcarray(:)
 type(time_type), intent(in), optional :: time
 character(len=*), intent(in), optional :: label

write(errstring, '(a)') "cannot call send_to() in the single process case"
      call error_handler(E_ERR,'send_to', errstring, source)

end subroutine send_to


!-----------------------------------------------------------------------------

!> Receive the dstarray from the source task id.
!> This communication style cannot be easily simulated correctly with one task.
!> If called, always throw an error.

subroutine receive_from(src_id, destarray, time, label)
 integer, intent(in) :: src_id
 real(r8), intent(inout) :: destarray(:)    ! really only out, but avoid compiler warnings
 type(time_type), intent(out), optional :: time
 character(len=*), intent(in), optional :: label

write(errstring, '(a)') "cannot call receive_from() in the single process case"
      call error_handler(E_ERR,'receive_from', errstring, source)

end subroutine receive_from



!-----------------------------------------------------------------------------

!> The array already has the values, nothing to do.  Not an error to call.

subroutine array_broadcast(array, root, icount)
 real(r8), intent(inout) :: array(:)
 integer, intent(in) :: root
 integer, intent(in), optional :: icount

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source)
endif

! Data is already in array, so you can return here.

end subroutine array_broadcast


!-----------------------------------------------------------------------------
! DART-specific cover utilities
!-----------------------------------------------------------------------------

!> Return .TRUE. if my local task id is 0, .FALSE. otherwise.
!> (Task numbers in MPI start at 0, contrary to the rules of polite fortran.)
!> This version always returns .TRUE. since there is only a single task ever.

function iam_task0()

logical :: iam_task0

if ( .not. module_initialized ) call initialize_mpi_utilities()

iam_task0 = (myrank == 0)

end function iam_task0

!-----------------------------------------------------------------------------

!> Returns with nothing to do.  Does validate the 'from' task id.
!> Not an error to call.

subroutine broadcast_send(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)
 integer, intent(in) :: from
! arrays are really only intent(in) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:)
 real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if (from /= myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "must be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_send', errstring, source)
endif

! this does nothing, because the array already has the data.
! the subroutine does validate 'from' to be sure it's a valid task id.
call array_broadcast(array1, from)

end subroutine broadcast_send

!-----------------------------------------------------------------------------

!> Returns with nothing to do.  Does validate the 'from' task id.
!> Not an error to call.

subroutine broadcast_recv(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)
 integer, intent(in) :: from
! arrays are really only intent(out) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:)
 real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if (from == myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "cannot be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_recv', errstring, source)
endif

! this does nothing, because the array already has the data.
! the subroutine does validate 'from' to be sure it's a valid task id.
call array_broadcast(array1, from)

end subroutine broadcast_recv

!-----------------------------------------------------------------------------
!> return sum for various input types/kinds

subroutine sum_across_tasks_int4(addend, sum)
 integer, intent(in) :: addend
 integer, intent(out) :: sum

sum = addend

end subroutine sum_across_tasks_int4

!-----------------------------------------------------------------------------

subroutine sum_across_tasks_int8(addend, sum)
 integer(i8), intent(in)  :: addend
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

!> Sum array items across all tasks and send
!> results in an array of same size to one task.

subroutine send_sum_to(local_val, task, global_val)

real(r8), intent(in)  :: local_val(:)  !! addend vals on each task
integer,  intent(in)  :: task          !! task to collect on
real(r8), intent(out) :: global_val(:) !! results returned only on given task

global_val(:) = local_val(:) ! only one task.

end subroutine send_sum_to

!-----------------------------------------------------------------------------

!> Collect min and max on task.

subroutine send_minmax_to(minmax, task, global_val)

real(r8), intent(in)  :: minmax(2)     !! min max on each task
integer,  intent(in)  :: task          !! task to collect on
real(r8), intent(out) :: global_val(2) !! results returned only on given task

global_val(:) = minmax(:) ! only one task.

end subroutine send_minmax_to

!-----------------------------------------------------------------------------

!> cover routine which is deprecated.  when all user code replaces this
!> with broadcast_minmax(), remove this.

subroutine all_reduce_min_max(min_var, max_var, num_elements)

integer,  intent(in)    :: num_elements
real(r8), intent(inout) :: min_var(num_elements)
real(r8), intent(inout) :: max_var(num_elements)

call broadcast_minmax(min_var, max_var, num_elements)

end subroutine all_reduce_min_max

!-----------------------------------------------------------------------------

!> Find min and max of each element of an array across tasks, put the result on every task.
!> For this null_mpi_version min_var and max_var are unchanged because there is
!> only 1 task.

subroutine broadcast_minmax(min_var, max_var, num_elements)

integer,  intent(in)    :: num_elements
real(r8), intent(inout) :: min_var(num_elements)
real(r8), intent(inout) :: max_var(num_elements)

end subroutine broadcast_minmax

!-----------------------------------------------------------------------------

!> Broadcast logical

subroutine broadcast_flag(flag, root)

logical, intent(inout) :: flag
integer, intent(in)    :: root !! relative to get_dart_mpi_comm()

! does nothing because data is already there

end subroutine broadcast_flag


!-----------------------------------------------------------------------------
! pipe-related utilities
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

subroutine block_task()

end subroutine block_task

!-----------------------------------------------------------------------------

subroutine restart_task()

end subroutine restart_task


!-----------------------------------------------------------------------------
! general system util wrappers.
!-----------------------------------------------------------------------------

!> Use the system() command to execute a command string.
!> Will wait for the command to complete and returns an
!> error code unless you end the command with & to put
!> it into background.   Function which returns the rc
!> of the command, 0 being all is ok.

function shell_execute(execute_string, serialize)
 character(len=*), intent(in) :: execute_string
 logical, intent(in), optional :: serialize
 integer :: shell_execute

!DEBUG: print *, "in-string is: ", trim(execute_string)

call do_system(execute_string, shell_execute)

!DEBUG: print *, "execution returns, rc = ", shell_execute
    
end function shell_execute

!-----------------------------------------------------------------------------

!> wrapper so you only have to make this work in a single place
!> 'shell_name' is a namelist item and normally is the null string.
!> on at least on cray system, the compute nodes only had one type
!> of shell and you had to specify it.

subroutine do_system(execute, rc)

character(len=*), intent(in)  :: execute
integer,          intent(out) :: rc

! !!NAG_BLOCK_EDIT START COMMENTED_OUT
!  call system(trim(shell_name)//' '//trim(execute)//' '//char(0), errno=rc)
! !!NAG_BLOCK_EDIT END COMMENTED_OUT
! !!OTHER_BLOCK_EDIT START COMMENTED_IN
    rc = system(trim(shell_name)//' '//trim(execute)//' '//char(0))
! !!OTHER_BLOCK_EDIT END COMMENTED_IN

end subroutine do_system

!-----------------------------------------------------------------------------

!> Wrapper for the sleep command.  Argument is a real
!> in seconds.  Different systems have different lower
!> resolutions for the minimum time it will sleep.
!> Subroutine, no return value.

subroutine sleep_seconds(naplength)

real(r8), intent(in) :: naplength

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
!> elapsed seconds.  this returns whatever units the fortran
!> intrinsic system_clock() returns.
!>
!> usage:
!>  real(digits12) :: base, time_elapsed
!>
!>  call start_mpi_timer(base)
!>  time_elapsed = read_mpi_timer(base)

subroutine start_mpi_timer(base)

real(digits12), intent(out) :: base

integer(i8) :: temp

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
integer(i8) :: temp

call system_clock(temp)
now = real(temp, digits12)

read_mpi_timer = now - base

end function read_mpi_timer

!-----------------------------------------------------------------------------
!> Return the communicator number.

function get_dart_mpi_comm()

integer :: get_dart_mpi_comm

get_dart_mpi_comm = my_local_comm

end function get_dart_mpi_comm

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! One sided communication

subroutine get_from_mean(owner, window, mindex, x)

integer,  intent(in)  :: owner  ! task in the window that owns the memory
integer,  intent(in)  :: window ! window object
integer,  intent(in)  :: mindex ! index in the tasks memory
real(r8), intent(out) :: x      ! result

call error_handler(E_ERR,'get_from_mean', 'cannot be used in serial mode', source)

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

call error_handler(E_ERR,'get_from_fwd', 'cannot be used in serial mode', source)

! NOT REACHED
x(:) = 0.0_r8

end subroutine get_from_fwd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

!> Collect global max values on each task.

subroutine get_global_max(max)

real(r8), intent(inout)  :: max       !> global max over tasks

! Nothing to do with only one task.

end subroutine get_global_max


end module mpi_utilities_mod

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!> NOTE: non-module code, so this subroutine can be called from the
!>  utilities module, which this module uses (and cannot have circular refs)
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!> Call exit with the specified code.  NOT PART of the mpi_utilities_mod, so
!> this can be called from any code in the system.

subroutine exit_all(exit_code)
! !!NAG_BLOCK_EDIT START COMMENTED_OUT
! use F90_unix_proc, only : exit
! !!NAG_BLOCK_EDIT END COMMENTED_OUT
 integer, intent(in) :: exit_code

   call exit(exit_code)

end subroutine exit_all
