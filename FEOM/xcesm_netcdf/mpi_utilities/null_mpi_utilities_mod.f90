! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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

use        types_mod, only : r8, digits12
use    utilities_mod, only : register_module, error_handler,          & 
                             E_ERR, E_WARN, E_MSG,                    &
                             initialize_utilities, finalize_utilities
use time_manager_mod, only : time_type

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


!   ---- private data for mpi_utilities ----

integer :: myrank        = -1  ! my mpi number
integer :: total_tasks   = -1  ! total mpi tasks/procs
integer :: my_local_comm =  0  ! duplicate communicator private to this file

public :: initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, get_dart_mpi_comm

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical :: module_initialized = .false.

character(len = 256) :: saved_progname = ''
character(len = 256) :: shell_name = ''   ! if needed, add ksh, tcsh, bash, etc

character(len = 256) :: errstring

! Namelist input - placeholder for now; no options yet in this module.
!namelist /mpi_utilities_nml/ x


contains

!-----------------------------------------------------------------------------
! mpi cover routines
!-----------------------------------------------------------------------------

!> Initialize the utilities module, and print out a message including the 
!> program name.

subroutine initialize_mpi_utilities(progname, alternatename)

character(len=*), intent(in), optional :: progname
character(len=*), intent(in), optional :: alternatename

if ( module_initialized ) then
   ! return without calling the code below multiple times.  Print out a warning each
   ! time this is called again because it may indicate an error in logic.  In a well-
   ! constructed program the initialize routine will only be called once.
   write(errstring, *) 'initialize_mpi_utilities has already been called'
   call error_handler(E_WARN,'initialize_mpi_utilities', errstring, source, revision, revdate)
   return
endif

call initialize_utilities(progname, alternatename)

if (present(progname)) then
   if (len_trim(progname) <= len(saved_progname)) then
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
my_local_comm = 0

if (r8 /= digits12) then
      write(errstring, *) "Using real * 4 for datasize of r8"
      call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring,source,revision,revdate)
   endif

! non-MPI successfully initialized.
call error_handler(E_MSG,'initialize_mpi_utilities: ','Running single process', &
                   source, revision, revdate)

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

integer,          intent(in)           :: dest_id
real(r8),         intent(in)           :: srcarray(:)
type(time_type),  intent(in), optional :: time
character(len=*), intent(in), optional :: label

write(errstring, '(a)') "cannot call send_to() in the single process case"
      call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)

end subroutine send_to


!-----------------------------------------------------------------------------

!> Receive the dstarray from the source task id.
!> This communication style cannot be easily simulated correctly with one task.
!> If called, always throw an error.

subroutine receive_from(src_id, destarray, time, label)

integer,          intent(in)            :: src_id
real(r8),         intent(out)           :: destarray(:)
type(time_type),  intent(out), optional :: time
character(len=*), intent(in),  optional :: label

write(errstring, '(a)') "cannot call receive_from() in the single process case"
      call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)

end subroutine receive_from



!-----------------------------------------------------------------------------

!> The array already has the values, nothing to do.  Not an error to call.

subroutine array_broadcast(array, root)

real(r8), intent(inout) :: array(:)
integer,  intent(in)    :: root

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
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

integer,  intent(in)              :: from
! arrays are really only intent(in) here, but must match array_broadcast() call.
real(r8), intent(inout)           :: array1(:)
real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if (from /= myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "must be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! Data is already in arrays, so you can return here.

end subroutine broadcast_send

!-----------------------------------------------------------------------------

!> Returns with nothing to do.  Does validate the 'from' task id.
!> Not an error to call.

subroutine broadcast_recv(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)

integer,  intent(in)              :: from
! arrays are really only intent(out) here, but must match array_broadcast() call.
real(r8), intent(inout)           :: array1(:)
real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

if ( .not. module_initialized ) call initialize_mpi_utilities()

! simple idiotproofing
if (from == myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "cannot be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! Data is already in arrays, so you can return here.

end subroutine broadcast_recv

!-----------------------------------------------------------------------------
   
!> Returns addend in sum for a single task.

subroutine sum_across_tasks(addend, sum)

integer, intent(in) :: addend
integer, intent(out) :: sum

if ( .not. module_initialized ) call initialize_mpi_utilities()

sum = addend

end subroutine sum_across_tasks


!-----------------------------------------------------------------------------
! pipe-related utilities - satisfy the subroutine names, but nothing to do.
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

character(len=*), intent(in)           :: execute_string
logical,          intent(in), optional :: serialize
integer                                :: shell_execute

character(len=255) :: doit

!DEBUG: print *, "in-string is: ", trim(execute_string)

   call do_system(shell_name, execute_string, shell_execute)
       
!DEBUG: print *, "execution returns, rc = ", shell_execute

end function shell_execute

!-----------------------------------------------------------------------------

!> wrapper so you only have to make this work in a single place

subroutine do_system(shell, execute, rc)

character(len=*), intent(in)  :: shell
character(len=*), intent(in)  :: execute
integer,          intent(out) :: rc

!#ifdef __NAG__
!  call system(trim(shell)//' '//trim(execute)//' '//char(0), errno=rc)
!#else
   rc = system(trim(shell)//' '//trim(execute)//' '//char(0))
!#endif

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

!> Return 0 for the communicator number.

function get_dart_mpi_comm()

integer :: get_dart_mpi_comm

get_dart_mpi_comm = my_local_comm

end function get_dart_mpi_comm

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

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

integer, intent(in) :: exit_code

call exit(exit_code)

end subroutine exit_all

!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
