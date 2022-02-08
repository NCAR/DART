! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


!> A collection of interfaces to the MPI (Message Passing Interface)
!> multi-processor communication library routines.
!>
!> This file and its companion file, null_mpi_utilities_mod.f90, 
!> are the only modules in the DART system that should be calling
!> the MPI library routines directly.  This isolates the MPI calls
!> and allows programs to swap in the null version to compile the
!> same source files into a serial program.
!>
!> The names of these routines were intentionally picked to be
!> more descriptive to someone who doesn't the MPI interfaces.
!> e.g. MPI_AllReduce() may not immediately tell a user what
!> it does, but broadcast_minmax() is hopefully more helpful.
!>
!> If you add any routines or change any arguments in this file
!> you must make the same changes in the null version.  These two
!> modules have the same module name and must have identical
!> public routines and calling formats.
!>
!> All MPI routines are called from here.  There is a companion
!> file which has the same module name and entry points but all
!> routines are stubs.  This allows a single-task version of the
!> code to be compiled without the MPI libraries.


module mpi_utilities_mod

use types_mod, only :  i4, i8, r8, digits12
use utilities_mod, only : error_handler, & 
                          E_ERR, E_WARN, E_MSG, E_DBG, get_unit, close_file, &
                          set_output, set_tasknum, initialize_utilities,     &
                          finalize_utilities,                                &
                          nmlfileunit, do_nml_file, do_nml_term,             &
                          find_namelist_in_file, check_namelist_read

use time_manager_mod, only : time_type, get_time, set_time

! BUILD TIP 1
! Many MPI installations have an MPI module; if one is present, use it.
! ('use mpi')
! If not, there will be an MPI include file which defines the parameters.
! ('include mpif.h')
! Use one but not both.   The 'use' line must be before the 'implicit none' 
! and 'private' lines, 'include' must come after.  Go figure.
! For more help on compiling a module which uses MPI see the 
! $DART/developer_tests/mpi_utilities/tests/README

use mpi


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

!include "mpif.h"


! BUILD TIP 2
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

integer :: myrank        = -1  ! my mpi number
integer :: total_tasks   = -1  ! total mpi tasks/procs
integer :: my_local_comm =  0  ! duplicate communicator private to this file
integer :: datasize      =  0  ! which MPI type corresponds to our r8 definition
integer :: longinttype   =  0  ! create an MPI type corresponding to our i8 definition



public :: initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, get_dart_mpi_comm, datasize, send_minmax_to,     &
          get_from_fwd, get_from_mean, broadcast_minmax, broadcast_flag,     &
          start_mpi_timer, read_mpi_timer, send_sum_to, get_global_max,      &
          all_reduce_min_max  ! deprecated, replace by broadcast_minmax

character(len=*), parameter :: source = 'mpi_utilities_mod.f90'

logical :: module_initialized   = .false.

character(len = 256) :: saved_progname = ''
character(len = 128) :: shell_name = ''   ! if needed, add ksh, tcsh, bash, etc
integer :: head_task = 0         ! def 0, but N-1 if reverse_task_layout true
logical :: print4status = .true. ! minimal messages for async4 handshake

logical :: given_communicator = .false.   ! if communicator passed in, use it

character(len = 256) :: errstring, errstring1

! for broadcasts, pack small messages into larger ones.  remember that the
! byte size will be this count * 8 because we only communicate r8s.  (unless
! the code is compiled with r8 redefined as r4, in which case it's * 4).
integer, parameter :: PACKLIMIT = 512

! also for broadcasts, make sure message size is not too large.  if so,
! split a single request into two or more broadcasts.  i know 2G is really
! 2 * 1024 * 1024 * 1024, but err on the conservative side here.
integer, parameter :: BCAST_MAXSIZE = 2 * 1000 * 1000 * 1000

! option for simple send/recvs to limit max single message size.
! split a single request into two or more broadcasts.  i know 2G is really
! 2 * 1024 * 1024 * 1024, but err on the conservative side here.
integer, parameter :: SNDRCV_MAXSIZE = 2 * 1000 * 1000 * 1000

! this turns on trace messages for most MPI communications
logical :: verbose        = .false.   ! very very very verbose, use with care
logical :: async2_verbose = .false.   ! messages only for system() in async2
logical :: async4_verbose = .false.   ! messages only for block/restart async4

! if your batch system does the task layout backwards, set this to true
! so the last task will communicate with the script in async 4 mode.
! as of now, mpich and mvapich do it forward, openmpi does it backwards.
logical :: reverse_task_layout  = .false.   ! task 0 on head node; task N-1 if .true.
logical :: separate_node_sync   = .false.   ! true if tasks & script do not share nodes
logical :: create_local_comm    = .true.    ! make a private communicator

! for large numbers of MPI tasks, you will get replicated messages, one
! per task, if this is set to true.  however, for debugging if you need
! messages from tasks which aren't 0, this will elicit them.  error messages
! from any task will print regardless of this setting.
logical :: all_tasks_print      = .false.   ! by default only msgs from 0 print

! make local array copy for send/recv/bcast.  was needed on an old, buggy version
! of the mpi libs but seems unneeded now. 
logical :: make_copy_before_sendrecv  = .false.   ! should not be needed; .true. is very slow
logical :: make_copy_before_broadcast = .false.   ! should not be needed; .true. is very slow

! NAMELIST: change the following from .false. to .true. to enable
! the reading of this namelist.  This is the only place you need
! to make this change.
logical :: read_namelist = .false.

namelist /mpi_utilities_nml/ reverse_task_layout, all_tasks_print, &
                             verbose, async2_verbose, async4_verbose, &
                             shell_name, separate_node_sync, create_local_comm, &
                             make_copy_before_sendrecv, make_copy_before_broadcast 

contains

!-----------------------------------------------------------------------------
! mpi cover routines
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

!> Initialize MPI and query it for global information.  Make a duplicate
!> communicator so that any user code which wants to call MPI will not 
!> interfere with any outstanding asynchronous requests, accidental tag
!> matches, etc.  This routine must be called before any other routine in
!> this file, and it should not be called more than once (but it does have
!> defensive code in case that happens.)

subroutine initialize_mpi_utilities(progname, alternatename, communicator)

character(len=*), intent(in), optional :: progname
character(len=*), intent(in), optional :: alternatename
integer,          intent(in), optional :: communicator

integer :: errcode, iunit
logical :: already

if ( module_initialized ) then
   ! return without calling the code below multiple times
   write(errstring, *) 'initialize_mpi_utilities has already been called'
   call error_handler(E_WARN,'initialize_mpi_utilities', errstring, source)
   return
endif

! prevent any other code from calling into this init
! routine and causing overlapping code execution
module_initialized = .true.

! some implementations of mpich need this to happen before any I/O is done.
! this makes the startup sequence very tricky. 
! still, allow for the possibility that the user code has already initialized
! the MPI libs and do not try to call initialize again.
errcode = -999 
call MPI_Initialized(already, errcode)
if (errcode /= MPI_SUCCESS) then
   write(*, *) 'MPI_Initialized returned error code ', errcode
   call exit(-99)
endif
if (.not.already) then
   call MPI_Init(errcode)
   if (errcode /= MPI_SUCCESS) then
      write(*, *) 'MPI_Init returned error code ', errcode
      call exit(-99)
   endif
endif

if (.not. present(communicator)) then
   ! give this a temporary initial value, in case we call the abort code.
   ! later we will dup the world comm and use a private comm for our comms.
   my_local_comm = MPI_COMM_WORLD
else
   my_local_comm = communicator
   given_communicator = .true.
endif

call MPI_Comm_rank(my_local_comm, myrank, errcode)
if (errcode /= MPI_SUCCESS) then
   write(*, *) 'MPI_Comm_rank returned error code ', errcode
   call exit(-99)
endif

! pass the arguments through so the utilities can log the program name
! only PE0 gets to output, whenever possible.
if (myrank == 0) then
   call initialize_utilities(progname, alternatename, .true.)
else
   call initialize_utilities(progname, alternatename, .false.)
endif

! save a copy of the initial program name for use in finalize.
! make sure it is not too long to fit into our variable.
if (present(progname)) then
   if (len_trim(progname) <= len(saved_progname)) then
      saved_progname = trim(progname)
   else
      saved_progname = progname(1:len(saved_progname))
   endif
endif

! if logging, add this info to the log
! (must come after regular utils are initialized)

! this must come AFTER the standard utils are initialized.
! Read the DART namelist for the mpi_utilities.
if (read_namelist) then
   call find_namelist_in_file('input.nml', 'mpi_utilities_nml', iunit)
   read(iunit, nml = mpi_utilities_nml, iostat = errcode)
   call check_namelist_read(iunit, errcode, "mpi_utilities_nml")
else
   errstring = ' !must edit assimilation_code/modules/utilities/mpi_utilities_mod.f90 to enable this namelist'
   if (do_nml_file()) write(nmlfileunit, '(A)') trim(errstring)
   if (do_nml_term()) write(     *     , '(A)') trim(errstring)
endif

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=mpi_utilities_nml)
if (do_nml_term()) write(     *     , nml=mpi_utilities_nml)


! duplicate the world communicator to isolate us from any other user
! calls to MPI.  All subsequent mpi calls here will use the local communicator
! and not the global world comm.
if (.not. given_communicator .and. create_local_comm) then
   call MPI_Comm_dup(MPI_COMM_WORLD, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Comm_dup returned error code ', errcode
      call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source)
   endif
endif

! find out who we are (0 to N-1).
call MPI_Comm_rank(my_local_comm, myrank, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_rank returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source)
endif

! number of tasks (if 10, returns 10.  task id numbers go from 0-9)
call MPI_Comm_size(my_local_comm, total_tasks, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_size returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source)
endif

! tell the utilities module what task number we are.
call set_tasknum(myrank)

! Turn off non-critical log messages from all but task 0, for performance
! and for sanity (e.g. only one copy of informational messages).  can be
! overridden by namelist if desired.
if (all_tasks_print) then
   call set_output(.true.)                     ! everyone gets to talk
else
   if (myrank /= 0) call set_output(.false.)   ! process 0 speaks for all
endif

! Users have the option of redefining the DART r8 kind to be the same size
! as r4.  But when we call the MPI routines we have to know what MPI type
! to use.  The internal dart kind 'digits12' is always defined to be *8, so
! if the parameters r8 and digits12 are the same value, use *8.  Otherwise
! assume r4 and r8 were defined to be the same and use *4.  (Some users
! do this to cut down on the space requirements.)

if (r8 == digits12) then
   datasize = MPI_REAL8
   !print *, "using real * 8 for datasize of r8"
else
   datasize = MPI_REAL4
   if (myrank == 0) then
      write(errstring, *) "Using real * 4 for datasize of r8"
      call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring, source)
   endif
endif

! create a type we can use for integer(i8) calls
longinttype = MPI_INTEGER8
! or:
!call MPI_Type_Create_F90_Integer(15, longinttype, errcode)
!if (errcode /= MPI_SUCCESS) then
!   write(errstring, '(a,i8)') 'MPI_Type_Create_F90_Integer returned error code ', errcode
!   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source)
!endif


! in async 4 mode, where the controlling job (usually filter) and the 
! model are both mpi tasks and they handshake via named pipes, the tasks
! writing and reading the named pipes must be on the same node as the 
! one running the start script.  many MPI implementations (mpich, mvapich) 
! put task 0 on the first node, with the script.  some (openmpi)
! lay them out in reverse order, and task N-1 is on the same node
! as the script.  if this is wrong, things will hang when the
! filter task tries to advance the model - it won't be able to
! write the 'go ahead' message to the pipe.  set this via namelist
! to match what works on your system.
if (reverse_task_layout) then
   head_task = total_tasks - 1
else
   head_task = 0   ! normal case
endif

! MPI successfully initialized.  Log for the record how many tasks.
if (verbose) write(*,*) "PE", myrank, ": MPI successfully initialized"

if (myrank == 0) then
   write(errstring, *) 'Running with ', total_tasks, ' MPI processes.'
   call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring, source)
endif

end subroutine initialize_mpi_utilities

!-----------------------------------------------------------------------------

!> Shut down MPI cleanly.  This must be done before the program exits; on
!> some implementations of MPI the final I/O flushes are not done until this
!> is called.  The optional argument can prevent us from calling MPI_Finalize,
!> so that user code can continue to use MPI after this returns.  Calling other
!> routines in this file after calling finalize will invalidate your warranty.

subroutine finalize_mpi_utilities(callfinalize, async)
 logical, intent(in), optional :: callfinalize
 integer, intent(in), optional :: async

integer :: errcode
logical :: dofinalize

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source)
endif

! give the async=4 case a chance to tell the script to shut down.
if (present(async)) then
   call finished_task(async)
endif

! close the log files and write a timestamp
if (saved_progname /= '') then
   call finalize_utilities(saved_progname)
else
   call finalize_utilities()
endif

! For the SGI implementation of MPI in particular, sync before shutting
! down MPI.  Once MPI is shut down, no other output can be written; it causes
! tasks to hang if they try.  Make sure all tasks are here and ready to
! close down at the same time.
call task_sync()

if (.not. given_communicator) then
   ! Release the private communicator we created at init time.
   if (my_local_comm /= MPI_COMM_WORLD) then
      call MPI_Comm_free(my_local_comm, errcode)
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Comm_free returned error code ', errcode
         call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source)
      endif
   endif
   my_local_comm = MPI_COMM_WORLD
endif

! If the optional argument is not given, or is given and is true, 
! shut down mpi.  Only if the argument is specified and is false do we
! skip the finalization.
if (.not.present(callfinalize)) then
   dofinalize = .TRUE.
else
   dofinalize = callfinalize
endif

! Normally we shut down MPI here.  If the user tells us not to shut down MPI
! they must call this routine from their own code before exiting.
if (dofinalize) then
   if (verbose) write(*,*) "PE", myrank, ": MPI finalize being called now"
   call MPI_Finalize(errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Finalize returned error code ', errcode
      call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source)
   endif
endif

! NO I/O after calling MPI_Finalize.  on some MPI implementations
! this can hang the job.

end subroutine finalize_mpi_utilities


!-----------------------------------------------------------------------------

!> Return the total number of MPI tasks.  e.g. if the number of tasks is 4,
!> it returns 4.  (The actual task numbers are 0-3.)

function task_count()

integer :: task_count

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_count', errstring, source)
endif

task_count = total_tasks

end function task_count


!-----------------------------------------------------------------------------

!> Return my unique task id.  Values run from 0 to N-1 (where N is the
!> total number of MPI tasks.

function my_task_id()

integer :: my_task_id

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'my_task_id', errstring, source)
endif

my_task_id = myrank

end function my_task_id


!-----------------------------------------------------------------------------

!> Synchronize all tasks.  This subroutine does not return until all tasks
!> execute this line of code.

subroutine task_sync()

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_sync', errstring, source)
endif

if (verbose) write(*,*) "PE", myrank, ": waiting at MPI Barrier"
call MPI_Barrier(my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Barrier returned error code ', errcode
   call error_handler(E_ERR,'task_sync', errstring, source)
endif

if (verbose) write(*,*) "PE", myrank, ": MPI Barrier released"

end subroutine task_sync


!-----------------------------------------------------------------------------

!> Send the srcarray to the destination id.
!> If time is specified, it is also sent in a separate communications call.  
!> This is a synchronous call; it will not return until the destination has 
!> called receive to accept the data.  If the send_to/receive_from calls are 
!> not paired correctly the code will hang.

subroutine send_to(dest_id, srcarray, time, label)
 integer, intent(in) :: dest_id
 real(r8), intent(in) :: srcarray(:)
 type(time_type), intent(in), optional :: time
 character(len=*), intent(in), optional :: label

integer :: tag, errcode
integer :: itime(2)
integer(i8) :: itemcount, offset, nextsize
real(r8), allocatable :: tmpdata(:)

if (verbose) write(*,*) "PE", myrank, ": start of send_to "

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'send_to', errstring, source)
endif

! simple idiotproofing
if ((dest_id < 0) .or. (dest_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "destination task id ", dest_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'send_to', errstring, source)
endif

itemcount = size(srcarray,KIND=i8)

if (present(label)) then
   write(*,*) trim(label)//" PE", myrank, ": send_to itemsize ", itemcount, " dest ", dest_id
else if (verbose) then
   write(*,*) "PE", myrank, ": send_to itemsize ", itemcount, " dest ", dest_id
endif

! use my task id as the tag; unused at this point.
tag = myrank

if (make_copy_before_sendrecv) allocate(tmpdata(min(itemcount, SNDRCV_MAXSIZE)))


if (itemcount <= SNDRCV_MAXSIZE) then

   if (verbose) write(*,*) "PE", myrank, ": send_to ", itemcount, " dest ", dest_id

   if (.not. make_copy_before_sendrecv) then
      call MPI_Ssend(srcarray, int(itemcount,i4), datasize, dest_id, tag, &
                    my_local_comm, errcode)
   else
      ! this copy should be unneeded, but on the intel fortran 9.0 compiler and mpich
      ! on one platform, calling this subroutine with an array section resulted in some
      ! apparently random memory corruption.  making a copy of the data into a local,
      ! contiguous buffer before send and receive fixed the corruption.  this shouldn't
      ! have been needed, and is a performance/memory sink.

      tmpdata = srcarray
      call MPI_Ssend(tmpdata, int(itemcount,i4), datasize, dest_id, tag, &
                    my_local_comm, errcode)
   endif
else
   ! there are a few places in the code where we send/receive a full state vector.
   ! as these get really large, they may exceed the limits of the MPI library.
   ! break large messages up into smaller chunks if needed.  
   offset = 1
   do while (offset <= itemcount)
      if (itemcount-offset >= SNDRCV_MAXSIZE) then
         nextsize = SNDRCV_MAXSIZE 
      else if (itemcount-offset == 0) then
         nextsize = 1
      else
         nextsize = itemcount - offset
      endif

      if (verbose) write(*,*) 'sending array items ', offset, ' thru ' , offset + nextsize - 1

      if (.not. make_copy_before_sendrecv) then
         call MPI_Ssend(srcarray(offset:offset+nextsize-1), int(nextsize,i4), datasize, dest_id, tag, &
                        my_local_comm, errcode)
      else
         tmpdata = srcarray(offset:offset+nextsize-1)
         call MPI_Ssend(tmpdata, int(nextsize,i4), datasize, dest_id, tag, &
                        my_local_comm, errcode)
      endif
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
         call error_handler(E_ERR,'send_to', errstring, source)
      endif
      offset = offset + nextsize
   enddo
endif

if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
   call error_handler(E_ERR,'send_to', errstring, source)
endif

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   if (verbose) write(*,*) "PE", myrank, ": time present"
   call get_time(time, itime(1), itime(2))
   call MPI_Ssend(itime, 2, MPI_INTEGER, dest_id, tag*2, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
      call error_handler(E_ERR,'send_to', errstring, source)
   endif
   if (verbose) write(*,*) "PE", myrank, ": sent time to ", dest_id
endif

if (make_copy_before_sendrecv) deallocate(tmpdata)


if (verbose) write(*,*) "PE", myrank, ": end of send_to "

end subroutine send_to


!-----------------------------------------------------------------------------

!> Receive data into the destination array from the src task.
!> If time is specified, it is received in a separate communications call.  
!> This is a synchronous call; it will not return until the source has 
!> sent the data.  If the send_to/receive_from calls are not paired correctly 
!> the code will hang.

subroutine receive_from(src_id, destarray, time, label)
 integer, intent(in) :: src_id
 real(r8), intent(inout) :: destarray(:)
 type(time_type), intent(out), optional :: time
 character(len=*), intent(in), optional :: label

integer :: tag, errcode
integer :: itime(2)
integer :: status(MPI_STATUS_SIZE)
integer(i8) :: itemcount, offset, nextsize
real(r8), allocatable :: tmpdata(:)

if (verbose) write(*,*) "PE", myrank, ": start of receive_from "

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'receive_from', errstring, source)
endif

! simple idiotproofing
if ((src_id < 0) .or. (src_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "source task id ", src_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'receive_from', errstring, source)
endif

itemcount = size(destarray,KIND=i8)

if (present(label)) then
   write(*,*) trim(label)//" PE", myrank, ": receive_from itemsize ", itemcount, " src ", src_id
else if (verbose) then
   write(*,*) "PE", myrank, ": receive_from itemsize ", itemcount, " src ", src_id
endif

! send_to uses its own id as the tag.
tag = src_id

if (verbose) write(*,*) "PE", myrank, ": receive ", itemcount, " src ", src_id

if (make_copy_before_sendrecv) allocate(tmpdata(min(itemcount,SNDRCV_MAXSIZE)))


if (itemcount <= SNDRCV_MAXSIZE) then
   if (.not. make_copy_before_sendrecv) then
      call MPI_Recv(destarray, int(itemcount,i4), datasize, src_id, MPI_ANY_TAG, &
                 my_local_comm, status, errcode)
   else
      ! this copy should be unneeded, but on the intel fortran 9.0 compiler and mpich
      ! on one platform, calling this subroutine with an array section resulted in some
      ! apparently random memory corruption.  making a copy of the data into a local,
      ! contiguous buffer before send and receive fixed the corruption.  this shouldn't
      ! have been needed, and is a performance/memory sink.

      call MPI_Recv(tmpdata, int(itemcount,i4), datasize, src_id, MPI_ANY_TAG, &
                    my_local_comm, status, errcode)
      destarray = tmpdata
   endif
else
   ! there are a few places in the code where we send/receive a full state vector.
   ! as these get really large, they may exceed the limits of the MPI library.
   ! break large messages up into smaller chunks if needed.  
   offset = 1
   do while (offset <= itemcount)
      if (itemcount-offset >= SNDRCV_MAXSIZE) then
         nextsize = SNDRCV_MAXSIZE 
      else if (itemcount-offset == 0) then
         nextsize = 1
      else
         nextsize = itemcount - offset
      endif

      if (verbose) write(*,*) 'recving array items ', offset, ' thru ' , offset + nextsize - 1

      if (.not. make_copy_before_sendrecv) then
         call MPI_Recv(destarray(offset:offset+nextsize-1), int(nextsize,i4), datasize, src_id, MPI_ANY_TAG, &
                       my_local_comm, status, errcode)
      else
         call MPI_Recv(tmpdata, int(nextsize,i4), datasize, src_id, MPI_ANY_TAG, &
                       my_local_comm, status, errcode)
         destarray(offset:offset+nextsize-1) = tmpdata(1:nextsize)
      endif
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
         call error_handler(E_ERR,'receive_from', errstring, source)
      endif
      offset = offset + nextsize
   end do
endif

if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
   call error_handler(E_ERR,'receive_from', errstring, source)
endif

if (verbose) write(*,*) "PE", myrank, ": received from ", src_id

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   if (verbose) write(*,*) "PE", myrank, ": time present"
   call MPI_Recv(itime, 2, MPI_INTEGER, src_id, tag*2, &
                 my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'receive_from', errstring, source)
   endif
   if (itime(2) < 0) then
      write(errstring, '(a,i8)') 'seconds in settime were < 0; using 0'
      call error_handler(E_MSG,'receive_from', errstring, source)
      time = set_time(itime(1), 0)
   else
      time = set_time(itime(1), itime(2))
   endif
   if (verbose) write(*,*) "PE", myrank, ": received time from ", src_id
endif

if (make_copy_before_sendrecv) deallocate(tmpdata)

if (verbose) write(*,*) "PE", myrank, ": end of receive_from "


end subroutine receive_from


!-----------------------------------------------------------------------------
! NOTE: so far we only seem to be sending real data, not integer, and 1D arrays.
!       this could be overloaded to send 2D arrays, and ints if needed.

!> The data array values on the root task will be broadcast to every other
!> task.  When this routine returns, all tasks will have the contents of the
!> root array in their own arrays.  Thus 'array' is intent(in) on root, and
!> intent(out) on all other tasks.

subroutine array_broadcast(array, root, icount)
 real(r8), intent(inout) :: array(:)
 integer, intent(in) :: root
 integer, intent(in), optional :: icount

integer :: itemcount, errcode, offset, nextsize
real(r8), allocatable :: tmpdata(:)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'array_broadcast', errstring, source)
endif

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source)
endif

! if the actual data to be sent is shorter than the size of 'array',
! there are performance advantages to sending only the actual data
! and not the full array size.  (performance tested on an ibm machine.)
! calling code must determine if this is the case and pass in a length
! shorter than the array size.
if (present(icount)) then
   if (icount > size(array)) then
      write(errstring,  '(a,i12)') "number of items to broadcast: ", icount
      write(errstring1, '(a,i12)') "cannot be larger than the array size: ", size(array)
      call error_handler(E_ERR,'array_broadcast', errstring, source, &
                         text2=errstring1)
   endif
   itemcount = icount
else
   itemcount = size(array)
endif

if (verbose .and. myrank == root) write(*,*) "PE", myrank, ": bcast itemsize from here ", itemcount

! comment this in only if you really have to.  it will flood the output with
! messages from every task except one.
!if (verbose .and. myrank /= root) write(*,*) "PE", myrank, ": bcast itemsize ", itemcount, " root ", root

if (make_copy_before_broadcast) allocate(tmpdata(min(itemcount,BCAST_MAXSIZE)))

! at least one user has run into a limit in the MPI implementation where you
! cannot broadcast too large an array.  if the size of this array is too large,
! broadcast it in chunks until all the data has been processed
if (itemcount <= BCAST_MAXSIZE) then
   if (.not. make_copy_before_broadcast) then
      call MPI_Bcast(array, itemcount, datasize, root, my_local_comm, errcode)
   else
      if (my_task_id() == root) tmpdata = array
      call MPI_Bcast(tmpdata, itemcount, datasize, root, my_local_comm, errcode)
      if (my_task_id() /= root) array = tmpdata
   endif
else
   offset = 1
   do while (offset <= itemcount)
      if (itemcount-offset >= BCAST_MAXSIZE) then
         nextsize = BCAST_MAXSIZE 
      else if (itemcount-offset == 0) then
         nextsize = 1
      else
         nextsize = itemcount - offset 
      endif

      if (verbose) write(*,*) 'bcasting array items ', offset, ' thru ' , offset + nextsize - 1

      ! test this - are array sections going to cause array temps to be created?
      if (.not. make_copy_before_broadcast) then
         call MPI_Bcast(array(offset:offset+nextsize-1), nextsize, datasize, root, my_local_comm, errcode)
      else
         if (my_task_id() == root) tmpdata = array(offset:offset+nextsize-1)
         call MPI_Bcast(tmpdata, nextsize, datasize, root, my_local_comm, errcode)
         if (my_task_id() /= root) array(offset:offset+nextsize-1) = tmpdata
      endif
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Bcast returned error code ', errcode
         call error_handler(E_ERR,'array_broadcast', errstring, source)
      endif
      offset = offset + nextsize
   end do
endif

if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Bcast returned error code ', errcode
   call error_handler(E_ERR,'array_broadcast', errstring, source)
endif

if (make_copy_before_broadcast) deallocate(tmpdata)

end subroutine array_broadcast


!-----------------------------------------------------------------------------
! DART-specific cover utilities
!-----------------------------------------------------------------------------

!> Return .TRUE. if my local task id is 0, .FALSE. otherwise.
!> (Task numbers in MPI start at 0, contrary to the rules of polite fortran.)

function iam_task0()

logical :: iam_task0

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'iam_task0', errstring, source)
endif

iam_task0 = (myrank == 0)

end function iam_task0

!-----------------------------------------------------------------------------

!> this must be paired with the same number of broadcast_recv()s on all 
!> other tasks.  it will not return until all tasks in the communications 
!> group have made the call.
!>
!> cover routine for array broadcast.  one additional sanity check -- make 
!> sure the 'from' matches my local task id.  also, these arrays are
!> intent(in) here, but they call a routine which is intent(inout) so they
!> must be the same here.

subroutine broadcast_send(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)
 integer, intent(in) :: from
! arrays are really only intent(in) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:)
 real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

real(r8) :: packbuf(PACKLIMIT)
real(r8) :: local(5)
logical  :: doscalar, morethanone
integer  :: itemcount

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_send', errstring, source)
endif

! simple idiotproofing
if (from /= myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "must be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_send', errstring, source)
endif

! for relatively small array sizes, pack them into a single send/recv pair.
call countup(array1, array2, array3, array4, array5, &
             scalar1, scalar2, scalar3, scalar4, scalar5, &
             itemcount, morethanone, doscalar)

if (itemcount <= PACKLIMIT .and. morethanone) then

   call packit(packbuf, array1, array2, array3, array4, array5, doscalar, &
                         scalar1, scalar2, scalar3, scalar4, scalar5)

   call array_broadcast(packbuf, from, itemcount)

else

   call array_broadcast(array1, from)

   if (morethanone) then
      if (present(array2)) call array_broadcast(array2, from)
      if (present(array3)) call array_broadcast(array3, from)
      if (present(array4)) call array_broadcast(array4, from)
      if (present(array5)) call array_broadcast(array5, from)
      
      if (doscalar) then
         call packscalar(local, scalar1, scalar2, scalar3, scalar4, scalar5)
         call array_broadcast(local, from)
      endif

   endif
endif


end subroutine broadcast_send

!-----------------------------------------------------------------------------

!> this must be paired with a single broadcast_send() on one other task, and
!> broadcast_recv() on all other tasks, and it must match exactly the number 
!> of args in the sending call.
!> it will not return until all tasks in the communications group have
!> made the call.
!>
!> cover routine for array broadcast.  one additional sanity check -- make 
!> sure the 'from' is not the same as my local task id.  these arrays are
!> intent(out) here, but they call a routine which is intent(inout) so they
!> must be the same here.

subroutine broadcast_recv(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)
 integer, intent(in) :: from
! arrays are really only intent(out) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:)
 real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

real(r8) :: packbuf(PACKLIMIT)
real(r8) :: local(5)
logical :: doscalar, morethanone
integer :: itemcount

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_recv', errstring, source)
endif

! simple idiotproofing
if (from == myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "cannot be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_recv', errstring, source)
endif

! for relatively small array sizes, pack them into a single send/recv pair.
call countup(array1, array2, array3, array4, array5, &
             scalar1, scalar2, scalar3, scalar4, scalar5, &
             itemcount, morethanone, doscalar)

if (itemcount <= PACKLIMIT .and. morethanone) then

   call array_broadcast(packbuf, from, itemcount)

   call unpackit(packbuf, array1, array2, array3, array4, array5, doscalar, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)

else

   call array_broadcast(array1, from)

   if (morethanone) then
      if (present(array2)) call array_broadcast(array2, from)
      if (present(array3)) call array_broadcast(array3, from)
      if (present(array4)) call array_broadcast(array4, from)
      if (present(array5)) call array_broadcast(array5, from)
   
      if (doscalar) then
         call array_broadcast(local, from)
         call unpackscalar(local, scalar1, scalar2, scalar3, &
                           scalar4, scalar5)
      endif

   endif

endif

end subroutine broadcast_recv

!-----------------------------------------------------------------------------

!> figure out how many items are in the specified arrays, total.
!> also note if there's more than a single array (array1) to send,
!> and if there are any scalars specified.

subroutine countup(array1, array2, array3, array4, array5, &
                   scalar1, scalar2, scalar3, scalar4, scalar5, &
                   numitems, morethanone, doscalar)
 real(r8), intent(in)           :: array1(:)
 real(r8), intent(in), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(in), optional :: scalar1, scalar2, scalar3, scalar4, scalar5
 integer,  intent(out)          :: numitems
 logical,  intent(out)          :: morethanone, doscalar

morethanone = .false.
numitems = size(array1)

if (present(array2)) then
   numitems = numitems + size(array2)
   morethanone = .true.
endif
if (present(array3)) then
   numitems = numitems + size(array3)
   morethanone = .true.
endif
if (present(array4)) then
   numitems = numitems + size(array4)
   morethanone = .true.
endif
if (present(array5)) then
   numitems = numitems + size(array5)
   morethanone = .true.
endif
if (present(scalar1)) then 
   numitems = numitems + 1
   morethanone = .true.
   doscalar = .true.
endif
if (present(scalar2)) then
   numitems = numitems + 1
   morethanone = .true.
   doscalar = .true.
endif
if (present(scalar3)) then
   numitems = numitems + 1
   morethanone = .true.
   doscalar = .true.
endif
if (present(scalar4)) then
   numitems = numitems + 1
   morethanone = .true.
   doscalar = .true.
endif
if (present(scalar5)) then
   numitems = numitems + 1
   morethanone = .true.
   doscalar = .true.
endif

end subroutine countup

!-----------------------------------------------------------------------------

!> pack multiple small arrays into a single buffer before sending.

subroutine packit(buf, array1, array2, array3, array4, array5, doscalar, &
                       scalar1, scalar2, scalar3, scalar4, scalar5)
 real(r8), intent(out)          :: buf(:)
 real(r8), intent(in)           :: array1(:)
 real(r8), intent(in), optional :: array2(:), array3(:), array4(:), array5(:)
 logical,  intent(in)           :: doscalar
 real(r8), intent(in), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

integer :: sindex, eindex

sindex = 1
eindex = sindex + size(array1) - 1
buf(sindex:eindex) = array1(:)
sindex = eindex+1

if (present(array2)) then
   eindex = sindex + size(array2) - 1
   buf(sindex:eindex) = array2(:)
   sindex = eindex+1
endif

if (present(array3)) then
   eindex = sindex + size(array3) - 1
   buf(sindex:eindex) = array3(:)
   sindex = eindex+1
endif

if (present(array4)) then
   eindex = sindex + size(array4) - 1
   buf(sindex:eindex) = array4(:)
   sindex = eindex+1
endif

if (present(array5)) then
   eindex = sindex + size(array5) - 1
   buf(sindex:eindex) = array5(:)
   sindex = eindex+1
endif

if (doscalar) then
   if (present(scalar1)) then
      buf(sindex) = scalar1
      sindex = sindex+1
   endif

   if (present(scalar2)) then
      buf(sindex) = scalar2
      sindex = sindex+1
   endif

   if (present(scalar3)) then
      buf(sindex) = scalar3
      sindex = sindex+1
   endif

   if (present(scalar4)) then
      buf(sindex) = scalar4
      sindex = sindex+1
   endif

   if (present(scalar5)) then
      buf(sindex) = scalar5
      sindex = sindex+1
   endif
endif

end subroutine packit

!-----------------------------------------------------------------------------

!> unpack multiple small arrays from a single buffer after receiving.

subroutine unpackit(buf, array1, array2, array3, array4, array5, doscalar, &
                         scalar1, scalar2, scalar3, scalar4, scalar5)
 real(r8), intent(in)            :: buf(:)
 real(r8), intent(out)           :: array1(:)
 real(r8), intent(out), optional :: array2(:), array3(:), array4(:), array5(:)
 logical,  intent(in)            :: doscalar
 real(r8), intent(out), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

integer :: sindex, eindex

sindex = 1
eindex = sindex + size(array1) - 1
array1(:) = buf(sindex:eindex)
sindex = eindex+1

if (present(array2)) then
   eindex = sindex + size(array2) - 1
   array2(:) = buf(sindex:eindex)
   sindex = eindex+1
endif

if (present(array3)) then
   eindex = sindex + size(array3) - 1
   array3(:) = buf(sindex:eindex) 
   sindex = eindex+1
endif

if (present(array4)) then
   eindex = sindex + size(array4) - 1
   array4(:) = buf(sindex:eindex)
   sindex = eindex+1
endif

if (present(array5)) then
   eindex = sindex + size(array5) - 1
   array5(:) = buf(sindex:eindex)
   sindex = eindex+1
endif

if (doscalar) then
   if (present(scalar1)) then
      scalar1 = buf(sindex)
      sindex = sindex+1
   endif
   
   if (present(scalar2)) then
      scalar2 = buf(sindex)
      sindex = sindex+1
   endif
   
   if (present(scalar3)) then
      scalar3 = buf(sindex)
      sindex = sindex+1
   endif
   
   if (present(scalar4)) then
      scalar4 = buf(sindex)
      sindex = sindex+1
   endif
   
   if (present(scalar5)) then
      scalar5 = buf(sindex)
      sindex = sindex+1
   endif
endif

end subroutine unpackit

!-----------------------------------------------------------------------------

!> for any values specified, pack into a single array

subroutine packscalar(local, scalar1, scalar2, scalar3, scalar4, scalar5)
 real(r8), intent(out)          :: local(5) 
 real(r8), intent(in), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

local = 0.0_r8
      
if (present(scalar1)) local(1) = scalar1
if (present(scalar2)) local(2) = scalar2
if (present(scalar3)) local(3) = scalar3
if (present(scalar4)) local(4) = scalar4
if (present(scalar5)) local(5) = scalar5

end subroutine packscalar
   
!-----------------------------------------------------------------------------

!> for any values specified, unpack from a single array

subroutine unpackscalar(local, scalar1, scalar2, scalar3, scalar4, scalar5)
 real(r8), intent(in)            :: local(5) 
 real(r8), intent(out), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

if (present(scalar1)) scalar1 = local(1)
if (present(scalar2)) scalar2 = local(2)
if (present(scalar3)) scalar3 = local(3)
if (present(scalar4)) scalar4 = local(4)
if (present(scalar5)) scalar5 = local(5)

end subroutine unpackscalar
   
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! overloaded global reduce routines

! The external32 representations of the datatypes returned by 
! MPI_TYPE_CREATE_F90_REAL/COMPLEX/INTEGER are given by the following rules.
! For MPI_TYPE_CREATE_F90_REAL:
! 
!    if      (p > 33) or (r > 4931) then  external32 representation 
!                                         is undefined   
!    else if (p > 15) or (r >  307) then  external32_size = 16 
!    else if (p >  6) or (r >   37) then  external32_size =  8 
!    else                                 external32_size =  4 
! 
! For MPI_TYPE_CREATE_F90_COMPLEX: twice the size as for MPI_TYPE_CREATE_F90_REAL.
! For MPI_TYPE_CREATE_F90_INTEGER:
! 
!    if      (r > 38) then  external32 representation is undefined 
!    else if (r > 18) then  external32_size =  16  
!    else if (r >  9) then  external32_size =  8  
!    else if (r >  4) then  external32_size =  4 
!    else if (r >  2) then  external32_size =  2  
!    else                   external32_size =  1  
!
!
!-----------------------------------------------------------------------------

!> take values from each task, add them, and return
!> the sum to all tasks.  integer version

subroutine sum_across_tasks_int4(addend, sum)
 integer, intent(in) :: addend
 integer, intent(out) :: sum

 integer :: errcode
 integer :: localaddend(1), localsum(1)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'sum_across_tasks', errstring, source)
endif

localaddend(1) = addend

if (verbose) write(*,*) "PE", myrank, ": Allreduce called"

call MPI_Allreduce(localaddend, localsum, 1, MPI_INTEGER, MPI_SUM, &
                   my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Allreduce returned error code ', errcode
   call error_handler(E_ERR,'sum_across_tasks', errstring, source)
endif

sum = localsum(1)

end subroutine sum_across_tasks_int4

!-----------------------------------------------------------------------------

!> take values from each task, add them, and return
!> the sum to all tasks. long integer version.

subroutine sum_across_tasks_int8(addend, sum)
 integer(i8), intent(in)  :: addend
 integer(i8), intent(out) :: sum

 integer :: errcode
 integer(i8) :: localaddend(1), localsum(1)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'sum_across_tasks', errstring, source)
endif

localaddend(1) = addend

if (verbose) write(*,*) "PE", myrank, ": Allreduce called"

call MPI_Allreduce(localaddend, localsum, 1, longinttype, MPI_SUM, &
                   my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Allreduce returned error code ', errcode
   call error_handler(E_ERR,'sum_across_tasks', errstring, source)
endif

sum = localsum(1)

end subroutine sum_across_tasks_int8
!-----------------------------------------------------------------------------

!> take values from each task, add them, and return
!> the sum to all tasks. real version.

subroutine sum_across_tasks_real(addend, sum)
 real(r8), intent(in) :: addend
 real(r8), intent(out) :: sum

 integer :: errcode
 real(r8) :: localaddend(1), localsum(1)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'sum_across_tasks', errstring, source)
endif

localaddend(1) = addend

if (verbose) write(*,*) "PE", myrank, ": Allreduce called"

call MPI_Allreduce(localaddend, localsum, 1, datasize, MPI_SUM, &
                   my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Allreduce returned error code ', errcode
   call error_handler(E_ERR,'sum_across_tasks', errstring, source)
endif

sum = localsum(1)

end subroutine sum_across_tasks_real

!-----------------------------------------------------------------------------

!> Sum array items across all tasks and send
!> results in an array of same size to one task.

subroutine send_sum_to(local_val, task, global_val)

real(r8), intent(in)  :: local_val(:)  !! addend vals on each task
integer,  intent(in)  :: task          !! task to collect on
real(r8), intent(out) :: global_val(:) !! results returned only on given task

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'send_sum_to', errstring, source)
endif

! collect values on a single given task 
call mpi_reduce(local_val(:), global_val(:), size(global_val), datasize, MPI_SUM, &
                task, get_dart_mpi_comm(), errcode)

end subroutine send_sum_to

!-----------------------------------------------------------------------------

!> Collect global min and max values on one task.

subroutine send_minmax_to(minmax, task, global_val)

real(r8), intent(in)  :: minmax(2)     !! min max on each task
integer,  intent(in)  :: task          !! task to collect on
real(r8), intent(out) :: global_val(2) !! results returned only on given task

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'send_minmax_to', errstring, source)
endif

! collect values on a single given task 
call mpi_reduce(minmax(1:1), global_val(1:1), 1, datasize, MPI_MIN, task, get_dart_mpi_comm(), errcode)
call mpi_reduce(minmax(2:2), global_val(2:2), 1, datasize, MPI_MAX, task, get_dart_mpi_comm(), errcode)

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

!> Find min and max of each element of an array, put the result on every task.
!> Overwrites arrays min_var, max_var with the minimum and maximum for each 
!> element across all tasks.

subroutine broadcast_minmax(min_var, max_var, num_elements)

integer,  intent(in)    :: num_elements
real(r8), intent(inout) :: min_var(num_elements)
real(r8), intent(inout) :: max_var(num_elements)

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_minmax', errstring, source)
endif

call mpi_allreduce(MPI_IN_PLACE, min_var, num_elements, datasize, MPI_MIN, get_dart_mpi_comm(), errcode)
call mpi_allreduce(MPI_IN_PLACE, max_var, num_elements, datasize, MPI_MAX, get_dart_mpi_comm(), errcode)

end subroutine broadcast_minmax

!-----------------------------------------------------------------------------
!> Broadcast logical

subroutine broadcast_flag(flag, root)

logical, intent(inout) :: flag
integer, intent(in)    :: root !! relative to get_dart_mpi_comm()

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_flag', errstring, source)
endif

call MPI_Bcast(flag, 1, MPI_LOGICAL, root, my_local_comm, errcode)

end subroutine broadcast_flag

!-----------------------------------------------------------------------------
! pipe-related utilities - used in 'async4' handshakes between mpi jobs
! and scripting to allow filter and an mpi model to alternate execution.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

!> block by reading a named pipe file until some other task
!> writes a string into it.  this ensures the task is not
!> spinning and using CPU cycles, but is asleep waiting in
!> the kernel.   one subtlety with this approach is that even
!> though named pipes are created in the filesystem, they are
!> implemented in the kernel, so on a multiprocessor machine
!> the write into the pipe file must occur on the same PE as
!> the reader is waiting.  see the 'wakeup_filter' program for
!> the MPI job which spreads out on all the PEs for this job
!> and writes into the file from the correct PE.

subroutine block_task()

character(len = 32) :: fifo_name, filter_to_model, model_to_filter, non_pipe
integer :: rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'block_task', errstring, source)
endif

! FIXME: this should be mpi or a string other than filter (this is generic 
! mpi wrapper code, callable from programs other than filter.)
filter_to_model = 'filter_to_model.lock'
model_to_filter = 'model_to_filter.lock'
non_pipe = 'filter_to_model.file'

! the i5.5 format below will not handle task counts larger than this.
if (total_tasks > 99999) then
   write(errstring, *) 'cannot handle task counts > 99999'
   call error_handler(E_ERR,'block_task', errstring, source)
endif

! make it so we only have to test 1 or 2 things here instead of 3
! when deciding whether to print status messages.
if (verbose) async4_verbose = .TRUE.

if ((myrank == head_task) .and. separate_node_sync) then

   if (async4_verbose) then
      write(*,*)  'checking master task host'
      rc = shell_execute('echo master task running on host `hostname`')
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose .or. print4status) write(*,*) 'MPI job telling script to advance model'
   rc = shell_execute('echo advance > '//trim(non_pipe))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

endif

if ((myrank == head_task) .and. .not. separate_node_sync) then

   if (async4_verbose) then
      write(*,*)  'checking master task host'
      rc = shell_execute('echo master task running on host `hostname`')
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose .or. print4status) write(*,*) 'MPI job telling script to advance model'
   rc = shell_execute('echo advance > '//trim(filter_to_model))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'MPI job now waiting to read from lock file'
   rc = shell_execute('cat < '//trim(model_to_filter)//'> /dev/null')
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

else

! if you change this in any way, change the corresponding string in 
! restart_task() below.
   ! FIXME: this should be 'task_lock', since it's generic code beyond filter.
   write(fifo_name, '(a, i5.5)') "filter_lock", myrank
   
   if (async4_verbose) then
      write(*,*)  'checking slave task host'
      rc = shell_execute('echo '//trim(fifo_name)//' accessed from host `hostname`')
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose) write(*,*) 'removing any previous lock file: '//trim(fifo_name)
   rc = shell_execute('rm -f '//trim(fifo_name))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'made fifo, named: '//trim(fifo_name)
   rc = shell_execute('mkfifo '//trim(fifo_name))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'ready to read from lock file: '//trim(fifo_name)
   rc = shell_execute('cat < '//trim(fifo_name)//'> /dev/null')
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'got response, removing lock file: '//trim(fifo_name)
   rc = shell_execute('rm -f '//trim(fifo_name))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

endif

! make sure all tasks get here before any get to proceed further.
! this hides a multitude of sins.  it could also cause
! tasks to hang forever, but right now it also makes some
! systems able to run async 4.  maybe this should be namelist
! selectable.  something named 'horrible_mpi_hack = .true.'
! if tasks are hanging out here instead of reading the lock files,
! they are burning cpu cycles and competing with the model advances.
! but, it works, as opposed to not working.
call task_sync()

end subroutine block_task

!-----------------------------------------------------------------------------

!> companion to block_task.  must be called by a different executable
!> and it writes into the named pipes to restart the waiting task.

subroutine restart_task()

character(len = 32) :: fifo_name, model_to_filter
integer :: rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'restart_task', errstring, source)
endif

! FIXME: ditto previous comment about using the string 'filter' here.
model_to_filter = 'model_to_filter.lock'

! the i5.5 format below will not handle task counts larger than this.
if (total_tasks > 99999) then
   write(errstring, *) 'cannot handle task counts > 99999'
   call error_handler(E_ERR,'block_task', errstring, source)
endif

! make it so we only have to test 1 or 2 things here instead of 3
! when deciding whether to print status messages.
if (verbose) async4_verbose = .TRUE.

! process 0 (or N-1) is handled differently in the code.
if ((myrank == head_task) .and. .not. separate_node_sync) then

   if (async4_verbose) then
      rc = shell_execute('echo master task running on host `hostname`')
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose .or. print4status) write(*,*) 'script telling MPI job ok to restart'
   rc = shell_execute('echo restart > '//trim(model_to_filter))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

else

   if (async4_verbose) write(*,*) 'waking up task id ', myrank

   ! FIXME: this should be 'task_lock', since it's generic code beyond filter.
   write(fifo_name,"(a,i5.5)") "filter_lock", myrank

   if (async4_verbose) then
      rc = shell_execute('echo '//trim(fifo_name)//' accessed from host `hostname`')
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose) write(*,*) 'ready to write to lock file: '//trim(fifo_name)
   rc = shell_execute('echo restart > '//trim(fifo_name))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

endif

end subroutine restart_task

!-----------------------------------------------------------------------------

!> must be called when filter is exiting so calling script knows the job is over.

subroutine finished_task(async)
 integer, intent(in) :: async

character(len = 32) :: filter_to_model, non_pipe
integer :: rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'restart_task', errstring, source)
endif

! only in the async=4 case does this matter.
if (async /= 4) return

! FIXME: ditto previous comment about using the string 'filter' here.
filter_to_model = 'filter_to_model.lock'
non_pipe = 'filter_to_model.file'


! only process 0 (or N-1) needs to do anything.
if (myrank == head_task) then

   if (print4status .or. verbose) write(*,*) 'MPI task telling script we are done'
   if (separate_node_sync) then
      rc = shell_execute('echo finished > '//trim(non_pipe))
   else
      rc = shell_execute('echo finished > '//trim(filter_to_model))
   endif

   
endif

end subroutine finished_task

!-----------------------------------------------------------------------------
! general system util wrappers.
!-----------------------------------------------------------------------------

!> Use the system() command to execute a command string.
!> Will wait for the command to complete and returns an
!> error code unless you end the command with & to put
!> it into background.   Function which returns the rc
!> of the command, 0 being all is ok.
!>
!> allow code to test the theory that maybe the system call is
!> not reentrant on some platforms.  if serialize is set and
!> is true, do each call serially.

function shell_execute(execute_string, serialize)
 character(len=*), intent(in) :: execute_string
 logical, intent(in), optional :: serialize
 integer :: shell_execute

logical :: all_at_once
integer :: errcode, dummy(1)
integer :: status(MPI_STATUS_SIZE)

if (verbose) async2_verbose = .true.

! default to everyone running concurrently, but if set and not true,
! serialize the calls to system() so they do not step on each other.
if (present(serialize)) then
   all_at_once = .not. serialize
else
   all_at_once = .TRUE.
endif

if (async2_verbose) write(*,*) "PE", myrank, ": system string is: ", trim(execute_string)
shell_execute = -1

! this is the normal (default) case
if (all_at_once) then

   ! all tasks call system at the same time
   call do_system(execute_string, shell_execute)
   if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

   return
endif

! only one task at a time calls system, and all wait their turn by
! making each task wait for a message from the (N-1)th task.

! this is used only to signal; the value it contains is unused.
dummy = 0

if (myrank == 0) then

   ! my turn to execute
   call do_system(execute_string, shell_execute)
   if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

   if (total_tasks > 1) then
      ! tell next task it can continue
      call MPI_Send(dummy, 1, MPI_INTEGER, 1, 1, my_local_comm, errcode)
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Send returned error code ', &
                                     errcode
         call error_handler(E_ERR,'shell_execute', errstring, source)
      endif
   endif

else if (myrank /= (total_tasks-1)) then
   ! wait for (me-1) to tell me it is my turn
   call MPI_Recv(dummy, 1, MPI_INTEGER, myrank-1, myrank, &
                 my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'shell_execute', errstring, source)
   endif

   ! my turn to execute
   call do_system(execute_string, shell_execute)
   if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

   ! and now tell (me+1) to go
   call MPI_Send(dummy, 1, MPI_INTEGER, myrank+1, myrank+1, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Send returned error code ', &
                                  errcode
      call error_handler(E_ERR,'shell_execute', errstring, source)
   endif
else
   ! last task, no one else to send to.
   call MPI_Recv(dummy, 1, MPI_INTEGER, myrank-1, myrank, &
                 my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'shell_execute', errstring, source)
   endif

   ! my turn to execute
   call do_system(execute_string, shell_execute)
   if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

endif

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

base = mpi_wtime()

end subroutine start_mpi_timer

!-----------------------------------------------------------------------------

!> return the time since the last call to start_timer().
!> can call multiple times to get running times.
!> call with a different base for nested timers.

function read_mpi_timer(base)

real(digits12), intent(in) :: base
real(digits12) :: read_mpi_timer

real(digits12) :: now

now = mpi_wtime()

read_mpi_timer = now - base

end function read_mpi_timer


!-----------------------------------------------------------------------------
!> return our communicator

function get_dart_mpi_comm()
 integer :: get_dart_mpi_comm

 get_dart_mpi_comm = my_local_comm

end function get_dart_mpi_comm

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! One sided communication

subroutine get_from_mean(owner, window, mindex, x)

integer,  intent(in)  :: owner  ! task in the window that owns the memory
integer,  intent(in)  :: window ! window object
integer,  intent(in)  :: mindex ! index in the tasks memory
real(r8), intent(out) :: x      ! result

integer(KIND=MPI_ADDRESS_KIND) :: target_disp
integer :: errcode

! Note to programmer: The data transfer is not guaranteed
! to have occured until the call to mpi_win_unlock. 
! => Don't do anything with x in between mpi_get and mpi_win_lock

! Note to programmer: openmpi 1.10.0 does not
! allow scalars in mpi calls. openmpi 1.10.1 fixes this.

target_disp = (mindex - 1)
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, window, errcode)
call mpi_get(x, 1, datasize, owner, target_disp, 1, datasize, window, errcode)
call mpi_win_unlock(owner, window, errcode)

end subroutine get_from_mean

!-----------------------------------------------------------------------------

subroutine get_from_fwd(owner, window, mindex, num_rows, x)

integer,  intent(in)  :: owner    ! task in the window that owns the memory
integer,  intent(in)  :: window   ! window object
integer,  intent(in)  :: mindex   ! index in the tasks memory
integer,  intent(in)  :: num_rows ! number of rows in the window
real(r8), intent(out) :: x(:)     ! result

integer(KIND=MPI_ADDRESS_KIND) :: target_disp
integer :: errcode

! Note to programmer: The data transfer is not guaranteed
! to have occured until the call to mpi_win_unlock. 
! => Don't do anything with x in between mpi_get and mpi_win_lock

target_disp = (mindex - 1)*num_rows
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, window, errcode)
call mpi_get(x, num_rows, datasize, owner, target_disp, num_rows, datasize, window, errcode)
call mpi_win_unlock(owner, window, errcode)

end subroutine get_from_fwd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

!> Collect global max values on each task.

subroutine get_global_max(max)

real(r8), intent(inout)  :: max       !> global max over tasks
integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'get_global_max', errstring, source)
endif

! collect max values over al tasks
call mpi_allreduce(MPI_IN_PLACE, max, 1, datasize, MPI_MAX, my_local_comm, errcode)

end subroutine get_global_max

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


end module mpi_utilities_mod

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!> NOTE: non-module code, so this subroutine can be called from the
!>  utilities module, which this module uses (and cannot have circular refs)
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!> In case of error, call this instead of the fortran intrinsic exit().
!> It will signal the other MPI tasks that something bad happened and they
!> should also exit.

subroutine exit_all(exit_code)
 use mpi_utilities_mod, only : get_dart_mpi_comm

 integer, intent(in) :: exit_code

integer :: ierror

! call abort on our communicator

!print *, 'calling abort on comm ', get_dart_mpi_comm()
call MPI_Abort(get_dart_mpi_comm(),  exit_code, ierror)

! execution should never get here

end subroutine exit_all

!-----------------------------------------------------------------------------

