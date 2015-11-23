! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module mpi_utilities_mod

!-----------------------------------------------------------------------------
!
!   A collection of interfaces to the MPI (Message Passing Interface)
!   multi-processor communication library routines.
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
!    # exit_all()         Subroutine that substitutes for the intrinsic exit.
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
!    # broadcast_send()   Subroutine which takes up to 5 real arrays and up to 
!                         5 scalar reals and broadcasts them to all other tasks.
!                         One array is required; the rest are optional.
!                         If sender ID is not the same as the local task ID
!                         it is an error.  This call does not return until
!                         all other tasks have called broadcast_recv() to pick
!                         up the data.  The arg lists for the data arrays and
!                         scalar values must match exactly between send/recv.
!
!    # broadcast_recv()   Subroutine which takes up to 5 real arrays and up to 
!                         5 scalar reals and receives them from a sending task.
!                         One array is required; the rest are optional.
!                         If sender ID is the same as the local task ID
!                         it is an error.  This call does not return until
!                         all other tasks have called broadcast_recv() to pick
!                         up the data.  The arg lists for the data arrays and
!                         scalar values must match exactly between send/recv.
!
!    # sum_across_tasks() Subroutine which takes a single integer argument
!                         from each task, and returns the sum of all integers
!                         across all tasks back to all tasks.  All tasks must
!                         call this routine before it can compute and return
!                         the value.
!
!   # get_dart_mpi_comm() This code creates a private communicator for use by
!                         the routines in this file.  This function returns it.
!
!   Lower level utility routines which interact with the utilities_mod.f90 
!   code to open a named pipe per MPI task, read and write from them, and
!   close and/or remove them.
!
!    # block_task()       Create a named pipe (fifo) and read from it to
!                         block the process in such a way that it consumes
!                         no CPU time.  NOTE that once you put yourself to
!                         sleep, you cannot wake yourself up.  Some other MPI
!                         task must call restart_task(), on the same set of
!                         processors the original program was distributed over.
!
!    # restart_task()     Write into the pipe to restart the reading task.
!                         Note that this must be an entirely separate executable
!                         from the one which called block_task(), because it is
!                         asleep like Sleeping Beauty and cannot wake itself.
!                         See filter and wakeup_filter for examples
!                         of a program pair which uses these calls.
! 
!    # finished_task()    Called from finalize_mpi_utilities, if async=4 
!                         writes out a string to the model pipe to tell it 
!                         the main executable is exiting.
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
!  *** read_pipe()        The character string is read from the pipe.
!                         (Can be overloaded to read ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine blocks until data is available.
!
!  *** write_pipe()       The character string is written to the pipe.
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

use types_mod, only : r8, digits12
use utilities_mod, only : register_module, error_handler, & 
                          E_ERR, E_WARN, E_MSG, E_DBG, get_unit, close_file, &
                          set_output, set_tasknum, initialize_utilities,     &
                          finalize_utilities,                                &
                          nmlfileunit, do_output, do_nml_file, do_nml_term,  &
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
! $DART/doc/mpi directory. 

use mpi

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
 ! interface block for getting return code back from system() routine
 interface
  function system(string)
   character(len=*) :: string
   integer :: system
  end function system
 end interface
 ! end block
 !!SYSTEM_BLOCK_EDIT END COMMENTED_IN


!   ---- private data for mpi_utilities ----

integer :: myrank        = -1  ! my mpi number
integer :: total_tasks   = -1  ! total mpi tasks/procs
integer :: my_local_comm =  0  ! duplicate communicator private to this file
integer :: datasize      =  0  ! which MPI type corresponds to our r8 definition


public :: initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, get_dart_mpi_comm, datasize, reduce_min_max,     &
          get_from_fwd, get_from_mean

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical :: module_initialized   = .false.

character(len = 129) :: saved_progname = ''
character(len = 129) :: shell_name = ''   ! if needed, add ksh, tcsh, bash, etc
integer :: head_task = 0         ! def 0, but N-1 if reverse_task_layout true
logical :: print4status = .true. ! minimal messages for async4 handshake

character(len = 256) :: errstring, errstring1

! for broadcasts, pack small messages into larger ones.  remember that the
! byte size will be this count * 8 because we only communicate r8s.  (unless
! the code is compiled with r8 redefined as r4, in which case it's * 4).
integer, parameter :: PACKLIMIT1 = 8
integer, parameter :: PACKLIMIT2 = 512

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

! make local copy for send/recv.  was needed on an old, buggy version
! of the mpi libs but seems unneeded now. 
logical :: make_copy_before_sendrecv = .false.   ! should not be needed

! NAMELIST: change the following from .false. to .true. to enable
! the reading of this namelist.  This is the only place you need
! to make this change.
logical :: use_namelist = .false.

namelist /mpi_utilities_nml/ reverse_task_layout, all_tasks_print, &
                             verbose, async2_verbose, async4_verbose, &
                             shell_name, separate_node_sync, create_local_comm, &
                             make_copy_before_sendrecv

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

integer :: errcode, iunit
logical :: already

if ( module_initialized ) then
   ! return without calling the code below multiple times
   write(errstring, *) 'initialize_mpi_utilities has already been called'
   call error_handler(E_WARN,'initialize_mpi_utilities', errstring, source, revision, revdate)
   return
endif

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

call MPI_Comm_rank(MPI_COMM_WORLD, myrank, errcode)
if (errcode /= MPI_SUCCESS) then
   write(*, *) 'MPI_Comm_rank returned error code ', errcode
   call exit(-99)
endif

! give this a temporary initial value, in case we call the abort code.
! later down, we will dup the world comm and use a private comm for
! our communication.
my_local_comm = MPI_COMM_WORLD

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

if ( .not. module_initialized ) then
   ! Initialize the module with utilities
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

! this must come AFTER the standard utils are initialized.
! Read the DART namelist for the mpi_utilities.
if (use_namelist) then
   call find_namelist_in_file('input.nml', 'mpi_utilities_nml', iunit)
   read(iunit, nml = mpi_utilities_nml, iostat = errcode)
   call check_namelist_read(iunit, errcode, "mpi_utilities_nml")
else
   errstring = ' !must edit mpi_utilities/mpi_utilities_mod.f90 to enable this namelist'
   if (do_nml_file()) write(nmlfileunit, '(A)') trim(errstring)
   if (do_nml_term()) write(     *     , '(A)') trim(errstring)
endif

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=mpi_utilities_nml)
if (do_nml_term()) write(     *     , nml=mpi_utilities_nml)


! duplicate the world communicator to isolate us from any other user
! calls to MPI.  All subsequent mpi calls here will use the local communicator
! and not the global world comm.
if (create_local_comm) then
   call MPI_Comm_dup(MPI_COMM_WORLD, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Comm_dup returned error code ', errcode
      call error_handler(E_ERR,'initialize_mpi_utilities', errstring, &
                         source, revision, revdate)
   endif
endif

! find out who we are (0 to N-1).
call MPI_Comm_rank(my_local_comm, myrank, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_rank returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! number of tasks (if 10, returns 10.  task id numbers go from 0-9)
call MPI_Comm_size(my_local_comm, total_tasks, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_size returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
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
      call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring,source,revision,revdate)
   endif
endif

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
   call error_handler(E_MSG,'initialize_mpi_utilities: ',errstring,source,revision,revdate)
endif

end subroutine initialize_mpi_utilities

!-----------------------------------------------------------------------------

subroutine finalize_mpi_utilities(callfinalize, async)
 logical, intent(in), optional :: callfinalize
 integer, intent(in), optional :: async

! Shut down MPI cleanly.  This must be done before the program exits; on
! some implementations of MPI the final I/O flushes are not done until this
! is called.  The optional argument can prevent us from calling MPI_Finalize,
! so that user code can continue to use MPI after this returns.  Calling other
! routines in this file after calling finalize will invalidate your warranty.

integer :: errcode
logical :: dofinalize

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
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

! Release the private communicator we created at init time.
if (my_local_comm /= MPI_COMM_WORLD) then
   call MPI_Comm_free(my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Comm_free returned error code ', errcode
      call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
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
      call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
   endif
endif

! NO I/O after calling MPI_Finalize.  on some MPI implementations
! this can hang the job.

end subroutine finalize_mpi_utilities


!-----------------------------------------------------------------------------

function task_count()

! Return the total number of MPI tasks.  e.g. if the number of tasks is 4,
! it returns 4.  (The actual task numbers are 0-3.)

integer :: task_count

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_count', errstring, source, revision, revdate)
endif

task_count = total_tasks

end function task_count


!-----------------------------------------------------------------------------

function my_task_id()

! Return my unique task id.  Values run from 0 to N-1 (where N is the
! total number of MPI tasks.

integer :: my_task_id

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'my_task_id', errstring, source, revision, revdate)
endif

my_task_id = myrank

end function my_task_id


!-----------------------------------------------------------------------------

subroutine task_sync()

! Synchronize all tasks.  This subroutine does not return until all tasks
! execute this line of code.

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_sync', errstring, source, revision, revdate)
endif

if (verbose) write(*,*) "PE", myrank, ": waiting at MPI Barrier"
call MPI_Barrier(my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Barrier returned error code ', errcode
   call error_handler(E_ERR,'task_sync', errstring, source, revision, revdate)
endif

if (verbose) write(*,*) "PE", myrank, ": MPI Barrier released"

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

integer :: i, tag, errcode
integer :: itime(2)
integer :: itemcount, offset, nextsize
real(r8), allocatable :: tmpdata(:)

if (verbose) write(*,*) "PE", myrank, ": start of send_to "

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((dest_id < 0) .or. (dest_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "destination task id ", dest_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

itemcount = size(srcarray)

if (present(label)) then
   write(*,*) trim(label)//" PE", myrank, ": send_to itemsize ", itemcount, " dest ", dest_id
else if (verbose) then
   write(*,*) "PE", myrank, ": send_to itemsize ", itemcount, " dest ", dest_id
endif

! use my task id as the tag; unused at this point.
tag = myrank

if (make_copy_before_sendrecv) then
   if (itemcount > SNDRCV_MAXSIZE) then
      write(errstring, '(a,i12,a,i12)') "MPI msg contains ", itemcount, &
            " items; above single msg limit of ", SNDRCV_MAXSIZE
      write(errstring1, '(a)') "cannot make local copy; change 'make_copy_before_sendrecv' to .false. so msg can be sent in multiple chunks"
      call error_handler(E_ERR,'send_to', errstring, source, revision, revdate, &
                         text2=errstring1)
   endif

   allocate(tmpdata(itemcount))

   ! this copy should be unneeded, but on the intel fortran 9.0 compiler and mpich
   ! on one platform, calling this routine with an array section resulted in some
   ! apparently random memory corruption.  making a copy of the data into a local,
   ! contiguous buffer before send and receive fixed the corruption.  this should
   ! be examined at some later time for any performance impact.
   tmpdata = srcarray

   if (verbose) write(*,*) "PE", myrank, ": send_to alloctmp ", itemcount, " dest ", dest_id
   call MPI_Ssend(tmpdata, itemcount, datasize, dest_id, tag, &
                 my_local_comm, errcode)
else
   ! there are a few places in the code where we send/receive a full state vector.
   ! as these get really large, they may exceed the limits of the MPI library.
   ! break large messages up into smaller chunks if needed.  
   if (itemcount > SNDRCV_MAXSIZE) then
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
   
         call MPI_Ssend(srcarray(offset:offset+nextsize-1), nextsize, datasize, dest_id, tag, &
                       my_local_comm, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
            call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
         endif
         offset = offset + nextsize
      end do
   else
      if (verbose) write(*,*) "PE", myrank, ": send_to ", itemcount, " dest ", dest_id
      call MPI_Ssend(srcarray, size(srcarray), datasize, dest_id, tag, &
                    my_local_comm, errcode)
   endif
endif

if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   if (verbose) write(*,*) "PE", myrank, ": time present"
   call get_time(time, itime(1), itime(2))
   call MPI_Ssend(itime, 2, MPI_INTEGER, dest_id, tag*2, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
      call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
   endif
   if (verbose) write(*,*) "PE", myrank, ": sent time to ", dest_id
endif

if (make_copy_before_sendrecv) deallocate(tmpdata)


if (verbose) write(*,*) "PE", myrank, ": end of send_to "

end subroutine send_to


!-----------------------------------------------------------------------------

subroutine receive_from(src_id, destarray, time, label)
 integer, intent(in) :: src_id
 real(r8), intent(inout) :: destarray(:)
 type(time_type), intent(out), optional :: time
 character(len=*), intent(in), optional :: label

! Receive data into the destination array from the src task.
! If time is specified, it is received in a separate communications call.  
! This is a synchronous call; it will not return until the source has 
! sent the data.  If the send_to/receive_from calls are not paired correctly 
! the code will hang.

integer :: i, tag, errcode
integer :: itime(2)
integer :: status(MPI_STATUS_SIZE)
integer :: itemcount, offset, nextsize
real(r8), allocatable :: tmpdata(:)

if (verbose) write(*,*) "PE", myrank, ": start of receive_from "

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((src_id < 0) .or. (src_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "source task id ", src_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

itemcount = size(destarray)

if (present(label)) then
   write(*,*) trim(label)//" PE", myrank, ": receive_from itemsize ", itemcount, " src ", src_id
else if (verbose) then
   write(*,*) "PE", myrank, ": receive_from itemsize ", itemcount, " src ", src_id
endif

! send_to uses its own id as the tag.
tag = src_id


if (make_copy_before_sendrecv) then
   if (itemcount > SNDRCV_MAXSIZE) then
      write(errstring, '(a,i12,a,i12)') "MPI msg contains ", itemcount, &
            " items; above single msg limit of ", SNDRCV_MAXSIZE
      write(errstring1, '(a)') "cannot make local copy; change 'make_copy_before_sendrecv' to .false. so msg can be received in multiple chunks"
      call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate, &
                         text2=errstring1)
   endif
   allocate(tmpdata(itemcount))

   if (verbose) write(*,*) "PE", myrank, ": receive_from alloctmp ", itemcount, " src ", src_id

   ! call MPI to receive the data from the remote task
   call MPI_Recv(tmpdata, itemcount, datasize, src_id, MPI_ANY_TAG, &
                 my_local_comm, status, errcode)
else
   ! there are a few places in the code where we send/receive a full state vector.
   ! as these get really large, they may exceed the limits of the MPI library.
   ! break large messages up into smaller chunks if needed.  
   if (itemcount > SNDRCV_MAXSIZE) then
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
   
         call MPI_Recv(destarray(offset:offset+nextsize-1), nextsize, datasize, src_id, MPI_ANY_TAG, &
                       my_local_comm, status, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
            call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
         endif
         offset = offset + nextsize
      end do
   else
      if (verbose) write(*,*) "PE", myrank, ": receive_from ", itemcount,  " src ", src_id
      call MPI_Recv(destarray, itemcount, datasize, src_id, MPI_ANY_TAG, &
                 my_local_comm, status, errcode)
   endif

endif

if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

if (verbose) write(*,*) "PE", myrank, ": received from ", src_id

if (make_copy_before_sendrecv) then
   ! see comment in send_to() about why this code receives into a temp array
   ! and then does a contents copy into the actual destination.
   destarray = tmpdata
   if (verbose) write(*,*) "PE", myrank, ": copied to final dataarray"
endif

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   if (verbose) write(*,*) "PE", myrank, ": time present"
   call MPI_Recv(itime, 2, MPI_INTEGER, src_id, tag*2, &
                 my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
   endif
   if (itime(2) < 0) then
      write(errstring, '(a,i8)') 'seconds in settime were < 0; using 0'
      call error_handler(E_MSG,'receive_from', errstring, source, revision, revdate)
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
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 1D, 2D, 3D inputs?


subroutine transpose_array

! not implemented here yet.  will have arguments -- several of them.

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)
endif

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

integer :: itemcount, errcode, offset, nextsize

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

itemcount = size(array)

!if (verbose) then
!   if (myrank == root) write(*,*) "PE", myrank, ": bcast itemsize here ", itemcount
!endif
!if (verbose) write(*,*) "PE", myrank, ": bcast itemsize ", itemcount, " root ", root

! at least one user has run into a limit in the MPI implementation where you
! cannot broadcast too large an array.  if the size of this array is too large,
! broadcast it in chunks until all the data has been processed
if (itemcount > BCAST_MAXSIZE) then
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
      call MPI_Bcast(array(offset:offset+nextsize-1), nextsize, datasize, root, my_local_comm, errcode)
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Bcast returned error code ', errcode
         call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
      endif
      offset = offset + nextsize
   end do
else
   call MPI_Bcast(array, itemcount, datasize, root, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Bcast returned error code ', errcode
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif
endif

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

real(r8), allocatable :: localchunk(:)
integer :: srccount, leftover
integer :: i, tag, errcode
logical :: iamroot
integer :: status(MPI_STATUS_SIZE)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'array_distribute', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

iamroot = (root == myrank)
tag = 1

srccount = size(srcarray)

! TODO: right now this code does contig chunks only
! TODO: it should select on the 'how' argument
dstcount = srccount / total_tasks
leftover = srccount - (dstcount * total_tasks)
if (myrank == total_tasks-1) dstcount = dstcount + leftover


! idiotproofing, continued...
if (size(dstarray) < dstcount) then
   write(errstring, '(a,i8,a,i8)') "size of dstarray is", size(dstarray), & 
                      " but must be >= ", dstcount
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif
if (size(which) < dstcount) then
   write(errstring, '(a,i8,a,i8)') "size of which is", size(which), & 
                      " but must be >= ", dstcount
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

! TODO: this code is separate from the 'dstcount' computation because we
! need to test to be sure the user has passed us in arrays large enough to
! hold the data, but then this section needs to have a select (how) and set
! the corresponding index numbers accordingly.
which(1:dstcount) = (/ (i, i= myrank *dstcount, (myrank+1)*dstcount - 1) /)
if (size(which) > dstcount) which(dstcount+1:) = -1
   

if (.not.iamroot) then

   ! my task is receiving data.
   call MPI_Recv(dstarray, dstcount, datasize, root, MPI_ANY_TAG, &
                 my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif

else
   ! my task must send to everyone else and copy to myself.
   allocate(localchunk(dstcount), stat=errcode)  
   if (errcode /= 0) then
      write(errstring, *) 'allocation error of allocatable array'
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif

   do i=0, total_tasks-1
      ! copy correct data from srcarray to localchunk for each destination
      if (i == myrank) then
         ! this is my task, so do a copy from localchunk to dstarray
         dstarray(1:dstcount) = localchunk(1:dstcount)
      else
         ! call MPI to send the data to the remote task
         call MPI_Ssend(localchunk, dstcount, datasize, i, tag, &
                        my_local_comm, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
            call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
         endif
      endif
      tag = tag + 1
   enddo

   deallocate(localchunk, stat=errcode)  
   if (errcode /= 0) then
      write(errstring, *) 'deallocation error of allocatable array'
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif
endif
   
! set any additional space which wasn't filled with zeros.
if (size(dstarray) > dstcount) dstarray(dstcount+1:) = 0.0

end subroutine array_distribute

!-----------------------------------------------------------------------------
! DART-specific cover utilities
!-----------------------------------------------------------------------------

function iam_task0()

! Return .TRUE. if my local task id is 0, .FALSE. otherwise.
! (Task numbers in MPI start at 0, contrary to the rules of polite fortran.)

logical :: iam_task0

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'iam_task0', errstring, source, revision, revdate)
endif

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

! this must be paired with the same number of broadcast_recv()s on all 
! other tasks.  it will not return until all tasks in the communications 
! group have made the call.

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' matches my local task id.  also, these arrays are
! intent(in) here, but they call a routine which is intent(inout) so they
! must be the same here.
 
real(r8) :: packbuf1(PACKLIMIT1), packbuf2(PACKLIMIT2)
real(r8) :: local(5)
logical  :: doscalar, morethanone
integer  :: itemcount, sindex, eindex

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! simple idiotproofing
if (from /= myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "must be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! for relatively small array sizes, pack them into a single send/recv pair.
call countup(array1, array2, array3, array4, array5, &
             scalar1, scalar2, scalar3, scalar4, scalar5, &
             itemcount, morethanone, doscalar)

! try to use the smallest buffer possible
if (itemcount <= PACKLIMIT1 .and. morethanone) then

   call packit(packbuf1, array1, array2, array3, array4, array5, doscalar, &
                         scalar1, scalar2, scalar3, scalar4, scalar5)

   call array_broadcast(packbuf1, from)

else if (itemcount <= PACKLIMIT2 .and. morethanone) then

   call packit(packbuf2, array1, array2, array3, array4, array5, doscalar, &
                         scalar1, scalar2, scalar3, scalar4, scalar5)

   call array_broadcast(packbuf2, from)

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
subroutine broadcast_recv(from, array1, array2, array3, array4, array5, &
                          scalar1, scalar2, scalar3, scalar4, scalar5)
 integer, intent(in) :: from
 ! really only intent(out) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:)
 real(r8), intent(inout), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(inout), optional :: scalar1, scalar2, scalar3, scalar4, scalar5

! this must be paired with broadcast_send() on all other tasks, and it must 
! match exactly the number of args in the sending call.
! it will not return until all tasks in the communications group have
! made the call.

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' is not the same as my local task id.  these arrays are
! intent(out) here, but they call a routine which is intent(inout) so they
! must be the same here.

real(r8) :: packbuf1(PACKLIMIT1), packbuf2(PACKLIMIT2)
real(r8) :: local(5)
logical :: doscalar, morethanone
integer :: itemcount, sindex, eindex

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! simple idiotproofing
if (from == myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "cannot be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! for relatively small array sizes, pack them into a single send/recv pair.
call countup(array1, array2, array3, array4, array5, &
             scalar1, scalar2, scalar3, scalar4, scalar5, &
             itemcount, morethanone, doscalar)

if (itemcount <= PACKLIMIT1 .and. morethanone) then

   call array_broadcast(packbuf1, from)

   call unpackit(packbuf1, array1, array2, array3, array4, array5, doscalar, &
                           scalar1, scalar2, scalar3, scalar4, scalar5)

else if (itemcount <= PACKLIMIT2 .and. morethanone) then

   call array_broadcast(packbuf2, from)

   call unpackit(packbuf2, array1, array2, array3, array4, array5, doscalar, &
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
subroutine countup(array1, array2, array3, array4, array5, &
                   scalar1, scalar2, scalar3, scalar4, scalar5, &
                   numitems, morethanone, doscalar)
 real(r8), intent(in)           :: array1(:)
 real(r8), intent(in), optional :: array2(:), array3(:), array4(:), array5(:)
 real(r8), intent(in), optional :: scalar1, scalar2, scalar3, scalar4, scalar5
 integer,  intent(out)          :: numitems
 logical,  intent(out)          :: morethanone, doscalar

! figure out how many items are in the specified arrays, total.
! also note if there's more than a single array (array1) to send,
! and if there are any scalars specified.

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
subroutine sum_across_tasks(addend, sum)
 integer, intent(in) :: addend
 integer, intent(out) :: sum

 integer :: errcode
 integer :: localaddend(1), localsum(1)

! cover routine for MPI all-reduce

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'sum_across_tasks', errstring, source, revision, revdate)
endif

localaddend(1) = addend

if (verbose) write(*,*) "PE", myrank, ": Allreduce called"

call MPI_Allreduce(localaddend, localsum, 1, MPI_INTEGER, MPI_SUM, &
                   my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Allreduce returned error code ', errcode
   call error_handler(E_ERR,'sum_across_tasks', errstring, source, revision, revdate)
endif

sum = localsum(1)

end subroutine sum_across_tasks


!-----------------------------------------------------------------------------
! pipe-related utilities
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine block_task()

! block by reading a named pipe file until some other task
! writes a string into it.  this ensures the task is not
! spinning and using CPU cycles, but is asleep waiting in
! the kernel.   one subtlety with this approach is that even
! though named pipes are created in the filesystem, they are
! implemented in the kernel, so on a multiprocessor machine
! the write into the pipe file must occur on the same PE as
! the reader is waiting.  see the 'wakeup_filter' program for
! the MPI job which spreads out on all the PEs for this job
! and writes into the file from the correct PE.

character(len = 32) :: fifo_name, filter_to_model, model_to_filter, non_pipe
integer :: rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'block_task', errstring, source, revision, revdate)
endif

! FIXME: this should be mpi or a string other than filter (this is generic 
! mpi wrapper code, callable from programs other than filter.)
filter_to_model = 'filter_to_model.lock'
model_to_filter = 'model_to_filter.lock'
non_pipe = 'filter_to_model.file'

! the i5.5 format below will not handle task counts larger than this.
if (total_tasks > 99999) then
   write(errstring, *) 'cannot handle task counts > 99999'
   call error_handler(E_ERR,'block_task', errstring, source, revision, revdate)
endif

! make it so we only have to test 1 or 2 things here instead of 3
! when deciding whether to print status messages.
if (verbose) async4_verbose = .TRUE.

if ((myrank == head_task) .and. separate_node_sync) then

   if (async4_verbose) then
      write(*,*)  'checking master task host'
      rc = system('echo master task running on host `hostname`'//' '//char(0))
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose .or. print4status) write(*,*) 'MPI job telling script to advance model'
   rc = system('echo advance > '//trim(non_pipe)//' '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

endif

if ((myrank == head_task) .and. .not. separate_node_sync) then

   if (async4_verbose) then
      write(*,*)  'checking master task host'
      rc = system('echo master task running on host `hostname`'//' '//char(0))
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose .or. print4status) write(*,*) 'MPI job telling script to advance model'
   rc = system('echo advance > '//trim(filter_to_model)//' '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'MPI job now waiting to read from lock file'
   rc = system('cat < '//trim(model_to_filter)//'> /dev/null '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

else

! if you change this in any way, change the corresponding string in 
! restart_task() below.
   ! FIXME: this should be 'task_lock', since it's generic code beyond filter.
   write(fifo_name, '(a, i5.5)') "filter_lock", myrank
   
   if (async4_verbose) then
      write(*,*)  'checking slave task host'
      rc = system('echo '//trim(fifo_name)//' accessed from host `hostname`'//' '//char(0))
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose) write(*,*) 'removing any previous lock file: '//trim(fifo_name)
   rc = system('rm -f '//trim(fifo_name)//' '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'made fifo, named: '//trim(fifo_name)
   rc = system('mkfifo '//trim(fifo_name)//' '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'ready to read from lock file: '//trim(fifo_name)
   rc = system('cat < '//trim(fifo_name)//'> /dev/null '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

   if (async4_verbose) write(*,*) 'got response, removing lock file: '//trim(fifo_name)
   rc = system('rm -f '//trim(fifo_name)//' '//char(0))
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
subroutine restart_task()

! companion to block_task.  must be called by a different executable
! and it writes into the named pipes to restart the waiting task.

character(len = 32) :: fifo_name, model_to_filter
integer :: rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'restart_task', errstring, source, revision, revdate)
endif

! FIXME: ditto previous comment about using the string 'filter' here.
model_to_filter = 'model_to_filter.lock'

! the i5.5 format below will not handle task counts larger than this.
if (total_tasks > 99999) then
   write(errstring, *) 'cannot handle task counts > 99999'
   call error_handler(E_ERR,'block_task', errstring, source, revision, revdate)
endif

! make it so we only have to test 1 or 2 things here instead of 3
! when deciding whether to print status messages.
if (verbose) async4_verbose = .TRUE.

! process 0 (or N-1) is handled differently in the code.
if ((myrank == head_task) .and. .not. separate_node_sync) then

   if (async4_verbose) then
      rc = system('echo master task running on host `hostname`'//' '//char(0))
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose .or. print4status) write(*,*) 'script telling MPI job ok to restart'
   rc = system('echo restart > '//trim(model_to_filter)//' '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

else

   if (async4_verbose) write(*,*) 'waking up task id ', myrank

   ! FIXME: this should be 'task_lock', since it's generic code beyond filter.
   write(fifo_name,"(a,i5.5)") "filter_lock", myrank

   if (async4_verbose) then
      rc = system('echo '//trim(fifo_name)//' accessed from host `hostname`'//' '//char(0))
      if (rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc
   endif

   if (async4_verbose) write(*,*) 'ready to write to lock file: '//trim(fifo_name)
   rc = system('echo restart > '//trim(fifo_name)//' '//char(0))
   if (async4_verbose .and. rc /= 0) write(*, *) 'system command returned nonzero rc, ', rc

endif

end subroutine restart_task

!-----------------------------------------------------------------------------
subroutine finished_task(async)
 integer, intent(in) :: async

character(len = 32) :: fifo_name, filter_to_model, non_pipe
integer :: rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'restart_task', errstring, source, revision, revdate)
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
      rc = system('echo finished > '//trim(non_pipe)//' '//char(0))
   else
      rc = system('echo finished > '//trim(filter_to_model)//' '//char(0))
   endif

   
endif

end subroutine finished_task


!-----------------------------------------------------------------------------
!    * make_pipe()        Function that creates a named pipe (fifo), opens it,
!                         and returns the unit number.  Ok to call if the pipe
!                         already exists or is already open; it will skip
!                         those steps and just return the unit number.  The 
!                         name argument is used as a base and a filename
!                         in the form 'base.NNNN' is generated, where the N's
!                         are the MPI rank number, 0 padded.
!
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

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'make_pipe', errstring, source, revision, revdate)
endif

write(fname, "(a,i4.4)") trim(pipename)//".", myrank
!print *, "fname now = ", trim(fname)

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

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'destroy_pipe', errstring, source, revision, revdate)
endif

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
 
 integer :: ios

 read(iunit, *, iostat=ios) chardata
 if (ios /= 0) then
   write(errstring, '(a,i4)') 'failure trying to read from unit ', iunit
   call error_handler(E_ERR,'read_pipe', errstring, source, revision, revdate)
endif

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
 
 integer :: ios

 write(iunit, *, iostat=ios) chardata
 if (ios /= 0) then
   write(errstring, '(a,i4)') 'failure trying to write from unit ', iunit
   call error_handler(E_ERR,'write_pipe', errstring, source, revision, revdate)
endif


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

! allow code to test the theory that maybe the system call is
! not reentrant on some platforms.  if serialize is set and
! is true, do each call serially.

logical :: all_at_once
integer :: i, errcode, dummy(1)
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
      !shell_execute = system(trim(execute_string)//' '//char(0))
      shell_execute = system(trim(shell_name)//' '//trim(execute_string)//' '//char(0))
      if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

      return
   endif

   ! only one task at a time calls system, and all wait their turn by
   ! making each task wait for a message from the (N-1)th task.
   
   ! this is used only to signal; the value it contains is unused.
   dummy = 0

   if (myrank == 0) then

      ! my turn to execute
      shell_execute = system(trim(shell_name)//' '//trim(execute_string)//' '//char(0))
      if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

      if (total_tasks > 1) then
         ! tell next task it can continue
         call MPI_Send(dummy, 1, MPI_INTEGER, 1, 1, my_local_comm, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Send returned error code ', &
                                        errcode
            call error_handler(E_ERR,'shell_execute', errstring, source, &
                               revision, revdate)
         endif
      endif

   else if (myrank /= (total_tasks-1)) then
      ! wait for (me-1) to tell me it is my turn
      call MPI_Recv(dummy, 1, MPI_INTEGER, myrank-1, myrank, &
                    my_local_comm, status, errcode)
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
         call error_handler(E_ERR,'shell_execute', errstring, source, &
                            revision, revdate)
      endif

      ! my turn to execute
      shell_execute = system(trim(shell_name)//' '//trim(execute_string)//' '//char(0))
      if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

      ! and now tell (me+1) to go
      call MPI_Send(dummy, 1, MPI_INTEGER, myrank+1, myrank+1, my_local_comm, errcode)
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Send returned error code ', &
                                     errcode
         call error_handler(E_ERR,'shell_execute', errstring, source, &
                            revision, revdate)
      endif
   else
      ! last task, no one else to send to.
      call MPI_Recv(dummy, 1, MPI_INTEGER, myrank-1, myrank, &
                    my_local_comm, status, errcode)
      if (errcode /= MPI_SUCCESS) then
         write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
         call error_handler(E_ERR,'shell_execute', errstring, source, &
                            revision, revdate)
      endif

      ! my turn to execute
      shell_execute = system(trim(shell_name)//' '//trim(execute_string)//' '//char(0))
      if (async2_verbose) write(*,*) "PE", myrank, ": execution returns, rc = ", shell_execute

   endif
       

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
function get_dart_mpi_comm()
 integer :: get_dart_mpi_comm

! return our private communicator (or world, if no private created)
 get_dart_mpi_comm = my_local_comm

end function get_dart_mpi_comm

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Collect min and max on task. This is for adaptive_inflate_mod
subroutine reduce_min_max(minmax, task, global_val)

real(r8), intent(in)  :: minmax(2) ! min max on each task
integer,  intent(in)  :: task ! task to collect on
real(r8), intent(out) :: global_val(2) ! only concerned with this on task collecting result

integer :: errcode

! collect on task 
call mpi_reduce(minmax(1:1), global_val(1:1), 1, datasize, MPI_MIN,task, get_dart_mpi_comm(), errcode)
call mpi_reduce(minmax(2:2), global_val(2:2), 1, datasize, MPI_MAX, task, get_dart_mpi_comm(), errcode)

end subroutine reduce_min_max

!-----------------------------------------------------------------------------
! One sided communication

subroutine get_from_mean(owner, window, index, x)

integer,  intent(in)  :: owner  ! task in the window that owns the memory
integer,  intent(in)  :: window ! window object
integer,  intent(in)  :: index  ! index in the tasks memory
real(r8), intent(out) :: x ! result

integer(KIND=MPI_ADDRESS_KIND) :: target_disp
integer :: errcode

! Note to programmer: The data transfer is not guaranteed
! to have occured until the call to mpi_win_unlock. 
! => Don't do anything with x in between mpi_get and mpi_win_lock

! Note to programmer: openmpi 1.10.0 does not
! allow scalars in mpi calls. openmpi 1.10.1 fixes
! this.
target_disp = (index - 1)
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, window, errcode)
call mpi_get(x, 1, datasize, owner, target_disp, 1, datasize, window, errcode)
call mpi_win_unlock(owner, window, errcode)

end subroutine get_from_mean

!-----------------------------------------------------------------------------

subroutine get_from_fwd(owner, window, index, num_rows, x)

integer,  intent(in)  :: owner    ! task in the window that owns the memory
integer,  intent(in)  :: window   ! window object
integer,  intent(in)  :: index    ! index in the tasks memory
integer,  intent(in)  :: num_rows ! number of rows in the window
real(r8), intent(out) :: x(:)     ! result

integer(KIND=MPI_ADDRESS_KIND) :: target_disp
integer :: errcode

! Note to programmer: The data transfer is not guaranteed
! to have occured until the call to mpi_win_unlock. 
! => Don't do anything with x in between mpi_get and mpi_win_lock


target_disp = (index - 1)*num_rows
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, window, errcode)
call mpi_get(x, num_rows, datasize, owner, target_disp, num_rows, datasize, window, errcode)
call mpi_win_unlock(owner, window, errcode)

end subroutine get_from_fwd

!-----------------------------------------------------------------------------


end module mpi_utilities_mod

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! NOTE -- non-module code, so this subroutine can be called from the
!  utilities module, which this module uses (and cannot have circular refs)
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine exit_all(exit_code)
 use mpi_utilities_mod, only : get_dart_mpi_comm

 integer, intent(in) :: exit_code

! In case of error, call this instead of the fortran intrinsic exit().
! It will signal the other MPI tasks that something bad happened and they
! should also exit.

integer :: ierror

! if we made a local communicator, call abort on it.
! otherwise call abort on the world comm.

!print *, 'calling abort on comm ', get_dart_mpi_comm()
call MPI_Abort(get_dart_mpi_comm(),  exit_code, ierror)

! execution should never get here

end subroutine exit_all

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
