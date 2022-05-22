MODULE mpi_utilities_mod
========================

Overview
--------

This module provides subroutines which utilize the MPI (Message Passing Interface) parallel communications library. DART
does **not** require MPI; to compile without using MPI substitute the ``null_mpi_utilities_mod.f90`` file for this one.
That file contains the same module name and public entry points as this one but implements a serial version of all the
routines. However, to be able to run most larger models with a reasonable number of ensemble members (e.g. 30-100) MPI
will be needed.

Several DART executables can be compiled and run as either a serial program or a parallel program. Most work
directories in the DART distribution source tree have a ``quickbuild.sh`` script. To build DART without MPI use
``./quickbuild.sh nompi``.  No source code changes are required to switch between using mpi or no mpi.

A parallel program generally runs faster and requires less memory per CPU than the serial code. It requires an
implementation of the MPI library and run-time system to pass data between different nodes on a parallel cluster or
supercomputer. There is a lot of information about MPI on the web. See here for `an intro to MPI and parallel
programming <https://computing.llnl.gov/tutorials/mpi/>`__, and here for `downloads and technical
help <http://www.open-mpi.org>`__.

Most of the larger models need to be compiled and run with MPI because of limitations on total memory accessible by a
single executable. The smaller models (e.g. any of the Lorenz models) can generally be run as a serial program without
needing MPI.

The MPI distributions usually include a module named ``mpi`` which defines the public entry points and the types and
names of the routine arguments. However there are build-time options and older distributions which only supply an
``mpi.h`` include file. If you get a compile-time error about the mpi module being missing, edit the source code in
``mpi_utilities/mpi_utilities_mod.f90`` and comment out the ``use mpi`` line and comment in the ``include 'mpi.h'``
line. The 'use' line must be before the 'contains' line, while the 'include' line must be after, so do not move the
existing lines. Just comment them in or out depending on which one you need to use.

To preserve backwards compatibility this code does not require a namelist. However there is a namelist defined in the
source file which contains some useful run-time options. To enable it edit the source file in
``mpi_utilities/mpi_utilities_mod.f90`` and set ``use_namelist`` to .TRUE. and recompile. The code will then read the
namelist described below. Messages printed to the nml output log file will confirm whether the defaults are being used
or if the namelist is being read in.

Namelist
--------

The source code defines a namelist, but for backwards compatibility it is not read in unless the source code in
``mpi_utilities/mpi_utilities_mod.f90`` is edited, the module global variable ``use_namelist`` is changed from .FALSE.
to .TRUE., and then all executables are recompiled.

If enabled, this namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with
a slash '/'. Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.

::

   &mpi_utilities_nml
       reverse_task_layout        = .false.
       all_tasks_print            = .false.
       verbose                    = .false.
       async2_verbose             = .false.
       async4_verbose             = .false.
       shell_name                 = ''
       separate_node_sync         = .false.
       create_local_comm          = .true.
       make_copy_before_sendrecv  = .false.
      /

| 

.. container::

   +---------------------------+--------------------+-------------------------------------------------------------------+
   | Item                      | Type               | Description                                                       |
   +===========================+====================+===================================================================+
   | reverse_task_layout       | logical            | The synchronizing mechanism between the job script and the        |
   |                           |                    | parallel filter in async=4 mode relies on the script and task 0   |
   |                           |                    | running on the same node (in the same memory space if the nodes   |
   |                           |                    | have multiple processors). Some MPI implementations (OpenMPI      |
   |                           |                    | being the most commonly used one) lay the tasks out so that the   |
   |                           |                    | last task is on the same node as the script. If the async 4 model |
   |                           |                    | advance never starts but there are no error messages, try setting |
   |                           |                    | this to .TRUE. before running. See also the 'async4_verbose' flag |
   |                           |                    | below.                                                            |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | all_tasks_print           | logical            | In the parallel filter, informational messages only print from    |
   |                           |                    | task 0 to avoid N copies of the same messages. Error messages and |
   |                           |                    | warnings print no matter which task they occur in. If this        |
   |                           |                    | variable is set to true, even messages will print from all tasks. |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | verbose                   | logical            | USE WITH CAUTION! This flag enables debugging print messages for  |
   |                           |                    | every MPI call - sends, receives, barriers - and is very, very    |
   |                           |                    | verbose. In most cases the size of the output file will exceed    |
   |                           |                    | the filesystem limits or will cause the executable to run so      |
   |                           |                    | slowly that it will not be useful. However in small testcases     |
   |                           |                    | this can be useful to trace problems.                             |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | async2_verbose            | logical            | Print out messages about the handshaking between filter and the   |
   |                           |                    | advance model scripts when running in async=2 mode. Not anywhere  |
   |                           |                    | as verbose as the flag above; in most cases the output volume is  |
   |                           |                    | reasonable.                                                       |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | async4_verbose            | logical            | Print out messages about the handshaking between filter and the   |
   |                           |                    | run script when running in async=4 mode. Not anywhere as verbose  |
   |                           |                    | as the flag above; in most cases the output volume is reasonable. |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | shell_name                | character(len=129) | If running on compute nodes which do not have the expected        |
   |                           |                    | default shell for async=2 or async=4 mode, specify the full       |
   |                           |                    | pathname of the shell to execute the script. Not normally needed  |
   |                           |                    | on most systems we run on. (However, at least one type of Cray    |
   |                           |                    | system has this need.)                                            |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | separate_node_sync        | logical            | Not supported yet. Will use files to handshake between the filter |
   |                           |                    | executable and the run script in async=4 mode when the launch     |
   |                           |                    | script is not running on any of the same nodes as the filter      |
   |                           |                    | tasks.                                                            |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | create_local_comm         | logical            | The DART MPI routines normally create a separate local MPI        |
   |                           |                    | communicator instead of using MPI_COMM_WORLD. This keeps DART     |
   |                           |                    | communications separate from any other user code. To use the      |
   |                           |                    | default world communicator set this to .FALSE. . Normal use       |
   |                           |                    | should leave this true.                                           |
   +---------------------------+--------------------+-------------------------------------------------------------------+
   | make_copy_before_sendrecv | logical            | Workaround for old MPI bug. Should be .false.                     |
   +---------------------------+--------------------+-------------------------------------------------------------------+

| 

Other modules used
------------------

::

   types_mod
   utilities_mod
   time_manager_mod
   mpi  (or mpif.h if mpi module not available)

Public interfaces
-----------------

=============================== ========================
*use mpi_utilities_mod, only :* initialize_mpi_utilities
\                               finalize_mpi_utilities
\                               task_count
\                               my_task_id
\                               task_sync
\                               block_task
\                               restart_task
\                               array_broadcast
\                               send_to
\                               receive_from
\                               iam_task0
\                               broadcast_send
\                               broadcast_recv
\                               shell_execute
\                               sleep_seconds
\                               sum_across_tasks
\                               get_dart_mpi_comm
\                               exit_all
=============================== ========================

| 

.. container:: routine

   *call initialize_mpi_utilities( [progname] [, alternatename])*
   ::

      character(len=*), intent(in), optional :: progname
      character(len=*), intent(in), optional :: alternatename

.. container:: indent1

   Initializes the MPI library, creates a private communicator, stores the total number of tasks and the local task
   number for later use, and registers this module. This routine calls ``initialize_utilities()`` internally before
   returning, so the calling program need only call this one routine to initialize the DART internals.

   On some implementations of MPI (in particular some variants of MPICH) it is best to initialize MPI before any I/O is
   done from any of the parallel tasks, so this routine should be called as close to the process startup as possible.

   It is not an error to try to initialize the MPI library more than once. It is still necessary to call this routine
   even if the application itself has already initialized the MPI library. Thise routine creates a private communicator
   so internal communications are shielded from any other communication called outside the DART libraries.

   It is an error to call any of the other routines in this file before calling this routine.

   ================= ================================================================================
   ``progname``      If given, written to the log file to document which program is being started.
   ``alternatename`` If given, use this name as the log file instead of the default ``dart_log.out``.
   ================= ================================================================================

| 

.. container:: routine

   *call finalize_mpi_utilities( [callfinalize] [, async])*
   ::

      logical, intent(in), optional  :: callfinalize
      integer, intent(in), optional  :: async

.. container:: indent1

   Frees the local communicator, and shuts down the MPI library unless ``callfinalize`` is specified and is ``.FALSE.``.
   On some hardware platforms it is problematic to try to call print or write from the parallel tasks after finalize has
   been executed, so this should only be called immediately before the process is ready to exit. This routine does an
   ``MPI_Barrier()`` call before calling ``MPI_Finalize()`` to ensure all tasks are finished writing.

   If the application itself is using MPI the ``callfinalize`` argument can be used to defer closing the MPI library
   until the application does it itself. This routine does close the DART log file and releases the local communicator
   even if not calling MPI_Finalize, so no other DART routines which might generate output can be used after calling
   this routine.

   It is an error to call any of the other routines in this file after calling this routine.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``callfinalize`` | If false, do not call the ``MPI_Finalize()`` routine.                                            |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``async``        | If the model advance mode (selected by the async namelist value in the filter_nml section)       |
   |                  | requires any synchronization or actions at shutdown, this is done. Currently async=4 requires an |
   |                  | additional set of actions at shutdown time.                                                      |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = task_count()*
   ::

      integer         :: task_count

.. container:: indent1

   Returns the total number of MPI tasks this job was started with. Note that MPI task numbers start at 0, but this is a
   count. So a 4-task job will return 4 here, but the actual task numbers will be from 0 to 3.

   ======= ======================================
   ``var`` Total number of MPI tasks in this job.
   ======= ======================================

| 

.. container:: routine

   *var = my_task_id()*
   ::

      integer         :: my_task_id

.. container:: indent1

   Returns the local MPI task number. This is one of the routines in which all tasks can make the same function call but
   each returns a different value. The return can be useful in creating unique filenames or otherwise distinguishing
   resources which are not shared amongst tasks. MPI task numbers start at 0, so valid task id numbers for a 4-task job
   will be 0 to 3.

   ======= =============================
   ``var`` My unique MPI task id number.
   ======= =============================

| 

.. container:: routine

   *call task_sync()*

.. container:: indent1

   Synchronize tasks. This call does not return until all tasks have called this routine. This ensures all tasks have
   reached the same place in the code before proceeding. All tasks must make this call or the program will hang.

| 

.. container:: routine

   *call send_to(dest_id, srcarray [, time])*
   ::

      integer,                   intent(in) :: dest_id
      real(r8), dimension(:),    intent(in) :: srcarray
      type(time_type), optional, intent(in) :: time

.. container:: indent1

   Use the MPI library to send a copy of an array of data from one task to another task. The sending task makes this
   call; the receiving task must make a corresponding call to ``receive_from()``.

   If ``time`` is specified, it is also sent to the receiving task. The receiving call must match this sending call
   regarding this argument; if ``time`` is specified here it must also be specified in the receive; if not given here it
   cannot be given in the receive.

   The current implementation uses ``MPI_Ssend()`` which does a synchronous send. That means this routine will not
   return until the receiving task has called the receive routine to accept the data. This may be subject to change; MPI
   has several other non-blocking options for send and receive.

   ============ ======================================
   ``dest_id``  The MPI task id of the receiver.
   ``srcarray`` The data to be copied to the receiver.
   ``time``     If specified, send the time as well.
   ============ ======================================

   The send and receive subroutines must be used with care. These calls must be used in pairs; the sending task and the
   receiving task must make corresponding calls or the tasks will hang. Calling them with different array sizes will
   result in either a run-time error or a core dump. The optional time argument must either be given in both calls or in
   neither or one of the tasks will hang. (Executive summary: There are lots of ways to go wrong here.)

| 

.. container:: routine

   *call receive_from(src_id, destarray [, time])*
   ::

      integer, intent(in)                    :: src_id
      real(r8), dimension(:), intent(out)    :: destarray
      type(time_type), intent(out), optional :: time

.. container:: indent1

   Use the MPI library to receive a copy of an array of data from another task. The receiving task makes this call; the
   sending task must make a corresponding call to ``send_to()``. Unpaired calls to these routines will result in the
   tasks hanging.

   If ``time`` is specified, it is also received from the sending task. The sending call must match this receiving call
   regarding this argument; if ``time`` is specified here it must also be specified in the send; if not given here it
   cannot be given in the send.

   The current implementation uses ``MPI_Recv()`` which does a synchronous receive. That means this routine will not
   return until the data has arrived in this task. This may be subject to change; MPI has several other non-blocking
   options for send and receive.

   ============= ============================================================
   ``src_id``    The MPI task id of the sender.
   ``destarray`` The location where the data from the sender is to be placed.
   ``time``      If specified, receive the time as well.
   ============= ============================================================

   See the notes section of ``send_to()``.

| 

.. container:: routine

   *call exit_all(exit_code)*
   ::

      integer, intent(in)   :: exit_code

.. container:: indent1

   A replacement for calling the Fortran intrinsic ``exit``. This routine calls ``MPI_Abort()`` to kill all MPI tasks
   associated with this job. This ensures one task does not exit silently and leave the rest hanging. This is not the
   same as calling ``finalize_mpi_utilities()`` which waits for the other tasks to finish, flushes all messages, closes
   log files cleanly, etc. This call immediately and abruptly halts all tasks associated with this job.

   Depending on the MPI implementation and job control system, the exit code may or may not be passed back to the
   calling job script.

   ============= ====================
   ``exit_code`` A numeric exit code.
   ============= ====================

   This routine is now called from the standard error handler. To avoid circular references this is NOT a module
   routine. Programs which are compiled without the mpi code must now compile with the ``null_mpi_utilities_mod.f90``
   file to satisfy the call to this routine in the error handler.

| 

.. container:: routine

   *call array_broadcast(array, root)*
   ::

      real(r8), dimension(:), intent(inout) :: array
      integer, intent(in)                   :: root

.. container:: indent1

   All tasks must make this call together, but the behavior in each task differs depending on whether it is the ``root``
   or not. On the task which has a task id equal to ``root`` the contents of the array will be sent to all other tasks.
   On any task which has a task id *not* equal to ``root`` the array is the location where the data is to be received
   into. Thus ``array`` is intent(in) on root, and intent(out) on all other tasks.

   When this routine returns, all tasks will have the contents of the root array in their own arrays.

   ========= ===========================================================================================
   ``array`` Array containing data to send to all other tasks, or the location in which to receive data.
   ``root``  Task ID which will be the data source. All others are destinations.
   ========= ===========================================================================================

   This is another of the routines which must be called by all tasks. The MPI call used here is synchronous, so all
   tasks block here until everyone has called this routine.

| 

.. container:: routine

   *var = iam_task0()*
   ::

      logical                        :: iam_task0

.. container:: indent1

   Returns ``.TRUE.`` if called from the task with MPI task id 0. Returns ``.FALSE.`` in all other tasks. It is
   frequently the case that some code should execute only on a single task. This allows one to easily write a block
   surrounded by ``if (iam_task0()) then ...`` .

   ======= ===========================================================================
   ``var`` Convenience function to easily test and execute code blocks on task 0 only.
   ======= ===========================================================================

| 

.. container:: routine

   *call broadcast_send(from, array1 [, array2] [, array3] [, array4] [, array5] [, scalar1] [, scalar2] [, scalar3] [,
   scalar4] [, scalar5] )*
   ::

      integer, intent(in)                   :: from
      real(r8), dimension(:), intent(inout) :: array1
      real(r8), dimension(:), intent(inout), optional :: array2
      real(r8), dimension(:), intent(inout), optional :: array3
      real(r8), dimension(:), intent(inout), optional :: array4
      real(r8), dimension(:), intent(inout), optional :: array5
      real(r8), intent(inout), optional :: scalar1
      real(r8), intent(inout), optional :: scalar2
      real(r8), intent(inout), optional :: scalar3
      real(r8), intent(inout), optional :: scalar4
      real(r8), intent(inout), optional :: scalar5

.. container:: indent1

   Cover routine for ``array_broadcast()``. This call must be matched with the companion call ``broadcast_recv()``. This
   routine should only be called on the task which is the root of the broadcast; it will be the data source. All other
   tasks must call ``broadcast_recv()``. This routine sends up to 5 data arrays and 5 scalars in a single call. A common
   pattern in the DART filter code is sending 2 arrays, but other combinations exist. This routine ensures that ``from``
   is the same as the current task ID. The arguments to this call must be matched exactly in number and type with the
   companion call to ``broadcast_recv()`` or an error (or hang) will occur.

   In reality the data here are ``intent(in)`` only but this routine will be calling ``array_broadcast()`` internally
   and so must be ``intent(inout)`` to match.

   ========== ======================================================
   ``from``   Current task ID; the root task for the data broadcast.
   ``array1`` First data array to be broadcast.
   *array2*   If given, second data array to be broadcast.
   *array3*   If given, third data array to be broadcast.
   *array4*   If given, fourth data array to be broadcast.
   *array5*   If given, fifth data array to be broadcast.
   *scalar1*  If given, first data scalar to be broadcast.
   *scalar2*  If given, second data scalar to be broadcast.
   *scalar3*  If given, third data scalar to be broadcast.
   *scalar4*  If given, fourth data scalar to be broadcast.
   *scalar5*  If given, fifth data scalar to be broadcast.
   ========== ======================================================

   This is another of the routines which must be called consistently; only one task makes this call and all other tasks
   call the companion ``broadcast_recv`` routine. The MPI call used here is synchronous, so all tasks block until
   everyone has called one of these two routines.

| 

.. container:: routine

   *call broadcast_recv(from, array1 [, array2] [, array3] [, array4] [, array5] [, scalar1] [, scalar2] [, scalar3] [,
   scalar4] [, scalar5] )*
   ::

      integer, intent(in)                   :: from
      real(r8), dimension(:), intent(inout) :: array1
      real(r8), dimension(:), intent(inout), optional :: array2
      real(r8), dimension(:), intent(inout), optional :: array3
      real(r8), dimension(:), intent(inout), optional :: array4
      real(r8), dimension(:), intent(inout), optional :: array5
      real(r8), intent(inout), optional :: scalar1
      real(r8), intent(inout), optional :: scalar2
      real(r8), intent(inout), optional :: scalar3
      real(r8), intent(inout), optional :: scalar4
      real(r8), intent(inout), optional :: scalar5

.. container:: indent1

   Cover routine for ``array_broadcast()``. This call must be matched with the companion call ``broadcast_send()``. This
   routine must be called on all tasks which are *not* the root of the broadcast; the arguments specify the location in
   which to receive data from the root. (The root task should call ``broadcast_send()``.) This routine receives up to 5
   data arrays and 5 scalars in a single call. A common pattern in the DART filter code is receiving 2 arrays, but other
   combinations exist. This routine ensures that ``from`` is *not* the same as the current task ID. The arguments to
   this call must be matched exactly in number and type with the companion call to ``broadcast_send()`` or an error (or
   hang) will occur.

   In reality the data arrays here are ``intent(out)`` only but this routine will be calling ``array_broadcast()``
   internally and so must be ``intent(inout)`` to match.

   ========== ==================================================
   ``from``   The task ID for the data broadcast source.
   ``array1`` First array location to receive data into.
   *array2*   If given, second data array to receive data into.
   *array3*   If given, third data array to receive data into.
   *array4*   If given, fourth data array to receive data into.
   *array5*   If given, fifth data array to receive data into.
   *scalar1*  If given, first data scalar to receive data into.
   *scalar2*  If given, second data scalar to receive data into.
   *scalar3*  If given, third data scalar to receive data into.
   *scalar4*  If given, fourth data scalar to receive data into.
   *scalar5*  If given, fifth data scalar to receive data into.
   ========== ==================================================

   This is another of the routines which must be called consistently; all tasks but one make this call and exactly one
   other task calls the companion ``broadcast_send`` routine. The MPI call used here is synchronous, so all tasks block
   until everyone has called one of these two routines.

| 

.. container:: routine

   *call sum_across_tasks(addend, sum)*
   ::

      integer, intent(in)                   :: addend
      integer, intent(out)                  :: sum

.. container:: indent1

   All tasks call this routine, each with their own different ``addend``. The returned value in ``sum`` is the total of
   the values summed across all tasks, and is the same for each task.

   ========== ============================================
   ``addend`` Single input value per task to be summed up.
   ``sum``    The sum.
   ========== ============================================

   This is another of those calls which must be made from each task, and the calls block until this is so.

| 

.. container:: routine

   *call block_task()*

.. container:: indent1

   Create a named pipe (fifo) and read from it to block the process in such a way that it consumes no CPU time. Beware
   that once you put yourself to sleep you cannot wake yourself up. Some other MPI program must call restart_task() on
   the same set of processors the original program was distributed over.

   Even though fifos appear to be files, in reality they are implemented in the kernel. The write into the fifo must be
   executed on the same node as the read is pending on. See the man pages for the mkfifo(1) command for more details.

| 

.. container:: routine

   *call restart_task()*

.. container:: indent1

   Write into the pipe to restart the reading task. Note that this must be an entirely separate executable from the one
   which called block_task(), because it is asleep like Sleeping Beauty and cannot wake itself. See filter and
   wakeup_filter for examples of a program pair which uses these calls in async=4 mode.

   Even though fifos appear to be files, in reality they are implemented in the kernel. The write into the fifo must be
   executed on the same node as the read is pending on. See the man pages for the mkfifo(1) command for more details.

| 

.. container:: routine

   *call finished_task(async)*
   ::

      integer, intent(in) :: async

.. container:: indent1

   For async=4 and task id = 0, write into the main filter-to-script fifo to tell the run script that filter is exiting.
   Does nothing else otherwise.

   Even though fifos appear to be files, in reality they are implemented in the kernel. The write into the fifo must be
   executed on the same node as the read is pending on. See the man pages for the mkfifo(1) command for more details.

| 

.. container:: routine

   *rc = shell_execute()*
   ::

      integer                       :: shell_execute
      character(len=*), intent(in)  :: execute_string
      logical, intent(in), optional :: serialize

.. container:: indent1

   Wrapper routine around the system() library function to execute shell level commands from inside the Fortran program.
   Will wait for the command to execute and will return the error code. 0 means ok, any other number indicates error.

   +--------------------+------------------------------------------------------------------------------------------------+
   | ``rc``             | Return code from the shell exit after the command has been executed.                           |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``execute_string`` | Command to be executed by the shell.                                                           |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``serialize``      | If specified and if .TRUE. run the command from each PE in turn, waiting for each to complete  |
   |                    | before beginning the next. The default is .FALSE. and does not require that all tasks call     |
   |                    | this routine. If given and .TRUE. then all tasks must make this call.                          |
   +--------------------+------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call sleep_seconds(naplength)*
   ::

      real(r8), intent(in) :: naplength

.. container:: indent1

   Wrapper routine for the sleep command. Argument is a real in seconds. Some systems have different lower resolutions
   for the minimum time it will sleep. This routine can round up to even seconds if a smaller than 1.0 time is given.

   ============= ===========================================
   ``naplength`` Number of seconds to sleep as a real value.
   ============= ===========================================

   The amount of time this routine will sleep is not precise and might be in units of whole seconds on some platforms.

| 

.. container:: routine

   *comm = get_dart_mpi_comm()*
   ::

      integer    :: get_dart_mpi_comm

.. container:: indent1

   This code creates a private communicator for DART MPI calls, in case other code in the executable is using the world
   communicator. This routine returns the private communicator. If it is called before the internal setup work is
   completed it returns MPI_COMM_WORLD. If it is called before MPI is initialized, it returns 0.

   ======== ==============================
   ``comm`` The private DART communicator.
   ======== ==============================

| 

Files
-----

-  mpi module or
-  mpif.h

Depending on the implementation of MPI, the library routines are either defined in an include file (``mpif.h``) or in a
proper Fortran 90 module (``use mpi``). If it is available the module is preferred; it allows for better argument
checking and optional arguments support in the MPI library calls.

References
----------

-  MPI: The Complete Reference; Snir, Otto, Huss-Lederman, Walker, Dongarra; MIT Press, 1996, ISBN 0-262-69184-1
-  `http://www-unix.mcs.anl.gov/mpi/ <http://www-unix.mcs.anl.gov/mpi/>`__

Private components
------------------

N/A
