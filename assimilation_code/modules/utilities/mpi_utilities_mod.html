<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module mpi_utilities_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE mpi_utilities_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>This module provides subroutines which utilize the MPI (Message
Passing Interface) parallel communications library. DART does
<strong>not</strong> require MPI; to compile without using MPI
substitute the <em class="file">null_mpi_utilities_mod.f90</em>
file for this one. That file contains the same module name and
public entry points as this one but implements a serial version of
all the routines. However, to be able to run most larger models
with a reasonable number of ensemble members (e.g. 30-100) MPI will
be needed.</p>
<p>The main DART executable <em class="file">filter</em> can be
compiled and run as either a serial program or a parallel program.
Most work directories in the DART distribution source tree have a
<em class="file">quickbuild.csh</em> script which can take a
<em class="code">-mpi</em> or a <em class="code">-nompi</em> flag.
This flag changes the list of files to be compiled to use either
the module which uses the MPI library or the one which makes no MPI
calls. No source code changes are required to switch between the
two options.</p>
<p>A parallel program generally runs faster and requires less
memory per CPU than the serial code. It requires an implementation
of the MPI library and run-time system to pass data between
different nodes on a parallel cluster or supercomputer. There is a
lot of information about MPI on the web. See here for <a href=
"https://computing.llnl.gov/tutorials/mpi/">an intro to MPI and
parallel programming</a>, and here for <a href=
"http://www.open-mpi.org">downloads and technical help</a>.</p>
<p>Most of the larger models need to be compiled and run with MPI
because of limitations on total memory accessible by a single
executable. The smaller models (e.g. any of the Lorenz models) can
generally be run as a serial program without needing MPI.</p>
<p>The MPI distributions usually include a module named <em class=
"file">mpi</em> which defines the public entry points and the types
and names of the routine arguments. However there are build-time
options and older distributions which only supply an <em class=
"file">mpi.h</em> include file. If you get a compile-time error
about the mpi module being missing, edit the source code in
<em class="code">mpi_utilities/mpi_utilities_mod.f90</em> and
comment out the <em class="code">use mpi</em> line and comment in
the <em class="code">include 'mpi.h'</em> line. The 'use' line must
be before the 'contains' line, while the 'include' line must be
after, so do not move the existing lines. Just comment them in or
out depending on which one you need to use.</p>
<p>To preserve backwards compatibility this code does not require a
namelist. However there is a namelist defined in the source file
which contains some useful run-time options. To enable it edit the
source file in <em class=
"code">mpi_utilities/mpi_utilities_mod.f90</em> and set <em class=
"code">use_namelist</em> to .TRUE. and recompile. The code will
then read the namelist described below. Messages printed to the nml
output log file will confirm whether the defaults are being used or
if the namelist is being read in.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>The source code defines a namelist, but for backwards
compatibility it is not read in unless the source code in
<em class="code">mpi_utilities/mpi_utilities_mod.f90</em> is
edited, the module global variable <em class=
"code">use_namelist</em> is changed from .FALSE. to .TRUE., and
then all executables are recompiled.</p>
<p>If enabled, this namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;mpi_utilities_nml
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
</pre></div>
<br>
<br>
<div>
<table border="0" cellpadding="10" width="100%" summary=
'namelist description'>
<thead align="left">
<tr>
<th>Item</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>reverse_task_layout</td>
<td>logical</td>
<td>The synchronizing mechanism between the job script and the
parallel filter in async=4 mode relies on the script and task 0
running on the same node (in the same memory space if the nodes
have multiple processors). Some MPI implementations (OpenMPI being
the most commonly used one) lay the tasks out so that the last task
is on the same node as the script. If the async 4 model advance
never starts but there are no error messages, try setting this to
.TRUE. before running. See also the 'async4_verbose' flag
below.</td>
</tr>
<tr>
<td>all_tasks_print</td>
<td>logical</td>
<td>In the parallel filter, informational messages only print from
task 0 to avoid N copies of the same messages. Error messages and
warnings print no matter which task they occur in. If this variable
is set to true, even messages will print from all tasks.</td>
</tr>
<tr>
<td>verbose</td>
<td>logical</td>
<td>USE WITH CAUTION! This flag enables debugging print messages
for every MPI call - sends, receives, barriers - and is very, very
verbose. In most cases the size of the output file will exceed the
filesystem limits or will cause the executable to run so slowly
that it will not be useful. However in small testcases this can be
useful to trace problems.</td>
</tr>
<tr>
<td>async2_verbose</td>
<td>logical</td>
<td>Print out messages about the handshaking between filter and the
advance model scripts when running in async=2 mode. Not anywhere as
verbose as the flag above; in most cases the output volume is
reasonable.</td>
</tr>
<tr>
<td>async4_verbose</td>
<td>logical</td>
<td>Print out messages about the handshaking between filter and the
run script when running in async=4 mode. Not anywhere as verbose as
the flag above; in most cases the output volume is reasonable.</td>
</tr>
<tr>
<td>shell_name</td>
<td>character(len=129)</td>
<td>If running on compute nodes which do not have the expected
default shell for async=2 or async=4 mode, specify the full
pathname of the shell to execute the script. Not normally needed on
most systems we run on. (However, at least one type of Cray system
has this need.)</td>
</tr>
<tr>
<td>separate_node_sync</td>
<td>logical</td>
<td>Not supported yet. Will use files to handshake between the
filter executable and the run script in async=4 mode when the
launch script is not running on any of the same nodes as the filter
tasks.</td>
</tr>
<tr>
<td>create_local_comm</td>
<td>logical</td>
<td>The DART MPI routines normally create a separate local MPI
communicator instead of using MPI_COMM_WORLD. This keeps DART
communications separate from any other user code. To use the
default world communicator set this to .FALSE. . Normal use should
leave this true.</td>
</tr>
<tr>
<td>make_copy_before_sendrecv</td>
<td>logical</td>
<td>Workaround for old MPI bug. Should be .false.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
<!--==================================================================-->
<!--==================================================================-->
 <a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
time_manager_mod
mpi  (or mpif.h if mpi module not available)
</pre>
<!--==================================================================-->
<!-- Public entities                                                  -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table>
<tr>
<td><em class="call">use mpi_utilities_mod, only :</em></td>
<td><a href=
"#initialize_mpi_utilities">initialize_mpi_utilities</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#finalize_mpi_utilities">finalize_mpi_utilities</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#task_count">task_count</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#my_task_id">my_task_id</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#task_sync">task_sync</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#block_task">block_task</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#restart_task">restart_task</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#array_broadcast">array_broadcast</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#send_to">send_to</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#receive_from">receive_from</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#iam_task0">iam_task0</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#broadcast_send">broadcast_send</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#broadcast_recv">broadcast_recv</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#shell_execute">shell_execute</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#sleep_seconds">sleep_seconds</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#sum_across_tasks">sum_across_tasks</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_dart_mpi_comm">get_dart_mpi_comm</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#exit_all">exit_all</a></td>
</tr>
</table>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
<a name="initialize_mpi_utilities" id=
"initialize_mpi_utilities"></a><br>
<div class="routine"><em class="call">call
initialize_mpi_utilities( <em class="optionalcode">[progname]</em>
<em class="optionalcode">[, alternatename]</em>)</em>
<pre>
character(len=*), intent(in), optional :: <em class=
"code">progname</em>
character(len=*), intent(in), optional :: <em class=
"code">alternatename</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initializes the MPI library, creates a private communicator,
stores the total number of tasks and the local task number for
later use, and registers this module. This routine calls <em class=
"code">initialize_utilities()</em> internally before returning, so
the calling program need only call this one routine to initialize
the DART internals.</p>
<p>On some implementations of MPI (in particular some variants of
MPICH) it is best to initialize MPI before any I/O is done from any
of the parallel tasks, so this routine should be called as close to
the process startup as possible.</p>
<p>It is not an error to try to initialize the MPI library more
than once. It is still necessary to call this routine even if the
application itself has already initialized the MPI library. Thise
routine creates a private communicator so internal communications
are shielded from any other communication called outside the DART
libraries.</p>
<p>It is an error to call any of the other routines in this file
before calling this routine.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">progname  </em></td>
<td>If given, written to the log file to document which program is
being started.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">alternatename  </em></td>
<td>If given, use this name as the log file instead of the default
<em class="code">dart_log.out</em>.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="finalize_mpi_utilities" id=
"finalize_mpi_utilities"></a><br>
<div class="routine"><em class="call">call finalize_mpi_utilities(
<em class="optionalcode">[callfinalize]</em> <em class=
"optionalcode">[, async]</em>)</em>
<pre>
logical, intent(in), optional  :: <em class=
"code">callfinalize</em>
integer, intent(in), optional  :: <em class="code">async</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Frees the local communicator, and shuts down the MPI library
unless <em class="code">callfinalize</em> is specified and is
<em class="code">.FALSE.</em>. On some hardware platforms it is
problematic to try to call print or write from the parallel tasks
after finalize has been executed, so this should only be called
immediately before the process is ready to exit. This routine does
an <em class="code">MPI_Barrier()</em> call before calling
<em class="code">MPI_Finalize()</em> to ensure all tasks are
finished writing.</p>
<p>If the application itself is using MPI the <em class=
"code">callfinalize</em> argument can be used to defer closing the
MPI library until the application does it itself. This routine does
close the DART log file and releases the local communicator even if
not calling MPI_Finalize, so no other DART routines which might
generate output can be used after calling this routine.</p>
<p>It is an error to call any of the other routines in this file
after calling this routine.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">callfinalize  </em></td>
<td>If false, do not call the <em class="code">MPI_Finalize()</em>
routine.</td>
</tr>
<tr>
<td valign="top"><em class="code">async  </em></td>
<td>If the model advance mode (selected by the async namelist value
in the filter_nml section) requires any synchronization or actions
at shutdown, this is done. Currently async=4 requires an additional
set of actions at shutdown time.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF FUNCTION =====================-->
 <a name="task_count" id="task_count"></a><br>
<div class="routine"><em class="call">var = task_count()</em>
<pre>
integer         :: <em class="code">task_count</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the total number of MPI tasks this job was started with.
Note that MPI task numbers start at 0, but this is a count. So a
4-task job will return 4 here, but the actual task numbers will be
from 0 to 3.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Total number of MPI tasks in this job.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF FUNCTION =====================-->
 <a name="my_task_id" id="my_task_id"></a><br>
<div class="routine"><em class="call">var = my_task_id()</em>
<pre>
integer         :: <em class="code">my_task_id</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the local MPI task number. This is one of the routines
in which all tasks can make the same function call but each returns
a different value. The return can be useful in creating unique
filenames or otherwise distinguishing resources which are not
shared amongst tasks. MPI task numbers start at 0, so valid task id
numbers for a 4-task job will be 0 to 3.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var   </em></td>
<td>My unique MPI task id number.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="task_sync" id="task_sync"></a><br>
<div class="routine"><em class="call">call task_sync()</em></div>
<div class="indent1"><!-- Description -->
<p>Synchronize tasks. This call does not return until all tasks
have called this routine. This ensures all tasks have reached the
same place in the code before proceeding. All tasks must make this
call or the program will hang.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="send_to" id="send_to"></a><br>
<div class="routine"><em class="call">call send_to(dest_id,
srcarray <em class="optionalcode">[, time]</em>)</em>
<pre>
integer,                   intent(in) :: <em class=
"code">dest_id</em>
real(r8), dimension(:),    intent(in) :: <em class=
"code">srcarray</em>
type(time_type), optional, intent(in) :: <em class=
"optionalcode">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Use the MPI library to send a copy of an array of data from one
task to another task. The sending task makes this call; the
receiving task must make a corresponding call to <em class=
"code"><a href="#receive_from">receive_from()</a></em>.</p>
<p>If <em class="code">time</em> is specified, it is also sent to
the receiving task. The receiving call must match this sending call
regarding this argument; if <em class="code">time</em> is specified
here it must also be specified in the receive; if not given here it
cannot be given in the receive.</p>
<p>The current implementation uses <em class=
"code">MPI_Ssend()</em> which does a synchronous send. That means
this routine will not return until the receiving task has called
the receive routine to accept the data. This may be subject to
change; MPI has several other non-blocking options for send and
receive.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">dest_id</em></td>
<td>The MPI task id of the receiver.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">srcarray   </em></td>
<td>The data to be copied to the receiver.</td>
</tr>
<tr>
<td valign="top"><em class="code">time</em></td>
<td>If specified, send the time as well.</td>
</tr>
</table>
<br>
<p>The send and receive subroutines must be used with care. These
calls must be used in pairs; the sending task and the receiving
task must make corresponding calls or the tasks will hang. Calling
them with different array sizes will result in either a run-time
error or a core dump. The optional time argument must either be
given in both calls or in neither or one of the tasks will hang.
(Executive summary: There are lots of ways to go wrong here.)</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="receive_from" id="receive_from"></a><br>
<div class="routine"><em class="call">call receive_from(src_id,
destarray <em class="optionalcode">[, time]</em>)</em>
<pre>
integer, intent(in)                    :: <em class=
"code">src_id</em>
real(r8), dimension(:), intent(out)    :: <em class=
"code">destarray</em>
type(time_type), intent(out), optional :: <em class=
"code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Use the MPI library to receive a copy of an array of data from
another task. The receiving task makes this call; the sending task
must make a corresponding call to <em class="code"><a href=
"#send_to">send_to()</a></em>. Unpaired calls to these routines
will result in the tasks hanging.</p>
<p>If <em class="code">time</em> is specified, it is also received
from the sending task. The sending call must match this receiving
call regarding this argument; if <em class="code">time</em> is
specified here it must also be specified in the send; if not given
here it cannot be given in the send.</p>
<p>The current implementation uses <em class="code">MPI_Recv()</em>
which does a synchronous receive. That means this routine will not
return until the data has arrived in this task. This may be subject
to change; MPI has several other non-blocking options for send and
receive.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">src_id   </em></td>
<td>The MPI task id of the sender.</td>
</tr>
<tr>
<td valign="top"><em class="code">destarray   </em></td>
<td>The location where the data from the sender is to be
placed.</td>
</tr>
<tr>
<td valign="top"><em class="code">time   </em></td>
<td>If specified, receive the time as well.</td>
</tr>
</table>
<br>
<p>See the notes section of <em class="code"><a href=
"#send_to">send_to()</a></em>.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="exit_all" id="exit_all"></a><br>
<div class="routine"><em class="call">call exit_all(exit_code)</em>
<pre>
integer, intent(in)   :: <em class="code">exit_code</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>A replacement for calling the Fortran intrinsic <em class=
"code">exit</em>. This routine calls <em class=
"code">MPI_Abort()</em> to kill all MPI tasks associated with this
job. This ensures one task does not exit silently and leave the
rest hanging. This is not the same as calling <em class=
"code"><a href=
"#finalize_mpi_utilities">finalize_mpi_utilities()</a></em> which
waits for the other tasks to finish, flushes all messages, closes
log files cleanly, etc. This call immediately and abruptly halts
all tasks associated with this job.</p>
<p>Depending on the MPI implementation and job control system, the
exit code may or may not be passed back to the calling job
script.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">exit_code   </em></td>
<td>A numeric exit code.</td>
</tr>
</table>
<br>
<p>This routine is now called from the standard error handler. To
avoid circular references this is NOT a module routine. Programs
which are compiled without the mpi code must now compile with the
<em class="file">null_mpi_utilities_mod.f90</em> file to satisfy
the call to this routine in the error handler.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="array_broadcast" id="array_broadcast"></a><br>
<div class="routine"><em class="call">call array_broadcast(array,
root)</em>
<pre>
real(r8), dimension(:), intent(inout) :: <em class=
"code">array</em>
integer, intent(in)                   :: <em class="code">root</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>All tasks must make this call together, but the behavior in each
task differs depending on whether it is the <em class=
"code">root</em> or not. On the task which has a task id equal to
<em class="code">root</em> the contents of the array will be sent
to all other tasks. On any task which has a task id <em>not</em>
equal to <em class="code">root</em> the array is the location where
the data is to be received into. Thus <em class="code">array</em>
is intent(in) on root, and intent(out) on all other tasks.</p>
<p>When this routine returns, all tasks will have the contents of
the root array in their own arrays.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">array   </em></td>
<td>Array containing data to send to all other tasks, or the
location in which to receive data.</td>
</tr>
<tr>
<td valign="top"><em class="code">root   </em></td>
<td>Task ID which will be the data source. All others are
destinations.</td>
</tr>
</table>
<br>
<p>This is another of the routines which must be called by all
tasks. The MPI call used here is synchronous, so all tasks block
here until everyone has called this routine.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="iam_task0" id="iam_task0"></a><br>
<div class="routine"><em class="call">var = iam_task0()</em>
<pre>
logical                        :: <em class="code">iam_task0</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns <em class="code">.TRUE.</em> if called from the task
with MPI task id 0. Returns <em class="code">.FALSE.</em> in all
other tasks. It is frequently the case that some code should
execute only on a single task. This allows one to easily write a
block surrounded by <em class="code">if (iam_task0()) then ...</em>
.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var   </em></td>
<td>Convenience function to easily test and execute code blocks on
task 0 only.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="broadcast_send" id="broadcast_send"></a><br>
<div class="routine"><em class="call">call broadcast_send(from,
array1 <em class="optionalcode">[, array2]</em> <em class=
"optionalcode">[, array3]</em> <em class="optionalcode">[,
array4]</em> <em class="optionalcode">[, array5]</em> <em class=
"optionalcode">[, scalar1]</em> <em class="optionalcode">[,
scalar2]</em> <em class="optionalcode">[, scalar3]</em> <em class=
"optionalcode">[, scalar4]</em> <em class="optionalcode">[,
scalar5]</em> )</em>
<pre>
integer, intent(in)                   :: <em class="code">from</em>
real(r8), dimension(:), intent(inout) :: <em class=
"code">array1</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array2</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array3</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array4</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array5</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar1</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar2</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar3</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar4</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar5</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Cover routine for <em class="code"><a href=
"#array_broadcast">array_broadcast()</a></em>. This call must be
matched with the companion call <em class="code"><a href=
"#broadcast_recv">broadcast_recv()</a></em>. This routine should
only be called on the task which is the root of the broadcast; it
will be the data source. All other tasks must call <em class=
"code">broadcast_recv()</em>. This routine sends up to 5 data
arrays and 5 scalars in a single call. A common pattern in the DART
filter code is sending 2 arrays, but other combinations exist. This
routine ensures that <em class="code">from</em> is the same as the
current task ID. The arguments to this call must be matched exactly
in number and type with the companion call to <em class=
"code"><a href="#broadcast_recv">broadcast_recv()</a></em> or an
error (or hang) will occur.</p>
<p>In reality the data here are <em class="code">intent(in)</em>
only but this routine will be calling <em class=
"code">array_broadcast()</em> internally and so must be <em class=
"code">intent(inout)</em> to match.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">from   </em></td>
<td>Current task ID; the root task for the data broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="code">array1   </em></td>
<td>First data array to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array2
  </em></td>
<td>If given, second data array to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array3
  </em></td>
<td>If given, third data array to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array4
  </em></td>
<td>If given, fourth data array to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array5
  </em></td>
<td>If given, fifth data array to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar1
  </em></td>
<td>If given, first data scalar to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar2
  </em></td>
<td>If given, second data scalar to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar3
  </em></td>
<td>If given, third data scalar to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar4
  </em></td>
<td>If given, fourth data scalar to be broadcast.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar5
  </em></td>
<td>If given, fifth data scalar to be broadcast.</td>
</tr>
</table>
<p>This is another of the routines which must be called
consistently; only one task makes this call and all other tasks
call the companion <em class="code">broadcast_recv</em> routine.
The MPI call used here is synchronous, so all tasks block until
everyone has called one of these two routines.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="broadcast_recv" id="broadcast_recv"></a><br>
<div class="routine"><em class="call">call broadcast_recv(from,
array1 <em class="optionalcode">[, array2]</em> <em class=
"optionalcode">[, array3]</em> <em class="optionalcode">[,
array4]</em> <em class="optionalcode">[, array5]</em> <em class=
"optionalcode">[, scalar1]</em> <em class="optionalcode">[,
scalar2]</em> <em class="optionalcode">[, scalar3]</em> <em class=
"optionalcode">[, scalar4]</em> <em class="optionalcode">[,
scalar5]</em> )</em>
<pre>
integer, intent(in)                   :: <em class="code">from</em>
real(r8), dimension(:), intent(inout) :: <em class=
"code">array1</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array2</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array3</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array4</em>
real(r8), dimension(:), intent(inout), optional :: <em class=
"optionalcode">array5</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar1</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar2</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar3</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar4</em>
real(r8), intent(inout), optional :: <em class=
"optionalcode">scalar5</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Cover routine for <em class="code"><a href=
"#array_broadcast">array_broadcast()</a></em>. This call must be
matched with the companion call <em class="code"><a href=
"#broadcast_send">broadcast_send()</a></em>. This routine must be
called on all tasks which are <em>not</em> the root of the
broadcast; the arguments specify the location in which to receive
data from the root. (The root task should call <em class=
"code">broadcast_send()</em>.) This routine receives up to 5 data
arrays and 5 scalars in a single call. A common pattern in the DART
filter code is receiving 2 arrays, but other combinations exist.
This routine ensures that <em class="code">from</em> is
<em>not</em> the same as the current task ID. The arguments to this
call must be matched exactly in number and type with the companion
call to <em class="code"><a href=
"#broadcast_recv">broadcast_send()</a></em> or an error (or hang)
will occur.</p>
<p>In reality the data arrays here are <em class=
"code">intent(out)</em> only but this routine will be calling
<em class="code">array_broadcast()</em> internally and so must be
<em class="code">intent(inout)</em> to match.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">from   </em></td>
<td>The task ID for the data broadcast source.</td>
</tr>
<tr>
<td valign="top"><em class="code">array1   </em></td>
<td>First array location to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array2
  </em></td>
<td>If given, second data array to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array3
  </em></td>
<td>If given, third data array to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array4
  </em></td>
<td>If given, fourth data array to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">array5
  </em></td>
<td>If given, fifth data array to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar1
  </em></td>
<td>If given, first data scalar to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar2
  </em></td>
<td>If given, second data scalar to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar3
  </em></td>
<td>If given, third data scalar to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar4
  </em></td>
<td>If given, fourth data scalar to receive data into.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">scalar5
  </em></td>
<td>If given, fifth data scalar to receive data into.</td>
</tr>
</table>
<br>
<p>This is another of the routines which must be called
consistently; all tasks but one make this call and exactly one
other task calls the companion <em class="code">broadcast_send</em>
routine. The MPI call used here is synchronous, so all tasks block
until everyone has called one of these two routines.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="sum_across_tasks" id="sum_across_tasks"></a><br>
<div class="routine"><em class="call">call sum_across_tasks(addend,
sum)</em>
<pre>
integer, intent(in)                   :: <em class=
"code">addend</em>
integer, intent(out)                  :: <em class="code">sum</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>All tasks call this routine, each with their own different
<em class="code">addend</em>. The returned value in <em class=
"code">sum</em> is the total of the values summed across all tasks,
and is the same for each task.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">addend   </em></td>
<td>Single input value per task to be summed up.</td>
</tr>
<tr>
<td valign="top"><em class="code">sum   </em></td>
<td>The sum.</td>
</tr>
</table>
<br>
<p>This is another of those calls which must be made from each
task, and the calls block until this is so.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="block_task" id="block_task"></a><br>
<div class="routine"><em class="call">call block_task()</em></div>
<div class="indent1"><!-- Description -->
<p>Create a named pipe (fifo) and read from it to block the process
in such a way that it consumes no CPU time. Beware that once you
put yourself to sleep you cannot wake yourself up. Some other MPI
program must call restart_task() on the same set of processors the
original program was distributed over.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
<br>
<p>Even though fifos appear to be files, in reality they are
implemented in the kernel. The write into the fifo must be executed
on the same node as the read is pending on. See the man pages for
the mkfifo(1) command for more details.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="restart_task" id="restart_task"></a><br>
<div class="routine"><em class="call">call
restart_task()</em></div>
<div class="indent1"><!-- Description -->
<p>Write into the pipe to restart the reading task. Note that this
must be an entirely separate executable from the one which called
block_task(), because it is asleep like Sleeping Beauty and cannot
wake itself. See filter and wakeup_filter for examples of a program
pair which uses these calls in async=4 mode.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
<br>
<p>Even though fifos appear to be files, in reality they are
implemented in the kernel. The write into the fifo must be executed
on the same node as the read is pending on. See the man pages for
the mkfifo(1) command for more details.</p>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="finished_task" id="finished_task"></a><br>
<div class="routine"><em class="call">call
finished_task(async)</em>
<pre>
integer, intent(in) :: <em class="code">async</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>For async=4 and task id = 0, write into the main
filter-to-script fifo to tell the run script that filter is
exiting. Does nothing else otherwise.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
<br>
<p>Even though fifos appear to be files, in reality they are
implemented in the kernel. The write into the fifo must be executed
on the same node as the read is pending on. See the man pages for
the mkfifo(1) command for more details.</p>
</div>
<br>
<!--===================== DESCRIPTION OF FUNCTION =====================-->
 <a name="shell_execute" id="shell_execute"></a><br>
<div class="routine"><em class="call">rc = shell_execute()</em>
<pre>
integer                       :: <em class=
"code">shell_execute</em>
character(len=*), intent(in)  :: <em class=
"code">execute_string</em>
logical, intent(in), optional :: <em class="code">serialize</em>

</pre></div>
<div class="indent1"><!-- Description -->
<p>Wrapper routine around the system() library function to execute
shell level commands from inside the Fortran program. Will wait for
the command to execute and will return the error code. 0 means ok,
any other number indicates error.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">rc   </em></td>
<td>Return code from the shell exit after the command has been
executed.</td>
</tr>
<tr>
<td valign="top"><em class="code">execute_string
  </em></td>
<td>Command to be executed by the shell.</td>
</tr>
<tr>
<td valign="top"><em class="code">serialize   </em></td>
<td>If specified and if .TRUE. run the command from each PE in
turn, waiting for each to complete before beginning the next. The
default is .FALSE. and does not require that all tasks call this
routine. If given and .TRUE. then all tasks must make this
call.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF SUBROUTINE =====================-->
 <a name="sleep_seconds" id="sleep_seconds"></a><br>
<div class="routine"><em class="call">call
sleep_seconds(naplength)</em>
<pre>
real(r8), intent(in) :: <em class="code">naplength</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Wrapper routine for the sleep command. Argument is a real in
seconds. Some systems have different lower resolutions for the
minimum time it will sleep. This routine can round up to even
seconds if a smaller than 1.0 time is given.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">naplength   </em></td>
<td>Number of seconds to sleep as a real value.</td>
</tr>
</table>
<br>
<p>The amount of time this routine will sleep is not precise and
might be in units of whole seconds on some platforms.</p>
</div>
<br>
<!--===================== DESCRIPTION OF FUNCTION =====================-->
 <a name="get_dart_mpi_comm" id="get_dart_mpi_comm"></a><br>
<div class="routine"><em class="call">comm =
get_dart_mpi_comm()</em>
<pre>
integer    :: <em class="code">get_dart_mpi_comm</em>

</pre></div>
<div class="indent1"><!-- Description -->
<p>This code creates a private communicator for DART MPI calls, in
case other code in the executable is using the world communicator.
This routine returns the private communicator. If it is called
before the internal setup work is completed it returns
MPI_COMM_WORLD. If it is called before MPI is initialized, it
returns 0.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">comm   </em></td>
<td>The private DART communicator.</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>mpi module or</li>
<li>mpif.h</li>
</ul>
<p>Depending on the implementation of MPI, the library routines are
either defined in an include file (<em class="code">mpif.h</em>) or
in a proper Fortran 90 module (<em class="code">use mpi</em>). If
it is available the module is preferred; it allows for better
argument checking and optional arguments support in the MPI library
calls.</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>MPI: The Complete Reference; Snir, Otto, Huss-Lederman, Walker,
Dongarra; MIT Press, 1996, ISBN 0-262-69184-1</li>
<li><a href="http://www-unix.mcs.anl.gov/mpi/"><em class=
"code">http://www-unix.mcs.anl.gov/mpi/</em></a></li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<p>If MPI returns an error, the DART error handler is called with
the numeric error code it received from MPI. See any of the MPI
references for an up-to-date list of error codes.</p>
<p>After printing to the standard output and log files, the DART
error handler calls the <em class="code">exit_all()</em> routine
which calls <em class="code">MPI_Abort()</em> to make sure all
tasks exit and the entire job does not hang if only one task has an
error.</p>
<!--==================================================================-->
<!-- Describe all known bugs                                          -->
<!--==================================================================-->
<a name="KnownBugs" id="KnownBugs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
