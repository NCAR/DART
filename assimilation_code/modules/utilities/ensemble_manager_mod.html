<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module ensemble_manager_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE ensemble_manager_mod</h1>
<table border="0" summary="dart logo" cellpadding="5">
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
<p>Manages storage and a number of operations for multiple copies
of a vector. The most obvious use is to manage ensembles of model
state vectors. In this case, the number of copies stored for each
state vector element is the ensemble size plus one or more
additional copies like the mean, variance, associated inflation
values, etc. The ensemble_manager provides routines to compute the
mean and variance of a subset of the copies, to track the time
associated with the copies, and to write and read restart files.
Most importantly, it provides a capability to do transposes between
two storage representations of an ensemble. In one representation,
each process stores all copies of a subset of the state variables
while in the other, each process stores all of the state variables
for a subset of copies. The ensemble manager is also used to manage
ensembles of observation priors and quality control and ensembles
of forward observation operator error status.</p>
<p>The ensemble manager interacts strongly with the multiple
process capability of the Message Passing Interface (MPI)
libraries. It is used to partition the data so each MPI process
stores only a subset of the copies and variables, dividing the data
as evenly as possible across the processes. At no time during the
execution does any one process have to store the entire dataset for
all ensemble members (unless running in serial mode without MPI, or
if running with 1 MPI task).</p>
<p>The ensemble manager is set of general purpose data management
routines. For run-time efficiency, the derived type information is
not marked private which means other modules can directly
manipulate the data arrays. However it means much care must be
taken to access the most recently updated representation of the
data, either the copies or variables arrays.</p>
<p>A set of sanity check routines have been added to track the last
modified version of the data: the copies array or the vars array.
Before directly reading or writing these arrays call one of the
'prepare' routines to indicate what kind of data access you are
about to make. If the most recently updated data is not as expected
an error message will occur. After the direct access if the
following operations detect that the data they are operating on is
not the most recently updated they will print an error message.
Routines inside the ensemble manager that alter the copies or vars
will set the state automatically so these routines are only
necessary to call if you are directly accessing the copies or vars
arrays from outside the ensemble manager.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;ensemble_manager_nml
   layout                      = 1
   tasks_per_node              = 1
   communication_configuration = 1
   debug                       = .false.
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
<td>layout</td>
<td>integer</td>
<td>Determines the logical process (PE) layout across MPI tasks. 1
is PE = MPI task. 2 is a round-robin layout around the nodes.
Layout 2 results in a more even usage of memory across nodes. This
may allow you to run with a larger state vector without hitting the
memory limit of the node. It may give a slight (5%) increase in
performance, but this is machine dependent. It has no effect on
serial runs.</td>
</tr>
<tr>
<td>tasks_per_node</td>
<td>integer</td>
<td>The number of MPI tasks per hardware node is generally fixed
when a batch job is submitted. This namelist item tells the
ensemble manager what the user selected at that time. Once a
program is running the code has no control to change how MPI tasks
are assigned to physical CPUs. This number is used only if layout =
2, and it allows the code spread high-memory-use PEs to different
hardware nodes by assigning them in a round-robin order. The job
will still run if this number does not match the real
"tasks_per_node" at the hardware level, but it may run out of
memory if the mismatch causes multiple high-memory-use tasks to be
run on the same node.</td>
</tr>
<tr>
<td>communication_configuration</td>
<td>integer</td>
<td>For most users, the default value of 1 is the best choice.
However there are multiple strategies for the internal MPI
communication patterns (see *Note below). Values from 1 to 4 select
different options; try the various options to see if one might be
faster than the others.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>If true print debugging information.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<p><i>*Note about MPI communication flags:</i><br>
The communication_configuration flags select various combinations
of the internal settings for use_copy2var_send_loop and
use_var2copy_rec_loop. These flags change the order of the MPI send
and MPI receives in the the routines all_copies_to_all_vars and
all_vars_to_all_copies. The figures below show the data transferred
between tasks for an 80 member ensemble. The left figure is using
96 tasks, the right figure is using 512 tasks. As the number of
tasks increases, the 'all to all' data transfer becomes a 'some to
all, all to some' transfer and the order of MPI send and MPI
receives becomes increasingly important. The default values give a
performance advantage as the number of tasks becomes much greater
than the the ensemble size. However, for small numbers of tasks,
i.e. less than the ensemble size, changing the default values may
improve performance.</p>
<div>
<table width="100%" summary='communication patterns'>
<tr>
<td><a href="../../../docs/images/comm_pattern96.png"><img src=
"../../../docs/images/comm_pattern96.png" width="400" alt=
"communication pattern"></a></td>
<td><a href="../../../docs/images/comm_pattern512.png"><img src=
"../../../docs/images/comm_pattern512.png" width="400" alt=
"communication pattern 2"></a></td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
 <a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
assim_model_mod
time_manager_mod
random_seq_mod
mpi_utilities_mod
sort_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table summary='public routine list'>
<tr>
<td><em class="call">use ensemble_manager_mod, only :</em></td>
<td><a href="#init_ensemble_manager">init_ensemble_manager</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_ensemble_restart">read_ensemble_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#write_ensemble_restart">write_ensemble_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_copy">get_copy</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#put_copy">put_copy</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#broadcast_copy">broadcast_copy</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_ensemble_time">set_ensemble_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_ensemble_time">get_ensemble_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#end_ensemble_manager">end_ensemble_manager</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#duplicate_ens">duplicate_ens</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_my_num_copies">get_my_num_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_my_copies">get_my_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_my_num_vars">get_my_num_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_my_vars">get_my_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_copy_owner_index">get_copy_owner_index</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_var_owner_index">get_var_owner_index</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#all_vars_to_all_copies">all_vars_to_all_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#all_copies_to_all_vars">all_copies_to_all_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#compute_copy_mean">compute_copy_mean</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#compute_copy_mean_sd">compute_copy_mean_sd</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#compute_copy_mean_var">compute_copy_mean_var</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#prepare_to_write_to_vars">prepare_to_write_to_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#prepare_to_write_to_copies">prepare_to_write_to_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#prepare_to_read_from_vars">prepare_to_read_from_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#prepare_to_read_from_copies">prepare_to_read_from_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#prepare_to_update_vars">prepare_to_update_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#prepare_to_update_copies">prepare_to_update_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#print_ens_handle">print_ens_handle</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#map_pe_to_task">map_pe_to_task</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#map_task_to_pe">map_task_to_pe</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--=================== DESCRIPTION OF A LOCAL TYPE ==================-->
<a name="ensemble_type" id="ensemble_type"></a><br>
<div class="type">
<pre>
<em class="call">type ensemble_type</em>
   !DIRECT ACCESS INTO STORAGE IS ALLOWED; BE CAREFUL
   integer :: num_copies
   integer :: num_vars
   integer :: my_num_copies
   integer :: my_num_vars
   integer, pointer :: my_copies(:)
   integer, pointer :: my_vars(:)
   ! Storage in next line is to be used when each PE has all copies of subset of vars
   real(r8), pointer :: copies(:, :)  ! Dimensioned (num_copies, my_num_vars)
   ! Storage on next line is used when each PE has subset of copies of all vars
   real(r8), pointer :: vars(:, :)    ! Dimensioned (num_vars, my_num_copies)
   ! Time is only related to var complete
   type(time_type), pointer :: time(:)
   integer :: distribution_type
   integer :: valid     ! copies modified last, vars modified last, both same
   integer :: id_num
   integer, allocatable :: task_to_pe_list(:) ! List of tasks
   integer, allocatable :: pe_to_task_list(:) ! List of tasks
   ! Flexible my_pe, layout_type which allows different task layouts for different ensemble handles
   integer :: my_pe
   integer :: layout_type
end type ensemble_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Provides a handle for an ensemble that manages copies of a
vector. For efficiency, the type internals are not private and
direct access to the storage arrays is used throughout DART.</p>
<table border="0" cellpadding="3" width="100%" summary=
'derived type description'>
<thead align="left">
<tr>
<th>Component</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>num_copies</td>
<td>Global number of copies of the vector.</td>
</tr>
<tr>
<td>num_vars</td>
<td>Global number of elements (variables) in the vector.</td>
</tr>
<tr>
<td>my_num_copies</td>
<td>Number of copies stored by this process.</td>
</tr>
<tr>
<td>my_num_vars</td>
<td>Number of variables stored by this process.</td>
</tr>
<tr>
<td>my_copies</td>
<td>Dimensioned to size my_num_copies. Contains a list of the
global indices of copies stored by this process.</td>
</tr>
<tr>
<td>my_vars</td>
<td>Dimensioned to size my_num_vars. Contains a list of the global
indices of variables stored by this process.</td>
</tr>
<tr>
<td>copies</td>
<td>Dimensioned (num_copies, my_num_vars). Storage for all copies
of variables stored by this process.</td>
</tr>
<tr>
<td>vars</td>
<td>Dimensioned (num_vars, my_num_copies). Storage for all
variables of copies stored by this process.</td>
</tr>
<tr>
<td>time</td>
<td>Dimensioned my_num_copies. A time_type that stores time
associated with a given copy of the vector.</td>
</tr>
<tr>
<td>distribution_type</td>
<td>Does nothing at present. Can be used for future releases to
control the layout of different copies and variables in
storage.</td>
</tr>
<tr>
<td>valid</td>
<td>Flag to track whether the copies array has the most recently
updated data, the vars array is most recently modified, or if both
the arrays have identical data, like after a transpose.</td>
</tr>
<tr>
<td>id_num</td>
<td>Internal number unique to each ensemble handle, used for
debugging purposes.</td>
</tr>
<tr>
<td>task_to_pe_list</td>
<td>Mapping from MPI task number to logical Processing Element (PE)
number. Enables different assignment of MPI tasks to PEs. If the
number of MPI tasks is larger than the number of copies of the
vector, when the ensemble is var complete then the first N MPI
tasks have allocated 'vars' arrays and the remaining ones do not.
Assigning the MPI tasks round-robin to multi-processor nodes can
make the memory usage more uniform across nodes, which may allow
more MPI tasks per node than the standard layout.</td>
</tr>
<tr>
<td>pe_to_task_list</td>
<td>Logical PE to MPI task mapping. See above for more
description.</td>
</tr>
<tr>
<td>my_pe</td>
<td>The logical PE number for the MPI task.</td>
</tr>
<tr>
<td>layout_type</td>
<td>Controls the mapping type between MPI tasks and PEs. Currently
type 1 is the standard layout (one-to-one mapping) and type 2 is a
round-robin mapping where each node gets a task in turn before
assigning a second task to each node, until all tasks are
assigned.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_ensemble_manager" id=
"init_ensemble_manager"></a><br>
<div class="routine"><em class="call">call
init_ensemble_manager(ens_handle, num_copies, num_vars <em class=
"optionalcode">[, distribution_type_in]</em> <em class=
"optionalcode">[, layout_type]</em>)</em>
<pre>
type(ensemble_type), intent(out) :: <em class=
"code">ens_handle</em>
integer,             intent(in)  :: <em class=
"code">num_copies</em>
integer,             intent(in)  :: <em class="code">num_vars</em>
integer, optional,   intent(in)  :: <em class=
"optionalcode">distribution_type_in</em>
integer, optional,   intent(in)  :: <em class=
"optionalcode">layout_type</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initializes an instance of an ensemble. Storage is allocated and
the size descriptions in the ensemble_type are initialized.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for the ensemble being initialized</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies</em></td>
<td>Number of copies of vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_vars</em></td>
<td>Number of variables in the vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">distribution_type_in</em></td>
<td>Controls layout of storage on PEs. Currently only option 1 is
supported.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">layout_type</em></td>
<td>Controls layout of MPI tasks on PEs. Type 1 is the default,
where MPI tasks are assigned to PEs on a one-to-one basis. Type 2
is a round-robin assignment where each node gets one task before
the nodes are assigned a second task. If running with more MPI
tasks than <em class="code">num_copies</em>, this can result in a
more uniform usage of memory across the nodes.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_ensemble_restart" id=
"read_ensemble_restart"></a><br>
<div class="routine"><em class="call">call
read_ensemble_restart(ens_handle, start_copy, end_copy,
start_from_restart, file_name <em class=
"optionalcode">[, init_time]</em> <em class=
"optionalcode">[, force_single_file]</em>)</em>
<pre>
type(ensemble_type),       intent(inout) :: <em class=
"code">ens_handle</em>
integer,                   intent(in)    :: <em class=
"code">start_copy</em>
integer,                   intent(in)    :: <em class=
"code">end_copy</em>
logical,                   intent(in)    :: <em class=
"code">start_from_restart</em>
character(len=*),          intent(in)    :: <em class=
"code">file_name</em>
type(time_type), optional, intent(in)    :: <em class=
"optionalcode">init_time</em>
logical, optional,         intent(in)    :: <em class=
"optionalcode">force_single_file</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Read in a set of copies of a vector from file <em class=
"code">file_name</em>. The copies read are place into global copies
start_copy:end_copy in the ens_handle. If start_from_restart is
false, then only a single copy of the vector is read from the file
and then it is perturbed using routines in assim_model_mod to
generate the required number of copies. The read can be from a
single file that contains all needed copies or from a different
file for each copy. This choice is controlled by the namelist entry
single_restart_file_in. However, the optional argument
force_single_file forces the read to be from a single file if it is
present and true. This is used for ensembles that contain the
inflation values for state space inflation. If multiple files are
to be read, the file names are generated by appending integers to
the input file_name. If the input is a single file all reads are
done sequentially by process 0 and then shipped to the PE that
stores that copy. If the input is multiple files each MPI task
reads the copies it stores directly and independently.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle of ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">start_copy</em></td>
<td>Global index of first of continguous set of copies to be
read.</td>
</tr>
<tr>
<td valign="top"><em class="code">end_copy</em></td>
<td>Global index of last of contiguous set of copies to be read,
copies(start_copy:end_copy).</td>
</tr>
<tr>
<td valign="top"><em class="code">start_from_restart</em></td>
<td>If true, read all copies from file. If false, read one copy and
perturb to get required number.</td>
</tr>
<tr>
<td valign="top"><em class="code">file_name</em></td>
<td>Name of file from which to read.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">init_time</em></td>
<td>If present, set time of all copies read to this value.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">force_single_file</em></td>
<td>If present and true, force the read to be from a single file
which contains all copies.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_ensemble_restart" id=
"write_ensemble_restart"></a><br>
<div class="routine"><em class="call">call
write_ensemble_restart(ens_handle, file_name, start_copy, end_copy
<em class="optionalcode">[, force_single_file]</em>)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
character(len=*),    intent(in)    :: <em class=
"code">file_name</em>
integer,             intent(in)    :: <em class=
"code">start_copy</em>
integer,             intent(in)    :: <em class=
"code">end_copy</em>
logical, optional,   intent(in)    :: <em class=
"optionalcode">force_single_file</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes a set of copies of a vector to file file_name. The copies
written are from global copies start_copy:end_copy in the
ens_handle. The write can be to a single file or to a different
file for each copy. This choice is controlled by the namelist entry
single_restart_file_out. However, the optional argument
force_single_file forces the write to be to a single file if it is
present and true. This is used for ensembles that contain the
inflation values for state space inflation. If multiple files are
to be written, the file names are generated by appending integers
to the input file_name. If the output is a single file all copies
are shipped from the PE that stores that copy to process 0, and
then written out sequentially. If the output is to multiple files
each MPI task writes the copies it stores directly and
independently.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td><em class="code">ens_handle</em></td>
<td>Handle of ensemble.</td>
</tr>
<tr>
<td><em class="code">file_name</em></td>
<td>Name of file from which to read.</td>
</tr>
<tr>
<td valign="top"><em class="code">start_copy</em></td>
<td>Global index of first of continguous set of copies to be
written.</td>
</tr>
<tr>
<td valign="top"><em class="code">end_copy</em></td>
<td>Global index of last of contiguous set of copies to be written,
copies(start_copy:end_copy).</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">force_single_file</em></td>
<td>If present and true, force the write to be to a single file
which contains all copies.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_copy" id="get_copy"></a><br>
<div class="routine"><em class="call">call get_copy(receiving_pe,
ens_handle, copy, vars <em class=
"optionalcode">[, mtime]</em>)</em>
<pre>
integer,                   intent(in)  :: <em class=
"code">receiving_pe</em>
type(ensemble_type),       intent(in)  :: <em class=
"code">ens_handle</em>
integer,                   intent(in)  :: <em class=
"code">copy</em>
real(r8), dimension(:),    intent(out) :: <em class=
"code">vars</em>
type(time_type), optional, intent(out) :: <em class=
"optionalcode">mtime</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Retrieves a copy of the state vector, indexed by the global
index copy. The process that is to receive the copy is receiving_pe
and the copy is returned in the one dimensional array vars. The
time of the copy is also returned if mtime is present. This is
generally used for operations, like IO, that require a single
processor to do things with the entire state vector. Data is only
returned in vars on the receiving PE; vars on all other PEs is
unset.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">receiving_pe</em></td>
<td>This process ends up with the requested copy of the state
vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">copy</em></td>
<td>The global index of the copy of the state vector that is to be
retrieved.</td>
</tr>
<tr>
<td valign="top"><em class="code">vars</em></td>
<td>One dimensional array in which the requested copy of the state
vector is returned. Data is only returned in vars on the receiving
PE; vars on all other PEs is unset.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">mtime</em></td>
<td>If present returns the time of the requested copy.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="put_copy" id="put_copy"></a><br>
<div class="routine"><em class="call">call put_copy(sending_pe,
ens_handle, copy, vars <em class=
"optionalcode">[, mtime]</em>)</em>
<pre>
integer,                   intent(in)    :: <em class=
"code">sending_pe</em>
type(ensemble_type),       intent(inout) :: <em class=
"code">ens_handle</em>
integer,                   intent(in)    :: <em class=
"code">copy</em>
real(r8), dimension(:),    intent(in)    :: <em class=
"code">vars</em>
type(time_type), optional, intent(in)    :: <em class=
"optionalcode">mtime</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sends a state vector, in vars, from the given process to the
process storing the global index copy. The time of the copy is also
sent if mtime is present. This is generally used for operations,
like IO, that require a single processor to do things with the
entire state vector. For instance, if a single process reads in a
state vector, it can be shipped to the storing process by this
subroutine. Only the data in vars on the sending PE is processed;
vars on all other PEs is ignored.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">sending_pe</em></td>
<td>This process sends the copy of the state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">copy</em></td>
<td>The global index of the copy of the state vector that is to be
sent.</td>
</tr>
<tr>
<td valign="top"><em class="code">vars</em></td>
<td>One dimensional array in which the requested copy of the state
vector is located. Only the data in vars on the sending PE is
processed; vars on all other PEs is ignored.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">mtime</em></td>
<td>If present send the time of the copy.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="broadcast_copy" id="broadcast_copy"></a><br>
<div class="routine"><em class="call">call
broadcast_copy(ens_handle, copy, arraydata)</em>
<pre>
type(ensemble_type),    intent(in)   :: <em class=
"code">ens_handle</em>
integer,                intent(in)   :: <em class="code">copy</em>
real(r8), dimension(:), intent(out)  :: <em class=
"code">arraydata</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Finds which PE has the global index copy and broadcasts that
copy to all PEs. <em class="code">arraydata</em> is an output on
all PEs, even on the PE which is the owner if it is separate
storage from the vars array in the ensemble handle. This is a
collective routine, which means it must be called by all processes
in the job.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">copy</em></td>
<td>The global index of the copy of the state vector that is to be
sent.</td>
</tr>
<tr>
<td valign="top"><em class="code">arraydata</em></td>
<td>One dimensional array into which the requested copy of the
state vector will be copied on all PEs, including the sending
PE.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_ensemble_time" id="set_ensemble_time"></a><br>
<div class="routine"><em class="call">call
set_ensemble_time(ens_handle, indx, mtime)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer,             intent(in)    :: <em class="code">indx</em>
type(time_type),     intent(in)    :: <em class="code">mtime</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the time of a copy to the given value. <em class=
"code">indx</em> in this case is the local copy number for a
specific task. <!-- Contrast this with 
<a href="#set_copy_time">set_copy_time()</a> which uses global copy numbers. -->
<a href="#get_copy_owner_index">get_copy_owner_index()</a> can be
called to see if you are the owning task for a given global copy
number, and to get the local index number for that copy.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">indx</em></td>
<td>The local index of the copy of the state vector that is to be
set.</td>
</tr>
<tr>
<td valign="top"><em class="code">mtime</em></td>
<td>The time to set for this copy.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_ensemble_time" id="get_ensemble_time"></a><br>
<div class="routine"><em class="call">call
get_ensemble_time(ens_handle, indx, mtime)</em>
<pre>
type(ensemble_type), intent(in)   :: <em class=
"code">ens_handle</em>
integer,             intent(in)   :: <em class="code">indx</em>
type(time_type),     intent(out)  :: <em class="code">mtime</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Get the time associated with a copy. <em class="code">indx</em>
in this case is the local copy number for a specific task. <a href=
"#get_copy_owner_index">get_copy_owner_index()</a> can be called to
see if you are the owning task for a given global copy number, and
to get the local index number for that copy. 
<!-- Contrast this with 
<a href="#get_copy_time">get_copy_time()</a> which uses global copy numbers. --></p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">indx</em></td>
<td>The local index of the copy to retrieve the time from.</td>
</tr>
<tr>
<td valign="top"><em class="code">mtime</em></td>
<td>The returned time value.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_ensemble_manager" id="end_ensemble_manager"></a><br>
<div class="routine"><em class="call">call
end_ensemble_manager(ens_handle)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Frees up storage associated with an ensemble.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="duplicate_ens" id="duplicate_ens"></a><br>
<div class="routine"><em class="call">call duplicate_ens(ens1,
ens2, duplicate_time)</em>
<pre>
type(ensemble_type), intent(in)    :: <em class="code">ens1</em>
type(ensemble_type), intent(inout) :: <em class="code">ens2</em>
logical, intent(in)                :: <em class=
"code">duplicate_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Copies the contents of the vars array from ens1 into ens2. If
the num_copies and num_vars are not consistent or if the
distribution_type is not consistent, fails with an error. If
duplicate_time is true, the times from ens1 are copied over the
times of ens2. Only the vars array data is copied from the source
to the destination. Transpose the data after duplication if you
want to access the copies.</p>
<table border="0" cellpadding="3" width="100%" summary=''>
<tr>
<td valign="top"><em class="code">ens1</em></td>
<td>Ensemble handle of ensemble to be copies into ens2. Data from
the vars array will be replicated.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens2</em></td>
<td>Ensemble handle of ensemble into which ens1 vars data will be
copied.</td>
</tr>
<tr>
<td valign="top"><em class="code">duplicate_time  </em></td>
<td>If true, copy the times from ens1 into ens2, else leave ens2
times unchanged.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_my_num_copies" id="get_my_num_copies"></a><br>
<div class="routine"><em class="call">var =
get_my_num_copies(ens_handle)</em>
<pre>
integer                          :: <em class=
"code">get_my_num_copies</em>
type(ensemble_type), intent(in)  :: <em class=
"code">ens_handle</em>

</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns number of copies stored by this process when storing all
variables for a subset of copies. Same as num_copies if running
with only a single process.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns the number of copies stored by this process when
storing all variables for a subset of copies.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_my_num_vars" id="get_my_num_vars"></a><br>
<div class="routine"><em class="call">var =
get_my_num_vars(ens_handle)</em>
<pre>
integer                         :: <em class=
"code">get_my_num_vars</em>
type(ensemble_type), intent(in) :: <em class="code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns number of variables stored by this process when storing
all copies of a subset of variables. Same as num_vars if running
with only a single process.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns the number of vars stored by this process when storing
all copies of a subset of variables.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_my_copies" id="get_my_copies"></a><br>
<div class="routine"><em class="call">call
get_my_copies(ens_handle, copies)</em>
<pre>
type(ensemble_type), intent(in) :: <em class="code">ens_handle</em>
integer, intent(out)            :: <em class="code">copies(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a list of the global copy numbers stored on this process
when storing subset of copies of all variables.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">copies</em></td>
<td>List of all copies stored by this process when storing subset
of copies of all variables.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_my_vars" id="get_my_vars"></a><br>
<div class="routine"><em class="call">call get_my_vars(ens_handle,
vars)</em>
<pre>
type(ensemble_type), intent(in) :: <em class="code">ens_handle</em>
integer, intent(out)            :: <em class="code">vars(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a list of the global variable numbers stored on this
process when storing all copies of a subset of variables.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">vars</em></td>
<td>List of all variables stored on this process when storing all
copies of a subset of variables.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_copy_owner_index" id="get_copy_owner_index"></a><br>
<div class="routine"><em class="call">call
get_copy_owner_index(copy_number, owner, owners_index)</em>
<pre>
integer, intent(in)  :: <em class="code">copy_number</em>
integer, intent(out) :: <em class="code">owner</em>
integer, intent(out) :: <em class="code">owners_index</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given the global index of a copy number, returns the PE that
stores this copy when all variables of a subset of copies are
stored and the local storage index for this copy on that
process.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">copy_number</em></td>
<td>Global index of a copy from an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">owner</em></td>
<td>Process Element (PE) that stores this copy when each has all
variables of a subset of copies.</td>
</tr>
<tr>
<td valign="top"><em class="code">owners_index  </em></td>
<td>Local storage index for this copy on the owning process.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_var_owner_index" id="get_var_owner_index"></a><br>
<div class="routine"><em class="call">call
get_var_owner_index(var_number, owner, owners_index)</em>
<pre>
integer, intent(in)  :: <em class="code">var_number</em>
integer, intent(out) :: <em class="code">owner</em>
integer, intent(out) :: <em class="code">owners_index</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given the global index of a variable in the vector, returns the
PE that stores this variable when all copies of a subset of
variables are stored and the local storage index for this variable
on that process.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">var_number</em></td>
<td>Global index of a variable in the vector from an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">owner</em></td>
<td>Process Element (PE) that stores this variable when each has
all copies of subset of variables.</td>
</tr>
<tr>
<td valign="top"><em class="code">owners_index  </em></td>
<td>Local storage index for this variable on the owning
process.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="all_vars_to_all_copies" id=
"all_vars_to_all_copies"></a><br>
<div class="routine"><em class="call">call
all_vars_to_all_copies(ens_handle, label)</em>
<pre>
type(ensemble_type), intent(inout)        :: <em class=
"code">ens_handle</em>
character(len=*),    intent(in), optional :: <em class=
"code">label</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Transposes data from a representation in which each PE has a
subset of copies of all variables to one in which each has all
copies of a subset of variables. In the current implementation,
storage is not released so both representations are always
available. However, one representation may be current while the
other is out of date.</p>
<p>Different different numbers of copies, different lengths of the
vectors, different numbers of PEs and different implementations of
the MPI parallel libraries can have very different performance
characteristics. The namelist item <em class=
"code">communication_configuration</em> controls one of four
possible combinations of the operation order during the transposes.
If performance is an issue the various settings on this namelist
item can be explored. See the <a href="#Namelist">namelist
section</a> for more details.</p>
<p>The transpose routines make both representations of the data
equivalent until the next update to either the copies or the vars
arrays, so either can be used as a data source.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>The handle of the ensemble being transposed.</td>
</tr>
<tr>
<td valign="top"><em class="code">label</em></td>
<td>A character string label. If present, a timestamp with this
label is printed at the start and end of the transpose.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="all_copies_to_all_vars" id=
"all_copies_to_all_vars"></a><br>
<div class="routine"><em class="call">call
all_copies_to_all_vars(ens_handle, label)</em>
<pre>
type(ensemble_type), intent(inout)        :: <em class=
"code">ens_handle</em>
character(len=*),    intent(in), optional :: <em class=
"code">label</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Transposes data from a representation in which each processor
has all copies of a subset of variables to one in which each has a
subset of copies of all variables. In the current implementation,
storage is not released so both representations are always
available. However, one representation may be current while the
other is out of date.</p>
<p>Different different numbers of copies, different lengths of the
vectors, different numbers of PEs and different implementations of
the MPI parallel libraries can have very different performance
characteristics. The namelist item <em class=
"code">communication_configuration</em> controls one of four
possible combinations of the operation order during the transposes.
If performance is an issue the various settings on this namelist
item can be explored. See the <a href="#Namelist">namelist
section</a> for more details.</p>
<p>The transpose routines make both representations of the data
equivalent until the next update to either the copies or the vars
arrays, so either can be used as a data source.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>The handle of the ensemble being transposed.</td>
</tr>
<tr>
<td valign="top"><em class="code">label</em></td>
<td>A character string label. If present, a timestamp with this
label is printed at the start and end of the transpose.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="compute_copy_mean" id="compute_copy_mean"></a><br>
<div class="routine"><em class="call">call
compute_copy_mean(ens_handle, start_copy, end_copy, mean_copy)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer,             intent(in)    :: <em class=
"code">start_copy</em>
integer,             intent(in)    :: <em class=
"code">end_copy</em>
integer,             intent(in)    :: <em class=
"code">mean_copy</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes the mean of a contiguous subset of copies starting with
global index start_copy and ending with global index ens_copy. Mean
is written to global index mean_copy.</p>
<p>When this routine is called the ensemble must have all copies of
a subset of the vars. It updates the copies array with the mean, so
after this call the copies array data is more current and the vars
data is stale.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">start_copy</em></td>
<td>Global index of first copy in mean and sd computation.</td>
</tr>
<tr>
<td valign="top"><em class="code">end_copy</em></td>
<td>Global index of last copy in mean and sd computation.</td>
</tr>
<tr>
<td valign="top"><em class="code">mean_copy</em></td>
<td>Global index of copy into which mean is written.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="compute_copy_mean_sd" id="compute_copy_mean_sd"></a><br>
<div class="routine"><em class="call">call
compute_copy_mean_sd(ens_handle, start_copy, end_copy, mean_copy,
sd_copy)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer,             intent(in)    :: <em class=
"code">start_copy</em>
integer,             intent(in)    :: <em class=
"code">end_copy</em>
integer,             intent(in)    :: <em class=
"code">mean_copy</em>
integer,             intent(in)    :: <em class="code">sd_copy</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes the mean and standard deviation of a contiguous subset
of copies starting with global index start_copy and ending with
global index ens_copy. Mean is written to index mean_copy and
standard deviation to index sd_copy.</p>
<p>When this routine is called the ensemble must have all copies of
a subset of the vars. It updates the copies arrays with the mean
and sd, so after this call the copies array data is more current
and the vars data is stale.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">start_copy</em></td>
<td>Global index of first copy in mean and sd computation.</td>
</tr>
<tr>
<td valign="top"><em class="code">end_copy</em></td>
<td>Global index of last copy in mean and sd computation.</td>
</tr>
<tr>
<td valign="top"><em class="code">mean_copy</em></td>
<td>Global index of copy into which mean is written.</td>
</tr>
<tr>
<td valign="top"><em class="code">sd_copy</em></td>
<td>Global index of copy into which standard deviation is
written.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="compute_copy_mean_var" id=
"compute_copy_mean_var"></a><br>
<div class="routine"><em class="call">call
compute_copy_mean_var(ens_handle, start_copy, end_copy, mean_copy,
var_copy)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer,             intent(in)  :: <em class=
"code">start_copy</em>
integer,             intent(in)  :: <em class="code">end_copy</em>
integer,             intent(in)  :: <em class="code">mean_copy</em>
integer,             intent(in)  :: <em class="code">var_copy</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes the mean and variance of a contiguous subset of copies
starting with global index start_copy and ending with global index
ens_copy. Mean is written to index mean_copy and variance to index
var_copy.</p>
<p>When this routine is called the ensemble must have all copies of
a subset of the vars. It updates the copies arrays with the mean
and variance, so after this call the copies array data is more
current and the vars data is stale.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">start_copy</em></td>
<td>Global index of first copy in mean and sd computation.</td>
</tr>
<tr>
<td valign="top"><em class="code">end_copy</em></td>
<td>Global index of last copy in mean and sd computation.</td>
</tr>
<tr>
<td valign="top"><em class="code">mean_copy</em></td>
<td>Global index of copy into which mean is written.</td>
</tr>
<tr>
<td valign="top"><em class="code">var_copy</em></td>
<td>Global index of copy into which variance is written.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="prepare_to_update_vars" id=
"prepare_to_update_vars"></a><br>
<div class="routine"><em class="call">call
prepare_to_update_vars(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)  :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Call this routine before directly accessing the <em class=
"code">ens_handle%vars</em> array when the data is going to be
updated, and the incoming vars array should have the most current
data representation.</p>
<p>Internally the ensemble manager tracks which of the copies or
vars arrays, or both, have the most recently updated representation
of the data. For example, before a transpose (<em class=
"code">all_vars_to_all_copies()</em> or <em class=
"code">all_copies_to_all_vars()</em>) the code checks to be sure
the source array has the most recently updated representation
before it does the operation. After a transpose both
representations have the same update time and are both valid.</p>
<p>For efficiency reasons we allow the copies and vars arrays to be
accessed directly from other code without going through a routine
in the ensemble manager. The "prepare" routines verify that the
desired array has the most recently updated representation of the
data, and if needed marks which one has been updated so the
internal consistency checks have an accurate accounting of the
representations.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for the ensemble being accessed directly.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="prepare_to_update_copies" id=
"prepare_to_update_copies"></a><br>
<div class="routine"><em class="call">call
prepare_to_update_copies(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)  :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Call this routine before directly accessing the <em class=
"code">ens_handle%copies</em> array when the data is going to be
updated, and the incoming copies array should have the most current
data representation.</p>
<p>Internally the ensemble manager tracks which of the copies or
vars arrays, or both, have the most recently updated representation
of the data. For example, before a transpose (<em class=
"code">all_vars_to_all_copies()</em> or <em class=
"code">all_copies_to_all_vars()</em>) the code checks to be sure
the source array has the most recently updated representation
before it does the operation. After a transpose both
representations have the same update time and are both valid.</p>
<p>For efficiency reasons we allow the copies and vars arrays to be
accessed directly from other code without going through a routine
in the ensemble manager. The "prepare" routines verify that the
desired array has the most recently updated representation of the
data, and if needed marks which one has been updated so the
internal consistency checks have an accurate accounting of the
representations.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for the ensemble being accessed directly.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="prepare_to_read_from_vars" id=
"prepare_to_read_from_vars"></a><br>
<div class="routine"><em class="call">call
prepare_to_read_from_vars(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)  :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Call this routine before directly accessing the <em class=
"code">ens_handle%vars</em> array for reading only, when the
incoming vars array should have the most current data
representation.</p>
<p>Internally the ensemble manager tracks which of the copies or
vars arrays, or both, have the most recently updated representation
of the data. For example, before a transpose (<em class=
"code">all_vars_to_all_copies()</em> or <em class=
"code">all_copies_to_all_vars()</em>) the code checks to be sure
the source array has the most recently updated representation
before it does the operation. After a transpose both
representations have the same update time and are both valid.</p>
<p>For efficiency reasons we allow the copies and vars arrays to be
accessed directly from other code without going through a routine
in the ensemble manager. The "prepare" routines verify that the
desired array has the most recently updated representation of the
data, and if needed marks which one has been updated so the
internal consistency checks have an accurate accounting of the
representations.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for the ensemble being accessed directly.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="prepare_to_read_from_copies" id=
"prepare_to_read_from_copies"></a><br>
<div class="routine"><em class="call">call
prepare_to_read_from_copies(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)  :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Call this routine before directly accessing the <em class=
"code">ens_handle%copies</em> array for reading only, when the
incoming copies array should have the most current data
representation.</p>
<p>Internally the ensemble manager tracks which of the copies or
vars arrays, or both, have the most recently updated representation
of the data. For example, before a transpose (<em class=
"code">all_vars_to_all_copies()</em> or <em class=
"code">all_copies_to_all_vars()</em>) the code checks to be sure
the source array has the most recently updated representation
before it does the operation. After a transpose both
representations have the same update time and are both valid.</p>
<p>For efficiency reasons we allow the copies and vars arrays to be
accessed directly from other code without going through a routine
in the ensemble manager. The "prepare" routines verify that the
desired array has the most recently updated representation of the
data, and if needed marks which one has been updated so the
internal consistency checks have an accurate accounting of the
representations.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for the ensemble being accessed directly.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="prepare_to_write_to_vars" id=
"prepare_to_write_to_vars"></a><br>
<div class="routine"><em class="call">call
prepare_to_write_to_vars(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)  :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Call this routine before directly accessing the <em class=
"code">ens_handle%vars</em> array for writing. This routine differs
from the 'update' version in that it doesn't care what the original
data state is. This routine might be used in the case where an
array is being filled for the first time and consistency with the
data in the copies array is not an issue.</p>
<p>Internally the ensemble manager tracks which of the copies or
vars arrays, or both, have the most recently updated representation
of the data. For example, before a transpose (<em class=
"code">all_vars_to_all_copies()</em> or <em class=
"code">all_copies_to_all_vars()</em>) the code checks to be sure
the source array has the most recently updated representation
before it does the operation. After a transpose both
representations have the same update time and are both valid.</p>
<p>For efficiency reasons we allow the copies and vars arrays to be
accessed directly from other code without going through a routine
in the ensemble manager. The "prepare" routines verify that the
desired array has the most recently updated representation of the
data, and if needed marks which one has been updated so the
internal consistency checks have an accurate accounting of the
representations.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for the ensemble being accessed directly.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="prepare_to_write_to_copies" id=
"prepare_to_write_to_copies"></a><br>
<div class="routine"><em class="call">call
prepare_to_write_to_copies(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)  :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Call this routine before directly accessing the <em class=
"code">ens_handle%copies</em> array for writing. This routine
differs from the 'update' version in that it doesn't care what the
original data state is. This routine might be used in the case
where an array is being filled for the first time and consistency
with the data in the vars array is not an issue.</p>
<p>Internally the ensemble manager tracks which of the copies or
vars arrays, or both, have the most recently updated representation
of the data. For example, before a transpose (<em class=
"code">all_vars_to_all_copies()</em> or <em class=
"code">all_copies_to_all_vars()</em>) the code checks to be sure
the source array has the most recently updated representation
before it does the operation. After a transpose both
representations have the same update time and are both valid.</p>
<p>For efficiency reasons we allow the copies and vars arrays to be
accessed directly from other code without going through a routine
in the ensemble manager. The "prepare" routines verify that the
desired array has the most recently updated representation of the
data, and if needed marks which one has been updated so the
internal consistency checks have an accurate accounting of the
representations.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for the ensemble being accessed directly.</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PRIVATE INTERFACES</h2>
<table summary='private routine list'>
<tr>
<td> </td>
<td><a href="#assign_tasks_to_pes">assign_tasks_to_pes</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#calc_tasks_on_each_node">calc_tasks_on_each_node</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#create_pe_to_task_list">create_pe_to_task_list</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_copy_list">get_copy_list</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_max_num_copies">get_max_num_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_max_num_vars">get_max_num_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_var_list">get_var_list</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#round_robin">round_robin</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#set_up_ens_distribution">set_up_ens_distribution</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#simple_layout">simple_layout</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#sort_task_list">sort_task_list</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#timestamp_message">timestamp_message</a></td>
</tr>
</table>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_max_num_copies" id="get_max_num_copies"></a><br>
<div class="routine"><em class="call">var =
get_max_num_copies(num_copies)</em>
<pre>
integer              :: <em class="code">get_max_num_copies</em>
integer, intent(in)  :: <em class="code">num_copies</em>

</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the largest number of copies that are on any pe when var
complete. Depends on distribution_type with only option 1 currently
implemented. Used to get size for creating storage to receive a
list of the copies on a PE.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns the largest number of copies any an individual PE when
var complete.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies</em></td>
<td>Total number of copies in the ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_max_num_vars" id="get_max_num_vars"></a><br>
<div class="routine"><em class="call">var =
get_max_num_vars(num_vars)</em>
<pre>
integer              :: <em class="code">get_max_num_vars</em>
integer, intent(in)  :: <em class="code">num_vars</em>

</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the largest number of vars that are on any pe when copy
complete. Depends on distribution_type with only option 1 currently
implemented. Used to get size for creating storage to receive a
list of the vars on a PE.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns the largest number of vars any an individual PE when
copy complete.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies</em></td>
<td>Total number of vars in an ensemble vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_up_ens_distribution" id=
"set_up_ens_distribution"></a><br>
<div class="routine"><em class="call">call
set_up_ens_distribution(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Figures out how to lay out the copy complete and vars complete
distributions. The distribution_type identifies different options.
Only distribution_type 1 is implemented. This puts every Nth var or
copy on a given processor where N is the total number of
processes.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_var_list" id="get_var_list"></a><br>
<div class="routine"><em class="call">call get_var_list(num_vars,
pe, var_list, pes_num_vars)</em>
<pre>
integer,   intent(in)     :: <em class="code">num_vars</em>
integer,   intent(in)     :: <em class="code">pe</em>
integer,   intent(out)    :: <em class="code">var_list(:)</em>
integer,   intent(out)    :: <em class="code">pes_num_vars</em>
</pre></div>

<p>Returns a list of the vars stored by process pe when copy
complete and the number of these vars. var_list must be dimensioned
large enough to hold all vars. Depends on distribution_type with
only option 1 currently implemented.</p>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_copy_list" id="get_copy_list"></a><br>
<div class="routine"><em class="call">call
get_copy_list(num_copies, pe, copy_list, pes_num_copies)</em>
<pre>
integer,   intent(in)     :: <em class="code">num_copies</em>
integer,   intent(in)     :: <em class="code">pe</em>
integer,   intent(out)    :: <em class="code">copy_list(:)</em>
integer,   intent(out)    :: <em class="code">pes_num_copies</em>
</pre></div>

<p>Returns a list of the copies stored by process pe when var
complete and the number of these copies. copy_list must be
dimensioned large enough to hold all copies. Depends on
distribution_type with only option 1 currently implemented.</p>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="timestamp_message" id="timestamp_message"></a><br>
<div class="routine"><em class="call">call timestamp_message(msg
<em class="optionalcode">[, sync]</em> <em class=
"optionalcode">[, alltasks]</em>)</em>
<pre>
character(len=*), intent(in)           :: <em class="code">msg</em>
logical,          intent(in), optional :: <em class=
"optionalcode">sync</em>
logical,          intent(in), optional :: <em class=
"optionalcode">alltasks</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Write current time and message to stdout and log file. If sync
is present and true, sync mpi jobs before printing time. If
alltasks is present and true, all tasks print the time. The default
is only task 0 prints a timestamp.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">msg</em></td>
<td>character string to prepend to the time info</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">sync</em></td>
<td>if present and true, execute an MPI_Barrier() to sync all MPI
tasks before printing the time. this means the time will be the
value of the slowest of the tasks to reach this point.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">alltasks</em></td>
<td>if present and true, have all tasks print out a timestamp. the
default is for just task 0 to print. the usual combination is
either sync=true and alltasks=false, or sync=false and
alltasks=true.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="print_ens_handle" id="print_ens_handle"></a><br>
<div class="routine"><em class="call">call
print_ens_handle(ens_handle, force, label)</em>
<pre>
type(ensemble_type),        intent(in) :: <em class=
"code">ens_handle</em>
logical,          optional, intent(in) :: <em class=
"code">force</em>
character(len=*), optional, intent(in) :: <em class=
"code">label</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>For debugging use, dump the contents of an ensemble handle
derived type. If the <em class="code">debug</em> namelist item is
true, this will print in any case. If <em class="code">debug</em>
is false, set <em class="code">force</em> to true to force
printing. The optional string label can help provide context for
the output.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>The derived type to print information about.</td>
</tr>
<tr>
<td valign="top"><em class="code">force</em></td>
<td>If the <em class="code">debug</em> namelist item is false, set
this to true to enable printing.</td>
</tr>
<tr>
<td valign="top"><em class="code">label</em></td>
<td>Optional string label to print to provide context for the
output.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="assign_tasks_to_pes" id="assign_tasks_to_pes"></a><br>
<div class="routine"><em class="call">call
assign_tasks_to_pes(ens_handle, nEns_members, layout_type)</em>
<pre>
type(ensemble_type), intent(inout)    :: <em class=
"code">ens_handle</em>
integer,             intent(in)       :: <em class=
"code">nEns_members</em>
integer,             intent(inout)    :: <em class=
"code">layout_type</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Calulate the task layout based on the tasks per node and the
total number of tasks. Allows the user to spread out the ensemble
members as much as possible to balance memory usage between nodes.
Possible options: 1. Standard task layout - first n tasks have the
ensemble members my_pe = my_task_id() 2. Round-robin on the
nodes</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"></td>
<td></td>
</tr>
<tr>
<td valign="top"></td>
<td></td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="round_robin" id="round_robin"></a><br>
<div class="routine"><em class="call">call
round_robin(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)    :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Round-robin MPI task layout starting at the first node. Starting
on the first node forces pe 0 = task 0. The smoother code assumes
task 0 has an ensemble member. If you want to break the assumption
that pe 0 = task 0, this routine is a good place to start. Test
with the smoother.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="create_pe_to_task_list" id=
"create_pe_to_task_list"></a><br>
<div class="routine"><em class="call">call
create_pe_to_task_list(ens_handle)</em>
<pre>
type(ensemble_type), intent(inout)    :: <em class=
"code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Creates the <em class="code">ens_handle%pe_to_task_list</em>.
<em class="code">ens_handle%task_to_pe_list</em> must have been
assigned first, otherwise this routine will just return
nonsense.</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="calc_tasks_on_each_node" id=
"calc_tasks_on_each_node"></a><br>
<div class="routine"><em class="call">call
calc_tasks_on_each_node(nodes, last_node_task_number)</em>
<pre>
integer, intent(out)  :: <em class=
"code">last_node_task_number</em>
integer, intent(out)  :: <em class="code">nodes</em>
</pre></div>

<p>Finds the of number nodes and how many tasks are on the last
node, given the number of tasks and the tasks_per_node (ptile). The
total number of tasks is num_pes = task_count() The last node may
have fewer tasks, for example, if ptile = 16 and the number of mpi
tasks = 17</p>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="simple_layout" id="simple_layout"></a><br>
<div class="routine"><em class="call">call
simple_layout(ens_handle, n)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer,             intent(in)    :: <em class="code">n</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>assigns the arrays task_to_pe_list and pe_to_task list for the
simple layout where my_pe = my_task_id()</p>

<dl>
   <dt>ens_handle</dt>
   <dd>Handle for an ensemble.</dd>
   <dt>n</dt>
   <dd>size</dd>
</dl>

<a name="sort_task_list" id="sort_task_list"></a><br>
<div class="routine"><em class="call">call sort_task_list(i, idx,
n)</em>
<pre>
integer, intent(in)    :: <em class="code">n</em>
integer, intent(inout) :: <em class=
"code">x(n)</em>   ! array to be sorted
integer, intent(out)   :: <em class=
"code">idx(n)</em> ! index of sorted array
</pre></div>

<p>sorts an array and returns the sorted array, and the index of
the original array</p>

<dl>
   <dt>n</dt>
   <dd>size</dd>
   <dt>x(n)</dt>
   <dd>array to be sorted</dd>
   <dt>idx(n)</dt>
   <dd>index of sorted array</dd>
</dl>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="map_pe_to_task" id="map_pe_to_task"></a><br>
<div class="routine"><em class="call">call
map_pe_to_task(ens_handle, p)</em>
<pre>
type(ensemble_type), intent(in) :: <em class="code">ens_handle</em>
integer,             intent(in) :: <em class="code">p</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Return the physical task for my_pe</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">p</em></td>
<td>The MPI task corresponding to the given PE number</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="map_task_to_pe" id="map_task_to_pe"></a><br>
<div class="routine"><em class="call">call
map_task_to_pe(ens_handle, t)</em>
<pre>
type(ensemble_type), intent(in) :: <em class="code">ens_handle</em>
integer,             intent(in) :: <em class="code">t</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Return my_pe corresponding to the physical task</p>
<table width="100%" border="0" summary="argument description"
cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Handle for an ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">t</em></td>
<td>Return the PE corresponding to the given MPI task number.</td>
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
<li>input.nml</li>
<li>State vector restart files, either one for all copies or one
per copy.</li>
<li>State vector output files, either one for all copies or one per
copy.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ol>
<li>none</li>
</ol>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%"
summary='error list'>
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">init_ensemble_manager</td>
<!-- message -->
<td valign="top">only distribution type 1 is implemented</td>
<!-- comment -->
<td valign="top">For now, can't request option other than 1 for
layout</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_ensemble_restart</td>
<!-- message -->
<td valign="top">start_from_restart in filter_nml and
single_restart_file_in in ensemble_manager_nml cannot both be
false</td>
<!-- comment -->
<td valign="top">Doesn't make sense to specify both of these
options.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_copy</td>
<!-- message -->
<td valign="top">Requested copy is &gt; maximum_copy</td>
<!-- comment -->
<td valign="top">Can't ask for a copy that is greater than the
maximum.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_copy</td>
<!-- message -->
<td valign="top">Size of vars ### Must be at least ###</td>
<!-- comment -->
<td valign="top">The vars array is not big enough to hold the
returned copy of the vector.</td>
</tr>
<tr><!-- routine -->
<td valign="top">put_copy</td>
<!-- message -->
<td valign="top">Requested copy: ### is &gt; maximum_copy: ###</td>
<!-- comment -->
<td valign="top">Can't ask for a copy that is greater than
maximum.</td>
</tr>
<tr><!-- routine -->
<td valign="top">put_copy</td>
<!-- message -->
<td valign="top">Size of vars: ### Must be at least ###</td>
<!-- comment -->
<td valign="top">The vars array is not big enough to hold the state
vector.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_ensemble_time</td>
<!-- message -->
<td valign="top">indx ### cannot exceed ###</td>
<!-- comment -->
<td valign="top">The index of the requested copy must be no greater
than the maximum number of copies.</td>
</tr>
<tr><!-- routine -->
<td valign="top">duplicate_ens</td>
<!-- message -->
<td valign="top">num_copies ### and ### must be equal</td>
<!-- comment -->
<td valign="top">Number of copies in ensembles being copied must be
the same.</td>
</tr>
<tr><!-- routine -->
<td valign="top">duplicate_ens</td>
<!-- message -->
<td valign="top">num_vars ### and ### must be equal</td>
<!-- comment -->
<td valign="top">Number of variables in ensembles being copied must
be the same.</td>
</tr>
<tr><!-- routine -->
<td valign="top">duplicate_ens</td>
<!-- message -->
<td valign="top">distribution_type ### and ### must be equal.</td>
<!-- comment -->
<td valign="top">Distribution types of ensembles being copies must
be the same.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_my_copies</td>
<!-- message -->
<td valign="top">Array copies only has size ### but must be at
least ###</td>
<!-- comment -->
<td valign="top">The copies array must be large enough to hold all
copies of the state vector.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_my_vars</td>
<!-- message -->
<td valign="top">Array vars only has size ### but must be at least
###</td>
<!-- comment -->
<td valign="top">The vars array must be large enough to hold all
variables of the state vector.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>Additional options for the layout of ensemble storage may lead
to improved performance for different problem sizes on different
architectures.</p>
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
