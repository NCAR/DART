<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program perfect_model_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>program <em class="program">perfect_model_obs</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>Main program for creating synthetic observation sequences given
a model for use in filter assimilations. Reads in an observation
sequence file which has only observation definitions and generates
synthetic observation values for an output observation sequence
file. The execution of perfect_model_obs is controlled by the input
observation sequence file and the model time-stepping capabilities
in a manner analogous to that used by the filter program.</p>
<!--================================================================-->
<!--============== DESCRIPTION OF A NAMELIST ========================-->
<!--================================================================-->
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
&amp;perfect_model_obs_nml
   single_file_in             = .false.,
   read_input_state_from_file = .false.,
   input_state_files          = "",
   init_time_days             = 0,
   init_time_seconds          = 0,

   single_file_out            = .false.,
   output_state_files         = "",
   write_output_state_to_file = .false.,
   output_interval            = 1,

   distributed_state          = .false.,
   async                      = 0,
   adv_ens_command            = "./advance_model.csh",
   tasks_per_model_advance    = 1,

   obs_seq_in_file_name       = "obs_seq.in",
   obs_seq_out_file_name      = "obs_seq.out",
   first_obs_days             = -1,
   first_obs_seconds          = -1,
   last_obs_days              = -1,
   last_obs_seconds           = -1,
   obs_window_days            = -1,
   obs_window_seconds         = -1,

   trace_execution            = .false.,
   output_timestamps          = .false.,
   print_every_nth_obs        = 0,
   output_forward_op_errors   = .false.,
   silence                    = .false.,
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
<td>read_input_state_from_file</td>
<td>logical</td>
<td>If false, model_mod must provide the input state.</td>
</tr>
<tr>
<td>single_file_in</td>
<td>logical</td>
<td>Get all states from a single file.</td>
</tr>
<tr>
<td>input_state_files</td>
<td>character(len=256) dimension(MAX_NUM_DOMS)</td>
<td>A list of files, one per domain. Each file must be a text file
containing the name of the NetCDF file to open.</td>
</tr>
<tr>
<td>write_output_state_to_file</td>
<td>logical</td>
<td>If false, state is not written out.</td>
</tr>
<tr>
<td>single_file_out</td>
<td>logical</td>
<td>Write all states to a single file.</td>
</tr>
<tr>
<td>output_state_files</td>
<td>character(len=256) dimension(MAX_NUM_DOMS)</td>
<td>A list of files, one per domain. Each file must be a text file
containing the names of the NetCDF file to open.</td>
</tr>
<tr>
<td>init_time_days</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, override the initial
data time read from restart file.</td>
</tr>
<tr>
<td>init_time_seconds</td>
<td>integer</td>
<td>If negative don't use. If non-negative, override the initial
data time read from restart file.</td>
</tr>
<tr>
<td>output_interval</td>
<td>integer</td>
<td>Output state and observation diagnostics every nth assimilation
time, n is output_interval.</td>
</tr>
<tr>
<td>distributed_state</td>
<td>logical</td>
<td>True means the ensemble data is distributed across all tasks as
it is read in, so a single task never has to have enough memory to
store the data for an ensemble member. Large models should always
set this to .true., while for small models it may be faster to set
this to .false.</td>
</tr>
<tr>
<td>async</td>
<td>integer</td>
<td>Controls method for advancing model:
<ul style="list-style: none;">
<li>0 = subroutine call</li>
<li>2 = shell command, single task model</li>
<li>4 = shell command, parallel model</li>
</ul>
</td>
</tr>
<tr>
<td>adv_ens_command</td>
<td>character(len=129)</td>
<td>Command sent to shell if async == 2 or 4.</td>
</tr>
<tr>
<td>tasks_per_model_advance</td>
<td>integer</td>
<td>Number of tasks to use while advancing the model.</td>
</tr>
<tr>
<td>obs_seq_in_file_name</td>
<td>character(len=256)</td>
<td>File name from which to read an observation sequence.</td>
</tr>
<tr>
<td>obs_seq_out_file_name</td>
<td>character(len=256)</td>
<td>File name to which to write output observation sequence.</td>
</tr>
<tr>
<td>first_obs_days</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore any
observations before this time.</td>
</tr>
<tr>
<td>first_obs_seconds</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore any
observations before this time.</td>
</tr>
<tr>
<td>last_obs_days</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore any
observations after this time.</td>
</tr>
<tr>
<td>last_obs_seconds</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore any
observations after this time.</td>
</tr>
<tr>
<td>obs_window_days</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, reserved for future
use.</td>
</tr>
<tr>
<td>obs_window_seconds</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, reserved for future
use.</td>
</tr>
<tr>
<td>trace_execution</td>
<td>logical</td>
<td>True means output very detailed messages about what routines
are being called in the main loop. Useful if a job hangs or
otherwise doesn't execute as expected.</td>
</tr>
<tr>
<td>output_timestamps</td>
<td>logical</td>
<td>True means output timestamps before and after the model advance
and the forward observation computation phases.</td>
</tr>
<tr>
<td>print_every_nth_obs</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, print a message noting
the processing of every Nth observation.</td>
</tr>
<tr>
<td>output_forward_op_errors</td>
<td>logical</td>
<td>True means output errors from forward observation operators.
This is the 'istatus' error return code from the model interpolate
routine. An ascii text file 'forward_op_errors' will be created in
the current directory. Each line will contain an observation key
number, and the istatus return code.</td>
</tr>
<tr>
<td>silence</td>
<td>logical</td>
<td>True means output almost no runtime messages. Not recommended
for general use, but can speed test programs if the execution time
becomes dominated by the volume of output.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
time_manager_mod
obs_sequence_mod
obs_def_mod
obs_model_mod
assim_model_mod
mpi_utilities_mod
random_seq_mod
ensemble_manager_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>observation sequence input file; name comes from
obs_seq_in_file_name</li>
<li>observation sequence output file; name comes from
obs_seq_out_file_name</li>
<li>input state vector file; name comes from
restart_in_file_name</li>
<li>output state vector file; name comes from
restart_out_file_name</li>
<li>perfect_model_mod.nml in input.nml</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>none</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">perfect_main</td>
<!-- message -->
<td valign="top">Only use one mpi process here: ### were
requested</td>
<!-- comment -->
<td valign="top">Don't use mpi for this.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Descibe Future Plans.                                            -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none</p>
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
