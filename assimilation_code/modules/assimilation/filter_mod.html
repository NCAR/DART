<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module filter_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE filter_mod</h1>
<table border="0" summary="dart header" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / 
<!-- <A HREF="#References">REFERENCES</A> / -->
<a href="#Errors">ERRORS</a> / <a href="#FuturePlans">PLANS</a> /
<a href="#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Main module for driving ensemble filter assimilations. Used by
filter.f90, perfect_model_obs.f90, model_mod_check.f90, and a
variety of test programs. See the <a href=
"../../programs/filter/filter.html">filter description</a> for a
general description of filter capabilities and controls.</p>
<p><em class="program">filter_mod</em> is a Fortran 90 module, and
provides a large number of options for controlling execution
behavior and parameter configuration that are driven from its
namelist. See the <a href="#Namelist">namelist</a> section below
for more details. The number of assimilation steps to be done is
controlled by the input observation sequence and by the
time-stepping capabilities of the model being used in the
assimilation.</p>
<p>See the <a href="http://www.image.ucar.edu/DAReS/DART">DART web
site</a> for more documentation, including a discussion of the
capabilities of the assimilation system, a diagram of the entire
execution cycle, the options and features.</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
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
&amp;filter_nml
   single_file_in               = .false.,
   input_state_files            = '',
   input_state_file_list        = '',
   init_time_days               = 0,
   init_time_seconds            = 0,
   perturb_from_single_instance = .false.,
   perturbation_amplitude       = 0.2,

   stages_to_write              = 'output'

   single_file_out              = .false.,
   output_state_files           = '',
   output_state_file_list       = '',
   output_interval              = 1,
   output_members               = .true.,
   num_output_state_members     = 0,
   output_mean                  = .true.,
   output_sd                    = .true.,
   write_all_stages_at_end      = .false.,
   compute_posterior            = .true.

   ens_size                     = 20,
   num_groups                   = 1,
   distributed_state            = .true.,

   async                        = 0,
   adv_ens_command              = "./advance_model.csh",
   tasks_per_model_advance      = 1,

   obs_sequence_in_name         = "obs_seq.out",
   obs_sequence_out_name        = "obs_seq.final",
   num_output_obs_members       = 0,
   first_obs_days               = -1,
   first_obs_seconds            = -1,
   last_obs_days                = -1,
   last_obs_seconds             = -1,
   obs_window_days              = -1,
   obs_window_seconds           = -1,

   inf_flavor                   = 0,                       0,
   inf_initial_from_restart     = .false.,                 .false.,
   inf_sd_initial_from_restart  = .false.,                 .false.,
   inf_deterministic            = .true.,                  .true.,
   inf_initial                  = 1.0,                     1.0,
   inf_lower_bound              = 1.0,                     1.0,
   inf_upper_bound              = 1000000.0,               1000000.0,
   inf_damping                  = 1.0,                     1.0,
   inf_sd_initial               = 0.0,                     0.0,
   inf_sd_lower_bound           = 0.0,                     0.0,
   inf_sd_max_change            = 1.05,                    1.05,

   trace_execution              = .false.,
   output_timestamps            = .false.,
   output_forward_op_errors     = .false.,
   write_obs_every_cycle        = .false.,
   silence                      = .false.,
 /
</pre></div>
<br>
<br>
<p>Particular options to be aware of are: async, ens_size, cutoff
(localization radius), inflation flavor, outlier_threshold, restart
filenames (including inflation), obs_sequence_in_name,
horiz_dist_only, binary or ascii controls for observation sequence
file formats. Some of these important items are located in other
namelists, but all are in the same input.nml file.</p>
<p>The inflation control variables are all dimensioned 2, the first
value controls the prior inflation and the second controls the
posterior inflation.</p>
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
<td>single_file_in</td>
<td>logical</td>
<td>True means all ensemble members are read from a single NetCDF
file. False means each member is in a separate file. NOT SUPPORTED
as of March, 2017 only multiple files can be used.</td>
</tr>
<tr>
<td>input_state_files</td>
<td>character(len=256) dimension(MAXFILES)</td>
<td>A list of the NetCDF files to open to read the state vectors.
Models using multiple domains must put the domain and ensemble
numbers in the file names. The order and format of those is to be
determined. NOT SUPPORTED as of March, 2017.</td>
</tr>
<tr>
<td>input_state_file_list</td>
<td>character(len=256) dimension(MAXFILES)</td>
<td>A list of files, one per domain. Each file must be a text file
containing the names of the NetCDF files to open, one per ensemble
member, one per line.</td>
</tr>
<tr>
<td>init_time_days</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, override the initial
days read from state data restart files.</td>
</tr>
<tr>
<td>init_time_seconds</td>
<td>integer</td>
<td>If negative don't use. If non-negative, override the initial
seconds read from state data restart files.</td>
</tr>
<tr>
<td>perturb_from_single_instance</td>
<td>logical</td>
<td>False means start from an ensemble-sized set of restart files.
True means perturb a single state vector from one restart file.
This may be done by model_mod, if model_mod provides subroutine
<em class="code">pert_model_copies</em>.</td>
</tr>
<tr>
<td>perturbation_amplitude</td>
<td>real(r8)</td>
<td>Standard deviation for the gaussian noise added when generating
perturbed ensemble members. Ignored if <em class=
"code">perturb_from_single_instance = .false.</em> or the
perturbed ensemble is created in model_mod.<br>
Random noise values drawn from a gaussian distribution with this
standard deviation will be added to the data in a single initial
ensemble member to generate the rest of the members.<br>
This option is more frequently used in the low order models and
less frequently used in large models. This is in part due to the
different scales of real geophysical variable values, and the
resulting inconsistencies between related field values. A more
successful initial condition generation strategy is to generate
climatological distributions from long model runs which have
internally consistent structures and values and then use
observations with a 'spin-up' period of assimilation to shape the
initial states into a set of members with enough spread and which
match the current set of observations.</td>
</tr>
<tr>
<td>stages_to_write</td>
<td>character(len=10) dimension(4)</td>
<td>Controls diagnostic and restart output. Valid values are
'input', 'preassim', 'postassim', 'output', and 'null'.</td>
</tr>
<tr>
<td>single_file_out</td>
<td>logical</td>
<td>True means all ensemble members are written to a single NetCDF
file. False means each member is output in a separate file. NOT
SUPPORTED as of March, 2017 - only multiple files can be used.</td>
</tr>
<tr>
<td>output_state_files</td>
<td>character(len=256) dimension(MAXFILES)</td>
<td>A list of the NetCDF files to open for writing updated state
vectors. Models using multiple domains must put the domain and
ensemble numbers in the file names. The order and format of those
is to be determined. NOT SUPPORTED as of March, 2017.</td>
</tr>
<tr>
<td>output_state_file_list</td>
<td>character(len=256) dimension(MAXFILES)</td>
<td>A list of files, one per domain. Each file must be a text file
containing the names of the NetCDF files to open, one per ensemble
member, one per line.</td>
</tr>
<tr>
<td>output_interval</td>
<td>integer</td>
<td>Output state and observation diagnostics every 'N'th
assimilation time, N is output_interval.</td>
</tr>
<tr>
<td>output_members</td>
<td>logical</td>
<td>True means output the ensemble members in any stage that is
enabled.</td>
</tr>
<tr>
<td>num_output_state_members</td>
<td>integer</td>
<td>Number of ensemble members to be included in the state
diagnostic output for stages 'preassim' and 'postassim'.
output_members must be TRUE.</td>
</tr>
<tr>
<td>output_mean</td>
<td>logical</td>
<td>True means output the ensemble mean in any stage that is
enabled.</td>
</tr>
<tr>
<td>output_sd</td>
<td>logical</td>
<td>True means output the ensemble standard deviation (spread) in
any stage that is enabled.</td>
</tr>
<tr>
<td>write_all_stages_at_end</td>
<td>logical</td>
<td>For most cases this should be .false. and data will be output
as it is generated for the 'preassim', 'postassim' diagnostics, and
then restart data will be output at the end. However, if I/O time
dominates the runtime, setting this to .true. will store the data
and it can all be written in parallel at the end of the execution.
This will require slightly more memory at runtime, but can lower
the cost of the job significantly in some cases.</td>
</tr>
<tr>
<td>compute_posterior</td>
<td>logical</td>
<td>If .false., skip computing posterior forward operators and do
not write posterior values in the obs_seq.final file. Saves time
and memory. Cannot enable posterior inflation and skip computing
the posteriors. For backwards compatibility the default for this is
.true.</td>
</tr>
<tr>
<td>ens_size</td>
<td>integer</td>
<td>Size of ensemble.</td>
</tr>
<tr>
<td>num_groups</td>
<td>integer</td>
<td>Number of groups for hierarchical filter. It should evenly
divide ens_size.</td>
</tr>
<tr>
<td>distributed_state</td>
<td>logical</td>
<td>True means the ensemble data is distributed across all tasks as
it is read in, so a single task never has to have enough memory to
store the data for an ensemble member. Large models should always
set this to .true., while for small models it may be faster to set
this to .false. This is different from <em>&amp;assim_tools_mod ::
distributed_mean</em> .</td>
</tr>
<tr>
<td>async</td>
<td>integer</td>
<td>Controls method for advancing model:
<ul>
<li>0 is subroutine call</li>
<li>2 is shell command</li>
<li>4 is mpi-job script</li>
</ul>
Ignored if filter is not controlling the model advance, e.g. in
CESM assimilations.</td>
</tr>
<tr>
<td>adv_ens_command</td>
<td>character(len=256)</td>
<td>Command sent to shell if async is 2.</td>
</tr>
<tr>
<td>tasks_per_model_advance</td>
<td>integer</td>
<td>Number of tasks to assign to each ensemble member advance.</td>
</tr>
<tr>
<td>obs_sequence_in_name</td>
<td>character(len=256)</td>
<td>File name from which to read an observation sequence.</td>
</tr>
<tr>
<td>obs_sequence_out_name</td>
<td>character(len=256)</td>
<td>File name to which to write output observation sequence.</td>
</tr>
<tr>
<td>num_output_obs_members</td>
<td>integer</td>
<td>Number of ensemble members to be included in the output
observation sequence file.</td>
</tr>
<tr>
<td>first_obs_days</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore all
observations before this time.</td>
</tr>
<tr>
<td>first_obs_seconds</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore all
observations before this time.</td>
</tr>
<tr>
<td>last_obs_days</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore all
observations after this time.</td>
</tr>
<tr>
<td>last_obs_seconds</td>
<td>integer</td>
<td>If negative, don't use. If non-negative, ignore all
observations after this time.</td>
</tr>
<tr>
<td>obs_window_days</td>
<td>integer</td>
<td>Assimilation window days; defaults to model timestep size.</td>
</tr>
<tr>
<td>obs_window_seconds</td>
<td>integer</td>
<td>Assimilation window seconds; defaults to model timestep
size.</td>
</tr>
<tr>
<td colspan="3">All variables named <code>inf_*</code> are arrays of length
2.<br>
The first element controls the prior inflation, the second element
controls the posterior inflation. See <a href=
"../../programs/filter/filter.html#Inflation">filter.html</a> for a
discussion of inflation and effective strategies.</td>
</tr>
<tr>
<td>inf_flavor</td>
<td>integer array dimension(2)</td>
<td>Inflation flavor for [prior, posterior]
<ul>
<li>0 = none</li>
<li>2 = spatially-varying state-space (gaussian)</li>
<li>3 = spatially-fixed state-space (gaussian)</li>
<li>4 = Relaxation To Prior Spread (Posterior inflation only)</li>
<li>5 = enhanced spatially-varying state-space (inverse gamma)</li>
</ul>
(See inf_sd_initial below for how to set the time evolution
options.)</td>
</tr>
<tr>
<td>inf_initial_from_restart</td>
<td>logical array dimension(2)</td>
<td>If true, get initial mean values for inflation from restart
file. If false, use the corresponding namelist value <em class=
"code">inf_initial</em>.</td>
</tr>
<tr>
<td>inf_sd_initial_from_restart</td>
<td>logical array dimension(2)</td>
<td>If true, get initial standard deviation values for inflation
from restart file. If false, use the corresponding namelist value
<em class="code">inf_sd_initial</em>.</td>
</tr>
<tr>
<td>inf_deterministic</td>
<td>logical array dimension(2)</td>
<td>True means deterministic inflation, false means
stochastic.</td>
</tr>
<tr>
<td>inf_initial</td>
<td>real(r8) dimension(2)</td>
<td>Initial value of inflation if not read from restart file.</td>
</tr>
<tr>
<td>inf_lower_bound</td>
<td>real(r8) dimension(2)</td>
<td>Lower bound for inflation value.</td>
</tr>
<tr>
<td>inf_upper_bound</td>
<td>real(r8) dimension(2)</td>
<td>Upper bound for inflation value.</td>
</tr>
<tr>
<td>inf_damping</td>
<td>real(r8) dimension(2)</td>
<td>Damping factor for inflation mean values. The difference
between the current inflation value and 1.0 is multiplied by this
factor before the next assimilation cycle. The value should be
between 0.0 and 1.0. Setting a value of 0.0 is full damping, which
in fact turns all inflation off by fixing the inflation value at
1.0. A value of 1.0 turns inflation damping off leaving the
original inflation value unchanged.</td>
</tr>
<tr>
<td>inf_sd_initial</td>
<td>real(r8) dimension(2)</td>
<td>Initial value of inflation standard deviation if not read from
restart file. If ≤ 0, do not update the inflation values, so they
are time-constant. If positive, the inflation values will adapt
through time, so they are time-varying.</td>
</tr>
<tr>
<td>inf_sd_lower_bound</td>
<td>real(r8) dimension(2)</td>
<td>Lower bound for inflation standard deviation. If using a
negative value for <em class="code">sd_initial</em> this should
also be negative to preserve the setting.</td>
</tr>
<tr>
<td>inf_sd_max_change</td>
<td>real(r8) dimension(2)</td>
<td>For inflation type 5 (enhanced inflation), controls the maximum
change of the inflation standard deviation when adapting for the
next assimilation cycle. The value should be between 1.0 and 2.0.
1.0 prevents any changes, while 2.0 allows 100% change. For the
enhanced inflation option, if the standard deviation initial value
is equal to the standard deviation lower bound the standard
deviation will not adapt in time. See <a href=
"../../programs/filter/filter.html#Inflation">this section</a> for
a discussion of how the standard deviation adapts based on
different types of inflation.</td>
</tr>
<tr>
<td>trace_execution</td>
<td>logical</td>
<td>True means output very detailed messages about what routines
are being called in the main filter loop. Useful if a job hangs or
otherwise doesn't execute as expected.</td>
</tr>
<tr>
<td>output_timestamps</td>
<td>logical</td>
<td>True means write timing information to the log before and after
the model advance and the observation assimilation phases.</td>
</tr>
<tr>
<td>output_forward_op_errors</td>
<td>logical</td>
<td>True means output errors from forward observation operators.
This is the 'istatus' error return code from the model_interpolate
routine. An ascii text file <em class=
"file">prior_forward_op_errors</em> and/or <em class=
"file">post_forward_op_errors</em> will be created in the current
directory. For each ensemble member which returns a non-zero return
code, a line will be written to this file. Each line will have
three values listed: the observation number, the ensemble member
number, and the istatus return code. Be cautious when turning this
option on. The number of lines in this file can be up to the number
of observations times the number of ensemble members times the
number of assimilation cycles performed. This option is generally
most useful when run with a small observation sequence file and a
small number of ensemble members to diagnose forward operator
problems.</td>
</tr>
<tr>
<td>write_obs_every_cycle</td>
<td>logical</td>
<td>For debug use; this option can significantly slow the execution
of filter. True means to write the entire output observation
sequence diagnostic file each time through the main filter loop
even though only observations with times up to and including the
current model time will have been assimilated. Unassimilated
observations have the value -888888.0 (the DART "missing value").
If filter crashes before finishing it may help to see the forward
operator values of observations that have been assimilated so
far.</td>
</tr>
<tr>
<td>silence</td>
<td>logical</td>
<td>True means output almost no runtime messages. Not recommended
for general use, but can speed long runs of the lower order models
if the execution time becomes dominated by the volume of
output.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe the modules used by this program.                       -->
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
obs_sequence_mod
obs_def_mod
obs_def_utilities_mod
time_manager_mod
utilities_mod
assim_model_mod
assim_tools_mod
obs_model_mod
ensemble_manager_mod
adaptive_inflate_mod
mpi_utilities_mod
smoother_mod
random_seq_mod
state_vector_io_mod
io_filenames_mod
forward_operator_mod
quality_control_mod
</pre>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
See the <a href=
"../../programs/filter/filter.html#FilesUsed">filter overview</a>
for the list of files.
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%"
summary='error table'>
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">ens_size in namelist is ###: Must be &gt; 1</td>
<!-- comment -->
<td valign="top">Ensemble size must be at least 2.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">inf_flavor= ### Must be 0, 2, 3.</td>
<!-- comment -->
<td valign="top">Observation Inflation is no longer supported (i.e
flavor 1).</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">Posterior observation space inflation (type 1) not
supported.</td>
<!-- comment -->
<td valign="top">Posterior observation space inflation doesn't
work.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">Number of processes &gt; model size.</td>
<!-- comment -->
<td valign="top">Number of processes can't exceed model size for
now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_generate_copy_meta_data</td>
<!-- message -->
<td valign="top">output metadata in filter needs state ensemble
size &lt; 10000, not ###.</td>
<!-- comment -->
<td valign="top">Only up to 10000 ensemble members with state
output for now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_generate_copy_meta_data</td>
<!-- message -->
<td valign="top">output metadata in filter needs obs ensemble size
&lt; 10000, not ###.</td>
<!-- comment -->
<td valign="top">Only up to 10000 ensemble members with obs space
output for now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_setup_obs_sequence</td>
<!-- message -->
<td valign="top">input obs_seq file has ### qc fields; must be &lt;
2.</td>
<!-- comment -->
<td valign="top">Only 0 or 1 qc fields in input obs sequence for
now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_obs_copy_index</td>
<!-- message -->
<td valign="top">Did not find observation copy with metadata
observation.</td>
<!-- comment -->
<td valign="top">Only 0 or 1 qc fields in input obs sequence for
now.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>Many. New assimilation algorithms, support for new observation
types, support for additional models, better performance on higher
numbers of MPI tasks... The list is long. Send email to <a href=
"mailto:dart@ucar.edu">dart@ucar.edu</a> if you are interested in
additional functionality in DART.</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
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
