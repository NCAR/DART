<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>RMA notes</title>
<link rel="stylesheet" type="text/css" href="doc.css">
<link href="../images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>RMA notes</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../images/Dartboard7.png" alt=
"DART project logo" height="70"></td>
<td>Jump to <a href="../index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
In the RMA version of DART, the state vector is not required to be
stored completely on any process. This is achieved using Remote
Memory Access (RMA). The RMA programing model allows processes to
read (and write) memory on other processors asynchronously. RMA
DART supported models:
<ul>
<li>9var</li>
<li>bgrid_solo</li>
<li>cam-fv</li>
<li>cice</li>
<li>cm1</li>
<li>lorenz_04</li>
<li>lorenz_63</li>
<li>lorenz_84</li>
<li>lorenz_96</li>
<li>mpas_atm (NetCDF overwrite not supported for
update_u_from_reconstruct = .true. )</li>
<li>POP</li>
<li>ROMS</li>
<li>wrf</li>
</ul>
There are six major differences between Lanai and RMA DART:
<ul>
<li>Read and write NetCDF restarts</li>
<li>Calculation of forward operators</li>
<li>Vertical conversion of observation locations</li>
<li>Diagnostic file changes</li>
<li><a href="state_structure.html">State structure module</a></li>
<li>Perturbation of the state</li>
</ul>
<p>Before bitwise testing with Lanai please read <a href=
"bitwise_considerations.html">bitwise considerations</a></p>
<h2>NetCDF Restarts</h2>
The programs filter and perfect_model_obs now read/write directly
from NetCDF files, rather than having to run converters
(<code>model_to_dart</code> and <code>dart_to_model</code>). To
facilitate this, there is a new required call
<code>add_domain</code> which must be called during
<code>static_init_model</code>. It can be called multiple times in
static_model_mod, e.g. once for each NetCDF file that contains
state variables. There are three ways to add a domain:
<ul>
<li><b>From Blank</b> : This is for small models such as lorenz_96
and no NetCDF restarts
<ul>
<li><code>dom_id = add_domain(model_size)</code></li>
</ul>
</li>
<li><b>From File</b> : This is for models which have NetCDF restart
files
<ul>
<li><code>dom_id = add_domain(template_file, num_vars, var_names,
... )</code></li>
</ul>
</li>
<li><b>From Spec</b> : Creates a skeleton structure for a domain (
currently only used in bgrid_solo )
<ul>
<li><code>dom_id = add_domain(num_vars, var_names, ... )</code><br>
<code>call add_dimension_to_variable(dom_id, var_id, dim_nam,
dim_size)</code><br>
<code>call finished_adding_domain</code></li>
</ul>
</li>
</ul>
<p>For models without NetCDF restarts, use
<code>add_domain(model_size)</code>. This is the minimum amount of
information needed by DART to create a netdcf file. For models with
NetCDF restarts use <code>add_domain(info_file, num_vars,
var_names)</code> which lets DART read the NetCDF dimensions for a
list of variables from a file (<code>info_file</code>). There are
several routines that can be used together to create a domain from
a description: <code>add_domain, add_dimension_to_variable,
finished_adding_domain</code>. This can be used in models such as
bgrid_solo where the model is spun up in perfect_model_obs, but the
model itself has variable structure (3D variables with names). See
<a href="#namelist_changes">Additions/Changes to existing namelists
for how to use NetCDF IO</a>.</p>
<p><b>Note</b> when using NetCDF restarts, inflation files are
NetCDF also. The inflation mean and inflation standard deviation
are in separate files when you use NetCDF restarts. See <a href=
"netcdf_inflation_files.html">NetCDF inflation files</a> for
details.</p>
<h2>Calculation of forward operators</h2>
<p>The forward operator code in model_mod now operates on an array
of state values. See <a href="forward_operator.html">forward
operator</a> for more detail about distributed vs. non-distributed
forward operators.</p>
<p>In distributed mode the forward operators for all ensemble
members are calculated in the same <code>model_interpolate</code>
call. In non-distributed mode, the forward oeprators for all
ensemble members a task owns (1-ens_size) are calculated at
once.</p>
<h2>Vertical conversion of observation locations</h2>
<p>The vertical conversion of observation locations is done before
the assimilation. In Lanai this calculation is done in the
assimilation as part of <code>get_close_obs</code> if a model_mod
does vertical conversion. See <a href=
"vertical_conversion.html">vertical conversion</a> for details
about this change. Note that not all models do vertical conversion
or even have a concept of vertical location, but every model_mod
must have the following routines:</p>
<ul>
<li><code>query_vert_localization_coord</code> - returns the
vertical localization coordiate (or does nothing if there is no
vertical conversion)</li>
<li><code>vert_convert</code> - converts location to required
vertical (or does nothing if there is no vertical conversion)</li>
</ul>
<h2>Diagnostic file changes</h2>
<p>For large models DART format diagnostic files (Prior_Diag.nc and
Posterior_Diag.nc) have been replaced with separate files for each
copy that would have gone into Prior_Diag.nc and
Posterior_Diag.nc.</p>
<p>For Prior_Diag.nc:</p>
<ul>
<li><b>Mean and standard deviation</b>:<br>
  preassim_mean.nc<br>
  preassim_sd.nc</li>
<li><b>Inflation mean and standard deviation</b> (if state space
inflation is used):<br>
  preassim_priorinf_mean.nc<br>
  preassim_priorinf_sd.nc</li>
<li><b>The number of ensemble members specifed</b> in filter_nml
(num_output_state_members):<br>
  preassim_member_####.nc</li>
</ul>
<p>For Posterior_Diag.nc:</p>
<ul>
<li><b>Mean and standard deviation</b>:<br>
  postassim_mean.nc<br>
  postassim_sd.nc</li>
<li><b>Inflation mean and standard deviation</b> (if state space
inflation is used):<br>
  postassim_priorinf_mean.nc<br>
  postassim_priorinf_sd.nc</li>
<li><b>The number of ensemble members specifed</b> in filter_nml
(num_output_state_members):<br>
  postassim_member_####.nc</li>
</ul>
<p>The <code>num_output_state_members</code> are not written
separately from the restarts. Note that restarts will have been
clamped if any clamping is applied (given as an arguement to
add_domain). This is <em>different</em> to Posterior_Diag.nc which
contains unclamped values. Note also that there are 2 more <a href=
"#STAGES_TO_WRITE">"stages"</a> which might be output, in addition
to the preassim and postassim discussed here.</p>
<p>For models with multiple domains the filenames above are
appended with the domain number, e.g. preassim_mean.nc becomes
preassim_mean_d01.nc, preassim_mean_d02.nc, etc.</p>
<h3>Changes to nc_write_model_atts</h3>
<code>nc_write_model_atts</code> has an new argument
'<code>model_mod_writes_state_variables</code>'. This is used to
communicate to DART whether the model will create and write state
variables in Prior_Diag.nc and Posterior_Diag.nc. If
<code>model_model_writes_state_variables = .false.</code> DART will
define and write state variables to the new diagnostic files. If
<code>model_model_writes_state_variables = .true.,
nc_write_model_vars</code> is called as normal.
<h2>Perturbations</h2>
The option to perturb one ensemble member to produce an ensemble is
in filter_nml:<code>perturb_from_single_instance</code>. The
model_mod interface is now <code>pert_model_copies</code> not
<code>pert_model_state</code>. Each task perturbs every ensemble
member for its own subsection of state. This is more complicated
than the Lanai routine <code>pert_model_state</code>, where a whole
state vector is available. If a model_mod does not provide a
perturb interface, filter will do the perturbing with an amplitude
set in filter_nml:perturbation_amplitude. Note the perturb namelist
options have been removed from ensemble_manager_nml
<h2>state_vector_io_nml</h2>
<div class="namelist">
<pre>
&amp;state_vector_io_nml
   buffer_state_io         = .false.,
   single_precision_output = .false.,
/
</pre></div>
<p>When <code>buffer_state_io</code> is <code>.false.</code> the
entire state is read into memory at once if .true. variables are
read one at a time. If your model can not fit into memory at once
this must be set to <code>.true.</code> .</p>
<p><code>single_precision_output</code> allows you to run filter in
double precision but write NetCDF files in single precision.</p>
<h2>quality_control_nml</h2>
These namelist options used to be in filter_nml, now they are in
quality_control_nml.
<div class="namelist">
<pre>
&amp;quality_control_nml
   input_qc_threshold          = 3,
   outlier_threshold           = 4,
   enable_special_outlier_code = .false.
/
</pre></div>
<a name="namelist_changes" id="namelist_changes"></a>
<h2>Additions/Changes to existing namelists</h2>
New namelist variables
<h3>filter_nml</h3>
<div class="namelist">
<pre>
&amp;filter_nml
   single_file_in               = .false.,
   single_file_out              = .false.,

   input_state_file_list        = 'null',
   output_state_file_list       = 'null',
   input_state_files            = 'null',
   output_state_files           = 'null',

   stages_to_write              = 'output'
   write_all_stages_at_end      = .false.
   output_restarts              = .true.
   output_mean                  = .true.
   output_sd                    = .true.

   perturb_from_single_instance = .false.,
   perturbation_amplitude       = 0.2_r8,

   distributed_state            = .true.
/
</pre></div>
<br>
<div><a name="STAGES_TO_WRITE" id="STAGES_TO_WRITE"></a>
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
<td>True means that all of the restart and inflation information is
read from a single NetCDF file. False means that you must specify
an input_state_file_list and DART will be expecting
input_{priorinf,postinf}_{mean,sd}.nc files for inflation.</td>
</tr>
<tr>
<td>single_file_out</td>
<td>logical</td>
<td>True means that all of the restart and inflation information is
written to a single NetCDF file. False means that you must specify
a output_state_file_list and DART will be output files specified in
the list. Inflation files will be written in the form
input_{priorinf,postinf}_{mean,sd}.nc.</td>
</tr>
<tr>
<td>input_restart_files</td>
<td>character array</td>
<td>This is used for single file input for low order models. For
multiple domains you can specify a file for each domain. When
specifying a list single_file_in, single_file_out must be set to
.true.</td>
</tr>
<tr>
<td>output_restart_files</td>
<td>character array</td>
<td>This is used for single file input for low order models. For
multiple domains you can specify a file for each domain. When
specifying a list single_file_in, single_file_out must be set to
.true.</td>
</tr>
<tr>
<td>input_state_file_list</td>
<td>character array</td>
<td>A list of files containing input model restarts. For multiple
domains you can specify a file for each domain. When specifying a
list single_file_in, single_file_out must be set to .false.</td>
</tr>
<tr>
<td>output_state_file_list</td>
<td>character array</td>
<td>A list of files containing output model restarts. For multiple
domains you can specify a file for each domain. When specifying a
list single_file_in, single_file_out must be set to .false.</td>
</tr>
</tbody>
<tbody valign="top">
<tr>
<td>stages_to_write</td>
<td>character array</td>
<td>Controls which stages to write. Currently there are four
options:
<ul>
<li><code>input</code> -- writes input mean and sd only</li>
<li><code>preassim</code> -- before assimilation, before prior
inflation is applied</li>
<li><code>postassim</code> -- after assimilation, before posterior
inflation is applied</li>
<li><code>output</code> -- final output for filter which includes
clamping and inflation</li>
</ul>
</td>
</tr>
<tr>
<td>write_all_stages_at_end</td>
<td>logical</td>
<td>True means output all stages at the end of filter. This is more
memory intensive but requires less time. For larger models IO
begins to dominate the overall cost of the assimilation, so
writting all stages at the end writes more files in parallel,
reducing the IO time. <code>output_state_file_list</code>.</td>
</tr>
<tr>
<td>output_restarts</td>
<td>logical</td>
<td>True means output a restart file(s). Filenames are defined in
<code>output_state_file_list</code>.</td>
</tr>
<tr>
<td>output_mean</td>
<td>logical</td>
<td>True means output a restart file which contains the ensemble
mean for the stages that have been turned on in
<code>stages_to_write</code>. The file name will have the stage
with <em class="code">_mean</em> appended.</td>
</tr>
<tr>
<td>output_sd</td>
<td>logical</td>
<td>True means output a restart file which contains the ensemble
standard deviation for the stages that have been turned on in
<code>stages_to_write</code>. The file name will have the stage
with <em class="code">_sd</em> appended.</td>
</tr>
<tr>
<td>perturb_from_single_instance</td>
<td>logical</td>
<td>Read a single file and perturb this to create an ensemble</td>
</tr>
<tr>
<td>perturbation_amplitude</td>
<td>float</td>
<td>Perturbation amplitude</td>
</tr>
<tr>
<td>distribute_state</td>
<td>logical</td>
<td>True keeps the state distributed across all tasks throughout
the entire execution of filter.</td>
</tr>
</tbody>
</table>
</div>
<p><b>For NetCDF reads and writes</b></p>
<p>For <b>input</b> file names:</p>
<ul>
<li>give <code>input_state_file_list</code> a file for each domain,
each of which contains a list of restart files.</li>
<li>if no <code>input_state_file_list</code> is provided then
default filenames will be used e.g. input_member_000*.nc,
input_priorinf_mean.nc, input_priorinf_sd.nc</li>
</ul>
<p>For <b>output</b> file names:</p>
<ul>
<li>give <code>output_state_file_list</code> a file for each
domain, each of which contains a list of restart files.</li>
<li>if no <code>output_state_file_list</code> is provided then
default filenames will be used e.g. output_member_000*.nc,
output_priorinf_mean.nc, output_priorinf_sd.nc</li>
</ul>
For small models you may want to use <code>single_file_in</code>,
<code>single_file_out</code> which contains all copies needed to
run filter. <!--
<p>
<b>For Lanai style dart input and output restart files:</b> <br>
<ul>
<li> output_state_file_list and restart_out_file_name are used in the same way as Lanai.
</ul>
-->
<h3>assim_tools_nml</h3>
<div class="namelist">
<pre>
&amp;assim_tools_nml
   distribute_mean  = .true.
/
</pre></div>
<p>In previous DART releases, each processor gets a copy of the
mean (in ens_mean_for_model). In RMA DART, the mean is distributed
across all processors. However, a user can choose to have a copy of
the mean on each processor by setting <code>distribute_mean =
.false.</code> . Note that the mean state is accessed through
<code>get_state</code> whether distribute_mean is
<code>.true.</code> or <code>.false.</code></p>
<h3>Removed from existing namelists</h3>
<pre>
&amp;filter_nml
   input_qc_threshold          = 3,
   outlier_threshold           = 4,
   enable_special_outlier_code = .false.
   start_from_restart          = .false.
   output_inflation            = .true.
   output_restart              = .true.
   /
</pre>
NOTE : <code>output_restart</code> has been renamed to
<code>output_restarts</code>. <b><code>output_inflation</code> is
no longer supported</b> and only writes inflation files if
<code>inf_flavor &gt; 1</code>
<pre>
&amp;ensemble_manager_nml
   single_restart_file_out = .true.
   perturbation_amplitude  = 0.2,
   /
</pre>
<pre>
&amp;assim_manager_nml
   write_binary_restart_files = .true.,
   netCDF_large_file_support  = .false.
   /
</pre>
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
