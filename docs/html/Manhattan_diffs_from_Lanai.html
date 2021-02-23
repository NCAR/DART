<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>DART Manhattan Differences from Lanai Release Notes</title>
<link rel="stylesheet" type="text/css" href="doc.css">
<link href="../images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>DART Manhattan Differences from Lanai Release Notes</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../images/Dartboard7.png" alt=
"DART project logo" height="70"></td>
<td>Jump to <a href="../index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Overview">Overview</a> / 
<!-- <a href="#NewModels">New Models or Changes to Existing Models</a> / -->
<a href="#NetcdfRestarts">NetCDF Restart Files</a> / <a href=
"#ForwardOps">Forward Operators</a> / <a href=
"#ObsConvert">Vertical Conversion</a> / <a href="#Diagnostics">DART
Diagnostics Changes</a> / <a href="#ModelMod">model_mod.f90
Interface Changes</a> / <a href="#ObsQuantity">Observation
Quantities</a> / <a href="#NameListChanges">Additions/Changes to
Existing Namelists</a> / <a href="#Perturbations">Perturbations</a>
/ <!-- <a href="#Misc">Tutorial, Scripting, Setup, Builds</a> / -->
<a href="#Legalese">Terms of Use</a> 
<!--==================================================================-->
 <a name="Overview" id="Overview"></a>
<h2>Overview</h2>
<p>This document includes an overview of the changes in the DART
system since the Lanai release. For further details on any of these
items look at the HTML documentation for that specific part of the
system.</p>
<p>The two most significant changes in the Manhattan version of
DART are it can support running models with a state vector larger
than the memory of a single task, removing a limit from the Lanai
version of DART. It also reads and writes NetCDF files directly
instead of requiring a conversion from one file to another. There
are many other smaller changes, detailed below.</p>
<p>Manhattan supported models:</p>
<ul>
<li>9var</li>
<li>bgrid_solo</li>
<li>cam-fv</li>
<li>cice</li>
<li>clm</li>
<li>cm1</li>
<li>forced_lorenz_96</li>
<li>ikeda</li>
<li>lorenz_63</li>
<li>lorenz_84</li>
<li>lorenz_96</li>
<li>lorenz_96_2scale</li>
<li>lorenz_04</li>
<li>mpas_atm (NetCDF overwrite not supported for
update_u_from_reconstruct = .true. )</li>
<li>null_model</li>
<li>POP</li>
<li>ROMS</li>
<li>simple_advection</li>
<li>wrf</li>
</ul>
<p>If your model of interest is not on the list consider checking
out the 'Classic' release of DART, which is Lanai plus bug fixes
and minor enhancements. All models previously supported by Lanai
are still in DART 'Classic'.</p>
<p>These are the major differences between the Lanai/Classic and
Manhattan releases of DART:</p>
<ul>
<li>Read and write NetCDF restarts</li>
<li>Calculation of forward operators</li>
<li>Vertical conversion of observation locations</li>
<li>Diagnostic file changes</li>
<li><a href="state_structure.html">State structure module</a></li>
<li>model_mod interface changes</li>
<li>Observation Quantity replaces Kind</li>
<li>Perturbation of the state</li>
</ul>
<!--==================================================================-->
<a name="NetcdfRestarts" id="NetcdfRestarts"></a>
<h2>NetCDF Restart Files</h2>
The programs filter and perfect_model_obs now read/write directly
from NetCDF files rather than having to run converters
(<code>model_to_dart</code> and <code>dart_to_model</code>). To
facilitate this there is a new required call
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
details. 
<!--==================================================================-->
 <a name="ForwardOps" id="ForwardOps"></a></p>
<h2>Calculation of Forward Operators</h2>
The forward operator code in model_mod now operates on an array of
state values. See <a href="forward_operator.html">forward
operator</a> for more detail about distributed vs. non-distributed
forward operators. In distributed mode the forward operators for
all ensemble members are calculated in the same
<code>model_interpolate</code> call. In non-distributed mode, the
forward operators for all ensemble members a task owns (1-ens_size)
are calculated at once. 
<!--==================================================================-->
 <a name="ObsConvert" id="ObsConvert"></a>
<h2>Vertical Conversion of Observation and State Locations</h2>
<p>The vertical conversion of observation locations is done before
the assimilation by default. This can be changed by namelist
options.</p>
<p>In Lanai this calculation is done in the assimilation as part of
<code>get_close_obs</code> if a model_mod does vertical conversion.
See <a href="vertical_conversion.html">vertical conversion</a> for
details about this change. Note that not all models do vertical
conversion or even have a concept of vertical location, but every
model_mod must have the following routines:</p>
<pre>
call set_vertical_localization_coord(vert_localization_coord)

call convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &amp;
                          which_vert, status)

call convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &amp;
                            which_vert, istatus)
</pre>
<p>If there are NOT multiple choices for a vertical coordinate
(e.g. cartesian, one dimensional), all these routines can be
no-ops.</p>
<p>If there are multiple types of vertical coordinates, the convert
routines must be able to convert between them. The
'set_vertical_localization_coord()' routine should be called from
'static_init_model()' to set what localization coordinate type is
being requested.</p>
<p>The three routines related to vertical coordinates/localization
choices are:</p>
<ul>
<li><code>set_vert_localization_coord</code> - sets the vertical
localization coordiate (not required if there is no vertical
conversion)</li>
<li><code>convert_vertical_obs</code> - converts observation
location to required vertical type (does nothing if there is no
vertical conversion)</li>
<li><code>convert_vertical_state</code> - converts state vector
location to required vertical type (does nothing if there is no
vertical conversion)</li>
</ul>
<p>
<!--==================================================================-->
 <a name="Diagnostics" id="Diagnostics"></a></p>
<h2>DART Diagnostic file changes</h2>
<p>For large models DART format diagnostic files
(<code>Prior_Diag.nc</code> and <code>Posterior_Diag.nc</code>)
have been replaced with separate files for each copy that would
have gone into Prior_Diag.nc and Posterior_Diag.nc.</p>
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
<p><code>nc_write_model_atts</code> now has 2 arguments:</p>
<ul>
<li>ncid - open netcdf file identifier</li>
<li>domain_id - domain number being written</li>
</ul>
<p>The calling code will write the model state, so this routine
should only add attributes and optionally, non-state information
like grid arrays.</p>
<p>This routine will only be called if DART is creating an output
NetCDF file from scratch. This may include any of the preassim,
postassim, or output files.</p>
<h3>Changes to nc_write_model_vars</h3>
<p><code>nc_write_model_vars</code> is currently unused (and in
fact uncalled). It remains for possible future expansion. 
<!--==================================================================-->
 <a name="ModelMod" id="ModelMod"></a></p>
<h2>model_mod.f90 Interface Changes</h2>
<p>The model_mod.f90 file contains all code that is specific to any
particular model. The code in this file is highly constrained since
these routines are *called by* other code in the DART system. All
routine interfaces -- the names, number of arguments, and the names
of those arguments -- must match the prescribed interfaces exactly.
Since not all required interfaces are needed for every model there
are default routines provided that can be referenced from a 'use'
statement and then the routine name can be put in the module
'public' list without any code for that routine having to be
written in the model_mod.f90 file.</p>
<p>The following 18 routines are required:</p>
<ul>
<li>static_init_model</li>
<li>get_model_size</li>
<li>get_state_meta_data</li>
<li>shortest_time_between_assimilations</li>
<li>model_interpolate</li>
<li>end_model</li>
<li>nc_write_model_atts</li>
<li>nc_write_model_vars</li>
<li>init_time</li>
<li>init_conditions</li>
<li>adv_1step</li>
<li>pert_model_copies</li>
<li>get_close_obs</li>
<li>get_close_state</li>
<li>convert_vertical_obs</li>
<li>convert_vertical_state</li>
<li>read_model_time</li>
<li>write_model_time</li>
</ul>
<!-- FIXME:  add more info about  exactly what these routines
need to do here. -->
<p>Here is an example of code from the top of a model_mod file,
including the modules where the default routines live and the
required public list.</p>
<pre>
<code>
use     location_mod, only : location_type, get_close_type, &amp;
                             get_close_obs, get_close_state, &amp;
                             convert_vertical_obs, convert_vertical_state, &amp;
                             set_location, set_location_missing, &amp;
                             set_vertical_localization_coord
use    utilities_mod, only : register_module, error_handler, &amp;
                             E_ERR, E_MSG
                             ! nmlfileunit, do_output, do_nml_file, do_nml_term,  &amp;
                             ! find_namelist_in_file, check_namelist_read
use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &amp;
                                 nc_add_global_creation_time, &amp;
                                 nc_begin_define_mode, nc_end_define_mode
use state_structure_mod, only : add_domain
use ensemble_manager_mod, only : ensemble_type
use dart_time_io_mod, only  : read_model_time, write_model_time
use default_model_mod, only : pert_model_copies, nc_write_model_vars

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: static_init_model,      &amp;
          get_model_size,         &amp;
          get_state_meta_data,    &amp;
          shortest_time_between_assimilations, &amp;
          model_interpolate,      &amp;
          end_model,              &amp;
          nc_write_model_atts,    &amp;
          adv_1step,              &amp;
          init_time,              &amp;
          init_conditions

! public but in another module
public :: nc_write_model_vars,    &amp;
          pert_model_copies,      &amp;
          get_close_obs,          &amp;
          get_close_state,        &amp;
          convert_vertical_obs,   &amp;
          convert_vertical_state, &amp;
          read_model_time,        &amp;
          write_model_time
</code>
</pre>
<!--==================================================================-->
<a name="ObsQuantity" id="ObsQuantity"></a>
<h2>Observation Quantity replaces Kinds</h2>
<p>Historically there has been confusion about the terms for
specific observation types (which often include the name of the
instrument collecting the data) and the generic quantity that is
being measured (e.g. temperature). The previous terms for these
were 'types' and 'kinds', respectively.</p>
<p>Starting with the Manhattan release we have tried to clarify the
terminology and make the interfaces consistent. The following table
lists the original names from the Lanai/Classic release and the
replacement routines in Manhattan.</p>
<p>All code that is part of the DART code repository has been
updated to use the replacment routines, but if you have your own
utilities written using this code, you will need to update your
code. Contact us ( <a href="mailto:dart@ucar.edu">dart@ucar.edu</a>
) for help if you have any questions.</p>
<!-- FIXME:  make this into a table and add some more info -->
<pre>
<code>
public subroutines, existing name on left, replacement on right:

assimilate_this_obs_kind()     =&gt;     assimilate_this_type_of_obs(type_index)
evaluate_this_obs_kind()       =&gt;       evaluate_this_type_of_obs(type_index)
use_ext_prior_this_obs_kind()  =&gt;  use_ext_prior_this_type_of_obs(type_index)

get_num_obs_kinds()      =&gt;  get_num_types_of_obs()
get_num_raw_obs_kinds()  =&gt;  get_num_quantities()

get_obs_kind_index()     =&gt; get_index_for_type_of_obs(type_name)
get_obs_kind_name()      =&gt; get_name_for_type_of_obs(type_index)

get_raw_obs_kind_index()  =&gt;  get_index_for_quantity(quant_name)
get_raw_obs_kind_name()   =&gt;  get_name_for_quantity(quant_index)

get_obs_kind_var_type()  =&gt;  get_quantity_for_type_of_obs(type_index)

get_obs_kind()      =&gt;  get_obs_def_type_of_obs(obs_def)
set_obs_def_kind()  =&gt;  set_obs_def_type_of_obs(obs_def)

get_kind_from_menu()      =&gt;  get_type_of_obs_from_menu()

read_obs_kind()     =&gt;   read_type_of_obs_table(file_unit, file_format)
write_obs_kind()    =&gt;  write_type_of_obs_table(file_unit, file_format)

maps obs_seq nums to specific type nums, only used in read_obs_seq:
map_def_index()  =&gt; map_type_of_obs_table()  

removed.  apparently unused, and simply calls get_obs_kind_name():
get_obs_name()

apparently unused anywhere, removed:
add_wind_names()
do_obs_form_pair()

public integer parameter constants and subroutine formal argument names,
old on left, new on right:

KIND_ =&gt; QTY_
kind =&gt; quantity

TYPE_ =&gt; TYPE_
type =&gt; type_of_obs

integer parameters:
max_obs_generic  =&gt;  max_defined_quantities  (not currently public, leave private)
max_obs_kinds    =&gt;  max_defined_types_of_obs 
</code>
</pre>
<!--==================================================================-->
<a name="NamelistChanges" id="NamelistChanges"></a>
<h2>Additions/Changes to Existing Namelists</h2>
<h3>quality_control_nml</h3>
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
a output_state_files and DART will be output files specified in the
list. Inflation files will be written in the form
input_{priorinf,postinf}_{mean,sd}.nc.</td>
</tr>
<tr>
<td>input_state_files</td>
<td>character array</td>
<td>This is used for single file input for low order models. For
multiple domains you can specify a file for each domain. When
specifying a list single_file_in, single_file_out must be set to
.true.</td>
</tr>
<tr>
<td>output_state_files</td>
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
reducing the IO time. Filenames are defined in
<code>output_state_files</code>.</td>
</tr>
<tr>
<td>output_restarts</td>
<td>logical</td>
<td>True means output a restart file(s). Filenames are defined in
<code>output_state_files</code>.</td>
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
<p><b>NetCDF reads and writes:</b></p>
<p>For <b>input</b> file names:</p>
<ul>
<li>give <code>input_state_file_list</code> a file for each domain,
each of which contains a list of restart files. An example of an
'input_list.txt' might look something like :<br>
<div class="namelist">
<pre>
advance_temp1/wrfinput_d01
advance_temp2/wrfinput_d01
advance_temp3/wrfinput_d01
advance_temp4/wrfinput_d01
advance_temp5/wrfinput_d01
....
</pre></div>
<br></li>
<li>if no <code>input_state_file_list</code> is provided then
default filenames will be used e.g. input_member_####.nc,
input_priorinf_mean.nc, input_priorinf_sd.nc</li>
</ul>
<p>For <b>output</b> file names:</p>
<ul>
<li>give <code>output_state_file_list</code> a file for each
domain, each of which contains a list of restart files. An example
of an 'input_list.txt' might for WRF might look something like
:<br>
<div class="namelist">
<pre>
wrf_out_d01.0001.nc
wrf_out_d01.0002.nc
wrf_out_d01.0003.nc
wrf_out_d01.0004.nc
wrf_out_d01.0005.nc
....
</pre></div>
<br>
if you would like to simply like to overwrite your previous data
input_list.txt = output_list.txt</li>
<li>if no <code>output_state_files</code> is provided then default
filenames will be used e.g. output_member_####.nc,
output_priorinf_mean.nc, output_priorinf_sd.nc</li>
</ul>
For small models you may want to use <code>single_file_in</code>,
<code>single_file_out</code> which contains all copies needed to
run filter. <!--
<p>
<b>For Lanai style dart input and output restart files:</b> <br>
<ul>
<li> output_state_files and restart_out_file_name are used in the same way as Lanai.
</ul>
-->
<h3>state_vector_io_nml</h3>
<div class="namelist">
<pre>
&amp;state_vector_io_nml
   buffer_state_io          = .false.,
   single_precision_output  = .false.,
/
</pre></div>
<p>When <code>buffer_state_io</code> is <code>.false.</code> the
entire state is read into memory at once if .true. variables are
read one at a time. If your model can not fit into memory at once
this must be set to <code>.true.</code> .</p>
<p><code>single_precision_output</code> allows you to run filter in
double precision but write NetCDF files in single presision</p>
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
<a name="Perturbations" id="Perturbations"></a>
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
