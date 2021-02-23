<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module assim_model_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE assim_model_mod</h1>
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
<p>This module acts as an intermediary between DART compliant
models and the filter. At one time the assim_model_type, which
combines a state vector and a time_type, was envisioned as being
fundamental to how DART views model states. This paradigm is
gradually being abandoned so that model state vectors and times are
handled as separate data types. It is important to call
static_init_assim_model before using routines in assim_model_mod.
Interfaces to work with model time stepping, restart files, and
computations about the locations of model state variables and the
distance between observations and state variables. Many of the
interfaces are passed through nearly directly to the model_mod.</p>
<h3 class="indent1">NOTES</h3>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This module does not have a namelist.</p>
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
location_mod (model dependent choice)
time_manager_mod
utilities_mod
model_mod
netcdf
typeSizes (part of netcdf)
</pre>
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table>
<tr>
<td><em class="call">use assim_model_mod, only :</em></td>
</tr>
<tr>
<td> </td>
<td><a href="#adv_1step">adv_1step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#aoutput_diagnostics">aoutput_diagnostics</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#aread_state_restart">aread_state_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#assim_model_type">assim_model_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#awrite_state_restart">awrite_state_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#close_restart">close_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#copy_assim_model">copy_assim_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#end_assim_model">end_assim_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ens_mean_for_model">ens_mean_for_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#finalize_diag_output">finalize_diag_output</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_close_maxdist_init">get_close_maxdist_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs">get_close_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs_init">get_close_obs_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_closest_state_time_to">get_closest_state_time_to</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_diag_input_copy_meta_data">get_diag_input_copy_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_initial_conditions">get_initial_condition</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_model_size">get_model_size</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_model_state_vector">get_model_state_vector</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_model_time">get_model_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_model_time_step">get_model_time_step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_state_meta_data">get_state_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_assim_model">init_assim_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_diag_input">init_diag_input</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_diag_output">init_diag_output</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#input_diagnostics">input_diagnostics</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#interpolate">interpolate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_append_time">nc_append_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_get_tindex">nc_get_tindex</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#nc_write_calendar_atts">nc_write_calendar_atts</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#netcdf_file_type">netcdf_file_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#open_restart_read">open_restart_read</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#open_restart_write">open_restart_write</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#output_diagnostics">output_diagnostics</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#pert_model_state">pert_model_state</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_state_restart">read_state_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#set_model_state_vector">set_model_state_vector</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_model_time">set_model_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#static_init_assim_model">static_init_assim_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_state_restart">write_state_restart</a></td>
</tr>
</table>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="assim_model_type" id="assim_model_type"></a><br>
<div class="routine">
<pre>
<em class="call">type assim_model_type</em>
   private
   real(r8), pointer   :: state_vector(:) 
   type(time_type)     :: time
   integer             :: model_size
   integer             :: copyID
end type assim_model_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>This type is used to represent both the state and time of a
state from a model.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">state_vector</td>
<td>A one dimensional representation of the model state
vector.</td>
</tr>
<tr>
<td valign="top">time</td>
<td>The time of the model state.</td>
</tr>
<tr>
<td valign="top">model_s</td>
<td>Size of the model state vector.</td>
</tr>
<tr>
<td valign="top">copyID</td>
<td>Not used in present implementation.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="netcdf_file_type" id="netcdf_file_type"></a><br>
<div class="routine">
<pre>
<em class="call">type netcdf_file_type</em>
   integer             :: ncid
   integer             :: Ntimes
   integer             :: NtimesMAX
   real(r8), pointer   :: rtimes(:)
   type(time_type), pointer :: times(:)
   character(len = 80)      :: fname
end type netcdf_file_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Basically, we want to keep a local mirror of the unlimited
dimension coordinate variable (i.e. time) because dynamically
querying it causes unacceptable performance degradation over "long"
integrations.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">ncid</td>
<td>The netcdf file unit id.</td>
</tr>
<tr>
<td valign="top">Ntimes</td>
<td>The current working length.</td>
</tr>
<tr>
<td valign="top">NtimesMAX</td>
<td>Allocated length.</td>
</tr>
<tr>
<td valign="top">rtimes</td>
<td>Times as real (r8).</td>
</tr>
<tr>
<td valign="top">times</td>
<td>Times as time_types.</td>
</tr>
<tr>
<td valign="top">fname</td>
<td>Netcdf file name.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="static_init_assim_model" id=
"static_init_assim_model"></a><br>
<div class="routine"><em class="call">call
static_init_assim_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Initializes the assim_model class. Must be called before any
other assim_model_mod interfaces are used. Also calls the static
initialization for the underlying model. There are no
arguments.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_diag_output" id="init_diag_output"></a><br>
<div class="routine"><em class="call">ncFileID =
init_diag_output(FileName, global_meta_data,
copies_of_field_per_time, meta_data_per_copy <em class=
"optionalcode">[, lagID]</em>)</em>
<pre>
type(netcdf_file_type)          :: <em class=
"code">init_diag_output </em>
character (len = *), intent(in) :: <em class="code">FileName </em>
character (len = *), intent(in) :: <em class=
"code">global_meta_data </em>
integer, intent(in)             :: <em class=
"code">copies_of_field_per_time </em>
character (len = *), intent(in) :: <em class=
"code">meta_data_per_copy(copies_of_field_per_time) </em>
integer, optional, intent(in)   :: <em class=
"optionalcode">lagID </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initializes a netCDF file for output of state space diagnostics.
A handle to the channel on which the file is opened is
returned.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>Identifier for the netcdf file is returned. This is not an
integer unit number, but a derived type containing additional
information about the opened file.</td>
</tr>
<tr>
<td valign="top"><em class="code">FileName</em></td>
<td>Name of file to open.</td>
</tr>
<tr>
<td valign="top"><em class="code">global_meta_data</em></td>
<td>Global metadata that describes the contents of this file.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">copies_of_field_per_time   </em></td>
<td>Number of copies of data to be written at each time. For
instance, these could be the prior ensemble members, prior ensemble
mean, prior ensemble spread, posterior ensemble members, posterior
spread and mean, etc..</td>
</tr>
<tr>
<td valign="top"><em class="code">meta_data_per_copy</em></td>
<td>Metadata describing each of the copies.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">lagID</em></td>
<td>If using the smoother, which lag number this output is
for.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_model_size" id="get_model_size"></a><br>
<div class="routine"><em class="call">var = get_model_size()</em>
<pre>
integer :: <em class="code">get_model_size </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the size of the model state vector. This is a direct
pass through to the model_mod.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_closest_state_time_to" id=
"get_closest_state_time_to"></a><br>
<div class="routine"><em class="call">var =
get_closest_state_time_to(model_time, time)</em>
<pre>
type(time_type)              :: <em class=
"code"> get_closest_state_time_to </em>
type(time_type), intent(in)  :: <em class="code"> model_time </em>
type(time_type), intent(in)  :: <em class="code"> time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the closest time that a model is capable of advancing a
given state to a specified time. For instance, what is the closest
time to 12GMT 01 January, 2004 that a model state at 00GMT 01
January, 2004 can be advanced? If the model time is past the time,
the model time is returned (new feature in releases after
Hawaii).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>The closest time to which the model can be advanced is
returned.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">model_time   </em></td>
<td>The time of a model state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">time</em></td>
<td>A time that one would like to get close to with the model.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_state_meta_data" id="get_state_meta_data"></a><br>
<div class="routine"><em class="call">call
get_state_meta_data()</em></div>
<div class="indent1"><!-- Description -->
<p>Pass through to model_mod. See model_mod documentation for
arguments and description.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_model_time" id="get_model_time"></a><br>
<div class="routine"><em class="call">var =
get_model_time(assim_model)</em>
<pre>
type(time_type)                    :: <em class=
"code">get_model_time</em>
type(assim_model_type), intent(in) :: <em class=
"code">assim_model</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns time from an assim_model type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returned time from assim_model</td>
</tr>
<tr>
<td valign="top"><em class=
"code">assim_model   </em></td>
<td>Assim_model type from which to extract time</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_model_state_vector" id=
"get_model_state_vector"></a><br>
<div class="routine"><em class="call">var =
get_model_state_vector(assim_model)</em>
<pre>
real(r8)                           :: <em class=
"code">get_model_state_vector(model_size)</em>
type(assim_model_type), intent(in) :: <em class=
"code">assim_model</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the state vector component from an assim_model_type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returned state vector</td>
</tr>
<tr>
<td valign="top"><em class=
"code">assim_model   </em></td>
<td>Input assim_model_type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="copy_assim_model" id="copy_assim_model"></a><br>
<div class="routine"><em class="call">call
copy_assim_model(model_out, model_in)</em>
<pre>
type(assim_model_type), intent(out) :: <em class=
"code">model_out</em>
type(assim_model_type), intent(in)  :: <em class=
"code">model_in</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Copies one assim_model_type to another.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">model_out   </em></td>
<td>Copy.</td>
</tr>
<tr>
<td valign="top"><em class="code">model_in</em></td>
<td>Data to be copied.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="interpolate" id="interpolate"></a><br>
<div class="routine"><em class="call">call interpolate(x, location,
loctype, obs_vals, istatus)</em>
<pre>
real(r8),            intent(in)  :: <em class="code">x(:)</em>
type(location_type), intent(in)  :: <em class="code">location</em>
integer,             intent(in)  :: <em class="code">loctype</em>
real(r8),            intent(out) :: <em class="code">obs_vals</em>
integer,             intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Interpolates a given model state variable type to a location
given the model state vector. Nearly direct call to
model_interpolate in model_mod. See model_mod for the error return
values in istatus.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>Model state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">location   </em></td>
<td>Location to which to interpolate.</td>
</tr>
<tr>
<td valign="top"><em class="code">loctype</em></td>
<td>Type of variable to interpolate.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_vals</em></td>
<td>Returned interpolated value.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Returned as 0 if all is well, else various errors.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_model_time" id="set_model_time"></a><br>
<div class="routine"><em class="call">call
set_model_time(assim_model, time)</em>
<pre>
type(assim_model_type), intent(inout) :: <em class=
"code">assim_model</em>
type(time_type), intent(in)           :: <em class="code">time</em>
</pre></div>
<!-- Description -->
<div class="indent1">
<p>Sets the time in an assim_model_type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">assim_model   </em></td>
<td>Set the time in this assim_model_type.</td>
</tr>
<tr>
<td valign="top"><em class="code">time</em></td>
<td>Set to this time</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_model_state_vector" id=
"set_model_state_vector"></a><br>
<div class="routine"><em class="call">call
set_model_state_vector(assim_model, state)</em>
<pre>
type(assim_model_type), intent(inout) :: <em class=
"code">assim_model</em>
real(r8), intent(in)                  :: <em class=
"code">state(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the state in an assim_model_type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">assim_model   </em></td>
<td>Set the state vector in this assim_model_type.</td>
</tr>
<tr>
<td valign="top"><em class="code">state</em></td>
<td>The state vector to be inserted.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_state_restart" id="write_state_restart"></a><br>
<div class="routine"><em class="call">call
write_state_restart(assim_model, funit <em class=
"optionalcode">[, target_time]</em>)</em>
<pre>
type(assim_model_type),    intent(in) :: <em class=
"code">assim_model</em>
integer,                   intent(in) :: <em class=
"code">funit</em>
type(time_type), optional, intent(in) :: <em class=
"optionalcode">target_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes a restart from an assim_model_type with an optional
target_time.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">assim_model   </em></td>
<td>Write a restart from this assim_model_type.</td>
</tr>
<tr>
<td valign="top"><em class="code">funit</em></td>
<td>Integer file unit id open for output of restart files.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">target_time</em></td>
<td>If present, put this target time at the front of the restart
file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_state_restart" id="read_state_restart"></a><br>
<div class="routine"><em class="call">call
read_state_restart(assim_model, funit <em class=
"optionalcode">[, target_time]</em>)</em>
<pre>
type(assim_model_type),    intent(out) :: <em class=
"code">assim_model</em>
integer,                   intent(in)  :: <em class=
"code">funit</em>
type(time_type), optional, intent(out) :: <em class=
"optionalcode">target_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Read a state restart file into assim_model_type. Optionally read
a prepended target time.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">assim_model   </em></td>
<td>Read the time and state vector from restart into this.</td>
</tr>
<tr>
<td valign="top"><em class="code">funit</em></td>
<td>File id that has been opened for reading restart files.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">target_time</em></td>
<td>If present, read a target time from the front of the file into
this.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="output_diagnostics" id="output_diagnostics"></a><br>
<div class="routine"><em class="call">call
output_diagnostics(ndFileID, state <em class=
"optionalcode">[, copy_index]</em>)</em>
<pre>
type(netcdf_file_type), intent(inout) :: <em class=
"code">ndFileID</em>
type(assim_model_type), intent(in)    :: <em class=
"code">state</em>
integer, optional,      intent(in)    :: <em class=
"optionalcode">copy_index</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes one copy of the state time and vector to a netCDF
file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ndFileID</em></td>
<td>An identifier for a netCDF file</td>
</tr>
<tr>
<td valign="top"><em class="code">state</em></td>
<td>State vector and time</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">copy_index   </em></td>
<td>Which copy of state is to be output</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_assim_model" id="end_assim_model"></a><br>
<div class="routine"><em class="call">call
end_assim_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Called to clean-up at end of assim_model use. For now just
passes through to model_mod.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="input_diagnostics" id="input_diagnostics"></a><br>
<div class="routine"><em class="call">call
input_diagnostics(file_id, state, copy_index)</em>
<pre>
integer,                intent(in)    :: <em class=
"code">file_id</em>
type(assim_model_type), intent(inout) :: <em class=
"code">state</em>
integer,                intent(out)   :: <em class=
"code">copy_index</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Used to read in a particular copy of the state vector from an
open state diagnostics file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">file_id</em></td>
<td>Integer descriptor (channel number) for a diagnostics file
being read.</td>
</tr>
<tr>
<td valign="top"><em class="code">state</em></td>
<td>Assim_model_type to read in data.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">copy_index   </em></td>
<td>Which copy of state to be read.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_diag_input" id="init_diag_input"></a><br>
<div class="routine"><em class="call">var =
init_diag_input(file_name, global_meta_data, model_size,
copies_of_field_per_time)</em>
<pre>
integer                       :: <em class=
"code">init_diag_input</em>
character(len=*), intent(in)  :: <em class="code">file_name</em>
character(len=*), intent(out) :: <em class=
"code">global_meta_data</em>
integer,          intent(out) :: <em class="code">model_size</em>
integer,          intent(out) :: <em class=
"code">copies_of_field_per_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Opens a state diagnostic file and reads the global meta data,
model size, and number of data copies.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns the unit number on which the file is open.</td>
</tr>
<tr>
<td valign="top"><em class="code">file_name</em></td>
<td>File name of state diagnostic file.</td>
</tr>
<tr>
<td valign="top"><em class="code">global_meta_data</em></td>
<td>Global metadata string from file.</td>
</tr>
<tr>
<td valign="top"><em class="code">model_size</em></td>
<td>Size of model.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">copies_of_field_per_time   </em></td>
<td>Number of copies of the state vector at each time.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_assim_model" id="init_assim_model"></a><br>
<div class="routine"><em class="call">call
init_assim_model(state)</em>
<pre>
type(assim_model_type), intent(inout) :: <em class=
"code">state</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Creates storage for an assim_model_type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state   </em></td>
<td>An assim_model_type that needs storage created.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_diag_input_copy_meta_data" id=
"get_diag_input_copy_meta_data"></a><br>
<div class="routine"><em class="call">call
get_diag_input_copy_meta_data(file_id, model_size_out, num_copies,
location, meta_data_per_copy)</em>
<pre>
integer,             intent(in)  :: <em class="code">file_id</em>
integer,             intent(in)  :: <em class=
"code">model_size_out</em>
integer,             intent(in)  :: <em class=
"code">num_copies</em>
type(location_type), intent(out) :: <em class=
"code">location(model_size_out)</em>
character(len = *)               :: <em class=
"code">meta_data_per_copy(num_copies)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads meta-data describing state vectors in a state diagnostics
file. Given the file, the model_size, and the number of copies,
returns the locations of each state variable and the text
description of each copy.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">file_id</em></td>
<td>Integer channel open to state diagostic file being read</td>
</tr>
<tr>
<td valign="top"><em class="code">Model_size_out</em></td>
<td>model size</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies</em></td>
<td>Number of copies of state in file</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Returned locations for state vector</td>
</tr>
<tr>
<td valign="top"><em class=
"code">meta_data_per_copy   </em></td>
<td>Meta data describing what is in each copy of state vector</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="finalize_diag_output" id="finalize_diag_output"></a><br>
<div class="routine"><em class="call">var =
finalize_diag_output(ncFileID)</em>
<pre>
integer                               :: <em class=
"code">finalize_diag_output</em>
type(netcdf_file_type), intent(inout) :: <em class=
"code">ncFileID</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Used to complete writing on and open netcdf file. An error
return is provided for passing to the netcdf error handling
routines.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns an error value.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ncFileID   </em></td>
<td>Netcdf file id of an open file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="aread_state_restart" id="aread_state_restart"></a><br>
<div class="routine"><em class="call">call
aread_state_restart(model_time, model_state, funit <em class=
"optionalcode">[, target_time]</em>)</em>
<pre>
type(time_type),           intent(out) :: <em class=
"code">model_time</em>
real(r8),                  intent(out) :: <em class=
"code">model_state(:)</em>
integer,                   intent(in)  :: <em class=
"code">funit</em>
type(time_type), optional, intent(out) :: <em class=
"optionalcode">target_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads a model time and state, and optionally a prepended target
time, from a state restart file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">model_time</em></td>
<td>Returned time of model state</td>
</tr>
<tr>
<td valign="top"><em class=
"code">model_state   </em></td>
<td>Returned model state.</td>
</tr>
<tr>
<td valign="top"><em class="code">funit</em></td>
<td>Channel open for reading a state restart file.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">target_time</em></td>
<td>If present, this time is read from the front of the restart
file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="aoutput_diagnostics" id="aoutput_diagnostics"></a><br>
<div class="routine"><em class="call">call
aoutput_diagnostics(ncFileID, model_time, model_state <em class=
"optionalcode">[, copy_index]</em>)</em>
<pre>
type(netcdf_file_type), intent(inout) :: <em class=
"code">ncFileID</em>
type(time_type),        intent(in)    :: <em class=
"code">model_time</em>
real(r8),               intent(in)    :: <em class=
"code">model_state(:)</em>
integer, optional,      intent(in)    :: <em class=
"optionalcode">copy_index</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Write a state vector to a state diagnostics netcdf file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>Unit for a state vector netcdf file open for output.</td>
</tr>
<tr>
<td valign="top"><em class="code">model_time</em></td>
<td>The time of the state to be output</td>
</tr>
<tr>
<td valign="top"><em class=
"code">model_state   </em></td>
<td>A model state vector to be output.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">copy_index</em></td>
<td>Which copy of state vector is to be written, default is copy
1</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="awrite_state_restart" id="awrite_state_restart"></a><br>
<div class="routine"><em class="call">call
awrite_state_restart(model_time, model_state, funit <em class=
"optionalcode">[, target_time]</em>)</em>
<pre>
type(time_type),           intent(in) :: <em class=
"code">model_time</em>
real(r8),                  intent(in) :: <em class=
"code">model_state(:)</em>
integer,                   intent(in) :: <em class=
"code">funit</em>
type(time_type), optional, intent(in) :: <em class=
"optionalcode">target_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes a model time and state vector to a restart file and
optionally prepends a target time.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">model_time</em></td>
<td>Time of model state.</td>
</tr>
<tr>
<td valign="top"><em class="code">model_state</em></td>
<td>Model state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">funit</em></td>
<td>Channel of file open for restart output.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">target_time   </em></td>
<td>If present, time to be prepended to state time / vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="pert_model_state" id="pert_model_state"></a><br>
<div class="routine"><em class="call">call
pert_model_state()</em></div>
<div class="indent1"><!-- Description -->
<p>Passes through to pert_model_state in model_mod. See model_mod
documentation for arguments and details.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_append_time" id="nc_append_time"></a><br>
<div class="routine"><em class="call">var =
nc_append_time(ncFileID, time)</em>
<pre>
integer                               :: <em class=
"code">nc_append_time</em>
type(netcdf_file_type), intent(inout) :: <em class=
"code">ncFileID</em>
type(time_type),        intent(in)    :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Appends the time to the time coordinate variable of the netcdf
file. The new length of the time variable is returned. Requires
that time is a coordinate variable AND it is the unlimited
dimension.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns new length of time variable.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ncFileID   </em></td>
<td>Points to open netcdf file.</td>
</tr>
<tr>
<td valign="top"><em class="code">time</em></td>
<td>The next time to be added to the file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_calendar_atts" id=
"nc_write_calendar_atts"></a><br>
<div class="routine"><em class="call">var =
nc_write_calendar_atts(ncFileID, TimeVarID)</em>
<pre>
integer                            :: <em class=
"code">nc_write_calendar_atts</em>
type(netcdf_file_type), intent(in) :: <em class=
"code">ncFileID</em>
integer,                intent(in) :: <em class=
"code">TimeVarID</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets up the metadata for the appropriate calendar being used in
the time manager an writes it to a netcdf file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns a netcdf error code.</td>
</tr>
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>Netcdf file id pointing to a file open for writing.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">TimeVarID   </em></td>
<td>The index of the time variable in the netcdf file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_get_tindex" id="nc_get_tindex"></a><br>
<div class="routine"><em class="call">var = nc_get_tindex(ncFileID,
statetime)</em>
<pre>
integer                               :: <em class=
"code">nc_get_tindex</em>
type(netcdf_file_type), intent(inout) :: <em class=
"code">ncFileID</em>
type(time_type),        intent(in)    :: <em class=
"code">statetime</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the index of a time from the time variable in a netcdf
file. This function has been replaced with more efficient
approaches and may be deleted from future releases.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>The index of the time in the netcdf file.</td>
</tr>
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>File id for an open netcdf file.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">statetime   </em></td>
<td>The time to be found in the netcdf file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_model_time_step" id="get_model_time_step"></a><br>
<div class="routine"><em class="call">var =
get_model_time_step()</em>
<pre>
type(time_type) :: <em class="code">get_model_time_step</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This passes through to model_mod. See model_mod documentation
for arguments and details.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var   </em></td>
<td>Returns time step of model.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="open_restart_read" id="open_restart_read"></a><br>
<div class="routine"><em class="call">var =
open_restart_read(file_name)</em>
<pre>
integer                      :: <em class=
"code">open_restart_read</em>
character(len=*), intent(in) :: <em class="code">file_name</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Opens a restart file for readig.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns a file descriptor (channel number).</td>
</tr>
<tr>
<td valign="top"><em class=
"code">file_name   </em></td>
<td>Name of restart file to be open for reading.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="open_restart_write" id="open_restart_write"></a><br>
<div class="routine"><em class="call">var =
open_restart_write(file_name)</em>
<pre>
integer                      :: <em class=
"code">open_restart_write</em>
character(len=*), intent(in) :: <em class="code">file_name</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Open a restart file for writing.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Returns a file descriptor (channel) for a restart file.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">file_name   </em></td>
<td>File name of restart file to be opened.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="close_restart" id="close_restart"></a><br>
<div class="routine"><em class="call">call
close_restart(file_unit)</em>
<pre>
integer, intent(in) :: <em class="code">file_unit</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Closes a restart file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">file_unit   </em></td>
<td>File descriptor (channel number) of open restart file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="adv_1step" id="adv_1step"></a><br>
<div class="routine"><em class="call">call adv_1step()</em></div>
<div class="indent1"><!-- Description -->
<p>Advances a model by one step. Pass through to model_mod. See
model_mod documentation for arguments and details.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_initial_conditions" id=
"get_initial_conditions"></a><br>
<div class="routine"><em class="call">call
get_initial_condition(time, x)</em>
<pre>
type(time_type), intent(out) :: <em class="code">time</em>
real(r8),        intent(out) :: <em class="code">x</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Obtains an initial condition from models that support this
option.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">time   </em></td>
<td>the valid time of the model state</td>
</tr>
<tr>
<td valign="top"><em class="code">x</em></td>
<td>the initial model state</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="ens_mean_for_model" id="ens_mean_for_model"></a><br>
<div class="routine"><em class="call">call
ens_mean_for_model(ens_mean)</em>
<pre>
type(r8), intent(in) :: <em class="code">ens_mean(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>An array of length model_size containing the ensemble means.
This is a direct pass through to the model_mod.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">ens_mean   </em></td>
<td>Array of length model_size containing the mean for each entry
in the state vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_maxdist_init" id=
"get_close_maxdist_init"></a><br>
<div class="routine"><em class="call">call
get_close_maxdist_init(gc, maxdist)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
type(r8), intent(in)                :: <em class=
"code">maxdist</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the threshold distance. Anything closer than this is deemed
to be close. This is a direct pass through to the model_mod, which
in turn can pass through to the location_mod.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Data for efficiently finding close locations.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">maxdist   </em></td>
<td>Anything closer than this distance is a close location.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs" id="get_close_obs"></a><br>
<div class="routine"><em class="call">call get_close_obs(gc,
base_obs_loc, base_obs_kind, obs, obs_kind, num_close, close_ind
<em class="optionalcode">[, dist]</em>)</em>
<pre>
type(get_close_type), intent(in)  :: <em class="code">gc</em>
type(location_type),  intent(in)  :: <em class=
"code">base_obs_loc</em>
integer,              intent(in)  :: <em class=
"code">base_obs_kind</em>
type(location_type),  intent(in)  :: <em class="code">obs(:)</em>
integer,              intent(in)  :: <em class=
"code">obs_kind(:)</em>
integer,              intent(out) :: <em class=
"code">num_close</em>
integer,              intent(out) :: <em class=
"code">close_ind(:)</em>
real(r8),  optional,  intent(out) :: <em class=
"optionalcode">dist(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a single location and a list of other locations, returns
the indices of all the locations close to the single one along with
the number of these and the distances for the close ones. The
observation kinds are passed in to allow more sophisticated
distance computations to be done if needed. This is a direct pass
through to the model_mod, which in turn can pass through to the
location_mod.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Data for efficiently finding close locations.</td>
</tr>
<tr>
<td valign="top"><em class="code">base_obs_loc</em></td>
<td>Single given location.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">base_obs_kind   </em></td>
<td>Kind of the single location.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs</em></td>
<td>List of observations from which close ones are to be
found.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_kind</em></td>
<td>Kind associated with observations in obs list.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_close</em></td>
<td>Number of observations close to the given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">close_ind</em></td>
<td>Indices of those locations that are close.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">dist</em></td>
<td>Distance between given location and the close ones identified
in close_ind.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs_init" id="get_close_obs_init"></a><br>
<div class="routine"><em class="call">call get_close_obs_init(gc,
num, obs)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
integer,              intent(in)    :: <em class="code">num</em>
type(location_type),  intent(in)    :: <em class="code">obs(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initialize storage for efficient identification of locations
close to a given location. Allocates storage for keeping track of
which 'box' each observation in the list is in. This is a direct
pass through to the model_mod, which in turn can pass through to
the location_mod.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Data for efficiently finding close locations.</td>
</tr>
<tr>
<td valign="top"><em class="code">num</em></td>
<td>The number of locations in the list.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs   </em></td>
<td>The location of each element in the list, not used in 1D
implementation.</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 
<!--kdr; Should these files be listed as "used by this module"? -->
<!--     What about diagnostics files? -->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<table border="0">
<tr>
<th>filename   </th>
<th>purpose/comment</th>
</tr>
<tr>
<td>filter_restart</td>
<td>specified in &amp;filter_nml:restart_in_filename</td>
</tr>
<tr>
<td>filter_restart</td>
<td>specified in &amp;filter_nml:restart_out_filename</td>
</tr>
<tr>
<td>input.nml</td>
<td>to read namelists</td>
</tr>
</table>
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
<td valign="top">init_diag_output</td>
<!-- message -->
<td valign="top">Compiler does not support required kinds of
variables.</td>
<!-- comment -->
<td valign="top">NetCDF-f90 interface function byteSizeOK returned
FALSE</td>
</tr>
<tr><!-- routine -->
<td valign="top">init_diag_output and various nc_XXX</td>
<!-- message -->
<td valign="top">various NetCDF-f90 messages</td>
<!-- comment -->
<td valign="top">Returned by one of the NetCDF calls in this
subroutine. Consult the NetCDF manual.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_diag_input_copy_meta_data</td>
<!-- message -->
<td valign="top">expected to read "locat" got ...</td>
<!-- comment -->
<td valign="top">The header of the metadata for the copies of the
data in diagnostic input file is not = 'locat'</td>
</tr>
<tr><!-- routine -->
<td valign="top">set_model_state_vector</td>
<!-- message -->
<td valign="top">state vector has length # model size (#) does not
match.</td>
<!-- comment -->
<td valign="top">Check your model resolution and fields included in
the state vector.</td>
</tr>
<tr><!-- routine -->
<td valign="top">aread_state_restart</td>
<!-- message -->
<td valign="top">read error is : #</td>
<!-- comment -->
<td valign="top">Unable to read model state from
assim_model_state_ic# file. # is error condition retured by read
statement.</td>
</tr>
<tr><!-- routine -->
<td valign="top">open_restart_read</td>
<!-- message -->
<td valign="top">OPEN status was #</td>
<!-- comment -->
<td valign="top">Failed to open file listed for reason #.</td>
</tr>
<tr><!-- routine -->
<td valign="top">aoutput_diagnostics</td>
<!-- message -->
<td valign="top">model time (d,s) (#,#) is index # in ncFileID
#</td>
<!-- comment -->
<td valign="top">Time index for file listed is &lt; 0</td>
</tr>
<tr><!-- routine -->
<td valign="top">ainput_diagnostics</td>
<!-- message -->
<td valign="top">expected "copy", got _____'</td>
<!-- comment -->
<td valign="top">Trying to read diagnostic state output
header.</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_append_time</td>
<!-- message -->
<td valign="top">"time" expected to be rank-1</td>
<!-- comment -->
<td valign="top">ndims /= 1<br>
The time array of the NetCDF file should be 1-dimensional</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_append_time</td>
<!-- message -->
<td valign="top">unlimited dimension expected to be
slowest-moving</td>
<!-- comment -->
<td valign="top">dimids(1) /= unlimitedDimID</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_append_time</td>
<!-- message -->
<td valign="top">time mirror and netcdf file time dimension
out-of-sync</td>
<!-- comment -->
<td valign="top">lngth /= ncFileId%Ntimes</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_append_time</td>
<!-- message -->
<td valign="top">various NetCDF-f90 error messages</td>
<!-- comment -->
<td valign="top">Returned from one of the NetCDF calls in this
subroutine. Consult the NetCDF manual.</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_get_tindex</td>
<!-- message -->
<td valign="top">trouble deep ... can go no farther. Stopping.</td>
<!-- comment -->
<td valign="top">timeindex &lt; -1</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_get_tindex</td>
<!-- message -->
<td valign="top">Model time preceeds earliest netCDF time.</td>
<!-- comment -->
<td valign="top">Time of current assim_model is earlier than all
the times on the NetCDF file to which the state is to be written by
aoutput_diagnostics.</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_get_tindex</td>
<!-- message -->
<td valign="top">subsequent netCDF time (days, seconds) # #</td>
<!-- comment -->
<td valign="top">Time of current assim_model is in the midst of the
times on the NetCDF file to which the state is to be written by
aoutput_diagnostics, but doesn't match any of them. Very bad.</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_get_tindex</td>
<!-- message -->
<td valign="top">various NetCDF-f90 error messages</td>
<!-- comment -->
<td valign="top">Returned from one of the NetCDF calls in this
subroutine. Consult the NetCDF manual.</td>
</tr>
<tr><!-- routine -->
<td valign="top">nc_write_calendar_atts</td>
<!-- message -->
<td valign="top">various NetCDF-f90 error messages</td>
<!-- comment -->
<td valign="top">Returned from one of the NetCDF calls in this
subroutine. Consult the NetCDF manual.</td>
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
