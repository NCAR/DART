<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module model_mod (SQG)</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>SQG</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>This is a uniform PV two-surface QG+1 spectral model contributed
by Rahul Majahan.</p>
<p>The underlying model is described in: Hakim, Gregory J., 2000:
Role of Nonmodal Growth and Nonlinearity in Cyclogenesis
Initial-Value Problems. J. Atmos. Sci., 57, 2951-2967. doi:
10.1175/1520-0469(2000)057&lt;2951:RONGAN&gt;2.0.CO;2</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
time_manager_mod
threed_sphere/location_mod
utilities_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table>
<tr>
<td><em class="call">use model_mod, only :</em></td>
<td><a href="#get_model_size">get_model_size</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#adv_1step">adv_1step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_state_meta_data">get_state_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#model_interpolate">model_interpolate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_model_time_step">get_model_time_step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#static_init_model">static_init_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#end_model">end_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_time">init_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_conditions">init_conditions</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_write_model_atts">nc_write_model_atts</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_write_model_vars">nc_write_model_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#pert_model_state">pert_model_state</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_close_maxdist_init">get_close_maxdist_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs_init">get_close_obs_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs">get_close_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ens_mean_for_model">ens_mean_for_model</a></td>
</tr>
</table>
<p>Optional namelist interface <a href="#Namelist"><em class=
"code">&amp;model_nml</em></a> may be read from file <em class=
"file">input.nml</em>.</p>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_model_size" id="get_model_size"></a><br>
<div class="routine"><em class="call">model_size = get_model_size(
)</em>
<pre>
integer :: <em class="code">get_model_size</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the length of the model state vector.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">model_size</em></td>
<td>The length of the model state vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="adv_1step" id="adv_1step"></a><br>
<div class="routine"><em class="call">call adv_1step(x, time)</em>
<pre>
real(r8), dimension(:), intent(inout) :: <em class="code">x</em>
type(time_type),        intent(in)    :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Advances the model for a single time step. The time associated
with the initial model state is also input although it is not used
for the computation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>State vector of length model_size.</td>
</tr>
<tr>
<td valign="top"><em class="code">time   </em></td>
<td>Specifies time of the initial model state.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_state_meta_data" id="get_state_meta_data"></a><br>
<div class="routine"><em class="call">call get_state_meta_data
(index_in, location, <em class=
"optionalcode">[, var_type]</em> )</em>
<pre>
integer,             intent(in)  :: <em class="code">index_in</em>
type(location_type), intent(out) :: <em class="code">location</em>
integer, optional,   intent(out) :: <em class=
"optionalcode"> var_type </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns metadata about a given element, indexed by index_in, in
the model state vector. The location defines where the state
variable is located.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">index_in   </em></td>
<td>Index of state vector element about which information is
requested.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>The location of state variable element.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">var_type</em></td>
<td>Returns the type (always 1) of the indexed state variable as an
optional argument.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="model_interpolate" id="model_interpolate"></a><br>
<div class="routine"><em class="call">call model_interpolate(x,
location, itype, obs_val, istatus)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">x</em>
type(location_type),    intent(in)  :: <em class=
"code">location</em>
integer,                intent(in)  :: <em class="code">itype</em>
real(r8),               intent(out) :: <em class=
"code">obs_val</em>
integer,                intent(out) :: <em class=
"code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given model state, returns the value interpolated to a given
location.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">location   </em></td>
<td>Location to which to interpolate.</td>
</tr>
<tr>
<td valign="top"><em class="code">itype</em></td>
<td>Not used.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_val</em></td>
<td>The interpolated value from the model.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Quality control information, always returned 0.</td>
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
<p>Returns the time step (forecast length) of the model;</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var   </em></td>
<td>Smallest time step of model.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="static_init_model" id="static_init_model"></a><br>
<div class="routine"><em class="call">call
static_init_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Used for runtime initialization of model; reads namelist,
initializes model parameters, etc. This is the first call made to
the model by any DART-compliant assimilation routine.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_model" id="end_model"></a><br>
<div class="routine"><em class="call">call end_model()</em></div>
<div class="indent1"><!-- Description -->
<p>A stub.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_time" id="init_time"></a><br>
<div class="routine"><em class="call">call init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the time at which the model will start if no input
initial conditions are to be used. This is used to spin-up the
model from rest.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">time   </em></td>
<td>Initial model time.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_conditions" id="init_conditions"></a><br>
<div class="routine"><em class="call">call init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class="code">x</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns default initial conditions for the model; generally used
for spinning up initial model states.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x   </em></td>
<td>Initial conditions for state vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_atts" id="nc_write_model_atts"></a><br>
<div class="routine"><em class="call">ierr =
nc_write_model_atts(ncFileID)</em>
<pre>
integer             :: <em class="code">nc_write_model_atts</em>
integer, intent(in) :: <em class="code">ncFileID</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Function to write model specific attributes to a netCDF file. At
present, DART is using the NetCDF format to output diagnostic
information. This is not a requirement, and models could choose to
provide output in other formats. This function writes the metadata
associated with the model to a NetCDF file opened to a file
identified by ncFileID.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">ncFileID   </em></td>
<td>Integer file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returns a 0 for successful completion.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_vars" id="nc_write_model_vars"></a><br>
<div class="routine"><em class="call">ierr =
nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)</em>
<pre>
integer                            :: <em class=
"code">nc_write_model_vars</em>
integer,                intent(in) :: <em class=
"code">ncFileID</em>
real(r8), dimension(:), intent(in) :: <em class=
"code">statevec</em>
integer,                intent(in) :: <em class=
"code">copyindex</em>
integer,                intent(in) :: <em class=
"code">timeindex</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes a copy of the state variables to a netCDF file. Multiple
copies of the state for a given time are supported, allowing, for
instance, a single file to include multiple ensemble estimates of
the state.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">statevec</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">copyindex   </em></td>
<td>Integer index of copy to be written.</td>
</tr>
<tr>
<td valign="top"><em class="code">timeindex</em></td>
<td>The timestep counter for the given state.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returns 0 for normal completion.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="pert_model_state" id="pert_model_state"></a><br>
<div class="routine"><em class="call">call pert_model_state(state,
pert_state, interf_provided)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">state</em>
real(r8), dimension(:), intent(out) :: <em class=
"code">pert_state</em>
logical,                intent(out) :: <em class=
"code">interf_provided</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a model state, produces a perturbed model state.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state</em></td>
<td>State vector to be perturbed.</td>
</tr>
<tr>
<td valign="top"><em class="code">pert_state</em></td>
<td>Perturbed state vector: NOT returned.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">interf_provided   </em></td>
<td>Returned false; interface is not implemented.</td>
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
real(r8),             intent(in)    :: <em class=
"code">maxdist</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3D Sphere locations module. See <a href=
"../../location/threed_sphere/location_mod.html#get_close_maxdist_init">
get_close_maxdist_init()</a> for the documentation of this
subroutine.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs_init" id="get_close_obs_init"></a><br>
<div class="routine"><em class="call">call get_close_obs_init(gc,
num, obs)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
integer,              intent(in)    :: <em class="code">num</em>
type(location_type),  intent(in)    :: <em class=
"code">obs(num)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3D Sphere locations module. See <a href=
"../../location/threed_sphere/location_mod.html#get_close_obs_init">
get_close_obs_init()</a> for the documentation of this
subroutine.</p>
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
real(r8), optional,   intent(out) :: <em class=
"optionalcode">dist(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3D Sphere locations module. See <a href=
"../../location/threed_sphere/location_mod.html#get_close_obs">get_close_obs()</a>
for the documentation of this subroutine.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="ens_mean_for_model" id="ens_mean_for_model"></a><br>
<div class="routine"><em class="call">call
ens_mean_for_model(ens_mean)</em>
<pre>
real(r8), dimension(:), intent(in) :: <em class=
"code">ens_mean</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>A NULL INTERFACE in this model.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">ens_mean   </em></td>
<td>State vector containing the ensemble mean.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
 <a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>We adhere to the F90 standard of starting a namelist with an
ampersand '&amp;' and terminating with a slash '/' for all our
namelist input.</p>
<div class="namelist">
<pre>
&amp;model_nml 
  output_state_vector = .false.
  channel_center = 45.0
  channel_width = 40.0
  assimilation_period_days = 0
  assimilation_period_seconds = 21600
  debug = .false.
/
</pre></div>
<div class="indent1"><!-- Description -->
<p>This namelist is read in a file called <em class=
"file">input.nml</em></p>
<table border="0" cellpadding="3" width="100%">
<thead align="left">
<tr>
<th>Contents</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>output_state_vector</td>
<td>logical</td>
<td>If .true. write state vector as a 1D array to the diagnostic
output file. If .false. break state vector up into fields before
writing to the outputfile.</td>
</tr>
<tr>
<td>channel_center</td>
<td>real(r8)</td>
<td>Channel center</td>
</tr>
<tr>
<td>channel_width</td>
<td>real(r8)</td>
<td>Channel width</td>
</tr>
<tr>
<td>assimilation_period_days</td>
<td>integer</td>
<td>Number of days for timestep</td>
</tr>
<tr>
<td>assimilation_period_seconds</td>
<td>integer</td>
<td>Number of seconds for timestep</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>Set to .true. for more output</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<hr>
<h2>FILES</h2>
<table border="0">
<tr>
<th align="left">filename</th>
<th align="left">purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read the model_mod namelist</td>
</tr>
<tr>
<td>preassim.nc</td>
<td>the time-history of the model state before assimilation</td>
</tr>
<tr>
<td>analysis.nc </td>
<td>the time-history of the model state after assimilation</td>
</tr>
<tr>
<td>dart_log.out [default name]</td>
<td>the run-time diagnostic output</td>
</tr>
<tr>
<td>dart_log.nml [default name]</td>
<td>the record of all the namelists actually USED - contains the
default values</td>
</tr>
</table>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<hr>
<h2>REFERENCES</h2>
<p>The underlying model is described in:<br>
Hakim, Gregory J., 2000: Role of Nonmodal Growth and Nonlinearity
in Cyclogenesis Initial-Value Problems. J. Atmos. Sci., 57,
2951-2967. doi:
10.1175/1520-0469(2000)057&lt;2951:RONGAN&gt;2.0.CO;2</p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
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
<td valign="top">nc_write_model_atts<br>
nc_write_model_vars</td>
<!-- message -->
<td valign="top">Various netCDF-f90 interface error messages</td>
<!-- comment -->
<td valign="top">From one of the netCDF calls in the named
routine</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
