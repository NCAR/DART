<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module model_mod (pe2lyr)</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>pe2lyr</h1>
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
<p>DART standard interfaces for a two-layer isentropic primitive
equation model.<br>
<br>
The 16 public interfaces are standardized for all DART compliant
models. These interfaces allow DART to advance the model, get the
model state and metadata describing this state, find state
variables that are close to a given location, and do spatial
interpolation for model state variables.<br>
<br>
This model is a 2-layer, isentropic, primitive equation model on a
sphere. TODO: add more detail here, including equations, etc.<br>
<br>
Contact: Jeffrey.S.Whitaker@noaa.gov</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
time_manager_mod
utilities_mod
random_seq_mod
threed_sphere/location_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
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
<p>Returns the size of the model as an integer. For this model the
default grid size is 96 (lon) by 48 (lat) by 2 levels, and 3
variables (U, V, Z) at each grid location, for a total size of
27,648. There are alternative include files which, if included at
compile time instead of the default file, defines a grid at twice
and 4 times this resolution. They have corresponding truncation
values of T63 and T127 (the default grid uses T31).</p>
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
<p>Returns metadata about a given element, indexed by <em class=
"code">index_in</em>, in the model state vector. The <em class=
"code">location</em> defines where the state variable is
located.<br>
<br>
For this model, the default grid is a global lat/lon grid, 96 (lon)
by 48 (lat) by 2 levels. The variable types are U, V, and Z:</p>
<ul>
<li>1 = TYPE_u</li>
<li>2 = TYPE_v</li>
<li>901 = TYPE_z</li>
</ul>
<p>Grids at twice and 4 times the resolution can be compiled in
instead by using one of the alternative header files (see
<em class="code">resolt31.h</em> (the default), <em class=
"code">resolt63.h</em>, and <em class="code">resolt127.h</em>).</p>
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
<td>The type of the state variable element.</td>
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
<p>Given a state vector, a location, and a model state variable
type, interpolates the state variable field to that location and
returns the value in obs_val. The istatus variable is always
returned as 0 (OK).</p>
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
<td>Type of state field to be interpolated.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_val</em></td>
<td>The interpolated value from the model.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Integer value returning 0 for successful, other values can be
defined for various failures.</td>
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
<p>Returns the the time step of the model; the smallest increment
in time that the model is capable of advancing the state in a given
implementation. For this model the default value is 20 minutes
(1200 seconds), but also comes with header files with times steps
of 10 and 5 minutes (for higher grid resolution and truncation
constants).</p>
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
<p>Used for runtime initialization of a model, for instance
calculating storage requirements, initializing model parameters,
etc. This is the first call made to a model by any DART compliant
assimilation routines.<br>
<br>
In this model, it allocates space for the grid, and initializes the
grid locations, data values, and various parameters, including
spherical harmonic weights.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_model" id="end_model"></a><br>
<div class="routine"><em class="call">call end_model()</em></div>
<div class="indent1"><!-- Description -->
<p>A stub since the pe2lyr model does no cleanup.</p>
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
initial conditions are to be used. This model sets the time to
0.</p>
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
<p>Returns default initial conditions for model; generally used for
spinning up initial model states. This model sets the default state
vector based on the initialized fields in the model. (TODO: which
are what?)</p>
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
<p>This routine writes the model-specific attributes to a netCDF
file. This includes coordinate variables and any metadata, but NOT
the model state vector. This model writes out the data as U, V, and
Z arrays on a lat/lon/height grid, so the attributes are organized
in the same way.</p>
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
<p>This routine writes the model-specific state vector (data) to a
netCDF file. This model writes out the data as U, V, and Z arrays
on a lat/lon/height grid.</p>
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
<p>Given a model state vector, perturbs this vector. Used to
generate initial conditions for spinning up ensembles. This model
has no code to generate these values, so it returns <em class=
"code">interf_provided</em> as .false. and the default algorithms
in filter are then used by the calling code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state</em></td>
<td>State vector to be perturbed.</td>
</tr>
<tr>
<td valign="top"><em class="code">pert_state</em></td>
<td>Perturbed state vector</td>
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
<p>In distance computations any two locations closer than the given
<em class="code">maxdist</em> will be considered close by the
<em class="code">get_close_obs()</em> routine. Pass-through to the
3-D sphere locations module. See <a href=
"../../location/threed_sphere/location_mod.html#get_close_maxdist_init">
get_close_maxdist_init()</a> for the documentation of this
subroutine.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc  </em></td>
<td>The get_close_type which stores precomputed information about
the locations to speed up searching</td>
</tr>
<tr>
<td valign="top"><em class="code">maxdist  </em></td>
<td>Anything closer than this will be considered close.</td>
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
type(location_type),  intent(in)    :: <em class=
"code">obs(num)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3-D sphere locations module. See <a href=
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
<p>Given a location and kind, compute the distances to all other
locations in the <em class="code">obs</em> list. The return values
are the number of items which are within maxdist of the base, the
index numbers in the original obs list, and optionally the
distances. The <em class="code">gc</em> contains precomputed
information to speed the computations.<br>
<br>
Pass-through to the 3-D sphere locations module. See <a href=
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
<p>Stub only. Not needed by this model.</p>
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
<div class="top">[<a href="#">top</a>]</div>
<hr>
<p>This model currently has no values settable by namelist.</p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>The model source is in pe2lyr_mod.f90, and the spherical
harmonic code is in spharmt_mod.f90. The various resolution
settings are in resolt31.h, resolt63.h, and resolt127.h.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<p>Zou, X., Barcilon, A., Navon, I.M., Whitaker, J., Cacuci, D.G..
1993: An Adjoint Sensitivity Study of Blocking in a Two-Layer
Isentropic Model. Monthly Weather Review: Vol. 121, No. 10, pp.
2833-2857.</p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<p>N/A</p>
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
