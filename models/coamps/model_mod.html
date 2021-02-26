<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module model_mod (COAMPS)</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>COAMPS</h1>
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
<p>DART interface module for the Coupled Ocean / Atmosphere
Mesoscale Prediction (COAMPS ®) model. The 16 public interfaces
listed here are standardized for all DART compliant models. These
interfaces allow DART to advance the model, get the model state and
metadata describing this state, find state variables that are close
to a given location, and do spatial interpolation for a variety of
variables required in observational operators.<br>
<br>
The following model description is taken from the <a href=
"http://www.nrlmry.navy.mil/coamps-web/web/view">COAMPS overview
web page:</a></p>
<blockquote>"The Coupled Ocean/Atmosphere Mesoscale Prediction
System (COAMPS) has been developed by the Marine Meteorology
Division (MMD) of the Naval Research Laboratory (NRL). The
atmospheric components of COAMPS, described below, are used
operationally by the U.S. Navy for short-term numerical weather
prediction for various regions around the world.<br>
<br>
The atmospheric portion of COAMPS represents a complete
three-dimensional data assimilation system comprised of data
quality control, analysis, initialization, and forecast model
components. Features include a globally relocatable grid,
user-defined grid resolutions and dimensions, nested grids, an
option for idealized or real-time simulations, and code that allows
for portability between mainframes and workstations. The
nonhydrostatic atmospheric model includes predictive equations for
the momentum, the non-dimensional pressure perturbation, the
potential temperature, the turbulent kinetic energy, and the mixing
ratios of water vapor, clouds, rain, ice, grauple, and snow, and
contains advanced parameterizations for boundary layer processes,
precipitation, and radiation.<br>
<br>
The distributed version of the COAMPS code that can be downloaded
from the web site has been designed to use the message-passing
interface (MPI), OpenMP directives, and horizontal domain
decomposition to achieve parallelism. The code is capable of
executing efficiently across vector, parallel, or symmetric
muti-processor (SMP) machines by simply changing run-time
options."</blockquote>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
time_manager_mod
threed_sphere/location_mod
utilities_mod
obs_kind_mod
random_seq_mod
netcdf
typesizes
coamps_grid_mod
coamps_interp_mod
coamps_restart_mod
coamps_util_mod
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
<tr>
<td> </td>
<td><a href="#adv_1step">adv_1step</a></td>
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
</table>
<p>The last 4 interfaces are only required for low-order models
where advancing the model can be done by a call to a subroutine.
The COAMPS model only advances by executing the coamps program.
Thus the last 4 interfaces only appear as stubs in this module.</p>
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
<p>Returns the length of the model state vector as an integer. This
includes all nested domains.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">model_size</em></td>
<td>The length of the model state vector.</td>
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
variable is located while the type of the variable (for instance
temperature, or u wind component) is returned by var_type. The
integer values used to indicate different variable types in
var_type are themselves defined as public interfaces to model_mod
if required.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">index_in   </em></td>
<td>Index of state vector element about which information is
requested.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Returns location of indexed state variable. The location should
use a location_mod that is appropriate for the model domain. For
realistic atmospheric models, for instance, a three-dimensional
spherical location module that can represent height in a variety of
ways is provided.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">var_type</em></td>
<td>Returns the type of the indexed state variable as an optional
argument.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="model_interpolate" id="model_interpolate"></a><br>
<div class="routine"><em class="call">call model_interpolate(x,
location, obs_kind, obs_val, istatus)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">x</em>
type(location_type),    intent(in)  :: <em class=
"code">location</em>
integer,                  intent(in)  :: <em class=
"code"> obs_kind </em>
real(r8),               intent(out) :: <em class=
"code">obs_val</em>
integer,                intent(out) :: <em class=
"code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given model state, returns the value of observation type
interpolated to a given location by a method of the model's
choosing. All observation kinds defined in obs_kind_mod are
supported. In the case where the observational operator is not
defined at the given location (e.g. the observation is below the
model surface or outside the domain), obs_val is returned as
-888888.0 and istatus = 1. Otherwise, istatus = 0. The
interpolation is performed in the domain with the highest
resolution containing the observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Location to which to interpolate.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_kind</em></td>
<td>Integer indexing which type of observation is to be
interpolated.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_val</em></td>
<td>The interpolated value from the model.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Integer flag indicating the result of the interpolation.</td>
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
<p>Returns the model base time step as a time_type. For now this is
set to 1 minute.</p>
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
<p>Used for runtime initialization of the model. This is the first
call made to the model by any DART compliant assimilation routine.
It reads the model namelist parameters, initializes the pressure
levels for the state vector, and generates the location data for
each member of the state.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_atts" id="nc_write_model_atts"></a><br>
<div class="routine"><em class="call">ierr =
nc_write_model_atts(ncFileId)</em>
<pre>
integer             :: <em class="code"> nc_write_model_atts</em>
integer, intent(in) :: <em class="code"> ncFileId </em>
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
<td valign="top"><em class="code">ncFileId   </em></td>
<td>Integer file descriptor opened to NetCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returned error code.</td>
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
"code"> nc_write_model_vars</em>
integer,                intent(in) :: <em class=
"code"> ncFileID </em>
real(r8), dimension(:), intent(in) :: <em class=
"code"> statevec </em>
integer,                intent(in) :: <em class=
"code"> copyindex</em>
integer,                intent(in) :: <em class=
"code"> timeindex </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes a copy of the state variables to a NetCDF file. Multiple
copies of the state for a given time are supported, allowing, for
instance, a single file to include multiple ensemble estimates of
the state.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID   </em></td>
<td>Integer file descriptor opened to NetCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">statevec</em></td>
<td>State vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">copyindex</em></td>
<td>Integer index to which copy is to be written.</td>
</tr>
<tr>
<td valign="top"><em class="code">timeindex</em></td>
<td>Integer index of which time in the file is being written.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returned error code.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="pert_model_state" id="pert_model_state"></a><br>
<div class="routine"><em class="call">call pert_model_state(state,
pert_state, interf_provided)</em>
<pre>
real(r8), dimension(:),   intent(in)    :: <em class=
"code"> state </em>
real(r8), dimension(:),   intent(out)   :: <em class=
"code"> pert_state </em>
logical,                  intent(out)   :: <em class=
"code"> interf_provided
</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a model state, produces a perturbed model state. This is
used to generate initial ensemble conditions perturbed around some
control trajectory state when one is preparing to spin-up
ensembles. In the COAMPS interface, this can be done three
different ways:</p>
<ul>
<li>No perturbation</li>
<li>Uniform perturbation - each element of the field has the same
additive perturbation</li>
<li>Individual perturbation - each element of the field has a
different additive perturbation The perturbation magnitude and
option are supplied out of the dynamic restart vector definition -
this allows us to supply a variance appropriate for each type of
variable at each level.</li>
</ul>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state</em></td>
<td>State vector to be perturbed.</td>
</tr>
<tr>
<td valign="top"><em class="code">pert_state</em></td>
<td>Perturbed state vector is returned.</td>
</tr>
<tr>
<td valign="top"><em class="code">interf_provided</em></td>
<td>Returns .true. for this model.</td>
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
<p>Pass-through to the 3-D sphere locations module. See <a href=
"../../location/threed-sphere/location_mod.html#get_close_maxdist_init">
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
<p>Pass-through to the 3-D sphere locations module. See <a href=
"../../location/threed-sphere/location_mod.html#get_close_obs_init">
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
<p>Pass-through to the 3-D sphere locations module. See <a href=
"../../location/threed-sphere/location_mod.html#get_close_obs">get_close_obs()</a>
for the documentation of this subroutine.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="ens_mean_for_model" id="ens_mean_for_model"></a><br>
<div class="routine"><em class="call">call
ens_mean_for_model(ens_mean)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class=
"code">ens_mean</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>A local copy is available here for use during other computations
in the model_mod code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_mean  </em></td>
<td>Ensemble mean state vector</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="adv_1step" id="adv_1step"></a><br>
<div class="routine"><em class="call">call adv_1step(x, time)</em>
<pre>
real(r8), dimension(:),   intent(inout) :: <em class=
"code"> x </em>
type(time_type),          intent(in)    :: <em class=
"code"> time </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This operation is not defined for the COAMPS model. This
interface is only required if `synchronous' model state advance is
supported (the model is called directly as a Fortran90 subroutine
from the assimilation programs). This is generally not the
preferred method for large models and a stub for this interface is
provided for the COAMPS model.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>State vector of length model_size.</td>
</tr>
<tr>
<td valign="top"><em class="code">time</em></td>
<td>Gives time of the initial model state. Needed for models that
have real time state requirements, for instance the computation of
radiational parameters. Note that DART provides a time_manager_mod
module that is used to support time computations throughout the
facility.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_model" id="end_model"></a><br>
<div class="routine"><em class="call">call end_model( )</em></div>
<div class="indent1"><!-- Description -->
<p>Called when use of a model is completed to clean up storage,
etc. A stub is provided for the COAMPS model.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_time" id="init_time"></a><br>
<div class="routine"><em class="call">call init_time(i_time)</em>
<pre>
type(time_type),        intent(in)  :: <em class=
"code"> i_time </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the time at which the model will start if no input
initial conditions are to be used. This is frequently used to
spin-up models from rest, but is not meaningfully supported for the
COAMPS model.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_conditions" id="init_conditions"></a><br>
<div class="routine"><em class="call">call init_conditions( x
)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class="code"> x </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns default initial conditions for model; generally used for
spinning up initial model states. For the COAMPS model just return
0's since initial state is always to be provided from input
files.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>Model state vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
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
&amp;model_nml
  cdtg = '2006072500',
  y_bound_skip = 3,
  x_bound_skip = 3,
  need_mean = .false.,
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
<td>cdtg</td>
<td>character(len=10)</td>
<td>Date/time group.</td>
</tr>
<tr>
<td>x_bound_skip, y_bound_skip</td>
<td>integer</td>
<td>Number of x and y boundary points to skip when perturbing the
model state.</td>
</tr>
<tr>
<td>need_mean</td>
<td>logical</td>
<td>Does the forward operator computation need the ensemble
mean?</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
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
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<p>The COAMPS registration web site is <a href=
"http://www.nrlmry.navy.mil/coamps-web/web/home">http://www.nrlmry.navy.mil/coamps-web/web/home</a>
and COAMPS is a registered trademark of the Naval Research
Laboratory.</p>
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
<td valign="top">nc_write_model_atts</td>
<!-- message -->
<td valign="top">Time dimension ID # must match Unlimited Dimension
ID #</td>
<!-- comment -->
<td valign="top">NetCDF file writing error</td>
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
