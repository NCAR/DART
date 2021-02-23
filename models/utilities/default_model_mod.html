<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module model_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE model_mod</h1>
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
<p>Every model that is DART compliant must provide an set of
interfaces that will be called by DART code. For models which have
no special code for some of these routines, they can pass through
the call to this default module, which satisfies the call but does
no work. To use these routines in a <em class=
"file">model_mod.f90</em>, add at the top:<br></p>
<pre>
use default_model_mod, only : xxx, yyy
</pre>
and then leave them in the public list. 
<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>The default routines have no namelist.</p>
<!--==================================================================-->
<p><!-- useless spacer to make the top behave correctly --></p>
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
time_manager_mod
location_mod
utilities_mod
netcdf_utilities_mod
ensemble_manager_mod
dart_time_io_mod
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
<td><a href=
"#shortest_time_between_assimilations">shortest_time_between_assimilations</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#static_init_model">static_init_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_time">init_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#fail_init_time">fail_init_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_conditions">init_conditions</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#fail_init_conditions">fail_init_conditions</a></td>
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
<td><a href="#pert_model_copies">pert_model_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs">get_close_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_state">get_close_state</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#convert_vertical_obs">convert_vertical_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#convert_vertical_state">convert_vertical_state</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_model_time">read_model_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_model_time">write_model_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#end_model">end_model</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_model_size" id="get_model_size"></a><br>
<div class="routine"><em class="call">model_size = get_model_size(
)</em>
<pre>
integer(i8) :: <em class="code">get_model_size</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the length of the model state vector as 1. Probably not
what you want. The model_mod should set this to the right size and
not use this routine.</p>
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
<p>Throws a fatal error. If the model_mod can advance the model it
should provide a real routine. This default routine is intended for
use by models which cannot advance themselves from inside
filter.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>State vector of length model_size.</td>
</tr>
<tr>
<td valign="top"><em class="code">time</em></td>
<td>Current time of model state.</td>
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
<p>Sets the location to missing and the variable type to 0. The
model_mod should provide a routine that sets a real location and a
state vector type for the requested item in the state vector.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">index_in</em></td>
<td>Index of state vector element about which information is
requested.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>The location of state variable element.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">var_type</em></td>
<td>The generic quantity of the state variable element.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="model_interpolate" id="model_interpolate"></a><br>
<div class="routine"><em class="call">call
model_interpolate(state_handle, ens_size, location, obs_quantity,
expected_obs, istatus)</em>
<pre>
type(ensemble_type),    intent(in)  :: <em class=
"code">state_handle</em>
integer,                intent(in)  :: <em class=
"code">ens_size</em>
type(location_type),    intent(in)  :: <em class=
"code">location</em>
integer,                intent(in)  :: <em class=
"code">obs_quantity</em>
real(r8),               intent(out) :: <em class=
"code">expected_obs(ens_size)</em>
integer,                intent(out) :: <em class=
"code">istatus(ens_size)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the expected obs to missing and returns an error code for
all obs. This routine should be supplied by the model_mod.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state_handle</em></td>
<td>The handle to the state structure containing information about
the state vector about which information is requested.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size</em></td>
<td>The ensemble size.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Location to which to interpolate.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_quantity</em></td>
<td>Quantity of state field to be interpolated.</td>
</tr>
<tr>
<td valign="top"><em class="code">expected_obs</em></td>
<td>The interpolated values from the model.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Integer values return 0 for success. Other positive values can
be defined for various failures.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="shortest_time_between_assimilations" id=
"shortest_time_between_assimilations"></a><br>
<div class="routine"><em class="call">var =
shortest_time_between_assimilations()</em>
<pre>
type(time_type) :: <em class=
"code">shortest_time_between_assimilations</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns 1 day.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Smallest advance time of the model.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="static_init_model" id="static_init_model"></a><br>
<div class="routine"><em class="call">call
static_init_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Does nothing.</p>
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
<p>Returns a time of 0.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">time</em></td>
<td>Initial model time.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="fail_init_time" id="fail_init_time"></a><br>
<div class="routine"><em class="call">call
fail_init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Throws a fatal error. This is appropriate for models that cannot
start from arbitrary initial conditions.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">time</em></td>
<td>NOT SET. Initial model time.</td>
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
<p>Returns x(:) = 0.0</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>Initial conditions for state vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="fail_init_conditions" id="fail_init_conditions"></a><br>
<div class="routine"><em class="call">call
fail_init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class="code">x</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Throws a fatal error. This is appropriate for models that cannot
start from arbitrary initial conditions.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>NOT SET: Initial conditions for state vector.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_atts" id="nc_write_model_atts"></a><br>
<div class="routine"><em class="call">call
nc_write_model_atts(ncFileID, domain_id)</em>
<pre>
integer, intent(in) :: <em class="code">ncFileID</em>
integer, intent(in) :: <em class="code">domain_id</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Does nothing.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>Integer file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">domain_id</em></td>
<td>integer describing the domain (which can be a nesting level, a
component model ...) Models with nested grids are decomposed into
'domains' in DART. The concept is extended to refer to 'coupled'
models where one model component may be the atmosphere, another
component may be the ocean, or land, or ionosphere ... these would
be referenced as different domains.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_vars" id="nc_write_model_vars"></a><br>
<div class="routine"><em class="call">call
nc_write_model_vars(ncFileID, domain_id, state_ens_handle
<em class="optionalcode">[, memberindex]</em> <em class=
"optionalcode">[, timeindex]</em>)</em>
<pre>
integer,             intent(in) :: <em class="code">ncFileID</em>
integer,             intent(in) :: <em class="code">domain_id</em>
type(ensemble_type), intent(in) :: <em class=
"code">state_ens_handle</em>
integer, optional,   intent(in) :: <em class=
"optionalcode">memberindex</em>
integer, optional,   intent(in) :: <em class=
"optionalcode">timeindex</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Does nothing</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">domain_id</em></td>
<td>integer describing the domain (which can be a nesting level, a
component model ...)</td>
</tr>
<tr>
<td valign="top"><em class="code">state_ens_handle</em></td>
<td>The handle to the state structure containing information about
the state vector about which information is requested.</td>
</tr>
<tr>
<td valign="top"><em class="code">memberindex</em></td>
<td>Integer index of ensemble member to be written.</td>
</tr>
<tr>
<td valign="top"><em class="code">timeindex</em></td>
<td>The timestep counter for the given state.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="pert_model_copies" id="pert_model_copies"></a><br>
<div class="routine"><em class="call">call
pert_model_copies(state_ens_handle, ens_size, pert_amp,
interf_provided)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">state_ens_handle</em>
integer,             intent(in)    :: <em class=
"code">ens_size</em>
real(r8),            intent(in)    :: <em class=
"code">pert_amp</em>
logical,             intent(out)   :: <em class=
"code">interf_provided</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns 'interface provided' flag as false, so the default
perturb routine in DART will add small amounts of gaussian noise to
all parts of the state vector.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state_ens_handle</em></td>
<td>The handle containing an ensemble of state vectors to be
perturbed.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size</em></td>
<td>The number of ensemble members to perturb.</td>
</tr>
<tr>
<td valign="top"><em class="code">pert_amp</em></td>
<td>the amplitude of the perturbations. The interpretation is based
on the model-specific implementation.</td>
</tr>
<tr>
<td valign="top"><em class="code">interf_provided</em></td>
<td>Returns false if model_mod cannot do this, else true.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 
<!-- these definitions came from the threed_sphere/location_mod.f90 -->
 <a name="get_close_obs" id="get_close_obs"></a><br>
<div class="routine"><em class="call">call get_close_obs(gc,
base_loc, base_type, locs, loc_qtys, loc_types, num_close,
close_ind <em class="optionalcode">[, dist]
[, state_handle</em>)</em>
<pre>
type(get_close_type),          intent(in)  :: <em class=
"code">gc</em>
type(location_type),           intent(in)  :: <em class=
"code">base_loc</em>
integer,                       intent(in)  :: <em class=
"code">base_type</em>
type(location_type),           intent(in)  :: <em class=
"code">locs(:)</em>
integer,                       intent(in)  :: <em class=
"code">loc_qtys(:)</em>
integer,                       intent(in)  :: <em class=
"code">loc_types(:)</em>
integer,                       intent(out) :: <em class=
"code">num_close</em>
integer,                       intent(out) :: <em class=
"code">close_ind(:)</em>
real(r8),            optional, intent(out) :: <em class=
"optionalcode">dist(:)</em>
type(ensemble_type), optional, intent(in)  :: <em class=
"optionalcode">state_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Passes the call through to the location module code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>The get_close_type which stores precomputed information about
the locations to speed up searching</td>
</tr>
<tr>
<td valign="top"><em class="code">base_loc</em></td>
<td>Reference location. The distances will be computed between this
location and every other location in the obs list</td>
</tr>
<tr>
<td valign="top"><em class="code">base_type</em></td>
<td>The DART quantity at the <em class="code">base_loc</em></td>
</tr>
<tr>
<td valign="top"><em class="code">locs(:)</em></td>
<td>Compute the distance between the <em class="code">base_loc</em>
and each of the locations in this list</td>
</tr>
<tr>
<td valign="top"><em class="code">loc_qtys(:)</em></td>
<td>The corresponding quantity of each item in the <em class=
"code">locs</em> list</td>
</tr>
<tr>
<td valign="top"><em class="code">loc_types(:)</em></td>
<td>The corresponding type of each item in the <em class=
"code">locs</em> list. This is not available in the default
implementation but may be used in custom implementations.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_close</em></td>
<td>The number of items from the <em class="code">locs</em> list
which are within maxdist of the base location</td>
</tr>
<tr>
<td valign="top"><em class="code">close_ind(:)</em></td>
<td>The list of index numbers from the <em class="code">locs</em>
list which are within maxdist of the base location</td>
</tr>
<tr>
<td valign="top"><em class="code">dist(:)</em></td>
<td>If present, return the distance between each entry in the
close_ind list and the base location. If not present, all items in
the obs list which are closer than maxdist will be added to the
list but the overhead of computing the exact distances will be
skipped.</td>
</tr>
<tr>
<td valign="top"><em class="code">state_handle</em></td>
<td>The handle to the state structure containing information about
the state vector about which information is requested.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_state" id="get_close_state"></a><br>
<div class="routine"><em class="call">call get_close_state(gc,
base_loc, base_type, state_loc, state_qtys, state_indx, num_close,
close_ind, dist, state_handle</em>)
<pre>
type(get_close_type), intent(in)    :: <em class="code">gc</em>
type(location_type),  intent(inout) :: <em class=
"code">base_loc</em>
integer,              intent(in)    :: <em class=
"code">base_type</em>
type(location_type),  intent(inout) :: <em class=
"code">state_loc(:)</em>
integer,              intent(in)    :: <em class=
"code">state_qtys(:)</em>
integer(i8),          intent(in)    :: <em class=
"code">state_indx(:)</em>
integer,              intent(out)   :: <em class=
"code">num_close</em>
integer,              intent(out)   :: <em class=
"code">close_ind(:)</em>
real(r8),             intent(out)   :: <em class=
"code">dist(:)</em>
type(ensemble_type),  intent(in)    :: <em class=
"code">state_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Passes the call through to the location module code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>The get_close_type which stores precomputed information about
the locations to speed up searching</td>
</tr>
<tr>
<td valign="top"><em class="code">base_loc</em></td>
<td>Reference location. The distances will be computed between this
location and every other location in the obs list</td>
</tr>
<tr>
<td valign="top"><em class="code">base_type</em></td>
<td>The DART quantity at the <em class="code">base_loc</em></td>
</tr>
<tr>
<td valign="top"><em class="code">state_loc(:)</em></td>
<td>Compute the distance between the <em class="code">base_loc</em>
and each of the locations in this list</td>
</tr>
<tr>
<td valign="top"><em class="code">state_qtys(:)</em></td>
<td>The corresponding quantity of each item in the <em class=
"code">state_loc</em> list</td>
</tr>
<tr>
<td valign="top"><em class="code">state_indx(:)</em></td>
<td>The corresponding DART index of each item in the <em class=
"code">state_loc</em> list. This is not available in the default
implementation but may be used in custom implementations.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_close</em></td>
<td>The number of items from the <em class="code">state_loc</em>
list which are within maxdist of the base location</td>
</tr>
<tr>
<td valign="top"><em class="code">close_ind(:)</em></td>
<td>The list of index numbers from the <em class=
"code">state_loc</em> list which are within maxdist of the base
location</td>
</tr>
<tr>
<td valign="top"><em class="code">dist(:)</em></td>
<td>If present, return the distance between each entry in the
<em class="code">close_ind</em> list and the base location. If not
present, all items in the <em class="code">state_loc</em> list
which are closer than maxdist will be added to the list but the
overhead of computing the exact distances will be skipped.</td>
</tr>
<tr>
<td valign="top"><em class="code">state_handle</em></td>
<td>The handle to the state structure containing information about
the state vector about which information is requested.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="convert_vertical_obs" id="convert_vertical_obs"></a><br>
<div class="routine"><em class="call">call
convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types,
which_vert, status)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=
"code">state_handle</em>
integer,             intent(in)  :: <em class="code">num</em>
type(location_type), intent(in)  :: <em class="code">locs(:)</em>
integer,             intent(in)  :: <em class=
"code">loc_qtys(:)</em>
integer,             intent(in)  :: <em class=
"code">loc_types(:)</em>
integer,             intent(in)  :: <em class=
"code">which_vert</em>
integer,             intent(out) :: <em class="code">status(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Passes the call through to the location module code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state_handle</em></td>
<td>The handle to the state.</td>
</tr>
<tr>
<td valign="top"><em class="code">num</em></td>
<td>the number of observation locations</td>
</tr>
<tr>
<td valign="top"><em class="code">locs</em></td>
<td>the array of observation locations</td>
</tr>
<tr>
<td valign="top"><em class="code">loc_qtys</em></td>
<td>the array of observation quantities.</td>
</tr>
<tr>
<td valign="top"><em class="code">loc_types</em></td>
<td>the array of observation types.</td>
</tr>
<tr>
<td valign="top"><em class="code">which_vert</em></td>
<td>the desired vertical coordinate system. There is a table in the
<em class="file">location_mod.f90</em> that relates integers to
vertical coordinate systems.</td>
</tr>
<tr>
<td valign="top"><em class="code">status</em></td>
<td>Success or failure of the vertical conversion. If <em class=
"code">istatus = 0</em>, the conversion was a success.
Any other value is a failure.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="convert_vertical_state" id=
"convert_vertical_state"></a><br>
<div class="routine"><em class="call">call
convert_vertical_state(state_handle, num, locs, loc_qtys,
loc_types, which_vert, status)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=
"code">state_handle</em>
integer,             intent(in)  :: <em class="code">num</em>
type(location_type), intent(in)  :: <em class="code">locs(:)</em>
integer,             intent(in)  :: <em class=
"code">loc_qtys(:)</em>
integer,             intent(in)  :: <em class=
"code">loc_types(:)</em>
integer,             intent(in)  :: <em class=
"code">which_vert</em>
integer,             intent(out) :: <em class="code">status(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Passes the call through to the location module code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state_handle</em></td>
<td>The handle to the state.</td>
</tr>
<tr>
<td valign="top"><em class="code">num</em></td>
<td>the number of state locations</td>
</tr>
<tr>
<td valign="top"><em class="code">locs</em></td>
<td>the array of state locations</td>
</tr>
<tr>
<td valign="top"><em class="code">loc_qtys</em></td>
<td>the array of state quantities.</td>
</tr>
<tr>
<td valign="top"><em class="code">loc_types</em></td>
<td>the array of state types.</td>
</tr>
<tr>
<td valign="top"><em class="code">which_vert</em></td>
<td>the desired vertical coordinate system. There is a table in the
<em class="file">location_mod.f90</em> that relates integers to
vertical coordinate systems.</td>
</tr>
<tr>
<td valign="top"><em class="code">status</em></td>
<td>Success or failure of the vertical conversion. If <em class=
"code">istatus = 0</em>, the conversion was a success.
Any other value is a failure.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_model_time" id="read_model_time"></a><br>
<div class="routine"><em class="call">model_time =
read_model_time(filename)</em>
<pre>
character(len=*), intent(in) :: <em class="code">filename</em>
type(time_type)              :: <em class="code">model_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Passes the call through to the dart_time_io module code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">filename</em></td>
<td>netCDF file name</td>
</tr>
<tr>
<td valign="top"><em class="code">model_time</em></td>
<td>The current time of the model state.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_model_time" id="write_model_time"></a><br>
<div class="routine"><em class="call">call write_model_time(ncid,
dart_time)</em>
<pre>
integer,          intent(in) :: <em class="code">ncid</em>
type(time_type),  intent(in) :: <em class="code">dart_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Passes the call through to the dart_time_io module code.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncid</em></td>
<td>handle to an open netCDF file</td>
</tr>
<tr>
<td valign="top"><em class="code">dart_time</em></td>
<td>The current time of the model state.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_model" id="end_model"></a><br>
<div class="routine"><em class="call">call end_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Does nothing.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<p><!-- useless spacer to make the top behave correctly --></p>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<a name="FilesUsed" id="FilesUsed"></a>
<h2>FILES</h2>
<p>none</p>
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
<p>Standard errors.</p>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none at this time.</p>
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
