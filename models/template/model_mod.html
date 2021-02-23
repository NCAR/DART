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
<p>Every model that is DART compliant must provide an interface as
documented here. The file <em class=
"file">models/template/model_mod.f90</em> provides the fortran
interfaces for a minimal implementation meeting these requirements.
When adding a new model to DART you can either start by modifying a
<em class="file">model_mod.f90</em> file from a similar model
already in DART or start with the template file. Either way, the
supplied interface must match these descriptions exactly; no
details of the underlying model can impact the interface.<br>
<br>
Several of the routines listed below are allowed to be a NULL
INTERFACE. This means the subroutine or function name must exist in
this file, but it is ok if it contains no executable code.<br>
<br>
A few of the routines listed below are allowed to be a PASS-THROUGH
INTERFACE. This means the subroutine or function name can be listed
on the 'use' line from the <em class="code">location_mod</em>, and
no subroutine or function with that name is supplied in this file.
Alternatively, this file can provide an implementation which calls
the underlying routines from the <em class="code">location_mod</em>
and then alters or augments the results based on model-specific
requirements.<br>
<br>
The system comes with several types of location modules for
computing distances appropriately. Two of the ones most commonly
used are for data in a 1D system and for data in a 3D spherical
coordinate system. Make the selection by listing the appropriate
choice from <em class="file">location/*/location_mod.f90</em> in
the corresponding <em class="file">path_names_*</em> file at
compilation time.</p>
<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->
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
 /
</pre></div>
<br>
<br>
<p>Models are free to include a model namelist which can be read
when <em class="code">static_init_model</em> is called. A good
example can be found in the lorenz_96 <em class=
"file">model_mod.f90</em>.</p>
<!--==================================================================-->
<p><!-- useless spacer to make the top behave correctly --></p>
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
time_manager_mod
location_mod (multiple choices here)
utilities_mod
POSSIBLY MANY OTHERS DEPENDING ON MODEL DETAILS
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
<p>A namelist interface <a href="#Namelist"><em class=
"code">&amp;model_nml</em></a> may be defined by the module, in
which case it will be read from file <em class=
"file">input.nml</em>. The details of the namelist are always
model-specific (there are no generic namelist values).</p>
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
<p>Returns the length of the model state vector. Required.</p>
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
<p>Does a single timestep advance of the model. The input value of
the vector x is the starting condition and x must be updated to
reflect the changed state after a timestep. The time argument is
intent in and is used for models that need to know the date/time to
compute a timestep, for instance for radiation computations. This
interface is only called if the namelist parameter async is set to
0 in <em class="program">perfect_model_obs</em> or <em class=
"program">filter</em> or if the program <em class=
"program">integrate_model</em> is to be used to advance the model
state as a separate executable. If one of these options is not
going to be used (the model will <em>only</em> be advanced as a
separate model-specific executable), this can be a NULL INTERFACE.
(The subroutine name must still exist, but it can contain no code
and it will not be called.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>State vector of length model_size.</td>
</tr>
<tr>
<td valign="top"><em class="code">time</em></td>
<td>Current time of the model state.</td>
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
<p>Given an integer index into the state vector, returns the
associated location. An optional argument returns the generic
quantity of this item, e.g. QTY_TEMPERATURE, QTY_DENSITY,
QTY_SALINITY, QTY_U_WIND_COMPONENT. This interface is required to
be functional for all applications.</p>
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
<p>Given a handle containing information for a state vector, an
ensemble size, a location, and a model state variable quantity
interpolates the state variable field to that location and returns
an ensemble-sized array of values in <em class=
"code">expected_obs(:)</em>. The <em class="code">istatus(:)</em>
array should be 0 for successful ensemble members and a positive
value for failures. The <em class="code">obs_quantity</em> variable
is one of the quantity (QTY) parameters defined in the <a href=
"../../assimilation_code/modules/observations/obs_kind_mod.html">obs_kind_mod.f90</a>
file and defines the quantity to interpolate. In low-order models
that have no notion of kinds of variables this argument may be
ignored. For applications in which only perfect model experiments
with identity observations (i.e. only the value of a particular
state variable is observed), this can be a NULL INTERFACE.
Otherwise it is required (which is the most common case).</p>
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
<p>Returns the smallest increment in time that the model is capable
of advancing the state in a given implementation. The actual value
may be set by the model_mod namelist (depends on the model). This
interface is required for all applications.</p>
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
<p>Called to do one time initialization of the model. As examples,
might define information about the model size or model timestep,
read in grid information, read a namelist, set options, etc. In
models that require pre-computed static data, for instance
spherical harmonic weights, these would also be computed here. Can
be a NULL INTERFACE for the simplest models.</p>
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
<p>Companion interface to init_conditions. Returns a time that is
somehow appropriate for starting up a long integration of the
model. At present, this is only used if the <em class=
"program">perfect_model_obs</em> namelist parameter <em class=
"code">read_input_state_from_file = .false.</em> If this
option should not be used in <em class=
"program">perfect_model_obs</em>, calling this routine should issue
a fatal error.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">time</em></td>
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
<p>Returns a model state vector, x, that is some sort of
appropriate initial condition for starting up a long integration of
the model. At present, this is only used if the <em class=
"program">perfect_model_obs</em> namelist parameter <em class=
"code">read_input_state_from_file = .false.</em> If this
option should not be used in <em class=
"program">perfect_model_obs</em>, calling this routine should issue
a fatal error.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>Initial conditions for state vector.</td>
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
<p>This routine writes the model-specific attributes to netCDF
files that DART creates. This includes coordinate variables and any
metadata, but NOT the actual model state vector. <em class=
"code">models/template/model_mod.f90</em> contains code that can be
used for any model as-is.<br>
<br>
The typical sequence for adding new dimensions, variables,
attributes:</p>
<pre>
NF90_OPEN             ! open existing netCDF dataset               
   NF90_redef         ! put into define mode                       
   NF90_def_dim       ! define additional dimensions (if any)     
   NF90_def_var       ! define variables: from name, kind, and dims
   NF90_put_att       ! assign attribute values                    
NF90_ENDDEF           ! end definitions: leave define mode         
   NF90_put_var       ! provide values for variable                
NF90_CLOSE            ! close: save updated netCDF dataset        
</pre>
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
<p>This routine may be used to write the model-specific state
vector (data) to a netCDF file. Only used if <em class=
"code">model_mod_writes_state_variables = .true.</em><br>
<br>
Typical sequence for adding new
dimensions,variables,attributes:</p>
<pre>
NF90_OPEN             ! open existing netCDF dataset               
   NF90_redef         ! put into define mode                       
   NF90_def_dim       ! define additional dimensions (if any)      
   NF90_def_var       ! define variables: from name, kind, and dims
   NF90_put_att       ! assign attribute values                    
NF90_ENDDEF           ! end definitions: leave define mode         
   NF90_put_var       ! provide values for variable                
NF90_CLOSE            ! close: save updated netCDF dataset         
</pre>
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
<p>Given an ensemble handle, the ensemble size, and a perturbation
amplitude; perturb the ensemble. Used to generate initial
conditions for spinning up ensembles. If the <em class=
"code">model_mod</em> does not want to do this, instead allowing
the default algorithms in <em class="program">filter</em> to take
effect, <em class=
"code">interf_provided =&amp;nbps;.false.</em> and the routine
can be trivial. Otherwise, <em class="code">interf_provided</em>
must be returned as <em class="code">.true.</em></p>
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
<p>Given a location and quantity, compute the distances to all
other locations in the <em class="code">obs</em> list. The return
values are the number of items which are within maxdist of the
base, the index numbers in the original obs list, and optionally
the distances. The <em class="code">gc</em> contains precomputed
information to speed the computations.<br>
<br>
In general this is a PASS-THROUGH ROUTINE. It is listed on the use
line for the locations_mod, and in the public list for this module,
but has no subroutine declaration and no other code in this
module:</p>
<pre>
use location_mod, only: get_close_obs

public :: get_close_obs
</pre>
<p>However, if the model needs to alter the values or wants to
supply an alternative implementation it can intercept the call like
so:</p>
<pre>
use location_mod, only: &amp;
        lm_get_close_obs =&gt; get_close_obs
        
public :: get_close_obs
</pre>
<p>In this case a local <em class="code">get_close_obs()</em>
routine must be supplied. To call the original code in the location
module use:</p>
<pre>
call lm_get_close_obs(gc, base_loc, ...)
</pre>
<p>This subroutine will be called after <em class=
"code">get_close_maxdist_init</em> and <em class=
"code">get_close_obs_init</em>.<br>
<br>
In most cases the PASS-THROUGH ROUTINE will be used, but some
models need to alter the actual distances depending on the
observation or state vector kind, or based on the observation or
state vector location. It is reasonable in this case to leave
<em class="code">get_close_maxdist_init()</em> and <em class=
"code">get_close_obs_init()</em> as pass-through routines and
intercept only <em class="code">get_close_obs()</em>. The local
<em class="code">get_close_obs()</em> can first call the location
mod routine and let it return a list of values, and then inspect
the list and alter or remove any entries as needed. See the CAM and
WRF model_mod files for examples of this use.</p>
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
close_ind <em class=
"optionalcode">[, dist, state_handle]</em>)</em>
<pre>
type(get_close_type),          intent(in)    :: <em class=
"code">gc</em>
type(location_type),           intent(inout) :: <em class=
"code">base_loc</em>
integer,                       intent(in)    :: <em class=
"code">base_type</em>
type(location_type),           intent(inout) :: <em class=
"code">state_loc(:)</em>
integer,                       intent(in)    :: <em class=
"code">state_qtys(:)</em>
integer(i8),                   intent(in)    :: <em class=
"code">state_indx(:)</em>
integer,                       intent(out)   :: <em class=
"code">num_close</em>
integer,                       intent(out)   :: <em class=
"code">close_ind(:)</em>
real(r8),            optional, intent(out)   :: <em class=
"code">dist(:)</em>
type(ensemble_type), optional, intent(in)    :: <em class=
"code">state_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a location and quantity, compute the distances to all
other locations in the <em class="code">state_loc</em> list. The
return values are the number of items which are within maxdist of
the base, the index numbers in the original state_loc list, and
optionally the distances. The <em class="code">gc</em> contains
precomputed information to speed the computations.<br>
<br>
In general this is a PASS-THROUGH ROUTINE. It is listed on the use
line for the locations_mod, and in the public list for this module,
but has no subroutine declaration and no other code in this
module:</p>
<pre>
use location_mod, only: get_close_state

public :: get_close_state
</pre>
<p>However, if the model needs to alter the values or wants to
supply an alternative implementation it can intercept the call like
so:</p>
<pre>
use location_mod, only: &amp;
        lm_get_close_state =&gt; get_close_state
        
public :: get_close_state
</pre>
<p>In this case a local <em class="code">get_close_state()</em>
routine must be supplied. To call the original code in the location
module use:</p>
<pre>
call loc_get_close_state(gc, base_loc, ...)
</pre>
<p>This subroutine will be called after <em class=
"code">get_close_maxdist_init</em> and <em class=
"code">get_close_state_init</em>.<br>
<br>
In most cases the PASS-THROUGH ROUTINE will be used, but some
models need to alter the actual distances depending on the
observation or state vector kind, or based on the observation or
state vector location. It is reasonable in this case to leave
<em class="code">get_close_maxdist_init()</em> and <em class=
"code">get_close_state_init()</em> as pass-through routines and
intercept only <em class="code">get_close_state()</em>. The local
<em class="code">get_close_state()</em> can first call the location
mod routine and let it return a list of values, and then inspect
the list and alter or remove any entries as needed. See the CAM and
WRF model_mod files for examples of this use.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>The get_close_type which stores precomputed information about
the locations to speed up searching</td>
</tr>
<tr>
<td valign="top"><em class="code">base_loc</em></td>
<td>Reference location. The distances will be computed between this
location and every other location in the list</td>
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
<p>Converts the observations to the desired vertical localization
coordinate system. Some models (toy models with no 'real'
observations) will not need this. Most (real) models have
observations in one or more coordinate systems (pressure, height)
and the model is generally represented in only one coordinate
system. To be able to interpolate the model state to the
observation location, or to compute the true distance between the
state and the observation, it is necessary to convert everything to
a single coodinate system.</p>
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
<p>Converts the state to the desired vertical localization
coordinate system. Some models (toy models with no 'real'
observations) will not need this. To compute the true distance
between the state and the observation, it is necessary to convert
everything to a single coodinate system.</p>
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
<p>Reads the valid time of the model state in a netCDF file. There
is a default routine in <em class=
"file">assimilation_code/modules/io/dart_time_io_mod.f90</em> that
can be used as a pass-through. That routine will read the
<strong>last</strong> timestep of a 'time' variable - which is the
same strategy used for reading netCDF files that have multiple
timesteps in them. If your model has some other representation of
time (i.e. it does not use a netCDF variable named 'time') - you
will have to write this routine.</p>
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
 <a name="write_model_time" id="write_model_time"></a><br>
<div class="routine"><em class="call">call write_model_time(ncid,
dart_time)</em>
<pre>
integer,          intent(in) :: <em class="code">ncid</em>
type(time_type),  intent(in) :: <em class="code">dart_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes the assimilation time to a netCDF file. There is a
default routine in <em class=
"file">assimilation_code/modules/io/dart_time_io_mod.f90</em> that
can be used as a pass-through. If your model has some other
representation of time (i.e. it does not use a netCDF variable
named 'time') - you will have to write this routine.</p>
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
<p>Does any shutdown and clean-up needed for model. Can be a NULL
INTERFACE if the model has no need to clean up storage, etc.</p>
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
<ul>
<li>Models are free to read and write files as they see fit.</li>
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
<ul>
<li>Models are free to issue calls to the error handler as they see
fit. No standard error handler calls are mandated.</li>
</ul>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>It is likely that a number of additional optional interfaces
will be added to the model_mod structure. For instance, hints about
how to divide the state vector into regions for parallel
assimilation will need to be obtained from the model. It is planned
that the interf_provided mechanism used in pert_model_copies will
allow those who do not wish to support enhanced interfaces to add
NULL interfaces by simply pasting in an interface block.</p>
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
