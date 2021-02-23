<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_def_gps_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE <em class="program">obs_def_gps_mod</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Interface">INTERFACES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>DART GPS Radio Occultation observation module, including the
observation operators for both local and non-local refractivity
computations.</p>
<p>Author and Contact information:</p>
<ul>
<li>GPS Science: Hui Liu, hliu at ncar.edu</li>
<li>DART Code: Nancy Collins, nancy at ucar.edu</li>
</ul>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is now enabled by default. The maximum number of
GPS observations is settable at runtime by changing the value in
the namelist. If you get an error about a missing namelist add
<em class="code">&amp;obs_def_gps_nml</em> using the example below
to your <em class="file">input.nml</em> namelist file and rerun. No
recompiling is needed.</p>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;obs_def_gps_nml
  max_gpsro_obs = 100000,
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
<td>max_gpsro_obs</td>
<td>integer</td>
<td>The maximum number of GPS refractivity observations supported
for a single execution. Generally the default will be sufficient
for a single run of filter, but not enough for a long diagnostics
run to produce a time series.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
location_mod (threed_sphere)
assim_model_mod
obs_kind_mod
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
<td><em class="call">use obs_def_gps_mod, only :</em></td>
<td><a href="#read_gpsro_ref">read_gpsro_ref</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_gpsro_ref">write_gpsro_ref</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_expected_gpsro_ref">get_expected_gpsro_ref</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#interactive_gpsro_ref">interactive_gpsro_ref</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_gpsro_ref">set_gpsro_ref</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_gpsro_ref">get_gpsro_ref</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="read_gpsro_ref" id="read_gpsro_ref"></a><br>
<div class="routine"><em class="call">call read_gpsro_ref(gpskey,
ifile, <em class="optionalcode">[, fform]</em>)</em>
<pre>
integer,          intent(out)          :: <em class=
"code">gpskey</em>
integer,          intent(in)           :: <em class=
"code">ifile</em>
character(len=*), intent(in), optional :: <em class=
"optionalcode">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Refractivity observations have several items of auxiliary data
to read or write. This routine reads in the data for the next
observation and returns the private GPS key index number that
identifies the auxiliary data for this observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gpskey  </em></td>
<td>GPS key number returned to the caller.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>Open file unit number to read from.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">fform</em></td>
<td>If specified, indicate whether the file was opened formatted or
unformatted. Default is 'formatted'.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_gpsro_ref" id="write_gpsro_ref"></a><br>
<div class="routine"><em class="call">call write_gpsro_ref(gpskey,
ifile, <em class="optionalcode">[, fform]</em>)</em>
<pre>
integer,          intent(in)           :: <em class=
"code">gpskey</em>
integer,          intent(in)           :: <em class=
"code">ifile</em>
character(len=*), intent(in), optional :: <em class=
"optionalcode">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Refractivity observations have several items of auxiliary data
to read or write. This routine writes out the auxiliary data for
the specified observation to the file unit given.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gpskey  </em></td>
<td>GPS key number identifying which observation to write aux data
for.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>Open file unit number to write to.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">fform</em></td>
<td>If specified, indicate whether the file was opened formatted or
unformatted. Default is 'formatted'.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_expected_gpsro_ref" id=
"get_expected_gpsro_ref"></a><br>
<div class="routine"><em class="call">call
get_expected_gpsro_ref(state_vector, location, gpskey, ro_ref,
istatus)</em>
<pre>
real(r8),            intent(in)  :: <em class=
"code">state_vector(:)</em>
type(location_type), intent(in)  :: <em class="code">location</em>
integer,             intent(in)  :: <em class="code">gpskey</em>
real(r8),            intent(out) :: <em class="code">ro_ref</em>
integer,             intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a location and the state vector from one of the ensemble
members, compute the model-predicted GPS refractivity that would be
observed at that location. There are two types of operators:
modeled <em>local</em> refractivity (N-1)*1.0e6 or
<em>non_local</em> refractivity (excess phase, m) The type is
indicated in the auxiliary information for each observation.<br>
<br></p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">state_vector  </em></td>
<td>A one dimensional representation of the model state vector</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Location of this observation</td>
</tr>
<tr>
<td valign="top"><em class="code">gpskey</em></td>
<td>Integer key identifying which GPS observation this is, so the
correct corresponding auxiliary information can be accessed.</td>
</tr>
<tr>
<td valign="top"><em class="code">ro_ref</em></td>
<td>The returned GPS refractivity value</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Returned integer status code describing problems with applying
forward operator. 0 is a good value; any positive value indicates
an error; negative values are reserved for internal DART use
only.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="interactive_gpsro_ref" id=
"interactive_gpsro_ref"></a><br>
<div class="routine"><em class="call">call
interactive_gpsro_ref(gpskey)</em>
<pre>
integer, intent(out) :: <em class="code">gpskey</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Prompts the user for the auxiliary information needed for a GPS
refractivity observation, and returns the new key associated with
this data.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gpskey  </em></td>
<td>Unique identifier associated with this GPS refractivity
observation. In this code it is an integer index into module local
arrays which hold the additional data. This routine returns the
incremented value associated with this data.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_gpsro_ref" id="set_gpsro_ref"></a><br>
<div class="routine"><em class="call">call set_gpsro_ref(gpskey,
nx, ny, nz, rfict0, ds, htop, subset0)</em>
<pre>
integer,          intent(out) :: <em class="code">gpskey</em>
real(r8),         intent(in)  :: <em class="code">nx</em>
real(r8),         intent(in)  :: <em class="code">ny</em>
real(r8),         intent(in)  :: <em class="code">nz</em>
real(r8),         intent(in)  :: <em class="code">rfict0</em>
real(r8),         intent(in)  :: <em class="code">ds</em>
real(r8),         intent(in)  :: <em class="code">htop</em>
character(len=6), intent(in)  :: <em class="code">subset0</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the auxiliary information associated with a GPS
refractivity observation. This routine increments and returns the
new key associated with these values.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gpskey</em></td>
<td>Unique identifier associated with this GPS refractivity
observation. In this code it is an integer index into module local
arrays which hold the additional data. This routine returns the
incremented value associated with this data.</td>
</tr>
<tr>
<td valign="top"><em class="code">nx</em></td>
<td>X component of direction of ray between the LEO (detector)
satellite and the GPS transmitter satellite at the tangent
point.</td>
</tr>
<tr>
<td valign="top"><em class="code">ny</em></td>
<td>Y component of tangent ray.</td>
</tr>
<tr>
<td valign="top"><em class="code">nz</em></td>
<td>Z component of tangent ray.</td>
</tr>
<tr>
<td valign="top"><em class="code">rfict0</em></td>
<td>Local curvature radius (meters).</td>
</tr>
<tr>
<td valign="top"><em class="code">ds</em></td>
<td>Delta S, increment to move along the ray in each direction when
integrating the non-local operator (meters).</td>
</tr>
<tr>
<td valign="top"><em class="code">htop</em></td>
<td>Elevation (in meters) where integration stops along the
ray.</td>
</tr>
<tr>
<td valign="top"><em class="code">subset0</em></td>
<td>The string 'GPSREF' for the local operator (refractivity
computed only at the tangent point), or 'GPSEXC' for the non-local
operator which computes excess phase along the ray.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_gpsro_ref" id="get_gpsro_ref"></a><br>
<div class="routine"><em class="call">call get_gpsro_ref(gpskey,
nx, ny, nz, rfict0, ds, htop, subset0)</em>
<pre>
integer,          intent(in)  :: <em class="code">gpskey</em>
real(r8),         intent(out) :: <em class="code">nx</em>
real(r8),         intent(out) :: <em class="code">ny</em>
real(r8),         intent(out) :: <em class="code">nz</em>
real(r8),         intent(out) :: <em class="code">rfict0</em>
real(r8),         intent(out) :: <em class="code">ds</em>
real(r8),         intent(out) :: <em class="code">htop</em>
character(len=6), intent(out) :: <em class="code">subset0</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Gets the auxiliary information associated with a GPS
refractivity observation, based on the GPS key number
specified.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gpskey</em></td>
<td>Unique identifier associated with this GPS refractivity
observation. In this code it is an integer index into module local
arrays which hold the additional data. The value specified selects
which observation to return data for.</td>
</tr>
<tr>
<td valign="top"><em class="code">nx</em></td>
<td>X component of direction of ray between the LEO (detector)
satellite and the GPS transmitter satellite at the tangent
point.</td>
</tr>
<tr>
<td valign="top"><em class="code">ny</em></td>
<td>Y component of tangent ray.</td>
</tr>
<tr>
<td valign="top"><em class="code">nz</em></td>
<td>Z component of tangent ray.</td>
</tr>
<tr>
<td valign="top"><em class="code">rfict0</em></td>
<td>Local curvature radius (meters).</td>
</tr>
<tr>
<td valign="top"><em class="code">ds</em></td>
<td>Delta S, increment to move along the ray in each direction when
integrating the non-local operator (meters).</td>
</tr>
<tr>
<td valign="top"><em class="code">htop</em></td>
<td>Elevation (in meters) where integration stops along the
ray.</td>
</tr>
<tr>
<td valign="top"><em class="code">subset0</em></td>
<td>The string 'GPSREF' for the local operator (refractivity
computed only at the tangent point), or 'GPSEXC' for the non-local
operator which computes excess phase along the ray.</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>A DART observation sequence file containing GPS obs.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>Assimilation of GPS Radio Occultation Data for Numerical
Weather Prediction, Kuo,Y.H., Sokolovskiy,S.V., Anthes,R.A.,
Vendenberghe,F., Terrestrial Atm and Ocn Sciences, Vol 11,
pp157-186, 2000.</li>
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
<td valign="top">initialize_module</td>
<!-- message -->
<td valign="top">initial allocation failed for gps observation
data, itemcount = (max_gpsro_obs)</td>
<!-- comment -->
<td valign="top">Need to increase max_gpsro_obs count in
namelist</td>
</tr>
<tr><!-- routine -->
<td valign="top">gpskey_out_of_range</td>
<!-- message -->
<td valign="top">gpskey (key#) exceeds max_radial_gps_obs
(maxval)</td>
<!-- comment -->
<td valign="top">The number of GPS observations exceeds the array
size allocated in the module. Need to increase max_gpsro_obs count
in namelist.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_gpsro_ref</td>
<!-- message -->
<td valign="top">Expected header 'gpsroref' in input file</td>
<!-- comment -->
<td valign="top">The format of the input obs_seq file is not
consistent.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_expected_gpsro_ref</td>
<!-- message -->
<td valign="top">vertical location must be height; gps obs key
#</td>
<!-- comment -->
<td valign="top">GPS observations must have vertical coordinates of
height</td>
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
<p>The current code first bins the very densely-sampled vertical
profile into 200 bins, and then interpolates the requested vertical
location from that. The original profiles have been plotted and are
smooth; there appears to be no need to pre-bin the ata.</p>
<p>The local operator needs no additional auxiliary data. The
observation files would be much smaller if the local operator
observation was a separate type without aux data, and only the
non-local operator observation types would need the ray direction,
the curvature, etc.</p>
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
