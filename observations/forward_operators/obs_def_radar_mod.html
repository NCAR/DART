<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_def_radar_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE <em class="program">obs_def_radar_mod</em></h1>
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
<p>DART radar observation module, including the observation
operators for the two primary radar-observation types -- Doppler
velocity and reflectivity -- plus other utility subroutines and
functions. A number of simplifications are employed for the
observation operators. Most notably, the model state is mapped to a
"point" observation, whereas a real radar observation is a
volumetric sample. The implications of this approximation have not
been investigated fully, so in the future it might be worth
developing and testing more sophisticated observation operators
that produce volumetric power- weighted samples.<br>
<br>
This module is able to compute reflectivity and precipitation fall
speed (needed for computing Doppler radial velocity) from the
prognostic model fields only for simple single-moment microphysics
schemes such as the Kessler and Lin schemes. If a more complicated
microphysics scheme is used, then reflectivity and fall speed must
be accessible instead as diagnostic fields in the model state.<br>
<br>
Author and Contact information:</p>
<ul>
<li>Radar Science: David Dowell, david.dowell at noaa.gov, Glen
Romine, romine at ucar.edu</li>
<li>DART Code: Nancy Collins, nancy at ucar.edu</li>
<li>Original DART/Radar work: Alain Caya</li>
</ul>
<h3>Backward compatibility note:</h3>
<p>For users of previous versions of the radar obs_def code, here
are a list of changes beginning with subversion revision 3616 which
are not backward compatible:</p>
<ul>
<li>The namelist has changed quite a bit; some items were removed,
some added, and some renamed. See the <a href="#Namelist">namelist
documention</a> in this file for the current item names and default
values.</li>
<li>Some constants which depend on the microphysics scheme have
been added to the namelist to make it easier to change the values
for different schemes, but the defaults have also changed. Verify
they are appropriate for the scheme being used.</li>
<li>The interactive create routine prompts for the beam direction
differently now. It takes azimuth and elevation, and then does the
trigonometry to compute the three internal values which are stored
in the file. The previous version prompted for the internal values
directly.</li>
<li>The get_expected routines try to call the model interpolate
routine for <em class="code">QTY_POWER_WEIGHTED_FALL_SPEED</em> and
<em class="code">QTY_RADAR_REFLECTIVITY</em> values. If they are
not available then the code calls the model interpolate routines
for several other quantities and computes these quantities.
However, this requires that the model_mod interpolate code returns
gracefully if the quantity is unknown or unsupported. The previous
version of the WRF model_mod code used to print an error message
and stop if the quantity was unknown. The updated version in the
repository which went in with this radar code has been changed to
return an error status code but continue if the quantity is
unknown.</li>
<li>The value for gravity is currently hardcoded in this module.
Previous versions of this code used the gravity constant in the
DART types_mod.f90 code, but in reality the code should be using
whatever value of gravity is being used in the model code. For now,
the value is at least separated so users can change the value in
this code if necessary.</li>
</ul>
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
<td><em class="call">use obs_def_radar_mod, only :</em></td>
<td><a href="#read_radar_ref">read_radar_ref</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_expected_radar_ref">get_expected_radar_ref</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_radial_vel">read_radial_vel</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_radial_vel">write_radial_vel</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#interactive_radial_vel">interactive_radial_vel</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_expected_radial_vel">get_expected_radial_vel</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_obs_def_radial_vel">get_obs_def_radial_vel</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_radial_vel">set_radial_vel</a></td>
</tr>
</table>
<p>Namelist interface <a href="#Namelist"><em class=
"code">&amp;obs_def_radar_mod_nml</em></a> is read from file
<em class="file">input.nml</em>.</p>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="read_radar_ref" id="read_radar_ref"></a><br>
<div class="routine"><em class="call">call read_radar_ref(obsvalue,
refkey)</em>
<pre>
real(r8),                   intent(inout) :: <em class=
"code">obsvalue</em>
integer,                    intent(out)   :: <em class=
"code">refkey</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reflectivity observations have no auxiliary data to read or
write, but there are namelist options that can alter the
observation value at runtime. This routine tests the observation
value and alters it if required.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obsvalue  </em></td>
<td>Observation value.</td>
</tr>
<tr>
<td valign="top"><em class="code">refkey</em></td>
<td>Set to 0 to avoid uninitialized values, but otherwise
unused.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_expected_radar_ref" id=
"get_expected_radar_ref"></a><br>
<div class="routine"><em class="call">call
get_expected_radar_ref(state_vector, location, ref, istatus)</em>
<pre>
real(r8),            intent(in)  :: <em class=
"code">state_vector(:)</em>
type(location_type), intent(in)  :: <em class="code">location</em>
real(r8),            intent(out) :: <em class="code">ref</em>
integer,             intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a location and the state vector from one of the ensemble
members, compute the model-predicted radar reflectivity that would
be observed at that location. The returned value is in dBZ.<br>
<br>
If <em class="code">apply_ref_limit_to_fwd_op</em> is .TRUE. in the
namelist, reflectivity values less than <em class=
"code">reflectivity_limit_fwd_op</em> will be set to <em class=
"code">lowest_reflectivity_fwd_op</em>.</p>
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
<td valign="top"><em class="code">ref</em></td>
<td>The returned radar reflectivity value</td>
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
 <a name="read_radial_vel" id="read_radial_vel"></a><br>
<div class="routine"><em class="call">call read_radial_vel(velkey,
ifile <em class="optionalcode">[, fform]</em>)</em>
<pre>
integer,                    intent(out) :: <em class=
"code">velkey</em>
integer,                    intent(in)  :: <em class=
"code">ifile</em>
character(len=*), optional, intent(in)  :: <em class=
"optionalcode">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads the additional auxiliary information associated with a
radial velocity observation. This includes the location of the
radar source, the beam direction, and the nyquist velocity.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">velkey  </em></td>
<td>Unique identifier associated with this radial velocity
observation. In this code it is an integer index into module local
arrays which hold the additional data. This routine increments it
and returns the new value.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>File unit descriptor for input file</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">fform</em></td>
<td>File format specifier: FORMATTED or UNFORMATTED; default
FORMATTED</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_radial_vel" id="write_radial_vel"></a><br>
<div class="routine"><em class="call">call write_radial_vel(velkey,
ifile <em class="optionalcode">[, fform]</em>)</em>
<pre>
integer,                    intent(in) :: <em class=
"code">velkey</em>
integer,                    intent(in) :: <em class=
"code">ifile</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes the additional auxiliary information associated with a
radial velocity observation. This includes the location of the
radar source, the beam direction, and the nyquist velocity.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">velkey  </em></td>
<td>Unique identifier associated with this radial velocity
observation. In this code it is an integer index into module local
arrays which hold the additional data. This routine uses the value
to select the appropriate data to write for this observation.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>File unit descriptor for output file</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">fform</em></td>
<td>File format specifier: FORMATTED or UNFORMATTED; default
FORMATTED</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_obs_def_radial_vel" id=
"get_obs_def_radial_vel"></a><br>
<div class="routine"><em class="call">call
get_obs_def_radial_vel(velkey, radar_location, beam_direction,
nyquist_velocity)</em>
<pre>
integer,             intent(in)  :: <em class="code">velkey</em>
type(location_type), intent(out) :: <em class=
"code">radar_location</em>
real(r8),            intent(out) :: <em class=
"code">beam_direction(3)</em>
real(r8),            intent(out) :: <em class=
"code">nyquist_velocity</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the auxiliary information associated with a given radial
velocity observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">velkey</em></td>
<td>Unique identifier associated with this radial velocity
observation. In this code it is an integer index into module local
arrays which hold the additional data. This routine uses the value
to select the appropriate data to return.</td>
</tr>
<tr>
<td valign="top"><em class="code">radar_location</em></td>
<td>Location of the radar.</td>
</tr>
<tr>
<td valign="top"><em class="code">beam_orientation</em></td>
<td>Orientation of the radar beam at the observation location. The
three values are: sin(azimuth)*cos(elevation),
cos(azimuth)*cos(elevation), and sin(elevation).</td>
</tr>
<tr>
<td valign="top"><em class=
"code">nyquist_velocity  </em></td>
<td>Nyquist velocity at the observation point in
meters/second.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_radial_vel" id="set_radial_vel"></a><br>
<div class="routine"><em class="call">call set_radial_vel(velkey,
radar_location, beam_direction, nyquist_velocity)</em>
<pre>
integer,             intent(out) :: <em class="code">velkey</em>
type(location_type), intent(in)  :: <em class=
"code">radar_location</em>
real(r8),            intent(in)  :: <em class=
"code">beam_direction(3)</em>
real(r8),            intent(in)  :: <em class=
"code">nyquist_velocity</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the auxiliary information associated with a radial velocity
observation. This routine increments and returns the new key
associated with these values.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">velkey</em></td>
<td>Unique identifier associated with this radial velocity
observation. In this code it is an integer index into module local
arrays which hold the additional data. This routine returns the
incremented value associated with this data.</td>
</tr>
<tr>
<td valign="top"><em class="code">radar_location</em></td>
<td>Location of the radar.</td>
</tr>
<tr>
<td valign="top"><em class="code">beam_orientation</em></td>
<td>Orientation of the radar beam at the observation location. The
three values are: sin(azimuth)*cos(elevation),
cos(azimuth)*cos(elevation), and sin(elevation).</td>
</tr>
<tr>
<td valign="top"><em class=
"code">nyquist_velocity  </em></td>
<td>Nyquist velocity at the observation point in
meters/second.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="interactive_radial_vel" id=
"interactive_radial_vel"></a><br>
<div class="routine"><em class="call">call
interactive_radial_vel(velkey)</em>
<pre>
integer, intent(out) :: <em class="code">velkey</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Prompts the user for the auxiliary information needed for a
radial velocity observation, and returns the new key associated
with this data.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">velkey  </em></td>
<td>Unique identifier associated with this radial velocity
observation. In this code it is an integer index into module local
arrays which hold the additional data. This routine returns the
incremented value associated with this data.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_expected_radial_vel" id=
"get_expected_radial_vel"></a><br>
<div class="routine"><em class="call">call
get_expected_radial_vel(state_vector, location, velkey, radial_vel,
istatus)</em>
<pre>
real(r8),            intent(in)  :: <em class=
"code">state_vector(:)</em>
type(location_type), intent(in)  :: <em class="code">location</em>
integer,             intent(in)  :: <em class="code">velkey</em>
real(r8),            intent(out) :: <em class=
"code">radial_vel</em>
integer,             intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a location and the state vector from one of the ensemble
members, compute the model-predicted radial velocity in
meters/second that would be observed at that location. <em class=
"code">velkey</em> is the unique index for this particular radial
velocity observation. The value is returned in <em class=
"code">radial_vel</em>, <em class="code">istatus</em> is the return
code.<br>
<br>
The along-beam component of the 3-d air velocity is computed from
the u, v, and w fields plus the beam_direction. The along-beam
component of power-weighted precipitation fall velocity is added to
the result.</p>
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
<td valign="top"><em class="code">velkey</em></td>
<td>Unique identifier associated with this radial velocity
observation</td>
</tr>
<tr>
<td valign="top"><em class="code">radial_vel</em></td>
<td>The returned radial velocity value in meters/second</td>
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
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
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
&amp;obs_def_radar_mod_nml
   apply_ref_limit_to_obs      =   .false.,
   reflectivity_limit_obs      =     -10.0,
   lowest_reflectivity_obs     =     -10.0,
   apply_ref_limit_to_fwd_op   =   .false.,
   reflectivity_limit_fwd_op   =     -10.0,
   lowest_reflectivity_fwd_op  =     -10.0,
   max_radial_vel_obs          =   1000000,
   allow_wet_graupel           =   .false.,
   microphysics_type           =       2  ,
   allow_dbztowt_conv          =   .false.,
   dielectric_factor           =     0.224,
   n0_rain                     =     8.0e6,
   n0_graupel                  =     4.0e6,
   n0_snow                     =     3.0e6,
   rho_rain                    =    1000.0,
   rho_graupel                 =     400.0,
   rho_snow                    =     100.0 
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
<td>apply_ref_limit_to_obs</td>
<td>logical</td>
<td>If .TRUE. replace all reflectivity values less than
"reflectivity_limit_obs" with "lowest_reflectivity_obs" value. If
.FALSE. leave all values as-is.</td>
</tr>
<tr>
<td>reflectivity_limit_obs</td>
<td>real(r8)</td>
<td>The threshold value. Observed reflectivity values less than
this threshold will be set to the "lowest_reflectivity_obs" value.
Units are dBZ.</td>
</tr>
<tr>
<td>lowest_reflectivity_obs</td>
<td>real(r8)</td>
<td>The 'set-to' value. Observed reflectivity values less than the
threshold will be set to this value. Units are dBZ.</td>
</tr>
<tr>
<td>apply_ref_limit_to_fwd_op</td>
<td>logical</td>
<td>Same as "apply_ref_limit_to_obs", but for the forward
operator.</td>
</tr>
<tr>
<td>reflectivity_limit_fwd_op</td>
<td>real(r8)</td>
<td>Same as "reflectivity_limit_obs", but for the forward operator
values.</td>
</tr>
<tr>
<td>lowest_reflectivity_fwd_op</td>
<td>real(r8)</td>
<td>Same as "lowest_reflectivity_obs", but for the forward operator
values.</td>
</tr>
<tr>
<td>max_radial_vel_obs</td>
<td>integer</td>
<td>Maximum number of observations of this type to support at run
time. This is combined total of all obs_seq files, for example the
observation diagnostic program potentially opens multiple
obs_seq.final files, or the obs merge program can also open
multiple obs files.</td>
</tr>
<tr>
<td>allow_wet_graupel</td>
<td>logical</td>
<td>It is difficult to predict/diagnose whether graupel/hail has a
wet or dry surface. Even when the temperature is above freezing,
evaporation and/or absorption can still result in a dry surface.
This issue is important because the reflectivity from graupel with
a wet surface is significantly greater than that from graupel with
a dry surface. Currently, the user has two options for how to
compute graupel reflectivity. If allow_wet_graupel is .false. (the
default), then graupel is always assumed to be dry. If
allow_wet_graupel is .true., then graupel is assumed to be wet
(dry) when the temperature is above (below) freezing. A consequence
is that a sharp gradient in reflectivity will be produced at the
freezing level. In the future, it might be better to provide the
option of having a transition layer.</td>
</tr>
<tr>
<td>microphysics_type</td>
<td>integer</td>
<td>If the state vector contains the reflectivity or the power
weighted fall speed, interpolate directly from those regardless of
the setting of this item. If the state vector does not contain the
fields, this value should be set to be compatible with whatever
microphysical scheme is being used by the model. If the model is
using a different microphysical scheme but has compatible fields to
the ones listed below, setting this value will select the scheme to
use.
<ul style="list-style: none;">
<li>1 = Kessler scheme.</li>
<li>2 = Lin et al. microphysics</li>
<li>3 = User selected scheme where 10 cm reflectivity and power
weighted fall velocity are expected in the state vector (failure if
not found)</li>
<li>4 = User selected scheme where only power weighted fall
velocity is expected (failure if not found)</li>
<li>5 = User selected scheme where only reflectivity is expected
(failure if not found)</li>
<li>-1 = ASSUME FALL VELOCITY IS ZERO, allows over-riding the
failure modes above if reflectivity and/or fall velocity are not
available but a result is desired for testing purposes only.</li>
</ul>
</td>
</tr>
<tr>
<td>allow_dbztowt_conv</td>
<td>logical</td>
<td>Flag to enable use of the dbztowt routine where reflectivity is
available, but not the power-weighted fall velocity. This scheme
uses emperical relations between reflectivity and fall velocity,
with poor accuracy for highly reflective, low density particles
(such as water coated snow aggregates). Expect questionable
accuracy in radial velocity from the forward operator with high
elevation angles where ice is present in the model state.</td>
</tr>
<tr>
<td>dielectric_factor</td>
<td>real(r8)</td>
<td>According to Smith (1984), there are two choices for the
dielectric factor depending on how the snow particle sizes are
specified. If melted raindrop diameters are used, then the factor
is 0.224. If equivalent ice sphere diameters are used, then the
factor is 0.189. The default is set to use the common convention of
melted raindrop diameters.</td>
</tr>
<tr>
<td>n0_rain</td>
<td>real(r8)</td>
<td>Intercept parameters (m^-4) for size distributions of each
hydrometeor. The default of 8.0e6 is for the Lin et al.
microphysics scheme with the Hobbs settings for graupel/hail. (The
Hobbs graupel settings are also the default for the Lin scheme in
WRF 2.2 and 3.0.)</td>
</tr>
<tr>
<td>n0_graupel</td>
<td>real(r8)</td>
<td>Intercept parameters (m^-4) for size distributions of each
hydrometeor. The default of 4.0e6 is for the Lin et al.
microphysics scheme with the Hobbs settings for graupel/hail. (The
Hobbs graupel settings are also the default for the Lin scheme in
WRF 2.2 and 3.0.)</td>
</tr>
<tr>
<td>n0_snow</td>
<td>real(r8)</td>
<td>Intercept parameters (m^-4) for size distributions of each
hydrometeor. The default of 3.0e6 is for the Lin et al.
microphysics scheme with the Hobbs settings for graupel/hail. (The
Hobbs graupel settings are also the default for the Lin scheme in
WRF 2.2 and 3.0.)</td>
</tr>
<tr>
<td>rho_rain</td>
<td>real(r8)</td>
<td>Density (kg m^-3) of each hydrometeor type. The default of
1000.0 is for the Lin et al. microphysics scheme with the Hobbs
setting for graupel/hail.</td>
</tr>
<tr>
<td>rho_graupel</td>
<td>real(r8)</td>
<td>Density (kg m^-3) of each hydrometeor type. The default of
400.0 is for the Lin et al. microphysics scheme with the Hobbs
setting for graupel/hail.</td>
</tr>
<tr>
<td>rho_snow</td>
<td>real(r8)</td>
<td>Density (kg m^-3) of each hydrometeor type. The default of
100.0 is for the Lin et al. microphysics scheme with the Hobbs
setting for graupel/hail.</td>
</tr>
<!-- FIXME: add max values, remove dbz, missing r8 probably not good val? --></tbody>
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
<ul>
<li>A DART observation sequence file containing Radar obs.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>Battan, L. J., 1973: <i>Radar Observation of the
Atmosphere.</i> Univ. of Chicago Press, 324 pp.</li>
<li>Caya, A. <i>Radar Observations in Dart.</i> DART Subversion
repository.</li>
<li>Doviak, R. J., and D. S. Zrnic, 1993: <i>Doppler Radar and
Weather Observations.</i> Academic Press, 562 pp.</li>
<li>Ferrier, B. S., 1994: A double-moment multiple-phase four-class
bulk ice scheme. Part I: Description. <i>J. Atmos. Sci.</i>,
<b>51</b>, 249-280.</li>
<li>Lin, Y.-L., Farley R. D., and H. D. Orville, 1983: Bulk
parameterization of the snow field in a cloud model. <i>J. Climate
Appl. Meteor.</i>, <b>22</b>, 1065-1092.</li>
<li>Smith, P. L. Jr., 1984: Equivalent radar reflectivity factors
for snow and ice particles. <i>J. Climate Appl. Meteor.</i>, 23,
1258-1260.</li>
<li>Smith, P. L. Jr., Myers C. G., and H. D. Orville, 1975: Radar
reflectivity factor calculations in numerical cloud models using
bulk parameterization of precipitation. <i>J. Appl. Meteor.</i>,
<b>14</b>, 1156-1165.</li>
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
<td valign="top">initial allocation failed for radial vel obs data,
itemcount = (max_radial_vel_obs)</td>
<!-- comment -->
<td valign="top">Need to increase max_radial_vel_obs count in
namelist</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_radial_vel</td>
<!-- message -->
<td valign="top">Expected location header "platform" in input
file</td>
<!-- comment -->
<td valign="top">The format of the input file is not
consistent.</td>
</tr>
<tr><!-- routine -->
<td valign="top">velkey_out_of_range</td>
<!-- message -->
<td valign="top">velkey (val) exceeds max_radial_vel_obs
(maxval)</td>
<!-- comment -->
<td valign="top">The number of radial velocity observations exceeds
the array size allocated in the module. Need to increase
max_radial_vel_obs count in namelist.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_nyquist_velocity</td>
<!-- message -->
<td valign="top">bad value for nyquist velocity</td>
<!-- comment -->
<td valign="top">The format of the input obs_seq file is not
consistent.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_beam_direction</td>
<!-- message -->
<td valign="top">beam_direction value must be between -1 and 1, got
()</td>
<!-- comment -->
<td valign="top">The format of the input obs_seq file is not
consistent.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_beam_direction</td>
<!-- message -->
<td valign="top">Expected orientation header "dir3d" in input
file</td>
<!-- comment -->
<td valign="top">The format of the input obs_seq file is not
consistent.</td>
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
<table>
<tr>
<td><em class="call">use obs_def_radar_mod, only :</em></td>
<td><a href="#initialize_module">initialize_module</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_beam_direction">read_beam_direction</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_nyquist_velocity">read_nyquist_velocity</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_beam_direction">write_beam_direction</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#write_nyquist_velocity">write_nyquist_velocity</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#interactive_beam_direction">interactive_beam_direction</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#interactive_nyquist_velocity">interactive_nyquist_velocity</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_reflectivity">get_reflectivity</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_precip_fall_speed">get_precip_fall_speed</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#initialize_constants">initialize_constants</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#print_constants">print_constants</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#pr_con">pr_con</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#velkey_out_of_range">velkey_out_of_range</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#check_namelist_limits">check_namelist_limits</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ascii_file_format">ascii_file_format</a></td>
</tr>
</table>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="initialize_module" id="initialize_module"></a><br>
<div class="routine"><em class="call">call
initialize_module()</em></div>
<div class="indent1"><!-- Description -->
<p>Reads the namelist, allocates space for the auxiliary data
associated wtih radial velocity observations, initializes the
constants used in subsequent computations (possibly altered by
values in the namelist), and prints out the list of constants and
the values in use. These may need to change depending on which
microphysics scheme is being used.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_beam_direction" id="read_beam_direction"></a><br>
<div class="routine"><em class="call">beam_direction =
read_beam_direction(ifile, is_asciiformat)</em>
<pre>
real(r8), dimension(3)            :: <em class=
"code">read_beam_direction</em>
integer,               intent(in) :: <em class="code">ifile</em>
logical,               intent(in) :: <em class=
"code">is_asciiformat</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads the beam direction at the observation location. Auxiliary
data for doppler radial velocity observations.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">read_beam_direction  </em></td>
<td>Returns three real values for the radar beam orientation</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>File unit descriptor for input file</td>
</tr>
<tr>
<td valign="top"><em class="code">is_asciiformat</em></td>
<td>File format specifier: .TRUE. if file is formatted/ascii, or
.FALSE. if unformatted/binary. Default .TRUE.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_nyquist_velocity" id=
"read_nyquist_velocity"></a><br>
<div class="routine"><em class="call">nyquist_velocity =
read_nyquist_velocity(ifile, is_asciiformat)</em>
<pre>
real(r8),            :: <em class="code">read_nyquist_velocity</em>
integer,  intent(in) :: <em class="code">ifile</em>
logical,  intent(in) :: <em class="code">is_asciiformat</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads nyquist velocity for a doppler radial velocity
observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">read_nyquist_velocity  </em></td>
<td>Returns a real value for the nyquist velocity value</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>File unit descriptor for input file</td>
</tr>
<tr>
<td valign="top"><em class="code">is_asciiformat</em></td>
<td>File format specifier: .TRUE. if file is formatted/ascii, or
.FALSE. if unformatted/binary. Default .TRUE.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_beam_direction" id="write_beam_direction"></a><br>
<div class="routine"><em class="call">call
write_beam_direction(ifile, beam_direction, is_asciiformat)</em>
<pre>
integer,                intent(in) :: <em class="code">ifile</em>
real(r8), dimension(3), intent(in) :: <em class=
"code">beam_direction</em>
logical,                intent(in) :: <em class=
"code">is_asciiformat</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes the beam direction at the observation location. Auxiliary
data for doppler radial velocity observations.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>File unit descriptor for output file</td>
</tr>
<tr>
<td valign="top"><em class=
"code">beam_direction  </em></td>
<td>Three components of the radar beam orientation</td>
</tr>
<tr>
<td valign="top"><em class="code">is_asciiformat</em></td>
<td>File format specifier: .TRUE. if file is formatted/ascii, or
.FALSE. if unformatted/binary. Default .TRUE.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_nyquist_velocity" id=
"write_nyquist_velocity"></a><br>
<div class="routine"><em class="call">call
write_nyquist_velocity(ifile, nyquist_velocity,
is_asciiformat)</em>
<pre>
integer,  intent(in) :: <em class="code">ifile</em>
real(r8), intent(in) :: <em class="code">nyquist_velocity</em>
logical,  intent(in) :: <em class="code">is_asciiformat</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes nyquist velocity for a doppler radial velocity
observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ifile</em></td>
<td>File unit descriptor for output file</td>
</tr>
<tr>
<td valign="top"><em class=
"code">nyquist_velocity  </em></td>
<td>The nyquist velocity value for this observation</td>
</tr>
<tr>
<td valign="top"><em class="code">is_asciiformat</em></td>
<td>File format specifier: .TRUE. if file is formatted/ascii, or
.FALSE. if unformatted/binary. Default .TRUE.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="interactive_beam_direction" id=
"interactive_beam_direction"></a><br>
<div class="routine"><em class="call">call
interactive_beam_direction(beam_direction)</em>
<pre>
real(r8), dimension(3), intent(out) :: <em class=
"code">beam_direction</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Prompts the user for input for the azimuth and elevation of the
radar beam at the observation location. Will be converted to the
three values actually stored in the observation sequence file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">beam_direction  </em></td>
<td>Three components of the radar beam orientation</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="interactive_nyquist_velocity" id=
"interactive_nyquist_velocity"></a><br>
<div class="routine"><em class="call">call
interactive_nyquist_velocity(nyquist_velocity)</em>
<pre>
real(r8), intent(out) :: <em class="code">nyquist_velocity</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Prompts the user for input for the nyquist velocity value
associated with a doppler radial velocity observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">nyquist_velocity  </em></td>
<td>Nyquist velocity value for the observation.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_reflectivity" id="get_reflectivity"></a><br>
<div class="routine"><em class="call">call get_reflectivity(qr, qg,
qs, rho, temp, ref)</em>
<pre>
real(r8), intent(in)  :: <em class="code">qr</em>
real(r8), intent(in)  :: <em class="code">qg</em>
real(r8), intent(in)  :: <em class="code">qs</em>
real(r8), intent(in)  :: <em class="code">rho</em>
real(r8), intent(in)  :: <em class="code">temp</em>
real(r8), intent(out) :: <em class="code">ref</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes the equivalent radar reflectivity factor in
mm<sup>6</sup> m<sup>-3</sup> for simple single-moment microphysics
schemes such as Kessler and Lin, et al. See the <a href=
"#References">references</a> for more details.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">qr</em></td>
<td>Rain water content (kg kg<sup>-1</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">qg</em></td>
<td>Graupel/hail content (kg kg<sup>-1</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">qs</em></td>
<td>Snow content (kg kg<sup>-1</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">rho</em></td>
<td>Air density (kg m<sup>-3</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">temp  </em></td>
<td>Air temperature (K)</td>
</tr>
<tr>
<td valign="top"><em class="code">ref</em></td>
<td>The returned radar reflectivity value</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_precip_fall_speed" id=
"get_precip_fall_speed"></a><br>
<div class="routine"><em class="call">call
get_precip_fall_speed(qr, qg, qs, rho, temp,
precip_fall_speed)</em>
<pre>
real(r8), intent(in)  :: <em class="code">qr</em>
real(r8), intent(in)  :: <em class="code">qg</em>
real(r8), intent(in)  :: <em class="code">qs</em>
real(r8), intent(in)  :: <em class="code">rho</em>
real(r8), intent(in)  :: <em class="code">temp</em>
real(r8), intent(out) :: <em class="code">precip_fall_speed</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes power-weighted precipitation fall speed in m
s<sup>-1</sup> for simple single-moment microphysics schemes such
as Kessler and Lin, et al. See the <a href=
"#References">references</a> for more details.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">qr</em></td>
<td>Rain water content (kg kg<sup>-1</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">qg</em></td>
<td>Graupel/hail content (kg kg<sup>-1</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">qs</em></td>
<td>Snow content (kg kg<sup>-1</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">rho</em></td>
<td>Air density (kg m<sup>-3</sup>)</td>
</tr>
<tr>
<td valign="top"><em class="code">temp</em></td>
<td>Air temperature (K)</td>
</tr>
<tr>
<td valign="top"><em class=
"code">precip_fall_speed  </em></td>
<td>The returned precipitation vall speed</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="initialize_constants" id="initialize_constants"></a><br>
<div class="routine"><em class="call">call
initialize_constants()</em></div>
<div class="indent1"><!-- Description -->
<p>Set values for a collection of constants used throughout the
module during the various calculations. These are set once in this
routine and are unchanged throughout the rest of the execution.
They cannot be true Fortran <em class="code">parameters</em>
because some of the values can be overwritten by namelist entries,
but once they are set they are treated as read-only parameters.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="print_constants" id="print_constants"></a><br>
<div class="routine"><em class="call">call
print_constants()</em></div>
<div class="indent1"><!-- Description -->
<p>Print out the names and values of all constant parameters used
by this module. The error handler message facility is used to print
the message, which by default goes to both the DART log file and to
the standard output of the program.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="pr_con" id="pr_con"></a><br>
<div class="routine"><em class="call">call pr_con(c_val,
c_str)</em>
<pre>
real(r8),         intent(in)  :: <em class="code">c_val</em>
character(len=*), intent(in)  :: <em class="code">c_str</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Calls the DART error handler routine to print out a string label
and a real value to both the log file and to the standard
output.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">Value of
constant  </em></td>
<td>A real value.</td>
</tr>
<tr>
<td valign="top"><em class="code">Name of constant</em></td>
<td>A character string.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="velkey_out_of_range" id="velkey_out_of_range"></a><br>
<div class="routine"><em class="call">call
velkey_out_of_range(velkey)</em>
<pre>
integer, intent(in)  :: <em class="code">velkey</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Range check key and trigger a fatal error if larger than the
allocated array for observation auxiliary data.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">velkey  </em></td>
<td>Integer key into a local array of auxiliary observation
data.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="check_namelist_limits" id=
"check_namelist_limits"></a><br>
<div class="routine"><em class="call">call
check_namelist_limits(apply_ref_limit_to_obs,
reflectivity_limit_obs, lowest_reflectivity_obs,
apply_ref_limit_to_fwd_op, reflectivity_limit_fwd_op,
lowest_reflectivity_fwd_op)</em>
<pre>
logical,  intent(in) :: <em class=
"code">apply_ref_limit_to_obs</em>
real(r8), intent(in) :: <em class=
"code">reflectivity_limit_obs</em>
real(r8), intent(in) :: <em class=
"code">lowest_reflectivity_obs</em>
logical,  intent(in) :: <em class=
"code">apply_ref_limit_to_fwd_op</em>
real(r8), intent(in) :: <em class=
"code">reflectivity_limit_fwd_op</em>
real(r8), intent(in) :: <em class=
"code">lowest_reflectivity_fwd_op</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Check the values set in the namelist for consistency. Print out
a message if the limits and set-to values are different; this may
be intentional but is not generally expected to be the case. In all
cases below, see the namelist documentation for a fuller
explanation of each value.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">apply_ref_limit_to_obs</em></td>
<td>Logical. See <a href="#Namelist">namelist</a>.</td>
</tr>
<tr>
<td valign="top"><em class="code">reflectivity_limit_obs</em></td>
<td>Real value. See <a href="#Namelist">namelist</a>.</td>
</tr>
<tr>
<td valign="top"><em class="code">lowest_reflectivity_obs</em></td>
<td>Real value. See <a href="#Namelist">namelist</a>.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">apply_ref_limit_to_fwd_op</em></td>
<td>Logical. See <a href="#Namelist">namelist</a>.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">reflectivity_limit_fwd_op</em></td>
<td>Real value. See <a href="#Namelist">namelist</a>.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">lowest_reflectivity_fwd_op  </em></td>
<td>Real value. See <a href="#Namelist">namelist</a>.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="ascii_file_format" id="ascii_file_format"></a><br>
<div class="routine"><em class="call">is_asciifile =
ascii_file_format(fform)</em>
<pre>
logical                                :: <em class=
"code">ascii_file_format</em>
character(len=*), intent(in), optional :: <em class=
"code">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Should be moved to DART utility module at some point. Returns
.TRUE. if the optional argument is missing or if it is not one of
the following values: <em class="code">"unformatted",
"UNFORMATTED", "unf", "UNF"</em>.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">ascii_file_format  </em></td>
<td>Return value. Logical. Default is .TRUE.</td>
</tr>
<tr>
<td valign="top"><em class="code">fform</em></td>
<td>Character string file format.</td>
</tr>
</table>
</div>
<br>
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
