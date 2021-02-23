<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module cov_cutoff_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE cov_cutoff_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#Interface">INTERFACES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Computes the weight with which an observation should impact a
state variable that is separated by a given distance. The distance
is in units determined by the location module being used.</p>
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
&amp;cov_cutoff_nml
   select_localization = 1  
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
<td>select_localization</td>
<td>integer</td>
<td>Selects the localization function.
<ul style="list-style: none;">
<li>1 = Gaspari-Cohn 5th order polynomial with halfwidth c.</li>
<li>2 = Boxcar with halfwidth c (goes to 0 for z_in &lt; 2c).</li>
<li>3 = Ramped Boxcar. Has value 1 for z_in &lt; c and then reduces
linearly to 0 at z_in = 2c.</li>
</ul>
</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--==================================================================-->
 <a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
location_mod
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
<table summary='modules used'>
<tr>
<td><em class="call">use cov_factor_mod, only :</em></td>
<td><a href="#comp_cov_factor">comp_cov_factor</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE ========================-->
<a name="comp_cov_factor" id="comp_cov_factor"></a><br>
<div class="routine"><em class="call">var = comp_cov_factor(z_in, c
<em class="optionalcode">[, obs_loc] [, obs_type]
[, target_loc] [, target_kind]
[, localization_override]</em>)</em>
<pre>
real(r8)                                  :: <em class=
"code">comp_cov_factor</em>
real(r8), intent(in)                      :: <em class=
"code">z_in</em>
real(r8), intent(in)                      :: <em class=
"code">c</em>
type(location_type), optional, intent(in) :: <em class=
"optionalcode">obs_loc</em>
integer, optional, intent(in)             :: <em class=
"optionalcode">obs_type</em>
type(location_type), optional, intent(in) :: <em class=
"optionalcode">target_loc</em>
integer, optional, intent(in)             :: <em class=
"optionalcode">target_kind</em>
integer, optional, intent(in)             :: <em class=
"optionalcode">localization_override</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a weighting factor for observation and a target variable
(state or observation) separated by distance z_in and with a
half-width distance, c. Three options are provided and controlled
by a namelist parameter. The optional argument
localization_override controls the type of localization function if
present. The optional arguments obs_loc, obs_type and target_loc,
target_kind are not used in the default code. They are made
available for users who may want to design more sophisticated
localization functions.</p>
<table width="100%" border="0" summary="function arguments"
cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Weighting factor.</td>
</tr>
<tr>
<td valign="top"><em class="code">z_in</em></td>
<td>The distance between observation and target.</td>
</tr>
<tr>
<td valign="top"><em class="code">c</em></td>
<td>Factor that describes localization function. Describes
half-width of functions used here.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">obs_loc</em></td>
<td>Location of the observation.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">obs_type</em></td>
<td>Observation specific type.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">target_loc</em></td>
<td>Location of target.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">target_kind</em></td>
<td>Generic kind of target.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">localization_override</em></td>
<td>Controls localization type if present. Same values as for
namelist control.</td>
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
<table border="0" summary="files used">
<tr>
<th>filename</th>
<th>purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read <em class="code">cov_cutoff_nml</em></td>
</tr>
</table>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ol>
<li>Gaspari and Cohn, 1999, QJRMS, <b>125</b>, 723-757. (eqn.
4.10)</li>
</ol>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%"
summary="errors">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">comp_cov_factor</td>
<!-- message -->
<td valign="top">Illegal value of "select_localization" in
cov_cutoff_mod namelist</td>
<!-- comment -->
<td valign="top">Only values 1 through 3 select a localization
function.</td>
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
