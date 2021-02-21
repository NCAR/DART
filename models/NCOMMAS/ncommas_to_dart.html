<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>program ncommas_to_dart</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">ncommas_to_dart</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70" /></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<p><em class="program">ncommas_to_dart</em> is the program that
reads a ncommas restart file (usually <em class=
"file">ncommas_restart.nc</em>) and creates a DART state vector
file (e.g. <em class="file">perfect_ics, filter_ics, ...</em>
).<br />
<br />
The list of variables used to create the DART state vector are
specified in the <em class="file">ncommas_vars.nml</em> file.<br />
<br />
Conditions required for successful execution of <em class=
"program">ncommas_to_dart</em>:</p>
<ul>
<li>a valid <em class="file">input.nml</em> namelist file for
DART</li>
<li>a valid <em class="file">ncommas_vars.nml</em> namelist file
for ncommas</li>
<li>the ncommas restart file mentioned in the <em class=
"file">input.nml&amp;model_nml:ncommas_restart_filename</em>
variable.</li>
</ul>
<p>Since this program is called repeatedly for every ensemble
member, we have found it convenient to link the ncommas restart
files to the default input filename (<em class=
"file">ncommas_restart.nc</em>). The default DART state vector
filename is <em class="file">dart_ics</em> - this may be moved or
linked as necessary.</p>
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;model_nml
   ncommas_restart_filename     = 'ncommas_restart.nc';
   assimilation_period_days     = 1,
   assimilation_period_seconds  = 0,
   model_perturbation_amplitude = 0.2,
   output_state_vector          = .true.,
   calendar                     = 'Gregorian',
   debug                        = 0
/
</pre></div>
<br />
<div class="namelist">
<pre>
&amp;ncommas_to_dart_nml
   ncommas_to_dart_output_file = 'dart_ics'  
/
</pre></div>
<br />
<br />
<p><em class="code">ncommas_to_dart_nml</em> and <em class=
"code">model_nml</em> are always read from a file called <em class=
"file">input.nml</em>. The full description of the <em class=
"code">model_nml</em> namelist is documented in the <a href=
"model_mod.html#Namelist">NCOMMAS model_mod</a>.</p>
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
<td>ncommas_to_dart_output_file</td>
<td>character(len=128)</td>
<td>The name of the DART file which contains the updated model
state info that should be written into the NCOMMAS file.</td>
</tr>
</tbody>
</table>
</div>
<br />
<br />
<p><em class="code">ncommas_vars_nml</em> is always read from a
file called <em class="file">ncommas_vars.nml</em>.</p>
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
<td>ncommas_state_variables</td>
<td>character(len=NF90_MAX_NAME) ::<br />
dimension(160)</td>
<td>The list of variable names in the NCOMMAS restart file to use
to create the DART state vector and their corresponding DART
kind.</td>
</tr>
</tbody>
</table>
</div>
<br />
<div class="namelist">
<pre>
&amp;ncommas_vars_nml
   ncommas_state_variables = 'U',   'QTY_U_WIND_COMPONENT',
                             'V',   'QTY_V_WIND_COMPONENT',
                             'W',   'QTY_VERTICAL_VELOCITY',
                             'TH',  'QTY_POTENTIAL_TEMPERATURE',
                             'DBZ', 'QTY_RADAR_REFLECTIVITY',
                             'WZ',  'QTY_VERTICAL_VORTICITY',
                             'PI',  'QTY_EXNER_FUNCTION',
                             'QV',  'QTY_VAPOR_MIXING_RATIO',
                             'QC',  'QTY_CLOUDWATER_MIXING_RATIO',
                             'QR',  'QTY_RAINWATER_MIXING_RATIO',
                             'QI',  'QTY_ICE_MIXING_RATIO',
                             'QS',  'QTY_SNOW_MIXING_RATIO',
                             'QH',  'QTY_GRAUPEL_MIXING_RATIO'
  /
</pre></div>
<br />
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>MODULES USED</h2>
<pre>
assim_model_mod
location_mod
model_mod
null_mpi_utilities_mod
obs_kind_mod
random_seq_mod
time_manager_mod
types_mod
utilities_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FILES Read</h2>
<ul>
<li>ncommas restart file; <em class=
"file">ncommas_restart.nc</em></li>
<li>DART namelist files; <em class="file">input.nml</em> and
<em class="file">ncommas_vars.nml</em></li>
</ul>
<h2>FILES Written</h2>
<ul>
<li>DART state vector file; e.g. <em class=
"file">dart_ics</em></li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>REFERENCES</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>ERROR CODES and CONDITIONS</h2>
<p>none - all error messages come from modules that have their own
documentation.</p>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FUTURE PLANS</h2>
<p>None.</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
