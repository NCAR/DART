<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>program dart_to_ncommas</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">dart_to_ncommas</em></h1>
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
<p><em class="program">dart_to_ncommas</em> is the program that
<strong>updates</strong> a ncommas netCDF-format restart file
(usually <em class="file">ncommas_restart.nc</em>) with the state
information contained in a DART output/restart file (e.g.
<em class="file">perfect_ics, filter_ics, ...</em> ). Only the
CURRENT values in the ncommas restart file will be updated. The
DART model time is compared to the time in the ncommas restart
file. If the last time in the restart file does not match the DART
model time, the program issues an error message and aborts.<br />
<br />
From the user perspective, most of the time <em class=
"program">dart_to_ncommas</em> will be used on DART files that have
a header containing one time stamp followed by the model
state.<br />
<br />
The <a href="#Namelist">dart_to_ncommas_nml</a> namelist allows
<em class="program">dart_to_ncommas</em> to read the <em class=
"file">assim_model_state_ic</em> files that have <em class=
"italic">two</em> timestamps in the header. These files are
temporarily generated when DART is used to advance the model. One
timestamp is the 'advance_to' time, the other is the 'valid_time'
of the model state. In this case, a namelist for ncommas (called
<em class="file">ncommas_in.DART</em>) is written that contains the
<em class="code">&amp;time_manager_nml</em> settings appropriate to
advance ncommas to the time requested by DART. The repository
version of the <em class="program">advance_model.csh</em> script
has a section to ensure the proper DART namelist settings for this
case.<br />
<br />
Conditions required for successful execution of <em class=
"program">dart_to_ncommas</em>:</p>
<ul>
<li>a valid <em class="file">input.nml</em> namelist file for
DART</li>
<li>a valid <em class="file">ncommas_vars.nml</em> namelist file
for ncommas - the same one used to create the DART state vector,
naturally,</li>
<li>a DART file (typically <em class=
"file">filter_restart.xxxx</em> or <em class=
"file">filter_ics.xxxx</em>)</li>
<li>a ncommas restart file (typically <em class=
"file">ncommas_restart.nc</em>).</li>
</ul>
<p>Since this program is called repeatedly for every ensemble
member, we have found it convenient to link the DART input file to
the default input filename (<em class="file">dart_restart</em>).
The same thing goes true for the ncommas output filename <em class=
"file">ncommas_restart.nc</em>.</p>
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
&amp;dart_to_ncommas_nml
   dart_to_ncommas_input_file = 'dart_restart',
   advance_time_present   = .false.  
/
</pre></div>
<br />
<br />
<p><em class="code">dart_to_ncommas_nml</em> and <em class=
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
<td>dart_to_ncommas_input_file</td>
<td>character(len=128)</td>
<td>The name of the DART file containing the model state to insert
into the ncommas restart file.</td>
</tr>
<tr>
<td>advance_time_present</td>
<td>logical</td>
<td>If you are converting a DART initial conditions or restart file
this should be <em class="code">.false.</em>; these files have a
single timestamp describing the valid time of the model state. If
<em class="code">.true.</em> TWO timestamps are expected to be the
DART file header. In this case, a namelist for ncommas (called
<em class="file">ncommas_in.DART</em>) is created that contains the
<em class="code">&amp;time_manager_nml</em> settings appropriate to
advance ncommas to the time requested by DART.</td>
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
<li>DART initial conditions/restart file; e.g. <em class=
"file">filter_ic</em></li>
<li>DART namelist file; <em class="file">input.nml</em></li>
<li>ncommas namelist file; <em class=
"file">ncommas_vars.nml</em></li>
<li>ncommas restart file <em class=
"file">ncommas_restart.nc</em></li>
</ul>
<h2>FILES Written</h2>
<ul>
<li>ncommas restart file; <em class=
"file">ncommas_restart.nc</em></li>
<li>ncommas namelist file; <em class=
"file">ncommas_in.DART</em></li>
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
