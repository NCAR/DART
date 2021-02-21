<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program compute_error</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">compute_error</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>Utility program to compute the time-mean ensemble error and
spread in the same manner that the DART MATLAB diagnostic routine
'plot_total_err' does. It runs from the command line, opens no
windows, and outputs several types of numerical results on standard
output. Grep for 'Total' to get the 2 lines with total error and
total spread. Intended for scripts where only the numeric results
are wanted instead of a time-series plot. This routine does not do
any weighted computations.</p>
<p>The default is to compare a True_State.nc file output from
perfect_model_obs to a Prior_Diag.nc file output from filter. Other
filenames can be specified in the namelist. These files must have
at least one overlapping value in the 'time' array. The statistics
will be done on the overlapping time region only.</p>
<p>The output includes the min and max error and spread values, and
the time index and time value where that occurs. There is also an
option to recompute the time mean ensemble error and spread after
skipping the first N times. This can be useful to skip an initial
error spike while the model is spinning up which can result in a
larger than expected total error.</p>
<p>Namelist interface <a href="#Namelist"><em class=
"code">&amp;compute_error_nml</em></a> is read from file <em class=
"file">input.nml</em>.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
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
&amp;compute_error_nml
   truth_file_name   = 'true_state.nc'
   diag_file_name    = 'preassim.nc'
   skip_first_ntimes = 0
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
<td>truth_file_name</td>
<td>character(len=256)</td>
<td>State-space diagnostic file from the 'perfect_model_obs'
program.</td>
</tr>
<tr>
<td>diag_file_name</td>
<td>character(len=256)</td>
<td>State space diagnostic file output from the 'filter'
program.</td>
</tr>
<tr>
<td>skip_first_ntimes</td>
<td>integer</td>
<td>If set to a value greater than 0, the error values will be
recomputed a second time, skipping the first N times. This can be
useful when running an experiment that has an initial error spike
as the model spins up and then decays down to a more steady
state.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>DART diagnosic files (True_State.nc, Prior_Diag.nc)</li>
<li>compute_error.nml</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>none</li>
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
<td valign="top">time dimension error</td>
<!-- message -->
<td valign="top">files must have overlapping time series</td>
<!-- comment -->
<td valign="top">The unlimited 'time' dimension must values in
common between both files.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>The matlab script has an option for doing weighted statistics.
This code does not.</p>
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
