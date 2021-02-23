<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program trans_pv_sv</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<center><a href="#Modules">MODULES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a></center>
<h1>PROGRAM <em class="program">trans_pv_sv</em></h1>
<!-- version tag follows, do not edit -->
<p></p>
<p><em class="program">trans_pv_sv</em> is responsible for
converting the ocean model 'snapshot' files to a DART 'initial
conditions' file. In order to do that, the valid time for the
snapshot files must be calculated from several pieces of
information: the filename contains a timestep index, the <em class=
"file">data</em><em class="unix">&amp;PARM03</em> namelist contains
information about the amount of time per timestep, and the
<em class="file">data.cal</em><em class="unix">&amp;CAL_NML</em>
namelist contains the start time. Additionally, the grid
characteristics must be read from <em class=
"file">data</em><em class="unix">&amp;PARM04</em>. Consequently,
the files <em class="file">data</em>, and <em class=
"file">data.cal</em> as well as the general <em class=
"file">input.nml</em> are needed in addition to the snapshot
files.<br>
<br>
This program has a number of options that are driven from namelists
and <strong>one</strong> piece of input read from STDIN: the
integer representing the timestep index of the snapshot file
set.</p>
<h2>Usage</h2>
<p>The output filename is hardwired to that expected by <em class=
"program">filter</em>. This example creates an output file named
<em class="file">assim_model_state_ud</em> from the following files
in the local directory:<br>
<br>
<em class="file">S.0000000096.data</em><br>
<em class="file">T.0000000096.data</em><br>
<em class="file">U.0000000096.data</em><br>
<em class="file">V.0000000096.data</em><br>
<em class="file">Eta.0000000096.data</em></p>
<div class="unix">./trans_pv_sv &lt; 96</div>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
model_mod
assim_model_mod
time_manager_mod
</pre>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>This program has no namelist of its own, but some of the
underlying modules require namelists. To avoid duplication and,
possibly, some inconsistency in the documentation, only a list of
the required namelists is provided here, with a hyperlink to the
full documentation for each namelist.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Namelist</th>
<th align="left">Primary Purpose</th>
</tr>
<tr>
<td><a href=
"../../assimilation_code/modules/utilities/utilities_mod.html#Namelist">
utilities_nml</a></td>
<td>set the termination level and file name for the run-time
log</td>
</tr>
<tr>
<td><a href=
"../../assimilation_code/modules/assimilation/assim_model_mod.html#Namelist">
assim_model_mod_nml</a></td>
<td>write DART restart files in binary or ASCII</td>
</tr>
<tr>
<td><a href="model_mod.html#Namelist">model_nml</a></td>
<td>write netCDF files with prognostic variables</td>
</tr>
<tr>
<td><a href="model_mod.html#namelist_cal_nml">CAL_NML</a></td>
<td>determine start time of the ocean model</td>
</tr>
<tr>
<td><a href="model_mod.html#namelist_parm03">PARM03</a></td>
<td>the amount of time per model timestep for deciphering snapshot
filenames</td>
</tr>
<tr>
<td><a href="model_mod.html#namelist_parm04">PARM04</a></td>
<td>ocean model grid parameters</td>
</tr>
</table>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<hr>
<h2>FILES</h2>
<ul>
<li>input namelist files: <em class="file">data, data.cal,
input.nml</em></li>
<li>input snapshot files: <em class=
"file">[S,T,U,V,Eta].nnnnnnnnnn.[data[,.meta]]</em></li>
<li>output initial conditions file: <em class=
"file">assim_model_state_ud</em></li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>none</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<p>The most common problem is trying to read the Fortran
direct-access big-endian snapshot files on a little-endian
architecture. This can manifest itself in very misleading ways.
Make sure you have the right compiler settings to be able to read
these files. There is no one error message that indicates the read
was unsuccessful.<br>
<br>
The read takes place in <a href=
"model_mod.html#read_snapshot">model_mod:read_snapshot()</a>.</p>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">trans_sv_pv</td>
<!-- message -->
<td valign="top">unable to read timestep from stdin.</td>
<!-- comment -->
<td valign="top">look at the example in the 'Usage' section.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>There are no known bugs.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<p>None at this time. Feel free to suggest improvements.</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
