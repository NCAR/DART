<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>program cam_to_dart</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">cam_to_dart</em></h1>
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
<h2>Overview</h2>
<p><em class="program">cam_to_dart</em> is the program that reads a
CAM restart file (usually <em class="file">caminput.nc</em>) and
creates a single DART output/restart file (e.g. <em class=
"file">perfect_ics, filter_ics, ...</em> ). If you have multiple
input files, you will need to rename the output files as you create
them.<br />
<br />
The list of variables extracted from the CAM netCDF file and
conveyed to DART is controlled by the set of <em class=
"file">input.nml</em> <em class=
"code">&amp;model_nml:state_names_*</em> variables. The <em class=
"code">date</em> and <em class="code">datesec</em> variables in the
CAM netcdf file are used to specify the valid time of the state
vector. The time may be changed with the <a href=
"../../assimilation_code/programs/restart_file_tool/restart_file_tool.html">
restart_file_tool</a> if desired.<br />
<br />
Some CAM restart files are from climatological runs and have a
valid time that predates the use of the Gregorian calendar. In such
instances, the year component of the original date is changed to be
a valid Gregorian year (by adding 1601). A warning is issued to the
screen and to the logfile. Please use the <a href=
"../../assimilation_code/programs/restart_file_tool/restart_file_tool.html">
restart_file_tool</a> to change this time.<br />
<br />
Conditions required for successful execution of <em class=
"program">cam_to_dart</em>:</p>
<ul>
<li>a valid <em class="file">input.nml</em> namelist file for
DART</li>
<li>a CAM 'phis' netCDF file [default: <em class=
"file">cam_phis.nc</em>]</li>
<li>a CAM restart file [default: <em class=
"file">caminput.nc</em>].</li>
</ul>
<p>Since this program is called repeatedly for every ensemble
member, we have found it convenient to link the CAM restart files
to the default input filename (<em class="file">caminput.nc</em>).
The default DART output filename is <em class="file">dart_ics</em>
- this may be moved or linked as necessary.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<!--==================================================================-->
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
&amp;cam_to_dart_nml
   cam_to_dart_input_file  = 'caminput.nc',
   cam_to_dart_output_file = 'dart_ics', 
   /
</pre></div>
<br />
<br />
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
<td>cam_to_dart_input_file</td>
<td>character(len=128)</td>
<td>The name of the DART file containing the CAM state.</td>
</tr>
<tr>
<td>cam_to_dart_output_file</td>
<td>character(len=128)</td>
<td>The name of the DART file containing the model state derived
from the CAM restart file.</td>
</tr>
</tbody>
</table>
</div>
<br />
<br />
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>MODULES USED</h2>
<pre>
assim_model_mod.f90
types_mod.f90
threed_sphere/location_mod.f90
model_mod.f90
null_mpi_utilities_mod.f90
obs_kind_mod.f90
random_seq_mod.f90
time_manager_mod.f90
utilities_mod.f90
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FILES Read</h2>
<ul>
<li>DART namelist file; <em class="file">input.nml</em></li>
<li>CAM restart file; <em class="file">caminput.nc</em></li>
<li>CAM "phis" file specified in <em class=
"code">&amp;model_nml::cam_phis</em> (normally <em class=
"file">cam_phis.nc</em>)</li>
</ul>
<h2>FILES Written</h2>
<ul>
<li>DART initial conditions/restart file; e.g. <em class=
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
