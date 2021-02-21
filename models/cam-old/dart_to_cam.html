<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>program dart_to_cam</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">dart_to_cam</em></h1>
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
<p><em class="program">dart_to_cam</em> is the program that reads a
DART restart or model advance file (e.g. <em class=
"file">perfect_ics, filter_ics, assim_model_state_id ...</em> ).
and overwrites the part of the CAM data in a single CAM restart
file (usually <em class="file">caminput.nc</em>) which is in the
DART state vector. If you have multiple input files, you will need
to rename the output files as you create them.</p>
<p>The list of variables extracted from the DART state vector and
exported to the CAM netCDF file is controlled by the set of
<em class="file">input.nml</em> <em class=
"code">&amp;model_nml:state_names_*</em> variables.</p>
<p>If the input is a model advance file, containing 2 timestamps
(the current model time and a future time the model should run
until), this program also creates a separate file named <em class=
"file">times</em> that contains three lines: the advance-to time,
the current model time, and the number of hours to advance. These
will need to be extracted and inserted in a CAM namelist to
indicate to CAM how long to run.</p>
<p>This program also updates the <em class="code">date</em> and
<em class="code">datesec</em> variables in the CAM netcdf file.
Generally these are identical times since the assimilation doesn't
change the time of the data, but in case the original file had a
different time that was overwritten in the state vector, it will
update the time for consistency.</p>
<p>Conditions required for successful execution of <em class=
"program">dart_to_cam</em>:</p>
<ul>
<li>a valid <em class="file">input.nml</em> namelist file for
DART</li>
<li>a CAM 'phis' netCDF file [default: <em class=
"file">cam_phis.nc</em>]</li>
<li>a DART restart file [default: <em class="file">dart_ics</em>]
(read)</li>
<li>a CAM restart file [default: <em class="file">caminput.nc</em>]
(read and written)</li>
</ul>
<p>Since this program is called repeatedly for every ensemble
member, we have found it convenient to link the DART input and CAM
restart files to the default filenames <em class=
"file">dart_ics</em> and <em class="file">caminput.nc</em>). The
output files may be moved or relinked as necessary.</p>
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
&amp;dart_to_cam_nml
   dart_to_cam_input_file  = 'dart_ics',
   dart_to_cam_output_file = 'caminput.nc',
   advance_time_present    = .true.,
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
<td>dart_to_cam_input_file</td>
<td>character(len=128)</td>
<td>The name of the DART restart file containing the CAM
state.</td>
</tr>
<tr>
<td>dart_to_cam_output_file</td>
<td>character(len=128)</td>
<td>The name of the CAM restart netcdf file.</td>
</tr>
<tr>
<td>advance_time_present</td>
<td>logical</td>
<td>Set to .false. for DART initial condition and restart files.
Use the .true. setting for the files written by filter during a
model advance.</td>
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
<li>DART initial conditions/restart file; e.g. <em class=
"file">dart_ics</em> (read)</li>
<li>CAM restart file; <em class="file">caminput.nc</em> (read and
written)</li>
<li>CAM "phis" file specified in <em class=
"code">&amp;model_nml::cam_phis</em> (normally <em class=
"file">cam_phis.nc</em>)</li>
</ul>
<h2>FILES Written</h2>
<ul>
<li>CAM restart file; <em class="file">caminput.nc</em> (read and
written)</li>
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
<a name="Bugs" id="Bugs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
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
