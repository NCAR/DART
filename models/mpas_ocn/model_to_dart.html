<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>program model_to_dart (MPAS OCN)</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">model_to_dart</em> for MPAS
OCN</h1>
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
<p><em class="program">model_to_dart</em> is the program that reads
an MPAS OCN analysis file (nominally named <em class=
"file">mpas_restart.nc</em>) and creates a DART state vector file
(e.g. <em class="file">perfect_ics, filter_ics, ...</em> ). The
MPAS analysis files have a <strong>Time</strong> UNLIMITED
Dimension, which indicates there may (at some point) be more than
one timestep in the file. The DART routines are currently designed
to use the LAST timestep. If the Time dimension of length 3, we use
the third timestep. A warning message is issued and indicates
exactly the time being used.<br />
<br />
<em class="file">input.nml</em><em class=
"code">&amp;mpas_vars_nml</em> defines the list of MPAS variables
used to build the DART state vector. This namelist is more fully
described in the <a href="model_mod.html">MPAS model_mod.html</a>
documentation. For example:</p>
<pre>&amp;mpas_vars_nml
   mpas_state_variables = 'temperature',  'QTY_TEMPERATURE',
                          'salinity',     'QTY_SALINITY',
                          'rho',          'QTY_DENSITY',
                          'u',            'QTY_EDGE_NORMAL_SPEED',
                          'h',            'QTY_SEA_SURFACE_HEIGHT'
                          'tracer1',      'QTY_TRACER_CONCENTRATION'
   /
</pre>
<p>Conditions required for successful execution of <em class=
"program">model_to_dart</em> are:</p>
<ul>
<li>a valid <em class="file">input.nml</em> namelist file for DART
which contains</li>
<li>a MPAS OCN analysis file (nominally named <em class=
"file">mpas_analysis.nc</em>).</li>
</ul>
<p>Since this program is called repeatedly for every ensemble
member, we have found it convenient to link the MPAS OCN analysis
files to a static input filename (e.g. <em class=
"file">mpas_analysis.nc</em>). The default DART filename is
<em class="file">dart_ics</em> - this may be moved or linked as
necessary.</p>
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
&amp;model_to_dart_nml
   model_to_dart_output_file = 'dart_ics'
   /
</pre></div>
<br />
<div class="namelist">
<pre>
&amp;model_nml
   model_analysis_filename  = 'mpas_analysis.nc'
   /

(partial namelist)
</pre></div>
<br />
<div class="namelist">
<pre>
&amp;mpas_vars_nml
   mpas_state_variables = '',
   mpas_state_bounds = '',
   /
</pre></div>
<br />
<br />
<p>The model_to_dart namelist includes:</p>
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
<td>model_to_dart_output_file</td>
<td>character(len=128)</td>
<td>The name of the DART file containing the model state derived
from the MPAS analysis file.</td>
</tr>
</tbody>
</table>
</div>
<br />
<p>Two more namelists need to be mentioned. The <a href=
"model_mod.html#Namelist">model_nml</a> namelist specifies the MPAS
analysis file to be used as the source. The <a href=
"model_mod.html#mpas_vars_nml">mpas_vars_nml</a> namelist specifies
the MPAS variables that will comprise the DART state vector.</p>
<p>For example:</p>
<pre>
&amp;mpas_vars_nml
   mpas_state_variables = 'temperature',  'QTY_TEMPERATURE',
                          'salinity',     'QTY_SALINITY',
                          'rho',          'QTY_DENSITY',
                          'u',            'QTY_EDGE_NORMAL_SPEED',
                          'h',            'QTY_SEA_SURFACE_HEIGHT'
                          'tracer1',      'QTY_TRACER_CONCENTRATION'
   /
</pre>
<br />
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>MODULES USED</h2>
<pre>
assim_model_mod.f90
types_mod.f90
location_mod.f90
model_to_dart.f90
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
<li>MPAS analysis file; <em class="file">mpas_analysis.nc</em></li>
<li>DART namelist file; <a href="work/input.nml">input.nml</a></li>
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
