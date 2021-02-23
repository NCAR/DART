<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>program gitm_blocks_to_netcdf</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">gitm_blocks_to_netcdf</em></h1>
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
<p>The <a href=
"http://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=GITM">Global
Ionosphere Thermosphere Model (GITM)</a> is a 3-dimensional
spherical code that models the Earth's thermosphere and ionosphere
system using a stretched grid in latitude and altitude. For a
fuller description of using GITM within DART, please see the
<a href="model_mod.html">DART GITM model documentation</a>.<br />
<br />
<em class="program">gitm_blocks_to_netcdf</em> is the program that
reads GITM restart files (i.e. <em class="file">b?????.rst</em>)
and creates a DART output/restart file (e.g. <em class=
"file">perfect_ics, filter_ics, ... </em>).<br />
<br />
The list of variables used to create the DART state vector are
specified in the <em class="file">input.nml</em> file.<br />
<br />
Conditions required for successful execution of <em class=
"program">gitm_blocks_to_netcdf</em>:</p>
<ul>
<li>a valid <em class="file">input.nml</em> namelist file for
DART</li>
<li>a valid <em class="file">UAM.in</em> control file for GITM</li>
<li>a set of <em class="file">b?????.rst</em> data files for
GITM</li>
<li>a <em class="file">header.rst</em> file for GITM</li>
<li>the DART/GITM interfaces must be compiled in a manner
consistent with the GITM data and control files. The following GITM
source files are required to build <em>any</em> DART interface:
<ul>
<li>models/gitm/GITM2/src/ModConstants.f90</li>
<li>models/gitm/GITM2/src/ModEarth.f90</li>
<li>models/gitm/GITM2/src/ModKind.f90</li>
<li>models/gitm/GITM2/src/ModSize.f90</li>
<li>models/gitm/GITM2/src/ModTime.f90</li>
<li>models/gitm/GITM2/src/time_routines.f90</li>
</ul>
Versions of these are included in the DART release. <em class=
"file">ModSize.f90</em>, in particular, must match what was used to
create the <em class="file">b????.rst</em> files.</li>
</ul>
<p>The individual model instances are run in unique directories.
This is also where the converter routines <em class=
"program">gitm_blocks_to_netcdf</em> and <em class=
"program">dart_to_gitm</em> are run. This makes it easy to use a
single 'static' name for the input and output filenames. <em class=
"program">advance_model.csh</em> is responsibile for linking the
appropriate files to these static filenames.</p>
<p>The simplest way to test the converter is to compile GITM and
run a single model state forward using <em class=
"program">work/clean.sh</em>. To build GITM ... download GITM and
unpack the code into <em class="file">DART/models/gitm/GITM2</em>
and follow these instructions:</p>
<div class="unix">
<pre>
cd models/gitm/GITM2
./Config.pl -install -compiler=ifortmpif90 -earth
make
cd ../work
./clean.sh 1 1 0 150.0 170.0 1.0 
</pre></div>
<p><!-- makes the 'top' link line up correctly --></p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>NAMELIST</h2>
<p>We adhere to the F90 standard of starting a namelist with an
ampersand '&amp;' and terminating with a slash '/' for all our
namelist input. Character strings that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the
namelist.</p>
<div class="namelist">
<pre>
&amp;gitm_blocks_to_netcdf_nml
   gitm_blocks_to_netcdf_output_file = 'dart_ics',
   /

&amp;model_nml
   gitm_restart_dirname         = 'advance_temp_e1/UA/restartOUT',
   assimilation_period_days     = 0,
   assimilation_period_seconds  = 1800,
   model_perturbation_amplitude = 0.2,
   output_state_vector          = .false.,
   calendar                     = 'Gregorian',
   debug                        = 0,
   gitm_state_variables = 'Temperature',            'QTY_TEMPERATURE',
                          'eTemperature',           'QTY_TEMPERATURE_ELECTRON',
                          'ITemperature',           'QTY_TEMPERATURE_ION',
                          'iO_3P_NDensityS',        'QTY_DENSITY_NEUTRAL_O3P',
   ...
</pre></div>
<table border="0" cellpadding="3" width="100%" summary=
'gitm_blocks_to_netcdf namelist description'>
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">
gitm_blocks_to_netcdf_output_file   </td>
<!--  type  -->
<td valign="top">character(len=128)  </td>
<!--descript-->
<td>The name of the DART file containing the model state derived
from the GITM restart files.</td>
</tr>
</table>
<br />
<p>The full description of the <em class="code">model_nml</em>
namelist is documented in the <a href=
"model_mod.html#Namelist">gitm model_mod</a>, but the most
important variable for <em class=
"program">gitm_blocks_to_netcdf</em> is repeated here.</p>
<table border="0" cellpadding="3" width="100%" summary=
'partial model_nml namelist description'>
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">gitm_restart_dirname  </td>
<!--  type  -->
<td valign="top">character(len=256)  </td>
<!--descript-->
<td>The name of the directory containing the GITM restart files and
runtime control information.</td>
</tr>
<tr><!--contents-->
<td valign="top">gitm_state_variables </td>
<!--  type  -->
<td valign="top">character(len=32),  <br />
dimension(2,80)</td>
<!--descript-->
<td>The list of variable names in the gitm restart file to use to
create the DART state vector and their corresponding DART kind. The
default list is specified in <a href=
"model_mod.nml">model_mod.nml</a></td>
</tr>
</table>
<p><!-- makes the 'top' link line up correctly --></p>
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>MODULES USED</h2>
<pre>
obs_def_upper_atm_mod.f90
assim_model_mod.f90
types_mod.f90
location/threed_sphere/location_mod.f90
models/gitm/GITM2/src/ModConstants.f90
models/gitm/GITM2/src/ModEarth.f90
models/gitm/GITM2/src/ModKind.f90
models/gitm/GITM2/src/ModSize.f90
models/gitm/GITM2/src/ModTime.f90
models/gitm/GITM2/src/time_routines.f90
models/gitm/dart_gitm_mod.f90
models/gitm/gitm_blocks_to_netcdf.f90
models/gitm/model_mod.f90
null_mpi_utilities_mod.f90
obs_kind_mod.f90
random_seq_mod.f90
time_manager_mod.f90
utilities_mod.f90
</pre>
<p><!-- makes the 'top' link line up correctly --></p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FILES Read</h2>
<ul>
<li>gitm restart files: <em class="file">b????.rst</em></li>
<li>gitm control files: <em class="file">header.rst</em></li>
<li>gitm control files: <em class="file">UAM.in.rst</em></li>
<li>DART namelist file: <em class="file">input.nml</em></li>
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
<ul>
<li>The official <em class="program">GITM</em> site is: can be
found at <a href=
"http://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=GITM">ccmc.gsfc.nasa.gov/models/modelinfo.php?model=GITM</a></li>
</ul>
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
