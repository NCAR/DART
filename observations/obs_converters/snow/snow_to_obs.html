<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>Snow observations</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">snow_to_obs</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#DataSources">DATA SOURCES</a> / <a href=
"#Programs">PROGRAMS</a> / <a href="#Namelist">NAMELIST</a> /
<a href="#Errors">ERRORS</a> / <a href="#FuturePlans">PLANS</a> /
<a href="#Legalese">TERMS OF USE</a>
<h2>MODIS Snowcover Fraction Observation Converter</h2>
<h4>Overview</h4>
<p>There are several satellite sources for snow observations.
Generally the data is distributed in HDF-EOS format. The converter
code in this directory DOES NOT READ HDF FILES as input. It expects
the files to have been preprocessed to contain text, one line per
observation, with northern hemisphere data only.</p>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<hr>
<h2>DATA SOURCES</h2>
<p>not sure.</p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<hr>
<h2>PROGRAMS</h2>
<p>The <em class="file">snow_to_obs.f90</em> file is the source for
the main converter program.</p>
<p>To compile and test, go into the work subdirectory and run the
<em class="file">quickbuild.csh</em> script to build the converter
and a couple of general purpose utilities. <em class=
"file">advance_time</em> helps with calendar and time computations,
and the <em class="file">obs_sequence_tool</em> manipulates DART
observation files once they have been created.</p>
<p>This converter creates observations of the
"MODIS_SNOWCOVER_FRAC" type.</p>
<p>There is another program in this directory called <em class=
"file">snow_to_obs_netcdf.f90</em> which is a prototype for reading
netcdf files that contain some metadata and presumably have been
converted from the original HDF. THIS HAS NOT BEEN TESTED but if
you have such data, please contact <a href=
'mailto:dart@ucar.edu'>dart@ucar.edu</a> for more assistance. If
you write something that reads the HDF-EOS MODIS files directly,
please, please contact us! Thanks.</p>
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;snow_to_obs_nml
  longrid         = 360,
  latgrid         = 90, 
  year            = 2000, 
  doy             = 1,
  snow_input_file = 'snowdata.input', 
  missing_value   = -20.0, 
  debug           = .false.
/
</pre></div>
<table border="0" cellpadding="10" width="100%" summary=
'snow_to_obs namelist description'>
<thead align="left">
<tr>
<th>Item</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>longrid</td>
<td>integer</td>
<td>The number of divisions in the longitude dimension.</td>
</tr>
<tr>
<td>latgrid</td>
<td>integer</td>
<td>The number of divisions in the latitude dimension. This
converter assumes the data is for the northern hemisphere only. A
namelist item could be added to select northern verses southern
hemisphere if needed.</td>
</tr>
<tr>
<td>year</td>
<td>integer</td>
<td>The year number of the data.</td>
</tr>
<tr>
<td>doy</td>
<td>integer</td>
<td>The day number in the year. Valid range 1 to 365 in a non-leap
year, 1 to 366 in a leap year.</td>
</tr>
<tr>
<td>snow_input_file</td>
<td>character(len=128)</td>
<td>The name of the input file.</td>
</tr>
<tr>
<td>missing_value</td>
<td>real(r8)</td>
<td>The value used to mark missing data.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>If set to .true. the converter will print out more information
as it does the conversion.</td>
</tr>
</tbody>
</table>
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
<a name="KnownBugs" id="KnownBugs"></a>
<hr>
<h2>KNOWN BUGS</h2>
<p>This program is hardcoded to read only northern hemisphere data.
It should handle global values.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<p>This program should use the HDF-EOS libraries to read the native
MODIS granule files. Right now the ascii intermediate files contain
no metadata, so if the namelist values don't match the actual
division of the globe, bad things will happen.</p>
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
