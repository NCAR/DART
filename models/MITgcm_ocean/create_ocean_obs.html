<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program create_ocean_obs</title>
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
<h1>PROGRAM <em class="program">create_ocean_obs</em></h1>
<!-- version tag follows, do not edit -->
<p></p>
<p><em class="program">create_ocean_obs</em> is responsible for
converting an interim ASCII file of ocean observations into a DART
observation sequence file. The interim ASCII file is a simple
'whitespace separated' table where each row is an observation and
each column is specific information about the observation.</p>
<table border="2" cellpadding="3">
<tr>
<th align="left">column number</th>
<th align="left" width="30%">quantity</th>
<th align="left">description</th>
</tr>
<tr>
<td>1</td>
<td>longitude (in degrees)</td>
<td>longitude of the observation</td>
</tr>
<tr>
<td>2</td>
<td>latitude (in degrees)</td>
<td>latitude of the observation</td>
</tr>
<tr>
<td>3</td>
<td>depth (in meters)</td>
<td>depth of the observation</td>
</tr>
<tr>
<td>4</td>
<td>observation value</td>
<td>such as it is ...</td>
</tr>
<tr>
<td>5</td>
<td>vertical coordinate flag</td>
<td>see <a href=
"../../assimilation_code/location/threed_sphere/location_mod.html#location_type">
location_mod:location_type</a> for a full explanation. The short
explanation is that <em>surface == -1</em>, and <em>depth == 3</em>
There is a pathological difference between a surface observation
and an observation with a depth of zero.</td>
</tr>
<tr>
<td>6</td>
<td>observation variance</td>
<td>good luck here ...</td>
</tr>
<tr>
<td>7</td>
<td>Quality Control flag</td>
<td>integer value passed through to DART. There is a namelist
parameter for <em class="program">filter</em> to ignore any
observation with a QC value &lt;= <a href=
"../../assimilation_code/programs/filter/filter.html#Namelist">input_qc_threshold</a></td>
</tr>
<tr>
<td>8</td>
<td>obs_kind_name</td>
<td>a character string that must match a string in <a href=
"../../observations/forward_operators/obs_def_MITgcm_ocean_model_mod.html">
obs_def/obs_def_MITgcm_ocean_mod.f90</a></td>
</tr>
<tr>
<td>9</td>
<td>startDate_1</td>
<td>the year-month-date of the observation (YYYYMMDD format)</td>
</tr>
<tr>
<td>10</td>
<td>startDate_2</td>
<td>the hour-minute-second of the observation (HHMMSS format)</td>
</tr>
</table>
<p>For example:</p>
<pre>
273.7500 21.3500 -2.5018 28.0441  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 21.4500 -2.5018 28.1524  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 21.5500 -2.5018 28.0808  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 21.6500 -2.5018 28.0143  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 21.7500 -2.5018 28.0242  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 21.8500 -2.5018 28.0160  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 21.9500 -2.5018 28.0077  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 22.0500 -2.5018 28.3399  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 22.1500 -2.5018 27.8852  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
273.7500 22.2500 -2.5018 27.8145  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
...
</pre>
<p>It is always possible to combine observation sequence files with
the program <a href=
"../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">
obs_sequence_tool</a>, so it was simply convenient to generate a
separate file for each observation platform and type ('GLIDER' and
'TEMPERATURE'), however it is by no means required.</p>
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<hr>
<h2>MODULES USED</h2>
<p>Some of these modules use modules ... <strong>those</strong>
modules and namelists are not discussed here. probably should be
...</p>
<pre>
types_mod
utilities_mod
dart_MITocean_mod
obs_sequence_mod
</pre>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>This program has a namelist of its own, and some of the
underlying modules require namelists. To avoid duplication and,
possibly, some inconsistency in the documentation; only a list of
the required namelists is provided - with a hyperlink to the full
documentation for each namelist.</p>
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
"../../assimilation_code/modules/observations/obs_sequence_mod.html#Namelist">
obs_sequence_nml</a></td>
<td>write binary or ASCII observation sequence files</td>
</tr>
</table>
<p>We adhere to the F90 standard of starting a namelist with an
ampersand '&amp;' and terminating with a slash '/'. Consider
yourself forewarned that filenames that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the
namelist.</p>
<div class="namelist">
<pre>
<em class=
"call">namelist /create_ocean_obs_nml/ </em> year, month, day, &amp;
         tot_days, max_num, fname, output_name, lon1, lon2, lat1, lat2
</pre></div>
<div class="indent1"><!-- Description -->
<p>This namelist is read in a file called <em class=
"file">input.nml</em></p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">year</td>
<!--  type  -->
<td valign="top">integer <em class="unit">[default: 1996]</em></td>
<!--descript-->
<td>The first year of interest.</td>
</tr>
<tr><!--contents-->
<td valign="top">month</td>
<!--  type  -->
<td valign="top">integer <em class="unit">[default: 1]</em></td>
<!--descript-->
<td>The first month of interest.</td>
</tr>
<tr><!--contents-->
<td valign="top">day</td>
<!--  type  -->
<td valign="top">integer <em class="unit">[default: 1]</em></td>
<!--descript-->
<td>The first day of interest.</td>
</tr>
<tr><!--contents-->
<td valign="top">tot_days</td>
<!--  type  -->
<td valign="top">integer <em class="unit">[default: 31]</em></td>
<!--descript-->
<td>Stop processing after this many days.</td>
</tr>
<tr><!--contents-->
<td valign="top">max_num</td>
<!--  type  -->
<td valign="top">integer <em class="unit">[default:
800000]</em></td>
<!--descript-->
<td>The maximum number of observations to read/write.</td>
</tr>
<tr><!--contents-->
<td valign="top">fname</td>
<!--  type  -->
<td valign="top">character(len=129)<br>
<em class="unit">[default: 'raw_ocean_obs.txt']</em></td>
<!--descript-->
<td>The name of the interim ASCII file of observations.</td>
</tr>
<tr><!--contents-->
<td valign="top">output_name</td>
<!--  type  -->
<td valign="top">character(len=129)<br>
<em class="unit">[default: 'raw_ocean_obs_seq.out']</em></td>
<!--descript-->
<td>The output file name.</td>
</tr>
<tr><!--contents-->
<td valign="top">lon1</td>
<!--  type  -->
<td valign="top">real <em class="unit">[default: 0.0]</em></td>
<!--descript-->
<td>The leftmost longitude of interest.</td>
</tr>
<tr><!--contents-->
<td valign="top">lon2</td>
<!--  type  -->
<td valign="top">real <em class="unit">[default: 360.0]</em></td>
<!--descript-->
<td>The rightmost longitude of interest.</td>
</tr>
<tr><!--contents-->
<td valign="top">lat1</td>
<!--  type  -->
<td valign="top">real <em class="unit">[default: -90.0]</em></td>
<!--descript-->
<td>The most southern latitude of interest.</td>
</tr>
<tr><!--contents-->
<td valign="top">lat2</td>
<!--  type  -->
<td valign="top">real <em class="unit">[default: 90.0]</em></td>
<!--descript-->
<td>The most northern latitude of interest.</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<hr>
<h2>FILES</h2>
<ul>
<li>input namelist file: <em class="file">input.nml</em></li>
<li>input data file: as listed by <em class=
"file">input.nml</em><em class=
"code">&amp;create_ocean_obs_nml:fname</em></li>
<li>output data file: as listed by <em class=
"file">input.nml</em><em class=
"code">&amp;create_ocean_obs_nml:output_name</em></li>
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
<p>There are no error conditions specific to <em class=
"program">create_ocean_obs</em>.</p>
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
