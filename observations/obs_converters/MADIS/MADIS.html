<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>MADIS Data</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MADIS Data Ingest System</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#DataSources">DATA SOURCES</a> / <a href=
"#Programs">PROGRAMS</a> / <a href="#Modules">MODULES</a> /
<a href="#Namelist">NAMELIST</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">FUTURE PLANS</a> / <a href="#Legalese">TERMS
OF USE</a>
<h2>Overview</h2>
<p>The <a href="http://madis.noaa.gov/">MADIS</a> (Meteorological
Assimilation Data Ingest System) service provides access to
real-time and archived data of a variety of types, with added
Quality Control (QC) and integration of data from a variety of
sources.</p>
<p>To convert a series of MADIS data files (where different types
of observations are distributed in separate files), one high level
view of the workflow is:</p>
<ol>
<li>convert each madis file, by platform type, into an obs_seq
file. one file in, one file out. no time changes. use the
<em class="file">shell_scripts/madis_conv.csh</em> script. there
are script options for hourly output files, or a single daily
output file.</li>
<li>if you aren't using the wrf preprocessing program, you're ready
to go.</li>
<li>if you do want to do subsequent wrf preprocessing, you need to:
<ol>
<li>decide on the windowing. each platform has a different
convention and if you're going to put them into the wrf
preprocessing you'll need to have the windowing match. use the
<em class="file">shell_scripts/windowing.csh</em> script.</li>
<li>the wrf preprocessing takes a list of files and assumes they
will all be assimilated at the same time, for superob'ing purposes,
so it should match the expected assimilation window when running
filter.</li>
</ol>
</li>
</ol>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<hr>
<h2>DATA SOURCES</h2>
<p><a href="http://madis.noaa.gov/">http://madis.noaa.gov</a></p>
<p>There are two satellite wind converter programs; the one in this
directory and one in the <a href="../SSEC/SSEC.html">SSEC</a>
directory. The observations distributed here come from <a href=
"http://www.nesdis.noaa.gov">NESDIS</a>. The SSEC observations are
processed by SSEC itself and will differ from the observations
converted here.</p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<hr>
<h2>PROGRAMS</h2>
<p>The programs in the <em class=
"file">DART/observations/MADIS/</em> directory extract data from
the distribution files and create DART observation sequence
(obs_seq) files. Build them in the <em class="file">work</em>
directory by running the <em class="program">./quickbuild.csh</em>
script. In addition to the converters, the <em class=
"file">advance_time</em> and <em class=
"file">obs_sequence_tool</em> utilities will be built.</p>
<p>There are currently converters for these data types:</p>
<table border="0" cellpadding="3" summary="">
<tbody valign="top">
<tr>
<td>ACARS aircraft T,U,V,Q data</td>
<td>convert_madis_acars</td>
</tr>
<tr>
<td>Marine surface data</td>
<td>convert_madis_marine</td>
</tr>
<tr>
<td>Mesonet surface data</td>
<td>convert_madis_mesonet</td>
</tr>
<tr>
<td>Metar data</td>
<td>convert_madis_metar</td>
</tr>
<tr>
<td>Wind Profiler data</td>
<td>convert_madis_profiler</td>
</tr>
<tr>
<td>Rawinsonde/Radiosonde data</td>
<td>convert_madis_rawin</td>
</tr>
<tr>
<td>Satellite Wind data</td>
<td>convert_madis_satwnd</td>
</tr>
</tbody>
</table>
<p>Example data files are in the <em class="file">data</em>
directory. Example scripts for converting batches of these files
are in the <em class="file">shell_scripts</em> directory. These are
<em>NOT</em> intended to be turnkey scripts; they will certainly
need to be customized for your use. There are comments at the top
of the scripts saying what options they include, and should be
commented enough to indicate where changes will be likely to need
to be made.</p>
<p>Several converters have compile-time choices for outputting
various types of moist variables. Check the source code for more
details. Some converters also read multiple T/F strings from the
console (standard input) to control at run-time what types of
observations to convert. Again, check the source code for more
details.</p>
<p>Each converter has hard-coded input and output filenames:</p>
<table border="0" cellpadding="5" summary="input/output filenames">
<tbody valign="top">
<tr>
<td>convert_madis_acars:</td>
<td>acars_input.nc</td>
<td>obs_seq.acars</td>
</tr>
<tr>
<td>convert_madis_marine:</td>
<td>marine_input.nc</td>
<td>obs_seq.marine</td>
</tr>
<tr>
<td>convert_madis_mesonet:</td>
<td>mesonet_input.nc</td>
<td>obs_seq.mesonet</td>
</tr>
<tr>
<td>convert_madis_metar:</td>
<td>metar_input.nc</td>
<td>obs_seq.metar</td>
</tr>
<tr>
<td>convert_madis_profiler:</td>
<td>profiler_input.nc</td>
<td>obs_seq.profiler</td>
</tr>
<tr>
<td>convert_madis_rawin:</td>
<td>rawin_input.nc</td>
<td>obs_seq.rawin</td>
</tr>
<tr>
<td>convert_madis_satwnd:</td>
<td>satwnd_input.nc</td>
<td>obs_seq.satwnd</td>
</tr>
</tbody>
</table>
<p>The expected usage pattern is that a script will copy, rename,
or make a symbolic link from the actual input file (which often
contains a timestamp in the name) to the fixed input name before
conversion, and move the output file to an appropriate filename
before the next invocation of the converter. If an existing
observation sequence file of the same output name is found when the
converter is run again, it will open that file and append the next
set of observations to it.</p>
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
<a name="KnownBugs" id="KnownBugs"></a>
<hr>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<p>none</p>
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
