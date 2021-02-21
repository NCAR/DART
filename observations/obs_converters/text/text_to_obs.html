<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program text_to_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">text_to_obs</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#DataSources">DATA SOURCES</a> / <a href=
"#Programs">PROGRAMS</a> / <a href="#Decisions">DECISIONS</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Text File to DART Converter</h2>
<h3>Overview</h3>
<p>If you have observations in spreadsheet or column format, in
text, with a single line per observation, then the files this
directory are a template for how to convert these observations into
a format suitable for DART use.</p>
<p>The workflow is usually:</p>
<ul>
<li>read in the needed information about each observation -
location, time, data value, observation type - from a data source
(usually a file)</li>
<li>call a series of DART library routines to construct a derived
type that contains all the information about a single
observation</li>
<li>call another set of DART library routines to put it into a
time-sorted series</li>
<li>repeat the last 2 steps until all observations are
processed</li>
<li>finally, call a write subroutine that writes out the entire
series to a file in a format that DART can read in</li>
</ul>
<p>It is not recommended that you try to mimic the ascii file
format by other means; the format is subject to change and the
library routines will continue to be supported even if the physical
format changes.</p>
<p>If your input data is in some kind of format like netCDF or HDF,
then one of the other converters (e.g. the MADIS ones for netCDF)
might be a better starting place for adapting code.</p>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<hr>
<h2>DATA SOURCES</h2>
<p>This part is up to you. For each observation you will need a
location, a data value, a type, a time, and some kind of error
estimate. The error estimate can be hardcoded in the converter if
they are not available in the input data. See below for more
details on selecting an appropriate error value.</p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<hr>
<h2>PROGRAMS</h2>
<p>The <em class="file">text_to_obs.f90</em> file is the source for
the main converter program. Look at the source code where it reads
the example data file. You will almost certainly need to change the
"read" statement to match your data format. The example code reads
each text line into a character buffer and then reads from that
buffer to parse up the data items.</p>
<p>To compile and test, go into the work subdirectory and run the
<em class="file">quickbuild.csh</em> script to build the converter
and a couple of general purpose utilities. <em class=
"file">advance_time</em> helps with calendar and time computations,
and the <em class="file">obs_sequence_tool</em> manipulates DART
observation files once they have been created.</p>
<p>To change the observation types, look in the <em class=
"file">DART/obs_def</em> directory. If you can find an
obs_def_XXX_mod.f90 file with an appropriate set of observation
types, change the 'use' lines in the converter source to include
those types. Then add that filename in the <em class=
"file">input.nml</em> namelist file to the &amp;preprocess_nml
namelist, the 'input_files' variable. Multiple files can be listed.
Then run quickbuild.csh again. It remakes the table of supported
observation types before trying to recompile the source code.</p>
<p>An example script for converting batches of files is in the
<em class="file">shell_scripts</em> directory. A tiny example data
file is in the <em class="file">data</em> directory. These are
<em>NOT</em> intended to be turnkey scripts; they will certainly
need to be customized for your use. There are comments at the top
of the script saying what options they include, and should be
commented enough to indicate where changes will be likely to need
to be made.</p>
<!--==================================================================-->
<a name="Decisions" id="Decisions"></a>
<hr>
<h2>DECISIONS YOU MIGHT NEED TO MAKE</h2>
<p>See the discussion in the <a href=
"../README.md#Decisions">obs_converters/README.md</a> page about
what options are available for the things you need to specify.
These include setting a time, specifying an expected error, setting
a location, and an observation type.</p>
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
