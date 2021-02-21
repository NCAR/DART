<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program dwl_to_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">dwl_to_obs</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#DataSources">DATA SOURCES</a> / <a href=
"#Programs">PROGRAMS</a> / <a href="#References">REFERENCES</a> /
<a href="#Errors">ERRORS</a> / <a href="#FuturePlans">PLANS</a> /
<a href="#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<h3>DWL to DART Converter</h3>
<p>These are Doppler Wind Lidar measurements which have previously
been extracted from the incoming format and output in ascii format,
one pair of wind component observations per line. This converter
reads in the ascii file and outputs the data in DART observation
sequence (obs_seq) format.</p>
<p>This is OSSE data from a satellite which is expected to be
launched in 2015. Information on the satellite mission is here at
<a href=
"http://en.wikipedia.org/wiki/ADM-Aeolus">http://en.wikipedia.org/wiki/ADM-Aeolus</a>.</p>
<p>The workflow is:</p>
<ul>
<li>read in the needed information about each observation -
location, time, observation values, obs errors - from an ascii
file</li>
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
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<hr>
<h2>DATA SOURCES</h2>
<p>Matic Savli at University of Ljubljana has programs which read
the expected instrument formats, do the proper conversions, and
write out ascii lines, one per wind observation.</p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<hr>
<h2>PROGRAMS</h2>
<p>The <em class="file">dwl_to_obs.f90</em> file is the source for
the main converter program. There is a sample data file in the
"data" directory. The converter reads each text line into a
character buffer and then reads from that buffer to parse up the
data items.</p>
<p>To compile and test, go into the work subdirectory and run the
<em class="file">quickbuild.csh</em> script to build the converter
and a couple of general purpose utilities. <em class=
"file">advance_time</em> helps with calendar and time computations,
and the <em class="file">obs_sequence_tool</em> manipulates DART
observation files once they have been created.</p>
<p>The observation types are defined in <em class=
"file">DART/obs_def/obs_def_dwl_mod.f90</em>. That filename must be
added to the <em class="file">input.nml</em> namelist file, to the
&amp;preprocess_nml namelist, the 'input_files' variable before
compiling any program that uses these observation types. Multiple
files can be listed. Then run quickbuild.csh again. It remakes the
table of supported observation types before trying to recompile the
source code.</p>
<p>An example script for converting batches of files is in the
<em class="file">shell_scripts</em> directory. It will need
customization before being used.</p>
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
