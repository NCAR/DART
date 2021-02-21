<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program tc_to_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">tc_to_obs</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#DataSources">DATA SOURCES</a> / <a href=
"#Programs">PROGRAMS</a> / <a href="#ObservationalErrors">EXPECTED
ERROR</a> / <a href="#Namelist">NAMELIST</a> / <a href=
"#KnownBugs">KNOWN BUGS</a> / <a href="#FuturePlans">PLANS</a> /
<a href="#Legalese">TERMS OF USE</a>
<h1>Tropical Cyclone ATCF File to DART Converter</h1>
<h2>Overview</h2>
<p>Tropical Cyclone data created by the 'Automated Tropical Cyclone
Forecast (ATCF) System' can be converted into DART observations of
the storm center location, minimum sea level pressure, and maximum
wind speed. Several of the options can be customized at runtime by
setting values in a Fortran namelist. See the <a href=
"#Namelist">namelist</a> section below for more details. In the
current release of DART only the <a href=
"../../../models/wrf/model_mod.html">WRF model</a> has forward
operator code to generate expected obs values for these vortex
observations.</p>
<p><a href=
"http://www.ral.ucar.edu/hurricanes/realtime/index.php#about_atcf_data_files">
This webpage</a> documents many things about the ATCF system and
the various file formats that are used for storm track data and
other characteristics.</p>
<p>The converter in this directory is only configured to read the
packed "b-deck" format (as described on the webpage referenced
above). There are sections in the fortran code which can be filled
in to read other format variants. This should mostly be a matter of
changing the read format string to match the data in the file.</p>
<p><br></p>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>DATA SOURCES</h2>
<p>A collection of past storm ATCF information can be found
<a href="http://www.ral.ucar.edu/hurricanes/repository">here</a>.
For each observation you will need a location, a data value, a
type, a time, and some kind of error estimate. The error estimates
will need to be hardcoded or computed in the converter since they
are not available in the input data. See below for more details on
selecting an appropriate error value.</p>
<p><br></p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PROGRAMS</h2>
<p>The <em class="file">tc_to_obs.f90</em> file is the source for
the main converter program. Look at the source code where it reads
the example data file. Given the variety of formatting details in
different files, you may quite possibly need to change the "read"
statement to match your data format. There is a 'select case'
section which is intended to let you add more formats and select
them at runtime via namelist.</p>
<p>To compile and test, go into the work subdirectory and run the
<em class="file">quickbuild.csh</em> script to build the converter
and a couple of general purpose utilities. <em class=
"file">advance_time</em> helps with calendar and time computations,
and the <em class="file">obs_sequence_tool</em> manipulates DART
observation files once they have been created.</p>
<p>This converter creates observation types defined in the
<em class=
"file">DART/observations/forward_operators/obs_def_vortex_mod.f90</em>
file. This file must be listed in the <em class=
"file">input.nml</em> namelist file, in the <em class=
"file">&amp;preprocess_nml</em> namelist, in the 'input_files'
variable, for any programs which are going to process these
observations. If you have to change the <em class=
"file">&amp;preprocess_nml</em> namelist you will have to run
<em class="file">quickbuild.csh</em> again to build and execute the
<em class="program">preprocess</em> program before compiling other
executables. It remakes the table of supported observation types
before trying to recompile other source code.</p>
<p>There is an example b-deck data file in the <em class=
"file">data</em> directory. This format is what is supported in the
code as distributed. There are other variants of this format which
have more spaces so the columns line up, and variants which have
many more fields than what is read here.</p>
<p><br></p>
<!--==================================================================-->
<a name="ObservationalErrors" id="ObservationalErrors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>SPECIFYING EXPECTED ERROR</h2>
<p>The ATCF files DO NOT include any estimated error values. The
source code currently has hardcoded values for location, sea level
pressure, and max wind errors. These may need to be adjusted as
needed if they do not give the expected results.</p>
<p><br></p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;tc_to_obs_nml
   input_atcf_file         = 'input.txt'
   fileformat              = 'b-deck'
   obs_out_file            = 'obs_seq.out'
   append_to_existing_file = .false.
   debug                   = .false.
 /
</pre></div>
<br>
<br>
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
<td>input_atcf_file</td>
<td>character(len=256)</td>
<td>Name of the input ascii text file in ATCF format.</td>
</tr>
<tr>
<td>fileformat</td>
<td>character(len=128)</td>
<td>Currently only supports 'b-deck' but if other format strings
are added, can switch at runtime between reading different
varieties of ATCF file formats.</td>
</tr>
<tr>
<td>obs_out_file</td>
<td>character(len=256)</td>
<td>Name of the output observation sequence file to create.</td>
</tr>
<tr>
<td>append_to_existing_file</td>
<td>logical</td>
<td>If .false., this program will overwrite an existing file. If
.true. and if a file already exists with the same name the newly
converted observations will be appended to that file. Useful if you
have multiple small input files that you want to concatenate into a
single output file. However, there is no code to check for
duplicated observations. If this is .true. and you run the
converter twice you will get duplicate observations in the file
which is bad. (It will affect the quality of your assimilation
results.) Use with care.<br>
You can concatenate multiple obs sequence files as a postprocessing
step with the <a href=
"../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">
observation sequence tool</a> which comes with DART and is in fact
built by the quickbuild.csh script in the TC converter work
directory.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>Set to .true. to print out more details during the conversion
process.</td>
</tr>
</tbody>
</table>
</div>
<p><br></p>
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
<a name="KnownBugs" id="KnownBugs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>KNOWN BUGS</h2>
<p>none</p>
<p><br></p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>If users add support for some of the other format variants, the
DART group would be happy to accept code contributions.</p>
<p><br></p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
