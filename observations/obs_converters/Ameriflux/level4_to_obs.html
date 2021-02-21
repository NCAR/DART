<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program level4_to_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">level4_to_obs</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#DataSources">DATA
SOURCES</a> / <a href="#Programs">PROGRAMS</a> / <a href=
"#Decisions">DECISIONS</a> / <a href="#References">REFERENCES</a> /
<a href="#Errors">ERRORS</a> / <a href="#FuturePlans">PLANS</a> /
<a href="#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<h3>AmeriFlux Level 4 data to DART Observation Sequence
Converter</h3>
<p>This routine is designed to convert the flux tower Level 4 data
from the <a href="http://ameriflux.lbl.gov">AmeriFlux</a> network
of observations from micrometeorological tower sites. AmeriFlux is
part of <a href="http://fluxnet.ornl.gov">FLUXNET</a> and the
converter is hoped to be a suitable starting point for the
conversion of observations from FLUXNET. As of May 2012, I have not
yet tried to work with any other observations from FLUXNET.<br>
<br>
The AmeriFlux Level 4 products are recorded using the local time.
DART observation sequence files use GMT. For more information about
AmeriFlux data products, go to <a href=
"http://ameriflux.lbl.gov">http://ameriflux.lbl.gov</a>.</p>
<div class="warning">
<div class="title">
<p>Warning</p>
</div>
<p>There was a pretty severe bug in the converter that swapped
latent heat flux and sensible heat flux. The bug was present
through revision 7200. It has been corrected in all subsequent
versions.</p>
</div>
<p>The workflow is usually:</p>
<ol>
<li>download the Level 4 data for the towers and years in question
(<a href="#DataSources">see DATA SOURCES below</a>)</li>
<li>record the TIME ZONE, latitude, longitude, and elevation for
each tower</li>
<li>build the DART executables with support for the tower
observations. This is done by running <em class=
"program">preprocess</em> with <em class=
"file">obs_def_tower_mod.f90</em> in the list of <em class=
"code">input_files</em> for <em class=
"code">preprocess_nml</em>.</li>
<li>provide basic tower information via the <em class=
"program">level4_to_obs_nml</em> namelist since this information is
not contained in the Level 4 data file</li>
<li>convert each Level 4 data file individually using <em class=
"program">level4_to_obs</em></li>
<li>combine all output files for the region and timeframe of
interest into one file using <a href=
"../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">
obs_sequence_tool</a></li>
</ol>
<p>For some models (CLM, for example), it is required to reorganize
the observation sequence files into a series of files that contains
ONLY the observations for each assimilation. This can be achieved
with the <a href="makedaily.sh">makedaily.sh</a> script.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
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
&amp;level4_to_obs_nml
   text_input_file = 'textdata.input',
   obs_out_file    = 'obs_seq.out',
   year            = -1,
   timezoneoffset  = -1,
   latitude        = -1.0,
   longitude       = -1.0,
   elevation       = -1.0,
   flux_height     = -1.0,
   maxgoodqc       = 3,
   verbose         = .false.
   /
</pre></div>
<div>
<table border="0" cellspacing="10" width="100%" summary=
'level4_to_obs namelist description'>
<thead align="left">
<tr>
<th>Contents</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr><!--contents-->
<td>text_input_file</td>
<!--  type  -->
<td>character(len=128)</td>
<!--descript-->
<td>Name of the Level 4 ASCII file of comma-separated values. This
may be a relative or absolute filename.</td>
</tr>
<tr><!--contents-->
<td>obs_out_file</td>
<!--  type  -->
<td>character(len=128)</td>
<!--descript-->
<td>Name of the output observation sequence file.</td>
</tr>
<tr><!--contents-->
<td>year</td>
<!--  type  -->
<td>integer</td>
<!--descript-->
<td>The year of the observations in the Level 4 text file.</td>
</tr>
<tr><!--contents-->
<td>timezoneoffset</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>the time zone offset (in hours) of the station. The tower
observation times are local time, we need to convert them to
GMT.</td>
</tr>
<tr><!--contents-->
<td>latitude</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>Latitude (in degrees N) of the tower.</td>
</tr>
<tr><!--contents-->
<td>longitude</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>Longitude (in degrees E) of the tower. For internal
consistency, DART uses longitudes in the range [0,360]. An input
value of -90 will be converted to 270, for example.</td>
</tr>
<tr><!--contents-->
<td>elevation</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>surface elevation (in meters) of the tower.</td>
</tr>
<tr><!--contents-->
<td>flux_height</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>height (in meters) of the flux instrument on the tower.</td>
</tr>
<tr><!--contents-->
<td>maxgoodqc</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>maximum value of any observation quality control flag to pass
through to the output observation sequence. Keep in mind that
<em class="program">filter</em> has the ability to discriminate on
the value, so there is really little to be gained by rejecting them
during the conversion.</td>
</tr>
<tr><!--contents-->
<td>verbose</td>
<!--  type  -->
<td>logical</td>
<!--descript-->
<td>Print extra information during the <em class=
"program">level4_to_obs</em> execution.</td>
</tr>
</tbody>
</table>
</div>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<hr>
<h2>DATA SOURCES</h2>
<p>The data was acquired from <a href=
"http://cdiac.ornl.gov/ftp/ameriflux/data/Level4/Sites_ByName">http://cdiac.ornl.gov/ftp/ameriflux/data/Level4/Sites_ByName</a><br>

and have names like <em class="file">USBar2004_L4_h.txt,
USHa12004_L4_h.txt, USNR12004_L4_h.txt, USSP32004_L4_h.txt,
USSRM2004_L4_h.txt, USWCr2004_L4_h.txt, USWrc2004_L4_h.txt,
...</em><br>
<br>
The Level 4 products in question are ASCII files of comma-separated
values taken every 30 minutes for an entire year. The first line is
a comma-separated list of column descriptors, all subsequent lines
are comma-separated numerical values. The converter presently
searches for the columns pertaining to <em>NEE_or_fMDS</em>,
<em>H_f</em>, <em>LE_f</em>, their corresponding quality control
fields, and those columns pertaining to the time of the
observation. These values are mapped as follows:</p>
<table width="100%" cellpadding="10" summary=
'data products summary'>
<tr>
<th align="left">Level 4 units</th>
<th align="left">Level 4 variable</th>
<th align="left">description</th>
<th align="left">DART type</th>
<th align="left">DART kind</th>
<th align="left">DART units</th>
</tr>
<tr>
<td colspan="6">
<hr></td>
</tr>
<tr><!-- >Level 4 units   <-->
<td>W/m^2</td>
<!-- >Level 4 variable<-->
<td>LE_f</td>
<!-- >description     <-->
<td>Latent Heat Flux</td>
<!-- >DART type       <-->
<td>TOWER_LATENT_HEAT_FLUX</td>
<!-- >DART kind       <-->
<td>QTY_LATENT_HEAT_FLUX</td>
<!-- >DART units      <-->
<td>W/m^2</td>
</tr>
<tr><!-- >Level 4 units   <-->
<td>[0-3]</td>
<!-- >Level 4 variable<-->
<td>LE_fqc</td>
<!-- >description     <-->
<td>QC for LE_f</td>
<!-- >DART type       <-->
<td>N/A</td>
<!-- >DART kind       <-->
<td>N/A</td>
<!-- >DART units      <-->
<td>same</td>
</tr>
<tr>
<td colspan="6">
<hr></td>
</tr>
<tr><!-- >Level 4 units   <-->
<td>W/m^2</td>
<!-- >Level 4 variable<-->
<td>H_f</td>
<!-- >description     <-->
<td>Sensible Heat Flux</td>
<!-- >DART type       <-->
<td>TOWER_SENSIBLE_HEAT_FLUX</td>
<!-- >DART kind       <-->
<td>QTY_SENSIBLE_HEAT_FLUX</td>
<!-- >DART units      <-->
<td>W/m^2</td>
</tr>
<tr><!-- >Level 4 units   <-->
<td>[0-3]</td>
<!-- >Level 4 variable<-->
<td>H_fqc</td>
<!-- >description     <-->
<td>QC for H_f</td>
<!-- >DART type       <-->
<td>N/A</td>
<!-- >DART kind       <-->
<td>N/A</td>
<!-- >DART units      <-->
<td>same</td>
</tr>
<tr>
<td colspan="6">
<hr></td>
</tr>
<tr><!-- >Level 4 units   <-->
<td>umolCO2/m^2/s</td>
<!-- >Level 4 variable<-->
<td>NEE_or_fMDS</td>
<!-- >description     <-->
<td>Net Ecosystem Production</td>
<!-- >DART type       <-->
<td>TOWER_NETC_ECO_EXCHANGE</td>
<!-- >DART kind       <-->
<td>QTY_NET_CARBON_PRODUCTION</td>
<!-- >DART units      <-->
<td>gC/m^2/s</td>
</tr>
<tr><!-- >Level 4 units   <-->
<td>[0-3]</td>
<!-- >Level 4 variable<-->
<td>NEE_or_fMDSqc</td>
<!-- >description     <-->
<td>QC for NEE_or_fMDS</td>
<!-- >DART type       <-->
<td>N/A</td>
<!-- >DART kind       <-->
<td>N/A</td>
<!-- >DART units      <-->
<td>same</td>
</tr>
</table>
<p>The <em class="code">LE_fqc</em>, <em class="code">H_fqc</em>,
and <em class="code">NEE_or_fMDSqc</em> variables use the following
convention:</p>
<blockquote>0 = original, 1 = category A (most reliable), 2 =
category B (medium), 3 = category C (least reliable). (Refer to
Reichstein et al. 2005 Global Change Biology for more
information)</blockquote>
<br>
<p>I am repeating the AmeriFlux <a href=
"http://ameriflux.lbl.gov/Data/Pages/DataUsagePolicy.aspx">Data
Fair-Use Policy</a> because I believe it is important to be a good
scientific citizen:</p>
<blockquote>"The AmeriFlux data provided on this site are freely
available and were furnished by individual AmeriFlux scientists who
encourage their use.<br>
<br>
Please kindly inform in writing (or e-mail) the appropriate
AmeriFlux scientist(s) of how you intend to use the data and of any
publication plans. It is also important to contact the AmeriFlux
investigator to assure you are downloading the latest revision of
the data and to prevent potential misuse or misinterpretation of
the data.<br>
<br>
Please acknowledge the data source as a citation or in the
acknowledgments if no citation is available. If the AmeriFlux
Principal Investigators (PIs) feel that they should be acknowledged
or offered participation as authors, they will let you know and we
assume that an agreement on such matters will be reached before
publishing and/or use of the data for publication.<br>
<br>
If your work directly competes with the PI's analysis they may ask
that they have the opportunity to submit a manuscript before you
submit one that uses unpublished data. In addition, when publishing
please acknowledge the agency that supported the research.<br>
<br>
Lastly, we kindly request that those publishing papers using
AmeriFlux data provide reprints to the PIs providing the data and
to the AmeriFlux archive via ameriflux.lbl.gov."</blockquote>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<hr>
<h2>PROGRAMS</h2>
<p>The <em class="file">level4_to_obs.f90</em> file is the source
for the main converter program. Look at the source code where it
reads the example data file. You will almost certainly need to
change the "read" statement to match your data format. The example
code reads each text line into a character buffer and then reads
from that buffer to parse up the data items.</p>
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
