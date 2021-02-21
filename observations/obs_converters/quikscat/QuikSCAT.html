<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>QuikSCAT SeaWinds Data</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>QuikSCAT SeaWinds Data</h1>
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
<p>NASA's QuikSCAT mission is described in <a href=
"http://winds.jpl.nasa.gov/missions/quikscat/index.cfm">http://winds.jpl.nasa.gov/missions/quikscat/</a>.
"QuikSCAT" refers to the satellite, "SeaWinds" refers to the
instrument that provides near-surface wind speeds and directions
over large bodies of water. QuikSCAT has an orbit of about 100
minutes, and the SeaWinds microwave radar covers a swath under the
satellite. The swath is comprised of successive scans (or rows) and
each scan has many wind-vector-cells (WVCs). For the purpose of
this document, we will focus only the <b>Level 2B</b> product at
25km resolution. If you go to the official JPL data distribution
site <a href=
"http://podaac.jpl.nasa.gov/DATA_CATALOG/quikscatinfo.html">http://podaac.jpl.nasa.gov/DATA_CATALOG/quikscatinfo.html</a>
, we are using the product labelled <b>L2B OWV 25km Swath</b>. Each
orbit consists of (potentially) 76 WVCs in each of 1624 rows or
scans. The azimuthal diversity of the radar returns affects the
error characteristics of the retrieved wind speeds and directions,
as does rain, interference of land in the radar footprint, and very
low wind speeds. Hence, not all wind retrievals are created
equal.</p>
<p>The algorithm that converts the 'sigma naughts' (the measure of
radar backscatter) into wind speeds and directions has multiple
solutions. Each candidate solution is called an 'ambiguity', and
there are several ways of choosing 'the best' ambiguity. Beauty is
in the eye of the beholder. At present, the routine to convert the
original L2B data files (one per orbit) in HDF format into the DART
observation sequence file makes several assumptions:</p>
<ol>
<li>All retrievals are labelled with a 10m height, in accordance
with the retrieval algorithm.</li>
<li>Only the highest-ranked (by the MLE method) solution is
desired.</li>
<li>Only the WVCs with a wvc_quality_flag of <b>zero</b> are
desired.</li>
<li>The mission specification of a wind speed rms error of 2 ms
(for winds less than 20 m/s) and 10% for windspeeds between 20 and
30 m/s can be extended to all winds with a qc flag of zero.</li>
<li>The mission specification of an error in direction of 20
degrees rms is applicable to all retrieved directions.</li>
<li>All retrievals with wind speeds less than 1.0 are not
used.</li>
<li>The above error characterstics can be simplified when deriving
the horizontal wind components (i.e. U,V). <b>Note :</b> this may
or may not be a good assumption, and efforts to assimilate the
speed and direction directly are under way.</li>
</ol>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<hr>
<h2>DATA SOURCES</h2>
<p>The NASA Jet Propulsion Laboratory (JPL) <a href=
"http://winds.jpl.nasa.gov/imagesAnim/quikscat.cfm">data
repository</a> has a collection of animations and data sets from
this instrument. In keeping with NASA tradition, these data are in
HDF format (specifically, HDF4), so if you want to read these files
directly, you will need to install the HDF4 libraries (which can be
downloaded from <a href=
"http://www.hdfgroup.org/products/hdf4/">http://www.hdfgroup.org/products/hdf4/</a>)</p>
<p>If you go to the official JPL data distribution site <a href=
"http://podaac.jpl.nasa.gov/DATA_CATALOG/quikscatinfo.html">http://podaac.jpl.nasa.gov/DATA_CATALOG/quikscatinfo.html</a>,
we are using the product labelled <b>L2B OWV 25km Swath</b>. They
are organized in folders by day ... with each orbit (each
revolution) in one compressed file. There are 14 revolutions per
day. The conversion to DART observation sequence format is done on
each revolution, multiple revolutions may be combined 'after the
fact' by any <em class="program">obs_sequence_tool</em> in the
<em class="file">work</em> directory of any model.</p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<hr>
<h2>PROGRAMS</h2>
<p>There are several programs that are distributed from the JPL
www-site, <a href=
"ftp://podaac.jpl.nasa.gov/pub/ocean_wind/quikscat/L2B/sw/">ftp://podaac.jpl.nasa.gov/pub/ocean_wind/quikscat/L2B/sw/</a>;
we specifically started from the Fortran file <a href=
"ftp://podaac.jpl.nasa.gov/pub/ocean_wind/quikscat/L2B/sw/FORTRAN/read_qscat2b.f">
read_qscat2b.f</a> and modified it to be called as a subroutine to
make it more similar to the rest of the DART framework. The
original <em class="file">Makefile</em> and <em class=
"file">read_qscat2b.f</em> are included in the DART distribution in
the <em class="file">DART/observations/quikscat</em> directory. You
will have to modify the <em class="file">Makefile</em> to build the
executable.</p>
<h3>convert_L2b.f90</h3>
<p><em class="program">convert_L2b</em> is the executable that
reads the HDF files distributed by JPL. <em class=
"file">DART/observations/quikscat/work</em> has the expected
<em class="file">mkmf_convert_L2b</em> and <em class=
"file">path_names_convert_L2b</em> files and compiles the
executable in the typical DART fashion - with one exception. The
location of the HDF (and possible dependencies) installation must
be conveyed to the <em class="program">mkmf</em> build mechanism.
Since this information is not required by the rest of DART, it made
sense (to me) to isolate it in the <em class=
"file">mkmf_convert_L2b</em> script. <b>It will be necessary to
modify the <em class="file">mkmf_convert_L2b</em> script to be able
to build <em class="program">convert_L2b</em></b>. In particular,
you will have to change the two lines specifying the location of
the HDF (and probably the JPG) libraries. The rest of the script
should require little, if any, modification.</p>
<div class="routine">set JPGDIR = <em class=
"input">/contrib/jpeg-6b_gnu-4.1.2-64</em><br>
set HDFDIR = <em class=
"input">/contrib/hdf-4.2r4_gnu-4.1.2-64</em><br></div>
<p>There are a lot of observations in every QuikSCAT orbit.
Consequently, the observation sequence files are pretty large -
particularly if you use the ASCII format. Using the binary format
(i.e. <em class="input">obs_sequence_nml:write_binary_obs_sequence
= .true.</em>) will result in observation sequence files that are
about <em>half</em> the size of the ASCII format.</p>
<p>Since there are about 14 QuikSCAT orbits per day, it may be
useful to convert individual orbits to an observation sequence file
and then concatenate multiple observation sequence files into one
file per day. This may be trivially accomplished with the
<em class="program">obs_sequence_tool</em> program in any
<em class="file">model/xxxx/work</em> directory. Be sure to include
the <em class=
"code">'../../../obs_def/obs_def_QuikSCAT_mod.f90'</em> string in
<em class="code">input.nml&amp;preprocess_nml:input_files</em> when
you run <em class="program">preprocess</em>.</p>
<h3>obs_to_table.f90, plot_wind_vectors.m</h3>
<p><em class=
"file">DART/diagnostics/threed_sphere/obs_to_table.f90</em> is a
potentially useful tool. You can run the observation sequence files
through this filter to come up with a 'XYZ'-like file that can be
readily plotted with <em class=
"file">DART/diagnostics/matlab/plot_wind_vectors.m</em>.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. We adhere to the F90 standard of starting a
namelist with an ampersand '&amp;' and terminating with a slash '/'
for all our namelist input. Character strings that contain a '/'
must be enclosed in quotes to prevent them from prematurely
terminating the namelist. The following values are the defaults for
these namelist items.</p>
<div class="namelist">
<pre>
&amp;convert_L2b_nml
   l2b_file = '',
   datadir = '.',
   outputdir = '.',
   lon1 = 0.0, 
   lon2 = 360.0, 
   lat1 = -90.0, 
   lat2 = 90.0,
   along_track_thin = 0,
   cross_track_thin = 0
 /
</pre></div>
<br>
<div>
<p>It is possible to restrict the output observation sequence to
contain data from a region of interest throught the use of the
namelist parameters. If you need a region that spans the Prime
Meridian lon1 can be a larger number than lon2, for example, a
region from 300 E to 40 E and 60 S to 30 S (some of the South
Atlantic), would be <em class="input">lon1 = 300, lon2 = 40, lat1 =
-60, lat2 = -30</em>.</p>
<table border="0" cellpadding="3" width="100%" summary=
'QuikSCAT namelist description'>
<thead align="left">
<tr>
<th>Contents</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>l2b_file</td>
<td>character(len=128)</td>
<td>name of the HDF file to read - NOT including the directory,
e.g. QS_S2B44444.20080021548</td>
</tr>
<tr>
<td>datadir</td>
<td>character(len=128)</td>
<td>the directory containing the HDF files</td>
</tr>
<tr>
<td>outputdir</td>
<td>character(len=128)</td>
<td>the directory for the output observation sequence files.</td>
</tr>
<tr>
<td>lon1</td>
<td>real(r4)</td>
<td>the West-most longitude of interest in degrees. [0.0, 360]</td>
</tr>
<tr>
<td>lon2</td>
<td>real(r4)</td>
<td>the East-most longitude of interest in degrees. [0.0, 360]</td>
</tr>
<tr>
<td>lat1</td>
<td>real(r4)</td>
<td>the South-most latitude of interest in degrees. [-90.0,
90.0]</td>
</tr>
<tr>
<td>lat2</td>
<td>real(r8)</td>
<td>the North-most latitude of interest in degrees. [-90.0,
90.0]</td>
</tr>
<tr>
<td>along_track_thin</td>
<td>integer</td>
<td>provides ability to thin the data by keeping only every Nth
row. e.g. 3 == keep every 3rd row.</td>
</tr>
<tr>
<td>cross_track_thin</td>
<td>integer</td>
<td>provides ability to thin the data by keeping only every Nth
wind vector cell in a particular row. e.g. 5 == keep every 5th
cell.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
 <a name="KnownBugs" id="KnownBugs"></a>
<hr>
<h2>KNOWN BUGS</h2>
<p>There are no known bugs at this time.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<ol>
<li>There is one bit of error-checking that did not survive the
conversion from F77 to F90. I need to restore the check that the
HDF file being read is a 'Level 2B' product.</li>
<li>There is a lot of error-checking that is not being done. I need
to bulletproof the code more.</li>
<li>We need namelist options to select something other than the
highest-ranked ambiguity.</li>
<li>We need namelist options to select more QC flags - not just the
ones with the 'perfect' QC value of 0</li>
<li>Add an option to leave the observations as speed and direction
instead of converting them to U,V components. This is a natural
implementation of the instrument error characteristics. However, it
would require writing a specialized forward operator in order to
assimilate obs of this type (speed, direction), and there is still
a numerical problem with trying to do the statistics required
during the assimilation of a cyclic direction value.</li>
</ol>
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
