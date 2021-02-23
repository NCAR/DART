<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>GTSPP Observations</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>GTSPP Observations</h1>
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
<a href="#Modules">MODULES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">FUTURE PLANS</a> / <a href="#Legalese">TERMS
OF USE</a>
<h2>Overview</h2>
<p>GTSPP (Global Temperature-Salinity Profile Program) data
measures vertical profiles of ocean temperature and salinity. The
<a href="http://www.nodc.noaa.gov/GTSPP/index.html">GTPSS home
page</a> has detailed information about the repository,
observations, and datasets. The programs in this directory convert
from the netcdf files found in the repository into DART observation
sequence (obs_seq) file format.</p>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>DATA SOURCES</h2>
<p>Data from the GTSPP can be downloaded interactively from
<a href="http://www.nodc.noaa.gov/cgi-bin/gtspp/gtsppform01.cgi">here</a>.
It is delivered in <a href=
"http://www.unidata.ucar.edu/software/netcdf">netCDF</a> file
format, one vertical profile per netCDF file.</p>
<p>Currently each vertical profile is stored in a separate file, so
converting a months's worth of observations involves downloading
many individual files. The converter program can take a list of
input files, so it is easy to collect a month of observations
together into a single output file with one execution of the
converter program.</p>
<p>The units in the source file are degrees C for temperature, g/kg
for salinity, and so far we have not found any error information
(not quality control, but observation instrument error values).
There is probably instrument source information encoded in these
files, but so far we don't have the key. The quality control values
are read and only those with a QC of 1 are retained.</p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PROGRAMS</h2>
<p>The data is distributed in <a href=
"http://www.unidata.ucar.edu/software/netcdf">netCDF</a> file
format. DART requires all observations to be in a proprietary
format often called DART "obs_seq" format. The files in this
directory, a combination of C shell scripts and a Fortran source
executable, do this data conversion.</p>
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
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
&amp;gtspp_to_obs_nml
   gtspp_netcdf_file     = '1234567.nc'
   gtspp_netcdf_filelist = 'gtspp_to_obs_filelist'
   gtspp_out_file        = 'obs_seq.gtspp'
   avg_obs_per_file      = 500
   debug                 = .false.
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
<td>gtspp_netcdf_file</td>
<td>character(len=128)</td>
<td>The input filename when converting a single profile. Only one
of the two file or filelist items can have a valid value, so to use
the single filename set the list name 'gtspp_netcdf_filelist' to
the empty string (' ').</td>
</tr>
<tr>
<td>gtspp_netcdf_filelist</td>
<td>character(len=128)</td>
<td>To convert a series of profiles in a single execution create a
text file which contains each input file, in ascii, one filename
per line. Set this item to the name of that file, and set
'gtspp_netcdf_file' to the empty string (' ').</td>
</tr>
<tr>
<td>gtspp_out_file</td>
<td>character(len=128)</td>
<td>The output file to be created. To be compatible with earlier
versions of this program, if this file already exists it will be
read in and the new data will be inserted into that file.</td>
</tr>
<tr>
<td>avg_obs_per_file</td>
<td>integer</td>
<td>The code needs an upper limit on the number of observations
generated by this program. It can be larger than the actual number
of observations converted. The total number of obs is computed by
multiplying this number by the number of input files. If you get an
error because there is no more room to add observations to the
output file, increase this number.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>If true, output more debugging messages.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
time_manager_mod
utilities_mod
location_mod
obs_sequence_mod
obs_def_mod
obs_def_ocean_mod
obs_kind_mod
netcdf
</pre>
<!--==================================================================-->
<a name="KnownBugs" id="KnownBugs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>KNOWN BUGS</h2>
<p>Does not have correct code for setting observation error
variance yet. Also, not sure if the incoming data qc is strict
enough.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none</p>
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
