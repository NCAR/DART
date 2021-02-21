<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>DART sst_to_obs, oi_sst_to_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">sst_to_obs, oi_sst_to_obs</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Decisions">DECISIONS</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>There are two gridded SST observation converters in this
directory, one for data from PODAAC, and one from NOAA/NCDC.
<em class="program">sst_to_obs</em> converts data from PODAAC and
has been used by Romain Escudier for regional studies with ROMS.
<em class="program">oi_sst_to_obs</em> converts data from NOAA/NCDC
and has been used by Fred Castruccio for global studies with
POP.</p>
<h3>sst_to_obs -- GHRSST to DART Observation Sequence
Converter</h3>
<p>These routines are designed to convert the <a href=
"https://podaac.jpl.nasa.gov/dataset/AVHRR_OI-NCEI-L4-GLOB-v2.0">GHRSST
Level 4 AVHRR_OI Global Blended Sea Surface Temperature Analysis
(GDS version 2) from NCEI data</a> distributed by the <a href=
"http://podaac.jpl.nasa.gov">Physical Oceanography Distributed
Active Archive Center</a>. Please remember to cite the data in your
publications, <a href=
"https://podaac.jpl.nasa.gov/dataset/AVHRR_OI-NCEI-L4-GLOB-v2.0">specific
instructions from PODAAC are available here.</a> This is an
example:</p>
<blockquote>National Centers for Environmental Information. 2016.
GHRSST Level 4 AVHRR_OI Global Blended Sea Surface Temperature
Analysis (GDS version 2) from NCEI. Ver. 2.0. PO.DAAC, CA, USA.
Dataset accessed [YYYY-MM-DD] at
http://dx.doi.org/10.5067/GHAAO-4BC02.</blockquote>
<p><strong>Many thanks to Romain Escudier (then at Rutgers) who did
the bulk of the work and graciously contributed his efforts to the
DART project.</strong> Romain gave us scripts and source code to
download the data from the PODAAC site, subset the global files to
a region of interest, and convert that subsetted file to a DART
observation sequence file. Those scripts and programs have been
only lightly modified to work with the Manhattan version of DART
and contain a bit more documentation.</p>
<p>The workflow is usually:</p>
<ol>
<li>compile the converters by running <em class=
"program">work/quickbuild.csh</em> in the usual way.</li>
<li>customize the <em class=
"file">shell_scripts/parameters_SST</em> resource file to specify
variables used by the rest of the scripting.</li>
<li>run <em class="program">shell_scripts/get_sst_ftp.sh</em> to
download the data from PODAAC.</li>
<li>provide a mask for the desired study area.</li>
<li>run <em class="program">shell_scripts/Prepare_SST.sh</em> to
subset the PODAAC data and create the DART observation sequence
files. Be aware that the <em class="program">Prepare_SST.sh</em>
modifies the <em class="file">shell_scripts/input.nml.template</em>
file and generates its own <em class="file">input.nml</em>.
<em class="file">work/input.nml</em> is not used.</li>
<li>combine all output files for the region and timeframe of
interest into one file using the <a href=
"../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">
obs_sequence_tool</a></li>
</ol>
<h3>Example:</h3>
<p>It is worth describing a small example. If you configure
<em class="program">get_sst_ftp.sh</em> to download the last two
days of 2010 and then specify the mask to subset for the
NorthWestAtlantic (NWA) and run <em class=
"program">Prepare_SST.sh</em> your directory structure should look
like the following:</p>
<pre>
<code>
0[1234] cheyenne6:/&lt;6&gt;obs_converters/SST
.
|-- ObsData
|   `-- SST
|       |-- ncfile
|       |   `-- 2010
|       |       |-- 20101230120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
|       |       `-- 20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
|       `-- nwaSST
|           `-- 2010
|               |-- 20101230120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc
|               `-- 20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc
|-- oi_sst_to_obs.f90
|-- oi_sst_to_obs.nml
|-- sst_to_obs.f90
|-- sst_to_obs.nml
|-- shell_scripts
|   |-- Prepare_SST.sh
|   |-- functions.sh
|   |-- get_sst_ftp.sh
|   |-- input.nml
|   |-- input.nml.template
|   |-- my_log.txt
|   |-- parameters_SST
|   `-- prepare_SST_file_NWA.sh
|-- masks
|   |-- Mask_NWA-NCDC-L4LRblend-GLOB-v01-fv02_0-AVHRR_OI.nc
|   `-- Mask_NWA120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
`-- work
    |-- Makefile
    |-- advance_time
    |-- input.nml
    |-- mkmf_advance_time
    |-- mkmf_obs_sequence_tool
    |-- mkmf_oi_sst_to_obs
    |-- mkmf_preprocess
    |-- mkmf_sst_to_obs
    |-- obs_sequence_tool
    |-- oi_sst_to_obs
    |-- path_names_advance_time
    |-- path_names_obs_sequence_tool
    |-- path_names_oi_sst_to_obs
    |-- path_names_preprocess
    |-- path_names_sst_to_obs
    |-- preprocess
    |-- quickbuild.csh
    `-- sst_to_obs
</code>
</pre>
<p>The location of the DART observation sequence files is specified
by <em class="file">parameter_SST</em>:<em class=
"code">DIR_OUT_DART</em>. That directory should contain the
following two files:</p>
<pre>
<strong>0[1236] cheyenne6:/&lt;6&gt;v2/Err30 &gt; ls -l</strong>
'total 7104
-rw-r--r-- 1 thoar p86850054 3626065 Jan 10 11:08 obs_seq.sst.20101230
-rw-r--r-- 1 thoar p86850054 3626065 Jan 10 11:08 obs_seq.sst.20101231
</pre>
<h2>oi_sst_to_obs -- NOAA/NCDC to DART Observation Sequence
Converter</h2>
<p><em class="program">oi_sst_to_obs</em> is designed to convert
the <a href=
"https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html">
NOAA High-resolution Blended Analysis: Daily Values using AVHRR
only</a> data. The global metadata of a typical file is shown
here:</p>
<pre>
<code>
:Conventions = "CF-1.5" ;
:title = "NOAA High-resolution Blended Analysis: Daily Values using AVHRR only" ;
:institution = "NOAA/NCDC" ;
:source = "NOAA/NCDC  ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/" ;
:comment = "Reynolds, et al., 2007:
     Daily High-Resolution-Blended Analyses for Sea Surface Temperature.
     J. Climate, 20, 5473-5496.
     Climatology is based on 1971-2000 OI.v2 SST, 
     Satellite data: Navy NOAA17 NOAA18 AVHRR, Ice data: NCEP ice." ;
:history = "Thu Aug 24 13:46:51 2017: ncatted -O -a References,global,d,, sst.day.mean.2004.v2.nc\n",
        "Version 1.0" ;
:references = "https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html" ;
:dataset_title = "NOAA Daily Optimum Interpolation Sea Surface Temperature" ;
</code>
</pre>
<p>The workflow is usually:</p>
<ol>
<li>compile the converters by running <em class=
"program">work/quickbuild.csh</em> in the usual way.</li>
<li><a href=
"https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html">
download the desired data.</a></li>
<li>customize the <em class="file">work/input.nml</em> file.</li>
<li>run <em class="program">work/oi_sst_to_obs</em> to create a
single DART observation sequence file.</li>
<li>combine all output files for the region and timeframe of
interest into one file using the <a href=
"../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">
obs_sequence_tool</a></li>
</ol>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>sst_to_obs NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;sst_to_obs_nml
   sst_netcdf_file     = '1234567.nc'
   sst_netcdf_filelist = 'sst_to_obs_filelist'
   sst_out_file        = 'obs_seq.sst'
   subsample_intv      = 1
   sst_rep_error       = 0.3
   debug               = .false.
   /
</pre></div>
<div>
<table border="0" cellspacing="10" width="100%" summary=
'sst_to_obs namelist description'>
<thead align="left">
<tr>
<th>Contents</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr><!--contents-->
<td>sst_netcdf_file</td>
<!--  type  -->
<td>character(len=256)</td>
<!--descript-->
<td>Name of the (usually subsetted) netcdf data file. This may be a
relative or absolute filename. If you run the scripts 'as is', this
will be something like:<br>
<em class=
"file">../ObsData/SST/nwaSST/2010/20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc</em></td>
</tr>
<tr><!--contents-->
<td>sst_netcdf_filelist</td>
<!--  type  -->
<td>character(len=256)</td>
<!--descript-->
<td>Name of the file that contains a list of (usually subsetted)
data files, one per line. <strong>You may not specify both
sst_netcdf_file AND sst_netcdf_filelist.</strong> One of them must
be empty.</td>
</tr>
<tr><!--contents-->
<td>sst_out_file</td>
<!--  type  -->
<td>character(len=256)</td>
<!--descript-->
<td>Name of the output observation sequence file.</td>
</tr>
<tr><!--contents-->
<td>subsample_intv</td>
<!--  type  -->
<td>integer</td>
<!--descript-->
<td>It is possible to 'thin' the observations. <em class=
"code">subsample_intv</em> allows one to take every Nth
observation.</td>
</tr>
<tr><!--contents-->
<td>sst_rep_error</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>In DART the observation error variance can be thought of as
having two components, an instrument error and a representativeness
error. In <em class="program">sst_to_obs</em> the instrument error
is specified in the netCDF file by the variable <em class=
"code">analysis_error</em>. The representativeness error is
specified by <em class="code">sst_rep_error</em>, which is
specified as a standard deviation. These two values are added
together and squared and used as the observation error variance.
<strong>Note:</strong> This algorithm maintains backwards
compatibility, but is technically not the right way to combine
these two quantities. If they both specified variance, adding them
together and then taking the square root would correctly specify a
standard deviation. Variances add, standard deviations do not.
Since the true observation error variance (in general) is not
known, we are content to live with an algorithm that produces
useful observation error variances. If your research comes to a
more definitive conclusion, please let us know.</td>
</tr>
<tr><!--contents-->
<td>debug</td>
<!--  type  -->
<td>logical</td>
<!--descript-->
<td>Print extra information during the <em class=
"program">sst_to_obs</em> execution.</td>
</tr>
</tbody>
</table>
</div>
<h2>oi_sst_to_obs NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;oi_sst_to_obs_nml
   input_file       = '1234567.nc'
   output_file_base = 'obs_seq.sst'
   subsample_intv   = 1
   sst_error_std    = 0.3
   debug            = .false.
   /
</pre></div>
<div>
<table border="0" cellspacing="10" width="100%" summary=
'oi_sst_to_obs namelist description'>
<thead align="left">
<tr>
<th>Contents</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr><!--contents-->
<td>input_file</td>
<!--  type  -->
<td>character(len=256)</td>
<!--descript-->
<td>Name of the input netcdf data file. This may be a relative or
absolute filename. If you run the scripts 'as is', this will be
something like:<br>
<em class=
"file">../ObsData/SST/nwaSST/2010/20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc</em></td>
</tr>
<tr><!--contents-->
<td>output_file_base</td>
<!--  type  -->
<td>character(len=256)</td>
<!--descript-->
<td>Partial filename for the output file. The date and time are
appended to <em class="code">output_file_base</em> to construct a
unique filename reflecting the time of the observations in the
file.</td>
</tr>
<tr><!--contents-->
<td>subsample_intv</td>
<!--  type  -->
<td>integer</td>
<!--descript-->
<td>It is possible to 'thin' the observations. <em class=
"code">subsample_intv</em> allows one to take every Nth
observation.</td>
</tr>
<tr><!--contents-->
<td>sst_error_std</td>
<!--  type  -->
<td>real</td>
<!--descript-->
<td>This is the total observation error standard deviation.</td>
</tr>
<tr><!--contents-->
<td>debug</td>
<!--  type  -->
<td>logical</td>
<!--descript-->
<td>Print extra information during the <em class=
"program">oi_sst_to_obs</em> execution.</td>
</tr>
</tbody>
</table>
</div>
<!--==================================================================-->
<a name="Decisions" id="Decisions"></a>
<hr>
<h2>DECISIONS YOU MIGHT NEED TO MAKE</h2>
<p>See the general discussion in the <a href=
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
<p>I do not believe <em class="program">sst_to_obs</em> will work
correctly if given multiple files in <em class=
"code">sst_netcdf_filelist</em>. The number of observation used to
declare the length of the output observation sequence is based on a
single file ... yet seems to be used by many. I have not tested
this configuration, since the scripting does not use the <em class=
"code">sst_netcdf_filelist</em> mechanism.</p>
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
