<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>WOD Observations</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>WOD Observations</h1>
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
<p>The WOD (World Ocean Database) data is a collection of data from
various sources, combined into a single format with uniform
treatment. The <a href=
"http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html">WOD 2009
page</a> has detailed information about the repository,
observations, and datasets. The programs in this directory convert
from the packed ASCII files found in the repository into DART
observation sequence (obs_seq) file format.<br>
<br>
There are 2 sets of available files - the raw observations, and the
observations binned onto standard levels. The recommended datasets
are the ones on standard levels. The raw data can be very dense in
the vertical and are not truly independent observations. This leads
to too much certainty in the updated values during the
assimilation.</p>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>DATA SOURCES</h2>
<p>Data from the WOD09 can be downloaded interactively from links
on <a href="http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html">this
page</a>. One suggestion is to pick the 'sorted by year link' and
download the files (you can select multiple years at once) for each
data type for the years of interest. Make sure to select the
standard level versions of each dataset.</p>
<p>UCAR/NCAR users with access to the DSS data repository can
download WOD09 files from <a href=
"http://dss.ucar.edu/datazone/dsszone/ds285.0/#WOD09">here</a>. A
UCAR DSS userid is required to access this page. The files to use
are named "yearly_*_STD.tar".</p>
<p>Requested citation if you use this data:</p>
<pre>
Johnson, D.R., T.P. Boyer, H.E. Garcia, R.A. Locarnini, O.K. Baranova, and M.M. Zweng, 
2009. World Ocean Database 2009 Documentation. Edited by Sydney Levitus. NODC 
Internal Report 20, NOAA Printing Office, Silver Spring, MD, 175 pp.  
Available at http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html. 
</pre>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PROGRAMS</h2>
<p>The data is distributed in a specialized packed ASCII format. In
this directory is a program called <em class="file">wodFOR.f</em>
which is an example reader program to print out data values from
the files. The program <em class="file">wod_to_obs</em> converts
these packed ASCII files into DART obs_sequence files.</p>
<p>As with most other DART directories, the <em class=
"file">work</em> directory contains a <em class=
"file">quickbuild.csh</em> script to build all necessary
executables.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
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
&amp;wod_to_obs_nml
   wod_input_file       =  'XBTS2005',
   wod_input_filelist   =  '',
   wod_out_file         =  'obs_seq.wod',
   avg_obs_per_file     =  500000,
   debug                =  .false.,
   timedebug            =  .false.,
   print_qc_summary     =  .true.,
   max_casts            =  -1,
   no_output_file       =  .false.,
   print_every_nth_cast =  -1,
   temperature_error    =  0.5,
   salinity_error       =  0.5, 
 /
! temperature error is in degrees C, salinity error in g/kg.
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
<td>wod_input_file</td>
<td>character(len=128)</td>
<td>The input filename when converting a single file. Only one of
the two namelist items that specify input files can have a valid
value, so to use a single filename set the list name
'wod_input_filelist' to the empty string (' ').</td>
</tr>
<tr>
<td>wod_input_filelist</td>
<td>character(len=128)</td>
<td>To convert one or more files in a single execution create a
text file which contains each input filename, in ascii, one
filename per line. Set this item to the name of that file, and set
'wod_input_file' to the empty string (' ').</td>
</tr>
<tr>
<td>wod_out_file</td>
<td>character(len=128)</td>
<td>The output file to be created. Note that unlike earlier
versions of some converters, this program will overwrite an
existing output file instead of appending to it. The risk of
replicated observations, which are difficult to detect since most
of the contents are floating point numbers, outweighed the possible
utility.</td>
</tr>
<tr>
<td>avg_obs_per_file</td>
<td>integer</td>
<td>The code needs an upper limit on the number of observations
generated by this program. It can be larger than the actual number
of observations converted. The total number of obs is computed by
multiplying this number by the number of input files. If you get an
error because there is no more room to add observations to the
output file, increase this number. Do not make this an unreasonably
huge number, however, since the code does preallocate space and
will be slow if the number of obs becomes very large.</td>
</tr>
<tr>
<td>print_every_nth_cast</td>
<td>integer</td>
<td>If a value greater than 0, the program will print a message
after processing every N casts. This allows the user to monitor the
progress of the conversion.</td>
</tr>
<tr>
<td>print_qc_summary</td>
<td>logical</td>
<td>If .TRUE. the program will print out a summary of the number of
casts which had a non-zero quality control values (current files
appear to use values of 1-9).</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>If .TRUE. the program will print out debugging
information.</td>
</tr>
<tr>
<td>timedebug</td>
<td>logical</td>
<td>If .TRUE. the program will print out specialized time-related
debugging information.</td>
</tr>
<tr>
<td>max_casts</td>
<td>integer</td>
<td>If a value greater than 0 the program will only convert at most
this number of casts from each input file. Generally only expected
to be useful for debugging. A negative value will convert all data
from the input file.</td>
</tr>
<tr>
<td>no_output_file</td>
<td>logical</td>
<td>If .TRUE. the converter will do all the work needed to convert
the observations, count the number of each category of QC values,
etc, but will not create the final obs_seq file. Can be useful if
checking an input file for problems, or for getting QC statistics
without waiting for a full output file to be constructed, which can
be slow for large numbers of obs. Only expected to be useful for
debugging.</td>
</tr>
<tr>
<td>temperature_error</td>
<td>real(r8)</td>
<td>The combined expected error of temperature observations from
all sources, including instrument error, model bias, and
representativeness error (e.g. larger or smaller grid box sizes
affecting expected accuracy), in degrees Centigrade. Values in
output file are error variance, which will be this value
squared.</td>
</tr>
<tr>
<td>salinity_error</td>
<td>real(r8)</td>
<td>The combined expected error of salinity observations from all
sources, including instrument error, model bias, and
representativeness error (e.g. larger or smaller grid box sizes
affecting expected accuracy) in g/kg (psu). Values in output file
are error variance, and use units of msu (kg/kg), so the numbers
will be this value / 1000.0, squared.</td>
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
</pre>
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERRORS and KNOWN BUGS</h2>
<p>The code for setting observation error variances is using fixed
values, and we are not certain if they are correct. Incoming QC
values larger than 0 are suspect, but it is not clear if they
really signal unusable values or whether there are some codes we
should accept.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>This converter is currently being used on WOD09 data, but the
standard files generally stop with early 2009 data. There are
subsequent additional new obs files available from the download
site.</p>
<p>The fractional-time field, and sometimes the day-of-month field
in a small percentage of the obs have bad values. The program
currently discards these obs, but it may be possible to recover the
original good day number and/or time of day. There is a subroutine
at the end of the <em class="file">wod_to_obs.f90</em> file which
contains all the reject/accept/correction information for the year,
month, day, time fields. To accept or correct the times on more
obs, edit this subroutine and make the necessary changes.</p>
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
