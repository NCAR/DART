<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program MIDAS_to_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">MIDAS_to_obs</em></h1>
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
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<h3>MIDAS netCDF File to DART observation converter</h3>
<p>Alex Chartier (University of Bath, UK) is the point-of-contact
for this effort.</p>
<blockquote>"MIDAS runs in Matlab. The raw observations come from
GPS receivers as RINEX files, but we can't use them directly just
yet ... Currently, the 'slant' (satellite-to-receiver path)
observations are inverted by MIDAS to make vertical,
column-integrated 'observations' of plasma density."</blockquote>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>DATA SOURCES</h2>
<p>Contact Alex for MIDAS observations.<br>
<br>
Alex writes out netCDF files that may be converted to DART
observation sequence files. The netCDF files have a pretty simple
format.</p>
<pre>
netcdf Test {
dimensions:
        latitude = 5 ;
        longitude = 6 ;
        height = 30 ;
        time = UNLIMITED ; // (1 currently)
variables:
        double latitude(latitude) ;
                latitude:units = "degrees_north" ;
                latitude:long_name = "latitude" ;
                latitude:standard_name = "latitude" ;
        double longitude(longitude) ;
                longitude:units = "degrees_east" ;
                longitude:long_name = "longitude" ;
                longitude:standard_name = "longitude" ;
        double height(height) ;
                height:units = "metres" ;
                height:long_name = "height" ;
                height:standard_name = "height" ;
        double time(time) ;
                time:units = "Days since 1601-01-01" ;
                time:long_name = "Time (UT)" ;
                time:standard_name = "Time" ;
        double Ne(height, latitude, longitude) ;
                Ne:grid_mapping = "standard" ;
                Ne:units = "1E11 e/m^3" ;
                Ne:long_name = "electron density" ;
                Ne:coordinates = "latitude longitude" ;
        double TEC(time, latitude, longitude) ;
                TEC:grid_mapping = "standard" ;
                TEC:units = "1E16 e/m^2" ;
                TEC:long_name = "total electron content" ;
                TEC:coordinates = "latitude longitude" ;
        double Variance(time, latitude, longitude) ;
                Variance:grid_mapping = "standard" ;
                Variance:units = "1E16 e/m^2" ;
                Variance:long_name = "Variance of total electron content" ;
                Variance:coordinates = "latitude longitude" ;
                Variance:standard_name = "TEC variance" ;
// global attributes:
                :Conventions = "CF-1.5" ;
}</pre>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PROGRAMS</h2>
<p>The <em class="file">MIDAS_to_obs.f90</em> file is the source
code for the main converter program.<br>
<br>
To compile and test, go into the <em class="file">MIDAS/work</em>
subdirectory and run the <em class="program">quickbuild.csh</em>
script to build the converter and a couple of general purpose
utilities. The <a href=
"../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">
obs_sequence_tool</a> manipulates (i.e. combines, subsets) DART
observation files once they have been created. The default
observations supported are those defined in <a href=
"../../forward_operators/obs_def_upper_atm_mod.f90">observations/forward_operators/obs_def_upper_atm_mod.f90</a>.
If you need additional observation types, you will have to add the
appropriate <em class="file">obs_def_XXX_mod.f90</em> file to the
<em class="file">input.nml</em> <em class=
"code">&amp;preprocess_nml:input_files</em> variable and run
<em class="program">quickbuild.csh</em> again. It rebuilds the
table of supported observation types before compiling the source
code.</p>
<p> </p>
<!-- needed to make 'top' align correctly -->
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
&amp;MIDAS_to_obs_nml
   input_file    = 'infile.nc'
   obs_out_file  = 'obs_seq.out',
   verbose       = .false.
   /
</pre></div>
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
<td>input_file</td>
<td>character(len=256)</td>
<td>Name of the input netCDF MIDAS file to read.</td>
</tr>
<tr>
<td>obs_out_file</td>
<td>character(len=256)</td>
<td>Name of the output observation sequence file that is
created.</td>
</tr>
<tr>
<td>verbose</td>
<td>logical</td>
<td>Controls how much informational output is printed during a
conversion. <em class="code">.true.</em> the most amount of output.
<em class="code">.false.</em> the least amount of output.</td>
</tr>
</tbody>
</table>
<h3 class="indent1">Example</h3>
<pre>
&amp;MIDAS_to_obs_nml
   input_file    = '../data/Test.nc',
   obs_out_file  = 'obs_seq.out',
   verbose       = .TRUE.,
</pre>
<!--==================================================================-->
<!-- References.                                                      -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>References</h2>
<!-- ul>
<li><a href="http://MIDAS.hwr.arizona.edu">The MIDAS web page.</a></li>

</ul -->
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
<a name="KnownBugs" id="KnownBugs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>KNOWN BUGS</h2>
<p>none</p>
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
