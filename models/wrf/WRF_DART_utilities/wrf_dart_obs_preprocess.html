<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program wrf_dart_obs_preprocess</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">wrf_dart_obs_preprocess</em></h1>
<table border="0" summary="dart header" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>Program to preprocess observations, with specific knowledge of
the WRF domain.</p>
<p>This program will exclude all observations outside of the given
WRF domain. There are options to exclude or increase the error
values of obs close to the domain boundaries. The program can
superob (average) aircraft and satellite wind obs if they are too
dense.</p>
<p>This program can read up to 9 additional obs_seq files and merge
their data in with the basic obs_sequence file which is the main
input.</p>
<p>This program can reject surface observations if the elevation
encoded in the observation is too different from the wrf surface
elevation.</p>
<p>This program can exclude observations above a specified height
or pressure.</p>
<p>This program can overwrite the incoming Data QC value with
another.</p>
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
&amp;wrf_obs_preproc_nml

  file_name_input          = 'obs_seq.old'
  file_name_output         = 'obs_seq.new'
  
  sonde_extra              = 'obs_seq.rawin'
  land_sfc_extra           = 'obs_seq.land_sfc'
  metar_extra              = 'obs_seq.metar'
  marine_sfc_extra         = 'obs_seq.marine'
  sat_wind_extra           = 'obs_seq.satwnd'
  profiler_extra           = 'obs_seq.profiler'
  gpsro_extra              = 'obs_seq.gpsro'
  acars_extra              = 'obs_seq.acars'
  trop_cyclone_extra       = 'obs_seq.tc'
  
  overwrite_obs_time       = .false.  
  
  obs_boundary             = 0.0
  increase_bdy_error       = .false.  
  maxobsfac                = 2.5   
  obsdistbdy               = 15.0  
  
  sfc_elevation_check      = .false.  
  sfc_elevation_tol        = 300.0  
  obs_pressure_top         = 0.0  
  obs_height_top           = 2.0e10  
  
  include_sig_data         = .true.   
  tc_sonde_radii           = -1.0  
  
  superob_aircraft         = .false.  
  aircraft_horiz_int       = 36.0  
  aircraft_pres_int        = 2500.0  
  
  superob_sat_winds        = .false.    
  sat_wind_horiz_int       = 100.0   
  sat_wind_pres_int        = 2500.0  
  
  overwrite_ncep_satwnd_qc = .false.    
  overwrite_ncep_sfc_qc    = .false.  
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
<td colspan="3"><strong>Generic parameters:</strong></td>
</tr>
<tr>
<td>file_name_input</td>
<td>character(len=129)</td>
<td>The input obs_seq file.</td>
</tr>
<tr>
<td>file_name_output</td>
<td>character(len=129)</td>
<td>The output obs_seq file.</td>
</tr>
<tr>
<td>sonde_extra, land_sfc_extra, metar_extra, marine_sfc_extra,
marine_sfc_extra, sat_wind_extra, profiler_extra, gpsro_extra,
acars_extra, trop_cyclone_extra</td>
<td>character(len=129)</td>
<td>The names of additional input obs_seq files, which if they
exist, will be merged in with the obs from the <em class=
"file">file_name_input</em> obs_seq file. If the files do not
exist, they are silently ignored without error.</td>
</tr>
<tr>
<td>overwrite_obs_time</td>
<td>logical</td>
<td>If true, replace the incoming observation time with the
analysis time. Not recommended.</td>
</tr>
<tr>
<td colspan="3"><strong>Boundary-specific parameters:</strong></td>
</tr>
<tr>
<td>obs_boundary</td>
<td>real(r8)</td>
<td>Number of grid points around domain boundary which will be
considered the new extent of the domain. Observations outside this
smaller area will be excluded.</td>
</tr>
<tr>
<td>increase_bdy_error</td>
<td>logical</td>
<td>If true, observations near the domain boundary will have their
observation error increased by <em class=
"code">maxobsfac</em>.</td>
</tr>
<tr>
<td>maxobsfac</td>
<td>real(r8)</td>
<td>If <em class="code">increase_bdy_error</em> is true, multiply
the error by a ramped factor. This item sets the maximum
error.</td>
</tr>
<tr>
<td>obsdistbdy</td>
<td>real(r8)</td>
<td>If <em class="code">increase_bdy_error</em> is true, this
defines the region around the boundary (in number of grid points)
where the observation error values will be altered. This is ramped,
so when you reach the innermost points the change in observation
error is 0.0.</td>
</tr>
<tr>
<td colspan="3"><strong>Parameters to reduce observation count
:</strong></td>
</tr>
<tr>
<td>sfc_elevation_check</td>
<td>logical</td>
<td>If true, check the height of surface observations against the
surface height in the model.</td>
</tr>
<tr>
<td>sfc_elevation_tol</td>
<td>real(r8)</td>
<td>If <em class="code">sfc_elevation_check</em> is true, the
maximum difference between the elevation of a surface observation
and the model surface height, in meters. If the difference is
larger than this value, the observation is excluded.</td>
</tr>
<tr>
<td>obs_pressure_top</td>
<td>real(r8)</td>
<td>Observations with a vertical coordinate in pressure which are
located above this pressure level (i.e. the obs vertical value is
smaller than the given pressure) will be excluded.</td>
</tr>
<tr>
<td>obs_height_top</td>
<td>real(r8)</td>
<td>Observations with a vertical coordinate in height which are
located above this height value (i.e. the obs vertical value is
larger than the given height) will be excluded.</td>
</tr>
<tr>
<td colspan="3"><strong>Radio/Rawinsonde-specific parameters
:</strong></td>
</tr>
<tr>
<td>include_sig_data</td>
<td>logical</td>
<td>If true, include significant level data from radiosondes.</td>
</tr>
<tr>
<td>tc_sonde_radii</td>
<td>real(r8)</td>
<td>If greater than 0.0 remove any sonde observations closer than
this distance in Kilometers to the center of a Tropical
Cyclone.</td>
</tr>
<tr>
<td colspan="3"><strong>Aircraft-specific parameters
:</strong></td>
</tr>
<tr>
<td>superob_aircraft</td>
<td>logical</td>
<td>If true, average all aircraft observations within the given
radius and output only a single observation. Any observation that
is used in computing a superob observation is removed from the list
and is not used in any other superob computation.</td>
</tr>
<tr>
<td>aircraft_horiz_int</td>
<td>real(r8)</td>
<td>If <em class="code">superob_aircraft</em> is true, the
horizontal distance in Kilometers which defines the superob area.
All other unused aircraft observations within this radius will be
averaged with the current observation.</td>
</tr>
<tr>
<td>aircraft_vert_int</td>
<td>real(r8)</td>
<td>If <em class="code">superob_aircraft</em> is true, the vertical
distance in Pascals which defines the maximum separation for
including an observation in the superob computation.</td>
</tr>
<tr>
<td colspan="3"><strong>Satellite Wind-specific parameters
:</strong></td>
</tr>
<tr>
<td>superob_sat_winds</td>
<td>logical</td>
<td>If true, average all sat_wind observations within the given
radius and output only a single observation. Any observation that
is used in computing a superob observation is removed from the list
and is not used in any other superob computation.</td>
</tr>
<tr>
<td>sat_wind_horiz_int</td>
<td>real(r8)</td>
<td>If <em class="code">superob_sat_winds</em> is true, the
horizontal distance in Kilometers which defines the superob area.
All other unused sat_wind observations within this radius will be
averaged with the current observation.</td>
</tr>
<tr>
<td>sat_wind_vert_int</td>
<td>real(r8)</td>
<td>If <em class="code">superob_sat_winds</em> is true, the
vertical distance in Pascals which defines the maximum separation
for including an observation in the superob computation.</td>
</tr>
<tr>
<td>overwrite_ncep_satwnd_qc</td>
<td>logical</td>
<td>If true, replace the incoming Data QC value in satellite wind
observations with 2.0.</td>
</tr>
<tr>
<td colspan="3"><strong>Surface Observation-specific parameters
:</strong></td>
</tr>
<tr>
<td>overwrite_ncep_sfc_qc</td>
<td>logical</td>
<td>If true, replace the incoming Data QC value in surface
observations with 2.0.</td>
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
obs_sequence_mod
utilities_mod
obs_kind_mod
time_manager_mod
model_mod
netcdf
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>Input namelist ; <em class="file">input.nml</em></li>
<li>Input WRF state netCDF files; <em class="file">wrfinput_d01,
wrfinput_d02, ...</em></li>
<li>Input obs_seq files (as specified in namelist)</li>
<li>Output obs_seq file (as specified in namelist)</li>
</ul>
<h3>File formats</h3>
<p>This utility can read one or more obs_seq files and combine them
while doing the rest of the processing. It uses the standard DART
observation sequence file format.</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>Generously contributed by Ryan Torn.</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<thead>
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>wrf_dart_obs_preprocess</td>
<td></td>
<td></td>
</tr>
<tr>
<td>wrf_dart_obs_preprocess</td>
<td></td>
<td></td>
</tr>
<tr>
<td>wrf_dart_obs_preprocess</td>
<td></td>
<td></td>
</tr>
<tr>
<td>wrf_dart_obs_preprocess</td>
<td></td>
<td></td>
</tr>
</tbody>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none.</p>
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
