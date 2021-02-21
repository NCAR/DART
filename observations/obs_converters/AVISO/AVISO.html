<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>Aviso+/CMEMS Observations</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>Aviso+/CMEMS Observations</h1>
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
<h3><em class="program">convert_aviso</em> was contributed by
Frederic Castruccio - Thanks Fred!</h3>
<p>This short description of the <em class=
"code">SEALEVEL_GLO_SLA_L3_REP_OBSERVATIONS_008_018</em> product is
repeated from the <strong>INFORMATION</strong> tab from the
<a href="http://marine.copernicus.eu/about-us/about-your-copernicus-marine-service/">
Copernicus Marine Environment Monitoring Service</a> online
catalogue (in April 2017).</p>
<blockquote>For the Global Ocean- Mono altimeter satellite
along-track sea surface heights computed with respect to a
twenty-year mean. Previously distributed by Aviso+, no change in
the scientific content. All the missions are homogenized with
respect to a reference mission which is currently Jason-2. This
product is computed with an optimal and centered computation time
window (6 weeks before and after the date). Two kinds of datasets
are proposed: filtered (nominal dataset) and
unfiltered.</blockquote>
The <em class="program">convert_aviso.f90</em> program is designed
to read a <a href=
"http://www.unidata.ucar.edu/software/netcdf">netCDF</a> file
containing the (Level 3) sea surface anomalies from any of the
following platforms: "<em class="code">Jason-1</em>", "<em class=
"code">Envisat</em>", or "<em class="code">Geosat Follow On</em>".
One of those platforms must be listed in the netCDF file global
attribute: <em class="code">platform</em><br>
<br>
The data files have names like:<br>
<em class=
"file">dt_global_j1_sla_vfec_20080101_20140106.nc</em>,<br>
<em class="file">dt_global_en_sla_vfec_20080101_20140106.nc</em>,
or<br>
<em class="file">dt_global_g2_sla_vfec_20080101_20140106.nc</em>;
corresponding to the <em class="code">Jason-1</em>, <em class=
"code">Envisat</em>, and the <em class="code">Geosat Follow On</em>
platforms.<br>
<br>
The DART observation TYPE corresponding to each of these platforms
are <em class="code">J1_SEA_SURFACE_ANOMALY</em>, <em class=
"code">EN_SEA_SURFACE_ANOMALY</em>, and <em class=
"code">GFO_SEA_SURFACE_ANOMALY</em>, respectively and are defined
in <a href=
"../../forward_operators/obs_def_ocean_mod.f90">obs_def_ocean_mod.f90</a>.<br>

<br>
Fred wrote a python script (<em class=
"program">shell_scripts/convert_aviso.py</em>) to repeatedly call
<em class="program">convert_aviso</em> and decided it was easiest
to simply provide the input file name as a command line argument
and always have the output file have the name <em class=
"file">obs_seq.aviso</em>. As such, there is no input namelist
specifically for these parameters, but other DART modules still
require run-time crontrol specified by <em class=
"file">input.nml</em>.
<p>After creating a large number of output observation sequence
files, it is usually necessary to consolidate the files and subset
them into files containing just the timeframe required for a single
assimilation. <strong>NOTE</strong>: the <em class=
"program">obs_sequence_tool</em> is constructed for just this
purpose.<br>
<br>
The <em class="program">shell_scripts/makedaily.sh</em> script
attempts to consolidate all the SLA observations and those that may
have been (separately) converted from the World Ocean Database into
24-hour segments centered at midnight GMT. You will have to modify
the <em class="program">makedaily.sh</em> script to suit your
filesystem and naming convention. It is provided as a starting
point.</p>
<p><strong>Reminder</strong>: (according to the data providers): In
order to compute Absolute Dynamic Topography, the Mean Dynamic
Topography (MDT) can be added. It is distributed by Aviso+ (
<a href=
"http://www.aviso.altimetry.fr/en/data/products/auxiliary-products/mdt.html">
http://www.aviso.altimetry.fr/en/data/products/auxiliary-products/mdt.html</a>
). Fred was using this product in assimilations with POP, so he
chose a different source for MDT - consistent with POP's
behavior.</p>
<p><!-- dummy space --></p>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>DATA SOURCES</h2>
<p>The Copernicus Marine and Environment Monitoring Service (CMEMS)
has taken over the processing and distribution of the Ssalto/Duacs
multimission altimeter products formerly administered by Aviso+.
After a registration process, the along-track sea level anomalies
(SLA) may be downloaded from <a href=
"http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&amp;view=details&amp;product_id=SEALEVEL_GLO_SLA_L3_REP_OBSERVATIONS_008_018">
http://marine.copernicus.eu/services-portfolio/access-to-products/</a>
- search for the <em class=
"code">SEALEVEL_GLO_SLA_L3_REP_OBSERVATIONS_008_018</em> if it does
not come up directly.</p>
<p><!-- dummy space --></p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PROGRAMS</h2>
<table border="0" cellpadding="10" width="100%" summary=
'program description'>
<tr>
<td valign="top"><em class="program">convert_aviso.f90</em></td>
<td valign="top">does the actual conversion from netCDF to a DART
observation sequence file, which may be ASCII or binary.</td>
</tr>
<tr>
<td valign="top"><em class=
"program">shell_scripts/convert_aviso.py</em></td>
<td valign="top">python script to convert a series of input files
and datestamp the output files.</td>
</tr>
<tr>
<td valign="top"><em class=
"program">shell_scripts/makedaily.sh</em></td>
<td valign="top">shell script to repeatedly call <em class=
"program">obs_sequence_tool</em> to consolidate multiple
observation sequence files into an observation sequence file that
has ALL the observations from ALL platforms in a single file.
<em class="program">makedaily.sh</em> is capable of looping over
time ranges and creating observation sequences for each time
range.</td>
</tr>
</table>
<p><!-- do nothing spacer --></p>
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>There is no namelist for <em class="program">convert_aviso</em>,
but other namelists control aspects of the execution, namely
<em class=
"code">&amp;obs_sequence_nml:write_binary_obs_sequence</em>. see
<a href=
"../../../assimilation_code/modules/observations/obs_sequence_mod.html">
obs_sequence_mod.html</a>.</p>
<p><!-- do nothing spacer --></p>
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
assimilation_code/location/threed_sphere/location_mod.f90
assimilation_code/modules/assimilation/assim_model_mod.f90
assimilation_code/modules/io/dart_time_io_mod.f90
assimilation_code/modules/observations/obs_kind_mod.f90
assimilation_code/modules/observations/obs_sequence_mod.f90
assimilation_code/modules/utilities/ensemble_manager_mod.f90
assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
assimilation_code/modules/utilities/random_seq_mod.f90
assimilation_code/modules/utilities/sort_mod.f90
assimilation_code/modules/utilities/time_manager_mod.f90
assimilation_code/modules/utilities/types_mod.f90
assimilation_code/modules/utilities/utilities_mod.f90
models/template/model_mod.f90
observations/forward_operators/obs_def_mod.f90
observations/obs_converters/AVISO/convert_aviso.f90
observations/obs_converters/utilities/obs_utilities_mod.f90
</pre>
<p><!-- dummy space --></p>
<!--==================================================================-->
<a name="Errors" id="Errors"></a> <a name="KnownBugs" id=
"KnownBugs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>KNOWN BUGS</h2>
<p>none</p>
<p><!-- dummy space --></p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none</p>
<p><!-- dummy space --></p>
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
