<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>program create_real_obs</title>
<link rel="stylesheet" type="text/css" href=
"../../../../docs/html/doc.css" />
<link href="../../../../docs/images/dart.ico" rel=
"shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM create_real_obs</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src=
"../../../../docs/images/Dartboard7.png" alt="DART project logo"
height="70" /></td>
<td>Jump to <a href="../../../../docs/index.html">DART
Documentation Main Index</a></td>
</tr>
</table>
<a href="#Instructions">INSTRUCTIONS</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Translating NCEP BUFR files into DART obs_seq.out files (input
file to filter) is a 2 stage process. The first stage uses NCEP
software to translate the BUFR file into an "intermediate" text
file. This is described in <a href=
"../prep_bufr/prep_bufr.html">prep_bufr</a>. The second step is to
translate the intermediate files into an <em class=
"file">obs_seq.out</em> files, which is done by <em class=
"program">create_real_obs</em>, as described in this document.</p>
<p>This program provides a number of options to select several
observation types (radiosonde, aircraft, and satellite data, etc.)
and the DART observation variables (U, V, T, Q, Ps) which are
specified in its optional namelist interface <a href=
"#Namelist"><em class="code">&amp;ncepobs_nml</em></a> which may be
read from file <em class="file">input.nml</em>.</p>
<!--==================================================================-->
<a name="Instructions" id="Instructions"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>INSTRUCTIONS</h2>
<ul>
<li>Go to DART/observations/NCEP/ascii_to_obs/work</li>
<li>Use <em class="program">quickbuild.csh</em> to compile all
executable programs in the directory. To rebuild just one program:
<ul>
<li>Use <em class="program">mkmf_create_real_obs</em> to generate
the makefile to compile <em class=
"file">create_real_obs.f90</em>.</li>
<li>Type <em class="unix">make</em> to get the executable.</li>
</ul>
</li>
<li>Make appropriate changes to the <em class=
"code">&amp;ncep_obs_nml</em> <a href="#Namelist">namelist</a> in
<em class="file">input.nml</em>, as follows.</li>
<li>run <em class="program">create_real_obs</em>.</li>
</ul>
<p>The selection of any combinations of the specific observation
fields (T, Q, U/V, and surface pressure) and types (radiosonde,
aircraft reports, or satellite wind, etc.) is made in the namelist
<em class="code">&amp;ncepobs_nml</em>. All the available
combinations of fields X types (i.e. ADPUPA and obs_U) will be
written to the obs_seq file. (You will be able to select which of
those to use during an assimilation in another namelist (<em class=
"code">assimilate_these_obs</em>, in <em class=
"code">&amp;obs_kind_nml</em>), so be sure to include all the
fields and types you might want.) You should change <em class=
"code">Obsbase</em> to the pathname of the decoded PREPBUFR text
data files. Be sure that <em class="code">daily_file</em> is set to
.TRUE. to create a single 24 hour file; .FALSE. converts input
files one-for-one with output files. The default action is to tag
each observation with the exact time it was taken and is the
recommended setting. However, if you want to bin the observations
in time, for example to do additional post-processing, the time on
all observations in the window can be overwritten and set to the
nearest synoptic time (e.g. 0Z, 6Z, 12Z, or 18Z), by setting
<em class="code">obs_time</em> to false.</p>
<p>Generally you will want to customize the namelist for your own
use. For example, here is a sample namelist:</p>
<pre>
&amp;ncepobs_nml
  year = 2007, 
  month = 3,
  day = 1,
  tot_days = 31,
  max_num = 700000,
  ObsBase = '../prep_bufr/work/temp_obs.'
  select_obs  = 1,
  ADPUPA = .true., 
  AIRCAR = .false.,  
  AIRCFT = .true., 
  SATEMP = .false., 
  SFCSHP = .false.,
  ADPSFC = .false.,  
  SATWND = .true., 
  obs_U  = .true., 
  obs_V  = .true.,
  obs_T  = .true.,
  obs_PS = .false.,
  obs_QV = .false.,
  daily_file = .true.
  obs_time = .true.,
/

&amp;obs_sequence_nml
  write_binary_obs_sequence = .false.  
/
</pre>
<p>This will produce daily observation sequence files for the
period of March 2007, which have the selected observation types and
fields; T, U, and V from radiosondes (ADPUPA) and aircraft
(AIRCFT). No surface pressure or specific humidity would appear in
the obs_seq files, nor observations from ACARS, satellites, and
surface stations. The output files look like "obs_seq200703dd",
with dd = 1,...,31.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;ncepobs_nml
   year       = 2003,
   month      = 1,
   day        = 1,
   tot_days   = 31,
   max_num    = 800000,
   select_obs = 0,
   ObsBase    = 'temp_obs.',
   ADPUPA     = .false., 
   AIRCAR     = .false., 
   AIRCFT     = .false., 
   SATEMP     = .false., 
   SFCSHP     = .false., 
   ADPSFC     = .false., 
   SATWND     = .false.,
   obs_U      = .false., 
   obs_V      = .false., 
   obs_T      = .false.,
   obs_PS     = .false.,
   obs_QV     = .false.,
   daily_file = .true.,
   obs_time   = .true.,
   lon1       =   0.0,
   lon2       = 360.0,
   lat1       = -90.0,
   lat2       =  90.0  
/
</pre></div>
<br />
<br />
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
<td>year, month, day</td>
<td>integer</td>
<td>Beginning year, month, day of the observation period.</td>
</tr>
<tr>
<td>tot_days</td>
<td>integer</td>
<td>Total days in the observation period. The converter cannot
cross month boundaries.</td>
</tr>
<tr>
<td>max_num</td>
<td>integer</td>
<td>Maximum observation number for the current one day files.</td>
</tr>
<tr>
<td>select_obs</td>
<td>integer</td>
<td>Controls whether to select a subset of observations from the
NCEP BUFR decoded daily ascii files.
<ul style="list-style: none;">
<li>0 = All observations are selected.</li>
<li>1 = Select observations using the logical parameters
below.</li>
</ul>
</td>
</tr>
<tr>
<td>daily_file</td>
<td>logical</td>
<td>Controls timespan of observations in each obs_seq file:
<ul>
<li>true = 24 hour spans (3:01Z to 3:00Z of the next day).
Filenames have the form obs_seqYYYYMMDD.</li>
<li>false = 6 hour spans (3:01Z to 9:00Z, 9:01Z to 15:00Z, 15:01Z
to 21:00Z, and 21:01Z to 3:00Z of the next day. Filenames have the
form obs_seqYYYYMMDDHH, where HH is 06, 12, 18, and 24.</li>
</ul>
</td>
</tr>
<tr>
<td>ObsBase</td>
<td>character(len=129)</td>
<td>Path that contains the decoded NCEP BUFR daily observation
files. To work with the example scripts this should be 'temp_obs.',
or if it includes a pathname then it should end with a
'/temp_obs.'</td>
</tr>
<tr>
<td>include_specific_humidity, include_relative_humidity,
include_dewpoint</td>
<td>logical</td>
<td>Controls which moisture observations are created. The default
is to create only specific humidity obs, but any, all, or none can
be requested. Set to .TRUE. to output that obs type, .FALSE. skips
it.</td>
</tr>
<tr>
<td>ADPUPA</td>
<td>logical</td>
<td>Select the NCEP type ADPUPA observations which includes land
and ship launched radiosondes and pibals as well as a few profile
dropsonde. This involves, at 00Z and 12Z, about 650 - 1000
stations, and at 06Z and 18Z (which are mostly pibals), about 150 -
400 stations.</td>
</tr>
<tr>
<td>AIRCFT</td>
<td>logical</td>
<td>Select the NCEP type AIRCFT observations, which includes
commercial, some military and reconnaissance reports. They are
flight level reports.</td>
</tr>
<tr>
<td>AIRCAR</td>
<td>logical</td>
<td>Select the NCEP type AIRCAR observations, which includes data
from aircraft takeoff and landings. Sometimes referred to as ACARS
obs.</td>
</tr>
<tr>
<td>SATEMP</td>
<td>logical</td>
<td>Select the NCEP type SATEMP observations, which includes NESDIS
ATOVS virtual temperature soundings.</td>
</tr>
<tr>
<td>SFCSHP</td>
<td>logical</td>
<td>Select the NCEP type SFCSHP observations, which includes
surface marine (ship, buoy, c-man) reports.</td>
</tr>
<tr>
<td>ADPSFC</td>
<td>logical</td>
<td>Select the NCEP type ADPSFC observations, which includes
surface land synoptic station reports.</td>
</tr>
<tr>
<td>SATWND</td>
<td>logical</td>
<td>Select the NCEP type SATWND observations, which includes winds
derived from satellite cloud drift analysis.</td>
</tr>
<tr>
<td>obs_U</td>
<td>logical</td>
<td>Select u-component of wind observations.</td>
</tr>
<tr>
<td>obs_V</td>
<td>logical</td>
<td>Select v-component of wind observations.</td>
</tr>
<tr>
<td>obs_T</td>
<td>logical</td>
<td>Select temperature observations.</td>
</tr>
<tr>
<td>obs_PS</td>
<td>logical</td>
<td>Select surface pressure observations.</td>
</tr>
<tr>
<td>obs_QV</td>
<td>logical</td>
<td>Select specific humidity observations.</td>
</tr>
<tr>
<td>lon1</td>
<td>real</td>
<td>Western longitude bound of observations to keep.</td>
</tr>
<tr>
<td>lon2</td>
<td>real</td>
<td>Eastern longitude bound of observations to keep. Can be less
than lon1 if region crosses prime meridian.</td>
</tr>
<tr>
<td>lat1</td>
<td>real</td>
<td>Lower latitude bound of observations to keep.</td>
</tr>
<tr>
<td>lat2</td>
<td>real</td>
<td>upper latitude bound of observations to keep.</td>
</tr>
<tr>
<td>obs_time</td>
<td>logical</td>
<td>If .true. use the full time in the input data. To force all
observation times in the output to the synoptic time (e.g. 0Z, 6Z,
12Z, or 18Z) set this to .false. (not recommended).</td>
</tr>
</tbody>
</table>
</div>
<br />
<br />
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
obs_utilities_mod
obs_sequence_mod
obs_kind_mod
obs_def_mod
assim_model_mod
model_mod
cov_cutoff_mod
location_mod
random_seq_mod
time_manager_mod
null_mpi_utilities_mod
real_obs_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FILES</h2>
<ul>
<li>path_names_create_real_obs; the list of modules used in the
compilation of create_real_obs.</li>
<li>temp_obs.yyyymmdd; (input) NCEP BUFR (decoded/intermediate)
observation file(s) Each one has 00Z of the next day on it.</li>
<li>input.nml; the namelist file used by create_real_obs.</li>
<li>obs_seqYYYYMMDD[HH]; (output) the obs_seq files used by
DART.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>REFERENCES</h2>
<ul>
<li>.../DART/observations/NCEP/prep_bufr/docs/* (NCEP text files
describing the BUFR files)</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top"> </td>
<!-- message -->
<td valign="top"> </td>
<!-- comment -->
<td valign="top"> </td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FUTURE PLANS</h2>
<p>Further development to get observations directly from
original<br />
(undecoded) NCEP BUFR files.</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
