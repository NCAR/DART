<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program littler_tf_dart</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">littler_tf_dart</em></h1>
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
<p>Programs to convert littler data files into DART observation
sequence files, and vice versa. The capability of the program is
limited to wind and temperature from radiosondes.</p>
<p>The littler data files do not contain observation errors. The
observation errors are in a separate file called <em class=
"file">obserr.txt</em>. The littler file generated here has to be
preprocessed by the program <em class="file">3dvar_obs.exe</em>
before beeing ingested in the WRF 3D-Var system.</p>
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
obs_sequence_mod
obs_def_mod
obs_kind_mod
location/threed_sphere/location_mod
time_manager_mod
utilities_mod
</pre>
<h2>MODULES INDIRECTLY USED</h2>
<pre>
assim_model_mod
models/wrf/model_mod
models/wrf/module_map_utils
random_seq_mod
</pre>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a><br>
<hr>
<br>
<h2>NAMELIST</h2>
<p>The program does not have its own namelist. However, an
<em class="file">input.nml</em> file is required for the modules
used by the program.</p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<hr>
<h2>FILES</h2>
<ul>
<li>input namelist ; <em class="file">input.nml</em></li>
<li>Input - output observation files; <em class=
"file">obs_seq.out</em> and <em class="file">little-r.dat</em></li>
<li>Input - output littler observation error files ; <em class=
"file">obserr.txt</em></li>
</ul>
<h3>File formats</h3>
<p>If there are no observation error at a particular pressure
level, the default value of -1 is written in <em class=
"file">obserr.txt</em>.</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<hr>
<h2>REFERENCES</h2>
<ul>
<li><a href="http://www.mmm.ucar.edu/wrf/WG4/">3DVAR GROUP
PAGE</a></li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">littler_tf_dart</td>
<!-- message -->
<td valign="top">Did not get all obs</td>
<!-- comment -->
<td valign="top">There is observation before dart time (0, 0) or
beyond day 200000 in the Gregorian calendar in the <em class=
"file">obs_seq.out</em> file.</td>
</tr>
<tr><!-- routine -->
<td valign="top">littler_tf_dart</td>
<!-- message -->
<td valign="top">No vertical coordinate.</td>
<!-- comment -->
<td valign="top">Only vertical pressure coordinate is
supported.</td>
</tr>
<tr><!-- routine -->
<td valign="top">littler_tf_dart</td>
<!-- message -->
<td valign="top">Error when reading or writing header of
sounding.</td>
<!-- comment -->
<td valign="top">Problem with littler file.</td>
</tr>
<tr><!-- routine -->
<td valign="top">littler_tf_dart</td>
<!-- message -->
<td valign="top">Error when reading or writing sounding.</td>
<!-- comment -->
<td valign="top">Problem with littler file.</td>
</tr>
<tr><!-- routine -->
<td valign="top">littler_tf_dart</td>
<!-- message -->
<td valign="top">Error when reading or writing footer of the
sounding.</td>
<!-- comment -->
<td valign="top">Problem with littler file.</td>
</tr>
<tr><!-- routine -->
<td valign="top">intplin, intplog</td>
<!-- message -->
<td valign="top">arrays xx and yy must have same size</td>
<!-- comment -->
<td valign="top">Each x value needs a corresponding y value.</td>
</tr>
<tr><!-- routine -->
<td valign="top">intplin, intplog</td>
<!-- message -->
<td valign="top">bad value in yy</td>
<!-- comment -->
<td valign="top">It is assumed that the input yy contains rms
errors and should not be negative. That should be traced back to
the file <em class="file">obserr.txt</em>.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<ol>
<li>Develop the program to support all observations in DART and WRF
3D-Var.</li>
<li>Use the preprocessor to include the observation list provided
by obs_kind_mod.</li>
</ol>
<!--==================================================================-->
<!-- Declare all private entities.                                    -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">call set_str_date(timestring,
dart_time)</em>
<pre>
type(time_type),   intent(in)  :: <em class="code"> dart_time </em>
character(len=20), intent(out) :: <em class=
"code"> timestring </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a dart_time (seconds, days), returns date as
bbbbbbyyyymmddhhmmss, where b is a blank space.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">call set_dart_time(tstring,
dart_time)</em>
<pre>
character(len=20), intent(in)  :: <em class="code"> tstring </em>
type(time_type),   intent(out) :: <em class="code"> dart_time </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a date as bbbbbbyyyymmddhhmmss, where b is a blank space,
returns the dart_time (seconds, days).</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">call StoreObsErr(obs_err_var,
pres, plevel, nlev, obs_err_std)</em>
<pre>
integer,  intent(in)    :: <em class="code"> nlev, pres </em>
real(r8), intent(in)    :: <em class="code"> obs_err_var </em>
integer,  intent(in)    :: <em class="code"> plevel(nlev) </em>
real(r8), intent(inout) :: <em class=
"code"> obs_err_std(nlev) </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>If the incoming pres corresponds exactly to a pressure level in
plevel, then transfers the incoming obs_err_var into the array
obs_err_std at the corresponding level.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">level_index =
GetClosestLevel(ilev, vlev, nlev)</em>
<pre>
integer,  intent(in) :: <em class="code"> nlev, ilev </em>
integer,  intent(in) :: <em class="code"> vlev(nlev) </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the index of the closest level in vlev to the incoming
ilev.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">call READ_OBSERR(filein,
platform, sensor_name, err, nlevels)</em>
<pre>
CHARACTER (LEN=80), intent(in)  :: <em class="code"> filein </em>
CHARACTER (LEN=80), intent(in)  :: <em class="code"> platform </em>
CHARACTER (LEN=80), intent(in   :: <em class=
"code"> sensor_name </em>
INTEGER,            intent(in)  :: <em class="code"> nlevels </em>
REAL(r8),           intent(out) :: <em class=
"code"> err(nlevels) </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Read observational error on pressure levels (in hPa) from the
incoming filein and store the result in the array err. It is
assumed that filein has the same format as WRF 3D-Var <em class=
"file">obserr.txt</em> file. It reads observational error for a
specific platform (e.g. RAOBS) and a specific sensor (e.g. WIND
SENSOR ERRORS).</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">f_obstype =
obstype(line)</em>
<pre>
CHARACTER (LEN= 80), intent(in) :: <em class="code"> line </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Read in a line the string present after keyword 'BOGUS', which
should be the sensor name.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">f_sensor = sensor(line)</em>
<pre>
CHARACTER (LEN= 80), intent(in) :: <em class="code"> line </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Read in a line the string present after numbers, which should be
the platform name.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">val = intplin(x,xx,yy)</em>
<pre>
INTEGER,  DIMENSION (:), intent(in) :: <em class="code"> xx </em>
REAL(r8), DIMENSION (:), intent(in) :: <em class="code"> yy </em>
REAL(r8),                intent(in) :: <em class="code"> x </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Do a linear interpolation.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">val = intplog(x,xx,yy)</em>
<pre>
INTEGER,  DIMENSION (:), intent(in) :: <em class="code"> xx </em>
REAL(r8), DIMENSION (:), intent(in) :: <em class="code"> yy </em>
REAL(r8),                intent(in) :: <em class="code"> x </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Do a log-linear interpolation.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<br>
<div class="routine"><em class="call">index = locate(x,xx)</em>
<pre>
INTEGER, DIMENSION (:), intent(in) :: <em class="code"> xx </em>
REAL(r8),               intent(in) :: <em class="code"> x </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Return the index in xx such that xx(index) &lt; x &lt;
xx(index+1).</p>
</div>
<br>
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
