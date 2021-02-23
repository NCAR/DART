<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>SSUSI Data</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>SSUSI F16 EDR-DSK format to observation sequence
converters</h1>
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
<p>The Special Sensor Ultraviolet Spectrographic Imager <a href=
"http://http://ssusi.jhuapl.edu/">SSUSI</a> is designed to remotely
sense the ionosphere and thermosphere. The following is repeated
from the SSUSI home page:</p>
<blockquote><em class="italic">Overview<br>
Beginning in 2003, the Defense Meteorological Satellite Program
(DMSP) satellites began carrying the SSUSI instrument - a
combination of spectrographic imaging and photometric systems
designed to remotely sense the ionosphere and thermosphere.<br>
<br>
The long term focus of the SSUSI program is to provide data
concerning the upper atmospheric response to the sun over the
changing conditions of the solar cycle. Data collected by SSUSI
instrument can help identify structure in the equatorial and polar
regions.<br>
<br>
Mission<br>
SSUSI was designed for the DMSP Block 5D-3 satellites. These
satellites are placed into nearly polar, sun-synchronous orbits at
an altitude of about 850 km. SSUSI is a remote-sensing instrument
which measures ultraviolet (UV) emissions in five different
wavelength bands from the Earth's upper atmosphere. SSUSI is
mounted on a nadir-looking panel of the satellite. The multicolor
images from SSUSI cover the visible Earth disk from horizon to
horizon and the anti-sunward limb up to an altitude of
approximately 520 km.<br>
<br>
The UV images and the derived environmental data provide the Air
Force Weather Agency (Offutt Air Force Base, Bellevue, NE) with
near real-time information that can be utilized in a number of
applications, such as maintenance of high frequency (HF)
communication links and related systems and assessment of the
environmental hazard to astronauts on the Space
Station.</em></blockquote>
<p><em class="program">convert_f16_edr_dsk.f90</em> will extract
the ON2 observations from the F16 "edr-dsk" format files and create
DART observation sequence files. There is one additional
preprocessing step before the edr-dsk files may be converted.<br>
<br>
The ON2_UNCERTAINTY variable in the netcdf files have IEEE NaN
values, but none of the required metadata to interpret them
correctly. These 2 lines will add the required attributes so that
NaNs are replaced with a fill value that can be queried and
manipulated. Since the ON2_UNCERTAINTY is a standard deviation, it
is sufficient to make the fill value negative. <a href-=
"#KnownBugs">See the section on Known Bugs</a></p>
<div class="unix">
<pre>
ncatted -a _FillValue,ON2_UNCERTAINTY,o,f,NaN        input_file.nc
ncatted -a _FillValue,ON2_UNCERTAINTY,m,f,-1.0       input_file.nc</pre></div>
<!--==================================================================-->
<a name="DataSources" id="DataSources"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>DATA SOURCES</h2>
<p><a href=
"http://ssusi.jhuapl.edu/data_products">http://ssusi.jhuapl.edu/data_products</a></p>
<p>Please read their <a href=
"http://ssusi.jhuapl.edu/home_data_usage">data usage</a>
policy.</p>
<!--==================================================================-->
<a name="Programs" id="Programs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PROGRAMS</h2>
<p><em class=
"file">DART/observations/SSUSI/convert_f16_edr_dsk.f90</em> will
extract ON2 data from the distribution files and create DART
observation sequence (obs_seq) files. Build it in the <em class=
"file">SSUSI/work</em> directory by running the <em class=
"program">./quickbuild.csh</em> script located there. In addition
to the converters, the <em class="file">advance_time</em> and
<em class="file">obs_sequence_tool</em> utilities will be
built.</p>
<p>An example data file is in the <em class="file">data</em>
directory. An example scripts for adding the required metadata to
the ON2_UNCERTAINTY variable in the <em class=
"file">shell_scripts</em> directory. These are <em>NOT</em>
intended to be turnkey scripts; they will certainly need to be
customized for your use. There are comments at the top of the
scripts saying what options they include, and should be commented
enough to indicate where changes will be likely to need to be
made.</p>
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERRORS</h2>
<p>The code for setting observation error variances is using fixed
values, and we are not certain if they are correct. Incoming QC
values larger than 0 are suspect, but it is not clear if they
really signal unusable values or whether there are some codes we
should accept.</p>
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
<a name="KnownBugs" id="KnownBugs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>KNOWN BUGS</h2>
<p>The netCDF files - as distributed - have NaN values to indicate
"MISSING". This makes it exceptionally hard to read or work with,
as almost everything will core dump when trying to perform any math
with NaNs. <em class="file">convert_f16_edr_dsk.f90</em> tries to
count how many values are missing. If the NaN has not been replaced
with a numerically valid MISSING value, the following FATAL ERROR
is generated (by the Intel compiler, with debug and traceback
enabled):</p>
<pre>
 set_nml_output Echo NML values to log file only
 Trying to open namelist log dart_log.nml
forrtl: error (65): floating invalid
Image              PC                Routine            Line        Source             
convert_f16_edr_d  000000000051717D  MAIN__                    143  convert_f16_edr_dsk.f90
convert_f16_edr_d  0000000000409B3C  Unknown               Unknown  Unknown
libc.so.6          0000003101E1ED5D  Unknown               Unknown  Unknown
convert_f16_edr_d  0000000000409A39  Unknown               Unknown  Unknown
Abort (core dumped)
</pre>
<p>The solution is to replace the NaN values with a viable MISSING
value using the <em class=
"program">shell_scripts/netcdf_manip.csh</em> script. It relies on
the netCDF Operators, freely available from  <a href=
"http://nco.sourceforge.net">http://nco.sourceforge.net</a></p>
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
