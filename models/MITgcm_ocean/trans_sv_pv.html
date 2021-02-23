<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program trans_sv_pv</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<center><a href="#Modules">MODULES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a></center>
<h1>PROGRAM <em class="program">trans_sv_pv</em></h1>
<!-- version tag follows, do not edit -->
<p></p>
<p><em class="program">trans_sv_pv</em> is responsible for
converting a DART 'initial conditions' file to a set of model
'snapshot' files and appropriate namelist files: <em class=
"file">data.cal</em> and <em class="file">data</em>. This is easier
than the reverse process because the DART initial conditions file
have a header that contains the valid time for the accompanying
state. This same header also has the 'advance-to' time. <em class=
"program">trans_sv_pv</em> uses this information to write out
appropriate <em class="unix">&amp;CAL_NML</em> and <em class=
"unix">&amp;PARM03</em> namelists in <em class=
"file">data.cal.DART</em> and <em class="file">data.DART</em>,
respectively. The rest of the information in <em class=
"file">data</em> is preserved, so it is possible to simply replace
<em class="file">data</em> with the new <em class=
"file">data.DART</em>.<br>
<br>
The input filename is hardwired to that expected by <em class=
"program">filter</em> and the output filenames are able to be
renamed into those defined by the <em class=
"file">data</em><em class="unix">&amp;PARM05</em> namelist
specifying the filenames to use to cold-start the ocean model. The
output filename is comprised of 4 parts: the variable name, the
startDate_1 component (YYYYMMDD), the startDate_2 component
(HHMMSS), and the extension (.data for the data and .meta for the
metadata). The startDate_1 and startDate_2 pieces are identical in
format to that used by identically named variables in the
<em class="file">data.cal</em><em class="unix">&amp;CAL_NML</em>
namelist.</p>
<h2>Usage</h2>
<p>There must be several input files in the current working
directory; most of these are required by the <em class=
"unix">model_mod</em> interface. The input filename is hardwired to
<em class="file">assim_model_state_ic</em>. Assuming the time tag
in the input file is set to 06Z 23 July 1996, this
example creates output files named<br>
<em class="file">S.19960723.060000.[data,meta]</em><br>
<em class="file">T.19960723.060000.[data,meta]</em><br>
<em class="file">U.19960723.060000.[data,meta]</em><br>
<em class="file">V.19960723.060000.[data,meta]</em><br>
<em class="file">Eta.19960723.060000.[data,meta]</em><br>
<em class="file">data.cal.DART</em>, and<br>
<em class="file">data.DART</em></p>
<div class="unix">mv some_DART_ics_input_file
assim_model_state_ic<br>
./trans_sv_pv<br>
cp data.cal.DART data.cal<br>
cp data.DART data</div>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
model_mod
assim_model_mod
time_manager_mod
</pre>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>This program has no namelist of its own, but some of the
underlying modules require namelists to be read, even if the values
are not used. To avoid duplication and, possibly, some
inconsistency in the documentation; only a list of the required
namelists is provided - with a hyperlink to the full documentation
for each namelist.</p>
<table border="0" cellpadding="10" width="100%">
<tr>
<th align="left">Namelist</th>
<th align="left">Primary Purpose</th>
</tr>
<tr>
<td valign="top"><a href=
"../../assimilation_code/modules/utilities/utilities_mod.html#Namelist">
utilities_nml</a></td>
<td>set the termination level and file name for the run-time
log</td>
</tr>
<tr>
<td valign="top"><a href=
"model_mod.html#namelist_cal_nml">CAL_NML</a></td>
<td>must be read, values are not used. The <em class=
"file">data.cal.DART</em> file has an updated namelist to be used
for the model advance.</td>
</tr>
<tr>
<td valign="top"><a href=
"model_mod.html#namelist_parm03">PARM03</a></td>
<td>must be read, values are not used, The <em class=
"file">data.DART</em> is an 'identical' version of <em class=
"file">data</em> with the exception of the <em class=
"unix">PARM03</em> namelist. The parameters <em class=
"unix">endTime</em>, <em class="unix">dumpFreq</em>, and <em class=
"unix">taveFreq</em> reflect the amount of time needed to advance
the model. The parameter <em class="unix">startTime</em> is set to
0.0, which is required to force the model to read the startup files
specified by <em class="unix">PARM05</em></td>
</tr>
<tr>
<td valign="top"><a href=
"model_mod.html#namelist_parm04">PARM04</a></td>
<td>ocean model grid parameters, read - never changed.</td>
</tr>
</table>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<hr>
<h2>FILES</h2>
<ul>
<li>input namelist files: <em class="file">data, data.cal,
input.nml</em></li>
<li>output namelist files: <em class="file">data.cal.DART,
data.DART</em></li>
<li>input data file: <em class=
"file">assim_model_state_ic</em></li>
<li>output data files: <em class=
"file">[S,T,U,V,Eta].YYYYMMDD.HHMMSS.[data,meta]</em></li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>none</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<p>There are no error conditions specific to <em class=
"program">trans_sv_pv</em>.</p>
<h2>KNOWN BUGS</h2>
<p>There are no known bugs.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<p>I may put in a routine to change the <em class=
"unix">PARM05</em> namelist parameters defining the input file
names. The hack in <em class="program">advance_model.csh</em> to
grep out the filenames from the namelist and rename the files at
the shell level is ugly.<br>
<br>
Feel free to suggest improvements.</p>
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
