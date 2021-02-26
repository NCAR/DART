<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program obs_loop</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>program <em class="program">obs_loop</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Usage%20Notes">USAGE</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#Discussion">DISCUSSION</a> /
<a href="#Building">BUILDING</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>This program is a template that is intended to be modified by
the user to do any desired operations on an observation sequence
file.</p>
<!--=====================================================================-->
<!--===================== USAGE NOTES ===================================-->
<!--=====================================================================-->
<a name="Usage Notes"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Usage</h2>
<p>This program is intended to be used as a template to read in
observations from one obs_seq file and write them, optionally
modified in some way, to another obs_seq file. It can be compiled
and run as-is, but it simply makes an exact copy of the input
file.</p>
<p>There are comments in the code (search for <em class=
"code">MODIFY HERE</em> ) where you can test values, types, times,
error values, and either modify them or skip copying that
observation to the output.</p>
<p>There are build files in <em class=
"file">observations/utilities/oned</em> and <em class=
"file">observations/utilities/threed_sphere</em> to build the
<em class="program">obs_loop</em> program.</p>
<br>
<br>
<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->
 <a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Namelist</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;obs_loop_nml
   filename_in  = ''
   filename_out = '' 
   print_only   = .false.
   calendar     = 'Gregorian'
   /
</pre></div>
<br>
<br>
<p>Items in this namelist set the input and output files.</p>
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
<td>filename_in</td>
<td>character(len=256)</td>
<td>Observation sequence file to read</td>
</tr>
<tr>
<td>filename_out</td>
<td>character(len=256)</td>
<td>Observation sequence file to create and write. If this file
exists it will be overwritten.</td>
</tr>
<tr>
<td>print_only</td>
<td>logical</td>
<td>If .TRUE. then do the work but only print out information about
what would be written as output without actually creating the
output file.</td>
</tr>
<tr>
<td>calendar</td>
<td>character(len=32)</td>
<td>The string name of a valid DART calendar type. (See the
<a href="../../modules/utilities/time_manager_mod.html">time_manager_mod</a>
documentation for a list of valid types.) The setting here does not
change what is written to the output file; it only changes how the
date information is printed to the screen in the informational
messages.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--=====================================================================-->
<!--===================== DISCUSSION ====================================-->
<!--=====================================================================-->
 <a name="Discussion" id="Discussion"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Discussion</h2>
<p>See the documentation in the obs_kind and obs_def modules for
things you can query about an observation, and how to set
(overwrite) existing values.</p>
<br>
<br>
<!--==================================================================-->
 <a name="Building" id="Building"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>BUILDING</h2>
<p>There are build files in <em class=
"file">observations/utilities/oned</em> and <em class=
"file">observations/utilities/threed_sphere</em> to build the
<em class="program">obs_loop</em> program.</p>
<p>The <em class="file">preprocess</em> program must be built and
run first to define what set of observation types will be
supported. See the <a href=
"../../../assimilation_code/programs/preprocess/preprocess.html">preprocess
documentation</a> for more details on how to define the list and
run it. The <em class="code">&amp;preprocess_nml</em> namelist in
the <em class="file">input.nml</em> file must contain files with
definitions for the combined set of all observation types which
will be encountered over all input obs_seq files.</p>
<p>If you have observation types which are not part of the default
list in the &amp;preprocess_nml namelist, add them to the input.nml
file and then either run quickbuild.csh or make and run preprocess
and then make the obs_loop tool.</p>
<p>Usually the directories where executables are built will include
a <em class="file">quickbuild.csh</em> script which builds and runs
preprocess and then builds the rest of the executables by executing
all files with names starting with <em class="file">mkmf_</em>.</p>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Files</h2>
<table border="0">
<tr>
<th>filename</th>
<th>purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read the &amp;obs_loop_nml namelist</td>
</tr>
</table>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>References</h2>
<ol>
<li>none</li>
</ol>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Error Codes and Conditions</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">obs_loop</td>
<!-- message -->
<td valign="top"></td>
<!-- comment -->
<td valign="top"></td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_loop</td>
<!-- message -->
<td valign="top"></td>
<!-- comment -->
<td valign="top"></td>
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
<hr>
<h2>Future Plans</h2>
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
