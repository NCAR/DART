<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program compare_states</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">compare_states</em></h1>
<table border="0" summary="" cellpadding="5">
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
<p>Utility program to compare fields in two NetCDF files and print
out the min and max values from each file and the min and max of
the differences between the two fields. The default is to compare
all numeric variables in the files, but specific variables can be
specified in the namelist or in a separate file. The two input
NetCDF filenames are read from the console or can be echo'd into
the standard input of the program.</p>
<p>If you want to restrict the comparison to only specific
variables in the files, specify the list of field names to compare
either in the namelist, or put a list of fields, one per line, in a
text file and specify the name of the text file. Only data arrays
can be compared, not character arrays, strings, or attribute
values.</p>
<p>Namelist interface <a href="#Namelist"><em class=
"code">&amp;compare_states_nml</em></a> must be read from file
<em class="file">input.nml</em>.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
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
&amp;compare_states_nml
   do_all_numeric_fields   = .true.
   fieldnames              = ''
   fieldlist_file          = ''
   fail_on_missing_field   = .true.
   only_report_differences = .true.
   debug                   = .false.
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
<td>do_all_numeric_fields</td>
<td>logical</td>
<td>If .true., all integer, float, and double variables in the
NetCDF files will have their values compared. If .false. the list
of specific variables to be compared must be given either directly
in the namelist in the <em class="code">fieldnames</em> item, or
else the field names must be listed in an ASCII file, one name per
line, and the name of that file is specified in <em class=
"code">fieldlist_file</em>.</td>
</tr>
<tr>
<td>fieldnames</td>
<td>character list</td>
<td>One or more names of arrays in the NetCDF files to be compared.
Only read if <em class="code">do_all_numeric_fields</em> is
.false.</td>
</tr>
<tr>
<td>fieldlist_file</td>
<td>character</td>
<td>Name of a text file containing the fieldnames, one per line. It
is an error to specify both the fieldnames namelist item and this
one. Only read if <em class="code">do_all_numeric_fields</em> is
.false.</td>
</tr>
<tr>
<td>fail_on_missing_field</td>
<td>logical</td>
<td>If .true. and any one of the field names is not found in both
files it is a fatal error. If .false. a message is printed about
the missing field but execution continues.</td>
</tr>
<tr>
<td>only_report_differences</td>
<td>logical</td>
<td>If .true. only print the name of the variable being tested;
skip printing the variable value min and max if the two files are
identical. If .false. print more details about both variables which
differ and varibles with the same values.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>If true print out debugging info.</td>
<td></td>
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
utilities_mod
parse_args_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>two NetCDF input files</li>
<li>compare_states.nml</li>
<li>field names text file (optionally)</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>none</li>
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
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">compare_states</td>
<!-- message -->
<td valign="top">Only use single process</td>
<!-- comment -->
<td valign="top">Only a single mpi process can be used with this
program</td>
</tr>
<tr><!-- routine -->
<td valign="top">compare_states</td>
<!-- message -->
<td valign="top">must specify data days and times</td>
<!-- comment -->
<td valign="top">If overwrite_data_time is true, the namelist must
include the new day and time.</td>
</tr>
<tr><!-- routine -->
<td valign="top">compare_states</td>
<!-- message -->
<td valign="top">output_is_model_advance_file must be true to set
advance time</td>
<!-- comment -->
<td valign="top">If overwrite_advance_time is true,
output_is_model_advance_file must also be true.</td>
</tr>
<tr><!-- routine -->
<td valign="top">compare_states</td>
<!-- message -->
<td valign="top">must specify advance days and times</td>
<!-- comment -->
<td valign="top">If overwrite_advance_time is true, the namelist
must include the new day and time.</td>
</tr>
<tr><!-- routine -->
<td valign="top">compare_states</td>
<!-- message -->
<td valign="top">overwrite_advance_time must be true if output file
has advance time</td>
<!-- comment -->
<td valign="top">If the incoming file does not have a model advance
time, the output cannot have one unless the user gives one in the
namelist, and sets overwrite_advance_time to true.</td>
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
