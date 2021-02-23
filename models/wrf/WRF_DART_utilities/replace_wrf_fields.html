<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program replace_wrf_fields</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">replace_wrf_fields</em></h1>
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
<p>Program to copy various fields from one WRF netCDF file to
another.</p>
<p>There are many existing utilities to process netCDF files, i.e.
the NCO operators and NCL scripts, which have more functionality
than this program. The only purpose for having this one is that it
is a standalone program with no prerequisites or dependencies other
than the netCDF libraries. If you already have other tools
available they can do the same functions that this program
does.</p>
<p>This program copies the given data fields from the input file to
the output file, failing if their sizes, shapes, or data types do
not match exactly. The expected use is to copy fields which are
updated by the WRF program but are not part of the DART state
vector, for example, sea surface temperature or soil fields. After
DART has updated the WRF restart <em class="file">wrfinput_d01</em>
file, this program can be used to update other fields in the file
before running the model.</p>
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
&amp;replace_wrf_fields_nml
   fieldnames = 'SST',
   fieldlist_file = '',
   fail_on_missing_field = .true.
   debug = .false.,
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
<td>fieldnames</td>
<td>character(len=129) (:)</td>
<td>An array of ASCII field names to be copied from the input
netCDF file to the output netCDF file. The names must match
exactly, and the size and shape of the data must be the same in the
input and output files for the data to be copied. If the field
names are set here, the fieldlist_file item must be ' '.</td>
</tr>
<tr>
<td>fieldlist_file</td>
<td>character(len=129)</td>
<td>An alternative to an explicit list of field names to copy. This
is a single string, the name of a file which contains a single
field name, one per line. If this option is set, the fieldnames
namelist item must be ' '.</td>
</tr>
<tr>
<td>fail_on_missing_field</td>
<td>logical</td>
<td>If any fields in the input list are not found in either the
input or output netcdf files, fail if this is set to true. If
false, a warning message will be printed but execution will
continue.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>If true, print out debugging messages about which fields are
found in the input and output files.</td>
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
<li>input namelist ; <em class="file">input.nml</em></li>
<li>Input - output WRF state netCDF files; <em class=
"file">wrfinput_d01, wrfinput_d02, ...</em></li>
<li>fieldlist_file (if specified in namelist)</li>
</ul>
<h3>File formats</h3>
<p>This utility works on any pair of netCDF files, doing a simple
read and copy from one to the other.</p>
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
<td valign="top">replace_wrf_fields</td>
<!-- message -->
<td valign="top">Usage: echo infile.nc outfile.nc |
./replace_wrf_fields</td>
<!-- comment -->
<td valign="top">The program did not read 2 filenames from the
console.</td>
</tr>
<tr><!-- routine -->
<td valign="top">replace_wrf_fields</td>
<!-- message -->
<td valign="top">cannot specify both fieldnames and
fieldlist_file</td>
<!-- comment -->
<td valign="top">In the namelist you must either specify an
explicit list of fieldnames to copy between the files, or give a
single filename which contains the list of field names. You cannot
specify both.</td>
</tr>
<tr><!-- routine -->
<td valign="top">replace_wrf_fields</td>
<!-- message -->
<td valign="top"><i>field</i> not found in input/output file</td>
<!-- comment -->
<td valign="top">If 'fail_on_missing_field' is true in the namelist
and a field is not found in either the input or output file.</td>
</tr>
<tr><!-- routine -->
<td valign="top">replace_wrf_fields</td>
<!-- message -->
<td valign="top"><i>field</i> does not match</td>
<!-- comment -->
<td valign="top">If the input and output files have different
sizes, number of dimensions, or data types, the program cannot copy
the data.</td>
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
