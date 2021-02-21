<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program closest_member_tool</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">closest_member_tool</em></h1>
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
<p>Utility program to compare the ensemble mean to an ensemble of
restart files, which can now be run in parallel. The program prints
out a sorted order of which members are 'closest' to the mean,
where the method used to determine 'close' is selectable by
namelist option. It also creates a file with a single number or
character string in it, for ease in scripting, which identifies the
closest member.</p>
<p>The ensemble mean is computed from the input ensemble. The
difference is computed point by point across the ensemble members.
There is an option to restrict the computation to just a subset of
the entire state vector by listing one or more generic quantities.
In this case, only state vector items matching one of these
quantities will contribute to the total difference value.</p>
<p>Available methods are:</p>
<dl>
<dt>1 - simple absolute difference:</dt>
<dd>The absolute value of the difference between each item in the
mean vector and the corresponding item in each ensemble member,
accumulated over the entire state vector.</dd>
<dt>2 - normalized absolute difference:</dt>
<dd>The absolute value of the difference between each item in the
mean vector and the corresponding item in each ensemble member
normalized by the mean value, accumulated over the entire state
vector.</dd>
<dt>3 - simple RMS difference:</dt>
<dd>The square root of the accumulated sum of the square of the
difference between each item in the mean vector and the
corresponding item in each ensemble member.</dd>
<dt>4 - normalized RMS difference:</dt>
<dd>The square root of the accumulated sum of the square of the
normalized difference between each item in the mean vector and the
corresponding item in each ensemble member.</dd>
</dl>
<p>This program could be used to select one or more ensemble
members to run a free model forecast forward in time after the
assimilation is finished. Each member is an equally likely
representation of the model state. Using the ensemble mean may not
be the best choice since the mean may not have self-consistent
fine-scale structures in the data.</p>
<p>In addition to printing out data about all members to both the
console and to the dart log file, this program creates a single
output file containing information about the closest member. If the
input restart data is in a single file, the output file
'closest_restart' contains a single number which is the ensemble
member number. If the input restart data is in separate files, the
output file contains the full filename of the closest member, e.g.
'filter_restart.0004' if member 4 is closest. For scripting the
contents of this file can be used to copy the corresponding member
data and convert it to the model input format for a free forecast,
for example.</p>
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
&amp;closest_member_tool_nml
   input_restart_files     = ''
   input_restart_file_list = ''      
   output_file_name        = 'closest_restart'
   ens_size                = 1
   difference_method       = 4      
   use_only_qtys           = ''
   single_restart_file_in  = .false.      
  /
</pre></div>
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
<td>input_restart_files</td>
<td>character(len=256),dimension(ens_size x num_domains)</td>
<td>An array of filenames each containing a list DART restart
data.</td>
</tr>
<tr>
<td>input_restart_file_list</td>
<td>character(len=256),dimension(num_domains)</td>
<td>A file containing a list of filenames for DART restart data,
one for each domain.</td>
</tr>
<tr>
<td>output_file_name</td>
<td>character(len=256)</td>
<td>This is a file containing the member number that is closest to
the ensemble mean.</td>
</tr>
<tr>
<td>ens_size</td>
<td>integer</td>
<td>Total number of ensemble members.</td>
</tr>
<tr>
<td>difference_method</td>
<td>integer</td>
<td>Select which method is used to compute 'distance' from mean:
<ul style="list-style: none;">
<li>1 = simple absolute difference</li>
<li>2 = absolute difference normalized by the mean</li>
<li>3 = simple RMS difference</li>
<li>4 = RMS of the normalized difference</li>
</ul>
</td>
</tr>
<tr>
<td>use_only_quantities</td>
<td>character(len=32)</td>
<td>If unspecified, all items in the state vector contribute to the
total difference. If one or more quantities are listed here, only
items in the state vector of these quantities contribute to the
total difference. These are the generic quantities, such as
QTY_TEMPERATURE, QTY_U_WIND_COMPONENT, QTY_DENSITY, etc. and not
specific types like RADIOSONDE_TEMPERATURE. Consult the model
interface code to determine which possible quantities are returned
by the <a href=
"../../../models/template/model_mod.html#get_state_meta_data">get_state_meta_data()</a>
routine.</td>
</tr>
<tr>
<td>single_restart_file_in</td>
<td>logical</td>
<td><strong>Not supported yet.</strong> Contact dart@ucar.edu if
you are interested in using this tool with files that contain all
ensemble members in a single file.</td>
</tr>
</tbody>
</table>
</div>
<p>Below is an example of a typical namelist for the
closest_member_tool.</p>
<div class="namelist">
<pre>
&amp;closest_member_tool_nml
   input_restart_files     = ''
   input_restart_file_list = 'restart_list.txt'      
   output_file_name        = 'closest_restart.txt'
   ens_size                = 3
   single_restart_file_in  = .false.      
   difference_method       = 4      
   use_only_qtys           = ''
  /
</pre></div>
<p>where <em class="file">restart_list.txt</em> contains</p>
<pre>
cam_restart_0001.nc
cam_restart_0002.nc
cam_restart_0003.nc
</pre>
<p>Currently <em class="code">single_restart_file_in</em> is not
supported. This is typically used for simpler models that have
built in model advances such as <em class=
"program">lorenz_96</em>.</p>
<br>
<br>
<br>
<br>
<!--==================================================================-->
 
<!-- A NAME="Modules"></A> @>todo modules have too many dependencies to make this meaninful
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>MODULES USED</H2>
<PRE>
types_mod
time_manager_mod
utilities_mod
sort_mod
assim_model_mod
mpi_utilities_mod
</PRE== >

<!====================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>inputfile.####.nc (list of restarts to find closest member)
-or-</li>
<li><em class="file">restart_list.txt</em> (a file containing a
list of restart files) and,</li>
<li><em class="file">input.nml</em></li>
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
<td valign="top">closest_member_tool</td>
<!-- message -->
<td valign="top">Invalid method number</td>
<!-- comment -->
<td valign="top">Values 1-4 are supported</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_variables</td>
<!-- message -->
<td valign="top">&lt;some variable name&gt;: NetCDF: Start+count
exceeds dimension bound</td>
<!-- comment -->
<td valign="top">The variable in the model definition is not
conformable with the variable in the restart file.</td>
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
<p>Add check to see that the model template variables are
conformable with the variables in the files being read.</p>
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
