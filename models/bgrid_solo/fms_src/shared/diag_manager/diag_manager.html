<!--

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

-->
<!DOCTYPE html>
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>Module diag_manager_mod</title>
<link type="text/css" href=
"http://www.gfdl.noaa.gov/~fms/style/doc.css" rel="stylesheet">
<style type="text/css">
          .fixed {
            font-size:medium;
            font-family:monospace;
            border-style:none;
            border-width:0.1em;
            padding:0.1em;
            color:#663366;
          }
</style>
</head>
<body>
<a name="TOP" id="TOP"></a><font class="header" size="1"><a href=
"#PUBLIC%20INTERFACE">PUBLIC INTERFACE</a> ~ <a href=
"#PUBLIC%20DATA">PUBLIC DATA</a> ~ <a href=
"#PUBLIC%20ROUTINES">PUBLIC ROUTINES</a> ~ <a href=
"#NAMELIST">NAMELIST</a> ~ <a href=
"#DIAGNOSTIC%20FIELDS">DIAGNOSTIC FIELDS</a> ~ <a href=
"#ERROR%20MESSAGES">ERROR MESSAGES</a> ~ <a href=
"#REFERENCES">REFERENCES</a> ~ <a href="#NOTES">NOTES</a></font>
<hr>
<h1>module diag_manager_mod</h1>
<a name="OVERVIEW" id="OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text">&lt;TT&gt;diag_manager_mod&lt;/TT&gt; is a set of
simple calls for parallel diagnostics on distributed systems. It is
geared toward the writing of data in netCDF format.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<div>&lt;TT&gt;diag_manager_mod&lt;/TT&gt; provides a convenient
set of interfaces for writing data to disk. It is built upon the
parallel I/O interface <a href=
"http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/models/bgrid_solo/fms_src/shared/mpp/mpp_io.html">
&lt;TT&gt;mpp_io&lt;/TT&gt;</a>. A single group of calls to the
&lt;TT&gt;diag_manager_mod&lt;/TT&gt; interfaces provides data to
disk at any number of sampling and/or averaging intervals specified
at run-time. Run-time specification of diagnostics are input
through the diagnostics table, which is described in the <a href=
"models/bgrid_solo/fms_src/shared/diag_manager/diag_table_tk.html">diag_table_tk</a>
documentation.<br>
<br>
&lt;B&gt;Features of &lt;TT&gt;diag_manager_mod&lt;/TT&gt;
include:&lt;/B&gt; Simple, minimal API.<br>
Run-time choice of diagnostics.<br>
Self-describing files: comprehensive header information (metadata)
in the file itself.<br>
Strong parallel write performance.<br>
Integrated netCDF capability: <a href=
"http://www.unidata.ucar.edu/packages/netcdf/">netCDF</a> is a data
format widely used in the climate/weather modeling community.
netCDF is considered the principal medium of data storage for
&lt;TT&gt;diag_manager_mod&lt;/TT&gt;. Raw unformatted fortran I/O
capability is also available.<br>
Requires off-line post-processing: a tool for this purpose,
&lt;TT&gt;mppnccombine&lt;/TT&gt;, is available. GFDL users may use
&lt;TT&gt;~hnv/pub/mppnccombine&lt;/TT&gt;. Outside users may
obtain the source <a href=
"ftp://ftp.gfdl.gov/perm/hnv/mpp/mppnccombine.c">here</a>. It can
be compiled on any C compiler and linked with the netCDF library.
The program is free and is covered by the <a href=
"ftp://ftp.gfdl.gov/perm/hnv/mpp/LICENSE">GPL license</a>.</div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>
time_manager_mod<br>      mpp_io_mod<br>         fms_mod<br>   diag_axis_mod<br> diag_output_mod</pre></div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h2>PUBLIC INTERFACE</h2>
<div>
<pre>
use diag_manager_mod [, only:  send_data,<br>                               register_diag_field,<br>                               register_static_field,<br>                               diag_manager_end,<br>                               diag_manager_init,<br>                               get_base_time,<br>                               get_base_date,<br>                               need_data ]</pre>
<dl>
<dt><a href="#send_data">send_data</a>:</dt>
<dd>Send data over to output fields.</dd>
<dt><a href="#register_diag_field">register_diag_field</a>:</dt>
<dd>Register Diagnostic Field.</dd>
<dt><a href=
"#register_static_field">register_static_field</a>:</dt>
<dd>Register Static Field.</dd>
<dt><a href="#diag_manager_end">diag_manager_end</a>:</dt>
<dd>Exit Diagnostics Manager.</dd>
<dt><a href="#diag_manager_init">diag_manager_init</a>:</dt>
<dd>Initialize Diagnostics Manager.</dd>
<dt><a href="#get_base_time">get_base_time</a>:</dt>
<dd>Return base time for diagnostics.</dd>
<dt><a href="#get_base_date">get_base_date</a>:</dt>
<dd>Return base date for diagnostics.</dd>
<dt><a href="#need_data">need_data</a>:</dt>
<dd>Determine whether data is needed for the current model time
step.</dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h2>PUBLIC DATA</h2>
<!-- BEGIN PUBLIC DATA -->
<div>None.<br>
<br></div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li><a name="send_data" id="send_data"></a>
<h2>send_data</h2>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>send_data is overloaded for 1 to 3-d arrays. diag_field_id
corresponds to the id returned from a previous call to
register_diag_field. The field array is restricted to the
computational range of the array. Optional argument is_in can be
used to update sub-arrays of the entire field. Additionally, an
optional logical or real mask can be used to apply missing values
to the array. For the real mask, the mask is applied if the mask
value is less than 0.5. The weight array is currently not
implemented.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>diag_field_id   </tt></td>
<td><br>
   <span class="type">[integer]</span><br>
   <span class="type">[integer]</span><br>
   <span class="type">[integer]</span><br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>field   </tt></td>
<td><br>
   <span class="type">[real]</span><br>
   <span class="type">[real,
dimension(:)]</span><br>
   <span class="type">[real,
dimension(:,:)]</span><br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time   </tt></td>
<td><br>
   <span class="type">[time_type]</span><br>
   <span class="type">[time_type]</span><br>
   <span class="type">[time_type]</span><br>
   <span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="register_diag_field" id="register_diag_field"></a>
<h2>register_diag_field</h2>
<pre> 
<b>register_diag_field</b> (module_name, field_name, axes, init_time, & long_name, units, missing_value, range)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Return field index for subsequent calls to <a href=
"#send_data">send_data</a></dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>module_name   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>field_name   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>axes   </tt></td>
<td><br>
   <span class="type">[integer,
dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>init_time   </tt></td>
<td><br>
   <span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>long_name   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>units   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>missing_value   </tt></td>
<td><br>
   <span class="type">[real]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>range   </tt></td>
<td><br>
   <span class="type">[real,
dimension(2)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="register_static_field" id="register_static_field"></a>
<h2>register_static_field</h2>
<pre> 
<b>register_static_field</b> (module_name, field_name, axes, & long_name, units, missing_value, range, require)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Return field index for subsequent call to send_data.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>module_name   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>field_name   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>axes   </tt></td>
<td><br>
   <span class="type">[integer,
dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>long_name   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>units   </tt></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>missing_value   </tt></td>
<td><br>
   <span class="type">[real]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>range   </tt></td>
<td><br>
   <span class="type">[real,
dimension(2)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="diag_manager_end" id="diag_manager_end"></a>
<h2>diag_manager_end</h2>
<pre>
<b>call diag_manager_end </b>(time)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Flushes diagnostic buffers where necessary. Close diagnostics
files.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>TIME   </tt></td>
<td><br>
   <span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="diag_manager_init" id="diag_manager_init"></a>
<h2>diag_manager_init</h2>
<pre>
<b>call diag_manager_init </b>()</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Open and read diag_table. Select fields and files for
diagnostic output.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="get_base_time" id="get_base_time"></a>
<h2>get_base_time</h2>
<pre>
<b>call get_base_time </b>()</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Return base time for diagnostics (note: base time must be &gt;=
model time).</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="get_base_date" id="get_base_date"></a>
<h2>get_base_date</h2>
<pre>
<b>call get_base_date </b>(year, month, day, hour, minute, second)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Return date information for diagnostic reference time.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="need_data" id="need_data"></a>
<h2>need_data</h2>
<pre> 
<b>need_data</b> (diag_field_id,next_model_time)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Determine whether data is needed for the current model time
step. Since diagnostic data are buffered, the "next" model time is
passed instead of the current model time. This call can be used to
minimize overhead for complicated diagnostics.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>inext_model_time   </tt></td>
<td>next_model_time = current model time + model time_step<br>
   <span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>diag_field_id   </tt></td>
<td><br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="NAMELIST" id="NAMELIST"></a> <!-- BEGIN NAMELIST -->
<!-- END NAMELIST --><a name="DIAGNOSTIC FIELDS"></a> 
<!-- BEGIN DIAGNOSTIC FIELDS -->
<!-- END DIAGNOSTIC FIELDS --><a name="DATA SETS"></a> 
<!-- BEGIN DATA SETS -->
<hr>
<h2>DATA SETS</h2>
<div>None.<br>
<br></div>
<!-- END DATA SETS -->
<a name="ERROR MESSAGES"></a> <!-- BEGIN ERROR MESSAGES -->
<hr>
<h2>ERROR MESSAGES</h2>
<div>None.<br>
<br></div>
<!-- END ERROR MESSAGES -->
<a name="REFERENCES" id="REFERENCES"></a>
<hr>
<h2>REFERENCES</h2>
<!-- BEGIN REFERENCES -->
<div>None.</div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h2>COMPILER SPECIFICS</h2>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
<dl>
<dt>COMPILING AND LINKING SOURCE</dt>
<dd>Any module or program unit using
&lt;TT&gt;diag_manager_mod&lt;/TT&gt; must contain the line
<pre>   use diag_manager_mod</pre>
If netCDF output is desired, the cpp flag
&lt;TT&gt;-Duse_netCDF&lt;/TT&gt; must be turned on. The loader
step requires an explicit link to the netCDF library (typically
something like &lt;TT&gt;-L/usr/local/lib -lnetcdf&lt;/TT&gt;,
depending on the path to the netCDF library). <a href=
"http://www.unidata.ucar.edu/packages/netcdf/guidef">netCDF release
3 for fortran</a> is required.</dd>
</dl>
</div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h2>PRECOMPILER OPTIONS</h2>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
<dl>
<dt>PORTABILITY</dt>
<dd>&lt;TT&gt;diag_manager_mod&lt;/TT&gt; uses standard f90.</dd>
</dl>
</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h2>LOADER OPTIONS</h2>
<!-- BEGIN LOADER -->
<div>
<p>GFDL users can checkout diag_manager_mod using the cvs command
&lt;TT&gt;setenv CVSROOT '/home/fms/cvs';cvs co
diag_manager&lt;/TT&gt;.</p>
<pre>        ACQUIRING SOURCE</pre></div>
<!-- END LOADER OPTIONS -->
<a name="TEST PROGRAM"></a>
<hr>
<h2>TEST PROGRAM</h2>
<!-- BEGIN TEST PROGRAM -->
<div>None.<br></div>
<br>
<!-- END TEST PROGRAM -->
<a name="KNOWN BUGS"></a>
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN KNOWN BUGS -->
<div>None.</div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES" id="NOTES"></a>
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<div>None.<br></div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN FUTURE PLANS -->
<div>None.</div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right"><a href="#TOP"><font size=
"-1">top</font></a></div>
</body>
</html>
