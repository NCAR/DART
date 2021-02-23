<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module schedule_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE schedule_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Provides a set of routines to generate a regular pattern of time
windows. This module is only used for converting observation
sequences files to netCDF format. If it stands the test of time, it
will likely be used to create an assimilation schedule independent
of the observation sequence file. Wouldn't that be nice ...</p>
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
&amp;schedule_nml
   first_bin_start      =  1601,  1,  1,  0,  0,  0
   first_bin_end        =  2999,  1,  1,  0,  0,  0
   last_bin_end         =  2999,  1,  1,  0,  0,  0
   bin_interval_days    = 1000000
   bin_interval_seconds = 0
   max_num_bins         = 1000
   calendar             = 'Gregorian'
   print_table          = .true.
  /
</pre></div>
<br>
<br>
<p>Controls various aspects of filter. The inflation control
variables are all dimensioned 2, the first value being for the
prior inflation and the second being for the posterior
inflation.</p>
<p>The default values will cause (pretty much) all possible
observations to be put into one output file.</p>
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
<td>first_bin_start</td>
<td>integer, dimension(6)</td>
<td>Date/time specification for starting time of first bin.</td>
</tr>
<tr>
<td>first_bin_end</td>
<td>integer, dimension(6)</td>
<td>Date/time specification for ending time of first bin. Sets the
bin width.</td>
</tr>
<tr>
<td>last_bin_end</td>
<td>integer, dimension(6)</td>
<td>Date/time specification for ending time of last bin. Sets the
length of the overall time of the schedule.</td>
</tr>
<tr>
<td>bin_interval_days</td>
<td>integer</td>
<td>Sets the time between bins. Must be larger or equal to the bin
width.</td>
</tr>
<tr>
<td>bin_interval_seconds</td>
<td>integer</td>
<td>Sets the time between bins. Must be larger or equal to the bin
width.</td>
</tr>
<tr>
<td>max_num_bins</td>
<td>integer</td>
<td>Upper limit on the number of bins.</td>
</tr>
<tr>
<td>calendar</td>
<td>character(len=32)</td>
<td>String calendar type. Valid types are listed in the <a href=
"time_manager_mod.html#cal_type">time_manager_mod</a> file.</td>
</tr>
<tr>
<td>print_table</td>
<td>logical</td>
<td>If .TRUE., print out information about the schedule each time
set_regular_schedule() is called.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
time_manager_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table>
<tr>
<td><em class="call">use schedule_mod, only :</em></td>
<td><a href="#schedule_type">schedule_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_regular_schedule">set_regular_schedule</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_time_from_schedule">get_time_from_schedule</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_schedule_length">get_schedule_length</a></td>
</tr>
</table>
<p>Namelist <a href="#Namelist"><em class=
"code">&amp;schedule_mod_nml</em></a> may be read from file
<em class="file">input.nml</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="set_regular_schedule" id="set_regular_schedule"></a><br>
<div class="routine"><em class="call">call
set_regular_schedule(schedule)</em>
<pre>
type(schedule_type), intent(out) :: <em class="code">schedule</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Uses the namelist information to compute and fill a
schedule_type variable.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">schedule</em></td>
<td>Fills this derived type with the information needed to generate
a series of regularly spaced time windows.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_time_from_schedule" id=
"get_time_from_schedule"></a><br>
<div class="routine"><em class="call">call
get_time_from_schedule(mytime, schedule, iepoch <em class=
"optionalcode">[, edge]</em>)</em>
<pre>
type(time_type),     intent(out) :: <em class="code">mytime</em>
 or
real(digits12),      intent(out) :: <em class="code">mytime</em>
type(schedule_type), intent(in)  :: <em class="code">schedule</em>
integer,             intent(in)  :: <em class="code">iepoch</em>
integer, optional,   intent(in)  :: <em class=
"optionalcode">edge</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns either the leading or trailing time for the specified
bin/epoch number for the given schedule. The time can be returned
in one of two formats, depending on the variable type specified for
the first argument: either a DART derived time_type, or a real of
kind digits12 (defined in the types_mod).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">mytime</em></td>
<td>Return value with the leading or trailing edge time for the
requested bin. There are two supported return formats, either as a
standard DART time_type, or as a real value which will contain the
number of days plus any fraction.</td>
</tr>
<tr>
<td valign="top"><em class="code">schedule</em></td>
<td>Schedule type to extract information from.</td>
</tr>
<tr>
<td valign="top"><em class="code">iepoch</em></td>
<td>The bin number, or epoch number, to return a time for. Unless
edge is specified and requests the ending time, the time returned
is the starting time for this bin.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">edge</em></td>
<td>If specified, and if edge is larger than 1, the trailing edge
time of the bin is returned. Any other value, or if this argument
is not specified, returns the leading edge time of the bin.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_schedule_length" id="get_schedule_length"></a><br>
<div class="routine"><em class="call">var =
get_schedule_length()</em>
<pre>
integer                             :: <em class=
"code">get_schedule_length</em>
type(schedule_type), intent(in)     :: <em class=
"code">schedule</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Return the total number of intervals/bins/epochs defined by this
schedule.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">schedule</em></td>
<td>Return number of time intervals in this schedule.</td>
</tr>
</table>
</div>
<br>
<!--=================== DESCRIPTION OF A LOCAL TYPE ===================-->
 <a name="schedule_type" id="schedule_type"></a><br>
<div class="type">
<pre>
<em class="call">type schedule_type</em>
   private
   integer :: num_bins
   integer :: current_bin
   logical :: last_bin
   integer :: calendar
   character(len=32) :: calendarstring
   type(time_type)          :: binwidth
   type(time_type)          :: bininterval
   type(time_type), pointer :: binstart(   :) =&gt; NULL()
   type(time_type), pointer :: binend(     :) =&gt; NULL()
   real(digits12),  pointer :: epoch_start(:) =&gt; NULL()
   real(digits12),  pointer :: epoch_end(  :) =&gt; NULL()
end type schedule_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>This type is used to define a schedule.</p>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<th align="left">filename</th>
<th align="left">purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read the schedule_mod namelist</td>
</tr>
</table>
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
<td valign="top"></td>
<!-- message -->
<td valign="top"></td>
<!-- comment -->
<td valign="top"></td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>Setting the schedule from the namelist values means you can only
have a single schedule object in the entire program. We also need a
subroutine to initialize a schedule type by giving explicit
arguments.</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
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
