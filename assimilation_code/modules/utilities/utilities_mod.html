<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module utilities_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE utilities_mod</h1>
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
<p>Provides a number of tools used by most DART modules including
tools for file IO, diagnostic tools for registering modules and
recording namelist arguments, and an error handler.</p>
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
&amp;utilities_nml
   TERMLEVEL      = 2,
   logfilename    = 'dart_log.out',
   nmlfilename    = 'dart_log.nml',
   module_details = .true.,
   print_debug    = .false.,
   write_nml      = 'file'
/
</pre></div>
<br>
<br>
<p>The namelist controls how the logging, namelist, messages, and
general utility routines behave.</p>
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
<td>TERMLEVEL</td>
<td>integer</td>
<td>Level at which calls to error manager terminate program. The
default setting is warnings and errors terminate the program.
Setting this to 2 (E_ERR) means only errors terminate. Setting this
to 3 means even errors do not cause an exit (which is not a good
idea).</td>
</tr>
<tr>
<td>logfilename</td>
<td>character(len=129)</td>
<td>File to which the log messages are written.</td>
</tr>
<tr>
<td>nmlfilename</td>
<td>character(len=129)</td>
<td>File to which the namelist output is written. Can be the same
name as the logfile.</td>
</tr>
<tr>
<td>module_details</td>
<td>logical</td>
<td>Each source code module can write out the repository version
number and filename to the logfile. Verbose, but useful for knowing
what version of the code was used during the run.</td>
</tr>
<tr>
<td>print_debug</td>
<td>logical</td>
<td>Setting this to .true. causes additional debug messages to
print. These can be very verbose and by default are turned
off.</td>
</tr>
<tr>
<td>write_nml</td>
<td>character(len=32)</td>
<td>String which controls where to write the namelist values that
are being used for this execution. Valid values are: 'none',
'file', 'terminal', 'both'. 'none' turns off this write. 'file'
writes a copy only to the <em class="code">nmlfilename</em>. Writes
are always in append mode, so the most recent information will be
at the end of an existing file. 'terminal' will write to the job's
standard output. 'both' will write both to the nml file and the
standard output unit.</td>
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
netCDF
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
<td><em class="call">use utilities, only :</em></td>
<td><a href="#file_exist">file_exist</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_unit">get_unit</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#open_file">open_file</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#close_file">close_file</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#timestamp">timestamp</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#register_module">register_module</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#error_handler">error_handler</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#to_upper">to_upper</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_check">nc_check</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#logfileunit">logfileunit</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nmlfileunit">nmlfileunit</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#initialize_utilities">initialize_utilities</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#finalize_utilities">finalize_utilities</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#dump_unit_attributes">dump_unit_attributes</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#find_namelist_in_file">find_namelist_in_file</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#check_namelist_read">check_namelist_read</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#find_textfile_dims">find_textfile_dims</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#file_to_text">file_to_text</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#is_longitude_between">is_longitude_between</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_next_filename">get_next_filename</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_filename_list">set_filename_list</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_tasknum">set_tasknum</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_output">set_output</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#do_output">do_output</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ERROR_LEVELS">E_DBG, DEBUG</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ERROR_LEVELS">E_MSG, MESSAGE</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ERROR_LEVELS">E_WARN, WARNING</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ERROR_LEVELS">E_ERR, FATAL</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="file_exist" id="file_exist"></a><br>
<div class="routine"><em class="call">var =
file_exist(file_name)</em>
<pre>
logical                      :: <em class="code">file_exist</em>
character(len=*), intent(in) :: <em class="code">file_name</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if file_name exists in the working directory, else
false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>True if file_name exists in working directory.</td>
</tr>
<tr>
<td valign="top"><em class="code">file_name  </em></td>
<td>Name of file to look for.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_unit" id="get_unit"></a><br>
<div class="routine"><em class="call">var = get_unit()</em>
<pre>
integer :: <em class="code">get_unit</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns an unused unit number for IO.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>An unused unit number.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="open_file" id="open_file"></a><br>
<div class="routine"><em class="call">var = open_file(fname
<em class="optionalcode">[, form, action]</em>)</em>
<pre>
integer                                :: <em class=
"code">open_file</em>
character(len=*), intent(in)           :: <em class=
"code">fname</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">form</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">action</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a unit number that is opened to the file fname. If form
is not present or if form is "formatted" or "FORMATTED", file is
opened for formatted IO. Otherwise, it is unformatted. The action
string is the standard action string for Fortran IO (see F90
language description).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Unit number opened to file fname.</td>
</tr>
<tr>
<td valign="top"><em class="code">fname  </em></td>
<td>Name of file to be opened.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">form  </em></td>
<td>Format: 'formatted' or 'FORMATTED' give formatted, anything
else is unformatted. Default is formatted.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">action  </em></td>
<td>Standard fortran string description of requested file open
action.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="timestamp" id="timestamp"></a><br>
<div class="routine"><em class="call">call timestamp(<em class=
"optionalcode">[string1, string2, string3,]  </em>pos)</em>
<pre>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">string1</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">string2</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">string3</em>
character(len=*), intent(in)           :: <em class="code">pos</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Prints the message 'Time is YYYY MM DD HH MM SS' to the logfile
along with three optional message strings. If the pos argument is
'end', the message printed is 'Finished... at YYYY MM DD HH MM SS'
and the logfile is closed.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"optionalcode">string1  </em></td>
<td>An optional message to be printed.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">string2  </em></td>
<td>An optional message to be printed.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">string3  </em></td>
<td>An optional message to be printed.</td>
</tr>
<tr>
<td valign="top"><em class="code">pos  </em></td>
<td>If 'end' terminates log_file output.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="close_file" id="close_file"></a><br>
<div class="routine"><em class="call">call close_file(iunit)</em>
<pre>
integer, intent(in) :: <em class="code">iunit</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Closes the given unit number. If the unit is not open, nothing
happens.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">iunit  </em></td>
<td>File unit to be closed.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="register_module" id="register_module"></a><br>
<div class="routine"><em class="call">call register_module(src,
rev, rdate)</em>
<pre>
character(len=*), intent(in) :: <em class="code">src</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">rev</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">rdate</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes the source name to both the logfileunit and to standard
out. The rev and revdate are deprecated as they are unsupported by
git.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">src  </em></td>
<td>source file name.</td>
</tr>
<tr>
<td valign="top"><em class="code">rev  </em></td>
<td>ignored</td>
</tr>
<tr>
<td valign="top"><em class="code">rdate  </em></td>
<td>ignored</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="error_handler" id="error_handler"></a><br>
<div class="routine"><em class="call">call error_handler(level,
routine, text, src, rev, rdate <em class="optionalcode">[, aut,
text2, text3]</em>)</em>
<pre>
integer, intent(in)                    :: <em class=
"code">level</em>
character(len=*), intent(in)           :: <em class=
"code">routine</em>
character(len=*), intent(in)           :: <em class=
"code">text</em>
character(len=*), intent(in)           :: <em class="code">src</em>
character(len=*), intent(in)           :: <em class="code">rev</em>
character(len=*), intent(in)           :: <em class=
"code">rdate</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">aut</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">text2</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">text3</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Prints an error message to standard out and to the logfileunit.
The message contains the routine name, an error message, the source
file, revision and revision date, and optionally the author. The
level of severity is message, debug, warning, or error. If the
level is greater than or equal to the TERMLEVEL (set in the
namelist), execution is terminated. The default TERMLEVEL only
stops for ERRORS.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">level  </em></td>
<td>Error severity (message, debug, warning, error). See below for
specific ations.</td>
</tr>
<tr>
<td valign="top"><em class="code">routine  </em></td>
<td>Name of routine generating error.</td>
</tr>
<tr>
<td valign="top"><em class="code">text  </em></td>
<td>Error message.</td>
</tr>
<tr>
<td valign="top"><em class="code">src  </em></td>
<td>Source file containing routine generating message.</td>
</tr>
<tr>
<td valign="top"><em class="code">rev  </em></td>
<td>Revision number of source file.</td>
</tr>
<tr>
<td valign="top"><em class="code">rdate  </em></td>
<td>Revision date of source file.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">aut  </em></td>
<td>Author of routine.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">text2  </em></td>
<td>If specified, the second line of text for the error
message.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">text3  </em></td>
<td>If specified, the third line of text for the error
message.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="find_namelist_in_file" id=
"find_namelist_in_file"></a><br>
<div class="routine"><em class="call">call
find_namelist_in_file(namelist_file_name, nml_name, iunit,
<em class="optionalcode">[,write_to_logfile_in]</em>)</em>
<pre>
character(len=*),  intent(in)          :: <em class=
"code">namelist_file_name</em>
character(len=*),  intent(in)          :: <em class=
"code">nml_name</em>
integer,           intent(out)         :: <em class=
"code">iunit</em>
logical, optional, intent(in)          :: <em class=
"optionalcode">write_to_logfile_in</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Opens the file namelist_file_name if it exists on unit iunit. A
fatal error occurs if the file does not exist (DART requires an
input.nml to be available, even if it contains no values). Searches
through the file for a line containing ONLY the string
&amp;nml_name (for instance &amp;filter_nml if nml_name is
"filter_nml"). If this line is found, the file is rewound and the
routine returns. Otherwise, a fatal error message is issued.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">namelist  </em></td>
<td>Name of file assumed to hold the namelist.</td>
</tr>
<tr>
<td valign="top"><em class="code">nml_name  </em></td>
<td>Name of the namelist to be searched for in the file, for
instance, filter_nml.</td>
</tr>
<tr>
<td valign="top"><em class="code">iunit  </em></td>
<td>Channel number on which file is opened.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">write_to_logfile_in  </em></td>
<td>When the namelist for the utilities module is read, the logfile
has not yet been open because its name is in the namelist. If
errors are found, have to write to standard out. So, when utilities
module calls this internally, this optional argument is set to
false. For all other applications, it is normally not used (default
is false).</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="check_namelist_read" id="check_namelist_read"></a><br>
<div class="routine"><em class="call">call
check_namelist_read(iunit, iostat_in, nml_name, <em class=
"optionalcode">[, write_to_logfile_in]</em>)</em>
<pre>
integer, intent(in)                    :: <em class=
"code">iunit</em>
integer, intent(in)                    :: <em class=
"code">iostat_in</em>
character(len=*), intent(in)           :: <em class=
"code">nml_name</em>
logical, optional, intent(in)          :: <em class=
"optionalcode">write_to_logfile_in</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Once a namelist has been read from an opened namelist file, this
routine checks for possible errors in the read. If the namelist
read was successful, the file opened on iunit is closed and the
routine returns. If iostat is not zero, an attempt is made to
rewind the file on iunit and read the last line that was
successfully read. If this can be done, this last line is printed
with the preamble "INVALID NAMELIST ENTRY". If the attempt to read
the line after rewinding fails, it is assumed that the original
read (before the call to this subroutine) failed by reaching the
end of the file. An error message stating that the namelist started
but was never terminated is issued.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">iunit  </em></td>
<td>Channel number on which file is opened.</td>
</tr>
<tr>
<td valign="top"><em class="code">iostat_in  </em></td>
<td>Error status return from an attempted read of a namelist from
this file.</td>
</tr>
<tr>
<td valign="top"><em class="code">nml_name  </em></td>
<td>The name of the namelist that is being read (for instance
filter_nml).</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">write_to_logfile_in  </em></td>
<td>When the namelist for the utilities module is read, the logfile
has not yet been open because its name is in the namelist. If
errors are found, have to write to standard out. So, when utilities
module calls this internally, this optional argument is set to
false. For all other applications, it is normally not used (default
is false).</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="find_textfile_dims" id="find_textfile_dims"></a><br>
<div class="routine"><em class="call">call find_textfile_dims
(fname, nlines, linelen)</em>
<pre>
character(len=*), intent (IN)  :: <em class="code">fname</em>
integer,          intent (OUT) :: <em class="code">nlines</em>
integer,          intent (OUT) :: <em class="code">linelen</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Determines the number of lines and maximum line length of an
ASCII text file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">fname</em></td>
<td>input, character string file name</td>
</tr>
<tr>
<td valign="top"><em class="code">nlines</em></td>
<td>output, number of lines in the file</td>
</tr>
<tr>
<td valign="top"><em class="code">linelen</em></td>
<td>output, length of longest line in the file</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="file_to_text" id="file_to_text"></a><br>
<div class="routine"><em class="call">call file_to_text (fname,
textblock)</em>
<pre>
character(len=*),               intent (IN)  :: <em class=
"code">fname</em>
character(len=*), dimension(:), intent (OUT) :: <em class=
"code">textblock</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Opens the given filename and reads ASCII text lines into a
character array.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">fname</em></td>
<td>input, character string file name</td>
</tr>
<tr>
<td valign="top"><em class="code">textblock</em></td>
<td>output, character array of text in the file</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="is_longitude_between" id="is_longitude_between"></a><br>
<div class="routine"><em class="call">var =
is_longitude_between(lon, minlon, maxlon <em class=
"optionalcode">[, doradians]</em>)</em>
<pre>
real(r8), intent(in)           :: <em class="code">lon</em>
real(r8), intent(in)           :: <em class="code">minlon</em>
real(r8), intent(in)           :: <em class="code">maxlon</em>
logical,  intent(in), optional :: <em class=
"optionalcode">doradians</em>
logical                        :: <em class=
"code">is_longitude_between</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Uniform way to test longitude ranges, in degrees, on a globe.
Returns true if lon is between min and max, starting at min and
going EAST until reaching max. Wraps across 0 longitude. If min
equals max, all points are inside. Includes endpoints. If optional
arg doradians is true, do computation in radians between 0 and 2*PI
instead of default 360. There is no rejection of input values based
on range; they are all converted to a known range by calling
modulo() first.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>True if lon is between min and max.</td>
</tr>
<tr>
<td valign="top"><em class="code">lon</em></td>
<td>Location to test.</td>
</tr>
<tr>
<td valign="top"><em class="code">minlon</em></td>
<td>Minimum longitude. Region will start here and go east.</td>
</tr>
<tr>
<td valign="top"><em class="code">maxlon</em></td>
<td>Maximum longitude. Region will end here.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">doradians</em></td>
<td>Optional argument. Default computations are in degrees. If this
argument is specified and is .true., do the computation in radians,
and wrap across the globe at 2 * PI. All inputs must then be
specified in radians.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_next_filename" id="get_next_filename"></a><br>
<div class="routine"><em class="call">var = get_next_filename(
listname, lineindex )</em>
<pre>
character(len=*),  intent(in) :: <em class="code">listname</em>
integer,           intent(in) :: <em class="code">lineindex</em>
character(len=128)            :: <em class=
"code">get_next_filename</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the specified line of a text file, given a filename and
a line number. It returns an empty string when the line number is
larger than the number of lines in a file.</p>
<p>Intended as an easy way to process a list of files. Use a
command like 'ls &gt; out' to create a file containing the list, in
order, of files to be processed. Then call this function with an
increasing index number until the return value is empty.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>An ascii string, up to 128 characters long, containing the
contents of line <em class="code">lineindex</em> of the input
file.</td>
</tr>
<tr>
<td valign="top"><em class="code">listname</em></td>
<td>The filename to open and read lines from.</td>
</tr>
<tr>
<td valign="top"><em class="code">lineindex</em></td>
<td>Integer line number, starting at 1. If larger than the number
of lines in the file, the empty string '' will be returned.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_filename_list" id="set_filename_list"></a><br>
<div class="routine"><em class="call">var = set_filename_list(
name_array, listname, caller_name )</em>
<pre>
character(len=*),  intent(inout) :: <em class=
"code">name_array</em>
character(len=*),  intent(in)    :: <em class="code">listname</em>
character(len=*),  intent(in)    :: <em class=
"code">caller_name</em>
integer                          :: <em class="code">var</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the count of filenames specified. Verifies that one of
either the name_array or the listname was specified but not both.
If the input was a listname copy the names into the name_array so
when this routine returns all the filenames are in name_array().
Verifies that no more than the allowed number of names was
specified if the input was a listname file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>The count of input files specified.</td>
</tr>
<tr>
<td valign="top"><em class="code">name_array</em></td>
<td>Array of input filename strings. Either this item or the
listname must be specified, but not both.</td>
</tr>
<tr>
<td valign="top"><em class="code">listname</em></td>
<td>The filename to open and read filenames from, one per line.
Either this item or the name_array must be specified but not
both.</td>
</tr>
<tr>
<td valign="top"><em class="code">caller_name</em></td>
<td>Calling subroutine name, used for error messages.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="to_upper" id="to_upper"></a><br>
<div class="routine"><em class="call">call to_upper(string)</em>
<pre>
character(len=*), intent (INOUT) :: <em class="code">string</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Converts the character string to UPPERCASE - in place. The input
string <strong>is</strong> modified.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top" width="15%"><em class="code">string</em></td>
<td>any character string</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_check" id="nc_check"></a><br>
<div class="routine"><em class="call">call nc_check(istatus,
subr_name <em class="optionalcode">[, context]</em>)</em>
<pre>
integer, intent(in)                    :: <em class=
"code">istatus</em>
character(len=*), intent(in)           :: <em class=
"code">subr_name</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">context</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Check the return code from a netcdf call. If no error, return
without taking any action. If an error is indicated (in the
<em class="code">istatus</em> argument) then call the error handler
with the subroutine name and any additional context information
(e.g. which file or which variable was being processed at the time
of the error). All errors are currently hardcoded to be <em class=
"code">FATAL</em> and this routine will not return.</p>
This routine calls a netCDF library routine to construct the text
error message corresponding to the error code in the first
argument. An example use of this routine is:
<pre>
call nc_check(nf90_create(path = trim(ncFileID%fname), cmode = nf90_share, ncid = ncFileID%ncid), &amp;
             'init_diag_output', 'create '//trim(ncFileID%fname))

</pre>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">istatus  </em></td>
<td>The return value from any netCDF call.</td>
</tr>
<tr>
<td valign="top"><em class="code">subr_name  </em></td>
<td>String name of the current subroutine, used in case of
error.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">context  </em></td>
<td>Additional text to be used in the error message, for example to
indicate which file or which variable is being processed.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_tasknum" id="set_tasknum"></a><br>
<div class="routine"><em class="call">call
set_tasknum(tasknum)</em>
<pre>
integer, intent(in)               :: <em class="code">tasknum</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Intended to be used in the MPI multi-task case. Sets the local
task number, which is then prepended to subsequent messages.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">tasknum  </em></td>
<td>Task number returned from MPI_Comm_Rank(). MPI task numbers are
0 based, so for a 4-task job these numbers are 0-3.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_output" id="set_output"></a><br>
<div class="routine"><em class="call">call set_output(doflag)</em>
<pre>
logical, intent(in)               :: <em class="code">doflag</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the status of output. Can be set on a per-task basis if you
are running with multiple tasks. If set to false only warnings and
fatal errors will write to the log. The default in the multi-task
case is controlled by the MPI module initialization code, which
sets task 0 to .TRUE. and all other tasks to .FALSE.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">doflag  </em></td>
<td>Sets, on a per-task basis, whether messages are to be written
to the logfile or standard output. Warnings and errors are always
output.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="do_output" id="do_output"></a><br>
<div class="routine"><em class="call">var = do_output()</em>
<pre>
logical                      :: <em class="code">do_output</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if this task should write to the log, false
otherwise. Set by the <em class="code">set_output()</em> routine.
Defaults to true for the single task case. Can be used in code like
so:</p>
<pre>
if (do_output()) then
 write(*,*) 'At this point in the code'
endif
</pre>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>True if this task should write output.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="initialize_utilities" id="initialize_utilities"></a><br>
<div class="routine"><em class="call">call initialize_utilities(
<em class="optionalcode">[progname]</em> <em class=
"optionalcode">[, alternatename]</em> )</em>
<pre>
character(len=*), intent(in), optional :: <em class=
"optionalcode">progname</em>
character(len=*), intent(in), optional :: <em class=
"optionalcode">alternatename</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads the namelist and opens the logfile. Records the values of
the namelist and registers this module.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"optionalcode">progname  </em></td>
<td>If given, use in the timestamp message in the log file to say
which program is being started.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">alternatename  </em></td>
<td>If given, log filename to use instead of the value in the
namelist. This permits, for example, different programs sharing the
same input.nml file to have different logs. If not given here and
no value is specified in the namelist, this defaults to
dart_log.out</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="finalize_utilities" id="finalize_utilities"></a><br>
<div class="routine"><em class="call">call
finalize_utilities()</em></div>
<div class="indent1"><!-- Description -->
<p>Closes the logfile; using utilities after this call is a bad
idea.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="dump_unit_attributes" id="dump_unit_attributes"></a><br>
<div class="routine"><em class="call">call
dump_unit_attributes(iunit)</em>
<pre>
integer, intent(in) :: <em class="code">iunit</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes all information about the status of the IO unit to the
error handler with error level message.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">iunit  </em></td>
<td>Unit about which information is requested.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A CONSTANT =====================-->
 <a name="ERROR_LEVELS" id="ERROR_LEVELS"></a><br>
<div class="routine">
<pre>
integer :: <em class="code">E_DBG, DEBUG</em>
integer :: <em class="code">E_MSG, MESSAGE</em>
integer :: <em class="code">E_WARN, WARNING</em>
integer :: <em class="code">E_ERR, FATAL</em>
</pre></div>
<div class="indent1"><!-- Description -->
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"></td>
<td>Severity levels to be passed to error handler. Levels are
debug, message, warning and fatal. The namelist parameter TERMLEVEL
can be used to control at which level program termination should
occur.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A CONSTANT =====================-->
 <a name="logfileunit" id="logfileunit"></a><br>
<div class="routine">
<pre>
integer :: <em class="code">logfileunit</em>
</pre></div>
<div class="indent1"><!-- Description -->
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top" width="15%"><em class="code">logfileunit</em></td>
<td>Unit opened to file for diagnostic output.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A CONSTANT =====================-->
 <a name="nmlfileunit" id="nmlfileunit"></a><br>
<div class="routine">
<pre>
integer :: <em class="code">nmlfileunit</em>
</pre></div>
<div class="indent1"><!-- Description -->
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top" width="15%"><em class="code">nmlfileunit</em></td>
<td>Unit opened to file for diagnostic output of namelist files.
Defaults to same as <em class="code">logfileunit</em>. Provides the
flexibility to log namelists to a separate file, reducing the
clutter in the log files and perhaps increasing readability.</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>assim_model_mod.nml in input.nml</li>
<li>logfile, name specified in namelist</li>
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
<td valign="top">get_unit</td>
<!-- message -->
<td valign="top">No available units</td>
<!-- comment -->
<td valign="top">Unable to open enough IO channels</td>
</tr>
<tr><!-- routine -->
<td valign="top">check_nml_error</td>
<!-- message -->
<td valign="top">while reading namelist _____</td>
<!-- comment -->
<td valign="top">Fatal error reading namelist. This could be caused
by having an entry in the namelist input file that is not in the
namelist, by having illegal values for namelist variables, or by a
variety of other compiler dependent problems.</td>
</tr>
<tr><!-- routine -->
<td valign="top">find_namelist_in_file</td>
<!-- message -->
<td valign="top">Namelist entry &amp;____ must exist in
namelist_nml.</td>
<!-- comment -->
<td valign="top">There must be an entry for the required namelist,
for instance &amp;filter_nml, in the input.nml namelist file. Even
if no values are to be changed from the default, an entry like
&amp;filter_nml followed by a line containing only / is
required.</td>
</tr>
<tr><!-- routine -->
<td valign="top">find_namelist_in_file</td>
<!-- message -->
<td valign="top">Namelist input file: input.nml must exist</td>
<!-- comment -->
<td valign="top">The namelist input file (usually input.nml) must
exist.</td>
</tr>
<tr><!-- routine -->
<td valign="top">check_namelist_read</td>
<!-- message -->
<td valign="top">INVALID NAMELIST ENTRY: ___ in namelist ____</td>
<!-- comment -->
<td valign="top">While reading the namelist, either a bad entry was
found or an end of file was encountered. The most confusing case is
when a namelist is being read successfully but is not appropriately
terminated with a /. The line printed out by the error message will
be the start of the next namelist in the input.nml file in this
case.</td>
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
<p>none</p>
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
