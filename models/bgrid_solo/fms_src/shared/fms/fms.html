<!DOCTYPE html>
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>Module fms_mod</title>
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
<h1>Module fms_mod</h1>
<a name="HEADER" id="HEADER"></a> <a name="OVERVIEW" id=
"OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text">The fms module provides routines that are commonly
used by most FMS modules.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<div>Here is a summary of the functions performed by routines in
the fms module.<br>
<br>
1. Output module version numbers to a common (<tt>log</tt>) file
using a common format.<br>
2. Open specific types of files common to many FMS modules. These
include namelist files, restart files, and 32-bit IEEE data files.
There also is a matching interface to close the files. If other
file types are needed the <tt>mpp_open</tt> and <tt>mpp_close</tt>
interfaces in module <a href=
"http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/models/bgrid_solo/fms_src/shared/mpp/mpp_io.html">
</a>must be used.<br>
3. Read and write distributed data to simple native unformatted
files. This type of file (called a restart file) is used to
checkpoint model integrations for a subsequent restart of the
run.<br>
4. Time sections of code (using a wrapper for <a href=
"http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/models/bgrid_solo/fms_src/shared/mpp/mpp.html">
mpp_mod</a>).<br>
5. For convenience there are several routines published from the
<a href=
"http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/models/bgrid_solo/fms_src/shared/mpp/mpp.html">
</a> module. These are routines for getting processor numbers,
commonly used I/O unit numbers, and error handling.</div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>
mpp_mod<br>mpp_domains_mod<br>     mpp_io_mod<br>     fms_io_mod</pre></div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h2>PUBLIC INTERFACE</h2>
<div>
<dl>
<dt><a href="#fms_init">fms_init</a>:</dt>
<dd>Initializes the FMS module and also calls the initialization
routines for all modules in the MPP package. Will be called
automatically if the user does not call it.</dd>
<dt><a href="#fms_end">fms_end</a>:</dt>
<dd>Calls the termination routines for all modules in the MPP
package.</dd>
<dt><a href="#file_exist">file_exist</a>:</dt>
<dd>Checks the existence of a given file name.</dd>
<dt><a href="#error_mesg">error_mesg</a>:</dt>
<dd>Print notes, warnings and error messages; terminates program
for warning and error messages. (use error levels
NOTE,WARNING,FATAL, see example below)</dd>
<dt><a href="#check_nml_error">check_nml_error</a>:</dt>
<dd>Checks the iostat argument that is returned after reading a
namelist and determines if the error code is valid.</dd>
<dt><a href="#write_version_number">write_version_number</a>:</dt>
<dd>Prints to the log file (or a specified unit) the (cvs) version
id string and (cvs) tag name.</dd>
<dt><a href="#mpp_clock_init">mpp_clock_init</a>:</dt>
<dd>Returns an identifier for performance timing a section of code
(similar to mpp_clock_id).</dd>
<dt><a href="#lowercase">lowercase</a>:</dt>
<dd>Convert character strings to all lower case.</dd>
<dt><a href="#uppercase">uppercase</a>:</dt>
<dd>Convert character strings to all upper case.</dd>
<dt><a href="#string_array_index">string_array_index</a>:</dt>
<dd>match the input character string to a string in an array/list
of character strings</dd>
<dt><a href="#monotonic_array">monotonic_array</a>:</dt>
<dd>Determines if a real input array has monotonically increasing
or decreasing values.</dd>
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
<li><a name="fms_init" id="fms_init"></a>
<h2>fms_init</h2>
<pre>
<b>call fms_init </b>( )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Initialization routine for the fms module. It also calls
initialization routines for the mpp, mpp_domains, and mpp_io
modules. Although this routine will be called automatically by
other fms_mod routines, users should explicitly call fms_init. If
this routine is called more than once it will return silently.
There are no arguments.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="fms_end" id="fms_end"></a>
<h2>fms_end</h2>
<pre>
<b>call fms_end </b>( )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Termination routine for the fms module. It also calls
destructor routines for the mpp, mpp_domains, and mpp_io modules.
If this routine is called more than once it will return silently.
There are no arguments.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="file_exist" id="file_exist"></a>
<h2>file_exist</h2>
<pre> 
<b>file_exist</b> ( file_name )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Checks the existence of the given file name. If the file_name
string has zero length or the first character is blank return a
false result.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>file_name   </tt></td>
<td>A file name (or path name) that is checked for existence.<br>
   <span class="type">[character]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">   </td>
<td>This function returns a logical result. If file_name exists the
result is true, otherwise false is returned. If the length of
character string "file_name" is zero or the first character is
blank, then the returned value will be false. When reading a file,
this function is often used in conjunction with routine
open_file.<br>
   <span class="type">[]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="error_mesg" id="error_mesg"></a>
<h2>error_mesg</h2>
<pre>
<b>call error_mesg </b>( routine, message, level )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Print notes, warnings and error messages; and terminates the
program for error messages. This routine is a wrapper around
mpp_error, and is provided for backward compatibility. This module
also publishes mpp_error, <b>users should try to use the mpp_error
interface</b>.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>routine   </tt></td>
<td>Routine name where the warning or error has occurred.<br>
   <span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>message   </tt></td>
<td>Warning or error message to be printed.<br>
   <span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>level   </tt></td>
<td>Level of severity; set to NOTE, WARNING, or FATAL Termination
always occurs for FATAL, never for NOTE, and is settable for
WARNING (see namelist).<br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>NOTE</b></dt>
<dd>Examples:
<pre>        use fms_mod, only: error_mesg, FATAL, NOTE
        call error_mesg ('fms_mod', 'initialization not called', FATAL)
        call error_mesg ('fms_mod', 'fms_mod message', NOTE)</pre></dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="check_nml_error" id="check_nml_error"></a>
<h2>check_nml_error</h2>
<pre> 
<b>check_nml_error</b> ( iostat, nml_name )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>The FMS allows multiple namelist records to reside in the same
file. Use this interface to check the iostat argument that is
returned after reading a record from the namelist file. If an
invalid iostat value is detected this routine will produce a fatal
error. See the NOTE below.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>iostat   </tt></td>
<td>The iostat value returned when reading a namelist record.<br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>nml_name   </tt></td>
<td>The name of the namelist. This name will be printed if an error
is encountered, otherwise the name is not used.<br>
   <span class="type">[character]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">   </td>
<td>This function returns the input iostat value (integer) if it is
an allowable error code. If the iostat error code is not allowable,
an error message is printed and the program terminated.<br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>NOTE</b></dt>
<dd>Some compilers will return non-zero iostat values when reading
through files with multiple namelist. This routine will try skip
these errors and only terminate for true namelist errors.<br>
<br>
Examples<br>
<br>
The following example checks if a file exists, reads a namelist
input from that file, and checks for errors in that namelist. When
the correct namelist is read and it has no errors the routine
check_nml_error will return zero and the while loop will exit. This
code segment should be used to read namelist files.
<pre>          integer :: unit, ierr, io

          if ( file_exist('input.nml') ) then
              unit = open_namelist_file ( )
              ierr=1
              do while (ierr /= 0)
                read  (unit, nml=moist_processes_nml, iostat=io, end=10)
                ierr = check_nml_error(io,'moist_processes_nml')
              enddo
        10    call close_file (unit)
          endif</pre></dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="write_version_number" id="write_version_number"></a>
<h2>write_version_number</h2>
<pre>
<b>call write_version_number </b>( version [, tag, unit] )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Prints to the log file (stdlog) or a specified unit the (cvs)
version id string and (cvs) tag name.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>version   </tt></td>
<td>string that contains routine name and version number.<br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>tag   </tt></td>
<td>The tag/name string, this is usually the Name string returned
by CVS when checking out the code.<br>
   <span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>unit   </tt></td>
<td>The Fortran unit number of an open formatted file. If this unit
number is not supplied the log file unit number is used
(stdlog).<br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_clock_init" id="mpp_clock_init"></a>
<h2>mpp_clock_init</h2>
<pre>id = <b>mpp_clock_init</b> ( name, level [, flags] )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Returns an identifier for performance timing sections of code.
Should be used in conjunction with mpp_clock_begin and
mpp_clock_end. For more details see the documentation for the MPP
module and look at the example below.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>name   </tt></td>
<td>A unique name string given to the code segment to be timed. The
length should not exceed 32 characters.<br>
   <span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>level   </tt></td>
<td>Level of timing. When level &gt; timing_level, which is set by
namelist <a href="#NAMELIST">&amp;fms_nml</a>, an identifier of
zero is returned. This will turn off performance timing for the
code section.<br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>flags   </tt></td>
<td>Use the flags published via the mpp_mod to control whether
synchronization or extra detail is desired. (flags =
MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED)<br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>id   </tt></td>
<td>The identification index returned by mpp_clocks_id. A zero
value is returned (turning clocks off) when input argument level
&gt; namelist variable timing_level.<br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>NOTE</b></dt>
<dd>1.The MPP_CLOCK_SYNC flag should be used whenever possible.
This flag causes mpp_sync to be called at the begin of a code
segment, resulting in more accurate performance timings. <b>Do not
use the MPP_CLOCK_SYNC flag for code sections that may not be
called on all processors.</b><br>
2.There is some amount of coordination required throughout an
entire program for consistency of the "timing levels". As a
guideline the following levels may be used, with higher levels
added as desired to specific component models.
<pre>                level 
                      example code section
                 1 
                      main program
                 2 
                      components models
                 3 
                      atmosphere dynamics or physics</pre>
Examples:<br>
<br>
The mpp_clock_init interface should be used in conjunction with the
mpp_mod interfaces mpp_clock_begin and mpp_clock_end. For example:
<pre>
use fms_mod, only: mpp_clock_init, mpp_clock_begin, &
                             mpp_clock_end. MPP_CLOCK_SYNC
          integer :: id_mycode
          integer :: timing_level = 5

          id_mycode = mpp_clock_init ('mycode loop', timing_level, &
                                      flags=MPP_CLOCK_SYNC)
          call mpp_clock_begin (id_mycode)
                        :
                        :
           ~~ this code will be timed ~~ 
                        :
                        :
          call mpp_clock_end (id_mycode)</pre></dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="lowercase" id="lowercase"></a>
<h2>lowercase</h2>
<pre>string = <b>lowercase</b> ( cs )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Converts a character string to all lower case letters. The
characters "A-Z" are converted to "a-z", all other characters are
left unchanged.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>cs   </tt></td>
<td>Character string that may contain upper case letters.<br>
   <span class="type">[character(len=*),
scalar]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>string   </tt></td>
<td>Character string that contains all lower case letters. The
length of this string must be the same as the input string.<br>
   <span class="type">[character(len=len(cs)),
scalar]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="uppercase" id="uppercase"></a>
<h2>uppercase</h2>
<pre>string = <b>uppercase</b> ( cs )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Converts a character string to all upper case letters. The
characters "a-z" are converted to "A-Z", all other characters are
left unchanged.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>cs   </tt></td>
<td>Character string that may contain lower case letters.<br>
   <span class="type">[character(len=*),
scalar]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>string   </tt></td>
<td>Character string that contains all upper case letters. The
length of this string must be the same as the input string.<br>
   <span class="type">[character(len=len(cs)),
scalar]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="string_array_index" id="string_array_index"></a>
<h2>string_array_index</h2>
<pre> 
<b>string_array_index</b> ( string, string_array [, index] )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Tries to find a match for a character string in a list of
character strings. The match is case sensitive and disregards blank
characters to the right of the string.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>string   </tt></td>
<td>Character string of arbitrary length.<br>
   <span class="type">[character(len=*),
scalar]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>string_array   </tt></td>
<td>Array/list of character strings.<br>
   <span class="type">[character(len=*),
dimension(:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>index   </tt></td>
<td>The index of string_array where the first match was found. If
no match was found then index = 0.<br>
   <span class="type">[]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>found   </tt></td>
<td>If an exact match was found then TRUE is returned, otherwise
FALSE is returned.<br>
   <span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>NOTE</b></dt>
<dd>Examples
<pre>       string = "def"
       string_array = (/ "abcd", "def ", "fghi" /)
       string_array_index ( string, string_array, index )
       Returns: TRUE, index = 2</pre></dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="monotonic_array" id="monotonic_array"></a>
<h2>monotonic_array</h2>
<pre> 
<b>monotonic_array</b> ( array [, direction] )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Determines if the real input array has monotonically increasing
or decreasing values.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>array   </tt></td>
<td>An array of real values. If the size(array) &lt; 2 this
function assumes the array is not monotonic, no fatal error will
occur.<br>
   <span class="type">[real,
dimension(:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>direction   </tt></td>
<td>If the input array is: &gt;&gt; monotonic (small to large) then
direction = +1. &gt;&gt; monotonic (large to small) then direction
= -1. &gt;&gt; not monotonic then direction = 0.<br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left">   </td>
<td>If the input array of real values either increases or decreases
monotonically then TRUE is returned, otherwise FALSE is
returned.<br>
   <span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="PUBLIC TYPES"></a> <!-- BEGIN PUBLIC TYPES -->
<!-- END PUBLIC TYPES --><a name="NAMELIST" id="NAMELIST"></a> 
<!-- BEGIN NAMELIST -->
<hr>
<h2>NAMELIST</h2>
<div><b>&amp;fms_nml</b><br>
<br>
<div>
<dl>
<dt><tt>timing_level</tt></dt>
<dd>
<dl>
<dt>The level of performance timing. If calls to the performance
timing routines have been inserted into the code then code sections
with a level &lt;= timing_level will be timed. The resulting output
will be printed to STDOUT. See the MPP module or <a href=
"#mpp_clock_init">mpp_clock_init</a> for more details.<br>
<span class="type">[integer, default: 0]</span></dt>
</dl>
</dd>
<dt><tt>read_all_pe</tt></dt>
<dd>
<dl>
<dt>Read global data on all processors extracting local part needed
(TRUE) or read global data on PE0 and broadcast to all PEs
(FALSE).<br>
<span class="type">[logical, default: true]</span></dt>
</dl>
</dd>
<dt><tt>warning_level</tt></dt>
<dd>
<dl>
<dt>Sets the termination condition for the WARNING flag to
interfaces error_mesg/mpp_error. set warning_level = 'fatal'
(program crashes for warning messages) or 'warning' (prints warning
message and continues).<br>
<span class="type">[character, default: 'warning']</span></dt>
</dl>
</dd>
<dt><tt>iospec_ieee32</tt></dt>
<dd>
<dl>
<dt>iospec flag used with the open_ieee32_file interface.<br>
<span class="type">[character, default: '-F
f77,cachea:48:1']</span></dt>
</dl>
</dd>
<dt><tt>stack_size</tt></dt>
<dd>
<dl>
<dt>The size in words of the MPP user stack. If stack_size &gt; 0,
the following MPP routine is called: call mpp_set_stack_size
(stack_size). If stack_size = 0 (default) then the default size set
by mpp_mod is used.<br>
<span class="type">[integer, default: 0]</span></dt>
</dl>
</dd>
<dt><tt>domains_stack_size</tt></dt>
<dd>
<dl>
<dt>The size in words of the MPP_DOMAINS user stack. If
domains_stack_size &gt; 0, the following MPP_DOMAINS routine is
called: call mpp_domains_set_stack_size (domains_stack_size). If
domains_stack_size = 0 (default) then the default size set by
mpp_domains_mod is used.<br>
<span class="type">[integer, default: 0]</span></dt>
</dl>
</dd>
</dl>
</div>
</div>
<br>
<!-- END NAMELIST -->
<a name="DIAGNOSTIC FIELDS"></a> <!-- BEGIN DIAGNOSTIC FIELDS -->
<!-- END DIAGNOSTIC FIELDS --><a name="DATA SETS"></a> 
<!-- BEGIN DATA SETS -->
<hr>
<h2>DATA SETS</h2>
<div>None.<br>
<br></div>
<!-- END DATA SETS -->
<a name="PUBLIC CODE"></a> <!-- BEGIN PUBLIC CODE -->
<!-- END PUBLIC CODE --><a name="ERROR MESSAGES"></a> 
<!-- BEGIN ERROR MESSAGES -->
<hr>
<h2>ERROR MESSAGES</h2>
<div>
<dl>
<dt><b>FATAL in fms_init</b></dt>
<dd><span class="errmsg">invalid entry for namelist variable
warning_level</span></dd>
<dd>The namelist variable warning_level must be either 'fatal' or
'warning' (case-insensitive).</dd>
<dt><b>FATAL in file_exist</b></dt>
<dd><span class="errmsg">set_domain not called</span></dd>
<dd>Before calling write_data you must first call set_domain with
domain2d data type associated with the distributed data you are
writing.</dd>
<dt><b>FATAL in check_nml_error</b></dt>
<dd><span class="errmsg">while reading namelist ...., iostat =
####</span></dd>
<dd>There was an error message reading the namelist specified.
Carefully examine all namelist variables for misspellings of type
mismatches (e.g., integer vs. real).</dd>
</dl>
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
<div>None.</div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h2>PRECOMPILER OPTIONS</h2>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>None.</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h2>LOADER OPTIONS</h2>
<!-- BEGIN LOADER -->
<div>None.<br>
<br></div>
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
<div>
<p>Namelist error checking may not work correctly with some
compilers.<br>
<br>
Users should beware when mixing Fortran reads and read_data calls.
If a Fortran read follows read_data and namelist variable
read_all_pe = FALSE (not the default), then the code will fail. It
is safest if Fortran reads precede calls to read_data.</p>
</div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES" id="NOTES"></a>
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<div>1) If the <b>MPP</b> or <b>MPP_DOMAINS</b> stack size is
exceeded the program will terminate after printing the required
size.<br>
<br>
2) When running on a very small number of processors or for high
resolution models the default domains_stack_size will probably be
insufficient.</div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN FUTURE PLANS -->
<div>
<p>NetCDF facilities for reading and writing restart files and
(IEEE32) data files.</p>
<p>May possible split the FMS module into two modules.<br>
<br>
i.general utilities (FMS_MOD)<br>
ii.I/O utilities (FMS_IO_MOD)</p>
</div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right"><a href="#TOP"><font size=
"-1">top</font></a></div>
</body>
</html>
