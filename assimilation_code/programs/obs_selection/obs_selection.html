<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program obs_selection</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>program <em class="program">obs_selection</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Building">BUILDING</a>
/ <a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>This specialized tool selects a subset of input observations
from an observation sequence file. For a more general purpose
observation sequence file tool, see the <a href=
"../obs_sequence_tool/obs_sequence_tool.html">obs_sequence_tool</a>.
This tool takes a selected list of observation types, times, and
locations, and extracts only the matching observations out of one
or more obs_sequence files. The tool which creates the input
selection file is usually <a href=
"../obs_seq_coverage/obs_seq_coverage.html">obs_seq_coverage</a>.
Alternatively, the selection file can be a full observation
sequence file, in which case the types, times, and locations of
those observations are used as the selection criteria.</p>
<p>This tool processes each observation sequence file listed in the
input namelist <em class="code">filename_seq</em> or <em class=
"code">filename_seq_list</em>. If the observation type, time and
location matches an entry in the selection file, it is copied
through to the output. Otherwise it is ignored.</p>
<p>The actions of the <em class="code">obs_selection</em> program
are controlled by a Fortran namelist, read from a file named
<em class="file">input.nml</em> in the current directory. A
detailed description of each namelist item is described in the
<a href="#Namelist">namelist section</a> of this document. The
names used in this discussion refer to these namelist items.</p>
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
&amp;obs_selection_nml
   filename_seq          = ''
   filename_seq_list     = ''
   filename_out          = 'obs_seq.processed'
   num_input_files       = 0
   selections_file       = 'obsdef_mask.txt'
   selections_is_obs_seq = .false.
   latlon_tolerance      = 0.000001
   match_vertical        = .false.
   surface_tolerance     = 0.0001
   pressure_tolerance    = 0.001
   height_tolerance      = 0.0001
   scaleheight_tolerance = 0.001
   level_tolerance       = 0.00001
   print_only            = .false.
   partial_write         = .false.
   print_timestamps      = .false.
   calendar              = 'Gregorian'
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
<td>filename_seq</td>
<td>character(len=256), dimension(500)</td>
<td>The array of names of the observation sequence files to
process, up to a max count of 500 files. (Specify only the actual
number of input files. It is not necessary to specify 500
entries.)</td>
</tr>
<tr>
<td>filename_seq_list</td>
<td>character(len=256)</td>
<td>An alternative way to specify the list of input files. The name
of a text file which contains, one per line, the names of the
observation sequence files to process. You can only specify one of
filename_seq OR filename_seq_list, not both.</td>
</tr>
<tr>
<td>num_input_files</td>
<td>integer</td>
<td>Optional. The number of observation sequence files to process.
Maximum of 500. If 0, the length is set by the number of input
files given. If non-zero, must match the given input file list
length. (Can be used to verify the right number of input files were
processed.)</td>
</tr>
<tr>
<td>filename_out</td>
<td>character(len=256)</td>
<td>The name of the resulting output observation sequence file.
There is only a single output file from this tool. If the input
specifies multiple obs_seq input files, the results are
concatinated into a single output file.</td>
</tr>
<tr>
<td>selections_file</td>
<td>character(len=256)</td>
<td>The name of the input file containing the mask of observation
definitions (the textfile output of <a href=
"../obs_seq_coverage/obs_seq_coverage.html">obs_seq_coverage</a>).
Alternatively, this can be the name of a full observation sequence
file. In this case, the types, times, and locations are extracted
from this file and then used in the same manner as a mask file from
the coverage tool.</td>
</tr>
<tr>
<td>selections_is_obs_seq</td>
<td>logical</td>
<td>If .TRUE. the filename given for the "selections_file" is a
full obs_sequence file and not a text file from the coverage
tool.</td>
</tr>
<tr>
<td>latlon_tolerance</td>
<td>real(r8)</td>
<td>Specified in degrees. For observations to match in the
horizontal the difference in degrees for each of latitude and
longitude must be less than this threshold. If less than or equal
to 0, the values must match exactly.</td>
</tr>
<tr>
<td>match_vertical</td>
<td>logical</td>
<td>If .TRUE. the locations of the observations in the input files
have to match the selection list not only the horizontal but also
in the vertical.</td>
</tr>
<tr>
<td>surface_tolerance</td>
<td>real(r8)</td>
<td>Specified in meters. If "match_vertical" is .FALSE. this value
is ignored. If "match_vertical" is .TRUE., this applies to
observations with a vertical type of VERTISSURFACE. For
observations which match in the horizontal, the vertical surface
elevation difference must be less than this to be considered the
same.</td>
</tr>
<tr>
<td>pressure_tolerance</td>
<td>real(r8)</td>
<td>Specified in pascals. If "match_vertical" is .FALSE. this value
is ignored. If "match_vertical" is .TRUE., this applies to
observations with a vertical type of VERTISPRESSURE. For
observations which match in the horizontal, the vertical difference
must be less than this to be considered the same.</td>
</tr>
<tr>
<td>height_tolerance</td>
<td>real(r8)</td>
<td>Specified in meters. If "match_vertical" is .FALSE. this value
is ignored. If "match_vertical" is .TRUE., this applies to
observations with a vertical type of VERTISHEIGHT. For observations
which match in the horizontal, the vertical difference must be less
than this to be considered the same.</td>
</tr>
<tr>
<td>scaleheight_tolerance</td>
<td>real(r8)</td>
<td>Specified in unitless values. If "match_vertical" is .FALSE.
this value is ignored. If "match_vertical" is .TRUE., this applies
to observations with a vertical type of VERTISSCALEHEIGHT. For
observations which match in the horizontal, the vertical difference
must be less than this to be considered the same.</td>
</tr>
<tr>
<td>level_tolerance</td>
<td>real(r8)</td>
<td>Specified in fractional model levels. If "match_vertical" is
.FALSE. this value is ignored. If "match_vertical" is .TRUE., this
applies to observations with a vertical type of VERTISLEVEL. For
observations which match in the horizontal, the vertical difference
must be less than this to be considered the same. Note that some
models only support integer level values, but others support
fractional levels. The vertical value in an observation is a
floating point/real value, so fractional levels are possible to
specify for an observation.</td>
</tr>
<tr>
<td>print_only</td>
<td>logical</td>
<td>If .TRUE. do not create an output file, but print a summary of
the number and types of each observation in each input file, and
then the number of observations and types which would have been
created in an output file.</td>
</tr>
<tr>
<td>partial_write</td>
<td>logical</td>
<td>Generally only used for debugging problems. After each input
obs_seq file is processed, this flag, if .TRUE., causes the code to
write out the partial results to the output file. The default is to
process all input files (if more than a single file is specified)
and write the output file only at the end of the processing.</td>
</tr>
<tr>
<td>print_timestamps</td>
<td>logical</td>
<td>Generally only used for debugging very slow execution runs.
This flag, if .TRUE., causes the code to output timestamps (wall
clock time) at various locations during the processing phases. It
may help isolate where particularly slow execution times are
occurring. For very large input files, or long lists of input
files, it can also help to estimate what the eventual run time of
the job will be.</td>
</tr>
<tr>
<td>calendar</td>
<td>character(len=32)</td>
<td>Set to the name of the calendar; only controls the printed
output for the dates of the first and last observations in the
file. Set this to "no_calendar" if the observations are not using
any calendar.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="Building" id="Building"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>BUILDING</h2>
<p>Most <em class="file">$DART/models/*/work</em> directories
contain files needed to build this tool along with the other
executable programs. It is also possible to build this tool in the
<em class="file">$DART/observations/utilities</em> directory. In
either case the <em class="file">preprocess</em> program must be
built and run first to define what set of observation types will be
supported. See the <a href=
"../../../assimilation_code/programs/preprocess/preprocess.html">preprocess
documentation</a> for more details on how to define the list and
run it. The <em class="code">&amp;preprocess_nml</em> namelist in
the <em class="file">input.nml</em> file must contain files with
definitions for the combined set of all observation types which
will be encountered over all input obs_seq files. The other
important choice when building the tool is to include a compatible
locations module in the <em class=
"file">path_names_obs_selection</em> file. For the low-order models
the <em class="file">oned</em> module should be used; for
real-world observations the <em class="file">threed_sphere</em>
module should be used.</p>
<p>Usually the directories where executables are built will include
a <em class="file">quickbuild.csh</em> script which builds and runs
preprocess and then builds the rest of the executables by executing
all files with names starting with <em class="file">mkmf_</em>. If
the obs_selection tool is not built because there is no <em class=
"file">mkmf_obs_selection</em> and <em class=
"file">path_names_obs_selection</em> file in the current directory
they can be copied from another model. The <em class=
"file">path_names_obs_selection</em> file will need to be edited to
be consistent with the model you are building.</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
time_manager_mod
obs_def_mod
obs_sequence_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li><em class="file">input.nml</em></li>
<li>The input files specified in the <em class=
"code">filename_seq</em> namelist variable.</li>
<li>The output file specified in the <em class=
"code">filename_out</em> namelist variable.</li>
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
<td valign="top">obs_selection</td>
<!-- message -->
<td valign="top">num_input_files &gt; max_num_input_files. change
max_num_input_files in source file</td>
<!-- comment -->
<td valign="top">The default is 500 files.</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_selection</td>
<!-- message -->
<td valign="top">num_input_files and filename_seq mismatch</td>
<!-- comment -->
<td valign="top">The number of filenames does not match the
filename count.</td>
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
<p>Long laundry list of things this tool could do, including:</p>
<ul>
<li>Remove duplicates. Lots of clever bookkeeping will be needed to
make sure this doesn't run excruciatingly slow.</li>
<li>General superob functions could be added - all observations
within a given tolerance could be combined. Hard questions are how
to specify the error of the combined observation, and how to
combine obs near or over the poles.</li>
</ul>
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
