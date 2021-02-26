<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program obs_common_subset</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>program <em class="program">obs_common_subset</em></h1>
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
<p>This specialized tool allows you to select subsets of
observations from two or more observation sequence files output
from <em class="file">filter</em>. It creates a new set of output
observation sequence files containing only the observations which
were successfully assimilated in all experiments.</p>
<p>Experiments using the same input observation sequence file but
with different configurations (e.g. different inflation values,
different localization radii, etc) can assimilate different numbers
of the available observations. In that case there will be
differences in the diagnostic plots which are not directly
relatable to the differences in the quality of the assimilation. If
this tool is run on the <em class="file">obs_seq.final</em> files
from all the experiments and then the diagnostics are generated,
only the observations which were assimilated in all experiments
will contribute to the summary statistics. A more direct comparison
can be made and improvements can be correctly attributed to the
differences in the experimental parameters.</p>
<p>This tool is intended to be used when comparing the results from
a group of related experiments in which <strong>the exact same
input observation sequence file</strong> is used for all runs. The
tool cannot process observation sequence files which differ in
anything other than whether an observation was successfully
assimilated/evaluated or not. Note that it is fine to add or remove
observation types from the <em class=
"code">assimilate_these_obs_types</em> or <em class=
"code">evaluate_these_obs_types</em> namelist items for different
experiments. The output observation sequence files will still
contain an identical list of observations, with some marked with a
DART QC indicating 'not assimilated because of namelist
control'.</p>
<p>See the <a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Documentation.php#obs_diagnostics">
two experiment diagnostic plot</a> documentation for Matlab scripts
supplied with DART to directly compare the observation diagnostic
output from multiple experiments (it does more than two, the script
has a poor name).</p>
<p>This is one of a set of tools which operate on observation
sequence files. For a more general purpose tool see the <a href=
"../obs_sequence_tool/obs_sequence_tool.html">obs_sequence_tool</a>,
and for a more flexible selection tool see the obs_selection_tool.</p>
<h4>Creating an Input Filelist</h4>
<p>One of the inputs to this tool is a list of filenames to
compare. The filenames can be directly in the namelist file, or
they can be in a set of separate text files. The latter may be
easier when there are more than just a few files to compare.</p>
<p>For experiments where there are multiple job steps, and so
multiple output observation sequence files per experiment, the
input to this tool would then be a list of lists of filenames. Each
set of names must be put into a text file with each filename on a
separate line.</p>
<p>If each experiment was run in a different set of directories,
and if a list of observation sequence filenames was made with the
<em class="file">ls</em> command:</p>
<pre>
&gt; ls exp1/*/obs_seq.final &gt; exp1list
&gt; cat exp1list
exp1/step1/obs_seq.final
exp1/step2/obs_seq.final
exp1/step3/obs_seq.final
exp1/step4/obs_seq.final
&gt; ls exp2/*/obs_seq.final &gt; exp2list
&gt; cat exp2list
exp2/step1/obs_seq.final
exp2/step2/obs_seq.final
exp2/step3/obs_seq.final
exp2/step4/obs_seq.final
&gt; ls exp3/*/obs_seq.final &gt; exp3list
&gt; cat exp2list
exp3/step1/obs_seq.final
exp3/step2/obs_seq.final
exp3/step3/obs_seq.final
exp3/step4/obs_seq.final
</pre>
<p>Then the namelist entries would be:</p>
<pre>
 filename_seq = ''
 filename_seq_list = 'exp1list', 'exp2list', exp3list'
 num_to_compare_at_once = 3
</pre>
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
&amp;obs_common_subset_nml
 num_to_compare_at_once = 2,
 filename_seq           = '',
 filename_seq_list      = '',
 filename_out_suffix    = '.common' ,
 print_every            = 10000,
 dart_qc_threshold      = 3,
 calendar               = 'Gregorian',
 print_only             = .false.,
 eval_and_assim_can_match = .false.,
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
<td>num_to_compare_at_once</td>
<td>integer</td>
<td>Number of observation sequence files to compare together at a
time. Most commonly the value is 2, but can be any number. If more
than this number of files are listed as inputs, the tool will loop
over the list N files at a time.</td>
</tr>
<tr>
<td>filename_seq</td>
<td>character(len=256), dimension(5000)</td>
<td>The array of names of the observation sequence files to
process. If more than N files (where N is num_to_compare_at_once)
are listed, they should be ordered so the first N files are
compared together, followed by the next set of N files, etc. You
can only specify one of filename_seq OR filename_seq_list, not
both.</td>
</tr>
<tr>
<td>filename_seq_list</td>
<td>character(len=256), dimension(100)</td>
<td>An alternative way to specify the list of input observation
sequence files. Give a list of N filenames which contain, one per
line, the names of the observation sequence files to process. There
should be N files specified (where N is num_to_compare_at_once),
and the first observation sequence filename listed in each file
will be compared together, then the second, until the lists are
exhausted. You can only specify one of filename_seq OR
filename_seq_list, not both.</td>
</tr>
<tr>
<td>filename_out_suffix</td>
<td>character(len=32)</td>
<td>A string to be appended to each of the input observation
sequence file names to create the output filenames.</td>
<td></td>
</tr>
<tr>
<td>print_every</td>
<td>integer</td>
<td>To indicate progress, a count of the successfully processed
observations is printed every Nth set of obs. To decrease the
output volume set this to a larger number. To disable this output
completely set this to -1.</td>
</tr>
<tr>
<td>dart_qc_threshold</td>
<td>integer</td>
<td>Observations with a DART QC value larger than this threshold
will be discarded. Note that this is the QC value set by <em class=
"file">filter</em> to indicate the outcome of trying to assimilate
an observation. This is not related to the incoming data QC. For an
observation which was successfully assimilated or evaluated in both
the Prior and Posterior this should be set to 1. To also include
observations which were successfully processed in the Prior but not
the Posterior, set to 3. To ignore the magnitude of the DART QC
values and keep observations only if the DART QCs match, set this
to any value higher than 7.</td>
</tr>
<tr>
<td>calendar</td>
<td>character(len=32)</td>
<td>Set to the name of the calendar; only controls the printed
output for the dates of the first and last observations in the
file. Set this to "no_calendar" if the observations are not using
any calendar.</td>
</tr>
<tr>
<td>print_only</td>
<td>logical</td>
<td>If .TRUE. do not create the output files, but print a summary
of the number and types of each observation in each of the input
and output files.</td>
</tr>
<tr>
<td>eval_and_assim_can_match</td>
<td>logical</td>
<td>Normally .FALSE. . If .TRUE. then observations which were
either successfully evaluated OR assimilated will match and are
kept.</td>
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
<p>Most <tt>$DART/models/*/work</tt> directories will build the
tool along with other executable programs. It is also possible to
build the tool in the <tt>$DART/observations/utilities</tt>
directory. The <tt>preprocess</tt> program must be built and run
first, to define what set of observation types will be supported.
See the <a href=
"../../../assimilation_code/programs/preprocess/preprocess.html">preprocess
documentation</a> for more details on how to define the list and
run it. The combined list of all observation types which will be
encountered over all input files must be in the preprocess input
list. The other important choice when building the tool is to
include a compatible locations module. For the low-order models,
the <tt>oned</tt> module should be used; for real-world
observations, the <tt>threed_sphere</tt> module should be used.</p>
<p>Generally the directories where executables are built will
include a "quickbuild.csh" script which will build and run
preprocess and then build the rest of the executables. The
"input.nml" namelists will need to be edited to include all the
required observation types first.</p>
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
"code">filename_seq</em> or <em class="code">filename_seq_list</em>
namelist variable.</li>
<li>The output files are specified by appending the string from the
<em class="code">filename_out_suffix</em> namelist item to the
input filenames.</li>
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
<td valign="top">obs_common_subset</td>
<!-- message -->
<td valign="top">num_input_files &gt; max_num_input_files.</td>
<!-- comment -->
<td valign="top">The default is 5000 total files. To process more,
change max_num_input_files in source code</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_common_subset</td>
<!-- message -->
<td valign="top">num_to_compare_at_once and filename_seq length
mismatch</td>
<!-- comment -->
<td valign="top">The number of filenames is not an even multiple of
the count.</td>
</tr>
<tr><!-- routine -->
<td valign="top">handle_filenames</td>
<!-- message -->
<td valign="top">cannot specify both filename_seq and
filename_seq_list</td>
<!-- comment -->
<td valign="top">You can either specify the files directly in the
namelist, or give a filename that contains the list of input files,
but not both.</td>
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
<p>none at this time.</p>
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
