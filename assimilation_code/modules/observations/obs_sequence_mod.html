<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_sequence_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE obs_sequence_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Interface">INTERFACES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Provides interfaces to the observation type and observation
sequence type. An observation contains everything there is to know
about an observation including all metadata contained in the
observation definition and any number of copies of data associated
with the observation (for instance an actual observation, an
ensemble of first guess values, etc). An observation sequence is a
time-ordered set of observations that is defined by a linked list
so that observations can be easily added or deleted. A number of
commands to extract observations depending on the times at which
they were taken are provided. For now, the observations are only
ordered by time, but the ability to add extra sort keys could be
added.</p>
<p>These routines are commonly used in conversion programs which
read observation data from various formats and create a DART
observation sequence in memory, and then write it out to a file.
See the <a href=
"../../../observations/obs_converters/README.md">observations</a>
directory for examples of programs which create and manipulate
observations using this routines.</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a><br>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
location_mod (depends on model_choice)
obs_def_mod
time_manager_mod
utilities_mod
obs_kind_mod
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
<td><em class="call">use obs_sequence_mod, only :</em></td>
<td><a href="#obs_sequence_type">obs_sequence_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_obs_sequence">init_obs_sequence</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#interactive_obs_sequence">interactive_obs_sequence</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_num_copies">get_num_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_num_qc">get_num_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_num_obs">get_num_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_max_num_obs">get_max_num_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_copy_meta_data">get_copy_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_qc_meta_data">get_qc_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_next_obs">get_next_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_prev_obs">get_prev_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_next_obs_from_key">get_next_obs_from_key</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_prev_obs_from_key">get_prev_obs_from_key</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#insert_obs_in_seq">insert_obs_in_seq</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#delete_obs_from_seq">delete_obs_from_seq</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_copy_meta_data">set_copy_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_qc_meta_data">set_qc_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_first_obs">get_first_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_last_obs">get_last_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#add_copies">add_copies</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#add_qc">add_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_obs_seq">write_obs_seq</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_obs_seq">read_obs_seq</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#append_obs_to_seq">append_obs_to_seq</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_obs_from_key">get_obs_from_key</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_obs_time_range">get_obs_time_range</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_obs">set_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_time_range_keys">get_time_range_keys</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_num_times">get_num_times</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#static_init_obs_sequence">static_init_obs_sequence</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#destroy_obs_sequence">destroy_obs_sequence</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_obs_seq_header">read_obs_seq_header</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_expected_obs">get_expected_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#delete_seq_head">delete_seq_head</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#delete_seq_tail">delete_seq_tail</a></td>
</tr>
<tr>
<td> </td>
<td></td>
</tr>
<tr>
<td> </td>
<td>LINKS BELOW FOR OBS_TYPE INTERFACES</td>
</tr>
<tr>
<td> </td>
<td></td>
</tr>
<tr>
<td> </td>
<td><a href="#obs_type">obs_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_obs">init_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#destroy_obs">destroy_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_obs_def">get_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_obs_def">set_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_obs_values">get_obs_values</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_obs_values">set_obs_values</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#replace_obs_values">replace_obs_values</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_qc">get_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_qc">set_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#replace_qc">replace_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_obs">write_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_obs">read_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#interactive_obs">interactive_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#copy_obs">copy_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#copy_obs">assignment(=)</a></td>
</tr>
</table>
<!--=================== DESCRIPTION OF A LOCAL TYPE ==================-->
<a name="obs_sequence_type" id="obs_sequence_type"></a><br>
<div class="type">
<pre>
<em class="call">type obs_sequence_type</em>
   private
   integer                       :: num_copies
   integer                       :: num_qc
   integer                       :: num_obs
   integer                       :: max_num_obs
   character(len=64), pointer    :: copy_meta_data(:)
   character(len=64), pointer    :: qc_meta_data(:)
   integer                       :: first_time
   integer                       :: last_time
   type(obs_type), pointer       :: obs(:)
end type obs_sequence_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>The obs_sequence type represents a series of observations
including multiple copies of data and quality control fields and
complete metadata about the observations. The sequence is organized
as an integer pointer linked list using a fixed array of storage
for obs (type obs_type). Each observation points to the previous
and next observation in time order (additional sort keys could be
added if needed) and has a unique integer key (see obs_type below).
The maximum number of observations in the sequence is represented
in the type as max_num_obs, the current number of observations is
in num_obs. The number of quality control (qc) fields per
observation is num_qc and the number of data values associated with
each observation is num_copies. Metadata for each copy of the data
is in copy_meta_data and metadata for the qc fields is in
qc_meta_data. The first and last pointers into the time linked list
are in first_time and last_time. A capability to write and read an
obs_sequence structure to disk is available. At present, the entire
observation sequence is read in to core memory. An on-disk
implementation may be necessary for very large observational
datasets.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">num_copies</td>
<td>Number of data values associated with each observation.</td>
</tr>
<tr>
<td valign="top">num_qc</td>
<td>Number of qc fields associated with each observation.</td>
</tr>
<tr>
<td valign="top">num_obs</td>
<td>Number of observations currently in sequence.</td>
</tr>
<tr>
<td valign="top">max_num_obs</td>
<td>Upper bounds on number of observations in sequence.</td>
</tr>
<tr>
<td valign="top">copy_meta_data</td>
<td>Text describing each copy of data associated with
observations.</td>
</tr>
<tr>
<td valign="top">qc_meta_data</td>
<td>Text describing each quality control field.</td>
</tr>
<tr>
<td valign="top">first_time</td>
<td>Location of first observation in sequence.</td>
</tr>
<tr>
<td valign="top">last_time</td>
<td>Location of last observation in sequence.</td>
</tr>
<tr>
<td valign="top">obs</td>
<td>Storage for all of the observations in the sequence.</td>
</tr>
</table>
</div>
<br>
<!--=================== DESCRIPTION OF A LOCAL TYPE ==================-->
 <a name="obs_type" id="obs_type"></a><br>
<div class="type">
<pre>
<em class="call">type obs_type</em>
   private
   integer            :: key
   type(obs_def_type) :: def
   real(r8), pointer  :: values(:)
   real(r8), pointer  :: qc(:)
   integer            :: prev_time
   integer            :: next_time
   integer            :: cov_group
end type obs_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Structure to represent everything known about a given
observation and to help with storing the observation in the
observation sequence structure (see above). The prev_time and
next_time are integer pointers that allow a linked list sorted on
time to be constructed. If needed, other sort keys could be
introduced (for instance by time available?). Each observation in a
sequence has a unique key and each observation has an obs_def_type
that contains all the definition and metadata for the observation.
A set of values is associated with the observation along with a set
of qc fields. The cov_group is not yet implemented but will allow
non-diagonal observation error covariances in a future release.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">key</td>
<td>Unique integer key when in an obs_sequence.</td>
</tr>
<tr>
<td valign="top">def</td>
<td>The definition of the observation (see obs_def_mod).</td>
</tr>
<tr>
<td valign="top">values</td>
<td>Values associated with the observation.</td>
</tr>
<tr>
<td valign="top">qc</td>
<td>Quality control fields associated with the observation.</td>
</tr>
<tr>
<td valign="top">prev_time</td>
<td>When in an obs_sequence, points to previous time sorted
observation.</td>
</tr>
<tr>
<td valign="top">next_time</td>
<td>When in an obs_sequence, points to next time sorted
observation.</td>
</tr>
<tr>
<td valign="top">cov_group</td>
<td>Not currently implemented.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="init_obs_sequence" id="init_obs_sequence"></a><br>
<div class="routine"><em class="call">call init_obs_sequence(seq,
num_copies, num_qc, expected_max_num_obs)</em>
<pre>
type(obs_sequence_type), intent(out) :: <em class="code">seq</em>
integer,                 intent(in)  :: <em class=
"code">num_copies</em>
integer,                 intent(in)  :: <em class=
"code">num_qc</em>
integer,                 intent(in)  :: <em class=
"code">expected_max_num_obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Constructor to create a variable of obs_sequence_type. This
routine must be called before using an obs_sequence_type. The
number of copies of the data to be associated with each observation
(for instance the observation from an instrument, an ensemble of
prior guesses, etc.) and the number of quality control fields
associated with each observation must be specified. Also, an
estimated upper bound on the number of observations to be stored in
the sequence is helpful in making creation of the sequence
efficient.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>The observation sequence being constructed</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies  </em></td>
<td>Number of copies of data to be associated with each
observation</td>
</tr>
<tr>
<td valign="top"><em class="code">num_qc  </em></td>
<td>Number of quality control fields associated with each
observation</td>
</tr>
<tr>
<td valign="top"><em class=
"code">expected_max_num_obs  </em></td>
<td>An estimate of the largest number of observations the sequence
might contain</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="interactive_obs_sequence" id=
"interactive_obs_sequence"></a><br>
<div class="routine"><em class="call">var =
interactive_obs_sequence()</em>
<pre>
type(obs_sequence_type) :: <em class=
"code">interactive_obs_sequence</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Uses input from standard in to create an observation sequence.
Initialization of the sequence is handled by the function.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>An observation sequence created from standard input</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_num_copies" id="get_num_copies"></a><br>
<div class="routine"><em class="call">var =
get_num_copies(seq)</em>
<pre>
integer                             :: <em class=
"code">get_num_copies</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns number of copies of data associated with each
observation in an observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns number of copies of data associated with each
observation in sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_num_qc" id="get_num_qc"></a><br>
<div class="routine"><em class="call">var = get_num_qc(seq)</em>
<pre>
integer                             :: <em class=
"code">get_num_qc</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns number of quality control fields associated with each
observation in an observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns number of quality control fields associated with each
observation in sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_num_obs" id="get_num_obs"></a><br>
<div class="routine"><em class="call">var = get_num_obs(seq)</em>
<pre>
integer                             :: <em class=
"code">get_num_obs</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns number of observations currently in an observation
sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns number of observations currently in an observation
sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_max_num_obs" id="get_max_num_obs"></a><br>
<div class="routine"><em class="call">var =
get_max_num_obs(seq)</em>
<pre>
integer                             :: <em class=
"code">get_max_num_obs</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns maximum number of observations an observation sequence
can hold.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns maximum number of observations an observation sequence
can hold</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_copy_meta_data" id="get_copy_meta_data"></a><br>
<div class="routine"><em class="call">var = get_copy_meta_data(seq,
copy_num)</em>
<pre>
character(len=64)                   :: <em class=
"code">get_copy_meta_data</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
integer,                 intent(in) :: <em class=
"code">copy_num</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns metadata associated with a given copy of data in an
observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns metadata associated with a copy of data in observation
sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">copy_num  </em></td>
<td>Return metadata for this copy</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_qc_meta_data" id="get_qc_meta_data"></a><br>
<div class="routine"><em class="call">var =
get_qc_meta_data(seq,qc_num)</em>
<pre>
character(len=64)                   :: <em class=
"code">get_qc_meta_data</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
integer,                 intent(in) :: <em class="code">qc_num</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns metadata associated with a given copy of quality control
fields associated with observations in an observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns metadata associated with a given qc copy</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">qc_num  </em></td>
<td>Return metadata for this copy</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_next_obs" id="get_next_obs"></a><br>
<div class="routine"><em class="call">call get_next_obs(seq, obs,
next_obs, is_this_last)</em>
<pre>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
type(obs_type),          intent(in)  :: <em class="code">obs</em>
type(obs_type),          intent(out) :: <em class=
"code">next_obs</em>
logical,                 intent(out) :: <em class=
"code">is_this_last</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an observation in a sequence, returns the next observation
in the sequence. If there is no next observation, is_this_last is
set to true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Find the next observation after this one</td>
</tr>
<tr>
<td valign="top"><em class="code">next_obs  </em></td>
<td>Return the next observation here</td>
</tr>
<tr>
<td valign="top"><em class=
"code">is_this_last  </em></td>
<td>True if obs is the last obs in sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_prev_obs" id="get_prev_obs"></a><br>
<div class="routine"><em class="call">call get_prev_obs(seq, obs,
prev_obs, is_this_first)</em>
<pre>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
type(obs_type),          intent(in)  :: <em class="code">obs</em>
type(obs_type),          intent(out) :: <em class=
"code">prev_obs</em>
logical,                 intent(out) :: <em class=
"code">is_this_first</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an observation in a sequence, returns the previous
observation in the sequence. If there is no previous observation,
is_this_first is set to true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Find the previous observation before this one</td>
</tr>
<tr>
<td valign="top"><em class="code">prev_obs  </em></td>
<td>Return the previous observation here</td>
</tr>
<tr>
<td valign="top"><em class=
"code">is_this_first  </em></td>
<td>True if obs is the first obs in sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_next_obs_from_key" id=
"get_next_obs_from_key"></a><br>
<div class="routine"><em class="call">call
get_next_obs_from_key(seq, last_key_used, next_obs,
is_this_last)</em>
<pre>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
integer,                 intent(in)  :: <em class=
"code">last_key_used</em>
type(obs_type),          intent(out) :: <em class=
"code">next_obs</em>
logical,                 intent(out) :: <em class=
"code">is_this_last</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given the last key used in a sequence, returns the next
observation in the sequence. If there is no next observation,
is_this_last is set to true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class=
"code">last_key_used  </em></td>
<td>Find the next observation after this key</td>
</tr>
<tr>
<td valign="top"><em class="code">next_obs  </em></td>
<td>Return the next observation here</td>
</tr>
<tr>
<td valign="top"><em class=
"code">is_this_last  </em></td>
<td>True if obs is the last obs in sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_prev_obs_from_key" id=
"get_prev_obs_from_key"></a><br>
<div class="routine"><em class="call">call
get_prev_obs_from_key(seq, last_key_used, prev_obs,
is_this_first)</em>
<pre>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
integer,                 intent(in)  :: <em class=
"code">last_key_used</em>
type(obs_type),          intent(out) :: <em class=
"code">prev_obs</em>
logical,                 intent(out) :: <em class=
"code">is_this_first</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given the last key used in a sequence, returns the previous
observation in the sequence. If there is no previous observation,
is_this_first is set to true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class=
"code">last_key_used  </em></td>
<td>Find the previous observation before this key</td>
</tr>
<tr>
<td valign="top"><em class="code">prev_obs  </em></td>
<td>Return the previous observation here</td>
</tr>
<tr>
<td valign="top"><em class=
"code">is_this_first  </em></td>
<td>True if obs is the first obs in sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_obs_from_key" id="get_obs_from_key"></a><br>
<div class="routine"><em class="call">call get_obs_from_key(seq,
key, obs)</em>
<pre>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
integer,                 intent(in)  :: <em class="code">key</em>
type(obs_type),          intent(out) :: <em class="code">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Each entry in an observation sequence has a unique integer key.
This subroutine returns the observation given an integer key.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>Return the observation with this key</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>The returned observation</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="insert_obs_in_seq" id="insert_obs_in_seq"></a><br>
<div class="routine"><em class="call">call insert_obs_in_seq(seq,
obs <em class="optionalcode">[, prev_obs]</em>)</em>
<pre>
type(obs_sequence_type),  intent(inout) :: <em class=
"code">seq</em>
type(obs_type),           intent(inout) :: <em class=
"code">obs</em>
type(obs_type), optional, intent(in)    :: <em class=
"optionalcode">prev_obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Inserts an observation in a sequence in appropriate time order.
If the optional argument prev_obs is present, the new observation
is inserted directly after the prev_obs. If an incorrect prev_obs
is provided so that the sequence is no longer time ordered, bad
things will happen.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>An observation to be inserted in the sequence</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">prev_obs  </em></td>
<td>If present, says the new observation belongs immediately after
this one</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="delete_obs_from_seq" id="delete_obs_from_seq"></a><br>
<div class="routine"><em class="call">call delete_obs_from_seq(seq,
obs)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
type(obs_type),          intent(inout) :: <em class="code">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an observation and a sequence, removes the observation
with the same key from the observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>The observation to be deleted from the sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="set_copy_meta_data" id="set_copy_meta_data"></a><br>
<div class="routine"><em class="call">call set_copy_meta_data(seq,
copy_num, meta_data)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
integer,                 intent(in)    :: <em class=
"code">copy_num</em>
character(len=64),       intent(in)    :: <em class=
"code">meta_data</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the copy metadata for this copy of the observations in an
observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">copy_num  </em></td>
<td>Set metadata for this copy of data</td>
</tr>
<tr>
<td valign="top"><em class="code">meta_data  </em></td>
<td>The metadata</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="set_qc_meta_data" id="set_qc_meta_data"></a><br>
<div class="routine"><em class="call">call set_qc_meta_data(seq,
qc_num, meta_data)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
integer,                 intent(in)    :: <em class=
"code">qc_num</em>
character(len=64),       intent(in)    :: <em class=
"code">meta_data</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the quality control metadata for this copy of the qc in an
observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">qc_num  </em></td>
<td>Set metadata for this quality control field</td>
</tr>
<tr>
<td valign="top"><em class="code">meta_data  </em></td>
<td>The metadata</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_first_obs" id="get_first_obs"></a><br>
<div class="routine"><em class="call">var = get_first_obs(seq,
obs)</em>
<pre>
logical                              :: <em class=
"code">get_first_obs</em>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
type(obs_type),          intent(out) :: <em class="code">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the first observation in a sequence. If there are no
observations in the sequence, the function returns false, else
true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns false if there are no obs in sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>The first observation in the sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_last_obs" id="get_last_obs"></a><br>
<div class="routine"><em class="call">var = get_last_obs(seq,
obs)</em>
<pre>
logical                              :: <em class=
"code">get_last_obs</em>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
type(obs_type),          intent(out) :: <em class="code">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the last observation in a sequence. If there are no
observations in the sequence, the function returns false, else
true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns false if there are no obs in sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>The last observation in the sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="add_copies" id="add_copies"></a><br>
<div class="routine"><em class="call">call add_copies(seq,
num_to_add)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
integer,                 intent(in)    :: <em class=
"code">num_to_add</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Increases the number of copies of data associated with each
observation by num_to_add. The current implementation re-creates
the entire observation sequence by deallocating and reallocating
each entry with a larger size.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">num_to_add  </em></td>
<td>Number of copies of data to add</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="add_qc" id="add_qc"></a><br>
<div class="routine"><em class="call">call add_qc(seq,
num_to_add)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
integer,                 intent(in)    :: <em class=
"code">num_to_add</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Increases the number of quality control fields associated with
each observation by num_to_add. The current implementation
re-creates the entire observation sequence by deallocating and
reallocating each entry with a larger size.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">num_to_add  </em></td>
<td>Number of quality control fields to add</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="read_obs_seq" id="read_obs_seq"></a><br>
<div class="routine"><em class="call">call read_obs_seq(file_name,
add_copies, add_qc, add_obs, seq)</em>
<pre>
character(len=*),        intent(in)  :: <em class=
"code">file_name</em>
integer,                 intent(in)  :: <em class=
"code">add_copies</em>
integer,                 intent(in)  :: <em class=
"code">add_qc</em>
integer,                 intent(in)  :: <em class=
"code">add_obs</em>
type(obs_sequence_type), intent(out) :: <em class="code">seq</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Read an observation sequence from <em class=
"code">file_name</em>. The sequence will have enough space for the
number of observations in the file plus any additional space
requested by the "add_xx" args. It is more efficient to allocate
the additional space at create time rather than try to add it in
later. The arguments can specify that the caller wants to add
additional data copies associated with each observation, or to add
additional quality control fields, or to add space for additional
observations. The format of the file (<tt>formatted</tt> vs.
<tt>unformatted</tt>) has been automatically detected since the I
release. The obs_sequence file format with I and later releases has
a header that associates observation type strings with an integer
which was not present in previous versions. I format files are no
longer supported.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">file_name  </em></td>
<td>Read from this file</td>
</tr>
<tr>
<td valign="top"><em class="code">add_copies  </em></td>
<td>Add this number of copies of data to the obs_sequence on
file</td>
</tr>
<tr>
<td valign="top"><em class="code">add_qc  </em></td>
<td>Add this number of qc fields to the obs_sequence on file</td>
</tr>
<tr>
<td valign="top"><em class="code">add_obs  </em></td>
<td>Add space for this number of additional observations to the
obs_sequence on file</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>The observation sequence read in with any additional space</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="write_obs_seq" id="write_obs_seq"></a><br>
<div class="routine"><em class="call">call write_obs_seq(seq,
file_name)</em>
<pre>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
character(len=*),        intent(in) :: <em class=
"code">file_name</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Write the observation sequence to file file_name. The format is
controlled by the namelist parameter write_binary_obs_sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">file_name  </em></td>
<td>Write the sequence to this file</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="set_obs" id="set_obs"></a><br>
<div class="routine"><em class="call">call set_obs(seq,obs
<em class="optionalcode">[, key_in]</em>)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
type(obs_type),          intent(in)    :: <em class="code">obs</em>
integer, optional,       intent(in)    :: <em class=
"optionalcode">key_in</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an observation, copies this observation into the
observation sequence using the key specified in the observation. If
the optional key_in argument is present, the observation is instead
copied into this element of the observation sequence (and the key
is changed to be key_in).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Observation to be put in sequence</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">key_in  </em></td>
<td>If present, the obs is copied into this key of the
sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="append_obs_to_seq" id="append_obs_to_seq"></a><br>
<div class="routine"><em class="call">call append_obs_to_seq(seq,
obs)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
type(obs_type),          intent(inout) :: <em class="code">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Append an observation to an observation sequence. An error
results if the time of the observation is not equal to or later
than the time of the last observation currently in the
sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Append this observation to the sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_obs_time_range" id="get_obs_time_range"></a><br>
<div class="routine"><em class="call">call get_obs_time_range(seq,
time1, time2, key_bounds, num_keys, out_of_range <em class=
"optionalcode">[, obs]</em>)</em>
<pre>
type(obs_sequence_type),  intent(in)  :: <em class="code">seq</em>
type(time_type),          intent(in)  :: <em class=
"code">time1</em>
type(time_type),          intent(in)  :: <em class=
"code">time2</em>
integer, dimension(2),    intent(out) :: <em class=
"code">key_bounds</em>
integer,                  intent(out) :: <em class=
"code">num_keys</em>
logical,                  intent(out) :: <em class=
"code">out_of_range</em>
type(obs_type), optional, intent(in)  :: <em class=
"optionalcode">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a time range specified by a beginning and ending time,
find the keys that bound all observations in this time range and
the number of observations in the time range. The routine
get_time_range_keys can then be used to get a list of all the keys
in the range if desired. The logical out_of_range is returned as
true if the beginning time of the time range is after the time of
the latest observation in the sequence. The optional argument obs
can increase the efficiency of the search through the sequence by
indicating that all observations before obs are definitely at times
before the start of the time range.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">time1  </em></td>
<td>Lower time bound</td>
</tr>
<tr>
<td valign="top"><em class="code">time2  </em></td>
<td>Upper time bound</td>
</tr>
<tr>
<td valign="top"><em class="code">key_bounds  </em></td>
<td>Lower and upper bounds on keys that are in the time range</td>
</tr>
<tr>
<td valign="top"><em class="code">num_keys  </em></td>
<td>Number of keys in the time range</td>
</tr>
<tr>
<td valign="top"><em class=
"code">out_of_range  </em></td>
<td>Returns true if the time range is entirely past the time of the
last obs in sequence</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">obs  </em></td>
<td>If present, can start search for time range from this
observation</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_time_range_keys" id="get_time_range_keys"></a><br>
<div class="routine"><em class="call">call get_time_range_keys(seq,
key_bounds, num_keys, keys)</em>
<pre>
type(obs_sequence_type),      intent(in)  :: <em class=
"code">seq</em>
integer, dimension(2),        intent(in)  :: <em class=
"code">key_bounds</em>
integer,                      intent(in)  :: <em class=
"code">num_keys</em>
integer, dimension(num_keys), intent(out) :: <em class=
"code">keys</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given the keys of the observations at the start and end of a
time range and the number of observations in the time range (these
are returned by <em class="code">get_obs_time_range()</em>), return
a list of the keys of all observations in the time range. Combining
the two routines allows one to get a list of all observations in
any time range by key. The <em class="code">keys</em> array must be
at least <em class="code">num_keys</em> long to hold the return
values.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">key_bounds  </em></td>
<td>Keys of first and last observation in a time range</td>
</tr>
<tr>
<td valign="top"><em class="code">num_keys  </em></td>
<td>Number of obs in the time range</td>
</tr>
<tr>
<td valign="top"><em class="code">keys  </em></td>
<td>Output list of keys of all obs in the time range</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_num_times" id="get_num_times"></a><br>
<div class="routine"><em class="call">var = get_num_times(seq)</em>
<pre>
integer                             :: <em class=
"code">get_num_times</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the number of unique times associated with observations
in an observation sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Number of unique times for observations in a sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_num_key_range" id="get_num_key_range"></a><br>
<div class="routine"><em class="call">var = get_num_key_range(seq,
key1, key2)</em>
<pre>
integer                             :: <em class=
"code">get_num_key_range</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
integer, optional,       intent(in) :: <em class=
"code">key1, key2</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the number of observations between the two given keys.
The default key numbers are the first and last in the sequence
file. This routine can be used to count the actual number of
observations in a sequence and will be accurate even if the
sequence has been trimmed with delete_seq_head() or
delete_seq_tail().</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Number of unique times for observations in a sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">key1  </em></td>
<td>The starting key number. Defaults to the first observation in
the sequence.</td>
</tr>
<tr>
<td valign="top"><em class="code">key2  </em></td>
<td>The ending key number. Defaults to the last observation in the
sequence.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="static_init_obs_sequence" id=
"static_init_obs_sequence"></a><br>
<div class="routine"><em class="call">call
static_init_obs_sequence()</em></div>
<div class="indent1"><!-- Description -->
<p>Initializes the obs_sequence module and reads namelists. This
MUST BE CALLED BEFORE USING ANY OTHER INTERFACES.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="destroy_obs_sequence" id="destroy_obs_sequence"></a><br>
<div class="routine"><em class="call">call
destroy_obs_sequence(seq)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Releases all allocated storage associated with an observation
sequence.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="read_obs_seq_header" id="read_obs_seq_header"></a><br>
<div class="routine"><em class="call">call
read_obs_seq_header(file_name, num_copies, num_qc, num_obs,
max_num_obs, file_id, read_format, pre_I_format <em class=
"optionalcode">[, close_the_file]</em>)</em>
<pre>
character(len=*),   intent(in)  :: <em class="code">file_name</em>
integer,            intent(out) :: <em class="code">num_copies</em>
integer,            intent(out) :: <em class="code">num_qc</em>
integer,            intent(out) :: <em class="code">num_obs</em>
integer,            intent(out) :: <em class=
"code">max_num_obs</em>
integer,            intent(out) :: <em class="code">file_id</em>
character(len=*),   intent(out) :: <em class=
"code">read_format</em>
logical,            intent(out) :: <em class=
"code">pre_I_format</em>
logical, optional,  intent(in)  :: <em class=
"optionalcode">close_the_file</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Allows one to see the global metadata associated with an
observation sequence that has been written to a file without
reading the whole file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">file_name  </em></td>
<td>File contatining an obs_sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies  </em></td>
<td>Number of copies of data associated with each observation</td>
</tr>
<tr>
<td valign="top"><em class="code">num_qc  </em></td>
<td>Number of quality control fields associated with each
observation</td>
</tr>
<tr>
<td valign="top"><em class="code">num_obs  </em></td>
<td>Number of observations in sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">max_num_obs  </em></td>
<td>Maximum number of observations sequence could hold</td>
</tr>
<tr>
<td valign="top"><em class="code">file_id  </em></td>
<td>File channel/descriptor returned from opening the file</td>
</tr>
<tr>
<td valign="top"><em class="code">read_format  </em></td>
<td>Either the string <tt>'unformatted'</tt> or
<tt>'formatted'</tt></td>
</tr>
<tr>
<td valign="top"><em class=
"code">pre_I_format  </em></td>
<td>Returns .true. if the file was written before the observation
type string/index number table was added to the standard header
starting with the I release.</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">close_the_file  </em></td>
<td>If specified and .TRUE. close the file after the header has
been read. The default is to leave the file open.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="init_obs" id="init_obs"></a><br>
<div class="routine"><em class="call">call init_obs(obs,
num_copies, num_qc)</em>
<pre>
type(obs_type), intent(out) :: <em class="code">obs</em>
integer,        intent(in)  :: <em class="code">num_copies</em>
integer,        intent(in)  :: <em class="code">num_qc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initializes an obs_type variable. This allocates storage for the
observation type and creates the appropriate obs_def_type and
related structures. IT IS ESSENTIAL THAT OBS_TYPE VARIABLES BE
INITIALIZED BEFORE USE.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>An obs_type data structure to be initialized</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies  </em></td>
<td>Number of copies of data associated with observation</td>
</tr>
<tr>
<td valign="top"><em class="code">num_qc  </em></td>
<td>Number of qc fields associated with observation</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="destroy_obs" id="destroy_obs"></a><br>
<div class="routine"><em class="call">call destroy_obs(obs)</em>
<pre>
type(obs_type), intent(inout) :: <em class="code">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Destroys an observation variable by releasing all associated
storage.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>An observation variable to be destroyed</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_obs_def" id="get_obs_def"></a><br>
<div class="routine"><em class="call">call get_obs_def(obs,
obs_def)</em>
<pre>
type(obs_type),     intent(in)  :: <em class="code">obs</em>
type(obs_def_type), intent(out) :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Extracts the definition portion of an observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>An observation</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>The definition portion of the observation</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="set_obs_def" id="set_obs_def"></a><br>
<div class="routine"><em class="call">call set_obs_def(obs,
obs_def)</em>
<pre>
type(obs_type),     intent(out) :: <em class="code">obs</em>
type(obs_def_type), intent(in)  :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an observation and an observation definition, insert the
definition in the observation structure.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>An observation whose definition portion will be updated</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>The observation definition that will be inserted in obs</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_obs_values" id="get_obs_values"></a><br>
<div class="routine"><em class="call">call get_obs_values(obs,
values <em class="optionalcode">[, copy_indx]</em>)</em>
<pre>
type(obs_type),         intent(in)  :: <em class="code">obs</em>
real(r8), dimension(:), intent(out) :: <em class="code">values</em>
integer, optional,      intent(in)  :: <em class=
"optionalcode">copy_indx</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Extract copies of the data from an observation. If <em class=
"optionalcode">copy_indx</em> is present extract a single value
indexed by <em class="optionalcode">copy_indx</em> into <em class=
"code">values(1)</em>.  <em class=
"optionalcode">copy_indx</em> must be between 1 and
<tt>num_copies</tt>, inclusive. If <em class=
"optionalcode">copy_indx</em> is not present extract all copies of
data into the <em class="code">values</em> array which must be
<tt>num_copies</tt> long (See <a href="#get_num_copies"><em class=
"code">get_num_copies</em></a>.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Observation from which to extract values</td>
</tr>
<tr>
<td valign="top"><em class="code">values  </em></td>
<td>The values extracted</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">copy_indx  </em></td>
<td>If present extract only this copy, otherwise extract all
copies</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_qc" id="get_qc"></a><br>
<div class="routine"><em class="call">call get_qc(obs, qc
<em class="optionalcode">[, qc_indx]</em>)</em>
<pre>
type(obs_type),         intent(in)  :: <em class="code">obs</em>
real(r8), dimension(:), intent(out) :: <em class="code">qc</em>
integer, optional,      intent(in)  :: <em class=
"optionalcode">qc_indx</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Extract quality control fields from an observation. If
<em class="optionalcode">qc_indx</em> is present extract a single
field indexed by <em class="optionalcode">qc_indx</em> into
<em class="code">qc(1)</em>.  <em class=
"optionalcode">qc_indx</em> must be between 1 and <tt>num_qc</tt>,
inclusive. If <em class="optionalcode">qc_indx</em> is not present
extract all quality control fields into the <em class=
"code">qc</em> array which must be <tt>num_qc</tt> long (See
<a href="#get_num_qc"><em class="code">get_num_qc</em></a>.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Observation from which to extract qc field(s)</td>
</tr>
<tr>
<td valign="top"><em class="code">qc  </em></td>
<td>Extracted qc fields</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">qc_indx  </em></td>
<td>If present extract only this field, otherwise extract all qc
fields</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="set_obs_values" id="set_obs_values"></a><br>
<div class="routine"><em class="call">call set_obs_values(obs,
values <em class="optionalcode">[, copy_indx]</em>)</em>
<pre>
type(obs_type),         intent(out) :: <em class="code">obs</em>
real(r8), dimension(:), intent(in)  :: <em class="code">values</em>
integer, optional,      intent(in)  :: <em class=
"optionalcode">copy_indx</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set value(s) of data in this observation. If <em class=
"optionalcode">copy_indx</em> is present set the single value
indexed by <em class="optionalcode">copy_indx</em> to <em class=
"code">values(1)</em>.  <em class=
"optionalcode">copy_indx</em> must be between 1 and
<tt>num_copies</tt>, inclusive. If <em class=
"optionalcode">copy_indx</em> is not present set all copies of data
from the <em class="code">values</em> array which must be
<tt>num_copies</tt> long (See <a href="#get_num_copies"><em class=
"code">get_num_copies</em></a>.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Observation whose values are being set</td>
</tr>
<tr>
<td valign="top"><em class="code">values  </em></td>
<td>Array of value(s) to be set</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">copy_indx  </em></td>
<td>If present set only this copy of data, otherwise set all
copies</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="replace_obs_values" id="replace_obs_values"></a><br>
<div class="routine"><em class="call">call replace_obs_values(seq,
key, values <em class="optionalcode">[, copy_indx]</em>)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
integer,                 intent(in)    :: <em class="code">key</em>
real(r8), dimension(:),  intent(in)    :: <em class=
"code">values</em>
integer, optional,       intent(in)    :: <em class=
"optionalcode">copy_indx</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set value(s) of data in the observation from a sequence with the
given <em class="code">key</em>. If <em class=
"optionalcode">copy_indx</em> is present set the single value
indexed by <em class="optionalcode">copy_indx</em> to <em class=
"code">values(1)</em>.  <em class=
"optionalcode">copy_indx</em> must be between 1 and
<tt>num_copies</tt>, inclusive. If <em class=
"optionalcode">copy_indx</em> is not present set all copies of data
from the <em class="code">values</em> array which must be
<tt>num_copies</tt> long (See <a href="#get_num_copies"><em class=
"code">get_num_copies</em></a>.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>Sequence which contains observation to update</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>Key to select which observation</td>
</tr>
<tr>
<td valign="top"><em class="code">values  </em></td>
<td>Array of value(s) to be set</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">copy_indx  </em></td>
<td>If present set only this copy of data, otherwise set all
copies</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="set_qc" id="set_qc"></a><br>
<div class="routine"><em class="call">call set_qc(obs, qc
<em class="optionalcode">[, qc_indx]</em>)</em>
<pre>
type(obs_type),         intent(out) :: <em class="code">obs</em>
real(r8), dimension(:), intent(in)  :: <em class="code">qc</em>
integer, optional,      intent(in)  :: <em class=
"optionalcode">qc_indx</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the quality control fields in an observation. If <em class=
"optionalcode">qc_indx</em> is present set a single field indexed
by <em class="optionalcode">qc_indx</em> to <em class=
"code">qc(1)</em>.  <em class="optionalcode">qc_indx</em> must
be between 1 and <tt>num_qc</tt>, inclusive. If <em class=
"optionalcode">qc_indx</em> is not present set all quality control
fields from the <em class="code">qc</em> array which must be
<tt>num_qc</tt> long (See <a href="#get_num_qc"><em class=
"code">get_num_qc</em></a>.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Observation having its qc fields set</td>
</tr>
<tr>
<td valign="top"><em class="code">qc  </em></td>
<td>Input values of qc fields</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">qc_indx  </em></td>
<td>If present update only this field, otherwise update all qc
fields</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="replace_qc" id="replace_qc"></a><br>
<div class="routine"><em class="call">call replace_qc(seq, key, qc
<em class="optionalcode">[, qc_indx]</em>)</em>
<pre>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
integer,                 intent(in)    :: <em class="code">key</em>
real(r8), dimension(:),  intent(in)    :: <em class="code">qc</em>
integer, optional,       intent(in)    :: <em class=
"optionalcode">qc_indx</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set value(s) of the quality control fields in the observation
from a sequence with the given <em class="code">key</em>. If
<em class="optionalcode">qc_indx</em> is present set the single
value indexed by <em class="optionalcode">qc_indx</em> to
<em class="code">qc(1)</em>.  <em class=
"optionalcode">qc_indx</em> must be between 1 and <tt>num_qc</tt>,
inclusive. If <em class="optionalcode">qc_indx</em> is not present
set all quality control fields from the <em class="code">qc</em>
array which must be <tt>num_qc</tt> long (See <a href=
"#get_num_qc"><em class="code">get_num_qc</em></a>.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>Observation sequence containing observation to update</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>Key to select which observation</td>
</tr>
<tr>
<td valign="top"><em class="code">qc  </em></td>
<td>Input values of qc fields</td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">qc_indx  </em></td>
<td>If present, only update single qc field, else update all qc
fields</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="write_obs" id="write_obs"></a><br>
<div class="routine"><em class="call">call write_obs(obs, file_id,
num_copies, num_qc)</em>
<pre>
type(obs_type), intent(in) :: <em class="code">obs</em>
integer,        intent(in) :: <em class="code">file_id</em>
integer,        intent(in) :: <em class="code">num_copies</em>
integer,        intent(in) :: <em class="code">num_qc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes an observation and all its associated metadata to a disk
file that has been opened with a format consistent with the
namelist parameter <tt>write_binary_obs_sequence</tt>.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Observation to be written to file</td>
</tr>
<tr>
<td valign="top"><em class="code">file_id  </em></td>
<td>Channel open to file for writing</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies  </em></td>
<td>The number of copies of data associated with the observation to
be output</td>
</tr>
<tr>
<td valign="top"><em class="code">num_qc  </em></td>
<td>The number of qc fields associated with the observation to be
output</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="read_obs" id="read_obs"></a><br>
<div class="routine"><em class="call">call read_obs(file_id,
num_copies, add_copies, num_qc, add_qc, key, obs, read_format
<em class="optionalcode">[, max_obs]</em>)</em>
<pre>
integer,            intent(in)    :: <em class="code">file_id</em>
integer,            intent(in)    :: <em class=
"code">num_copies</em>
integer,            intent(in)    :: <em class=
"code">add_copies</em>
integer,            intent(in)    :: <em class="code">num_qc</em>
integer,            intent(in)    :: <em class="code">add_qc</em>
integer,            intent(in)    :: <em class="code">key</em>
type(obs_type),     intent(inout) :: <em class="code">obs</em>
character(len=*),   intent(in)    :: <em class=
"code">read_format</em>
integer, optional,  intent(in)    :: <em class=
"optionalcode">max_obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads an observation from an obs_sequence file. The number of
copies of data and the number of qc values associated with each
observation must be provided. If additional copies of data or
additional qc fields are needed, arguments allow them to be added.
WARNING: The key argument is no longer used and should be
removed.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">file_id  </em></td>
<td>Channel open to file from which to read</td>
</tr>
<tr>
<td valign="top"><em class="code">num_copies  </em></td>
<td>Number of copies of data associated with observation in
file</td>
</tr>
<tr>
<td valign="top"><em class="code">add_copies  </em></td>
<td>Number of additional copies of observation to be added</td>
</tr>
<tr>
<td valign="top"><em class="code">num_qc  </em></td>
<td>Number of qc fields associated with observation in file</td>
</tr>
<tr>
<td valign="top"><em class="code">add_qc  </em></td>
<td>Number of additional qc fields to be added</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>No longer used, should be deleted</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>The observation being read in</td>
</tr>
<tr>
<td valign="top"><em class="code">read_format  </em></td>
<td>Either the string <tt>'formatted'</tt> or
<tt>'unformatted'</tt></td>
</tr>
<tr>
<td valign="top"><em class=
"optionalcode">max_obs  </em></td>
<td>If present, specifies the largest observation key number in the
sequence. This is used only for additional error checks on the next
and previous obs linked list values.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="interactive_obs" id="interactive_obs"></a><br>
<div class="routine"><em class="call">call
interactive_obs(num_copies, num_qc, obs, key)</em>
<pre>
integer,        intent(in)    :: <em class="code">num_copies</em>
integer,        intent(in)    :: <em class="code">num_qc</em>
type(obs_type), intent(inout) :: <em class="code">obs</em>
integer,        intent(in)    :: <em class="code">key</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Use standard input to create an observation. The number of
values, number of qc fields, and an observation type-specific key
associated with the observation are input. (Note that the key here
is not the same as the key in an observation sequence.)</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">num_copies  </em></td>
<td>Number of copies of data to be associated with observation</td>
</tr>
<tr>
<td valign="top"><em class="code">num_qc  </em></td>
<td>Number of qc fields to be associated with observation</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>Observation created via standard input</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>An observation type-specific key can be associated with each
observation for use by the obs_def code.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="copy_obs" id="copy_obs"></a><br>
<div class="routine"><em class="call">call copy_obs(obs1,
obs2)</em>
<pre>
type(obs_type), intent(out) :: <em class="code">obs1</em>
type(obs_type), intent(in)  :: <em class="code">obs2</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Copies the observation type obs2 to obs1. If the sizes of obs
fields are not compatible, the space in obs1 is deallocated and
reallocated with the appropriate size. This is overloaded to
assignment(=).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs1  </em></td>
<td>Copy obs2 to here (destination)</td>
</tr>
<tr>
<td valign="top"><em class="code">obs2  </em></td>
<td>Copy into obs1 (source)</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="get_expected_obs" id="get_expected_obs"></a><br>
<div class="routine"><em class="call">call
get_expected_obs_from_def_distrib_state(state_handle, ens_size,
copy_indices, key, &amp; obs_def, obs_kind_ind, state_time,
isprior, assimilate_this_ob, evaluate_this_ob, expected_obs, &amp;
istatus)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=
"code">state_handle</em>
integer,             intent(in)  :: <em class="code">ens_size</em>
integer,             intent(in)  :: <em class=
"code">copy_indices(ens_size)</em>
integer,             intent(in)  :: <em class="code">key</em>
type(obs_def_type),  intent(in)  :: <em class="code">obs_def</em>
integer,             intent(in)  :: <em class=
"code">obs_kind_ind</em>
type(time_type),     intent(in)  :: <em class=
"code">state_time</em>
logical,             intent(in)  :: <em class="code">isprior</em>
integer,             intent(out) :: <em class=
"code">istatus(ens_size)</em>
logical,             intent(out) :: <em class=
"code">assimilate_this_ob, evaluate_this_ob</em>
real(r8),            intent(out) :: <em class=
"code">expected_obs(ens_size)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Used to compute the expected value of a set of observations in
an observation sequence given a model state vector. Also returns a
status variable that reports on problems taking forward operators.
This version returns forward operator values for the entire
ensemble in a single call.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state_handle</em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">keys</em></td>
<td>List of integer keys that specify observations in seq</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_index</em></td>
<td>The ensemble number for this state vector</td>
</tr>
<tr>
<td valign="top"><em class="code">state</em></td>
<td>Model state vector</td>
</tr>
<tr>
<td valign="top"><em class="code">state_time</em></td>
<td>The time of the state data</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_vals</em></td>
<td>Returned expected values of the observations</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Integer error code for use in quality control (0 means no
error)</td>
</tr>
<tr>
<td valign="top"><em class="code">assimilate_this_ob</em></td>
<td>Returns true if this observation type is being assimilated</td>
</tr>
<tr>
<td valign="top"><em class="code">evaluate_this_ob</em></td>
<td>Returns true if this observation type is being evaluated but
not assimilated</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="delete_seq_head" id="delete_seq_head"></a><br>
<div class="routine"><em class="call">call
delete_seq_head(first_time, seq, all_gone)</em>
<pre>
type(time_type),         intent(in)    :: <em class=
"code">first_time</em>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
logical,                 intent(out)   :: <em class=
"code">all_gone</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Deletes all observations in the sequence with times before
first_time. If no observations remain, return all_gone as .true. If
no observations fall into the time window (e.g. all before
first_time or empty sequence to begin with), no deletions are done
and all_gone is simply returned as .true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">first_time  </em></td>
<td>Delete all observations with times before this</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">all_gone  </em></td>
<td>Returns true if there are no valid observations remaining in
the sequence after first_time</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE ======================-->
 <a name="delete_seq_tail" id="delete_seq_tail"></a><br>
<div class="routine"><em class="call">call
delete_seq_tail(last_time, seq, all_gone)</em>
<pre>
type(time_type),         intent(in)    :: <em class=
"code">last_time</em>
type(obs_sequence_type), intent(inout) :: <em class="code">seq</em>
logical,                 intent(out)   :: <em class=
"code">all_gone</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Deletes all observations in the sequence with times after
last_time. If no observations remain, return all_gone as .true. If
no observations fall into the time window (e.g. all after last_time
or empty sequence to begin with), no deletions are done and
all_gone is simply returned as .true.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">last_time  </em></td>
<td>Delete all observations with times after this</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">all_gone  </em></td>
<td>Returns true if there are no valid observations remaining in
the sequence before last_time</td>
</tr>
</table>
</div>
<br>
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
&amp;obs_sequence_nml
   write_binary_obs_sequence = .false.
   read_binary_file_format   = 'native'
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
<td valign="top">write_binary_obs_sequence</td>
<td>logical</td>
<td>If true, write binary obs_sequence files. If false, write ascii
obs_sequence files.</td>
</tr>
<tr>
<td valign="top">read_binary_file_format</td>
<td>character(len=32)</td>
<td>The 'endian'ness of binary obs_sequence files. May be 'native'
(endianness matches hardware default), 'big-endian',
'little-endian', and possibly 'cray'. Ignored if observation
sequence files are ASCII.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>obs_sequence_mod.nml in input.nml</li>
<li>Files for reading and writing obs_sequences and obs specified
in filter_nml.</li>
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
<td valign="top">insert_obs_in_seq</td>
<!-- message -->
<td valign="top">ran out of room, num_obs # &gt; max_num_obs #</td>
<!-- comment -->
<td valign="top">Overflowed number of obs in sequence. Called from
many public entries.</td>
</tr>
<tr><!-- routine -->
<td valign="top">append_obs_to_seq</td>
<!-- message -->
<td valign="top">tried to append an obs to sequence with bad
time</td>
<!-- comment -->
<td valign="top">Tried to append an obs with earlier time than last
obs in sequence.</td>
</tr>
<tr><!-- routine -->
<td valign="top">append_obs_to_seq</td>
<!-- message -->
<td valign="top">ran out of room, max_num_obs = #</td>
<!-- comment -->
<td valign="top">Overflowed the obs sequence.</td>
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
<p>Future versions should automate file reading and only the write
namelist parameter should remain.</p>
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
