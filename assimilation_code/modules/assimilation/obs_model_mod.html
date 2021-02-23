<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_model_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE obs_model_mod</h1>
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
<p>The code in this module computes the assimilation windows, and
decides if the model needs to run in order for the data to be at
the appropriate time to assimilate the next available observations.
It also has the code to write out the current states, advance the
model (in a variety of ways) and then read back in the updated
states.</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
assim_model_mod
obs_sequence_mod
obs_def_mod
time_manager_mod
ensemble_manager_mod
mpi_utilities_mod
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
<td><em class="call">use obs_model_mod, only :</em></td>
<td><a href="#advance_state">advance_state</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#move_ahead">move_ahead</a></td>
</tr>
</table>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="move_ahead" id="move_ahead"></a><br>
<div class="routine"><em class="call">call move_ahead(ens_handle,
ens_size, seq, last_key_used, window_time, key_bounds,
num_obs_in_set, curr_ens_time, next_ens_time, trace_messages)</em>
<pre>
type(ensemble_type),     intent(in)  :: <em class=
"code">ens_handle</em>
integer,                 intent(in)  :: <em class=
"code">ens_size</em>
type(obs_sequence_type), intent(in)  :: <em class="code">seq</em>
integer,                 intent(in)  :: <em class=
"code">last_key_used</em>
type(time_type),         intent(in)  :: <em class=
"code">window_time</em>
integer, dimension(2),   intent(out) :: <em class=
"code">key_bounds</em>
integer,                 intent(out) :: <em class=
"code">num_obs_in_set</em>
type(time_type),         intent(out) :: <em class=
"code">curr_ens_time</em>
type(time_type),         intent(out) :: <em class=
"code">next_ens_time</em>
logical, optional,       intent(in)  :: <em class=
"code">trace_messages</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an observation sequence and an ensemble, determines how to
advance the model so that the next set of observations can be
assimilated. Also returns the first and last keys and the number of
observations to be assimilated at this time. The algorithm
implemented here (one might want to have other variants) first
finds the time of the next observation that has not been
assimilated at a previous time. It also determines the time of the
ensemble state vectors. It then uses information about the model's
time stepping capabilities to determine the time to which the model
can be advanced that is CLOSEST to the time of the next
observation. For now, this algorithm assumes that the model's
timestep is a constant. A window of width equal to the model
timestep is centered around the closest model time to the next
observation and all observations in this window are added to the
set to be assimilated.<br>
<br>
Previous versions of this routine also made the call which actually
advanced the model before returning. This is no longer true. The
routine only determines the time stepping and number of
observations. The calling code must then call advance_state() if
indeed the next observation to be assimilated is not within the
current window. This is determined by comparing the current
ensemble time with the next ensemble time. If equal no advance is
needed. Otherwise, next ensemble time is the target time for
advance_state().</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Identifies the model state ensemble</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size</em></td>
<td>Number of ensemble members</td>
</tr>
<tr>
<td valign="top"><em class="code">seq</em></td>
<td>An observation sequence</td>
</tr>
<tr>
<td valign="top"><em class="code">last_key_used</em></td>
<td>Identifies the last observation from the sequence that has been
used</td>
</tr>
<tr>
<td valign="top"><em class="code">window_time</em></td>
<td>Reserved for future use.</td>
</tr>
<tr>
<td valign="top"><em class="code">key_bounds</em></td>
<td>Returned lower and upper bound on observations to be used at
this time</td>
</tr>
<tr>
<td valign="top"><em class="code">num_obs_in_set</em></td>
<td>Number of observations to be used at this time</td>
</tr>
<tr>
<td valign="top"><em class="code">curr_ens_time</em></td>
<td>The time of the ensemble data passed into this routine.</td>
</tr>
<tr>
<td valign="top"><em class="code">next_ens_time</em></td>
<td>The time the ensemble data should be advanced to. If equal to
curr_ens_time, the model does not need to advance to assimilate the
next observation.</td>
</tr>
<tr>
<td valign="top"><em class="code">trace_messages</em></td>
<td>Optional argument. By default, detailed time trace messages are
disabled but can be turned on by passing this in as .True. . The
messages will print the current window times, data time, next
observation time, next window time, next data time, etc.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="advance_state" id="advance_state"></a><br>
<div class="routine"><em class="call">call
advance_state(ens_handle, ens_size, target_time, async,
adv_ens_command, tasks_per_model_advance)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer, intent(in)                :: <em class=
"code">ens_size</em>
type(time_type), intent(in)        :: <em class=
"code">target_time</em>
integer, intent(in)                :: <em class="code">async</em>
character(len=*), intent(in)       :: <em class=
"code">adv_ens_command</em>
integer, intent(in)                :: <em class=
"code">tasks_per_model_advance</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Advances all ensemble size copies of an ensemble stored in
ens_handle to the target_time. If async=0 this is done by repeated
calls to the <tt>adv_1step()</tt> subroutine. If async=2, a call to
the shell with the command <em class="code">adv_ens_command</em> is
used. If async=4, the filter program synchronizes with the MPI job
shell script using the <tt>block_task()</tt> and
<tt>restart_task()</tt> routines to suspend execution until all
model advances have completed. The script can start the model
advances using MPI and have it execute in parallel in this
mode.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Structure for holding ensemble information and data</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size</em></td>
<td>Ensemble size.</td>
</tr>
<tr>
<td valign="top"><em class="code">target_time</em></td>
<td>Time to which model is to be advanced.</td>
</tr>
<tr>
<td valign="top"><em class="code">async</em></td>
<td>How to advance model:
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td>0 = subroutine adv_1step</td>
</tr>
<tr>
<td>2 = shell executes adv_ens_command</td>
</tr>
<tr>
<td>4 = MPI job script advances models and syncs with filter
task</td>
</tr>
</table>
</td>
</tr>
<tr>
<td valign="top"><em class="code">adv_ens_command</em></td>
<td>Command to be issued to shell to advance model if async=2.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">tasks_per_model_advance   </em></td>
<td>Reserved for future use.</td>
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
<p>This module does not have a namelist.</p>
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
<td valign="top">assim_model_state_ic<em>####</em></td>
<td>a binary representation of the state vector prepended by a
small header consisting of the 'advance-to' time and the
'valid-time' of the state vector. The <em>####</em> represents the
ensemble member number if <em class=
"code">&amp;ensemble_manager_nml</em>: <em class=
"code">single_restart_file_out = .true.</em>.</td>
</tr>
<tr>
<td valign="top">
assim_model_state_ud<em>####   </em></td>
<td>a binary representation of the state vector prepended by a
small header consisting of the 'valid-time' of the state vector.
This is the 'updated' model state (after the model has advanced the
state to the desired 'advance-to' time).</td>
</tr>
<tr>
<td valign="top">filter_control<em>####</em></td>
<td>a text file containing information needed to advance the
ensemble members; i.e., the ensemble member number, the input state
vector file, the output state vector file - that sort of
thing.</td>
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
<td valign="top">move_ahead</td>
<!-- message -->
<td valign="top">next obs time not in model time window</td>
<!-- comment -->
<td valign="top">Error in algorithm to compute observation
window</td>
</tr>
<tr><!-- routine -->
<td valign="top">advance_state</td>
<!-- message -->
<td valign="top">target time ###,### is before model_time
###,###</td>
<!-- comment -->
<td valign="top">Target time must not be before current model
time.</td>
</tr>
<tr><!-- routine -->
<td valign="top">advance_state</td>
<!-- message -->
<td valign="top">Trying to use ### model states -- too many. Use
less than 10000 member ensemble.</td>
<!-- comment -->
<td valign="top">Maximum of 9999 ensemble members is allowed.</td>
</tr>
<tr><!-- routine -->
<td valign="top">advance_state</td>
<!-- message -->
<td valign="top">Can only have 10000 processes.</td>
<!-- comment -->
<td valign="top">No more than 9999 processes can run.</td>
</tr>
<tr><!-- routine -->
<td valign="top">advance_state</td>
<!-- message -->
<td valign="top">input.nml - async is #, must be 0, or 2.</td>
<!-- comment -->
<td valign="top">Only 0 or 2 work for async.</td>
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
<p>none at this time</p>
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
