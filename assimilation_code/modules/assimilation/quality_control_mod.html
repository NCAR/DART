<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module quality_control_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE quality_control_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Usage%20Notes">USAGE</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#Discussion">DISCUSSION</a> /
<a href="#Interface">INTERFACES</a> / <a href=
"#FilesUsed">FILES</a> / <a href="#References">REFERENCES</a> /
<a href="#Errors">ERRORS</a> / <a href="#FuturePlans">PLANS</a> /
<a href="#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Routines in this module deal with two different types of quality
control (QC) related functions. The first is to support
interpretation of the <em>incoming</em> data quality, to reject
observations at assimilation time which are marked as poor quality.
The second is to document how DART disposed of each observation;
whether it was successfully assimilated or rejected, and if
rejected, for which reason.</p>
<!--=====================================================================-->
<!--===================== USAGE NOTES ===================================-->
<!--=====================================================================-->
<a name="Usage Notes"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Usage</h2>
<h4>Incoming Data Quality Control</h4>
<p>DART currently supports a single incoming quality control scheme
compatible with NCEP usage. Lower values are considered better and
higher values are considered poorer. A single namelist item,
<em class="code">input_qc_threshold</em> sets the boundary between
accepted and rejected observations. Values <em>larger</em> than
this value are rejected; values equal to or lower are accepted.
Note that observations could be subsequently rejected for other
reasons, including failing the outlier threshold test or all
observations of this type being excluded by namelist control. See
the <a href=
"../observations/obs_kind_mod.html#Namelist">obs_kind_mod</a>
namelist documentation for more details on how to enable or disable
assimilation by observation type at runtime.</p>
<p>The incoming quality control value is set when an observation
sequence file is created. If the data provider user a different
scheme the values must be translated into NCEP-consistent values.
Generally we use the value 3 for most runs.</p>
<p>Observations can also be rejected by the assimilation if the
observation value is too far from the mean of the ensemble of
expected values (the forward operator results). This is controlled
by the <em class="code">outlier_threshold</em> namelist item.</p>
<p>Specifically, the outlier test computes the difference between
the observation value and the prior ensemble mean. It then computes
a standard deviation by taking the square root of the sum of the
observation error variance and the prior ensemble variance for the
observation. If the difference between the ensemble mean and the
observation value is more than the specified number of standard
deviations then the observation is not used. This can be an
effective way to discard clearly erroneous observation values. A
commonly used value is 3. -1 disables this check.</p>
<p>There is an option to add code to this module to specialize the
outlier threshold routine. For example, it is possible to allow all
observations of one type to be assimilated regardless of the
outlier value, and enforce the outlier threshold only on other
types of observations. To enable this capability requires two
actions: setting the <em class=
"code">enable_special_outlier_code</em> namelist to <em class=
"code">.TRUE.</em>, and adding your custom code to the <em class=
"code">failed_outlier()</em> subroutine in this module.</p>
<h4>DART Outgoing Quality Control</h4>
<p>As DART assimilates each observation it adds a <em>DART Quality
Control</em> value to the output observation sequence (frequently
written to a file named <em class="file">obs_seq.final)</em>. This
flag indicates how the observation was used during the
assimilation. The flag is a numeric value with the following
meanings:</p>
<table border="0" cellpadding="3" width="100%" summary=
'dart quality control flags'>
<tr>
<td>0:</td>
<td>Observation was assimilated successfully</td>
</tr>
<tr>
<td>1:</td>
<td>Observation was evaluated only so not used in the
assimilation</td>
</tr>
<tr>
<td>2:</td>
<td>The observation was used but one or more of the posterior
forward observation operators failed</td>
</tr>
<tr>
<td>3:</td>
<td>The observation was evaluated only so not used AND one or more
of the posterior forward observation operators failed</td>
</tr>
<tr>
<td>4:</td>
<td>One or more prior forward observation operators failed so the
observation was not used</td>
</tr>
<tr>
<td>5:</td>
<td>The observation was not used because it was not selected in the
namelist to be assimilated or evaluated</td>
</tr>
<tr>
<td>6:</td>
<td>The incoming quality control value was larger than the
threshold so the observation was not used</td>
</tr>
<tr>
<td>7:</td>
<td>Outlier threshold test failed (as described above)</td>
</tr>
<tr>
<td>8:</td>
<td>The location conversion to the vertical localization unit
failed so the observation was not used</td>
</tr>
</table>
<br>
<br>
<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->
 <a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Namelist</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;quality_control_nml
   input_qc_threshold          = 3
   outlier_threshold           = -1
   enable_special_outlier_code = .false.
  /
</pre></div>
<br>
<br>
<p>Items in this namelist control whether an observation is
assimilated or not.</p>
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
<td>input_qc_threshold</td>
<td>integer</td>
<td>Numeric value indicating whether this observation is considered
"good quality" and should be assimilated, or whether it is suspect
because of previous quality control processes. This value would
have been set when the observation was created and added to the
observation sequence file. Observations with an incoming QC value
larger than this threshold are rejected and not assimilated.</td>
</tr>
<tr>
<td>outlier threshold</td>
<td>integer</td>
<td>This numeric value defines the maximum number of standard
deviations an observation value can be away from the ensemble
forward operator mean and still be assimilated. Setting it to the
value -1 disables this check.</td>
</tr>
<tr>
<td>enable_special_outlier_code</td>
<td>logical</td>
<td>Setting this value to .TRUE. will call a subroutine <em class=
"code">failed_outlier()</em> instead of using the default code. The
user can then customize the tests in this subroutine, for example
to accept all observations of a particular type, or use different
numerical thresholds for different observation types or
locations.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--=====================================================================-->
<!--===================== DISCUSSION ====================================-->
<!--=====================================================================-->
 <a name="Discussion" id="Discussion"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Discussion</h2>
<h4>Small ensemble spread</h4>
<p>If an ensemble is spun up from a single state the ensemble
spread may be very small to begin and many observations may be
rejected by the <em class="code">outlier_threshold</em>. But as the
ensemble spread increases the assimilation should be able to
assimilate more and more observations as the model trajectory
becomes consistent with those observations.</p>
<!--==================================================================-->
<!--=====================================================================-->
<!--===================== INTERFACES ====================================-->
<!--=====================================================================-->
<!-- this section should be replaced by auto-generated code. -->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
random_seq_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<hr>
<h2>Public Interfaces</h2>
<table>
<tr>
<td><em class="code">use quality_control_mod, only :</em></td>
<td><a href="#initialize_qc">initialize_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#input_qc_ok">input_qc_ok</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_dart_qc">get_dart_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#check_outlier_threshold">check_outlier_threshold</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#good_dart_qc">good_dart_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_input_qc">set_input_qc</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#dart_flags">dart_flags</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="check_outlier_threshold" id=
"check_outlier_threshold"></a><br>
<div class="routine"><em class="call">call
check_outlier_threshold(obs_prior_mean, obs_prior_var, obs_val,
obs_err_var, &amp; obs_seq, this_obs_key, dart_qc)</em>
<pre>
real(r8),                intent(in)    :: obs_prior_mean !&gt;  prior observation mean
real(r8),                intent(in)    :: obs_prior_var  !&gt;  prior observation variance
real(r8),                intent(in)    :: obs_val        !&gt;  observation value
real(r8),                intent(in)    :: obs_err_var    !&gt;  observation error variance
type(obs_sequence_type), intent(in)    :: obs_seq        !&gt;  observation sequence
integer,                 intent(in)    :: this_obs_key   !&gt;  index for this observation
integer,                 intent(inout) :: dart_qc        !&gt;  possibly modified DART QC
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes whether this observation failed the outlier threshold
test and if so, updates the DART QC.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="input_qc_ok" id="input_qc_ok"></a><br>
<div class="routine"><em class="call">var = input_qc_ok(input_qc,
qc_to_use)</em>
<pre>
real(r8), intent(in)  :: input_qc    !&gt; incoming QC data value
integer,  intent(out) :: qc_to_use   !&gt; resulting DART QC
logical               :: input_qc_ok !&gt; true if input_qc is good
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if the input qc indicates this observation is good
to use.</p>
</div>
<br>
<!--============= DESCRIPTION OF A PUBLIC CONSTANT =================-->
 <a name="DART QC VALUES"></a><br>
<div class="routine">
<pre>
! Dart quality control variables
integer, parameter :: DARTQC_ASSIM_GOOD_FOP        = 0
integer, parameter :: DARTQC_EVAL_GOOD_FOP         = 1
integer, parameter :: DARTQC_ASSIM_FAILED_POST_FOP = 2
integer, parameter :: DARTQC_EVAL_FAILED_POST_FOP  = 3
integer, parameter :: DARTQC_FAILED_FOP            = 4
integer, parameter :: DARTQC_NOT_IN_NAMELIST       = 5
integer, parameter :: DARTQC_BAD_INCOMING_QC       = 6
integer, parameter :: DARTQC_FAILED_OUTLIER_TEST   = 7
integer, parameter :: DARTQC_FAILED_VERT_CONVERT   = 8
!!integer, parameter :: DARTQC_OUTSIDE_DOMAIN        = 9  ! we have no way (yet) for the model_mod to signal this
</pre></div>
<div class="indent1"><!-- Description -->
<p>These are public constants for use in other parts of the DART
code.</p>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Files</h2>
<table border="0">
<tr>
<th>filename</th>
<th>purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read the quality_control_mod namelist</td>
</tr>
</table>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>References</h2>
<ol>
<li>none</li>
</ol>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Error Codes and Conditions</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">routine name</td>
<!-- message -->
<td valign="top">output string</td>
<!-- comment -->
<td valign="top">description or comment</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Future Plans</h2>
<p>Should support different incoming data QC schemes.</p>
<p>It would be nice to have a different DART QC flag for
observations which fail the forward operator because they are
simply outside the model domain. The diagnosic routines may
indicate a large number of failed forward operators which make it
confusing to identify observations where the forward operator
should have been computed and can skew the statistics.
Unfortunately, this requires adding an additional requirement on
the model-dependent <em classs="file">model_mod.f90</em> code in
the <em class="code">model_interpolate()</em> routine. The current
interface defines a return status code of 0 as success, any
positive value as failure, and negative numbers are reserved for
other uses. To identify obs outside the domain would require
reserving another value that the interpolate routine could
return.</p>
<p>At this time the best suggestion is to cull out-of-domain obs
from the input observation sequence file by a preprocessing program
before assimilation.</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Private Components</h2>
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
