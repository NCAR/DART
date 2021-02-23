<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_def_1d_state_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE <em class="program">obs_def_1d_state_mod</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Interface">INTERFACES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>The list of observation types to be supported by the DART
executables is defined at compile time. The observations DART
supports can be changed at any time by adding or removing items
from the preprocess namelist and rerunning
<em>quickbuild.csh</em>.</p>
<p><em class="unix">Preprocess</em> takes observation specific code
sections from special obs_def files to generate <em class=
"file">obs_def_mod.f90</em> and <em class=
"file">obs_kind_mod.f90</em> which are then compiled into
<em class="unix">filter</em> and other DART programs. One of the
motivations behind creating <em class=
"file">obs_def_1d_state_mod.f90</em> was to provide a prototype for
people developing more complicated specialized observation
definition modules.<br>
<br>
<em class="file">Obs_def_1d_state_mod.f90</em> is an extended
format Fortran 90 module that provides the definition for
observation types designed for use with idealized low-order models
that use the 1D location module and can be thought of as having a
state vector that is equally spaced on a 1D cyclic domain.
Observation types include:</p>
<ul>
<li>RAW_STATE_VARIABLE - A straight linear interpolation to a point
on a [0,1] domain.</li>
<li>RAW_STATE_VAR_POWER - The interpolated RAW_STATE_VARIABLE
raised to a real-valued power.</li>
<li>RAW_STATE_1D_INTEGRAL - An area-weighted 'integral' of the
state variable over some part of the cyclic 1D domain.</li>
</ul>
RAW_STATE_VAR_POWER is convenient for studying non-gaussian,
non-linear assimilation problems. RAW_STATE_VAR_POWER can be used
to do idealized studies related to remote sensing observations that
are best thought of as weighted integrals of some quantity over a
finite volume.<br>
<br>
The RAW_STATE_1D_INTEGRAL has an associated half_width and
localization type (see the <a href=
"../../assimilation_code/modules/assimilation/cov_cutoff_mod.html">cov_cutoff
module</a> documentation) and a number of points at which to
compute the associated integral by quadrature. The location of the
observation defines the center of mass of the integral. The
integral is centered around the location and extends outward on
each side to 2*half_width. The weight associated with the integral
is defined by the weight of the localization function (for instance
Gaspari Cohn) using the same localization options as defined by the
cov_cutoff module. The number of points are used to equally divide
the range for computing the integral by quadrature.<br>
<br>
Special observation modules like <em class=
"file">obs_def_1d_state_mod.f90</em> contain Fortran 90 code
<em>and</em> additional specially formatted commented code that is
used to guide the preprocess program in constructing
obs_def_mod.f90 and obs_kind_mod.f90. The specially formatted
comments are most conveniently placed at the beginning of the
module and comprise seven sections, each beginning and ending with
a special F90 comment line that must be included
<em>verbatim</em>.<br>
<br>
The seven sections and their specific instances for the
1d_raw_state_mod are:
<ol>
<li>A list of all observation types defined by this module and
their associated generic quantities (see <a href=
"../../assimilation_code/programs/preprocess/preprocess.html">preprocess</a>
for details on quantity files). The header line is followed by
lines that have the observation type name (an all caps Fortran 90
identifier) and their associated generic quantity identifier. If
there is no special processing needed for an observation type, and
no additional data needed beyond the standard contents of an
observation then a third word on the line, <em class=
"unix">COMMON_CODE</em>, will instruct the preprocess program to
automatically generate all stubs and code needed for this type. For
observation types needing special code or additional data, this
word should not be specified and the user must supply the code
manually.
<pre>
! BEGIN DART PREPROCESS KIND LIST
! RAW_STATE_VARIABLE,    QTY_STATE_VARIABLE,   COMMON_CODE
! RAW_STATE_1D_INTEGRAL, QTY_1D_INTEGRAL
! END DART PREPROCESS KIND LIST
</pre>
<br>
<br></li>
<li>A list of all the use statements that the completed
obs_def_mod.f90 must have in order to use the public interfaces
provided by this special obs_def module. This section is optional
if there are no external interfaces.
<pre>
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral,  &amp;
!                                    interactive_1d_integral, get_expected_1d_integral, &amp;
!                                    set_1d_integral
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
</pre>
<br>
<br></li>
<li>Case statement entries for each observation type defined by
this special obs_def module stating how to compute the forward
observation operator. There must be a case statement entry for each
type of observation, <em>except</em> for observation types defined
with COMMON_CODE.
<pre>
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RAW_STATE_1D_INTEGRAL)
!            call get_expected_1d_integral(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
</pre>
<br>
<br></li>
<li>Case statement entries for each observation type defined by
this special obs_def module stating how to read any extra required
information from an obs sequence file. There must be a case
statement entry for each type of observation, <em>except</em> for
observation types defined with COMMON_CODE. If no special action is
required put a <code>continue</code> statement as the body of the
case instead of a subroutine call.
<pre>
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call read_1d_integral(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF
</pre>
<br>
<br></li>
<li>Case statement entries for each observation type defined by
this special obs_def module stating how to write any extra required
information to an obs sequence file. There must be a case statement
entry for each type of observation, <em>except</em> for observation
types defined with COMMON_CODE. If no special action is required
put a <code>continue</code> statement as the body of the case
instead of a subroutine call.
<pre>
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call write_1d_integral(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF
</pre>
<br>
<br></li>
<li>Case statement entries for each observation type defined by
this special obs_def module stating how to interactively create any
extra required information. There must be a case statement entry
for each type of observation, <em>except</em> for observation types
defined with COMMON_CODE. If no special action is required put a
<code>continue</code> statement as the body of the case instead of
a subroutine call.
<pre>
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call interactive_1d_integral(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
</pre>
<br>
<br></li>
<li>Any executable F90 module code must be tagged with the
following comments. All lines between these markers will be copied,
verbatim, to obs_def_mod.f90. This section is not required if there
are no observation-specific subroutines.
<pre>
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_1d_state_mod

... (module executable code)

end module obs_def_1d_state_mod
! END DART PREPROCESS MODULE CODE
</pre>
<br>
<br></li>
</ol>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
location_mod (1d_location_mod_only)
time_manager_mod
assim_model_mod
cov_cutoff_mod
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
<td><em class="call">use obs_def_mod, only :</em></td>
<td><a href="#write_1d_integral">write_1d_integral</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_1d_integral">read_1d_integral</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#interactive_1d_integral">interactive_1d_integral</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_expected_1d_integral">get_expected_1d_integral</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_1d_integral">set_1d_integral</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_power">write_power</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_power">read_power</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#interactive_power">interactive_power</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_expected_power">get_expected_power</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_power">set_power</a></td>
</tr>
</table>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
<a name="write_1d_integral" id="write_1d_integral"></a><br>
<div class="routine"><em class="call">call
write_1d_integral(igrkey, ifile, fform)</em>
<pre>
integer,          intent(in) :: <em class="code">igrkey</em>
integer,          intent(in) :: <em class="code">ifile</em>
character(len=*), intent(in) :: <em class="code">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes out the extra information for observation with unique
identifier key for a 1d_integral observation type. This includes
the half-width, localization type and number of quadrature points
for this observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">igrkey  </em></td>
<td>Unique integer key associated with the 1d integral observation
being processed. This is not the same as the key that all types of
observations have and uniquely distinguishes all observations from
each other; this is a key that is only set and retrieved by this
code for 1d integral observations. It is stored in the obs_def
derived type, not in the main obs_type definition.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile  </em></td>
<td>Unit number on which observation sequence file is open</td>
</tr>
<tr>
<td valign="top"><em class="code">fform  </em></td>
<td>String noting whether file is opened for 'formatted' or
'unformatted' IO.</td>
</tr>
</table>
</div>
<br>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
 <a name="read_1d_integral" id="read_1d_integral"></a><br>
<div class="routine"><em class="call">call read_1d_integral(igrkey,
ifile, fform)</em>
<pre>
integer,          intent(out) :: <em class="code">igrkey</em>
integer,          intent(in)  :: <em class="code">ifile</em>
character(len=*), intent(in)  :: <em class="code">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads the extra information for observation with unique
identifier key for a 1d_integral observation type. This information
includes the half-width, localization type and number of quadrature
points for this observation. The key that is returned is uniquely
associated with the definition that has been created and is used by
this module to keep track of the associated parameters for this
observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">igrkey  </em></td>
<td>Unique integer key associated with the observation being
processed.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile  </em></td>
<td>Unit number on which observation sequence file is open</td>
</tr>
<tr>
<td valign="top"><em class="code">fform  </em></td>
<td>String noting whether file is opened for 'formatted' or
'unformatted' IO.</td>
</tr>
</table>
</div>
<br>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
 <a name="interactive_1d_integral" id=
"interactive_1d_integral"></a><br>
<div class="routine"><em class="call">call
interactive_1d_integral(igrkey)</em>
<pre>
integer, intent(out) :: <em class="code">igrkey</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Uses input from standard in to define the characteristics of a
1D integral observation. The key that is returned is uniquely
associated with the definition that has been created and can be
used by this module to keep track of the associated parameters
(half_width, localization option, number of quadrature points) for
this key.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">igrkey  </em></td>
<td>Unique identifier associated with the created observation
definition in the obs sequence.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_expected_1d_integral" id=
"get_expected_1d_integral"></a><br>
<div class="routine"><em class="call">call
get_expected_1d_integral(state, location, igrkey, val,
istatus)</em>
<pre>
real(r8), intent(in)            :: <em class="code">state</em>
type(location_type), intent(in) :: <em class="code">location</em>
integer, intent(in)             :: <em class="code">igrkey</em>
real(r8), intent(out)           :: <em class="code">val</em>
integer, intent(out)            :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes the forward observation operator for a 1d integral
observation. Calls the <em class="code">interpolate()</em> routine
multiple times to invoke the forward operator code in whatever
model this has been compiled with.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state  </em></td>
<td>Model state vector (or extended state vector).</td>
</tr>
<tr>
<td valign="top"><em class="code">location  </em></td>
<td>Location of this observation.</td>
</tr>
<tr>
<td valign="top"><em class="code">igrkey  </em></td>
<td>Unique integer key associated with this observation.</td>
</tr>
<tr>
<td valign="top"><em class="code">val  </em></td>
<td>Returned value of forward observation operator.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus  </em></td>
<td>Returns 0 if forward operator was successfully computed, else
returns a positive value. (Negative values are reserved for system
use.)</td>
</tr>
</table>
</div>
<br>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
 <a name="set_1d_integral" id="set_1d_integral"></a><br>
<div class="routine"><em class="call">call
set_1d_integral(integral_half_width, num_eval_pts, localize_type,
igrkey, istatus)</em>
<pre>
real(r8), intent(in)  :: <em class="code">integral_half_width</em>
integer,  intent(in)  :: <em class="code">num_eval_pts</em>
integer,  intent(in)  :: <em class="code">localize_type</em>
integer,  intent(out) :: <em class="code">igrkey</em>
integer,  intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Available for use by programs that create observations to set
the additional metadata for these observation types. This
information includes the integral half-width, localization type and
number of quadrature points for this observation. The key that is
returned is uniquely associated with the definition that has been
created and should be set in the obs_def structure by calling
<em class="code">set_obs_def_key()</em>. This key is different from
the main observation key which all observation types have. This key
is unique to this observation type and is used when reading in the
observation sequence to match the corresponding metadata with each
observation of this type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">integral_half_width  </em></td>
<td>Real value setting the half-width of the integral.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">num_eval_pts  </em></td>
<td>Integer, number of evaluation points. 5-20 recommended.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">localize_type  </em></td>
<td>Integer localization type: 1=Gaspari-Cohn; 2=Boxcar; 3=Ramped
Boxcar</td>
</tr>
<tr>
<td valign="top"><em class="code">igrkey  </em></td>
<td>Unique integer key associated with the observation being
processed.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus  </em></td>
<td>Return code. 0 means success, any other value is an error</td>
</tr>
</table>
</div>
<br>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
 <a name="write_power" id="write_power"></a><br>
<div class="routine"><em class="call">call write_power(powkey,
ifile, fform)</em>
<pre>
integer,          intent(in) :: <em class="code">powkey</em>
integer,          intent(in) :: <em class="code">ifile</em>
character(len=*), intent(in) :: <em class="code">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes out the extra information, the power, for observation
with unique identifier key for a power observation type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">powkey  </em></td>
<td>Unique integer key associated with the power observation being
processed. This is not the same as the key that all types of
observations have and uniquely distinguishes all observations from
each other; this is a key that is only set and retrieved by this
code for power observations. It is stored in the obs_def derived
type, not in the main obs_type definition.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile  </em></td>
<td>Unit number on which observation sequence file is open</td>
</tr>
<tr>
<td valign="top"><em class="code">fform  </em></td>
<td>String noting whether file is opened for 'formatted' or
'unformatted' IO.</td>
</tr>
</table>
</div>
<br>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
 <a name="read_power" id="read_power"></a><br>
<div class="routine"><em class="call">call read_power(powkey,
ifile, fform)</em>
<pre>
integer,          intent(out) :: <em class="code">powkey</em>
integer,          intent(in)  :: <em class="code">ifile</em>
character(len=*), intent(in)  :: <em class="code">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads the extra information, the power, for observation with
unique identifier key for a power observation type. The key that is
returned is uniquely associated with the definition that has been
created and is used by this module to keep track of the associated
parameters for this observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">powkey  </em></td>
<td>Unique integer key associated with the observation being
processed.</td>
</tr>
<tr>
<td valign="top"><em class="code">ifile  </em></td>
<td>Unit number on which observation sequence file is open</td>
</tr>
<tr>
<td valign="top"><em class="code">fform  </em></td>
<td>String noting whether file is opened for 'formatted' or
'unformatted' IO.</td>
</tr>
</table>
</div>
<br>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
 <a name="interactive_power" id="interactive_power"></a><br>
<div class="routine"><em class="call">call
interactive_power(powkey)</em>
<pre>
integer, intent(out) :: <em class="code">powkey</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Uses input from standard in to define the characteristics of a
power observation. The key that is returned is uniquely associated
with the definition that has been created and can be used by this
module to keep track of the associated parameter, the power, for
this key.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">powkey  </em></td>
<td>Unique identifier associated with the created observation
definition in the obs sequence.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_expected_power" id="get_expected_power"></a><br>
<div class="routine"><em class="call">call
get_expected_power(state, location, powkey, val, istatus)</em>
<pre>
real(r8), intent(in)            :: <em class="code">state</em>
type(location_type), intent(in) :: <em class="code">location</em>
integer, intent(in)             :: <em class="code">powkey</em>
real(r8), intent(out)           :: <em class="code">val</em>
integer, intent(out)            :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes the forward observation operator for a power
observation. Calls the <em class="code">interpolate()</em> routine
to invoke the forward operator code in whatever model this has been
compiled with, then raises the result to the specified power
associated with this powkey.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state  </em></td>
<td>Model state vector (or extended state vector).</td>
</tr>
<tr>
<td valign="top"><em class="code">location  </em></td>
<td>Location of this observation.</td>
</tr>
<tr>
<td valign="top"><em class="code">powkey  </em></td>
<td>Unique integer key associated with this observation.</td>
</tr>
<tr>
<td valign="top"><em class="code">val  </em></td>
<td>Returned value of forward observation operator.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus  </em></td>
<td>Returns 0 if forward operator was successfully computed, else
returns a positive value. (Negative values are reserved for system
use.)</td>
</tr>
</table>
</div>
<br>
<!--============= DESCRIPTION OF A SUBROUTINE =======================-->
 <a name="set_power" id="set_power"></a><br>
<div class="routine"><em class="call">call set_power(power_in,
powkey, istatus)</em>
<pre>
real(r8), intent(in)  :: <em class="code">power_in</em>
integer,  intent(out) :: <em class="code">powkey</em>
integer,  intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Available for use by programs that create observations to set
the additional metadata for these observation types. This
information includes the power to which to raise the state
variable. The key that is returned is uniquely associated with the
definition that has been created and should be set in the obs_def
structure by calling <em class="code">set_obs_def_key()</em>. This
key is different from the main observation key which all
observation types have. This key is unique to this observation type
and is used when reading in the observation sequence to match the
corresponding metadata with each observation of this type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">power_in  </em></td>
<td>Real value setting the power.</td>
</tr>
<tr>
<td valign="top"><em class="code">powkey  </em></td>
<td>Unique integer key associated with the observation being
processed.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus  </em></td>
<td>Return code. 0 means success, any other value is an error</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Namelist for this module.                           -->
<!--==================================================================-->
 <a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This module has no namelist.</p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>NONE</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ol>
<li>none</li>
</ol>
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
<td valign="top">interactive_1d_integral</td>
<!-- message -->
<td valign="top">Out of space, max_1d_integral_obs limit NNNN
(currently 1000).</td>
<!-- comment -->
<td valign="top">There is only room for a fixed number of 1d
integral observations. The max number is defined by
max_1d_integral_obs. Set this to a larger value if more are
needed.</td>
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
