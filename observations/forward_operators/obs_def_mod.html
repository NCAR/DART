<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_def_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE obs_def_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Interface">INTERFACES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>The DART Fortran90 derived type <em class="unix">obs_def</em>
provide an abstraction of the definition of an observation. An
observation sequence <em class="unix">obs_seq</em> at a higher
level is composed of observation definitions associated with
observed values. For now, the basic operations required to
implement an observation definition are an ability to compute a
forward operator given the model state vector, the ability to
read/write the observation definition from/to a file, and a
capability to do a standard input driven interactive definition of
the observation definition.<br>
<br>
DART makes a distinction between specific <em class=
"unix">observation types</em> and generic <em class=
"unix">observation quantities</em>. The role of the various obs_def
input files is to define the mapping between the types and
quantities, and optionally to provide type-specific processing
routines.<br>
<br>
A single obs_def output module is created by the program <em class=
"unix">preprocess</em> from two kinds of input files. First, a
DEFAULT obs_def module (normally called <em class=
"file">DEFAULT_obs_def_mod.F90</em> and documented in this
directory) is used as a template into which the preprocessor
incorporates information from zero or more special obs_def modules
(such as <em class="file">obs_def_1d_state_mod.f90</em> or
<em class="file">obs_def_reanalysis_bufr_mod.f90</em>, also
documented in this directory). If no special obs_def files are
included in the preprocessor namelist, a minimal <em class=
"file">obs_def_mod.f90</em> is created which can only support
identity forward observation operators.<br>
<br>
To add a new observation type which does not fit into any of the
already-defined obs_def files, a new file should be created in the
<em class="file">obs_def</em> directory. These files are usually
named according the the pattern <em class=
"file">obs_def_</em>X<em class="file">_mod.f90</em>, where the X is
either an instrument name, a data source, or a class of
observations. See the existing filenames in that directory for
ideas. Then this new filename must be listed in the <em class=
"file">input.nml</em> namelist for the model, in the <em class=
"unix">&amp;preprocess_nml</em> section, in the <em class=
"unix">obs_type_files</em> variable. This variable is a string list
type which can contain multiple filenames. Running the <em class=
"unix">preprocess</em> program will then use the contents of the
new file to generate the needed output files for use in linking to
the rest of the DART system.</p>
<h3>Simple observations</h3>
<p>If the new observation type can be directly interpolated by a
model_mod interpolation routine, and has no additional
observation-specific code for reading, writing, or initializing the
observation, then the entire contents of the new file is:</p>
<pre>
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! <em>type</em>, <em>quantity</em>, COMMON_CODE
! <em>(repeat lines for each type)</em>
! END DART PREPROCESS TYPE DEFINITIONS
</pre>
<p>DART will automatically generate all interface code needed for
these new observation types. For example, here is a real list:</p>
<pre>
! BEGIN DART PREPROCESS TYPE DEFINITIONS
!VELOCITY,                     QTY_VELOCITY,              COMMON_CODE
!TRACER_CONCENTRATION,         QTY_TRACER_CONCENTRATION,  COMMON_CODE
!TRACER_SOURCE,                QTY_TRACER_SOURCE,         COMMON_CODE
!MEAN_SOURCE,                  QTY_MEAN_SOURCE,           COMMON_CODE
!SOURCE_PHASE,                 QTY_SOURCE_PHASE,          COMMON_CODE
! END DART PREPROCESS TYPE DEFINITIONS
</pre>
<p>The first column is the specific observation <em class=
"unix">type</em> and should be unique. The second column is the
generic observation <em class="unix">quantity</em>. The quantities
available to DART are defined at compile time by
<em>preprocess</em> via the option 'quantity_files' in the
<em>preprocess_nml</em> namelist. The third column must be the
keyword <em class="unix">COMMON_CODE</em> which tells the
<em class="unix">preprocess</em> program to automatically generate
all necessary interface code for this type.</p>
<h3>Observations needing special handling</h3>
<p>For observation types which have observation-specific routines,
must interpolate using a combination of other generic quantities,
or require additional observation-specific data to be stored, the
following format is used:</p>
<pre>
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! <em>type</em>, <em>quantity</em>
! <em>(repeat lines for each type/quantity pair)</em>
! END DART PREPROCESS TYPE DEFINITIONS
</pre>
<p>DART will need user-supplied interface code for each of the
listed types. For example, here is a real list:</p>
<pre>
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! DOPPLER_RADIAL_VELOCITY, QTY_VELOCITY
! RADAR_REFLECTIVITY,      QTY_RADAR_REFLECTIVITY
! END DART PREPROCESS TYPE DEFINITIONS
</pre>
<p>In this case, DART needs additional information for how to
process these types. They include code sections delimited by
precisely formatted comments, and possibly module code
sections:</p>
<ol>
<li>
<pre>
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
</pre>
Any fortran use statements for public subroutines or variables from
other modules should be placed between these lines, with comment
characters in the first column.<br>
<br>
For example, if the forward operator code includes a module with
public routines then a "use" statement like:
<pre>
use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &amp;
                                 interactive_1d_integral, get_expected_1d_integral
</pre>
needs to be added to the obs_def_mod so the listed subroutines are
available to be called. This would look like:
<pre>
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
! use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &amp;
!                                  interactive_1d_integral, get_expected_1d_integral
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
</pre></li>
<li>
<pre>
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
</pre>
These comments must enclose a case statement for each defined type
that returns the expected observation value based on the current
values of the state vector. The code must be in comments, with the
comment character in the first column.<br>
<br>
The variables available to be passed to subroutines or used in this
section of code are:
<table>
<tr>
<td><em class="code">state</em></td>
<td>the entire model state vector</td>
</tr>
<tr>
<td><em class="code">state_time</em></td>
<td>the time of the state data</td>
</tr>
<tr>
<td><em class="code">ens_index</em></td>
<td>the ensemble member number</td>
</tr>
<tr>
<td><em class="code">location</em></td>
<td>the observation location</td>
</tr>
<tr>
<td><em class="code">obs_kind_ind </em></td>
<td>the index of the specific observation type</td>
</tr>
<tr>
<td><em class="code">obs_time</em></td>
<td>the time of the observation</td>
</tr>
<tr>
<td><em class="code">error_val</em></td>
<td>the observation error variance</td>
</tr>
</table>
<br>
The routine must fill in the values of these variables:
<table>
<tr>
<td><em class="code">obs_val </em></td>
<td>the computed forward operator value</td>
</tr>
<tr>
<td><em class="code">istatus</em></td>
<td>return code: 0=ok, &gt;0 is error, &lt;0 reserved for system
use</td>
</tr>
</table>
<br>
To call a model_mod interpolate routine directly, the argument list
must match exactly:
<pre>
interpolate(state, location, QTY_xxx, obs_val, istatus)
</pre>
This can be useful if the forward operator needs to retrieve values
for fields which are typically found in a model and then compute a
derived value from them.</li>
<li>
<pre>
! BEGIN DART PREPROCESS READ_OBS_DEF
! END DART PREPROCESS READ_OBS_DEF
</pre>
These comments must enclose a case statement for each defined type
that reads any additional data associated with a single
observation. If there is no information beyond that for the basic
obs_def type, the case statement must still be provided, but the
code can simply be <em class="unix">continue</em>. The code must be
in comments, with the comment character in the first column.<br>
<br>
The variables available to be passed to subroutines or used in this
section of code are:
<table>
<tr>
<td><em class="code">ifile</em></td>
<td>the open unit number positioned ready to read, read-only</td>
</tr>
<tr>
<td><em class="code">obs_def </em></td>
<td>the rest of the obs_def derived type for this obs,
read-write</td>
</tr>
<tr>
<td><em class="code">key</em></td>
<td>the index observation number in this sequence, read-only</td>
</tr>
<tr>
<td><em class="code">obs_val</em></td>
<td>the observation value, if needed. in general should not be
changed</td>
</tr>
<tr>
<td><em class="code">is_ascii</em></td>
<td>logical to indicate how the file was opened, formatted or
unformatted</td>
</tr>
</table>
<br>
The usual use of this routine is to read in additional metadata per
observation and to set the private key in the <em class=
"code">obs_def</em> to indicate which index to use for this
observation to look up the corresponding metadata in arrays or
derived types. Do not confuse the key in the obs_def with the key
argument to this routine; the latter is the global observation
sequence number for this observation.</li>
<li>
<pre>
! BEGIN DART PREPROCESS WRITE_OBS_DEF
! END DART PREPROCESS WRITE_OBS_DEF
</pre>
These comments must enclose a case statement for each defined type
that writes any additional data associated with a single
observation. If there is no information beyond that for the basic
obs_def type, the case statement must still be provided, but the
code can simply be <em class="unix">continue</em>. The code must be
in comments, with the comment character in the first column.<br>
<br>
The variables available to be passed to subroutines or used in this
section of code are:
<table>
<tr>
<td><em class="code">ifile</em></td>
<td>the open unit number positioned ready to write, read-only</td>
</tr>
<tr>
<td><em class="code">obs_def </em></td>
<td>the rest of the obs_def derived type for this obs,
read-only</td>
</tr>
<tr>
<td><em class="code">key</em></td>
<td>the index observation number in this sequence, read-only</td>
</tr>
<tr>
<td><em class="code">is_ascii</em></td>
<td>logical to indicate how the file was opened, formatted or
unformatted</td>
</tr>
</table>
<br>
The usual use of this routine is to write the additional metadata
for this observation based on the private key in the <em class=
"code">obs_def</em>. Do not confuse this with the key in the
subroutine call which is the observation number relative to the
entire observation sequence file.</li>
<li>
<pre>
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! END DART PREPROCESS INTERACTIVE_OBS_DEF
</pre>
These comments must enclose a case statement for each defined type
that prompts the user for any additional data associated with a
single observation. If there is no information beyond that for the
basic obs_def type, the case statement must still be provided, but
the code can simply be <em class="unix">continue</em>. The code
must be in comments, with the comment character in the first
column.<br>
<br>
The variables available to be passed to subroutines or used in this
section of code are:
<table>
<tr>
<td><em class="code">obs_def </em></td>
<td>the rest of the obs_def derived type for this obs,
read-write</td>
</tr>
<tr>
<td><em class="code">key</em></td>
<td>the index observation number in this sequence, read-only</td>
</tr>
</table>
<br>
The DART code will prompt for the rest of the obs_def values
(location, type, value, error) but any additional metadata needed
by this observation type should be prompted to, and read from, the
console (e.g. <em class="code">write(*,*)</em>, and <em class=
"code">read(*, *)</em>). The code will generally set the <em class=
"code">obs_def%key</em> value as part of setting the metadata.</li>
<li>
<pre>
! BEGIN DART PREPROCESS MODULE CODE
! END DART PREPROCESS MODULE CODE
</pre>
If the code to process this observation requires module data and/or
subroutines, then these comments must surround the module
definitions. Unlike all the other sections, this comment pair is
optional, and if used, the code must not be in comments; it will be
copied verbatim over to the output file.<br>
<br>
Generally the code for a forward operator should be defined inside
a module, to keep module variables and other private subroutines
from colliding with unrelated routines and variables in other
forward operator files.</li>
</ol>
<p>It is possible to mix automatic code types and user-supplied
code types in the same list. Simply add the COMMON_CODE keyword on
the lines which need no special data or interfaces. For example,
here is an extract from the 1d state obs_def module, where the raw
state variable needs only autogenerated code, but the 1d integral
has user-supplied processing code:</p>
<pre>
! BEGIN DART PREPROCESS TYPE LIST
! RAW_STATE_VARIABLE,    QTY_STATE_VARIABLE, COMMON_CODE
! RAW_STATE_1D_INTEGRAL, QTY_1D_INTEGRAL
! END DART PREPROCESS TYPE LIST


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &amp;
!                                    interactive_1d_integral, get_expected_1d_integral
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(RAW_STATE_1D_INTEGRAL)
!            call get_expected_1d_integral(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call read_1d_integral(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call write_1d_integral(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(RAW_STATE_1D_INTEGRAL)
!         call interactive_1d_integral(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_1d_state_mod

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location
use  assim_model_mod, only : interpolate
use   cov_cutoff_mod, only : comp_cov_factor

implicit none

public :: write_1d_integral, read_1d_integral, interactive_1d_integral, &amp;
          get_expected_1d_integral

...  (module code here)

end module obs_def_1d_state_mod
! END DART PREPROCESS MODULE CODE
</pre>
<p>See the <a href=
"obs_def_1d_state_mod.html">obs_def_1d_state_mod</a> documentation
for more details and examples of each section. Also see <a href=
"obs_def_wind_speed_mod.f90">obs_def_wind_speed_mod.f90</a> for an
example of a 3D geophysical forward operator.<br>
<br>
In addition to collecting and managing any additional observation
type-specific code, this module provides the definition of the
obs_def_type derived type, and a collection of subroutines for
creating, accessing, and updating this type. The remainder of this
document describes the subroutines provided by this module.</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
location_mod (depends on model choice)
time_manager_mod
assim_model_mod
obs_kind_mod
Other special obs_def_kind modules as required
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
<td><a href="#obs_def_type">obs_def_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_obs_def">init_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_obs_def_location">get_obs_def_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_obs_def_type_of_obs">get_obs_def_type_of_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_obs_def_time">get_obs_def_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_obs_def_error_variance">get_obs_def_error_variance</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_obs_def_key">get_obs_def_key</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_obs_def_location">set_obs_def_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#set_obs_def_type_of_obs">set_obs_def_type_of_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_obs_def_time">set_obs_def_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#set_obs_def_error_variance">set_obs_def_error_variance</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_obs_def_key">set_obs_def_key</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#interactive_obs_def">interactive_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_obs_def">write_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_obs_def">read_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_expected_obs_from_def">get_expected_obs_from_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#destroy_obs_def">destroy_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#copy_obs_def">copy_obs_def</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#copy_obs_def">assignment(=)</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_name_for_type_of_obs">get_name_for_type_of_obs</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--=================== DESCRIPTION OF A LOCAL TYPE ==================-->
<a name="obs_def_type" id="obs_def_type"></a><br>
<div class="routine">
<pre>
<em class="call">type obs_def_type</em>
   private
   type(location_type)  :: location
   integer              :: kind
   type(time_type)      :: time
   real(r8)             :: error_variance
   integer              :: key
end type obs_def_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Models all that is known about an observation except for actual
values. Includes a location, type, time and error variance.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">location</td>
<td>Location of the observation.</td>
</tr>
<tr>
<td valign="top">kind</td>
<td>Despite the name, the specific type of the observation.</td>
</tr>
<tr>
<td valign="top">time</td>
<td>Time of the observation.</td>
</tr>
<tr>
<td valign="top">error_variance</td>
<td>Error variance of the observation.</td>
</tr>
<tr>
<td valign="top">key</td>
<td>Unique identifier for observations of a particular type.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_obs_def" id="init_obs_def"></a><br>
<div class="routine"><em class="call">call init_obs_def(obs_def,
location, kind, time, error_variance)</em>
<pre>
type(obs_def_type),  intent(out) :: <em class="code">obs_def</em>
type(location_type), intent(in)  :: <em class="code">location</em>
integer,             intent(in)  :: <em class="code">kind</em>
type(time_type),     intent(in)  :: <em class="code">time</em>
real(r8),            intent(in)  :: <em class=
"code">error_variance</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Creates an obs_def type with location, type, time and
error_variance specified.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>The obs_def that is created</td>
</tr>
<tr>
<td valign="top"><em class="code">location  </em></td>
<td>Location for this obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">kind  </em></td>
<td>Observation type for obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">time  </em></td>
<td>Time for obs_def</td>
</tr>
<tr>
<td valign="top"><em class=
"code">error_variance  </em></td>
<td>Error variance of this observation</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="copy_obs_def" id="copy_obs_def"></a><br>
<div class="routine"><em class="call">call copy_obs_def(obs_def1,
obs_def2)</em>
<pre>
type(obs_def_type), intent(out) :: <em class="code">obs_def1</em>
type(obs_def_type), intent(in)  :: <em class="code">obs_def2</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Copies obs_def2 to obs_def1, overloaded as assignment (=).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def1  </em></td>
<td>obs_def to be copied into</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def2  </em></td>
<td>obs_def to be copied from</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_obs_def_key" id="get_obs_def_key"></a><br>
<div class="routine"><em class="call">var =
get_obs_def_key(obs_def)</em>
<pre>
integer                        :: <em class=
"code">get_obs_def_key</em>
type(obs_def_type), intent(in) :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns key from an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns key from an obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_obs_def_error_variance" id=
"get_obs_def_error_variance"></a><br>
<div class="routine"><em class="call">var =
get_obs_def_error_variance(obs_def)</em>
<pre>
real(r8)                       :: <em class=
"code">get_obs_def_error_variance</em>
type(obs_def_type), intent(in) :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns error variance from an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Error variance from an obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_obs_def_location" id="get_obs_def_location"></a><br>
<div class="routine"><em class="call">var =
get_obs_def_location(obs_def)</em>
<pre>
type(location_type)              :: <em class=
"code">get_obs_def_location</em>
type(obs_def_type), intent(in)   :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the location from an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns location from an obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_obs_def_type_of_obs" id=
"get_obs_def_type_of_obs"></a><br>
<div class="routine"><em class="call">var =
get_obs_def_type_of_obs(obs_def)</em>
<pre>
integer                         :: <em class=
"code">get_obs_def_type_of_obs</em>
type(obs_def_type),  intent(in) :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns an observation type from an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns the observation type from an obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_obs_def_time" id="get_obs_def_time"></a><br>
<div class="routine"><em class="call">var =
get_obs_def_time(obs_def)</em>
<pre>
type(time_type)                :: <em class=
"code">get_obs_def_time</em>
type(obs_def_type), intent(in) :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns time from an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns time from an obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_name_for_type_of_obs" id=
"get_name_for_type_of_obs"></a><br>
<div class="routine"><em class="call">obs_name =
get_name_for_type_of_obs(obs_kind_ind)</em>
<pre>
character(len = 32)            :: <em class=
"code">get_name_for_type_of_obs</em>
integer, intent(in)            :: <em class=
"code">obs_kind_ind</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns an observation name from an observation type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var  </em></td>
<td>Returns name from an observation type</td>
</tr>
<tr>
<td valign="top"><em class=
"code">obs_kind_ind  </em></td>
<td>An observation type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_obs_def_location" id="set_obs_def_location"></a><br>
<div class="routine"><em class="call">call
set_obs_def_location(obs_def, location)</em>
<pre>
type(obs_def_type),  intent(inout) :: <em class="code">obs_def</em>
type(location_type), intent(in)    :: <em class=
"code">location</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the location in an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">location  </em></td>
<td>A location</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_obs_def_error_variance" id=
"set_obs_def_error_variance"></a><br>
<div class="routine"><em class="call">call
set_obs_def_error_variance(obs_def, error_variance)</em>
<pre>
type(obs_def_type), intent(inout) :: <em class="code">obs_def</em>
real(r8), intent(in)              :: <em class=
"code">error_variance</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set error variance for an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
<tr>
<td valign="top"><em class=
"code">error_variance  </em></td>
<td>Error variance</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_obs_def_key" id="set_obs_def_key"></a><br>
<div class="routine"><em class="call">call set_obs_def_key(obs_def,
key)</em>
<pre>
type(obs_def_type), intent(inout) :: <em class="code">obs_def</em>
integer,            intent(in)    :: <em class="code">key</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the key for an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>Unique identifier for this observation</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_obs_def_type_of_obs" id=
"set_obs_def_type_of_obs"></a><br>
<div class="routine"><em class="call">call
set_obs_def_type_of_obs(obs_def, kind)</em>
<pre>
type(obs_def_type), intent(inout) :: <em class="code">obs_def</em>
integer,            intent(in)    :: <em class="code">kind</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the type of observation in an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">kind  </em></td>
<td>An integer observation type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_obs_def_time" id="set_obs_def_time"></a><br>
<div class="routine"><em class="call">call
set_obs_def_time(obs_def, time)</em>
<pre>
type(obs_def_type), intent(inout) :: <em class="code">obs_def</em>
type(time_type), intent(in)       :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets time for an observation definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def</td>
</tr>
<tr>
<td valign="top"><em class="code">time  </em></td>
<td>Time to set</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_expected_obs_from_def" id=
"get_expected_obs_from_def"></a><br>
<div class="routine"><em class="call">call
get_expected_obs_from_def(key, obs_def, obs_kind_ind, ens_index,
state, state_time, obs_val, istatus, assimilate_this_ob,
evaluate_this_ob)</em>
<pre>
integer,            intent(in)  :: <em class="code">key</em>
type(obs_def_type), intent(in)  :: <em class="code">obs_def</em>
integer,            intent(in)  :: <em class=
"code">obs_kind_ind</em>
integer,            intent(in)  :: <em class="code">ens_index</em>
real(r8),           intent(in)  :: <em class="code">state(:)</em>
type(time_type),    intent(in)  :: <em class="code">state_time</em>
real(r8),           intent(out) :: <em class="code">obs_val</em>
integer,            intent(out) :: <em class="code">istatus</em>
logical,            intent(out) :: <em class=
"code">assimilate_this_ob</em>
logical,            intent(out) :: <em class=
"code">evaluate_this_ob</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Compute the observation (forward) operator for a particular obs
definition.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>descriptor for observation type</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>The input obs_def</td>
</tr>
<tr>
<td valign="top"><em class=
"code">obs_kind_ind  </em></td>
<td>The obs type</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_index  </em></td>
<td>The ensemble member number of this state vector</td>
</tr>
<tr>
<td valign="top"><em class="code">state  </em></td>
<td>Model state vector</td>
</tr>
<tr>
<td valign="top"><em class="code">state_time  </em></td>
<td>Time of the data in the model state vector</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus  </em></td>
<td>Returned integer describing problems with applying forward
operator (0 == OK, &gt;0 == error, &lt;0 reserved for sys
use).</td>
</tr>
<tr>
<td valign="top"><em class=
"code">assimilate_this_ob  </em></td>
<td>Indicates whether to assimilate this obs or not</td>
</tr>
<tr>
<td valign="top"><em class=
"code">evaluate_this_ob  </em></td>
<td>Indicates whether to evaluate this obs or not</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_obs_def" id="read_obs_def"></a><br>
<div class="routine"><em class="call">call read_obs_def(ifile,
obs_def, key, obs_val <em class="optionalcode">[,fform]</em>)</em>
<pre>
integer,                    intent(in)    :: <em class=
"code">ifile</em>
type(obs_def_type),         intent(inout) :: <em class=
"code">obs_def</em>
integer,                    intent(in)    :: <em class=
"code">key</em>
real(r8),                   intent(inout) :: <em class=
"code">obs_val</em>
character(len=*), optional, intent(in)    :: <em class=
"code">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads an obs_def from file open on channel ifile. Uses format
specified in fform or FORMATTED if fform is not present.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ifile  </em></td>
<td>File unit open to output file</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>Observation definition to be read</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>Present if unique identifier key is needed by some obs type.
Unused by default code.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_val  </em></td>
<td>Present if needed to perform operations based on value. Unused
by default code.</td>
</tr>
<tr>
<td valign="top"><em class="code">fform  </em></td>
<td>File format specifier: FORMATTED or UNFORMATTED; default
FORMATTED (FORMATTED in this case is the human readable/text option
as opposed to UNFORMATTED which is binary.)</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="interactive_obs_def" id="interactive_obs_def"></a><br>
<div class="routine"><em class="call">call
interactive_obs_def(obs_def, key)</em>
<pre>
type(obs_def_type), intent(inout) :: <em class="code">obs_def</em>
integer,            intent(in)    :: <em class="code">key</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Creates an obs_def via input from standard in.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def to be created</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>Present if unique identifier key is needed by some obs type.
Unused by default code.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_obs_def" id="write_obs_def"></a><br>
<div class="routine"><em class="call">call write_obs_def(ifile,
obs_def, key <em class="optionalcode">[,fform]</em>)</em>
<pre>
integer,                    intent(in) :: <em class=
"code">ifile</em>
type(obs_def_type),         intent(in) :: <em class=
"code">obs_def</em>
integer,                    intent(in) :: <em class="code">key</em>
character(len=*), optional, intent(in) :: <em class=
"code">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes an obs_def to file open on channel ifile. Uses format
specified in fform or FORMATTED if fform is not present.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ifile  </em></td>
<td>File unit open to output file</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>Observation definition to be written</td>
</tr>
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>Present if unique identifier key is needed by some obs type.
Unused by default code.</td>
</tr>
<tr>
<td valign="top"><em class="code">fform  </em></td>
<td>File format specifier: FORMATTED or UNFORMATTED; default
FORMATTED</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="destroy_obs_def" id="destroy_obs_def"></a><br>
<div class="routine"><em class="call">call
destroy_obs_def(obs_def)</em>
<pre>
type(obs_def_type), intent(inout) :: <em class="code">obs_def</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Releases all storage associated with an obs_def and its
subcomponents.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">obs_def  </em></td>
<td>An obs_def to be released.</td>
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
<li>The read_obs_def() and write_obs_def() routines are passed an
already-opened file channel/descriptor and read to or write from
it.</li>
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
<td valign="top">get_expected_obs_from_def</td>
<!-- message -->
<td valign="top">Attempt to evaluate undefined observation
type</td>
<!-- comment -->
<td valign="top">An observation type for which no forward operator
has been defined is an error.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_obs_def</td>
<!-- message -->
<td valign="top">Expected header "obdef" in input file</td>
<!-- comment -->
<td valign="top">The format of the input file is not
consistent.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_obs_def</td>
<!-- message -->
<td valign="top">Expected kind header "kind " in input file</td>
<!-- comment -->
<td valign="top">The format of the input file is not
consistent.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_obs_def</td>
<!-- message -->
<td valign="top">Attempt to read for undefined obs_kind index</td>
<!-- comment -->
<td valign="top">Reading for an observation type for which no
forward operator has been defined is an error.</td>
</tr>
<tr><!-- routine -->
<td valign="top">write_obs_def</td>
<!-- message -->
<td valign="top">Attempt to write for undefined obs_kind index</td>
<!-- comment -->
<td valign="top">Writing for an observation type for which no
forward operator has been defined is an error.</td>
</tr>
<tr><!-- routine -->
<td valign="top">interactive_obs_def</td>
<!-- message -->
<td valign="top">Attempt to interactively create undefined obs_kind
index</td>
<!-- comment -->
<td valign="top">Creating an observation type for which no forward
operator has been defined is an error.</td>
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
