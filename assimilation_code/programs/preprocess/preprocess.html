<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program preprocess</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM preprocess</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>Preprocess is a DART-supplied preprocessor program. Preprocess
is used to insert observation specific code into DART at compile
time.</p>
<p>In DART, forward operators are not specific to any one model. To
achieve this separation between models and forward operators DART
makes a distinction between an observation <em>type</em> and a
physical <em>quantity</em>. For example, a radiosonde used to
measure windspeed would be a <em>type</em> of observation. Zonal
wind and meridional wind are <em>quantities</em> used to calculate
windspeed. Specifying many observation types allows DART to be able
to evaluate some observations and assimilate others even if the
instruments measure the same quantity.</p>
<p>Preprocess takes user supplied observation and quantity files
and combines them with template files to produce code for DART. Use
the namelist option 'obs_type_files' to specify the input
observation files and the namelist option 'quantity_files' to
specify the input quantity files.</p>
<ul>
<li>If no quantity files are given, a default list of quantities is
used.</li>
<li>If no obs_type_files are given, only identity observations can
be used in the filter (i.e. the state variable values are directly
observed; forward operator is an identity)</li>
</ul>
<p>The template files <em class="file">DEFAULT_obs_def_mod.F90</em>
and <em class="file">DEFAULT_obs_kind_mod.F90</em> contain
specially formatted comment lines. These comment lines are used as
markers to insert observation specific information. Prepreocess
relies these comment lines being used <em>verbatim</em>.</p>
<p>There is no need to to alter <em class=
"file">DEFAULT_obs_def_mod.F90</em> or <em class=
"file">DEFAULT_obs_kind_mod.F90</em>. Detailed instructions for
adding new observation types can be found in <a href=
"../../../observations/forward_operators/obs_def_mod.html">obs_def_mod.html</a>.
New quantities should be added to a quantity file, for example a new atmosphere
quantity should be added to <em class="file">atmosphere_quantities_mod.f90</em>.</p>

<p>Every line in a quantity file between the start and end markers
must be a comment or a quantity definition (QTY_string). Multiple
name-value pairs can be specified for a quantity but are not
required. For example, temperature may be defined: <code>!
QTY_TEMPERATURE units="K" minval=0.0</code>. Comments are allowed
between quantity definitions or on the same line as the definition.
The code snippet below shows acceptable formats for quantity
definitions</p>
<p><code>! BEGIN DART PREPROCESS QUANTITY DEFINITIONS<br>
!<br>
! Formats accepted:<br>
!<br>
! QTY_string<br>
! QTY_string name=value<br>
! QTY_string name=value name2=value2<br>
!<br>
! QTY_string ! comments<br>
!<br>
! ! comment<br>
!<br>
! END DART PREPROCESS QUANTITY DEFINITIONS<br></code></p>
<p>The output files produced by preprocess are named <em class=
"file">assimilation_code/modules/observations/obs_kind_mod.f90</em>
and <em class=
"file">observations/forward_operators/obs_def_mod.f90</em>, but can
be renamed by namelist control if needed. Be aware that if you
change the name of these output files, you will need to change the
path_names files for DART executables.<br>
<br>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
 <a name="Namelist" id="Namelist"></a></p>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>When you run preprocess, the namelist is read from the file
<em class="file">input.nml</em> in the directory where preprocess
is run.</p>
<p>Namelists start with an ampersand '&amp;' and terminate with a
slash '/'. Character strings that contain a '/' must be enclosed in
quotes to prevent them from prematurely terminating the
namelist.</p>
<div class="namelist">
<pre>
&amp;preprocess_nml
  overwrite_output        = .true.,
  input_obs_def_mod_file  = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90',
  output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90',
  input_obs_qty_mod_file  = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90',
  output_obs_qty_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90',
  quantity_files          = '../../../assimilation_code/modules/observations/atmosphere_quantities_mod.f90',
  obs_type_files          = '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90',
                            '../../../observations/forward_operators/obs_def_rel_humidity_mod.f90',
                            '../../../observations/forward_operators/obs_def_altimeter_mod.f90'
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
<td>input_obs_def_mod_file</td>
<td>character(len=256)<br></td>
<td>Path name of the template obs def module to be preprocessed.
The default is <em class=
"file">../../../observations/forward_operators/DEFAULT_obs_def_mod.F90</em>.
This file must have the appropriate commented lines indicating
where the different parts of the input special obs definition
modules are to be inserted.</td>
</tr>
<tr>
<td>output_obs_def_mod_file</td>
<td>character(len=256)<br></td>
<td>Path name of output obs def module to be created by preprocess.
The default is <em class=
"file">../../../observations/forward_operators/obs_def_mod.f90</em>.</td>
</tr>
<tr>
<td>input_obs_qty_mod_file</td>
<td>character(len=256)<br></td>
<td>Path name of input obs quantity file to be preprocessed. The
default path name is <em class=
"file">../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90</em>.
This file must have the appropriate commented lines indicating
where the different quantity modules are to be inserted.</td>
</tr>
<tr>
<td>output_obs_qty_mod_file</td>
<td>character(len=256)<br></td>
<td>Path name of output obs quantity module to be created by
preprocess. The default is <em class=
"file">../../../assimilation_code/modules/observations/obs_kind_mod.f90</em>.</td>
</tr>
<tr>
<td>obs_type_files</td>
<td>character(len=256)(:)<br></td>
<td>A list of files containing observation definitions for the type
of observations you want to use with DART. The maximum number of
files is limited to MAX_OBS_TYPE_FILES = 1000. The DART obs_def
files are in <em class=
"file">observations/forward_operators/obs_def_*.mod.f90</em>.</td>
</tr>
<tr>
<td>overwrite_output</td>
<td>logical<br></td>
<td>By defualt, preprocess will overwrite the existing
obs_kind_mod.f90 and obs_def_mod.f90 files. Set overwrite_output =
.false. if you want to preprocess to not overwrite existing
files.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
parse_arges_mod
types_mod
utilities_mod
</pre>
<p>Namelist interface <em class="code">&amp;preprocess_nml</em>
must be read from file <em class="file">input.nml</em>.</p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>input_obs_def_mod_file, specified by namelist; usually
<em class="file">DEFAULT_obs_def_mod.F90</em>.</li>
<li>output_obs_def_mod_file, specified by namelist; usually
<em class="file">obs_def_mod.f90</em>.</li>
<li>input_obs_qty_mod_file, specified by namelist; usually
<em class="file">DEFAULT_obs_kind_mod.F90</em>.</li>
<li>output_obs_qty_mod_file, specified by namelist; usually
<em class="file">obs_kind_mod.f90</em>.</li>
<li>obs_type_files, specified by namelist; usually files like
<em class="file">obs_def_reanalysis_bufr_mod.f90</em>.</li>
<li>quantity_files, specified by namelist; usually files like
<em class="file">atmosphere_quantities_mod.f90</em>.</li>
<li>namelistfile</li>
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
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">file ____ does not exist (and must)</td>
<!-- comment -->
<td valign="top">The input obs_type_files and qty_files must
exist.</td>
</tr>
<tr><!-- routine -->
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">file _____ does NOT contain ! BEGIN DART
PREPROCESS QUANTITY LIST</td>
<!-- comment -->
<td valign="top">Each special obs_def input file must contain this
comment string.</td>
</tr>
<tr><!-- routine -->
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">file _____ does NOT contain " END DART PREPROCESS
KIND LIST</td>
<!-- comment -->
<td valign="top">Each special obs_def input file must contain this
comment string.</td>
</tr>
<tr><!-- routine -->
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">Input DEFAULT obs_kind file ended
unexpectedly.</td>
<!-- comment -->
<td valign="top">Did not find strings indicating where to insert
special obs_def sections in the input obs_kind module.</td>
</tr>
<tr><!-- routine -->
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">Input DEFAULT obs_def file ended
unexpectedly.</td>
<!-- comment -->
<td valign="top">Did not find strings indicating where to insert
special obs_def sections in the input obs_def module.</td>
</tr>
<tr><!-- routine -->
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">file _____ does NOT contain ! BEGIN DART
PREPROCESS.</td>
<!-- comment -->
<td valign="top">Input special obs_def file must contain this
comment string.</td>
</tr>
<tr><!-- routine -->
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">file _____ does NOT contain ! END DART
PREPROCESS.</td>
<!-- comment -->
<td valign="top">Input special obs_def file must contain this
comment string.</td>
</tr>
<tr><!-- routine -->
<td valign="top">preprocess</td>
<!-- message -->
<td valign="top">'Incompatible duplicate entry detected'</td>
<!-- comment -->
<td valign="top">A quantity has been defined more than once but
with conflicting metadata in each definition.</td>
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
<p>none</p>
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
