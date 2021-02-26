<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module model_mod (CLM)</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>CLM</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>This is the DART interface to the Community Land Model (CLM). It
is run as part of the <a href=
"http://www.cesm.ucar.edu/models/cesm1.1/">Community Earth System
Model (CESM)</a> framework. It is <strong>strongly</strong>
recommended that you become familiar with running a multi-instance
experiment in CESM <strong>before</strong> you try to run DART/CLM.
The DART/CLM facility uses language and concepts that should be
familiar to CESM users. The DART/CLM capability is entirely
dependent on the multi-instance capability of CESM, first supported
in its entirety in CESM1.1.1. Consequently, this version or newer
is required to run CLM/DART. The <a href=
"http://www.cesm.ucar.edu/models/cesm1.1/clm/models/lnd/clm/doc/UsersGuide/clm_ug.pdf">
CLM User's Guide</a> is an excellent reference for CLM. <em class=
"new">As of (V7195) 3 October 2014, CESM1.2.1 is also
supported.</em><br>
<br>
DART uses the multi-instance capability of CESM, which means that
DART is not responsible for advancing the model. This GREATLY
simplifies the traditional DART workflow, but it means <em>CESM has
to stop and write out a restart file every time an assimilation is
required</em>. The multi-instance capability is very new to CESM
and we are in close collaboration with the CESM developers to make
using DART with CESM as easy as possible. While we strive to keep
DART requirements out of the model code, there are a few SourceMods
needed to run DART from within CESM. Appropriate SourceMods for
each CESM version are available at <a href=
"http://www.image.ucar.edu/pub/DART/CESM">http://www.image.ucar.edu/pub/DART/CESM</a>
and should be unpacked into your HOME directory. They will create a
<em class="file">~/cesm_?_?_?</em> directory with the appropriate
SourceMods structure. The ensuing scripts require these SourceMods
and expect them to be in your HOME directory.<br>
<br>
Our notes on how to set up, configure, build, and run CESM for an
assimilation experiment evolved into scripts. These scripts are not
intended to be a 'black box'; you will have to read and understand
them and modify them to your own purpose. They are heavily
commented -- in keeping with their origins as a set of notes. If
you would like to offer suggestions on how to improve those notes -
please send them to dart@ucar.edu - we'd love to hear them.</p>
<table border="0" cellpadding="5" width="100%" summary=
'script description'>
<thead align="left">
<tr>
<th>Script</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td><a href=
"shell_scripts/CESM1_1_1_setup_pmo">CESM1_1_1_setup_pmo</a></td>
<td>runs a single instance of CLM to harvest synthetic observations
for an OSSE or "perfect model" experiment. It requires a single CLM
state from a previous experiment and uses a specified DATM stream
for forcing. This parallels an assimilation experiment in that in
the multi-instance setting each CLM instance may use (should use?)
a unique DATM forcing. This script has almost nothing to do with
DART. There is one (trivial) section that records some
configuration information in the DART setup script, but that's
about it. This script should initially be run without DART to
ensure a working CESM environment.<br>
<br>
As of (V7195) 3 October 2014, this script demonstrates how to
create 'vector'-based CLM history files (which requires a bugfix)
and has an option to use a bugfixed snow grain-size code.<br>
<a href=
'http://bugs.cgd.ucar.edu/show_bug.cgi?id=1730'>http://bugs.cgd.ucar.edu/show_bug.cgi?id=1730</a><br>

<a href=
'http://bugs.cgd.ucar.edu/show_bug.cgi?id=1934'>http://bugs.cgd.ucar.edu/show_bug.cgi?id=1934</a></td>
</tr>
<tr>
<td><a href=
"shell_scripts/CESM1_2_1_setup_pmo">CESM1_2_1_setup_pmo</a></td>
<td>Is functionally identical to <em class=
"program">CESM1_1_1_setup_pmo</em> but is appropriate for the the
CESM 1_2_1 release, which supports both CLM 4 and CLM 4.5.</td>
</tr>
<tr>
<td><a href=
"shell_scripts/CESM1_1_1_setup_hybrid">CESM1_1_1_setup_hybrid</a></td>
<td>runs a multi-instance CLM experiment and can be used to perform
a free run or 'open loop' experiment. By default, each CLM instance
uses a unique DATM forcing. This script also has almost nothing to
do with DART. There is one (trivial) section that records some
configuration information in the DART setup script, but that's
about it. This script should initially be run without DART to
ensure a working CESM.<br>
<br>
As of (V7195) 3 October 2014, this script demonstrates how to
create 'vector'-based CLM history files (which requires a bugfix)
and has an option to use a bugfixed snow grain-size code.<br>
<a href=
'http://bugs.cgd.ucar.edu/show_bug.cgi?id=1730'>http://bugs.cgd.ucar.edu/show_bug.cgi?id=1730</a><br>

<a href=
'http://bugs.cgd.ucar.edu/show_bug.cgi?id=1934'>http://bugs.cgd.ucar.edu/show_bug.cgi?id=1934</a></td>
</tr>
<tr>
<td><a href=
"shell_scripts/CESM1_2_1_setup_hybrid">CESM1_2_1_setup_hybrid</a></td>
<td>Is functionally identical to <em class=
"program">CESM1_1_1_setup_hybrid</em> but is appropriate for the
the CESM 1_2_1 release, which supports both CLM 4 and CLM 4.5.</td>
</tr>
<tr>
<td><a href=
"shell_scripts/CESM_DART_config">CESM_DART_config</a></td>
<td>augments a CESM case with the bits and pieces required to run
DART. When either <em class="program">CESM1_?_1_setup_pmo</em> or
<em class="program">CESM1_?_1_setup_hybrid</em> gets executed,
<em class="program">CESM_DART_config</em> gets copied to the CESM
"caseroot" directory. It is designed such that you can execute it
at any time during a CESM experiment. When you do execute it, it
will build the DART executables and copy them into the CESM "bld"
directory, stage the run-time configurable <em class=
"file">input.nml</em> in the "caseroot" directory, etc. and also
<em>modifies</em> the CESM <em class="program">case.run</em> script
to call the DART scripts for assimilation or to harvest synthetic
observations.</td>
</tr>
</tbody>
</table>
<p>In addition to the script above, there are a couple scripts that
will either perform an assimilation (<a href=
"shell_scripts/assimilate.csh">assimilate.csh</a>) or harvest
observations for a perfect model experiment (<a href=
"shell_scripts/perfect_model.csh">perfect_model.csh</a>). These
scripts are designed to work on several compute platforms although
they require configuration, mainly to indicate the location of the
DART observation sequence files on your system.</p>
<h2>Pertinent details of the CLM gridcell.</h2>
<table border="0" cellpadding="5" width="100%" summary=
'CLM gridcell'>
<tbody valign="top">
<tr>
<td><a href=
"http://www.cesm.ucar.edu/models/clm/surface.heterogeneity.html"><img src="../../docs/images/clm_landcover.jpg"
alt="CLM gridcell breakdown" height="250"></a></td>
<td>"The land surface is represented by 5 primary sub-grid land
cover types (landunits: glacier, lake, wetland, urban, vegetated)
in each grid cell. The vegetated portion of a grid cell is further
divided into patches of plant functional types, each with its own
leaf and stem area index and canopy height. Each subgrid land cover
type and PFT patch is a separate column for energy and water
calculations." -- CLM documentation.<br>
<br>
The only location information available is at the gridcell level.
All landunits, columns, and PFTs in that gridcell have the same
location. This has ramifications for the forward observation
operators. If the observation metadata has information about land
use/land cover, it can be used to select only those patches that
are appropriate. Otherwise, an area-weighted average of ALL patches
in the gridcell is used to calculate the observation value for that
location.</td>
</tr>
</tbody>
</table>
<h2>A word about forward observation operators.</h2>
<p>"Simple" observations like snowcover fraction come directly from
the DART state. It is possible to configure the CLM history files
to contain the CLM estimates of some quantities (mostly flux tower
observations e.g, net ecosystem production, sensible heat flux,
latent heat flux) that are very complicated combinations of
portions of the CLM state. The forward observation operators for
these flux tower observations read these quantities from the CLM
<em class="file">.h1.</em> history file. The smaller the CLM
gridcell, the more likely it seems that these values will agree
with point observations.<br>
<br>
The prior and posterior values for these will naturally be
identical as the history file is unchanged by the assimilation.
Configuring the CLM user_nl_clm files to output the desired
quantities must be done at the first execution of CLM. As soon as
CONTINUE_RUN=TRUE, the namelist values for history file generation
are ignored. Because the history file creation is very flexible,
some additional information must be passed to DART to construct the
filename of the <em class="file">.h1.</em> history file needed for
any particular time.</p>
<h2>Major changes as of (V7195) 3 October 2014</h2>
<p>The DART state vector may be constructed in a much more flexible
way. Variables from two different CLM history files may also be
incorporated directly into the DART state - which should GREATLY
speed up the forward observation operators - and allow the
observation operators to be constructed in a more flexible manner
so that they can be used by any model capable of providing required
inputs. It is now possible to read some variables from the restart
file, some variables from a traditional history file, and some from
a 'vector-based' history file that has the same structure
(gridcell/landunit/column/pft) as the restart file. This should
allow more accurate forward observation operators since the
quantities are not gridcell-averaged a priori.<br>
<br>
Another namelist item has been added <em class=
"code">clm_vector_history_filename</em> to support the concept that
two history files can be supported. My intent was to have the
original history file (required for grid metadata) and another for
support of vector-based quantities in support of forward
observation operators. Upon reflection, I'm not sure I need two
different history files - BUT - I'm sure there will be a situation
where it comes in handy.<br>
<br>
The new namelist specification of what goes into the DART state
vector includes the ability to specify if the quantity should have
a lower bound, upper bound, or both, what file the variable should
be read from, and if the variable should be modified by the
assimilation or not. <strong>Only variables in the CLM restart file
will be candidates for updating.</strong> No CLM history files are
modified. <strong>It is important to know that the variables in the
DART diagnostic files <em class="file">preassim.nc</em> and
<em class="file">analysis.nc</em> will contain the unbounded
versions of ALL the variables specied in <em class=
"code">clm_variables</em>.</strong><br>
<br>
The example <em class="file">input.nml</em> <em class=
"code">model_nml</em> demonstrates how to construct the DART state
vector. The following table explains in detail each entry for
<em class="code">clm_variables</em>:</p>
<div>
<table border="0" cellpadding="4" width="100%" summary=
'clm_variables description'>
<thead align="left">
<tr>
<th>Column 1</th>
<th>Column 2</th>
<th>Column 3</th>
<th>Column 4</th>
<th>Column 5</th>
<th>Column 6</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>Variable name</td>
<td>DART KIND</td>
<td>minimum</td>
<td>maximum</td>
<td>filename</td>
<td>update</td>
</tr>
</tbody>
</table>
<br>
<br>
<table border="0" cellpadding="4" width="100%" summary=
'clm_variables description'>
<tbody valign="top">
<tr>
<td><strong>Column 1</strong></td>
<td>Variable name  </td>
<td>This is the CLM variable name as it appears in the CLM netCDF
file.</td>
</tr>
<tr>
<td><strong>Column 2</strong></td>
<td>DART KIND</td>
<td>This is the character string of the corresponding DART
KIND.</td>
</tr>
<tr>
<td><strong>Column 3</strong></td>
<td>minimum</td>
<td>If the variable is to be updated in the CLM restart file, this
specifies the minimum value. If set to 'NA', there is no minimum
value.</td>
</tr>
<tr>
<td><strong>Column 4</strong></td>
<td>maximum</td>
<td>If the variable is to be updated in the CLM restart file, this
specifies the maximum value. If set to 'NA', there is no maximum
value.</td>
</tr>
<tr>
<td><strong>Column 5</strong></td>
<td>filename</td>
<td>This specifies which file should be used to obtain the
variable.<br>
<em class=
"code">'restart'</em> =&gt; clm_restart_filename<br>
<em class=
"code">'history'</em> =&gt; clm_history_filename<br>
<em class=
"code">'vector' </em> =&gt; clm_vector_history_filename</td>
</tr>
<tr>
<td><strong>Column 6</strong></td>
<td>update</td>
<td>If the variable comes from the restart file, it may be updated
after the assimilation.<br>
<em class=
"code">'UPDATE'      </em> =&gt;
the variable in the restart file is updated.<br>
<em class="code">'NO_COPY_BACK'</em> =&gt; the variable in the
restart file remains unchanged.<br></td>
</tr>
</tbody>
</table>
</div>
<p>The following are only meant to be examples - they are not
scientifically validated. Some of these that are UPDATED are
probably diagnostic quantities, Some of these that should be
updated may be marked NO_COPY_BACK. There are multiple choices for
some DART kinds. This list is by no means complete.</p>
<pre>
       'livecrootc',  'QTY_ROOT_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'deadcrootc',  'QTY_ROOT_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'livestemc',   'QTY_STEM_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'deadstemc',   'QTY_STEM_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'livecrootn',  'QTY_ROOT_NITROGEN',          'NA', 'NA', 'restart', 'UPDATE',
       'deadcrootn',  'QTY_ROOT_NITROGEN',          'NA', 'NA', 'restart', 'UPDATE',
       'livestemn',   'QTY_STEM_NITROGEN',          'NA', 'NA', 'restart', 'UPDATE',
       'deadstemn',   'QTY_STEM_NITROGEN',          'NA', 'NA', 'restart', 'UPDATE',
       'litr1c',      'QTY_LEAF_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'litr2c',      'QTY_LEAF_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'litr3c',      'QTY_LEAF_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'soil1c',      'QTY_SOIL_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'soil2c',      'QTY_SOIL_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'soil3c',      'QTY_SOIL_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'soil4c',      'QTY_SOIL_CARBON',            'NA', 'NA', 'restart', 'UPDATE',
       'fabd',        'QTY_FPAR_DIRECT',            'NA', 'NA', 'restart', 'UPDATE',
       'fabi',        'QTY_FPAR_DIFFUSE',           'NA', 'NA', 'restart', 'UPDATE',
       'T_VEG',       'QTY_VEGETATION_TEMPERATURE', 'NA', 'NA', 'restart', 'UPDATE',
       'fabd_sun_z',  'QTY_FPAR_SUNLIT_DIRECT',     'NA', 'NA', 'restart', 'UPDATE',
       'fabd_sha_z',  'QTY_FPAR_SUNLIT_DIFFUSE',    'NA', 'NA', 'restart', 'UPDATE',
       'fabi_sun_z',  'QTY_FPAR_SHADED_DIRECT',     'NA', 'NA', 'restart', 'UPDATE',
       'fabi_sha_z',  'QTY_FPAR_SHADED_DIFFUSE',    'NA', 'NA', 'restart', 'UPDATE',
       'elai',        'QTY_LEAF_AREA_INDEX',        'NA', 'NA', 'restart', 'UPDATE',
</pre>
<p><strong>Only the first variable for a DART kind in the
clm_variables list will be used for the forward observation
operator.</strong> The following is perfectly legal (for CLM4, at
least):</p>
<pre>
clm_variables = 'LAIP_VALUE', 'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'restart' , 'UPDATE',
                'tlai',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'restart' , 'UPDATE',
                'elai',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'restart' , 'UPDATE',
                'ELAI',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                'LAISHA',     'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                'LAISUN',     'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                'TLAI',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                'TLAI',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'vector'  , 'NO_COPY_BACK'
   /
</pre>
<p>however, only LAIP_VALUE will be used to calculate the LAI when
an observation of LAI is encountered. All the other LAI variables
in the DART state will be modified by the assimilation based on the
relationship of LAIP_VALUE and the observation. Those coming from
the restart file and marked 'UPDATE' <strong>will</strong> be
updated in the CLM restart file.</p>
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>These namelists are read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;model_nml 
  clm_restart_filename         = 'clm_restart.nc',
  clm_history_filename         = 'clm_history.nc',
  clm_vector_history_filename  = 'clm_vector_history.nc',
  output_state_vector          = .false.,
  assimilation_period_days     = 2,
  assimilation_period_seconds  = 0,
  model_perturbation_amplitude = 0.2,
  calendar                     = 'Gregorian',
  debug                        = 0
  clm_variables  = 'frac_sno',    'QTY_SNOWCOVER_FRAC',         'NA' , 'NA', 'restart' , 'NO_COPY_BACK',
                   'H2OSNO',      'QTY_SNOW_WATER',             '0.0', 'NA', 'restart' , 'UPDATE',
                   'H2OSOI_LIQ',  'QTY_SOIL_MOISTURE',          '0.0', 'NA', 'restart' , 'UPDATE',
                   'H2OSOI_ICE',  'QTY_ICE',                    '0.0', 'NA', 'restart' , 'UPDATE',
                   'T_SOISNO',    'QTY_SOIL_TEMPERATURE',       'NA' , 'NA', 'restart' , 'UPDATE',
                   'SNOWDP',      'QTY_SNOW_THICKNESS',         'NA' , 'NA', 'restart' , 'UPDATE',
                   'LAIP_VALUE',  'QTY_LEAF_AREA_INDEX',        'NA' , 'NA', 'restart' , 'NO_COPY_BACK',
                   'cpool',       'QTY_CARBON',                 '0.0', 'NA', 'restart' , 'UPDATE',
                   'frootc',      'QTY_ROOT_CARBON',            '0.0', 'NA', 'restart' , 'UPDATE',
                   'leafc',       'QTY_LEAF_CARBON',            '0.0', 'NA', 'restart' , 'UPDATE',
                   'leafn',       'QTY_LEAF_NITROGEN',          '0.0', 'NA', 'restart' , 'UPDATE',
                   'NEP',         'QTY_NET_CARBON_PRODUCTION',  'NA' , 'NA', 'history' , 'NO_COPY_BACK',
                   'TV',          'QTY_VEGETATION_TEMPERATURE', 'NA' , 'NA', 'vector'  , 'NO_COPY_BACK',
                   'RH2M_R',      'QTY_SPECIFIC_HUMIDITY',      'NA' , 'NA', 'vector'  , 'NO_COPY_BACK',
                   'PBOT',        'QTY_SURFACE_PRESSURE',       'NA' , 'NA', 'vector'  , 'NO_COPY_BACK',
                   'TBOT',        'QTY_TEMPERATURE',            'NA' , 'NA', 'vector'  , 'NO_COPY_BACK'
   /
</pre></div>
<div>
<table border="0" cellpadding="3" width="100%" summary=
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
<td>clm_restart_filename</td>
<td>character(len=256)</td>
<td>this is the filename of the CLM restart file. The DART scripts
resolve linking the specific CLM restart file to this generic name.
This file provides the elements used to make up the DART state
vector. The variables are in their original landunit, column, and
PFT-based representations.</td>
</tr>
<tr>
<td>clm_history_filename</td>
<td>character(len=256)</td>
<td>this is the filename of the CLM <em class="file">.h0.</em>
history file. The DART scripts resolve linking the specific CLM
history file to this generic name. Some of the metadata needed for
the DART/CLM interfaces is contained only in this history file, so
it is needed for all DART routines.</td>
</tr>
<tr>
<td>clm_vector_history_filename</td>
<td>character(len=256)</td>
<td>this is the filename of a second CLM history file. The DART
scripts resolve linking the specific CLM history file to this
generic name. The default setup scripts actually create 3 separate
CLM history files, the <em class="file">.h2.</em> ones are linked
to this filename. It is possible to create this history file at the
same resolution as the restart file, which should make for better
forward operators. It is only needed if some of the variables
specified in <em class="code">clm_variables</em> come from this
file.</td>
</tr>
<tr>
<td>output_state_vector</td>
<td>logical</td>
<td>If .true. write state vector as a 1D array to the DART
diagnostic output files. If .false. break state vector up into
variables before writing to the output files.</td>
</tr>
<tr>
<td>assimilation_period_days,<br>
assimilation_period_seconds</td>
<td>integer</td>
<td>Combined, these specify the width of the assimilation window.
The current model time is used as the center time of the
assimilation window. All observations in the assimilation window
are assimilated. BEWARE: if you put observations that occur before
the beginning of the assimilation_period, DART will error out
because it cannot move the model 'back in time' to process these
observations.</td>
</tr>
<tr>
<td>model_perturbation_amplitude</td>
<td>real(r8)</td>
<td>Required by the DART interfaces, but not used by CLM.</td>
</tr>
<tr>
<td>calendar</td>
<td>character(len=32)</td>
<td>string specifying the calendar to use with DART. The CLM dates
will be interpreted with this same calendar. For assimilations with
real observations, this should be 'Gregorian'.</td>
</tr>
<tr>
<td>debug</td>
<td>integer</td>
<td>Set to 0 (zero) for minimal output. Successively higher values
generate successively more output. Not all values are important,
however. It seems I've only used values [3,6,7,8]. Go figure.</td>
</tr>
<tr>
<td><em class="removed">clm_state_variables</em><br>
clm_variables</td>
<td>character(:,6)</td>
<td>Strings that identify the CLM variables, their DART kind, the
min &amp; max values, what file to read from, and whether or not
the file should be updated after the assimilation. <em class=
"removed">Only CLM variable names in the CLM restart file are
valid.</em> The DART kind must be one found in the <em class=
"file">DART/assimilation_code/modules/observations/obs_kind_mod.f90</em>
AFTER it gets built by <em class="program">preprocess</em>. Most of
the land observation kinds are specified by <em class=
"file">DART/observations/forward_operators/obs_def_land_mod.f90</em>
and <em class=
"file">DART/observations/forward_operators/obs_def_tower_mod.f90</em>,
so they should be specified in the preprocess_nml:input_files
variable.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<div class="namelist">
<pre>
&amp;obs_def_tower_nml
   casename    = '../clm_dart',
   hist_nhtfrq = -24,
   debug       = .false.
   /
</pre></div>
<div>
<table border="0" cellpadding="3" width="100%" summary=
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
<td>casename</td>
<td>character(len=256)</td>
<td>this is the name of the CESM case. It is used by the forward
observation operators to help construct the filename of the CLM
<em class="file">.h1.</em> history files for the flux tower
observations. When the <em class="file">input.nml</em> gets staged
in the CASEROOT directory by <em class=
"program">CESM_DART_config</em>, the appropriate value should
automatically be inserted.</td>
</tr>
<tr>
<td>hist_nhtfrq</td>
<td>integer</td>
<td>this is the same value as in the CLM documentation. A negative
value indicates the number of hours contained in the <em class=
"file">.h1.</em> file. This value is needed to constuct the right
<em class="file">.h1.</em> filename. When the <em class=
"file">input.nml</em> gets staged in the CASEROOT directory by
<em class="program">CESM_DART_config</em>, the appropriate value
should automatically be inserted. Due to the large number of ways
of specifying the CLM history file information, the correct value
here is very dependent on how the case was configured. You would be
wise to check it.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>Set to .false. for minimal output.</td>
</tr>
</tbody>
</table>
</div>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED (directly)</h2>
<pre>
types_mod
time_manager_mod
threed_sphere/location_mod
utilities_mod
obs_kind_mod
obs_def_land_mod
obs_def_tower_mod
random_seq_mod
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
<h2>PUBLIC INTERFACES - Required</h2>
<table summary="list of all required public interfaces">
<tr>
<td><em class="call">use model_mod, only :</em></td>
<td><a href="#get_model_size">get_model_size</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#adv_1step">adv_1step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_state_meta_data">get_state_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#model_interpolate">model_interpolate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_model_time_step">get_model_time_step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#static_init_model">static_init_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#end_model">end_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_time">init_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_conditions">init_conditions</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_write_model_atts">nc_write_model_atts</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_write_model_vars">nc_write_model_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#pert_model_state">pert_model_state</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_close_maxdist_init">get_close_maxdist_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs_init">get_close_obs_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs">get_close_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ens_mean_for_model">ens_mean_for_model</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_model_size" id="get_model_size"></a><br>
<div class="routine"><em class="call">model_size = get_model_size(
)</em>
<pre>
integer :: <em class="code">get_model_size</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the length of the model state vector.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">model_size</em></td>
<td>The length of the model state vector.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="adv_1step" id="adv_1step"></a><br>
<div class="routine"><em class="call">call adv_1step(x, time)</em>
<pre>
real(r8), dimension(:), intent(inout) :: <em class="code">x</em>
type(time_type),        intent(in)    :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Advances the model for a single time step. The time associated
with the initial model state is also input although it is not used
for the computation.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">x</em></td>
<td>State vector of length model_size.</td>
</tr>
<tr>
<td><em class="code">time   </em></td>
<td>Specifies time of the initial model state.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_state_meta_data" id="get_state_meta_data"></a><br>
<div class="routine"><em class="call">call get_state_meta_data
(index_in, location, <em class=
"optionalcode">[, var_type]</em> )</em>
<pre>
integer,             intent(in)  :: <em class="code">index_in</em>
type(location_type), intent(out) :: <em class="code">location</em>
integer, optional,   intent(out) :: <em class=
"optionalcode"> var_type </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns metadata about a given element, indexed by index_in, in
the model state vector. The location defines where the state
variable is located.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">index_in   </em></td>
<td>Index of state vector element about which information is
requested.</td>
</tr>
<tr>
<td><em class="code">location</em></td>
<td>The location of state variable element.</td>
</tr>
<tr>
<td><em class="optionalcode">var_type</em></td>
<td>The generic DART kind of the state variable element.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="model_interpolate" id="model_interpolate"></a><br>
<div class="routine"><em class="call">call model_interpolate(x,
location, itype, obs_val, istatus)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">x</em>
type(location_type),    intent(in)  :: <em class=
"code">location</em>
integer,                intent(in)  :: <em class="code">itype</em>
real(r8),               intent(out) :: <em class=
"code">obs_val</em>
integer,                intent(out) :: <em class=
"code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given model state, returns the value interpolated to a given
location.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">x</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td><em class="code">location   </em></td>
<td>Location to which to interpolate.</td>
</tr>
<tr>
<td><em class="code">itype</em></td>
<td>Not used.</td>
</tr>
<tr>
<td><em class="code">obs_val</em></td>
<td>The interpolated value from the model.</td>
</tr>
<tr>
<td><em class="code">istatus</em></td>
<td>If the interpolation was successful <em class="code">istatus =
0</em>. If <em class="code">istatus /= 0</em> the interpolation
failed. Values less than zero are reserved for DART.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_model_time_step" id="get_model_time_step"></a><br>
<div class="routine"><em class="call">var =
get_model_time_step()</em>
<pre>
type(time_type) :: <em class="code">get_model_time_step</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the time step (forecast length) of the model;</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">var   </em></td>
<td>Smallest time step of model.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="static_init_model" id="static_init_model"></a><br>
<div class="routine"><em class="call">call
static_init_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Used for runtime initialization of model; reads namelist,
initializes model parameters, etc. This is the first call made to
the model by any DART-compliant assimilation routine.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary=""></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_model" id="end_model"></a><br>
<div class="routine"><em class="call">call end_model()</em></div>
<div class="indent1"><!-- Description -->
<p>A stub.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary=""></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_time" id="init_time"></a><br>
<div class="routine"><em class="call">call init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the time at which the model will start if no input
initial conditions are to be used. This is used to spin-up the
model from rest.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">time   </em></td>
<td>Initial model time.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_conditions" id="init_conditions"></a><br>
<div class="routine"><em class="call">call init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class="code">x</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns default initial conditions for the model; generally used
for spinning up initial model states.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">x   </em></td>
<td>Initial conditions for state vector.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_atts" id="nc_write_model_atts"></a><br>
<div class="routine"><em class="call">ierr =
nc_write_model_atts(ncFileID)</em>
<pre>
integer             :: <em class="code">nc_write_model_atts</em>
integer, intent(in) :: <em class="code">ncFileID</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Function to write model specific attributes to a netCDF file. At
present, DART is using the NetCDF format to output diagnostic
information. This is not a requirement, and models could choose to
provide output in other formats. This function writes the metadata
associated with the model to a NetCDF file opened to a file
identified by ncFileID.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">ncFileID   </em></td>
<td>Integer file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td><em class="code">ierr</em></td>
<td>Returns a 0 for successful completion.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_vars" id="nc_write_model_vars"></a><br>
<div class="routine"><em class="call">ierr =
nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)</em>
<pre>
integer                            :: <em class=
"code">nc_write_model_vars</em>
integer,                intent(in) :: <em class=
"code">ncFileID</em>
real(r8), dimension(:), intent(in) :: <em class=
"code">statevec</em>
integer,                intent(in) :: <em class=
"code">copyindex</em>
integer,                intent(in) :: <em class=
"code">timeindex</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes a copy of the state variables to a netCDF file. Multiple
copies of the state for a given time are supported, allowing, for
instance, a single file to include multiple ensemble estimates of
the state.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">ncFileID</em></td>
<td>file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td><em class="code">statevec</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td><em class="code">copyindex   </em></td>
<td>Integer index of copy to be written.</td>
</tr>
<tr>
<td><em class="code">timeindex</em></td>
<td>The timestep counter for the given state.</td>
</tr>
<tr>
<td><em class="code">ierr</em></td>
<td>Returns 0 for normal completion.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="pert_model_state" id="pert_model_state"></a><br>
<div class="routine"><em class="call">call pert_model_state(state,
pert_state, interf_provided)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">state</em>
real(r8), dimension(:), intent(out) :: <em class=
"code">pert_state</em>
logical,                intent(out) :: <em class=
"code">interf_provided</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a model state, produces a perturbed model state.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">state</em></td>
<td>State vector to be perturbed.</td>
</tr>
<tr>
<td><em class="code">pert_state</em></td>
<td>Perturbed state vector: NOT returned.</td>
</tr>
<tr>
<td><em class="code">interf_provided   </em></td>
<td>Returned false; interface is not implemented.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_maxdist_init" id=
"get_close_maxdist_init"></a><br>
<div class="routine"><em class="call">call
get_close_maxdist_init(gc, maxdist)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
real(r8),             intent(in)    :: <em class=
"code">maxdist</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>In distance computations any two locations closer than the given
<em class="code">maxdist</em> will be considered close by the
<em class="code">get_close_obs()</em> routine. Pass-through to the
3D Sphere locations module. See <a href=
"../../assimilation_code/location/threed_sphere/location_mod.html#get_close_maxdist_init">
get_close_maxdist_init()</a> for the documentation of this
subroutine.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs_init" id="get_close_obs_init"></a><br>
<div class="routine"><em class="call">call get_close_obs_init(gc,
num, obs)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
integer,              intent(in)    :: <em class="code">num</em>
type(location_type),  intent(in)    :: <em class=
"code">obs(num)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3D Sphere locations module. See <a href=
"../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs_init">
get_close_obs_init()</a> for the documentation of this
subroutine.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs" id="get_close_obs"></a><br>
<div class="routine"><em class="call">call get_close_obs(gc,
base_obs_loc, base_obs_kind, obs, obs_kind, num_close, close_ind
<em class="optionalcode">[, dist]</em>)</em>
<pre>
type(get_close_type), intent(in)  :: <em class="code">gc</em>
type(location_type),  intent(in)  :: <em class=
"code">base_obs_loc</em>
integer,              intent(in)  :: <em class=
"code">base_obs_kind</em>
type(location_type),  intent(in)  :: <em class="code">obs(:)</em>
integer,              intent(in)  :: <em class=
"code">obs_kind(:)</em>
integer,              intent(out) :: <em class=
"code">num_close</em>
integer,              intent(out) :: <em class=
"code">close_ind(:)</em>
real(r8), optional,   intent(out) :: <em class=
"optionalcode">dist(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3D Sphere locations module. See <a href=
"../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs">
get_close_obs()</a> for the documentation of this subroutine.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="ens_mean_for_model" id="ens_mean_for_model"></a><br>
<div class="routine"><em class="call">call
ens_mean_for_model(ens_mean)</em>
<pre>
real(r8), dimension(:), intent(in) :: <em class=
"code">ens_mean</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>A NULL INTERFACE in this model.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">ens_mean   </em></td>
<td>State vector containing the ensemble mean.</td>
</tr>
</tbody>
</table>
</div>
<!--==================================================================-->
<!-- list of all optional public interfaces                           -->
<!--==================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES - Optional</h2>
<table summary="list of all optional public interfaces">
<tr>
<td><em class="call">use model_mod, only :</em></td>
<td><a href="#get_gridsize">get_gridsize</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#clm_to_dart_state_vector">clm_to_dart_state_vector</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#sv_to_restart_file">sv_to_restart_file</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_clm_restart_filename">get_clm_restart_filename</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_state_time">get_state_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_grid_vertval">get_grid_vertval</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#compute_gridcell_value">compute_gridcell_value</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#gridcell_components">gridcell_components</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#DART_get_var">DART_get_var</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_model_time">get_model_time</a></td>
</tr>
</table>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_gridsize" id="get_gridsize"></a><br>
<div class="routine"><em class="call">call get_gridsize(num_lon,
num_lat, num_lev)</em>
<pre>
integer, intent(out) :: <em class=
"code">num_lon, num_lat, num_lev</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the number of longitudes, latitudes, and total number of
levels in the CLM state.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">num_lon</em></td>
<td>The number of longitude grid cells in the CLM state. This comes
from the CLM history file.</td>
</tr>
<tr>
<td><em class="code">num_lat</em></td>
<td>The number of latitude grid cells in the CLM state. This comes
from the CLM history file.</td>
</tr>
<tr>
<td><em class="code">num_lev</em></td>
<td>The number of levels grid cells in the CLM state. This comes
from 'nlevtot' in the CLM restart file.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="clm_to_dart_state_vector" id=
"clm_to_dart_state_vector"></a><br>
<div class="routine"><em class="call">call
clm_to_dart_state_vector(state_vector, restart_time)</em>
<pre>
real(r8),         intent(inout) :: <em class=
"code">state_vector(:)</em>
type(time_type),  intent(out)   :: <em class=
"code">restart_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads the current time and state variables from CLM netCDF
file(s) and packs them into a DART state vector. This MUST happen
in the same fashion as the metadata arrays are built. The variables
are specified by <em class="code">model_nml:clm_variables</em>.
Each variable specifies its own file of origin. If there are
multiple times in the file of origin, only the time that matches
the restart file are used.<br>
<br></p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">state_vector</em></td>
<td>The DART state vector.</td>
</tr>
<tr>
<td><em class="code">restart_time</em></td>
<td>The valid time of the CLM state.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="sv_to_restart_file" id="sv_to_restart_file"></a><br>
<div class="routine"><em class="call">call
sv_to_restart_file(state_vector, filename, dart_time)</em>
<pre>
real(r8),         intent(in) :: <em class=
"code">state_vector(:)</em>
character(len=*), intent(in) :: <em class="code">filename</em>
type(time_type),  intent(in) :: <em class="code">dart_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This routine updates the CLM restart file with the posterior
state from the assimilation. Some CLM variables that are useful to
include in the DART state (frac_sno, for example) are diagnostic
quantities and are not used for subsequent model advances. The
known diagnostic variables are NOT updated. If the values created
by the assimilation are outside physical bounds, or if the original
CLM value was 'missing', the <em class=
"program">vector_to_prog_var()</em> subroutine ensures that the
values in the original CLM restart file are <strong>not
updated</strong>.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">state_vector</em></td>
<td>The DART state vector containing the state modified by the
assimilation.</td>
</tr>
<tr>
<td><em class="code">filename</em></td>
<td>The name of the CLM restart file. <strong>The contents of some
of the variables will be overwritten with new values.</strong></td>
</tr>
<tr>
<td><em class="code">dart_time</em></td>
<td>The valid time of the DART state. This has to match the time in
the CLM restart file.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_clm_restart_filename" id=
"get_clm_restart_filename"></a><br>
<div class="routine"><em class="call">call
get_clm_restart_filename( filename )</em>
<pre>
character(len=*), intent(out) :: <em class="code">filename</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>provides access to the name of the CLM restart file to routines
outside the scope of this module.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">filename</em></td>
<td>The name of the CLM restart file.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_state_time" id="get_state_time"></a><br>
<div class="routine"><em class="call">time =
get_state_time(file_handle)</em>
<pre>
integer,          intent(in) :: <em class="code">file_handle</em> 
character(len=*), intent(in) :: <em class="code">file_handle</em> 
type(time_type)              :: <em class=
"code">get_state_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This routine has two interfaces - one for an integer input, one
for a filename. They both return the valid time of the model state
contained in the file. The file referenced is the CLM restart file
in netCDF format.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">file_handle</em></td>
<td>If specified as an integer, it must be the netCDF file
identifier from nf90_open(). If specified as a filename, the name
of the netCDF file.</td>
</tr>
<tr>
<td><em class="code">time</em></td>
<td>A DART time-type that contains the valid time of the model
state in the CLM restart file.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_grid_vertval" id="get_grid_vertval"></a><br>
<div class="routine"><em class="call">call get_grid_vertval(x,
location, varstring, interp_val, istatus)</em>
<pre>
real(r8),            intent(in)  :: <em class="code">x(:)</em>
type(location_type), intent(in)  :: <em class="code">location</em>
character(len=*),    intent(in)  :: <em class="code">varstring</em>
real(r8),            intent(out) :: <em class=
"code">interp_val</em>
integer,             intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Calculate the value of quantity at depth. The gridcell value at
the levels above and below the depth of interest are calculated and
then the value for the desired depth is linearly interpolated. Each
gridcell value is an area-weighted value of an unknown number of
column- or pft-based quantities. This is one of the workhorse
routines for <em class="program">model_interpolate()</em>.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">x</em></td>
<td>The DART state vector.</td>
</tr>
<tr>
<td><em class="code">location</em></td>
<td>The location of the desired quantity.</td>
</tr>
<tr>
<td><em class="code">varstring</em></td>
<td>The CLM variable of interest - this must be part of the DART
state. e.g, T_SOISNO, H2OSOI_LIQ, H2OSOI_ICE ...</td>
</tr>
<tr>
<td><em class="code">interp_val</em></td>
<td>The quantity at the location of interest.</td>
</tr>
<tr>
<td><em class="code">istatus</em></td>
<td>error code. 0 (zero) indicates a successful interpolation.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="compute_gridcell_value" id=
"compute_gridcell_value"></a><br>
<div class="routine"><em class="call">call
compute_gridcell_value(x, location, varstring, interp_val,
istatus)</em>
<pre>
real(r8),            intent(in)  :: <em class="code">x(:)</em>
type(location_type), intent(in)  :: <em class="code">location</em>
character(len=*),    intent(in)  :: <em class="code">varstring</em>
real(r8),            intent(out) :: <em class=
"code">interp_val</em>
integer,             intent(out) :: <em class="code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Calculate the value of a CLM variable in the DART state vector
given a location. Since the CLM location information is only
available at the gridcell level, all the columns in a gridcell are
area-weighted to derive the value for the location. This is one of
the workhorse routines for <em class=
"program">model_interpolate()</em>, and only select CLM variables
are currently supported. Only CLM variables that have no vertical
levels may use this routine.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">x</em></td>
<td>The DART state vector.</td>
</tr>
<tr>
<td><em class="code">location</em></td>
<td>The location of the desired quantity.</td>
</tr>
<tr>
<td><em class="code">varstring</em></td>
<td>The CLM variable of interest - this must be part of the DART
state. e.g, frac_sno, leafc, ZWT ...</td>
</tr>
<tr>
<td><em class="code">interp_val</em></td>
<td>The quantity at the location of interest.</td>
</tr>
<tr>
<td><em class="code">istatus</em></td>
<td>error code. 0 (zero) indicates a successful interpolation.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="gridcell_components" id="gridcell_components"></a><br>
<div class="routine"><em class="call">call
gridcell_components(varstring)</em>
<pre>
character(len=*), intent(in) :: <em class="code">varstring</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a utility routine that helps identify how many land
units,columns, or PFTs are in each gridcell for a particular
variable. Helps answer exploratory questions about which gridcells
are appropriate to test code. The CLM variable is read from the CLM
restart file.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">varstring</em></td>
<td>The CLM variable name of interest.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="DART_get_var" id="DART_get_var"></a><br>
<div class="routine"><em class="call">call DART_get_var(ncid,
varname, datmat)</em>
<pre>
integer,                  intent(in)  :: <em class="code">ncid</em>
character(len=*),         intent(in)  :: <em class=
"code">varname</em>
real(r8), dimension(:),   intent(out) :: <em class=
"code">datmat</em>
real(r8), dimension(:,:), intent(out) :: <em class=
"code">datmat</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads a 1D or 2D variable of 'any' type from a netCDF file and
processes and applies the offset/scale/FillValue attributes
correctly.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">ncid</em></td>
<td>The netCDF file identifier to an open file. ncid is the output
from a nf90_open() call.</td>
</tr>
<tr>
<td><em class="code">varname</em></td>
<td>The name of the netCDF variable of interest. The variables can
be integers, floats, or doubles.</td>
</tr>
<tr>
<td><em class="code">datmat</em></td>
<td>The shape of datmat must match the shape of the netCDF
variable. Only 1D or 2D variables are currently supported.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_model_time" id="get_model_time"></a><br>
<div class="routine"><em class="call">model_time = get_model_time(
)</em>
<pre>
integer :: <em class="code">get_model_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the valid time of the model state vector.</p>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<td><em class="code">model_time</em></td>
<td>The valid time of the model state vector.</td>
</tr>
</tbody>
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
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<th align="left">filename</th>
<th align="left">purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read the model_mod namelist</td>
</tr>
<tr>
<td>clm_restart.nc</td>
<td>both read and modified by the CLM model_mod</td>
</tr>
<tr>
<td>clm_history.nc</td>
<td>read by the CLM model_mod for metadata purposes.</td>
</tr>
<tr>
<td>*.h1.* history files</td>
<td>may be read by the obs_def_tower_mod for observation operator
purposes.</td>
</tr>
<tr>
<td>dart_log.out</td>
<td>the run-time diagnostic output</td>
</tr>
<tr>
<td>dart_log.nml</td>
<td>the record of all the namelists actually USED - contains the
default values</td>
</tr>
</tbody>
</table>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<p><a href=
"http://www.cesm.ucar.edu/models/cesm1.1/clm/models/lnd/clm/doc/UsersGuide/clm_ug.pdf">
CLM User's Guide</a> is an excellent reference for CLM.</p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="3" width="100%"
summary="error situations">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td>nc_write_model_atts<br>
nc_write_model_vars</td>
<!-- message -->
<td>Various netCDF-f90 interface error messages</td>
<!-- comment -->
<td>From one of the netCDF calls in the named routine</td>
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
<p>Almost too many to list.</p>
<ol>
<li>Implement a robust update_snow() routine that takes the
modified SWE and repartitions it into the respective snow layers in
a manner that works with both CLM4 and CLM4.5. This may mean
modifying the clm_variables list to contain SNOWDP, H2OSOI_LIQ,
H2OSOI_ICE, T_SOISNO, and others that may not be in the UPDATE
list.</li>
<li>Implement a fast way to get the quantities needed for the
calculation of radiative transfer models - needs a whole column of
CLM variables, redundant if multiple frequencies are used.</li>
<li>Figure out what to do when one or more of the ensemble members
does not have snow/leaves/etc. when the observation indicates there
should be. Ditto for removing snow/leaves/etc. when the observation
indicates otherwise.</li>
<li>Right now, the soil moisture observation operator is used by
the COSMOS code to calculate the expected neutron intensity counts.
This is the right idea, however, the COSMOS forward operator uses
m3/m3 and the CLM units are kg/m2 ... I have not checked to see if
they are, in fact, identical. This brings up a bigger issue in that
the soil moisture observation operator would also be used to
calculate whatever a TDT probe or ??? would measure. What units are
they in? Can one operator support both?</li>
<li>...</li>
</ol>
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
