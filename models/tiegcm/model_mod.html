<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
"http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module model_mod (TIEGCM)</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP"></a>
<h1>TIEGCM</h1>
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
/ <a href="#FuturePlans">PLANS</a> / / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>This is the DART interface to the Thermosphere Ionosphere
Electrodynamic General Circulation Model (<a href=
"http://www.hao.ucar.edu/modeling/tgcm/tie.php">TIEGCM</a>), which
is a community model developed at the NCAR High Altitude
Observatory. TIEGCM is widely used by the space physics and
aeronomy community and is one of the most well-validated models of
the Earth's upper atmosphere. DART/TIEGCM has been used to
assimilate neutral mass density retrieved from satellite-borne
accelerometers and electon density obtained from ground-based and
space-based GNSS signals. Unlike other ionospheric data
assimilation applications, this approach allows simultaneous
assimilation of thermospheric and ionospheric parameters by taking
advantage of the coupling of plasma and neutral constituents
described in TIEGCM. DART/TIEGCM's demonstrated capability to infer
under-observed thermospheric parameters from abundant electron
density observations has important implications for the future of
upper atmosphere research.<br>
<br>
DART is designed so that the TIEGCM source code can be used with no
modifications, as DART runs TIEGCM as a completely separate
executable. The TIEGCM source code and restart files are
<strong>not</strong> included in DART, so you must obtain them from
the NCAR High Altitude Observatory (<a href=
"http://www.hao.ucar.edu/modeling/tgcm/download.php">download
website</a>). It is <strong>strongly</strong> recommended that you
become familiar with running TIEGCM <strong>before</strong> you try
to run DART/TIEGCM (See the <a href=
"http://www.hao.ucar.edu/modeling/tgcm/doc/userguide/html">TIEGCM
User's Guide</a>). Some assumptions are made about the mannner in
which TIEGCM is run: (1) There can only be 1 each of the TIEGCM
primary (restart) and secondary NetCDF history files. The TIEGCM
primary history files contain the prognostic variables necessary to
restart the model, while the secondary history files contain
diagnostic variables; (2) The last timestep in the restart file is
the only timestep which is converted to a DART state vector, and
only the last timestep in the TIEGCM primary file is ever modified
by DART. The TIEGCM variables to be included in a DART state
vector, and possibly updated by the assimilation, are specified in
the DART namelist. (Some of the TIEGCM variables used to compute
observation priors need not to be updated.) It is required to
associate the TIEGCM variable name with a 'generic' DART
counterpart (e.g., <em class="code">NE</em> is <em class=
"code">QTY_ELECTRON_DENSITY</em>). The composition of the DART
state vector and which variables get updated in the TIEGCM primary
file are under complete user control.<br>
<br>
In the course of a filtering experiment, it is necessary to make a
short forecast with TIEGCM. DART writes out an ancillary file with
the information necessary to advance TIEGCM to the required time.
The DART script <em class="program">advance_model.csh</em> reads
this information and modifies the TIEGCM namelist <em class=
"file">tiegcm.nml</em> such that TIEGCM runs upto the requested
time when DART assimilates the next set of observations. The run
scripts <em class="program">run_filter.csh</em> and <em class=
"program">run_perfect_model_obs.csh</em> are configured to run
under the LSF queueing system. The scripting examples exploit an
'embarassingly-simple' parallel paradigm in that each TIEGCM
instance is a single-threaded executable and all ensemble members
may run simultaneously. To use these run scripts, the TIECGM
executable needs to be compiled with no MPI option. As such, there
is an advantage to matching the ensemble size to the number of
tasks. Requesting more tasks than the number of ensemble members
may speed up the DART portion of an assimilation (i.e., <em class=
"program">filter</em>) but will not make the model advance faster.
The <em class="program">filter</em> may be compiled with MPI and
can exploit all available tasks.</p>
<h2>Quickstart guide to running</h2>
<p>It is important to understand basic DART nomenclature and
mechanisms. Please take the time to read and run the <a href=
"../../tutorial/index.pdf">DART tutorial</a>.<br>
<br>
Both <em class="program">run_filter.csh</em> and <em class=
"program">run_perfect_model_obs.csh</em> are heavily internally
commented. Please read and understand the scripts. The overall
process is to</p>
<ol>
<li>Specify resources (wall-clock time, number of nodes, tasks that
sort of thing).</li>
<li>Set shell variables to identify the location of the DART
exectuables, the TIEGCM executables, initial ensemble, etc.</li>
<li>Establish a temporary working directory for the
experiment.</li>
<li>Populate that directory with the initial ensemble and required
namelists.</li>
<li>Convert each TIEGCM ensemble member to a DART initial
conditions file.</li>
<li>Run either <em class="program">filter</em> or <em class=
"program">run_perfect_model_obs.csh</em>.</li>
</ol>
<ol>
<li style="list-style: none"><em class=
"program">perfect_model_obs</em> will</li>
<li>Check for any desired observations at the current time of the
model state and create the synthetic observations for all
observation times in the specified assimilation window. If the
model needs to be advanced, it then</li>
<li>creates a unique run-time directory for the model advance,</li>
<li>copies the required information into that directory,</li>
<li>conveys the desired forecast stopping time to TIEGCM via the
<em class="file">tiegcm.nml</em> and</li>
<li>runs a single executable of TIEGCM.</li>
<li>Steps 1-5 are repeated until the input DART observation
sequence file has been exhausted.</li>
</ol>
<ol>
<li style="list-style: none"><em class="program">filter</em>
will</li>
<li>Check for any desired observations at the current time of the
model state and assimilates all the observations in the specified
assimilation window. If the model needs to be advanced, it
then</li>
<li>creates a set of run-time directories, one for each task. A
single task may be responsible for advancing more than one TIEGCM
instance. If so, each instance is done serially, one after another.
See the documentation for <a href=
"../../docs/html/filter_async_modes.html">options for running DART
in parallel</a>.</li>
<li>Copy the required information into that directory.</li>
<li>Update the TIEGCM restart file with the most current
DART-modified state and convey the desired forecast stopping time
to TIEGCM via the unique <em class="file">tiegcm.nml</em> for this
ensemble member.</li>
<li>Runs a single executable of TIEGCM.</li>
<li>Steps 1-5 are repeated until the input DART observation
sequence file</li>
</ol>
<h3>What to check when things go wrong.</h3>
<p>The scripts are designed to send email to the user that contains
the run-time output from the script. Check that first. If that does
not provide the information needed, go to the run directory (i.e.
CENTRALDIR) and check the <em class="file">dart_log.out</em>. It
usually provides the same information as the email, but sometimes
it can help. If that does not help, go to any of the
CENTRALDIR/<em class="file">advance_temp<i>nnnn</i></em>
directories and read the <em class=
"file">log_advance.<i>nnnn</i>.txt</em> file.</p>
<p>
<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->
 <a name="Namelist"></a></p>
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
&amp;model_nml 
   output_state_vector         = .false.
   tiegcm_restart_file_name    = 'tiegcm_restart_p.nc'
   tiegcm_secondary_file_name  = 'tiegcm_s.nc'
   tiegcm_namelist_file_name   = 'tiegcm.nml'
   assimilation_period_seconds = 3600
   estimate_f10_7              = .false.
   debug                       = 1
   variables = 'NE',    'QTY_ELECTRON_DENSITY',          '1000.0',  'NA',      'restart',    'UPDATE'
               'OP',    'QTY_DENSITY_ION_OP',            'NA',      'NA',      'restart',    'UPDATE',
               'TI',    'QTY_TEMPERATURE_ION',           'NA',      'NA',      'restart',    'UPDATE',
               'TE',    'QTY_TEMPERATURE_ELECTRON',      'NA',      'NA',      'restart',    'UPDATE',
               'OP_NM', 'QTY_DENSITY_ION_OP',            'NA',      'NA',      'restart',    'UPDATE',
               'O1',    'QTY_ATOMIC_OXYGEN_MIXING_RATIO','0.00001', '0.99999', 'secondary',  'NO_COPY_BACK',
               'O2',    'QTY_MOLEC_OXYGEN_MIXING_RATIO', '0.00001', '0.99999', 'secondary',  'NO_COPY_BACK',
               'TN',    'QTY_TEMPERATURE',               '0.0',     '6000.0',  'secondary',  'NO_COPY_BACK',
               'ZG',    'QTY_GEOMETRIC_HEIGHT',          'NA',      'NA',      'secondary',  'NO_COPY_BACK',
               'VTEC',  'QTY_VERTICAL_TEC',              'NA',      'NA',      'calculate',  'NO_COPY_BACK'
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
<td>output_state_vector</td>
<td>logical</td>
<td>If .true. write state vector as a 1D array to the DART
diagnostic output files. If .false. break state vector up into
variables before writing to the output files.</td>
</tr>
<tr>
<td>tiegcm_restart_file_name</td>
<td>character(len=256)</td>
<td>The TIEGCM restart file name.</td>
</tr>
<tr>
<td>tiegcm_secondary_file_name</td>
<td>character(len=256)</td>
<td>The TIEGCM secondary file name.</td>
</tr>
<tr>
<td>tiegcm_namelist_file_name</td>
<td>character(len=256)</td>
<td>The TIEGCM namelist file name.</td>
</tr>
<tr>
<td>assimilation_period_seconds</td>
<td>integer</td>
<td>This specifies the width of the assimilation window. The
current model time is used as the center time of the assimilation
window. All observations in the assimilation window are
assimilated. BEWARE: if you put observations that occur before the
beginning of the assimilation_period, DART will error out because
it cannot move the model 'back in time' to process these
observations. <em class="code">assimilation_period_seconds</em>
must be an integer number of TIEGCM dynamical timesteps (as
specified by tiegcm.nml:STEP) AND be able to be expressed by
tiegcm.nml:STOP. Since STOP has three components: day-of-year,
hour, and minute, the <em class=
"code">assimilation_period_seconds</em> must be an integer number
of minutes.</td>
</tr>
<tr>
<td>estimate_f10_7</td>
<td>logical</td>
<td>Switch to specify that the f10.7 index should be estimated by
augmenting the DART state vector with a scalar. The location of the
f10.7 index is taken to be longitude of local noon and latitude
zero. WARNING: this is provided with no guarantees. Please read the
comments in <em class="file">model_mod.f90</em> and act
accordingly.</td>
</tr>
<tr>
<td>debug</td>
<td>integer</td>
<td>Set to 0 (zero) for minimal output. Successively larger values
generate successively more output.</td>
</tr>
<tr>
<td>variables</td>
<td>character(:,6)</td>
<td>Strings that identify the TIEGCM variables, their DART kind,
the min &amp; max values, what file to read from, and whether or
not the file should be updated after the assimilation. The DART
kind must be one found in the <em class=
"file">DART/assimilation_code/modules/observations/obs_kind_mod.f90</em>
AFTER it gets built by <em class="program">preprocess</em>. Most of
the upper atmosphere observation kinds are specified by <em class=
"file">DART/observations/forward_operators/obs_def_upper_atm_mod.f90</em>,
so it should be specified in the <em class=
"file">preprocess_nml</em>:<em class="code">input_files</em>
variable. Since TIEGCM has an entire class of variables (all the
variables that end in <em class="code">_NM</em>) that are simply 1
dynamical timestep behind the variables at the output time, it is
<strong>imperative</strong> that these variables be specified to
occur AFTER their counterparts in the DART namelist. This will
ensure that the most current variables are used in the calculation
of the forward observation operators.<br>
<table border="0" cellpadding="3" width="100%" summary=
'variable description'>
<tr>
<td valign="top"><em class="code">variables(:,1)</em></td>
<td valign="top">Specifies the TIEGCM variable name in the netCDF
file.</td>
</tr>
<tr>
<td valign="top"><em class="code">variables(:,2)</em></td>
<td valign="top">Specifies the DART kind for that variable.</td>
</tr>
<tr>
<td valign="top"><em class="code">variables(:,3)</em></td>
<td valign="top">Specifies a minimum bound (if any) for that
variable.</td>
</tr>
<tr>
<td valign="top"><em class="code">variables(:,4)</em></td>
<td valign="top">Specifies a maximum bound (if any) for that
variable.</td>
</tr>
<tr>
<td valign="top"><em class="code">variables(:,5)</em></td>
<td valign="top">Specifies what file the variable should come from.
The only valid possibilies are "restart", "secondary", or
"calculate". "restart" will read from whatever file is specified by
<em class="code">tiegcm_restart_file_name</em>. "secondary" will
read from whatever file is specified by <em class=
"code">tiegcm_secondary_file_name</em>. "calculate" will call a
variable-dependent function -- see <em class=
"file">model_mod.f90</em>:<em class=
"code">tiegcm_to_dart_vector()</em> for the <em class=
"code">create_vtec()</em> example.</td>
</tr>
<tr>
<td valign="top"><em class="code">variables(:,6)</em></td>
<td valign="top">Specifies if the variable should be updated in the
TIEGCM restart file. The value may be "UPDATE" or anything else. If
<strong>and only if</strong> the variable comes from the restart
file <strong>and</strong> <em class="code">variables(:,6)</em> ==
"UPDATE" will the variable be modified in the TIEGCM restart file.
No variables in the secondary file are EVER modified.</td>
</tr>
</table>
</td>
</tr>
</tbody>
</table>
</div>
<!--==================================================================-->
<a name="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
adaptive_inflate_mod.f90
assim_model_mod.f90
assim_tools_mod.f90
types_mod.f90
cov_cutoff_mod.f90
ensemble_manager_mod.f90
filter.f90
location/threed_sphere/location_mod.f90
[null_,]mpi_utilities_mod.f90
obs_def_mod.f90
obs_kind_mod.f90
obs_model_mod.f90
obs_sequence_mod.f90
random_seq_mod.f90
reg_factor_mod.f90
smoother_mod.f90
sort_mod.f90
time_manager_mod.f90
utilities_mod.f90
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES - REQUIRED</h2>
<table>
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
<h2>PUBLIC INTERFACES - OPTIONAL</h2>
<table>
<tr>
<td><em class="call">use model_mod, only :</em></td>
<td><a href="#tiegcm_to_dart_vector">tiegcm_to_dart_vector</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#dart_vector_to_tiegcm">dart_vector_to_tiegcm</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_f107_value">get_f107_value</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#test_interpolate">test_interpolate</a></td>
</tr>
</table>
<p>A namelist interface <a href="#Namelist"><em class=
"code">&amp;model_nml</em></a> is defined by the module, and is
read from file <em class="file">input.nml</em>.</p>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_model_size"></a><br>
<div class="routine"><em class="call">model_size = get_model_size(
)</em>
<pre>
integer :: <em class="code">get_model_size</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the length of the model state vector. Required.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">model_size</em></td>
<td>The length of the model state vector.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="adv_1step"></a><br>
<div class="routine"><em class="call">call adv_1step(x, time)</em>
<pre>
real(r8), dimension(:), intent(inout) :: <em class="code">x</em>
type(time_type),        intent(in)    :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Since TIEGCM is not called as a subroutine, this is a NULL
interface. TIEGCM is advanced as a separate executable - i.e.
<em class="code">async == 2</em>. <em href="program">adv_1step</em>
only gets called if <em class="code">async == 0</em>. The
subroutine must still exist, but contains no code and will not be
called. An error message is issued if an unsupported value of
<em class="file">filter,perfect_model_obs</em>:<em class=
"code">async</em> is used.</p>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_state_meta_data"></a><br>
<div class="routine"><em class="call">call get_state_meta_data
(index_in, location, <em class=
"optionalcode">[, var_kind]</em> )</em>
<pre>
integer,             intent(in)  :: <em class="code">index_in</em>
type(location_type), intent(out) :: <em class="code">location</em>
integer, optional,   intent(out) :: <em class=
"optionalcode"> var_kind </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an integer index into the state vector structure, returns
the associated location. A second intent(out) optional argument
returns the generic kind of this item, e.g.
QTY_MOLEC_OXYGEN_MIXING_RATIO, QTY_ELECTRON_DENSITY, ... This
interface is required to be functional for all applications.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">index_in</em></td>
<td>Index of state vector element about which information is
requested.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>The location of state variable element.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">var_kind</em></td>
<td>The generic kind of the state variable element.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="model_interpolate"></a><br>
<div class="routine"><em class="call">call model_interpolate(x,
location, ikind, obs_val, istatus)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">x</em>
type(location_type),    intent(in)  :: <em class=
"code">location</em>
integer,                intent(in)  :: <em class="code">ikind</em>
real(r8),               intent(out) :: <em class=
"code">obs_val</em>
integer,                intent(out) :: <em class=
"code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a state vector, a location, and a model state variable
kind interpolates the state variable field to that location and
returns the value in obs_val. The istatus variable should be
returned as 0 unless there is some problem in computing the
interpolation in which case a positive value should be returned.
The ikind variable is one of the KIND parameters defined in the
<a href=
"../../assimilation_code/modules/observations/obs_kind_mod.html">obs_kind_mod.f90</a>
file and defines which generic kind of item is being
interpolated.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Location to which to interpolate.</td>
</tr>
<tr>
<td valign="top"><em class="code">itype</em></td>
<td>Kind of state field to be interpolated.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_val</em></td>
<td>The interpolated value from the model.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Integer value returning 0 for success. Other values can be
defined for various failures.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_model_time_step"></a><br>
<div class="routine"><em class="call">var =
get_model_time_step()</em>
<pre>
type(time_type) :: <em class="code">get_model_time_step</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the smallest useful forecast length (time step) of the
model. This is set by <em class="file">input.nml</em>:<em class=
"code">assimilation_period_seconds</em> and must be an integer
number of TIEGCM dynamical timesteps (as specified by <em class=
"file">tiegcm.nml</em>:<em class="code">STEP</em>) AND be able to
be expressed by <em class="file">tiegcm.nml</em>:<em class=
"code">STOP</em>. Since <em class="code">STOP</em> has three
components: day-of-year, hour, and minute, the <em class=
"code">assimilation_period_seconds</em> must be an integer number
of minutes.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var</em></td>
<td>Smallest forecast step of model.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="static_init_model"></a><br>
<div class="routine"><em class="call">call
static_init_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Called to do one-time initialization of the model. There are no
input arguments. <em class="code">static_init_model</em> reads the
DART and TIEGCM namelists and reads the grid geometry and
constructs the shape of the DART vector given the TIEGCM variables
specified in the DART namelist.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="end_model"></a><br>
<div class="routine"><em class="call">call end_model()</em></div>
<div class="indent1"><!-- Description -->
<p>Does all required shutdown and clean-up needed.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="init_time"></a><br>
<div class="routine"><em class="call">call init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a NULL INTERFACE for TIEGCM. If <em class=
"file">input.nml</em>:<em class="code">start_from_restart ==
.FALSE.</em>, this routine is called and will generate a fatal
error.</p>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="init_conditions"></a><br>
<div class="routine"><em class="call">call init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class="code">x</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a NULL INTERFACE for TIEGCM. If <em class=
"file">input.nml</em>:<em class="code">start_from_restart ==
.FALSE.</em>, this routine is called and will generate a fatal
error.</p>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="nc_write_model_atts"></a><br>
<div class="routine"><em class="call">ierr =
nc_write_model_atts(ncFileID)</em>
<pre>
integer             :: <em class="code">nc_write_model_atts</em>
integer, intent(in) :: <em class="code">ncFileID</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This routine writes the model-specific attributes to a netCDF
file. This includes the coordinate variables and any metadata, but
NOT the model state vector. We do have to allocate SPACE for the
model state vector, but that variable gets filled as the model
advances. If <em class="file">input.nml</em>:<em class=
"code">model_nml:output_state_vector == .TRUE.</em>, the DART state
vector is written as one long vector. If <em class=
"file">input.nml</em>:<em class=
"code">model_nml:output_state_vector == .FALSE.</em>, the DART
state vector is reshaped into the original TIEGCM variables and
those variables are written.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>Integer file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returns a 0 for successful completion.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="nc_write_model_vars"></a><br>
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
<p>This routine writes the DART state vector to a netCDF file. If
<em class="file">input.nml</em>:<em class=
"code">model_nml:output_state_vector == .TRUE.</em>, the DART state
vector is written as one long vector. If <em class=
"file">input.nml</em>:<em class=
"code">model_nml:output_state_vector == .FALSE.</em>, the DART
state vector is reshaped into the original TIEGCM variables and
those variables are written.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">statevec</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">copyindex</em></td>
<td>Integer index of copy to be written.</td>
</tr>
<tr>
<td valign="top"><em class="code">timeindex</em></td>
<td>The timestep counter for the given state.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returns 0 for normal completion.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="pert_model_state"></a><br>
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
<p><em class="code">pert_model_state</em> is intended to take a
single model state vector and perturbs it in some way to generate
initial conditions for spinning up ensembles. TIEGCM does this is a
manner that is different than most other models. The F10_7
parameter must be included in the DART state vector as a
QTY_1D_PARAMETER and gaussian noise is added to it. That value must
be conveyed to the tiegcm namelist and used to advance the
model.<br>
<br>
Most other models simply add noise with certain characteristics to
the model state.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state</em></td>
<td>State vector to be perturbed.</td>
</tr>
<tr>
<td valign="top"><em class="code">pert_state</em></td>
<td>Perturbed state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">interf_provided</em></td>
<td>This is returned as .TRUE. since the routine exists. A value of
.FALSE. would indicate that the default DART routine should just
add noise to every element of state.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_close_maxdist_init"></a><br>
<div class="routine"><em class="call">call
get_close_maxdist_init(gc, maxdist)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
real(r8),             intent(in)    :: <em class=
"code">maxdist</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a PASS-THROUGH routine, the actual routine is the
default one in <em class="file">location_mod</em>. In distance
computations any two locations closer than the given <em class=
"code">maxdist</em> will be considered close by the <em class=
"code">get_close_obs()</em> routine. <em class=
"code">get_close_maxdist_init</em> is listed on the <em class=
"code">use</em> line for the locations_mod, and in the public list
for this module, but has no subroutine declaration and no other
code in this module.</p>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_close_obs_init"></a><br>
<div class="routine"><em class="call">call get_close_obs_init(gc,
num, obs)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
integer,              intent(in)    :: <em class="code">num</em>
type(location_type),  intent(in)    :: <em class=
"code">obs(num)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a PASS-THROUGH routine. The default routine in the
location module precomputes information to accelerate the distance
computations done by <em class="code">get_close_obs()</em>. Like
the other PASS-THROUGH ROUTINES it is listed on the use line for
the locations_mod, and in the public list for this module, but has
no subroutine declaration and no other code in this module:</p>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_close_obs"></a><br>
<div class="routine"><em class="call">call get_close_obs(gc,
base_obs_loc, base_obs_kind, obs_loc, obs_kind, num_close,
close_ind <em class="optionalcode">[, dist]</em>)</em>
<pre>
type(get_close_type), intent(in)  :: <em class="code">gc</em>
type(location_type),  intent(in)  :: <em class=
"code">base_obs_loc</em>
integer,              intent(in)  :: <em class=
"code">base_obs_kind</em>
type(location_type),  intent(in)  :: <em class=
"code">obs_loc(:)</em>
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
<p>Given a location and kind, compute the distances to all other
locations in the <em class="code">obs_loc</em> list. The return
values are the number of items which are within maxdist of the
base, the index numbers in the original obs_loc list, and
optionally the distances. The <em class="code">gc</em> contains
precomputed information to speed the computations.<br>
<br>
This is different than the default <em class=
"code">location_mod:get_close_obs()</em> in that it is possible to
modify the 'distance' based on the DART 'kind'. This allows one to
apply specialized localizations.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>The get_close_type which stores precomputed information about
the locations to speed up searching</td>
</tr>
<tr>
<td valign="top"><em class="code">base_obs_loc</em></td>
<td>Reference location. The distances will be computed between this
location and every other location in the obs list</td>
</tr>
<tr>
<td valign="top"><em class="code">base_obs_kind</em></td>
<td>The kind of base_obs_loc</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_loc</em></td>
<td>Compute the distance between the base_obs_loc and each of the
locations in this list</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_kind</em></td>
<td>The corresponding kind of each item in the obs list</td>
</tr>
<tr>
<td valign="top"><em class="code">num_close</em></td>
<td>The number of items from the obs_loc list which are within
maxdist of the base location</td>
</tr>
<tr>
<td valign="top"><em class="code">close_ind</em></td>
<td>The list of index numbers from the obs_loc list which are
within maxdist of the base location</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">dist</em></td>
<td>If present, return the distance between each entry in the
close_ind list and the base location. If not present, all items in
the obs_loc list which are closer than maxdist will be added to the
list but the overhead of computing the exact distances will be
skipped.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="ens_mean_for_model"></a><br>
<div class="routine"><em class="call">call
ens_mean_for_model(ens_mean)</em>
<pre>
real(r8), dimension(:), intent(in) :: <em class=
"code">ens_mean</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>A model-size vector with the means of the ensembles for each of
the state vector items. The model should save a local copy of this
data if it needs to use it later to compute distances or other
values. This routine is called after each model advance and
contains the updated means.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_mean</em></td>
<td>State vector containing the ensemble mean.</td>
</tr>
</table>
</div>
<!--==================================================================-->
<h3>TIEGCM public routines</h3>
<!--==================================================================-->
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="tiegcm_to_dart_vector"></a><br>
<div class="routine"><em class="call">call
tiegcm_to_dart_vector(statevec, model_time)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class=
"code">statevec</em>
type(time_type),        intent(out) :: <em class=
"code">model_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Read TIEGCM fields from the TIEGCM restart file and/or TIEGCM
secondary file and pack them into a DART vector.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">statevec</em></td>
<td>variable that contains the DART state vector</td>
</tr>
<tr>
<td valign="top"><em class="code">model_time</em></td>
<td>variable that contains the LAST TIME in the TIEGCM restart
file.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="dart_vector_to_tiegcm"></a><br>
<div class="routine"><em class="call">call
dart_vector_to_tiegcm(statevec, dart_time)</em>
<pre>
real(r8), dimension(:), intent(in) :: <em class=
"code">statevec</em>
type(time_type),        intent(in) :: <em class=
"code">dart_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Unpacks a DART vector and updates the TIEGCM restart file
variables. Only those variables designated as 'UPDATE' are put into
the TIEGCM restart file. All variables are written to the DART
diagnostic files <strong>prior</strong> to the application of any
"clamping". The variables <strong>are "clamped"</strong> before
being written to the TIEGCM restart file. The clamping limits are
specified in columns 3 and 4 of <em class=
"code">&amp;model_nml:variables</em>.<br>
<br>
The time of the DART state is compared to the time in the restart
file to ensure that we are not improperly updating a restart
file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">statevec</em></td>
<td>Variable containing the DART state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">dart_time</em></td>
<td>Variable containing the time of the DART state vector.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_f107_value"></a><br>
<div class="routine"><em class="call">var = get_f107_value(x)</em>
<pre>
real(r8)                           :: <em class=
"code">get_f107_value</em>
real(r8), dimension(:), intent(in) :: <em class="code">x</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>If the F10_7 value is part of the DART state, return that value.
If it is not part of the DART state, just return the F10_7 value
from the TIEGCM namelist.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>Variable containing the DART state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">var</em></td>
<td>The f10_7 value.</td>
</tr>
</table>
</div>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="test_interpolate"></a><br>
<div class="routine"><em class="call">call test_interpolate(x,
locarray)</em>
<pre>
real(r8), dimension(:), intent(in) :: <em class="code">x</em>
real(r8), dimension(3), intent(in) :: <em class=
"code">locarray</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This function is <strong>only</strong> used by <a href=
"../../assimilation_code/programs/model_mod_check/model_mod_check.html%20models/POP/model_mod_check.html">
model_mod_check.f90</a> and can be modified to suit your needs.
<em class="program">test_interpolate()</em> exercises <em class=
"program">model_interpolate()</em>, <em class=
"program">get_state_meta_data()</em>, <em class=
"program">static_init_model()</em> and a host of supporting
routines.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>variable containing the DART state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">locarray</em></td>
<td>variable containing the location of interest.<br>
locarray(1) is the longitude (in degrees East)<br>
locarray(2) is the latitude (in degrees North)<br>
locarray(3) is the height (in meters).</td>
</tr>
</table>
</div>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<table width="100%" border="0" cellspacing="10" cellpadding="3"
summary="">
<tbody valign="top">
<tr>
<th align="left"><em class="file">filename</em></th>
<th align="left">purpose</th>
</tr>
<tr>
<td><em class="file">tiegcm.nml</em></td>
<td>TIEGCM control file modified to control starting and
stopping.</td>
</tr>
<tr>
<td><em class="file">input.nml</em></td>
<td>to read the model_mod namelist</td>
</tr>
<tr>
<td><em class="file">tiegcm_restart_p.nc</em></td>
<td>both read and modified by the TIEGCM model_mod</td>
</tr>
<tr>
<td><em class="file">tiegcm_s.nc</em></td>
<td>read by the GCOM model_mod for metadata purposes.</td>
</tr>
<tr>
<td><em class="file">namelist_update</em></td>
<td>DART file containing information useful for starting and
stopping TIEGCM. <em class="program">advance_model.csh</em> uses
this to update the TIEGCM file <em class=
"program">tiegcm.nml</em></td>
</tr>
<tr>
<td><em class="file">dart_log.out</em></td>
<td>the run-time diagnostic output</td>
</tr>
<tr>
<td><em class="file">dart_log.nml</em></td>
<td>the record of all the namelists (and their values) actually
USED</td>
</tr>
<tr>
<td><em class="file">log_advance.<i>nnnn</i>.txt</em></td>
<td>the run-time output of everything that happens in <em class=
"program">advance_model.csh</em>. This file will be in the
<em class="file">advance_temp<i>nnnn</i></em> directory.</td>
</tr>
</tbody>
</table>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>Matsuo, T., and E. A. Araujo-Pradere (2011),<br>
Role of thermosphere-ionosphere coupling in a global ionosphere
specification,<br>
<i>Radio Science</i>, <b>46</b>, RS0D23, <a href=
"http://dx.doi.org/doi:10.1029/2010RS004576">doi:10.1029/2010RS004576</a></li>
<li style="list-style: none"><!-- whitespace --><br></li>
<li>Lee, I. T., T, Matsuo, A. D. Richmond, J. Y. Liu, W. Wang, C.
H. Lin, J. L. Anderson, and M. Q. Chen (2012),<br>
Assimilation of FORMOSAT-3/COSMIC electron density profiles into
thermosphere/Ionosphere coupling model by using ensemble Kalman
filter,<br>
<i>Journal of Geophysical Research</i>, <b>117</b>, A10318,
<a href="http://dx.doi.org/doi:10.1029/2012JA017700">doi:10.1029/2012JA017700</a></li>
<li style="list-style: none"><!-- whitespace --><br></li>
<li>Matsuo, T., I. T. Lee, and J. L. Anderson (2013),<br>
Thermospheric mass density specification using an ensemble Kalman
filter,<br>
<i>Journal of Geophysical Research</i>, <b>118</b>, 1339-1350,
<a href=
"http://dx.doi.org/doi:10.1002/jgra.50162">doi:10.1002/jgra.50162</a></li>
<li style="list-style: none"><!-- whitespace --><br></li>
<li>Lee, I. T., H. F. Tsai, J. Y. Liu, Matsuo, T., and L. C. Chang
(2013),<br>
Modeling impact of FORMOSAT-7/COSMIC-2 mission on ionospheric space
weather monitoring,<br>
<i>Journal of Geophysical Research</i>, <b>118</b>, 6518-6523,
<a href=
"http://dx.doi.org/doi:10.1002/jgra.50538">doi:10.1002/jgra.50538</a></li>
<li style="list-style: none"><!-- whitespace --><br></li>
<li>Matsuo, T. (2014),<br>
Upper atmosphere data assimilation with an ensemble Kalman filter,
in Modeling the Ionosphere-Thermosphere System,<br>
<i>Geophys. Monogr. Ser.</i>, vol. 201, edited by J. Huba, R.
Schunk, and G. Khazanov, pp. 273-282, John Wiley &amp; Sons, Ltd,
Chichester, UK, <a href=
"http://dx.doi.org/doi:10.1002/9781118704417">doi:10.1002/9781118704417</a></li>
<li style="list-style: none"><!-- whitespace --><br></li>
<li>Hsu, C.-H., T. Matsuo, W. Wang, and J. Y. Liu (2014),<br>
Effects of inferring unobserved thermospheric and ionospheric state
variables by using an ensemble Kalman filter on global ionospheric
specification and forecasting,<br>
<i>Journal of Geophysical Research</i>, <b>119</b>, 9256-9267,
<a href=
"http://dx.doi.org/doi:10.1002/2014JA020390">doi:10.1002/2014JA020390</a></li>
<li style="list-style: none"><!-- whitespace --><br></li>
<li>Chartier, A., T. Matsuo, J. L. Anderson, G. Lu, T. Hoar, N.
Collins, A. Coster, C. Mitchell, L. Paxton, G. Bust (2015),<br>
Ionospheric Data Assimilation and Forecasting During Storms,<br>
<i>Journal of Geophysical Research</i>, under review</li>
<li style="list-style: none; display: inline">
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<ul>
<li>Models are free to issue calls to the error handler as they see
fit. No standard error handler calls are mandated.</li>
</ul>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>It is likely that a number of additional optional interfaces
will be added to the model_mod structure. For instance, hints about
how to divide the state vector into regions for parallel
assimilation will need to be obtained from the model. It is planned
that the interf_provided mechanism used in pert_model_state will
allow those who do not wish to support enhanced interfaces to add
NULL interfaces by simply pasting in an interface block.</p>

<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================--></li>
</ul>
</body>
</html>
