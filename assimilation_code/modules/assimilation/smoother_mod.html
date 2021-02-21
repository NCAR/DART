<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module smoother_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE smoother_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Implements a fixed lag ensemble smoother as part of the filter.
For now, this is done inefficiently with a separate call to
<em class="program">assim_tools_mod:filter_assim()</em> for each
lag.<br>
<br>
To enable the smoother, set the number of lags (num_lags) to
something larger than 0 in the <em class=
"program">smoother_nml</em> section of your <em class=
"file">input.nml</em> file and run <em class="program">filter</em>
as before.</p>
<div class="routine">
<pre>
&amp;smoother_nml
   num_lags              = <em class="changed">10</em>,
   start_from_restart    = .false.,
   output_restart        = .true.,
   restart_in_file_name  = "ics",
   restart_out_file_name = "restart"  /
</pre></div>
<p>In the low order models, 10 is a plausible number.<br>
<br>
In addition to generating <em class="file">preassim.nc</em> and
<em class="file">analysis.nc</em> files, files of the form
<em class="file">Lag_NNNNN_Diag.nc</em> will be generated. Each of
these has N fewer timesteps than the lag=0 run, starting at the
same time but ending N timesteps sooner. The <em class=
"file">obs_seq.final</em> file and the <em class=
"file">preassim.nc</em> and <em class="file">analysis.nc</em> files
will be the same as the non-lagged version; the new output will be
in each of the <em class="file">Lag_NNNNN_Diag.nc</em> files.</p>
<a name="Example" id="Example"></a>
<h2>EXAMPLE</h2>
<p>If you have a <em class="file">true_state.nc</em> file and want
to use the <em class="program">plot_total_err</em> matlab function
to plot the error, you must do the following steps to generate
analogs of lagged <em class="file">true_state.nc</em> files to use
as a comparison. (The logic is not currently implemented in the
matlab scripts to be able to compare netCDF files with unequal time
coordinates.)<br>
<br>
Make N separate versions of the true_state.nc with the last N
timesteps removed. Using the netCDF NCO operator program 'ncks' is
one way. If the true_state.nc file has 1000 time steps, then this
command removes the last one:</p>
<div class="unix">ncks -d time,0,998 true_state.nc
True_Lag01.nc</div>
<p>Note that the first time is at index 0, so the last timestep is
index 999 in the full file, and 998 in the truncated file. Repeat
this step for all N lags. Here are NCO commands to generate 10
truth files for num_lags = 10, 1000 time steps in
true_state.nc:</p>
<div class="unix">ncks -d time,0,998 true_state.nc
True_Lag01.nc<br>
ncks -d time,0,997 true_state.nc True_Lag02.nc<br>
ncks -d time,0,996 true_state.nc True_Lag03.nc<br>
ncks -d time,0,995 true_state.nc True_Lag04.nc<br>
ncks -d time,0,994 true_state.nc True_Lag05.nc<br>
ncks -d time,0,993 true_state.nc True_Lag06.nc<br>
ncks -d time,0,992 true_state.nc True_Lag07.nc<br>
ncks -d time,0,991 true_state.nc True_Lag08.nc<br>
ncks -d time,0,990 true_state.nc True_Lag09.nc<br>
ncks -d time,0,989 true_state.nc True_Lag10.nc<br></div>
<p>Here is an example matlab session which plots the lag=0 results
and then odd numbered lags from 1 to 9. It uses the <em class=
"program">plot_total_err</em> function from the $DART/matlab
directory:</p>
<pre>
datadir    = '.';
truth_file = fullfile(datadir,'true_state.nc');
diagn_file = fullfile(datadir,'preassim.nc');
plot_total_err
reply = input('original data.  hit enter to continue ');

truth_file = fullfile(datadir,'True_Lag01.nc');
diagn_file = fullfile(datadir,'Lag_00001_Diag.nc');
plot_total_err
reply = input('Lag 01.  hit enter to continue ');

truth_file = fullfile(datadir,'True_Lag03.nc');
diagn_file = fullfile(datadir,'Lag_00003_Diag.nc');
plot_total_err
reply = input('Lag 03.  hit enter to continue ');

truth_file = fullfile(datadir,'True_Lag05.nc');
diagn_file = fullfile(datadir,'Lag_00005_Diag.nc');
plot_total_err
reply = input('Lag 05.  hit enter to continue ');

truth_file = fullfile(datadir,'True_Lag07.nc');
diagn_file = fullfile(datadir,'Lag_00007_Diag.nc');
plot_total_err
reply = input('Lag 07.  hit enter to continue ');

truth_file = fullfile(datadir,'True_Lag09.nc');
diagn_file = fullfile(datadir,'Lag_00009_Diag.nc');
plot_total_err
reply = input('Lag 09.  hit enter to continue ');
</pre>
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
&amp;smoother_nml
   num_lags              = 0,
   start_from_restart    = .false.,
   output_restart        = .false.,
   restart_in_file_name  = 'ics',
   restart_out_file_name = 'restart'  
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
<td>num_lags</td>
<td>integer</td>
<td>Number of smoother lags; &lt; 1 means no smoother.</td>
</tr>
<tr>
<td>start_from_restart</td>
<td>logical</td>
<td>True if smoother states are to come from restart file(s). False
if they are to be spun up from scratch.</td>
</tr>
<tr>
<td>output_restart</td>
<td>logical</td>
<td>True if restart file(s) are to be written, else false.</td>
</tr>
<tr>
<td>restart_in_file_name</td>
<td>character(len=129)</td>
<td>String used to construct the file name from which to read
restart data. <code>Lag_NNNNN_</code> will be prepended to the specified value
to create the actual filename. If each ensemble is to be read from
a separate file, the .NNNN ensemble number will also be appended.
e.g. specifying 'ics' here results in 'Lag_00001_ics' if all
ensemble members are read from a single file, 'Lag_00001_ics.0001',
'Lag_00001_ics.0002', etc for multiples.</td>
</tr>
<tr>
<td>restart_out_file_name   </td>
<td>character(len=129)</td>
<td>String used to construct the file name to which to write
restart data. <code>Lag_NNNNN_</code> will be prepended to the specified value
to create the actual filename. If each ensemble is to be written to
a separate file, the .NNNN ensemble number will also be appended.
e.g. specifying 'restart' here results in 'Lag_00001_restart' if
all ensemble members are written to a single file,
'Lag_00001_restart.0001', 'Lag_00001_restart.0002', etc for
multiples.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
mpi_utilities_mod
utilities_mod
ensemble_manager_mod
time_manager_mod
assim_model_mod
assim_tools_mod
obs_sequence_mod
adaptive_inflate_mod
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
<td><em class="call">use smoother_mod, only :</em></td>
<td><a href="#smoother_read_restart">smoother_read_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#advance_smoother">advance_smoother</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#smoother_gen_copy_meta_data">smoother_gen_copy_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#smoother_write_restart">smoother_write_restart</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_smoother">init_smoother</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#do_smoothing">do_smoothing</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#smoother_mean_spread">smoother_mean_spread</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#smoother_assim">smoother_assim</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#filter_state_space_diagnostics">filter_state_space_diagnostics</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#smoother_ss_diagnostics">smoother_ss_diagnostics</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#smoother_end">smoother_end</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="smoother_read_restart" id="smoother_read_restart"></a><br>
<div class="routine"><em class="call">call
smoother_read_restart(ens_handle, ens_size, model_size, time1,
init_time_days)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer, intent(in)                :: <em class=
"code">ens_size</em>
integer, intent(in)                :: <em class=
"code">model_size</em>
type(time_type), intent(inout)     :: <em class="code">time1</em>
integer, intent(in)                :: <em class=
"code">init_time_days</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads in ensemble of states for all lag estimates from a restart
file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle  </em></td>
<td>Handle of ensemble manager structure of single state; copied
into all lags for startup.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size  </em></td>
<td>Size of the ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">model_size  </em></td>
<td>Size of the model state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">time1  </em></td>
<td>Overwrite the time in the restart file with this value if
init_time_days is non-negative.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">init_time_days  </em></td>
<td>If non-negative, use time1 instead of time in restart
file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="advance_smoother" id="advance_smoother"></a><br>
<div class="routine"><em class="call">call
advance_smoother(ens_handle)</em>
<pre>
type(ensemble_type), intent(in) :: <em class="code">ens_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Advances smoother state estimates at all lags forward in time.
This entails copying the most recent smoother state, contained in
ens_handle, into the lag 1 smoother state and pushing back all
other lags by 1 (i.e. lag 1 becomes lag 2, etc.).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle  </em></td>
<td>Ensemble handle with most recent filtered state.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="smoother_gen_copy_meta_data" id=
"smoother_gen_copy_meta_data"></a><br>
<div class="routine"><em class="call">call
smoother_gen_copy_meta_data(num_output_state_members,
output_inflation)</em>
<pre>
integer, intent(in) :: <em class=
"code">num_output_state_members</em>
logical, intent(in) :: <em class="code">output_inflation</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initializes the metadata required for the smoother state space
diagnostic files.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">num_output_state_members  </em></td>
<td>Number of copies of smoother state vector that should be in
state space diagnostic output.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_inflation  </em></td>
<td>True if smoother state space output should include inflation
values.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="smoother_write_restart" id=
"smoother_write_restart"></a><br>
<div class="routine"><em class="call">call
smoother_write_restart(start_copy, end_copy)</em>
<pre>
integer, intent(in) :: <em class="code">start_copy</em>
integer, intent(in) :: <em class="code">end_copy</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Outputs restart files for all lags of smoother state. Integer
arguments specify the start and end global indices of a continguous
set of copies that contain the ensemble members.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">start_copy  </em></td>
<td>Global index of ensemble copy that starts the actual ensemble
members for smoother.</td>
</tr>
<tr>
<td valign="top"><em class="code">end_copy  </em></td>
<td>Global index of ensemble copy that ends the actual ensemble
members for smoother.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_smoother" id="init_smoother"></a><br>
<div class="routine"><em class="call">call
init_smoother(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=
"code">ens_handle</em>
integer, intent(in)                :: <em class=
"code">POST_INF_COPY</em>
integer, intent(in)                :: <em class=
"code">POST_INF_SD_COPY</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initializes the storage needed for a smoother. Also initializes
an adaptive inflation type that does NO inflation (not currently
supported for smoothers).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle  </em></td>
<td>An ensemble handle for the filter that contains information
about ensemble and model size.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">POST_INF_COPY  </em></td>
<td>Global index of ensemble copy that holds posterior state space
inflation values.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">POST_INF_SD_COPY  </em></td>
<td>Global index of ensemble copy that holds posterior inflation
standard deviation values.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="do_smoothing" id="do_smoothing"></a><br>
<div class="routine"><em class="call">var = do_smoothing()</em>
<pre>
logical, intent(out) :: <em class="code">do_smoothing</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if smoothing is to be done, else false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">do_smoothing  </em></td>
<td>Returns true if smoothing is to be done.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="smoother_mean_spread" id="smoother_mean_spread"></a><br>
<div class="routine"><em class="call">call
smoother_mean_spread(ens_size,ENS_MEAN_COPY,ENS_SD_COPY,
output_state_ens_mean,output_state_ens_spread)</em>
<pre>
integer, intent(in) :: <em class="code">ens_size</em>
integer, intent(in) :: <em class="code">ENS_MEAN_COPY</em>
integer, intent(in) :: <em class="code">ENS_SD_COPY</em>
logical, intent(in) :: <em class="code">output_state_ens_mean</em>
logical, intent(in) :: <em class=
"code">output_state_ens_spread</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Computes the ensemble mean (and spread if required) of all state
variables for all lagged ensembles. Spread is only computed if it
is required for output.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_size  </em></td>
<td>Size of ensemble.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ENS_MEAN_COPY  </em></td>
<td>Global index of copy that stores ensemble mean.</td>
</tr>
<tr>
<td valign="top"><em class="code">ENS_SD_COPY  </em></td>
<td>Global index of copy that stores ensemble spread.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_state_ens_mean  </em></td>
<td>True if the ensemble mean is to be output to state diagnostic
file.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_state_ens_spread  </em></td>
<td>True if ensemble spread is to be output to state diagnostic
file.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="smoother_assim" id="smoother_assim"></a><br>
<div class="routine"><em class="call">call
smoother_assim(obs_ens_handle, seq, keys, ens_size, num_groups,
obs_val_index, ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY,
PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY,
OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START,
OBS_PRIOR_VAR_END)</em>
<pre>
type(ensemble_type), intent(inout)  :: <em class=
"code">obs_ens_handle</em>
type(obs_sequence_type), intent(in) :: <em class="code">seq</em>
integer, dimension(:), intent(in)   :: <em class="code">keys</em>
integer, intent(in)                 :: <em class=
"code">ens_size</em>
integer, intent(in)                 :: <em class=
"code">num_groups</em>
integer, intent(in)                 :: <em class=
"code">obs_val_index</em>
integer, intent(in)                 :: <em class=
"code">ENS_MEAN_COPY</em>
integer, intent(in)                 :: <em class=
"code">ENS_SD_COPY</em>
integer, intent(in)                 :: <em class=
"code">PRIOR_INF_COPY</em>
integer, intent(in)                 :: <em class=
"code">PRIOR_INF_SD_COPY</em>
integer, intent(in)                 :: <em class=
"code">OBS_KEY_COPY</em>
integer, intent(in)                 :: <em class=
"code">OBS_GLOBAL_QC_COPY</em>
integer, intent(in)                 :: <em class=
"code">OBS_PRIOR_MEAN_START</em>
integer, intent(in)                 :: <em class=
"code">OBS_PRIOR_MEAN_END</em>
integer, intent(in)                 :: <em class=
"code">OBS_PRIOR_VAR_START</em>
integer, intent(in)                 :: <em class=
"code">OBS_PRIOR_VAR_END</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Does assimilation of a set of observations for each smoother
lag.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">obs_ens_handle  </em></td>
<td>Handle for ensemble manager holding prior estimates of
observations.</td>
</tr>
<tr>
<td valign="top"><em class="code">seq  </em></td>
<td>Observation sequence being assimilated.</td>
</tr>
<tr>
<td valign="top"><em class="code">keys  </em></td>
<td>A one dimensional array containing indices in seq of
observations to as similate at current time.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size  </em></td>
<td>Ensemble size.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_groups  </em></td>
<td>Number of groups in filter.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">obs_val_index  </em></td>
<td>Integer index of copy of data in seq that contains the observed
value from instruments.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ENS_MEAN_COPY  </em></td>
<td>Global index in smoother's state ensemble that holds ensemble
mean.</td>
</tr>
<tr>
<td valign="top"><em class="code">ENS_SD_COPY  </em></td>
<td>Global index in smoother's state ensemble that holds ensemble
standard deviation.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">PRIOR_INF_COPY  </em></td>
<td>Global index in obs_ens_handle that holds inflation values (not
used for smoother).</td>
</tr>
<tr>
<td valign="top"><em class=
"code">PRIOR_INF_SD_COPY  </em></td>
<td>Global index in obs_ens_handle that holds inflation sd values
(not used for smoother).</td>
</tr>
<tr>
<td valign="top"><em class=
"code">OBS_KEY_COPY  </em></td>
<td>Global index in obs_ens_handle that holds the key for the
observation.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">OBS_GLOBAL_QC_COPY  </em></td>
<td>Global index in obs_ens_handle that holds the quality control
value.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">OBS_PRIOR_MEAN_START  </em></td>
<td>Global index in obs_ens_handle that holds the first group's
prior mean.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">OBS_PRIOR_MEAN_END  </em></td>
<td>Global index in obs_ens_handle that holds the last group's
prior mean.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">OBS_PRIOR_VAR_START  </em></td>
<td>Global index in obs_ens_handle that holds the first group's
prior variance.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">OBS_PRIOR_VAR_END  </em></td>
<td>Global index in obs_ens_handle that holds the last group's
prior variance.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="filter_state_space_diagnostics" id=
"filter_state_space_diagnostics"></a><br>
<div class="routine"><em class="call">call
filter_state_space_diagnostics(out_unit, ens_handle, model_size,
num_output_state_members, output_state_mean_index,
output_state_spread_index, output_inflation, temp_ens,
ENS_MEAN_COPY, ENS_SD_COPY, inflate, INF_COPY, INF_SD_COPY)</em>
<pre>
type(netcdf_file_type), intent(inout)   :: <em class=
"code">out_unit</em>
type(ensemble_type), intent(inout)      :: <em class=
"code">ens_handle</em>
integer, intent(in)                     :: <em class=
"code">model_size</em>
integer, intent(in)                     :: <em class=
"code">num_output_state_members</em>
integer, intent(in)                     :: <em class=
"code">output_state_mean_index</em>
integer, intent(in)                     :: <em class=
"code">output_state_spread_index</em>
logical, intent(in)                     :: <em class=
"code">output_inflation</em>
real(r8), intent(out)                   :: <em class=
"code">temp_ens(model_size)</em>
integer, intent(in)                     :: <em class=
"code">ENS_MEAN_COPY</em>
integer, intent(in)                     :: <em class=
"code">ENS_SD_COPY</em>
type(adaptive_inflate_type), intent(in) :: <em class=
"code">inflate</em>
integer, intent(in)                     :: <em class=
"code">INF_COPY</em>
integer, intent(in)                     :: <em class=
"code">INF_SD_COPY</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Writes state space diagnostic values including ensemble members,
mean and spread, and inflation mean and spread to a netcdf
file.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">out_unit  </em></td>
<td>Descriptor for the netcdf file being written.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_handle  </em></td>
<td>Ensemble handle whose state space values are to be
written.</td>
</tr>
<tr>
<td valign="top"><em class="code">model_size  </em></td>
<td>Size of the model state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">num_output_state_members  </em></td>
<td>Number of individual state members to be output.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_state_mean_index  </em></td>
<td>Index in netcdf file for ensemble mean.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_state_spread_index  </em></td>
<td>Index in netcdf file for ensemble spread.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_inflation  </em></td>
<td>True if the inflation values are to be output. Default is
.TRUE.</td>
</tr>
<tr>
<td valign="top"><em class="code">temp_ens  </em></td>
<td>Storage passed in to avoid having to allocate extra space.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ENS_MEAN_COPY  </em></td>
<td>Global index in ens_handle for ensemble mean.</td>
</tr>
<tr>
<td valign="top"><em class="code">ENS_SD_COPY  </em></td>
<td>Global index in ens_handle for ensemble spread.</td>
</tr>
<tr>
<td valign="top"><em class="code">inflate  </em></td>
<td>Contains description and values of state space inflation.</td>
</tr>
<tr>
<td valign="top"><em class="code">INF_COPY  </em></td>
<td>Global index in ens_handle of inflation values.</td>
</tr>
<tr>
<td valign="top"><em class="code">INF_SD_COPY  </em></td>
<td>Global index in ens_handle of inflation standard deviation
values.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="smoother_ss_diagnostics" id=
"smoother_ss_diagnostics"></a><br>
<div class="routine"><em class="call">call
smoother_ss_diagnostics(model_size, num_output_state_members,
output_inflation, temp_ens, ENS_MEAN_COPY, ENS_SD_COPY,
POST_INF_COPY, POST_INF_SD_COPY)</em>
<pre>
integer, intent(in)   :: <em class="code">model_size</em>
integer, intent(in)   :: <em class=
"code">num_output_state_members</em>
logical, intent(in)   :: <em class="code">output_inflation</em>
real(r8), intent(out) :: <em class="code">temp_ens(model_size)</em>
integer, intent(in)   :: <em class="code">ENS_MEAN_COPY</em>
integer, intent(in)   :: <em class="code">ENS_SD_COPY</em>
integer, intent(in)   :: <em class="code">POST_INF_COPY</em>
integer, intent(in)   :: <em class="code">POST_INF_SD_COPY</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Outputs state space diagnostics files for all smoother lags.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">model_size  </em></td>
<td>Size of the model state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">num_output_state_members  </em></td>
<td>Number of state copies to be output in the state space
diagnostics file.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_inflation  </em></td>
<td>True if the inflation values are to be output. Default is
.TRUE.</td>
</tr>
<tr>
<td valign="top"><em class="code">temp_ens  </em></td>
<td>Storage passed in to avoid having to allocate extra space.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ENS_MEAN_COPY  </em></td>
<td>Global index of the ensemble mean in the lag smoother ensemble
handles.</td>
</tr>
<tr>
<td valign="top"><em class="code">ENS_SD_COPY  </em></td>
<td>Global index of the ensemble spread in the lag smoother
ensemble handles.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">POST_INF_COPY  </em></td>
<td>Global index of the inflation value in the lag smoother
ensemble handles (not currently used).</td>
</tr>
<tr>
<td valign="top"><em class=
"code">POST_INF_SD_COPY  </em></td>
<td>Global index of the inflation spread in the lag smoother
ensemble handles (not currently used).</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="smoother_end" id="smoother_end"></a><br>
<div class="routine"><em class="call">call
smoother_end()</em></div>
<div class="indent1"><!-- Description -->
<p>Releases storage allocated for smoother.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="smoother_inc_lags" id="smoother_inc_lags"></a><br>
<div class="routine"><em class="call">call
smoother_inc_lags()</em></div>
<div class="indent1"><!-- Description -->
<p>Increments the number of lags that are in use for smoother. Used
when a smoother is being started up and there have not been enough
times to propagate the state to all requested lags.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
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
<li>input.nml</li>
<li>smoother initial condition files</li>
<li>smoother restart files</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<!-- FIXME 
should put a reference to any of Shree's papers here -->
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
<td valign="top">smoother_gen_copy_meta_data</td>
<!-- message -->
<td valign="top">output metadata in smoother needs ensemble size
&lt; 10000, not ###</td>
<!-- comment -->
<td valign="top">Can't output more than 9999 copies.</td>
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
