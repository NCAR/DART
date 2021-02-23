<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module adaptive_inflate_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE adaptive_inflate_mod</h1>
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
<p>This module implements a variety of hierarchical Bayesian
adaptive inflation algorithms for use with ensemble filters. It can
provide constant valued inflation in state or observation space,
consistent with previous DART releases. It can provide
spatially-constant, time-varying adaptive inflation. It can provide
spatially-varying, time-varying adaptive inflation and it can
provide temporally-varying observation space inflation. And
finally, it can provide adaptive damped inflation, which decreases
inflation through time when observation density varies. Diagnostic
output and restart files are available. Several papers on the NCAR
<a href="http://www.image.ucar.edu/DAReS">IMAGe/DAReS</a> web page
document the algorithms in detail. The <em class=
"file">DART/tutorial/section12</em> chapter has more
information.</p>
<p>Details on controlling the inflation options are contained in
the documentation for the filter. The filter_nml controls what
inflation options are used.</p>
<p>Inflation flavor 3 (spatially-constant state space) reads and
writes a restart file that is the full size of the state vector,
however it takes the first value in the array and replicates that
throughout the array. This allows one to switch between flavors 2
and 3. Going from inflation flavor 3 to 2 the initial value for all
items in the state vector will be a constant value and will then
start to adapt. Going from inflation flavor 2 to 3 whatever value
is in the array at index 1 will be replicated and used for the
entire rest of the state vector items.</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
random_seq_mod
time_manager_mod
ensemble_manager_mod
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
<td><em class="call">use adaptive_inflate_mod, only :</em></td>
<td><a href="#update_inflation">update_inflation</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#adaptive_inflate_end">adaptive_inflate_end</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#inflate_ens">inflate_ens</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#output_inflate_diagnostics">output_inflate_diagnostics</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#do_obs_inflate">do_obs_inflate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#do_single_ss_inflate">do_single_ss_inflate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#do_varying_ss_inflate">do_varying_ss_inflate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#adaptive_inflate_init">adaptive_inflate_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#adaptive_inflate_type">adaptive_inflate_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_inflate">get_inflate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_inflate">set_inflate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_sd">set_sd</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_sd">set_sd</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#deterministic_inflate">deterministic_inflate</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="update_inflation" id="update_inflation"></a><br>
<div class="routine"><em class="call">call
update_inflation(inflate_handle, inflate, inflate_sd, prior_mean,
prior_var, obs, obs_var, gamma)</em>
<pre>
type(adaptive_inflate_type), intent(in)    :: <em class=
"code">inflate_handle</em>
real(r8),                    intent(inout) :: <em class=
"code">inflate</em>
real(r8),                    intent(inout) :: <em class=
"code">inflate_sd</em>
real(r8),                    intent(in)    :: <em class=
"code">prior_mean</em>
real(r8),                    intent(in)    :: <em class=
"code">prior_var</em>
real(r8),                    intent(in)    :: <em class=
"code">obs</em>
real(r8),                    intent(in)    :: <em class=
"code">obs_var</em>
real(r8),                    intent(in)    :: <em class=
"code">gamma</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Updates the mean and standard deviation of an inflation
distribution given the prior values, the prior observation ensemble
mean and variance, and the observation and its error variance. The
factor gamma is the expected impact (0 to 1) of the state variable
corresponding to the inflation on the observation and is the
product of the ensemble correlation plus an additional localization
factor or group regression factors.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle to object that describes the inflation type and
values.</td>
</tr>
<tr>
<td valign="top"><em class="code">inflate  </em></td>
<td>Prior mean value of the inflation distribution.</td>
</tr>
<tr>
<td valign="top"><em class="code">inflate_sd  </em></td>
<td>Prior standard deviation of the inflation distribution.</td>
</tr>
<tr>
<td valign="top"><em class="code">prior_mean  </em></td>
<td>Mean of the prior observation ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">prior_var  </em></td>
<td>Variance of the prior observation ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs  </em></td>
<td>The observed value.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_var  </em></td>
<td>Observational error variance.</td>
</tr>
<tr>
<td valign="top"><em class="code">gamma  </em></td>
<td>Expected impact factor, product of correlation, localization,
regression factor.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="adaptive_inflate_end" id="adaptive_inflate_end"></a><br>
<div class="routine"><em class="call">call
adaptive_inflate_end(inflate_handle, ens_handle, ss_inflate_index,
ss_inflate_sd_index)</em>
<pre>
type(adaptive_inflate_type), intent(in)    :: <em class=
"code">inflate_handle</em>
type(ensemble_type),         intent(inout) :: <em class=
"code">ens_handle</em>
integer,                     intent(in)    :: <em class=
"code">ss_inflate_index</em>
integer,                     intent(in)    :: <em class=
"code">ss_inflate_sd_index</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Outputs the values of inflation to restart files using the
ensemble_manager for state space inflation and file output for
observation space inflation. Releases allocated storage in
inflate_handle.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for the details of the inflation being performed.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_handle  </em></td>
<td>Handle for ensemble storage that holds values of state space
inflation.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ss_inflate_index  </em></td>
<td>Index in ensemble storage copies for state space
inflation.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ss_inflate_sd_index  </em></td>
<td>Index in ensemble storage copies for state space inflation
standard deviation.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="inflate_ens" id="inflate_ens"></a><br>
<div class="routine"><em class="call">call
inflate_ens(inflate_handle, ens,mean, inflate <em class=
"optionalcode">[,var_in]</em>)</em>
<pre>
type(adaptive_inflate_type),               intent(in)  :: <em class="code">inflate_handle</em>
real(r8),                    dimension(:), intent(out) :: <em class="code">ens</em>
real(r8),                                  intent(in)  :: <em class="code">mean</em>
real(r8),                                  intent(in)  :: <em class="code">inflate</em>
real(r8),                    optional,     intent(in)  :: <em class="optionalcode">var_in</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an ensemble, its mean and the covarance inflation factor,
inflates the ensemble.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for the details of the inflation being performed.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens  </em></td>
<td>Values for the ensemble to be inflated</td>
</tr>
<tr>
<td valign="top"><em class="code">mean  </em></td>
<td>The mean of the ensemble.</td>
</tr>
<tr>
<td valign="top"><em class="code">inflate  </em></td>
<td>The covariance inflation factor.</td>
</tr>
<tr>
<td valign="top"><em class="code">var_in  </em></td>
<td>The variance of the ensemble.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="output_inflate_diagnostics" id=
"output_inflate_diagnostics"></a><br>
<div class="routine"><em class="call">call
output_inflate_diagnostics(inflate_handle, time)</em>
<pre>
type(adaptive_inflate_type), intent(in) :: <em class=
"code">inflate_handle</em>
type(time_type),             intent(in) :: <em class=
"code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Outputs diagnostic record of inflation for the observation space
of spatially constant state space inflation. Spatially varying
state space diagnostics are in the Posterior and Prior Diagnostic
netcdf files and are written with calls from filter.f90.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for the details of the inflation being performed.</td>
</tr>
<tr>
<td valign="top"><em class="code">time  </em></td>
<td>Time of this diagnostic info.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="do_obs_inflate" id="do_obs_inflate"></a><br>
<div class="routine"><em class="call">var =
do_obs_inflate(inflate_handle)</em>
<pre>
logical,               intent(out) :: <em class=
"code">do_obs_inflate</em>
adaptive_inflate_type, intent(in)  :: <em class=
"code">inflate_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if observation space inflation is being done by
this handle.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">do_obs_inflate  </em></td>
<td>True if obs space inflation is being done by this handle.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle to inflation details.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="do_varying_ss_inflate" id=
"do_varying_ss_inflate"></a><br>
<div class="routine"><em class="call">var =
do_varying_ss_inflate(inflate_handle)</em>
<pre>
logical,               intent(out) :: <em class=
"code">do_varying_ss_inflate</em>
adaptive_inflate_type, intent(in)  :: <em class=
"code">inflate_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if spatially varying state space inflation is being
done by this handle.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">do_varying_ss_inflate  </em></td>
<td>True if spatially varying state space inflation is being done
by this handle.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle to inflation details.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="do_single_ss_inflate" id="do_single_ss_inflate"></a><br>
<div class="routine"><em class="call">var =
do_single_ss_inflate(inflate_handle)</em>
<pre>
logical,               intent(out) :: <em class=
"code">do_single_ss_inflate</em>
adaptive_inflate_type, intent(in)  :: <em class=
"code">inflate_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if spatially fixed state space inflation is being
done by this handle.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">do_single_ss_inflate  </em></td>
<td>True if spatially fixed state space inflation is being done by
this handle.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle to inflation details.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="adaptive_inflate_init" id=
"adaptive_inflate_init"></a><br>
<div class="routine"><em class="call">call
adaptive_inflate_init(inflate_handle, inf_flavor,
mean_from_restart, sd_from_restart, output_restart, deterministic,
in_file_name, out_file_name, diag_file_name, inf_initial,
sd_initial, inf_lower_bound, inf_upper_bound, sd_lower_bound,
ens_handle, ss_inflate_index, ss_inflate_sd_index, label)</em>
<pre>
type(adaptive_inflate_type), intent(inout) :: <em class=
"code">inflate_handle</em>
integer, intent(in)                        :: <em class=
"code">inf_flavor</em>
logical, intent(in)                        :: <em class=
"code">mean_from_restart</em>
logical, intent(in)                        :: <em class=
"code">sd_from_restart</em>
logical, intent(in)                        :: <em class=
"code">output_restart</em>
logical, intent(in)                        :: <em class=
"code">deterministic</em>
character(len=*), intent(in)               :: <em class=
"code">in_file_name</em>
character(len=*), intent(in)               :: <em class=
"code">out_file_name</em>
character(len=*), intent(in)               :: <em class=
"code">diag_file_name</em>
real(r8), intent(in)                       :: <em class=
"code">inf_initial</em>
real(r8), intent(in)                       :: <em class=
"code">sd_initial</em>
real(r8), intent(in)                       :: <em class=
"code">inf_lower_bound</em>
real(r8), intent(in)                       :: <em class=
"code">inf_upper_bound</em>
real(r8), intent(in)                       :: <em class=
"code">sd_lower_bound</em>
type(ensemble_type), intent(inout)         :: <em class=
"code">ens_handle</em>
integer, intent(in)                        :: <em class=
"code">ss_inflate_index</em>
integer, intent(in)                        :: <em class=
"code">ss_inflate_sd_index</em>
character(len=*), intent(in)               :: <em class=
"code">label</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initializes a descriptor of an inflation object.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for the inflation descriptor being initialized.</td>
</tr>
<tr>
<td valign="top"><em class="code">inf_flavor  </em></td>
<td>Type of inflation, 1=obs_inflate, 2=varying_ss_inflate,
3=single_ss_inflate.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">mean_from_restart  </em></td>
<td>True if inflation mean values to be read from restart
file.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">sd_from_restart  </em></td>
<td>True if inflation standard deviation values to be read from
restart file.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">output_restart  </em></td>
<td>True if an inflation restart file is to be output.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">deterministic  </em></td>
<td>True if deterministic inflation is to be done.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">in_file_name  </em></td>
<td>File name from which to read restart.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">out_file_name  </em></td>
<td>File name to which to write restart.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">diag_file_name  </em></td>
<td>File name to which to write diagnostic output; obs space
inflation only .</td>
</tr>
<tr>
<td valign="top"><em class="code">inf_initial  </em></td>
<td>Initial value of inflation for start_from_restart=.false.</td>
</tr>
<tr>
<td valign="top"><em class="code">sd_initial  </em></td>
<td>Initial value of inflation standard deviation for
start_from_restart=.false.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inf_lower_bound  </em></td>
<td>Lower bound on inflation value.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inf_upper_bound  </em></td>
<td>Upper bound on inflation value.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">sd_lower_bound  </em></td>
<td>Lower bound on inflation standard deviation.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_handle  </em></td>
<td>Ensemble handle with storage for state space inflation.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ss_inflate_index  </em></td>
<td>Index op copy in ensemble storage for inflation value.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">ss_inflate_sd_index  </em></td>
<td>Index of copy in ensemble storage for inflation standard
deviation.</td>
</tr>
<tr>
<td valign="top"><em class="code">label  </em></td>
<td>Character label to be used in diagnostic output (e.g. 'Prior',
'Posterior').</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_sd" id="get_sd"></a><br>
<div class="routine"><em class="call">var =
get_sd(inflate_handle)</em>
<pre>
real(r8), intent(out)                   :: <em class=
"code">get_sd</em>
type(adaptive_inflate_type), intent(in) :: <em class=
"code">inflate_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns value of observation space inflation standard
deviation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">get_sd  </em></td>
<td>Returns the value of observation space inflation.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for inflation descriptor.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_inflate" id="get_inflate"></a><br>
<div class="routine"><em class="call">var =
get_inflate(inflate_handle)</em>
<pre>
real(r8), intent(out)                   :: <em class=
"code">get_inflate</em>
type(adaptive_inflate_type), intent(in) :: <em class=
"code">inflate_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns value of observation space inflation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">get_inflate  </em></td>
<td>Returns the value of observation space inflation.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for inflation descriptor.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_inflate" id="set_inflate"></a><br>
<div class="routine"><em class="call">call
set_inflate(inflate_handle, inflate)</em>
<pre>
type(adaptive_inflate_type), intent(inout) :: <em class=
"code">inflate_handle</em>
real(r8), intent(in)                       :: <em class=
"code">inflate</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the value of observation space inflation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for inflation descriptor.</td>
</tr>
<tr>
<td valign="top"><em class="code">inflate  </em></td>
<td>Set observation space inflation to this value.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_sd" id="set_sd"></a><br>
<div class="routine"><em class="call">call set_sd(inflate_handle,
sd)</em>
<pre>
type(adaptive_inflate_type), intent(inout) :: <em class=
"code">inflate_handle</em>
real(r8), intent(in)                       :: <em class=
"code">sd</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Set the value of observation space inflation standard
deviation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for inflation descriptor.</td>
</tr>
<tr>
<td valign="top"><em class="code">sd  </em></td>
<td>Set observation space inflation standard deviation to this
value.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="deterministic_inflate" id=
"deterministic_inflate"></a><br>
<div class="routine"><em class="call">var =
deterministic_inflate(inflate_handle)</em>
<pre>
logical, intent(out)                    :: <em class=
"code">deterministic_inflate</em>
type(adaptive_inflate_type), intent(in) :: <em class=
"code">inflate_handle</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if deterministic inflation is being done.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">deterministic_inflate  </em></td>
<td>Returns true if deterministic inflation is being done.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">inflate_handle  </em></td>
<td>Handle for inflation descriptor.</td>
</tr>
</table>
</div>
<br>
<!--=================== DESCRIPTION OF A LOCAL TYPE ==================-->
 <a name="adaptive_inflate_type" id=
"adaptive_inflate_type"></a><br>
<div class="type">
<pre>
<em class="call">type adaptive_inflate_type</em>
   private
   integer :: inflation_flavor
   integer :: obs_diag_unit
   logical :: start_from_restart
   logical :: output_restart
   logical :: deterministic
   character(len = 129) :: in_file_name
   character(len = 129) :: out_file_name
   character(len = 129) :: diag_file_name
   real(r8) :: inflate
   real(r8) :: sd
   real(r8) :: sd_lower_bound
   real(r8) :: inf_lower_bound
   real(r8) :: inf_upper_bound
   type(random_seq_type) :: ran_seq
end type adaptive_inflate_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Provides a handle for a descriptor of inflation. Includes type
of inflation, values controlling it, input and output file names,
an output file descriptor for observation space inflation
diagnotics, and a random sequence for doing reproducible
non-determinstic inflation. There are 2 instances of this type, one
for Prior and one for Posterior inflation.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">inflation_flavor</td>
<td>Type of inflation; 0=none, 1=obs. space, 2=spatially varying,
3=spatially-fixed.</td>
</tr>
<tr>
<td valign="top">obs_diag_unit</td>
<td>Unit descriptor for output diagnostic file.</td>
</tr>
<tr>
<td valign="top">start_from_restart</td>
<td>True if initial inflation to be read from file.</td>
</tr>
<tr>
<td valign="top">output_restart</td>
<td>True if final inflation values to be written to file.</td>
</tr>
<tr>
<td valign="top">deterministic</td>
<td>True if inflation is to be done be deterministic
algorithm.</td>
</tr>
<tr>
<td valign="top">in_file_name</td>
<td>File name containing restart.</td>
</tr>
<tr>
<td valign="top">out_file_name</td>
<td>File to contain output restart.</td>
</tr>
<tr>
<td valign="top">diag_file_name</td>
<td>File to hold observation space diagnostics.</td>
</tr>
<tr>
<td valign="top">inflate</td>
<td>Initial value of inflation for all types; current value for
obs. space.</td>
</tr>
<tr>
<td valign="top">sd</td>
<td>Initial value of sd for all types; current value for obs.
space.</td>
</tr>
<tr>
<td valign="top">sd_lower_bound</td>
<td>Don't allow standard deviation to get smaller than this.</td>
</tr>
<tr>
<td valign="top">inf_lower_bound</td>
<td>Don't let inflation get smaller than this.</td>
</tr>
<tr>
<td valign="top">inf_upper_bound</td>
<td>Don't let inflation get larger than this.</td>
</tr>
<tr>
<td valign="top">ran_seq</td>
<td>Handle to random number sequence to allow reproducing
non-deterministic inflate.</td>
</tr>
</table>
</div>
<br>
<!--=================================================================-->
<!--============== DESCRIPTION OF A NAMELIST ========================-->
<!--=================================================================-->
 <a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>The adaptive_inflate module no longer has a namelist. Control
has been moved to <a href=
"filter_mod.html#Namelist">&amp;filter_nml</a> in filter.</p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<p>Three files are opened from this module, but all names are
passed in from the filter_nml now, and there are 2 values for each
name: one for the prior and one for the posterior inflation.</p>
<ul>
<li>inf_in_file_name<br>
Mean and standard deviation values read in restart file
format.</li>
<li>inf_out_file_name<br>
Mean and standard deviation values written in restart file
format.</li>
<li>inf_diag_file_name<br>
Contains diagnostic history of inflation values for obs space and
spatially-fixed state space inflation. Diagnostics for
spatially-varying state space inflation are extra fields on the
Posterior and Prior diagnostic netcdf files created in
filter.f90.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>Anderson, J. L., 2007: An adaptive covariance inflation error
correction algorithm for ensemble filters. <span style=
"font-style: italic;">Tellus A</span>, <span style=
"font-weight: bold;">59</span>, 210-224.<br>
<a href="http://dx.doi.org/10.1111/j.1600-0870.2006.00216.x"
target="_blank">doi: 10.1111/j.1600-0870.2006.00216.x</a><br></li>
<li>Anderson, J. L., 2009: Spatially and temporally varying
adaptive covariance inflation for ensemble filters. <span style=
"font-style: italic;">Tellus A</span>, <span style=
"font-weight: bold;">61</span>, 72-83.<br>
<a href="http://dx.doi.org/10.1111/j.1600-0870.2008.00361.x"
target="_blank">doi: 10.1111/j.1600-0870.2008.00361.x</a><br></li>
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
<td valign="top">adaptive_inflate_init</td>
<!-- message -->
<td valign="top">Cannot have non-deterministic inflation and
inf_lower_bound &lt; 1.</td>
<!-- comment -->
<td valign="top">Algorithm can't work in this case.<br></td>
</tr>
<tr><!-- routine -->
<td valign="top">adaptive_inflate_init</td>
<!-- message -->
<td valign="top">ss_inflate_index = ### and ss_inflate_sd_index =
### must be contiguous.</td>
<!-- comment -->
<td valign="top">Storage for these two must be contiguous in
ensemble_manager.</td>
</tr>
<tr><!-- routine -->
<td valign="top">adaptive_inflate_end</td>
<!-- message -->
<td valign="top">ss_inflate_index = ### and ss_inflate_sd_index =
### must be contiguous.</td>
<!-- comment -->
<td valign="top">Storage for these two must be contiguous in
ensemble_manager.</td>
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
<p>no discussion</p>
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
