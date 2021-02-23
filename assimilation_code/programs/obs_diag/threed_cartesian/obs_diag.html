<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program obs_diag (for 3D XYZ observations)</title>
<link rel="stylesheet" type="text/css" href=
"../../../../docs/html/doc.css">
<link href="../../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">obs_diag</em> (for observations
that use the threed_cartesian location module)</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src=
"../../../../docs/images/Dartboard7.png" alt="DART project logo"
height="70"></td>
<td>Jump to <a href="../../../../docs/index.html">DART
Documentation Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Modules">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#Usage">USAGE</a> / <a href="#References">REFERENCES</a> /
<a href="#Errors">ERRORS</a> / <a href="#FuturePlans">PLANS</a> /
<a href="#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>Main program for evaluating filter performance in observation
space. Primarily, the prior or posterior ensemble mean (and spread)
are compared to the observation and several quantities are
calculated. These quantities are then saved in a netCDF file that
has all the metadata to create meaningful figures.</p>
<p>Each <em class="file">obs_seq.final</em> file contains an
observation sequence that has multiple 'copies' of the observation.
One copy is the actual observation, another copy is the prior
ensemble mean estimate of the observation, one is the spread of the
prior ensemble estimate, one may be the prior estimate from
ensemble member 1, ... etc. If the original observation sequence is
the result of a 'perfect model' experiment, there is an additional
copy called the 'truth' - the noise-free expected observation given
the true model state. Since this copy does not, in general, exist
for the high-order models, all comparisons are made with the copy
labelled 'observation'. There is also a namelist variable
(<em class="code">use_zero_error_obs</em>) to compare against the
'truth' instead; the observation error variance is then
automatically set to zero.</p>
<p>Each ensemble member applies a forward observation operator to
the state to compute the "expected" value of an observation. Please
note: the forward observation operator is applied
<strong>AFTER</strong> any prior inflation has taken place!
Similarly, the forward observation operator is applied AFTER any
posterior inflation. This has always been the case. For a detailed
look at the relationship between the observation operators and
inflation, please look at the <a href=
"../../filter/filter.html#DetailedProgramFlow">Detailed Program Execution
Flow</a> section of <a href=
"../../filter/filter.html">filter.html</a>.<br>
<br>
Given multiple estimates of the observation, several quantities can
be calculated. It is possible to compute the expected observations
from the state vector before assimilating (the "guess", "forecast",
or "prior") or after the assimilation (the "analysis", or
"posterior").</p>
<p>Even with <em class="file">input.nml</em>:<em class=
"code">filter_nml:num_output_obs_members</em> set to <em class=
"code">0</em>; the full [prior,posterior] ensemble mean and
[prior,posterior] ensemble spread are preserved in the <em class=
"file">obs_seq.final</em> file. Consequently, the ensemble means
and spreads are used to calculate the diagnostics. If the
<em class="file">input.nml</em>:<em class=
"code">filter_nml:num_output_obs_members</em> is set to <em class=
"code">80</em> (for example); the first 80 ensemble members prior
and posterior "expected" values of the observation are also
included. In this case, the <em class="file">obs_seq.final</em>
file contains enough information to calculate a rank histograms,
verify forecasts, etc. The ensemble means are still used for many
other calculations.</p>
<table width="100%">
<tr>
<td rowspan="2"><a href=
"../../../../docs/images/obs_diag_evolution_example.png"><img src=
"../../../../docs/images/obs_diag_evolution_example.png" width=
"300"></a></td>
<td rowspan="2"><a href=
"../../../../docs/images/obs_diag_profile_example.png"><img src=
"../../../../docs/images/obs_diag_profile_example.png" width=
"300"></a></td>
<td><a href=
"../../../../docs/images/RankHistogram_ncview.png"><img src=
"../../../../docs/images/RankHistogram_ncview.png" width=
"200"></a></td>
</tr>
<tr>
<td><a href=
"../../../../docs/images/RankHistogram_matlab.png"><img src=
"../../../../docs/images/RankHistogram_matlab.png" width=
"200"></a></td>
</tr>
</table>
<p>Since this program is fundamentally interested in the response
as a function of region, there are three versions of this program;
one for each of the <em class="file">oned, threed_sphere, or
threed_cartesian</em> location modules (<em class=
"file">location_mod.f90</em>). It did not make sense to ask the
<em class="program">lorenz_96</em> model what part of North America
you'd like to investigate or how you would like to bin in the
vertical. The low-order models write out similar netCDF files and
the Matlab scripts have been updated accordingly. The oned
observations have locations conceptualized as being on a unit
circle, so only the namelist input variables pertaining to
longitude are used.</p>
<p>Identity observations (only possible from "perfect model
experiments") are already explored with state-space diagnostics, so
<em class="program">obs_diag</em> simply skips them.</p>
<p><em class="program">obs_diag</em> is designed to explore the
effect of the assimilation in three ways; 1) as a function of time
for a particular variable and level (this is the figure on the
left), 2) as a time-averaged vertical profile (figure in the
middle), and sometimes 3) in terms of a rank histogram - "Where
does the actual observation rank relative to the rest of the
ensemble?" (figures on the right). The figures on the left and
center were created by several Matlab® scripts that query the
<em class="file">obs_diag_output.nc</em> file: <em class=
"file">DART/diagnostics/matlab/<a href=
"../../../../diagnostics/matlab/plot_evolution.m">plot_evolution.m</a></em>
and <a href=
"../../../../diagnostics/matlab/plot_profile.m">plot_profile.m</a>.
Both of these takes as input a file name and a 'quantity' to plot
('rmse','spread','totalspread', ...) and exhaustively plots
the quantity (for every variable, every level, every region) in a
single matlab figure window - and creates a series of .ps files
with multiple pages for each of the figures. The directory gets
cluttered with them. The rank histogram information can easily be
plotted with <a href=
"http://meteora.ucsd.edu/~pierce/ncview_home_page.html">ncview</a>,
a free third-party piece of software or with <a href=
"../../../../diagnostics/matlab/plot_rank_histogram.m">plot_rank_histogram.m</a>.</p>
<p><em class="program">obs_diag</em> can be configured to compare
the ensemble estimates against the 'observation' copy or the
'truth' copy based on the setting of the <em class=
"code">use_zero_error_obs</em> namelist variable.</p>
<p>The observation sequence files contain only the time of the
observation, nothing of the assimilation interval, etc. - so it
requires user guidance to declare what sort of temporal binning for
the temporal evolution plots. I do a 'bunch' of arithmetic on the
namelist times to convert them to a series of temporal bin edges
that are used when traversing the observation sequence. The actual
algorithm is that the user input for the start date and bin width
set up a sequence that ends in one of two ways ... the last time is
reached or the number of bins has been reached.</p>
<p><em class="program">obs_diag</em> reads <em class=
"file">obs_seq.final</em> files and calculates the following
quantities (in no particular order) for an arbitrary number of
regions and levels. <em class="program">obs_diag</em> creates a
netCDF file called <em class="file">obs_diag_output.nc</em>. It is
necessary to query the <em class="code">CopyMetaData</em> variable
to determine the storage order (i.e. "which copy is what?") if you
want to use your own plotting routines.</p>
<div class="unix">ncdump -f F -v CopyMetaData
obs_diag_output.nc</div>
<br>
<br>
<table width="90%">
<tr>
<td valign="top"><b>Nposs</b></td>
<td>The number of observations available to be assimilated.</td>
</tr>
<tr>
<td valign="top"><b>Nused</b></td>
<td>The number of observations that were assimilated.</td>
</tr>
<tr>
<td valign="top"><b>NbadUV</b></td>
<td>the number of velocity observations that had a matching
component that was not assimilated;</td>
</tr>
<tr>
<td valign="top"><b>NbadLV</b></td>
<td>the number of observations that were above or below the highest
or lowest model level, respectively;</td>
</tr>
<tr>
<td valign="top"><b>rmse</b></td>
<td>The root-mean-squared error (the horizontal wind components are
also used to calculate the vector wind velocity and its RMS
error).</td>
</tr>
<tr>
<td valign="top"><b>bias</b></td>
<td>The simple sum of forecast - observation. The bias of the
horizontal wind speed (not velocity) is also computed.</td>
</tr>
<tr>
<td valign="top"><b>spread</b></td>
<td>The standard deviation of the univariate obs. DART does not
exploit the bivariate nature of U,V winds and so the spread of the
horizontal wind is defined as the sum of the spreads of the U and V
components.</td>
</tr>
<tr>
<td valign="top"><b>totalspread   </b></td>
<td>The total standard deviation of the estimate. We pool the
ensemble variance of the observation plus the observation error
variance and take the square root.</td>
</tr>
<tr>
<td valign="top"><b>NbadDARTQC   </b></td>
<td>the number of observations that had a DART QC value (&gt; 1 for
a prior, &gt; 3 for a posterior)</td>
</tr>
<tr>
<td valign="top"><b>observation</b></td>
<td>the mean of the observation values</td>
</tr>
<tr>
<td valign="top"><b>ens_mean</b></td>
<td>the ensemble mean of the model estimates of the observation
values</td>
</tr>
<tr>
<td valign="top"><b>N_trusted</b></td>
<td>the number of implicitly trusted observations, regardless of
DART QC</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_0</b></td>
<td>the number of observations that had a DART QC value of 0</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_1</b></td>
<td>the number of observations that had a DART QC value of 1</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_2</b></td>
<td>the number of observations that had a DART QC value of 2</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_3</b></td>
<td>the number of observations that had a DART QC value of 3</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_4</b></td>
<td>the number of observations that had a DART QC value of 4</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_5</b></td>
<td>the number of observations that had a DART QC value of 5</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_6</b></td>
<td>the number of observations that had a DART QC value of 6</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_7</b></td>
<td>the number of observations that had a DART QC value of 7</td>
</tr>
<tr>
<td valign="top"><b>N_DARTqc_8</b></td>
<td>the number of observations that had a DART QC value of 8</td>
</tr>
</table>
<p>The temporal evolution of the above quantities for every
observation type (RADIOSONDE_U_WIND_COMPONENT,
AIRCRAFT_SPECIFIC_HUMIDITY, ...) is recorded in the output netCDF
file - <em class="file"><b>obs_diag_output.nc</b></em>. This netCDF
file can then be loaded and displayed using the Matlab® scripts in
<em class="file">..../DART/diagnostics/matlab</em>. (which may
depend on functions in <em class="file">..../DART/matlab</em>). The
temporal, geographic, and vertical binning are under namelist
control. Temporal averages of the above quantities are also stored
in the netCDF file. Normally, it is useful to skip the 'burn-in'
period - the amount of time to skip is under namelist control.</p>
<p>The DART QC flag is intended to provide information about
whether the observation was assimilated, evaluated only, whether
the assimilation resulted in a 'good' observation, etc. <em class=
"green">DART QC values lower than <strong>2</strong> indicate the
prior and posteriors are OK.</em> DART QC values higher than
<strong>3</strong> were <strong>not</strong> assimilated or
evaluated. Here is the table that should explain things more
fully:</p>
<table width="80%">
<tr>
<th align="left">DART QC flag value</th>
<th align="left">meaning</th>
</tr>
<tr>
<td>0</td>
<td>observation assimilated</td>
</tr>
<tr>
<td>1</td>
<td>observation evaluated only (because of namelist settings)</td>
</tr>
<tr>
<td colspan="2">
<hr></td>
</tr>
<tr>
<td>2</td>
<td>assimilated, but the posterior forward operator failed</td>
</tr>
<tr>
<td>3</td>
<td>evaluated only, but the posterior forward operator failed</td>
</tr>
<tr>
<td colspan="2">
<hr></td>
</tr>
<tr>
<td>4</td>
<td>prior forward operator failed</td>
</tr>
<tr>
<td>5</td>
<td>not used because observation type not listed in namelist</td>
</tr>
<tr>
<td>6</td>
<td>rejected because incoming observation QC too large</td>
</tr>
<tr>
<td>7</td>
<td>rejected because of a failed outlier threshold test</td>
</tr>
<tr>
<td><em class="changed">8</em></td>
<td><em class="changed">vertical conversion failed</em></td>
</tr>
<tr>
<td>9+</td>
<td>reserved for future use</td>
</tr>
</table>
<!--==================================================================-->
<p><!-- spacer for top --></p>
<a name="New" id="New"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>What is new in the Manhattan Release.</h2>
<ol>
<li>Support for DART QC = 8 (failed vertical conversion).<br></li>
<li>Simplified input file specification.<br></li>
<li>Replace namelist integer variable <em class="code">debug</em>
with logical variable <em class="code">verbose</em> to control
amount of run-time output.<br></li>
<li>Removed <em class="code">rat_cri</em> and <em class=
"code">input_qc_threshold</em> from the namelists. They had been
deprecated for quite some time.<br></li>
<li>Some of the internal variable names have been changed to make
it easier to distinguish between variances and standard
deviations.</li>
</ol>
<h2>What is new in the Lanai Release.</h2>
<p><em class="program">obs_diag</em> has several improvements:</p>
<ol>
<li>Improved vertical specification. Namelist variables <em class=
"code">[h,p,m]level_edges</em> allow fine-grained control over the
vertical binning. It is not allowed to specify both the edges and
midpoints for the vertical bins.<br></li>
<li>Improved error-checking for input specification, particularly
the vertical bins. Repeated values are squeezed out.<br></li>
<li>Support for 'trusted' observations. Trusted observation types
may be specified in the namelist and all observations of that type
will be counted in the statistics despite the DART QC code (as long
as the forward observation operator succeeds). See namelist
variable <em class="code">trusted_obs</em>. For more details, see
the section on <a href="#Trusted">Trusted
observations</a>.<br></li>
<li>Support for 'true' observations (i.e. from an OSSE). If the
'truth' copy of an observation is desired for comparison (instead
of the default copy) the observation error variance is set to 0.0
and the statistics are calculated relative to the 'truth' copy (as
opposed to the normal 'noisy' or 'observation' copy). See namelist
variable <em class="code">use_zero_error_obs</em>.<br></li>
<li>discontinued the use of <em class="code">rat_cri</em> and
<em class="code">input_qc_threshold</em> namelist variables. Their
functionality was replaced by the DART QC mechanism long
ago.<br></li>
<li>The creation of the rank histogram (if possible) is now
namelist-controlled by namelist variable <em class=
"code">create_rank_histogram</em>.</li>
</ol>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<p><!-- spacer for top --></p>
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
&amp;obs_diag_nml
   obs_sequence_name     = ''
   obs_sequence_list     = ''
   first_bin_center      =  2003, 1, 1, 0, 0, 0
   last_bin_center       =  2003, 1, 2, 0, 0, 0
   bin_separation        =     0, 0, 0, 6, 0, 0
   bin_width             =     0, 0, 0, 6, 0, 0
   time_to_skip          =     0, 0, 0, 6, 0, 0
   max_num_bins          = 1000
   hlevel                = -888888.0
   hlevel_edges          = -888888.0
   Nregions              = 0
   xlim1                 = -888888.0
   xlim2                 = -888888.0
   ylim1                 = -888888.0
   ylim2                 = -888888.0
   reg_names             = 'null'
   trusted_obs           = 'null'
   create_rank_histogram = .true.
   outliers_in_histogram = .false.
   use_zero_error_obs    = .false.
   verbose               = .false.
   /
</pre>
<p>These values may be more useful:</p>
<pre>
   obs_sequence_list     = 'input_file_list.txt'
   hlevel_edges          = -1.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0, 32000.0, 64000.0
   Nregions              = 1
   reg_names             = 'everywhere - maybe'
   xlim1                 = -1.0
   xlim2                 = 1000000.0
   ylim1                 = -1.0
   ylim2                 = 1000000.0
</pre></div>
<br>
<br>
<p>The date-time integer arrays in this namelist have the form
(YYYY, MM, DY, HR, MIN, SEC).<br>
The allowable ranges for the region boundaries are: latitude
[-90.,90], longitude [0.,Inf.]</p>
<p>You can only specify <strong>either</strong> <em class=
"code">obs_sequence_name</em> <strong>or</strong> <em class=
"code">obs_sequence_list</em> -- not both. One of them has to be an
empty string ... i.e. <em class="code">''</em>.</p>
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
<td>obs_sequence_name</td>
<td>character(len=256), dimension(100)</td>
<td>An array of names of observation sequence files. These may be
relative or absolute filenames. If this is set, <em class=
"code">obs_sequence_list</em> must be set to ' ' (empty
string).</td>
</tr>
<tr>
<td>obs_sequence_list</td>
<td>character(len=256)</td>
<td>Name of an ascii text file which contains a list of one or more
observation sequence files, one per line. If this is specified,
<em class="code">obs_sequence_name</em> must be set to ' '.
Can be created by any method, including sending the output of the
'ls' command to a file, a text editor, or another program. If this
is set, <em class="code">obs_sequence_name</em> must be set to
' ' (empty string).</td>
</tr>
<tr>
<td>first_bin_center</td>
<td>integer, dimension(6)</td>
<td>first timeslot of the first obs_seq.final file to process. The
six integers are: year, month, day, hour, hour, minute, second, in
that order. <em class="program">obs_diag</em> has improved run-time
output that reports the time and date of the first and last
observations in every observation sequence file. Look for the
string 'First observation date' in the logfile. If the <em class=
"code">verbose</em> is 'true', it is also written to the
screen.</td>
</tr>
<tr>
<td>last_bin_center</td>
<td>integer, dimension(6)</td>
<td>last timeslot of interest. (reminder: the last timeslot of day
1 is hour 0 of day 2) The six integers are: year, month, day, hour,
hour, minute, second, in that order. This does not need to be
exact, the values from <em class="code">first_bin_center</em> and
<em class="code">bin_separation</em> are used to populate the time
array and stop on or before the time defined by <em class=
"code">last_bin_center</em>. See also <em class=
"code">max_num_bins</em>.</td>
</tr>
<tr>
<td>bin_separation</td>
<td>integer, dimension(6)</td>
<td>Time between bin centers. The year and month values
<em>must</em> be zero.</td>
</tr>
<tr>
<td>bin_width</td>
<td>integer, dimension(6)</td>
<td>Time span around bin centers in which obs will be compared. The
year and month values <em>must</em> be zero. Frequently, but not
required to be, the same as the values for bin_separation. 0</td>
</tr>
<tr>
<td>time_to_skip</td>
<td>integer, dimension(6)</td>
<td>Time span at the beginning to skip when calculating vertical
profiles of rms error and bias. The year and month values
<em>must</em> be zero. Useful because it takes some time for the
assimilation to settle down from the climatological spread at the
start. <em class="code">time_to_skip</em> is an amount of time
AFTER the first edge of the first bin.</td>
</tr>
<tr>
<td>max_num_bins</td>
<td>integer</td>
<td>This provides an alternative way to declare the <em class=
"code">last_bin_center</em>. If <em class="code">max_num_bins</em>
is set to '10', only 10 timesteps will be output - provided
<em class="code">last_bin_center</em> is set to some later
date.</td>
</tr>
<tr>
<td>hlevel</td>
<td>real, dimension(50)</td>
<td>Same, but for observations that have height(m) or depth(m) as
the vertical coordinate. An example of defining the midpoints
is:<br>
<em class="code">hlevel = 1000, 2000, 3000, 4000, 5000, 6000, 7000,
8000, 9000, 10000, 11000,</em></td>
</tr>
<tr>
<td>hlevel_edges</td>
<td>real, dimension(51)</td>
<td>The edges defining the height (or depth) levels for the
vertical binning. You may specify either <em class=
"code">hlevel</em> or <em class="code">hlevel_edges</em>, but not
both. An example of defining the edges is:<br>
<em class="code">hlevel_edges = 0, 1500, 2500, 3500, 4500, 5500,
6500,</em></td>
</tr>
<tr>
<td>Nregions</td>
<td>integer</td>
<td>Number of regions of the globe for which obs space diagnostics
are computed separately. Must be between [1,50]. If 50 is not
enough, increase <em class="file">obs_diag.f90</em><em class=
"code">MaxRegions</em> and recompile.</td>
</tr>
<tr>
<td>xlim1</td>
<td>real, dimension(50)</td>
<td>western extent of each of the regions.</td>
</tr>
<tr>
<td>xlim2</td>
<td>real, dimension(50)</td>
<td>eastern extent of each of the regions.</td>
</tr>
<tr>
<td>ylim1</td>
<td>real, dimension(50)</td>
<td>southern extent of the regions.</td>
</tr>
<tr>
<td>ylim2</td>
<td>real, dimension(50)</td>
<td>northern extent of the regions.</td>
</tr>
<tr>
<td>reg_names</td>
<td>character(len=129), dimension(50)</td>
<td>Array of names for the regions to be analyzed. Will be used for
plot titles.</td>
</tr>
<tr>
<td>trusted_obs</td>
<td>character(len=32), dimension(50)</td>
<td>list of observation types that <strong>must</strong>
participate in the calculation of the statistics, regardless of the
DART QC (provided that the forward observation operator can still
be applied without failure). e.g. 'RADIOSONDE_TEMPERATURE', ... For
more details, see the section on <a href="#Trusted">Trusted
observations</a>.</td>
</tr>
<tr>
<td>use_zero_error_obs</td>
<td>logical</td>
<td>if <em class="code">.true.</em>, the observation copy used for
the statistics calculations will be 'truth'. Only 'perfect'
observations (from <em class="program">perfect_model_obs</em>) have
this copy. The observation error variance will be set to zero.</td>
</tr>
<tr>
<td>create_rank_histogram</td>
<td>logical</td>
<td>if <em class="code">.true.</em> and there are actual ensemble
estimates of the observations in the <em class=
"file">obs_seq.final</em> (i.e. <em class=
"code">filter_nml:num_output_obs_members</em> is larger than zero),
a rank histogram will be created.</td>
</tr>
<tr>
<td>outliers_in_histogram</td>
<td>logical</td>
<td>if <em class="code">.true.</em> the observations that have been
rejected by the outlier threshhold mechanism will be
<em>included</em> in the calculation of the rank histogram.</td>
</tr>
<tr>
<td>verbose</td>
<td>logical</td>
<td>switch controlling amount of run-time output.</td>
</tr>
</tbody>
</table>
</div>
<!--==================================================================-->
<p><!-- spacer for top --></p>
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
obs_sequence_mod
obs_kind_mod
obs_def_mod (and possibly other obs_def_xxx mods)
assim_model_mod
random_seq_mod
model_mod
location_mod
types_mod
time_manager_mod
utilities_mod
sort_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<p><!-- spacer for top --></p>
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li><em class="file">input.nml</em> is used for <em class=
"code">obs_diag_nml</em></li>
<li><em class="file">obs_diag_output.nc</em> is the netCDF output
file</li>
<li><em class="file">dart_log.out</em> list directed output from
the obs_diag.</li>
<li><em class="file">LargeInnov.txt</em> contains the distance
ratio histogram -- useful for estimating the distribution of the
magnitudes of the innovations.</li>
</ul>
<p>Obs_diag may require a model input file from which to get grid
information, metadata, and links to modules providing the models
expected observations. It all depends on what's needed by the
<em class="file">model_mod.f90</em></p>
<h3 class="indent1">Discussion of obs_diag_output.nc</h3>
<p>Every observation type encountered in the observation sequence
file is tracked separately, and aggregated into temporal and 3D
spatial bins. There are two main efforts to this program. One is to
track the temporal evolution of any of the quantities available in
the netCDF file for any possible observation type:</p>
<div class="unix">ncdump -v CopyMetaData,ObservationTypes
obs_diag_output.nc</div>
<p>The other is to explore the vertical profile of a particular
observation kind. By default, each observation kind has a
'guess/prior' value and an 'analysis/posterior' value - which shed
some insight into the innovations.</p>
<hr width="40%" align="left">
<h4 class="indent1">temporal evolution</h4>
<p>The <em class="file">obs_diag_output.nc</em> output file has all
the metadata I could think of, as well as separate variables for
every observation type in the observation sequence file.
Furthermore, there is a separate variable for the 'guess/prior' and
'analysis/posterior' estimate of the observation. To distinguish
between the two, a suffix is appended to the variable name. An
example seems appropriate:</p>
<pre>
  ...
  char CopyMetaData(copy, stringlength) ;
          CopyMetaData:long_name = "quantity names" ;
  char ObservationTypes(obstypes, stringlength) ;
          ObservationTypes:long_name = "DART observation types" ;
          ObservationTypes:comment = "table relating integer to observation type string" ;
  float RADIOSONDE_U_WIND_COMPONENT_guess(time, copy, hlevel, region) ;
          RADIOSONDE_U_WIND_COMPONENT_guess:_FillValue = -888888.f ;
          RADIOSONDE_U_WIND_COMPONENT_guess:missing_value = -888888.f ;
  float RADIOSONDE_V_WIND_COMPONENT_guess(time, copy, hlevel, region) ;
          RADIOSONDE_V_WIND_COMPONENT_guess:_FillValue = -888888.f ;
          RADIOSONDE_V_WIND_COMPONENT_guess:missing_value = -888888.f ;
  ...
  float MARINE_SFC_ALTIMETER_guess(time, copy, surface, region) ;
          MARINE_SFC_ALTIMETER_guess:_FillValue = -888888.f ;
          MARINE_SFC_ALTIMETER_guess:missing_value = -888888.f ;
  ...
  float RADIOSONDE_WIND_VELOCITY_guess(time, copy, hlevel, region) ;
          RADIOSONDE_WIND_VELOCITY_guess:_FillValue = -888888.f ;
          RADIOSONDE_WIND_VELOCITY_guess:missing_value = -888888.f ;
  ...
  float RADIOSONDE_U_WIND_COMPONENT_analy(time, copy, hlevel, region) ;
          RADIOSONDE_U_WIND_COMPONENT_analy:_FillValue = -888888.f ;
          RADIOSONDE_U_WIND_COMPONENT_analy:missing_value = -888888.f ;
  float RADIOSONDE_V_WIND_COMPONENT_analy(time, copy, hlevel, region) ;
          RADIOSONDE_V_WIND_COMPONENT_analy:_FillValue = -888888.f ;
          RADIOSONDE_V_WIND_COMPONENT_analy:missing_value = -888888.f ;
  ...
</pre>
<p>There are several things to note:</p>
<ol>
<li>the 'WIND_VELOCITY' component is nowhere 'near' the
corresponding U,V components.</li>
<li>all of the 'guess' variables come before the matching 'analy'
variables.</li>
<li>surface variables (i.e. <tt>MARINE_SFC_ALTIMETER</tt> have a
coordinate called 'surface' as opposed to 'hlevel' for the others
in this example).</li>
</ol>
<hr width="40%" align="left">
<h4 class="indent1">vertical profiles</h4>
<p>Believe it or not, there are another set of netCDF variables
specifically for the vertical profiles, essentially duplicating the
previous variables but <b>without the 'time' dimension</b>. These
are distinguished by the suffix added to the observation kind -
'VPguess' and 'VPanaly' - 'VP' for Vertical Profile.</p>
<pre>
  ...
  float SAT_WIND_VELOCITY_VPguess(copy, hlevel, region) ;
          SAT_WIND_VELOCITY_VPguess:_FillValue = -888888.f ;
          SAT_WIND_VELOCITY_VPguess:missing_value = -888888.f ;
  ...
  float RADIOSONDE_U_WIND_COMPONENT_VPanaly(copy, hlevel, region) ;
          RADIOSONDE_U_WIND_COMPONENT_VPanaly:_FillValue = -888888.f ;
          RADIOSONDE_U_WIND_COMPONENT_VPanaly:missing_value = -888888.f ;
  ...
</pre>
<p>Observations flagged as 'surface' do not participate in the
vertical profiles (Because surface variables cannot exist on any
other level, there's not much to plot!). Observations on the lowest
level DO participate. There's a difference!</p>
<hr width="40%" align="left">
<h4 class="indent1">rank histograms</h4>
<p>If it is possible to calculate a rank histogram, there will also
be :</p>
<pre>
   ...
   int RADIOSONDE_U_WIND_COMPONENT_guess_RankHi(time, rank_bins, hlevel, region) ;
   ...
   int RADIOSONDE_V_WIND_COMPONENT_guess_RankHi(time, rank_bins, hlevel, region) ;
   ...
   int MARINE_SFC_ALTIMETER_guess_RankHist(time, rank_bins, surface, region) ;
   ...
</pre>
<p>as well as some global attributes. The attributes reflect the
namelist settings and can be used by plotting routines to provide
additional annotation for the histogram.</p>
<pre>
                :DART_QCs_in_histogram = 0, 1, 2, 3, 7 ;
                :outliers_in_histogram = "TRUE" ;
</pre>
<p>Please note:</p>
<ol>
<li>netCDF restricts variable names to 40 characters, so
'_Rank_Hist' may be truncated.</li>
<li>It is sufficiently vague to try to calculate a rank histogram
for a velocity derived from the assimilation of U,V components such
that NO rank histogram is created for velocity. A run-time log
message will inform as to which variables are NOT having a rank
histogram variable preserved in the <em class=
"file">obs_diag_output.nc</em> file - IFF it is possible to
calculate a rank histogram in the first place.</li>
</ol>
<table width="100%">
<tr>
<td><a href=
"../../../../docs/images/RankHistogram_ncview.png"><img src=
"../../../../docs/images/RankHistogram_ncview.png" width=
"200"></a></td>
<td><a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Documentation.php#ncview_histogram">
Instructions for viewing the rank histogram with ncview</a>.</td>
</tr>
<tr>
<td><a href=
"../../../../docs/images/RankHistogram_matlab.png"><img src=
"../../../../docs/images/RankHistogram_matlab.png" width=
"200"></a></td>
<td><a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Documentation.php#mat_obs">
Instructions for viewing the rank histogram with Matlab</a>.</td>
</tr>
</table>
<hr width="40%" align="left">
<a name="Trusted" id="Trusted"></a>
<h4 class="indent1">"Trusted" observation types</h4>
<p>This needs to be stated up front: <em class=
"program">obs_diag</em> is a post-processor; it cannot influence
the assimilation. One interpretation of a TRUSTED observation is
that the assimilation should <strong>always</strong> use the
observation, even if it is far from the ensemble. At present (23
Feb 2015), the filter in DART does not forcibly assimilate any one
observation and selectively assimilate the others. Still, it is
useful to explore the results using a set of 'trusted type'
observations, whether they were assimilated, evaluated, or rejected
by the outlier threshhold. This is the important distinction. The
diagnostics can be calculated differently for each <em>observation
type</em>.</p>
<p>The normal diagnostics calculate the metrics (rmse, bias, etc.)
only for the 'good' observations - those that were assimilated or
evaluated. The <em class="code">outlier_threshold</em> essentially
defines what observations are considered too far from the ensemble
<strong>prior</strong> to be useful. These observations get a DART
QC of 7 and are not assimilated. The observations with a DART QC of
7 do not contribute the the metrics being calculated. Similarly, if
the forward observation operator fails, these observations cannot
contribute. When the operator fails, the 'expected' observation
value is 'MISSING', and there is no ensemble mean or spread.</p>
<p>'Trusted type' observation metrics are calculated using all the
observations that were assimilated or evaluated
<strong>AND</strong> the observations that were rejected by the
outlier threshhold. <em class="program">obs_diag</em> can
post-process the DART QC and calculate the metrics appropriately
for <strong>observation types</strong> listed in the <em class=
"code">trusted_obs</em> namelist variable. If there are trusted
observation types specified for <em class="program">obs_diag</em>,
the <em class="file">obs_diag_output.nc</em> has global metadata to
indicate that a different set of criteria were used to calculate
the metrics. The individual variables also have an extra attribute.
In the following output, <em class=
"file">input.nml:obs_diag_nml:trusted_obs</em> was set: <em class=
"mono">trusted_obs = 'RADIOSONDE_TEMPERATURE', 'RADIOSONDE_U_WIND_COMPONENT'</em></p>
<pre>
  ...
        float RADIOSONDE_U_WIND_COMPONENT_guess(time, copy, hlevel, region) ;
                RADIOSONDE_U_WIND_COMPONENT_guess:_FillValue = -888888.f ;
                RADIOSONDE_U_WIND_COMPONENT_guess:missing_value = -888888.f ;
                RADIOSONDE_U_WIND_COMPONENT_guess:<em>TRUSTED = "TRUE" ;</em>
        float RADIOSONDE_V_WIND_COMPONENT_guess(time, copy, hlevel, region) ;
                RADIOSONDE_V_WIND_COMPONENT_guess:_FillValue = -888888.f ;
                RADIOSONDE_V_WIND_COMPONENT_guess:missing_value = -888888.f ;
  ...
// global attributes:
  ...
                :<em>trusted_obs_01 = "RADIOSONDE_TEMPERATURE" ;</em>
                :<em>trusted_obs_02 = "RADIOSONDE_U_WIND_COMPONENT" ;</em>
                :obs_seq_file_001 = "cam_obs_seq.1978-01-01-00000.final" ;
                :obs_seq_file_002 = "cam_obs_seq.1978-01-02-00000.final" ;
                :obs_seq_file_003 = "cam_obs_seq.1978-01-03-00000.final" ;
  ...
                :MARINE_SFC_ALTIMETER = 7 ;
                :LAND_SFC_ALTIMETER = 8 ;
                :RADIOSONDE_U_WIND_COMPONENT<em>--TRUSTED</em> = 10 ;
                :RADIOSONDE_V_WIND_COMPONENT = 11 ;
                :RADIOSONDE_TEMPERATURE<em>--TRUSTED</em> = 14 ;
                :RADIOSONDE_SPECIFIC_HUMIDITY = 15 ;
                :AIRCRAFT_U_WIND_COMPONENT = 21 ;
  ...
</pre>
<table width="100%">
<tr>
<td width="45%">The Matlab scripts try to ensure that the trusted
observation graphics clarify that the metrics plotted are somehow
'different' than the normal processing stream. Some text is added
to indicate that the values include the outlying observations.
<strong>IMPORTANT:</strong> The interpretation of the number of
observations 'possible' and 'used' still reflects what was used
<strong>in the assimilation!</strong> The number of observations
rejected by the outlier threshhold is not explicilty plotted. To
reinforce this, the text for the observation axis on all graphics
has been changed to <em class=
"mono">"o=possible, *=assimilated"</em>. In short, the
distance between the number of observations possible and the number
assimilated still reflects the number of observations rejected by
the outlier threshhold and the number of failed forward observation
operators.</td>
<td>  <a href=
"../../../../docs/images/RAD_T_trusted_bias_evolution.png"><img src="../../../../docs/images/RAD_T_trusted_bias_evolution.png"
width="600"></a></td>
</tr>
</table>
<p>There is ONE ambiguous case for trusted observations. There may
be instances in which the observation fails the outlier threshhold
test (which is based on the prior) and the posterior forward
operator fails. DART does not have a QC that explicilty covers this
case. The current logic in <em class="program">obs_diag</em>
correctly handles these cases <strong>except</strong> when trying
to use 'trusted' observations. There is a section of code in
<em class="program">obs_diag</em> that may be enabled if you are
encountering this ambiguous case. As <em class=
"program">obs_diag</em> runs, a warning message is issued and a
summary count is printed if the ambiguous case is encountered. What
normally happens is that if that specific observation type is
trusted, the posterior values include a MISSING value in the
calculation which makes them inaccurate. If the block of code is
enabled, the DART QC is recast as the PRIOR forward observation
operator fails. This is technically incorrect, but for the case of
trusted observations, it results in only calculating statistics for
trusted observations that have a useful prior and posterior.
<strong>This should not be used unless you are willing to
intentionally disregard 'trusted' observations that were rejected
by the outlier threshhold.</strong> Since the whole point of a
trusted observation is to <em>include</em> observations potentially
rejected by the outlier threshhold, you see the problem. Some
people like to compare the posteriors. <em>THAT</em> can be the
problem.</p>
<pre>
if ((qc_integer == 7) .and. (abs(posterior_mean(1) - MISSING_R8) &lt; 1.0_r8)) then
            write(string1,*)'WARNING ambiguous case for obs index ',obsindex
            string2 = 'obs failed outlier threshhold AND posterior operator failed.'
            string3 = 'Counting as a Prior QC == 7, Posterior QC == 4.'
            if (trusted) then
! COMMENT      string3 = 'WARNING changing DART QC from 7 to 4'
! COMMENT      qc_integer = 4
            endif
            call error_handler(E_MSG,'obs_diag',string1,text2=string2,text3=string3)
            num_ambiguous = num_ambiguous + 1
         endif</pre>
<!--==================================================================-->
<!-- Discuss  typical usage of obs_diag.                              -->
<!--==================================================================-->
<p><!-- spacer for top --></p>
<a name="Usage" id="Usage"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>USAGE</h2>
<p><em class="program">obs_diag</em> is built in
.../DART/models/<em class="italic">your_model</em>/work, in the
same way as the other DART components.</p>
<h3 class="indent1">multiple observation sequence files</h3>
<p>There are two ways to specify input files for <em class=
"program">obs_diag</em>. You can either specify the name of a file
containing a list of files (in <em class=
"code">obs_sequence_list</em>), or you may specify a list of files
via <em class="code">obs_sequence_name</em>.</p>
<h3 class="indent1">Example: observation sequence files
spanning 30 days.</h3>
<table width="95%">
<tr>
<td>In this example, we will be accumulating metrics for 30 days.
The <em class="file">obs_diag_output.nc</em> file will have exactly
ONE timestep in it (so it won't be much use for the <em class=
"program">plot_evolution</em> functions) - but the <em class=
"program">plot_profile</em> functions and the <em class=
"program">plot_rank_histogram</em> function will be used to explore
the assimilation. By way of an example, we will NOT be using
outlier observations in the rank histogram. Lets presume that all
your <em class="file">obs_seq.final</em> files are in
alphabetically-nice directories:</td>
<td>  <img src=
"../../../../docs/images/RankHistogram_matlab.png" width=
"200"></td>
</tr>
</table>
<pre>
/Exp1/Dir01/obs_seq.final
/Exp1/Dir02/obs_seq.final
/Exp1/Dir03/obs_seq.final
...
/Exp1/Dir99/obs_seq.final
</pre>
<p>The first step is to create a file containing the list of
observation sequence files you want to use. This can be done with
the unix command 'ls' with the -1 option (that's a number one) to
put one file per line.</p>
<div class="unix">ls -1 /Exp1/Dir*/obs_seq.final &gt;
obs_file_list.txt</div>
<p>It is necessary to turn on the verbose option to check the
first/last times that will be used for the histogram. Then, the
namelist settings for 2008 07 31 12Z through 2008 08 30 12Z
are:</p>
<div class="routine">
<pre>
&amp;obs_diag_nml
   obs_sequence_name     = ''
   obs_sequence_list     = <em class=
"input">'obs_file_list.txt'</em>
   first_bin_center      = <em class=
"input"> 2008, 8,15,12, 0, 0</em>
   last_bin_center       = <em class=
"input"> 2008, 8,15,12, 0, 0</em>
   bin_separation        = <em class=
"input">    0, 0,30, 0, 0, 0</em>
   bin_width             = <em class=
"input">    0, 0,30, 0, 0, 0</em>
   time_to_skip          = <em class=
"input">    0, 0, 0, 0, 0, 0</em>
   max_num_bins          = <em class="input">1000</em>
   Nregions              = <em class="input">1</em>
   xlim1                 = <em class="input">-1.0</em>
   xlim2                 = <em class="input">1000000.0</em>
   ylim1                 = <em class="input">-1.0</em>
   ylim2                 = <em class="input">1000000.0</em>
   reg_names             = <em class="input">'Entire Domain'</em>
   create_rank_histogram = <em class="input">.true.</em>
   outliers_in_histogram = <em class="input">.false.</em>
   verbose               = <em class="input">.true.</em>
   /
</pre></div>
<p>then, simply run <em class="program">obs_diag</em> in the usual
manner - you may want to save the run-time output to a file. Here
is a portion of the run-time output:</p>
<div class="unix">
<pre>
...
Region  1 Entire Domain                    (WESN):     0.0000   360.0000   -90.0000    90.0000
 Requesting            1  assimilation periods.

epoch      1  start day=148865, sec=43201
epoch      1 center day=148880, sec=43200
epoch      1    end day=148895, sec=43200
epoch      1  start 2008 Jul 31 12:00:01
epoch      1 center 2008 Aug 15 12:00:00
epoch      1    end 2008 Aug 30 12:00:00
...
MARINE_SFC_HORIZONTAL_WIND_guess_RankHis has            0 "rank"able observations.
SAT_HORIZONTAL_WIND_guess_RankHist       has            0 "rank"able observations.
...
</pre></div>
<p>Discussion: It should be pretty clear that there is exactly 1
assimilation period, it may surprise you that the start is 1 second
past 12Z. This is deliberate and reflects the DART convention of
starting intervals 1 second after the end of the previous interval.
The times in the netCDF variables reflect the defined start/stop of
the period, regardless of the time of the first/last
observation.<br>
<br>
Please note that none of the 'horizontal_wind' variables will have
a rank histogram, so they are not written to the netCDF file. ANY
variable that does not have a rank histogram with some observations
will NOT have a rank histogram variable in the netCDF file.<br>
<br>
Now that you have the <em class="file">obs_diag_output.nc</em>, you
can explore it with <a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Documentation.php#mat_obs">
plot_profile.m, plot_bias_xxx_profile.m, or
plot_rmse_xxx_profile.m</a>, and look at the rank histograms with
<a href=
"http://meteora.ucsd.edu/~pierce/ncview_home_page.html">ncview</a>
or <em class="program">plot_rank_histogram.m</em>.</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<p><!-- spacer for top --></p>
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
<p><!-- spacer for top --></p>
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
<td valign="top">obs_diag</td>
<!-- message -->
<td valign="top">No first observation in sequence.</td>
<!-- comment -->
<td valign="top">get_first_obs couldn't find a "first obs" in the
obs_seq.final.</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_diag</td>
<!-- message -->
<td valign="top">No last observation in sequence</td>
<!-- comment -->
<td valign="top">get_last_obs couldn't find a "last obs" in the
obs_seq.final</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_diag</td>
<!-- message -->
<td valign="top">metadata incomplete</td>
<!-- comment -->
<td valign="top">Couldn't find the index for the observation value
in the observation sequence file. It is the only one that is
required.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_get_obs_info</td>
<!-- message -->
<td valign="top">Vertical coordinate not recognized</td>
<!-- comment -->
<td valign="top">It must be "surface", "pressure", or "height"</td>
</tr>
<tr><!-- routine -->
<td valign="top">Convert2Time</td>
<!-- message -->
<td valign="top">namelist parameter out-of-bounds. Fix and try
again</td>
<!-- comment -->
<td valign="top">bin_width, bin_separation, and time_to_skip must
have year = 0 and month = 0</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<p><!-- spacer for top --></p>
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>The RMSE of the vector wind velocity is being used. The bias is
actually the bias of the wind speed, with no regard to direction.
Seems like this is not consistent ...<br>
<br>
Use log(p) instead of pressure for binning/plotting.<br>
<br>
Hope to have a separate DART QC flag for observations outside the
model state. Right now, all observations that have a DARt forward
operator fail (extrapolate, mainly) get counted as observations
that are rejected. Logically, these observations are not possible
because most (all?) DART observation operators cannot
extrapolate.<br>
<br>
If the specified levels are not in the proper range for the
observations, no vertical profile is built. The 'verbose' option
should highlight the fact there are no observations in the
specified vertical bins - AND provide some sort of guidance about
the min and max values of the vertical values.</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<p><!-- spacer for top --></p>
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<p><!-- spacer for top --></p>
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
