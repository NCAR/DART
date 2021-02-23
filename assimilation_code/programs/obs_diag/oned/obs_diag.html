<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program obs_diag (for 1D observations)</title>
<link rel="stylesheet" type="text/css" href=
"../../../../docs/html/doc.css">
<link href="../../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">obs_diag</em> (for 1D
observations)</h1>
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
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview/Usage</h2>
<p>Main program for observation-space diagnostics for the models
with 1D locations. 18 quantities are calculated for each region for
each temporal bin specified by user input. The result of the code
is a netCDF file that contains the 18 quantities of the prior (aka
'guess') and posterior (aka 'analysis') estimates as a function of
time and region as well as all the metadata to create meaningful
figures. <strong>The 1D version of <em class=
"program">obs_diag</em> has defaults that automatically set the
first and last bin center based on the first and last observation
time in the set of observations being processed.</strong> This is
different behavior than the 3D versions.</p>
<p>Each <em class="file">obs_seq.final</em> file contains an
observation sequence that has multiple 'copies' of the observation.
One copy is the actual observation, another copy is the prior
ensemble mean estimate of the observation, one is the spread of the
prior ensemble estimate, one may be the prior estimate from
ensemble member 1, ... etc. The only observations for the 1D models
are generally the result of a 'perfect model' experiment, so there
is an additional copy called the 'truth' - the noise-free expected
observation given the true model state. Since this copy does not,
in general, exist for the high-order models, all comparisons are
made with the copy labelled 'observation'. There is also a namelist
variable (<em class="code">use_zero_error_obs</em>) to compare
against the 'truth' instead; the observation error variance is then
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
<td align="left"><a href=
"../../../../docs/images/lorenz_63_rmse_evolution.png"><img src=
"../../../../docs/images/lorenz_63_rmse_evolution.png" width=
"300"></a></td>
<td align="right"><a href=
"../../../../docs/images/lorenz_63_rank_histogram.png"><img src=
"../../../../docs/images/lorenz_63_rank_histogram.png" width=
"300"></a></td>
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
<p><em class="program">obs_diag</em> is designed to explore the
effect of the assimilation in two ways; 1) as a function of time
for a particular variable (this is the figure on the left), and
sometimes 2) in terms of a rank histogram - "Where does the actual
observation rank relative to the rest of the ensemble?" (figure on
the right). The figures were created by Matlab® scripts that query
the <em class="file">obs_diag_output.nc</em> file: <em class=
"file">DART/diagnostics/matlab/<a href=
"../../../../diagnostics/matlab/plot_evolution.m">plot_evolution.m</a></em>
and <a href=
"../../../../diagnostics/matlab/plot_rank_histogram.m">plot_rank_histogram.m</a>.
Both of these takes as input a file name and a 'quantity' to plot
('rmse','spread','totalspread', ...) and exhaustively plots
the quantity (for every variable, every region) in a single matlab
figure window - and creates a series of .ps files with multiple
pages for each of the figures. The directory gets cluttered with
them.</p>
<p>The observation sequence files contain only the time of the
observation, nothing of the assimilation interval, etc. - so it
requires user guidance to declare what sort of temporal binning for
the temporal evolution plots. I do a 'bunch' of arithmetic on the
namelist times to convert them to a series of temporal bin edges
that are used when traversing the observation sequence. The actual
algorithm is that the user input for the start date and bin width
set up a sequence that ends in one of two ways ... the last time is
reached or the number of bins has been reached.
<strong>NOTE:</strong> for the purpose of interpretability, the 1D
<em class="program">obs_diag</em> routines saves 'dates' as
GREGORIAN dates despite the fact these systems have no concept of a
calendar.</p>
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
<li>Support for DART QC = 8 (failed vertical conversion). This is
provided simply to make the netCDF files as consistent as needed
for plotting purposes.<br></li>
<li>Simplified input file specification.<br></li>
<li>Some of the internal variable names have been changed to make
it easier to distinguish between variances and standard
deviations.</li>
</ol>
<h2>What is new in the Lanai Release.</h2>
<p><em class="program">obs_diag</em> has several improvements:</p>
<ol>
<li>Support for 'trusted' observations. Trusted observation types
may be specified in the namelist and all observations of that type
will be counted in the statistics despite the DART QC code (as long
as the forward observation operator succeeds). See namelist
variable <em class="code">trusted_obs</em>.<br></li>
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
   bin_width_days        = -1
   bin_width_seconds     = -1
   init_skip_days        = 0
   init_skip_seconds     = 0
   max_num_bins          = 9999
   Nregions              = 3
   lonlim1               = 0.0, 0.0, 0.5
   lonlim2               = 1.0, 0.5, 1.0
   reg_names             = 'whole', 'yin', 'yang'
   trusted_obs           = 'null'
   use_zero_error_obs    = .false.
   create_rank_histogram = .true.
   outliers_in_histogram = .true.
   verbose               = .false.
   /
</pre>
<p>This value may be more useful:</p>
<pre>
   obs_sequence_name     = 'obs_seq.final'
</pre></div>
<br>
<br>
<p>The allowable ranges for the region boundaries are: lon [0.0,
1.0). The 1D locations are conceived as the distance around a unit
sphere. An observation with a location exactly ON a region boundary
cannot 'count' for both regions. The logic used to resolve this
is:</p>
<blockquote>if((lon ≥ lon1) .and. (lon &lt; lon2)) keeper =
.true.</blockquote>
<p>Consequently, if you want to include an observation precisely AT
1.0, (for example), you need to specify something a little larger
than 1.0.<br>
<br>
You can only specify <strong>either</strong> <em class=
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
<td>bin_width_days, bin_width_seconds</td>
<td>integer</td>
<td>Specifies the width of the analysis window. All observations
within a window centered at the observation time +/-
bin_width_[days,seconds] is used. If both values are 0, half the
separation between observation times as defined in the observation
sequence file is used for the bin width (i.e. all observations
used).</td>
</tr>
<tr>
<td>init_skip_days, init_skip_seconds</td>
<td>integer</td>
<td>Ignore all observations before this time. This allows one to
skip the 'spinup' or stabilization period of an assimilation.</td>
</tr>
<tr>
<td>max_num_bins</td>
<td>integer</td>
<td>This provides a way to restrict the number of temporal bins. If
<em class="code">max_num_bins</em> is set to '10', only 10
timesteps will be output, provided there are that many.</td>
</tr>
<tr>
<td>Nregions</td>
<td>integer</td>
<td>The number of regions for the unit circle for which you'd like
observation-space diagnostics. If 3 is not enough increase
<em class="code">MaxRegions</em> in <em class=
"file">obs_diag.f90</em> and recompile.</td>
</tr>
<tr>
<td>lonlim1</td>
<td>real(r8) array of length(Nregions)</td>
<td>starting value of coordinates defining 'regions'. A value of -1
indicates the start of 'no region'.</td>
</tr>
<tr>
<td>lonlim2</td>
<td>real(r8) array of length(Nregions)</td>
<td>ending value of coordinates defining 'regions'. A value of -1
indicates the end of 'no region'.</td>
</tr>
<tr>
<td>reg_names</td>
<td>character(len=6), dimension(Nregions)</td>
<td>Array of names for each of the regions. The default example has
the unit circle as a whole and divided into two equal parts, so
there are only three regions.</td>
</tr>
<tr>
<td>trusted_obs</td>
<td>character(len=32), dimension(5)</td>
<td>Array of names for observation TYPES that will be included in
the statistics if at all possible (i.e. the forward observation
operator succeeds). For 1D observations the only choices in the
code as distributed are 'RAW_STATE_VARIABLE' and/or
'RAW_STATE_1D_INTEGRAL'. (Additional 1D observation types can be
added by the user.)</td>
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
<h2>MODULES DIRECTLY USED</h2>
<pre>
types_mod
obs_sequence_mod
obs_def_mod
obs_kind_mod
location_mod
time_manager_mod
utilities_mod
sort_mod
random_seq_mod
</pre>
<h2>MODULES INDIRECTLY USED</h2>
<pre>
assim_model_mod
cov_cutoff_mod
model_mod
null_mpi_utilities_mod
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
<h3 class="indent1">Discussion of obs_diag_output.nc</h3>
<p>Every observation type encountered in the observation sequence
file is tracked separately, and aggregated into temporal and
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
  ...
  int rank_bins(rank_bins) ;
          rank_bins:long_name = "rank histogram bins" ;
          rank_bins:comment = "position of the observation among the sorted noisy ensemble members" ;
  float RAW_STATE_VARIABLE_guess(time, copy, region) ;
          RAW_STATE_VARIABLE_guess:_FillValue = -888888.f ;
          RAW_STATE_VARIABLE_guess:missing_value = -888888.f ;
  float RAW_STATE_VARIABLE_analy(time, copy, region) ;
          RAW_STATE_VARIABLE_analy:_FillValue = -888888.f ;
          RAW_STATE_VARIABLE_analy:missing_value = -888888.f ;
  ...
</pre>
<h4 class="indent1">rank histograms</h4>
<p>If it is possible to calculate a rank histogram, there will also
be :</p>
<pre>
   ...
  int RAW_STATE_VARIABLE_guess_RankHist(time, rank_bins, region) ;
   ...
</pre>
<p>as well as some global attributes. The attributes reflect the
namelist settings and can be used by plotting routines to provide
additional annotation for the histogram.</p>
<pre>
                :DART_QCs_in_histogram = 0, 1, 2, 3, 7 ;
                :outliers_in_histogram = "TRUE" ;
</pre>
<p>Please note:<br>
netCDF restricts variable names to 40 characters, so '_Rank_Hist'
may be truncated.</p>
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
<td valign="top">get_last_obs</td>
<!-- message -->
<td valign="top">No "last" observation in sequence.</td>
<!-- comment -->
<td valign="top">Generated by an incomplete observation sequence
file.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_first_obs</td>
<!-- message -->
<td valign="top">No Observations in sequence.</td>
<!-- comment -->
<td valign="top">Empty observation sequence file.</td>
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
<p>none.</p>
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
