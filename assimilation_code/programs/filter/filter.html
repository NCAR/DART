<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program filter</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">filter</em></h1>
<table border="0" summary="dart header" cellpadding="5">
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
<p>Main program for driving ensemble filter assimilations.</p>
<p><em class="program">filter</em> is a Fortran 90 program, and
provides a large number of options for controlling execution
behavior and parameter configuration that are driven from its
namelist. See the <a href="#Namelist">namelist</a> section below
for more details. The number of assimilation steps to be done is
controlled by the input observation sequence and by the
time-stepping capabilities of the model being used in the
assimilation.</p>
<p>This overview includes these subsections:</p>
<ul>
<li><a href="#ProgramFlow">Program Flow</a></li>
<li><a href="#FilterTypes">Filter Types</a></li>
<li><a href="#GettingStarted">Getting Started</a></li>
<li><a href="#FreeRun">Free Model Run after Assimilation</a></li>
<li><a href="#EvalOnly">Evaluate a Model State against
Observations</a></li>
<li><a href="#Verify">Compare Results with and without
Assimilation</a></li>
<li><a href="#QCVals">DART Quality Control Values on
Output</a></li>
<li><a href="#Inflation">Description of Inflation Options</a></li>
<li><a href="#DetailedProgramFlow">Detailed Program Flow</a></li>
</ul>
<p>See the <a href="http://www.image.ucar.edu/DAReS/DART">DART web
site</a> for more documentation, including a discussion of the
capabilities of the assimilation system, a diagram of the entire
execution cycle, the options and features.</p>
<a name="ProgramFlow" id="ProgramFlow"></a>
<h3>Program Flow</h3>
<p>The basic execution loop is:</p>
<ul>
<li>Read in model initial conditions, observations, set up and
initialize</li>
<li>Until out of observations:
<ul>
<li>Run multiple copies of the model to get forecasts of model
state</li>
<li>Assimilate all observations in the current time window</li>
<li>Repeat</li>
</ul>
</li>
<li>Write out diagnostic files, restart files, final observation
sequence file</li>
</ul>
<p>The time of the observations in the input observation sequence
file controls the length of execution of filter.</p>
<p>For large, parallel models, the execution loop is usually
wrapped in an external script which does these additional
steps:</p>
<ul>
<li>Link to an observation sequence file which contains only
observation times within the next assimilation window</li>
<li>Link any output inflation files from the previous step to be
the input files for this step</li>
<li>Run filter, which will exit after doing the assimilation
without trying to advance the model</li>
<li>Save the output diagnostic files for later</li>
<li>Advance the N copies of the model using the model scripts or
whatever method is appropriate</li>
<li>Repeat until all data is assimilated</li>
</ul>
<p>For large models filter is almost always compiled to be a
parallel MPI program, and most large models are themselves a
parallel program using OpenMP, MPI, or both. MPI programs usually
cannot start other MPI programs, so the external script submits
both the filter job and the N model advances to a batch system so
all run as independent parallel jobs.</p>
<p>The same source code is used for all applications of filter. The
code specific to the types of observations and the interface code
for the computational model is configured at compile time. The top
level directory has been simplified from previous versions to look
like :</p>
<ul>
<li><em class="file">README</em></li>
<li><em class="file">COPYRIGHT</em></li>
<li><em class="dir">assimilation_code</em></li>
<li><em class="dir">build_templates</em></li>
<li><em class="dir">diagnostics</em></li>
<li><em class="dir">documentation</em></li>
<li><em class="dir">models</em></li>
<li><em class="dir">observations</em></li>
</ul>
the <em class="dir">assimilation_code</em> contains all <em class=
"dir">module</em> and <em class="dir">program</em> source code for
all of the main programs including filter. Specifically in the
modules directory there is a <em class="file">filter_mod.f90</em>
which contains the source for the filter main program. Each model
has a separate directory under DART/models, and under each model is
a work directory where the code is compiled and can be run for
testing. Generally when a full-size experiment is done the
executables are copied to a different location - e.g. scratch space
on a large filesystem - since the data files for 10s to 100s of
copies of a model can get very large. A lightly pruned directory
tree can be browsed in the main <a href=
"../../../docs/index.html#Directories">index.html</a>. <a name=
"FilterTypes" id="FilterTypes"></a>
<h3>Types of Filters available</h3>
<p>The different types of assimilation algorithms (EAKF, ENKF,
Kernel filter, Particle filter, etc.) are determined by the
<em class="code">&amp;assim_tools_nml:filter_kind</em> entry,
described in <a href=
"../../modules/assimilation/assim_tools_mod.html">assim_tools_mod.html</a>.
Despite having 'filter' in the name, they are assimilation
algorithms and so are implemented in <em class=
"file">assim_tools_mod.f90</em>.</p>
<a name="GettingStarted" id="GettingStarted"></a>
<h3>Getting Started</h3>
<p>Running a successful assimilation takes careful diagnostic work
and experiment iterations to find the best settings for your
specific case. The basic Kalman filter can be coded in only a
handful of lines; the hard work is making the right choices to
compensate for sampling errors, model bias, observation error, lack
of model forecast divergence, variations in observation density in
space and time, random correlations, etc. There are tools built
into DART to deal with most of these problems but it takes careful
work to apply them correctly.</p>
<p>If you are adding a new model or a new observation type, we
suggest you assimilate exactly one observation, with no model
advance, with inflation turned off, with a large cutoff, and with
the outlier threshold off (see below for how to set these namelist
items). Run an assimilation. Look at the <em class=
"file">obs_seq.final</em> file to see what the forward operator
computed. Use ncdiff to difference the <em class=
"file">preassim_mean.nc</em> and <em class=
"file">postassim_mean.nc</em> (or <em class=
"file">output_mean.nc</em>) diagnostic NetCDF files and look at the
changes (the "innovations") in the various model fields. Is it in
the right location for that observation? Does it have a reasonable
value?</p>
<p>Then assimilate a group of observations and check the results
carefully. Run the observation diagnostics and look at the total
error and spread. Look carefully at the number of observations
being assimilated compared to how many are available. Assimilations
that are not working can give good looking statistics if they
reject all but the few observations that happen to match the
current state. The errors should grow as the model advances and
then shrink when new observations are assimilated, so a timeseries
plot of the RMSE should show a sawtooth pattern. The initial error
entirely depends on the match between the initial ensemble and the
observations and may be large but it should decrease and then reach
a roughly stable level. The ensemble spread should ultimately
remain relatively steady, at a value around the expected
observation error level. Once you believe you have a working
assimilation, this will be your baseline case. If the ensemble
spread is too small, several of the DART facilities described below
are intended to compensate for ensemble members getting too close
to each other. Then one by one enable or tune each of the items
below, checking each time to see what is the effect on the
results.</p>
<p>Suggestions for the most common namelist settings and features
built into DART for running a successful assimilation include:</p>

<h4>Ensemble Size</h4>
<p>In practice, ensemble sizes between 20 and 100 seem to work
best. Fewer than 20-30 members leads to statistical errors which
are too large. More than 100 members takes longer to run with very
little benefit, and eventually the results get worse again. Often
the limit on the number of members is based on the size of the
model since you have to run N copies of the model each time you
move forward in time. If you can, start with 50-60 members and then
experiment with fewer or more once you have a set of baseline
results to compare it with. The namelist setting for ensemble size
is <em class="code">&amp;filter_nml :: ens_size</em></p>

<h4>Localization</h4>
<p>There are two main advantages to using localization. One is it
avoids an observation impacting unrelated state variables because
of spurious correlations. The other is that, especially for large
models, it improves run-time performance because only points within
the localization radius need to be considered. Because of the way
the parallelization was implemented in DART, localization was easy
to add and using it usually results in a very large performance
gain. See <a href=
"../../modules/assimilation/assim_tools_mod.html#Localization">here</a>
for a discussion of localization-related namelist items.</p>

<h4>Inflation</h4>
<p>Since the filter is run with a number of members which is
usually small compared to the number of degrees of freedom of the
model (i.e. the size of the state vector or the number of EOFs
needed to characterize the variability), the model uncertainty is
under-represented. Other sources of error and uncertainty are not
represented at all. These factors lead to the ensemble being
'over-confident', or having too little spread. More observations
leads to more over-confidence. This characteristic can worsen with
time, leading to ensemble collapse to a single solution. Inflation
increases the spread of the members in a systematic way to overcome
this problem. There are several sophisticated options on inflation,
including spatial and temporal adaptive and damping options, which
help deal with observations which vary in density over time and
location. See <a href="#Inflation">here</a> for a discussion of
inflation-related namelist items.</p>

<h4>Outlier Rejection</h4>
<p>Outlier rejection can be used to avoid bad observations (ones
where the value was recorded in error or the processing has an
error and a non-physical value was generated). It also avoids
observations which have accurate values but the mean of the
ensemble members is so far from the observation value that
assimilating it would result in unacceptably large increments that
might destablize the model run. If the difference between the
observation and the prior ensemble mean is more than N standard
deviations from the square root of the sum of the prior ensemble
and observation error variance, the observation will be rejected.
The namelist setting for the number of standard deviations to
include is <em class="code">&amp;filter_nml ::
outlier_threshold</em> and we typically suggest starting with a
value of 3.0.</p>

<h4>Sampling Error</h4>
<p>For small ensemble sizes a table of expected statistical error
distributions can be generated before running DART. Corrections
accounting for these errors are applied during the assimilation to
increase the ensemble spread which can improve the assimilation
results. The namelist item to enable this option is <em class=
"code">&amp;assim_tools_nml :: sampling_error_correction</em>.
Additionally you will need to have the precomputed correction file
<em class="file">sampling_error_correction_table.nc</em>, in the
run directory. See the description of the namelist item in the
<a href=
"../../modules/assimilation/assim_tools_mod.html#Namelist">&amp;assim_tools_nml</a>
namelist, and <a href=
"../system_simulation/system_simulation.html">look here</a> for
instructions on where to find (or how to generate) the auxiliary
file needed by this code. See Anderson (2011).</p>

<h3>Free run/Forecast After Assimilation</h3>
<p>Separate scripting can be done to support forecasts starting
from the analyzed model states. After filter exits, the models can
be run freely (with no assimilated data) further forward in time
using one or more of the last updated model states from filter.
Since all ensemble members are equally likely a member can be
selected at random, or a member close to the mean can be chosen.
See the <a href=
"../../../assimilation_code/programs/closest_member_tool/closest_member_tool.html">
closest_member_tool</a> for one way to select a "close" member. The
ensemble mean is available to be used, but since it is a
combination of all the member states it may not have
self-consistent features, so using a single member is usually
preferred.</p>
<a name="EvalOnly" id="EvalOnly"></a>
<h3>Evaluating Observations Without Assimilation</h3>
<p>Filter can be used to evaluate the accuracy of a single model
state based on a set of available observations. Either copy or link
the model state file so there appear to be 2 separate ensemble
members (which are identical). Set the filter namelist ensemble
size to 2 by setting <em class="code">ens_size</em> to 2 in the
&amp;filter_nml namelist. Turn off the outlier threshold and both
Prior and Posterior inflation by setting <em class=
"code">outlier_threshold</em> to -1, and both the <em class=
"code">inf_flavor</em> values to 0 in the same &amp;filter_nml
namelist. Set all observation types to be 'evaluate-only' and have
no types in the 'assimilate' list by listing all types in the
<em class="code">evaluate_these_obs_types</em> list in the
<em class="code">&amp;obs_kind_nml</em> section of the namelist,
and none in the assimilation list. Run filter as usual, including
model advances if needed. Run observation diagnostics on the
resulting <em class="file">obs_seq.final</em> file to compute the
difference between the observed values and the predicted values
from this model state.</p>
<a name="Verify" id="Verify"></a>
<h3>Verification/Comparison With and Without Assimilation</h3>
<p>To compare results of an experiment with and without
assimilating data, do one run assimilating the observations. Then
do a second run where all the observation types are moved to the
<em class="code">evaluate_these_obs_types</em> list in the
<em class="code">&amp;obs_kind_nml</em> section of the namelist.
Also turn inflation off by setting both <em class=
"code">inf_flavor</em> values to 0 in the &amp;filter_nml namelist.
The forward operators will still be called, but they will have no
impact on the model state. Then the two sets of diagnostic state
space netcdf files can be compared to evaluate the impact of
assimilating the observations, and the observation diagnostic files
can also be compared.</p>
<a name="QCVals" id="QCVals"></a>
<h3>DART Quality Control Flag added to Output Observation Sequence
File</h3>
<p>The filter adds a quality control field with metadata 'DART
quality control' to the <em class="file">obs_seq.final</em> file.
At present, this field can have the following values:</p>
<table border="0" cellpadding="3" width="100%" summary=
'dart quality control flags'>
<tr>
<td>0:</td>
<td>Observation was assimilated successfully</td>
</tr>
<tr>
<td>1:</td>
<td>Observation was evaluated only but not used in the
assimilation</td>
</tr>
<tr>
<td>2:</td>
<td>The observation was used but one or more of the posterior
forward observation operators failed</td>
</tr>
<tr>
<td>3:</td>
<td>The observation was evaluated only but not used AND one or more
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
<td>The prior quality control value was too high so the observation
was not used.</td>
</tr>
<tr>
<td>7:</td>
<td>Outlier test failed (see below)</td>
</tr>
</table>
<p>The outlier test computes the difference between the observation
value and the prior ensemble mean. It then computes a standard
deviation by taking the square root of the sum of the observation
error variance and the prior ensemble variance for the observation.
If the difference between the ensemble mean and the observation
value is more than the specified number of standard deviations,
then the observation is not used and the DART quality control field
is set to 7.</p>
<a name="Inflation" id="Inflation"></a>
<h3>Discussion of Inflation Options</h3>
<p>In pre-Manhattan DART, there were two choices for the basic type
of inflation: observation-space or state-space. Observation-space
inflation is no longer supported. (If you are interested in
observation-space inflation, talk to Jeff first.) The rest of this
discussion applies to state-space inflation.</p>
<p>State-space inflation changes the spread of an ensemble without
changing the ensemble mean. The algorithm computes the ensemble
mean and standard deviation for each variable in the state vector
in turn, and then moves the member's values away from the mean in
such a way that the mean remains unchanged. The resulting standard
deviation is larger than before. It can be applied to the Prior
state, before observations are assimilated (the most frequently
used case), or it can be applied to the Posterior state, after
assimilation. See <a href="http://dx.doi.org/10.1175/JTECH2049.1"
target="_blank">Anderson (2007)</a>, <a href=
"http://dx.doi.org/10.1111/j.1600-0870.2008.00361.x" target=
"_blank">Anderson (2009)</a>.<br></p>
<p>Inflation values can vary in space and time, depending on the
specified namelist values. Even though we talk about a single
inflation value, the inflation has a gaussian distribution with a
mean and standard deviation. We use the mean value when we inflate,
and the standard deviation indicates how sure of the value we are.
Larger standard deviation values mean "less sure" and the inflation
value can increase more quickly with time. Smaller values mean
"more sure" and the time evolution will be slower since we are more
confident that the mean (inflation value) is correct.</p>
<p>The standard deviation of inflation allows inflation values to
increase with time, if required by increasing density or frequency
of observations, but it does not provide a mechanism to reduce the
inflation when the frequency or density of observations declines.
So there is also an option to damp inflation through time. In
practice with large geophysical models using damped inflation has
been a successful strategy.</p>
<p>The following namelist items which control inflation are found
in the <em class="file">input.nml</em> file, in the &amp;filter_nml
namelist. The detailed descriptions are in the <a href=
"../../modules/assimilation/filter_mod.html#Namelist">namelist</a>
page. Here we try to give some basic advice about commonly used
values and suggestions for where to start. Spatial variation is
controlled by <em class="code">inf_flavor</em>, which also controls
whether there's any inflation, <em class=
"code">inf_initial_from_restart</em>, and <em class=
"code">inf_initial</em>, as described below. Time variation is
controlled by <em class="code">inf_sd_initial_from_restart</em>,
<em class="code">inf_sd_initial</em>, <em class=
"code">inf_sd_lower_bound</em>, <em class="code">inf_damping</em>,
<em class="code">inf_lower_bound</em> and <em class=
"code">inf_upper_bound</em>.</p>
<p>In the namelist each entry has two values. The first is for
Prior inflation and the second is for Posterior inflation.</p>


<dl>
<dt><code>&amp;filter_nml :: inf_flavor</code></dt>
<dd>
  <details><i>valid values:</i>0, 2, 3, 4, 5</details>
  <p>Set the type of Prior and Posterior inflation applied to the
  state vector. Values mean:</p>
    
  <ul>
    <li><b>0:</b> No inflation (Prior and/or Posterior) and all other inflation
      variables are ignored</li>
    <li>[<b>1:</b> Deprecated: Observation space inflation]</li>
    <li><b>2:</b> Spatially-varying state space inflation (gaussian)</li>
    <li><b>3:</b> Spatially-uniform state space inflation (gaussian)</li>
    <li><b>4:</b> Relaxation To Prior Spread (Posterior inflation only)</li>
    <li><b>5:</b> Enhanced Spatially-varying state space inflation (inverse
      gamma)</li>
  </ul>

  <p>Spatially-varying state space inflation stores an array of inflation
  values, one for each item in the state vector. If time-evolution is enabled
  each value can evolve independently. Spatially-uniform state space inflation
  uses a single inflation value for all items in the state vector. If
  time-evolution is enabled that single value can evolve. See
  <code>inf_sd_*</code> below for control of the time-evolution behavior.
  Enhanced spatially-varying inflation uses an inverse-gamma distribution which
  allows the standard deviation of the inflation to increase or decrease
  through time and may produce better results. In practice we recommend
  starting with no inflation (both values 0). Then try inflation type 2 or 5
  prior inflation and no inflation (0) for posterior. WARNING: even if
  inf_flavor is not 0, inflation will be turned off if
  <em class="code">inf_damping</em> is set to 0.</p>
</dd>
<dt><code>&amp;filter_nml :: inf_initial_from_restart</code>
<dd>
  <details><i>valid values:</i> .true. or .false.</details>
  <p>If true, read the inflation values from an inflation restart
    file named <em class="file">input_{prior,post}inf_mean.nc.</em> An
    initial run could be done to let spatially-varying inflation values
    evolve in a spinup phase, and then the saved values can be read
    back in and used as fixed values in further runs. Or if
    time-varying inflation is used, then the restart file from the
    previous job step must be supplied as an input file for the next
    step.</p>
</dd>
<dt><code>&amp;filter_nml :: inf_initial</code></dt>
<dd>
  <details><i>valid values:</i> real numbers, usually 1.0 or slightly larger</details>
  If not reading in inflation values from a restart file, the
  initial value to set for the inflation. Generally we recommend
  starting with just slightly above 1.0, maybe 1.02, for a slight
  amount of initial inflation.
</dd>

<dt><code>&amp;filter_nml :: inf_lower_bound</code></dt>
<dd>
  <details><i>valid values:</i> real numbers, usually 1.0 or slightly larger</details>
  <p>If inflation is time-evolving (see <em class= "code">inf_sd_*</em> below),
  then this sets the lowest value the inflation can evolve to. Setting a number
  less than one allows for deflation but generally in a well-observed system
  the ensemble needs more spread and not less. We recommend a setting of 1.0.
</dd>

<dt><code>&amp;filter_nml :: inf_upper_bound</code></dt>
<dd>
  <details><i>valid values:</i> real numbers, larger than 1.0</details>
  <p>If inflation is time-evolving (see <code>inf_sd_*</code> below), then this
    sets the largest value the inflation can evolve to. We recommend a setting
    of 100.0, although if the inflation values reach those levels there is
    probably a problem with the assimilation.</p>
</dd>

<dt><code>&amp;filter_nml :: inf_damping</code></dt>
<dd>
  <details><i>valid values:</i> 0.0 to 1.0</details>
  <p>Applies to all state-space inflation types, but most frequently used with
    time-adaptive inflation variants. The difference between the current
    inflation value and 1.0 is multiplied by this factor before the next
    assimilation cycle. So the inflation values are pushed towards 1.0, from
    above or below (if inf_lower_bound allows inflation values less than 1.0).
    A value of 0.0 turns all inflation off by forcing the inflation value to
    1.0. A value of 1.0 turns damping off by leaving the original inflation
    value unchanged. We have had good results in large geophysical models
    using time- and space-adaptive state-space inflation and setting the
    damping to a value of 0.9, which damps slowly.</p>
</dd>

<dt><code>&amp;filter_nml :: inf_sd_initial_from_restart</code></dt>
<dd>
  <details><i>valid values:</i> .true. or .false.</details>
  <p>If true, read the inflation standard deviation values from an restart file
    named <code>input_{prior,post}inf_sd.nc.</code> See the comments above
    about <em class="code">inflation_initial_from_restart</em>.</p>
</dd>

<dt><code>&amp;filter_nml :: inf_sd_initial</code></dt>
<dd>
  <details><i>valid values:</i> ≤ 0.0 to disable evolution of inflation, &gt; 0.0
    otherwise</details>
  <p>The initial value to set for the inflation standard deviation,
    IF not reading in inflation standard deviation values from a file.
    This value (or these values) control whether the inflation values
    evolve with time or not. A negative value or 0.0 prevents the
    inflation values from being updated, so they are constant
    throughout the run. If positive, the inflation values evolve
    through time. We have had good results setting this and <em class=
    "code">inf_sd_lower_bound</em> to 0.6 for large geophysical
    models.</p>
</dd>

<dt><code>&amp;filter_nml :: inf_sd_lower_bound</code></dt>
<dd>
  <details><i>valid values:</i> ≤ 0.0 to disable evolution of inflation, &gt; 0.0
    otherwise</details>
  <p>If the setting of <em class="code">inf_sd_initial</em> is ≤ 0 (to disable
    time evolution of inflation) then set this to the same value.</p>
  <p>Otherwise, the standard deviation of the inflation cannot fall
  below this value. Smaller values will restrict the inflation to
  vary more slowly with time; larger values will allow the inflation
  to adapt more quickly. We have had good results setting this and
  <em class="code">inf_sd_initial</em> to 0.6 for large geophysical
  models. Since the <em class="code">inf_sd_lower_bound</em> is a
  scalar, it is not possible to set different lower bounds for
  different parts of the state vector.</p>
  <p>Time-varying inflation with flavor 2 generally results in the
  inflation standard deviation for all state variables shrinking to
  the lower bound and staying there. For flavor 5, the inflation
  standard deviation value is allowed to increase and decrease.</p>
</dd>

<dt><code>&amp;filter_nml :: inf_sd_max_change</code></dt>
<dd>
  <details><i>valid values: </i>1.0 to 2.0</details>
  <p>Used only with the Enhanced inflation (flavor 5). The Enhanced
  inflation algorithm allows the standard deviation to increase as
  well as decrease. The <em class="code">inf_sd_max_change</em>
  controls the maximum increase of the standard deviation in an
  assimilation cycle. A value of 1.0 means it will not increase, a
  value of 2.0 means it can double; a value inbetween sets the
  percentage it can increase, e.g. 1.05 is a limit of 5%. Suggested
  value is 1.05 (max increase of 5% per cycle).</p>
  <p>Because the standard deviation for original flavor 2 could never
  increase, setting the <em class="code">inf_sd_initial</em> value
  equal to the <em class="code">inf_sd_lower_bound</em> value
  effectively fixed the standard deviation at a constant value. To
  match the same behavior, if they are equal and Enhanced inflation
  (flavor 5) is used it will also use that fixed value for the
  standard deviation of the inflation. Otherwise the standard
  deviation will adapt as needed during each assimilation cycle.</p>
</dd>

<dt><code>&amp;filter_nml :: inf_deterministic</code></dt>
<dd>
  <details><i>valid values:</i> .true. or .false.</details>
  <p>Recommend always using <code>.true.</code>.</p>
</dd>
</dl>

<h3>Guidance regarding inflation</h3>

<p>The suggested procedure for testing inflation options is to
start without any (both <em class="code">inf_flavor</em> values set
to 0 and <em class="code">inf_damping</em> &gt; 0.). Then enable
Prior state space, spatially-varying inflation, with no Posterior
inflation (set <em class="code">inf_flavor</em> to [2, 0]). Then
try damped inflation (set <em class="code">inf_damping</em> to 0.9
and set <em class="code">inf_sd_initial</em> and <em class=
"code">inf_sd_lower_bound</em> to 0.6). The inflation values and
standard deviation are written out to files with <em class=
"file">_{prior,post}inf_{mean,sd}</em> in their names. These NetCDF
files can be viewed with common tools (we often use <a href=
"http://meteora.ucsd.edu/~pierce/ncview_home_page.html">ncview</a>
). Expected inflation values are generally in the 1 to 30 range; if
values grow much larger than this it usually indicates a problem
with the assimilation.</p>
<p>It is possible to set inflation values in an existing netCDF
file by using one of the standard NCO utilities like "<em class=
"program">ncap2</em>" on a copy of a restart file. Inflation mean
and sd values look exactly like restart values, arranged by
variable type like T, U, V, etc.</p>
<p>Here's an example of using ncap2 to set the T,U and V inf
values:</p>
<div class="unix">
<pre>
  ncap2 -s 'T=1.0;U=1.0;V=1.0' wrfinput_d01 input_priorinf_mean.nc
  ncap2 -s 'T=0.6;U=0.6;V=0.6' wrfinput_d01 input_priorinf_sd.nc
  -or-
  ncap2 -s 'T(:,:,:)=1.0;U(:,:,:)=1.0;V(:,:,:)=1.0' wrfinput_d01 input_priorinf_mean.nc
  ncap2 -s 'T(:,:,:)=0.6;U(:,:,:)=0.6;V(:,:,:)=0.6' wrfinput_d01 input_priorinf_sd.nc
</pre></div>
<p>Some versions of the NCO utilities change the full 3D arrays
into a single scalar. If that's your result (check your output with
<tt>ncdump -h</tt>) use the alternate syntax or a more recent
version of the NCO tools.</p>
<a name="WhereToModify" id="WhereToModify"></a>
<h3>Directories expected to be Modified</h3>
<p>DART is distributed as a toolkit/library/facility that can be
used as-is with the existing models and observations, but is also
designed so that users can add new models, new observation types
and forward operators, and new assimilation algorithms.</p>
<p>The locations in the DART <a href=
"../../../docs/index.html#Directories">code tree</a> which are
intended to be modified by users are:</p>
<dl>
<dt>New Models</dt>
<dd>Add a new directory in the <em class="file">models</em>
subdirectory. Copy (recursively, e.g. <em class="code">cp -r</em>)
the contents of the <em class="file">template</em> directory and
modify from there. Note that the <em class=
"file">model_mod.f90</em> file in the template dir is appropriate
for small models; for large geophysical models see the <em class=
"file">full_model_mod.f90</em> file and also examine other model
directories for ideas. See additional documentation in the <a href=
"../../../models/template/model_mod.html">model_mod</a>
documentation, and the <a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Documentation.php#adding_a_model">
DART web pages</a> on adding new models.</dd>
<dt>New Observation Platforms</dt>
<dd>To convert observations from other formats to DART format, add
a new directory in the <em class=
"file">observations/obs_converters</em> subdirectory and populate
it with converter code.</dd>
<dt>New Observation Types and Forward Operators</dt>
<dd>Define a new type (a measurement from an observing platform)
via a file in the <em class=
"file">observations/forward_operators</em> subdirectory. If the
forward operator is more complicated than directly interpolating a
field in the model state, this is where the code for that goes. See
additional documentation in the <a href=
"../../../observations/forward_operators/obs_def_mod.html">obs_def_mod</a>
documentation, and the <a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Observations.php#adding_types">
DART web pages</a> on adding new types. Adding a new type may
require adding a new <em class="code">generic kind</em>, which is
documented in <a href=
"../../modules/observations/obs_kind_mod.html">obs_def_mod</a>.</dd>
<dt>New Assimilation Algorithms</dt>
<dd>If you want to try out a different filter type modify the
filter code in the <em class="file">assim_tools_mod.f90</em> file.
See the <a href=
"../../modules/assimilation/assim_tools_mod.html">assim_tools_mod</a>
documentation.</dd>
</dl>
<a name="DetailedProgramFlow" id="DetailedProgramFlow"></a>
<h3>Detailed Program Execution Flow</h3>
<p>The Manhattan release of DART includes state space output
expanded from the previous two stages (Prior and Posterior) to up
to four (input, preassim, postassim, and output). This makes it
possible to examine the states with and without either kind of
inflation, as described below. In addition, the state space vectors
are each written to a separate NetCDF file: <em class=
"file">${stage}_mean.nc, ${stage}_sd.nc,
${stage}_member_####.nc</em> . The detailed execution flow inside
the filter program is:</p>
<ul>
<li>Read in observations.</li>
<li>Read in state vectors from model netcdf restart files.</li>
<li>Initialize inflation fields, possibly reading netcdf restart
files.</li>
<li>If requested, initialize and write to "input" netcdf diagnostic
files.</li>
<li>Trim off any observations if start/stop times specified.</li>
<li>Begin main assimilation loop:
<ul>
<li>Check model time vs observation times:
<ul>
<li>If current assimilation window is earlier than model time,
error.</li>
<li>If current assimilation window includes model time, begin
assimilating.</li>
<li>If current assimilation window is later than model time,
advance model:
<ul>
<li>Write out current state vectors for all ensemble members.</li>
<li>Advance the model by subroutine call or by shell script:
<ul>
<li>Tell the model to run up to the requested time.</li>
</ul>
</li>
<li>Read in new state vectors from netcdf files for all ensemble
members.</li>
</ul>
</li>
</ul>
</li>
<li>Apply prior inflation if requested.</li>
<li>Compute ensemble of prior observation values with forward
operators.</li>
<li>If requested, compute and write the "preassim" netcdf
diagnostic files. This is AFTER any prior inflation has been
applied.</li>
<li>Compute prior observation space diagnostics.</li>
<li>Assimilate all observations in this window:
<ul>
<li>Get all obs locations and kinds.</li>
<li>Get all state vector locations and kinds.</li>
<li>For each observation:
<ul>
<li>Compute the observation increments.</li>
<li>Find all other obs and states within localization radius.</li>
<li>Compute the covariance between obs and state variables.</li>
<li>Apply increments to state variables weighted by correlation
values.</li>
<li>Apply increments to any remaining unassimilated
observations.</li>
<li>Loop until all observations in window processed.</li>
</ul>
</li>
</ul>
</li>
<li>If requested, compute and write the "postassim" netcdf
diagnostic files (members, mean, spread). This is BEFORE any
posterior inflation has been applied.</li>
<li>Apply posterior inflation if requested.</li>
<li>Compute ensemble of posterior observation values with forward
operators.</li>
<li>Compute posterior observation space diagnostics.</li>
<li>If requested, compute and write out the "output" netcdf
diagnostic files (members, mean, spread). This is AFTER any
posterior inflation has been applied.</li>
<li>Loop until all observations in input file processed.</li>
</ul>
</li>
<li>Close diagnostic files.</li>
<li>Write out final observation sequence file.</li>
<li>Write out inflation restart files if requested.</li>
<li>Write out final state vectors to model restart files if
requested.</li>
<li>Release memory for state vector and observation ensemble
members.</li>
</ul>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>See the <a href=
"../../modules/assimilation/filter_mod.html#Namelist">filter
namelist</a> page for a detailed description of all <em class=
"code">&amp;filter_nml</em> variables. This namelist is read from
the file <em class="file">input.nml</em>.</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe the modules used by this program.                       -->
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
mpi_utilities_mod
filter_mod
</pre>
<p>Note that <a href=
"../../modules/assimilation/filter_mod.html#Modules">filter_mod.f90</a>
uses many more modules.</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<p>See <a href="#DetailedProgramFlow">Detailed Program Flow</a> for
a short description of DART's new 'stages'. In addition, the
Manhattan release simplifies some namelists by replacing many
user-settable file names with hardwired filenames. Files can then
be renamed in the run scripts to suit the user's needs.</p>
<ul>
<li>input ensemble member states; from <em>&amp;filter_nml ::
input_state_files</em> or <em>input_state_file_list</em></li>
<li>output ensemble member states; to <em>&amp;filter_nml ::
output_state_files</em> or <em>output_state_file_list</em></li>
<li>input observation sequence file; from <em class=
"code">&amp;filter_nml :: obs_sequence_in_name</em></li>
<li>output observation sequence file; from <em class=
"code">&amp;filter_nml :: obs_sequence_out_name</em></li>
<li>output state space diagnostics files; <em class=
"file">${stage}_mean.nc, ${stage}_sd.nc,</em> where stage =
{input,preassim,postassim,output}</li>
<li>input state space inflation data (if enabled); from <em class=
"file">input_{prior,post}inf_{mean,sd}.nc.</em></li>
<li>output state space inflation data (if enabled); to <em class=
"file">${stage}_{prior,post}inf_{mean,sd}.nc.</em>, where stage ≠
"input"</li>
<li>input.nml, to read &amp;filter_nml</li>
</ul>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>Anderson, J. L., 2001: An Ensemble Adjustment Kalman Filter for
Data Assimilation. <span style="font-style: italic;">Mon. Wea.
Rev.</span>, <span style="font-weight: bold;">129</span>,
2884-2903.<br>
<a href=
"http://dx.doi.org/10.1175/1520-0493%282001%29129%3C2884%3AAEAKFF%3E2.0.CO%3B2"
target="_blank">doi:
10.1175/1520-0493(2001)129&lt;2884:AEAKFF&gt;2.0.CO;2</a><br></li>
<li>Anderson, J. L., 2003: A Local Least Squares Framework for
Ensemble Filtering. <span style="font-style: italic;">Mon. Wea.
Rev.</span>, <span style="font-weight: bold;">131</span>,
634-642.<br>
<a href=
"http://dx.doi.org/10.1175/1520-0493%282003%29131%3C0634%3AALLSFF%3E2.0.CO%3B2"
target="_blank">doi:
10.1175/1520-0493(2003)131&lt;0634:ALLSFF&gt;2.0.CO;2</a><br></li>
<li>Anderson, J. L., 2007: An adaptive covariance inflation error
correction algorithm for ensemble filters. <span style=
"font-style: italic;">Tellus A</span>, <span style=
"font-weight: bold;">59</span>, 210-224.<br>
<a href="http://dx.doi.org/10.1111/j.1600-0870.2006.00216.x"
target="_blank">doi: 10.1111/j.1600-0870.2006.00216.x</a><br></li>
<li>Anderson, J. L., 2007: Exploring the need for localization in
ensemble data assimilation using a hierarchical ensemble filter.
<span style="font-style: italic;">Physica D</span>, <span style=
"font-weight: bold;">230</span>, 99-111.<br>
<a href="http://dx.doi.org/10.1016/j.physd.2006.02.011" target=
"_blank">doi:10.1016/j.physd.2006.02.011</a><br></li>
<li>Anderson, J., Collins, N., 2007: Scalable Implementations of
Ensemble Filter Algorithms for Data Assimilation. <span style=
"font-style: italic;">Journal of Atmospheric and Oceanic
Technology</span>, <span style="font-weight: bold;">24</span>,
1452-1463.<br>
<a href="http://dx.doi.org/10.1175/JTECH2049.1" target=
"_blank">doi: 10.1175/JTECH2049.1</a><br></li>
<li>Anderson, J. L., 2009: Spatially and temporally varying
adaptive covariance inflation for ensemble filters. <span style=
"font-style: italic;">Tellus A</span>, <span style=
"font-weight: bold;">61</span>, 72-83.<br>
<a href="http://dx.doi.org/10.1111/j.1600-0870.2008.00361.x"
target="_blank">doi: 10.1111/j.1600-0870.2008.00361.x</a><br></li>
<li>Anderson, J., T. Hoar, K. Raeder, H. Liu, N. Collins, R. Torn,
and A. Arellano, 2009: The Data Assimilation Research Testbed: A
Community Facility. <span style="font-style: italic;">Bull. Amer.
Meteor. Soc.</span>, <span style="font-weight: bold;">90</span>,
1283-1296.<br>
<a href="http://dx.doi.org/10.1175/2009BAMS2618.1" target=
"_blank">doi: 10.1175/2009BAMS2618.1</a><br></li>
<li>Anderson, J. L., 2010: A Non-Gaussian Ensemble Filter Update
for Data Assimilation. <span style="font-style: italic;">Mon. Wea.
Rev.</span>, <span style="font-weight: bold;">139</span>,
4186-4198.<br>
<a href="http://dx.doi.org/10.1175/2010MWR3253.1" target=
"_blank">doi: 10.1175/2010MWR3253.1</a><br></li>
<li>Anderson, J. L., 2011: Localization and Sampling Error
Correction in Ensemble Kalman Filter Data Assimilation. Submitted
for publication, Jan 2011. Contact author.</li>
</ul>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%"
summary='error table'>
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">ens_size in namelist is ###: Must be &gt; 1</td>
<!-- comment -->
<td valign="top">Ensemble size must be at least 2.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">inf_flavor= ### Must be 0, 2, 3.</td>
<!-- comment -->
<td valign="top">Observation Inflation is no longer supported (i.e
flavor 1).</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">Posterior observation space inflation (type 1) not
supported.</td>
<!-- comment -->
<td valign="top">Posterior observation space inflation doesn't
work.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_main</td>
<!-- message -->
<td valign="top">Number of processes &gt; model size.</td>
<!-- comment -->
<td valign="top">Number of processes can't exceed model size for
now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_generate_copy_meta_data</td>
<!-- message -->
<td valign="top">output metadata in filter needs state ensemble
size &lt; 10000, not ###.</td>
<!-- comment -->
<td valign="top">Only up to 10000 ensemble members with state
output for now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_generate_copy_meta_data</td>
<!-- message -->
<td valign="top">output metadata in filter needs obs ensemble size
&lt; 10000, not ###.</td>
<!-- comment -->
<td valign="top">Only up to 10000 ensemble members with obs space
output for now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">filter_setup_obs_sequence</td>
<!-- message -->
<td valign="top">input obs_seq file has ### qc fields; must be &lt;
2.</td>
<!-- comment -->
<td valign="top">Only 0 or 1 qc fields in input obs sequence for
now.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_obs_copy_index</td>
<!-- message -->
<td valign="top">Did not find observation copy with metadata
observation.</td>
<!-- comment -->
<td valign="top">Only 0 or 1 qc fields in input obs sequence for
now.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>Many. New assimilation algorithms, support for new observation
types, support for additional models, better performance on higher
numbers of MPI tasks... The list is long. Send email to <a href=
"mailto:dart@ucar.edu">dart@ucar.edu</a> if you are interested in
additional functionality in DART.</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
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
