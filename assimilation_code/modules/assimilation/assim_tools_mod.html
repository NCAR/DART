<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module assim_tools_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE assim_tools_mod</h1>
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
<p>This module provides subroutines that implement the parallel
versions of the sequential scalar filter algorithms. These include
the standard sequential filter as described in <a href=
"#References">Anderson 2001, 2003</a> along with systematic
correction algorithms for both mean and spread. In addition,
algorithms to do a variety of flavors of filters including the
EAKF, ENKF, particle filter, and kernel filters are included. The
parallel implementation that allows each observation to update all
state variables that are close to it at the same time is described
in <a href="#References">Anderson and Collins, 2007</a>.</p>
<a name="FilterTypes" id="FilterTypes"></a>
<h2>Filter Types</h2>
<p>Available observation space filter types include:</p>
<ul>
<li>1 = EAKF (Ensemble Adjustment Kalman Filter, see <a href=
"#References">Anderson 2001</a>)</li>
<li>2 = ENKF (Ensemble Kalman Filter)</li>
<li>3 = Kernel filter</li>
<li>4 = Observation Space Particle filter</li>
<li>5 = Random draw from posterior (contact dart@ucar.edu before
using)</li>
<li>6 = Deterministic draw from posterior with fixed kurtosis
(ditto)</li>
<li>7 = Boxcar kernel filter</li>
<li>8 = Rank Histogram filter (see <a href="#References">Anderson
2010</a>)</li>
<li>9 = Particle filter (see <a href="#References">Poterjoy
2016</a>)</li>
</ul>
<p>We recommend using type=1, the EAKF. Note that although the
algorithm is expressed in a slightly different form, the EAKF is
identical to the EnSRF (Ensemble Square Root Filter) described by
Whitaker and Hamill in 2002. Highly non-gaussian distributions may
get better results from type=8, Rank Histogram filter.</p>
<a name="Localization" id="Localization"></a>
<h2>Localization</h2>
<p><em>Localization</em> controls how far the impact of an
observation extends. The namelist items related to localization are
spread over several different individual namelists, so we have made
a single collected description of them here along with some
guidance on setting the values.</p>
<p>This discussion centers on the mechanics of how you control
localization in DART with the namelist items, and a little bit
about pragmatic approaches to picking the values. There is no
discussion about the theory behind localization - contact Jeff
Anderson for more details. Additionally, the discussion here
applies specifically to models using the 3d-sphere location module.
The same process takes place in 1d models but the details of the
location module namelist is different.</p>
<p>The following namelist items related to 3d-sphere localization
are all found in the <em class="file">input.nml</em> file:</p>

<dl>
<dt><code>&amp;assim_tools_nml :: cutoff</code></dt>
<dd>
  <details><i>valid values:</i> 0.0 to infinity</details>
  <p>This is the value, in radians, of the half-width of the
  localization radius (this follows the terminology of an early paper
  on localization). For each observation, a state vector item
  increment is computed based on the covariance values. Then a
  multiplier, based on the 'select_localization' setting (see below)
  decreases the increment as the distance between the obs and the
  state vector item increases. In all cases if the distance exceeds
  2*cutoff, the increment is 0.</p>
</dd>

<dt><code>&amp;cov_cutoff_nml :: select_localization</code></dt>
<dd>
  <details><i>valid values:</i> 1=Gaspari-Cohn; 2=Boxcar; 3=Ramped Boxcar</details>
  <p>Controls the shape of the multiplier function applied to the computed
  increment as the distance increases between the obs and the state vector item.
  Most users use type 1 localization.</p>

  <ul>
  <li>Type 1 (Gaspari-Cohn) has a value of 1 at 0 distance, 0 at
  2*cutoff, and decreases in an approximation of a gaussian in
  between.<br>
  <br></li>
  <li>Type 2 (Boxcar) is 1 from 0 to 2*cutoff, and then 0 beyond.<br>
  <br></li>
  <li>Type 3 (Ramped Boxcar) is 1 to cutoff and then ramps linearly
  down to 0 at 2*cutoff.</li>
  </ul>

  <img src= "../../../docs/images/cutoff_fig.png" alt="Shapes of Cutoff curves"
  width="100%" border="0">
</dd>

<dt><code>&amp;location_nml :: horiz_dist_only</code></dt>
<dd>
  <details><i>valid values:</i> .true., .false.</details>
  <p>If set to .true., then the vertical location of all items,
  observations and state vector both, are ignored when computing
  distances between pairs of locations. This has the effect that all
  items within a vertical-cylindrical area are considered the same
  distance away.</p>

  <p>If set to .false., then the full 3d separation is computed. Since
  the localization is computed in radians, the 2d distance is easy to
  compute but a scaling factor must be given for the vertical since
  vertical coordinates can be in meters, pressure, or model levels.
  See below for the 'vert_normalization_xxx' namelist items.</p>
</dd>

<dt><code>&amp;location_nml :: vert_normalization_{pressure,height,level,scale_height}</code></dt>
<dd>
  <details><i>valid values:</i> real numbers, in pascals, meters, index, and value
  respectively</details>

  <p>If 'horiz_dist_only' is set to .true., these are ignored. If
  set to .false., these are required. They are the amount of that
  quantity that is equivalent to 1 radian in the horizontal. If the
  model is an earth-based one, then one radian is roughly 6366
  kilometers, so if vert_normalization_height is set to 6366000
  meters, then the localization cutoff will be a perfect sphere. If
  you want to localize over a larger distance in the vertical than
  horizontal, use a larger value. If you want to localize more
  sharply in the vertical, use a smaller number. The type of
  localization used is set by which type of vertical coordinate the
  observations and state vector items have.</p>

  <p>If you have observations with different vertical coordinates (e.g.
  pressure and height), or if your observations have a different
  vertical coordinate than your state vector items, or if you want to
  localize in a different type of unit than your normal vertical
  coordinate (e.g. your model uses pressure in the vertical but you
  wish to localize in meters), then you will need to modify or add a
  <em class="code">get_close()</em> routine in your <em class=
  "file">model_mod.f90</em> file. See the discussion in the <a href=
  "../../location/threed_sphere/location_mod.html">location
  module</a> documentation for how to transform vertical coordinates
  before localization.</p>
</dd>

<dt><code>&amp;assim_tools_nml ::adaptive_localization_threshold</code></dt>
<dd>
  <details><i>valid values:</i> integer counts, or -1 to disable</details>
  <p>Used to dynamically shrink the localization cutoff in areas of
  dense observations. If set to something larger than 0, first the
  number of other observations within 2*cutoff is computed. If it is
  larger than this given threshold, the cutoff is decreased
  proportionally so if the observations were evenly distributed in
  space, the number of observations within 2*revised_cutoff would now
  be the threshold value. The cutoff value is computed for each
  observation as it is assimilated, so can be different for each
  one.</p>
</dd>

<dt><code>&amp;assim_tools_nml :: adaptive_cutoff_floor</code></dt>
<dd><details><i>valid values:</i> 0.0 to infinity, or -1 to disable</details>
  <p>If using adaptive localization (adaptive_localization_threshold
  set to a value greater than 0), then this value can be used to set
  a minimum cutoff distance below which the adaptive code will not
  shrink. Set to -1 to disable. Ignored if not using adaptive
  localization.</p>
</dd>

<dt><code>&amp;assim_tools_nml :: output_localization_diagnostics</code></dt>
<dd>
  <details><i>valid values:</i> .true., .false.</details>
  <p>If .true. and if adaptive localization is on, a single text
  line is printed to a file giving the original cutoff and number of
  observations, and the revised cutoff and new number of counts
  within this smaller cutoff for any observation which has nearby
  observations which exceed the adaptive threshold count.</p>
</dd>

<dt><code>&amp;assim_tools_nml :: localization_diagnostics_file</code></dt>
<dd>
  <details><i>valid values:</i> text string</details>
  <p>Name of the file where the adaptive localization diagnostic
  information is written.</p>
</dd>


<dt><code>&amp;assim_tools_nml :: special_localization_obs_types</code></dt>
<dd>
  <details><i>valid values:</i> list of 1 or more text strings</details>
  <p>The cutoff localization setting is less critical in DART than
  it might be in other situations since during the assimilation DART
  computes the covariances between observations and nearby state
  vector locations and that is the major factor in controlling the
  impact an observation has. For conventional observations
  fine-tuning the cutoff based on observation type is not recommended
  (it is possible to do more harm than good with it). But in certain
  special cases there may be valid reasons to want to change the
  localization cutoff distances drastically for certain kinds of
  observations. This and the following namelist items allow this.</p>

  <p>Optional list of observation types (e.g. "RADAR_REFLECTIVITY",
  "AIRS_TEMPERATURE") which will use a different cutoff distance. Any
  observation types not listed here will use the standard cutoff
  distance (set by the 'cutoff' namelist value). This is only
  implemented for the threed_sphere location module (the one used by
  most geophysical models.)</p>
</dd>

<dt><code>&amp;assim_tools_nml :: special_localization_cutoffs</code></dt>
<dd>
  <details><i>valid values:</i> list of 1 or more real values, 0.0 to infinity</details>
  <p>A list of real values, the same length as the list of
  observation types, to be used as the cutoff value for each of the
  given observation types. This is only implemented for the
  threed_sphere location module (the one used by most geophysical
  models.)</p>
</dd>
</dl>

<h3>Guidance regarding localization</h3>

<p>There are a large set of options for localization. Individual
cases may differ but in general the following guidelines might
help. Most users use the Gaspari-Cohn covariance cutoff type. The
value of the cutoff itself is the item most often changed in a
sensitivity run to pick a good general value, and then left as-is
for subsequent runs. Most localize in the vertical, but tend to use
large values so as to not disturb vertical structures. Users do not
generally use adaptive localization, unless their observations are
very dense in some areas and sparse in others.</p>

<p>The advice for setting good values for the cutoff value is to
err on the larger side - to estimate for all types of observations
under all conditions what the farthest feasible impact or
correlated structure size would be. The downsides of guessing too
large are 1) run time is slower, and 2) there can be spurious
correlations between state vector items and observations which
aren't physically related and noise can creep into the assimilation
results this way. The downside of guessing too small is that state
vector items that should get an impact from an observation won't.
This might disrupt organized features in a field and the model may
take more time to recover/reconstruct the feature.</p>

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
&amp;assim_tools_nml
   filter_kind                       = 1
   cutoff                            = 0.2
   distribute_mean                   = .false.
   sort_obs_inc                      = .false.
   spread_restoration                = .false.
   sampling_error_correction         = .false.
   adaptive_localization_threshold   = -1
   adaptive_cutoff_floor             = 0.0
   output_localization_diagnostics   = .false.
   localization_diagnostics_file     = "localization_diagnostics"
   print_every_nth_obs               = 0
   rectangular_quadrature            = .true.
   gaussian_likelihood_tails         = .false.
   close_obs_caching                 = .true.
   allow_missing_in_clm              = .false.
   adjust_obs_impact                 = .false.
   obs_impact_filename               = ""
   allow_any_impact_values           = .false.
   convert_all_obs_verticals_first   = .true.
   convert_all_state_verticals_first = .false.
   special_localization_obs_types    = 'null'
   special_localization_cutoffs      = -888888.0
  /

</pre></div>

<h3>Description of each namelist entry</h3>

<dl>
<dt><code>filter_kind</code></dt>
<dd>
<details><i>type:</i> integer</details>
<p>
  Selects the variant of filter to be used.
</p>
<ul>
  <li>1 = EAKF (Ensemble Adjustment Kalman Filter, see <a href=
  "#References">Anderson 2001</a>)</li>
  <li>2 = ENKF (Ensemble Kalman Filter)</li>
  <li>3 = Kernel filter</li>
  <li>4 = Observation Space Particle filter</li>
  <li>5 = Random draw from posterior (contact dart@ucar.edu before
  using)</li>
  <li>6 = Deterministic draw from posterior with fixed kurtosis
  (ditto)</li>
  <li>7 = Boxcar kernel filter</li>
  <li>8 = Rank Histogram filter (see <a href="#References">Anderson
  2010</a>)</li>
  <li>9 = Particle filter (see <a href="#References">Poterjoy
  2016</a>)
</ul>
<p>
  The EAKF is the most commonly used filter. Note that although the
  algorithm is expressed in a slightly different form, the EAKF is
  identical to the EnSRF (Ensemble Square Root Filter) described by
  Whitaker and Hamill in 2002.
</p>
<p>
  The Rank Histgram filter can be more successful for highly
  nongaussian distributions.
</p>
<p>
  Jon Poterjoy's Particle filter is included with this code release. To use,
  it, overwrite <code>assim_tools_mod.f90</code> with <code>assim_tools_mod.pf.f90</code>
  and rebuild filter.
</p>
<pre>
<code>
$ mv assimilation_code/modules/assimilation/assim_tools_mod.pf.f90 assimilation_code/modules/assimilation/assim_tools_mod.f90
</code>
</pre>
<p>
  There are additional namelist items in this version specific to the particle
  filter. Read the code for more details.
</p>
</dd>

<dt><code>cutoff</code></dt>
<dd>
<details><i>type:</i> real(r8)</details>
<p>
  Cutoff controls a distance dependent weight that modulates the
  impact of an observation on a state variable. The units depend both
  on the location module being used and on the covariance cutoff
  module options selected. As defined in the original paper, this is
  the half-width; the localization goes to 0 at 2 times this
  value.
</p>
</dd>

<dt><code>distribute_mean</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  If your model uses coordinates that have no options for
  different vertical coordinates then this setting has no effect on
  speed and should be .true. to use less memory. If your model has
  code to convert between different coordinate systems, for example
  Pressure, Height, Model Levels, etc, then setting this .false. will
  generally run much faster at assimilation time but will require
  more memory per MPI task. If you run out of memory, setting this to
  .true. may allow you to run but take longer.
</p>
</dd>

<dt><code>sort_obs_inc</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  If true, the final increments from obs_increment are sorted so
  that the mean increment value is as small as possible. This
  minimizes regression errors when non-deterministic filters or error
  correction algorithms are applied. HOWEVER, when using
  deterministic filters (filter_kind == 1 or 8) with no inflation or
  a combination of a determinstic filter and deterministic inflation
  (filter_nml:inf_deterministic = .TRUE.) sorting the increments is
  both unnecessary and expensive. A warning is printed to stdout and
  the log and the sorting is skipped.
</p>
</dd>

<dt><code>spread_restoration</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  True turns on algorithm to restore amount of spread that would
  be expected to be lost if underlying obs/state variable correlation
  were really 0.
</p>
</dd>

<dt><code>sampling_error_correction</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  If true, apply sampling error corrections to the correlation
  values based on the ensemble size. See <a href=
  "#References">Anderson 2012</a>. This option uses special input
  files generated by the gen_sampling_err_table tool in the
  assimilation_code/programs directory. The values are generated for
  a specific ensemble size and most common ensemble sizes have
  precomputed entries in the table. There is no dependence on which
  model is being used, only on the number of ensemble members. The
  input file must exist in the directory where the filter program is
  executing.
</p>
</dd>

<dt><code>adaptive_localization_threshold</code></dt>
<dd>
<details><i>type:</i> integer</details>
<p>
  Used to reduce the impact of observations in densely observed
  regions. If the number of observations close to a given observation
  is greater than the threshold number, the cutoff radius for
  localization is adjusted to try to make the number of observations
  close to the given observation be the threshold number. This should
  be dependent on the location module and is tuned for a
  three_dimensional spherical implementation for numerical weather
  prediction models at present.
</p>
</dd>

<dt><code>adaptive_cutoff_floor</code></dt>
<dd>
<details><i>type:</i> real</details>
<p>
  If adaptive localization is enabled and if this value is
  greater than 0, then the adaptive cutoff distance will be set to a
  value no smaller than the distance specified here. This guarentees
  a minimum cutoff value even in regions of very dense
  observations.
</p>
</dd>

<dt><code>output_localization_diagnostics</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  Setting this to <code>.true.</code> will output an additional text file
  that contains the obs key, the obs time, the obs location, the cutoff
  distance and the number of other obs which are within that radius. If
  adaptive localization is enabled, the output also contains the updated cutoff
  distance and the number of other obs within that new radius. Without adaptive
  localization there will be a text line for each observation, so this file
  could get very large. With adaptive localization enabled, there will only be
  one line per observation where the radius is changed, so the size of the file
  will depend on the number of changed cutoffs.
</p>
</dd>

<dt><code>localization_diagnostics_file</code></dt>
<dd>
<details><i>type:</i> character(len=129)</details>
<p>
  Filename for the localization diagnostics information. This
  file will be opened in append mode, so new information will be
  written at the end of any existing data.
</p>
</dd>

<dt><code>print_every_nth_obs</code></dt>
<dd>
<details><i>type:</i> integer</details>
<p>
  If set to a value <em class="code">N</em> greater than 0, the
  observation assimilation loop prints out a progress message every
  <em class="code">N</em>th observations. This can be useful to
  estimate the expected run time for a large observation file, or to
  verify progress is being made in cases with suspected
  problems.
</p>
</dd>

<dt><code>rectangular_quadrature</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  Only relevant for filter type 8 and recommended to leave <code>.true.</code>.
</p>
</dd>

<dt><code>gaussian_likelihood_tails</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  Only relevant for filter type 8 and recommended to leave <code>.false.</code>.
</p>
</dd>

<dt><code>close_obs_caching</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  Should remain .TRUE. unless you are using
  specialized_localization_cutoffs. In that case to get accurate
  results, set it to .FALSE.. This also needs to be .FALSE. if you
  have a get_close_obs() routine in your model_mod file that uses the
  types/kinds of the obs to adjust the distances.
</p>
</dd>

<dt><code>allow_missing_in_clm</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  If true, missing values (MISSING_R8 as defined in the
  types_mod.f90 file) are allowed in the state vector. Model
  interpolation routines must be written to recognize this value and
  fail the interpolation. During assimilation any state vector items
  where one or more of the ensemble members are missing will be
  skipped and their values will be unchanged by the assimilation. The
  system currently has limited support for this option; the CLM model
  has been tested and is known to work. If your model would benefit from
  setting missing values in the state vector, contact DAReS staff by emailing
  <a href="mailto:dart@ucar.edu">dart@ucar.edu</a>.
</p>
</dd>

<dt><code>adjust_obs_impact</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  If true, reads a table of observation quantities and types
  which should be artifically adjusted regardless of the actual
  correlation computed during assimilation. Setting the impact value
  to 0 prevents items from being adjusted by that class of
  observations. The input file can be constructed by the
  'obs_impact_tool' program, included in this release. See the
  documentation for more details.
</p>
</dd>

<dt><code>obs_impact_filename</code></dt>
<dd>
<details><i>type:</i> character(len=256)</details>
<p>
  If adjust_obs_impact is true, the name of the file with the
  observation types and quantities and state quantities that should
  have have an additional factor applied to the correlations during
  assimilation.
</p>
</dd>

<dt><code>allow_any_impact_values</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  If .false., then the impact values can only be zero or one (0.0
  or 1.0) - any other value will throw an error. .false. is the
  recommended setting.
</p>
</dd>

<dt><code>convert_all_obs_verticals_first</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  Should generally always be left .True.. For models without
  vertical conversion choices the setting of this item has no
  impact.
</p>
</dd>

<dt><code>convert_all_state_verticals_first</code></dt>
<dd>
<details><i>type:</i> logical</details>
<p>
  If the model has multiple choices for the vertical coordinate
  system during localization (e.g. pressure, height, etc) then this
  should be .true. if previous versions of get_state_meta_data() did
  a vertical conversion or if most of the state is going to be
  impacted by at least one observation. If only part of the state is
  going to be updated or if get_state_meta_data() never used to do
  vertical conversions, leave it .false.. The results should be the
  same but the run time may be impacted by doing unneeded conversions
  up front. For models without vertical conversion choices the
  setting of this item has no impact.
</p>
</dd>

<dt><code>special_localization_obs_types</code></dt>
<dd>
<details><i>type:</i> character(len=32), dimension(:)</details>
<p>
  Optional list of observation types (e.g. "RADAR_REFLECTIVITY",
  "RADIOSONDE_TEMPERATURE") which will use a different cutoff value
  other than the default specified by the 'cutoff' namelist. This is
  only implemented for the 'threed_sphere' locations module.
</p>
</dd>

<dt><code>special_localization_cutoffs</code></dt>
<dd>
<details><i>type:</i> real(r8), dimension(:)</details>
<p>
  Optional list of real values which must be the same length and
  in the same order as the observation types list given for the
  'special_localization_obs_types' item. These values will set a
  different cutoff distance for localization based on the type of the
  observation currently being assimilated. Any observation type not
  in the list will use the default cutoff value. This is only
  implemented for the 'threed_sphere' locations module.
</p>
</dd>
</dl>

<!--==================================================================-->
 <a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
sort_mod
random_seq_mod
obs_sequence_mod
obs_def_mod
cov_cutoff_mod
reg_factor_mod
location_mod (model dependent choice)
ensemble_manager_mod
mpi_utilities_mod
adaptive_inflate_mod
time_manager_mod
assim_model_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table summary='public interface list'>
<tr>
<td><em class="call">use assim_tools_mod, only :</em></td>
<td><a href="#filter_assim">filter_assim</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="filter_assim" id="filter_assim"></a><br>
<div class="routine"><em class="call">call filter_assim(ens_handle,
obs_ens_handle, obs_seq, keys, ens_size, num_groups, obs_val_index,
inflate, ens_mean_copy, ens_sd_copy, ens_inf_copy, ens_inf_sd_copy,
obs_key_copy, obs_global_qc_copy, obs_prior_mean_start,
obs_prior_mean_end, obs_prior_var_start, obs_prior_var_end,
inflate_only)</em>
<pre>
type(ensemble_type), intent(inout)         :: <em class=
"code">ens_handle</em>
type(ensemble_type), intent(inout)         :: <em class=
"code">obs_ens_handle</em>
type(obs_sequence_type), intent(in)        :: <em class=
"code">obs_seq</em>
integer, intent(in)                        :: <em class=
"code">keys(:)</em>
integer, intent(in)                        :: <em class=
"code">ens_size</em>
integer, intent(in)                        :: <em class=
"code">num_groups</em>
integer, intent(in)                        :: <em class=
"code">obs_val_index</em>
type(adaptive_inflate_type), intent(inout) :: <em class=
"code">inflate</em>
integer, intent(in)                        :: <em class=
"code">ens_mean_copy</em>
integer, intent(in)                        :: <em class=
"code">ens_sd_copy</em>
integer, intent(in)                        :: <em class=
"code">ens_inf_copy</em>
integer, intent(in)                        :: <em class=
"code">ens_inf_sd_copy</em>
integer, intent(in)                        :: <em class=
"code">obs_key_copy</em>
integer, intent(in)                        :: <em class=
"code">obs_global_qc_copy</em>
integer, intent(in)                        :: <em class=
"code">obs_prior_mean_start</em>
integer, intent(in)                        :: <em class=
"code">obs_prior_mean_end</em>
integer, intent(in)                        :: <em class=
"code">obs_prior_var_start</em>
integer, intent(in)                        :: <em class=
"code">obs_prior_var_end</em>
logical, intent(in)                        :: <em class=
"code">inflate_only</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Does assimilation and inflation for a set of observations that
is identified by having integer indices listed in keys. Only the
inflation is updated if inflation_only is true, otherwise the state
is also updated.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_handle</em></td>
<td>Contains state variable ensemble data and description.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_ens_handle</em></td>
<td>Contains observation prior variable ensemble and
description.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_seq</em></td>
<td>Contains the observation sequence including observed values and
error variances.</td>
</tr>
<tr>
<td valign="top"><em class="code">keys</em></td>
<td>A list of integer indices of observations in obs_seq that are
to be used at this time.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size</em></td>
<td>Number of ensemble members in state and observation prior
ensembles.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_groups</em></td>
<td>Number of groups being used in assimilation.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_val_index</em></td>
<td>Integer index of copy in obs_seq that contains the observed
value from instrument.</td>
</tr>
<tr>
<td valign="top"><em class="code">inflate</em></td>
<td>Contains inflation values and all information about inflation
to be used.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_mean_copy</em></td>
<td>Index of copy containing ensemble mean in ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_sd_copy</em></td>
<td>Index of copy containing ensemble standard deviation in
ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_inf_copy</em></td>
<td>Index of copy containing state space inflation in
ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_inf_sd_copy</em></td>
<td>Index of copy containing state space inflation standard
deviation in ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_key_copy</em></td>
<td>Index of copy containing unique key for observation in
obs_ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_global_qc_copy</em></td>
<td>Index of copy containing global quality control value in
obs_ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">obs_prior_mean_start   </em></td>
<td>Index of copy containing first group's prior mean in
obs_ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_prior_mean_end</em></td>
<td>Index of copy containing last group's prior mean in
obs_ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_prior_var_start</em></td>
<td>Index of copy containing first group's ensemble variance in
obs_ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_prior_var_end</em></td>
<td>Index of copy containing last group's ensemble variance in
obs_ens_handle.</td>
</tr>
<tr>
<td valign="top"><em class="code">inflate_only</em></td>
<td>True if only inflation is to be updated, and not state.</td>
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
<table border="0" summary='files used by this module'>
<tr>
<th>filename</th>
<th>purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read <em class="code">assim_tools_nml</em></td>
</tr>
</table>
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
10.1175/1520-0493(2001)129&lt;2884:AEAKFF&gt;2.0.CO;2</a><br>
<br></li>
<li>Anderson, J. L., 2003: A Local Least Squares Framework for
Ensemble Filtering. <span style="font-style: italic;">Mon. Wea.
Rev.</span>, <span style="font-weight: bold;">131</span>,
634-642.<br>
<a href=
"http://dx.doi.org/10.1175/1520-0493%282003%29131%3C0634%3AALLSFF%3E2.0.CO%3B2"
target="_blank">doi:
10.1175/1520-0493(2003)131&lt;0634:ALLSFF&gt;2.0.CO;2</a><br>
<br></li>
<li>Anderson, J., Collins, N., 2007: Scalable Implementations of
Ensemble Filter Algorithms for Data Assimilation. <span style=
"font-style: italic;">Journal of Atmospheric and Oceanic
Technology</span>, <span style="font-weight: bold;">24</span>,
1452-1463.<br>
<a href="http://dx.doi.org/10.1175/JTECH2049.1" target=
"_blank">doi: 10.1175/JTECH2049.1</a><br>
<br></li>
<li>Anderson, J. L., 2010: A Non-Gaussian Ensemble Filter Update
for Data Assimilation. <span style="font-style: italic;">Mon. Wea.
Rev.</span>, <span style="font-weight: bold;">139</span>,
4186-4198.<br>
<a href="http://dx.doi.org/10.1175/2010MWR3253.1" target=
"_blank">doi: 10.1175/2010MWR3253.1</a><br>
<br></li>
<li>Anderson, J. L., 2012:, Localization and Sampling Error
Correction in Ensemble Kalman Filter Data Assimilation.
<span style="font-style: italic;">Mon. Wea. Rev.</span>,
<span style="font-weight: bold;">140</span>, 2359-2371.<br>
<a href="http://dx.doi.org/10.1175/MWR-D-11-00013.1" target=
"_blank">doi: 10.1175/MWR-D-11-00013.1</a><br>
<br></li>
<li>Poterjoy, J., 2016:, A localized particle filter for
high-dimensional nonlinear systems. <span style=
"font-style: italic;">Mon. Wea. Rev.</span> <span style=
"font-weight: bold;">144</span> 59-76.<br>
<a href="http://dx.doi.org/10.1175/MWR-D-15-0163.1" target=
"_blank">doi:10.1175/MWR-D-15-0163.1</a><br>
<br></li>
</ul>
<br>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
 <a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%"
summary='error codes and conditions'>
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">assim_tools_init</td>
<!-- message -->
<td valign="top">cant combine spread_restoration and filter_kind
###</td>
<!-- comment -->
<td valign="top">Spread restoration only compatible with
filter_kind 1</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_increment</td>
<!-- message -->
<td valign="top">Illegal value of filter_kind in assim_tools
namelist [1-8 OK]</td>
<!-- comment -->
<td valign="top">Only 1-4,8 currently supported.</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_increment_eakf<br>
obs_increment_ran_kf<br>
obs_increment_det_kf</td>
<!-- message -->
<td valign="top">Both obs_var and prior_var are zero. This is
inconsistent</td>
<!-- comment -->
<td valign="top">Product of two delta functions doesn't work.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_correction_from_file</td>
<!-- message -->
<td valign="top">Use less than 10000 ensemble</td>
<!-- comment -->
<td valign="top">File only works up to 9999 members.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_correction_from_file</td>
<!-- message -->
<td valign="top">Correction file ______ does not exist.</td>
<!-- comment -->
<td valign="top">Couldn't find the correction file in the working
directory.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>The current version of the systematic error correction algorithm
does not work in a logical fashion with the parallel region
assimilation feature.</p>
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
