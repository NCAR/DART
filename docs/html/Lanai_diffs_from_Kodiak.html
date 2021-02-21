<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>DART Lanai Differences from Kodiak Release Notes</title>
<link rel="stylesheet" type="text/css" href="doc.css">
<link href="../images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>DART Lanai Differences from Kodiak Release Notes</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../images/Dartboard7.png" alt=
"DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Changes">Changes to Core DART routines</a> / <a href=
"#NewModels">New Models or Changes to Existing Models</a> /
<a href="#ForwardOps">New or Changed Forward Operators</a> /
<a href="#ObsConvert">Observation Converters</a> / <a href=
"#Diagnostics">New or Updated DART Diagnostics</a> / <a href=
"#Misc">Tutorial, Scripting, Setup, Builds</a> / <a href=
"#Legalese">Terms of Use</a> 
<!--==================================================================-->
<a name="Overview" id="Overview"></a>
<h2>Overview</h2>
<p>This document includes an overview of the changes in the DART
system since the Kodiak release. For further details on any of
these items look at the HTML documentation for that specific part
of the system.</p>
<p>There is a longer companion document for this release, the
<a href="Lanai_release.html">Lanai Release Notes</a>, which include
installation instructions, a walk-through of running one of the
low-order models, the diagnostics, and a description of
non-backward compatible changes. See the <a href=
"Lanai_release.html#CurrentUsers">Notes for Current Users</a>
section for additional information on changes in this release.</p>
<!--==================================================================-->
<a name="Changes" id="Changes"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Changes to Core DART routines</h2>
<p>This section describes changes in the basic DART library
routines since the Kodiak release.</p>
<ul>
<li>Added a completely new random number generator based on the
Mersenne Twister algorithm from the GNU scientific library. It
seems to have better behavior if reseeded frequently, which is a
possible usage pattern if perfect_model_obs is run for only single
steps and the model is advanced in an external script. As part of
this code update all random number code was moved into the
random_seq_mod and random_nr_mod is deprecated.</li>
<li>Perfect_model_obs calls a seed routine in the time manager now
that generates a consistent seed based on the current time of the
state. This makes subsequent runs give consistent results and yet
separate runs don't get identical error values.</li>
<li>Added random number generator seeds in several routines to try
to get consistent results no matter how many MPI tasks the code was
run with. This includes:
<ul>
<li>cam model_mod.f90, pert_model_state()</li>
<li>assim_tools_mod.f90, filter_assim(), filter kinds 2, 3, and
5</li>
<li>wrf model_mod.f90, pert_model_state()</li>
<li>adaptive_inflate_mod.f90, adaptive_inflate_init(),
non-deterministic inf</li>
</ul>
</li>
<li>There is a new &amp;filter_nml namelist item:
enable_special_outlier_code. If .true. the DART quality control
code will call a separate subroutine at the end of filter.f90 to
evaluate the outlier threshold. The user can add code to that
routine to change the threshold based on observation type or values
as they wish. If .false. the default filter outlier threshold code
will be called and the user routine ignored.</li>
<li>If your <em class="file">model_mod.f90</em> provides a
customized <em class="code">get_close_obs()</em> routine that makes
use of the types/kinds arguments for either the base location or
the close location list, there is an important change in this
release. The fifth argument to the <em class=
"code">get_close_obs()</em> call is now a list of generic kinds
corresponding to the location list. The fourth argument to the
<em class="code">get_dist()</em> routine is now also a generic kind
and not a specific type. In previous versions of the system the
list of close locations was sometimes a list of specific types and
other times a list of generic kinds. The system now always passes
generic kinds for the close locations list for consistency. The
base location and specific type remains the same as before. If you
have a <em class="code">get_close_obs()</em> routine in your
<em class="file">model_mod.f90</em> file and have questions about
usage, <a href="mailto:dart@ucar.edu">contact</a> the DART
development team.</li>
<li>Filter will call the end_model() subroutine in the model_mod
for the first time. It should have been called all along, but was
not.</li>
<li>Added a time sort routine in the time_manager_mod.</li>
<li>Avoid a pair of all-to-all transposes when setting the
inflation mean and sd from the namelist. The new code finds the
task which has the two copies and sets them directly without a
transpose. The log messages were also moved to the end of the
routine - if you read in the mean/sd values from a restart file the
log messages that printed out the min/max values needed to be after
the read from the file.</li>
<li>Reordered the send/receive loops in the all-to-all transposes
to scale better on yellowstone.</li>
<li>Remove a state-vector size array from the stack in
read_ensemble_restart(). The array is now allocated only if needed
and then deallocated. The ensemble write routine was changed before
the Kodiak release but the same code in read was apparently not
changed simply as an oversight.</li>
<li>If the ensemble mean is selected to be written out in dart
restart file format, the date might not have been updated
correctly. The code was fixed to ensure the ensemble mean date in
the file was correct.</li>
<li>filter writes the ensemble size into the log file.</li>
<li>Reorganized the code in the section of obs_model_mod that
prints out the time windows, with and without verbose details.
Should be clearer if the next observation is in or out of the
current assimilation window, and if the model needs to advance or
not.</li>
<li>Added a fill_inflation_restart utility which can write a file
with a fixed mean and sd, so the first step of a long assimilation
run can use the same 'start_from_restart_file' as subsequent
steps.</li>
<li>Added new location module options:
<ul>
<li>Channel coordinate system</li>
<li>[0-1] periodic 3D coordinate system</li>
<li>X,Y,Z 3D Cartesian coordinate system</li>
<li>2D annulus coordinate system</li>
</ul>
</li>
</ul>
<!--==================================================================-->
<br>
<br>
<a name="NewModels" id="NewModels"></a>
<hr>
<h2>New Models or Changes to Existing Models</h2>
<p>Several new models have been incorporated into DART. This
section details both changes to existing models and descriptions of
new models that have been added since the Kodiak release.</p>
<ul>
<li>Support for components under the CESM framework:
<ul>
<li>Added support for the Community Land Model (CLM).</li>
<li>Added support to run the Community Atmospheric Model (CAM)
under the CESM framework.</li>
<li>Added support for the CESM 1.1.1 release for CAM, POP, CLM:
includes experiment setup scripts and assimilation scripts.</li>
<li>CAM, POP, and/or CLM can be assimilated either individually or
in combination while running under the CESM framework. If
assimilating into multiple components, they are assimilated
sequentially with observations only affecting a single component
directly. Other components are indirectly affected through
interactions with the coupler.</li>
<li>Setup scripts are provided to configure a CESM experiment using
the multi-instance feature of CESM to support ensembles for
assimilation.</li>
<li>POP state vector contains potential temperature; observations
from the World Ocean Database are in-situ or sensible temperature.
The model_mod now corrects for this.
<ul>
<li>The state vector has all along contained potential temperature
and not in-situ (sensible) temperature. The observations from the
World Ocean Database are of sensible temperature. Changed the
specific kind in the model_mod to be
<tt>QTY_POTENTIAL_TEMPERATURE</tt> and added new code to convert
from potential to in-situ temperature. Differences for even the
deeper obs (4-5km) is still small ( ~ 0.2 degree). (in-situ or
sensible temperature is what you measure with a regular
thermometer.)</li>
</ul>
</li>
<li>Support for the SE core (HOMME) of CAM has been developed but
<strong>is not</strong> part of the current release. Contact the
DART group if you have an interest in running this configuration of
CAM.</li>
</ul>
</li>
<li>Changes to the WRF model_mod:
<ul>
<li>Allow advanced microphysics schemes (needed interpolation for 7
new kinds)</li>
<li>Interpolation in the vertical is done in log(p) instead of
linear pressure space. log(p) is the default, but a compile-time
variable can restore the linear interpolation.</li>
<li>Added support in the namelist to avoid writing updated fields
back into the wrf netcdf files. The fields are still updated during
the assimilation but the updated data is not written back to the
wrfinput file during the dart_to_wrf step.</li>
<li>Fixed an obscure bug in the vertical convert routine of the wrf
model_mod that would occasionally fail to convert an obs. This
would make tiny differences in the output as the number of mpi
tasks change. No quantitative differences in the results but they
were not bitwise compatible before and they are again now.</li>
</ul>
</li>
<li>Added support for the MPAS_ATM and MPAS_OCN models.
<ul>
<li>Added interpolation routines for the voroni-tesselation grid
(roughly hexagonal)</li>
<li>Includes vertical conversion routines for vertical
localization.</li>
<li>Added code to the mpas_atm model to interpolate specific
humidity and pressure, so we can assimilate GPS obs now.</li>
</ul>
</li>
<li>Added support for the 'SQG' uniform PV two-surface QC+1
spectral model.</li>
<li>Added support for a flux-transport solar dynamo model.</li>
<li>Added support for the GITM upper atmosphere model.</li>
<li>Added support for the NOAH land model.</li>
<li>Added support for the NAAPS model.</li>
<li>Added model_mod interface code for the NOGAPS model to the SVN
repository.</li>
<li>Simple advection model:
<ul>
<li>Fix where the random number seed is set in the
models/simple_advection model_mod - it needed to be sooner than it
was being called.</li>
</ul>
</li>
</ul>
<!--==================================================================-->
<a name="ForwardOps" id="ForwardOps"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>New or changed Forward Operators</h2>
<p>This section describes changes to the Foward Operators and new
Generic Kinds or Specific Types that have been added since the
Kodiak release.</p>
<ul>
<li>Many new kinds added to the DEFAULT_obs_kind_mod.f90:
<ul>
<li>QTY_CANOPY_WATER</li>
<li>QTY_CARBON</li>
<li>QTY_CLW_PATH</li>
<li>QTY_DIFFERENTIAL_REFLECTIVITY</li>
<li>QTY_DUST</li>
<li>QTY_EDGE_NORMAL_SPEED</li>
<li>QTY_FLASH_RATE_2D</li>
<li>QTY_GRAUPEL_VOLUME</li>
<li>QTY_GROUND_HEAT_FLUX QTY_HAIL_MIXING_RATIO</li>
<li>QTY_HAIL_NUMBER_CONCENTR</li>
<li>QTY_HAIL_VOLUME QTY_ICE QTY_INTEGRATED_AOD</li>
<li>QTY_INTEGRATED_DUST</li>
<li>QTY_INTEGRATED_SEASALT QTY_INTEGRATED_SMOKE</li>
<li>QTY_INTEGRATED_SULFATE</li>
<li>QTY_LATENT_HEAT_FLUX</li>
<li>QTY_LEAF_AREA_INDEX</li>
<li>QTY_LEAF_CARBON</li>
<li>QTY_LEAF_NITROGEN QTY_LIQUID_WATER</li>
<li>QTY_MICROWAVE_BRIGHT_TEMP</li>
<li>QTY_NET_CARBON_FLUX</li>
<li>QTY_NET_CARBON_PRODUCTION</li>
<li>QTY_NEUTRON_INTENSITY</li>
<li>QTY_NITROGEN QTY_RADIATION</li>
<li>QTY_ROOT_CARBON</li>
<li>QTY_ROOT_NITROGEN</li>
<li>QTY_SEASALT</li>
<li>QTY_SENSIBLE_HEAT_FLUX</li>
<li>QTY_SMOKE</li>
<li>QTY_SNOWCOVER_FRAC</li>
<li>QTY_SNOW_THICKNESS</li>
<li>QTY_SNOW_WATER</li>
<li>QTY_SO2</li>
<li>QTY_SOIL_CARBON</li>
<li>QTY_SOIL_NITROGEN</li>
<li>QTY_SPECIFIC_DIFFERENTIAL_PHASE</li>
<li>QTY_STEM_CARBON</li>
<li>QTY_STEM_NITROGEN</li>
<li>QTY_SULFATE</li>
<li>QTY_VORTEX_WMAX</li>
<li>QTY_WATER_TABLE_DEPTH</li>
<li>QTY_WIND_TURBINE_POWER</li>
<li>plus slots 151-250 reserved for Chemistry (specifically
WRF-Chem) kinds</li>
</ul>
</li>
<li>Added a forward operator for total precipitable water. It loops
over model levels so it can be used as an example of how to handle
this without having to hardcode the number of levels into the
operator.</li>
<li>Added a forward operator (and obs_seq file converter) for
COSMOS ground moisture observations.</li>
<li>Added a forward operator (and obs_seq file converter) for MIDAS
observations of Total Electron Count.</li>
<li>Added a 'set_1d_integral()' routine to the
obs_def_1d_state_mod.f90 forward operator for the low order models.
This subroutine isn't used by filter but it would be needed if
someone wanted to write a standalone program to generate obs of
this type. We use this file as an example of how to write an obs
type that has metadata, but we need to give an example of how to
set the metadata if you aren't using create_obs_sequence
interactively (e.g. your data is in netcdf and you have a separate
converter program.)</li>
</ul>
<!--==================================================================-->
<a name="ObsConvert" id="ObsConvert"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Observation Converters</h2>
<p>This section describes support for new observation types or
sources that have been added since the Kodiak release.</p>
<ul>
<li>Added an obs_sequence converter for wind profiler data from
MADIS.</li>
<li>Added an obs_sequence converter for Ameriflux land
observations(latent heat flux, sensible heat flux, net ecosystem
production).</li>
<li>Added an obs_sequence converter for MODIS snow coverage
measurements.</li>
<li>Added an obs_sequence converter for COSMOS ground moisture
observations.</li>
<li>Added an obs_sequence converter for MIDAS observations of Total
Electron Count.</li>
<li>Updated scripts for the GPS converter; added options to convert
data from multiple satellites.</li>
<li>More scripting support in the MADIS obs converters; more error
checks added to the rawin converter.</li>
<li>Added processing for wind profiler observation to the
wrf_dart_obs_preprocess program.</li>
<li>Fix BUG in airs converter - the humidity obs are accumulated
across the layers and so the best location for them is the layer
midpoint and not on the edges (levels) as the temperature obs are.
Also fixed off-by-one error where the converter would make one more
obs above the requested top level.</li>
<li>Made gts_to_dart converter create separate obs types for
surface dewpoint vs obs aloft because they have different vertical
coordinates.</li>
<li>Converted mss commands to hpss commands for a couple
observation converter shell scripts (inc AIRS).</li>
<li>New matlab code to generate evenly spaced observations on the
surface of a sphere (e.g. the globe).</li>
<li>Added obs_loop.f90 example file in obs_sequence directory;
example template for how to construct special purpose obs_sequence
tools.</li>
<li>Change the default in the script for the prepbufr converter so
it will swap bytes, since all machines except ibms will need this
now.</li>
<li>The 'wrf_dart_obs_preprocess' program now refuses to superob
observations that include the pole, since the simple averaging of
latitude and longitude that works everyplace else won't work there.
Also treats observations near the prime meridian more
correctly.</li>
</ul>
<!--==================================================================-->
<a name="Diagnostics" id="Diagnostics"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>New or updated DART Diagnostics</h2>
<p>This section describes new or updated diagnostic routines that
have been added since the Kodiak release.</p>
<ul>
<li>Handle empty epochs in the obs_seq_to_netcdf converter.</li>
<li>Added a matlab utility to show the output of a 'hop' test
(running a model for a continuous period vs. stopping and
restarting a run).</li>
<li>Improved the routine that computes axes tick values in plots
with multiple values plotted on the same plot.</li>
<li>The obs_common_subset program can select common observations
from up to 4 observation sequence files at a time.</li>
<li>Add code in obs_seq_verify to ensure that the ensemble members
are in the same order in all netcdf files.</li>
<li>Added support for the unstructured grids of mpas to our matlab
diagnostics.</li>
<li>Fix to writing of ReportTime in obs_seq_coverage.</li>
<li>Fixed logic in obs_seq_verify when determining the forecast
lat.</li>
<li>Fixed loops inside obs_seq_coverage which were using the wrong
limits on the loops. Fixed writing of 'ntimes' in output netcdf
variable.</li>
<li>The obs_common_subset tool supports comparing more than 2
obs_seq.final files at a time, and will loop over sets of
files.</li>
<li>Rewrote the algorithm in the obs_selection tool so it had
better scaling with large numbers of obs.</li>
<li>Several improvements to the 'obs_diag' program:
<ul>
<li>Added preliminary support for a list of 'trusted obs' in the
obs_diag program.</li>
<li>Can disable the rank histogram generation with a namelist
item.</li>
<li>Can define height_edges or heights in the namelist, but not
both.</li>
<li>The 'rat_cri' namelist item (critical ratio) has been
deprecated.</li>
</ul>
</li>
<li>Extend obs_seq_verify so it can be used for forecasts from a
single member. minor changes to obs_selection, obs_seq_coverage and
obs_seq_verify to support a single member.</li>
<li>Added Matlab script to read/print timestamps from binary dart
restart/ic files.</li>
<li>Default for obs_seq_to_netcdf in all the namelists is now 'one
big time bin' so you don't have to know the exact timespan of an
obs_seq.final file before converting to netCDF.</li>
</ul>
<!--==================================================================-->
<a name="Misc" id="Misc"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Tutorial, Scripting, Setup, Builds</h2>
<p>This section describes updates and changes to the tutorial
materials, scripting, setup, and build information since the Kodiak
release.</p>
<ul>
<li>The mkmf-generated Makefiles now take care of calling
'fixsystem' if needed so the mpi utilities code compiles without
further user intervention.</li>
<li>Make the default input.nml for the Lorenz 96 and Lorenz 63
model gives good assimilation results. Rename the original
input.nml to input.workshop.nml. The workshop_setup script renames
it back before doing anything else so this won't break the workshop
instructions. Simplify all the workshop_setup.csh scripts to do the
minimal work needed by the DART tutorial.</li>
<li>Updates to the models/template directory with the start of a
full 3d geophysical model template. Still under construction.</li>
<li>Move the pdf files in the tutorial directory up a level.
Removed framemaker source files because we no longer have access to
a working version of the Framemaker software. Moved routines that
generate figures and diagrams to a non-distributed directory of the
subversion repository.</li>
<li>Enable netCDF large file support in the work/input.nml for
models which are likely to have large state vectors.</li>
<li>Minor updates to the doc.css file, make pages look identical in
the safari and firefox browsers.</li>
<li>Added a utility that sorts and reformats namelists, culls all
comments to the bottom of the file. Useful for doing diffs and
finding duplicated namelists in a file.</li>
<li>Cleaned up mkmf files - removed files for obsolete platforms
and compilers, updated suggested default flags for intel.</li>
<li>Update the mkmf template for gfortran to allow fortran source
lines longer than 132 characters.</li>
</ul>
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
