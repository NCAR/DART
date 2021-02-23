<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>DART Manhattan Release Notes</title>
<link rel="stylesheet" type="text/css" href="doc.css">
<link href="../images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>Manhattan</h1>
<h2>DART Manhattan release documentation</h2>
<h2>DART Overview</h2>
<p>The Data Assimilation Research Testbed (DART) is designed to
facilitate the combination of assimilation algorithms, models, and
real (or synthetic) observations to allow increased understanding
of all three. The DART programs are highly portable, having been
compiled with many Fortran 90 compilers and run on linux
compute-servers, linux clusters, OSX laptops/desktops, SGI Altix
clusters, supercomputers running AIX, and more. Read the <a href=
"https://www.image.ucar.edu/DAReS/DART/DART2_Starting.php#customizations">
Customizations</a> section for help in building on new
platforms.</p>
<p>DART employs a modular programming approach to apply an Ensemble
Kalman Filter which adjusts model values toward a state that is
more consistent with information from a set of observations. Models
may be swapped in and out, as can different algorithms in the
Ensemble Kalman Filter. The method requires running multiple
instances of a model to generate an ensemble of states. A forward
operator appropriate for the type of observation being assimilated
is applied to each of the states to generate the model's estimate
of the observation. Comparing these estimates and their uncertainty
to the observation and its uncertainty ultimately results in the
adjustments to the model states. See the <a href=
"../DART_LAB/DART_LAB.html">DART_LAB</a> demos or read more
<a href="../tutorial/index.html">in the DART tutorial</a>.</p>
<p>DART diagnostic output can be written that contains the model
state before and after the adjustment, along with the ensemble mean
and standard deviation, and prior or posterior inflation values if
inflation is enabled. There is also a text file, <em class=
"file">obs_seq.final</em>, with the model estimates of the
observations. There is a suite of MATLAB® functions that facilitate
exploration of the results, but the netCDF files are inherently
portable and contain all the necessary metadata to interpret the
contents with other analysis programs such as NCL, R, etc.</p>
<p>To get started running with Lorenz 63 model refer to <a href=
"Manhattan_getting_started.html">Getting Started</a>.</p>

<!--==================================================================-->
<a name="CurrentUsers" id="CurrentUsers"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Notes for Current Users</h2>
<p>If you have been updating from the rma_trunk branch of the DART
subversion repository you will notice that the code tree has been
simplified to be more intuitive for users. The new top level
directory structure looks like :</p>
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
if you do try to do an 'svn update' on an existing directory, you
will encounter many 'tree conflicts'.
<p>We suggest that current users checkout a fresh version of
Manhattan in a new location. To see which files need to be moved,
run 'svn status' on your original checked out version. Anything
with an M or ? in the first column needs to be moved to the new
location in the new tree. Please <a href=
"mailto:dart@ucar.edu">contact</a> DART if you have any issues
migrating your existing code to the new tree structure.</p>
<p>There is a list of non-backwards compatible changes (<a href=
"#Nonbackward">see below</a>), and a list of new options and
functions.</p>
<p>The Manhattan release will continue to be updated for the next
few months as we continue to add features. Checking out the
Manhattan release branch and running 'svn update' from time to time
is the recommended way to update your DART tree.</p>
<!--==================================================================-->
<a name="Nonbackward" id="Nonbackward"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Non-backwards Compatible Changes</h2>
<p>Unlike previous releases of DART, this version contains more
non-backwards compatible changes than usual. Please examine the
following list carefully. We do suggest you check out the Manhattan
release into a new location and migrate any local changes from
previous versions as a second step.</p>
<p>Changes in the Manhattan release (15 May 2015) which are
<em>not</em> backwards compatible with the Lanai release (13 Dec
2013):</p>
<ol>
<li>We no longer require model data to be converted to DART format
restart files. We directly read and write NetCDF format only. To
specify the input and output files for filter, there are new
namelist items in the &amp;filter_nml namelist:
<code>'input_state_file_list'</code> and
<code>'output_state_file_list'</code> .</li>
<li>The information formerly in <em class="file">Prior_Diag.nc</em>
and <em class="file">Posterior_Diag.nc</em> has been moved. If you
are reading and writing ensemble members from different files, the
state information, the ensemble mean and standard deviation, and
the inflation mean and standard deviation will all be read and
written to separate files:
<ul>
<li><em class="file">[stage]_member_####.nc</em></li>
<li><em class="file">[stage]_mean.nc</em></li>
<li><em class="file">[stage]_sd.nc</em></li>
<li><em class="file">[stage]_priorinf_{mean,sd}.nc</em> (if prior
inflation is turned on)</li>
<li><em class="file">[stage]_postinf_{mean,sd}.nc</em> (if
posterior inflation is turned on)</li>
</ul>
<br>
If you are reading and writing ensemble members from a single file,
all this information will now be in a single NetCDF file but will
be stored in different variables inside that file:
<ul>
<li><em class="file">[var].nc</em></li>
<li><em class="file">[var]_mean.nc</em></li>
<li><em class="file">[var]_sd.nc</em></li>
<li><em class="file">[var]_priorinf_{mean,sd}.nc</em> (if prior
inflation is turned on)</li>
<li><em class="file">[var]_postinf_{mean,sd}.nc</em> (if posterior
inflation is turned on)</li>
</ul>
<br>
We also now have options for writing files at four stages of the
assimilation cycle: <code>'input', 'preassim', 'postassim',
'output'</code>. This is set in the &amp;filter_nml namelist with
stages_to_write.</li>
<li>New model_mod.f90 required routines:
<ul>
<li><em class="code">vert_convert()</em></li>
<li><em class="code">query_vert_localization_coord()</em></li>
<li><em class="code">pert_model_copies()</em></li>
<li><em class="code">read_model_time()</em></li>
<li><em class="code">write_model_time()</em></li>
</ul>
There are default version of these available to use if you have no
special requirements.</li>
<li>Several of the model_mod.f90 argument lists have changed
<ul>
<li><em class="code">model_interpolate()</em> now takes in the
<code>state_handle</code> as an argument rather than a state vector
array. It also return an array of <code>expected_obs</code> and
<code>istatus</code> for each of the ensemble members</li>
<li><em class="code">get_state_meta_data()</em> also requires the
<code>state_handle</code> as an argument rather than a state vector
array.</li>
<li><em class="code">nc_write_model_atts()</em> has an additional
argument <code>moel_mod_writes_state_variables</code>. If true then
the model_mod is expected to write out the state variables, if
false DART will write out the state variable (this is the prefered
method for adding new models, it requires less code from the model
developer)</li>
</ul>
</li>
<li>There are several namelist changes mainly in the
&amp;filter_nml and &amp;perfect_model_mod which are outlined in
detail in <a href=
"Manhattan_diffs_from_Lanai.html">Manhattan_diffs_from_Lanai</a></li>
<li>All modules have been moved to <em class=
"dir">DART/assimilation_code/modules/</em> directory. And similarly
all of the programs have moved to <em class=
"dir">DART/assimilation_code/programs/</em></li>
<li>The location modules which were stored in <em class=
"dir">locations</em> have moved to <em class=
"dir">DART/assimilation_code/location</em> directory</li>
<li>The observation converters which were stored in <em class=
"dir">observations</em> have moved to <em class=
"dir">DART/observations/obs_converters</em> directory</li>
<li>The forward operators have moved from <em class=
"dir">obs_def/obs_def_*_mod.f90</em> to <em class=
"dir">observations/forward_operators</em></li>
<li>The tutorial files have moved to <em class=
"dir">DART/docs/tutorial directory</em></li>
<li>The program <em class="file">fill_inflation_restart</em> is
OBSOLETE since DART inflation files are now in NetCDF format. Now
inflation files can be filled using <em class="program">ncap2</em>.
Here is an example using version 4.4.2 or later of the NCO tools:
<pre>
  ncap2 -s "T=1.0;U=1.0;V=1.0" wrfinput_d01 prior_inf.nc'
  ncap2 -s "T=0.6;U=0.6;V=0.6" wrfinput_d01 prior_sd.nc'
</pre></li>
<li>The default flags in the mkmf_template.XXX files have been
updated to be more consistent with current compiler versions.</li>
<li>If you enable the sampling error correction option, the
required data is now read from a single netcdf file which supports
multiple ensemble sizes. A program is provided to compute
additional ensemble sizes if they are not in the default file.</li>
<li>
<p>Our use of TYPES and KINDS has been very confusing in the past.
In Manhattan we have tried to make it clearer which things in DART
are generic quantities (QTY) - temperature, pressure, etc - and
which things are specific types of observations -
Radiosonde_temperature, Argo_salinity etc.</p>
<p>Below is a mapping between old and new subroutine names here for
reference. We have made these changes to all files distributed with
DART. If you have lots of code developed outside of the subversion
repository, please contact <a href="mailto:dart@ucar.edu">DART</a>
for a sed script to help automate the changes.</p>
Public subroutines, existing name on left, replacement on right:
<pre>
    
    assimilate_this_obs_kind()     =&gt;     assimilate_this_type_of_obs(type_index)
    evaluate_this_obs_kind()       =&gt;       evaluate_this_type_of_obs(type_index)
    use_ext_prior_this_obs_kind()  =&gt;  use_ext_prior_this_type_of_obs(type_index)
    
    get_num_obs_kinds()            =&gt;  get_num_types_of_obs()
    get_num_raw_obs_kinds()        =&gt;  get_num_quantities()
    
    get_obs_kind_index()           =&gt; get_index_for_type_of_obs(type_name)
    get_obs_kind_name()            =&gt; get_name_for_type_of_obs(type_index)
    
    get_raw_obs_kind_index()       =&gt;  get_index_for_quantity(qty_name)
    get_raw_obs_kind_name()        =&gt;  get_name_for_quantity(qty_index)
    
    get_obs_kind_var_type()        =&gt;  get_quantity_for_type_of_obs(type_index)
    
    get_obs_kind()                 =&gt;  get_obs_def_type_of_obs(obs_def)
    set_obs_def_kind()             =&gt;  set_obs_def_type_of_obs(obs_def)
    
    get_kind_from_menu()           =&gt;  get_type_of_obs_from_menu()
    
    read_obs_kind()                =&gt;   read_type_of_obs_table(file_unit, file_format)
    write_obs_kind()               =&gt;  write_type_of_obs_table(file_unit, file_format)
    
    maps obs_seq nums to specific type nums, only used in read_obs_seq:
    map_def_index()                =&gt; map_type_of_obs_table()
    
    removed this.  apparently unused, and simply calls get_obs_kind_name():
    get_obs_name()
    
    apparently unused anywhere, removed:
    add_wind_names()
    do_obs_form_pair()
</pre>
Public integer parameter constants and subroutine formal argument
names, old on left, new on right:
<pre>

   KIND_ =&gt; QTY_
   kind  =&gt; quantity
   
   TYPE_ =&gt; TYPE_
   type  =&gt; type_of_obs
   
   integer parameters:
   max_obs_generic  =&gt;  max_defined_quantities  (not currently public, stays private)
   max_obs_kinds    =&gt;  max_defined_types_of_obs 
</pre></li>
<li>For smaller models we support single file input and output.
These files contain all of the member information, mean, standard
deviation and inflation values for all of the state variables. This
can be run with cycling and all time steps will be appended to the
file.
<p>For <em class="program">perfect_model_obs</em> we provide a
<em class="file">perfect_input.cdl</em> file which contains a
single ensemble member which will be considered the 'truth' and
observations will be generated based on those values. The output
will contain all of the cycling timesteps all of the state
variables.</p>
<p>For <em class="program">filter</em> we provide a <em class=
"file">filter_input.cdl</em> file which contains all of the state
member variables and potentially inflation mean and standard
deviation values. The output will contain all of the cycling
timesteps all of the state variables. Additionally you have the
option to write out different stages during the assimilation in the
&amp;filter_nml <code>stages_to_write</code> mentioned above.</p>
<p>To generate a NetCDF file from a .cdl file run:</p>
<pre>
   ncgen -o perfect_input.nc perfect_input.cdl
   ncgen -o filter_input.nc filter_input.cdl
   </pre></li>
</ol>
<!--==================================================================-->
<a name="NewFeatures" id="NewFeatures"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>New Features</h2>
<ul>
<li>DART now reads and writes NetCDF files for the model state
information. If your model uses NetCDF file format, you no longer
need model_to_dart or dart_to_model to translate to a DART format
file. If your model does not use NetCDF, you can adapt your
model_to_dart and dart_to_model executables to read and write a
NetCDF file for DART to use. The read/write code is part of the core DART
routines so no code is needed in the model_mod model-specific module. There is
a new routine <a href="state_structure.html">add_domain()</a> that a
model_mod::static_init_model() can user to define which NetCDF
variables should be part of the model state, and what DART quantity
(formerly kind) they correspond to.</li>
<li>DART no longer limits the size of a model state to the size of
a single MPI task's memory. The state is read in variable by
variable and distributed across all MPI tasks, so the memory use is
much smaller than previous versions of DART. One-sided MPI
communication is used during the computation of forward operator
values to get required parts of the state from other tasks.</li>
<li>Many of the DART namelists have been simplified, and some items
have moved to a more specific namelist.</li>
<li>Observation sequence files can include externally computed
forward operator values which can be used in the assimilation
instead of calling a forward operator inside DART.</li>
<li>The DART directory structure has been reorganized to make it
easier to identify the various software tools, modules,
documentation and tutorials supplied with the system.</li>
<li>The MATLAB® diagnostic routines have been updated to not
require the MEXNC toolbox. These routines use the built-in NetCDF
support that comes with MATLAB®.</li>
<li>There is a new Particle Filter type. Please contact us if you
are interested in using it.</li>
<li>DART can now take subsets of observation types and restrict
them from impacting certain quantities in the state during the
assimilation. A tool to simplify constructing the table of
interactions is provided (obs_impact_tool).</li>
<li>State Structure
<ul>
<li>Contains information about dimensions and size of variables in
your state. There is a number of accessor functions to get variable
information such as <code>get_variable_size()</code>. See the
<a href="state_structure.html">state_structure.html</a> for more
details.</li>
</ul>
</li>
<li>The POP model_mod now can interpolate Sea Surface Anomaly
observations.</li>
</ul>
<!--==================================================================-->
<a name="SupportedModels" id="SupportedModels"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Supported Models</h2>
Currently we support the models listed below. There are several new
models that have been added that are not on the Lanai Release
including CM1, CICE, and ROMS. Any previously supported models not
on this list are still supported in DART <a href=
"http://www.image.ucar.edu/DAReS/DART/classic/index.html">classic</a>
<ul>
<li><b>9var</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/9var/model_mod.html">9var</a> model.</li>
</ul>
</li>
<li><b>bgrid_solo</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/bgrid_solo/model_mod.html">bgrid solo</a> model.</li>
</ul>
</li>
<li><b>cam-fv</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/cam-fv/model_mod.html">CAM finite volume</a> global
atmospheric model.</li>
<li>Documentation for the <a href=
"http://www.cesm.ucar.edu/models/atm-cam/">CAM model</a>.</li>
</ul>
</li>
<li><b>cice (NEW)</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/cice/model_mod.html">CICE</a> model.</li>
<li>Documentation for the <a href=
"http://www.cesm.ucar.edu/models/ccsm4.0/cice/">CICE
model</a>.</li>
</ul>
</li>
<li><b>cm1 (NEW)</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/cm1/model_mod.html">CM1 cloud-resolving
model</a>.</li>
<li>Documentation for the <a href=
"http://www2.mmm.ucar.edu/people/bryan/cm1/">CM1 model</a>.</li>
</ul>
</li>
<li><b>forced_lorenz_96</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/forced_lorenz_96/model_mod.html">forced lorenz_96</a>
model.</li>
</ul>
</li>
<li><b>lorenz_63</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/lorenz_63/model_mod.html">lorenz_96</a> model.</li>
</ul>
</li>
<li><b>lorenz_84</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/lorenz_84/model_mod.html">lorenz_84</a> model.</li>
</ul>
</li>
<li><b>lorenz_96</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/lorenz_96/model_mod.html">lorenz_96</a> model.</li>
</ul>
</li>
<li><b>lorenz_04</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/lorenz_04/model_mod.html">lorenz_04</a> model.</li>
</ul>
</li>
<li><b>mpas_atm</b> (NetCDF overwrite not supported for
update_u_from_reconstruct = .true. )
<ul>
<li>DART interface documentation for the <a href=
"../../models/mpas_atm/model_mod.html">MPAS atmosphere</a>
model.</li>
<li>Documentation for the <a href=
"https://mpas-dev.github.io/atmosphere/atmosphere.html">MPAS
model</a>.</li>
</ul>
</li>
<li><b>POP</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/POP/model_mod.html">POP</a> global ocean model.</li>
<li>Documentation for the <a href=
"http://www.cesm.ucar.edu/models/ccsm2.0/pop/">POP model</a>.</li>
</ul>
</li>
<li><b>ROMS (NEW)</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/ROMS/model_mod.html">ROMS</a> regional ocean
model.</li>
<li>Documentation for the <a href="https://www.myroms.org/">ROMS
model</a>.</li>
</ul>
</li>
<li><b>simple_advection</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/simple_advection/model_mod.html">simple advection</a>
model.</li>
</ul>
</li>
<li><b>wrf</b>
<ul>
<li>DART interface documentation for the <a href=
"../../models/wrf/model_mod.html">WRF</a> regional forecast
model.</li>
<li>Documentation for the <a href=
"http://www.wrf-model.org/index.php">WRF model</a>.</li>
</ul>
</li>
</ul>

<p>The <em class="file">DART/models/template</em> directory
contains sample files for adding a new model. See the <a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Documentation.php#adding_a_model">
Adding a Model</a> section of the DART web pages for more help on
adding a new model.</p>
<!--==================================================================-->
<a name="ChangedModels" id="ChangedModels"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Changed Models</h2>
<ul>
<li>WRF
<ul>
<li>Allow advanced microphysics schemes (needed interpolation for 7
new kinds)</li>
<li>Interpolation in the vertical is now done in log(p) instead of
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
<li>CAM
<ul>
<li>DART/CAM now runs under the CESM framework, so all options
available with the framework can be used.</li>
<li>Support for the SE core (HOMME) has been developed but is NOT
part of this release. Please contact the <a href=
"mailto:dart@ucar.edu">DART Development Group</a> if you have an
interest in this configuration of CAM.</li>
</ul>
</li>
<li>Simple Advection Model
<ul>
<li>Fixed a bug where the random number generator was being used
before being called with an initial seed.</li>
</ul>
</li>
</ul>
<!--==================================================================-->
<a name="NewFOs" id="NewFOs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>New Observation Types/Forward Operators</h2>
<ul>
<li>Many new observation types related to land and atmospheric
chemistry have been added. See the <a href=
"../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90">
obs_kind_mod.f90</a> for a list of the generic quantities now
available.</li>
<li>New forward operator for Sea Ice (cice) ice thickness
observations. See the <a href=
"../../observations/forward_operators/obs_def_cice_mod.f90">obs_def_cice_mod.f90</a>
file for details.</li>
<li>New forward operator for Carbon Monoxide (CO) Nadir
observations. See the <a href=
"../../observations/forward_operators/obs_def_CO_Nadir_mod.f90">obs_def_CO_Nadir_mod.f90</a>
file for details.</li>
<li>New forward operator for Total Cloud Water in a column
observations. See the <a href=
"../../observations/forward_operators/obs_def_cwp_mod.f90">obs_def_cwp_mod.f90</a>
file for details.</li>
</ul>
<!--==================================================================-->
<a name="NewObs" id="NewObs"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>New Observation Types/Sources</h2>
<ul>
<li>AVISO<br>
Added an observation converter for Sea Surface Height Anomaly
observations. Documentation in <a href=
"../../observations/obs_converters/AVISO/convert_aviso.f90">convert_aviso.f90</a>
(source).</li>
<li>cice<br>
Added an obs_sequence converter for Sea Ice observations.
Documentation in <a href=
"../../observations/obs_converters/cice/cice_to_obs.html">cice_to_obs.html</a>.</li>
<li>GPSPW<br>
Added an obs_sequence converter for GPS precipitable water
observations. Documentation in <a href=
"../../observations/obs_converters/GPSPW/convert_gpspw.f90">convert_gpspw.f90</a>
(source).</li>
<li>MODIS<br>
Added an obs_sequence converter for MODIS FPAR (Fraction of
Photosynthetically Active Radiation) and LAI (Leaf Area Index)
obseverations. Documentation in <a href=
"../../observations/obs_converters/MODIS/MOD15A2_to_obs.html">MOD15A2_to_obs.html</a>.</li>
<li>ok_mesonet<br>
Added an obs_sequence converter for the Oklahoma Mesonet
observations. Documentation in <a href=
"../../observations/obs_converters/ok_mesonet/ok_mesonet.html">ok_mesonet.html</a>.</li>
<li>ROMS<br>
Added an obs_sequence converter for ROMS ocean data. This converter
includes externally computed forward operators output from the ROMS
model using FGAT (First Guess At Time) during the model run.
Documentation in <a href=
"../../observations/obs_converters/ROMS/convert_roms_obs.f90">convert_roms_obs.f90</a>
(source).</li>
<li>SSUSI<br>
Added an obs_sequence converter for wind profiler observations.
Documentation in <a href=
"../../observations/obs_converters/SSUSI/convert_f16_edr_dsk.html">convert_f16_edr_dsk.html</a>.</li>
<li>tropical_cyclone<br>
Added an obs_sequence converter for ASCII format tropical cyclone
track observations. Documentation in <a href=
"../../observations/obs_converters/tropical_cyclone/tc_to_obs.html">
tc_to_obs.html</a>.</li>
</ul>
<!--==================================================================-->
<a name="NewDiagnostics" id="NewDiagnostics"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>New Diagnostics and Documentation</h2>
<p><strong>Better Web Pages.</strong> We've put a lot of effort
into expanding our documentation. For example, please check out
<a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Diagnostics.php#mat_obs">
the MATLAB diagnostics section</a> or the pages outlining the
<a href=
"http://www.image.ucar.edu/DAReS/DART/DART2_Observations.php#obs_seq_overview">
observation sequence file contents</a>.<br></p>
<ul>
<li>The MATLAB® diagnostic routines have been updated to remove the
dependency on third-party toolboxes. These routines use the
built-in netCDF support that comes with basic MATLAB® (no other
toolboxes needed).</li>
</ul>
<p>But there's always more to add. <strong>Please let us know where we are
lacking.</strong></p>
<!--==================================================================-->
<a name="NewUtilities" id="NewUtilities"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>New Utilities</h2>
<p>This section describes updates and changes to the tutorial
materials, scripting, setup, and build information since the Lanai
release.</p>
<ul>
<li><em class="program">obs_impact_tool</em> please refer to
<a href=
"https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/obs_impact_tool/obs_impact_tool.html">
Website</a> or <a href=
"../../assimilation_code/programs/obs_impact_tool/obs_impact_tool.html">
local file</a></li>
<li><em class="program">gen_sampling_error_table</em> now computes
sampling error correction tables for any ensemble size. <!--
       <a href="https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/gen_sampling_err_table/gen_sampling_err_table.html">Website</a>  
       or <a href="../../assimilation_code/programs/gen_sampling_err_table/gen_sampling_err_table.html">local file</a></li>
   --></li>
<li><em class="program">compute_error</em> <a href=
"https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/compute_error/compute_error.html">
Website</a> or <a href=
"../../assimilation_code/programs/compute_error/compute_error.html">
local file</a></li>
</ul>
<!--==================================================================-->
<a name="KnownProblems" id="KnownProblems"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Known Problems</h2>

<p>There are many changes in this release and more updates are
expected to come soon. We are not aware of any obvious bugs, but if
you encounter any unexpected behavior please contact us. Please
watch the dart-users email list for announcements of updates to the
release code, and be prepared to do an 'svn update' from time to
time to get updated files.</p>
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
