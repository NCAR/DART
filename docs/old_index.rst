DART Documentation Main Index
=============================

Overview
--------

The Data Assimilation Research Testbed (DART) is a public-domain, community facility for doing research on and applying
ensemble data assimilation techniques over a variety of models and observation types. It includes many example models
and support for common observation types and for different filter types. It also includes material for teaching and
learning the basic principles of data assimilation.

DART strives to implement general solutions which work over a range of models and observation types, and to develop
theory-based algorithms that solve the many real problems encountered when doing ensemble data assimilation. The
algorithms in DART are tested on both simple one-dimensional models (e.g. the Lorenz series of models) as well as
full-up 3D NWP (Numerical Weather Prediction) models and GCMs (Global Climate Models). The basic Kalman filter code can
be written in a few lines. In practice, however, there are a variety of difficulties resulting from sampling error,
model bias, observation error, lack of model divergence, variations in observation density in space and time, etc. There
are tools built into the DART framework to address many of these problems.

This release of DART includes many new features. DART will now read directly from NetCDF files. If your model uses
NetCDF file format the model_to_dart and dart_to_model steps are no longer needed. Given that many jobs spend a large
percentage of time doing file I/O, this can be a significant speedup in the overall assimilation cycle. DART now
distributes the ensemble of model states across the MPI tasks, removing the hard memory limit that a single ensemble
member's data fit into the memory of a single task. This removes the memory limit for models at high resolution or with
nested grids. DART can assimilate forward operators computed outside of the filter process. Users can provide a table of
observations and state vector items along with a factor that can easily prevent one class of observations from impacting
another class of state vector data, and vice versa. The DART directory structure has been revamped to help users more
easily find the various utilities, tools, supported models, and observation types that come with the DART distribution.
The Matlab diagnostic routines have been rewritten to no longer require the external MEXNC toolbox. They now use the
intrinsic Matlab NetCDF functions.

To get started, look here:

-  Our extensive `web pages <http://www.image.ucar.edu/DAReS/>`__
-  The :doc:`./html/Manhattan_release` which include installation hints, a walk-through of building and running a model,
   and an overview of the diagnostics.
-  Documentation for the main assimilation program :doc:`../assimilation_code/programs/filter/filter`
-  Documentation for the other related programs that come with DART
-  Documentation for :doc:`../observations/obs_converters/observations` supported by DART
-  Documentation for the Matlab diagnostic tools
-  Discussion of :doc:`../assimilation_code/modules/assimilation/assim_tools_mod`,
   `inflation <../assimilation_code/programs/filter/filter.html#Inflation>`__
-  The web pages specific to each model
-  The web pages specific to each observation type

The best way to get to know the DART software is to follow along while reading the tutorial documents in the
:doc:`./DART_LAB/DART_LAB` and then the :doc:`./tutorial/index` directory.

The latest official release is named "Manhattan". See the extensive release notes :doc:`./html/Manhattan_release` which
include installation help, a walk-through of building and running a model, and then examples of how to use the
diagnostics to evaluate the success of the assimilation. See :doc:`./html/Manhattan_diffs_from_Lanai` for a brief
summary of changes since Lanai, including new functionality, new models and tools, and any non-backwards-compatible
changes since the Lanai release.

Future releases of DART are expected to have substantial changes and will be less backwards compatible than has been
historically true with DART releases. New development will continue on a separate subversion branch.

Every source file in the DART system has a corresponding .html file that contains specifics for public interfaces in
each of the DART modules, for the executable programs which come with DART, and for how to interface new models into the
DART system.

The remainder of this page contains links to all the documentation for this DART release.

Software updates
----------------

Updates to the release are now summarized in a file at the top level called `CHANGELOG <../CHANGELOG>`__. The latest
updates are at the end of this file.

User documentation and tutorials
--------------------------------

Start `here <http://www.image.ucar.edu/DAReS/>`__ if you are looking for DART User-level HTML documentation. DART comes
with an extensive set of documentation including release notes for each version, a walk-through on-line tutorial, and a
full set of pdf and framemaker tutorial materials.

The Manhattan :doc:`./html/Manhattan_release` include installation hints, a walk-through of building and running a
model, and an overview of the diagnostics. It also includes a list of new or changed models, observation support,
diagnostics, and non-backwards compatible changes.

For a shorter summary document of only the changes since the last release see :doc:`./html/Manhattan_diffs_from_Lanai`.
This may be more helpful for current DART users who are looking for pointers to differences when they update.

Three tutorials, in PDF format, are available. The first is more introductory and interactive with
:doc:`./DART_LAB/DART_LAB`. The :doc:`./tutorial/index` is more of a workshop format, with multiple sections covering
various parts of DART with suggested exercises at the end of most sections.

All sections below this one are detailed information on the programming interfaces of all the DART modules, the namelist
details, the executable programs which are part of DART. For introductory materials, see the links above.

Programs
--------

DART contains many library functions and separate executable programs. The main DART executable is the ``filter``
program. Other programs generate data or format the diagnostic information.

The executable programs that come with DART include:

-  :doc:`../assimilation_code/programs/filter/filter` - the main assimilation code
-  :doc:`../assimilation_code/programs/perfect_model_obs/perfect_model_obs` - run a model to generate synthetic
   observation values
-  :doc:`../assimilation_code/programs/create_obs_sequence/create_obs_sequence` - interactive program to generate
   observations
-  :doc:`../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq` - repeat a set of observations
   at multiple times
-  :doc:`../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool` - general observation sequence file
   manipulation tool
-  :doc:`../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart` - [deprecated] - initialize
   inflation files
-  :doc:`../assimilation_code/programs/advance_time/advance_time` - increment calendar times, useful for scripting loops
-  :doc:`../assimilation_code/programs/closest_member_tool/closest_member_tool` - select DART restart file closest to
   mean
-  :doc:`../assimilation_code/programs/integrate_model/integrate_model` - wrapper for models called as subroutines
-  :doc:`../assimilation_code/programs/preprocess/preprocess` - used during compiling
-  :doc:`../build_templates/mkmf` - used to generate makefiles during compiling
-  :doc:`../assimilation_code/programs/wakeup_filter/wakeup_filter` - used when filter runs a parallel model advance
-  :doc:`../assimilation_code/programs/system_simulation/system_simulation` (sampling error correction) - generate the
   files used for Sampling Error Correction option

The diagnostic programs that process observations after being assimilated by DART include:

-  :doc:`../assimilation_code/programs/obs_diag/oned/obs_diag` - low order model diagnostics
-  :doc:`../assimilation_code/programs/obs_diag/threed_sphere/obs_diag` - full 3d model diagnostics
-  :doc:`../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf` - convert output obs sequence files into
   netcdf format
-  :doc:`../assimilation_code/programs/obs_common_subset/obs_common_subset` - select a common subset of obs from
   multiple files
-  :doc:`../assimilation_code/programs/obs_selection/obs_selection` - select a given set of obs from a longer sequence
-  :doc:`../assimilation_code/programs/obs_seq_coverage/obs_seq_coverage` - select a consistent set of obs through time
-  :doc:`../assimilation_code/programs/obs_seq_verify/obs_seq_verify` - convert obs to a netcdf file formatted for
   forecast verification
-  :doc:`../assimilation_code/programs/compare_states/compare_states` - compare fields within multiple restart files
-  :doc:`../assimilation_code/programs/model_mod_check/model_mod_check` - development and testing tool during interface
   development
-  :doc:`../developer_tests/utilities/PrecisionCheck` - compiler/platform check of Fortran real/integer precision

The executable programs that convert observations into DART format include:

-  :doc:`../observations/obs_converters/observations`
-   
-  :doc:`../observations/obs_converters/AIRS/AIRS`
-  `AURA temperature data <../observations/obs_converters/AURA/convert_aura.f90>`__ (source)
-  :doc:`../observations/obs_converters/AVISO/AVISO`
-  :doc:`../observations/obs_converters/Ameriflux/level4_to_obs`
-  :doc:`../observations/obs_converters/cice/cice_to_obs`
-  `CHAMP data <../observations/obs_converters/CHAMP/CHAMP_density_text_to_obs.f90>`__ (source)
-  `CNOFS data <../observations/obs_converters/CNOFS/CNOFS_text_to_obs.f90>`__ (source)
-  :doc:`../observations/obs_converters/COSMOS/COSMOS_development` (development format)
-  :doc:`../observations/obs_converters/COSMOS/COSMOS_to_obs`
-  :doc:`../observations/obs_converters/DWL/dwl_to_obs`
-  `Even Sphere data <../observations/obs_converters/even_sphere/even_sphere.m>`__ (source)
-  `GITM data <../observations/obs_converters/text_GITM/text_to_obs.f90>`__ (source)
-  :doc:`../observations/obs_converters/gps/gps`
-  `Ground GPS Vtec data <../observations/obs_converters/gnd_gps_vtec/gnd_gps_vtec_text_to_obs.f90>`__ (source)
-  `GSI data <../observations/obs_converters/GSI2DART/gsi_to_dart.f90>`__ (source)
-  :doc:`../observations/obs_converters/GTSPP/GTSPP`
-  :doc:`../observations/obs_converters/MADIS/MADIS`
-  :doc:`../observations/obs_converters/MIDAS/MIDAS_to_obs` (netcdf intermediate files)
-  :doc:`../observations/obs_converters/snow/snow_to_obs` (source)
-  :doc:`../observations/obs_converters/MODIS/MODIS_README` (ORNL DAAC)
-  :doc:`../observations/obs_converters/NCEP/prep_bufr/prep_bufr`
-  :doc:`../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs`
-  :doc:`../observations/obs_converters/ok_mesonet/ok_mesonet`
-  :doc:`../observations/obs_converters/quikscat/QuikSCAT`
-  `ROMS data <../observations/obs_converters/ROMS/convert_roms_obs.f90>`__ (source)
-  :doc:`../observations/obs_converters/radar/radar`
-  `SABER data <../observations/obs_converters/SABER/convert_saber_cdf.f90>`__ (source)
-  :doc:`../observations/obs_converters/SSEC/SSEC`
-  :doc:`../observations/obs_converters/SSUSI/convert_f16_edr_dsk`
-  :doc:`../observations/obs_converters/text/text_to_obs`
-  :doc:`../observations/obs_converters/tropical_cyclone/tc_to_obs` (source)
-  :doc:`../observations/obs_converters/tpw/tpw` (source)
-  :doc:`../observations/obs_converters/var/littler_tf_dart`
-  :doc:`../observations/obs_converters/var/rad_3dvar_to_dart`
-  :doc:`../observations/obs_converters/var/var`
-  :doc:`../observations/obs_converters/WOD/WOD`

Models
------

DART comes with several models which can be used to learn about data assimilation, to do actual experiments with real
observations, or to use as a template for adding additional models to DART.

All models in the DART project have individual documentation pages, which can be found here (if an html document is not
available, the link is to the .f90 source):

Currently Manhattan has support for many of our larger models such as WRF, POP, CAM, CICE, CLM, ROMS, MPAS_ATM, ... and
all lower models such as lorenz_96. Models previously available on Lanai can still be used with DART
`classic <https://svn-dares-dart.cgd.ucar.edu/DART/releases/classic/>`__.

**Supported in Manhattan**

-  :doc:`../models/9var/model_mod` (html)
-  :doc:`../models/bgrid_solo/model_mod` (html)
-  :doc:`../models/cam-fv/model_mod` (html)
-  `cice <../models/cice/model_mod.f90>`__ (source)
-  :doc:`../models/clm/model_mod` (html)
-  :doc:`../models/forced_lorenz_96/model_mod` (html)
-  :doc:`../models/lorenz_04/model_mod` (html)
-  :doc:`../models/lorenz_63/model_mod` (html)
-  :doc:`../models/lorenz_84/model_mod` (html)
-  :doc:`../models/lorenz_96/model_mod` (html)
-  `Lorenz_96_2scale <../models/lorenz_96_2scale/model_mod.f90>`__ (source)
-  :doc:`../models/mpas_atm/model_mod` (html)
-  :doc:`../models/null_model/model_mod` (html)
-  :doc:`../models/POP/model_mod` (html)
-  :doc:`../models/ROMS/model_mod` (html)
-  :doc:`../models/simple_advection/model_mod` (html)
-  :doc:`../models/template/model_mod`
-  :doc:`../models/wrf/model_mod` (html)

**Supported in Classic**

-  `AM2 <../models/am2/model_mod.f90>`__ (source)
-  :doc:`../models/coamps/model_mod` (html)
-  `COAMPS_nest <../models/coamps_nest/model_mod.f90>`__ (source)
-  `dynamo <../models/dynamo/model_mod.f90>`__ (source)
-  `forced_barot <../models/forced_barot/model_mod.f90>`__ (source)
-  :doc:`../models/gitm/model_mod` (html)
-  :doc:`../models/ikeda/model_mod` (html)
-  `MITgcm_annulus <../models/MITgcm_annulus/model_mod.f90>`__ (source)
-  :doc:`../models/MITgcm_ocean/model_mod` (html)
-  :doc:`../models/mpas_ocn/model_mod` (html)
-  `NAAPS <../models/NAAPS/model_mod.f90>`__ (source)
-  :doc:`../models/NCOMMAS/model_mod` (html)
-  :doc:`../models/noah/model_mod` (html)
-  :doc:`../models/pe2lyr/model_mod` (html)
-  `Rose <../models/rose/model_mod.f90>`__ (source)
-  :doc:`../models/sqg/model_mod` (html)
-  `TIEgcm <../models/tiegcm/model_mod.f90>`__ (source)

Namelists
---------

Generally read from the file ``input.nml``. We adhere to the F90 standard of starting a namelist with an ampersand '&'
and terminating with a slash '/'.

Namelists for Programs:

-  `&closest_member_tool_nml <../assimilation_code/programs/closest_member_tool/closest_member_tool.html#Namelist>`__
-  `&compare_states_nml <../assimilation_code/programs/compare_states/compare_states.html#Namelist>`__
-  `&filter_nml <../assimilation_code/programs/filter/filter.html#Namelist>`__
-  `&full_error_nml <../assimilation_code/programs/system_simulation/system_simulation.html#Namelist>`__ (system
   simulation)
-  `&model_mod_check_nml <../assimilation_code/programs/model_mod_check/model_mod_check.html#Namelist>`__
-  `&obs_common_subset_nml <../assimilation_code/programs/obs_common_subset/obs_common_subset.html#Namelist>`__
-  `&obs_diag_nml <../assimilation_code/programs/obs_diag/oned/obs_diag.html#Namelist>`__ (oned)
-  `&obs_diag_nml <../assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html#Namelist>`__ (threed_sphere)
-  `&obs_loop_nml <../assimilation_code/programs/obs_loop/obs_loop.nml>`__
-  `&obs_selection_nml <../assimilation_code/programs/obs_selection/obs_selection.html#Namelist>`__
-  `&obs_seq_coverage_nml <../assimilation_code/programs/obs_seq_coverage/obs_seq_coverage.html#Namelist>`__
-  `&obs_seq_to_netcdf_nml <../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html#Namelist>`__
-  `&obs_seq_verify_nml <../assimilation_code/programs/obs_seq_verify/obs_seq_verify.html#Namelist>`__
-  `&obs_sequence_tool_nml <../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html#Namelist>`__
-  `&perfect_model_obs_nml <../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html#Namelist>`__
-  `&preprocess_nml <../assimilation_code/programs/preprocess/preprocess.html#Namelist>`__

Namelists for Observation Conversion Programs:

-  `&convert_airs_L2_nml <../observations/obs_converters/AIRS/AIRS.html#Namelist>`__
-  `&convert_L2b_nml <../observations/obs_converters/quikscat/QuikSCAT.html#Namelist>`__
-  `&convert_tpw_nml <../observations/obs_converters/tpw/tpw.html#Namelist>`__
-  `&COSMOS_development_nml <../observations/obs_converters/COSMOS/COSMOS_development.html#Namelist>`__
-  `&COSMOS_to_obs_nml <../observations/obs_converters/COSMOS/COSMOS_to_obs.html#Namelist>`__
-  `&convert_cosmic_gps_nml <../observations/obs_converters/gps/gps.html#Namelist>`__
-  `&level4_to_obs_nml <../observations/obs_converters/Ameriflux/level4_to_obs.html#Namelist>`__
-  `&MIDAS_to_obs_nml <../observations/obs_converters/MIDAS/MIDAS_to_obs.html#Namelist>`__
-  `&MOD15A2_to_obs_nml <../observations/obs_converters/MODIS/MOD15A2_to_obs.html#Namelist>`__
-  `&ncepobs_nml <../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs.html#Namelist>`__
-  `&tc_to_obs_nml <../observations/obs_converters/tropical_cyclone/tc_to_obs.html#Namelist>`__
-  `&rad_3dvar_to_dart_nml <../observations/obs_converters/var/rad_3dvar_to_dart.html#Namelist>`__
-  `&wod_to_obs_nml <../observations/obs_converters/WOD/WOD.html#Namelist>`__

Namelists for Modules:

-  `&assim_model_mod_nml <../assimilation_code/modules/assimilation/assim_model_mod.html#Namelist>`__
-  `&assim_tools_mod_nml <../assimilation_code/modules/assimilation/assim_tools_mod.html#Namelist>`__
-  `&cov_cutoff_mod_nml <../assimilation_code/modules/assimilation/cov_cutoff_mod.html#Namelist>`__
-  `&ensemble_manager_mod_nml <../assimilation_code/modules/utilities/ensemble_manager_mod.html#Namelist>`__
-  `&location_mod_nml <../assimilation_code/location/channel/location_mod.html#Namelist>`__ (channel)
-  `&location_mod_nml <../assimilation_code/location/column/location_mod.nml>`__ (column)
-  `&location_mod_nml <../assimilation_code/location/threed_cartesian/location_mod.html#Namelist>`__ (threed_cartesian)
-  `&location_mod_nml <../assimilation_code/location/threed_sphere/location_mod.html#Namelist>`__ (threed_sphere)
-  `&mpi_utilities_mod_nml <../assimilation_code/modules/utilities/mpi_utilities_mod.html#Namelist>`__
-  `&obs_def_gps_mod_nml <../observations/forward_operators/obs_def_gps_mod.html#Namelist>`__
-  `&obs_def_ocean_mod_nml <../observations/forward_operators/obs_def_ocean_mod.nml>`__
-  `&obs_def_radar_mod_nml <../observations/forward_operators/obs_def_radar_mod.html#Namelist>`__
-  `&obs_def_tower_mod_nml <../observations/forward_operators/obs_def_tower_mod.nml>`__
-  `&obs_def_tpw_mod_nml <../observations/forward_operators/obs_def_tpw_mod.nml>`__
-  `&obs_kind_mod_nml <../assimilation_code/modules/observations/obs_kind_mod.html#Namelist>`__
-  `&obs_sequence_mod_nml <../assimilation_code/modules/observations/obs_sequence_mod.html#Namelist>`__
-  `&reg_factor_mod_nml <../assimilation_code/modules/assimilation/reg_factor_mod.html#Namelist>`__
-  `&smoother_mod_nml <../assimilation_code/modules/assimilation/smoother_mod.html#Namelist>`__
-  `&schedule_mod_nml <../assimilation_code/modules/utilities/schedule_mod.html#Namelist>`__
-  `&utilities_mod_nml <../assimilation_code/modules/utilities/utilities_mod.html#Namelist>`__

Namelists for Models:

-  9var `&model_nml <../models/9var/model_mod.html#Namelist>`__
-  bgrid_solo `&model_nml <../models/bgrid_solo/model_mod.html#Namelist>`__
-  cam `&model_nml <../models/cam-fv/model_mod.html#Namelist>`__
-  clm `&model_nml <../models/clm/model_mod.html#Namelist>`__
-  coamps `&model_nml <../models/coamps/model_mod.html#Namelist>`__
-  coamps_nest `&model_nml <../models/coamps_nest/model_mod.html#Namelist>`__
-  forced_lorenz_96 `&model_nml <../models/forced_lorenz_96/model_mod.html#Namelist>`__
-  ikeda `&model_nml <../models/ikeda/model_mod.html#Namelist>`__
-  lorenz_04 `&model_nml <../models/lorenz_04/model_mod.html#Namelist>`__
-  lorenz_63 `&model_nml <../models/lorenz_63/model_mod.html#Namelist>`__
-  lorenz_84 `&model_nml <../models/lorenz_84/model_mod.html#Namelist>`__
-  lorenz_96 `&model_nml <../models/lorenz_96/model_mod.html#Namelist>`__
-  lorenz_96_2scale `&model_nml <../models/lorenz_96_2scale/model_mod.html#Namelist>`__
-  MITgcm_ocean `&create_ocean_obs_nml <../models/MITgcm_ocean/create_ocean_obs.html#Namelist>`__
-  MITgcm_ocean `&model_nml <../models/MITgcm_ocean/model_mod.html#Namelist>`__
-  mpas_atm `&model_nml <../models/mpas_atm/model_mod.html#Namelist>`__
-  mpas_ocn `&model_nml <../models/mpas_ocn/model_mod.html#Namelist>`__
-  NAAPS `&model_nml <../models/NAAPS/model_mod.nml>`__
-  NAAPS `&model_mod_check_nml <../models/NAAPS/model_mod_check.nml>`__
-  NCOMMAS `&model_nml <../models/NCOMMAS/model_mod.html#Namelist>`__
-  NCOMMAS `&ncommas_vars_nml <../models/NCOMMAS/model_mod.html#Namelist>`__
-  noah `&model_nml <../models/noah/model_mod.html#Namelist>`__
-  null_model `&model_nml <../models/null_model/model_mod.html#Namelist>`__
-  POP `&model_nml <../models/POP/model_mod.html#Namelist>`__
-  ROMS `&model_nml <../models/ROMS/model_mod.html#Namelist>`__
-  simple_advection `&model_nml <../models/simple_advection/model_mod.html#Namelist>`__
-  sqg `&model_nml <../models/sqg/model_mod.html#Namelist>`__
-  template `&model_nml <../models/template/model_mod.html#Namelist>`__
-  wrf `&model_nml <../models/wrf/model_mod.html#Namelist>`__
-  wrf `&replace_wrf_fields_nml <../models/wrf/WRF_DART_utilities/replace_wrf_fields.html#Namelist>`__
-  wrf `&wrf_dart_obs_preprocess_nml <../models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess.html#Namelist>`__

Modules
-------

All modules in the DART project have individual documentation pages, which can be found here:

Assimilation Modules

-  :doc:`../assimilation_code/modules/assimilation/adaptive_inflate_mod`
-  :doc:`../assimilation_code/modules/assimilation/assim_tools_mod`
-  :doc:`../assimilation_code/modules/assimilation/assim_model_mod`
-  :doc:`../assimilation_code/modules/assimilation/assim_tools_mod`
-  :doc:`../assimilation_code/modules/assimilation/cov_cutoff_mod`
-  :doc:`../assimilation_code/modules/assimilation/filter_mod`
-  :doc:`../assimilation_code/modules/assimilation/obs_model_mod`
-  `assimilation_code/modules/assimilation/quality_control.f90 <../assimilation_code/modules/assimilation/quality_control_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/modules/assimilation/reg_factor_mod`
-  `assimilation_code/modules/assimilation/sampling_error_correction_mod.f90 <../assimilation_code/modules/assimilation/sampling_error_correction_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/modules/assimilation/smoother_mod`

Location Modules

-  `assimilation_code/location/annulus/location_mod.f90 <../assimilation_code/location/annulus/location_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/location/channel/location_mod`
-  `assimilation_code/location/column/location_mod.f90 <../assimilation_code/location/column/location_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/location/oned/location_mod`
-  `assimilation_code/location/threed/location_mod.f90 <../assimilation_code/location/threed/location_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/location/threed_cartesian/location_mod`
-  `assimilation_code/location/threed_cartesian/xyz_location_mod.f90 <../assimilation_code/location/threed_cartesian/xyz_location_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/location/threed_sphere/location_mod`
-  `assimilation_code/location/twod/location_mod.f90 <../assimilation_code/location/twod/location_mod.f90>`__ (source)
-  `assimilation_code/location/twod_annulus/location_mod.f90 <../assimilation_code/location/twod_annulus/location_mod.f90>`__
   (source)
-  `assimilation_code/location/twod_sphere/location_mod.f90 <../assimilation_code/location/twod_sphere/location_mod.f90>`__
   (source)

Observation Modules

-  :doc:`../assimilation_code/modules/observations/DEFAULT_obs_kind_mod`
-  `assimilation_code/modules/observations/forward_operator_mod.f90 <../assimilation_code/modules/observations/forward_operator_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/modules/observations/obs_kind_mod`
-  :doc:`../assimilation_code/modules/observations/obs_sequence_mod`

I/O Modules

-  `assimilation_code/modules/io/dart_time_io_mod.f90 <../assimilation_code/modules/io/dart_time_io_mod.f90>`__ (source)
-  `assimilation_code/modules/io/direct_netcdf_mod.f90 <../assimilation_code/modules/io/direct_netcdf_mod.f90>`__
   (source)
-  `assimilation_code/modules/io/io_filenames_mod.f90 <../assimilation_code/modules/io/io_filenames_mod.f90>`__ (source)
-  `assimilation_code/modules/io/single_file_io_mod.f90 <../assimilation_code/modules/io/single_file_io_mod.f90>`__
   (source)
-  `assimilation_code/modules/io/state_structure_mod.f90 <../assimilation_code/modules/io/state_structure_mod.f90>`__
   (source)
-  `assimilation_code/modules/io/state_vector_io_mod.f90 <../assimilation_code/modules/io/state_vector_io_mod.f90>`__
   (source)

Utilities Modules

-  `assimilation_code/modules/utilities/assert_mod.f90 <../assimilation_code/modules/utilities/assert_mod.f90>`__
   (source)
-  `assimilation_code/modules/utilities/cray_win_mod.f90 <../assimilation_code/modules/utilities/cray_win_mod.f90>`__
   (source)
-  `assimilation_code/modules/utilities/distributed_state_mod.f90 <../assimilation_code/modules/utilities/distributed_state_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/modules/utilities/ensemble_manager_mod`
-  `assimilation_code/modules/utilities/obs_impact_mod.f90 <../assimilation_code/modules/utilities/obs_impact_mod.f90>`__
   (source)
-  `assimilation_code/modules/utilities/parse_args_mod.f90 <../assimilation_code/modules/utilities/parse_args_mod.f90>`__
   (source)
-  :doc:`../assimilation_code/modules/utilities/mpi_utilities_mod`
-  :doc:`../assimilation_code/modules/utilities/random_seq_mod`
-  :doc:`../assimilation_code/modules/utilities/schedule_mod`
-  `assimilation_code/modules/utilities/sort_mod.f90 <../assimilation_code/modules/utilities/sort_mod.f90>`__ (source)
-  :doc:`../assimilation_code/modules/utilities/time_manager_mod`
-  :doc:`../assimilation_code/modules/utilities/types_mod`
-  :doc:`../assimilation_code/modules/utilities/utilities_mod`

Example Model Module

-  :doc:`../models/POP/dart_pop_mod`

Forward Operators Modules

-  :doc:`../observations/forward_operators/DEFAULT_obs_def_mod`
-  :doc:`../observations/forward_operators/DEFAULT_obs_def_mod`
-  :doc:`../observations/forward_operators/obs_def_1d_state_mod`
-  `observations/forward_operators/obs_def_AIRS_mod.f90 <../observations/forward_operators/obs_def_AIRS_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_altimeter_mod.f90 <../observations/forward_operators/obs_def_altimeter_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_AOD_mod.f90 <../observations/forward_operators/obs_def_AOD_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_AURA_mod.f90 <../observations/forward_operators/obs_def_AURA_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_cice_mod.f90 <../observations/forward_operators/obs_def_cice_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_cloud_mod.f90 <../observations/forward_operators/obs_def_cloud_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_CO_Nadir_mod.f90 <../observations/forward_operators/obs_def_CO_Nadir_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_COSMOS_mod.f90 <../observations/forward_operators/obs_def_COSMOS_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_cwp_mod.f90 <../observations/forward_operators/obs_def_cwp_mod.f90>`__
   (source)
-  :doc:`../observations/forward_operators/obs_def_dew_point_mod`
-  `observations/forward_operators/obs_def_dwl_mod.f90 <../observations/forward_operators/obs_def_dwl_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_eval_mod.f90 <../observations/forward_operators/obs_def_eval_mod.f90>`__
   (source)
-  :doc:`../observations/forward_operators/obs_def_gps_mod`
-  `observations/forward_operators/obs_def_gts_mod.f90 <../observations/forward_operators/obs_def_gts_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_GWD_mod.f90 <../observations/forward_operators/obs_def_GWD_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_metar_mod.f90 <../observations/forward_operators/obs_def_metar_mod.f90>`__
   (source)
-  :doc:`../observations/forward_operators/obs_def_mod`
-  :doc:`../observations/forward_operators/obs_def_ocean_mod`
-  `observations/forward_operators/obs_def_pe2lyr_mod.f90 <../observations/forward_operators/obs_def_pe2lyr_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_QuikSCAT_mod.f90 <../observations/forward_operators/obs_def_QuikSCAT_mod.f90>`__
   (source)
-  :doc:`../observations/forward_operators/obs_def_radar_mod`
-  `observations/forward_operators/obs_def_radiance_mod.f90 <../observations/forward_operators/obs_def_radiance_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_reanalysis_bufr_mod.f90 <../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_rel_humidity_mod.f90 <../observations/forward_operators/obs_def_rel_humidity_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_SABER_mod.f90 <../observations/forward_operators/obs_def_SABER_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_simple_advection_mod.f90 <../observations/forward_operators/obs_def_simple_advection_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_sqg_mod.f90 <../observations/forward_operators/obs_def_sqg_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_surface_mod.f90 <../observations/forward_operators/obs_def_surface_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_TES_nadir_mod.f90 <../observations/forward_operators/obs_def_TES_nadir_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_tower_mod.f90 <../observations/forward_operators/obs_def_tower_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_tpw_mod.f90 <../observations/forward_operators/obs_def_tpw_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_upper_atm_mod.f90 <../observations/forward_operators/obs_def_upper_atm_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_vortex_mod.f90 <../observations/forward_operators/obs_def_vortex_mod.f90>`__
   (source)
-  `observations/forward_operators/obs_def_wind_speed_mod.f90 <../observations/forward_operators/obs_def_wind_speed_mod.f90>`__
   (source)

Directory tree
--------------

NOTE: 'work', 'matlab', and 'shell_scripts' directory names have been removed from this list.

::

     |--assimilation_code
     |  |--location
     |  |  |--annulus
     |  |  |--channel
     |  |  |--column
     |  |  |--oned
     |  |  |--threed
     |  |  |--threed_cartesian
     |  |  |--threed_sphere
     |  |  |--twod
     |  |  |--twod_annulus
     |  |  |--twod_sphere
     |  |--modules
     |  |  |--assimilation
     |  |  |--io
     |  |  |--observations
     |  |  |--utilities
     |  |--programs
     |  |  |--advance_time
     |  |  |--closest_member_tool
     |  |  |--compare_states
     |  |  |  |--work
     |  |  |--compute_error
     |  |  |--create_fixed_network_seq
     |  |  |--create_obs_sequence
     |  |  |--fill_inflation_restart
     |  |  |--filter
     |  |  |--gen_sampling_err_table
     |  |  |  |--work
     |  |  |--integrate_model
     |  |  |--model_mod_check
     |  |  |--obs_common_subset
     |  |  |--obs_diag
     |  |  |  |--oned
     |  |  |  |--threed_cartesian
     |  |  |  |--threed_sphere
     |  |  |--obs_impact_tool
     |  |  |--obs_loop
     |  |  |--obs_selection
     |  |  |--obs_seq_coverage
     |  |  |--obs_seq_to_netcdf
     |  |  |--obs_sequence_tool
     |  |  |--obs_seq_verify
     |  |  |--perfect_model_obs
     |  |  |--preprocess
     |  |  |--system_simulation
     |  |  |  |--final_full_precomputed_tables
     |  |  |  |--work
     |  |  |--wakeup_filter
     |  |--scripts
     |--build_templates
     |--developer_tests
     |  |--forward_operators
     |  |--harnesses
     |  |  |--filename_harness
     |  |  |--read_transpose_write
     |  |--io
     |  |  |--work
     |  |--location
     |  |  |--annulus
     |  |  |  |--test
     |  |  |--channel
     |  |  |  |--test
     |  |  |--column
     |  |  |  |--test
     |  |  |--oned
     |  |  |  |--test
     |  |  |--threed
     |  |  |  |--test
     |  |  |--threed_cartesian
     |  |  |  |--test
     |  |  |--threed_sphere
     |  |  |  |--test
     |  |  |--twod
     |  |  |  |--test
     |  |  |--twod_annulus
     |  |  |  |--test
     |  |  |--twod_sphere
     |  |     |--test
     |  |--mpi_utilities
     |  |  |--tests
     |  |--obs_sequence
     |  |  |--data
     |  |  |--work
     |  |--random_seq
     |  |  |--test
     |  |--reg_factor
     |  |--time_manager
     |  |--utilities
     |     |--work
     |--diagnostics
     |  |--matlab
     |     |--deprecated
     |     |--private
     |--docs
     |  |--DART_LAB
     |  |  |--matlab
     |  |  |  |--private
     |  |  |--presentation
     |  |--Graphs
     |  |--html
     |  |  |--boilerplate
     |  |  |--design
     |  |  |--history
     |  |--images
     |  |--tutorial
     |--observations
        |--forward_operators
        |  |--test
        |--obs_converters
           |--AIRS
           |  |--data
           |  |--output
           |--Ameriflux
           |--AURA
           |  |--data
           |--AVISO
           |--CHAMP
           |--cice
           |  |--data
           |--CNOFS
           |--COSMOS
           |  |--data
           |--DWL
           |  |--data
           |--even_sphere
           |--gnd_gps_vtec
           |--gps
           |  |--cosmic
           |  |  |--20071001
           |  |--matlab
           |--GPSPW
           |  |--data
           |--GTSPP
           |  |--data
           |  |--matlab
           |--MADIS
           |  |--data
           |--MIDAS
           |  |--data
           |--MODIS
           |  |--data
           |--NCEP
           |  |--ascii_to_obs
           |  |--prep_bufr
           |     |--blk_ublk
           |     |--convert_bufr
           |     |--data
           |     |  |--201012
           |     |--docs
           |     |  |--Reason_codes
           |     |--exe
           |     |--lib
           |     |--src
           |--obs_error
           |--ok_mesonet
           |  |--data
           |--quikscat
           |  |--data
           |--radar
           |  |--examples
           |--ROMS
           |  |--data
           |--SABER
           |  |--data
           |  |--progs
           |--snow
           |  |--data
           |--SSEC
           |  |--data
           |--SSUSI
           |  |--data
           |--text
           |  |--data
           |--text_GITM
           |--tpw
           |  |--data
           |  |--doc
           |--tropical_cyclone
           |  |--data
           |--utilities
           |  |--oned
           |  |--threed_sphere
           |--var
           |  |--3DVAR_OBSPROC
           |  |--data
           |--WOD
              |--data
            

::

    
     |--models
        |--9var
        |--am2
        |--bgrid_solo
        |  |--fms_src
        |  |  |--atmos_bgrid
        |  |  |  |--driver
        |  |  |  |  |--solo
        |  |  |  |--model
        |  |  |  |--tools
        |  |  |--atmos_param
        |  |  |  |--hs_forcing
        |  |  |--atmos_shared
        |  |  |  |--tracer_driver
        |  |  |  |--vert_advection
        |  |  |--atmos_solo
        |  |  |--shared
        |  |     |--axis_utils
        |  |     |--constants
        |  |     |--diag_manager
        |  |     |--fft
        |  |     |--field_manager
        |  |     |--fms
        |  |     |--horiz_interp
        |  |     |--mpp
        |  |     |--platform
        |  |     |--sat_vapor_pres
        |  |     |--time_manager
        |  |     |--topography
        |  |     |--tracer_manager
        |  |     |--udunits
        |  |     |--utilities
        |  |--test
        |--cam-fv
        |  |--deprecated
        |  |--doc
        |  |--shell_scripts
        |     |--cesm1_5
        |     |--cesm2_0
        |--cam-old
        |  |--deprecated
        |  |--doc
        |  |--full_experiment
        |  |--perfect_model
        |--CESM
        |  |--doc
        |     |--CESM_DART_assim_modes
        |--cice
        |--clm
        |  |--datm
        |  |--docs
        |--cm1
        |--coamps
        |  |--doc
        |  |--externals
        |  |  |--obs_def
        |  |--templates
        |--coamps_nest
        |  |--doc
        |  |--externals
        |  |  |--obs_def
        |  |--shell_scripts
        |  |  |--COAMPS_RESTART_SCRIPTS
        |  |  |--TEMPLATES
        |  |--templates
        |  |  |--EXPERIMENT_EXAMPLE
        |--dynamo
        |  |--data
        |--ECHAM
        |--forced_barot
        |  |--obs
        |--forced_lorenz_96
        |--gitm
        |  |--GITM2
        |  |  |--src
        |  |--python
        |  |--testdata1
        |--ikeda
        |--LMDZ
        |--lorenz_04
        |--lorenz_63
        |--lorenz_84
        |--lorenz_96
        |  |--tests
        |--lorenz_96_2scale
        |--MITgcm_annulus
        |--MITgcm_ocean
        |  |--inputs
        |--model_mod_tools
        |--mpas_atm
        |  |--data
        |--mpas_ocn
        |  |--data
        |--NAAPS
        |--NCOMMAS
        |  |--docs
        |--noah
        |  |--ensemble_source
        |  |--forcing
        |  |--templates
        |--null_model
        |--PBL_1d
        |--pe2lyr
        |--POP
        |--ROMS
        |  |--data
        |  |--doc
        |--rose
        |--simple_advection
        |--sqg
        |--template
        |--tiegcm
        |--wrf
           |--experiments
           |  |--Radar
           |     |--IC
           |     |  |--sounding_perturbation
           |     |--obs
           |--namelist
           |--PERTURB
           |  |--3DVAR-COVAR
           |--regression
           |  |--CONUS-V2
           |  |--CONUS-V3
           |  |--Global-V3
           |  |--Radar
           |  |--WRF
           |--WRF_BC
           |--WRF_DART_utilities
         

Other documentation
-------------------

Additional documentation which didn't fit neatly into the other categories.

-  :doc:`./html/Manhattan_release`
-  :doc:`./html/Manhattan_diffs_from_Lanai`
-  :doc:`./html/mpi_intro`
-  :doc:`./html/filter_async_modes`
-  :doc:`../build_templates/mkmf`
-  :doc:`./tutorial/index`
-  :doc:`./DART_LAB/DART_LAB`

Complete list of all documentation
----------------------------------

The kitchen sink - quick links to all existing html docs plus all model_mod source files in the DART distribution tree:

-  :doc:`../models/9var/model_mod`
-  :doc:`../observations/obs_converters/AIRS/AIRS`
-  `models/am2/model_mod.f90 <../models/am2/model_mod.f90>`__
-  `convert_aura.f90 <../observations/obs_converters/AURA/convert_aura.f90>`__
-  :doc:`../observations/obs_converters/AVISO/AVISO`
-  :doc:`../observations/obs_converters/Ameriflux/level4_to_obs`
-  :doc:`../models/CESM/doc/setup_guidelines`
-  `CHAMP_density_text_to_obs.f90 <../observations/obs_converters/CHAMP/CHAMP_density_text_to_obs.f90>`__
-  `CNOFS_text_to_obs.f90 <../observations/obs_converters/CNOFS/CNOFS_text_to_obs.f90>`__
-  `models/coamps_nest/model_mod.f90 <../models/coamps_nest/model_mod.f90>`__
-  :doc:`../observations/obs_converters/COSMOS/COSMOS_development`
-  :doc:`../observations/obs_converters/COSMOS/COSMOS_to_obs`
-  :doc:`../docs/DART_LAB/DART_LAB`
-  :doc:`../observations/forward_operators/DEFAULT_obs_def_mod`
-  :doc:`../assimilation_code/modules/observations/DEFAULT_obs_kind_mod`
-  :doc:`../observations/obs_converters/DWL/dwl_to_obs`
-  `even_sphere.m <../observations/obs_converters/even_sphere/even_sphere.m>`__
-  `text_GITM/text_to_obs.f90 <../observations/obs_converters/text_GITM/text_to_obs.f90>`__
-  `gsi_to_dart.f90 <../observations/obs_converters/GSI2DART/gsi_to_dart.f90>`__
-  :doc:`../observations/obs_converters/GTSPP/GTSPP`
-  `gnd_gps_vtec_text_to_obs.f90 <../observations/obs_converters/gnd_gps_vtec/gnd_gps_vtec_text_to_obs.f90>`__
-  :doc:`../docs/html/Lanai_diffs_from_Kodiak`
-  :doc:`../docs/html/Lanai_release`
-  `models/lorenz_96_2scale/model_mod.f90 <../models/lorenz_96_2scale/model_mod.f90>`__
-  :doc:`../observations/obs_converters/MADIS/MADIS`
-  :doc:`../observations/obs_converters/MIDAS/MIDAS_to_obs`
-  `models/MITgcm_annulus/model_mod.f90 <../models/MITgcm_annulus/model_mod.f90>`__
-  :doc:`../models/MITgcm_ocean/create_ocean_obs`
-  :doc:`../models/MITgcm_ocean/model_mod`
-  :doc:`../models/MITgcm_ocean/trans_pv_sv`
-  :doc:`../models/MITgcm_ocean/trans_sv_pv`
-  :doc:`../observations/obs_converters/snow/snow_to_obs`
-  :doc:`../observations/obs_converters/MODIS/MOD15A2_to_obs`
-  :doc:`../observations/obs_converters/MODIS/MODIS_README`
-  :doc:`../docs/html/Manhattan_diffs_from_Lanai`
-  :doc:`../docs/html/Manhattan_getting_started`
-  :doc:`../docs/html/Manhattan_release`
-  :doc:`../observations/obs_converters/NCEP/prep_bufr/prep_bufr`
-  :doc:`../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs`
-  :doc:`../observations/obs_converters/NCEP/prep_bufr/prep_bufr`
-  :doc:`../models/NCOMMAS/dart_to_ncommas`
-  :doc:`../models/NCOMMAS/model_mod`
-  :doc:`../models/NCOMMAS/ncommas_to_dart`
-  :doc:`../models/POP/dart_pop_mod`
-  :doc:`../models/POP/model_mod`
-  :doc:`../models/POP/model_mod_check`
-  :doc:`../developer_tests/utilities/PrecisionCheck`
-  :doc:`../observations/obs_converters/ROMS/ROMS`
-  :doc:`../models/ROMS/model_mod`
-  `models/rose/model_mod.f90 <../models/rose/model_mod.f90>`__
-  `convert_saber_cdf.f90 <../observations/obs_converters/SABER/convert_saber_cdf.f90>`__
-  :doc:`../observations/obs_converters/SSEC/SSEC`
-  :doc:`../observations/obs_converters/SSUSI/convert_f16_edr_dsk`
-  `models/tiegcm/model_mod.f90 <../models/tiegcm/model_mod.f90>`__
-  :doc:`../observations/obs_converters/tpw/tpw`
-  :doc:`../observations/obs_converters/tropical_cyclone/tc_to_obs`
-  :doc:`./tutorial/index`
-  :doc:`../observations/obs_converters/WOD/WOD`
-  :doc:`../assimilation_code/modules/assimilation/adaptive_inflate_mod`
-  :doc:`../assimilation_code/programs/advance_time/advance_time`
-  `assert_mod.f90 <../assimilation_code/modules/utilities/assert_mod.f90>`__
-  :doc:`../assimilation_code/modules/assimilation/assim_model_mod`
-  :doc:`../assimilation_code/modules/assimilation/assim_tools_mod`
-  :doc:`../models/bgrid_solo/model_mod`
-  :doc:`../docs/html/bitwise_considerations`
-  :doc:`../docs/html/boilerplate/boilerplate`
-  :doc:`../models/cam-fv/model_mod`
-  :doc:`../models/cam-old/cam_to_dart`
-  :doc:`../models/cam-old/dart_to_cam`
-  :doc:`../models/cam-old/model_mod`
-  :doc:`../assimilation_code/location/channel/location_mod`
-  :doc:`../observations/obs_converters/cice/cice_to_obs`
-  `models/cice/model_mod.f90 <../models/cice/model_mod.f90>`__
-  :doc:`../models/clm/model_mod`
-  :doc:`../assimilation_code/programs/closest_member_tool/closest_member_tool`
-  :doc:`../models/cm1/model_mod`
-  :doc:`../models/coamps/model_mod`
-  :doc:`../assimilation_code/programs/compare_states/compare_states`
-  :doc:`../assimilation_code/programs/compute_error/compute_error`
-  :doc:`../assimilation_code/modules/assimilation/cov_cutoff_mod`
-  `cray_win_mod.f90 <../assimilation_code/modules/utilities/cray_win_mod.f90>`__
-  :doc:`../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq`
-  :doc:`../assimilation_code/programs/create_obs_sequence/create_obs_sequence`
-  `dart_time_io_mod.f90 <../assimilation_code/modules/io/dart_time_io_mod.f90>`__
-  `direct_netcdf_mod.f90 <../assimilation_code/modules/io/direct_netcdf_mod.f90>`__
-  :doc:`../docs/html/distributed_state`
-  `distributed_state_mod.f90 <../assimilation_code/modules/utilities/distributed_state_mod.f90>`__
-  `models/dynamo/model_mod.f90 <../models/dynamo/model_mod.f90>`__
-  :doc:`../assimilation_code/modules/utilities/ensemble_manager_mod`
-  `adaptive_inflate_mod.f90 <../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90>`__
-  :doc:`../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart`
-  :doc:`../assimilation_code/programs/filter/filter`
-  :doc:`../docs/html/filter_async_modes`
-  :doc:`../assimilation_code/modules/assimilation/filter_mod`
-  `models/forced_barot/model_mod.f90 <../models/forced_barot/model_mod.f90>`__
-  :doc:`../models/forced_lorenz_96/model_mod`
-  :doc:`../docs/html/forward_operator`
-  `forward_operator_mod.f90 <../assimilation_code/modules/observations/forward_operator_mod.f90>`__
-  :doc:`../assimilation_code/programs/gen_sampling_err_table/gen_sampling_err_table`
-  :doc:`../docs/html/generating_ensemble_ics`
-  :doc:`../docs/html/generating_obs_sequence`
-  :doc:`../models/gitm/dart_to_gitm`
-  :doc:`../models/gitm/gitm_to_dart`
-  :doc:`../models/gitm/model_mod`
-  :doc:`../observations/obs_converters/gps/gps`
-  :doc:`../docs/html/history/Fiji_release`
-  :doc:`../docs/html/history/Guam_release`
-  :doc:`../docs/html/history/I_diffs_from_workshop`
-  :doc:`../docs/html/history/Iceland_release`
-  :doc:`../docs/html/history/Jamaica_diffs_from_I`
-  :doc:`../docs/html/history/Jamaica_release`
-  :doc:`../docs/html/history/Kodiak_release`
-  :doc:`../docs/html/history/PostI_diffs_from_I`
-  :doc:`../docs/html/history/Post_Iceland_release`
-  :doc:`../docs/html/history/hawaii_release`
-  :doc:`../docs/html/history/pre_guam_release`
-  :doc:`../docs/html/history/pre_hawaii_release`
-  :doc:`../docs/html/history/pre_j_release`
-  :doc:`../models/ikeda/model_mod`
-  :doc:`../assimilation_code/programs/integrate_model/integrate_model`
-  `io_filenames_mod.f90 <../assimilation_code/modules/io/io_filenames_mod.f90>`__
-  :doc:`../assimilation_code/location/location_mod`
-  :doc:`../models/lorenz_04/model_mod`
-  :doc:`../models/lorenz_63/model_mod`
-  :doc:`../docs/html/lorenz_63_example`
-  :doc:`../models/lorenz_84/model_mod`
-  :doc:`../models/lorenz_96/model_mod`
-  :doc:`../build_templates/mkmf`
-  :doc:`../assimilation_code/programs/model_mod_check/model_mod_check`
-  :doc:`../models/mpas_atm/model_mod`
-  :doc:`../models/mpas_atm/mpas_dart_obs_preprocess`
-  :doc:`../models/mpas_ocn/model_mod`
-  :doc:`../models/mpas_ocn/model_to_dart`
-  :doc:`../docs/html/mpi_intro`
-  :doc:`../assimilation_code/modules/utilities/mpi_utilities_mod`
-  :doc:`../docs/html/netcdf_inflation_files`
-  :doc:`../models/noah/dart_to_noah`
-  :doc:`../models/noah/model_mod`
-  :doc:`../models/noah/noah_to_dart`
-  :doc:`../models/null_model/model_mod`
-  :doc:`../assimilation_code/programs/obs_common_subset/obs_common_subset`
-  :doc:`../observations/forward_operators/obs_def_1d_state_mod`
-  `obs_def_AIRS_mod.f90 <../observations/forward_operators/obs_def_AIRS_mod.f90>`__
-  `obs_def_AOD_mod.f90 <../observations/forward_operators/obs_def_AOD_mod.f90>`__
-  `obs_def_AURA_mod.f90 <../observations/forward_operators/obs_def_AURA_mod.f90>`__
-  `obs_def_COSMOS_mod.f90 <../observations/forward_operators/obs_def_COSMOS_mod.f90>`__
-  `obs_def_CO_Nadir_mod.f90 <../observations/forward_operators/obs_def_CO_Nadir_mod.f90>`__
-  `obs_def_GWD_mod.f90 <../observations/forward_operators/obs_def_GWD_mod.f90>`__
-  `obs_def_QuikSCAT_mod.f90 <../observations/forward_operators/obs_def_QuikSCAT_mod.f90>`__
-  `obs_def_SABER_mod.f90 <../observations/forward_operators/obs_def_SABER_mod.f90>`__
-  `obs_def_TES_nadir_mod.f90 <../observations/forward_operators/obs_def_TES_nadir_mod.f90>`__
-  `obs_def_altimeter_mod.f90 <../observations/forward_operators/obs_def_altimeter_mod.f90>`__
-  `obs_def_cice_mod.f90 <../observations/forward_operators/obs_def_cice_mod.f90>`__
-  `obs_def_cloud_mod.f90 <../observations/forward_operators/obs_def_cloud_mod.f90>`__
-  `obs_def_cwp_mod.f90 <../observations/forward_operators/obs_def_cwp_mod.f90>`__
-  :doc:`../observations/forward_operators/obs_def_dew_point_mod`
-  `obs_def_dwl_mod.f90 <../observations/forward_operators/obs_def_dwl_mod.f90>`__
-  `obs_def_eval_mod.f90 <../observations/forward_operators/obs_def_eval_mod.f90>`__
-  :doc:`../observations/forward_operators/obs_def_gps_mod`
-  `obs_def_gts_mod.f90 <../observations/forward_operators/obs_def_gts_mod.f90>`__
-  `obs_def_metar_mod.f90 <../observations/forward_operators/obs_def_metar_mod.f90>`__
-  :doc:`../observations/forward_operators/obs_def_mod`
-  :doc:`../observations/forward_operators/obs_def_ocean_mod`
-  `obs_def_pe2lyr_mod.f90 <../observations/forward_operators/obs_def_pe2lyr_mod.f90>`__
-  :doc:`../observations/forward_operators/obs_def_radar_mod`
-  `obs_def_radiance_mod.f90 <../observations/forward_operators/obs_def_radiance_mod.f90>`__
-  `obs_def_reanalysis_bufr_mod.f90 <../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90>`__
-  `obs_def_rel_humidity_mod.f90 <../observations/forward_operators/obs_def_rel_humidity_mod.f90>`__
-  `obs_def_simple_advection_mod.f90 <../observations/forward_operators/obs_def_simple_advection_mod.f90>`__
-  `obs_def_sqg_mod.f90 <../observations/forward_operators/obs_def_sqg_mod.f90>`__
-  `obs_def_surface_mod.f90 <../observations/forward_operators/obs_def_surface_mod.f90>`__
-  `obs_def_tower_mod.f90 <../observations/forward_operators/obs_def_tower_mod.f90>`__
-  `obs_def_tpw_mod.f90 <../observations/forward_operators/obs_def_tpw_mod.f90>`__
-  `obs_def_upper_atm_mod.f90 <../observations/forward_operators/obs_def_upper_atm_mod.f90>`__
-  `obs_def_vortex_mod.f90 <../observations/forward_operators/obs_def_vortex_mod.f90>`__
-  `obs_def_wind_speed_mod.f90 <../observations/forward_operators/obs_def_wind_speed_mod.f90>`__
-  :doc:`../assimilation_code/programs/obs_diag/oned/obs_diag`
-  :doc:`../assimilation_code/programs/obs_diag/threed_cartesian/obs_diag`
-  :doc:`../assimilation_code/programs/obs_diag/threed_sphere/obs_diag`
-  `obs_impact_mod.html <../assimilation_code/modules/utilities/obs_impact_mod.f90>`__
-  :doc:`../assimilation_code/programs/obs_impact_tool/obs_impact_tool`
-  :doc:`../assimilation_code/modules/observations/obs_kind_mod`
-  :doc:`../assimilation_code/modules/assimilation/obs_model_mod`
-  :doc:`../assimilation_code/programs/obs_selection/obs_selection`
-  :doc:`../assimilation_code/programs/obs_seq_coverage/obs_seq_coverage`
-  :doc:`../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf`
-  :doc:`../assimilation_code/programs/obs_seq_verify/obs_seq_verify`
-  :doc:`../assimilation_code/modules/observations/obs_sequence_mod`
-  :doc:`../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`
-  :doc:`../observations/obs_converters/observations`
-  :doc:`../observations/obs_converters/ok_mesonet/ok_mesonet`
-  :doc:`../assimilation_code/location/oned/location_mod`
-  `parse_args_mod.f90 <../assimilation_code/modules/utilities/parse_args_mod.f90>`__
-  :doc:`../models/pe2lyr/model_mod`
-  :doc:`../assimilation_code/programs/perfect_model_obs/perfect_model_obs`
-  :doc:`../assimilation_code/programs/preprocess/preprocess`
-  :doc:`../observations/obs_converters/quikscat/QuikSCAT`
-  :doc:`../observations/obs_converters/radar/radar`
-  :doc:`../assimilation_code/modules/utilities/random_seq_mod`
-  :doc:`../assimilation_code/modules/assimilation/reg_factor_mod`
-  :doc:`../assimilation_code/programs/restart_file_tool/restart_file_tool`
-  :doc:`../docs/html/running_lorenz_63`
-  `sampling_error_correction_mod.f90 <../assimilation_code/modules/assimilation/sampling_error_correction_mod.f90>`__
-  :doc:`../assimilation_code/modules/utilities/schedule_mod`
-  :doc:`../models/simple_advection/model_mod`
-  `single_file_io_mod.f90 <../assimilation_code/modules/io/single_file_io_mod.f90>`__
-  :doc:`../assimilation_code/modules/assimilation/smoother_mod`
-  :doc:`../observations/obs_converters/snow/snow_to_obs`
-  `sort_mod.f90 <../assimilation_code/modules/utilities/sort_mod.f90>`__
-  :doc:`../observations/obs_converters/text/text_to_obs`
-  :doc:`../models/sqg/model_mod`
-  :doc:`../docs/html/state_structure`
-  `state_structure_mod.f90 <../assimilation_code/modules/io/state_structure_mod.f90>`__
-  `state_vector_io_mod.f90 <../assimilation_code/modules/io/state_vector_io_mod.f90>`__
-  :doc:`../assimilation_code/programs/system_simulation/system_simulation`
-  :doc:`../models/template/model_mod`
-  :doc:`../docs/html/boilerplate/template`
-  :doc:`../assimilation_code/location/threed_cartesian/location_mod`
-  :doc:`../assimilation_code/location/threed_sphere/location_mod`
-  :doc:`../models/tiegcm/model_mod`
-  :doc:`../assimilation_code/modules/utilities/time_manager_mod`
-  :doc:`../observations/obs_converters/tpw/tpw`
-  :doc:`../observations/obs_converters/tropical_cyclone/tc_to_obs`
-  :doc:`../docs/tutorial/index`
-  :doc:`../assimilation_code/modules/utilities/types_mod`
-  :doc:`../assimilation_code/modules/utilities/utilities_mod`
-  :doc:`../observations/obs_converters/var/littler_tf_dart`
-  :doc:`../observations/obs_converters/var/rad_3dvar_to_dart`
-  :doc:`../observations/obs_converters/var/var`
-  :doc:`../docs/html/vertical_conversion`
-  :doc:`../assimilation_code/programs/wakeup_filter/wakeup_filter`
-  :doc:`../models/wrf/WRF_DART_utilities/dart_to_wrf`
-  :doc:`../models/wrf/WRF_DART_utilities/replace_wrf_fields`
-  :doc:`../models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess`
-  :doc:`../models/wrf/model_mod`
-  :doc:`../models/wrf/shell_scripts/advance_model`
-  `xyz_location_mod.html <../assimilation_code/location/threed_cartesian/xyz_location_mod.f90>`__
