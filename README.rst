Welcome to the Data Assimilation Research Testbed
=================================================

The Data Assimilation Research Testbed (DART) is an open-source, freely
available community facility for ensemble data assimilation (DA). [1]_ DART is
developed and maintained by the `Data Assimilation Research Section
(DAReS) <https://dart.ucar.edu/about/>`_ at the `National Center
for Atmospheric Research (NCAR) <https://ncar.ucar.edu>`_.

+------------------------------+------------------------------+
| |spaghetti_square|           | |assim_anim|                 |
+------------------------------+------------------------------+

Ensemble Data Assimilation
--------------------------

Ensemble DA is a technique for combining observations with numerical models
to estimate the state of a physical system.

It enables modelers, observational scientists, and geophysicists to:

- Generate initial conditions for forecasts.
- Create a retrospective estimate of the state of a system, a practice known as
  producing a *reanalysis*.
- Assess the relative value of specific observations on forecast skill, a
  practice known as conducting an *observing system experiment (OSE)*.
- Estimate the value of hypothetical observations in order to inform the design
  of an observing system, a practice known as conducting an *observing system
  simulation experiment (OSSE)*. 
- Determine a model's systematic bias in estimating the state of a system, a
  practice known as diagnosing *model error*.

The DART software environment makes it easy to explore a variety of data
assimilation methods and observations with different numerical models. It
provides powerful, flexible DA tools that are easy to use and customize to
support efficient and reliable DA applications. While DART is primarily
oriented for DA research, it has also been used in operational settings.

DART includes:

- A comprehensive tutorial introducing the concepts of ensemble DA.
- Extensive documentation of its source code.
- Interfaces to a variety of models and observation sets that can be used to
  introduce new users or graduate students to ensemble DA.

DART is also designed to facilitate the combination of assimilation algorithms,
models, and real or synthetic observations to allow increased
understanding of all three. It provides a framework for developing, testing,
and distributing advances in ensemble DA to a broad community of users by
removing the implementation-specific peculiarities of one-off DA systems.

These tools are intended for use by the full range of geosciencies community:
beginners and experts; students and teachers; national centers and university
research labs.

Organization of the documentation
---------------------------------

Because of DART's extensive scope, this documentation is detailed and
carefully organized, enabling you to easily find the information you need. If
you have any questions or suggestions for improvements, please contact DAReS
staff by emailing dart@ucar.edu.

The documentation is partitioned into three parts:

- a user guide that explains how to install DART and perform data assimilation
- source code documentation that provides a detailed description of the
  programs and modules in the repository
- a comprehensive description of data assimilation theory

Manhattan Release
-----------------

The Manhattan release is new and currently supports only a subset of the
models. DAReS staff will port over any requested model. Email dart@ucar.edu
if yours is not on the list.

For more information on this release, see :doc:`guide/Manhattan_release`.

Quick-start
-----------

1. fork the NCAR/DART repo
2. clone your (new) fork to your machine - this will set up a remote named
   ‘origin’
3. create a remote to point back to the NCAR/DART repo … convention dictates
   that this remote should be called ‘upstream’
4. check out the appropriate branch
5. Download one of the tar files (listed below) of ‘large’ files so you can test
   your DART installation.
6. If you want to issue a PR, create a feature branch and push that to your fork
   and issue the PR.

There are several large files that are needed to run some of the tests and
examples but are not included in order to keep the repository as small as
possible. If you are interested in running *bgrid_solo*, *cam-fv*, or testing
the *NCEP/prep_bufr* observation converter, you will need these files. These
files are available at:

+-------------------+------+----------------------------------------------------------------------------------------------------------------------------------+
| Release           | Size | Filename                                                                                                                         |
+===================+======+==================================================================================================================================+
| “Manhattan”       | 189M | `Manhattan_large_files.tar.gz <https://www.image.ucar.edu/pub/DART/Release_datasets/Manhattan_large_files.tar.gz>`__             |
+-------------------+------+----------------------------------------------------------------------------------------------------------------------------------+
| “wrf-chem.r13172” | 141M | `wrf-chem.r13172_large_files.tar.gz <https://www.image.ucar.edu/pub/DART/Release_datasets/wrf-chem.r13172_large_files.tar.gz>`__ |
+-------------------+------+----------------------------------------------------------------------------------------------------------------------------------+
| “Lanai”           | 158M | `Lanai_large_files.tar.gz <https://www.image.ucar.edu/pub/DART/Release_datasets/Lanai_large_files.tar.gz>`__                     |
+-------------------+------+----------------------------------------------------------------------------------------------------------------------------------+
| “Kodiak”          | 158M | `Kodiak_large_files.tar.gz <https://www.image.ucar.edu/pub/DART/Release_datasets/Kodiak_large_files.tar.gz>`__                   |
+-------------------+------+----------------------------------------------------------------------------------------------------------------------------------+
| “Jamaica”         | 32M  | `Jamaica_large_files.tar.gz <https://www.image.ucar.edu/pub/DART/Release_datasets/Jamaica_large_files.tar.gz>`__                 |
+-------------------+------+----------------------------------------------------------------------------------------------------------------------------------+
| “Hawaii”          | 32M  | `Hawaii_large_files.tar.gz <https://www.image.ucar.edu/pub/DART/Release_datasets/Hawaii_large_files.tar.gz>`__                   |
+-------------------+------+----------------------------------------------------------------------------------------------------------------------------------+

Download the appropriate tar file and untar it into your DART repository. Ignore
any warnings about ``tar: Ignoring unknown extended header keyword``.

Go into the ``build_templates`` directory and copy over the closest
``mkmf.template``._compiler.system\_ file into ``mkmf.template``.

Edit it to set the NETCDF directory location if not in ``/usr/local`` or comment
it out and set $NETCDF in your environment. *This NetCDF library must have been
compiled with the same compiler that you use to compile DART and must include
the F90 interfaces.*

Go into ``models/lorenz_63/work`` and run *quickbuild.csh*.

.. code-block::

   $ cd models/lorenz_63/work
   $ ./quickbuild.csh

If it compiles, run this series of commands to do a very basic test:

.. code-block::

   $ ./perfect_model_obs
   $ ./filter

If that runs and you have Matlab installed on your system add
``DART/diagnostics/matlab`` to your matlab search path and run the
``plot_total_err`` diagnostic script while in the ``models/lorenz_63/work``
directory. If the output plots and looks reasonable (error level stays around 2
and doesn’t grow unbounded) you have successfully installed DART and completed
your first assimilation with it.

If you are planning to run one of the larger models and want to use the Lorenz
63 model as a test, run ``./quickbuild.csh -mpi``. It will build filter and any
other MPI-capable executables with MPI.

.. important::

   The ``mpif90`` command you use must have been built with the same version of
   the compiler as you are using.

If any of these steps fail or you don’t know how to do them, go to the DART
project web page listed above for very detailed instructions that should get you
over any bumps in the process.

Citing DART
-----------

Cite DART using the following text:

   The Data Assimilation Research Testbed (Version X.Y.Z) [Software]. (2019).
   Boulder, Colorado: UCAR/NCAR/CISL/DAReS. http://doi.org/10.5065/D6WQ0202

Update the DART version and year as appropriate.

References
----------

.. [1] Anderson, J. L., T. Hoar, K. Raeder, H. Liu, N. Collins, R. Torn and A.
       Arellano, 2009 The Data Assimilation Research Testbed: A Community 
       Facility. *Bulletin of the American Meteorological Society*, **90**,
       1283-1296, `doi:10.1175/2009BAMS2618.1
       <http://dx.doi.org/10.1175/2009BAMS2618.1>`_

.. |spaghetti_square| image:: ./guide/images/DARTspaghettiSquare.gif
   :width: 100%

.. |assim_anim| image:: ./guide/images/science_nuggets/AssimAnim.gif
   :width: 100%

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   guide/system-requirements
   guide/fortran-compiler
   guide/locating-netcdf-library
   guide/downloading-dart
   guide/compiling-dart
   guide/verifying-installation

.. toctree::
   :maxdepth: 2
   :caption: Intermediate topics

   guide/assimilation-complex-model
   guide/mpi_intro
   guide/filters
   guide/inflation
   guide/Radiance_support
   guide/controlling-files-output

.. toctree::
   :maxdepth: 2
   :caption: Observations

   guide/preprocess-program
   guide/introduction-obs-seq-file
   guide/detailed-structure-obs-seq
   guide/creating-obs-seq-synthetic
   guide/creating-obs-seq-real
   guide/available-observation-converters
   guide/manipulating-with-obs-sequence-tool
   guide/difference-between-type-and-quantity
   guide/adding-support-new-type

.. toctree::
   :maxdepth: 2 
   :caption: Observation Converters

   observations/obs_converters/README


.. toctree::
   :hidden:

   observations/obs_converters/AIRS/AIRS
   observations/obs_converters/AIRS/README
   observations/obs_converters/AVISO/AVISO
   observations/obs_converters/Ameriflux/level4_to_obs
   observations/obs_converters/CHAMP/work/README
   observations/obs_converters/cice/cice_to_obs
   observations/obs_converters/CONAGUA/README
   observations/obs_converters/COSMOS/COSMOS_to_obs
   observations/obs_converters/COSMOS/COSMOS_development
   observations/obs_converters/DWL/dwl_to_obs
   observations/obs_converters/GMI/README
   observations/obs_converters/GOES/README
   observations/obs_converters/GPSPW/README
   observations/obs_converters/GSI2DART/readme
   observations/obs_converters/GTSPP/GTSPP
   observations/obs_converters/MADIS/MADIS
   observations/obs_converters/MIDAS/MIDAS_to_obs
   observations/obs_converters/MODIS/readme
   observations/obs_converters/MODIS/MOD15A2_to_obs
   observations/obs_converters/MPD/README
   observations/obs_converters/NCEP/prep_bufr/prep_bufr
   observations/obs_converters/NCEP/ascii_to_obs/create_real_obs
   observations/obs_converters/ROMS/ROMS
   observations/obs_converters/SSEC/SSEC
   observations/obs_converters/SST/SST
   observations/obs_converters/SSUSI/convert_f16_edr_dsk
   observations/obs_converters/WOD/WOD
   observations/obs_converters/gnd_gps_vtec/README
   observations/obs_converters/gps/gps
   observations/obs_converters/ok_mesonet/ok_mesonet
   observations/obs_converters/quikscat/QuikSCAT
   observations/obs_converters/even_sphere/README
   observations/obs_converters/obs_error/README
   observations/obs_converters/radar/radar
   observations/obs_converters/snow/snow_to_obs
   observations/obs_converters/text/text_to_obs
   observations/obs_converters/tpw/tpw
   observations/obs_converters/tropical_cyclone/tc_to_obs
   observations/obs_converters/var/littler_tf_dart
   observations/obs_converters/var/rad_3dvar_to_dart
   observations/obs_converters/var/var
   

.. toctree::
   :maxdepth: 2
   :caption: Diagnostics

   guide/checking-your-assimilation
   guide/computing-filter-increments
   guide/how-does-output-differ-from-input-increments
   guide/dart-missing-data-value
   guide/dart-quality-control
   guide/examining-obs-seq-final
   guide/matlab-observation-space

.. toctree::
   :maxdepth: 2
   :caption: Theory

   theory/readme
   theory/conditional-probability-bayes-theorem
   guide/DART_LAB/DART_LAB

.. toctree::
   :maxdepth: 2
   :caption: Releases

   guide/Manhattan_release
   guide/Lanai_release
   guide/history/Kodiak_release
   guide/history/Jamaica_release
   guide/history/Iceland_release
   guide/history/hawaii_release
   guide/history/Guam_release
   guide/history/Fiji_release

.. toctree::
   :maxdepth: 2
   :caption: Models

   models/9var/readme
   models/am2/readme
   models/bgrid_solo/readme
   models/cam-fv/readme
   models/CESM/readme
   models/cice/readme
   models/clm/readme
   models/cm1/readme
   models/coamps_nest/readme
   models/coamps/readme
   models/ECHAM/readme
   models/FESOM/readme
   models/gitm/readme
   models/ikeda/readme
   models/LMDZ/readme
   models/lorenz_04/readme
   models/lorenz_63/readme
   models/lorenz_84/readme
   models/lorenz_96/readme
   models/lorenz_96_2scale/readme
   models/forced_lorenz_96/readme
   models/MITgcm_ocean/readme
   models/mpas_atm/readme
   models/mpas_ocn/readme
   models/NCOMMAS/readme
   models/noah/readme
   models/PBL_1d/readme
   models/pe2lyr/readme
   models/POP/readme
   models/ROMS/readme
   models/rose/readme
   models/simple_advection/readme
   models/sqg/readme
   models/tiegcm/readme
   models/wrf_hydro/readme
   models/wrf/readme

.. toctree::
   :maxdepth: 2
   :caption: Contributing and Community

   guide/contributors-guide
   guide/requesting-features-reporting-bugs
   guide/mailing-list

.. toctree::
   :maxdepth: 2
   :caption: Guide

   guide/Manhattan_getting_started
   guide/rma
   guide/Manhattan_diffs_from_Lanai
   guide/forward_operator
   guide/boilerplate/boilerplate
   guide/boilerplate/template
   guide/vertical_conversion
   guide/bitwise_considerations
   guide/netcdf_inflation_files
   guide/state_structure
   guide/filter_async_modes
   guide/distributed_state

.. toctree::
   :maxdepth: 2
   :caption: History

   guide/Lanai_diffs_from_Kodiak   
   guide/history/Jamaica_diffs_from_I
   guide/history/pre_j_release
   guide/history/PostI_diffs_from_I
   guide/history/Post_Iceland_release
   guide/history/I_diffs_from_workshop
   guide/history/pre_hawaii_release
   guide/history/pre_guam_release

.. toctree::
   :maxdepth: 2
   :caption: Assimilation code

   assimilation_code/location/channel/location_mod
   assimilation_code/location/location_mod
   assimilation_code/location/oned/location_mod
   assimilation_code/location/threed_cartesian/location_mod
   assimilation_code/location/threed_sphere/location_mod
   assimilation_code/programs/obs_seq_verify/obs_seq_verify
   assimilation_code/programs/wakeup_filter/wakeup_filter
   assimilation_code/programs/compare_states/compare_states
   assimilation_code/programs/gen_sampling_err_table/gen_sampling_err_table
   assimilation_code/programs/perturb_single_instance/perturb_single_instance
   assimilation_code/programs/system_simulation/system_simulation
   assimilation_code/programs/compute_error/compute_error
   assimilation_code/programs/preprocess/preprocess
   assimilation_code/programs/obs_impact_tool/obs_impact_tool
   assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq
   assimilation_code/programs/obs_loop/obs_loop
   assimilation_code/programs/perfect_model_obs/perfect_model_obs
   assimilation_code/programs/obs_selection/obs_selection
   assimilation_code/programs/obs_sequence_tool/obs_sequence_tool
   assimilation_code/programs/integrate_model/integrate_model
   assimilation_code/programs/obs_diag/oned/obs_diag
   assimilation_code/programs/obs_diag/threed_cartesian/obs_diag
   assimilation_code/programs/obs_diag/threed_sphere/obs_diag
   assimilation_code/programs/fill_inflation_restart/fill_inflation_restart
   assimilation_code/programs/obs_seq_coverage/obs_seq_coverage
   assimilation_code/programs/advance_time/advance_time
   assimilation_code/programs/model_mod_check/model_mod_check
   assimilation_code/programs/closest_member_tool/closest_member_tool
   assimilation_code/programs/restart_file_tool/restart_file_tool
   assimilation_code/programs/filter/filter
   assimilation_code/programs/obs_keep_a_few/obs_keep_a_few
   assimilation_code/programs/create_obs_sequence/create_obs_sequence
   assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf
   assimilation_code/programs/obs_common_subset/obs_common_subset
   assimilation_code/modules/utilities/ensemble_manager_mod
   assimilation_code/modules/utilities/random_seq_mod
   assimilation_code/modules/utilities/mpi_utilities_mod
   assimilation_code/modules/utilities/time_manager_mod
   assimilation_code/modules/utilities/utilities_mod
   assimilation_code/modules/utilities/types_mod
   assimilation_code/modules/utilities/schedule_mod
   assimilation_code/modules/observations/obs_kind_mod
   assimilation_code/modules/observations/DEFAULT_obs_kind_mod
   assimilation_code/modules/observations/obs_sequence_mod
   assimilation_code/modules/assimilation/smoother_mod
   assimilation_code/modules/assimilation/assim_readme
   assimilation_code/modules/assimilation/assim_tools_mod
   assimilation_code/modules/assimilation/cov_cutoff_mod
   assimilation_code/modules/assimilation/obs_readme
   assimilation_code/modules/assimilation/reg_factor_mod
   assimilation_code/modules/assimilation/adaptive_inflate_mod
   assimilation_code/modules/assimilation/quality_control_mod
   assimilation_code/modules/assimilation/filter_mod

.. toctree::
   :maxdepth: 2
   :caption: Developer tests

   developer_tests/location/location_mod
   developer_tests/forward_operators/readme
   developer_tests/utilities/PrecisionCheck

.. toctree::
   :maxdepth: 2
   :caption: Observations

   observations/forward_operators/obs_def_gps_mod
   observations/forward_operators/obs_def_dew_point_mod
   observations/forward_operators/obs_def_ocean_mod
   observations/forward_operators/obs_def_1d_state_mod
   observations/forward_operators/obs_def_radar_mod
   observations/forward_operators/DEFAULT_obs_def_mod
   observations/forward_operators/obs_def_mod
   observations/forward_operators/obs_def_rttov_mod
   
.. toctree::
   :maxdepth: 2
   :caption: Misc

   models/null_model/readme
   models/NCOMMAS/dart_to_ncommas
   models/NCOMMAS/ncommas_to_dart
   models/POP/dart_pop_mod
   models/mpas_ocn/model_to_dart
   models/mpas_atm/mpas_dart_obs_preprocess
   models/CESM/doc/setup_guidelines
   models/template/readme
   models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_radon
   models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_sulfur_hex
   models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_driver
   models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_carbon_aerosol
   models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_utilities
   models/bgrid_solo/fms_src/atmos_shared/vert_advection/vert_advection
   models/bgrid_solo/fms_src/shared/time_manager/time_manager
   models/bgrid_solo/fms_src/shared/field_manager/field_manager
   models/bgrid_solo/fms_src/shared/horiz_interp/horiz_interp
   models/bgrid_solo/fms_src/shared/fms/fms
   models/bgrid_solo/fms_src/shared/constants/constants
   models/bgrid_solo/fms_src/shared/platform/platform
   models/bgrid_solo/fms_src/shared/utilities/utilities
   models/bgrid_solo/fms_src/shared/tracer_manager/tracer_manager
   models/bgrid_solo/fms_src/shared/mpp/mpp_domains
   models/bgrid_solo/fms_src/shared/mpp/mpp_io
   models/bgrid_solo/fms_src/shared/mpp/mpp
   models/bgrid_solo/fms_src/shared/fft/fft
   models/bgrid_solo/fms_src/shared/sat_vapor_pres/sat_vapor_pres
   models/bgrid_solo/fms_src/shared/topography/topography
   models/bgrid_solo/fms_src/shared/topography/gaussian_topog
   models/bgrid_solo/fms_src/shared/diag_manager/diag_manager
   models/bgrid_solo/fms_src/shared/diag_manager/diag_table_tk
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_polar_filter
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_halo
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_horiz
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_cold_start
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_prog_var
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_diagnostics
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_integrals
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_change_grid
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_masks
   models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_vert
   models/bgrid_solo/fms_src/atmos_bgrid/driver/solo/atmosphere
   models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core
   models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core_driver
   models/bgrid_solo/fms_src/atmos_param/hs_forcing/hs_forcing
   models/bgrid_solo/fms_src/atmos_solo/atmos_model
   models/wrf/WRF_DART_utilities/replace_wrf_fields
   models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess
   models/gitm/netcdf_to_gitm_blocks
   models/gitm/gitm_blocks_to_netcdf
   models/utilities/default_readme
   models/cam-old/cam_to_dart
   models/cam-old/readme
   models/cam-old/dart_to_cam
   models/MITgcm_ocean/trans_pv_sv
   models/MITgcm_ocean/create_ocean_obs
   models/MITgcm_ocean/trans_sv_pv

.. toctree::
   :maxdepth: 2
   :caption: Build templates

   build_templates/mkmf

.. toctree::
   :maxdepth: 2
   :caption: Root
   
   copyright
   changelog
   guide/404
