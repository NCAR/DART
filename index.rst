.. _Welcome page:

Welcome to the Data Assimilation Research Testbed
=================================================

The Data Assimilation Research Testbed (DART) is an open-source, freely
available community facility for ensemble data assimilation (DA). [1]_ DART is
developed and maintained by the `Data Assimilation Research Section
(DAReS) <https://dart.ucar.edu/about/>`_ at the NSF `National Center
for Atmospheric Research (NSF NCAR) <https://ncar.ucar.edu>`_.

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
- Nonlinear and Non-Gaussian DA Capabilities

DART is also designed to facilitate the combination of assimilation algorithms,
models, and real or synthetic observations to allow increased
understanding of all three. It provides a framework for developing, testing,
and distributing advances in ensemble DA to a broad community of users by
removing the implementation-specific peculiarities of one-off DA systems.

These tools are intended for use by the full range of geosciencies community:
beginners and experts; students and teachers; national centers and university
research labs.

Nonlinear and Non-Gaussian Data Assimilation Capabilities in DART
-----------------------------------------------------------------

One of the historical drawbacks of ensemble data assimilation techniques is the
assumption that the quantities being assimilated obey a normal distribution.
While this is often a safe assumption -- distributions of temperature and
pressure can be approximated using a normal distribution -- many quantities
such as precipitation, snow depth and tracer concentration, as well as many
model parameters aren't normally distributed.

Applying traditional ensemble data assimilation techniques in situations where
assumptions of gaussianity are invalid can lead to poor forecast skill and
inconclusive results.

To overcome these problems, DART now implements a novel data assimilation
technique that no longer requires quantities to be normally distributed. The
Quantile-Conserving Ensemble Filtering Framework :ref:`(QCEFF) <QCEFF>`
provides a general method of computing increments for the prior ensemble of an
observed quantity by allowing the use of arbitrary distributions for the prior
and the observation error. For a detailed description of the QCEFF, see
Anderson (2022). [2]_ 

While the QCEFF for computing observation increments can lead to significant improvements in 
analysis estimates for observed variables, those improvements can be lost when using standard 
linear regression of observation increments to update state variables. The QCEFF also 
implements a capability to do regression in a probit probability integral transformed space. 
Doing the regression of observation quantile increments in the transformed space guarantees 
that the posterior ensembles for state variables also retain the advantages of the observation space
quantile conserving posteriors. For example, if state variables are bounded, then posterior 
ensembles will respect the bounds. The posterior ensembles also respect other aspects of the 
continuous prior distributions. For a detailed description of this process, see
Anderson (2023) [3]_ and Anderson et al. (2023). [4]_

Inflation and localization, methods that improve the quality of ensemble DA, can also negate 
the advantages of the QCEFF methods. For this reason, both localization and inflation can be 
done in the probit-transformed quantile space as well.  Combining these new methods can 
significantly improve data assimilation for non-Gaussian quantities in Earth system models by 
extending the capabilities of ensemble DA to general non-Gaussian and nonlinear distributions. 
Transformative improvements in DA can result for many applications. The largest improvements are 
found for bounded variables like tracer concentrations, snow and ice depth, soil moisture, and 
similar quantities in other parts of the Earth system. Model parameters can also be estimated 
with DA and large improvements can occur for bounded parameters. Variables that have distinctly 
non-Gaussian prior distributions can also see large improvements. Examples can include atmospheric 
quantities like moisture and cloud amount in the presence of convection, and many land surface variables.

For instructions on how to use these tools, see :ref:`QCEFF`.

For step-by-step examples of the QCEFF tools, you can work through 
:ref:`examples with the Lorenz 96 tracer model <quantile tracer>`

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

The current version of DART is the Manhattan release. 
Email dart@ucar.edu for advice if you are interested in a model which has not been converted from the previous Lanai release.


Quick-start
-----------

DART is available through `GitHub <https://github.com/NCAR/DART>`__. To
download the latest version of DART, use:

.. code::

   git clone https://github.com/NCAR/DART.git


Go into the ``build_templates`` directory and copy over the closest
``mkmf.template``._compiler.system\_ file into ``mkmf.template``.

Edit it to set the NETCDF directory location if not in ``/usr/local`` or comment
it out and set $NETCDF in your environment. *This NetCDF library must have been
compiled with the same compiler that you use to compile DART and must include
the F90 interfaces.*

Go into ``models/lorenz_63/work`` and run *quickbuild.sh nompi*.

.. code-block::

   $ cd models/lorenz_63/work
   $ ./quickbuild.sh nompi

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
63 model as a test, run ``./quickbuild.sh``. It will build filter and any
other MPI-capable executables with MPI.

``./quickbuild.sh help`` will print out the quickbuild.sh usage.

.. important::

   The ``mpif90`` command you use must have been built with the same version of
   the compiler as you are using.

If any of these steps fail or you don’t know how to do them, go to the DART
project web page listed above for very detailed instructions that should get you
over any bumps in the process.

Quick-start for developers
^^^^^^^^^^^^^^^^^^^^^^^^^^

To create a fork of DART for your own development you will need
a `GitHub <https://github.com/>`__ account. 

1. fork the NCAR/DART repo on GitHub
2. clone your (new) fork to your machine - this will set up a remote named
   ‘origin’.

.. code::

   git clone https://github.com/USERNAME/DART.git

where `USERNAME` is your GitHub username. 

3. create a remote to point back to the NCAR/DART repo. Convention dictates
   that this remote should be called ‘upstream’

.. code::

   git remote add upstream https://github.com/NCAR/DART.git

Use ‘upstream’ to keep your fork up to date with NCAR/DART. GitHub has documentation
on `working with forks <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/working-with-forks>`__.

4. Download one of the tar files (listed below) of ‘large’ files so you can test
   your DART installation.
5. If you want to contribute your work back to the DART community, create a feature
   branch with your work, then issue a `pull request <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork>`__
   to propose changes to NCAR/DART.

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


Citing DART
-----------

Cite DART using the following text:

   The Data Assimilation Research Testbed (Version X.Y.Z) [Software]. (2024).
   Boulder, Colorado: UCAR/NSF NCAR/CISL/DAReS. http://doi.org/10.5065/D6WQ0202

Update the DART version and year as appropriate.

References
----------

.. [1] Anderson, J. L., T. Hoar, K. Raeder, H. Liu, N. Collins, R. Torn and A.
       Arellano, 2009 The Data Assimilation Research Testbed: A Community 
       Facility. *Bulletin of the American Meteorological Society*, **90**,
       1283-1296, `doi:10.1175/2009BAMS2618.1
       <http://dx.doi.org/10.1175/2009BAMS2618.1>`_
.. [2] Anderson, J. L., 2022: A Quantile-Conserving Ensemble Filter Framework.
       Part I: Updating an Observed Variable. *Monthly Weather Review*, **150**,
       1061–1074, `doi:10.1175/MWR-D-21-0229.1 <http://n2t.net/ark:/85065/d7mk6hm4>`_
.. [3] Anderson, J. L., 2023: A Quantile-Conserving Ensemble Filter Framework.
       Part II: Regression of Observation Increments in a Probit and
       Probability Integral Transformed Space. *Monthly Weather Review*,
       **151**, 2759–2777, `doi:10.1175/MWR-D-23-0065.1 <http://n2t.net/ark:/85065/d7nv9pbt>`_ 
.. [4] Anderson, J. L., Riedel, C., Wieringa, M., Ishraque, F., Smith, M., Kershaw, H.
       2023: A Quantile-Conserving
       Ensemble Filter Framework. Part III: Data Assimilation for Mixed Distributions
       with Application to a Low-Order Tracer Advection Model. *Monthly Weather Review*
       `[Manuscript submitted for publication] <_static/papers/QCEFF_3_submitted.pdf>`_

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
   :caption: What is data assimilation?

   guide/introduction-ensemble-da
   guide/lorenz-63-model
   guide/da-in-dart-with-lorenz-63

.. toctree::
   :maxdepth: 2
   :caption: What is DART?

   guide/what-is-dart
   guide/benefits-of-using-dart
   guide/brief-history-of-dart
   guide/high-level-da-workflows
   guide/dart-design-philosophy
   guide/important-capabilities-dart
   guide/qceff_probit

.. toctree::
   :maxdepth: 2
   :caption: Run DART with your model

   guide/advice-for-new-collaborators
   guide/instructions-for-porting-a-new-model-to-dart
   DART build system <guide/quickbuild.rst>
   guide/assimilation-complex-model
   guide/mpi_intro
   guide/filters
   guide/inflation
   guide/required-model-mod-routines
   guide/suggestions-for-a-simple-model
   guide/suggestions-for-a-complex-model
   guide/how-to-test-your-model-mod-routines
   guide/controlling-files-output
   guide/vertical-conversion
   guide/data-management-issues
   assimilation_code/programs/readme

.. toctree::
   :maxdepth: 2
   :caption: Observations

   guide/adding-your-observations-to-dart
   guide/preprocess-program
   guide/introduction-obs-seq-file
   guide/detailed-structure-obs-seq
   guide/creating-obs-seq-synthetic
   guide/creating-obs-seq-real
   guide/available-observation-converters
   guide/manipulating-with-obs-sequence-tool
   guide/difference-between-type-and-quantity
   guide/adding-support-new-type
   Radiances <guide/Radiance_support>

.. toctree::
   :maxdepth: 2 
   :caption: Observation Converters

   observations/obs_converters/README


.. toctree::
   :hidden:

   observations/obs_converters/AIRS/README
   observations/obs_converters/AIRS/convert_airs_L2
   observations/obs_converters/AIRS/convert_amsu_L1
   observations/obs_converters/AVISO/AVISO
   observations/obs_converters/Ameriflux/fluxnetfull_to_obs
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
   observations/obs_converters/NASA_Earthdata/README
   observations/obs_converters/NCEP/prep_bufr/prep_bufr
   observations/obs_converters/NCEP/ascii_to_obs/create_real_obs
   observations/obs_converters/NSIDC/SMAP_L2_to_obs
   observations/obs_converters/ROMS/ROMS
   observations/obs_converters/SIF/SIF_to_obs_netcdf
   observations/obs_converters/SSEC/SSEC
   observations/obs_converters/SST/SST
   observations/obs_converters/SSUSI/convert_f16_edr_dsk
   observations/obs_converters/WOD/WOD
   observations/obs_converters/gnd_gps_vtec/README
   observations/obs_converters/gps/gps
   observations/obs_converters/ok_mesonet/ok_mesonet
   observations/obs_converters/ocean_color/README
   observations/obs_converters/quikscat/QuikSCAT
   observations/obs_converters/even_sphere/README
   observations/obs_converters/obs_error/README
   observations/obs_converters/radar/README
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
   CLM-DART Tutorial <models/clm/tutorial/README>
   WRF-DART Tutorial <models/wrf/tutorial/README>
   guide/qceff-examples.rst
   
.. toctree::
   :maxdepth: 2
   :caption: Models

   models/README

.. toctree::
   :maxdepth: 2
   :caption: Models
   :hidden:

   models/9var/readme
   models/aether_lat-lon/readme
   models/am2/readme
   models/bgrid_solo/readme
   models/cam-fv/readme
   models/cam-fv/cam_dart_obs_preprocessor
   models/CESM/readme
   models/cice/readme
   models/clm/readme
   models/clm/clm_to_dart
   models/clm/dart_to_clm
   models/cm1/readme
   models/coamps_nest/readme
   models/coamps/readme
   models/ECHAM/readme
   models/FESOM/readme
   models/gitm/readme
   models/gitm/netcdf_to_gitm_blocks
   models/gitm/gitm_blocks_to_netcdf
   models/ikeda/readme
   models/LMDZ/readme
   models/lorenz_04/readme
   models/lorenz_63/readme
   models/lorenz_84/readme
   models/lorenz_96/readme
   models/lorenz_96_2scale/readme
   models/lorenz_96_tracer_advection/readme
   models/forced_lorenz_96/readme
   models/MITgcm_ocean/readme
   models/MOM6/readme
   models/mpas_atm/readme
   models/mpas_atm/mpas_dart_obs_preprocess
   models/mpas_ocn/readme
   models/mpas_ocn/model_to_dart
   models/NCOMMAS/readme
   models/noah/readme
   models/null_model/readme
   models/PBL_1d/readme
   models/pe2lyr/readme
   models/POP/readme
   models/POP/dart_pop_mod
   models/ROMS/readme
   models/rose/readme
   models/seir/readme
   models/simple_advection/readme
   models/sqg/readme
   models/template/new_model
   models/tiegcm/readme
   models/wrf_hydro/readme
   models/wrf/readme
   models/wrf/WRF_DART_utilities/replace_wrf_fields
   models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess
   models/template/readme
   models/utilities/default_model_mod


.. toctree::
   :maxdepth: 2
   :caption: Contributing and Community

   guide/contributors-guide
   guide/requesting-features-reporting-bugs
   guide/mailing-list

.. toctree::
   :maxdepth: 2
   :caption: Guide

   guide/Manhattan_diffs_from_Lanai
   guide/forward_operator
   guide/netcdf_inflation_files
   guide/state_structure
   guide/filter_async_modes
   guide/distributed_state


.. toctree::
   :maxdepth: 2
   :caption: Assimilation code
   :hidden:

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
   assimilation_code/modules/assimilation/assim_model_mod
   assimilation_code/modules/assimilation/assim_tools_mod
   assimilation_code/modules/assimilation/cov_cutoff_mod
   assimilation_code/modules/assimilation/obs_model_mod
   assimilation_code/modules/assimilation/reg_factor_mod
   assimilation_code/modules/assimilation/adaptive_inflate_mod
   assimilation_code/modules/assimilation/quality_control_mod
   assimilation_code/modules/assimilation/filter_mod

.. toctree::
   :maxdepth: 2
   :caption: Developer tests
   :hidden:

   developer_tests/location/location_mod
   developer_tests/forward_operators/readme
   developer_tests/utilities/PrecisionCheck

.. toctree::
   :maxdepth: 2
   :caption: Forward Operators
   :hidden:

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
   
   models/CESM/doc/setup_guidelines

   
.. toctree::   
   :caption: non-compiling models
   :hidden:
         
   models/MITgcm_ocean/create_ocean_obs
   models/NCOMMAS/dart_to_ncommas
   models/NCOMMAS/ncommas_to_dart


.. toctree::
   :maxdepth: 2
   :caption: Root
   
   copyright
   Changelog <CHANGELOG>

.. toctree::
   :hidden:

   guide/404
