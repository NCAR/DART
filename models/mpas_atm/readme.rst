.. index:: mpas atm: mpas atmosphere; model for prediction across scales

MPAS_ATM
========

Overview
--------

This document describes the DART interface module for the atmospheric component 
of the Model for Prediction Across Scales
`MPAS <https://ncar.ucar.edu/what-we-offer/models/model-prediction-across-scales-mpas>`__ 
(or briefly, MPAS-ATM) global model.

The DART interface has constants set to match MPAS v5.0 onwards as defined in 
`MPAS-Model/src/framework/mpas_constants.F <https://github.com/MPAS-Dev/MPAS-Model/blob/master/src/framework/mpas_constants.F>`__.
If you need to reproduce work with DART and MPAS v4 you will need to change the model_mod.f90 
parameters ``cp``, ``cv`` and ``rvord`` to match MPAS v4. 


The mpas-atm model uses an unstructured Voronoi grid mesh,
formally Spherical Centriodal Voronoi Tesselations (SCVTs). This allows for both
quasi-uniform discretization of the sphere and local refinement. The MPAS/DART
interface was built on the SCVT-dual mesh and does not regrid to regular lat/lon
grids. In the C-grid discretization, the normal component of velocity on cell
edges is prognosed; zonal and meridional wind components are diagnosed on the
cell centers. We provide several options to choose from in the assimilation of
wind observations as shown below.

The grid terminology used in MPAS is as shown in the figure below:

|MPAS_grid_structure|

The wind options during a DART assimilation are controlled by combinations of 4
different namelist values. The values determine which fields the forward
operator uses to compute expected observation values; how the horizontal
interpolation is computed in that forward operator; and how the assimilation
increments are applied to update the wind quantities in the state vector.
Preliminary results based on real data assimilation experiments indicate that
performance is better when the zonal and meridional winds are used as input to
the forward operator that uses Barycentric interpolation, and when the
prognostic *u* wind is updated by the incremental method described in the figure
below. However there remain scientific questions about how best to handle the
wind fields under different situations. Thus we have kept all implemented
options available for use in experimental comparisons. See the figure below for
a flow-chart representation of how the 4 namelist items interact:

|WindDA_options|

Cycling of MPAS/DART is run in a *restart* mode. As for all DART experiments,
the overall design for an experiment is this: the DART program ``filter`` will
read the initial condition file, the observation sequence file, and the DART
namelist to decide whether or not to advance the MPAS-ATM model. All of the
control of the execution of the MPAS model is done by DART directly. If the
model needs to be advanced, ``filter`` makes a call to the shell to execute the
script ``advance_model.csh``, which is ENTIRELY responsible for getting all the
input files, data files, namelists, etc. into a temporary directory, running the
model, and copying the results back to the parent directory (which we call
CENTRALDIR). The whole process hinges on setting the MPAS-ATM model namelist
values such that it is doing a restart for every model advance. Unlike MPAS-ATM
free forecast runs, the forecast step in MPAS/DART requires to set up one more
namelist parameter called ``config_do_DAcycling = .true.`` in ``&restart``
section of ``namelist.input`` to recouple the state vectors (updated by filter)
with the mass field for the restart mode. For more information, check the
``advance_model.csh`` script in ``./shell_scripts/`` directory.

Since DART is an ensemble algorithm, there are multiple analysis files for a
single analysis time: one for each ensemble member. Because MPAS/DART is run in
a restart mode, each member should keep its own MPAS restart file from the
previous cycle (rather than having a single template file in CENTRALDIR).
Creating the initial ensemble of states is an area of active research.

Namelist
--------

The two namelists that are part of the MPAS-DART interface are
``&model_nml`` and ``&mpas_vars_nml``. The ``&model_nml`` namelist options control the 
behavior of the MPAS-DART interface. The ``&mpas_vars_nml`` namelist is used to 
define the MPAS variables that make up the DART state vector. 

The namelist is read from the file *input.nml*. Namelists start with an
ampersand '&' and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

model_nml
^^^^^^^^^

.. code-block:: fortran

   &model_nml
      init_template_filename       = 'mpas_init.nc',
      vert_localization_coord      = 3,
      assimilation_period_days     = 0,
      assimilation_period_seconds  = 21600,
      model_perturbation_amplitude = 0.0001,
      log_p_vert_interp            = .true.,
      calendar                     = 'Gregorian',
      use_u_for_wind               = .false.,
      use_rbf_option               = 2,
      update_u_from_reconstruct    = .true.,
      use_increments_for_u_update  = .true.,
      highest_obs_pressure_mb      = 100.0,
      sfc_elev_max_diff            = -1.0,
      outside_grid_level_tolerance = -1.0,
      write_grid_to_diag_files     = .false.,
      no_normalization_of_scale_heights = .true.

   /

+---------------------------------------+---------------------------------------+-----------------------------------------+
| Item                                  | Type                                  | Description                             |
+=======================================+=======================================+=========================================+
| init_template_filename                | character(len=256)                    | The name of the MPAS analysis file to   |
|                                       | *[default: 'mpas_init.nc']*           | be read and/or written by the DART      |
|                                       |                                       | programs for the state data.            |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| highest_obs_pressure_mb               | real(r8)                              | Observations higher than this           |
|                                       | *[default: 100.0]*                    | pressure are ignored. Set to -1.0 to    |
|                                       |                                       | ignore this test. For models with a     |
|                                       |                                       | prescribed top boundary layer, trying   |
|                                       |                                       | to assimilate very high observations    |
|                                       |                                       | results in problems because the model   |
|                                       |                                       | damps out any changes the               |
|                                       |                                       | assimilation tries to make. With        |
|                                       |                                       | adaptive algorithms this results in     |
|                                       |                                       | larger and larger coefficients as the   |
|                                       |                                       | assimilation tries to effect state      |
|                                       |                                       | vector change.                          |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| assimilation_period_days              | integer *[default: 0]*                | The number of days to advance the       |
|                                       |                                       | model for each assimilation. Even if    |
|                                       |                                       | the model is being advanced outside     |
|                                       |                                       | of the DART filter program, the         |
|                                       |                                       | assimilation period should be set       |
|                                       |                                       | correctly. Only observations with a     |
|                                       |                                       | time within +/- 1/2 this window size    |
|                                       |                                       | will be assimilated.                    |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| assimilation_period_seconds           | integer *[default: 21600]*            | In addition to                          |
|                                       |                                       | ``assimilation_period_days``, the       |
|                                       |                                       | number of seconds to advance the        |
|                                       |                                       | model for each assimilation.            |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| vert_localization_coord               | integer *[default: 3]*                | Vertical coordinate for vertical        |
|                                       |                                       | localization.                           |
|                                       |                                       |                                         |
|                                       |                                       | -  1 = model level                      |
|                                       |                                       | -  2 = pressure (in pascals)            |
|                                       |                                       | -  3 = height (in meters)               |
|                                       |                                       | -  4 = scale height (unitless)          |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| sfc_elev_max_diff                     | real(r8)\ *[default: -1.0]*           | If > 0, the maximum difference, in      |
|                                       |                                       | meters, between an observation marked   |
|                                       |                                       | as a 'surface obs' as the vertical      |
|                                       |                                       | type (with the surface elevation, in    |
|                                       |                                       | meters, as the numerical vertical       |
|                                       |                                       | location), and the surface elevation    |
|                                       |                                       | as defined by the model. Observations   |
|                                       |                                       | further away from the surface than      |
|                                       |                                       | this threshold are rejected and not     |
|                                       |                                       | assimilated. If the value is            |
|                                       |                                       | negative, this test is skipped.         |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| log_p_vert_interp                     | logical *[default: .true.]*           | If ``.true.``, vertical interpolation   |
|                                       |                                       | is done in log-pressure. Otherwise,     |
|                                       |                                       | linear.                                 |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| use_u_for_wind                        | logical *[default: .false.]*          | If ``.false.``, zonal and meridional    |
|                                       |                                       | winds at cell centers are used for      |
|                                       |                                       | the wind observation operator           |
|                                       |                                       | [default]. In that case, triangular     |
|                                       |                                       | meshes are used for the barycentric     |
|                                       |                                       | (e.g., area-weighted) interpolation.    |
|                                       |                                       | If ``.true.``, wind vectors at an       |
|                                       |                                       | arbitrary (e.g., observation) point     |
|                                       |                                       | are reconstructed from the normal       |
|                                       |                                       | component of velocity on cell edges     |
|                                       |                                       | *(u)* using radial basis functions      |
|                                       |                                       | (RBFs) provided by the MPAS model.      |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| use_rbf_option                        | integer *[default: 2]*                | If ``use_u_for_wind = .true.``, this    |
|                                       |                                       | option controls how many points will    |
|                                       |                                       | be used in the RBF interpolation.       |
|                                       |                                       | Options are available as 0, 1, 2, and   |
|                                       |                                       | 3. All the edges available in N (=      |
|                                       |                                       | 0,1,2, or 3) neighboring cells go       |
|                                       |                                       | into the RBF reconstruction.            |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| update_u_from_reconstruct             | logical *[default: .true.]*           | When zonal and meridional winds at      |
|                                       |                                       | cell centers are used for the wind      |
|                                       |                                       | observation operator                    |
|                                       |                                       | (``use_u_for_wind = .false.``), this    |
|                                       |                                       | option decides if the normal            |
|                                       |                                       | component of velocity on cell edges     |
|                                       |                                       | (which is the only wind prognostic      |
|                                       |                                       | variable in MPAS-ATM) should be         |
|                                       |                                       | updated from the winds at cell          |
|                                       |                                       | centers. If ``.true.``,                 |
|                                       |                                       | ``use_increments_for_u_update``         |
|                                       |                                       | should be also decided.                 |
|                                       |                                       | If ``use_u_for_wind = .true.``          |
|                                       |                                       | and the normal component of             |
|                                       |                                       | velocity on cell edges is defined as    |
|                                       |                                       | a state vector, this option should be   |
|                                       |                                       | ``.false.`` so the edge winds can be    |
|                                       |                                       | directly updated by filter.             |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| use_increments_for_u_update           | logical *[default: .true.]*           | Only if ``update_u_from_reconstruct     |
|                                       |                                       | = .true.``, this option is used to      |
|                                       |                                       | decide if the edge winds are replaced   |
|                                       |                                       | by averaging from the analysis winds    |
|                                       |                                       | at cell centers (``.false.``), or       |
|                                       |                                       | just updated by the analysis            |
|                                       |                                       | increments at cell centers              |
|                                       |                                       | (``.true.``). If ``.true.``, all        |
|                                       |                                       | the wind components (e.g., both at      |
|                                       |                                       | cell centers and edges) are read from   |
|                                       |                                       | prior and used to compute the           |
|                                       |                                       | increments [Recommended].               |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| model_perturbation_amplitude          | real(r8) *[default: 0.0001]*          | The amplitude of random noise to add    |
|                                       |                                       | when trying to perturb a single state   |
|                                       |                                       | vector to create an ensemble. Only      |
|                                       |                                       | used when ``start_from_restart =        |
|                                       |                                       | .false.`` in the ``&filter_nml``        |
|                                       |                                       | namelist within ``input.nml``           |
|                                       |                                       | Multiplied by the state vector, it      |
|                                       |                                       | produces standard deviation of a        |
|                                       |                                       | gaussian distribution with the mean     |
|                                       |                                       | at the value of the state vector        |
|                                       |                                       | element.                                |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| calendar                              | character(len=32)                     | Character string specifying the         |
|                                       | *[default: 'Gregorian']*              | calendar being used by MPAS.            |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| outside_grid_level_tolerance          | real(r8) *[default: -1.0]*            | If greater than 0.0, amount of          |
|                                       |                                       | distance in fractional model levels     |
|                                       |                                       | that a vertical location can be above   |
|                                       |                                       | or below the top or bottom of the       |
|                                       |                                       | grid and still be evaluated without     |
|                                       |                                       | error. Since *extrapolate* is not       |
|                                       |                                       | implemented yet, the value of           |
|                                       |                                       | ``.false.`` will be assumed. In this    |
|                                       |                                       | case, vertical locations equivalent     |
|                                       |                                       | to level 1 or level N will be used.     |
|                                       |                                       | Eventually, if *extrapolate* is         |
|                                       |                                       | ``.true.``, extrapolate from the        |
|                                       |                                       | first or last model level. If           |
|                                       |                                       | *extrapolate* is ``.false.``, simply    |
|                                       |                                       | use the value at level 1 for low        |
|                                       |                                       | vertical locations, or at level N for   |
|                                       |                                       | high vertical locations.                |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| write_grid_to_diag_files              | logical *[default: .false.]*          | If ``.true.``, write the grid           |
|                                       |                                       | information to netcdf files created     |
|                                       |                                       | by DART. Results in larger files.       |
+---------------------------------------+---------------------------------------+-----------------------------------------+
| no_normalization_of_scale_heights     | logical *[default: .true.]*           | When converting to scale height for the |
|                                       |                                       | vertical, set this to .false. to use    | 
|                                       |                                       | the log of the pressure. To normalize   |
|                                       |                                       | by the surface pressure (backwards      |
|                                       |                                       | compatible with previous code),         |
|                                       |                                       | set to .true.                           |
+---------------------------------------+---------------------------------------+-----------------------------------------+




mpas_vars_nml
^^^^^^^^^^^^^

.. code-block:: fortran

   &mpas_vars_nml
      mpas_state_variables = 'theta',                 'QTY_POTENTIAL_TEMPERATURE',
                             'uReconstructZonal',     'QTY_U_WIND_COMPONENT',
                             'uReconstructMeridional','QTY_V_WIND_COMPONENT',
                             'qv',                    'QTY_VAPOR_MIXING_RATIO',
                             'qc',                    'QTY_CLOUDWATER_MIXING_RATIO',
                             'surface_pressure',      'QTY_SURFACE_PRESSURE'
      mpas_state_bounds    = 'qv','0.0','NULL'
                             'qc','0.0','NULL'
   /

**mpas_state_variables**

 - The first column must match the exact NetCDF field name in the MPAS file.
 - The second column must be a :ref:`physical quantity<obs_kind_mod>` recognized by DART.

**mpas_state_bounds**

 - The first column must match the exact NetCDF field name in the MPAS file.
 - The second and third columns are the minimum and maximum values for the field.


``mpas_state_bounds`` is for variables with fixed minimum or maximum limits. 
When writing back to MPAS NetCDF files, out-of-range values are adjusted to be in bounds.
Note the adjustment is done only when writing to MPAS NetCDF files, not during assimilation.
DART files, mean, sd, and inflation files are not adjusted.

Note that changing values at the edges of the distribution means it is no longer completely Gaussian.
In practice this technique has worked effectively, but if the assimilation is continually trying to 
move the values outside the permitted range the results may be of poor quality. 
Examine the diagnostics for these fields carefully when using bounds to restrict their values.
You may want to consider using the :ref:`QCEFF<qceff>` filter when working with bounded quantities.


Grid Information
----------------

As the forward operators use the unstructured grid meshes in MPAS-ATM, the
DART/MPAS interface needs to read static variables related to the grid structure
from the MPAS ATM netCDF file (specified in ``init_template_filename``).
The calculations to find the closest mesh cell to an observation location is performed in the
cartesian coordinate to avoid the polar issues.

References
----------

The Data Assimilation section in the MPAS documentation found at
http://mpas-dev.github.io.

.. |MPAS_grid_structure| image:: ../../guide/images/MPAS_grid_structure.png

.. |WindDA_options| image:: ../../guide/images/MPAS_WindDA_options.png
