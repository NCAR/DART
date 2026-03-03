pywatershed
==============

**pywatershed** is a python package for simulating hydrologic processes
motivated by the need to modernize important, legacy hydrologic models at the
USGS, particularly the Precipitation-Runoff Modeling System
(`PRMS <https://www.usgs.gov/software/precipitation-runoff-modeling-system-prms>`_).
pywatershed calculates explicit solutions of spatially distributed
hydrologic process representations including evaporation, transpiration, runoff,
infiltration, interflow, snowpack, soil moisture, conceptual groundwater storage,
and channel flow. These process representations simulate hydrologic response and
water budgets given inputs of spatially distributed weather variables and land
use change at temporal scales ranging from days to centuries.
For more information on the goals and status of pywatershed, please see the
`pywatershed docs <https://pywatershed.readthedocs.io/en/main/>`_.

Overview
--------
This software is aimed at building an interface between DART and pywatershed.
In its current form, it is designed to assimilate streamflow observations into
two distinct domains: [1] a Channel domain which includes streamflow, and [2] an
HRU (Hydrologic Response Unit) domain which includes groundwater. Groundwater
variables from more than one basin could contribute to a single streamflow
segment (i.e., stream or link). The localization code follows the
*Along-The-Stream Localization (ATS)* developed by
`El Gharamti et al. (2021) <https://doi.org/10.5194/hess-25-5315-2021>`_

Namelist
--------
The ``&model_nml`` variables and their default values are listed here:

.. code-block:: fortran

  &model_nml
      assimilation_period_days     = 1
      assimilation_period_seconds  = 0
      model_perturbation_amplitude = 0.05
      perturb_distribution         = 'lognormal'
      max_link_distance            = 100000.0
      domain_order                 = 'channel', 'hru'
      domain_shapefiles            = 'seg_inflow.nc', 'gwres_stor.nc'
      channel_config_file          = 'parameters_dis_seg_app.nc'
      hru_config_file              = 'parameters_dis_hru.nc'
      hydro_config_file            = 'parameters_PRMSChannel.nc'
      channel_variables            = 'seg_inflow', 'QTY_STREAM_FLOW', '0.0', 'NA', 'update',
      hru_variables                = 'gwres_stor', 'QTY_GROUNDWATER', '0.0', 'NA', 'update',
      debug                        = 0
      /

The namelist provides control over model variables, ensemble perturbation and
domain choices.

+-------------------------------------+--------------------+------------------------------------------------------------+
| Item                                | Type               | Description                                                |
+=====================================+====================+============================================================+
| assimilation_period_days            | integer            | Number of days to advance the model for each assimilation. |
+-------------------------------------+--------------------+------------------------------------------------------------+
| assimilation_period_seconds         | integer            | In addition to ``assimilation_period_days``, the number of |
|                                     |                    | seconds to advance the model for each assimilation.        |
+-------------------------------------+--------------------+------------------------------------------------------------+
| model_perturbation_amplitude        | real(r8)           | Perturbation parameter for initial ensemble generation.    |
+-------------------------------------+--------------------+------------------------------------------------------------+
| perturb_distribution                | character(len=256) | Distribution needed for initial ensemble generation,       |
|                                     |                    | options: lognormal, normal                                 |
+-------------------------------------+--------------------+------------------------------------------------------------+
| max_link_distance                   | real(r8)           | Maximum distance in meters along the stream for            |
|                                     |                    | localization, e.g., 10000 (10 km)                          |
+-------------------------------------+--------------------+------------------------------------------------------------+
| domain_order                        | character(256, 3)  | There are three possible domains to include: 'channel',    |
|                                     |                    | 'hru', 'parameters'. This variable specifies the order.    |
+-------------------------------------+--------------------+------------------------------------------------------------+
| domain_shapefiles                   | character(256, 3)  | Input (restart) files used to determine the shape of the   |
|                                     |                    | input variables and other geographic metadata. They must   |
|                                     |                    | be ordered as listed in ``domain_order``.                  |
+-------------------------------------+--------------------+------------------------------------------------------------+
| channel_config_file                 | character(256)     | Channel parameters are read from this file. These include  |
|                                     |                    | longitude, latitude, upstream and downstream segments, ..  |
+-------------------------------------+--------------------+------------------------------------------------------------+
| hru_config_file                     | character(256)     | HRU parameters are read from this file. These include      |
|                                     |                    | HRU ID, HRU-segment relations, elevation, ...              |
+-------------------------------------+--------------------+------------------------------------------------------------+
| hydro_config_file                   | character(256)     | The relation between HRUs and segments is read from this   |
|                                     |                    | file. This is used to construct groundwater masking.       |
+-------------------------------------+--------------------+------------------------------------------------------------+
| channel_variables                   | character(:, 5)    | The list of variables in the channel model such as         |
|                                     |                    | streamflow and their corresponding DART QTY and clamps.    |
+-------------------------------------+--------------------+------------------------------------------------------------+
| hru_variables                       | character(:, 5)    | The list of variables in the HRU such as groundwater.      |
+-------------------------------------+--------------------+------------------------------------------------------------+
| debug                               | integer            | The switch to specify the run-time verbosity.              |
|                                     |                    |                                                            |
|                                     |                    | - ``0`` is as quiet as it gets                             |
|                                     |                    | - ``>0`` provides more detailed run-time messages.         |
+-------------------------------------+--------------------+------------------------------------------------------------+

On top of the ``model_mod`` code, the other program in this directory is ``create_identity_streamflow_obs`` which is
based on the already available USGS observation converter. This converter provides a capability of evaluating identity
obs by adding in a desired list of gauges.    
