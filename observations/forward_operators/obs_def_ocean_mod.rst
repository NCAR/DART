MODULE obs_def_ocean_mod
========================

Overview
--------

| DART includes a flexible, powerful, and slightly complicated mechanism for incorporating new types of observations.
  The ``obs_def_ocean_mod`` module being described here is used by the program ``preprocess`` to insert appropriate
  definitions of ocean observations into the ``DEFAULT_obs_def_mod.f90`` template and generate the source files
  ``obs_def_mod.f90`` and ``obs_kind_mod.f90`` that are used by ``filter`` and other DART programs.

| Only ``HFRADAR_RADIAL_VELOCITY`` observations require a forward operator, as evidenced by the fact there is no
  ``COMMON_CODE`` in the third column of the type definitions table.
  All other observations types map to quantities that must be available in the model state;
  the observations types flagged with ``COMMON_CODE`` will use the ``model_interpolate()`` routine
  as the forward operator.

| The mandatory header line is followed by lines that have the observation type name (an all caps Fortran 90 identifier)
  and their associated generic quantity identifier from the obs_kind module. If there is no special processing needed
  for an observation type, and no additional data needed beyond the standard contents of an observation, then a third
  word on the line, the ``COMMON_CODE`` will instruct the preprocess program to automatically generate all stubs and
  code needed for this type. For observation types needing any special code or additional data, this word should not be
  specified and the user must supply the code manually. One of the future extensions of this module will be to support
  acoustic tomographic observations, which will necessitate specific support routines.

Ocean variable types and their corresponding quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::


   ! BEGIN DART PREPROCESS TYPE DEFINITIONS
   ! SALINITY,                      QTY_SALINITY,              COMMON_CODE
   ! TEMPERATURE,                   QTY_TEMPERATURE,           COMMON_CODE
   ! U_CURRENT_COMPONENT,           QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! V_CURRENT_COMPONENT,           QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! SEA_SURFACE_HEIGHT,            QTY_SEA_SURFACE_HEIGHT,    COMMON_CODE
   ! SEA_SURFACE_PRESSURE,          QTY_SEA_SURFACE_PRESSURE,  COMMON_CODE
   ! ARGO_U_CURRENT_COMPONENT,      QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! ARGO_V_CURRENT_COMPONENT,      QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! ARGO_SALINITY,                 QTY_SALINITY,              COMMON_CODE
   ! ARGO_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
   ! ADCP_U_CURRENT_COMPONENT,      QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! ADCP_V_CURRENT_COMPONENT,      QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! ADCP_SALINITY,                 QTY_SALINITY,              COMMON_CODE
   ! ADCP_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
   ! FLOAT_SALINITY,                QTY_SALINITY,              COMMON_CODE
   ! FLOAT_TEMPERATURE,             QTY_TEMPERATURE,           COMMON_CODE
   ! DRIFTER_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! DRIFTER_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! DRIFTER_SALINITY,              QTY_SALINITY,              COMMON_CODE
   ! DRIFTER_TEMPERATURE,           QTY_TEMPERATURE,           COMMON_CODE
   ! GLIDER_U_CURRENT_COMPONENT,    QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! GLIDER_V_CURRENT_COMPONENT,    QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! GLIDER_SALINITY,               QTY_SALINITY,              COMMON_CODE
   ! GLIDER_TEMPERATURE,            QTY_TEMPERATURE,           COMMON_CODE
   ! MOORING_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! MOORING_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! MOORING_SALINITY,              QTY_SALINITY,              COMMON_CODE
   ! MOORING_TEMPERATURE,           QTY_TEMPERATURE,           COMMON_CODE
   ! MOORING_PRESSURE,              QTY_PRESSURE,              COMMON_CODE
   ! BOTTLE_SALINITY,               QTY_SALINITY,              COMMON_CODE
   ! BOTTLE_TEMPERATURE,            QTY_TEMPERATURE,           COMMON_CODE
   ! CTD_SALINITY,                  QTY_SALINITY,              COMMON_CODE
   ! CTD_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
   ! TCTD_SALINITY,                 QTY_SALINITY,              COMMON_CODE
   ! TCTD_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
   ! STD_SALINITY,                  QTY_SALINITY,              COMMON_CODE
   ! STD_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
   ! XCTD_SALINITY,                 QTY_SALINITY,              COMMON_CODE
   ! XCTD_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
   ! MBT_SALINITY,                  QTY_SALINITY,              COMMON_CODE
   ! MBT_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
   ! XBT_SALINITY,                  QTY_SALINITY,              COMMON_CODE
   ! XBT_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
   ! DBT_SALINITY,                  QTY_SALINITY,              COMMON_CODE
   ! DBT_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
   ! APB_SALINITY,                  QTY_SALINITY,              COMMON_CODE
   ! APB_TEMPERATURE,               QTY_TEMPERATURE,           COMMON_CODE
   ! DOPPLER_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! DOPPLER_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! DOPPLER_W_CURRENT_COMPONENT,   QTY_W_CURRENT_COMPONENT,   COMMON_CODE
   ! SATELLITE_MICROWAVE_SST,       QTY_TEMPERATURE,           COMMON_CODE
   ! SATELLITE_INFRARED_SST,        QTY_TEMPERATURE,           COMMON_CODE
   ! SATELLITE_BLENDED_SST,         QTY_TEMPERATURE,           COMMON_CODE
   ! SATELLITE_SSH,                 QTY_SEA_SURFACE_HEIGHT,    COMMON_CODE
   ! SATELLITE_SSS,                 QTY_SALINITY,              COMMON_CODE
   ! J1_SEA_SURFACE_ANOMALY,        QTY_SEA_SURFACE_ANOMALY,   COMMON_CODE
   ! EN_SEA_SURFACE_ANOMALY,        QTY_SEA_SURFACE_ANOMALY,   COMMON_CODE
   ! GFO_SEA_SURFACE_ANOMALY,       QTY_SEA_SURFACE_ANOMALY,   COMMON_CODE
   ! DRY_LAND,                      QTY_DRY_LAND,              COMMON_CODE
   ! OI_SEA_SURFACE_TEMPERATURE,    QTY_TEMPERATURE,           COMMON_CODE
   ! HFRADAR_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
   ! HFRADAR_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
   ! HFRADAR_RADIAL_VELOCITY,       QTY_VELOCITY
   ! FERRYBOX_SALINITY,             QTY_SALINITY,              COMMON_CODE
   ! FERRYBOX_TEMPERATURE,          QTY_TEMPERATURE,           COMMON_CODE
   ! END DART PREPROCESS TYPE DEFINITIONS

New observation types may be added to this list with no loss of generality. Supporting the observations and actually
**assimilating** them are somewhat different and is controlled by the ``input.nml``\ ``&obs_kind_nml``
`assimilate_these_obs_types <../../assimilation_code/modules/observations/obs_kind_mod.html#Namelist>`__ variable. This
provides the flexibility to have an observation sequence file containing many different observation types and being able
to selectively choose what types will be assimilated.

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod (threed_sphere)
   assim_model_mod
   obs_kind_mod
   ensemble_manager_mod
   obs_def_utilities_mod


Public interfaces
-----------------

=============================== ==========================
*use obs_def_ocean_mod, only :* read_hf_radial_vel
\                               write_hf_radial_vel
\                               interactive_hf_radial_vel
\                               get_expected_hf_radial_vel
\                               get_obs_def_hf_radial_vel 
\                               set_hf_radial_vel
=============================== ==========================


| 

Namelist
--------

Namelist interface ``obs_def_ocean_nml`` is read from file ``input.nml``.
Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent 
them from prematurely terminating the namelist.

::

   &obs_def_ocean_nml
      max_radial_vel_obs = 1000000
      debug = .false.
      /

| 

.. container::

   +---------------------+------------+---------------------------------------------------------+
   | Item                | Type       | Description                                             |
   +=====================+============+=========================================================+
   | max_radial_vel_obs  | integer    | The maximum number of radial velocity observations      |
   |                     |            | to be read at one time. An error is thrown if more      |
   |                     |            | observations are encountered.                           |
   |                     |            | Increase value and rerun.                               |
   +---------------------+------------+---------------------------------------------------------+
   | debug               | logical    | Switch to control how much run-time output is created.  |
   |                     |            | ``.false.`` indicates less output,                      |
   |                     |            | ``.true.`` indicates more output.                       |
   +---------------------+------------+---------------------------------------------------------+

| 

Public components
-----------------

none

Files
-----

none

References
----------

none

Private components
------------------

N/A
