.. index:: wrf_unified, WRF, wrf_chem, WRF-Chem, WRF CHEM

.. _wrf_unified:

WRF WRF-Chem Unified Model Interface
====================================

DART interface module for the Weather Research and Forecasting
`(WRF) <https://www.mmm.ucar.edu/weather-research-and-forecasting-model>`__
model including the WRF-Chem extension.

The model interface code supports WRF configurations with multiple domains. Data
for all domains is read into the DART state vector. During the computation of
the forward operators (getting the estimated observation values from each
ensemble member), the search starts in the domain with the highest number, which
is generally the finest nest or one of multiple finer nests. The search stops as
soon as a domain contains the observation location, working its way from the
largest numbered domain to the smallest, ending with domain 1. For example, in a
4 domain case the data in the state vector that came from ``wrfinput_d04`` is
searched first, then ``wrfinput_d03``, ``wrfinput_d02``, and finally
``wrfinput_d01``. The forward operator is computed from the first (highest
resolution) domain that contains the lat/lon of the observation.

During the assimilation phase, when the state values are adjusted based on the
correlations and assimilation increments, all points in all domains that are
within the localization radius are adjusted, regardless of domain.

The fields from WRF that are read into the DART state vector are controlled by
namelist. See below for the documentation on the ``&model_nml`` entries. The
state vector should include all fields needed to restart a WRF run. There may be
additional fields needed depending on the microphysics scheme selected.

Input files
-----------

- Meteorological fields are read from ``wrfinput_d01`` ... ``wrfinput_d0N``
  (one file per domain, ``N`` = ``num_domains``, at most 9).
- Chemistry fields can either be in the same file as the meteorological fields,
  or in separate files. This is set by ``chemistry_separate_file``:

  - ``.false.`` (default): chemistry fields are in ``wrfinput_d0N``. List the
    chemistry variables and bounds in ``wrf_state_variables`` and
    ``wrf_state_bounds``, together with the meteorological variables.
    ``chem_state_variables`` and ``chem_state_bounds`` are not used.
  - ``.true.``: chemistry fields are read from ``wrfchem_d01`` ...
    ``wrfchem_d0N``. List the chemistry variables and bounds in
    ``chem_state_variables`` and ``chem_state_bounds``, and the meteorological
    variables in ``wrf_state_variables`` and ``wrf_state_bounds``. The
    ``wrfinput_d0N`` files are still needed for the grid and base state
    information.

- Variable names in the namelist must match the netCDF variable names exactly.

Filter input and output files with separate chemistry files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``chemistry_separate_file = .true.``, each WRF domain contributes two DART
state domains, so there are ``2 * num_domains`` in total. They are ordered with
all the meteorological domains first, then all the chemistry domains:

  #. ``wrfinput_d01`` ... ``wrfinput_d0N`` (meteorology)
  #. ``wrfchem_d01`` ... ``wrfchem_d0N`` (chemistry)

``filter`` needs a set of input and output files for each of these state
domains, in that order. For example, with 2 WRF domains:

.. code-block:: text

   &filter_nml
      input_state_file_list  = 'input_wrf_d01.txt',  'input_wrf_d02.txt',
                               'input_chem_d01.txt', 'input_chem_d02.txt'
      output_state_file_list = 'output_wrf_d01.txt',  'output_wrf_d02.txt',
                               'output_chem_d01.txt', 'output_chem_d02.txt'
   /

Each text file lists one file per ensemble member. Alternatively, for a single
member, use ``input_state_files`` and ``output_state_files`` with
``num_domains * 2`` file names, in the same order.

The model interface makes these assumptions:

- The chemistry file for domain ``n`` is on the same grid as ``wrfinput_d0n``.
  The grid, map projection and base state are read only from ``wrfinput_d0n``,
  so these files must be present in the run directory even if they are not in
  the state. The chemistry files are used only for the state variables.
- The meteorological and chemistry files are listed in the order above, with
  the same number of domains for each. A chemistry domain is needed for every
  WRF domain, even if you list no chemistry variables for it.
- A quantity should be listed in either the meteorological or the chemistry
  variables, not both. When a quantity is needed, the meteorological domain is
  searched first, then the chemistry domain.
- Chemistry variables use the same domain selection strings (e.g. ``'999'``,
  ``'12'``) as meteorological variables, but select which chemistry domains the
  variable is in.
- The same ``num_domains`` applies to both sets of files.

.. note::

   ``PHB`` (base state geopotential) should be included in the state vector,
   with the ``NO_COPY_BACK`` option. It is needed to compute model heights, and
   it is not changed by the assimilation, so it is not written back to the WRF
   netCDF file.

.. note::

   Some variables are assigned a DART quantity by the model interface
   regardless of the quantity given in the namelist: ``MU`` is ``QTY_PRESSURE``,
   ``PSFC`` is ``QTY_SURFACE_PRESSURE``, ``T2`` is ``QTY_2M_TEMPERATURE``,
   ``TH2`` is ``QTY_2M_POTENTIAL_TEMPERATURE`` and ``Q2`` is
   ``QTY_2M_SPECIFIC_HUMIDITY``. ``QTY_TEMPERATURE`` is converted to
   ``QTY_POTENTIAL_TEMPERATURE``. 

Hybrid vertical coordinate
~~~~~~~~~~~~~~~~~~~~~~~~~~

The interface supports the WRF hybrid vertical coordinate. It is detected from
the ``HYBRID_OPT`` global attribute of ``wrfinput_d0N`` (``HYBRID_OPT = 2``). If
the attribute is absent, the terrain following coordinate is assumed. 

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: text

   &model_nml
      wrf_state_variables  = 'U',     'QTY_U_WIND_COMPONENT',     'UPDATE',       '999',
                             'V',     'QTY_V_WIND_COMPONENT',     'UPDATE',       '999',
                             'W',     'QTY_VERTICAL_VELOCITY',    'UPDATE',       '999',
                             'PH',    'QTY_GEOPOTENTIAL_HEIGHT',  'UPDATE',       '999',
                             'T',     'QTY_POTENTIAL_TEMPERATURE','UPDATE',       '999',
                             'MU',    'QTY_PRESSURE',             'UPDATE',       '999',
                             'QVAPOR','QTY_VAPOR_MIXING_RATIO',   'UPDATE',       '999',
                             'PSFC',  'QTY_SURFACE_PRESSURE',     'UPDATE',       '999',
                             'PHB',   'QTY_BASE_STATE_GEOP',      'NO_COPY_BACK', '999'
      wrf_state_bounds     = 'QVAPOR','0.0','NULL',
                             'QRAIN', '0.0','NULL',
                             'QCLOUD','0.0','NULL'
      chemistry_separate_file = .true.
      chem_state_variables = 'o3', 'QTY_O3', 'UPDATE', '999',
                             'no', 'QTY_NO', 'UPDATE', '999'
      chem_state_bounds    = 'o3', '0.0', 'NULL'
      num_domains                 = 1
      calendar_type               = 3
      assimilation_period_seconds = 21600
      sfc_elev_max_diff           = -1.0
      vert_localization_coord     = 3
      allow_perturbed_ics         = .false.
      allow_obs_below_vol         = .false.
      log_vert_interp             = .true.
      log_horz_interpM            = .false.
      log_horz_interpQ            = .false.
      polar                       = .false.
      periodic_x                  = .false.
      periodic_y                  = .false.
   /


Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
    :header-rows: 1
    :widths: 25 20 55

    * - Item
      - Type and default
      - Description
    * - wrf_state_variables
      - character(:,4)

        default: 'NULL'
      - A list of strings, 4 per WRF field to be added to the DART state vector. The 4 strings are:

          #. WRF field name - must match the netCDF name exactly
          #. DART QTY name - must match a valid DART QTY_xxx exactly
          #. 'UPDATE' or 'NO_COPY_BACK'. If 'UPDATE', the data is written to the netCDF file after the assimilation. If 'NO_COPY_BACK', it is not.
          #. A numeric string listing the domain numbers this field is part of. The special string '999' means all domains. For example, '12' means domains 1 and 2, '13' means 1 and 3.

    * - wrf_state_bounds
      - character(:,3)

        default: 'NULL'
      - A list of strings, 3 per WRF field. When data is written to the WRF netCDF file, variables listed here have minimum and maximum values enforced. The 3 strings are:

          #. WRF field name - must match the netCDF name exactly
          #. Minimum - specified as a string but must be a numeric value (e.g. '0.1'). Can be 'NULL' to allow any minimum value.
          #. Maximum - specified as a string but must be a numeric value (e.g. '10.0'). Can be 'NULL' to allow any maximum value.

    * - chemistry_separate_file
      - logical

        default: .false.
      - If .false., chemistry fields are in the same netCDF file as the meteorological fields (``wrfinput_d0N``). ``wrf_state_variables`` and ``wrf_state_bounds`` then hold both the meteorological and the chemistry fields, and ``chem_state_variables`` and ``chem_state_bounds`` are ignored.
        If .true., chemistry fields are in separate files (``wrfchem_d0N``) and are listed in ``chem_state_variables`` and ``chem_state_bounds``; ``wrf_state_variables`` and ``wrf_state_bounds`` hold only the meteorological fields. The chemistry fields are separate DART domains, numbered after the ``num_domains`` meteorological domains.
    * - chem_state_variables
      - character(:,4)

        default: 'NULL'
      - Chemistry state variables, same format as ``wrf_state_variables``. Only used if ``chemistry_separate_file`` is .true.
    * - chem_state_bounds
      - character(:,3)

        default: 'NULL'
      - Chemistry state bounds, same format as ``wrf_state_bounds``. Only used if ``chemistry_separate_file`` is .true.
    * - num_domains
      - integer

        default: 1
      - Total number of WRF domains, including nested domains.
    * - calendar_type
      - integer

        default: 3
      - Calendar type. Should be 3 (GREGORIAN) for WRF.
    * - assimilation_period_seconds
      - integer

        default: 21600
      - The time (in seconds) between assimilations. This is modified if necessary to be an integer multiple of the underlying model timestep.
    * - sfc_elev_max_diff
      - real(r8)

        default: -1.0
      - If > 0, the maximum difference, in meters, between an observation marked as a 'surface obs' as the vertical type (with the surface elevation, in meters, as the numerical vertical location), and the surface elevation as defined by the model. Observations further away from the surface than this threshold are rejected and not assimilated. If the value is negative, this test is skipped.
    * - vert_localization_coord
      - integer

        default: 3
      - Vertical coordinate for vertical localization.

          -  1 = model level
          -  2 = pressure (in pascals)
          -  3 = height (in meters)
          -  4 = scale height (unitless)
    * - allow_perturbed_ics
      - logical

        default: .false.
      - Should not be used in most cases. Provided only for testing purposes to create a tiny ensemble for non-advancing tests.
    * - allow_obs_below_vol
      - logical

        default: .false.
      - If .false., pressure or height observations above the surface but below the lowest model level are rejected. If .true., the model values are extrapolated downward from the lowest levels so these observations can be used.
    * - log_vert_interp
      - logical

        default: .true.
      - If .true., vertical interpolation (and extrapolation) of pressure is done after taking the log of the pressure values. If .false., it is linear in pressure.
    * - log_horz_interpM
      - logical

        default: .false.
      - If .true., horizontal interpolation of pressure on the mass grid points is done after taking the log. If .false., it is linear in pressure.
    * - log_horz_interpQ
      - logical

        default: .false.
      - If .true., horizontal interpolation of the moisture (Q) fields is done after taking the log. If .false., it is linear.
    * - polar
      - logical

        default: .false.
      - Set to .true. if the WRF domain 1 is a polar (global) domain. Applies to domain 1 only.
    * - periodic_x
      - logical

        default: .false.
      - Set to .true. if the WRF domain 1 is periodic in the west-east direction (e.g. global domain). Applies to domain 1 only.
    * - periodic_y
      - logical

        default: .false.
      - Set to .true. if the WRF domain 1 is periodic in the south-north direction. Applies to domain 1 only.

References
----------

- `WRF user guide <https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/contents.html>`__
- `WRF-Chem <https://www2.acom.ucar.edu/wrf-chem>`__
