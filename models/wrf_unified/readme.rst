.. index:: wrf_chem, WRF-Chem, WRF CHEM

.. _wrf_unified: 

WRF-Chem
=================

DART interface module for the Weather Research and Forecasting
`(WRF) <https://www.mmm.ucar.edu/weather-research-and-forecasting-model>`__
model including the WRF-Chem extension.

The model interface code supports WRF configurations with multiple domains. Data
for all domains is read into the DART state vector. During the computation of
the forward operators (getting the estimated observation values from each
ensemble member), the search starts in the domain with the highest number, which
is generally the finest nest or one of multiple finer nests. The search stops as
soon as a domain contains the observation location, working its way from largest
number to smallest number domain, ending with domain 1. For example, in a 4
domain case the data in the state vector that came from ``wrfinput_d04`` is
searched first, then ``wrfinput_d03``, ``wrfinput_d02``, and finally 
``wrfinput_d01``.

The forward operator is computed from the first (highest resolution) domain  that contains the
lat/lon of the observation. During the assimilation phase, when the state values
are adjusted based on the correlations and assimilation increments, all points
in all domains that are within the localization radius are adjusted, regardless
of domain. 

The fields from WRF that are copied into the DART state vector are controlled by
namelist. See below for the documentation on the &model_nml entries. The state
vector should include all fields needed to restart a WRF run. There may be
additional fields needed depending on the microphysics scheme selected. 

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: text

   &model_nml
      wrf_state_variables         = 'NULL'
      wrf_state_bounds            = 'NULL'
      num_domains                 = 1
      calendar_type               = 3
      assimilation_period_seconds = 21600
      allow_obs_below_vol         = .false.
      vert_localization_coord     = 3
      sfc_elev_max_diff           = -1.0
      polar                       = .false.
      periodic_x                  = .false.
      periodic_y                  = .false.
      allow_perturbed_ics         = .false.   # testing purposes only
   /


Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+---------------------------------------+-------------------+---------------------------------------+
| Item                                  | Type              | Description                           |
+=======================================+===================+=======================================+
| wrf_state_variables                   | character(:,4)    | A 2D array of strings, 4 per wrf      |
|                                       |                   | field to be added to the dart state   |
|                                       |                   | vector. The 4 strings are:            |
|                                       |                   |                                       |
|                                       |                   | #. WRF field name - must match netcdf |
|                                       |                   |    name exactly                       |
|                                       |                   | #. DART KIND name - must match a      |
|                                       |                   |    valid DART QTY_xxx exactly         |
|                                       |                   | #. 'UPDATE' or 'NO_COPY_BACK'         |
|                                       |                   |    If 'UPDATE', the data is written   |
|                                       |                   |    to netcdf file after               |
|                                       |                   |    the assimilation. If 'NO_COPY_BACK'|
|                                       |                   |    the data is not copied back to the |
|                                       |                   |    wrf netcdf file after the          |
|                                       |                   |    assimilation.                      |
|                                       |                   | #. A numeric string listing the       |
|                                       |                   |    domain numbers this array is part  |
|                                       |                   |    of. The special string 999 means   |
|                                       |                   |    all domains. For example, '12'     |
|                                       |                   |    means domains 1 and 2, '13' means  |
|                                       |                   |    1 and 3.                           |
+---------------------------------------+-------------------+---------------------------------------+
| wrf_state_bounds                      | character(:,3)    | A 2D array of strings, 3 per wrf      |
|                                       |                   | array. During the writing of data to  |
|                                       |                   | from the wrf netcdf file, variables   |
|                                       |                   | listed here will have minimum and     |
|                                       |                   | maximum values enforced. The 4        |
|                                       |                   | strings are:                          |
|                                       |                   |                                       |
|                                       |                   | #. WRF field name - must match netcdf |
|                                       |                   |    name exactly                       |
|                                       |                   | #. Minimum -- specified as a string   |
|                                       |                   |    but must be a numeric value (e.g.  |
|                                       |                   |    '0.1') Can be 'NULL' to allow any  |
|                                       |                   |    minimum value.                     |
|                                       |                   | #. Maximum -- specified as a string   |
|                                       |                   |    but must be a numeric value (e.g.  |
|                                       |                   |    '0.1') Can be 'NULL' to allow any  |
|                                       |                   |    maximum value.                     |
+---------------------------------------+-------------------+---------------------------------------+
| num_domains                           | integer           | Total number of WRF domains,          |
|                                       |                   | including nested domains.             |
+---------------------------------------+-------------------+---------------------------------------+
| calendar_type                         | integer           | Calendar type. Should be 3            |
|                                       |                   | (GREGORIAN) for WRF.                  |
+---------------------------------------+-------------------+---------------------------------------+
| assimilation_period_seconds           | integer           | The time (in seconds) between         |
|                                       |                   | assimilations. This is modified if    |
|                                       |                   | necessary to be an integer multiple   |
|                                       |                   | of the underlying model timestep.     |
+---------------------------------------+-------------------+---------------------------------------+
| periodic_x                            | logical           | If *.true.*, the grid is periodic in  |
|                                       |                   | longitude, and points above the last  |
|                                       |                   | grid cell and points below the first  |
|                                       |                   | grid cell are wrapped. Note this is   |
|                                       |                   | not the same as a grid which crosses  |
|                                       |                   | the prime meridian. WRF handles that  |
|                                       |                   | with an offset in longitude and       |
|                                       |                   | points beyond the last grid index are |
|                                       |                   | outside the domain.                   |
+---------------------------------------+-------------------+---------------------------------------+
| periodic_y                            | logical           | Used for the Single Column Model to   |
|                                       |                   | make the grid wrap in Y (see scm      |
|                                       |                   | below). This is NOT the same as       |
|                                       |                   | wrapping in latitude (see polar       |
|                                       |                   | below).                               |
+---------------------------------------+-------------------+---------------------------------------+
| polar                                 | logical           | If *.true.*, points at the poles are  |
|                                       |                   | wrapped across the grid. It is not    |
|                                       |                   | clear this is a good idea since the   |
|                                       |                   | grid is degenerate here.              |
+---------------------------------------+-------------------+---------------------------------------+
| sfc_elev_max_diff                     | real(r8)          | If > 0, the maximum difference, in    |
|                                       |                   | meters, between an observation marked |
|                                       |                   | as a 'surface obs' as the vertical    |
|                                       |                   | type (with the surface elevation, in  |
|                                       |                   | meters, as the numerical vertical     |
|                                       |                   | location), and the surface elevation  |
|                                       |                   | as defined by the model. Observations |
|                                       |                   | further away from the surface than    |
|                                       |                   | this threshold are rejected and not   |
|                                       |                   | assimilated. If the value is          |
|                                       |                   | negative, this test is skipped.       |
+---------------------------------------+-------------------+---------------------------------------+
| allow_obs_below_vol                   | logical           | If *.false.* then if an observation   |
|                                       |                   | with a vertical coordinate of         |
|                                       |                   | pressure or height (i.e. not a        |
|                                       |                   | surface observation) is below the     |
|                                       |                   | lowest 3d sigma level, it is outside  |
|                                       |                   | the field volume and the              |
|                                       |                   | interpolation routine rejects it. If  |
|                                       |                   | this is set to *.true.* and the       |
|                                       |                   | observation is above the surface      |
|                                       |                   | elevation but below the lowest field  |
|                                       |                   | volume level, the code will           |
|                                       |                   | extrapolate downward from data values |
|                                       |                   | at levels 1 and 2.                    |
+---------------------------------------+-------------------+---------------------------------------+
| vert_localization_coord               | integer           | Vertical coordinate for vertical      |
|                                       |                   | localization.                         |
|                                       |                   |                                       |
|                                       |                   | -  1 = model level                    |
|                                       |                   | -  2 = pressure (in pascals)          |
|                                       |                   | -  3 = height (in meters)             |
|                                       |                   | -  4 = scale height (unitless)        |
+---------------------------------------+-------------------+---------------------------------------+
| allow_perturbed_ics                   | logical           | *allow_perturbed_ics* should not be   |
|                                       |                   | used in most cases. It is provided    |
|                                       |                   | only as a means to create a tiny      |
|                                       |                   | ensemble for non-advancing tests.     |
+---------------------------------------+-------------------+---------------------------------------+

References
----------

https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/contents.html
