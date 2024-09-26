WRF
===

Overview
--------


DART interface module for the Weather Research and Forecasting
`(WRF) <https://www.mmm.ucar.edu/weather-research-and-forecasting-model>`__
model. This page documents the details of the
module compiled into DART that interfaces with the WRF data in the state vector.
**The WRF-DART interface is compatible with WRF versions 4 and later, and is 
no longer backwards compatible with WRFv3.9 and earlier.**  
For more information on the interface changes required between 
different WRF versions see the WRF tutorial link in the next section.

WRF+DART Tutorial
-----------------

**There is additional overview and tutorial documentation for running a WRF/DART
assimilation in** :doc:`./tutorial/README`

Please work through the tutorial in order to learn how to run WRF and DART.

Items of Note
~~~~~~~~~~~~~

- The ``model_mod`` reads WRF netCDF files directly to acquire the model state
  data. The ``wrf_to_dart`` and ``dart_to_wrf`` programs are no longer
  necessary.
- A netCDF file named ``wrfinput_d01`` is required and must be at the same
  resolution and have the same surface elevation data as the files converted to
  create the DART initial conditions. No data will be read from this file, but
  the grid information must match exactly.

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

The forward operator is computed from the first domain grid that contains the
lat/lon of the observation. During the assimilation phase, when the state values
are adjusted based on the correlations and assimilation increments, all points
in all domains that are within the localization radius are adjusted, regardless
of domain. The impact of an observation on the state depends only on the
distance between the observation and the state vector point, and the regression
coefficient based on the correlation between the distributions of the ensemble
of state vector points and the ensemble of observation forward operator values.

The fields from WRF that are copied into the DART state vector are controlled by
namelist. See below for the documentation on the &model_nml entries. The state
vector should include all fields needed to restart a WRF run. There may be
additional fields needed depending on the microphysics scheme selected. See the
ascii file ``wrf_state_variables_table`` in the ``models/wrf`` directory for a
list of fields that are often included in the DART state.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block::

   &model_nml
      default_state_variables     = .true.
      wrf_state_variables         = 'NULL'
      wrf_state_bounds            = 'NULL'
      num_domains                 = 1
      calendar_type               = 3
      assimilation_period_seconds = 21600
      allow_obs_below_vol         = .false.
      vert_localization_coord     = 3
      center_search_half_length   = 500000.
      center_spline_grid_scale    = 10
      circulation_pres_level      = 80000.0
      circulation_radius          = 108000.0
      sfc_elev_max_diff           = -1.0
      polar                       = .false.
      periodic_x                  = .false.
      periodic_y                  = .false.
      scm                         = .false.  
      allow_perturbed_ics         = .false.   # testing purposes only
   /

      # Notes for model_nml:
      # (1) vert_localization_coord must be one of:
      #     1 = model level
      #     2 = pressure
      #     3 = height
      #     4 = scale height
      # (2) see bottom of this file for explanations of polar, periodic_x, 
      #     periodic_y, and scm
      # (3) calendar = 3 is GREGORIAN, which is what WRF uses.
      # (4) if 'default_state_variables' is .true. the model_mod.f90 code will
      #     fill the state variable table with the following wrf vars: 
      #        U, V, W, PH, T, MU
      #     you must set it to false before you change the value 
      #     of 'wrf_state_variables' and have it take effect.
      # (5) the format for 'wrf_state_variables' is an array of 5 strings:
      #     wrf netcdf variable name, dart QTY_xxx string, type string (must be 
      #     unique, will soon be obsolete, we hope), 'UPDATE', and '999' if the 
      #     array is part of all domains.  otherwise, it is a string with the domain
      #     numbers (e.g. '12' for domains 1 and 2, '13' for domains 1 and 3).
      #   example:
      # wrf_state_variables='U','QTY_U_WIND_COMPONENT','TYPE_U','UPDATE','999',
      #                     'V','QTY_V_WIND_COMPONENT','TYPE_V','UPDATE','999',
      #                     'W','QTY_VERTICAL_VELOCITY','TYPE_W','UPDATE','999',
      #                     'T','QTY_POTENTIAL_TEMPERATURE','TYPE_T','UPDATE','999',
      #                     'PH','QTY_GEOPOTENTIAL_HEIGHT','TYPE_GZ','UPDATE','999',
      #                     'MU','QTY_PRESSURE','TYPE_MU','UPDATE','999',
      #                     'QVAPOR','QTY_VAPOR_MIXING_RATIO','TYPE_QV','UPDATE','999',
      #                     'QCLOUD','QTY_CLOUD_LIQUID_WATER','TYPE_QC','UPDATE','999',
      #                     'QRAIN','QTY_RAINWATER_MIXING_RATIO','TYPE_QR','UPDATE','999',
      #                     'U10','QTY_U_WIND_COMPONENT','TYPE_U10','UPDATE','999',
      #                     'V10','QTY_V_WIND_COMPONENT','TYPE_V10','UPDATE','999',
      #                     'T2','QTY_TEMPERATURE','TYPE_T2','UPDATE','999',
      #                     'TH2','QTY_POTENTIAL_TEMPERATURE','TYPE_TH2','UPDATE','999',
      #                     'Q2','QTY_SPECIFIC_HUMIDITY','TYPE_Q2','UPDATE','999',
      #                     'PSFC','QTY_PRESSURE','TYPE_PS','UPDATE','999',
      # (6) the format for 'wrf_state_bounds' is an array of 4 strings:
      #     wrf netcdf variable name, minimum value, maximum value, and either
      #     FAIL or CLAMP.  FAIL will halt the program if an out of range value
      #     is detected.  CLAMP will set out of range values to the min or max.
      #     The special string 'NULL' will map to plus or minus infinity and will
      #     not change the values.  arrays not listed in this table will not
      #     be changed as they are read or written.
      #
      #
      # polar and periodic_x are used in global wrf.  if polar is true, the 
      # grid interpolation routines will wrap over the north and south poles.  
      # if periodic_x is true, when the east and west edges of the grid are
      # reached the interpolation will wrap.  note this is a separate issue
      # from regional models which cross the GMT line; those grids are marked
      # as having a negative offset and do not need to wrap; this flag controls
      # what happens when the edges of the grid are reached.

      # the scm flag is used for the 'single column model' version of WRF.
      # it needs the periodic_x and periodic_y flags set to true, in which
      # case the X and Y directions are periodic; no collapsing of the grid
      # into a single location like the 3d-spherical polar flag implies.

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+---------------------------------------+-------------------+---------------------------------------+
| Item                                  | Type              | Description                           |
+=======================================+===================+=======================================+
| default_state_variables               | logical           | If *.true.*, the dart state vector    |
|                                       |                   | contains the fields U, V, W, PH, T,   |
|                                       |                   | MU, in that order, and only those.    |
|                                       |                   | Any values listed in the              |
|                                       |                   | *wrf_state_variables* namelist item   |
|                                       |                   | will be ignored.                      |
+---------------------------------------+-------------------+---------------------------------------+
| wrf_state_variables                   | character(:, 5)   | A 2D array of strings, 5 per wrf      |
|                                       |                   | array to be added to the dart state   |
|                                       |                   | vector. If *default_state_variables*  |
|                                       |                   | is *.true.*, this is ignored. When    |
|                                       |                   | *.false.*, this list of array names   |
|                                       |                   | controls which arrays and the order   |
|                                       |                   | that they are added to the state      |
|                                       |                   | vector. The 5 strings are:            |
|                                       |                   |                                       |
|                                       |                   | #. WRF field name - must match netcdf |
|                                       |                   |    name exactly                       |
|                                       |                   | #. DART KIND name - must match a      |
|                                       |                   |    valid DART QTY_xxx exactly         |
|                                       |                   | #. TYPE_NN - will hopefully be        |
|                                       |                   |    obsolete, but for now NN should    |
|                                       |                   |    match the field name.              |
|                                       |                   | #. the string UPDATE. at some future  |
|                                       |                   |    point, non-updatable fields may    |
|                                       |                   |    become part of the state vector.   |
|                                       |                   | #. A numeric string listing the       |
|                                       |                   |    domain numbers this array is part  |
|                                       |                   |    of. The specical string 999 means  |
|                                       |                   |    all domains. For example, '12'     |
|                                       |                   |    means domains 1 and 2, '13' means  |
|                                       |                   |    1 and 3.                           |
+---------------------------------------+-------------------+---------------------------------------+
| wrf_state_bounds                      | character(:, 4)   | A 2D array of strings, 4 per wrf      |
|                                       |                   | array. During the copy of data to and |
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
|                                       |                   | #. Action -- valid strings are        |
|                                       |                   |    'CLAMP', 'FAIL'. 'FAIL' means if a |
|                                       |                   |    value is found outside the range,  |
|                                       |                   |    the code fails with an error.      |
|                                       |                   |    'CLAMP' simply sets the out of     |
|                                       |                   |    range values to the given minimum  |
|                                       |                   |    or maximum without error.          |
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
|                                       |                   | grid is degnerate here.               |
+---------------------------------------+-------------------+---------------------------------------+
| scm                                   | logical           | If *.true.* the Single Column Model   |
|                                       |                   | is assumed. The grid is a single      |
|                                       |                   | vertical column, and there are 9      |
|                                       |                   | cells arranged in a 3x3 grid. See the |
|                                       |                   | WRF documentation for more            |
|                                       |                   | information on this configuration.    |
|                                       |                   | *periodic_x* and *periodic_y* should  |
|                                       |                   | also be *.true.* in this case.        |
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
| center_search_half_length             | real(r8)          | The model_mod now contains two        |
|                                       |                   | schemes for searching for a vortex    |
|                                       |                   | center location. If the **old**       |
|                                       |                   | scheme is compiled in, then this and  |
|                                       |                   | the center_spline_grid_scale namelist |
|                                       |                   | items are used. (Search code for      |
|                                       |                   | 'use_old_vortex'.) Half length (in    |
|                                       |                   | meters) of a square box for searching |
|                                       |                   | the vortex center.                    |
+---------------------------------------+-------------------+---------------------------------------+
| center_spline_grid_scale              | integer           | The model_mod now contains two        |
|                                       |                   | schemes for searching for a vortex    |
|                                       |                   | center location. If the **old**       |
|                                       |                   | scheme is compiled in, then this and  |
|                                       |                   | the center_search_half_length         |
|                                       |                   | namelist items are used. (Search code |
|                                       |                   | for 'use_old_vortex'.) Ratio of       |
|                                       |                   | refining grid for                     |
|                                       |                   | spline-interpolation in determining   |
|                                       |                   | the vortex center.                    |
+---------------------------------------+-------------------+---------------------------------------+
| circulation_pres_level                | real(r8)          | The model_mod now contains two        |
|                                       |                   | schemes for searching for a vortex    |
|                                       |                   | center location. If the **new**       |
|                                       |                   | scheme is compiled in, then this and  |
|                                       |                   | the circulation_radius namelist items |
|                                       |                   | are used. (Search code for            |
|                                       |                   | 'use_old_vortex'.) Pressure, in       |
|                                       |                   | pascals, of the level at which the    |
|                                       |                   | circulation is computed when          |
|                                       |                   | searching for the vortex center.      |
+---------------------------------------+-------------------+---------------------------------------+
| circulation_radius                    | real(r8)          | The model_mod now contains two        |
|                                       |                   | schemes for searching for a vortex    |
|                                       |                   | center location. If the **new**       |
|                                       |                   | scheme is compiled in, then this and  |
|                                       |                   | the circulation_pres_level namelist   |
|                                       |                   | items are used. (Search code for      |
|                                       |                   | 'use_old_vortex'.) Radius, in meters, |
|                                       |                   | of the circle over which the          |
|                                       |                   | circulation calculation is done when  |
|                                       |                   | searching for the vortex center.      |
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
|                                       |                   | Creating an initial ensemble is       |
|                                       |                   | covered in :doc:`./tutorial/README`   |
+---------------------------------------+-------------------+---------------------------------------+


The following items used to be in the WRF namelist but have been removed. The
first 4 are no longer needed, and the last one was moved to the
``&dart_to_wrf_nml`` namelist in 2010. In the Lanai release having these values
in the namelist does not cause a fatal error, but more recent versions of the
code will fail if any of these values are specified. Remove them from your
namelist to avoid errors.

=================== ================= =========================================
Item                Type              Description
=================== ================= =========================================
``surf_obs``        logical           OBSOLETE -- now an error to specify this.
``soil_data``       logical           OBSOLETE -- now an error to specify this.
``h_diab``          logical           OBSOLETE -- now an error to specify this.
``num_moist_vars``  integer           OBSOLETE -- now an error to specify this.
``adv_mod_command`` character(len=32) OBSOLETE -- now an error to specify this.
=================== ================= =========================================

Files
-----

-  model_nml in input.nml
-  wrfinput_d01, wrfinput_d02, ... (one file for each domain)
-  netCDF output state diagnostics files

References
----------

https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/contents.html
