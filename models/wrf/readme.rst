WRF
===

Overview
--------

The following is a description of the DART interface module for the 
Weather Research and Forecasting model
`(WRF) <https://www.mmm.ucar.edu/weather-research-and-forecasting-model>`__
model. This page provides an overview of the module compiled into DART 
that interfaces with the WRF data in the state vector.
**The WRF-DART interface is compatible with WRF versions 4 and later, and is 
no longer backwards compatible with WRFv3.9 and earlier.**  
For more information on the interface changes required between 
different WRF versions, read through this documentation *and* the 
WRF-DART tutorial link in the next section.  

There have been several important updates to the WRF-DART interface starting
with `DARTv11.5.0. <https://github.com/NCAR/DART/releases/tag/v11.5.0>`__ 
Some important WRF-DART updates include:

- Version 11.4.1: Detects use of the Hybrid Vertical Coordinate system
  (terrain following at surface) and accounts for this in the forward
  operator calculations.

- Version 11.5.0: Improves compatibility with WRFv4+ versions where
  the prognostic 3D temperature variable is THM.  It is now mandatory to
  include THM instead of T in the ``wrf_state_variables`` namelist.



It is always recommended that you update your DART version to the 
`latest release <https://github.com/NCAR/DART/releases>`__ before beginning new research.

WRF-DART Tutorial
-----------------

This tutorial provides a real-world example of assimilating a wide variety of atmospheric
observations during an extreme storm event for the United States during April 2017.
**It is strongly recommended that you also review and perform the tutorial for 
running a WRF-DART assimilation** `here. <https://docs.dart.ucar.edu/en/latest/models/wrf/tutorial/README.html>`__


WRF Interface Overview
----------------------

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

.. Important::
   
   Although the model interface code is compatible with multiple domains, the 
   supporting `shell scripts <https://github.com/NCAR/DART/tree/main/models/wrf/shell_scripts>`__
   and WRF-DART tutorial are currently  designed for a single domain and will
   require some modifications for multiple (nested) domain functionality. If you
   need help with these modifications please contact DART support.


In summary, the forward operator is computed from the first domain grid (searching from
finest grid to coarsest grid) that contains the lat/lon of the observation. During the
assimilation phase, however,  when the state values are adjusted based on the correlations
and assimilation increments, all points in all domains that are within the 
localization radius are adjusted, regardless of domain. The impact of an observation 
on the state depends only on the distance between the observation and the state 
vector point, and the regression coefficient based on the correlation between the 
distributions of the ensemble of state vector points and the ensemble of observation 
forward operator values.

The fields from WRF that are copied into the DART state vector are controlled by
the namelist within ``input.nml``. See below for the documentation on the ``&model_nml`` entries within
``input.nml``. The state vector should include all fields needed to restart a WRF run.
There may be additional fields needed depending on the microphysics scheme selected. See the
ascii file `wrf_state_variables_table  <https://github.com/NCAR/DART/blob/main/models/wrf/wrf_state_variables_table>`__ 
that includes a list of fields that may be added to the DART state.

.. _wrfnamelist:

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

::

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


Namelist Description:
~~~~~~~~~~~~~~~~~~~~~

+-------------------------------+-------------------+---------------------------------------+
| Item                          | Type              | Description                           |
+===============================+===================+=======================================+
| default_state_variables       | logical           | If *.true.*, the dart state vector    |
|                               |                   | contains the fields U, V, W, PH, THM, |
|                               |                   | MU, in that order, and only those.    |
|                               |                   | Any values listed in the              |
|                               |                   | *wrf_state_variables* namelist item   |
|                               |                   | will be ignored.                      |
+-------------------------------+-------------------+---------------------------------------+
| wrf_state_variables           | character(:,5)    | A 2D array of strings, 5 per wrf      |
|                               |                   | array to be added to the dart state   |
|                               |                   | vector. If *default_state_variables*  |
|                               |                   | is *.true.*, this is ignored. When    |
|                               |                   | *.false.*, this list of array names   |
|                               |                   | controls which arrays and the order   |
|                               |                   | that they are added to the state      |
|                               |                   | vector. The 5 strings are:            |
|                               |                   |                                       |
|                               |                   | #. WRF field name - must match netcdf |
|                               |                   |    name exactly                       |
|                               |                   | #. DART Quantity name - must match a  |
|                               |                   |    valid DART QTY_xxx exactly         |
|                               |                   | #. WRF Type - supplements the quantity|
|                               |                   |    name to control the operation of   |
|                               |                   |    forward operator.                  |
|                               |                   | #. The string UPDATE or NO_COPY_BACK  |
|                               |                   |    Determines whether WRF state       |
|                               |                   |    is updated by DART                 |
|                               |                   | #. A numeric string listing the       |
|                               |                   |    domain(s) that include the WRF     |
|                               |                   |    state variable.                    |
|                               |                   |    The special string '999' means     |
|                               |                   |    all domains. For example, '12'     |
|                               |                   |    means domains 1 and 2, '13' means  |
|                               |                   |    1 and 3.                           |
+-------------------------------+-------------------+---------------------------------------+
| wrf_state_bounds              | character(:,4)    | A 2D array of strings, 4 per wrf      |
|                               |                   | array. During the copy of data to and |
|                               |                   | from the WRF (wrfinput*) file,        |
|                               |                   | variables listed here will have       |
|                               |                   | minimum and maximum values enforced.  |
|                               |                   | The 4 strings are:                    |
|                               |                   |                                       |
|                               |                   | #. WRF field name - must match        |
|                               |                   |    WRF variable name exactly          |
|                               |                   | #. Minimum -- specified as a string   |
|                               |                   |    but must be a numeric value (e.g.  |
|                               |                   |    '0.1') Can be 'NULL' to allow any  |
|                               |                   |    minimum value.                     |
|                               |                   | #. Maximum -- specified as a string   |
|                               |                   |    but must be a numeric value (e.g.  |
|                               |                   |    '0.1') Can be 'NULL' to allow any  |
|                               |                   |    maximum value.                     |
|                               |                   | #. Action -- valid strings are        |
|                               |                   |    'CLAMP' or 'FAIL'. Ignored by      |
|                               |                   |    filter. Filter will always clamp   |
|                               |                   |    if min and/or max is set.          |
+-------------------------------+-------------------+---------------------------------------+
| num_domains                   | integer           | Total number of WRF domains,          |
|                               |                   | including nested domains.             |
+-------------------------------+-------------------+---------------------------------------+
| calendar_type                 | integer           | Calendar type. Should be 3            |
|                               |                   | (GREGORIAN) for WRF.                  |
+-------------------------------+-------------------+---------------------------------------+
| assimilation_period_seconds   | integer           | The time (in seconds) between         |
|                               |                   | assimilations. This is modified if    |
|                               |                   | necessary to be an integer multiple   |
|                               |                   | of the underlying model timestep.     |
+-------------------------------+-------------------+---------------------------------------+
| periodic_x                    | logical           | If *.true.*, the grid is periodic in  |
|                               |                   | longitude, and points above the last  |
|                               |                   | grid cell and points below the first  |
|                               |                   | grid cell are wrapped. Note this is   |
|                               |                   | not the same as a grid which crosses  |
|                               |                   | the prime meridian. WRF handles that  |
|                               |                   | with an offset in longitude and       |
|                               |                   | points beyond the last grid index are |
|                               |                   | outside the domain.                   |
+-------------------------------+-------------------+---------------------------------------+
| periodic_y                    | logical           | Used for the WRF single column model  |
|                               |                   | to make the grid wrap in Y (see scm   |
|                               |                   | below). This is NOT the same as       |
|                               |                   | wrapping in latitude (see polar       |
|                               |                   | below).                               |
+-------------------------------+-------------------+---------------------------------------+
| polar                         | logical           | If *.true.*, points at the poles are  |
|                               |                   | wrapped across the grid. It is not    |
|                               |                   | clear this is a good idea because the |
|                               |                   | grid is degnerate here.               |
+-------------------------------+-------------------+---------------------------------------+
| scm                           | logical           | If *.true.* the single column model   |
|                               |                   | is assumed. The grid is a single      |
|                               |                   | vertical column, and there are 9      |
|                               |                   | cells arranged in a 3x3 grid. See the |
|                               |                   | WRF documentation for more            |
|                               |                   | information on this configuration.    |
|                               |                   | *periodic_x* and *periodic_y* should  |
|                               |                   | also be *.true.* in this case.        |
+-------------------------------+-------------------+---------------------------------------+
| sfc_elev_max_diff             | real(r8)          | The maximum elevation difference      |
|                               |                   | (in meters) between a 'surface'       |
|                               |                   | observation and the land surface      |
|                               |                   | elevation defined in WRF.             |
|                               |                   | If the value is > 0, that value is    |
|                               |                   | the threshold at which the surface    |
|                               |                   | observations are rejected. If the     |
|                               |                   | value is negative the test is skipped.|
+-------------------------------+-------------------+---------------------------------------+
| allow_obs_below_vol           | logical           | If *.false.* then if an observation   |
|                               |                   | with a vertical coordinate of         |
|                               |                   | pressure or height (i.e. not a        |
|                               |                   | surface observation) is below the     |
|                               |                   | lowest 3d sigma level, it is outside  |
|                               |                   | the field volume and the              |
|                               |                   | interpolation routine rejects it. If  |
|                               |                   | this is set to *.true.* and the       |
|                               |                   | observation is above the surface      |
|                               |                   | elevation but below the lowest field  |
|                               |                   | volume level, the code will           |
|                               |                   | extrapolate downward from data values |
|                               |                   | at levels 1 and 2.                    |
+-------------------------------+-------------------+---------------------------------------+
| center_search_half_length     | real(r8)          | A parameter in the 'use_old_vortex'   | 
|                               |                   | scheme used to search for a vortex    |
|                               |                   | center location. It is the half-length|   
|                               |                   | (meters) of a square box used during  |
|                               |                   | the vortex search. This value and the |
|                               |                   | 'center_spline_grid_scale' namelist   |
|                               |                   | items are required. To implement, set |
|                               |                   | ``use_old_vortex = .true.`` in        |
|                               |                   | ``model_mod.f90`` prior to compiling  |
|                               |                   | DART.                                 |
+-------------------------------+-------------------+---------------------------------------+
| center_spline_grid_scale      | integer           | A parameter in the 'use_old_vortex'   |
|                               |                   | scheme used to search for a vortex    |
|                               |                   | center location. It is the fine grid  |
|                               |                   | ratio for the spline interpolation    |
|                               |                   | used during the vortex search. This   |
|                               |                   | value and the                         | 
|                               |                   | 'center_search_half_length' namelist  |
|                               |                   | items are required. To implement, set |
|                               |                   | ``use_old_vortex = .true.`` in        |
|                               |                   | ``model_mod.f90`` prior to compiling  |
|                               |                   | DART.                                 |
+-------------------------------+-------------------+---------------------------------------+
| circulation_pres_level        | real(r8)          | A parameter in the 'circulation'      |
|                               |                   | scheme used to search for a vortex    |
|                               |                   | center location. It is the pressure   |
|                               |                   | (Pascals) at which the circulation is |
|                               |                   | computed during the vortex search.    |
|                               |                   | This value and the                    |
|                               |                   | 'circulation_radius' namelist items   |
|                               |                   | are required. To implement, set       |
|                               |                   | ``use_old_vortex = .false.`` in       |
|                               |                   | ``model_mod.f90`` prior to compiling  |
|                               |                   | DART.                                 |
+-------------------------------+-------------------+---------------------------------------+
| circulation_radius            | real(r8)          | A parameter in the 'circulation'      |
|                               |                   | scheme used to search for a vortex    |
|                               |                   | center location. It is the radius     |
|                               |                   | (meters) of the circle over which the |
|                               |                   | search for the vortex center is       |
|                               |                   | performed. This value and the         |
|                               |                   | 'circulation_pres_level' namelist     |
|                               |                   | items are required.  To implement,    |
|                               |                   | set ``use_old_vortex = .false.`` in   |
|                               |                   | ``model_mod.f90`` prior to compiling  |
|                               |                   | DART.                                 |
+-------------------------------+-------------------+---------------------------------------+
| vert_localization_coord       | integer           | Vertical coordinate for vertical      |
|                               |                   | localization.                         |
|                               |                   |                                       |
|                               |                   | -  1 = model level                    |
|                               |                   | -  2 = pressure (in pascals)          |
|                               |                   | -  3 = height (in meters)             |
|                               |                   | -  4 = scale height (unitless)        |
+-------------------------------+-------------------+---------------------------------------+
| allow_perturbed_ics           | logical           | *allow_perturbed_ics* should not be   |
|                               |                   | used in most cases. It is provided    |
|                               |                   | only as a means to create a tiny      |
|                               |                   | ensemble for non-advancing tests.     |
|                               |                   | Creating an initial ensemble is       |
|                               |                   | covered in :doc:`./tutorial/README`   |
+-------------------------------+-------------------+---------------------------------------+


Additional Namelist Information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- default_state_variables

You must set ``default_state_variables = .false.`` before changing the value
of ``wrf_state_variables`` to have it take effect.


- wrf_state_variables  

The format for ``wrf_state_variables`` is an array of 5 strings:
WRF output field, DART Quantity, WRF TYPE, 'UPDATE' or 'NO_COPY_BACK', and a numerical
string 'XXX'. If XXX=999 the variable is part of all domains, otherwise it is limited
to specific domains (e.g. '12' for domains 1 and 2, '13' for domains 1 and 3).
For example:

::

       wrf_state_variables='U','QTY_U_WIND_COMPONENT','TYPE_U','UPDATE','999',
                           'V','QTY_V_WIND_COMPONENT','TYPE_V','UPDATE','999',
                           'W','QTY_VERTICAL_VELOCITY','TYPE_W','UPDATE','999',
                           'THM','QTY_POTENTIAL_TEMPERATURE','TYPE_T','UPDATE','999',
                           'PH','QTY_GEOPOTENTIAL_HEIGHT','TYPE_GZ','UPDATE','999',
                           'MU','QTY_PRESSURE','TYPE_MU','UPDATE','999',
                           'QVAPOR','QTY_VAPOR_MIXING_RATIO','TYPE_QV','UPDATE','999',
                           'QCLOUD','QTY_CLOUD_LIQUID_WATER','TYPE_QC','UPDATE','999',
                           'QRAIN','QTY_RAINWATER_MIXING_RATIO','TYPE_QR','UPDATE','999',
                           'U10','QTY_U_WIND_COMPONENT','TYPE_U10','UPDATE','999',
                           'V10','QTY_V_WIND_COMPONENT','TYPE_V10','UPDATE','999',
                           'T2','QTY_TEMPERATURE','TYPE_T2','UPDATE','999',
                           'TH2','QTY_POTENTIAL_TEMPERATURE','TYPE_TH2','UPDATE','999',
                           'Q2','QTY_SPECIFIC_HUMIDITY','TYPE_Q2','UPDATE','999',
                           'PSFC','QTY_PRESSURE','TYPE_PS','UPDATE','999',


.. Important::

   It is mandatory to include ``THM`` instead of ``T`` as the ``TYPE_T`` WRF 
   temperature variable. This is because ``THM`` is the prognostic temperature variable
   that will impact the forecast when updated. 


- polar, periodic_x

The ``Polar`` and ``periodic_x`` namelist values are used in global WRF simulations.
If ``polar`` is true, the grid interpolation routines will wrap over the north and south poles.
If ``periodic_x`` is true, when the east and west edges of the grid are
reached the interpolation will wrap.  Note this is a separate issue
from regional models which cross the GMT line. Those grids are marked
as having a negative offset and do not need to wrap. This flag controls
what happens when the edges of the grid are reached.


- Single Column Model (scm)

The ``scm`` flag is used for the single column model version of WRF.
It needs the periodic_x and periodic_y flags set to true, in which
case the X and Y directions are periodic. There is no collapsing of the grid
into a single location like the 3d-spherical polar flag implies.

    
- sfc_elev_max_diff

The intent of the ``sfc_elev_max_diff`` quality control check is to eliminate
surface observations that are mismatched from the WRF model's surface elevation.
Mismatch can occur if the WRF land surface elevation is not finely resolved (coarse grid)
thus there is a significant representation mismatch between a point observation
and the WRF model. Assimilating surface observations with large mismatch can
deprecate assimilation forecast skill.
This check can only be applied to **surface observations** which are automatically
assigned to observations that use the ``VERTISSURFACE`` vertical coordinate
defined in the ``obs_seq.out`` file.   


- allow_obs_below_vol  

The ``allow_obs_below_vol`` enables vertical extrapolation in cases where the 
observation vertical location is below the lowest WRF model vertical layer, thus
used as an alternative for the standard vertical interpolation routine. 
The bottom WRF layer can vary based on total vertical levels, however, in general,
descends to (roughly) 10-50 meters above the surface and does not encompass common 
surface observations at 2 and 10 meters. This is not recommended given
(linear) extrapolation is a poor approximation of surface observations  at the
land-atmosphere boundary where energy and vapor exchange are controlled by 
similarity theory. When using  surface observations it is preferred
(and the default of the WRF ``model_mod.f90``) to operate on the WRF 2D 
surface output (e.g. T2, U10) instead of WRF 3D output (e.g. T, THM) to 
avoid the need for extrapolation.


- Vortex option

The vortex searching namelist options are only required during WRF simulations
where the spatial domain of interest is dynamic such as with a hurricane.




References
----------

https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/contents.html
