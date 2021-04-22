MODULE ``obs_def_gps_mod``
==========================

Overview
--------

DART GPS Radio Occultation observation module, including the observation operators 
for both local and non-local refractivity computations.

Author information:

-  Dr. Hui Liu

Namelist
--------

This namelist is now enabled by default. The maximum number of GPS observations is settable at runtime by changing the
value in the namelist. If you get an error about a missing namelist add ``&obs_def_gps_nml`` using the example below to
your ``input.nml`` namelist file and rerun. No recompiling is needed.

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_def_gps_nml
     max_gpsro_obs = 100000,
   / 

| 

.. container::

   +---------------+---------+------------------------------------------------------------------------------------------+
   | Item          | Type    | Description                                                                              |
   +===============+=========+==========================================================================================+
   | max_gpsro_obs | integer | The maximum number of GPS refractivity observations supported for a single execution.    |
   |               |         | Generally the default will be sufficient for a single run of filter, but not enough for  |
   |               |         | a long diagnostics run to produce a time series.                                         |
   +---------------+---------+------------------------------------------------------------------------------------------+

| 

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod (threed_sphere)
   assim_model_mod
   obs_kind_mod

Public interfaces
-----------------

============================= ======================
*use obs_def_gps_mod, only :* read_gpsro_ref
\                             write_gpsro_ref
\                             get_expected_gpsro_ref
\                             interactive_gpsro_ref
\                             set_gpsro_ref
\                             get_gpsro_ref
============================= ======================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *call read_gpsro_ref(gpskey, ifile, [, fform])*
   ::

      integer,          intent(out)          :: gpskey
      integer,          intent(in)           :: ifile
      character(len=*), intent(in), optional :: fform

.. container:: indent1

   Refractivity observations have several items of auxiliary data to read or write. This routine reads in the data for
   the next observation and returns the private GPS key index number that identifies the auxiliary data for this
   observation.

   ========== ====================================================================================================
   ``gpskey`` GPS key number returned to the caller.
   ``ifile``  Open file unit number to read from.
   *fform*    If specified, indicate whether the file was opened formatted or unformatted. Default is 'formatted'.
   ========== ====================================================================================================

| 

.. container:: routine

   *call write_gpsro_ref(gpskey, ifile, [, fform])*
   ::

      integer,          intent(in)           :: gpskey
      integer,          intent(in)           :: ifile
      character(len=*), intent(in), optional :: fform

.. container:: indent1

   Refractivity observations have several items of auxiliary data to read or write. This routine writes out the
   auxiliary data for the specified observation to the file unit given.

   ========== ====================================================================================================
   ``gpskey`` GPS key number identifying which observation to write aux data for.
   ``ifile``  Open file unit number to write to.
   *fform*    If specified, indicate whether the file was opened formatted or unformatted. Default is 'formatted'.
   ========== ====================================================================================================

| 

.. container:: routine

   *call get_expected_gpsro_ref(state_vector, location, gpskey, ro_ref, istatus)*
   ::

      real(r8),            intent(in)  :: state_vector(:)
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: gpskey
      real(r8),            intent(out) :: ro_ref
      integer,             intent(out) :: istatus

.. container:: indent1

   | Given a location and the state vector from one of the ensemble members, compute the model-predicted GPS
     refractivity that would be observed at that location. There are two types of operators: modeled *local*
     refractivity (N-1)*1.0e6 or *non_local* refractivity (excess phase, m) The type is indicated in the auxiliary
     information for each observation.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_vector`` | A one dimensional representation of the model state vector                                       |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``location``     | Location of this observation                                                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``gpskey``       | Integer key identifying which GPS observation this is, so the correct corresponding auxiliary    |
   |                  | information can be accessed.                                                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``ro_ref``       | The returned GPS refractivity value                                                              |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``istatus``      | Returned integer status code describing problems with applying forward operator. 0 is a good     |
   |                  | value; any positive value indicates an error; negative values are reserved for internal DART use |
   |                  | only.                                                                                            |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call interactive_gpsro_ref(gpskey)*
   ::

      integer, intent(out) :: gpskey

.. container:: indent1

   Prompts the user for the auxiliary information needed for a GPS refractivity observation, and returns the new key
   associated with this data.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``gpskey`` | Unique identifier associated with this GPS refractivity observation. In this code it is an integer     |
   |            | index into module local arrays which hold the additional data. This routine returns the incremented    |
   |            | value associated with this data.                                                                       |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call set_gpsro_ref(gpskey, nx, ny, nz, rfict0, ds, htop, subset0)*
   ::

      integer,          intent(out) :: gpskey
      real(r8),         intent(in)  :: nx
      real(r8),         intent(in)  :: ny
      real(r8),         intent(in)  :: nz
      real(r8),         intent(in)  :: rfict0
      real(r8),         intent(in)  :: ds
      real(r8),         intent(in)  :: htop
      character(len=6), intent(in)  :: subset0

.. container:: indent1

   Sets the auxiliary information associated with a GPS refractivity observation. This routine increments and returns
   the new key associated with these values.

   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``gpskey``  | Unique identifier associated with this GPS refractivity observation. In this code it is an integer    |
   |             | index into module local arrays which hold the additional data. This routine returns the incremented   |
   |             | value associated with this data.                                                                      |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``nx``      | X component of direction of ray between the LEO (detector) satellite and the GPS transmitter          |
   |             | satellite at the tangent point.                                                                       |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``ny``      | Y component of tangent ray.                                                                           |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``nz``      | Z component of tangent ray.                                                                           |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``rfict0``  | Local curvature radius (meters).                                                                      |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``ds``      | Delta S, increment to move along the ray in each direction when integrating the non-local operator    |
   |             | (meters).                                                                                             |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``htop``    | Elevation (in meters) where integration stops along the ray.                                          |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``subset0`` | The string 'GPSREF' for the local operator (refractivity computed only at the tangent point), or      |
   |             | 'GPSEXC' for the non-local operator which computes excess phase along the ray.                        |
   +-------------+-------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_gpsro_ref(gpskey, nx, ny, nz, rfict0, ds, htop, subset0)*
   ::

      integer,          intent(in)  :: gpskey
      real(r8),         intent(out) :: nx
      real(r8),         intent(out) :: ny
      real(r8),         intent(out) :: nz
      real(r8),         intent(out) :: rfict0
      real(r8),         intent(out) :: ds
      real(r8),         intent(out) :: htop
      character(len=6), intent(out) :: subset0

.. container:: indent1

   Gets the auxiliary information associated with a GPS refractivity observation, based on the GPS key number specified.

   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``gpskey``  | Unique identifier associated with this GPS refractivity observation. In this code it is an integer    |
   |             | index into module local arrays which hold the additional data. The value specified selects which      |
   |             | observation to return data for.                                                                       |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``nx``      | X component of direction of ray between the LEO (detector) satellite and the GPS transmitter          |
   |             | satellite at the tangent point.                                                                       |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``ny``      | Y component of tangent ray.                                                                           |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``nz``      | Z component of tangent ray.                                                                           |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``rfict0``  | Local curvature radius (meters).                                                                      |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``ds``      | Delta S, increment to move along the ray in each direction when integrating the non-local operator    |
   |             | (meters).                                                                                             |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``htop``    | Elevation (in meters) where integration stops along the ray.                                          |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``subset0`` | The string 'GPSREF' for the local operator (refractivity computed only at the tangent point), or      |
   |             | 'GPSEXC' for the non-local operator which computes excess phase along the ray.                        |
   +-------------+-------------------------------------------------------------------------------------------------------+

| 

Files
-----

-  A DART observation sequence file containing GPS obs.

References
----------

-  Assimilation of GPS Radio Occultation Data for Numerical Weather Prediction, Kuo,Y.H., Sokolovskiy,S.V., Anthes,R.A.,
   Vendenberghe,F., Terrestrial Atm and Ocn Sciences, Vol 11, pp157-186, 2000.


Error codes and conditions
--------------------------


+------------------------+---------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|         Routine        |                                     Message                                     |                                                              Comment                                                             |
+========================+=================================================================================+==================================================================================================================================+
| initialize_module      | initial allocation failed for gps observation data, itemcount = (max_gpsro_obs) | Need to increase max_gpsro_obs count in namelist                                                                                 |
+------------------------+---------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| gpskey_out_of_range    | gpskey (key#) exceeds max_radial_gps_obs (maxval)                               | The number of GPS observations exceeds the array size allocated in the module. Need to increase max_gpsro_obs count in namelist. |
+------------------------+---------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| read_gpsro_ref         | Expected header 'gpsroref' in input file                                        | The format of the input obs_seq file is not consistent.                                                                          |
+------------------------+---------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| get_expected_gpsro_ref | vertical location must be height; gps obs key #                                 | GPS observations must have vertical coordinates of height                                                                        |
+------------------------+---------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+

Future Plans
------------

- The current code first bins the very densely-sampled vertical profile into 200 bins, and then interpolates 
  the requested vertical location from that. The original profiles have been plotted and are smooth; 
  there appears to be no need to pre-bin the ata.

- The local operator needs no additional auxiliary data. The observation files would be much smaller if the
  local operator observation was a separate type without aux data, and only the non-local operator observation
  types would need the ray direction, the curvature, etc.
