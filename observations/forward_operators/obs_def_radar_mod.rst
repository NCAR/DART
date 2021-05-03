MODULE ``obs_def_radar_mod``
============================

Overview
--------

| DART radar observation module, including the observation operators for the two primary radar-observation types --
  Doppler velocity and reflectivity -- plus other utility subroutines and functions. A number of simplifications are
  employed for the observation operators. Most notably, the model state is mapped to a "point" observation, whereas a
  real radar observation is a volumetric sample. The implications of this approximation have not been investigated
  fully, so in the future it might be worth developing and testing more sophisticated observation operators that produce
  volumetric power- weighted samples.
| This module is able to compute reflectivity and precipitation fall speed (needed for computing Doppler radial
  velocity) from the prognostic model fields only for simple single-moment microphysics schemes such as the Kessler and
  Lin schemes. If a more complicated microphysics scheme is used, then reflectivity and fall speed must be accessible
  instead as diagnostic fields in the model state.
| Author and Contact information:

-  Radar Science: David Dowell, david.dowell at noaa.gov, Glen Romine, romine at ucar.edu
-  DART Code: Nancy Collins, nancy at ucar.edu
-  Original DART/Radar work: Alain Caya

Backward compatibility note
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For users of previous versions of the radar obs_def code, here are a list of changes beginning with subversion revision
3616 which are not backward compatible:

-  The namelist has changed quite a bit; some items were removed, some added, and some renamed. See the namelist
   documention in this file for the current item names and default values.
-  Some constants which depend on the microphysics scheme have been added to the namelist to make it easier to change
   the values for different schemes, but the defaults have also changed. Verify they are appropriate for the scheme
   being used.
-  The interactive create routine prompts for the beam direction differently now. It takes azimuth and elevation, and
   then does the trigonometry to compute the three internal values which are stored in the file. The previous version
   prompted for the internal values directly.
-  The get_expected routines try to call the model interpolate routine for ``QTY_POWER_WEIGHTED_FALL_SPEED`` and
   ``QTY_RADAR_REFLECTIVITY`` values. If they are not available then the code calls the model interpolate routines for
   several other quantities and computes these quantities. However, this requires that the model_mod interpolate code
   returns gracefully if the quantity is unknown or unsupported. The previous version of the WRF model_mod code used to
   print an error message and stop if the quantity was unknown. The updated version in the repository which went in with
   this radar code has been changed to return an error status code but continue if the quantity is unknown.
-  The value for gravity is currently hardcoded in this module. Previous versions of this code used the gravity constant
   in the DART types_mod.f90 code, but in reality the code should be using whatever value of gravity is being used in
   the model code. For now, the value is at least separated so users can change the value in this code if necessary.

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

=============================== =======================
*use obs_def_radar_mod, only :* read_radar_ref
\                               get_expected_radar_ref
\                               read_radial_vel
\                               write_radial_vel
\                               interactive_radial_vel
\                               get_expected_radial_vel
\                               get_obs_def_radial_vel
\                               set_radial_vel
=============================== =======================

Namelist interface ``&obs_def_radar_mod_nml`` is read from file ``input.nml``.

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *call read_radar_ref(obsvalue, refkey)*
   ::

      real(r8),                   intent(inout) :: obsvalue
      integer,                    intent(out)   :: refkey

.. container:: indent1

   Reflectivity observations have no auxiliary data to read or write, but there are namelist options that can alter the
   observation value at runtime. This routine tests the observation value and alters it if required.

   ============ =============================================================
   ``obsvalue`` Observation value.
   ``refkey``   Set to 0 to avoid uninitialized values, but otherwise unused.
   ============ =============================================================

| 

.. container:: routine

   *call get_expected_radar_ref(state_vector, location, ref, istatus)*
   ::

      real(r8),            intent(in)  :: state_vector(:)
      type(location_type), intent(in)  :: location
      real(r8),            intent(out) :: ref
      integer,             intent(out) :: istatus

.. container:: indent1

   | Given a location and the state vector from one of the ensemble members, compute the model-predicted radar
     reflectivity that would be observed at that location. The returned value is in dBZ.
   | If ``apply_ref_limit_to_fwd_op`` is .TRUE. in the namelist, reflectivity values less than
     ``reflectivity_limit_fwd_op`` will be set to ``lowest_reflectivity_fwd_op``.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_vector`` | A one dimensional representation of the model state vector                                       |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``location``     | Location of this observation                                                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``ref``          | The returned radar reflectivity value                                                            |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``istatus``      | Returned integer status code describing problems with applying forward operator. 0 is a good     |
   |                  | value; any positive value indicates an error; negative values are reserved for internal DART use |
   |                  | only.                                                                                            |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call read_radial_vel(velkey, ifile [, fform])*
   ::

      integer,                    intent(out) :: velkey
      integer,                    intent(in)  :: ifile
      character(len=*), optional, intent(in)  :: fform

.. container:: indent1

   Reads the additional auxiliary information associated with a radial velocity observation. This includes the location
   of the radar source, the beam direction, and the nyquist velocity.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``velkey`` | Unique identifier associated with this radial velocity observation. In this code it is an integer      |
   |            | index into module local arrays which hold the additional data. This routine increments it and returns  |
   |            | the new value.                                                                                         |
   +------------+--------------------------------------------------------------------------------------------------------+
   | ``ifile``  | File unit descriptor for input file                                                                    |
   +------------+--------------------------------------------------------------------------------------------------------+
   | *fform*    | File format specifier: FORMATTED or UNFORMATTED; default FORMATTED                                     |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call write_radial_vel(velkey, ifile [, fform])*
   ::

      integer,                    intent(in) :: velkey
      integer,                    intent(in) :: ifile
      character(len=*), optional, intent(in) :: fform

.. container:: indent1

   Writes the additional auxiliary information associated with a radial velocity observation. This includes the location
   of the radar source, the beam direction, and the nyquist velocity.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``velkey`` | Unique identifier associated with this radial velocity observation. In this code it is an integer      |
   |            | index into module local arrays which hold the additional data. This routine uses the value to select   |
   |            | the appropriate data to write for this observation.                                                    |
   +------------+--------------------------------------------------------------------------------------------------------+
   | ``ifile``  | File unit descriptor for output file                                                                   |
   +------------+--------------------------------------------------------------------------------------------------------+
   | *fform*    | File format specifier: FORMATTED or UNFORMATTED; default FORMATTED                                     |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_obs_def_radial_vel(velkey, radar_location, beam_direction, nyquist_velocity)*
   ::

      integer,             intent(in)  :: velkey
      type(location_type), intent(out) :: radar_location
      real(r8),            intent(out) :: beam_direction(3)
      real(r8),            intent(out) :: nyquist_velocity

.. container:: indent1

   Returns the auxiliary information associated with a given radial velocity observation.

   +----------------------+----------------------------------------------------------------------------------------------+
   | ``velkey``           | Unique identifier associated with this radial velocity observation. In this code it is an    |
   |                      | integer index into module local arrays which hold the additional data. This routine uses the |
   |                      | value to select the appropriate data to return.                                              |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``radar_location``   | Location of the radar.                                                                       |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``beam_orientation`` | Orientation of the radar beam at the observation location. The three values are:             |
   |                      | sin(azimuth)*cos(elevation), cos(azimuth)*cos(elevation), and sin(elevation).                |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``nyquist_velocity`` | Nyquist velocity at the observation point in meters/second.                                  |
   +----------------------+----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call set_radial_vel(velkey, radar_location, beam_direction, nyquist_velocity)*
   ::

      integer,             intent(out) :: velkey
      type(location_type), intent(in)  :: radar_location
      real(r8),            intent(in)  :: beam_direction(3)
      real(r8),            intent(in)  :: nyquist_velocity

.. container:: indent1

   Sets the auxiliary information associated with a radial velocity observation. This routine increments and returns the
   new key associated with these values.

   +----------------------+----------------------------------------------------------------------------------------------+
   | ``velkey``           | Unique identifier associated with this radial velocity observation. In this code it is an    |
   |                      | integer index into module local arrays which hold the additional data. This routine returns  |
   |                      | the incremented value associated with this data.                                             |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``radar_location``   | Location of the radar.                                                                       |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``beam_orientation`` | Orientation of the radar beam at the observation location. The three values are:             |
   |                      | sin(azimuth)*cos(elevation), cos(azimuth)*cos(elevation), and sin(elevation).                |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``nyquist_velocity`` | Nyquist velocity at the observation point in meters/second.                                  |
   +----------------------+----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call interactive_radial_vel(velkey)*
   ::

      integer, intent(out) :: velkey

.. container:: indent1

   Prompts the user for the auxiliary information needed for a radial velocity observation, and returns the new key
   associated with this data.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``velkey`` | Unique identifier associated with this radial velocity observation. In this code it is an integer      |
   |            | index into module local arrays which hold the additional data. This routine returns the incremented    |
   |            | value associated with this data.                                                                       |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_expected_radial_vel(state_vector, location, velkey, radial_vel, istatus)*
   ::

      real(r8),            intent(in)  :: state_vector(:)
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: velkey
      real(r8),            intent(out) :: radial_vel
      integer,             intent(out) :: istatus

.. container:: indent1

   | Given a location and the state vector from one of the ensemble members, compute the model-predicted radial velocity
     in meters/second that would be observed at that location. ``velkey`` is the unique index for this particular radial
     velocity observation. The value is returned in ``radial_vel``, ``istatus`` is the return code.
   | The along-beam component of the 3-d air velocity is computed from the u, v, and w fields plus the beam_direction.
     The along-beam component of power-weighted precipitation fall velocity is added to the result.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_vector`` | A one dimensional representation of the model state vector                                       |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``location``     | Location of this observation                                                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``velkey``       | Unique identifier associated with this radial velocity observation                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``radial_vel``   | The returned radial velocity value in meters/second                                              |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``istatus``      | Returned integer status code describing problems with applying forward operator. 0 is a good     |
   |                  | value; any positive value indicates an error; negative values are reserved for internal DART use |
   |                  | only.                                                                                            |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_def_radar_mod_nml
      apply_ref_limit_to_obs      =   .false.,
      reflectivity_limit_obs      =     -10.0,
      lowest_reflectivity_obs     =     -10.0,
      apply_ref_limit_to_fwd_op   =   .false.,
      reflectivity_limit_fwd_op   =     -10.0,
      lowest_reflectivity_fwd_op  =     -10.0,
      max_radial_vel_obs          =   1000000,
      allow_wet_graupel           =   .false.,
      microphysics_type           =       2  ,
      allow_dbztowt_conv          =   .false.,
      dielectric_factor           =     0.224,
      n0_rain                     =     8.0e6,
      n0_graupel                  =     4.0e6,
      n0_snow                     =     3.0e6,
      rho_rain                    =    1000.0,
      rho_graupel                 =     400.0,
      rho_snow                    =     100.0 
      /

| 

.. container::

   +-----------------------------+------------+---------------------------------------+
   | Item                        | Type       | Description                           |
   +=============================+============+=======================================+
   | apply_ref_limit_to_obs      | logical    | If .TRUE. replace all reflectivity    |
   |                             |            | values less than                      |
   |                             |            | "reflectivity_limit_obs" with         |
   |                             |            | "lowest_reflectivity_obs" value. If   |
   |                             |            | .FALSE. leave all values as-is.       |
   +-----------------------------+------------+---------------------------------------+
   | reflectivity_limit_obs      | real(r8)   | The threshold value. Observed         |
   |                             |            | reflectivity values less than this    |
   |                             |            | threshold will be set to the          |
   |                             |            | "lowest_reflectivity_obs" value.      |
   |                             |            | Units are dBZ.                        |
   +-----------------------------+------------+---------------------------------------+
   | lowest_reflectivity_obs     | real(r8)   | The 'set-to' value. Observed          |
   |                             |            | reflectivity values less than the     |
   |                             |            | threshold will be set to this value.  |
   |                             |            | Units are dBZ.                        |
   +-----------------------------+------------+---------------------------------------+
   | apply_ref_limit_to_fwd_op   | logical    | Same as "apply_ref_limit_to_obs", but |
   |                             |            | for the forward operator.             |
   +-----------------------------+------------+---------------------------------------+
   | reflectivity_limit_fwd_op   | real(r8)   | Same as "reflectivity_limit_obs", but |
   |                             |            | for the forward operator values.      |
   +-----------------------------+------------+---------------------------------------+
   | lowest_reflectivity_fwd_op  | real(r8)   | Same as "lowest_reflectivity_obs",    |
   |                             |            | but for the forward operator values.  |
   +-----------------------------+------------+---------------------------------------+
   | max_radial_vel_obs          | integer    | Maximum number of observations of     |
   |                             |            | this type to support at run time.     |
   |                             |            | This is combined total of all obs_seq |
   |                             |            | files, for example the observation    |
   |                             |            | diagnostic program potentially opens  |
   |                             |            | multiple obs_seq.final files, or the  |
   |                             |            | obs merge program can also open       |
   |                             |            | multiple obs files.                   |
   +-----------------------------+------------+---------------------------------------+
   | allow_wet_graupel           | logical    | It is difficult to predict/diagnose   |
   |                             |            | whether graupel/hail has a wet or dry |
   |                             |            | surface. Even when the temperature is |
   |                             |            | above freezing, evaporation and/or    |
   |                             |            | absorption can still result in a dry  |
   |                             |            | surface. This issue is important      |
   |                             |            | because the reflectivity from graupel |
   |                             |            | with a wet surface is significantly   |
   |                             |            | greater than that from graupel with a |
   |                             |            | dry surface. Currently, the user has  |
   |                             |            | two options for how to compute        |
   |                             |            | graupel reflectivity. If              |
   |                             |            | allow_wet_graupel is .false. (the     |
   |                             |            | default), then graupel is always      |
   |                             |            | assumed to be dry. If                 |
   |                             |            | allow_wet_graupel is .true., then     |
   |                             |            | graupel is assumed to be wet (dry)    |
   |                             |            | when the temperature is above (below) |
   |                             |            | freezing. A consequence is that a     |
   |                             |            | sharp gradient in reflectivity will   |
   |                             |            | be produced at the freezing level. In |
   |                             |            | the future, it might be better to     |
   |                             |            | provide the option of having a        |
   |                             |            | transition layer.                     |
   +-----------------------------+------------+---------------------------------------+
   | microphysics_type           | integer    | If the state vector contains the      |
   |                             |            | reflectivity or the power weighted    |
   |                             |            | fall speed, interpolate directly from |
   |                             |            | those regardless of the setting of    |
   |                             |            | this item. If the state vector does   |
   |                             |            | not contain the fields, this value    |
   |                             |            | should be set to be compatible with   |
   |                             |            | whatever microphysical scheme is      |
   |                             |            | being used by the model. If the model |
   |                             |            | is using a different microphysical    |
   |                             |            | scheme but has compatible fields to   |
   |                             |            | the ones listed below, setting this   |
   |                             |            | value will select the scheme to use.  |
   |                             |            |                                       |
   |                             |            | -  1 = Kessler scheme.                |
   |                             |            | -  2 = Lin et al. microphysics        |
   |                             |            | -  3 = User selected scheme where 10  |
   |                             |            |    cm reflectivity and power weighted |
   |                             |            |    fall velocity are expected in the  |
   |                             |            |    state vector (failure if not       |
   |                             |            |    found)                             |
   |                             |            | -  4 = User selected scheme where     |
   |                             |            |    only power weighted fall velocity  |
   |                             |            |    is expected (failure if not found) |
   |                             |            | -  5 = User selected scheme where     |
   |                             |            |    only reflectivity is expected      |
   |                             |            |    (failure if not found)             |
   |                             |            | -  -1 = ASSUME FALL VELOCITY IS ZERO, |
   |                             |            |    allows over-riding the failure     |
   |                             |            |    modes above if reflectivity and/or |
   |                             |            |    fall velocity are not available    |
   |                             |            |    but a result is desired for        |
   |                             |            |    testing purposes only.             |
   +-----------------------------+------------+---------------------------------------+
   | allow_dbztowt_conv          | logical    | Flag to enable use of the dbztowt     |
   |                             |            | routine where reflectivity is         |
   |                             |            | available, but not the power-weighted |
   |                             |            | fall velocity. This scheme uses       |
   |                             |            | emperical relations between           |
   |                             |            | reflectivity and fall velocity, with  |
   |                             |            | poor accuracy for highly reflective,  |
   |                             |            | low density particles (such as water  |
   |                             |            | coated snow aggregates). Expect       |
   |                             |            | questionable accuracy in radial       |
   |                             |            | velocity from the forward operator    |
   |                             |            | with high elevation angles where ice  |
   |                             |            | is present in the model state.        |
   +-----------------------------+------------+---------------------------------------+
   | dielectric_factor           | real(r8)   | According to Smith (1984), there are  |
   |                             |            | two choices for the dielectric factor |
   |                             |            | depending on how the snow particle    |
   |                             |            | sizes are specified. If melted        |
   |                             |            | raindrop diameters are used, then the |
   |                             |            | factor is 0.224. If equivalent ice    |
   |                             |            | sphere diameters are used, then the   |
   |                             |            | factor is 0.189. The default is set   |
   |                             |            | to use the common convention of       |
   |                             |            | melted raindrop diameters.            |
   +-----------------------------+------------+---------------------------------------+
   | n0_rain                     | real(r8)   | Intercept parameters (m^-4) for size  |
   |                             |            | distributions of each hydrometeor.    |
   |                             |            | The default of 8.0e6 is for the Lin   |
   |                             |            | et al. microphysics scheme with the   |
   |                             |            | Hobbs settings for graupel/hail. (The |
   |                             |            | Hobbs graupel settings are also the   |
   |                             |            | default for the Lin scheme in WRF 2.2 |
   |                             |            | and 3.0.)                             |
   +-----------------------------+------------+---------------------------------------+
   | n0_graupel                  | real(r8)   | Intercept parameters (m^-4) for size  |
   |                             |            | distributions of each hydrometeor.    |
   |                             |            | The default of 4.0e6 is for the Lin   |
   |                             |            | et al. microphysics scheme with the   |
   |                             |            | Hobbs settings for graupel/hail. (The |
   |                             |            | Hobbs graupel settings are also the   |
   |                             |            | default for the Lin scheme in WRF 2.2 |
   |                             |            | and 3.0.)                             |
   +-----------------------------+------------+---------------------------------------+
   | n0_snow                     | real(r8)   | Intercept parameters (m^-4) for size  |
   |                             |            | distributions of each hydrometeor.    |
   |                             |            | The default of 3.0e6 is for the Lin   |
   |                             |            | et al. microphysics scheme with the   |
   |                             |            | Hobbs settings for graupel/hail. (The |
   |                             |            | Hobbs graupel settings are also the   |
   |                             |            | default for the Lin scheme in WRF 2.2 |
   |                             |            | and 3.0.)                             |
   +-----------------------------+------------+---------------------------------------+
   | rho_rain                    | real(r8)   | Density (kg m^-3) of each hydrometeor |
   |                             |            | type. The default of 1000.0 is for    |
   |                             |            | the Lin et al. microphysics scheme    |
   |                             |            | with the Hobbs setting for            |
   |                             |            | graupel/hail.                         |
   +-----------------------------+------------+---------------------------------------+
   | rho_graupel                 | real(r8)   | Density (kg m^-3) of each hydrometeor |
   |                             |            | type. The default of 400.0 is for the |
   |                             |            | Lin et al. microphysics scheme with   |
   |                             |            | the Hobbs setting for graupel/hail.   |
   +-----------------------------+------------+---------------------------------------+
   | rho_snow                    | real(r8)   | Density (kg m^-3) of each hydrometeor |
   |                             |            | type. The default of 100.0 is for the |
   |                             |            | Lin et al. microphysics scheme with   |
   |                             |            | the Hobbs setting for graupel/hail.   |
   +-----------------------------+------------+---------------------------------------+

| 

Files
-----

-  A DART observation sequence file containing Radar obs.

References
----------

-  Battan, L. J., 1973: *Radar Observation of the Atmosphere.* Univ. of Chicago Press, 324 pp.
-  Caya, A. *Radar Observations in Dart.* DART Subversion repository.
-  Doviak, R. J., and D. S. Zrnic, 1993: *Doppler Radar and Weather Observations.* Academic Press, 562 pp.
-  Ferrier, B. S., 1994: A double-moment multiple-phase four-class bulk ice scheme. Part I: Description. *J. Atmos.
   Sci.*, **51**, 249-280.
-  Lin, Y.-L., Farley R. D., and H. D. Orville, 1983: Bulk parameterization of the snow field in a cloud model. *J.
   Climate Appl. Meteor.*, **22**, 1065-1092.
-  Smith, P. L. Jr., 1984: Equivalent radar reflectivity factors for snow and ice particles. *J. Climate Appl. Meteor.*,
   23, 1258-1260.
-  Smith, P. L. Jr., Myers C. G., and H. D. Orville, 1975: Radar reflectivity factor calculations in numerical cloud
   models using bulk parameterization of precipitation. *J. Appl. Meteor.*, **14**, 1156-1165.


Error codes and conditions
--------------------------

+-----------------------+-------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
|        Routine        |                                       Message                                       |                                                                      Comment                                                                      |
+=======================+=====================================================================================+===================================================================================================================================================+
| initialize_module     | initial allocation failed for radial vel obs data, itemcount = (max_radial_vel_obs) | Need to increase max_radial_vel_obs count in namelist                                                                                             |
+-----------------------+-------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| read_radial_vel       | Expected location header "platform" in input file                                   | The format of the input file is not consistent.                                                                                                   |
+-----------------------+-------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| velkey_out_of_range   | velkey (val) exceeds max_radial_vel_obs (maxval)                                    | The number of radial velocity observations exceeds the array size allocated in the module. Need to increase max_radial_vel_obs count in namelist. |
+-----------------------+-------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| read_nyquist_velocity | bad value for nyquist velocity                                                      | The format of the input obs_seq file is not consistent.                                                                                           |
+-----------------------+-------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| read_beam_direction   | beam_direction value must be between -1 and 1, got ()                               | The format of the input obs_seq file is not consistent.                                                                                           |
+-----------------------+-------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| read_beam_direction   | Expected orientation header "dir3d" in input file                                   | The format of the input obs_seq file is not consistent.                                                                                           |
+-----------------------+-------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+


Private components
------------------

=============================== ============================
*use obs_def_radar_mod, only :* initialize_module
\                               read_beam_direction
\                               read_nyquist_velocity
\                               write_beam_direction
\                               write_nyquist_velocity
\                               interactive_beam_direction
\                               interactive_nyquist_velocity
\                               get_reflectivity
\                               get_precip_fall_speed
\                               initialize_constants
\                               print_constants
\                               pr_con
\                               velkey_out_of_range
\                               check_namelist_limits
\                               ascii_file_format
=============================== ============================

| 

.. container:: routine

   *call initialize_module()*

.. container:: indent1

   Reads the namelist, allocates space for the auxiliary data associated wtih radial velocity observations, initializes
   the constants used in subsequent computations (possibly altered by values in the namelist), and prints out the list
   of constants and the values in use. These may need to change depending on which microphysics scheme is being used.

| 

.. container:: routine

   *beam_direction = read_beam_direction(ifile, is_asciiformat)*
   ::

      real(r8), dimension(3)            :: read_beam_direction
      integer,               intent(in) :: ifile
      logical,               intent(in) :: is_asciiformat

.. container:: indent1

   Reads the beam direction at the observation location. Auxiliary data for doppler radial velocity observations.

   +-------------------------+-------------------------------------------------------------------------------------------+
   | ``read_beam_direction`` | Returns three real values for the radar beam orientation                                  |
   +-------------------------+-------------------------------------------------------------------------------------------+
   | ``ifile``               | File unit descriptor for input file                                                       |
   +-------------------------+-------------------------------------------------------------------------------------------+
   | ``is_asciiformat``      | File format specifier: .TRUE. if file is formatted/ascii, or .FALSE. if                   |
   |                         | unformatted/binary. Default .TRUE.                                                        |
   +-------------------------+-------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *nyquist_velocity = read_nyquist_velocity(ifile, is_asciiformat)*
   ::

      real(r8),            :: read_nyquist_velocity
      integer,  intent(in) :: ifile
      logical,  intent(in) :: is_asciiformat

.. container:: indent1

   Reads nyquist velocity for a doppler radial velocity observation.

   +---------------------------+-----------------------------------------------------------------------------------------+
   | ``read_nyquist_velocity`` | Returns a real value for the nyquist velocity value                                     |
   +---------------------------+-----------------------------------------------------------------------------------------+
   | ``ifile``                 | File unit descriptor for input file                                                     |
   +---------------------------+-----------------------------------------------------------------------------------------+
   | ``is_asciiformat``        | File format specifier: .TRUE. if file is formatted/ascii, or .FALSE. if                 |
   |                           | unformatted/binary. Default .TRUE.                                                      |
   +---------------------------+-----------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call write_beam_direction(ifile, beam_direction, is_asciiformat)*
   ::

      integer,                intent(in) :: ifile
      real(r8), dimension(3), intent(in) :: beam_direction
      logical,                intent(in) :: is_asciiformat

.. container:: indent1

   Writes the beam direction at the observation location. Auxiliary data for doppler radial velocity observations.

   +--------------------+------------------------------------------------------------------------------------------------+
   | ``ifile``          | File unit descriptor for output file                                                           |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``beam_direction`` | Three components of the radar beam orientation                                                 |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``is_asciiformat`` | File format specifier: .TRUE. if file is formatted/ascii, or .FALSE. if unformatted/binary.    |
   |                    | Default .TRUE.                                                                                 |
   +--------------------+------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call write_nyquist_velocity(ifile, nyquist_velocity, is_asciiformat)*
   ::

      integer,  intent(in) :: ifile
      real(r8), intent(in) :: nyquist_velocity
      logical,  intent(in) :: is_asciiformat

.. container:: indent1

   Writes nyquist velocity for a doppler radial velocity observation.

   +----------------------+----------------------------------------------------------------------------------------------+
   | ``ifile``            | File unit descriptor for output file                                                         |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``nyquist_velocity`` | The nyquist velocity value for this observation                                              |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``is_asciiformat``   | File format specifier: .TRUE. if file is formatted/ascii, or .FALSE. if unformatted/binary.  |
   |                      | Default .TRUE.                                                                               |
   +----------------------+----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call interactive_beam_direction(beam_direction)*
   ::

      real(r8), dimension(3), intent(out) :: beam_direction

.. container:: indent1

   Prompts the user for input for the azimuth and elevation of the radar beam at the observation location. Will be
   converted to the three values actually stored in the observation sequence file.

   ================== ==============================================
   ``beam_direction`` Three components of the radar beam orientation
   ================== ==============================================

| 

.. container:: routine

   *call interactive_nyquist_velocity(nyquist_velocity)*
   ::

      real(r8), intent(out) :: nyquist_velocity

.. container:: indent1

   Prompts the user for input for the nyquist velocity value associated with a doppler radial velocity observation.

   ==================== ===========================================
   ``nyquist_velocity`` Nyquist velocity value for the observation.
   ==================== ===========================================

| 

.. container:: routine

   *call get_reflectivity(qr, qg, qs, rho, temp, ref)*
   ::

      real(r8), intent(in)  :: qr
      real(r8), intent(in)  :: qg
      real(r8), intent(in)  :: qs
      real(r8), intent(in)  :: rho
      real(r8), intent(in)  :: temp
      real(r8), intent(out) :: ref

.. container:: indent1

   Computes the equivalent radar reflectivity factor in mm\ :sup:`6` m\ :sup:`-3` for simple single-moment microphysics
   schemes such as Kessler and Lin, et al. See the references for more details.

   ======== =======================================
   ``qr``   Rain water content (kg kg\ :sup:`-1`)
   ``qg``   Graupel/hail content (kg kg\ :sup:`-1`)
   ``qs``   Snow content (kg kg\ :sup:`-1`)
   ``rho``  Air density (kg m\ :sup:`-3`)
   ``temp`` Air temperature (K)
   ``ref``  The returned radar reflectivity value
   ======== =======================================

| 

.. container:: routine

   *call get_precip_fall_speed(qr, qg, qs, rho, temp, precip_fall_speed)*
   ::

      real(r8), intent(in)  :: qr
      real(r8), intent(in)  :: qg
      real(r8), intent(in)  :: qs
      real(r8), intent(in)  :: rho
      real(r8), intent(in)  :: temp
      real(r8), intent(out) :: precip_fall_speed

.. container:: indent1

   Computes power-weighted precipitation fall speed in m s\ :sup:`-1` for simple single-moment microphysics schemes such
   as Kessler and Lin, et al. See the references for more details.

   ===================== =======================================
   ``qr``                Rain water content (kg kg\ :sup:`-1`)
   ``qg``                Graupel/hail content (kg kg\ :sup:`-1`)
   ``qs``                Snow content (kg kg\ :sup:`-1`)
   ``rho``               Air density (kg m\ :sup:`-3`)
   ``temp``              Air temperature (K)
   ``precip_fall_speed`` The returned precipitation vall speed
   ===================== =======================================

| 

.. container:: routine

   *call initialize_constants()*

.. container:: indent1

   Set values for a collection of constants used throughout the module during the various calculations. These are set
   once in this routine and are unchanged throughout the rest of the execution. They cannot be true Fortran
   ``parameters`` because some of the values can be overwritten by namelist entries, but once they are set they are
   treated as read-only parameters.

| 

.. container:: routine

   *call print_constants()*

.. container:: indent1

   Print out the names and values of all constant parameters used by this module. The error handler message facility is
   used to print the message, which by default goes to both the DART log file and to the standard output of the program.

| 

.. container:: routine

   *call pr_con(c_val, c_str)*
   ::

      real(r8),         intent(in)  :: c_val
      character(len=*), intent(in)  :: c_str

.. container:: indent1

   Calls the DART error handler routine to print out a string label and a real value to both the log file and to the
   standard output.

   ===================== ===================
   ``Value of constant`` A real value.
   ``Name of constant``  A character string.
   ===================== ===================

| 

.. container:: routine

   *call velkey_out_of_range(velkey)*
   ::

      integer, intent(in)  :: velkey

.. container:: indent1

   Range check key and trigger a fatal error if larger than the allocated array for observation auxiliary data.

   ========== =============================================================
   ``velkey`` Integer key into a local array of auxiliary observation data.
   ========== =============================================================

| 

.. container:: routine

   *call check_namelist_limits(apply_ref_limit_to_obs, reflectivity_limit_obs, lowest_reflectivity_obs,
   apply_ref_limit_to_fwd_op, reflectivity_limit_fwd_op, lowest_reflectivity_fwd_op)*
   ::

      logical,  intent(in) :: apply_ref_limit_to_obs
      real(r8), intent(in) :: reflectivity_limit_obs
      real(r8), intent(in) :: lowest_reflectivity_obs
      logical,  intent(in) :: apply_ref_limit_to_fwd_op
      real(r8), intent(in) :: reflectivity_limit_fwd_op
      real(r8), intent(in) :: lowest_reflectivity_fwd_op

.. container:: indent1

   Check the values set in the namelist for consistency. Print out a message if the limits and set-to values are
   different; this may be intentional but is not generally expected to be the case. In all cases below, see the namelist
   documentation for a fuller explanation of each value.

   ============================== =========================
   ``apply_ref_limit_to_obs``     Logical. See namelist.
   ``reflectivity_limit_obs``     Real value. See namelist.
   ``lowest_reflectivity_obs``    Real value. See namelist.
   ``apply_ref_limit_to_fwd_op``  Logical. See namelist.
   ``reflectivity_limit_fwd_op``  Real value. See namelist.
   ``lowest_reflectivity_fwd_op`` Real value. See namelist.
   ============================== =========================

| 

.. container:: routine

   *is_asciifile = ascii_file_format(fform)*
   ::

      logical                                :: ascii_file_format
      character(len=*), intent(in), optional :: fform

.. container:: indent1

   Should be moved to DART utility module at some point. Returns .TRUE. if the optional argument is missing or if it is
   not one of the following values: ``"unformatted", "UNFORMATTED", "unf", "UNF"``.

   ===================== ========================================
   ``ascii_file_format`` Return value. Logical. Default is .TRUE.
   ``fform``             Character string file format.
   ===================== ========================================

| 
