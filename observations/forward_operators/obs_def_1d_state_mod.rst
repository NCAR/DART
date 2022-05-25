MODULE ``obs_def_1d_state_mod``
===============================

Overview
--------

The list of observation types to be supported by the DART executables is defined at compile time. The observations DART
supports can be changed at any time by adding or removing items from the preprocess namelist and rerunning
*quickbuild.sh*.

| ``Preprocess`` takes observation specific code sections from special obs_def files to generate ``obs_def_mod.f90`` and
  ``obs_kind_mod.f90`` which are then compiled into ``filter`` and other DART programs. One of the motivations behind
  creating ``obs_def_1d_state_mod.f90`` was to provide a prototype for people developing more complicated specialized
  observation definition modules.
| ``Obs_def_1d_state_mod.f90`` is an extended format Fortran 90 module that provides the definition for observation
  types designed for use with idealized low-order models that use the 1D location module and can be thought of as having
  a state vector that is equally spaced on a 1D cyclic domain. Observation types include:

-  RAW_STATE_VARIABLE - A straight linear interpolation to a point on a [0,1] domain.
-  RAW_STATE_VAR_POWER - The interpolated RAW_STATE_VARIABLE raised to a real-valued power.
-  RAW_STATE_1D_INTEGRAL - An area-weighted 'integral' of the state variable over some part of the cyclic 1D domain.

| RAW_STATE_VAR_POWER is convenient for studying non-gaussian, non-linear assimilation problems. RAW_STATE_VAR_POWER can
  be used to do idealized studies related to remote sensing observations that are best thought of as weighted integrals
  of some quantity over a finite volume.
| The RAW_STATE_1D_INTEGRAL has an associated half_width and localization type (see the
  :doc:`../../assimilation_code/modules/assimilation/cov_cutoff_mod` documentation) and a number of points at which to
  compute the associated integral by quadrature. The location of the observation defines the center of mass of the
  integral. The integral is centered around the location and extends outward on each side to 2*half_width. The weight
  associated with the integral is defined by the weight of the localization function (for instance Gaspari Cohn) using
  the same localization options as defined by the cov_cutoff module. The number of points are used to equally divide the
  range for computing the integral by quadrature.
| Special observation modules like ``obs_def_1d_state_mod.f90`` contain Fortran 90 code *and* additional specially
  formatted commented code that is used to guide the preprocess program in constructing obs_def_mod.f90 and
  obs_kind_mod.f90. The specially formatted comments are most conveniently placed at the beginning of the module and
  comprise seven sections, each beginning and ending with a special F90 comment line that must be included *verbatim*.
| The seven sections and their specific instances for the 1d_raw_state_mod are:

#. A list of all observation types defined by this module and their associated generic quantities (see
   :doc:`../../assimilation_code/programs/preprocess/preprocess` for details on quantity files). The header line is
   followed by lines that have the observation type name (an all caps Fortran 90 identifier) and their associated
   generic quantity identifier. If there is no special processing needed for an observation type, and no additional data
   needed beyond the standard contents of an observation then a third word on the line, ``COMMON_CODE``, will instruct
   the preprocess program to automatically generate all stubs and code needed for this type. For observation types
   needing special code or additional data, this word should not be specified and the user must supply the code
   manually.

   ::

      ! BEGIN DART PREPROCESS KIND LIST
      ! RAW_STATE_VARIABLE,    QTY_STATE_VARIABLE,   COMMON_CODE
      ! RAW_STATE_1D_INTEGRAL, QTY_1D_INTEGRAL
      ! END DART PREPROCESS KIND LIST

   | 

#. A list of all the use statements that the completed obs_def_mod.f90 must have in order to use the public interfaces
   provided by this special obs_def module. This section is optional if there are no external interfaces.

   ::

      ! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
      !   use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral,  &
      !                                    interactive_1d_integral, get_expected_1d_integral, &
      !                                    set_1d_integral
      ! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

   | 

#. Case statement entries for each observation type defined by this special obs_def module stating how to compute the
   forward observation operator. There must be a case statement entry for each type of observation, *except* for
   observation types defined with COMMON_CODE.

   ::

      ! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
      !         case(RAW_STATE_1D_INTEGRAL)
      !            call get_expected_1d_integral(state, location, obs_def%key, obs_val, istatus)
      ! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

   | 

#. Case statement entries for each observation type defined by this special obs_def module stating how to read any extra
   required information from an obs sequence file. There must be a case statement entry for each type of observation,
   *except* for observation types defined with COMMON_CODE. If no special action is required put a ``continue``
   statement as the body of the case instead of a subroutine call.

   ::

      ! BEGIN DART PREPROCESS READ_OBS_DEF
      !      case(RAW_STATE_1D_INTEGRAL)
      !         call read_1d_integral(obs_def%key, ifile, fform)
      ! END DART PREPROCESS READ_OBS_DEF

   | 

#. Case statement entries for each observation type defined by this special obs_def module stating how to write any
   extra required information to an obs sequence file. There must be a case statement entry for each type of
   observation, *except* for observation types defined with COMMON_CODE. If no special action is required put a
   ``continue`` statement as the body of the case instead of a subroutine call.

   ::

      ! BEGIN DART PREPROCESS WRITE_OBS_DEF
      !      case(RAW_STATE_1D_INTEGRAL)
      !         call write_1d_integral(obs_def%key, ifile, fform)
      ! END DART PREPROCESS WRITE_OBS_DEF

   | 

#. Case statement entries for each observation type defined by this special obs_def module stating how to interactively
   create any extra required information. There must be a case statement entry for each type of observation, *except*
   for observation types defined with COMMON_CODE. If no special action is required put a ``continue`` statement as the
   body of the case instead of a subroutine call.

   ::

      ! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
      !      case(RAW_STATE_1D_INTEGRAL)
      !         call interactive_1d_integral(obs_def%key)
      ! END DART PREPROCESS INTERACTIVE_OBS_DEF

   | 

#. Any executable F90 module code must be tagged with the following comments. All lines between these markers will be
   copied, verbatim, to obs_def_mod.f90. This section is not required if there are no observation-specific subroutines.

   ::

      ! BEGIN DART PREPROCESS MODULE CODE
      module obs_def_1d_state_mod

      ... (module executable code)

      end module obs_def_1d_state_mod
      ! END DART PREPROCESS MODULE CODE

   | 

| 

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod (1d_location_mod_only)
   time_manager_mod
   assim_model_mod
   cov_cutoff_mod

Public interfaces
-----------------

========================= ========================
*use obs_def_mod, only :* write_1d_integral
\                         read_1d_integral
\                         interactive_1d_integral
\                         get_expected_1d_integral
\                         set_1d_integral
\                         write_power
\                         read_power
\                         interactive_power
\                         get_expected_power
\                         set_power
========================= ========================

| 

.. container:: routine

   *call write_1d_integral(igrkey, ifile, fform)*
   ::

      integer,          intent(in) :: igrkey
      integer,          intent(in) :: ifile
      character(len=*), intent(in) :: fform

.. container:: indent1

   Writes out the extra information for observation with unique identifier key for a 1d_integral observation type. This
   includes the half-width, localization type and number of quadrature points for this observation.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``igrkey`` | Unique integer key associated with the 1d integral observation being processed. This is not the same   |
   |            | as the key that all types of observations have and uniquely distinguishes all observations from each   |
   |            | other; this is a key that is only set and retrieved by this code for 1d integral observations. It is   |
   |            | stored in the obs_def derived type, not in the main obs_type definition.                               |
   +------------+--------------------------------------------------------------------------------------------------------+
   | ``ifile``  | Unit number on which observation sequence file is open                                                 |
   +------------+--------------------------------------------------------------------------------------------------------+
   | ``fform``  | String noting whether file is opened for 'formatted' or 'unformatted' IO.                              |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call read_1d_integral(igrkey, ifile, fform)*
   ::

      integer,          intent(out) :: igrkey
      integer,          intent(in)  :: ifile
      character(len=*), intent(in)  :: fform

.. container:: indent1

   Reads the extra information for observation with unique identifier key for a 1d_integral observation type. This
   information includes the half-width, localization type and number of quadrature points for this observation. The key
   that is returned is uniquely associated with the definition that has been created and is used by this module to keep
   track of the associated parameters for this observation.

   ========== =========================================================================
   ``igrkey`` Unique integer key associated with the observation being processed.
   ``ifile``  Unit number on which observation sequence file is open
   ``fform``  String noting whether file is opened for 'formatted' or 'unformatted' IO.
   ========== =========================================================================

| 

.. container:: routine

   *call interactive_1d_integral(igrkey)*
   ::

      integer, intent(out) :: igrkey

.. container:: indent1

   Uses input from standard in to define the characteristics of a 1D integral observation. The key that is returned is
   uniquely associated with the definition that has been created and can be used by this module to keep track of the
   associated parameters (half_width, localization option, number of quadrature points) for this key.

   ========== =========================================================================================
   ``igrkey`` Unique identifier associated with the created observation definition in the obs sequence.
   ========== =========================================================================================

| 

.. container:: routine

   *call get_expected_1d_integral(state, location, igrkey, val, istatus)*
   ::

      real(r8), intent(in)            :: state
      type(location_type), intent(in) :: location
      integer, intent(in)             :: igrkey
      real(r8), intent(out)           :: val
      integer, intent(out)            :: istatus

.. container:: indent1

   Computes the forward observation operator for a 1d integral observation. Calls the ``interpolate()`` routine multiple
   times to invoke the forward operator code in whatever model this has been compiled with.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``state``    | Model state vector (or extended state vector).                                                       |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``location`` | Location of this observation.                                                                        |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``igrkey``   | Unique integer key associated with this observation.                                                 |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``val``      | Returned value of forward observation operator.                                                      |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``istatus``  | Returns 0 if forward operator was successfully computed, else returns a positive value. (Negative    |
   |              | values are reserved for system use.)                                                                 |
   +--------------+------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call set_1d_integral(integral_half_width, num_eval_pts, localize_type, igrkey, istatus)*
   ::

      real(r8), intent(in)  :: integral_half_width
      integer,  intent(in)  :: num_eval_pts
      integer,  intent(in)  :: localize_type
      integer,  intent(out) :: igrkey
      integer,  intent(out) :: istatus

.. container:: indent1

   Available for use by programs that create observations to set the additional metadata for these observation types.
   This information includes the integral half-width, localization type and number of quadrature points for this
   observation. The key that is returned is uniquely associated with the definition that has been created and should be
   set in the obs_def structure by calling ``set_obs_def_key()``. This key is different from the main observation key
   which all observation types have. This key is unique to this observation type and is used when reading in the
   observation sequence to match the corresponding metadata with each observation of this type.

   ======================= ====================================================================
   ``integral_half_width`` Real value setting the half-width of the integral.
   ``num_eval_pts``        Integer, number of evaluation points. 5-20 recommended.
   ``localize_type``       Integer localization type: 1=Gaspari-Cohn; 2=Boxcar; 3=Ramped Boxcar
   ``igrkey``              Unique integer key associated with the observation being processed.
   ``istatus``             Return code. 0 means success, any other value is an error
   ======================= ====================================================================

| 

.. container:: routine

   *call write_power(powkey, ifile, fform)*
   ::

      integer,          intent(in) :: powkey
      integer,          intent(in) :: ifile
      character(len=*), intent(in) :: fform

.. container:: indent1

   Writes out the extra information, the power, for observation with unique identifier key for a power observation type.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``powkey`` | Unique integer key associated with the power observation being processed. This is not the same as the  |
   |            | key that all types of observations have and uniquely distinguishes all observations from each other;   |
   |            | this is a key that is only set and retrieved by this code for power observations. It is stored in the  |
   |            | obs_def derived type, not in the main obs_type definition.                                             |
   +------------+--------------------------------------------------------------------------------------------------------+
   | ``ifile``  | Unit number on which observation sequence file is open                                                 |
   +------------+--------------------------------------------------------------------------------------------------------+
   | ``fform``  | String noting whether file is opened for 'formatted' or 'unformatted' IO.                              |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call read_power(powkey, ifile, fform)*
   ::

      integer,          intent(out) :: powkey
      integer,          intent(in)  :: ifile
      character(len=*), intent(in)  :: fform

.. container:: indent1

   Reads the extra information, the power, for observation with unique identifier key for a power observation type. The
   key that is returned is uniquely associated with the definition that has been created and is used by this module to
   keep track of the associated parameters for this observation.

   ========== =========================================================================
   ``powkey`` Unique integer key associated with the observation being processed.
   ``ifile``  Unit number on which observation sequence file is open
   ``fform``  String noting whether file is opened for 'formatted' or 'unformatted' IO.
   ========== =========================================================================

| 

.. container:: routine

   *call interactive_power(powkey)*
   ::

      integer, intent(out) :: powkey

.. container:: indent1

   Uses input from standard in to define the characteristics of a power observation. The key that is returned is
   uniquely associated with the definition that has been created and can be used by this module to keep track of the
   associated parameter, the power, for this key.

   ========== =========================================================================================
   ``powkey`` Unique identifier associated with the created observation definition in the obs sequence.
   ========== =========================================================================================

| 

.. container:: routine

   *call get_expected_power(state, location, powkey, val, istatus)*
   ::

      real(r8), intent(in)            :: state
      type(location_type), intent(in) :: location
      integer, intent(in)             :: powkey
      real(r8), intent(out)           :: val
      integer, intent(out)            :: istatus

.. container:: indent1

   Computes the forward observation operator for a power observation. Calls the ``interpolate()`` routine to invoke the
   forward operator code in whatever model this has been compiled with, then raises the result to the specified power
   associated with this powkey.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``state``    | Model state vector (or extended state vector).                                                       |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``location`` | Location of this observation.                                                                        |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``powkey``   | Unique integer key associated with this observation.                                                 |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``val``      | Returned value of forward observation operator.                                                      |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``istatus``  | Returns 0 if forward operator was successfully computed, else returns a positive value. (Negative    |
   |              | values are reserved for system use.)                                                                 |
   +--------------+------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call set_power(power_in, powkey, istatus)*
   ::

      real(r8), intent(in)  :: power_in
      integer,  intent(out) :: powkey
      integer,  intent(out) :: istatus

.. container:: indent1

   Available for use by programs that create observations to set the additional metadata for these observation types.
   This information includes the power to which to raise the state variable. The key that is returned is uniquely
   associated with the definition that has been created and should be set in the obs_def structure by calling
   ``set_obs_def_key()``. This key is different from the main observation key which all observation types have. This key
   is unique to this observation type and is used when reading in the observation sequence to match the corresponding
   metadata with each observation of this type.

   ============ ===================================================================
   ``power_in`` Real value setting the power.
   ``powkey``   Unique integer key associated with the observation being processed.
   ``istatus``  Return code. 0 means success, any other value is an error
   ============ ===================================================================

| 

Namelist
--------

This module has no namelist.

Files
-----

-  NONE

References
----------

#. none

Error codes and conditions
--------------------------

+-------------------------+----------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|         Routine         |                             Message                            |                                                                               Comment                                                                               |
+=========================+================================================================+=====================================================================================================================================================================+
| interactive_1d_integral | Out of space, max_1d_integral_obs limit NNNN (currently 1000). | There is only room for a fixed number of 1d integral observations. The max number is defined by max_1d_integral_obs. Set this to a larger value if more are needed. |
+-------------------------+----------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
