MODULE obs_def_mod
==================

Overview
--------

The DART Fortran90 derived type ``obs_def`` provide an abstraction of the definition of an observation. An observation
sequence ``obs_seq`` at a higher level is composed of observation definitions associated with observed values. For
now, the basic operations required to implement an observation definition are an ability to compute a forward operator
given the model state vector, the ability to read/write the observation definition from/to a file, and a capability to
do a standard input driven interactive definition of the observation definition.

DART makes a distinction between specific ``observation types`` and generic ``observation quantities``. The role of
the various obs_def input files is to define the mapping between the types and quantities, and optionally to provide
type-specific processing routines.

A single obs_def output module is created by the program ``preprocess`` from two kinds of input files. First, a
DEFAULT obs_def module (normally called ``DEFAULT_obs_def_mod.F90`` and documented in this directory) is used as a
template into which the preprocessor incorporates information from zero or more special obs_def modules (such as
``obs_def_1d_state_mod.f90`` or ``obs_def_reanalysis_bufr_mod.f90``, also documented in this directory). If no special
obs_def files are included in the preprocessor namelist, a minimal ``obs_def_mod.f90`` is created which can only
support identity forward observation operators.

New Observation Types
~~~~~~~~~~~~~~~~~~~~~

To add a new observation type which does not fit into any of the already-defined obs_def files, a new file should be
created in the ``obs_def`` directory. These files are usually named according the the pattern
``obs_def_``\ X\ ``_mod.f90``, where the X is either an instrument name, a data source, or a class of observations.
See the existing filenames in that directory for ideas. Then this new filename must be listed in the ``input.nml``
namelist for the model, in the ``&preprocess_nml`` section, in the ``obs_type_files`` variable. This variable is a
string list type which can contain multiple filenames. Running the ``preprocess`` program will then use the contents
of the new file to generate the needed output files for use in linking to the rest of the DART system.

Simple observations
~~~~~~~~~~~~~~~~~~~

If the new observation type can be directly interpolated by a model_mod interpolation routine, and has no additional
observation-specific code for reading, writing, or initializing the observation, then the entire contents of the new
file is:

::

   ! BEGIN DART PREPROCESS TYPE DEFINITIONS
   ! type, quantity, COMMON_CODE
   ! (repeat lines for each type)
   ! END DART PREPROCESS TYPE DEFINITIONS

DART will automatically generate all interface code needed for these new observation types. For example, here is a real
list:

::

   ! BEGIN DART PREPROCESS TYPE DEFINITIONS
   !VELOCITY,                     QTY_VELOCITY,              COMMON_CODE
   !TRACER_CONCENTRATION,         QTY_TRACER_CONCENTRATION,  COMMON_CODE
   !TRACER_SOURCE,                QTY_TRACER_SOURCE,         COMMON_CODE
   !MEAN_SOURCE,                  QTY_MEAN_SOURCE,           COMMON_CODE
   !SOURCE_PHASE,                 QTY_SOURCE_PHASE,          COMMON_CODE
   ! END DART PREPROCESS TYPE DEFINITIONS

The first column is the specific observation ``type`` and should be unique. The second column is the generic observation
``quantity``. The quantities available to DART are defined at compile time by *preprocess* via the option
'quantity_files' in the *preprocess_nml* namelist. The third column must be the keyword ``COMMON_CODE`` which tells the
``preprocess`` program to automatically generate all necessary interface code for this type.

Observations needing special handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For observation types which have observation-specific routines, must interpolate using a combination of other generic
quantities, or require additional observation-specific data to be stored, the following format is used:

::

   ! BEGIN DART PREPROCESS TYPE DEFINITIONS
   ! type, quantity
   ! (repeat lines for each type/quantity pair)
   ! END DART PREPROCESS TYPE DEFINITIONS

DART will need user-supplied interface code for each of the listed types. For example, here is a real list:

::

   ! BEGIN DART PREPROCESS TYPE DEFINITIONS
   ! DOPPLER_RADIAL_VELOCITY, QTY_VELOCITY
   ! RADAR_REFLECTIVITY,      QTY_RADAR_REFLECTIVITY
   ! END DART PREPROCESS TYPE DEFINITIONS

In this case, DART needs additional information for how to process these types. They include code sections delimited by
precisely formatted comments, and possibly module code sections:

#. ::

      ! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
      ! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

   | Any fortran use statements for public subroutines or variables from other modules should be placed between these
     lines, with comment characters in the first column.
   | For example, if the forward operator code includes a module with public routines then a "use" statement like:

   ::

      use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &
                                       interactive_1d_integral, get_expected_1d_integral

   needs to be added to the obs_def_mod so the listed subroutines are available to be called. This would look like:

   ::

      ! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
      ! use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &
      !                                  interactive_1d_integral, get_expected_1d_integral
      ! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

#. ::

      ! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
      ! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

   | These comments must enclose a case statement for each defined type that returns the expected observation value
     based on the current values of the state vector. The code must be in comments, with the comment character in the
     first column.
   | The variables available to be passed to subroutines or used in this section of code are:

   ================ ==========================================
   ``state``        the entire model state vector
   ``state_time``   the time of the state data
   ``ens_index``    the ensemble member number
   ``location``     the observation location
   ``obs_kind_ind`` the index of the specific observation type
   ``obs_time``     the time of the observation
   ``error_val``    the observation error variance
   ================ ==========================================

   | 
   | The routine must fill in the values of these variables:

   =========== ==========================================================
   ``obs_val`` the computed forward operator value
   ``istatus`` return code: 0=ok, >0 is error, <0 reserved for system use
   =========== ==========================================================

   | 
   | To call a model_mod interpolate routine directly, the argument list must match exactly:

   ::

      interpolate(state, location, QTY_xxx, obs_val, istatus)

   This can be useful if the forward operator needs to retrieve values for fields which are typically found in a model
   and then compute a derived value from them.

#. ::

      ! BEGIN DART PREPROCESS READ_OBS_DEF
      ! END DART PREPROCESS READ_OBS_DEF

   | These comments must enclose a case statement for each defined type that reads any additional data associated with a
     single observation. If there is no information beyond that for the basic obs_def type, the case statement must
     still be provided, but the code can simply be ``continue``. The code must be in comments, with the comment
     character in the first column.
   | The variables available to be passed to subroutines or used in this section of code are:

   ============ =====================================================================
   ``ifile``    the open unit number positioned ready to read, read-only
   ``obs_def``  the rest of the obs_def derived type for this obs, read-write
   ``key``      the index observation number in this sequence, read-only
   ``obs_val``  the observation value, if needed. in general should not be changed
   ``is_ascii`` logical to indicate how the file was opened, formatted or unformatted
   ============ =====================================================================

   | 
   | The usual use of this routine is to read in additional metadata per observation and to set the private key in the
     ``obs_def`` to indicate which index to use for this observation to look up the corresponding metadata in arrays or
     derived types. Do not confuse the key in the obs_def with the key argument to this routine; the latter is the
     global observation sequence number for this observation.

#. ::

      ! BEGIN DART PREPROCESS WRITE_OBS_DEF
      ! END DART PREPROCESS WRITE_OBS_DEF

   | These comments must enclose a case statement for each defined type that writes any additional data associated with
     a single observation. If there is no information beyond that for the basic obs_def type, the case statement must
     still be provided, but the code can simply be ``continue``. The code must be in comments, with the comment
     character in the first column.
   | The variables available to be passed to subroutines or used in this section of code are:

   ============ =====================================================================
   ``ifile``    the open unit number positioned ready to write, read-only
   ``obs_def``  the rest of the obs_def derived type for this obs, read-only
   ``key``      the index observation number in this sequence, read-only
   ``is_ascii`` logical to indicate how the file was opened, formatted or unformatted
   ============ =====================================================================

   | 
   | The usual use of this routine is to write the additional metadata for this observation based on the private key in
     the ``obs_def``. Do not confuse this with the key in the subroutine call which is the observation number relative
     to the entire observation sequence file.

#. ::

      ! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
      ! END DART PREPROCESS INTERACTIVE_OBS_DEF

   | These comments must enclose a case statement for each defined type that prompts the user for any additional data
     associated with a single observation. If there is no information beyond that for the basic obs_def type, the case
     statement must still be provided, but the code can simply be ``continue``. The code must be in comments, with the
     comment character in the first column.
   | The variables available to be passed to subroutines or used in this section of code are:

   =========== =============================================================
   ``obs_def`` the rest of the obs_def derived type for this obs, read-write
   ``key``     the index observation number in this sequence, read-only
   =========== =============================================================

   | 
   | The DART code will prompt for the rest of the obs_def values (location, type, value, error) but any additional
     metadata needed by this observation type should be prompted to, and read from, the console (e.g. ``write(*,*)``,
     and ``read(*, *)``). The code will generally set the ``obs_def%key`` value as part of setting the metadata.

#. ::

      ! BEGIN DART PREPROCESS MODULE CODE
      ! END DART PREPROCESS MODULE CODE

   | If the code to process this observation requires module data and/or subroutines, then these comments must surround
     the module definitions. Unlike all the other sections, this comment pair is optional, and if used, the code must
     not be in comments; it will be copied verbatim over to the output file.
   | Generally the code for a forward operator should be defined inside a module, to keep module variables and other
     private subroutines from colliding with unrelated routines and variables in other forward operator files.

It is possible to mix automatic code types and user-supplied code types in the same list. Simply add the COMMON_CODE
keyword on the lines which need no special data or interfaces. For example, here is an extract from the 1d state obs_def
module, where the raw state variable needs only autogenerated code, but the 1d integral has user-supplied processing
code:

::

   ! BEGIN DART PREPROCESS TYPE LIST
   ! RAW_STATE_VARIABLE,    QTY_STATE_VARIABLE, COMMON_CODE
   ! RAW_STATE_1D_INTEGRAL, QTY_1D_INTEGRAL
   ! END DART PREPROCESS TYPE LIST


   ! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
   !   use obs_def_1d_state_mod, only : write_1d_integral, read_1d_integral, &
   !                                    interactive_1d_integral, get_expected_1d_integral
   ! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

   ! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
   !         case(RAW_STATE_1D_INTEGRAL)
   !            call get_expected_1d_integral(state, location, obs_def%key, obs_val, istatus)
   ! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

   ! BEGIN DART PREPROCESS READ_OBS_DEF
   !      case(RAW_STATE_1D_INTEGRAL)
   !         call read_1d_integral(obs_def%key, ifile, fileformat)
   ! END DART PREPROCESS READ_OBS_DEF

   ! BEGIN DART PREPROCESS WRITE_OBS_DEF
   !      case(RAW_STATE_1D_INTEGRAL)
   !         call write_1d_integral(obs_def%key, ifile, fileformat)
   ! END DART PREPROCESS WRITE_OBS_DEF

   ! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
   !      case(RAW_STATE_1D_INTEGRAL)
   !         call interactive_1d_integral(obs_def%key)
   ! END DART PREPROCESS INTERACTIVE_OBS_DEF

   ! BEGIN DART PREPROCESS MODULE CODE
   module obs_def_1d_state_mod

   use        types_mod, only : r8
   use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
   use     location_mod, only : location_type, set_location, get_location
   use  assim_model_mod, only : interpolate
   use   cov_cutoff_mod, only : comp_cov_factor

   implicit none

   public :: write_1d_integral, read_1d_integral, interactive_1d_integral, &
             get_expected_1d_integral

   ...  (module code here)

   end module obs_def_1d_state_mod
   ! END DART PREPROCESS MODULE CODE

| See the :doc:`./obs_def_1d_state_mod` documentation for more details and examples of each section. Also see
  ``obs_def_wind_speed_mod.f90`` for an example of a 3D geophysical forward operator.
| In addition to collecting and managing any additional observation type-specific code, this module provides the
  definition of the obs_def_type derived type, and a collection of subroutines for creating, accessing, and updating
  this type. The remainder of this document describes the subroutines provided by this module.

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod (depends on model choice)
   time_manager_mod
   assim_model_mod
   obs_kind_mod
   Other special obs_def_kind modules as required

Public interfaces
-----------------

========================= ==========================
*use obs_def_mod, only :* obs_def_type
\                         init_obs_def
\                         get_obs_def_location
\                         get_obs_def_type_of_obs
\                         get_obs_def_time
\                         get_obs_def_error_variance
\                         get_obs_def_key
\                         set_obs_def_location
\                         set_obs_def_type_of_obs
\                         set_obs_def_time
\                         set_obs_def_error_variance
\                         set_obs_def_key
\                         interactive_obs_def
\                         write_obs_def
\                         read_obs_def
\                         get_expected_obs_from_def
\                         destroy_obs_def
\                         copy_obs_def
\                         assignment(=)
\                         get_name_for_type_of_obs
========================= ==========================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   ::

      type obs_def_type
         private
         type(location_type)  :: location
         integer              :: kind
         type(time_type)      :: time
         real(r8)             :: error_variance
         integer              :: key
      end type obs_def_type

.. container:: indent1

   Models all that is known about an observation except for actual values. Includes a location, type, time and error
   variance.

   ============== ========================================================
   Component      Description
   ============== ========================================================
   location       Location of the observation.
   kind           Despite the name, the specific type of the observation.
   time           Time of the observation.
   error_variance Error variance of the observation.
   key            Unique identifier for observations of a particular type.
   ============== ========================================================

| 

.. container:: routine

   *call init_obs_def(obs_def, location, kind, time, error_variance)*
   ::

      type(obs_def_type),  intent(out) :: obs_def
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: kind
      type(time_type),     intent(in)  :: time
      real(r8),            intent(in)  :: error_variance

.. container:: indent1

   Creates an obs_def type with location, type, time and error_variance specified.

   ================== ==================================
   ``obs_def``        The obs_def that is created
   ``location``       Location for this obs_def
   ``kind``           Observation type for obs_def
   ``time``           Time for obs_def
   ``error_variance`` Error variance of this observation
   ================== ==================================

| 

.. container:: routine

   *call copy_obs_def(obs_def1, obs_def2)*
   ::

      type(obs_def_type), intent(out) :: obs_def1
      type(obs_def_type), intent(in)  :: obs_def2

.. container:: indent1

   Copies obs_def2 to obs_def1, overloaded as assignment (=).

   ============ =========================
   ``obs_def1`` obs_def to be copied into
   ``obs_def2`` obs_def to be copied from
   ============ =========================

| 

.. container:: routine

   *var = get_obs_def_key(obs_def)*
   ::

      integer                        :: get_obs_def_key
      type(obs_def_type), intent(in) :: obs_def

.. container:: indent1

   Returns key from an observation definition.

   =========== ===========================
   ``var``     Returns key from an obs_def
   ``obs_def`` An obs_def
   =========== ===========================

| 

.. container:: routine

   *var = get_obs_def_error_variance(obs_def)*
   ::

      real(r8)                       :: get_obs_def_error_variance
      type(obs_def_type), intent(in) :: obs_def

.. container:: indent1

   Returns error variance from an observation definition.

   =========== ==============================
   ``var``     Error variance from an obs_def
   ``obs_def`` An obs_def
   =========== ==============================

| 

.. container:: routine

   *var = get_obs_def_location(obs_def)*
   ::

      type(location_type)              :: get_obs_def_location
      type(obs_def_type), intent(in)   :: obs_def

.. container:: indent1

   Returns the location from an observation definition.

   =========== ================================
   ``var``     Returns location from an obs_def
   ``obs_def`` An obs_def
   =========== ================================

| 

.. container:: routine

   *var = get_obs_def_type_of_obs(obs_def)*
   ::

      integer                         :: get_obs_def_type_of_obs
      type(obs_def_type),  intent(in) :: obs_def

.. container:: indent1

   Returns an observation type from an observation definition.

   =========== ============================================
   ``var``     Returns the observation type from an obs_def
   ``obs_def`` An obs_def
   =========== ============================================

| 

.. container:: routine

   *var = get_obs_def_time(obs_def)*
   ::

      type(time_type)                :: get_obs_def_time
      type(obs_def_type), intent(in) :: obs_def

.. container:: indent1

   Returns time from an observation definition.

   =========== ============================
   ``var``     Returns time from an obs_def
   ``obs_def`` An obs_def
   =========== ============================

| 

.. container:: routine

   *obs_name = get_name_for_type_of_obs(obs_kind_ind)*
   ::

      character(len = 32)            :: get_name_for_type_of_obs
      integer, intent(in)            :: obs_kind_ind

.. container:: indent1

   Returns an observation name from an observation type.

   ================ =====================================
   ``var``          Returns name from an observation type
   ``obs_kind_ind`` An observation type
   ================ =====================================

| 

.. container:: routine

   *call set_obs_def_location(obs_def, location)*
   ::

      type(obs_def_type),  intent(inout) :: obs_def
      type(location_type), intent(in)    :: location

.. container:: indent1

   Set the location in an observation definition.

   ============ ==========
   ``obs_def``  An obs_def
   ``location`` A location
   ============ ==========

| 

.. container:: routine

   *call set_obs_def_error_variance(obs_def, error_variance)*
   ::

      type(obs_def_type), intent(inout) :: obs_def
      real(r8), intent(in)              :: error_variance

.. container:: indent1

   Set error variance for an observation definition.

   ================== ==============
   ``obs_def``        An obs_def
   ``error_variance`` Error variance
   ================== ==============

| 

.. container:: routine

   *call set_obs_def_key(obs_def, key)*
   ::

      type(obs_def_type), intent(inout) :: obs_def
      integer,            intent(in)    :: key

.. container:: indent1

   Set the key for an observation definition.

   =========== ======================================
   ``obs_def`` An obs_def
   ``key``     Unique identifier for this observation
   =========== ======================================

| 

.. container:: routine

   *call set_obs_def_type_of_obs(obs_def, kind)*
   ::

      type(obs_def_type), intent(inout) :: obs_def
      integer,            intent(in)    :: kind

.. container:: indent1

   Set the type of observation in an observation definition.

   =========== ===========================
   ``obs_def`` An obs_def
   ``kind``    An integer observation type
   =========== ===========================

| 

.. container:: routine

   *call set_obs_def_time(obs_def, time)*
   ::

      type(obs_def_type), intent(inout) :: obs_def
      type(time_type), intent(in)       :: time

.. container:: indent1

   Sets time for an observation definition.

   =========== ===========
   ``obs_def`` An obs_def
   ``time``    Time to set
   =========== ===========

| 

.. container:: routine

   *call get_expected_obs_from_def(key, obs_def, obs_kind_ind, ens_index, state, state_time, obs_val, istatus,
   assimilate_this_ob, evaluate_this_ob)*
   ::

      integer,            intent(in)  :: key
      type(obs_def_type), intent(in)  :: obs_def
      integer,            intent(in)  :: obs_kind_ind
      integer,            intent(in)  :: ens_index
      real(r8),           intent(in)  :: state(:)
      type(time_type),    intent(in)  :: state_time
      real(r8),           intent(out) :: obs_val
      integer,            intent(out) :: istatus
      logical,            intent(out) :: assimilate_this_ob
      logical,            intent(out) :: evaluate_this_ob

.. container:: indent1

   Compute the observation (forward) operator for a particular obs definition.

   +------------------------+--------------------------------------------------------------------------------------------+
   | ``key``                | descriptor for observation type                                                            |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``obs_def``            | The input obs_def                                                                          |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``obs_kind_ind``       | The obs type                                                                               |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``ens_index``          | The ensemble member number of this state vector                                            |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``state``              | Model state vector                                                                         |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``state_time``         | Time of the data in the model state vector                                                 |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``istatus``            | Returned integer describing problems with applying forward operator (0 == OK, >0 == error, |
   |                        | <0 reserved for sys use).                                                                  |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``assimilate_this_ob`` | Indicates whether to assimilate this obs or not                                            |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``evaluate_this_ob``   | Indicates whether to evaluate this obs or not                                              |
   +------------------------+--------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call read_obs_def(ifile, obs_def, key, obs_val [,fform])*
   ::

      integer,                    intent(in)    :: ifile
      type(obs_def_type),         intent(inout) :: obs_def
      integer,                    intent(in)    :: key
      real(r8),                   intent(inout) :: obs_val
      character(len=*), optional, intent(in)    :: fform

.. container:: indent1

   Reads an obs_def from file open on channel ifile. Uses format specified in fform or FORMATTED if fform is not
   present.

   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``ifile``   | File unit open to output file                                                                         |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``obs_def`` | Observation definition to be read                                                                     |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``key``     | Present if unique identifier key is needed by some obs type. Unused by default code.                  |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``obs_val`` | Present if needed to perform operations based on value. Unused by default code.                       |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``fform``   | File format specifier: FORMATTED or UNFORMATTED; default FORMATTED (FORMATTED in this case is the     |
   |             | human readable/text option as opposed to UNFORMATTED which is binary.)                                |
   +-------------+-------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call interactive_obs_def(obs_def, key)*
   ::

      type(obs_def_type), intent(inout) :: obs_def
      integer,            intent(in)    :: key

.. container:: indent1

   Creates an obs_def via input from standard in.

   =========== ====================================================================================
   ``obs_def`` An obs_def to be created
   ``key``     Present if unique identifier key is needed by some obs type. Unused by default code.
   =========== ====================================================================================

| 

.. container:: routine

   *call write_obs_def(ifile, obs_def, key [,fform])*
   ::

      integer,                    intent(in) :: ifile
      type(obs_def_type),         intent(in) :: obs_def
      integer,                    intent(in) :: key
      character(len=*), optional, intent(in) :: fform

.. container:: indent1

   Writes an obs_def to file open on channel ifile. Uses format specified in fform or FORMATTED if fform is not present.

   =========== ====================================================================================
   ``ifile``   File unit open to output file
   ``obs_def`` Observation definition to be written
   ``key``     Present if unique identifier key is needed by some obs type. Unused by default code.
   ``fform``   File format specifier: FORMATTED or UNFORMATTED; default FORMATTED
   =========== ====================================================================================

| 

.. container:: routine

   *call destroy_obs_def(obs_def)*
   ::

      type(obs_def_type), intent(inout) :: obs_def

.. container:: indent1

   Releases all storage associated with an obs_def and its subcomponents.

   =========== ==========================
   ``obs_def`` An obs_def to be released.
   =========== ==========================

| 

Files
-----

-  The read_obs_def() and write_obs_def() routines are passed an already-opened file channel/descriptor and read to or
   write from it.

References
----------

-  none


Error codes and conditions
--------------------------

+---------------------------+----------------------------------------------------------+---------------------------------------------------------------------------------------------+
|          Routine          |                          Message                         |                                           Comment                                           |
+===========================+==========================================================+=============================================================================================+
| get_expected_obs_from_def | Attempt to evaluate undefined observation type           | An observation type for which no forward operator has been defined is an error.             |
+---------------------------+----------------------------------------------------------+---------------------------------------------------------------------------------------------+
| read_obs_def              | Expected header "obdef" in input file                    | The format of the input file is not consistent.                                             |
+---------------------------+----------------------------------------------------------+---------------------------------------------------------------------------------------------+
| read_obs_def              | Expected kind header "kind " in input file               | The format of the input file is not consistent.                                             |
+---------------------------+----------------------------------------------------------+---------------------------------------------------------------------------------------------+
| read_obs_def              | Attempt to read for undefined obs_kind index             | Reading for an observation type for which no forward operator has been defined is an error. |
+---------------------------+----------------------------------------------------------+---------------------------------------------------------------------------------------------+
| write_obs_def             | Attempt to write for undefined obs_kind index            | Writing for an observation type for which no forward operator has been defined is an error. |
+---------------------------+----------------------------------------------------------+---------------------------------------------------------------------------------------------+
| interactive_obs_def       | Attempt to interactively create undefined obs_kind index | Creating an observation type for which no forward operator has been defined is an error.    |
+---------------------------+----------------------------------------------------------+---------------------------------------------------------------------------------------------+


Private components
------------------

N/A
