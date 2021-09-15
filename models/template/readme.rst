MODULE model_mod
================

Overview
--------

| Every model that is DART compliant must provide an interface as documented here. The file
  ``models/template/model_mod.f90`` provides the fortran interfaces for a minimal implementation meeting these
  requirements. When adding a new model to DART you can either start by modifying a ``model_mod.f90`` file from a
  similar model already in DART or start with the template file. Either way, the supplied interface must match these
  descriptions exactly; no details of the underlying model can impact the interface.
| Several of the routines listed below are allowed to be a NULL INTERFACE. This means the subroutine or function name
  must exist in this file, but it is ok if it contains no executable code.
| A few of the routines listed below are allowed to be a PASS-THROUGH INTERFACE. This means the subroutine or function
  name can be listed on the 'use' line from the ``location_mod``, and no subroutine or function with that name is
  supplied in this file. Alternatively, this file can provide an implementation which calls the underlying routines from
  the ``location_mod`` and then alters or augments the results based on model-specific requirements.
| The system comes with several types of location modules for computing distances appropriately. Two of the ones most
  commonly used are for data in a 1D system and for data in a 3D spherical coordinate system. Make the selection by
  listing the appropriate choice from ``location/*/location_mod.f90`` in the corresponding ``path_names_*`` file at
  compilation time.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &model_nml 
    /

| 

Models are free to include a model namelist which can be read when ``static_init_model`` is called. A good example can
be found in the lorenz_96 ``model_mod.f90``.

Other modules used
------------------

::

   types_mod
   time_manager_mod
   location_mod (multiple choices here)
   utilities_mod
   POSSIBLY MANY OTHERS DEPENDING ON MODEL DETAILS

Public interfaces
-----------------

======================= ===================================
*use model_mod, only :* get_model_size
\                       adv_1step
\                       get_state_meta_data
\                       model_interpolate
\                       shortest_time_between_assimilations
\                       static_init_model
\                       init_time
\                       init_conditions
\                       nc_write_model_atts
\                       nc_write_model_vars
\                       pert_model_copies
\                       get_close_obs
\                       get_close_state
\                       convert_vertical_obs
\                       convert_vertical_state
\                       read_model_time
\                       write_model_time
\                       end_model
======================= ===================================

A namelist interface ``&model_nml`` may be defined by the module, in which case it will be read from file ``input.nml``.
The details of the namelist are always model-specific (there are no generic namelist values).

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *model_size = get_model_size( )*
   ::

      integer(i8) :: get_model_size

.. container:: indent1

   Returns the length of the model state vector. Required.

   ============== =====================================
   ``model_size`` The length of the model state vector.
   ============== =====================================

| 

.. container:: routine

   *call adv_1step(x, time)*
   ::

      real(r8), dimension(:), intent(inout) :: x
      type(time_type),        intent(in)    :: time

.. container:: indent1

   Does a single timestep advance of the model. The input value of the vector x is the starting condition and x must be
   updated to reflect the changed state after a timestep. The time argument is intent in and is used for models that
   need to know the date/time to compute a timestep, for instance for radiation computations. This interface is only
   called if the namelist parameter async is set to 0 in ``perfect_model_obs`` or ``filter`` or if the program
   ``integrate_model`` is to be used to advance the model state as a separate executable. If one of these options is not
   going to be used (the model will *only* be advanced as a separate model-specific executable), this can be a NULL
   INTERFACE. (The subroutine name must still exist, but it can contain no code and it will not be called.)

   ======== ==================================
   ``x``    State vector of length model_size.
   ``time`` Current time of the model state.
   ======== ==================================

| 

.. container:: routine

   *call get_state_meta_data (index_in, location, [, var_type] )*
   ::

      integer,             intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer, optional,   intent(out) ::  var_type 

.. container:: indent1

   Given an integer index into the state vector, returns the associated location. An optional argument returns the
   generic quantity of this item, e.g. QTY_TEMPERATURE, QTY_DENSITY, QTY_SALINITY, QTY_U_WIND_COMPONENT. This interface
   is required to be functional for all applications.

   ============ ===================================================================
   ``index_in`` Index of state vector element about which information is requested.
   ``location`` The location of state variable element.
   *var_type*   The generic quantity of the state variable element.
   ============ ===================================================================

| 

.. container:: routine

   *call model_interpolate(state_handle, ens_size, location, obs_quantity, expected_obs, istatus)*
   ::

      type(ensemble_type),    intent(in)  :: state_handle
      integer,                intent(in)  :: ens_size
      type(location_type),    intent(in)  :: location
      integer,                intent(in)  :: obs_quantity
      real(r8),               intent(out) :: expected_obs(ens_size)
      integer,                intent(out) :: istatus(ens_size)

.. container:: indent1

   Given a handle containing information for a state vector, an ensemble size, a location, and a model state variable
   quantity interpolates the state variable field to that location and returns an ensemble-sized array of values in
   ``expected_obs(:)``. The ``istatus(:)`` array should be 0 for successful ensemble members and a positive value for
   failures. The ``obs_quantity`` variable is one of the quantity (QTY) parameters defined in the
   :doc:`../../assimilation_code/modules/observations/obs_kind_mod` file and defines the quantity to interpolate. In
   low-order models that have no notion of kinds of variables this argument may be ignored. For applications in which
   only perfect model experiments with identity observations (i.e. only the value of a particular state variable is
   observed), this can be a NULL INTERFACE. Otherwise it is required (which is the most common case).

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state structure containing information about the state vector about which      |
   |                  | information is requested.                                                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``ens_size``     | The ensemble size.                                                                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``location``     | Location to which to interpolate.                                                                |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``obs_quantity`` | Quantity of state field to be interpolated.                                                      |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``expected_obs`` | The interpolated values from the model.                                                          |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``istatus``      | Integer values return 0 for success. Other positive values can be defined for various failures.  |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = shortest_time_between_assimilations()*
   ::

      type(time_type) :: shortest_time_between_assimilations

.. container:: indent1

   Returns the smallest increment in time that the model is capable of advancing the state in a given implementation.
   The actual value may be set by the model_mod namelist (depends on the model). This interface is required for all
   applications.

   ======= ===================================
   ``var`` Smallest advance time of the model.
   ======= ===================================

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   Called to do one time initialization of the model. As examples, might define information about the model size or
   model timestep, read in grid information, read a namelist, set options, etc. In models that require pre-computed
   static data, for instance spherical harmonic weights, these would also be computed here. Can be a NULL INTERFACE for
   the simplest models.

| 

.. container:: routine

   *call init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   Companion interface to init_conditions. Returns a time that is somehow appropriate for starting up a long integration
   of the model. At present, this is only used if the ``perfect_model_obs`` namelist parameter
   ``read_input_state_from_file = .false.`` If this option should not be used in ``perfect_model_obs``, calling this
   routine should issue a fatal error.

   ======== ===================
   ``time`` Initial model time.
   ======== ===================

| 

.. container:: routine

   *call init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   Returns a model state vector, x, that is some sort of appropriate initial condition for starting up a long
   integration of the model. At present, this is only used if the ``perfect_model_obs`` namelist parameter
   ``read_input_state_from_file = .false.`` If this option should not be used in ``perfect_model_obs``, calling this
   routine should issue a fatal error.

   ===== ====================================
   ``x`` Initial conditions for state vector.
   ===== ====================================

| 

.. container:: routine

   *call nc_write_model_atts(ncFileID, domain_id)*
   ::

      integer, intent(in) :: ncFileID
      integer, intent(in) :: domain_id

.. container:: indent1

   | This routine writes the model-specific attributes to netCDF files that DART creates. This includes coordinate
     variables and any metadata, but NOT the actual model state vector. ``models/template/model_mod.f90`` contains code
     that can be used for any model as-is.
   | The typical sequence for adding new dimensions, variables, attributes:

   ::

      NF90_OPEN             ! open existing netCDF dataset               
         NF90_redef         ! put into define mode                       
         NF90_def_dim       ! define additional dimensions (if any)     
         NF90_def_var       ! define variables: from name, kind, and dims
         NF90_put_att       ! assign attribute values                    
      NF90_ENDDEF           ! end definitions: leave define mode         
         NF90_put_var       ! provide values for variable                
      NF90_CLOSE            ! close: save updated netCDF dataset        

   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``ncFileID``  | Integer file descriptor to previously-opened netCDF file.                                           |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``domain_id`` | integer describing the domain (which can be a nesting level, a component model ...) Models with     |
   |               | nested grids are decomposed into 'domains' in DART. The concept is extended to refer to 'coupled'   |
   |               | models where one model component may be the atmosphere, another component may be the ocean, or      |
   |               | land, or ionosphere ... these would be referenced as different domains.                             |
   +---------------+-----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call nc_write_model_vars(ncFileID, domain_id, state_ens_handle [, memberindex] [, timeindex])*
   ::

      integer,             intent(in) :: ncFileID
      integer,             intent(in) :: domain_id
      type(ensemble_type), intent(in) :: state_ens_handle
      integer, optional,   intent(in) :: memberindex
      integer, optional,   intent(in) :: timeindex

.. container:: indent1

   | This routine may be used to write the model-specific state vector (data) to a netCDF file. Only used if
     ``model_mod_writes_state_variables = .true.``
   | Typical sequence for adding new dimensions,variables,attributes:

   ::

      NF90_OPEN             ! open existing netCDF dataset               
         NF90_redef         ! put into define mode                       
         NF90_def_dim       ! define additional dimensions (if any)      
         NF90_def_var       ! define variables: from name, kind, and dims
         NF90_put_att       ! assign attribute values                    
      NF90_ENDDEF           ! end definitions: leave define mode         
         NF90_put_var       ! provide values for variable                
      NF90_CLOSE            ! close: save updated netCDF dataset         

   +----------------------+----------------------------------------------------------------------------------------------+
   | ``ncFileID``         | file descriptor to previously-opened netCDF file.                                            |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``domain_id``        | integer describing the domain (which can be a nesting level, a component model ...)          |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``state_ens_handle`` | The handle to the state structure containing information about the state vector about which  |
   |                      | information is requested.                                                                    |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``memberindex``      | Integer index of ensemble member to be written.                                              |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``timeindex``        | The timestep counter for the given state.                                                    |
   +----------------------+----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)*
   ::

      type(ensemble_type), intent(inout) :: state_ens_handle
      integer,             intent(in)    :: ens_size
      real(r8),            intent(in)    :: pert_amp
      logical,             intent(out)   :: interf_provided

.. container:: indent1

   Given an ensemble handle, the ensemble size, and a perturbation amplitude; perturb the ensemble. Used to generate
   initial conditions for spinning up ensembles. If the ``model_mod`` does not want to do this, instead allowing the
   default algorithms in ``filter`` to take effect, ``interf_provided =&nbps;.false.`` and the routine can be trivial.
   Otherwise, ``interf_provided`` must be returned as ``.true.``

   +----------------------+----------------------------------------------------------------------------------------------+
   | ``state_ens_handle`` | The handle containing an ensemble of state vectors to be perturbed.                          |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``ens_size``         | The number of ensemble members to perturb.                                                   |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``pert_amp``         | the amplitude of the perturbations. The interpretation is based on the model-specific        |
   |                      | implementation.                                                                              |
   +----------------------+----------------------------------------------------------------------------------------------+
   | ``interf_provided``  | Returns false if model_mod cannot do this, else true.                                        |
   +----------------------+----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, num_close, close_ind [, dist] [,
   state_handle)*
   ::

      type(get_close_type),          intent(in)  :: gc
      type(location_type),           intent(in)  :: base_loc
      integer,                       intent(in)  :: base_type
      type(location_type),           intent(in)  :: locs(:)
      integer,                       intent(in)  :: loc_qtys(:)
      integer,                       intent(in)  :: loc_types(:)
      integer,                       intent(out) :: num_close
      integer,                       intent(out) :: close_ind(:)
      real(r8),            optional, intent(out) :: dist(:)
      type(ensemble_type), optional, intent(in)  :: state_handle

.. container:: indent1

   | Given a location and quantity, compute the distances to all other locations in the ``obs`` list. The return values
     are the number of items which are within maxdist of the base, the index numbers in the original obs list, and
     optionally the distances. The ``gc`` contains precomputed information to speed the computations.
   | In general this is a PASS-THROUGH ROUTINE. It is listed on the use line for the locations_mod, and in the public
     list for this module, but has no subroutine declaration and no other code in this module:

   ::

      use location_mod, only: get_close_obs

      public :: get_close_obs

   However, if the model needs to alter the values or wants to supply an alternative implementation it can intercept the
   call like so:

   ::

      use location_mod, only: &
              lm_get_close_obs => get_close_obs
              
      public :: get_close_obs

   In this case a local ``get_close_obs()`` routine must be supplied. To call the original code in the location module
   use:

   ::

      call lm_get_close_obs(gc, base_loc, ...)

   | This subroutine will be called after ``get_close_maxdist_init`` and ``get_close_obs_init``.
   | In most cases the PASS-THROUGH ROUTINE will be used, but some models need to alter the actual distances depending
     on the observation or state vector kind, or based on the observation or state vector location. It is reasonable in
     this case to leave ``get_close_maxdist_init()`` and ``get_close_obs_init()`` as pass-through routines and intercept
     only ``get_close_obs()``. The local ``get_close_obs()`` can first call the location mod routine and let it return a
     list of values, and then inspect the list and alter or remove any entries as needed. See the CAM and WRF model_mod
     files for examples of this use.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``gc``           | The get_close_type which stores precomputed information about the locations to speed up          |
   |                  | searching                                                                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``base_loc``     | Reference location. The distances will be computed between this location and every other         |
   |                  | location in the obs list                                                                         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``base_type``    | The DART quantity at the ``base_loc``                                                            |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``locs(:)``      | Compute the distance between the ``base_loc`` and each of the locations in this list             |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_qtys(:)``  | The corresponding quantity of each item in the ``locs`` list                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_types(:)`` | The corresponding type of each item in the ``locs`` list. This is not available in the default   |
   |                  | implementation but may be used in custom implementations.                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num_close``    | The number of items from the ``locs`` list which are within maxdist of the base location         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``close_ind(:)`` | The list of index numbers from the ``locs`` list which are within maxdist of the base location   |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``dist(:)``      | If present, return the distance between each entry in the close_ind list and the base location.  |
   |                  | If not present, all items in the obs list which are closer than maxdist will be added to the     |
   |                  | list but the overhead of computing the exact distances will be skipped.                          |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state structure containing information about the state vector about which      |
   |                  | information is requested.                                                                        |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_close_state(gc, base_loc, base_type, state_loc, state_qtys, state_indx, num_close, close_ind [, dist,
   state_handle])*
   ::

      type(get_close_type),          intent(in)    :: gc
      type(location_type),           intent(inout) :: base_loc
      integer,                       intent(in)    :: base_type
      type(location_type),           intent(inout) :: state_loc(:)
      integer,                       intent(in)    :: state_qtys(:)
      integer(i8),                   intent(in)    :: state_indx(:)
      integer,                       intent(out)   :: num_close
      integer,                       intent(out)   :: close_ind(:)
      real(r8),            optional, intent(out)   :: dist(:)
      type(ensemble_type), optional, intent(in)    :: state_handle

.. container:: indent1

   | Given a location and quantity, compute the distances to all other locations in the ``state_loc`` list. The return
     values are the number of items which are within maxdist of the base, the index numbers in the original state_loc
     list, and optionally the distances. The ``gc`` contains precomputed information to speed the computations.
   | In general this is a PASS-THROUGH ROUTINE. It is listed on the use line for the locations_mod, and in the public
     list for this module, but has no subroutine declaration and no other code in this module:

   ::

      use location_mod, only: get_close_state

      public :: get_close_state

   However, if the model needs to alter the values or wants to supply an alternative implementation it can intercept the
   call like so:

   ::

      use location_mod, only: &
              lm_get_close_state => get_close_state
              
      public :: get_close_state

   In this case a local ``get_close_state()`` routine must be supplied. To call the original code in the location module
   use:

   ::

      call loc_get_close_state(gc, base_loc, ...)

   | This subroutine will be called after ``get_close_maxdist_init`` and ``get_close_state_init``.
   | In most cases the PASS-THROUGH ROUTINE will be used, but some models need to alter the actual distances depending
     on the observation or state vector kind, or based on the observation or state vector location. It is reasonable in
     this case to leave ``get_close_maxdist_init()`` and ``get_close_state_init()`` as pass-through routines and
     intercept only ``get_close_state()``. The local ``get_close_state()`` can first call the location mod routine and
     let it return a list of values, and then inspect the list and alter or remove any entries as needed. See the CAM
     and WRF model_mod files for examples of this use.

   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``gc``            | The get_close_type which stores precomputed information about the locations to speed up         |
   |                   | searching                                                                                       |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``base_loc``      | Reference location. The distances will be computed between this location and every other        |
   |                   | location in the list                                                                            |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``base_type``     | The DART quantity at the ``base_loc``                                                           |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_loc(:)``  | Compute the distance between the ``base_loc`` and each of the locations in this list            |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_qtys(:)`` | The corresponding quantity of each item in the ``state_loc`` list                               |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_indx(:)`` | The corresponding DART index of each item in the ``state_loc`` list. This is not available in   |
   |                   | the default implementation but may be used in custom implementations.                           |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``num_close``     | The number of items from the ``state_loc`` list which are within maxdist of the base location   |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``close_ind(:)``  | The list of index numbers from the ``state_loc`` list which are within maxdist of the base      |
   |                   | location                                                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``dist(:)``       | If present, return the distance between each entry in the ``close_ind`` list and the base       |
   |                   | location. If not present, all items in the ``state_loc`` list which are closer than maxdist     |
   |                   | will be added to the list but the overhead of computing the exact distances will be skipped.    |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``state_handle``  | The handle to the state structure containing information about the state vector about which     |
   |                   | information is requested.                                                                       |
   +-------------------+-------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, which_vert, status)*
   ::

      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: num
      type(location_type), intent(in)  :: locs(:)
      integer,             intent(in)  :: loc_qtys(:)
      integer,             intent(in)  :: loc_types(:)
      integer,             intent(in)  :: which_vert
      integer,             intent(out) :: status(:)

.. container:: indent1

   Converts the observations to the desired vertical localization coordinate system. Some models (toy models with no
   'real' observations) will not need this. Most (real) models have observations in one or more coordinate systems
   (pressure, height) and the model is generally represented in only one coordinate system. To be able to interpolate
   the model state to the observation location, or to compute the true distance between the state and the observation,
   it is necessary to convert everything to a single coodinate system.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state.                                                                         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num``          | the number of observation locations                                                              |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``locs``         | the array of observation locations                                                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_qtys``     | the array of observation quantities.                                                             |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_types``    | the array of observation types.                                                                  |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``which_vert``   | the desired vertical coordinate system. There is a table in the ``location_mod.f90`` that        |
   |                  | relates integers to vertical coordinate systems.                                                 |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``status``       | Success or failure of the vertical conversion. If ``istatus = 0``, the conversion was a success. |
   |                  | Any other value is a failure.                                                                    |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call convert_vertical_state(state_handle, num, locs, loc_qtys, loc_types, which_vert, status)*
   ::

      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: num
      type(location_type), intent(in)  :: locs(:)
      integer,             intent(in)  :: loc_qtys(:)
      integer(i8),         intent(in)  :: loc_indx(:)
      integer,             intent(in)  :: which_vert
      integer,             intent(out) :: status(:)

.. container:: indent1

   Converts the state to the desired vertical localization coordinate system. Some models (toy models with no 'real'
   observations) will not need this. To compute the true distance between the state and the observation, it is necessary
   to convert everything to a single coodinate system.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The handle to the state.                                                                         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num``          | the number of state locations                                                                    |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``locs``         | the array of state locations                                                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_qtys``     | the array of state quantities.                                                                   |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``loc_indx``     | the array of state indices.                                                                      |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``which_vert``   | the desired vertical coordinate system. There is a table in the ``location_mod.f90`` that        |
   |                  | relates integers to vertical coordinate systems.                                                 |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``status``       | Success or failure of the vertical conversion. If ``istatus = 0``, the conversion was a success. |
   |                  | Any other value is a failure.                                                                    |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *model_time = read_model_time(filename)*
   ::

      character(len=*), intent(in) :: filename
      type(time_type)              :: model_time

.. container:: indent1

   Reads the valid time of the model state in a netCDF file. There is a default routine in
   ``assimilation_code/modules/io/dart_time_io_mod.f90`` that can be used as a pass-through. That routine will read the
   **last** timestep of a 'time' variable - which is the same strategy used for reading netCDF files that have multiple
   timesteps in them. If your model has some other representation of time (i.e. it does not use a netCDF variable named
   'time') - you will have to write this routine.

   ============= ====================================
   ``ncid``      handle to an open netCDF file
   ``dart_time`` The current time of the model state.
   ============= ====================================

| 

.. container:: routine

   *call write_model_time(ncid, dart_time)*
   ::

      integer,          intent(in) :: ncid
      type(time_type),  intent(in) :: dart_time

.. container:: indent1

   Writes the assimilation time to a netCDF file. There is a default routine in
   ``assimilation_code/modules/io/dart_time_io_mod.f90`` that can be used as a pass-through. If your model has some
   other representation of time (i.e. it does not use a netCDF variable named 'time') - you will have to write this
   routine.

   ============= ====================================
   ``ncid``      handle to an open netCDF file
   ``dart_time`` The current time of the model state.
   ============= ====================================

| 

.. container:: routine

   *call end_model()*

.. container:: indent1

   Does any shutdown and clean-up needed for model. Can be a NULL INTERFACE if the model has no need to clean up
   storage, etc.

Files
-----

-  Models are free to read and write files as they see fit.

References
----------

#. none

Private components
------------------

N/A
