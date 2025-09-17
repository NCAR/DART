.. index:: required model_mod routines, 
           required model mod routines,
           required routines
.. _Required model_mod routines:

Required model_mod routines
===========================

There are 18 Fortran subroutines required to interface a model with DART. 
Place these routines in your ``model_mod.f90`` within  ``DART/models/{your_model}``.
Many routines have sensible defaults in default_model_mod. 
For example, to use the default pert_model_copies routine add the following line in your model_mod.f90:

.. code-block:: fortran

   use default_model_mod, only : pert_model_copies

As in all Fortran programs, multiple routines can be listed after the colon.

.. important::

   Always use the **only** keyword to import just the routines you need from default_model_mod.


Most models will still require custom implementations for some or all routines. 
The table below lists all 18 routines and their defaults. 
See :doc:`suggestions-for-a-simple-model` (advanced by DART) or 
:doc:`suggestions-for-a-complex-model` (advanced externally) for guidance.

.. list-table:: Required model_mod routines
    :header-rows: 1
    :widths: 10 40 35

    * - Routine
      - Purpose
      - Default behavior
    * - **init_time**
      - Set the *initial time* if not read from the restart file.
      - Sets the initial time to 0 days, 0 seconds. To fail use ``init_time => fail_init_time``.
    * - **init_conditions**
      - For a “cold start” fill in an empty state vector with *initial conditions*. Many models cannot just make up values from thin air and thus choose to fail when this is requested.
      - Sets the initial state to 0. To fail use ``init_conditions => fail_init_conditions``.
    * - **get_model_size**
      - Return the *number of items in the state vector*.
      - Returns 1; i.e. there is only one item in the state.
    * - **static_init_model**
      - *Initialize* DART with information about the model that will be used by the remaining ``model_mod`` routines. The procedure for doing this will depend on how complex the model is.
      - Does nothing.
    * - **get_state_meta_data**
      - Takes an index into the state vector and returns the *location* corresponding to that value and optionally the *variable quantity*.
      - Sets a missing location and the default variable quantity.
    * - **end_model**
      - *Deallocate* any arrays allocated in **static_init_model**.
      - Does nothing.
    * - **adv_1step**
      - If possible, *advance the model* state from one time to another. Complex models will be unable to implement this method and should fail.
      - Call the error handler with the message “unable to advance model”.
    * - **shortest_time_ between_assimilations**
      - Return a namelist or a fixed value for the *minimum model advance time* between assimilations. Note that complex models will handle advancing the time externally.
      - Returns a time period of 1 day.
    * - **model_interpolate**
      - *Interpolate* a requested quantity to the given location to get an array of expected values for all ensemble members. \ *NOTE*: this is often the most time consuming method to implement.
      - Fail and set the expected observation to “missing.”
    * - **nc_write_model_atts**
      - Add any *additional information* to the netCDF output diagnostic files. *NOTE*: the state will already be output by other routines, so this method should not create or write the state variables.
      - Does nothing.
    * - **read_model_time**
      - *Read* the model time from a netCDF file.
      - Attempt to read the “time” variable from a state file in an intelligent way.
    * - **write_model_time**
      - *Write* the model time to a netCDF file.
      - Write the “time” variable from the file according to the DART calendar.
    * - **pert_model_copies**
      - *Perturb* a state vector in order to create an ensemble.
      - Add Gaussian noise with a specified amplitude to all parts of the state vector.
    * - **convert_vertical_obs**
      - Some 3D models have multiple vertical coordinates (e.g. pressure, height, or model level); this method *converts observations* between different vertical coordinate systems.
      - Do no conversion.
    * - **convert_vertical_state**
      - Some 3D models have multiple vertical coordinates (e.g. pressure, height, or model level); this method *converts state* between different vertical coordinate systems.
      - Do no conversion.
    * - **get_close_obs**
      - Calculate *which observations are “close”* to a given location and, optionally, the distance. This is used for localization to reduce sampling error.
      - Uses the default behavior for determining distance.
    * - **get_close_state**
      - Calculate *which state points are “close”* to a given location and, optionally, the distance. This is used for localization to reduce sampling error.
      - Uses the default behavior for determining distance.
    * - **nc_write_model_vars**
      - This method is not currently called, so just use the default routine.
      - Does nothing.
