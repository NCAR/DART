.. _Required model_mod routines:

Required model_mod routines
===========================

There are 18 Fortran subroutines necessary to implement in order to
successfully integrate a model in DART. You will place these routines in your
``model_mod.f90`` in a subdirectory with the name of your model in
``DART/models``. There is often a sensible default implementation that can be
used for each of these routines. For example, in the case of a model that
starts at a time of “0”, for the required routine **init_time()** the following
code will use this default implementation:

.. code-block:: fortran

   use default_model_mod,     only : init_time

As in all Fortran programs, a comma-separated list of routines can be listed
after the colon.

.. important

   Do not “use” the entire module without the keyword “only” in order to avoid
   including the default behavior for all subroutines contained in that module
   (in this example ``default_model_mod``).

The following table lists each of the 18 routines, their default modules
relative to ``DART``, and the default behavior. If the default behavior is not
desired, see the section :doc:`suggestions-for-a-simple-model` for a model that
DART can advance, or :doc:`suggestions-for-a-complex-model` for a model that is
advanced externally from DART.

+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| Routine # / name                                | Purpose                                                                        | Default module / directory                                                          | Default behavior                                  |
+=================================================+================================================================================+=====================================================================================+===================================================+
| 1. **init_time()**                              | Set the *initial time* if not read from the restart file.                      | ``default_model_mod`` / ``models/utilities``                                        | Sets the initial time to 0 days, 0 seconds        |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 2. **init_conditions()**                        | For a “cold start” fill in an empty state vector with *initial conditions*.    | ``default_model_mod`` / ``models/utilities``                                        | Sets the initial state to 0. To fail use          |
|                                                 | Many models cannot just make up values from thin air and thus choose to fail   |                                                                                     | ``init_conditions => fail_init_conditions``.      |
|                                                 | when this is requested.                                                        |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 3. **get_model_size()**                         | Return the *number of items in the state vector*.                              | ``default_model_mod`` / ``models/utilities``                                        | Returns 1; i.e. there is only one item in the     |
|                                                 |                                                                                |                                                                                     | state.                                            |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 4. **static_init_model()**                      | *Initialize* DART with information about the model that will be used by the    | ``default_model_mod`` / ``models/utilities``                                        | Does nothing.                                     |
|                                                 | remaining ``model_mod`` routines. The procedure for doing this will depend on  |                                                                                     |                                                   |
|                                                 | how complex the model is; see below for suggestions for implementation.        |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 5. **get_state_meta_data()**                    | Takes an index into the state vector and returns the *location* corresponding  | ``default_model_mod`` / ``models/utilities``                                        | Sets a missing location and the default variable  |
|                                                 | to that value and optionally the *variable type*. See below for suggestions on |                                                                                     | type.                                             |
|                                                 | implementation.                                                                |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 6. **end_model()**                              | *Deallocate* any arrays allocated in **static_init_model()**.                  | ``default_model_mod`` / ``models/utilities``                                        | Does nothing.                                     |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 7. **adv_1step()**                              | If possible, *advance the model* state from one time to another. Complex       | ``default_model_mod`` / ``models/utilities``                                        | Call the error handler with the message “unable   |
|                                                 | models will be unable to implement this method and should fail.                |                                                                                     | to advance model”.                                |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 8. **shortest_time_between_assimilations()**    | Return a namelist or a fixed value for the *minimum model advance time*        | ``default_model_mod`` / ``models/utilities``                                        | Returns a time period of 1 day.                   |
|                                                 | between assimilations. Note that complex models will handle advancing the time |                                                                                     |                                                   |
|                                                 | externally.                                                                    |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 9. **model_interpolate()**                      | *Interpolate* a requested quantity to the given location to get an array of    | ``default_model_mod`` / ``models/utilities``                                        | Fail and set the expected observation to          |
|                                                 | expected values for all ensemble members. \ *NOTE*: this is often the most     |                                                                                     | “missing.”                                        |
|                                                 | time consuming method to implement.                                            |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 10. **nc_write_model_atts()**                   | Add any *additional information* to the netCDF output diagnostic files.        | ``default_model_mod`` / ``models/utilities``                                        | Does nothing.                                     |
|                                                 | *NOTE*: the state will already be output by other routines, so this method     |                                                                                     |                                                   |
|                                                 | should not create or write the state variables.                                |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 11. **read_model_time()**                       | *Read* the model time from a state vector netCDF file.                         | ``dart_time_io`` / ``assimilation_code/io/utilities``                               | Attempt to read the “time” variable from a state  |
|                                                 |                                                                                |                                                                                     | file in an intelligent way.                       |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 12. **write_model_time()**                      | *Write* the model time to a state vector netCDF file.                          | ``dart_time_io`` / ``assimilation_code/io/utilities``                               | Write the “time” variable from the file according |
|                                                 |                                                                                |                                                                                     | to the DART calendar.                             |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 13. **pert_model_copies()**                     | *Perturb* a state vector in order to create an ensemble.                       | ``default_model_mod`` / ``models/utilities``                                        | Add Gaussian noise with a specified amplitude to  |
|                                                 |                                                                                |                                                                                     | all parts of the state vector.                    |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 14. **convert_vertical_obs()**                  | Some 3D models have multiple vertical coordinates (e.g. pressure, height, or   | ``location_mod/`` ``assimilation_code/`` ``location/XXX``                           | Do no conversion.                                 |
|                                                 | model level); this method *converts observations* between different vertical   |                                                                                     |                                                   |
|                                                 | coordinate systems.                                                            |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 15. **convert_vertical_state()**                | Some 3D models have multiple vertical coordinates (e.g. pressure, height, or   | ``location_mod/`` ``assimilation_code/`` ``location/XXX``                           | Do no conversion.                                 |
|                                                 | model level); this method *converts state* between different vertical          |                                                                                     |                                                   |
|                                                 | coordinate systems.                                                            |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 16. **get_close_obs()**                         | Calculate *which observations are “close”* to a given location and,            | ``location_mod/`` ``assimilation_code/`` ``location/XXX``                           | Uses the default behavior for determining         |
|                                                 | optionally, the distance. This is used for localization to reduce sampling     |                                                                                     | distance.                                         |
|                                                 | error.                                                                         |                                                                                     |                                                   |
|                                                 |                                                                                |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 17. **get_close_state()**                       | Calculate *which state points are “close”* to a given location and,            | ``location_mod/`` ``assimilation_code/`` ``location/XXX``                           | Uses the default behavior for determining         |
|                                                 | optionally, the distance. This is used for localization to reduce sampling     |                                                                                     | distance.                                         |
|                                                 | error.                                                                         |                                                                                     |                                                   |
|                                                 |                                                                                |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
| 18. **nc_write_model_vars()**                   | This method is not currently called, so just use the default routine for now.  | ``default_model_mod`` / ``models/utilities``                                        | Does nothing.                                     |
|                                                 | This method will be used in a future implementation.                           |                                                                                     |                                                   |
+-------------------------------------------------+--------------------------------------------------------------------------------+-------------------------------------------------------------------------------------+---------------------------------------------------+
