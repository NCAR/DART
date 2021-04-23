Suggestions for a "simple" model
================================

A “simple” model is one where DART can advance the model through a function call. As we saw above, Lorenz 63 falls into
this category and can be used as a reference. Here we provide some further advice on how to add this kind of model to
DART.

The primary consideration with a simple model is how you will store the state. If you have only a *single type of
variable* in your state vector (for example, the Lorenz 63 model), here are some hints on how to implement your
initialization and meta data routines:

+------------------------------+---------------------------------------------------------------------------------------+
| Routine # / name             | Suggested implementation                                                              |
+==============================+=======================================================================================+
| 4. **static_init_model()**   | Your model_size will likely be set by namelist, so read it, allocate an array of that |
|                              | size, and precompute all the locations for each state vector item. Call               |
|                              | **add_domain()** with the model size so DART knows how long the state vector is.      |
+------------------------------+---------------------------------------------------------------------------------------+
| 5. **get_state_meta_data()** | Return ``QTY_STATE_VARIABLE`` as the quantity, and return the location for that index |
|                              | through a look-up into the location array created during **static_init\_ model()**.   |
+------------------------------+---------------------------------------------------------------------------------------+

If you have *more than a single type of variable* in the state vector (for example, “concentration”, “wind”, etc. as in
the ``DART/models/simple_advection`` model):

+------------------------------+---------------------------------------------------------------------------------------+
| Routine # / name             | Suggested implementation                                                              |
+==============================+=======================================================================================+
| 4. **static_init_model()**   | Read from the namelist the number of fields to be used in the state vector. Use       |
|                              | **add_domain()** to indicate which netCDF vars should be read. Read in any auxiliary  |
|                              | data needed by interpolation code (for example, the grid topology). Cache the grid    |
|                              | locations of the state variables as appropriate, and use **get_domain_size()** to     |
|                              | compute the model_size.                                                               |
+------------------------------+---------------------------------------------------------------------------------------+
| 5. **get_state_meta_data()** | Call **get_model_variable_indices()** and **get_state_kind()** to figure out the      |
|                              | (*i*,\ *j*,\ *k*) indices and which variable this offset is. Use the                  |
|                              | (*i*,\ *j*,\ *k*) index to compute the grid location and return it along with the     |
|                              | quantity.                                                                             |
+------------------------------+---------------------------------------------------------------------------------------+

Now, for either type of simple model, the following applies:

+-----------------------------------------------+---------------------------------------------------------------------------------------+
| Routine # / name                              | Suggested implementation                                                              |
+===============================================+=======================================================================================+
| 6. **end_model()**                            | Deallocate any arrays allocated in **static_init_model()**                            |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| 7. **adv_1step()**                            | If possible, embed the code that computes **x**\ (*t*\ +1) = **F**\ (**x**\ (*t*)) or |
|                                               | call a separate subroutine to advance the model state from one time to another.       |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| 8. **shortest_time_between_assimilations()**  | Return a namelist or a fixed value for the minimum model advance time.                |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| 9. **model_interpolate()**                    | Find the (*i*,\ *j*,\ *k*) indices which enclose that location, or search for the     |
|                                               | cell number. For some models you can compute (*i*,\ *j*) directly from a regular      |
|                                               | lat/lon grid, and in others you may have to search over a deformed grid. Any model    |
|                                               | code or utilities available for this purpose may prove very helpful as a starting     |
|                                               | point. In the end, you will use **get_state()** to retrieve an ensemble-sized array   |
|                                               | of values for each offset into the state vector, and then do interpolation to get an  |
|                                               | array of expected values.                                                             |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| 10. **nc_write_model_atts()**                 | Optionally add any desired attributes to the output diagnostic files.                 |
+-----------------------------------------------+---------------------------------------------------------------------------------------+

The remaining routines can mostly use the defaults, except possibly for 11. **read_model_time()** and 12.
**write_model_time()**, which may need to be customized if using a model restart file that already stores time in a
particular format.

Note that there is often no need to convert vertical obs or states in a simple model without vertical coordinate
choices.
