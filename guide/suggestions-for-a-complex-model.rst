Suggestions for a "complex" model
=================================

A “complex” model is typically a large geophysical model where the model must be advanced outside of DART execution
control. Here we provide some advice on how to integrate this kind of model with DART.

First of all, the 4. **static_init_model**, 5. **get_state_meta_data()** and 6. **end_model()** suggestions will match
the multiple state variable in the previous section as complex models will typically have multiple fields.

An additional twist is that complex models may have different grid locations for different variables, (i.e. grid
staggering), but the above instructions still apply.

The 7. **adv_1step()** method for a complex model should fail, so the default behavior is sufficient.

The advice for the 8. **shortest_time_between_assimilations()** routine is similar to the advice for a simple model:
read the value from the namelist or return a fixed time as appropriate.

.. note::

   Since the model will not be advanced by DART, the value returned here is
   irrelevant except for user information purposes.

For the remaining routines, we give the following implementation suggestions:

+-----------------------------------+---------------------------------------------------------------------------------------+
| Routine # / name                  | Suggested implementation                                                              |
+===================================+=======================================================================================+
| 9. **model_interpolate()**        | Find the (*i*,\ *j*,\ *k*) indices which enclose that location, or search for the     |
|                                   | cell number. For some models you can compute (*i*,\ *j*) directly from a regular      |
|                                   | lat/lon grid, and in others you may have to search over a deformed grid. Any model    |
|                                   | code or utilities available for this purpose may prove very helpful as a starting     |
|                                   | point. In the end, you will use **get_state()** to retrieve an ensemble-sized array   |
|                                   | of values for each offset into the state vector, and then do interpolation to get an  |
|                                   | array of expected values.                                                             |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 10. **nc_write_model_atts()**     | It is very helpful (but optional) to add grid information to assist in plotting your  |
|                                   | results.                                                                              |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 11. **read_model_time()**         | (see **write_model_time()** below)                                                    |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 12. **write_model_time()**        | If the model time is stored in the netCDF files, supply routines that can read and    |
|                                   | write it in the correct format. The default routines will work if the model time      |
|                                   | matches what those routines expect: a time variable with an optional calendar         |
|                                   | variable. If no calendar is provided, the routine assumes fractional days. If the     |
|                                   | time variable is an array (i.e. more than one time step is stored in the file),       |
|                                   | read/write the last one.                                                              |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 13. **pert_model_copies()**       | The default of adding Gaussian noise to all state variables may be undesirable.       |
|                                   | Complex models often have a method to perturb a state according to a particular       |
|                                   | formula or method. Otherwise, it may be necessary to perturb each variable with       |
|                                   | separate noise levels, only perturb certain variables, etc.                           |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 14. **convert_vertical_obs()**    | (see **convert_vertical_state()** below)                                              |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 15. **convert_vertical_state()**  | Add code to convert between vertical coordinates (e.g. pressure, height, sigma        |
|                                   | levels, etc.) if appropriate. Code from the model or a model utility may be a very    |
|                                   | helpful starting point.                                                               |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 16. **get_close_obs()**           | (see **get_close_state()** below)                                                     |
+-----------------------------------+---------------------------------------------------------------------------------------+
| 17. **get_close_state()**         | If you want to change the localization impact based on something other than the type  |
|                                   | or kind, put code here. You should test for vertical type and do the conversion on    |
|                                   | demand if it hasn’t already been done.                                                |
+-----------------------------------+---------------------------------------------------------------------------------------+

As mentioned above, the most difficult routine to implement for a complex model is typically 9. **model_interpolate()**.
