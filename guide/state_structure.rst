State Stucture
==============

state_structure_mod is a module that holds all the domain, variable, dimension info about the model_mods in the state.
Note it stores only **metadata** about the state, not the actual state variables themselves.

It is the foundation for two parts of the code:

-  Read/write state variables from/to netcdf files
-  Calculate DART index from x,y,z variable indices and the inverse: x,y,z, variable from DART index.

Inside ``static_init_model`` a call is made to ``add_domain``. This call is *required* as it communicates to the state
structure that a new domain has been added to the state. The state structure keeps track of the number of domains in the
state. These may be multiple domains in one model_mod, e.g. nested domains in WRF, or multiple model_mods, e.g. POP
coupled with CAM. The minimum amount of information ``add_domain`` needs is model size which means vector of length
model size has been added to the state. This equivalent to Lanai where the only information filter has is that the model
is a vector of length model_size. For models with netcdf restart files you supply ``add_domain`` with:

-  a netcdf file
-  the number of variables to read from the file
-  the name of the variables
-  Optionally:

   -  the DART KINDS of the variables
   -  clamping upper and lower bounds
   -  update/not update this variable

For models that are spun up in perfect_model_obs you can manually describe the variables so you can create netcdf files
containing the varibles in the model state, e.g. Temperature, Surface Pressure, etc. There are 3 steps to this process:

#. Supply ``add_domain`` with almost the same arguments as you would for a netcdf file, but skip the first arguement
   (netcdf filename).
#. For each variable, loop around the required number of dimensions and call ``add_dimension_to_variable``
#. Call ``finished_adding_domain`` to let the state structure know that you have finished adding dimensions to
   variables.

DART index
^^^^^^^^^^

| To get the dart index for an i,j,k,variable in a domain use:
| ``get_dart_vector_index(i, j, k, dom_id, var_id)``

| To get the i,j,k, variable, domain from the dart index use:
| ``get_model_variable_indices(dart_index, i, j, k, var_id, dom_id)``

**Note** That (i,j,k) needs to be converted to (lon, lat, lev) or to whatever grid the variable is on. ``get_dim_name``
can be used to get the dimension name from i,j,k if needed.

Unlimited dimensions: io vs model_mod routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some model restart files have an unlimited dimension. For IO purposes, e.g. creating netcdf files, the unlimited
dimension is used. For state structure accessor functions called be the model_mod the unlimited dimension is ignored. So
if you have a variable TEMPERATURE in your netcdf file, with dimensions (lon, lat, level, time) the IO routines will see
a 4D variable, but ``get_num_dims`` used in model_mod will return 3D.
