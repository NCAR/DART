RMA notes
=========

In the RMA version of DART, the state vector is not required to be stored completely on any process. This is achieved
using Remote Memory Access (RMA). The RMA programing model allows processes to read (and write) memory on other
processors asynchronously. RMA DART supported models:

Before bitwise testing with Lanai please read :doc:`./bitwise_considerations`

NetCDF restarts
---------------

The programs filter and perfect_model_obs now read/write directly from NetCDF files, rather than having to run
converters (``model_to_dart`` and ``dart_to_model``). To facilitate this, there is a new required call ``add_domain``
which must be called during ``static_init_model``. It can be called multiple times in static_model_mod, e.g. once for
each NetCDF file that contains state variables. There are three ways to add a domain:

-  **From Blank** : This is for small models such as lorenz_96 and no NetCDF restarts

   -  ``dom_id = add_domain(model_size)``

-  **From File** : This is for models which have NetCDF restart files

   -  ``dom_id = add_domain(template_file, num_vars, var_names, ... )``

-  **From Spec** : Creates a skeleton structure for a domain ( currently only used in bgrid_solo )

   -  ``dom_id = add_domain(num_vars, var_names, ... )``
      ``call add_dimension_to_variable(dom_id, var_id, dim_nam, dim_size)``
      ``call finished_adding_domain``

For models without NetCDF restarts, use ``add_domain(model_size)``. This is the minimum amount of information needed by
DART to create a netdcf file. For models with NetCDF restarts use ``add_domain(info_file, num_vars, var_names)`` which
lets DART read the NetCDF dimensions for a list of variables from a file (``info_file``). There are several routines
that can be used together to create a domain from a description:
``add_domain, add_dimension_to_variable, finished_adding_domain``. This can be used in models such as bgrid_solo where
the model is spun up in perfect_model_obs, but the model itself has variable structure (3D variables with names). See
Additions/Changes to existing namelists for how to use NetCDF IO.

**Note** when using NetCDF restarts, inflation files are NetCDF also. The inflation mean and inflation standard
deviation are in separate files when you use NetCDF restarts. See :doc:`./netcdf_inflation_files` for details.

Calculation of forward operators
--------------------------------

The forward operator code in model_mod now operates on an array of state values. See :doc:`./forward_operator` for more
detail about distributed vs. non-distributed forward operators.

In distributed mode the forward operators for all ensemble members are calculated in the same ``model_interpolate``
call. In non-distributed mode, the forward oeprators for all ensemble members a task owns (1-ens_size) are calculated at
once.

Vertical conversion of observation locations
--------------------------------------------

The vertical conversion of observation locations is done before the assimilation. In Lanai this calculation is done in
the assimilation as part of ``get_close_obs`` if a model_mod does vertical conversion. See :doc:`./vertical_conversion`
for details about this change. Note that not all models do vertical conversion or even have a concept of vertical
location, but every model_mod must have the following routines:

-  ``query_vert_localization_coord`` - returns the vertical localization coordiate (or does nothing if there is no
   vertical conversion)
-  ``vert_convert`` - converts location to required vertical (or does nothing if there is no vertical conversion)

Diagnostic file changes
-----------------------

For large models DART format diagnostic files (Prior_Diag.nc and Posterior_Diag.nc) have been replaced with separate
files for each copy that would have gone into Prior_Diag.nc and Posterior_Diag.nc.

For Prior_Diag.nc:

-  **Mean and standard deviation**:
   preassim_mean.nc
   preassim_sd.nc
-  **Inflation mean and standard deviation** (if state space inflation is used):
   preassim_priorinf_mean.nc
   preassim_priorinf_sd.nc
-  **The number of ensemble members specifed** in filter_nml (num_output_state_members):
   preassim_member_####.nc

For Posterior_Diag.nc:

-  **Mean and standard deviation**:
   postassim_mean.nc
   postassim_sd.nc
-  **Inflation mean and standard deviation** (if state space inflation is used):
   postassim_priorinf_mean.nc
   postassim_priorinf_sd.nc
-  **The number of ensemble members specifed** in filter_nml (num_output_state_members):
   postassim_member_####.nc

The ``num_output_state_members`` are now written separately from the restarts. Note that restarts will have been clamped
if any clamping is applied (given as an arguement to add_domain). This is *different* to Posterior_Diag.nc which
contains unclamped values. Note also that there are 2 more "stages" which might be output, in addition to the preassim
and postassim discussed here.

For models with multiple domains the filenames above are appended with the domain number, e.g. preassim_mean.nc becomes
preassim_mean_d01.nc, preassim_mean_d02.nc, etc.

Changes to nc_write_model_atts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nc_write_model_atts`` has an new argument '``model_mod_writes_state_variables``'. This is used to communicate to DART
whether the model will create and write state variables in Prior_Diag.nc and Posterior_Diag.nc. If
``model_model_writes_state_variables = .false.`` DART will define and write state variables to the new diagnostic files.
If ``model_model_writes_state_variables = .true., nc_write_model_vars`` is called as normal.

Perturbations
-------------

The option to perturb one ensemble member to produce an ensemble is in filter_nml:``perturb_from_single_instance``. The
model_mod interface is now ``pert_model_copies`` not ``pert_model_state``. Each task perturbs every ensemble member for
its own subsection of state. This is more complicated than the Lanai routine ``pert_model_state``, where a whole state
vector is available. If a model_mod does not provide a perturb interface, filter will do the perturbing with an
amplitude set in filter_nml:perturbation_amplitude. Note the perturb namelist options have been removed from
ensemble_manager_nml

State_vector_io_nml
-------------------

::

   &state_vector_io_nml
      buffer_state_io         = .false.,
      single_precision_output = .false.,
   /

When ``buffer_state_io`` is ``.false.`` the entire state is read into memory at once if .true. variables are read one at
a time. If your model can not fit into memory at once this must be set to ``.true.`` .

``single_precision_output`` allows you to run filter in double precision but write NetCDF files in single precision.

Quality_control_nml
-------------------

These namelist options used to be in filter_nml, now they are in quality_control_nml.

::

   &quality_control_nml
      input_qc_threshold          = 3,
      outlier_threshold           = 4,
      enable_special_outlier_code = .false.
   /

Additions/changes to existing namelists
---------------------------------------

New namelist variables

filter_nml
~~~~~~~~~~

::

   &filter_nml
      single_file_in               = .false.,
      single_file_out              = .false.,

      input_state_file_list        = 'null',
      output_state_file_list       = 'null',
      input_state_files            = 'null',
      output_state_files           = 'null',

      stages_to_write              = 'output'
      write_all_stages_at_end      = .false.
      output_restarts              = .true.
      output_mean                  = .true.
      output_sd                    = .true.

      perturb_from_single_instance = .false.,
      perturbation_amplitude       = 0.2_r8,

      distributed_state            = .true.
   /

| 

.. container::

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Item                                  | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | single_file_in                        | logical                               | True means that all of the restart    |
   |                                       |                                       | and inflation information is read     |
   |                                       |                                       | from a single NetCDF file. False      |
   |                                       |                                       | means that you must specify an        |
   |                                       |                                       | input_state_file_list and DART will   |
   |                                       |                                       | be expecting                          |
   |                                       |                                       | input_{priorinf,postinf}_{mean,sd}.nc |
   |                                       |                                       | files for inflation.                  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | single_file_out                       | logical                               | True means that all of the restart    |
   |                                       |                                       | and inflation information is written  |
   |                                       |                                       | to a single NetCDF file. False means  |
   |                                       |                                       | that you must specify a               |
   |                                       |                                       | output_state_file_list and DART will  |
   |                                       |                                       | be output files specified in the      |
   |                                       |                                       | list. Inflation files will be written |
   |                                       |                                       | in the form                           |
   |                                       |                                       | i                                     |
   |                                       |                                       | nput_{priorinf,postinf}_{mean,sd}.nc. |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | input_restart_files                   | character array                       | This is used for single file input    |
   |                                       |                                       | for low order models. For multiple    |
   |                                       |                                       | domains you can specify a file for    |
   |                                       |                                       | each domain. When specifying a list   |
   |                                       |                                       | single_file_in, single_file_out must  |
   |                                       |                                       | be set to .true.                      |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_restart_files                  | character array                       | This is used for single file input    |
   |                                       |                                       | for low order models. For multiple    |
   |                                       |                                       | domains you can specify a file for    |
   |                                       |                                       | each domain. When specifying a list   |
   |                                       |                                       | single_file_in, single_file_out must  |
   |                                       |                                       | be set to .true.                      |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | input_state_file_list                 | character array                       | A list of files containing input      |
   |                                       |                                       | model restarts. For multiple domains  |
   |                                       |                                       | you can specify a file for each       |
   |                                       |                                       | domain. When specifying a list        |
   |                                       |                                       | single_file_in, single_file_out must  |
   |                                       |                                       | be set to .false.                     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_state_file_list                | character array                       | A list of files containing output     |
   |                                       |                                       | model restarts. For multiple domains  |
   |                                       |                                       | you can specify a file for each       |
   |                                       |                                       | domain. When specifying a list        |
   |                                       |                                       | single_file_in, single_file_out must  |
   |                                       |                                       | be set to .false.                     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | stages_to_write                       | character array                       | Controls which stages to write.       |
   |                                       |                                       | Case-insensitive input.               |
   |                                       |                                       | Currently there are six options:      |
   |                                       |                                       |                                       |
   |                                       |                                       | -  ``input`` -- writes input mean and |
   |                                       |                                       |    sd only                            |
   |                                       |                                       | -  ``forecast`` -- before             |
   |                                       |                                       |    assimilation, before prior         |
   |                                       |                                       |    inflation is applied               |
   |                                       |                                       | -  ``preassim`` -- before             |
   |                                       |                                       |    assimilation, before prior         |
   |                                       |                                       |    inflation is applied               |
   |                                       |                                       | -  ``postassim`` -- after             |
   |                                       |                                       |    assimilation, before posterior     |
   |                                       |                                       |    inflation is applied               |
   |                                       |                                       | -  ``analysis`` -- after              |
   |                                       |                                       |    assimilation, after posterior      |
   |                                       |                                       |    inflation is applied               |
   |                                       |                                       | -  ``output`` -- final output from    |
   |                                       |                                       |    filter which includes clamping and |
   |                                       |                                       |    inflation                          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | write_all_stages_at_end               | logical                               | True means output all stages at the   |
   |                                       |                                       | end of filter. This is more memory    |
   |                                       |                                       | intensive but requires less time. For |
   |                                       |                                       | larger models IO begins to dominate   |
   |                                       |                                       | the overall cost of the assimilation, |
   |                                       |                                       | so writting all stages at the end     |
   |                                       |                                       | writes more files in parallel,        |
   |                                       |                                       | reducing the IO time.                 |
   |                                       |                                       | ``output_state_file_list``.           |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_restarts                       | logical                               | True means output a restart file(s).  |
   |                                       |                                       | Filenames are defined in              |
   |                                       |                                       | ``output_state_file_list``.           |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_mean                           | logical                               | True means output a restart file      |
   |                                       |                                       | which contains the ensemble mean for  |
   |                                       |                                       | the stages that have been turned on   |
   |                                       |                                       | in ``stages_to_write``. The file name |
   |                                       |                                       | will have the stage with ``_mean``    |
   |                                       |                                       | appended.                             |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_sd                             | logical                               | True means output a restart file      |
   |                                       |                                       | which contains the ensemble standard  |
   |                                       |                                       | deviation for the stages that have    |
   |                                       |                                       | been turned on in                     |
   |                                       |                                       | ``stages_to_write``. The file name    |
   |                                       |                                       | will have the stage with ``_sd``      |
   |                                       |                                       | appended.                             |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | perturb_from_single_instance          | logical                               | Read a single file and perturb this   |
   |                                       |                                       | to create an ensemble                 |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | perturbation_amplitude                | float                                 | Perturbation amplitude                |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | distribute_state                      | logical                               | True keeps the state distributed      |
   |                                       |                                       | across all tasks throughout the       |
   |                                       |                                       | entire execution of filter.           |
   +---------------------------------------+---------------------------------------+---------------------------------------+

**For NetCDF reads and writes**

For **input** file names:

-  give ``input_state_file_list`` a file for each domain, each of which contains a list of restart files.
-  if no ``input_state_file_list`` is provided then default filenames will be used e.g. input_member_000*.nc,
   input_priorinf_mean.nc, input_priorinf_sd.nc

For **output** file names:

-  give ``output_state_file_list`` a file for each domain, each of which contains a list of restart files.
-  if no ``output_state_file_list`` is provided then default filenames will be used e.g. output_member_000*.nc,
   output_priorinf_mean.nc, output_priorinf_sd.nc

For small models you may want to use ``single_file_in``, ``single_file_out`` which contains all copies needed to run
filter.

Assim_tools_nml
~~~~~~~~~~~~~~~

::

   &assim_tools_nml
      distribute_mean  = .true.
   /

In previous DART releases, each processor gets a copy of the mean (in ens_mean_for_model). In RMA DART, the mean is
distributed across all processors. However, a user can choose to have a copy of the mean on each processor by setting
``distribute_mean = .false.`` . Note that the mean state is accessed through ``get_state`` whether distribute_mean is
``.true.`` or ``.false.``

Removed from existing namelists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   &filter_nml
      input_qc_threshold          = 3,
      outlier_threshold           = 4,
      enable_special_outlier_code = .false.
      start_from_restart          = .false.
      output_inflation            = .true.
      output_restart              = .true.
      /

NOTE : ``output_restart`` has been renamed to ``output_restarts``. **``output_inflation`` is no longer supported** and
only writes inflation files if ``inf_flavor > 1``

::

   &ensemble_manager_nml
      single_restart_file_out = .true.
      perturbation_amplitude  = 0.2,
      /

::

   &assim_manager_nml
      write_binary_restart_files = .true.,
      netCDF_large_file_support  = .false.
      /
