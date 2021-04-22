DART Manhattan Differences from Lanai Release Notes
===================================================

Overview
--------

This document includes an overview of the changes in the DART system since the Lanai release. For further details on any
of these items look at the HTML documentation for that specific part of the system.

The two most significant changes in the Manhattan version of DART are it can support running models with a state vector
larger than the memory of a single task, removing a limit from the Lanai version of DART. It also reads and writes
NetCDF files directly instead of requiring a conversion from one file to another. There are many other smaller changes,
detailed below.

Manhattan supported models:

-  9var
-  bgrid_solo
-  cam-fv
-  cice
-  clm
-  cm1
-  forced_lorenz_96
-  ikeda
-  lorenz_63
-  lorenz_84
-  lorenz_96
-  lorenz_96_2scale
-  lorenz_04
-  mpas_atm (NetCDF overwrite not supported for update_u_from_reconstruct = .true. )
-  null_model
-  POP
-  ROMS
-  simple_advection
-  wrf

If your model of interest is not on the list consider checking out the 'Classic' release of DART, which is Lanai plus
bug fixes and minor enhancements. All models previously supported by Lanai are still in DART 'Classic'.

These are the major differences between the Lanai/Classic and Manhattan releases of DART:

-  Read and write NetCDF restarts
-  Calculation of forward operators
-  Vertical conversion of observation locations
-  Diagnostic file changes
-  :doc:`./state_structure`
-  model_mod interface changes
-  Observation Quantity replaces Kind
-  Perturbation of the state

NetCDF restart files
--------------------

The programs filter and perfect_model_obs now read/write directly from NetCDF files rather than having to run converters
(``model_to_dart`` and ``dart_to_model``). To facilitate this there is a new required call ``add_domain`` which must be
called during ``static_init_model``. It can be called multiple times in static_model_mod, e.g. once for each NetCDF file
that contains state variables. There are three ways to add a domain:

-  **From File** : This is for models which have NetCDF restart files

   -  ``dom_id = add_domain(template_file, num_vars, var_names, ... )``

-  **From Spec** : Creates a skeleton structure for a domain ( currently only used in bgrid_solo )

   -  ``dom_id = add_domain(num_vars, var_names, ... )``
   -  ``call add_dimension_to_variable(dom_id, var_id, dim_nam, dim_size)``
   -  ``call finished_adding_domain``

-  **From Blank** : This is for small models such as lorenz_96 and no NetCDF restarts

   -  ``dom_id = add_domain(model_size)``

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
detail about distributed vs. non-distributed forward operators. In distributed mode the forward operators for all
ensemble members are calculated in the same ``model_interpolate`` call. In non-distributed mode, the forward operators
for all ensemble members a task owns (1-ens_size) are calculated at once.

Vertical conversion of observation and state locations
------------------------------------------------------

The vertical conversion of observation locations is done before the assimilation by default. This can be changed by
namelist options.

In Lanai this calculation is done in the assimilation as part of ``get_close_obs`` if a model_mod does vertical
conversion. Note that not all models do vertical
conversion or even have a concept of vertical location, but every model_mod must have the following routines:

::

   call set_vertical_localization_coord(vert_localization_coord)

   call convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                             which_vert, status)

   call convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                               which_vert, istatus)

If there are NOT multiple choices for a vertical coordinate (e.g. cartesian, one dimensional), all these routines can be
no-ops.

If there are multiple types of vertical coordinates, the convert routines must be able to convert between them. The
'set_vertical_localization_coord()' routine should be called from 'static_init_model()' to set what localization
coordinate type is being requested.

The three routines related to vertical coordinates/localization choices are:

-  ``set_vert_localization_coord`` - sets the vertical localization coordiate (not required if there is no vertical
   conversion)
-  ``convert_vertical_obs`` - converts observation location to required vertical type (does nothing if there is no
   vertical conversion)
-  ``convert_vertical_state`` - converts state vector location to required vertical type (does nothing if there is no
   vertical conversion)

DART diagnostic file changes
----------------------------

For large models DART format diagnostic files (``Prior_Diag.nc`` and ``Posterior_Diag.nc``) have been replaced with
separate files for each copy that would have gone into Prior_Diag.nc and Posterior_Diag.nc.

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

The ``num_output_state_members`` are not written separately from the restarts. Note that restarts will have been clamped
if any clamping is applied (given as an arguement to add_domain). This is *different* to Posterior_Diag.nc which
contains unclamped values. Note also that there are 2 more "stages" which might be output, in addition to the preassim
and postassim discussed here.

For models with multiple domains the filenames above are appended with the domain number, e.g. preassim_mean.nc becomes
preassim_mean_d01.nc, preassim_mean_d02.nc, etc.

Changes to nc_write_model_atts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nc_write_model_atts`` now has 2 arguments:

-  ncid - open netcdf file identifier
-  domain_id - domain number being written

The calling code will write the model state, so this routine should only add attributes and optionally, non-state
information like grid arrays.

This routine will only be called if DART is creating an output NetCDF file from scratch. This may include any of the
preassim, postassim, or output files.

Changes to nc_write_model_vars
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nc_write_model_vars`` is currently unused (and in fact uncalled). It remains for possible future expansion.

Model_mod.f90 interface changes
-------------------------------

The model_mod.f90 file contains all code that is specific to any particular model. The code in this file is highly
constrained since these routines are \*called by\* other code in the DART system. All routine interfaces -- the names,
number of arguments, and the names of those arguments -- must match the prescribed interfaces exactly. Since not all
required interfaces are needed for every model there are default routines provided that can be referenced from a 'use'
statement and then the routine name can be put in the module 'public' list without any code for that routine having to
be written in the model_mod.f90 file.

The following 18 routines are required:

-  static_init_model
-  get_model_size
-  get_state_meta_data
-  shortest_time_between_assimilations
-  model_interpolate
-  end_model
-  nc_write_model_atts
-  nc_write_model_vars
-  init_time
-  init_conditions
-  adv_1step
-  pert_model_copies
-  get_close_obs
-  get_close_state
-  convert_vertical_obs
-  convert_vertical_state
-  read_model_time
-  write_model_time

Here is an example of code from the top of a model_mod file, including the modules where the default routines live and
the required public list.

::


   use     location_mod, only : location_type, get_close_type, &
                                get_close_obs, get_close_state, &
                                convert_vertical_obs, convert_vertical_state, &
                                set_location, set_location_missing, &
                                set_vertical_localization_coord
   use    utilities_mod, only : register_module, error_handler, &
                                E_ERR, E_MSG
                                ! nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                                ! find_namelist_in_file, check_namelist_read
   use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                    nc_add_global_creation_time, &
                                    nc_begin_define_mode, nc_end_define_mode
   use state_structure_mod, only : add_domain
   use ensemble_manager_mod, only : ensemble_type
   use dart_time_io_mod, only  : read_model_time, write_model_time
   use default_model_mod, only : pert_model_copies, nc_write_model_vars

   implicit none
   private

   ! required by DART code - will be called from filter and other
   ! DART executables.  interfaces to these routines are fixed and
   ! cannot be changed in any way.
   public :: static_init_model,      &
             get_model_size,         &
             get_state_meta_data,    &
             shortest_time_between_assimilations, &
             model_interpolate,      &
             end_model,              &
             nc_write_model_atts,    &
             adv_1step,              &
             init_time,              &
             init_conditions

   ! public but in another module
   public :: nc_write_model_vars,    &
             pert_model_copies,      &
             get_close_obs,          &
             get_close_state,        &
             convert_vertical_obs,   &
             convert_vertical_state, &
             read_model_time,        &
             write_model_time

Observation quantity replaces kinds
-----------------------------------

Historically there has been confusion about the terms for specific observation types (which often include the name of
the instrument collecting the data) and the generic quantity that is being measured (e.g. temperature). The previous
terms for these were 'types' and 'kinds', respectively.

Starting with the Manhattan release we have tried to clarify the terminology and make the interfaces consistent. The
following table lists the original names from the Lanai/Classic release and the replacement routines in Manhattan.

All code that is part of the DART code repository has been updated to use the replacment routines, but if you have your
own utilities written using this code, you will need to update your code. Contact us ( dart@ucar.edu ) for help if you
have any questions.

::


   public subroutines, existing name on left, replacement on right:

   assimilate_this_obs_kind()     =>     assimilate_this_type_of_obs(type_index)
   evaluate_this_obs_kind()       =>       evaluate_this_type_of_obs(type_index)
   use_ext_prior_this_obs_kind()  =>  use_ext_prior_this_type_of_obs(type_index)

   get_num_obs_kinds()      =>  get_num_types_of_obs()
   get_num_raw_obs_kinds()  =>  get_num_quantities()

   get_obs_kind_index()     => get_index_for_type_of_obs(type_name)
   get_obs_kind_name()      => get_name_for_type_of_obs(type_index)

   get_raw_obs_kind_index()  =>  get_index_for_quantity(quant_name)
   get_raw_obs_kind_name()   =>  get_name_for_quantity(quant_index)

   get_obs_kind_var_type()  =>  get_quantity_for_type_of_obs(type_index)

   get_obs_kind()      =>  get_obs_def_type_of_obs(obs_def)
   set_obs_def_kind()  =>  set_obs_def_type_of_obs(obs_def)

   get_kind_from_menu()      =>  get_type_of_obs_from_menu()

   read_obs_kind()     =>   read_type_of_obs_table(file_unit, file_format)
   write_obs_kind()    =>  write_type_of_obs_table(file_unit, file_format)

   maps obs_seq nums to specific type nums, only used in read_obs_seq:
   map_def_index()  => map_type_of_obs_table()  

   removed.  apparently unused, and simply calls get_obs_kind_name():
   get_obs_name()

   apparently unused anywhere, removed:
   add_wind_names()
   do_obs_form_pair()

   public integer parameter constants and subroutine formal argument names,
   old on left, new on right:

   KIND_ => QTY_
   kind => quantity

   TYPE_ => TYPE_
   type => type_of_obs

   integer parameters:
   max_obs_generic  =>  max_defined_quantities  (not currently public, leave private)
   max_obs_kinds    =>  max_defined_types_of_obs 

Additions/changes to existing namelists
---------------------------------------

Quality_control_nml
~~~~~~~~~~~~~~~~~~~

These namelist options used to be in filter_nml, now they are in quality_control_nml.

::

   &quality_control_nml
      input_qc_threshold          = 3,
      outlier_threshold           = 4,
      enable_special_outlier_code = .false.
   /

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

   +----------------------------------+--------------------------+----------------------------------------+
   | Item                             | Type                     | Description                            |
   +==================================+==========================+========================================+
   | single_file_in                   | logical                  | True means that all of the restart     |
   |                                  |                          | and inflation information is read      |
   |                                  |                          | from a single NetCDF file. False       |
   |                                  |                          | means that you must specify an         |
   |                                  |                          | input_state_file_list and DART will    |
   |                                  |                          | be expecting                           |
   |                                  |                          | input_{priorinf,postinf}_{mean,sd}.nc  |
   |                                  |                          | files for inflation.                   |
   +----------------------------------+--------------------------+----------------------------------------+
   | single_file_out                  | logical                  | True means that all of the restart     |
   |                                  |                          | and inflation information is written   |
   |                                  |                          | to a single NetCDF file. False means   |
   |                                  |                          | that you must specify a                |
   |                                  |                          | output_state_files and DART will be    |
   |                                  |                          | output files specified in the list.    |
   |                                  |                          | Inflation files will be written in     |
   |                                  |                          | the form                               |
   |                                  |                          | input_{priorinf,postinf}_{mean,sd}.nc. |
   +----------------------------------+--------------------------+----------------------------------------+
   | input_state_files                | character array          | This is used for single file input     |
   |                                  |                          | for low order models. For multiple     |
   |                                  |                          | domains you can specify a file for     |
   |                                  |                          | each domain. When specifying a list    |
   |                                  |                          | single_file_in, single_file_out must   |
   |                                  |                          | be set to .true.                       |
   +----------------------------------+--------------------------+----------------------------------------+
   | output_state_files               | character array          | This is used for single file input     |
   |                                  |                          | for low order models. For multiple     |
   |                                  |                          | domains you can specify a file for     |
   |                                  |                          | each domain. When specifying a list    |
   |                                  |                          | single_file_in, single_file_out must   |
   |                                  |                          | be set to .true.                       |
   +----------------------------------+--------------------------+----------------------------------------+
   | input_state_file_list            | character array          | A list of files containing input       |
   |                                  |                          | model restarts. For multiple domains   |
   |                                  |                          | you can specify a file for each        |
   |                                  |                          | domain. When specifying a list         |
   |                                  |                          | single_file_in, single_file_out must   |
   |                                  |                          | be set to .false.                      |
   +----------------------------------+--------------------------+----------------------------------------+
   | output_state_file_list           | character array          | A list of files containing output      |
   |                                  |                          | model restarts. For multiple domains   |
   |                                  |                          | you can specify a file for each        |
   |                                  |                          | domain. When specifying a list         |
   |                                  |                          | single_file_in, single_file_out must   |
   |                                  |                          | be set to .false.                      |
   +----------------------------------+--------------------------+----------------------------------------+
   | stages_to_write                  | character array          | Controls which stages to write.        |
   |                                  |                          | Case-insensitive input.                |
   |                                  |                          | Currently there are six options:       |
   |                                  |                          |                                        |
   |                                  |                          | -  ``input`` -- writes input mean and  |
   |                                  |                          |    sd only                             |
   |                                  |                          | -  ``forecast`` -- before              |
   |                                  |                          |    assimilation, before prior          |
   |                                  |                          |    inflation is applied                |
   |                                  |                          | -  ``preassim`` -- before              |
   |                                  |                          |    assimilation, before prior          |
   |                                  |                          |    inflation is applied                |
   |                                  |                          | -  ``postassim`` -- after              |
   |                                  |                          |    assimilation, before posterior      |
   |                                  |                          |    inflation is applied                |
   |                                  |                          | -  ``analysis`` -- after               |
   |                                  |                          |    assimilation, after posterior       |
   |                                  |                          |    inflation is applied                |
   |                                  |                          | -  ``output`` -- final output from     |
   |                                  |                          |    filter which includes clamping and  |
   |                                  |                          |    inflation                           |
   +----------------------------------+--------------------------+----------------------------------------+
   | write_all_stages_at_end          | logical                  | True means output all stages at the    |
   |                                  |                          | end of filter. This is more memory     |
   |                                  |                          | intensive but requires less time. For  |
   |                                  |                          | larger models IO begins to dominate    |
   |                                  |                          | the overall cost of the assimilation,  |
   |                                  |                          | so writting all stages at the end      |
   |                                  |                          | writes more files in parallel,         |
   |                                  |                          | reducing the IO time. Filenames are    |
   |                                  |                          | defined in ``output_state_files``.     |
   +----------------------------------+--------------------------+----------------------------------------+
   | output_restarts                  | logical                  | True means output a restart file(s).   |
   |                                  |                          | Filenames are defined in               |
   |                                  |                          | ``output_state_files``.                |
   +----------------------------------+--------------------------+----------------------------------------+
   | output_mean                      | logical                  | True means output a restart file       |
   |                                  |                          | which contains the ensemble mean for   |
   |                                  |                          | the stages that have been turned on    |
   |                                  |                          | in ``stages_to_write``. The file name  |
   |                                  |                          | will have the stage with ``_mean``     |
   |                                  |                          | appended.                              |
   +----------------------------------+--------------------------+----------------------------------------+
   | output_sd                        | logical                  | True means output a restart file       |
   |                                  |                          | which contains the ensemble standard   |
   |                                  |                          | deviation for the stages that have     |
   |                                  |                          | been turned on in                      |
   |                                  |                          | ``stages_to_write``. The file name     |
   |                                  |                          | will have the stage with ``_sd``       |
   |                                  |                          | appended.                              |
   +----------------------------------+--------------------------+----------------------------------------+
   | perturb_from_single_instance     | logical                  | Read a single file and perturb this    |
   |                                  |                          | to create an ensemble                  |
   +----------------------------------+--------------------------+----------------------------------------+
   | perturbation_amplitude           | float                    | Perturbation amplitude                 |
   +----------------------------------+--------------------------+----------------------------------------+
   | distribute_state                 | logical                  | True keeps the state distributed       |
   |                                  |                          | across all tasks throughout the        |
   |                                  |                          | entire execution of filter.            |
   +----------------------------------+--------------------------+----------------------------------------+

**NetCDF reads and writes:**

For **input** file names:

-  | give ``input_state_file_list`` a file for each domain, each of which contains a list of restart files. An example
     of an 'input_list.txt' might look something like :

   ::

      advance_temp1/wrfinput_d01
      advance_temp2/wrfinput_d01
      advance_temp3/wrfinput_d01
      advance_temp4/wrfinput_d01
      advance_temp5/wrfinput_d01
      ....

   | 

-  if no ``input_state_file_list`` is provided then default filenames will be used e.g. input_member_####.nc,
   input_priorinf_mean.nc, input_priorinf_sd.nc

For **output** file names:

-  | give ``output_state_file_list`` a file for each domain, each of which contains a list of restart files. An example
     of an 'input_list.txt' might for WRF might look something like :

   ::

      wrf_out_d01.0001.nc
      wrf_out_d01.0002.nc
      wrf_out_d01.0003.nc
      wrf_out_d01.0004.nc
      wrf_out_d01.0005.nc
      ....

   | 
   | if you would like to simply like to overwrite your previous data input_list.txt = output_list.txt

-  if no ``output_state_files`` is provided then default filenames will be used e.g. output_member_####.nc,
   output_priorinf_mean.nc, output_priorinf_sd.nc

For small models you may want to use ``single_file_in``, ``single_file_out`` which contains all copies needed to run
filter.

State_vector_io_nml
~~~~~~~~~~~~~~~~~~~

::

   &state_vector_io_nml
      buffer_state_io          = .false.,
      single_precision_output  = .false.,
   /

When ``buffer_state_io`` is ``.false.`` the entire state is read into memory at once if .true. variables are read one at
a time. If your model can not fit into memory at once this must be set to ``.true.`` .

``single_precision_output`` allows you to run filter in double precision but write NetCDF files in single presision

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

Perturbations
-------------

The option to perturb one ensemble member to produce an ensemble is in filter_nml:``perturb_from_single_instance``. The
model_mod interface is now ``pert_model_copies`` not ``pert_model_state``. Each task perturbs every ensemble member for
its own subsection of state. This is more complicated than the Lanai routine ``pert_model_state``, where a whole state
vector is available. If a model_mod does not provide a perturb interface, filter will do the perturbing with an
amplitude set in filter_nml:perturbation_amplitude. Note the perturb namelist options have been removed from
ensemble_manager_nml
