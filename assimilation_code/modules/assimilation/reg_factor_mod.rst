MODULE reg_factor
=================

Overview
--------

Computes a weighting factor to reduce the impact of observations on state variables using information from groups of
ensembles. Can be run using groups or using archived summary information available from previous group filter
experiments.

Other modules used
------------------

::

   types_mod
   utilities_mod
   time_manager_mod

Public interfaces
-----------------

============================ ===============
*use reg_factor_mod, only :* comp_reg_factor
============================ ===============

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *var = comp_reg_factor(num_groups, regress, obs_index, state_index [, obs_state_ind] [, obs_state_max])*
   ::

      real(r8)                                    ::  comp_reg_factor 
      integer,                         intent(in) ::  num_groups 
      real(r8), dimension(num_groups), intent(in) ::  regress 
      integer,                         intent(in) ::  obs_index 
      integer,                         intent(in) ::  state_index 
      integer, optional,               intent(in) ::  obs_state_ind 
      integer, optional,               intent(in) ::  obs_state_max 

.. container:: indent1

   Returns a weighting factor given regression factors from each group of a group filter or retrieves a factor generated
   by previous group filter runs from a file.

   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``num_groups``  | Number of groups. Set to 1 when using information from previously run group filter from file.     |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``regress``     | Regression factor from each group for a given state variable and observation variable pair.       |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``obs_index``   | Integer index of the observation being processed. Not used in current implementation .            |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``state_index`` | Integer index of state variable being processed. Not used in current implementation.              |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | *obs_state_ind* | Index into file generated for Bgrid model which could be duplicated in other large models.        |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | *obs_state_max* | Maximum number of observation state variable pairs with non-zero impacts for a given model and    |
   |                 | observation sequence. Used for generating Bgrid statistic files.                                  |
   +-----------------+---------------------------------------------------------------------------------------------------+

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &reg_factor_nml
      select_regression    = 1,
      input_reg_file       = "time_mean_reg",
      save_reg_diagnostics = .false.,
      reg_diagnostics_file = "reg_diagnostics"  
   /

| 

.. container::

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Item                                  | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | select_regression                     | integer                               | Selects the method for computing      |
   |                                       |                                       | regression factor.                    |
   |                                       |                                       |                                       |
   |                                       |                                       | -  1 = compute using sampling theory  |
   |                                       |                                       |    for any ensemble size.             |
   |                                       |                                       | -  2 = low order model format. Works  |
   |                                       |                                       |    from archived time mean or time    |
   |                                       |                                       |    median regression files generated  |
   |                                       |                                       |    by low-order models like           |
   |                                       |                                       |    Lorenz-96.                         |
   |                                       |                                       | -  3 = selects bgrid archived file.   |
   |                                       |                                       |    This is not currently supported in |
   |                                       |                                       |    released versions.                 |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | input_reg_file                        | character(len=129)                    | File name from which statistics are   |
   |                                       |                                       | to be read for select_regression = 3. |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | save_reg_diagnostics                  | logical                               | True if regression diagnostics should |
   |                                       |                                       | be computed.                          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | reg_diagnostics_file                  | character(len=129)                    | File name to which to write           |
   |                                       |                                       | diagnostics.                          |
   +---------------------------------------+---------------------------------------+---------------------------------------+

| 

Files
-----

-  (optional) input regression file from namelist variable input_reg_file.
-  reg_factor_mod.nml in input.nml

================================================== ===============================
filename                                           purpose
================================================== ===============================
from ``input.nml``\ &reg_factor_mod:input_reg_file file of regression coefficients
================================================== ===============================

References
----------

-  none

Private components
------------------

N/A
