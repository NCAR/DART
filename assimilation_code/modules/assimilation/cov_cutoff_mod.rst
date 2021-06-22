MODULE cov_cutoff_mod
=====================

Overview
--------

Computes the weight with which an observation should impact a state variable that is separated by a given distance. The
distance is in units determined by the location module being used.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &cov_cutoff_nml
      select_localization = 1  
   /

| 

.. container::

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Item                                  | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | select_localization                   | integer                               | Selects the localization function.    |
   |                                       |                                       |                                       |
   |                                       |                                       | -  1 = Gaspari-Cohn 5th order         |
   |                                       |                                       |    polynomial with halfwidth c.       |
   |                                       |                                       | -  2 = Boxcar with halfwidth c (goes  |
   |                                       |                                       |    to 0 for z_in < 2c).               |
   |                                       |                                       | -  3 = Ramped Boxcar. Has value 1 for |
   |                                       |                                       |    z_in < c and then reduces linearly |
   |                                       |                                       |    to 0 at z_in = 2c.                 |
   +---------------------------------------+---------------------------------------+---------------------------------------+

| 

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod

Public interfaces
-----------------

============================ ===============
*use cov_factor_mod, only :* comp_cov_factor
============================ ===============

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *var = comp_cov_factor(z_in, c [, obs_loc] [, obs_type] [, target_loc] [, target_kind] [, localization_override])*
   ::

      real(r8)                                  :: comp_cov_factor
      real(r8), intent(in)                      :: z_in
      real(r8), intent(in)                      :: c
      type(location_type), optional, intent(in) :: obs_loc
      integer, optional, intent(in)             :: obs_type
      type(location_type), optional, intent(in) :: target_loc
      integer, optional, intent(in)             :: target_kind
      integer, optional, intent(in)             :: localization_override

.. container:: indent1

   Returns a weighting factor for observation and a target variable (state or observation) separated by distance z_in
   and with a half-width distance, c. Three options are provided and controlled by a namelist parameter. The optional
   argument localization_override controls the type of localization function if present. The optional arguments obs_loc,
   obs_type and target_loc, target_kind are not used in the default code. They are made available for users who may want
   to design more sophisticated localization functions.

   ======================= =========================================================================================
   ``var``                 Weighting factor.
   ``z_in``                The distance between observation and target.
   ``c``                   Factor that describes localization function. Describes half-width of functions used here.
   *obs_loc*               Location of the observation.
   *obs_type*              Observation specific type.
   *target_loc*            Location of target.
   *target_kind*           Generic kind of target.
   *localization_override* Controls localization type if present. Same values as for namelist control.
   ======================= =========================================================================================

| 

Files
-----

========= ==========================
filename  purpose
========= ==========================
input.nml to read ``cov_cutoff_nml``
========= ==========================

References
----------

#. Gaspari and Cohn, 1999, QJRMS, **125**, 723-757. (eqn. 4.10)


Error codes and conditions
--------------------------

+-----------------+-------------------------------------------------------------------+---------------------------------------------------------+
|     Routine     |                              Message                              |                         Comment                         |
+=================+===================================================================+=========================================================+
| comp_cov_factor | Illegal value of "select_localization" in cov_cutoff_mod namelist | Only values 1 through 3 select a localization function. |
+-----------------+-------------------------------------------------------------------+---------------------------------------------------------+


Private components
------------------

N/A
