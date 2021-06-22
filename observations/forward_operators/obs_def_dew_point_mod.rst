MODULE obs_def_dew_point_mod
============================

Overview
--------

Provides a subroutine to calculate the dew point temperature from model temperature, 
specific humidity, and pressure.

Revision 2801 (April 2007) implements a more robust method (based on Bolton's 
Approximation) for calculating dew point. This has been further revised to 
avoid a numerical instability that could lead to failed forward operators for 
dewpoints almost exactly 0 C.

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod (most likely threed_sphere)
   assim_model_mod
   obs_kind_mod

Public interfaces
-----------------

=================================== ======================
*use obs_def_dew_point_mod, only :* get_expected_dew_point
=================================== ======================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *call get_expected_dew_point(state_vector, location, key, td, istatus)*
   ::

      real(r8),            intent(in)  :: state_vector
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: key
      real(r8),            intent(out) :: td
      integer,             intent(out) :: istatus

.. container:: indent1

   Calculates the dew point temperature (Kelvin).

   ================ =========================================================================
   ``state_vector`` A one dimensional representation of the model state vector
   ``location``     Location for this obs
   ``key``          Controls whether upper levels (key = 1) or 2-meter (key = 2) is required.
   ``td``           The returned dew point temperature value
   ``istatus``      Returned integer describing problems with applying forward operator
   ================ =========================================================================

| 

Files
-----

-  NONE

References
----------

#. Bolton, David, 1980: The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108, 1046-1053.

Error codes and conditions
--------------------------

+------------------------+-----------------------------------------------------------+----------------------------------------+
|         Routine        |                          Message                          |                 Comment                |
+========================+===========================================================+========================================+
| get_expected_dew_point | 'key has to be 1 (upper levels) or 2 (2-meter), got ',key | The input value of key is not allowed. |
+------------------------+-----------------------------------------------------------+----------------------------------------+
