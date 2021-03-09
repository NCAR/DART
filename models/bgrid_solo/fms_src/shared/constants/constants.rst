module constants_mod
====================

Overview
--------

Defines useful constants for Earth in mks units.

.. container::

   Constants are defined as real parameters.They are accessed through the "use" statement. While the local name of
   constant may be changed, their values can not be redefined.

| 

Other modules used
------------------

.. container::

   ::

      fms_mod

Public interface
----------------

.. container::

   ::

      use constants_mod [, only:  constants_init ]

   constants_init:
      A optional initialization routine. The only purpose of this routine is to write the version and tag name
      information to the log file.

| 

Public data
-----------

.. container::

   ======== ==== ====================== ========== ======================================================
   Name     Type Value                  Units      Description
   ======== ==== ====================== ========== ======================================================
   RADIUS   real 6376.e3                meters     radius of the earth
   OMEGA    real 7.292e-5               1/sec      rotation rate of the planet (earth)
   GRAV     real 9.80                   m/s2       acceleration due to gravity
   RDGAS    real 287.04                 J/kg/deg   gas constant for dry air
   KAPPA    real 2./7.                             RDGAS / CP
   CP       real RDGAS/KAPPA            J/kg/deg   specific heat capacity of dry air at constant pressure
   RVGAS    real 461.50                 J/Kg/deg   gas constant for water vapor
   DENS_H2O real 1000.                  Kg/m3      density of liquid water
   HLV      real 2.500e6                J/Kg       latent heat of evaporation
   HLF      real 3.34e5                 J/kg       latent heat of fusion
   HLS      real 2.834e6                J/Kg       latent heat of sublimation
   TFREEZE  real 273.16                 deg K      temp where fresh water freezes
   STEFAN   real 5.6734e-8              (W/m2/deg4 Stefan-Boltzmann constant
   VONKARM  real 0.40                   ---        Von Karman constant
   PI       real 3.14159265358979323846 ---        is it enough?
   ======== ==== ====================== ========== ======================================================

Public routines
---------------

a. .. rubric:: Constants_init
      :name: constants_init

   ::

      call constants_init 

   **DESCRIPTION**
      The only purpose of this routine is to write the version and tag name information to the log file. This routine
      does not have to be called. If it is called more than once or called from other than the root PE it will return
      silently. There are no arguments.

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   None.

References
----------

.. container::

   None.

| 

Compiler specifics
------------------

.. container::

   None.

| 

Precompiler options
-------------------

.. container::

   None.

| 

Loader options
--------------

.. container::

   None.

Test PROGRAM
------------

.. container::

   None.

| 

Notes
-----

.. container::

   <B>NOTES ON USAGE:</B>
   All constants have been declared as type REAL, PARAMETER.
   The value a constant can not be changed in a users program. New constants can be defined in terms of values from the
   constants module using a parameter statement.<br><br>
   The name given to a particular constant may be changed.<br><br>
   Constants can be used on the right side on an assignment statement (their value can not be reassigned).
   <B>EXAMPLES:</B>
   ::

           use constants_mod, only:  TFREEZE, grav_new => GRAV
           real, parameter :: grav_inv = 1.0 / grav_new
           tempc(:,:,:) = tempk(:,:,:) - TFREEZE
           geopotential(:,:) = height(:,:) * grav_new

| 
