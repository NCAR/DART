module platform_mod
===================

Overview
--------

.. container::

   ``platform_mod`` is a module that provides public entities whose value may depend on the operating system and
   compiler.

.. container::

   | The combination of operating system and compiler is referred to here as the *platform*. ``platform_mod`` provides:

   ::

        integer, parameter :: r8_kind, r4_kind, &
                              c8_kind, c4_kind, &
                              l8_kind, l4_kind, &
                              i8_kind, i4_kind, i2_kind

   by use association. These are set to the appropriate KIND values that provide 4 and 8-byte versions of the fortran
   intrinsic types. The actual numerical value of the KIND may differ depending on platform.

   These should be used when the actually bytelength of the variable is important. To specify numerical precision,
   Fortran recommends the use of the intrinsics ``SELECTED_REAL_KIND`` and ``SELECTED_INT_KIND``.

Other other modules used
------------------------

.. container::

   None.

| 

Public interface
----------------

.. container::

   None.

| 

Public data
-----------

.. container::

   None.

| 

Public routines
---------------

.. container::

   None.

| 

Namelist
--------

.. container::

   None.

| 

Compiling and linking source
----------------------------

.. container::

   Any module or program unit using ``platform_mod`` must contain the line
   ::

      use platform_mod

Portability
-----------

.. container::

   None.

| 

Acquiring source
----------------

.. container::

   The ``platform_mod`` source consists of the main source file ``platform.F90`` and also requires the following include
   files:
   ::

      os.h

   GFDL users can check it out of the main CVS repository as part of the ``platform`` CVS module. The current public tag
   is ``fez``. Public access to the GFDL CVS repository will soon be made available.

| 

Notes
-----

.. container::

   None.

| 
