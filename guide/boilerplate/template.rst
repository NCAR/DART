PROGRM OR MODULE name_of_thing
==============================

Overview
--------

Explain in general terms what this is.

Namelist
--------

DART namelists are always read from file ``input.nml``.

We adhere to the F90 standard of starting a namelist with an ampersand '&' and terminating with a slash '/' for all our
namelist input. Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.

::

   &NAMELIST_NML 
      name=value,
      name=value, 
      name=value
   /

Any comments about the namelist as a whole.

| 

.. container::

   ==== ==== ==============================
   Item Type Description
   ==== ==== ==============================
   name type (often multi-line) description
   ==== ==== ==============================

| 

Modules used
------------

::

   types_mod
   utilities_mod
   random_seq_mod
   time_manager_mod
   ensemble_manager_mod

Public interfaces
-----------------

================================== ==================
*use this_module_name_mod, only :* subr/function name
\                                  name2
\                                  name3
================================== ==================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *subroutine subroutine1(arg1, [, arg2])*
   ::

      real(r8),           intent(in) :: arg1
      real(r8), optional, intent(in) :: arg2

.. container:: indent1

   describe what this subroutine does.

   ======== =======================
   ``arg1`` Describe arg1.
   *arg2*   Describe optional arg2.
   ======== =======================

| 

.. container:: routine

   *function function1(arg1)*
   ::

      logical,             :: function1
      integer, intent(in)  :: arg1

.. container:: indent1

   Describe function.

   ============= ===================================
   ``function1`` describe what this function returns
   ``arg1``      describe function argument
   ============= ===================================

| 

.. container:: type

   ::

      type bob_type
         private
         integer :: bob1
         integer :: bob2
      end type bob_type

.. container:: indent1

   describe bob

   ========= ==============
   Component Description
   ========= ==============
   bob1      Describe bob1.
   bob2      Describe bob2.
   ========= ==============

| 

Files
-----

describe files used by code

References
----------

-  author, title. publication, volume, pages.
   `doi: nn.yyyy/rest_of_number <http://dx.doi.org/nn.yyyy/rest_of_number>`__

Private components
------------------

no discussion
