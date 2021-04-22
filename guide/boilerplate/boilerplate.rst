PROGRAM OR MODULE name_of_thing
===============================

Overview
--------

This is the place for the general description of the module or program or whatever.

overview issues, how/when/why to use this thing, etc.

more stuff about code, usage, etc.

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

| 

.. container::

   ==== ==== ==============================
   Item Type Description
   ==== ==== ==============================
   name type (often multi-line) description
   ==== ==== ==============================

| 

Other modules used
------------------

::

   types_mod
   model_mod
   assim_tools_mod
   time_manager_mod
   fms_mod

Public interfaces
-----------------

===================== ===========
*use xxxxxxx, only :* yyypubtype1
\                     yyyroutine1
\                     yyyroutine2
\                     yyyroutine3
===================== ===========

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   ::

      type location_type
         private
         real(r8) :: x
      end type location_type

.. container:: indent1

   The location type is essential to the fabric of the universe. If we don't know where we are, how do we know which
   exit to take?

   ========= ====================================
   Component Description
   ========= ====================================
   x         is the nondimensional distance [0,1]
   ========= ====================================

| 

.. container:: routine

   *subroutine yyyroutine1( var1, var2, var3 [, global])*
   ::

      type(time_type),          intent(in)    ::  var1 
      real(r8), dimension(:),   intent(inout) ::  var2 
      real(r8), dimension(:,:), intent(out)   ::  var3 
      real(r8), optional,       intent(in)    ::  global 

.. container:: indent1

   | Returns the resolution of compute domain for either the current processor or the global domain.
   | This is the best thing since sliced bread. All you have to do is throw some arguments in the call and the
     subroutine automatically slices, dices and makes julienne fries. *But wait!* There's more!

   +----------+----------------------------------------------------------------------------------------------------------+
   | ``var1`` | is the number of spark plugs in a '67 Cuda.                                                              |
   +----------+----------------------------------------------------------------------------------------------------------+
   | ``var2`` | is used for both input and output. Upon input, this contains the furlongs per weekday for every vertical |
   |          | level. Upon exit, we now have the number of spotted gobies per square hectare at that level.             |
   +----------+----------------------------------------------------------------------------------------------------------+
   | ``var3`` | The leftmost dimension pertains to the number of feet in an orange-footed brush fowl, indiginous to      |
   |          | Australia. The next dimension is the number of feathers on said fowl, naturally.                         |
   +----------+----------------------------------------------------------------------------------------------------------+
   | *global* | is really, REALLY useful in certain situations.                                                          |
   +----------+----------------------------------------------------------------------------------------------------------+

   notes would go here

| 

.. container:: routine

   *function yyyroutine2( var1, var2, var3, [,bob] )*
   ::

      type(time_type),          intent(in)  ::  var1 
      real(r8), dimension(:),   intent(in)  ::  var2 
      real(r8), dimension(:,:), intent(in)  ::  var3 
      real(r8), optional,       intent(in)  ::  bob 
      integer, dimension(size(var2))        ::  yyyroutine2 

.. container:: indent1

   | Returns the resolution of compute domain for either the current processor or the global domain. All input variables
     are not changed. Otherwise, this would be a subroutine.
   | This is the second-best thing since sliced bread. All you have to do is throw some arguments in the call and the
     function automatically deep fries.

   +----------------+----------------------------------------------------------------------------------------------------+
   | ``var1``       | the first time you changed your oil.                                                               |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``var2``       | miles between every oil change you've ever done. Don't lie.                                        |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``var3``       | the distances you've ridden. Each row corresponds to the hour-of-day, each column is a different   |
   |                | day-of-the-week.                                                                                   |
   +----------------+----------------------------------------------------------------------------------------------------+
   | *bob*          | mean time between failures. in msec.                                                               |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``yyroutine2`` | number of gray hairs as a function of time. in kilohairs.                                          |
   +----------------+----------------------------------------------------------------------------------------------------+

   notes would go here

| 

Files
-----

This is the place to discuss the files that are associated with this module. They could be input files, output files,
data files, shell scripts ... anything.

=========================== ===========================================================================
filename                    purpose
=========================== ===========================================================================
inputfile1                  to read some input
input.nml                   to read namelists
preassim.nc                 the time-history of the model state before assimilation
analysis.nc                 the time-history of the model state after assimilation
dart_log.out [default name] the run-time diagnostic output
dart_log.nml [default name] the record of all the namelists actually USED - contains the default values
=========================== ===========================================================================

References
----------

-  Anderson, J., T. Hoar, K. Raeder, H. Liu, N. Collins, R. Torn, and A. Arellano, 2009:
   The Data Assimilation Research Testbed: A Community Facility. Bull. Amer. Meteor. Soc., 90, 1283-1296.
   `DOI: 10.1175/2009BAMS2618.1 <http://dx.doi.org/10.1175%2F2009BAMS2618.1>`__

-  none

Private components
------------------

N/A

Any routines or 'local' variables of interest may be discussed here. There are generally lots of 'internal' functions
that make life simpler, but you don't want to make them available outside the scope of the current module. This is the
place to point them out, if you like.

.. container:: routine

   ::

      type location_type
         private
         real(r8) :: x
      end type location_type
