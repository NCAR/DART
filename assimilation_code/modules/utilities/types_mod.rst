MODULE types_mod
================

Overview
--------

Provides some commonly used mathematical constants, and a set of Fortran integer and real kinds, to be used to select
the right variable size (e.g. 4 bytes, 8 bytes) to match the rest of the DART interfaces. (DART does not depend on
compiler flags to set precision, but explicitly specifies a kind for each variable in the public interfaces.)

Other modules used
------------------

::

   none

Public interfaces
-----------------

| This routine provides the following constants, but no routines of any kind.
| The constants defined here *may* or *may not* be declared the same as constants used in non-DART pieces of code. It
  would seem like a good idea to match the DART definition of 'gas_constant' to the WRF equivalent if you are going to
  be running WRF/DART experiments (for example).

======================= ==============
*use types_mod, only :* i4
\                       i8
\                       r4
\                       r8
\                       c4
\                       c8
\                       digits12
\                       PI
\                       DEG2RAD
\                       RAD2DEG
\                       SECPERDAY
\                       MISSING_R4
\                       MISSING_R8
\                       MISSING_I
\                       MISSING_DATA
\                       metadatalength
\                       obstypelength
\                       t_kelvin
\                       es_alpha
\                       es_beta
\                       es_gamma
\                       gas_constant_v
\                       gas_constant
\                       L_over_Rv
\                       ps0
\                       earth_radius
\                       gravity
======================= ==============

| 

.. container:: type

   ::

      integer, parameter :: i4
      integer, parameter :: i8
      integer, parameter :: r4
      integer, parameter :: r8
      integer, parameter :: c4
      integer, parameter :: c8
      integer, parameter :: digits12

.. container:: indent1

   These kinds are used when declaring variables, like:

   ::

      real(r8)    :: myvariable
      integer(i4) :: shortint

   | All DART public interfaces use types on the real values to ensure they are consistent across various compilers and
     compile-time options. The digits12 is generally only used for reals which require extra precision.
   | Some models are able to run with single precision real values, which saves both memory when executing and file
     space when writing and reading restart files. To accomplish this, the users edit this file, redefine r8 to equal
     r4, and then rebuild all of DART.

| 

.. container:: type

   ::

      real(KIND=R8), parameter :: PI
      real(KIND=R8), parameter :: DEG2RAD
      real(KIND=R8), parameter :: RAD2DEG
      real(KIND=R8), parameter :: SECPERDAY

.. container:: indent1

   Some commonly used math constants, defined here for convenience.

| 

.. container:: type

   ::

      real(KIND=R4), parameter :: MISSING_R4
      real(KIND=R8), parameter :: MISSING_R8
      integer,       parameter :: MISSING_I
      integer,       parameter :: MISSING_DATA

.. container:: indent1

   Numeric constants used in the DART code when a numeric value is required, but the data is invalid or missing. These
   are typically defined as negative and a series of 8's, so they are distinctive when scanning a list of values.

| 

.. container:: type

   ::

      integer, parameter :: metadatalength
      integer, parameter :: obstypelength

.. container:: indent1

   Some common string limits used system-wide by DART code. The obstypelength is limited by the Fortran-imposed maximum
   number of characters in a parameter; the metadatalength was selected to be long enough to allow descriptive names but
   short enough to keep printing to less than a single line.

| 

.. container:: type

   ::

      real(KIND=R8), parameter :: t_kevin
      real(KIND=R8), parameter :: es_alpha
      real(KIND=R8), parameter :: es_beta
      real(KIND=R8), parameter :: es_gamma
      real(KIND=R8), parameter :: gas_constant_v
      real(KIND=R8), parameter :: gas_constant
      real(KIND=R8), parameter :: L_over_Rv
      real(KIND=R8), parameter :: ps0
      real(KIND=R8), parameter :: earth_radius
      real(KIND=R8), parameter :: gravity

.. container:: indent1

   | A set of geophysical constants, which could be argued do not belong in a DART-supplied file since they are quite
     probably specific to a model or a particular forward operator.
   | Best case would be if we could engineer the code so these constants were provided by the model and then used when
     compiling the forward operator files. But given that Fortran use statements cannot be circular, this poses a
     problem. Perhaps we could work out how the obs_def code could define these constants and then they could be used by
     the model code. For now, they are defined here but it is up to the model and obs_def code writers whether to use
     these or not.

| 

Namelist
--------

There is no namelist for this module.

Files
-----

None.

References
----------

#. none

Private components
------------------

N/A
