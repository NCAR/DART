Module sat_vapor_pres_mod
=========================

Overview
--------

Routines for determining the saturation vapor pressure (``ES``) and the derivative of ``ES`` with respect to
temperature.

.. container::

   This module contains routines for determining the saturation vapor pressure (``ES``) from lookup tables constructed
   using equations given in the Smithsonian tables. The ``ES`` lookup tables are valid between -160C and +100C (approx
   113K to 373K). The values of ``ES`` are computed over ice from -160C to -20C, over water from 0C to 100C, and a
   blended value (over water and ice) from -20C to 0C. This version was written for non-vector machines. See the notes
   section for details on vectorization.

| 

Other modules used
------------------

.. container::

   ::

      constants_mod
            fms_mod

Public interface
----------------

.. container::

   Description summarizing public interface.

   lookup_es:
      For the given temperatures, returns the saturation vapor pressures.
   lookup_des:
      For the given temperatures, returns the derivative of saturation vapor pressure with respect to temperature.
   compute_es:
      For the given temperatures, computes the saturation vapor pressures.
   sat_vapor_pres_init:
      Initializes the lookup tables for saturation vapor pressure.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Lookup_es
      :name: lookup_es

   ::

      call lookup_es ( temp, esat )

   **DESCRIPTION**
      For the given temperatures these routines return the saturation vapor pressure (esat). The return values are
      derived from lookup tables (see notes below).
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``temp``                                                  | Temperature in degrees Kelvin.                            |
      |                                                           | [real, dimension(scalar)]                                 |
      |                                                           | [real, dimension(:)]                                      |
      |                                                           | [real, dimension(:,:)]                                    |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``esat``                                                  | Saturation vapor pressure in pascals. May be a scalar,    |
      |                                                           | 1d, 2d, or 3d array. Must have the same order and size as |
      |                                                           | temp.                                                     |
      |                                                           | [real, dimension(scalar)]                                 |
      |                                                           | [real, dimension(:)]                                      |
      |                                                           | [real, dimension(:,:)]                                    |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Lookup_des
      :name: lookup_des

   ::

      call lookup_des ( temp, desat )

   **DESCRIPTION**
      For the given temperatures these routines return the derivative of esat w.r.t. temperature (desat). The return
      values are derived from lookup tables (see notes below).
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``temp``                                                  | Temperature in degrees Kelvin.                            |
      |                                                           | [real, dimension(scalar)]                                 |
      |                                                           | [real, dimension(:)]                                      |
      |                                                           | [real, dimension(:,:)]                                    |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``desat``                                                 | Derivative of saturation vapor pressure w.r.t.            |
      |                                                           | temperature in pascals/degree. May be a scalar, 1d, 2d,   |
      |                                                           | or 3d array. Must have the same order and size as temp.   |
      |                                                           | [real, dimension(scalar)]                                 |
      |                                                           | [real, dimension(:)]                                      |
      |                                                           | [real, dimension(:,:)]                                    |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Compute_es
      :name: compute_es

   ::

      es = compute_es ( temp )

   **DESCRIPTION**
      Computes saturation vapor pressure for the given temperature using the equations given in the Smithsonian
      Meteorological Tables. Between -20C and 0C a blended value over ice and water is returned.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``temp``                                                  | Temperature in degrees Kelvin.                            |
      |                                                           | [real, dimension(:)]                                      |
      |                                                           | [real, dimension(scalar)]                                 |
      |                                                           | [real, dimension(:,:)]                                    |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``es``                                                    | Saturation vapor pressure in pascals. May be a scalar,    |
      |                                                           | 1d, 2d, or 3d array. Must have the same order and size as |
      |                                                           | temp.                                                     |
      |                                                           | [real, dimension(:)]                                      |
      |                                                           | [real, dimension(scalar)]                                 |
      |                                                           | [real, dimension(:,:)]                                    |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Sat_vapor_pres_init
      :name: sat_vapor_pres_init

   ::

      call sat_vapor_pres_init 

   **DESCRIPTION**
      Initializes the lookup tables for saturation vapor pressure. This routine will be called automatically the first
      time **lookup_es** or **lookup_des** is called, the user does not need to call this routine. There are no
      arguments.

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   **FATAL in lookup_es**
      table overflow, nbad=##
      Temperature(s) provided to the saturation vapor pressure lookup are outside the valid range of the lookup table
      (-160 to 100 deg C). This may be due to a numerical instability in the model. Information should have been printed
      to standard output to help determine where the instability may have occurred. If the lookup table needs a larger
      temperature range, then parameters in the module header must be modified.

References
----------

.. container::

   #. Smithsonian Meteorological Tables Page 350.

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

   test_sat_vapor_pres
      ::

         use sat_vapor_pres_mod
         implicit none

         integer, parameter :: ipts=500, jpts=100, kpts=50, nloop=1
         real, dimension(ipts,jpts,kpts) :: t,es,esn,des,desn
         integer :: n

          generate temperatures between 120K and 340K
           call random_number (t)
           t = 130. + t * 200.

          initialize the tables (optional)
           call sat_vapor_pres_init

          compute actual es and "almost" actual des
            es = compute_es  (t)
           des = compute_des (t)

         do n = 1, nloop
          es and des
           call lookup_es  (t, esn)
           call lookup_des (t,desn)
         enddo

          terminate, print deviation from actual
           print *, 'size=',ipts,jpts,kpts,nloop
           print *, 'err es  = ', sum((esn-es)**2)
           print *, 'err des = ', sum((desn-des)**2)

         contains

         ----------------------------------
          routine to estimate derivative

          function compute_des (tem) result (des)
          real, intent(in) :: tem(:,:,:)
          real, dimension(size(tem,1),size(tem,2),size(tem,3)) :: des,esp,esm
          real, parameter :: tdel = .01
             esp = compute_es (tem+tdel)
             esm = compute_es (tem-tdel)
             des = (esp-esm)/(2*tdel)
          end function compute_des
         ----------------------------------

         end program test_sat_vapor_pres

| 

Notes
-----

.. container::

   1. **Vectorization**
   To create a vector version the lookup routines need to be modified. The local variables: tmp, del, ind, should be
   changed to arrays with the same size and order as input array temp.
   2. **Construction of the ``ES`` tables**
   The tables are constructed using the saturation vapor pressure (``ES``) equations in the Smithsonian tables. The
   tables are valid between -160C to +100C with increments at 1/10 degree. Between -160C and -20C values of ``ES`` over
   ice are used, between 0C and 100C values of ``ES`` over water are used, between -20C and 0C blended values of ``ES``
   (over water and over ice) are used.
   There are three tables constructed: ``ES``, first derivative (``ES'``), and second derivative (``ES``''). The ES
   table is constructed directly from the equations in the Smithsonian tables. The ``ES``' table is constructed by
   bracketing temperature values at +/- 0.01 degrees. The ``ES``'' table is estimated by using centered differencing of
   the ``ES``' table.
   3. **Determination of ``es`` and ``es'`` from lookup tables**
   Values of the saturation vapor pressure (``es``) and the derivative (``es'``) are determined at temperature (T) from
   the lookup tables (``ES``, ``ES'``, ``ES''``) using the following formula.
   ::

          es (T) = ES(t) + ES'(t) * dt + 0.5 * ES''(t) * dt**2
          es'(T) = ES'(t) + ES''(t) * dt

          where     t = lookup table temperature closest to T
                   dt = T - t

   4. Internal (private) parameters
   These parameters can be modified to increase/decrease the size/range of the lookup tables.
   ::

      tcmin   The minimum temperature (in deg C) in the lookup tables.
                    [integer, default: tcmin = -160]

          tcmax   The maximum temperature (in deg C) in the lookup tables.
                    [integer, default: tcmin = +100]

| 
