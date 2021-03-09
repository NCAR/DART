module atmos_sulfur_hex_mod
===========================

Overview
--------

This code allows the implementation of sulfur hexafluoride tracer in the FMS framework.

| 

Other modules used
------------------

.. container::

   ::

      fms_mod
          time_manager_mod
          diag_manager_mod
        tracer_manager_mod
         field_manager_mod
      tracer_utilities_mod
             constants_mod

Public interface
----------------

.. container::

   ::

      use atmos_sulfur_hex_mod [, only:  atmos_sf6_sourcesink,
                                         atmos_sulfur_hex_init,
                                         sulfur_hex_end ]

   atmos_sf6_sourcesink:
      A routine to calculate the sources and sinks of sulfur hexafluoride.
   atmos_sulfur_hex_init:
      The constructor routine for the sulfur hexafluoride module.
   sulfur_hex_end:
      The destructor routine for the sulfur hexafluoride module.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Atmos_sf6_sourcesink
      :name: atmos_sf6_sourcesink

   ::

      call atmos_sf6_sourcesink (lon, lat, land, pwt, sf6, sf6_dt, Time, is, ie, js, je, kbot)

   **DESCRIPTION**
      A routine to calculate the sources and sinks of sulfur hexafluoride.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lon``                                                   | Longitude of the centre of the model gridcells.           |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lat``                                                   | Latitude of the centre of the model gridcells.            |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``land``                                                  | Land/sea mask.                                            |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pwt``                                                   | The pressure weighting array. = dP/grav                   |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``sf6``                                                   | The array of the sulfur hexafluoride mixing ratio.        |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``Time``                                                  | Model time.                                               |
      |                                                           | [type(time_type)]                                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``is, ie, js, je``                                        | Local domain boundaries.                                  |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``kbot``                                                  | Integer array describing which model layer intercepts the |
      |                                                           | surface.                                                  |
      |                                                           | [integer, optional, dimension(:,:)]                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``sf6_dt``                                                | The array of the tendency of the sulfur hexafluoride      |
      |                                                           | mixing ratio.                                             |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Atmos_sulfur_hex_init
      :name: atmos_sulfur_hex_init

   ::

      call atmos_sulfur_hex_init (lonb, latb, r, axes, Time, mask)

   **DESCRIPTION**
      A routine to initialize the sulfur hexafluoride module.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lonb``                                                  | The longitudes for the local domain.                      |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``latb``                                                  | The latitudes for the local domain.                       |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``mask``                                                  | optional mask (0. or 1.) that designates which grid       |
      |                                                           | points are above (=1.) or below (=0.) the ground          |
      |                                                           | dimensioned as (nlon,nlat,nlev).                          |
      |                                                           | [real, optional, dimension(:,:,:)]                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``Time``                                                  | Model time.                                               |
      |                                                           | [type(time_type)]                                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``axes``                                                  | The axes relating to the tracer array dimensioned as      |
      |                                                           | (nlon, nlat, nlev, ntime)                                 |
      |                                                           | [integer, dimension(4)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **INPUT/OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``r``                                                     | Tracer fields dimensioned as (nlon,nlat,nlev,ntrace).     |
      |                                                           | [real, dimension(:,:,:,:)]                                |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Sulfur_hex_end
      :name: sulfur_hex_end

   ::

      call sulfur_hex_end 

   **DESCRIPTION**
      This subroutine is the exit routine for the sulfur hexafluoride module.

Data sets
---------

.. container::

   Sulfur hexaflouride emissions
      Monthly.emissions contains the estimated global emission rate of SF6 in Gg/yr for 62 months between December 1988
      and January 1994, inclusive. These are based on the annual estimates of Levin and Hesshaimer (submitted), and have
      been linearly interpolated to monthly values. The last half of 1993 has been extrapolated using the trend for the
      previous 12 months.
      The dataset can be obtained from the contact person above.

Error messages
--------------

.. container::

   None.

References
----------

.. container::

   #. Levin, I. and V. Hessahimer: Refining of atmospheric transport model entries by the globally observed passive
      tracer distributions of 85Krypton and Sulfur Hexafluoride (SF6). Submitted to the Journal of Geophysical Research.

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

   None.

| 
