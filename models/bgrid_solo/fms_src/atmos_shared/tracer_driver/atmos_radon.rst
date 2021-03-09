module atmos_radon_mod
======================

Overview
--------

| This code allows the implementation of an extremely simplified radon tracer in the FMS framework.
| It should be taken as the implementation of a very simple tracer which bears some characteristics of radon.

.. container::

   This module presents an implementation of a tracer. It should be taken as representing radon only in a rudimentary
   manner.

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

Public interface
----------------

.. container::

   ::

      use atmos_radon_mod [, only:  atmos_radon_sourcesink,
                                    atmos_radon_init,
                                    atmos_radon_end ]

   atmos_radon_sourcesink:
      The routine that calculate the sources and sinks of radon.
   atmos_radon_init:
      The constructor routine for the radon module.
   atmos_radon_end:
      The destructor routine for the radon module.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Atmos_radon_sourcesink
      :name: atmos_radon_sourcesink

   ::

      call atmos_radon_sourcesink (lon, lat, land, pwt, radon, radon_dt, Time, kbot)

   **DESCRIPTION**
      This is a very rudimentary implementation of radon.
      It is assumed that the Rn222 flux is 3.69e-21 kg/m*m/sec over land for latitudes < 60N
      Between 60N and 70N the source = source \* .5
      Rn222 has a half-life time of 3.83 days, which corresponds to an e-folding time of 5.52 days.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lon``                                                   | Longitude of the centre of the model gridcells            |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lat``                                                   | Latitude of the centre of the model gridcells             |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``land``                                                  | Land/sea mask.                                            |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pwt``                                                   | The pressure weighting array. = dP/grav                   |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``radon``                                                 | The array of the radon mixing ratio.                      |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``Time``                                                  | Model time.                                               |
      |                                                           | [type(time_type)]                                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``kbot``                                                  | Integer array describing which model layer intercepts the |
      |                                                           | surface.                                                  |
      |                                                           | [integer, optional, dimension(:,:)]                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``radon_dt``                                              | The array of the tendency of the radon mixing ratio.      |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Atmos_radon_init
      :name: atmos_radon_init

   ::

      call atmos_radon_init 

   **DESCRIPTION**
      A routine to initialize the radon module.
   **INPUT**
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

c. .. rubric:: Atmos_radon_end
      :name: atmos_radon_end

   ::

      call atmos_radon_end 

   **DESCRIPTION**
      This subroutine writes the version name to logfile and exits.

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   None.

.. container::

   top
