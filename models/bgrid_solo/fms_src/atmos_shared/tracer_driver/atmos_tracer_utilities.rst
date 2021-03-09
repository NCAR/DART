module atmos_tracer_utilities_mod
=================================

Overview
--------

This code provides some utility routines for atmospheric tracers in the FMS framework.

.. container::

   This module gives utility routines which can be used to provide consistent removal mechanisms for atmospheric
   tracers.
   In particular it provides schemes for wet and dry deposiiton that can be easily utilized.

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
           constants_mod
        horiz_interp_mod

Public interface
----------------

.. container::

   ::

      use atmos_tracer_utilities_mod [, only:  atmos_tracer_utilities_init,
                                               dry_deposition,
                                               wet_deposition,
                                               interp_emiss,
                                               tracer_utilities_end ]

   atmos_tracer_utilities_init:
      This is a routine to create and register the dry and wet deposition fields of the tracers.
   dry_deposition:
      Routine to calculate the fraction of tracer to be removed by dry deposition.
   wet_deposition:
      Routine to calculate the fraction of tracer removed by wet deposition
   interp_emiss:
      A routine to interpolate emission fields of arbitrary resolution onto the resolution of the model.
   tracer_utilities_end:
      The destructor routine for the tracer utilities module.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Atmos_tracer_utilities_init
      :name: atmos_tracer_utilities_init

   ::

      call atmos_tracer_utilities_init (lonb,latb, mass_axes, Time)

   **DESCRIPTION**
      This routine creates diagnostic names for dry and wet deposition fields of the tracers. It takes the tracer name
      and appends "ddep" for the dry deposition field and "wdep" for the wet deposition field. This names can then be
      entered in the diag_table for diagnostic output of the tracer dry and wet deposition. The module name associated
      with these fields in "tracers". The units of the deposition fields are assumed to be kg/m2/s.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lonb``                                                  | The longitudes for the local domain.                      |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``latb``                                                  | The latitudes for the local domain.                       |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``mass_axes``                                             | The axes relating to the tracer array.                    |
      |                                                           | [integer, dimension(3)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``Time``                                                  | Model time.                                               |
      |                                                           | [type(time_type)]                                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Dry_deposition
      :name: dry_deposition

   ::

      call dry_deposition ( n, is, js, u, v, T, pwt, pfull, u_star, landmask, dsinku, tracer, Time)

   **DESCRIPTION**
      | There are two types of dry deposition coded.
      | 1) Wind driven derived dry deposition velocity.
      | 2) Fixed dry deposition velocity.
      | The theory behind the wind driven dry deposition velocity calculation assumes that the deposition can be modeled
        as a parallel resistance type problem.
      | Total resistance to HNO3-type dry deposition,

      ::

                R = Ra + Rb
           resisa = aerodynamic resistance
           resisb = surface resistance (laminar layer + uptake)
                  = 5/u*  [s/cm]        for neutral stability
               Vd = 1/R

      | For the fixed dry deposition velocity, there is no change in the deposition velocity but the variation of the
        depth of the surface layer implies that there is variation in the amount deposited.
      | To utilize this section of code add one of the following lines as a method for the tracer of interest in the
        field table.

      ::

          "dry_deposition","wind_driven","surfr=XXX"
              where XXX is the total resistance defined above.

          "dry_deposition","fixed","land=XXX, sea=YYY"
              where XXX is the dry deposition velocity (m/s) over land
                and YYY is the dry deposition velocity (m/s) over sea.

   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | The tracer number.                                        |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``is, js``                                                | Start indices for array (computational indices).          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``u``                                                     | U wind field.                                             |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``v``                                                     | V wind field.                                             |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``T``                                                     | Temperature.                                              |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pwt``                                                   | Pressure differential of half levels.                     |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pfull``                                                 | Full pressure levels.                                     |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``u_star``                                                | Friction velocity.                                        |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``landmask``                                              | Land - sea mask.                                          |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``dsinku``                                                | The amount of tracer in the surface layer which is dry    |
      |                                                           | deposited per second.                                     |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Wet_deposition
      :name: wet_deposition

   ::

      call wet_deposition (n, T, pfull, phalf, rain, snow, qdt, tracer, tracer_dt, Time, cloud_param, is, js)

   **DESCRIPTION**
      | Schemes allowed here are
      | 1) Deposition removed in the same fractional amount as the modeled precipitation rate is to a standardized
        precipitation rate. Basically this scheme assumes that a fractional area of the gridbox is affected by
        precipitation and that this precipitation rate is due to a cloud of standardized cloud liquid water content.
        Removal is constant throughout the column where precipitation is occuring.
      | 2) Removal according to Henry's Law. This law states that the ratio of the concentation in cloud water and the
        partial pressure in the interstitial air is a constant. In this instance, the units for Henry's constant are
        kg/L/Pa (normally it is M/L/Pa) Parameters for a large number of species can be found at
        http://www.mpch-mainz.mpg.de/~sander/res/henry.html To utilize this section of code add one of the following
        lines as a method for the tracer of interest in the field table.

      ::

          "wet_deposition","henry","henry=XXX, dependence=YYY"
              where XXX is the Henry's constant for the tracer in question
                and YYY is the temperature dependence of the Henry's Law constant.

          "wet_deposition","fraction","lslwc=XXX, convlwc=YYY"
              where XXX is the liquid water content of a standard large scale cloud
                and YYY is the liquid water content of a standard convective cloud.

   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Tracer number                                             |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``is, js``                                                | start indices for array (computational indices)           |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``T``                                                     | Temperature                                               |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pfull``                                                 | Full level pressure field                                 |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``phalf``                                                 | Half level pressure field                                 |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``rain``                                                  | Precipitation in the form of rain                         |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``snow``                                                  | Precipitation in the form of snow                         |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``qdt``                                                   | The tendency of the specific humidity due to the cloud    |
      |                                                           | parametrization                                           |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``tracer``                                                | The tracer field                                          |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``Time``                                                  | The time structure for submitting wet deposition as a     |
      |                                                           | diagnostic                                                |
      |                                                           | [type(time_type)]                                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``cloud_param``                                           | Is this a convective (convect) or large scale (lscale)    |
      |                                                           | cloud parametrization?                                    |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``tracer_dt``                                             | The tendency of the tracer field due to wet deposition.   |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Interp_emiss
      :name: interp_emiss

   ::

      call interp_emiss (global_source, start_lon, start_lat, & lon_resol, lat_resol, data_out)

   **DESCRIPTION**
      Routine to interpolate emission fields (or any 2D field) to the model resolution. The local section of the global
      field is returned to the local processor.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``global_source``                                         | Global emission field.                                    |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``start_lon``                                             | Longitude of starting point of emission field (in         |
      |                                                           | radians). This is the westernmost boundary of the global  |
      |                                                           | field.                                                    |
      |                                                           | [real]                                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``start_lat``                                             | Latitude of starting point of emission field (in          |
      |                                                           | radians). This is the southern boundary of the global     |
      |                                                           | field.                                                    |
      |                                                           | [real]                                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lon_resol``                                             | Longitudinal resolution of the emission data (in          |
      |                                                           | radians).                                                 |
      |                                                           | [real]                                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``lat_resol``                                             | Latitudinal resolution of the emission data (in radians). |
      |                                                           | [real]                                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``data_out``                                              | Interpolated emission field on the local PE.              |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

e. .. rubric:: Tracer_utilities_end
      :name: tracer_utilities_end

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
