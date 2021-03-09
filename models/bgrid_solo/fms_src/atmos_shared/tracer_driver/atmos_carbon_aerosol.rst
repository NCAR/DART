module atmos_carbon_aerosol_mod
===============================

Overview
--------

This code allows the implementation of black and organic carbon tracers in the FMS framework.

.. container::

   This module presents the method of Cooke et al. (1999, 2002) In its present implementation the black and organic
   carbon tracers are from the combustion of fossil fuel. While the code here should provide insights into the
   carbonaceous aerosol cycle it is provided here more as an example of how to implement a tracer module in the FMS
   infrastructure. The parameters of the model should be checked and set to values corresponding to previous works if a
   user wishes to try to reproduce those works.

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

      use atmos_carbon_aerosol_mod [, only:  atmos_blackc_sourcesink,
                                             atmos_organic_sourcesink,
                                             atmos_carbon_aerosol_init,
                                             atmos_carbon_aerosol_end ]

   atmos_blackc_sourcesink:
      A subroutine to calculate the source and sinks of black carbon aerosol.
   atmos_organic_sourcesink:
      A subroutine to calculate the source and sinks of organic carbon aerosol.
   atmos_carbon_aerosol_init:
      Subroutine to initialize the carbon aerosol module.
   atmos_carbon_aerosol_end:
      The destructor routine for the carbon aerosol module.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Atmos_blackc_sourcesink
      :name: atmos_blackc_sourcesink

   ::

      call atmos_blackc_sourcesink (lon, lat, land, pwt, & black_cphob, black_cphob_dt, & black_cphil, black_cphil_dt, & Time, is, ie, js, je, kbot)

   **DESCRIPTION**
      | This routine calculates the source and sink terms for black carbon. Simply put, the hydrophobic aerosol has
        sources from emissions and sinks from dry deposition and transformation into hydrophilic aerosol. The
        hydrophilic aerosol also has emission sources and has sinks of wet and dry deposition.
      | The following schematic shows how the black carbon scheme is implemented.

      ::

          +------------+  Trans-   +------------+
          |  Hydro-    | formation |  Hydro-    |
          |  phobic    |           |  philic    |
          |  black     |---------->|  black     |
          |  carbon    |           |  carbon    |
          |            |           |            |
          +------------+           +------------+
             ^      |                ^    |   |
             |      |                |    |   |
             |      =                |    =   =
           Source  Dry            Source Dry Wet
                   Dep.                  Dep Dep

      The transformation time used here is 1 day, which corresponds to an e-folding time of 1.44 days. This can be
      varied as necessary.

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
      | ``black_cphob``                                           | The array of the hydrophobic black carbon aerosol mixing  |
      |                                                           | ratio                                                     |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``black_cphil``                                           | The array of the hydrophilic black carbon aerosol mixing  |
      |                                                           | ratio                                                     |
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
      | ``black_cphob_dt``                                        | The array of the tendency of the hydrophobic black carbon |
      |                                                           | aerosol mixing ratio.                                     |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``black_cphil_dt``                                        | The array of the tendency of the hydrophilic black carbon |
      |                                                           | aerosol mixing ratio.                                     |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Atmos_organic_sourcesink
      :name: atmos_organic_sourcesink

   ::

      call atmos_organic_sourcesink (lon, lat, land, pwt, organic_carbon, organic_carbon_dt, & Time, is, ie, js, je, kbot)

   **DESCRIPTION**
      | This routine calculates the source and sink terms for organic carbon. Simply put, the hydrophobic aerosol has
        sources from emissions and sinks from dry deposition and transformation into hydrophilic aerosol. The
        hydrophilic aerosol also has emission sources and has sinks of wet and dry deposition.
      | The following schematic shows how the organic carbon scheme is implemented.

      ::

          +------------+  Trans-   +------------+
          |  Hydro-    | formation |  Hydro-    |
          |  phobic    |           |  philic    |
          |  organic   |---------->|  organic   |
          |  carbon    |           |  carbon    |
          |            |           |            |
          +------------+           +------------+
             ^      |                ^    |   |
             |      |                |    |   |
             |      =                |    =   =
           Source  Dry            Source Dry Wet
                   Dep.                  Dep Dep

      | The transformation time used here is 2 days, which corresponds to an e-folding time of 2.88 days. This can be
        varied as necessary.

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
      | ``organic_carbon``                                        | The array of the organic carbon aerosol mixing ratio      |
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
      | ``organic_carbon_dt``                                     | The array of the tendency of the organic carbon aerosol   |
      |                                                           | mixing ratio.                                             |
      |                                                           | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Atmos_carbon_aerosol_init
      :name: atmos_carbon_aerosol_init

   ::

      call atmos_carbon_aerosol_init (lonb, latb, r, axes, Time, mask)

   **DESCRIPTION**
      This subroutine querys the tracer manager to find the indices for the various carbonaceous aerosol tracers. It
      also registers the emission fields for diagnostic purposes.
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

d. .. rubric:: Atmos_carbon_aerosol_end
      :name: atmos_carbon_aerosol_end

   ::

      call atmos_carbon_aerosol_end 

   **DESCRIPTION**
      This subroutine writes the version name to logfile and exits.

Data sets
---------

.. container::

   Black carbon emissions
      The black carbon emission dataset is that derived in Cooke et al. (1999) The dataset can be obtained from the
      contact person above.
   Organic carbon emissions
      The organic carbon emission dataset is that derived in Cooke et al. (1999) The dataset can be obtained from the
      contact person above.

Error messages
--------------

.. container::

   None.

References
----------

.. container::

   #. Cooke, W. F. and J. J. N. Wilson, A global black carbon aerosol model, J. Geophys. Res., 101, 19395-19409, 1996.
   #. Cooke, W. F., C. Liousse, H. Cachier and J. Feichter, Construction of a 1 x 1 fossil fuel emission dataset for
      carbonaceous aerosol and implementation and radiative impact in the ECHAM-4 model, J. Geophys. Res., 104,
      22137-22162, 1999
   #. Cooke, W.F., V. Ramaswamy and P. Kasibathla, A GCM study of the global carbonaceous aerosol distribution. J.
      Geophys. Res., 107, accepted, 2002

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
