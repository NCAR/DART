MOM6
==============

A new ocean component model based on the Modular Ocean Model version 6 (MOM6) has been incorporated into
`CESM <https://www.cesm.ucar.edu/>`_ and is anticipated to replace POP2 as the default ocean component in CESM3.
An early functional release of the MOM6 ocean component has been made available to users beginning with CESM2.2.
Instructions for using MOM6 in CESM are available on the `MOM_interface GitHub Wiki
<https://github.com/ESCOMP/MOM_interface/wiki>`_.

This DART-MOM6 interface was developed for `MOM6 <https://github.com/NCAR/MOM6>`_ within the CESM framework.

MOM6 time
---------

The default in CESM is to run with no leap years.
To assimilate real observations, we need to switch to the Gregorian 
calendar to account for leap years.

.. code-block:: text

    ./xmlchange CALENDAR=GREGORIAN

To illustrate what happens if you do not set CALENDAR=GREGORIAN, here is
an example where the RUN_STARTDATE is set to 2015-02-01 and MOM6 is run for 10 days.

.. code-block:: text

    ./xmlchange RUN_STARTDATE=2015-02-01

The MOM6 restart file has the following meta data, where Time is days from year 1.

.. code-block:: text

    double Time(Time) ;
                    Time:long_name = "Time" ;
                    Time:units = "days" ;
                    Time:axis = "T" ;
    ...
    // global attributes:
    		:filename = "./c.T62_g16.ens3.mom6.r.2015-02-11-00000._0001.nc" ;
    data:
    
     Time = 735151 ;
    }


The absence of leap years gives you inconsistent time information when comparing 
to observation times in YYYY-MM-DD:

- Restart filename has the time 2015-02-11-00000.
- The Time variable is Time = 735151 days, which is 2013/10/11


MOM6 checksum of restart files
------------------------------

When reading in restart files, MOM6 verifies a checksum for each variable
in the restart file. Data assimilation updates the data in the MOM6 restart file,
which will cause the checksum verification to fail.  To use DART-MOM6 with CESM
turn off the checksum verification using the ``user_nml_mom`` namelist option:


.. code-block:: text

    RESTART_CHECKSUMS_REQUIRED = False

model_nml
---------

The namelist options for DART-MOM6 are as follows:

.. code-block:: text

    &model_nml
       template_file = 'mom6.r.nc',
       ocean_geometry = 'ocean_geometry.nc',
       static_file = 'c.e22.GMOM.T62_g16.nuopc.001.mom6.static.nc',
       model_state_variables        = 'Salt ', 'QTY_SALINITY             ', 'UPDATE',
                                      'Temp ', 'QTY_POTENTIAL_TEMPERATURE', 'UPDATE',
                                      'u    ', 'QTY_U_CURRENT_COMPONENT  ', 'UPDATE',
                                      'v    ', 'QTY_V_CURRENT_COMPONENT  ', 'UPDATE',
                                      'h    ', 'QTY_LAYER_THICKNESS      ', 'UPDATE',
       assimilation_period_days     = 1
       assimilation_period_seconds  = 1
       /

* ``template_file`` is a MOM6 restart file. The size and shape of the state variables will be read from this netCDF file.

* ``ocean_geometry`` is a MOM6 netCDF file containing the variables ``D``, the basin depth in meters; and ``wet``, whether a point is land or ocean at the Earth's surface.

* ``static_file`` is a MOM6 netCDF file containing the grid information. The following three grids are read into DART:

    .. code-block:: text

         geolon(:,:)   Longitude of tracer (T) points
         geolat(:,:)   Latitude of tracer (T) points
         
         geolon_u(:,:) Longitude of zonal velocity (Cu) points
         geolat_u(:,:) Latitude of zonal velocity (Cu) points
         
         geolon_v(:,:) Longitude of meridional velocity (Cv) points
         geolat_v(:,:) Latitude of meridional velocity (Cv) point


Vertical Coordinate
-------------------

The vertical coordinate in MOM6 is layer thickness which varies across the ensemble.
To get the depth in meters at a particular layer, you must sum the layer thicknesses.

.. code-block:: text

    double h(Time, Layer, lath, lonh) ;
           h:long_name = "Layer Thickness" ;
           h:units = "m" ;

.. Note

   Layer interface thickness maybe available from MOM6. But the restarts we have
   available have "Layer thickness" only.


Land in the state
------------------

The MOM6 grid is global, so land is included in the state. To avoid updating land,
``get_close_state`` forces the distance to be very large for dry land and locations
below the sea floor.

Identifying land/ocean points at the surface can be simply done using the the longitude, latitude
of the given point. ``get_state_meta_data`` is used to assign land points at the surface QTY_DRY_LAND.

The process to identify points below the sea floor requires the vertical location of the
point in meters. The conversion from model layer to depth in meters is done in
``convert_vertical_state``. The depth is then used to identify points below the
basin depth in ``get_close_state``.


.. code-block:: fortran
   :emphasize-lines: 5, 9
   :caption: snippet from get_close_state

    ! Put any land or sea floor points very far away
    ! so they are not updated by assimilation
    do ii = 1, num_close
    
      if(loc_qtys(close_ind(ii)) == QTY_DRY_LAND) dist = 1.0e9_r8
    
      lon_lat_vert = get_location(locs(close_ind(ii))) ! assuming VERTISHEIGHT
      call get_model_variable_indices(loc_indx(ii), i, j, k)
      if ( below_sea_floor(i,j,lon_lat_vert(3)) ) dist = 1.0e9_r8
    
    enddo



