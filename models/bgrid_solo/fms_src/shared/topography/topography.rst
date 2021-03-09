module topography_mod
=====================

Overview
--------

Routines for creating land surface topography fields and land-water masks for latitude-longitude grids.

.. container::

   This module generates realistic mountains and land-water masks on a specified latitude-longitude grid by
   interpolating from the 1/6 degree Navy mean topography and percent water data sets. The fields that can be generated
   are mean and standard deviation of topography within the specified grid boxes; and land-ocean (or water) mask and
   land-ocean (or water) fractional area.
   The interpolation scheme conserves the area-weighted average of the input data by using module horiz_interp.
   The interfaces get_gaussian_topog and gaussian_topog_init are documented in :doc:`./gaussian_topog`.

| 

Other modules used
------------------

.. container::

   ::

      gaussian_topog_mod
        horiz_interp_mod
                 fms_mod

Public interface
----------------

.. container::

   ::

      use topography_mod [, only:  get_topog_mean,
                                   get_topog_stdev,
                                   get_ocean_frac,
                                   get_ocean_mask,
                                   get_water_frac,
                                   get_water_mask ]

   get_topog_mean:
      Returns a "realistic" mean surface height field.
   get_topog_stdev:
      Returns a standard deviation of higher resolution topography with the given model grid boxes.
   get_ocean_frac:
      Returns fractional area covered by ocean in a grid box.
   get_ocean_mask:
      Returns a land-ocean mask in a grid box.
   get_water_frac:
      Returns fractional area covered by water.
   get_water_mask:
      Returns a land-water mask in a grid box.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Get_topog_mean
      :name: get_topog_mean

   ::

      flag = <B> get_topog_mean </B> ( blon, blat, zmean )

   **DESCRIPTION**
      Returns realistic mountains on a latitude-longtude grid. The returned field is the mean topography for the given
      grid boxes. Computed using a conserving area-weighted interpolation. The current input data set is the 1/6 degree
      Navy mean topography.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blon``                                                  | The longitude (in radians) at grid box boundaries.        |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blat``                                                  | The latitude (in radians) at grid box boundaries.         |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``zmean``                                                 | The mean surface height (meters). The size of this field  |
      |                                                           | must be size(blon)-1 by size(blat)-1.                     |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``get_topog_mean``                                        | A logical value of TRUE is returned if the surface height |
      |                                                           | field was successfully created. A value of FALSE may be   |
      |                                                           | returned if the input topography data set was not         |
      |                                                           | readable.                                                 |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Get_topog_stdev
      :name: get_topog_stdev

   ::

      flag = <B> get_topog_stdev </B> ( blon, blat, stdev )

   **DESCRIPTION**
      Returns the standard deviation of the "finer" input topography data set, currently the Navy 1/6 degree mean
      topography data, within the boundaries of the given input grid.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blon``                                                  | The longitude (in radians) at grid box boundaries.        |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blat``                                                  | The latitude (in radians) at grid box boundaries.         |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``stdev``                                                 | The standard deviation of surface height (in meters)      |
      |                                                           | within given input model grid boxes. The size of this     |
      |                                                           | field must be size(blon)-1 by size(blat)-1.               |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``get_topog_stdev``                                       | A logical value of TRUE is returned if the output field   |
      |                                                           | was successfully created. A value of FALSE may be         |
      |                                                           | returned if the input topography data set was not         |
      |                                                           | readable.                                                 |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Get_ocean_frac
      :name: get_ocean_frac

   ::

      flag = <B> get_ocean_frac </B> ( blon, blat, ocean_frac )

   **DESCRIPTION**
      Returns fractional area covered by ocean in the given model grid boxes.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blon``                                                  | The longitude (in radians) at grid box boundaries.        |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blat``                                                  | The latitude (in radians) at grid box boundaries.         |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``ocean_frac``                                            | The fractional amount (0 to 1) of ocean in a grid box.    |
      |                                                           | The size of this field must be size(blon)-1 by            |
      |                                                           | size(blat)-1.                                             |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``get_ocean_frac``                                        | A logical value of TRUE is returned if the output field   |
      |                                                           | was successfully created. A value of FALSE may be         |
      |                                                           | returned if the Navy 1/6 degree percent water data set    |
      |                                                           | was not readable.                                         |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Get_ocean_mask
      :name: get_ocean_mask

   ::

      flag = <B> get_ocean_mask </B> ( blon, blat, ocean_mask )

   **DESCRIPTION**
      Returns a land-ocean mask in the given model grid boxes.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blon``                                                  | The longitude (in radians) at grid box boundaries.        |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blat``                                                  | The latitude (in radians) at grid box boundaries.         |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``ocean_frac``                                            | The fractional amount (0 to 1) of ocean in a grid box.    |
      |                                                           | The size of this field must be size(blon)-1 by            |
      |                                                           | size(blat)-1.                                             |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``get_ocean_mask``                                        | A logical value of TRUE is returned if the output field   |
      |                                                           | was successfully created. A value of FALSE may be         |
      |                                                           | returned if the Navy 1/6 degree percent water data set    |
      |                                                           | was not readable.                                         |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

e. .. rubric:: Get_water_frac
      :name: get_water_frac

   ::

      flag = <B> get_water_frac </B> ( blon, blat, water_frac )

   **DESCRIPTION**
      Returns the percent of water in a grid box.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blon``                                                  | The longitude (in radians) at grid box boundaries.        |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blat``                                                  | The latitude (in radians) at grid box boundaries.         |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``water_frac``                                            | The fractional amount (0 to 1) of water in a grid box.    |
      |                                                           | The size of this field must be size(blon)-1 by            |
      |                                                           | size(blat)-1.                                             |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``get_water_frac``                                        | A logical value of TRUE is returned if the output field   |
      |                                                           | was successfully created. A value of FALSE may be         |
      |                                                           | returned if the Navy 1/6 degree percent water data set    |
      |                                                           | was not readable.                                         |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

f. .. rubric:: Get_water_mask
      :name: get_water_mask

   ::

      flag = <B> get_water_mask </B> ( blon, blat, water_mask )

   **DESCRIPTION**
      Returns a land-water mask in the given model grid boxes.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blon``                                                  | The longitude (in radians) at grid box boundaries.        |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``blat``                                                  | The latitude (in radians) at grid box boundaries.         |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``water_mask``                                            | A binary mask for water (true) or land (false). The size  |
      |                                                           | of this field must be size(blon)-1 by size(blat)-1.       |
      |                                                           | [real, dimension(:,:)]                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``get_water_mask``                                        | A logical value of TRUE is returned if the output field   |
      |                                                           | was successfully created. A value of FALSE may be         |
      |                                                           | returned if the Navy 1/6 degree percent water data set    |
      |                                                           | was not readable.                                         |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

Namelist
--------

.. container::

   **&topography_nml**
   ``topog_file``
   Name of topography file.
   [character, default: DATA/navy_topography.data]
   ``water_file``
   Name of percent water file.
   [character, default: DATA/navy_pctwater.data]

| 

Data sets
---------

.. container::

   This module uses the 1/6 degree U.S. Navy mean topography and percent water data sets.
   These data sets have been re-formatted to separate 32-bit IEEE files. The names of these files is specified by the
   namelist input.
   The format for both files is as follows:
   ::

           record = 1    nlon, nlat
           record = 2    blon, blat
           record = 3    data

   where:
   ::

           nlon, nlat = The number of longitude and latitude points
                        in the horizontal grid.  For the 1/6 degree
                        data sets this is 2160 x 1080. [integer]
           blon, blat = The longitude and latitude grid box boundaries in degrees.
                           [real :: blon(nlon+1), blat(nlat+1)]

           data       = The topography or percent water data.
                          [real :: data(nlon,nlat)]

Error messages
--------------

.. container::

   **FATAL in get_topog_mean**
      shape(zmean) is not equal to (/size(blon)-1,size(blat)-1/))
      Check the input grid size and output field size.
   **FATAL in get_water_frac**
      shape(water_frac) is not equal to (/size(blon)-1,size(blat)-1/))
      Check the input grid size and output field size.

References
----------

.. container::

   None.

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

   To run this program you will need the topography and percent water data sets and use the following namelist (in file
   input.nml).
   &gaussian_topog_nml height = 5000., 3000., 3000., 3000., olon = 90., 255., 285., 0., olat = 45., 45., -15., -90.,
   wlon = 15., 10., 5., 180., wlat = 15., 25., 25., 20., /
   program test
   test program for topography and gaussian_topog modules
   ::

        use topography_mod
        implicit none
        
        integer, parameter :: nlon=24, nlat=18
        real :: x(nlon), y(nlat), xb(nlon+1), yb(nlat+1), z(nlon,nlat)
        real :: hpi, rtd
        integer :: i,j
        logical :: a
        
       gaussian mountain parameters
        real, parameter :: ht=4000.
        real, parameter :: x0=90., y0=45. ! origin in degrees
        real, parameter :: xw=15., yw=15. ! half-width in degees
        real, parameter :: xr=30., yr= 0. ! ridge-width in degrees
        
       create lat/lon grid in radians
          hpi = acos(0.0)
          rtd = 90./hpi ! rad to deg
          do i=1,nlon
            xb(i) = 4.*hpi*real(i-1)/real(nlon)
          enddo
            xb(nlon+1) = xb(1)+4.*hpi
            yb(1) = -hpi
          do j=2,nlat
            yb(j) = yb(j-1) + 2.*hpi/real(nlat)
          enddo
            yb(nlat+1) = hpi
       mid-point of grid boxes
          x(1:nlon) = 0.5*(xb(1:nlon)+xb(2:nlon+1))
          y(1:nlat) = 0.5*(yb(1:nlat)+yb(2:nlat+1))
       test topography_mod routines
          a = get_topog_mean(xb,yb,z)
          call printz ('get_topog_mean')
        
          a = get_water_frac(xb,yb,z)
          z = z*100. ! in percent
          call printz ('get_water_frac')
        
          a = get_ocean_frac(xb,yb,z)
          z = z*100. ! in percent
          call printz ('get_ocean_frac')
        
       test gaussian_topog_mod routines
          a = .true.
          z = get_gaussian_topog(x,y,ht,x0,y0,xw,yw,xr,yr)
          call printz ('get_gaussian_topog')
        
          call gaussian_topog_init (x,y,z)
          call printz ('gaussian_topog_init')
        
        contains
        
       simple printout of topog/water array
          subroutine printz (lab)
          character(len=*), intent(in) :: lab
           if (a) then
              print '(/a)', trim(lab)
           else
              print '(/a)', 'no data available: '//trim(lab)
              return
           endif
       print full grid
              print '(3x,25i5)', (nint(x(i)*rtd),i=1,nlon)
            do j=nlat,1,-1
              print '(i3,25i5)', nint(y(j)*rtd), (nint(z(i,j)),i=1,nlon)
            enddo
          end subroutine printz
        
        end program test

| 

Notes
-----

.. container::

   None.

| 
