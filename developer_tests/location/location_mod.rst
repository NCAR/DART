MODULE location_mod
===================

Overview
--------

DART provides a selection of options for the coordinate system in which all observations and all model state vector
locations are described. All executables are built with a single choice from the available location modules. The names
of these modules are all ``location_mod``.

Introduction
------------

The core algorithms of DART work with many different models which have a variety of coordinate systems. This directory
provides code for creating, setting/getting, copying location information (coordinates) independently of the actual
specific coordinate information. It also contains distance routines needed by the DART algorithms.

Each of the different location_mod.f90 files provides the same set of interfaces and defines a 'module location_mod', so
by selecting the proper version in your path_names_xxx file you can compile your model code with the main DART routines.

-  :doc:`../../assimilation_code/location/threed_sphere/location_mod`:
   The most frequently used version for real-world 3d models. It uses latitude and longitude for horizontal coordinates,
   plus a vertical coordinate which can be meters, pressure, model level, surface, or no specific vertical location.
-  :doc:`../../assimilation_code/location/oned/location_mod`:
   The most frequently used for small models (e.g. the Lorenz family). It has a cyclic domain from 0 to 1.
-  :doc:`../../assimilation_code/location/threed_cartesian/location_mod`:
    A full 3D X,Y,Z coordinate system.
-  :doc:`../../assimilation_code/location/channel/location_mod`:
   a 3d domain periodic in x, limited in y, and unlimited z.
-  column: no x,y but 1d height, pressure, or model level for vertical.
-  annulus: a hollow 3d cylinder with azimuth, radius, and depth.
-  twod: a periodic 2d domain with x,y coordinates between 0 and 1.
-  twod_sphere: a 2d shell with latitude, longitude pairs.
-  threed: a periodic 3d domain with x,y,z coordinates between 0 and 1.

Other schemes can be added, as needed by the models. Possible ideas are a non-periodic version of the 1d, 2d cartesian
versions. Email `dart at ucar.edu <mailto:dart@ucar.edu>`__ if you have a different coordinate scheme which we might
want to support.

Namelist
--------

Each location module option has a different namelist. See the specific documentation for the location option of choice.

Files
-----

-  none

References
----------

-  none

Private components
------------------

N/A
