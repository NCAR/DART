Aether Cube Sphere 
==================

This document describes the DART interface to the Aether ionosphere-thermosphere model in its cube
sphere implementation.

In addition to the standard DART programs associated with a given model, this interface creates
three additional programs for the Aether cube sphere: ``aether_to_dart``, ``dart_to_aether``, and 
``create_geometry_file``.

Use `this link <https://www.image.ucar.edu/pub/DART/aether_cube_sphere_restart_files.zip>`_ to 
download a set of Aether cube sphere restart files.

aether_to_dart
--------------

This program takes six Aether restart files and combines them into a single filter input file for
DART. This program reads entries from two namelists in ``input.nml``:

In ``&directory_nml``:

- ``restart_directory`` specifies the path where the restart files are saved
- ``filter_directory`` specifies the path where the ``filter_input`` file will be created (this
  entry can be the same as ``restart_directory`` but an additional entry is included for
  flexibility)

In ``&transform_state_nml``:

- ``restart_file_prefix``, ``restart_file_middle``, and ``restart_file_suffix`` specify substrings
  of the restart filename that are combined to create the full
  filename. The ensemble member string is inserted in between ``restart_file_prefix`` and
  ``restart_file_middle``. The grid face is inserted in between ``restart_file_middle`` and 
  ``restart_file suffix``.
- ``filter_input_prefix`` and ``filter_input_suffix`` specify substrings of the filter input
  filename. The ensemble member string is inserted in between ``filter_input_prefix`` and
  ``filter_input_suffix``.

Using aether_to_dart
~~~~~~~~~~~~~~~~~~~~

When executing ``aether_to_dart`` the ensemble member must be specified as a command line argument:

.. code-block::

    ./aether_to_dart 0000

The program must be run once for each ensemble member, ``./aether_to_dart 0000``,
``./aether_to_dart 0001`` and so on. 

dart_to_aether
--------------

This program takes a single filter output file and inserts them back into the six initial Aether 
restart files. This program uses all of the same namelist entries as ``aether_to_dart`` except for
instead of using ``filter_input_prefix`` and ``filter_input_suffix`` in ``&transform_state_nml`` it 
uses:

- ``filter_output_prefix`` and ``filter_output_suffix`` specify substrings of the filter output
  filename. The ensemble member string is inserted in between ``filter_output_prefix`` and
  ``filter_output_suffix``.

Using dart_to_aether
~~~~~~~~~~~~~~~~~~~~

When executing ``aether_to_dart`` the ensemble member must be specified as a command line argument:

.. code-block::

    ./dart_to_aether 0000

The program must be run once for each ensemble member, ``./dart_to_aether 0000``,
``./dart_to_aether 0001`` and so on. 

create_geometry_file
--------------------

This program takes two sets of six Aether grid files:

- ``grid_corners_g000?.nc``, where ``?`` corresponds to the six grid corner files, numbered 0-5, and
- ``grid_g000?.nc``, where ``?`` corresponds to the six grid center files, numbered 0-5.

Use `this link <https://www.image.ucar.edu/pub/DART/aether_cube_sphere_grid_files.zip>`_ to 
download a set of Aether cube sphere grid files.

Aether's model fields are defined at the grid centers. The grid corners are at locations 
approximately equidistant from either three (if the grid corner is at at a cube vertex) or four 
(if the grid corner is on a cube face) grid centers.

The program iterates through the grid corners and finds the nearest three or four grid centers in 
order to create a single ``geometry_file.nc`` that saves the cube sphere grid topology so that the
center fields can be used for interpolated by DART.

The two key compontents of the ``geometry_file.nc`` are:

1. The locations of the eight grid corners that coincide with the cube vertices and the column
   indices of the nearest three grid centers to each of the cube vertices. This structure defines a
   triangle surrounding each cube vertex.
2. The locations of the remaining ``N`` grid corners that lie on the cube faces and the column
   indices of the nearest four grid centers to each of the cube face grid corners. This structure 
   defines a quadrilateral surrouding each of the grid corners that lie on a cube face.

.. note::
   
   The ``in_quad`` function in DART's ``quad_utilities_mod`` does not implement a method to
   determine whether a point is in a quad if the given quad encloses either the North or South Pole.
   Thus ``create_geometry_file`` also determines which two grid corners coincide with the North and
   South Poles and stores the index of these corners as global attributes in ``geometry_file.nc``
   that can be read and used by the model mod to account for the special case where a quad encloses
   a pole.

The header of the ``geometry_file.nc`` created by ``create_geometry_file`` should resemble the
following:

.. code-block::

   netcdf geometry_file {
   dimensions:
      center_altitudes = 44 ;
      vertex_columns = 8 ;
      vertex_neighbors = 3 ;
      quad_columns = 1938 ;
      quad_neighbors = 4 ;
      center_columns = 1944 ;
   variables:
	  float center_altitude(center_altitudes) ;
	  float vertex_longitude(vertex_columns) ;
	  float vertex_latitude(vertex_columns) ;
	  float quad_longitude(quad_columns) ;
	  float quad_latitude(quad_columns) ;
	  float center_longitude(center_columns) ;
	  float center_latitude(center_columns) ;
	  int vertex_neighbor_indices(vertex_neighbors, vertex_columns) ;
	  int quad_neighbor_indices(quad_neighbors, quad_columns) ;

   // global attributes:
		:index_of_north_pole_quad_column = 1760 ;
		:index_of_south_pole_quad_column = 1403 ;
   }

geometry_file.nc dimensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``vertex_columns`` and ``quad_columns`` correspond to all of the grid corners on the sphere. 
Each of the ``vertex_columns`` are enclosed by three grid centers, which are referred to as
``vertex_neighbors``. Each of the ``quad_columns`` are enclosed by four grid centers, which are
refered to as ``quad_neighbors``.

The dimensions of the grid centers are ``center_altitudes`` in the vertical and ``center_columns``
in the horizontal.

geometry_file.nc variables
~~~~~~~~~~~~~~~~~~~~~~~~~~

The longitudes and latitudes for each of grid corners corresponding to the eight cube vertices are
stored in the ``vertex_longitude`` and ``vertex_latitude`` fields, respectively.

The longitudes and latitudes for each of the grid corners on the cube faces enclosed by
quadrilaterals are stored in the ``quad_longitude`` and ``quad_latitude`` fields, respectively.

The altitudes, longitudes and latitudes for each of the grid centers are stored in the
``center_altitude``, ``center_longitude`` and ``center_latitude`` fields, respectively.

.. important::

   The key feature of the ``geometry_file.nc`` is the relationship between the grid corners and the
   grid centers. This relationship is stored in the ``vertex_neighbor_indices`` and
   ``quad_neighbor_indices`` integer fields.

``vertex_neighbor_indices`` is a 3x8 integer array where each of the 8 columns corresponds to the 
a grid corner coinciding with a cube vertex and each of the three rows corresponds to the indices of
the center columns that define a triangle enclosing the cube vertex. 
``quad_neighbor_indices`` is a 4xN integer array where each of the N columns corresponds to 
a grid corner on a cube face and each of the four rows corresponds to the indices of the center
columns that define a quad that encloses the grid corner.
