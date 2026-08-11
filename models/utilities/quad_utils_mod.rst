.. _quad_utils_mod:

======================
quad_utils_mod
======================

.. contents:: 
   :depth: 3
   :local:

Overview
--------

The ``quad_utils_mod`` module provides horizontal interpolation utilities for
longitude–latitude grids that are logically rectangular. It supports regular,
irregularly spaced, and fully curvilinear grids. Various forward operators use this 
interpolation routine to map between observation locations and model state variables.

The module implements both grid search and interpolation capabilities.

Supported Grid Types
--------------------

The module supports three grid types:

``GRID_QUAD_FULLY_REGULAR``
   A regular latitude–longitude grid with uniform spacing. Grid coordinates are
   defined analytically using starting values and constant increments. Example: 

.. code-block:: text

   lat ↑
        (0,2) --- (1,2) --- (2,2)
          |        |        |
        (0,1) --- (1,1) --- (2,1)
          |        |        |
        (0,0) --- (1,0) --- (2,0)
                        → lon

``GRID_QUAD_IRREG_SPACED_REGULAR``
   A logically rectangular grid with non-uniform spacing in longitude and/or
   latitude. Coordinates are defined using 1D arrays. Example: 

.. code-block:: text

   lat ↑
        (0,2) ---- (1.5,2) ------ (3,2)
          |          |            |
        (0,1) ---- (1.5,1) ------ (3,1)
          |          |            |
        (0,0) ---- (1.5,0) ------ (3,0)
                             → lon

``GRID_QUAD_FULLY_IRREGULAR``
   A curvilinear grid defined by 2D longitude and latitude arrays. Grid cells are
   general quadrilaterals and may be distorted. This type is used for dipole and 
   tripole grids. Example: 

.. code-block:: text

   lat ↑
            o--------o-------o
           /        /       /
          /        /       /
         o--------o-------o
         |      / |     /
         |    /   |   /
         o------o---o
                          → lon

All grid types are assumed to be logically rectangular.

General Workflow
----------------

Before using the interpolation utilities, the user must first identify the
structure of the model grid. In particular, determine whether the grid is:

- fully regular (uniform spacing),
- logically rectangular with non-uniform spacing, or
- fully irregular (curvilinear, defined by 2D coordinate arrays).

This classification determines how the interpolation handle is initialized and
how grid coordinates are provided.

The typical usage pattern is:

1. Initialize an interpolation handle using ``init_quad_interp`` with the
   appropriate grid type.
2. Provide grid coordinate information using one of the ``set_*`` routines
   corresponding to the grid type.
3. Locate the grid cell containing a given longitude/latitude using
   ``quad_lon_lat_locate``.
4. Interpolate values using ``quad_lon_lat_evaluate``.
5. Finalize and deallocate resources using ``finalize_quad_interp``.

Interpolation Methods
---------------------

Regular and Irregularly Spaced Grids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For regular and semi-regular grids, interpolation is performed using standard
bilinear interpolation in index space. The method:

- identifies the bounding grid indices,
- computes fractional distances in longitude and latitude, and 
- performs sequential linear interpolation.

Fully Irregular Grids
^^^^^^^^^^^^^^^^^^^^^

For fully irregular grids, interpolation is performed in physical space using
the four corners of the enclosing quadrilateral.

The interpolant is assumed a surface, defined as:

.. math::

   f(x, y) = a + b x + c y + d x y

The coefficients are determined from the values at the four corners of the
quadrilateral by solving a local 3x3 linear system. This system is solved 
using Gaussian elimination with partial pivoting.

An optional rotation of the quadrilateral can be applied prior to interpolation
to improve numerical behavior.

Search Algorithm
""""""""""""""""

For fully irregular grids, a two-stage search is used to efficiently locate the
containing quadrilateral:

1. A coarse regular latitude–longitude grid is constructed.
2. Each coarse grid cell stores a list of overlapping quadrilaterals.
3. At runtime, only the subset of candidate quadrilaterals associated with the
   relevant coarse cell is searched.

This significantly reduces the computational cost compared to a global search.

Interpolation Options
---------------------

The interpolation handle supports several configuration options:

``global_grid``
   Indicates that the grid is global. Enables cyclic longitude and pole wrapping.

``spans_lon_zero``
   Indicates that the grid crosses the 0/360 longitude discontinuity.

``pole_wrap``
   Enables wrapping across the poles.

``north_to_south`` 
   Indicates latitude ordering.

``grid_staggering``
   Specifies the location of grid points relative to the cell (center, edges,
   or corners).

These options influence both search and interpolation behavior.

Special Considerations
----------------------

Dipole and Tripole Grids
^^^^^^^^^^^^^^^^^^^^^^^^

Dipole and tripole grids are treated as fully irregular grids. No special
interpolation algorithm is used for these grids; instead, they are handled
through the general quadrilateral interpolation framework.

Mask Handling
^^^^^^^^^^^^^

An optional mask can be provided for fully irregular grids. If any of the four
corners of a quadrilateral cell are masked, interpolation within that cell is
not performed.

Pole Handling
^^^^^^^^^^^^^

The module includes limited support for grids that wrap over the poles. In some
cases, particularly for staggered grids, interpolation may fail near the pole.
These cases are identified and flagged during the location step. 

Limitations
-----------

- Only logically rectangular grids are supported.
- Fully unstructured meshes are not supported.
- Fully irregular grids require storage of 2D longitude and latitude arrays,
  which may increase memory usage for very high-resolution grids.
- Mask handling is conservative; partial-cell interpolation is not supported.
- Interpolation in highly distorted or nearly degenerate quadrilaterals may
  fail due to ill-conditioning of the local system.
- Multiple interpolation schemes are not currently implemented, although the
  framework allows for future extension.
