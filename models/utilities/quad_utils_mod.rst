MODULE quad_utils_mod
================

Overview
--------

! For any model which uses a rectangular grid or an irregular grid where each element of the grid is a quad (even if deformed), these routines compute which quad contains the location and the interpolated data value.


::

   use quad_util_mod


Namelist
--------

::

namelist /quad_interpolate_nml/ do_rotate, debug


Other modules used
------------------

::

   types_mod
   location_mod
   utilities_mod
   options_mod


Public interfaces
-----------------

======================= ===================================
*use quad_utils_mod, only :* quad_interp_handle,
\          init_quad_interp,                
\          finalize_quad_interp,            
\          set_quad_coords,                 
\          quad_lon_lat_locate,             
\          quad_lon_lat_evaluate,           
\          GRID_QUAD_FULLY_REGULAR,         
\          GRID_QUAD_IRREG_SPACED_REGULAR,  
\          GRID_QUAD_FULLY_IRREGULAR,       
\          GRID_QUAD_UNKNOWN_TYPE,          
\          QUAD_LOCATED_UNKNOWN,            
\          QUAD_LOCATED_CELL_CENTERS,       
\          QUAD_LOCATED_LON_EDGES,          
\          QUAD_LOCATED_LAT_EDGES,          
\          QUAD_LOCATED_CELL_CORNERS,       
\          get_quad_grid_size,              
\          get_quad_global,                 
\          print_quad_handle                
======================= ===================================

Description
-----------

Interpolation routines for longitude/latitude grids which are logically 
rectangular and either fully regular, partially regular or fully deformed.

This module includes initialization, search, interpolation, and finalization
routines.

The size of the grid is specified by a count of the longitudes and latitudes.

The actual coordinates of the grid can be specified in one of 3 ways:
           
   fully regular grid, evenly spaced and fully orthogonal:
    (origin, delta) each for lon and lat.

   fully orthogonal but possibly irregularly spaced along the axes:
    1D array(counts) each for lon and lat.

   logically rectangular but corners of quads are fully irregular:
    2D array(lon counts, lat counts) each for lon and lat.

To search for a given location and return the (i,j) indices of the
corners of the enclosing quad:

   for the fully regular grid the enclosing quad is found by computing
   the corresponding (i,j) indices.
      
   for the irregularly spaced grid, the enclosing quad is found by searching
   each of the 1D arrays for the enclosing (i,j) indices.
      
   For the fully irregular grid, the enclosing quad is found by this process:
    At setup time: 
    1) create a coarse fully regular grid.
    2) compute the intersection of each regular grid box with the target grid quads
    3) keep a list of which target grid quads overlap for each regular grid box.
    At search time:
    4) find the location in the fully regular grid (which can be done quickly)
    5) do an exhaustive search of the target grid quads which overlap that
       regular grid box, returning when one of the target grid quads encloses
       the given location. 

Examples of expected usage:

Interpolation of a single field value at a given location:

:: 

  call init_quad_interp()
  call set_quad_coords()
  do i=1, num_locations_to_interpolate
      call quad_lon_lat_locate()
      call quad_lon_lat_evaluate()
  enddo
  call finalize_quad_interp()
   
Interpolation of multiple field values at the same location:

:: 

  call init_quad_interp()
  call set_quad_coords()
  do i=1, num_locations_to_interpolate
      call quad_lon_lat_locate()
      do j=1, nfields_at_this_loc
          call quad_lon_lat_evaluate()
      enddo
   enddo
   call finalize_quad_interp()



Files
-----

none

References
----------

#. none

Private components
------------------

N/A
