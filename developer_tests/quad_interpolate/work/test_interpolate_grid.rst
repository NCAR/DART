
PROGRAM ``test_interpolate_grid``
========================

Overview
--------
 
| Test the quad interpolation routines by moving known data between different grids.

Usage
----

From ~nancy on cheyenne, copy caminput.nc into the work directory.

::

   test_interpolate_grid


Output
------

Creates several NetCDF files which can be viewed with ncview.


- field0.nc is the original data grid
- fieldT.nc is that data transferred to the cam grid
- field1.nc is that data transferred back to a copy of the original grid
- field2.nc is that data transferred back to a denser grid

To view the data, run:

::

 ncview field0.nc

Since the grid is lat/lon, ncview may try to add continent outlines.
To turn that off, use the 'opts' button and set Overlays to 'None'.  
Use the 'OK' button.  Closing any ncview windows with the red button
on the menu bar will crash the program.

To diff the data, run:

::

 ncdiff field0.nc field1.nc diff.nc
 ncview diff.nc


This is going from a regular grid to another regular grid and back.
A test case with an irregular MOM 6 grid is still to come.


