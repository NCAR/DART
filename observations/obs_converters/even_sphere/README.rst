===========
Even Sphere
===========

Generate a series of synthetic observations located at roughly
evenly distributed locations on a sphere.  At each location
generate a vertical column of observations.  This could mimic
a radiosonde observing network, for example.


This directory contains MATLAB scripts (even_sphere*.m) 
that generate input for the `create_obs_sequence 
<../../../assimilation_code/programs/create_obs_sequence/create_obs_sequence.html>`_ program.  
They takes a number of vertical levels and a total number of points,
and generate a roughly evenly distributed set of observations
across the entire globe.  Note that the number of obs
will be the number of points times the number of vertical levels.  
``Even_sphere.m`` also takes observation error variances 
and includes them in the observation templates.

Here there are also scripts (run_fixed_network_*.csh) which use the
output from create_obs_sequence and the `create_fixed_network_seq 
<../../../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq.html>`_ 
program to generate a series of obs_seq.in files.


The process, end to end:

MATLAB:

Edit even_sphere.m and set the number of levels, the
number of profiles, the vertical coordinate type, etc.     

Run it in MATLAB.  It will make a plot (which you can 
save from the menu) and it will create a text file 'even_create_input'.

DART:

Then you have a choice about building and running the ``create_obs_sequence``
and ``create_fixed_network_seq`` programs:

A. building them in the models/template/work directory 
B. using the ones which were built in models/your_model/work directory 
   by quick_build.csh. 

Choice A uses programs which have no model specific file dependencies,
but may involve more separate steps than B:

1. Build the programs in template/work
2. Link (or copy) these files to the directory 
   in which you want to create obs_seq files.

.. code-block:: text

   ./even_create_sequence 
   ./run_fixed_network_{seq or daily}.csh
   models/template/work/create_fixed_network_seq
   models/template/work/create_obs_sequence
   models/template/work/input.nml

3. In your obs_seq directory, run create_obs_sequence, 
   which creates a ``set_def.out`` file.

.. code-block:: text

   ./create_obs_sequence < even_create_input > /dev/null

4. Edit and run your choice of ``run_fixed_network_*.csh`` for the desired dates.
   These call create_fixed_network_seq, which creates an ``obs_seq.in`` file
   for each specified date.

\(B)
This choice may involve fewer steps, *if*\ there is a model specific script
which combines the steps in A).  
See `the cam-fv example <models/cam-fv/shell_scripts/synth_obs_locs_to_seqs.csh>`_.
If there is *not*\ a script like that for your model,
you can follow the steps in A), 
substituting your model name for the "template" in the pathnames. 
NOTE: you may need to link additional input files, which your model needs to start, 
into the directory where you will run the programs.
These typically contain grid information and are found in your_model/work.
For example, cam-fv needs a caminput.nc and cam_phis.nc.


DETAILS on generating points evenly distributed on a sphere

This is the algorithm that's being used:

.. code-block:: text

  N     := the number of profiles you want
  dlong := pi*(3-sqrt(5))  /* ~2.39996323 */
  dy    := 2.0/N
  long  := 0
  y     := 1 - dy/2

  for k := 0 .. N-1
      r       := sqrt(1-y*y)
      node[k] := (cos(long)*r, sin(long)*r, y)
      y       := y - dy
      long    := long + dlong

For the geometric and visually minded: 

#. Picture a unit sphere in cartesion space (x,y,z).
#. Choose a value -1 < y < 1, which defines an x-z plane.
   That plane intersects with the unit sphere to form a circle
   whose center is on the y axis.  (The circle radius is small 
   near y = +/-1 and is 1 at y=0.)
#. Choose an angle ("long" above, "phi" in the script) and draw a ray 
   from the center of the circle to a point on the circle using this angle 
   relative to the x positive direction.  Where the ray intersects the circle
   (and sphere) is one of the evenly distributed points on the sphere 
   which we want.  
#. Its x and z coordinates can then be combined
   with the already defined y coordinate to define the cartesian location 
   of the point.
#. The choice of the y and angle for each point is where the magic enters the algorithm.
   They are derived from the Fibonacci or Golden Spiral formula (derived elsewhere).

