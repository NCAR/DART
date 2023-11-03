===========
Even Sphere
===========

It is frequently useful to generate a series of synthetic observations 
located at roughly evenly-distributed locations [1]_ on a sphere.  

There are three methods described here.  

1.  A Matlab script and the standard DART observation generation utilities.  
2.  A csh script with all the parts of 1.  (Not available for all models).
3.  A stand-alone `Fortran`_ program. 

The Fortran program does not generate the nice plots that the Matlab process does, 
but it may be faster and easier to automate for generating a large number of obs.


Matlab Scripts Plus Standard DART Observation Executables
---------------------------------------------------------

This involves multiple steps:

1. determine how many locations are needed
2. determine the vertical levels needed
3. run the MATLAB function `even_sphere.m`_ to create the text file containing the input 
   for :code:`create_obs_sequence`
4. run :doc:`../../../assimilation_code/programs/create_obs_sequence/create_obs_sequence`
   to create an observation sequence (usually :code:`set_def.out`, although it is possible to 
   create :code:`obs_seq.out` files directly if you don't really care about the observation values).
5. if desired, run :doc:`../../../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq`
   to create a longer observation sequence file.
6. run :doc:`../../../assimilation_code/programs/perfect_model_obs/perfect_model_obs` 
   to harvest the synthetic observations from a chosen model.  

This directory contains a MATLAB function (:code:`even_sphere.m`) 
that generates input for the 
:doc:`../../../assimilation_code/programs/create_obs_sequence/create_obs_sequence` .  
After executing :code:`create_obs_sequence`, the resulting observation sequence file
will have a template for 'RADIOSONDE_TEMPERATURE','RADIOSONDE_U_WIND_COMPONENT',and 
'RADIOSONDE_V_WIND_COMPONENT' observations at specified pressure levels and roughly 
evenly-distributed locations across the entire globe. Optionally, bogus observation 
values may also inserted; which may be useful in certain circumstances.

even_sphere.m
-------------

:code:`even_sphere.m` has many optional arguments to tailor its behavior.
It has exactly 1 required argument - the number of horizontal locations desired.

- It will create a text file :code:`even_create_input` to be used as input to :code:`create_obs_sequence`
- The choice of pressure levels is described `here <Levels_>`_.
- The default number of pressure levels is 21.  Argument :code:`nlevels` specifies how many levels
  to use from the beginning of the levels list.  
- The default observation error variances for each observation type are level-dependent 
  and are consistent with 
  :code:`DART/observations/obs_converters/obs_error/ncep_obs_err_mod.f90`
- The default is to create 'empty' observation sequences - i.e. they have no actual 
  observation values and are suitable to be used with :code:`perfect_model_obs`
- The default date of the observations is 2017-12-25 00:00:00
- A plot of the locations will be created. The number of gridlines is configurable but 
  defaults to 288 in longitude and 192 in latitude.
- **All** the defaults can be changed by specifying 'variable-value' pairs of options, 
  as described below. Examples of some options are also available via the normal 
  MATLAB *help* facility. (Documenting all of them in the *help* makes the help page too long.)

Note that the number of observations will be the number of locations \* 
the number of vertical levels \* the number of variables (i.e. 3) 
:code:`even_sphere.m` also takes observation error variances 
and includes them in the observation sequences.

Optional Argument Variable-Value pairs
--------------------------------------

The optional variable-value pairs can appear in any order.

+-------------------+----------------------------+--------------------------------------------------+
| optional variable | example value              | Description                                      |
+===================+============================+==================================================+
| 'nlevels'         | 5                          | number of pressure levels to use.                |
|                   |                            | May be less than the length of the               |
|                   |                            | 'levels' array, but cannot be more.              |
+-------------------+----------------------------+--------------------------------------------------+
| levels            | [1000  500  300  200  100] | pressure levels desired.                         |
|                   |                            | see `Levels`_ section for discussion.            |
+-------------------+----------------------------+--------------------------------------------------+
| T_error_var       | [1.44 0.64 0.81 1.44 0.64] | level-specific                                   |
|                   |                            | Temperature error variances.                     |
|                   |                            | see `Levels`_ section for discussion.            |
+-------------------+----------------------------+--------------------------------------------------+
| W_error_var       | [1.96 4.41 9.00 7.29 4.41] | level-specific error variances                   |
|                   |                            | for both U, V wind components.                   |
|                   |                            | see `Levels`_ section for discussion.            |
+-------------------+----------------------------+--------------------------------------------------+
| 'YMD'             | '2017-12-25'               | Date required for :code:`create_obs_sequence`.   |
|                   |                            | If :code:`create_fixed_network_seq` is run, this |
|                   |                            | time is replaced.                                |
+-------------------+----------------------------+--------------------------------------------------+
| fill_obs          | false                      | 'true' inserts a bogus observation value of 1.0  |
|                   |                            | and a bogus QC value of 0.'false' does not insert|
|                   |                            | bogus values and essentially creates an empty    |
|                   |                            | obs sequence file (typically :code:`set_def.out`)|
+-------------------+----------------------------+--------------------------------------------------+
| 'nlon'            | 288                        | number of longitude grid lines in plot           |
+-------------------+----------------------------+--------------------------------------------------+
| 'nlat'            | 192                        | number of latitude grid lines in plot            |
+-------------------+----------------------------+--------------------------------------------------+

Examples
--------

1. 30 horizontal locations at 6 pressure levels:

.. code-block::

   nprofiles   = 30;
   levels      = [1000  850  500  300  200  100];
   T_error_var = [1.44 0.64 0.64 0.81 1.44 0.64];
   W_error_var = [1.96 2.25 4.41 9.00 7.29 4.41];
   even_sphere(nprofiles, 'levels', levels, ...
              'T_error_var', T_error_var, 'W_error_var', W_error_var)


2. 30 horizontal locations at 3 pressure levels. Note that the
   *nlevels* argument specifies that only the first 3 pressure levels
   are used even though there are 6 potential pressure levels. 
   Similarly, only the matching error variances are used.

.. code-block::

   nprofiles   = 30;
   nlevels     = 3 ;
   levels      = [1000  850  500  300  200  100];
   T_error_var = [1.44 0.64 0.64 0.81 1.44 0.64];
   W_error_var = [1.96 2.25 4.41 9.00 7.29 4.41];
   even_sphere(nprofiles, 'nlevels', nlevels, 'levels', levels, ...
              'T_error_var', T_error_var, 'W_error_var', W_error_var)

Levels
------
 
.. attention::

   If you need realistic error variances attached to your observations,
   be careful to align your levels and variances.

The default levels that this program generates are the *mandatory pressure levels* defined in the
`AMS glossary <https://glossary.ametsoc.org/wiki/Mandatory_level>`_.
The corresponding error variances are from ncep_obs_err_mod.  
See :doc:`../obs_error/README`.
Levels at the top can be excluded by setting *nprofiles* < 21 (size(levels)).

.. code::

   levels      = [1000  925  850  700  500  400  300   250  200  150  100   70   50   30   20   10    7    5    3    2    1];
   T_error_var = [1.44 1.00 0.64 0.64 0.64 0.64 0.81  1.44 1.44 1.00 0.64 0.64 0.81 1.00 1.69 2.25 2.25 2.25 2.25 2.25 2.25];
   W_error_var = [1.96 2.25 2.25 2.56 4.41 6.76 9.00 10.24 7.29 5.76 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41 4.41];

Here's an example of replacing the AMS levels with a set that has more levels near the surface
and none above 150 hPa.  Note that the error variances should change to be consistent
with the levels.

.. code::

   levels      = [1000  950  900  850  800  750  700  650  600  550  500  400  300  200  150];
   T_error_var = [1.44 1.21 0.81 0.64 0.64 0.64 0.64 0.64 0.64 0.64 0.64 0.64 0.81 1.44 1.00];
   W_error_var = [1.96 2.25 2.25 2.25 2.56 2.56 2.56 3.24 3.61 4.00 4.41 6.76 9.00 7.29 5.76];

Running Matlab in Batch Mode
----------------------------
 
If you would prefer to run :code:`even_sphere.m` in batch mode (i.e. from within a shell script),
here is an example syntax that worked for me. The script ran in the same directory
as :code:`even_sphere.m`. There are many ways to construct the input, naturally - but you don't have
to explicitly edit :code:`even_sphere.m` this way. 

.. code::

    #!/bin/csh

    \rm -rf matlab_input.m

    cat >> matlab_input.m << EndOfInput

       nprofiles   = 30;
       levels      = [1000  850  500  300  200  100];
       T_error_var = [1.44 0.64 0.64 0.81 1.44 0.64];
       W_error_var = [1.96 2.25 4.41 9.00 7.29 4.41];
       even_sphere(nprofiles, 'levels', levels, ...
                  'T_error_var', T_error_var, 'W_error_var', W_error_var)
       fname = sprintf('even_sphere_%d_profiles',nprofiles);
       orient landscape
       print(fname,'-dpdf')

    EndOfInput

    matlab -nosplash -nodesktop -r "try; cd $PWD; matlab_input; catch; end; exit";


Automation Scripts
------------------

Here there are also scripts (:code:`run_fixed_network_\*.csh`) which use the
output from :code:`create_obs_sequence` and the 
:doc:`../../../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq` 
to generate a series of observation sequence files.

run_fixed_network_seq.csh
~~~~~~~~~~~~~~~~~~~~~~~~~

Calls :code:`create_fixed_network_seq` to create a separate file for each time period.
By default, it makes 2 files/day, 12 hours apart, single time per file.
The intervals and dates can be changed by editing the script.
It assumes that :code:`create_fixed_network` has any model-specific files it needs in this directory.
It requires a :code:`set_def.out` file (usually created by :code:`create_obs_sequence`).

run_fixed_network_daily.csh
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calls :code:`create_fixed_network_seq` to create a separate file for each time period.
By default it makes 1 file/day, single time (noon) per file.
The dates and time can be changed by editing the script.
It assumes that :code:`create_fixed_network` has any model-specific files it needs in this directory.
It requires a :code:`set_def.out` file (usually created by :code:`create_obs_sequence`).

The process, end to end:

MATLAB:

Set the number of profiles, the levels, etc. and run :code:`even_sphere.m` in
MATLAB. It creates the necessary text file :code:`even_create_input` for the next step.
It will also make a plot - which you can save.

DART:

Then you have a choice about building and running the :code:`create_obs_sequence`
and :code:`create_fixed_network_seq` programs:

A. building them in the :code:`models/template/work` directory 
B. using the ones which were built in :code:`models/your_model/work` directory 
   by :code:`quickbuild.sh`.

Choice A uses programs which have no model specific file dependencies,
but may involve more separate steps than B.

A
~~~~~~

1. Build the programs in :code:`template/work`
2. Link (or copy) these files to the directory 
   in which you want to create obs_seq files.

.. code-block:: text

   ./even_create_sequence 
   ./run_fixed_network_{seq or daily}.csh
   models/template/work/create_fixed_network_seq
   models/template/work/create_obs_sequence
   models/template/work/input.nml

3. In your obs_seq directory, run create_obs_sequence, 
   which creates a :code:`set_def.out` file.

.. code-block:: text

   ./create_obs_sequence < even_create_input > /dev/null

4. Edit and run your choice of :code:`run_fixed_network_\*.csh` for the desired dates.
   These call create_fixed_network_seq, which creates an :code:`obs_seq.in` file
   for each specified date.

B
~~~~~~

This choice may involve fewer steps, *if* there is a model specific script
which combines the steps in A).  
See the cam-fv example (models/cam-fv/shell_scripts/synth_obs_locs_to_seqs.csh).
If there is *not* a script like that for your model,
you can follow the steps in A), 
substituting your model name for the "template" in the pathnames. 
NOTE: you may need to link any additional input files which your model requires
into the directory where you will run the programs.
These typically contain grid information and are found in :code:`your_model/work`.
For example, *cam-fv* needs a :code:`caminput.nc` and :code:`cam_phis.nc`.

.. _Fortran:

Fortran program for generating obs directly
-------------------------------------------


cd into the work directory and run ``quickbuild.sh``.

This builds the ``create_even_sphere`` executable.  Edit the ``input.nml``
to set the number of obs to generate and the date in the namelist.  Run
the program and the output file will be generated.



DETAILS of generating points evenly-distributed on a sphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the algorithm that's being used [1]_:

.. code-block:: text

  N     := the number of profiles you want
  dlong := pi*(3-sqrt(5))  /* ~2.39996323 */
  dy    := 2.0/N
  phi   := 0
  y     := 1 - dy/2

  for k := 0 .. N-1
      r       := sqrt(1-y*y)
      node[k] := (cos(phi)*r, sin(phi)*r, y)
      y       := y - dy
      phi     := phi + dlong

For the geometric and visually minded: 

#. Picture a unit sphere in cartesion space (x,y,z).
#. Choose a value -1 < y < 1, which defines an x-z plane.
   That plane intersects with the unit sphere to form a circle
   whose center is on the y axis.  (The circle radius is small 
   near y = +/-1 and is 1 at y=0.)
#. Choose an angle ("phi") and draw a ray 
   from the center of the circle to a point on the circle using this angle 
   relative to the x positive direction.  Where the ray intersects the circle
   (and sphere) is one of the evenly distributed points on the sphere 
   which we want.  
#. Its x and z coordinates can then be combined
   with the already defined y coordinate to define the cartesian location 
   of the point.
#. The choice of the y and angle for each point is where the magic enters the algorithm.
   They are derived from the Fibonacci or Golden Spiral formula (derived elsewhere).


.. [1] A python example of the Golden Section spiral algorithm can be found in
    https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    See the contribution from Fab von Bellinghousen.
