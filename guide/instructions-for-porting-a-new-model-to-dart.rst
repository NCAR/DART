.. _Porting new models:

Instructions for porting a new model to DART
============================================
To determine if your model is compatible with DART, see :ref:`Using new models`.

DART provides a script ``new_model.sh`` which will create the necessary files
for a new model interface.

Templates are chosen based on location module input. The currently supported
location templates are for 3D and 1D modules, with the possibility for more
in the future. At the moment, ``threed_sphere``, ``threed_cartesian``, and
``oned``  will produce model_mod.f90 code that compile will sucessfully with ``./quickbuild.sh``.

Enter ``./new_model.sh``, then the desired model name and location module separated
by spaces. This will create the necessary files to get started.

For example to create a model interface for a model called BOUMME which uses
the 3D sphere location module:

.. code-block::  text

   cd models
   ./new_model.sh BOUMME threed_sphere

This will create an BOUMME model directory with the following files:

.. code-block:: text

     BOUMME/
     ├── model_mod.f90
     ├── readme.rst
     └── work
         ├── input.nml
         └── quickbuild.sh

- ``model_mod.f90`` is where to add the :ref:`Required model_mod routines`.
- ``readme.rst``  is a stub to add documenation for your model interface.  
- ``quickbuild.sh`` is used to compile DART for your model.

Navigate to the work directory and enter ``./quickbuild.sh`` and everything should compile at this point.
Please note that you will need to run ``./quickbuild.sh`` again after making edits to ``model_mod.f90`` to recompile. 

The DAReS team recommends that you look at the existing supported models and reusing code from them if
possible when you write the code required for DART. Models with similar grid types 
or vertical coordinates are good candidates.

There is often a sensible default implementation that can be used for each of these routines as well. For
more information on what the default behavior for each routine is and how to use the default implementations,
see :ref:`Required model_mod routines`.

The required subroutines are these:

.. code-block:: text

   public :: get_model_size, &
             get_state_meta_data,  &
             model_interpolate, &
             shortest_time_between_assimilations, &
             static_init_model, &
             init_conditions,    &
             adv_1step, &
             nc_write_model_atts, &
             pert_model_copies, &
             nc_write_model_vars, &
             init_time, &
             get_close_obs, &
             get_close_state, &
             end_model, &
             convert_vertical_obs, &
             convert_vertical_state, &
             read_model_time, &
             write_model_time


If needed, model_mod can contain additional subroutines that are used
for any model-specific utility programs.  No routines other than
these will be called by programs in the DART distribution.

Edit the model_mod and fill in these routines:

#. ``static_init_model()`` - make it read in any grid information
   and the number of variables that will be in the state vector.
   Fill in the model_size variable.    Now ``get_model_size()`` and 
   ``get_model_time_step()`` from the template should be ok as-is.

#. ``get_state_meta_data()`` - given an index number into the state vector 
   return the location and kind.

#. ``model_interpolate()`` - given a location (lon/lat/vert in 3d, x in 1d)
   and a state QTY_xxx kind, return the interpolated value the field
   has at that location.   this is probably one of the routines that
   will take the most code to write.

For now, ignore these routines:

.. code-block:: text

   nc_write_model_vars()
   get_close_obs()
   get_close_state()
   end_model()
   convert_vertical_obs()
   convert_vertical_state()
   read_model_time()
   write_model_time()

If you have data in a initial condition/restart file, then you
can ignore these routines:

.. code-block:: text

   shortest_time_between_assimilations()
   init_conditions()

Otherwise, have them return an initial time and an initial default
ensemble state.

If your model is NOT subroutine callable, you can ignore this routine:

.. code-block:: text

   adv_1step()

Otherwise have it call the interface to your model and put the files
necessary to build your model into the models/YOUR_MODEL directory.

If you want to let filter add gaussian noise to a single state vector
to generate an ensemble, you can ignore this routine:

.. code-block:: text

   pert_model_copies()

Otherwise fill in code that does whatever perturbation makes sense
to have an initial ensemble of states.  in some cases that means
adding a different range of values to each different field in the
state vector.

At this point you should have enough code to start testing with
the ``model_mod_check`` program.  It is a stand-alone utility
that calls many of the model_mod interface routines and should
be easier to debug than some of the other DART programs.


Once you have that program working you should have enough code
to test and run simple experiments.


The general flow is:

#. ``./create_obs_sequence`` - make a file with a single observation in it

#. ``./perfect_model_obs`` - should interpolate a value for the obs

#. generate an ensemble of states, or set 'perturb_from_single_instance' to .true.

#. run ``./filter`` with the single observation 

#. Look at the preassim.nc and analysis.nc files
   Diff them with ``ncdiff``:

   .. code-block:: text

      ncdiff analysis.nc preassim.nc Innov.nc

   plot it, with ``ncview`` if possible:  

   .. code-block:: text

      ncview Innov.nc

   The difference between the two is the impact of that single observation
   see if it's at the right location and if the differences seem reasonable


If your model data cannot be output in NetCDF file format, or cannot
be directly converted to NetCDF file format with the ncgen program,
there are 2 additional steps:

* ``model_to_dart`` - read your native format and output data in NetCDF format

* ``dart_to_model`` - write the updated data back to the native file format


More details on each of these 5 steps follows. There is a more in-depth description of each individual program here: :ref:`DART programs`.

Running ``model_to_dart`` if needed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your model data is not stored in NetCDF file format, a program to
convert your data from the model to NetCDF is needed.  It needs to
read your model data in whatever format it uses and create NetCDF
variables with the field names, and appropriate dimensions if these
are multi-dimensional fields (e.g. 2d or 3d).  If the data is ASCII,
the generic NetCDF utility ncgen may be helpful.

Running ``create_obs_sequence``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can make a synthetic observation (or a series of them) with this
interactive program and use them for testing.  Before running make sure
the observation types you want to use are in the input.nml file in the
&obs_kind_nml section, either in the assimilate or evaluate lists.

Run the program.  Give the total number of obs you want to create
(start with 1).  Answer 0 to number of data items and 0 to number of
quality control items.  Answer 0 when it says enter -1 to quit.  You
will be prompted for an observation number to select what type of
observation you are going to test.  

Give it a location that should be inside your domain, someplace where
you can compute (by hand) what the correct value should be.  When it
asks for time, give it a time that is the same as the time on your
model data.

When it asks for error variance, at this point it doesn't matter.
give it something like 10% of the expected data value.  Later on
this is going to matter a lot, but for testing the interpolation of
a single synthetic obs, this will do.

For an output filename, it suggests 'set_def.out' but in this case
tell it 'obs_seq.in'.


Running ``perfect_model_obs``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure the NetCDF file with your input data matches the input name 
in the input.nml file, the &perfect_model_obs_nml namelist.  
Make sure the input obs_sequence is still set to 'obs_seq.in'.
run perfect_model_obs.  Something bad will happen, most likely.  Fix it.

Eventually it will run and you will get an 'obs_seq.out' file.  For these
tests, make sure &obs_sequence_nml : write_binary_obs_sequence = .false.
in the input.nml file.  The sequence files will be short and in ascii.
You can check to see what the interpolated value is.  if it's right, congratulations.
If not, debug the interpolation code in the model_mod.f90 file.


Using a single input state
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the &filter_nml namelist, set 'perturb_from_single_instance' to .true.
this tells filter that you have not generated N initial conditions,
that you are only going to supply one and it needs to perturb that
one to generate an initial ensemble.  Make sure the 'input_state_files' 
matches the name of the single state vector file you have.  You can
use the 'obs_seq.out' file from the perfect_model run because now
it has data for that observation.  Later on you will need to decide
on how to generate a real set of initial states, and then you will
set 'perturb_from_single_instance' back to .false. and 
supply N files instead of one.  You may need to set the 
&ensemble_manager_nml : perturbation_amplitude
down to something smaller than 0.2 for these tests - 0.00001 is a good
first guess for adding small perturbations to a state.


Running ``filter``
~~~~~~~~~~~~~~~~~~

Set the ens_size to something small for testing - between 4 and 10 is
usually a good range.  Make sure your observation type is in the
'assimilate_these_obs_types' list and not in the evaluate list.
run filter.  Find bugs and fix them until the output 'obs_seq.final' 
seems to have reasonable values.  Running filter will generate 
NetCDF diagnostic files.  The most useful for diagnosis will
be comparing preassim.nc and analysis.nc.


Diagnostics
~~~~~~~~~~~

Run 'ncdiff analysis.nc preassim.nc differences.nc' and use
your favorite netcdf plotting tool to see if there are any differences
between the 2 files.  For modules using a regular lat/lon grid 'ncview'
is a quick way to scan files.  For something on an irregular
grid a more complicated tool will have to be used.  If the files are
identical the assimilation didn't do anything.  Check to see if there
is a non-zero DART quality control value in the obs_seq.final file.
Check to see if there are errors in the dart_log.out file.  Figure out
why there's no change.  If there is a difference, it should be at
the location of the observation and extend out from it for a short
distance.  If it isn't in the right location, look at your get_state_meta_data()
code.  If it doesn't have a reasonable value, look at your model_interpolate() code.


Running ``dart_to_model`` if needed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After you have run filter, the files named in the 'output_state_files' namelist
item will contain the changed values.  If your model is reading NetCDF format
it can ingest these directly.  If not, an additional step is needed to copy
over the updated values for the next model run.

