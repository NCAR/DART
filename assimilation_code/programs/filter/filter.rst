PROGRAM ``filter``
==================

Overview
--------

Main program for driving ensemble filter assimilations.

``filter`` is a Fortran 90 program, and provides a large number of options for controlling execution behavior and
parameter configuration that are driven from its namelist. See the namelist section below for more details. The number
of assimilation steps to be done is controlled by the input observation sequence and by the time-stepping capabilities
of the model being used in the assimilation.

This overview includes these subsections:

-  Program Flow
-  Filter Types
-  Getting Started
-  Free Model Run after Assimilation
-  Evaluate a Model State against Observations
-  Compare Results with and without Assimilation
-  DART Quality Control Values on Output
-  Description of Inflation Options
-  Detailed Program Flow

See :ref:`Welcome page` for more documentation, including a discussion of the
capabilities of the assimilation system, a diagram of the entire execution cycle, the options and features.

Program flow
~~~~~~~~~~~~

The basic execution loop is:

-  Read in model initial conditions, observations, set up and initialize
-  Until out of observations:

   -  Run multiple copies of the model to get forecasts of model state
   -  Assimilate all observations in the current time window
   -  Repeat

-  Write out diagnostic files, restart files, final observation sequence file

The time of the observations in the input observation sequence file controls the length of execution of filter.

For large, parallel models, the execution loop is usually wrapped in an external script which does these additional
steps:

-  Link to an observation sequence file which contains only observation times within the next assimilation window
-  Link any output inflation files from the previous step to be the input files for this step
-  Run filter, which will exit after doing the assimilation without trying to advance the model
-  Save the output diagnostic files for later
-  Advance the N copies of the model using the model scripts or whatever method is appropriate
-  Repeat until all data is assimilated

For large models filter is almost always compiled to be a parallel MPI program, and most large models are themselves a
parallel program using OpenMP, MPI, or both. MPI programs usually cannot start other MPI programs, so the external
script submits both the filter job and the N model advances to a batch system so all run as independent parallel jobs.

The same source code is used for all applications of filter. The code specific to the types of observations and the
interface code for the computational model is configured at compile time. The top level directory has been simplified
from previous versions to look like :

-  ``README``
-  ``COPYRIGHT``
-  *assimilation_code*
-  *build_templates*
-  *diagnostics*
-  *documentation*
-  *models*
-  *observations*

the *assimilation_code* contains all *module* and *program* source code for all of the main programs including filter.
Specifically in the modules directory there is a ``filter_mod.f90`` which contains the source for the filter main
program. Each model has a separate directory under DART/models, and under each model is a work directory where the code
is compiled and can be run for testing. Generally when a full-size experiment is done the executables are copied to a
different location - e.g. scratch space on a large filesystem - since the data files for 10s to 100s of copies of a
model can get very large.

Directories expected to be modified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DART is distributed as a toolkit/library/facility that can be used as-is with the existing models and observations, but
is also designed so that users can add new models, new observation types and forward operators, and new assimilation
algorithms.

The locations in the DART `code tree <../../../docs/index.html#Directories>`__ which are intended to be modified by
users are:

New Models
   Add a new directory in the ``models`` subdirectory. Copy (recursively, e.g. ``cp -r``) the contents of the
   ``template`` directory and modify from there. Note that the ``model_mod.f90`` file in the template dir is appropriate
   for small models; for large geophysical models see the ``full_model_mod.f90`` file and also examine other model
   directories for ideas. See additional documentation in the :doc:`../../../models/template/readme` documentation,
   and the `DART web pages <http://www.image.ucar.edu/DAReS/DART/DART2_Documentation.php#adding_a_model>`__ on adding
   new models.
New Observation Platforms
   To convert observations from other formats to DART format, add a new directory in the ``observations/obs_converters``
   subdirectory and populate it with converter code.
New Observation Types and Forward Operators
   Define a new type (a measurement from an observing platform) via a file in the ``observations/forward_operators``
   subdirectory. If the forward operator is more complicated than directly interpolating a field in the model state,
   this is where the code for that goes. See additional documentation in the
   :doc:`../../../observations/forward_operators/obs_def_mod` documentation, and the `DART web
   pages <http://www.image.ucar.edu/DAReS/DART/DART2_Observations.php#adding_types>`__ on adding new types. Adding a new
   type may require adding a new ``generic kind``, which is documented in
   :doc:`../../modules/observations/obs_kind_mod`.
New Assimilation Algorithms
   If you want to try out a different filter type modify the filter code in the ``assim_tools_mod.f90`` file. See the
   :doc:`../../modules/assimilation/assim_tools_mod` documentation.

Detailed program execution flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Manhattan release of DART includes state space output expanded from the 
previous two stages (Prior and Posterior) to up to six (input, forecast, preassim, 
postassim, analysis, and output). This makes it possible to examine the states with 
and without either kind of inflation, as described below. In addition, the state 
space vectors are each written to a separate NetCDF file: 
``${stage}_mean.nc, ${stage}_sd.nc, ${stage}_member_####.nc`` .
The detailed execution flow inside the filter program is:

-  Read in observations.
-  Read in state vectors from model netcdf restart files.
-  Initialize inflation fields, possibly reading netcdf restart files.
-  If requested, initialize and write to "input" netcdf diagnostic files.
-  Trim off any observations if start/stop times specified.
-  Begin main assimilation loop:

   -  Check model time vs observation times:

      -  If current assimilation window is earlier than model time, error.
      -  If current assimilation window includes model time, begin assimilating.
      -  If current assimilation window is later than model time, advance model:

         -  Write out current state vectors for all ensemble members.
         -  Advance the model by subroutine call or by shell script:

            -  Tell the model to run up to the requested time.

         -  Read in new state vectors from netcdf files for all ensemble members.

   -  Apply prior inflation if requested.
   -  Compute ensemble of prior observation values with forward operators.
   -  If requested, compute and write the "preassim" netcdf diagnostic files. This is AFTER any prior inflation has been
      applied.
   -  Compute prior observation space diagnostics.
   -  Assimilate all observations in this window:

      -  Get all obs locations and kinds.
      -  Get all state vector locations and kinds.
      -  For each observation:

         -  Compute the observation increments.
         -  Find all other obs and states within localization radius.
         -  Compute the covariance between obs and state variables.
         -  Apply increments to state variables weighted by correlation values.
         -  Apply increments to any remaining unassimilated observations.
         -  Loop until all observations in window processed.

   -  If requested, compute and write the "postassim" netcdf diagnostic files (members, mean, spread). This is BEFORE
      any posterior inflation has been applied.
   -  Apply posterior inflation if requested.
   -  Compute ensemble of posterior observation values with forward operators.
   -  Compute posterior observation space diagnostics.
   -  If requested, compute and write out the "output" netcdf diagnostic files (members, mean, spread). This is AFTER
      any posterior inflation has been applied.
   -  Loop until all observations in input file processed.

-  Close diagnostic files.
-  Write out final observation sequence file.
-  Write out inflation restart files if requested.
-  Write out final state vectors to model restart files if requested.
-  Release memory for state vector and observation ensemble members.

Namelist
--------

See the `filter namelist <../../modules/assimilation/filter_mod.html#Namelist>`__ page for a detailed description of all
``&filter_nml`` variables. This namelist is read from the file ``input.nml``.

Modules used
------------

::

   mpi_utilities_mod
   filter_mod

Note that `filter_mod.f90 <../../modules/assimilation/filter_mod.html#Modules>`__ uses many more modules.

Files
-----

See Detailed Program Flow for a short description of DART's new 'stages'. In addition, the Manhattan release simplifies
some namelists by replacing many user-settable file names with hardwired filenames. Files can then be renamed in the run
scripts to suit the user's needs.

-  input ensemble member states; from *&filter_nml :: input_state_files* or *input_state_file_list*
-  output ensemble member states; to *&filter_nml :: output_state_files* or *output_state_file_list*
-  input observation sequence file; from ``&filter_nml :: obs_sequence_in_name``
-  output observation sequence file; from ``&filter_nml :: obs_sequence_out_name``
-  output state space diagnostics files; ``${stage}_mean.nc, ${stage}_sd.nc,`` where stage =
   {input,forecast,preassim,postassim,analysis,output}
-  input state space inflation data (if enabled); from ``input_{prior,post}inf_{mean,sd}.nc.``
-  output state space inflation data (if enabled); to ``${stage}_{prior,post}inf_{mean,sd}.nc.``, where stage ≠ "input"
-  input.nml, to read &filter_nml

References
----------

-  Anderson, J. L., 2001: An Ensemble Adjustment Kalman Filter for Data Assimilation. Mon. Wea. Rev., 129, 2884-2903.
   `doi:
   10.1175/1520-0493(2001)129<2884:AEAKFF>2.0.CO;2 <http://dx.doi.org/10.1175/1520-0493%282001%29129%3C2884%3AAEAKFF%3E2.0.CO%3B2>`__
-  Anderson, J. L., 2003: A Local Least Squares Framework for Ensemble Filtering. Mon. Wea. Rev., 131, 634-642.
   `doi:
   10.1175/1520-0493(2003)131<0634:ALLSFF>2.0.CO;2 <http://dx.doi.org/10.1175/1520-0493%282003%29131%3C0634%3AALLSFF%3E2.0.CO%3B2>`__
-  Anderson, J. L., 2007: An adaptive covariance inflation error correction algorithm for ensemble filters. Tellus A,
   59, 210-224.
   `doi: 10.1111/j.1600-0870.2006.00216.x <http://dx.doi.org/10.1111/j.1600-0870.2006.00216.x>`__
-  Anderson, J. L., 2007: Exploring the need for localization in ensemble data assimilation using a hierarchical
   ensemble filter. Physica D, 230, 99-111.
   `doi:10.1016/j.physd.2006.02.011 <http://dx.doi.org/10.1016/j.physd.2006.02.011>`__
-  Anderson, J., Collins, N., 2007: Scalable Implementations of Ensemble Filter Algorithms for Data Assimilation.
   Journal of Atmospheric and Oceanic Technology, 24, 1452-1463.
   `doi: 10.1175/JTECH2049.1 <http://dx.doi.org/10.1175/JTECH2049.1>`__
-  Anderson, J. L., 2009: Spatially and temporally varying adaptive covariance inflation for ensemble filters. Tellus A,
   61, 72-83.
   `doi: 10.1111/j.1600-0870.2008.00361.x <http://dx.doi.org/10.1111/j.1600-0870.2008.00361.x>`__
-  Anderson, J., T. Hoar, K. Raeder, H. Liu, N. Collins, R. Torn, and A. Arellano, 2009: The Data Assimilation Research
   Testbed: A Community Facility. Bull. Amer. Meteor. Soc., 90, 1283-1296.
   `doi: 10.1175/2009BAMS2618.1 <http://dx.doi.org/10.1175/2009BAMS2618.1>`__
-  Anderson, J. L., 2010: A Non-Gaussian Ensemble Filter Update for Data Assimilation. Mon. Wea. Rev., 139, 4186-4198.
   `doi: 10.1175/2010MWR3253.1 <http://dx.doi.org/10.1175/2010MWR3253.1>`__
-  Anderson, J. L., 2011: Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation.
   Submitted for publication, Jan 2011. Contact author.
