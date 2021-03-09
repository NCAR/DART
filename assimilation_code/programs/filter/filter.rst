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

See the `DART web site <http://www.image.ucar.edu/DAReS/DART>`__ for more documentation, including a discussion of the
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
model can get very large. A lightly pruned directory tree can be browsed in the main
`index.html <../../../docs/index.html#Directories>`__.

Types of filters available
~~~~~~~~~~~~~~~~~~~~~~~~~~

The different types of assimilation algorithms (EAKF, ENKF, Kernel filter, Particle filter, etc.) are determined by the
``&assim_tools_nml:filter_kind`` entry, described in :doc:`../../modules/assimilation/assim_tools_mod`. Despite having
'filter' in the name, they are assimilation algorithms and so are implemented in ``assim_tools_mod.f90``.

Getting started
~~~~~~~~~~~~~~~

Running a successful assimilation takes careful diagnostic work and experiment iterations to find the best settings for
your specific case. The basic Kalman filter can be coded in only a handful of lines; the hard work is making the right
choices to compensate for sampling errors, model bias, observation error, lack of model forecast divergence, variations
in observation density in space and time, random correlations, etc. There are tools built into DART to deal with most of
these problems but it takes careful work to apply them correctly.

If you are adding a new model or a new observation type, we suggest you assimilate exactly one observation, with no
model advance, with inflation turned off, with a large cutoff, and with the outlier threshold off (see below for how to
set these namelist items). Run an assimilation. Look at the ``obs_seq.final`` file to see what the forward operator
computed. Use ncdiff to difference the ``preassim_mean.nc`` and ``postassim_mean.nc`` (or ``output_mean.nc``) diagnostic
NetCDF files and look at the changes (the "innovations") in the various model fields. Is it in the right location for
that observation? Does it have a reasonable value?

Then assimilate a group of observations and check the results carefully. Run the observation diagnostics and look at the
total error and spread. Look carefully at the number of observations being assimilated compared to how many are
available. Assimilations that are not working can give good looking statistics if they reject all but the few
observations that happen to match the current state. The errors should grow as the model advances and then shrink when
new observations are assimilated, so a timeseries plot of the RMSE should show a sawtooth pattern. The initial error
entirely depends on the match between the initial ensemble and the observations and may be large but it should decrease
and then reach a roughly stable level. The ensemble spread should ultimately remain relatively steady, at a value around
the expected observation error level. Once you believe you have a working assimilation, this will be your baseline case.
If the ensemble spread is too small, several of the DART facilities described below are intended to compensate for
ensemble members getting too close to each other. Then one by one enable or tune each of the items below, checking each
time to see what is the effect on the results.

Suggestions for the most common namelist settings and features built into DART for running a successful assimilation
include:

Ensemble size
^^^^^^^^^^^^^

In practice, ensemble sizes between 20 and 100 seem to work best. Fewer than 20-30 members leads to statistical errors
which are too large. More than 100 members takes longer to run with very little benefit, and eventually the results get
worse again. Often the limit on the number of members is based on the size of the model since you have to run N copies
of the model each time you move forward in time. If you can, start with 50-60 members and then experiment with fewer or
more once you have a set of baseline results to compare it with. The namelist setting for ensemble size is
``&filter_nml :: ens_size``

Localization
^^^^^^^^^^^^

There are two main advantages to using localization. One is it avoids an observation impacting unrelated state variables
because of spurious correlations. The other is that, especially for large models, it improves run-time performance
because only points within the localization radius need to be considered. Because of the way the parallelization was
implemented in DART, localization was easy to add and using it usually results in a very large performance gain. See
`here <../../modules/assimilation/assim_tools_mod.html#Localization>`__ for a discussion of localization-related
namelist items.

Inflation
^^^^^^^^^

Since the filter is run with a number of members which is usually small compared to the number of degrees of freedom of
the model (i.e. the size of the state vector or the number of EOFs needed to characterize the variability), the model
uncertainty is under-represented. Other sources of error and uncertainty are not represented at all. These factors lead
to the ensemble being 'over-confident', or having too little spread. More observations leads to more over-confidence.
This characteristic can worsen with time, leading to ensemble collapse to a single solution. Inflation increases the
spread of the members in a systematic way to overcome this problem. There are several sophisticated options on
inflation, including spatial and temporal adaptive and damping options, which help deal with observations which vary in
density over time and location. See here for a discussion of inflation-related namelist items.

Outlier rejection
^^^^^^^^^^^^^^^^^

Outlier rejection can be used to avoid bad observations (ones where the value was recorded in error or the processing
has an error and a non-physical value was generated). It also avoids observations which have accurate values but the
mean of the ensemble members is so far from the observation value that assimilating it would result in unacceptably
large increments that might destablize the model run. If the difference between the observation and the prior ensemble
mean is more than N standard deviations from the square root of the sum of the prior ensemble and observation error
variance, the observation will be rejected. The namelist setting for the number of standard deviations to include is
``&filter_nml :: outlier_threshold`` and we typically suggest starting with a value of 3.0.

Sampling error
^^^^^^^^^^^^^^

For small ensemble sizes a table of expected statistical error distributions can be generated before running DART.
Corrections accounting for these errors are applied during the assimilation to increase the ensemble spread which can
improve the assimilation results. The namelist item to enable this option is
``&assim_tools_nml :: sampling_error_correction``. Additionally you will need to have the precomputed correction file
``sampling_error_correction_table.nc``, in the run directory. See the description of the namelist item in the
`&assim_tools_nml <../../modules/assimilation/assim_tools_mod.html#Namelist>`__ namelist, and
:doc:`../system_simulation/system_simulation` for instructions on where to find (or how to generate) the auxiliary file
needed by this code. See Anderson (2011).

Free run/forecast after assimilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Separate scripting can be done to support forecasts starting from the analyzed model states. After filter exits, the
models can be run freely (with no assimilated data) further forward in time using one or more of the last updated model
states from filter. Since all ensemble members are equally likely a member can be selected at random, or a member close
to the mean can be chosen. See the :doc:`../../../assimilation_code/programs/closest_member_tool/closest_member_tool`
for one way to select a "close" member. The ensemble mean is available to be used, but since it is a combination of all
the member states it may not have self-consistent features, so using a single member is usually preferred.

Evaluating observations without assimilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Filter can be used to evaluate the accuracy of a single model state based on a set of available observations. Either
copy or link the model state file so there appear to be 2 separate ensemble members (which are identical). Set the
filter namelist ensemble size to 2 by setting ``ens_size`` to 2 in the &filter_nml namelist. Turn off the outlier
threshold and both Prior and Posterior inflation by setting ``outlier_threshold`` to -1, and both the ``inf_flavor``
values to 0 in the same &filter_nml namelist. Set all observation types to be 'evaluate-only' and have no types in the
'assimilate' list by listing all types in the ``evaluate_these_obs_types`` list in the ``&obs_kind_nml`` section of the
namelist, and none in the assimilation list. Run filter as usual, including model advances if needed. Run observation
diagnostics on the resulting ``obs_seq.final`` file to compute the difference between the observed values and the
predicted values from this model state.

Verification/comparison with and without assimilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compare results of an experiment with and without assimilating data, do one run assimilating the observations. Then
do a second run where all the observation types are moved to the ``evaluate_these_obs_types`` list in the
``&obs_kind_nml`` section of the namelist. Also turn inflation off by setting both ``inf_flavor`` values to 0 in the
&filter_nml namelist. The forward operators will still be called, but they will have no impact on the model state. Then
the two sets of diagnostic state space netcdf files can be compared to evaluate the impact of assimilating the
observations, and the observation diagnostic files can also be compared.

DART quality control flag added to output observation sequence file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The filter adds a quality control field with metadata 'DART quality control' to the ``obs_seq.final`` file. At present,
this field can have the following values:

== =====================================================================================================================
0: Observation was assimilated successfully
1: Observation was evaluated only but not used in the assimilation
2: The observation was used but one or more of the posterior forward observation operators failed
3: The observation was evaluated only but not used AND one or more of the posterior forward observation operators failed
4: One or more prior forward observation operators failed so the observation was not used
5: The observation was not used because it was not selected in the namelist to be assimilated or evaluated
6: The prior quality control value was too high so the observation was not used.
7: Outlier test failed (see below)
== =====================================================================================================================

The outlier test computes the difference between the observation value and the prior ensemble mean. It then computes a
standard deviation by taking the square root of the sum of the observation error variance and the prior ensemble
variance for the observation. If the difference between the ensemble mean and the observation value is more than the
specified number of standard deviations, then the observation is not used and the DART quality control field is set to
7.

Discussion of inflation options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In pre-Manhattan DART, there were two choices for the basic type of inflation: observation-space or state-space.
Observation-space inflation is no longer supported. (If you are interested in observation-space inflation, talk to Jeff
first.) The rest of this discussion applies to state-space inflation.

| State-space inflation changes the spread of an ensemble without changing the ensemble mean. The algorithm computes the
  ensemble mean and standard deviation for each variable in the state vector in turn, and then moves the member's values
  away from the mean in such a way that the mean remains unchanged. The resulting standard deviation is larger than
  before. It can be applied to the Prior state, before observations are assimilated (the most frequently used case), or
  it can be applied to the Posterior state, after assimilation. See `Anderson
  (2007) <http://dx.doi.org/10.1175/JTECH2049.1>`__, `Anderson
  (2009) <http://dx.doi.org/10.1111/j.1600-0870.2008.00361.x>`__.

Inflation values can vary in space and time, depending on the specified namelist values. Even though we talk about a
single inflation value, the inflation has a gaussian distribution with a mean and standard deviation. We use the mean
value when we inflate, and the standard deviation indicates how sure of the value we are. Larger standard deviation
values mean "less sure" and the inflation value can increase more quickly with time. Smaller values mean "more sure" and
the time evolution will be slower since we are more confident that the mean (inflation value) is correct.

The standard deviation of inflation allows inflation values to increase with time, if required by increasing density or
frequency of observations, but it does not provide a mechanism to reduce the inflation when the frequency or density of
observations declines. So there is also an option to damp inflation through time. In practice with large geophysical
models using damped inflation has been a successful strategy.

The following namelist items which control inflation are found in the ``input.nml`` file, in the &filter_nml namelist.
The detailed descriptions are in the `namelist <../../modules/assimilation/filter_mod.html#Namelist>`__ page. Here we
try to give some basic advice about commonly used values and suggestions for where to start. Spatial variation is
controlled by ``inf_flavor``, which also controls whether there's any inflation, ``inf_initial_from_restart``, and
``inf_initial``, as described below. Time variation is controlled by ``inf_sd_initial_from_restart``,
``inf_sd_initial``, ``inf_sd_lower_bound``, ``inf_damping``, ``inf_lower_bound`` and ``inf_upper_bound``.

In the namelist each entry has two values. The first is for Prior inflation and the second is for Posterior inflation.

``&filter_nml :: inf_flavor``
   *valid values:*\ 0, 2, 3, 4, 5

   Set the type of Prior and Posterior inflation applied to the state vector. Values mean:

   -  **0:** No inflation (Prior and/or Posterior) and all other inflation variables are ignored
   -  [**1:** Deprecated: Observation space inflation]
   -  **2:** Spatially-varying state space inflation (gaussian)
   -  **3:** Spatially-uniform state space inflation (gaussian)
   -  **4:** Relaxation To Prior Spread (Posterior inflation only)
   -  **5:** Enhanced Spatially-varying state space inflation (inverse gamma)

   Spatially-varying state space inflation stores an array of inflation values, one for each item in the state vector.
   If time-evolution is enabled each value can evolve independently. Spatially-uniform state space inflation uses a
   single inflation value for all items in the state vector. If time-evolution is enabled that single value can evolve.
   See ``inf_sd_*`` below for control of the time-evolution behavior. Enhanced spatially-varying inflation uses an
   inverse-gamma distribution which allows the standard deviation of the inflation to increase or decrease through time
   and may produce better results. In practice we recommend starting with no inflation (both values 0). Then try
   inflation type 2 or 5 prior inflation and no inflation (0) for posterior. WARNING: even if inf_flavor is not 0,
   inflation will be turned off if ``inf_damping`` is set to 0.

``&filter_nml :: inf_initial_from_restart``
   *valid values:* .true. or .false.

   If true, read the inflation values from an inflation restart file named ``input_{prior,post}inf_mean.nc.`` An initial
   run could be done to let spatially-varying inflation values evolve in a spinup phase, and then the saved values can
   be read back in and used as fixed values in further runs. Or if time-varying inflation is used, then the restart file
   from the previous job step must be supplied as an input file for the next step.

``&filter_nml :: inf_initial``
   *valid values:* real numbers, usually 1.0 or slightly larger
   If not reading in inflation values from a restart file, the initial value to set for the inflation. Generally we
   recommend starting with just slightly above 1.0, maybe 1.02, for a slight amount of initial inflation.
``&filter_nml :: inf_lower_bound``
   *valid values:* real numbers, usually 1.0 or slightly larger

   If inflation is time-evolving (see ``inf_sd_*`` below), then this sets the lowest value the inflation can evolve to.
   Setting a number less than one allows for deflation but generally in a well-observed system the ensemble needs more
   spread and not less. We recommend a setting of 1.0.

``&filter_nml :: inf_upper_bound``
   *valid values:* real numbers, larger than 1.0

   If inflation is time-evolving (see ``inf_sd_*`` below), then this sets the largest value the inflation can evolve to.
   We recommend a setting of 100.0, although if the inflation values reach those levels there is probably a problem with
   the assimilation.

``&filter_nml :: inf_damping``
   *valid values:* 0.0 to 1.0

   Applies to all state-space inflation types, but most frequently used with time-adaptive inflation variants. The
   difference between the current inflation value and 1.0 is multiplied by this factor before the next assimilation
   cycle. So the inflation values are pushed towards 1.0, from above or below (if inf_lower_bound allows inflation
   values less than 1.0). A value of 0.0 turns all inflation off by forcing the inflation value to 1.0. A value of 1.0
   turns damping off by leaving the original inflation value unchanged. We have had good results in large geophysical
   models using time- and space-adaptive state-space inflation and setting the damping to a value of 0.9, which damps
   slowly.

``&filter_nml :: inf_sd_initial_from_restart``
   *valid values:* .true. or .false.

   If true, read the inflation standard deviation values from an restart file named ``input_{prior,post}inf_sd.nc.`` See
   the comments above about ``inflation_initial_from_restart``.

``&filter_nml :: inf_sd_initial``
   *valid values:* ≤ 0.0 to disable evolution of inflation, > 0.0 otherwise

   The initial value to set for the inflation standard deviation, IF not reading in inflation standard deviation values
   from a file. This value (or these values) control whether the inflation values evolve with time or not. A negative
   value or 0.0 prevents the inflation values from being updated, so they are constant throughout the run. If positive,
   the inflation values evolve through time. We have had good results setting this and ``inf_sd_lower_bound`` to 0.6 for
   large geophysical models.

``&filter_nml :: inf_sd_lower_bound``
   *valid values:* ≤ 0.0 to disable evolution of inflation, > 0.0 otherwise

   If the setting of ``inf_sd_initial`` is ≤ 0 (to disable time evolution of inflation) then set this to the same value.

   Otherwise, the standard deviation of the inflation cannot fall below this value. Smaller values will restrict the
   inflation to vary more slowly with time; larger values will allow the inflation to adapt more quickly. We have had
   good results setting this and ``inf_sd_initial`` to 0.6 for large geophysical models. Since the
   ``inf_sd_lower_bound`` is a scalar, it is not possible to set different lower bounds for different parts of the state
   vector.

   Time-varying inflation with flavor 2 generally results in the inflation standard deviation for all state variables
   shrinking to the lower bound and staying there. For flavor 5, the inflation standard deviation value is allowed to
   increase and decrease.

``&filter_nml :: inf_sd_max_change``
   *valid values:* 1.0 to 2.0

   Used only with the Enhanced inflation (flavor 5). The Enhanced inflation algorithm allows the standard deviation to
   increase as well as decrease. The ``inf_sd_max_change`` controls the maximum increase of the standard deviation in an
   assimilation cycle. A value of 1.0 means it will not increase, a value of 2.0 means it can double; a value inbetween
   sets the percentage it can increase, e.g. 1.05 is a limit of 5%. Suggested value is 1.05 (max increase of 5% per
   cycle).

   Because the standard deviation for original flavor 2 could never increase, setting the ``inf_sd_initial`` value equal
   to the ``inf_sd_lower_bound`` value effectively fixed the standard deviation at a constant value. To match the same
   behavior, if they are equal and Enhanced inflation (flavor 5) is used it will also use that fixed value for the
   standard deviation of the inflation. Otherwise the standard deviation will adapt as needed during each assimilation
   cycle.

``&filter_nml :: inf_deterministic``
   *valid values:* .true. or .false.

   Recommend always using ``.true.``.

Guidance regarding inflation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The suggested procedure for testing inflation options is to start without any (both ``inf_flavor`` values set to 0 and
``inf_damping`` > 0.). Then enable Prior state space, spatially-varying inflation, with no Posterior inflation (set
``inf_flavor`` to [2, 0]). Then try damped inflation (set ``inf_damping`` to 0.9 and set ``inf_sd_initial`` and
``inf_sd_lower_bound`` to 0.6). The inflation values and standard deviation are written out to files with
``_{prior,post}inf_{mean,sd}`` in their names. These NetCDF files can be viewed with common tools (we often use
`ncview <http://meteora.ucsd.edu/~pierce/ncview_home_page.html>`__ ). Expected inflation values are generally in the 1
to 30 range; if values grow much larger than this it usually indicates a problem with the assimilation.

It is possible to set inflation values in an existing netCDF file by using one of the standard NCO utilities like
"``ncap2``" on a copy of a restart file. Inflation mean and sd values look exactly like restart values, arranged by
variable type like T, U, V, etc.

Here's an example of using ncap2 to set the T,U and V inf values:

.. container:: unix

   ::

        ncap2 -s 'T=1.0;U=1.0;V=1.0' wrfinput_d01 input_priorinf_mean.nc
        ncap2 -s 'T=0.6;U=0.6;V=0.6' wrfinput_d01 input_priorinf_sd.nc
        -or-
        ncap2 -s 'T(:,:,:)=1.0;U(:,:,:)=1.0;V(:,:,:)=1.0' wrfinput_d01 input_priorinf_mean.nc
        ncap2 -s 'T(:,:,:)=0.6;U(:,:,:)=0.6;V(:,:,:)=0.6' wrfinput_d01 input_priorinf_sd.nc

Some versions of the NCO utilities change the full 3D arrays into a single scalar. If that's your result (check your
output with ``ncdump -h``) use the alternate syntax or a more recent version of the NCO tools.

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
   directories for ideas. See additional documentation in the :doc:`../../../models/template/model_mod` documentation,
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

The Manhattan release of DART includes state space output expanded from the previous two stages (Prior and Posterior) to
up to four (input, preassim, postassim, and output). This makes it possible to examine the states with and without
either kind of inflation, as described below. In addition, the state space vectors are each written to a separate NetCDF
file: ``${stage}_mean.nc, ${stage}_sd.nc, ${stage}_member_####.nc`` . The detailed execution flow inside the filter
program is:

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
   {input,preassim,postassim,output}
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
