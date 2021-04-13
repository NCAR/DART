Assimilation in a complex model
===============================

Introduction
------------

Running a successful assimilation takes careful diagnostic work and experiment
iterations to find the best settings for your specific case.

The basic Kalman filter can be coded in only a handful of lines. The difficulty
in getting an assimilation system working properly involves making the right
choices to compensate for sampling errors, model bias, observation error, lack
of model forecast divergence, variations in observation density in space and
time, random correlations, etc. There are tools built into DART to deal with
most of these problems but it takes careful work to apply them correctly.

This document guides you through the process of using DART with your model. It
uses a questionnaire to guide you through the questions you'll need to answer
in order to get a data assimilation system working.

Is your model appropriate for any kind of DA?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your model isn't chaotic, you don't need DA. You run the model, you look at
the difference between the prediction and the observations, and you invert the
equations inside the model to compute what different inputs would have
produced outputs closer to the observations.

Chaotic models don't have a simple relationship between inputs and outputs.
There are internal feedbacks and non-linear behaviors that make it difficult
to adjust the inputs to make the outputs better match the observations.

What is your model state?
~~~~~~~~~~~~~~~~~~~~~~~~~

"Model state" has a specific definition that can be the source of much
confusion if you are running a model and haven't though about DA before.
Formally, it is the minimal set of variables that must be saved when a model
stops so it can be restarted again exactly.

At first glance this means all the variables on the right side of
the equals sign for the governing equations of the system.  However
many models which have not been designed with DA in mind may have
no clear time when all parts of the model are at a consistent time.
For example, some variables may be 1/2 timestep ahead or behind others.
Or some derived variables may be expensive to compute and so are
precomputed and stored and not recomputed. If the DA process changes
the state variables all derived variables must be recomputed before
proceeding.

Restart files often store many more variables than the minimal set
needed to restart the model. Often other variables are used in 
diagnostic routines or are of interest on their own. Generally
these aren't considered part of the model state.

How is your model execution controlled?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In general, larger and more complex models have an environment they
are expecting to run within.  This environment includes scripts to control
the execution parameters or input parameter files, how many processors are
used in a parallel system, how the tasks are distributed over the hardware,
how an execution is expected to run in model time, and what variables are
written to the output files.

For DA, there must, at a minimum, be a way to control how long the model 
runs before it writes out the results and exits.  

For large models, the DA filter process is a large parallel program
generally requiring a multi-processor supercomputer or cluster.  Many
models themselves are large parallel programs, so there can be issues
with how the switch between model and DA process is done.

New or adjusted scripting is generally required to include the DA process
in the overall execution flow.

Are you able to start and stop your model at specific times?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DA process is generally a cycle of running the model for a certain 
amount of model time, then running the DA filter to adjust the model 
state before continuing.

These two steps happen over and over as observations are available to
guide the adjustments to the model state.

Models may be written with the assumption that startup costs are
only done once and then the model runs for a long period of time.  
When used with DA models are generally started and stopped after 
running a relatively short amount of model time.  If model startup 
time is long this can result in unacceptably slow performance.

A small amount of round-off error is often introduced when a model 
writes restart files before stopping.  So running a model N timesteps 
forward vs. running N/2, stopping, writing restart files, starting, 
reading restart files, and finishing the last N/2 timesteps will 
usually not result in identical values.

The goal is to minimize the differences.  This can require small or
large changes to make the model behave as expected with repeated 
starting and stopping.

Some models include external forcing, for example boundary conditions
from a separate model.  If cycling the forcing files may need to be
updated periodically outside of the DA system.

What coordinate system is used by your model?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coordinate systems use a series of numbers to describe the
relationship in space between parts of the model state and
where observations are located.  In Earth-system models,
often a latitude-longitude-vertical coordinate system
is used.  X,Y,Z Cartesian coordinates are also used to describe
3D space.  Other options include cyclindrical or spherical coordinates,
and unit-line, -square or -cube coordinates with cyclical boundaries.

Only a single coordinate system can be selected and it applies to
both the model state locations as well as the observations.

If the model coordinate system is based on some other space
it may be necessary to transform it into physical coordinates
before running DA.  For example, some models compute in spectral
space and the output must be translated into a physical space
before DA can be done.

What file format is used for model restart files?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DART reads and writes NetCDF file format.  Many earth-system models
already use this format.  If the model does not, converter programs
from the native format to NetCDF and back are needed.  NetCDF is a
self-describing format with metadata that allows DART to read and
process model data without additional configuration files.

What quantities are in the model state?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DART defines a "Quantity" as the fundamental physical object
a value is measuring.  Examples are Temperature, Pressure,
Salinity, etc.  Each value in a model state must be 
associated with a defined quantity.

What observations are you intending to assimilate?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any observation you intend to assimilate requires a method to
compute an "expected value" based on the model state.  Often
the observation is of the same quantity as exists in the model
state, so computing the expected value is a direct process.

Other times the expected value is a function of quantities in
the model state, and code called a "forward operator" uses
one or more quantities from the model state and computes the
expected value.

If the model state does not contain quantities that are needed
to compute an expected value, auxiliary data values can be read
and used to compute the expected value.  But if the expected value
cannot be computed or is not in some way a function of the model
state, the observations cannot be assimilated.

How are you going to generate your initial ensemble?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most models don't have an existing ensemble of states ready
for ingestion into an ensemble DA system. Options for generating
the initial ensemble include adding random perturbations to a 
single variable in a single state, perturbing forcing variables
differently for each ensemble member, or perturbing the entire state.

For models which have a lot of error growth, it may be enough to
add a very small amount of noise to a single variable in the state
to generate an ensemble of states and then run them forward in time
with the model to generate states which have sufficient differences.

For models with slower error growth, larger perturbations may be
needed, a longer model advance time before starting assimilation, 
or perturbations of forcing or boundary files may be needed.

The goal is to generate a set of model states which are different
but contain internally-consistent values.  

An ensemble of states without sufficient differences (spread) will
reject assimilating observations.

General advice
--------------

If you are adding a new model or a new observation type, you should assimilate
exactly one observation, with no model advance, with inflation turned off, with
a large cutoff, and with the outlier threshold off (see below for how to
set these namelist items).

Run an assimilation. Look at the ``obs_seq.final`` file to see what the forward
operator computed. Use ncdiff to difference the ``preassim_mean.nc`` and
``postassim_mean.nc`` (or ``output_mean.nc``) diagnostic NetCDF files and look
at the changes (the "innovations") in the various model fields. Is it in the
right location for that observation? Does it have a reasonable value?

Then assimilate a group of observations and check the results carefully. Run
the observation diagnostics and look at the total error and spread. Look
carefully at the number of observations being assimilated compared to how many
are available.

Assimilations that are not working can give good looking statistics if they
reject all but the few observations that happen to match the current state.
The errors should grow as the model advances and then shrink when new
observations are assimilated, so a timeseries plot of the RMSE should show a
sawtooth pattern. The initial error entirely depends on the match between the
initial ensemble and the observations and may be large but it should decrease
and then reach a roughly stable level. The ensemble spread should ultimately
remain relatively steady, at a value around the expected observation error
level. Once you believe you have a working assimilation, this will be your
baseline case.

If the ensemble spread is too small, several of the DART facilities described
below are intended to compensate for ensemble members getting too close to each
other. Then one by one enable or tune each of the items below, checking each
time to see what is the effect on the results.

Suggestions for the most common namelist settings and features built into DART
for running a successful assimilation include:

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
`here <../assimilation_code/modules/assimilation/assim_tools_mod.html#Localization>`__ for a discussion of localization-related
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
`&assim_tools_nml <../assimilation_code/modules/assimilation/assim_tools_mod.html#Namelist>`__ namelist, and
:doc:`../assimilation_code/programs/system_simulation/system_simulation` for instructions on where to find (or how to generate) the auxiliary file
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
1: Observation was evaluated (as specified in namelist) and not used in the assimilation
2: The observation was used but one or more of the posterior forward observation operators failed
3: The observation was evaluated AND one or more of the posterior forward observation operators failed
4: One or more prior forward observation operators failed so the observation was not used
5: The observation was not used because it was not selected in the namelist to be assimilated or evaluated
6: The prior quality control value was too high so the observation was not used.
7: Outlier test failed (see below)
8: Vertical conversion failed
== =====================================================================================================================

The outlier test computes the difference between the observation value and the prior ensemble mean. It then computes a
standard deviation by taking the square root of the sum of the observation error variance and the prior ensemble
variance for the observation. If the difference between the ensemble mean and the observation value is more than the
specified number of standard deviations, then the observation is not used and the DART quality control field is set to
7.
