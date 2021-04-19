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

Your first attempt
------------------

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

Next attempts
------------------

:doc:`High-level data assimilation workflows <high-level-da-workflows>`
gives an overview of a variety of complete assimilation experiments,
including the programs which need to be run and their input and output.

Important features of assimilations
-----------------------------------

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
density over time and location. 
See :doc:`Inflation <inflation>` for a discussion of inflation-related namelist items.

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
:doc:`../assimilation_code/programs/gen_sampling_err_table/gen_sampling_err_table` 
for instructions on where to find (or how to generate) the auxiliary file
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
