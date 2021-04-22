Inflation
=========

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

Historically, inflation was first introduced to address sampling errors (the fact
that we are limited to a small ensemble size).
Latest research, e.g. `El Gharamti et al. (2019) <https://doi.org/10.1175/MWR-D-18-0389.1>`__ 
suggests that prior and posterior inflation can be used to address different issues
in the filtering problem. Prior inflation is able to address issues in the forecast 
step such as model errors while posterior inflation can help mitigate sampling errors 
in the analysis step. 


Inflation values can vary in space and time, depending on the specified namelist values. Even though we talk about a
single inflation value, the inflation has a probability density with a mean and standard deviation. We use the mean
value when we inflate, and the standard deviation indicates how sure of the value we are. Larger standard deviation
values mean "less sure" and the inflation value can increase more quickly with time. Smaller values mean "more sure" and
the time evolution will be slower since we are more confident that the mean (inflation value) is correct.

The standard deviation of inflation allows inflation values to increase with time, if required by increasing density or
frequency of observations, but it does not provide a mechanism to reduce the inflation when the frequency or density of
observations declines. So there is also an option to damp inflation through time. In practice with large geophysical
models using damped inflation has been a successful strategy.

The following namelist items which control inflation are found in the ``input.nml`` file, in the &filter_nml namelist.
The detailed descriptions are in the `filter_mod <../assimilation_code/modules/assimilation/filter_mod.html#Namelist>`__ page. Here we
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
   If time-evolution is enabled, each value can evolve independently. Spatially-uniform state space inflation uses a
   single inflation value for all items in the state vector. If time-evolution is enabled, that single value can evolve.
   See ``inf_sd_*`` below for control of the time-evolution behavior. Enhanced spatially-varying inflation uses an
   inverse-gamma distribution which allows the standard deviation of the inflation to increase or decrease through time
   and may produce better results (see `El Gharamti (2018) <https://doi.org/10.1175/MWR-D-17-0187.1>`__). 
   In practice we recommend starting with no inflation (both values 0). Then try
   inflation type 2 or 5 prior inflation and no inflation (0) for posterior. WARNING: even if inf_flavor is not 0,
   inflation will be turned off if ``inf_damping`` is set to 0.

   .. important::
   
       Relaxation to prior spread (aka RTPS, i.e., ``inf_flavor=4``) is a 
       spatially varying **posterior** inflation algorithm.  


   When using RTPS you cannot set the prior inflation 
   flavor to 4. The code will exit with an error messge. Unlike all other flavors, RTPS does 
   not use files to handle inflation in time. So, if the user supplies ``input_postinf_{mean,sd}.nc``, 
   these will be **ignored**. The ONLY namelist option that RTPS uses (other than ``inf_flavor=4``)
   is the second entry of ``inf_initial``. This value is technically not the 
   posterior inflation value but rather a *weighting* factor (denoted by :math:`{\alpha}`; in 
   `Whitaker and Hamill (2012) <https://doi.org/10.1175/MWR-D-11-00276.1>`__)
   that is used to relax the posterior spread to the prior spread. For instance, if :math:`{\alpha}=0.3`
   then the inflated posterior spread is as follows: 70% of the analysis spread plus
   30% of the prior spread. If :math:`{\alpha}=1.0`, then the inflated posterior spread is simply set 
   to the prior spread. Using :math:`{\alpha}`, RTPS calculates the effective posterior inflation *under the hood*
   and writes out the inflation values to the user. These can be looked at for diagnostic purposes. 
   The algorithm disregards them for the next data assimilation cycle. In short, RTPS is 
   adaptive in time but unlike flavors 2, 3 and 5 it has no memory. 
   The recommendation is to set the second entry of ``inf_initial``
   to any number between 0.0 and 1.0.     

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

First and foremost, if you are using one of the temporally-varying inflation options, 
save the entire series of inflation files to explore how inflation evolves through time.
As part of the workflow, you have to take the output of one assimilation cycle and rename
it to be the input for the next assimilation cycle. That is the time to make a copy 
that has a unique name - usually with some sort of date or timestamp. This also makes
it possible to restart an experiment.

The suggested procedure for testing inflation options is to start without any (both ``inf_flavor`` values set to 0 and
``inf_damping`` > 0.). Then enable Prior state space, spatially-varying inflation, with no Posterior inflation (set
``inf_flavor`` to [2, 0]). Then try damped inflation (set ``inf_damping`` to 0.9 and set ``inf_sd_initial`` and
``inf_sd_lower_bound`` to 0.6). The inflation values and standard deviation are written out to files with
``_{prior,post}inf_{mean,sd}`` in their names. These NetCDF files can be viewed with common tools (we often use
`ncview <http://meteora.ucsd.edu/~pierce/ncview_home_page.html>`__ ). Expected inflation values are generally in the 1
to 30 range; if values grow much larger than this it usually indicates a problem with the assimilation.

:doc:`../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart` 
may be used to create netCDF files with initial values such that the 
input.nml settings for reading from file vs. reading from namelist can stay constant 
throughout the entire experiment.

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
