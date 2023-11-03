.. _Using new models:

Can I run my model with DART?
=============================

The DART team often collaborates with other groups to help write the interface
code to a new model. The most efficient way to get started is to meet with
DAReS staff either virtually or in person, to discuss what is involved in
supporting a different model.

If part of your team isn't familiar with data assimilation yet, you should
review the introductory material in this documentation and and also look at
work through the concepts in the :doc:`/theory/readme`.

Goals of using DART
-------------------

DART is the Data Assimilation Research Testbed.  It is a collection of 
tools, routines, and scripts that allow users to build custom solutions
and explore a variety of DA related efforts.  It is not a turnkey system;
it must be built before use and is often customized based on needs and goals.

DART is often used for the following types of projects:

- Learning about Data Assimilation (DA)
- Using DART with an existing model and supported observations
- Using DART with a new model: :ref:`Porting new models`
- Using new observations with DART in an existing model
- Using both a new model and new observations with DART
- Using DART to teach DA

You can view a list of models that are already supported at :ref:`Supported models`
and a list of supported observations at :ref:`programs`.

Everything on this "possible goals" list except adding support for a new model
can generally be done by a single user with minimal help from the DART team.
Therefore this discussion focuses only on adding a new model to DART.

Should I consider using DART with my model?
-------------------------------------------

DART is an ensemble-based DA system. It makes multiple runs of a model with
slightly different inputs and uses the statistical distribution of the results
to decide how to adjust the model state to be more consistent with the
observations.

The advantage of ensemble systems is that no changes to the model itself are
required. The disadvantage is that multiple runs of the model are needed and
this can be computationally expensive.

Simple models can be added to DART with a single person effort, but
larger, more complex models can require multiple person-months with
support from the DART team to add the interfaces and scripts needed 
to perform a large-scale DA experiment.

The DART code is in Fortran. The supporting scripts and tools are
a mix of shell scripts and python. The model can be written in any language;
it will only be run and the input and output files will be used by DART.

Things to discuss before beginning
----------------------------------

Is your model appropriate for any kind of DA?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your model isn't chaotic, you don't need data assimilation.
In non-chaotic models, you can improve your predictions by running the model, 
examining the difference between the prediction and the observations, inverting
the equations inside the model to compute how different inputs would have
produced outputs closer to the observations.

Chaotic models do not have a simple relationship between inputs and
outputs. There are internal feedbacks and non-linear behaviors that make
it difficult to adjust the inputs to make the outputs better match the
observations.  

What is your model state?
~~~~~~~~~~~~~~~~~~~~~~~~~

"Model state" has a very specific definition that can be the source
of much confusion if someone running a model has not thought about
DA before.  Formally it is the minimal set of variables that must be 
saved when a model stops so it can be restarted again exactly.

At first glance this means all the variables on the right side of
the equals sign for the governing equations of the system.  However
many models which have not been designed with DA in mind may have
no clear time when all parts of the model are at a consistent time.
e.g. some variables may be 1/2 timestep ahead or behind others.
Some derived variables may be expensive to compute and so are
precomputed and stored and not recomputed.  If the DA process changes
the state variables all derived variables must be recomputed before
proceeding.

Restart files often store many more variables than the minimal set
needed to restart the model.  Often other variables are used in 
diagnostic routines or are of interest on their own.  Generally
these are not considered part of the model state.

How is your model execution controlled?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generally larger and more complex models have an environment they
are expecting to run in/with.  e.g. scripts to control the execution
parameters, or input parameter files; how many processors are used in
a parallel system, how the tasks are distributed over the hardware;
how long does the execution run, in model time, and what variables are
written to the output files.

For DA, at a minimum there must be a way to control how long the model 
runs before it writes out the results and exits.  

For large models, the DA filter process is a large parallel program
generally requiring a multi-processor supercomputer or cluster.  Many
models themselves are large parallel programs, so there can be issues
with how the switch between model and DA process is done.

New or adjusted scripting is generally required to include the DA process
in the overall execution flow.

Cycling with a DA system
~~~~~~~~~~~~~~~~~~~~~~~~

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
may not result in identical values. Large changes suggest that the
model is not a good candidate for a cycling DA system.

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

For models which have a lot of error growth it may be enough to
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


What code is required to interface a model with DART?
-----------------------------------------------------

There is a single FORTRAN module that hides the model details from the
rest of the DART system.  Generally the routines which require the most
work are the interpolation routine, followed by the metadata routine
and the "get close" localization routines.

Interpolation
~~~~~~~~~~~~~

Given an observation quantity and location, the model interface routines
must return an array of values, one for each ensemble member.  The values
must be the best estimate of what a real instrument would return if the
real state of the system were each of the ensemble values.  

For a regular grid this can be computed fairly simply with routines
already provided in the DART system.  It involves locating the grid
values that enclose the observation location, and doing bi- or tri-linear
interpolation to the actual location.

However, many models have non-regular grid, especially in the vertical
coordinates for an earth-system-based model.  Or the grid can be an 
irregular mesh or deformed mesh.  It may take searching or transforms
to identify the closest values in the model state to use for interpolation.

Metadata
~~~~~~~~

Given an offset into the model state, the model interface routines
must return the location in the selected coordinate system, and the 
quantity at that offset.

There are routines provided which simplify this for regular or deformed
grids, so this generally is not too complex but may require additional
arrays for irregular grids or unstructured grids.

Localization
~~~~~~~~~~~~

DART bases the impact of observations on the model state on the
correlation between the array of predicted observation values, the
actual observation value and error, and the array of model state values.

In practice observations are only correlated with model state values
"close" to the observation.  Spurrious correlations can occur which
degrade the results after assimilation.  Also there are efficiency gains
if only parts of the model state which are "close" to the observation
are processed.  

DART includes routines which can compute what part of the state are
close to a given observation.  However some models have special considerations
for whether they want to control the impact of observations on parts
of the model state and this can be adjusted based on code added to the
model-specific parts of getting close observations and model state.

Vertical issues
~~~~~~~~~~~~~~~

Most Earth System models use Latitude and Longitude for horizontal
coordinates or can generate them if needed (e.g. spectral models
can transform their state into Lat/Lon coords).  But often vertical
coordinates pose additional complications.

If the model and the observations both use the same coordinates for
vertical, e.g. pressure or height, then there are no need for
conversion routines.  But some models use terrain-following
coordinates, or a mix of pressure and terrain coordinates.
Observation vertical locations can be reported in height or in
pressure.

Additionally, if vertical localization is to be done in a different
coordinate than the model or observations (e.g. scale height), then
conversion routines are needed.

The interface code may need to read in additional arrays from the
model in order to convert the vertical coordinates accurately.

During the run of filter there are two options for when vertical
conversion is done: all at the start, or on demand.  If the observations
to be assimilated are expected to impact all or almost all of the
state, doing all vertical conversion at the start is more efficient.
If the observations are expected to impact only a small percentage
of the state variables then doing it on demand is more efficient.
The options here are namelist selectable at runtime and the impact
on total runtime can be easily measured and compared.


