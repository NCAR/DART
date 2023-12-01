DART Observations
=================

Overview
--------

Real-world observations of earth-system data come from a variety of
sources, including radiosondes, satellites, ships, aircraft, weather
stations, etc. The files in this *observations* directory can be used to
convert data from a variety of native formats into a common DART
*observation sequence* format.

Synthetic observations are those not based on an actual instrument
reading of a system, but instead are fabricated to have a known value,
or have values computed by running a model, possibly with a fixed amount
of simulated noise added. These observations can be used for testing,
determining the sensitivity of the model to assimilation, and for
designing new observation systems. The DART system includes several ways
to create synthetic observations. See the :ref:`programs <programs>`
section below for more details.

The DART framework enforces a clean separation between observations and
the models they are assimilated into. The same observations can be used
in any model which understands how to generate a value for the requested
type of observation from its state space values.

In many cases a single, self-contained program can convert directly from
the observation location, time, value, and error into the DART format.
In other cases, especially those linking with a complicated external
library (e.g. BUFR), there is a two-step process with two programs and
an ASCII intermediate file. We are currently leaning towards single-step
conversions but either approach can be used for new programs.

Frequently the original datasets are in a standard scientific format
like netCDF, HDF, or BUFR, and library routines for those formats can be
used to read in the original observation data.

The DART software distribution includes Fortran subroutines and
functions to help create a sequence of observations in memory, and then
a call to the DART observation sequence write routine will create an
entire *obs_seq* file in the correct format.

The DART system comes with several types of location modules for
computing distances appropriately. Two of the ones most commonly used
are for data in a 1D system and for data in a 3D spherical coordinate
system. Most of the programs here assume the
*location/threed_sphere/location_mod.f90* 3D sphere location module is
being used.

There are currently some additional observation sources and types which
we are in the process of collecting information and conversion programs
for and which will eventually be added to this directory. In the
meantime, if you have converters for data or interest in something that
is not in the repository, please `email the DART
group <mailto:dart@ucar.edu>`__.


Data Sources and Formats
------------------------

See the various subdirectories here, which generally include information
on where the example data was obtained and in what format it is
distributed. Most data is available for download off the web. The Data
Support Section (DSS) at NSF NCAR has large data repositories, the MADIS
data center distributes observations in NetCDF format, GTS real-time
weather data is available from various sources. For new converters, if
you can find what format the data is distributed in you may be able to
adapt one of the existing converters here for your own use. Formats read
by the existing converters include NetCDF, HDF, little-r, text,
Prepbufr, amongst others.

**See the**  :ref:`programs <programs>` **section below for a list of the
current converter programs. It might save you from reinventing the
wheel.**

If you have looked and none of the existing converters are right for
your data, here are some suggestions for where to start creating a new
converter. Create a new subdirectory in the *observations* directory.
Copy with the recursive option (*cp -r*) one of the existing converters
and adapt to your needs. Our suggestions for which converter to start
from depends on the format of your input observations to be converted.
If your input data format is:

Start with the MADIS converters, and in particular try the
convert_madis_profiler.f90 file because it is the most straightforward.
Another good option is SST/oi_sst_to_obs.f90


+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| netCDF                            | Start with the MADIS converters, and in particular try the convert_madis_profiler.f90 file because it is the most straightforward. Another good option is SST/oi_sst_to_obs.f90 |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Comma separated text              | Start with the Ameriflux converter.                                                                                                                                             |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Generic text                      | Start with the text converter.                                                                                                                                                  |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| HDF-EOS                           | Start with the AIRS converter.                                                                                                                                                  |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| BUFR or prepBUFR                  | Start with the NCEP converter.                                                                                                                                                  |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Dense data, like Satellite swaths | Start with the tpw converter, which includes code that averages the raw data in space and time.                                                                                 |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Ray-path integrated data          | Start with the GPS converter, which includes code that traces a path and integrates values along the ray.                                                                       |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| World Ocean Database packed ASCII | Start with the WOD converter.                                                                                                                                                   |
+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Decisions You May Need to Make
------------------------------

Time
~~~~

Time enters into the assimilation system in 3 places: the time of
the state vector data (the current model time when this data was
produced), the time of each observation, and the assimilation window
length. The window length is set by the model-dependent routine
``shortest_time_between_assimilations()``.
The internal timestepping of the model is unrelated to any of these times
and is outside the scope of the assimilation system.

The basic time type in DART is a pair of integers; one for the day
number and one for the number of seconds. Generally the low order
models, which aren’t direct geophysical models, use time directly as a
sequence of days starting at 0 and incrementing in any appropriate
number of seconds or days. The observations assimilated into these
systems do not need to use a calendar.

Observations of a real-world system usually are distributed with a
year/month/day, hour/min/seconds timestamp. There are routines in DART
to convert back and forth between the (day-number/seconds) format and a
variety of (year/month/day) calendars. See `the time manager
documentation <../../assimilation_code/modules/utilities/time_manager_mod.html#time_type>`__
for more details on how DART stores time information and the types of
available calendars. Some climate models which do long runs (100s or
1000s of years) use a modified calendar for simplicity in computation,
e.g. months which always have 30 days, or no leap years. When trying to
assimilate real observations into these models there may be calendar
issues to solve.

The smallest resolvable unit of time in DART is a second. To model a
system which operates on sub-second time scales the time can be scaled
up by some factor. As long as the observation time, the state data time,
and the minimum model advance time are expressed in the same scaled time
units, there is no problem.

Error
~~~~~

Observations must specify an associated expected error. Each individual
observation stores its own error value, so it can be a constant value
for all observations of that type or it can vary by location, by height,
by magnitude of the observed value, etc. This value is the expected
instrument error plus the representativeness error of the model. The
model error includes deficiencies in the equations representing the
processes of the system as well as errors introduced by representing a
continuous system as a series of discrete points. While the instrument
error and the representativeness error could be specified separately,
they each have the same impact on the assimilation and can be difficult
to determine with any real accuracy. For simplicity, in DART (and most
current assimilation software) they are combined and specified as a
single value.

The instrument error is generally supplied by the instrument maker.
Sadly, it is frequently surprisingly difficult to find these values. For
the representativeness error, a set of artificial observations could be
generated with the
`perfect_model_obs <../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html>`__
program and an assimilation experiment could be run to generate an
estimate of the error in the model. In practice however most people make
an educated guess on the values of the error and then start with a
larger than expected value and decrease it based on the results of
running some test assimilations. For these tests the namelist for the
`outlier
threshold <../../assimilation_code/programs/filter/filter.html#Namelist>`__
should be disabled by setting it to -1 (the default value is 3). This
value controls whether the observation is rejected because the observed
value is too far from the ensemble mean.

If the diagnostics show that the difference between the mean of the
forward operators and the observed value is consistently smaller than
the specified observation error, then the error is probably too large. A
too-large error reduces the impact of an observation on the state. If
the specified observation error is too small it is likely the
observation will be rejected when the outlier threshold is enabled, and
the observation will not be assimilated. It is important to look at the
output observation sequence files after an assimilation to see how many
observations were assimilated or rejected, and also at the RMSE (`root
mean squared error <http://www.wikipedia.org/wiki/RMSE>`__) versus the
total spread. DART includes Matlab diagnostic routines to create these
types of plots. The observation RMSE and total spread should be roughly
commensurate. The total spread includes contributions from both the
ensemble variance and the observational error variance, so it can be
adjusted by changing the error values on the incoming observations.
There are other ways to adjust the ensemble spread, including
`inflation <../../assimilation_code/programs/filter/filter.html#Inflation>`__,
so the observation error is not the only factor to consider.

One last recommendation: if possible, the Prior forward operator values
should be compared against the observations after several assimilation
cycles. If you plot results using the Posterior values it is always
possible for the assimilation to overfit the observations and look good
on the diagnostic plots. But the actual test is to then advance the
model and look at how the forecast of the state compares to the
observations.

Types
~~~~~

All observations have to have a specific ‘type’. There are namelist
controls to turn on and off the assimilation of observations at run-time
by type, or to only evaluate the forward operator for an observation but
have no impact on the state. Several of the diagnostics also group
observations by type to give aggregate statistics after an assimilation.
Generally types are based on both the observing platform or instrument
as well as the kind of observation, e.g. RADIOSONDE_TEMPERATURE,
ARGO_SALINITY, etc. Each type is associated with a single underlying
generic ‘kind’, which controls what forward operator code is called
inside the model, e.g. QTY_TEMPERATURE, QTY_DENSITY, etc.

See `here <../forward_operators/obs_def_mod.html>`__ for more details on
how to use and add new DART types. The DART obs_kind_mod.f90 defines a
list of already defined observation kinds, and users can either use
existing observation types in ‘obs_def_xxx_mod.f90’ files, or define
their own.

Locations
~~~~~~~~~

The two most common choices for specifying the location of an
observation are the
`threed_sphere <../../assimilation_code/location/threed_sphere/location_mod.html>`__
and the
`oned <../../assimilation_code/location/oned/location_mod.html>`__
locations. For observations of a real-world system, the 3D Sphere is
generally the best choice. For low-order, 1D models, the 1D locations
are the most commonly used. The observation locations need to match the
type of locations used in the model.



Converting a series of observations
-----------------------------------

If you are running a series of assimilation steps you may
need a separate observation sequence (obs_seq) file per step.
The suggested process is to create the first few files by hand to check
the resulting obs_seq files and then write scripts (python, shell)
to automate the creation of the remainder of the files.
The following are some of the considerations to take
into account when creating scripts for a series of obs_seq files.

Looping in Time
~~~~~~~~~~~~~~~

Often observations are distributed in files that contain observations
from a particular time period, e.g. a file per day or per week.
The output obs_seq files need to include observations from the same
time period as the assimilation window; how often the assimilation
is stopped and the model is advanced in time.  The conversion process
can either convert all the observations from an input file into a single
output file and in a subsequent step break the file into the required
time ranges, or the conversion process can extract and convert only
the observations required for a single output file and loop multiple
times over the same input file.

Generally earth system models use calendar dates, including months,
days, years, hours, minutes and seconds.
The ``advance_time`` program is very useful in adding or subtracting time periods
from calendar dates taking into account changing months and years,
accounting for leap days, etc.

Observation conversion programs usually take one of two strategies
for their input and output filenames.

* Have fixed input and output filenames for the converter.
  Have the script make symbolic links from the actual filenames to the
  fixed names for the files for each conversion run.

* Have a Fortran namelist variable that sets the input and output
  filenames for the converter.  Have the script generate or edit the
  namelist file (e.g. with the `sed` stream editor) to set the actual
  filenames for each conversion run.

Generally it is a good idea to encode the date information in the
output filename so each file is guarenteed to be unique.
This can also make it simpler at filter runtime to generate the
required input observation sequence filenames using a program
like ``advance_time``.


Multiple Observation Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

It is common that an assimilation will want to use observations
from different sources.  Generally it is easier to convert observations
from each source separately and then merge them together with the
``obs_sequence_tool``.

Creating filenames and directory names which follow a pattern
that can be generated with the ``advance_time`` program makes this easier to do.

The ``obs_sequence_tool`` can read the input filenames from a separate ascii file.
This makes generating the filenames easy from a script; it can
simply concatinate the input filenames echo'd to an ascii file and
then run the obs_sequence_tool.  The output file can either be set
by using ``sed`` on the namelist, or a fixed output filename can be used
and then the file renamed after the tool has run.


Conversion Run Time for Large File Counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If 100s of files need to be generated and a supercomputer or other multiple-CPU
resource is available, batch files which convert multiple files at the same time
can be a large time savings.  Care must be taken that each conversion has its own
settings and unique filenames. Often a separate working directory from other
conversions running at the same time simplifies the scripting needed.


Verification
~~~~~~~~~~~~

Observations taken from real-world sources can have missing values, illegal
values, missing files, duplicated data, etc.  The list is as long as your
imagination.  It can be very useful to write or adapt programs like ``obs_info``
to print out the first and last obs times in a file, the count of each obs type, etc.
Especially for observations which are close to the start/end of a month or year,
it is easy to find truncated data files.

If converting a large number of files it is also common for computer system
failures to occur at random times.  File systems fill up, batch jobs exit early,
power glitches stop programs before they finish.  Look for anomolous observation
counts, unexpected first and last times of obs in a file, missing files, files
with many fewer bytes than others, and anything else you can think of.


Output Formats
~~~~~~~~~~~~~~

There are options to write output obs_seq files in binary, which are roughly
half the size of ascii files.  However it greatly increases the effort to
examine the contents of a file for problems.  Generally we have used the ascii
format. It is portable between systems of different "endians" (order of bytes
in a multi-byte number) and can be browsed much more easily.


.. _programs:

Converter programs
==================

The *DART/observations/obs_converters* directory contains a variety of
converter programs to read various external formats and convert the
observations into the format required by DART.

The current list of converters (some directories contain multiple
converters) include:

-  ``AIRS``: :doc:`./AIRS/README`
-  ``AURA``: See ``./AURA``
-  ``Aviso+/CMEMS``: :doc:`./AVISO/AVISO`
-  ``Ameriflux``: :doc:`./Ameriflux/level4_to_obs`
-  ``CHAMP``: :doc:`./CHAMP/work/README`
-  ``cice``: :doc:`./cice/cice_to_obs`
-  ``CNOFS``: See ``./CNOFS``
-  ``CONAGUA``: :doc:`./CONAGUA/README`
-  ``COSMOS``: :doc:`./COSMOS/COSMOS_to_obs`
-  ``DWL``: :doc:`./DWL/dwl_to_obs`
-  ``GMI``: :doc:`./GMI/README`
-  ``GOES``: :doc:`./GOES/README`
-  ``GPSPW``: :doc:`./GPSPW/README`
-  ``GRACE``: See ``./GRACE``
-  ``GSI2DART``: :doc:`./GSI2DART/readme`
-  ``GTSPP``: :doc:`./GTSPP/GTSPP`
-  ``MADIS``: :doc:`./MADIS/MADIS`
-  ``MIDAS``: :doc:`./MIDAS/MIDAS_to_obs`
-  ``MODIS``: :doc:`./MODIS/MOD15A2_to_obs`
-  ``MPD``: See ``./MPD``
-  ``NASA_Earthdata``:doc:`./NASA_Earthdata/README`
-  ``NCEP``: (prepbufr-> ascii) :doc:`./NCEP/prep_bufr/prep_bufr`
-  ``NCEP``: (ascii-> obs_seq) :doc:`./NCEP/ascii_to_obs/create_real_obs`
-  ``NSIDC``:doc:`./NSIDC/SMAP_L2_to_obs`
-  ``ROMS``: :doc:`./ROMS/ROMS`
-  ``SIF``: :doc:`./SIF/SIF_to_obs_netcdf`
-  ``SSEC``: :doc:`./SSEC/SSEC`
-  ``SST``: :doc:`./SST/SST`
-  ``ocean color``: :doc:`./ocean_color/README`
-  ``SSUSI``: :doc:`./SSUSI/convert_f16_edr_dsk`
-  ``WOD``: :doc:`./WOD/WOD`
-  ``gnd_gps_vtec``: :doc:`./gnd_gps_vtec/README`
-  ``GPS``: :doc:`./gps/gps`
-  ``ok_mesonet``: :doc:`./ok_mesonet/ok_mesonet`
-  ``QuikSCAT``: :doc:`./quikscat/QuikSCAT`
-  ``Radar``: :doc:`./radar/README`
-  ``snow``: :doc:`./snow/snow_to_obs`
-  ``Text``: :doc:`./text/text_to_obs`
-  ``text_GITM``: See ``./text_GITM``
-  ``tpw``: :doc:`./tpw/tpw`
-  ``Tropical Cyclones``: :doc:`./tropical_cyclone/tc_to_obs`
-  ``3DVAR/4DVAR``: :doc:`./var/var`
-  ``Var (little-r)``: :doc:`./var/littler_tf_dart`
-  ``Var (radar)``: :doc:`./var/rad_3dvar_to_dart`

In addition the following external program produces DART observation
sequence files:

-  `Observation Processing And Wind Synthesis
   (OPAWS) <http://code.google.com/p/opaws/>`__: OPAWS can process NSF NCAR
   Dorade (sweep) and NSF NCAR EOL Foray (netcdf) radar data. It analyzes
   (grids) data in either two-dimensions (on the conical surface of each
   sweep) or three-dimensions (Cartesian). Analyses are output in
   netcdf, Vis5d, and/or DART (Data Assimilation Research Testbed)
   formats.

For generating synthetic observations, see the
`create_obs_sequence <../../assimilation_code/programs/create_obs_sequence/create_obs_sequence.html>`__
program documentation. You can also generate observation files based on
text input. See the `text_to_obs <text/text_to_obs.html>`__ program
documentation and even_sphere. Or for simulating a large complex observing system, you
can use the DART library routines in a Fortran program to compute the
observation information and have the DART routines write the output
file.

There are a couple utilities of note:

-  `even_sphere <even_sphere/README.html>`__ - a utility for generating
   a text file of evenly-spaced observation locations that can then be used in a
   perfect model experiment.
-  `obs_error <obs_error/README.html>`__ - modules that specify observation
   errors based on what is used by ECMWF and NCEP


See the
`perfect_model <../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html>`__
program documentation on how to run a model with a set of observations
that have only locations, types, and times, and have the forward
operators compute the observation values.

Contact the `DART development group <mailto:dart@ucar.edu>`__ if you
have observations in a different format that you want to convert. We can
give you advice and pointers on how to approach writing the code.
