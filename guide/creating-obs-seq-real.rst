Creating an obs_seq file from real observations
===============================================

Real observations come in a mind-boggling diversity of formats. We have
converters for many formats in the ``DART/observations/obs_converters``
directory. The documentation for that directory is listed in
:doc:`../observations/obs_converters/README`.

The converters are designed to work on one input file format and create (or add
to) an output observation sequence. It may be desirable to post-process multiple
observation sequence files with the
:doc:`../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool` to
select for timeframe, geographic region, etc.

Many of the formats require their own libraries (like HDF), and require intimate
knowledge of the data format to extract the portions required for the :doc:`DART
observation sequence file <detailed-structure-obs-seq>`.

You should feel free to browse the converters and their companion
documentation. If you create a new observation coverter for a format that DART
doesn't already support, please follow the :doc:`contributors-guide` to add 
your code to DART. These types of contributions are greatly appreciated by 
DAReS staff and by the geoscience community!

The DART framework enforces a clean separation between observations and the
models used for assimilation. The same observations can be used in any model
which understands how to generate a value for the requested type of observation
from the models’ state-space values (i.e. the forward observation operator must
exist - DART provides many for the most common state variables).

In many cases, the original datasets are in a standard scientific format like
netCDF, HDF, or BUFR, and library routines for those formats can be used to read
in the original observation data. The DART software distribution includes
Fortran subroutines and functions to help create a sequence of observations in
memory, and then a call to the DART observation sequence write routine will
create an entire *obs_seq* file in the correct format.

In many cases, a single, self-contained program can convert directly from the
observation location, time, value, and error into the DART format. In other
cases, especially those linking with a complicated external library (e.g. BUFR),
there is a two-step process with two programs and an ASCII intermediate file. We
are currently leaning towards single-step conversions but either approach can be
used for new programs.

The DART system comes with several types of location modules for computing
distances appropriately. The two most commonly used are for data in a 1D system
and for data in a 3D spherical coordinate system. All the programs in the
``DART/observations`` directory assume the
``assimilation_code/location/threed_sphere/location_mod.f90`` 3D sphere location
module is being used.

With the myriad of observation file formats, HDF, Grib, BUFR, netCDF, … we
simply have not had the time nor need to support all of them. The converters
are a work in progress. There are currently about 10 other observation sources
and types which we are in the process of collecting information and conversion
programs for and which will eventually be added to this directory. In the
meantime, if you have converters for data or interest in something that is not
in the repository, please email the DART group. Your best bet is to contact our
group at *dart@ucar.edu* with a specific request and we can steer you to the
most similar process.

Overview
--------

Real-world observations of earth-system data come from a variety of sources,
including radiosondes, satellites, ships, aircraft, weather stations, etc. The
files in this ``observations`` directory can be used to convert data from a
variety of native formats into a common DART *observation sequence* format.

Synthetic observations are those not based on an actual instrument reading of a
system, but instead are fabricated to have a known value, or have values
computed by running a model, possibly with a fixed amount of simulated noise
added. These observations can be used for testing, determining the sensitivity
of the model to assimilation, and for designing new observation systems. The
DART system includes several ways to create synthetic observations. For more 
information, see :doc:`creating-obs-seq-synthetic`.

The DART framework enforces a clean separation between observations and the
models they are assimilated into. The same observations can be used in any model
which understands how to generate a value for the requested type of observation
from its state space values.

In many cases a single, self-contained program can convert directly from the
observation location, time, value, and error into the DART format. In other
cases, especially those linking with a complicated external library (e.g. BUFR),
there is a two-step process with two programs and an ASCII intermediate file. We
are currently leaning towards single-step conversions but either approach can be
used for new programs.

Frequently the original datasets are in a standard scientific format like
netCDF, HDF, or BUFR, and library routines for those formats can be used to read
in the original observation data.

The DART software distribution includes Fortran subroutines and functions to
help create a sequence of observations in memory, and then a call to the DART
observation sequence write routine will create an entire *obs_seq* file in the
correct format.

The DART system comes with several types of location modules for computing
distances appropriately. Two of the ones most commonly used are for data in a 1D
system and for data in a 3D spherical coordinate system. All the programs here
assume the ``location/threed_sphere/location_mod.f90`` 3D sphere location module
is being used.

There are currently some additional observation sources and types which we are
in the process of collecting information and conversion programs for and which
will eventually be added to this directory. In the meantime, if you have
converters for data or interest in something that is not in the repository,
please contact DAReS staff by emailing dart@ucar.edu.

Data sources and formats
------------------------

See the various subdirectories here, which generally include information on
where the example data was obtained and in what format it is distributed. Most
data is available for download off the web. The Data Support Section (DSS) at
NSF NCAR has large data repositories, the MADIS data center distributes observations
in netCDF format, GTS real-time weather data is available from various sources.
For new converters, if you can find what format the data is distributed in you
may be able to adapt one of the existing converters here for your own use.
Formats read by the existing converters include netCDF, HDF, little-r, text,
Prepbufr, amongst others.

See the current list of :doc:`converter programs <available-observation-converters>`

If you have looked and none of the existing converters are right for your data,
here are some suggestions for where to start creating a new converter. Create a
new subdirectory in the *observations* directory. Copy with the recursive option
(*cp -r*) one of the existing converters and adapt to your needs. Our
suggestions for which converter to start from depends on the format of your
input observations to be converted. If your input data format is:

+---------------------------------------+---------------------------------------+
| format                                | advice                                |
+=======================================+=======================================+
| netCDF                                | Start with the *MADIS* converters,    |
|                                       | and in particular try the             |
|                                       | ``convert_madis_profiler.f90`` file   |
|                                       | because it is the most                |
|                                       | straightforward. Another good option  |
|                                       | is ``SST/oi_sst_to_obs.f90``.         |
+---------------------------------------+---------------------------------------+
| Comma separated text                  | Start with the *Ameriflux* converter. |
+---------------------------------------+---------------------------------------+
| Generic text                          | Start with the *text* converter.      |
+---------------------------------------+---------------------------------------+
| HDF-EOS5                              | Start with the *AIRS* converter.      |
+---------------------------------------+---------------------------------------+
| BUFR or prepBUFR                      | Start with the *NCEP* converter.      |
+---------------------------------------+---------------------------------------+
| Dense data, like Satellite swaths     | Start with the *tpw* converter, which |
|                                       | includes code that averages the raw   |
|                                       | data in space and time.               |
+---------------------------------------+---------------------------------------+
| Ray-path integrated data              | Start with the *GPS* converter, which |
|                                       | includes code that traces a path and  |
|                                       | integrates values along the ray.      |
+---------------------------------------+---------------------------------------+
| World Ocean Database packed ASCII     | Start with the *WOD* converter.       |
+---------------------------------------+---------------------------------------+

.. raw:: html

   <!--
   The existing DART csv readers are:
   vi -R Ameriflux/level4_to_obs.f90 \
   CHAMP/CHAMP_density_text_to_obs.f90 \
   CNOFS/CNOFS_text_to_obs.f90 \
   COSMOS/COSMOS_development.f90 \
   COSMOS/COSMOS_to_obs.f90 \
   MODIS/MOD15A2_to_obs.f90 \
   ROMS/convert_roms_obs.f90 \
   gnd_gps_vtec/gnd_gps_vtec_text_to_obs.f90 \
   gps/convert_cosmic_gps_cdf.f90 \
   gps/convert_cosmic_ionosphere.f90 \
   quikscat/quikscat_JPL_mod.f90 \
   snow/snow_to_obs.f90 \
   text/text_to_obs.f90 \
   text_GITM/text_to_obs.f90   -->

Decisions you might need to make
--------------------------------

Time
~~~~

Time enters into the assimilation system in 3 places: the timestamp of the state
vector data (the current model time when this data was produced), the time of
each observation, and the minimum time period the model should be called to
advance (the assimilation window size). The internal timestepping of the model
is unrelated to any of these times and is outside the scope of the assimilation
system.

The basic time type in DART is a pair of integers; one for the day number and
one for the number of seconds. Generally the low order models, which aren’t
direct geophysical models, use time directly as a sequence of days starting at 0
and incrementing in any appropriate number of seconds or days. The observations
assimilated into these systems do not need to use a calendar.

Observations of a real-world system usually are distributed with a
year/month/day, hour/min/seconds timestamp. There are routines in DART to
convert back and forth between the (day-number/seconds) format and a variety of
(year/month/day) calendars. For more details on how DART stores time
information and the types of available calendars, see
:doc:`../assimilation_code/modules/utilities/time_manager_mod`.

Some climate models which do long runs (100s or 1000s of years) use a modified
calendar for simplicity in computation, e.g. months which always have 30 days,
or no leap years. When trying to assimilate real observations into these models
there may be calendar issues to solve.

The smallest resolvable unit of time in DART is a second. To model a system
which operates on sub-second time scales the time can be scaled up by some
factor. As long as the observation time, the state data time, and the minimum
model advance time are expressed in the same scaled time units, there is no
problem.

Error variances
~~~~~~~~~~~~~~~

Observations must specify an associated expected error variance. Each individual
observation stores its own error variance value, so it can be a constant value
for all observations of that type or it can vary by location, by height, by
magnitude of the observed value, etc. This value is the expected instrument
error variance plus the representativeness error variance of the model. The
model error variance includes deficiencies in the equations representing the
processes of the system as well as errors introduced by representing a
continuous system as a series of discrete points. While the instrument error and
the representativeness error could be specified separately, they each have the
same impact on the assimilation and can be difficult to determine with any real
accuracy. For simplicity, in DART (and most current assimilation software) they
are combined and specified as a single value, which we frequently call the
‘observation error’. Keep in mind we really mean ‘observation error variance’.

The instrument error is generally supplied by the instrument maker. Sadly, it is
frequently surprisingly difficult to find these values. For the
representativeness error, you can generate a set of artificial observations
with the
:doc:`../assimilation_code/programs/perfect_model_obs/perfect_model_obs`
and then run an assimilation experiment to generate an estimate of the error in
the model.

In practice, however, most people make an educated guess on the values of the
error and then start with a larger than expected value and decrease it based on
the results of running some test assimilations.

For these tests, the namelist for the outlier threshold in the ``filter_nml``
namelist of ``input.nml`` should be disabled by setting it to -1 (the default
value is 3). This value controls whether the observation is rejected because
the observed value is too far from the ensemble mean.

If the diagnostics show that the difference between the mean of the forward
operators and the observed value is consistently smaller than the specified
observation error, then the error is probably too large. A error that is too
large reduces the impact of an observation on the state. If the specified
observation error is too small it is likely the observation will be rejected
when the outlier threshold is enabled, and the observation will not be
assimilated. It is important to look at the output observation sequence files
after an assimilation to see how many observations were assimilated or
rejected, and also at the RMSE (`root mean squared
error <http://www.wikipedia.org/wiki/RMSE>`__) versus the total spread. DART
includes Matlab diagnostic routines to create these types of plots. The
observation RMSE and total spread should be roughly commensurate. The total
spread includes contributions from both the ensemble variance and the
observational error variance, so it can be adjusted by changing the error
values on the incoming observations.

There are other ways to adjust the ensemble spread, including :doc:`inflation`,
so the observation error is not the only factor to consider.

One last recommendation: if possible, the Prior forward operator values should
be compared against the observations after several assimilation cycles. If you
plot results using the Posterior values it is always possible for the
assimilation to overfit the observations and look good on the diagnostic plots.
But the actual test is to then advance the model and look at how the forecast of
the state compares to the observations.

Observation types
~~~~~~~~~~~~~~~~~

All observations have to have a specific ‘type’. There are namelist controls to
turn on and off the assimilation of observations at run-time by type, or to only
evaluate the forward operator for an observation but have no impact on the
state. Several of the diagnostics also group observations by type to give
aggregate statistics after an assimilation. Generally types are based on both
the observing platform or instrument as well as the ‘kind’ of observation,
e.g. RADIOSONDE_TEMPERATURE, ARGO_SALINITY, etc. Each type is associated with a
single underlying generic ‘kind’, which controls what forward operator code is
called inside the model, e.g. QTY_TEMPERATURE, QTY_DENSITY, etc.

For more details on how to use and add new DART types, see the 
:doc:`../observations/forward_operators/obs_def_mod`.

The DART ``obs_kind_mod.f90`` defines a list of already defined observation
types, and users can either use existing observation types in
‘obs_def_xxx_mod.f90’ files, or define their own. Be aware that
``obs_kind_mod.f90`` is autogenerated by the
:doc:`/assimilation_code/programs/preprocess/preprocess`, so until you
configure and run ``preprocess``, ``obs_kind_mod.f90`` will not exist.

Observation locations
~~~~~~~~~~~~~~~~~~~~~

The two most common choices for specifying the location of an observation are
the :doc:`../assimilation_code/location/threed_sphere/location_mod` and the
:doc:`../assimilation_code/location/oned/location_mod` locations.

For observations of a real-world system, the 3D Sphere is generally the best
choice. For low-order, 1D models, the 1D locations are the most commonly used.
The observation locations need to match the type of locations used in the model
in that you cannot read observations on a unit circle (1D) when using models
that require 3D Sphere locations.

The choice of the vertical coordinate system may also be important. For the 3D
Sphere, the vertical coordinate system choices are:

================= ============= ===============================================
string            integer value meaning
================= ============= ===============================================
VERTISUNDEF       -2            has no specific vertical location (undefined)
VERTISSURFACE     -1            surface value (value is surface elevation in m)
VERTISLEVEL       1             by model level
VERTISPRESSURE    2             by pressure (in pascals)
VERTISHEIGHT      3             by height (in meters)
VERTISSCALEHEIGHT 4             by scale height (unitless)
================= ============= ===============================================

The choice of the vertical coordinate system may have ramifications for vertical
localization, depending on your model’s ability to convert from one coordinate
system to another. ``VERTISUNDEF`` is typically used for column-integrated
quantities. ``VERTISLEVEL`` only makes sense for synthetic observations.

When observations are declared to be ``VERTISSURFACE`` or ``VERTISUNDEF``
it is not possible to compute a vertical distance between the observation and
anything else. Consequently, the distance between that observation and everything
else (state, other observations) is strictly a horizontal distance, and the observation
will impact the entire column (all levels) within the horizontal localization radius.

