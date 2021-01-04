
<span id="TOP" class="anchor"></span>

DART Observations
=================

---
![DARTlogo](https://github.com/NCAR/DART/blob/Manhattan/docs/images/Dartboard7.png)
---

[OVERVIEW](#Overview) / [DATA SOURCES](#datasourcesandformats) /
[DECISIONS](#Decisions) / [PROGRAMS](#Programs) / [KNOWN BUGS](#knownbugs) / 
[FUTURE PLANS](#futureplans) / [TERMS OF USE](#termsofuse)

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
to create synthetic observations. See the [Programs](#Programs) section
below for more details.

The DART framework enforces a clean separation between observations and
the models they are assimilated into. The same observations can be used
in any model which understands how to generate a value for the requested
type of observation from its state space values.

In many cases a single, self-contained program can convert directly from
the observation location, time, value, and error into the DART format.
In other cases, especially those linking with a complicated external
library (e.g. BUFR), there is a two-step process with two programs and
an ASCII intermediate file. We are currently leaning towards single-step
conversions but either approach can be used for new programs.

Frequently the original datasets are in a standard scientific format
like netCDF, HDF, or BUFR, and library routines for those formats can be
used to read in the original observation data.

The DART software distribution includes Fortran subroutines and
functions to help create a sequence of observations in memory, and then
a call to the DART observation sequence write routine will create an
entire *obs\_seq* file in the correct format.

The DART system comes with several types of location modules for
computing distances appropriately. Two of the ones most commonly used
are for data in a 1D system and for data in a 3D spherical coordinate
system. Most of the programs here assume the
*location/threed\_sphere/location\_mod.f90* 3D sphere location module is
being used.

There are currently some additional observation sources and types which
we are in the process of collecting information and conversion programs
for and which will eventually be added to this directory. In the
meantime, if you have converters for data or interest in something that
is not in the repository, please [email the DART group](mailto:dart@ucar.edu).

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

DATA SOURCES AND FORMATS
------------------------

See the various subdirectories here, which generally include information
on where the example data was obtained and in what format it is
distributed. Most data is available for download off the web. The Data
Support Section (DSS) at NCAR has large data repositories, the MADIS
data center distributes observations in NetCDF format, GTS real-time
weather data is available from various sources. For new converters, if
you can find what format the data is distributed in you may be able to
adapt one of the existing converters here for your own use. Formats read
by the existing converters include NetCDF, HDF, little-r, text,
Prepbufr, amongst others.

**See the [Programs](#Programs) section below for a list of the current
converter programs. It might save you from reinventing the wheel.**

If you have looked and none of the existing converters are right for
your data, here are some suggestions for where to start creating a new
converter. Create a new subdirectory in the *observations* directory.
Copy with the recursive option (*cp -r*) one of the existing converters
and adapt to your needs. Our suggestions for which converter to start
from depends on the format of your input observations to be converted.
If your input data format is:

<table>
<tr>
<td>netCDF</td>
<td>Start with the <strong>MADIS</strong> converters, and in particular try the 
        <strong>convert_madis_profiler.f90</strong> file because it is the most 
        straightforward. Another good option is <strong>SST/oi_sst_to_obs.f90</strong></td>
</tr><tr>
<td>Comma separated text</td>
<td>Start with the <strong>Ameriflux</strong> converter.</td>
</tr><tr>
<td>Generic text</td>
<td>Start with the <strong>text</strong> converter.</td>
</tr><tr>
<td> HDF-EOS</td>
<td>Start with the <strong>AIRS</strong> converter.</td>
</tr><tr>
<td>BUFR or prepBUFR</td>
<td>Start with the <strong>NCEP</strong> converter.</td>
</tr><tr>
<td>Dense data, like Satellite swaths</td>
<td>Start with the <strong>tpw</strong> converter, which includes code that 
    averages the raw data in space and time.</td>
</tr><tr>
<td>Ray-path integrated data</td>
<td>Start with the <strong>GPS</strong> converter, which includes code that 
traces a path and integrates values along the ray.</td>
</tr><tr>
<td>World Ocean Database packed ASCII</td><td>Start with the <strong>WOD</strong> converter.</td>
</tr>
</table>

<span id="Decisions" class="anchor"></span>

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

DECISIONS YOU MIGHT NEED TO MAKE
--------------------------------

#### Time

Time enters into the assimilation system in 3 places: the timestamp of
the state vector data (the current model time when this data was
produced), the time of each observation, and the minimum time period the
model should be called to advance (the assimilation window size). The
internal timestepping of the model is unrelated to any of these times
and is outside the scope of the assimilation system.

The basic time type in DART is a pair of integers; one for the day
number and one for the number of seconds. Generally the low order
models, which aren't direct geophysical models, use time directly as a
sequence of days starting at 0 and incrementing in any appropriate
number of seconds or days. The observations assimilated into these
systems do not need to use a calendar.

Observations of a real-world system usually are distributed with a
year/month/day, hour/min/seconds timestamp. There are routines in DART
to convert back and forth between the (day-number/seconds) format and a
variety of (year/month/day) calendars. See [the time manager
documentation](../../assimilation_code/modules/utilities/time_manager_mod.html#time_type)
for more details on how DART stores time information and the types of
available calendars. Some climate models which do long runs (100s or
1000s of years) use a modified calendar for simplicity in computation,
e.g. months which always have 30 days, or no leap years. When trying to
assimilate real observations into these models there may be calendar
issues to solve.

The smallest resolvable unit of time in DART is a second. To model a
system which operates on sub-second time scales the time can be scaled
up by some factor. As long as the observation time, the state data time,
and the minimum model advance time are expressed in the same scaled time
units, there is no problem.

#### Error

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
[perfect\_model\_obs](../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html)
program and an assimilation experiment could be run to generate an
estimate of the error in the model. In practice however most people make
an educated guess on the values of the error and then start with a
larger than expected value and decrease it based on the results of
running some test assimilations. For these tests the namelist for the
[outlier
threshold](../../assimilation_code/programs/filter/filter.html#Namelist)
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
observations were assimilated or rejected, and also at the RMSE ([root
mean squared error](http://www.wikipedia.org/wiki/RMSE)) versus the
total spread. DART includes Matlab diagnostic routines to create these
types of plots. The observation RMSE and total spread should be roughly
commensurate. The total spread includes contributions from both the
ensemble variance and the observational error variance, so it can be
adjusted by changing the error values on the incoming observations.
There are other ways to adjust the ensemble spread, including
[inflation](../../assimilation_code/programs/filter/filter.html#Inflation),
so the observation error is not the only factor to consider.

One last recommendation: if possible, the Prior forward operator values
should be compared against the observations after several assimilation
cycles. If you plot results using the Posterior values it is always
possible for the assimilation to overfit the observations and look good
on the diagnostic plots. But the actual test is to then advance the
model and look at how the forecast of the state compares to the
observations.

#### Types

All observations have to have a specific 'type'. There are namelist
controls to turn on and off the assimilation of observations at run-time
by type, or to only evaluate the forward operator for an observation but
have no impact on the state. Several of the diagnostics also group
observations by type to give aggregate statistics after an assimilation.
Generally types are based on both the observing platform or instrument
as well as the kind of observation, e.g. RADIOSONDE\_TEMPERATURE,
ARGO\_SALINITY, etc. Each type is associated with a single underlying
generic 'kind', which controls what forward operator code is called
inside the model, e.g. QTY\_TEMPERATURE, QTY\_DENSITY, etc.

See [here](../forward_operators/obs_def_mod.html) for more details on
how to use and add new DART types. The DART obs\_kind\_mod.f90 defines a
list of already defined observation kinds, and users can either use
existing observation types in 'obs\_def\_xxx\_mod.f90' files, or define
their own.

#### Locations

The two most common choices for specifying the location of an
observation are the
[threed\_sphere](../../assimilation_code/location/threed_sphere/location_mod.html)
and the [oned](../../assimilation_code/location/oned/location_mod.html)
locations. For observations of a real-world system, the 3D Sphere is
generally the best choice. For low-order, 1D models, the 1D locations
are the most commonly used. The observation locations need to match the
type of locations used in the model.

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

PROGRAMS
--------

The *DART/observations/obs\_converters* directory contains a variety of
converter programs to read various external formats and convert the
observations into the format required by DART.

The current list of converters (some directories contain multiple converters) include:

-   [AIRS](AIRS/README.md) [atmospheric variables](AIRS/AIRS.html) and [AMSUA radiances](AIRS/README.md)
-   AURA (uses a combination of IDL and Fortran)
-   [Aviso+/CMEMS](AVISO/AVISO.html)
-   [Ameriflux](Ameriflux/level4_to_obs.html)
-   [CHAMP](CHAMP/work/README)
-   CNOFS
-   [CONAGUA](README)
-   [COSMOS](COSMOS/COSMOS_to_obs.html)
-   [DWL](DWL/dwl_to_obs.html)
-   [GMI](GMI/README.md)
-   [GOES](GOES/README.md)
-   [GPSPW](GPSPW/README)
-   GRACE
-   [GSI2DART](GSI2DART/README)
-   [GTSPP](GTSPP/GTSPP.html)
-   [MADIS](MADIS/MADIS.html)
-   [MIDAS](MIDAS/MIDAS_to_obs.html)
-   [MODIS](MODIS/MOD15A2_to_obs.htm)
-   [MPD](MPD/README.md)
-   [NCEP (prepbufr-&gt;ascii)](NCEP/prep_bufr/prep_bufr.html)
-   [NCEP (ascii-&gt;obs\_seq)](NCEP/ascii_to_obs/create_real_obs.html)
-   [ROMS](ROMS/ROMS.htm)
-   [SSEC](SSEC/SSEC.html)
-   [SST](SST/SST.html)
-   [SSUSI](SSUSI/convert_f16_edr_dsk.html)
-   [WOD](WOD/WOD.html)
-   [cice](cice/cice_to_obs.html)
-   [gnd\_gps\_vtec](gnd_gps_vtec/README)
-   [GPS](gps/gps.html)
-   [ok\_mesonet](ok_mesonet/ok_mesonet.html)
-   [QuikSCAT](quikscat/QuikSCAT.html)
-   [Radar](radar/radar.html)
-   [snow](snow/snow_to_obs.html)
-   [Text](text/text_to_obs.html)
-   text_GITM
-   [tpw](tpw/tpw.html)
-   [Tropical Cyclones](tropical_cyclone/tc_to_obs.html)
-   [Var (little-r)](var/littler_tf_dart.html)
-   [Var (radar)](var/rad_3dvar_to_dart.html)

There are also a couple utilities of note:

-   [even\_sphere](even_sphere/README) - a utility for generating
    evenly-spaced observation locations that can then be used in a
    perfect model experiment.
-   [obs\_error](obs_error/README) - modules that specify observation
    errors based on what is used by ECMWF and NCEP

In addition the following external program produces DART observation
sequence files:

-   [Observation Processing And Wind Synthesis
    (OPAWS)](http://code.google.com/p/opaws/): OPAWS can process NCAR
    Dorade (sweep) and NCAR EOL Foray (netcdf) radar data. It analyzes
    (grids) data in either two-dimensions (on the conical surface of
    each sweep) or three-dimensions (Cartesian). Analyses are output in
    netcdf, Vis5d, and/or DART (Data Assimilation Research Testbed)
    formats.

For generating synthetic observations, see the
[create\_obs\_sequence](../../assimilation_code/programs/create_obs_sequence/create_obs_sequence.html)
program documentation. You can also generate observation files based on
text input. See the [text\_to\_obs](text/text_to_obs.html) program
documentation. Or for simulating a large complex observing system, you
can use the DART library routines in a Fortran program to compute the
observation information and have the DART routines write the output
file.

See the
[perfect\_model](../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html)
program documentation on how to run a model with a set of observations
that have only locations, types, and times, and have the forward
operators compute the observation values.


<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

KNOWN BUGS
----------

 

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

FUTURE PLANS
------------

Contact the [DART development group](mailto:dart@ucar.edu) if you have
observations in a different format that you want to convert. We can give
you advice and pointers on how to approach writing the code.

<div class="top">

\[[top](#)\]

</div>

------------------------------------------------------------------------

Terms of Use
------------

DART software - Copyright UCAR. This open source software is provided by
UCAR, "as is", without charge, subject to all terms of use at
<http://www.image.ucar.edu/DAReS/DART/DART_download>
