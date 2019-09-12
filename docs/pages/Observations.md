---
title: Observations
layout: default
---

# DART Observation support and processing

## Links to major sections of this document:
- [An overview of the DART 'preprocess' program](#preprocess)
- [An overview of the observation sequence](#obs_seq_overview)
    - [Creating observations and sequences](#obs_creation)
        - [synthetic observations](#obs_seq_osse)
        - [real observations](#obs_real)
- [Working with observation sequences](#obs_seq_manip)
- [The difference between observation TYPE and QUANTITY](#adding_types)  

~~The DART distribution includes a full set of documentation. Once you
download DART, you may view the documentation offline by opening the
*index.html* file in the top-level DART directory. If you want to
explore the documentation page without downloading DART, you may
[view the documentation for the Manhattan release](Manhattan/documentation/index.html).~~

<span id="preprocess" class="anchor"></span>

---

## An overview of the DART 'preprocess' program

~~First and foremost, check out
[preprocess.html](https://ncar.github.io/DART/api/v2.1.10/program/preprocess.html)
for detailed information.~~

**The *preprocess* program actually *builds* the source code that supports
observations.** That source code is then used by all
the remaining modules. It is **imperative** to actually **run**
*preprocess* before building any executables. This is how the same code
can assimilate synthetic 'observations' for the Lorenz_63 model and
real radar reflectivities for WRF without needing to specify a set of
radar operators for the Lorenz_63 model\!  

*preprocess* combines multiple 'obs_def' modules into one
`obs_def_mod.f90` that is then used by the rest of DART. Additionally,
a new `obs_kind_mod.f90` is built that will provide support for
associating the specific observation **TYPES** with corresponding
(generic) observation **QUANTITIES**. More on that later. The list of
source codes is contained in the *&preprocess_nml* namelist and they
ultimately determine what observations and operators are supported. If
you want to add another 'obs_def' module, you **must** rerun
*preprocess* and recompile the rest of your project. *preprocess* is
designed to abort if the files it is supposed to build already exist.
For this reason, it is necessary to remove a couple files (if they
exist) before you run the preprocessor. It is just a good habit to
develop and is done automatically if you
use the DART *quickbuild.csh* script.

~~~
\rm -f ../../../observations/forward_operators/obs_def_mod.f90
\rm -f ../../../observations/forward_operators/obs_kind_mod.f90
./preprocess
ls -l ../../../observations/forward_operators/obs_def_mod.f90
ls -l ../../../observations/forward_operators/obs_kind_mod.f90
~~~

For example, with a namelist that looks like:

~~~
&preprocess_nml
    input_obs_kind_mod_file  = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90'
    output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90'
    input_obs_def_mod_file   = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
    input_files              = '../../../observations/forward_operators/obs_def_gps_mod.f90',
                               '../../../observations/forward_operators/obs_def_QuikSCAT_mod.f90',
                               '../../../observations/forward_operators/obs_def_GWD_mod.f90',
                               '../../../observations/forward_operators/obs_def_altimeter_mod.f90',
                               '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90'
    output_obs_def_mod_file  = '../../../observations/forward_operators/obs_def_mod.f90'
    /
~~~

*preprocess* will combine `DEFAULT_obs_def_mod.F90, obs_def_gps_mod.f90,
obs_def_QuikSCAT_mod.f90, obs_def_GWD_mod.f90, obs_def_altimeter_mod.f90`,
and `obs_def_reanalysis_bufr_mod.f90`, into `obs_def_mod.f90` -
which can be used by the rest of the project.

### Building and Running 'preprocess'

*preprocess* is an executable, so it should come as no surprise that it
must be built in the normal DART fashion. The
`DART/build_templates/mkmf.template` must be correct for your
environment, and the *input.nml* must have your desired preprocess_nml
set correctly. Given that ...  

~~~
csh mkmf_preprocess  
make  
./preprocess
~~~

will build and run *preprocess*.  

The first command generates an appropriate `Makefile` and the
`input.nml.preprocess_default` file. The second command results in the
compilation of a series of Fortran90 modules which ultimately produces
an executable file: *preprocess*. The third command actually runs
preprocess - which builds the new `obs_kind_mod.f90` and
`obs_def_mod.f90` source code files. The rest of DART may now be built.

### The rationale for 'preprocess'

**IMPORTANT**: Since each 'observation quantity' may require different
amounts of metadata to be read or written; any routine to read or write
an observation sequence **must** be compiled with support for those
particular observations. This is the **whole point** of the
'preprocess' process\! The supported observations are listed in the
`input.nml&obs_kind_nml` block.

For example, radar observations need extra metadata to specify the
location of the radar in addition to the location of the observation,
radiosondes only require the observation location. GPS occultations need
the locations of the two satellites so the forward operator can
integrate along the raypath, cosmic ray soil moisture sensors (yes, they
exist) have forward operators that require site-specific calibration
parameters that are not part of the model and must be included in the
observation metadata. That sort of thing.

For this reason, we strongly recommend that you use the DART routines
to read and process DART observation sequence files.  

<span id="obs_seq_overview" class="anchor"></span>  

---

## An overview of the observation sequence

Observation sequences are complicated, there's just no better way to
describe it. Trying to automatically accomodate a myriad of observation
file formats, structure, and metadata is simply not an easy task. For
this reason, DART has its own format for observations and a set of
programs to convert observations from their original formats to DART's
format. There are definitely some things to know ...  

An `obs_seq.in` file actually contains no observation quantities.
It may be best thought of as a **perfectly**-laid-out notebook - just
waiting for an observer to fill in the actual observation quantities.
All the rows and columns are ready, labelled, and repeated for every
observation time and platform. This file is generally the start of a
"perfect model" experiment. Essentially, one instance of the model is
run through *perfect_model_obs* - which applies the appropriate
forward operators to the model state and 'writes them down' in our
notebook. The completed notebook is renamed `obs_seq.out`.  

An `obs_seq.out` file contains a linked list of observations -
potentially (and usually) observations from different platforms and of
different quantities - each with their own error characteristics and
metadata. These files arise from running *perfect_model_obs* **OR**
from any number of converter programs. The creation of observation
sequences from real observations is not automatic and an email to the
DART team asking for advice for your specific types of observations is
perfectly within reason.  

There is something called an `obs_seq.final` file - which contains
everything in the `obs_seq.out` file as well as a few additional
'copies' of the observation. Remember, DART is an ensemble algorithm.
Each ensemble member must compute its own estimate of the observation
for the algorithm. The `obs_seq.final` file *may* contain each of these
estimates (namelist controlled). Minimally, the mean and spread of the
ensemble estimates is recorded in the `obs_seq.final` file. The best
method of determining the performance of your 'real world' experiment is
to compare in *observation-space* since we can never know the model
state that perfectly represents the real world.  

**IMPORTANT**: Since each 'observation type' may require different
amounts of metadata to be read or written; any routine to read or write
an observation sequence **must** be compiled with support for those
particular observations. The supported observations are listed in the
`input.nml&obs_kind_nml` block. This is the **whole point** of the
'preprocess' process ...  


<table>
<tr>
<th>observation sequence file structure</th>
<th>obs_seq.out</th>
<th>obs_seq.final</th>
</tr>
<tr>
<td valign="top">There are extensible parts of the observation sequence
file; for example, the number of observation kinds contained in the
file, whether the locations have 1 or more components, how
many quality control values are available for each
observation, where those quality control values come from,
how many 'copies' of each observation there
are ... et cetera. The images to the right are links to
full-size images. **They are from entirely separate
experiments. They are just meant to show the flexibility
of the file format.**</td>
<td valign="top"><a href="../images/science_nuggets/obs_seq_out_diagram.png"><img
     alt="The structure of an obs_seq.out file"
     border="0" width="300"
     src="../images/science_nuggets/obs_seq_out_diagram.png" /></a></td>
<td valign="top"><a href="../images/science_nuggets/obs_seq_final_diagram.png"><img
     alt="The structure of an obs_seq.final file"
     border="0" width="300"
     src="../images/science_nuggets/obs_seq_final_diagram.png" /></a></td></tr>
</table>


<span id="obs_creation" class="anchor"></span>

---

## Creating observations and sequences.

**It is strongly encouraged that you use a single observation to test a
new model implementation.**  

Experience has shown that starting 'simple' is the fastest way to good
results. Starting with a **single** observation will exercise a
sufficient portion of the procedure and provide insight into where to
spend more effort. Starting with a single *synthetic* observation will
allow you to focus on the more interesting parts of the DART scheme
without getting bogged down in the world of observation data formats.

<span id="obs_seq_osse"></span>  

### Creating a synthetic observation sequence.

There are several steps to create an observation sequence file, which
follows directly from the modular nature of the DART programming
philosophy.

1.  Decide what observations you want to investigate and edit the
    `input.nml&obs_kind_nml` block.
2.  Build and run *preprocess* to create code that supports the
    observations you want.
3.  Build and run *create_obs_sequence* to define the specifics about
    the observation you want.
4.  Build and run *create_fixed_network_sequence* to replicate those
    specifics through time.
5.  Build and run *perfect_model_obs* to create an observation
    consistent with the model state and specified error distribution at
    the requested times and locations.


<span id="L63_obs_generation"></span>

#### Example: generating observations for the Lorenz '63 model.

1\) There are no 'real' observations for the Lorenz '63 model, so the
appropriate namelist settings are:

~~~
&obs_kind_nml
    assimilate_these_obs_types = 'RAW_STATE_VARIABLE'  /

&preprocess_nml
     input_obs_def_mod_file = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
    output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90'
    input_obs_kind_mod_file = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90'
   output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90'
                input_files = '../../../observations/forward_operators/obs_def_1d_state_mod.f90'
  /
~~~

2\) Run *preprocess* in the normal fashion.

3\) *create_obs_sequence* creates an *observation set definition*
(typically named `set_def.out`), the time-independent part of an
observation sequence. It may help to think of it as trying to define
what sorts of observations will be taken at one 'reading' ... you walk
out to the box and take temperature, humidity, and wind observations all
at the same time and place, for example. You can think of it as one page
in an observer's notebook, and only contains the *location, type,* and
*observational error characteristics* (normally just the diagonal
observational error variance) for a related set of observations. There
are no actual observation values, nor are there any times associated
with the definition. The program is interactive and queries the user for
the information it needs. Begin by creating a minimal observation set
definition in which each of the 3 state variables of L63 is directly
observed with an observational error variance of 1.0 for each
observation. To do this, use the following input sequence (the text
including and after \# is a comment and does not need to be entered):  

The following is a screenshot (much of the verbose logging has been left
off for clarity), the user input looks *like this*.

~~~
   [unixprompt]$ ./create_obs_sequence
    Starting program create_obs_sequence
    Initializing the utilities module.
    Trying to log to unit   10
    Trying to open file dart_log.out

    --------------------------------------
    Starting ... at YYYY MM DD HH MM SS =
                    2017  3 28 10 15 30
    Program create_obs_sequence
    --------------------------------------

    set_nml_output Echo NML values to log file only
    Trying to open namelist log dart_log.nml
    ------------------------------------------------------


    -------------- ASSIMILATE_THESE_OBS_TYPES --------------
    RAW_STATE_VARIABLE
    -------------- EVALUATE_THESE_OBS_TYPES --------------
    ------------------------------------------------------

    ---------- USE_PRECOMPUTED_FO_OBS_TYPES --------------
    ------------------------------------------------------

    Input upper bound on number of observations in sequence
   4
    Input number of copies of data (0 for just a definition)
   0
    Input number of quality control values per field (0 or greater)
   0
    input a -1 if there are no more obs
   0
         Input -1 * state variable index for identity observations
         OR input the name of the observation type from table below:
         OR input the integer index, BUT see documentation...
           1 RAW_STATE_VARIABLE
   -1
    input time in days and seconds
   0 0
    Input error variance for this observation definition
   1.0
    input a -1 if there are no more obs
   0

    { this gets repeated ... until you tell it to stop ... }

    input a -1 if there are no more obs
   -1
    Input filename for sequence (  set_def.out   usually works well)
    set_def.out
    write_obs_seq  opening formatted file set_def.out
    write_obs_seq  closed file set_def.out
~~~

Rest assured that if you requested to assimilate more realistic
observation types, you will be queried for appropriate information by
*create_obs_sequence*. Below is a table that explains all of the input
you should need to supply for observations of the L63 model state.

~~~
4            # upper bound on num of observations in sequence
0            # number of copies of data (0 for just a definition)
0            # number of quality control values per field (0 or greater)
0            # -1 to exit/end observation definitions

-1           # observe state variable 1
0   0        # time -- days, seconds
1.0          # observational variance
0            # -1 to exit/end observation definitions

-2           # observe state variable 2
0   0        # time -- days, seconds
1.0          # observational variance
0            # -1 to exit/end observation definitions

-3           # observe state variable 3
0   0        # time -- days, seconds
1.0          # observational variance
-1           # -1 to exit/end observation definitions

set_def.out  # Output file name
~~~


4\) *create_fixed_network_sequence* takes the observation set
definition and repeats it in time, essentially making multiple pages in
our notebook. Again, the program is interactive and queries the user for
information. You should be able to simply follow the prompts. The table
below represents the input needed for the L63
example:

~~~
set_def.out    # Input observation set definition file          
1              # Regular spaced observation interval in time      
1000           # 1000 observation times
0, 43200       # First observation after 12 hours (0 days, 12 * 3600 seconds)
0, 43200       # Observations every 12 hours
obs_seq.in     # Output file for observation sequence definition
~~~

5\) *perfect_model_obs* advances the model from the state defined by
the initial conditions file specified in the `input.nml` and 'applies
the forward operator' to harvest observations to fill in the observation
sequence specified in `obs_seq.in`. The observation sequence finally
has values for the observations and is saved in a file generally named
*obs_seq.out*. *perfect_model_obs* is namelist-driven, as opposed to
the previous two (whose input is a lot harder to specify in a namelist).
Take a look at (and modify if you like) the
 `input.nml&perfect_model_obs_nml` section of the namelist.  

The End. Not only should you have an observation sequence file (usually
`obs_seq.out`) , you also have a file containing the exact evolution of
the model consistent with those observations - the true state:
`perfect_output.nc`.

<span id="obs_real" class="anchor"></span>

---

## Real Observations - Converting to a DART-compatible format.

[OVERVIEW](#Overview) /
[DATA SOURCES](#DataSources) /
[DECISIONS](#Decisions) /
[PROGRAMS](#Programs) /
[FUTURE PLANS](#FuturePlans)

Real observations come in a mind-boggling diversity of formats. We have
converters for many formats in the `DART/observations/obs_converters`
directory. The documentation for that directory is listed in
[observations.html](Manhattan/observations/obs_converters/observations.html).  

The converters are designed to work on one input file format and create
(or add to) an output observation sequence. It may be desirable to
post-process multiple observation sequence files with the
[obs_sequence_tool](https://ncar.github.io/DART/api/v2.1.10/program/obs_sequence_tool.html)
... to select for timeframe, geographic region, etc.  

Many of the formats require their own libraries (like HDF), and require
intimate knowledge of the data format to extract the portions required
for the [DART observation sequence file](#obs_seq_overview).  Please
feel free to
browse the converters and their companion documentation. Feel free to
donate converters for formats we don't already support\! We like that
kind of stuff.  

The DART framework enforces a clean separation between observations and
the models used for assimilation. The same observations can be used in
any model which understands how to generate a value for the requested
type of observation from the models' state-space values (i.e. the
forward observation operator must exist - DART provides many for the
most common state variables).  

In many cases, the original datasets are in a standard scientific format
like netCDF, HDF, or BUFR, and library routines for those formats can be
used to read in the original observation data. The DART software
distribution includes Fortran subroutines and functions to help create a
sequence of observations in memory, and then a call to the DART
observation sequence write routine will create an entire *obs_seq* file
in the correct format.  

In many cases, a single, self-contained program can convert directly
from the observation location, time, value, and error into the DART
format. In other cases, especially those linking with a complicated
external library (e.g. BUFR), there is a two-step process with two
programs and an ASCII intermediate file. We are currently leaning
towards single-step conversions but either approach can be used for new
programs.  

The DART system comes with several types of location modules for
computing distances appropriately. The two most commonly used are for
data in a 1D system and for data in a 3D spherical coordinate system.
All the programs in the `DART/observations` directory assume the
`assimilation_code/location/threed_sphere/location_mod.f90` 3D sphere
location module is being used.  

With the myriad of observation file formats, HDF, Grib, BUFR, netCDF,
... we simply have not had the time nor need to support all of them. The
converters are a work in progress. There are currently about 10 other
observation sources and types which we are in the process of collecting
information and conversion programs for and which will eventually be
added to this directory. In the meantime, if you have converters for
data or interest in something that is not in the repository, please
email the DART group. Your best bet is to contact our group at
*dart@ucar.edu* with a specific request and we can steer you to the most
similar process.

### Overview

Real-world observations of earth-system data come from a variety of
sources, including radiosondes, satellites, ships, aircraft, weather
stations, etc. The files in this `observations` directory can be used to
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
entire *obs_seq* file in the correct format.

The DART system comes with several types of location modules for
computing distances appropriately. Two of the ones most commonly used
are for data in a 1D system and for data in a 3D spherical coordinate
system. All the programs here assume the
`location/threed_sphere/location_mod.f90` 3D sphere location module is
being used.

There are currently some additional observation sources and types which
we are in the process of collecting information and conversion programs
for and which will eventually be added to this directory. In the
meantime, if you have converters for data or interest in something that
is not in the repository,
please [email the DART group](mailto:dart@ucar.edu).

<span id="DataSources" class="anchor"></span>

-----

### DATA SOURCES AND FORMATS

See the various subdirectories here, which generally include information
on where the example data was obtained and in what format it is
distributed. Most data is available for download off the web. The Data
Support Section (DSS) at NCAR has large data repositories, the MADIS
data center distributes observations in netCDF format, GTS real-time
weather data is available from various sources. For new converters, if
you can find what format the data is distributed in you may be able to
adapt one of the existing converters here for your own use. Formats read
by the existing converters include netCDF, HDF, little-r, text,
Prepbufr, amongst others.

See the [Programs](#Programs) section below for a list of the current
converter programs.

If you have looked and none of the existing converters are right for
your data, here are some suggestions for where to start creating a new
converter. Create a new subdirectory in the *observations* directory.
Copy with the recursive option (*cp -r*) one of the existing converters
and adapt to your needs. Our suggestions for which converter to start
from depends on the format of your input observations to be converted.
If your input data format is:  

| format | advice |  
| :----- | :----- |  
| netCDF | Start with the *MADIS* converters, and in particular try the `convert_madis_profiler.f90` file because it is the most straightforward. Another good option is `SST/oi_sst_to_obs.f90`. |  
| Comma separated text | Start with the *Ameriflux* converter. |  
| Generic text | Start with the *text* converter. |  
| HDF-EOS | Start with the *AIRS* converter. |  
| BUFR or prepBUFR | Start with the *NCEP* converter. |  
| Dense data, like Satellite swaths | Start with the *tpw* converter, which includes code that averages the raw data in space and time. |  
| Ray-path integrated data | Start with the *GPS* converter, which includes code that traces a path and integrates values along the ray. |  
| World Ocean Database packed ASCII | Start with the *WOD* converter. |  


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


<span id="Decisions" class="anchor"></span>

-----

### DECISIONS YOU MIGHT NEED TO MAKE

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
variety of (year/month/day) calendars. See
[the time manager documentation](../../assimilation_code/modules/utilities/time_manager_mod.html#time_type)
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

#### Error Variances

Observations must specify an associated expected error variance. Each individual
observation stores its own error variance value, so it can be a constant value
for all observations of that type or it can vary by location, by height,
by magnitude of the observed value, etc. This value is the expected
instrument error variance plus the representativeness error variance of the model.
The model error variance includes deficiencies in the equations representing the
processes of the system as well as errors introduced by representing a
continuous system as a series of discrete points. While the instrument
error and the representativeness error could be specified separately,
they each have the same impact on the assimilation and can be difficult
to determine with any real accuracy. For simplicity, in DART (and most
current assimilation software) they are combined and specified as a
single value, which we frequently call the 'observation error'. Keep in
mind we really mean 'observation error variance'.

The instrument error is generally supplied by the instrument maker.
Sadly, it is frequently surprisingly difficult to find these values. For
the representativeness error, a set of artificial observations could be
generated with the
[perfect_model_obs](../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html)
program and an assimilation experiment could be run to generate an
estimate of the error in the model. In practice however most people make
an educated guess on the values of the error and then start with a
larger than expected value and decrease it based on the results of
running some test assimilations. For these tests the namelist for the
[outlier threshold](../../assimilation_code/programs/filter/filter.html#Namelist)
should be disabled by setting it to -1 (the default value is 3). This
value controls whether the observation is rejected because the observed
value is too far from the ensemble mean.

If the diagnostics show that the difference between the mean of the
forward operators and the observed value is consistently smaller than
the specified observation error, then the error is probably too large. A
error that is too large reduces the impact of an observation on the state. If
the specified observation error is too small it is likely the
observation will be rejected when the outlier threshold is enabled, and
the observation will not be assimilated. It is important to look at the
output observation sequence files after an assimilation to see how many
observations were assimilated or rejected, and also at the RMSE
([root mean squared error](http://www.wikipedia.org/wiki/RMSE)) versus the
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

#### Observation Types

All observations have to have a specific 'type'. There are namelist
controls to turn on and off the assimilation of observations at run-time
by type, or to only evaluate the forward operator for an observation but
have no impact on the state. Several of the diagnostics also group
observations by type to give aggregate statistics after an assimilation.
Generally types are based on both the observing platform or instrument
as well as the 'kind' of observation, e.g. RADIOSONDE_TEMPERATURE,
ARGO_SALINITY, etc. Each type is associated with a single underlying
generic 'kind', which controls what forward operator code is called
inside the model, e.g. QTY_TEMPERATURE, QTY_DENSITY, etc.

See the [obs_def_mod.html](../forward_operators/obs_def_mod.html) for more details on
how to use and add new DART types. The DART `obs_kind_mod.f90` defines a
list of already defined observation types, and users can either use
existing observation types in 'obs_def_xxx_mod.f90' files, or define
their own. Be aware that `obs_kind_mod.f90` is autogenerated by *preprocess*,
so until you configure and run *preprocess*, `obs_kind_mod.f90` will not exist.

#### Observation Locations

The two most common choices for specifying the location of an
observation are the
[threed_sphere](../../assimilation_code/location/threed_sphere/location_mod.html)
and the [oned](../../assimilation_code/location/oned/location_mod.html)
locations. For observations of a real-world system, the 3D Sphere is
generally the best choice. For low-order, 1D models, the 1D locations
are the most commonly used. The observation locations need to match the
type of locations used in the model in that you cannot read observations
on a unit circle (1D) when using models that require 3D Sphere locations.

The choice of the vertical coordinate system may also be important.
For the 3D Sphere, the vertical coordinate system choices are:  

| string            | integer value | meaning  |  
| :-----            | :------------ | :------  |  
| VERTISUNDEF       | -2            | has no specific vertical location (undefined) |  
| VERTISSURFACE     | -1            | surface value (value is surface elevation in m) |  
| VERTISLEVEL       |  1            | by model level |  
| VERTISPRESSURE    |  2            | by pressure (in pascals) |  
| VERTISHEIGHT      |  3            | by height (in meters) |  
| VERTISSCALEHEIGHT |  4            | by scale height (unitless) |  

The choice of the vertical coordinate system may have ramifications for vertical
localization, depending on your model's ability to convert from one coordinate
system to another. `VERTISUNDEF` is typically used for column-integrated quantities.
`VERTISLEVEL` only makes sense for synthetic observations.

<span id="Programs" class="anchor"></span>

\[[top](#)\]

-----

### PROGRAMS

The **DART/observations/obs_converters** directory contains a variety of
converter programs to read various external formats and convert the
observations into the format required by DART.

Each directory has at least one converter:

  - [AIRS](AIRS/AIRS.html) <!-- AURA -->
  - [Aviso+/CMEMS](AVISO/AVISO.html)
  - [Ameriflux](Ameriflux/level4_to_obs.html) <!-- CHAMP --> <!-- CNOFS -->
  - [COSMOS](COSMOS/COSMOS_to_obs.html)
  - [DWL](DWL/dwl_to_obs.html)
  - [GPSPW](GPSPW/README)
  - [GSI2DART](GSI2DART/README)
  - [GTSPP](GTSPP/GTSPP.html)
  - [MADIS](MADIS/MADIS.html)
  - [MIDAS](MIDAS/MIDAS_to_obs.html)
  - [MODIS](MODIS/MOD15A2_to_obs.htm)
  - [NCEP (prepbufr -\> ascii)](NCEP/prep_bufr/prep_bufr.html)
  - [NCEP (ascii -\> obs_seq)](NCEP/ascii_to_obs/create_real_obs.html)
  - [ROMS](ROMS/ROMS.htm) <!-- SABER -->
  - [SSEC](SSEC/SSEC.html)
  - [SST](SST/SST.html)
  - [SSUSI](SSUSI/convert_f16_edr_dsk.html)
  - [WOD](WOD/WOD.html)
  - [cice](cice/cice_to_obs.html)
  - [gnd_gps_vtec](gnd_gps_vtec/README)
  - [GPS](gps/gps.html)
  - [ok_mesonet](ok_mesonet/ok_mesonet.html)
  - [QuikSCAT](quikscat/QuikSCAT.html)
  - [Radar](radar/radar.html)
  - [snow](snow/snow_to_obs.html)
  - [Text](text/text_to_obs.html) <!-- text_GITM -->
  - [tpw](tpw/tpw.html)
  - [Tropical Cyclones](tropical_cyclone/tc_to_obs.html)
  - [Var (little-r)](var/littler_tf_dart.html)
  - [Var (radar)](var/rad_3dvar_to_dart.html)

There are also a couple utilities of note:

  - [even_sphere](even_sphere/README) - a utility for generating
    evenly-spaced observation locations that can then be used in a
    perfect model experiment.
  - [obs_error](obs_error/README) - modules that specify observation
    errors based on what is used by ECMWF and NCEP

In addition the following external program produces DART observation
sequence files:

  - [Observation Processing And Wind Synthesis (OPAWS)](http://code.google.com/p/opaws/):
    OPAWS can process NCAR Dorade (sweep) and NCAR EOL Foray (netCDF)
    radar data. It analyzes (grids) data in either two-dimensions
    (on the conical surface of each sweep) or three-dimensions (Cartesian).
    Analyses are output in netCDF, Vis5d, and/or DART
    (Data Assimilation Research Testbed) formats.

For generating synthetic observations, see the
[create_obs_sequence](../../assimilation_code/programs/create_obs_sequence/create_obs_sequence.html)
program documentation. You can also generate observation files based on
text input. See the [text_to_obs](text/text_to_obs.html) program
documentation. Or for simulating a large complex observing system, you
can use the DART library routines in a Fortran program to compute the
observation information and have the DART routines write the output
file.

See the
[perfect_model](../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html)
program documentation on how to run a model with a set of observations
that have only locations, types, and times, and have the forward
operators compute the observation values.

<span id="FuturePlans" class="anchor"></span>

-----

### FUTURE PLANS

Contact the [DART development group](mailto:dart@ucar.edu) if you have
observations in a different format that you want to convert. We can give
you advice and pointers on how to approach writing the code.

<span id="obs_seq_manip" class="anchor"></span>

\[[top](#)\]

---

# Working with observation sequences.

First and foremost, check out the
[obs_sequence_tool.html](https://ncar.github.io/DART/api/v2.1.10/program/obs_sequence_tool.html)
document for detailed information and examples.  

*obs_sequence_tool* is the primary tool for manipulating observation
sequence files. Observations sequence files are linked lists of
observations organized by time. That is to say, the observations may
appear in any order in the file, but traversing the linked list will
result in observations ordered by time. *obs_sequence_tool* can be
used to combine observation sequences, convert from ASCII to binary or
vice-versa, extract a subset of observations, etc.  

For testing, it is terribly useful to extract a small number of
observations (like ONE) from an existing observation sequence file.

<span id="adding_types" class="anchor"></span>

\[[top](#)\]

---

# The difference between observation TYPE and QUANTITY.

Broadly speaking, observation TYPES are specific instances of a generic
observation QUANTITY. The distinction is useful for several reasons, not
the least of which is to evaluate observation platforms. Zonal wind
observations from QuikSCAT vs. radiosondes, for example. They are both
observations of zonal winds (what we call QTY_U_WIND_COMPONENT), but
they are different observation TYPES; QKSWND_U_WIND_COMPONENT, and
RADIOSONDE_U_WIND_COMPONENT, respectively. The forward observation
operators are implemented based on observation QUANTITY. When requested,
the model generates a QTY_U_WIND_COMPONENT, it doesn't need to know
that it will be compared to a QuikSCAT value or a radiosonde value.  

**However**, it is usually scientifically very interesting to be able to
compare the assimilations one TYPE of observation vs. another. One
observation sequence file can have lots of types of observations; DART
has the capability to assimilate (or evaluate) any combination of
observation types without getting bogged down in dataset management. The
same observation sequence can be used for experiments that
include/exclude certain observation types - ensuring that you are
performing the experiment you THINK you are performing
...

# Adding support for a new observation TYPE.

[DART/observations/forward_operators/obs_def_mod.html](https://ncar.github.io/DART/api/v2.1.10/module/obs_def_mod.html)
is the source for detailed information.

---
