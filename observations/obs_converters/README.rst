DART Observations
=================

.. toctree::
   :maxdepth: 2
   :caption: Contents

   overview
   converter_programs
   adding_observations
   converter_automation
   decisions
   data_sources
   converting_series


Overview
========

Real-world observations of earth-system data come from a variety of sources, including radiosondes, satellites, ships, aircraft, weather stations, etc. The files in this *observations* directory can be used to convert data from native formats into a common DART *observation sequence* format.

Synthetic observations are those not based on an actual instrument reading of a system but instead are fabricated to have a known value or have values computed by running a model, possibly with a fixed amount of simulated noise added. These observations can be used for testing, determining the model's sensitivity to assimilation, and designing new observation systems. The DART system includes several ways to create synthetic observations. See the :ref:`programs <programs>` section below for more details.

The DART framework enforces a clean separation between observations and the models into which they are assimilated. The same observations can be used in any model that understands how to generate a value for the requested type of observation from its state space values.

In many cases, a single, self-contained program can convert directly from the observation location, time, value, and error into the DART format. In other cases, especially those linking with a complicated external library (e.g., BUFR), there is a two-step process with two programs and an ASCII intermediate file. We are currently leaning towards single-step conversions, but either approach can be used for new programs.

Frequently the original datasets are in a standard scientific format like netCDF, HDF, or BUFR, and library routines for those formats can be used to read in the original observation data.

The DART software distribution includes Fortran subroutines and functions to help create a sequence of observations in memory, and then a call to the DART observation sequence write routine will create an entire *obs_seq* file in the correct format.

The DART system has several location modules for appropriately computing distances. Two of the ones most commonly used are for data in a 1D system and for data in a 3D spherical coordinate system. Most of the programs here assume the *location/threed_sphere/location_mod.f90* 3D sphere location module is being used.

We are currently collecting information and conversion programs for some additional observation sources and types, which will eventually be added to this directory. In the meantime, if you have converters for data or are interested in something that is not in the repository, please email the DART group.


Converter Programs
==================

The *DART/observations/obs_converters* directory contains a variety of converter programs to read various external formats and convert the observations into the format required by DART.

The current list of converters (some directories contain multiple converters) include:

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
-  ``text_GITM``: See ``./text_GITM``
-  ``tpw``: :doc:`./tpw/tpw`
-  ``Tropical Cyclones``: :doc:`./tropical_cyclone/tc_to_obs`
-  ``3DVAR/4DVAR``: :doc:`./var/var`
-  ``Var (little-r)``: :doc:`./var/littler_tf_dart`
-  ``Var (radar)``: :doc:`./var/rad_3dvar_to_dart`

In addition, the following external program produces DART observation sequence files:

-  `Observation Processing And Wind Synthesis (OPAWS) <http://code.google.com/p/opaws/>`__: OPAWS can process NSF NCAR Dorade (sweep) and NSF NCAR EOL Foray (netcdf) radar data. It analyzes (grids) data in either two-dimensions (on the conical surface of each sweep) or three-dimensions (Cartesian). Analyses are output in netcdf, Vis5d, and/or DART (Data Assimilation Research Testbed) formats.

For generating synthetic observations, see the `create_obs_sequence <../../assimilation_code/programs/create_obs_sequence/create_obs_sequence.html>`__ program documentation. You can also generate observation files based on text input. See the `text_to_obs <text/text_to_obs.html>`__ program documentation and even_sphere. Or, for simulating a large complex observing system, you can use the DART library routines in a Fortran program to compute the observation information and have the DART routines write the output file.

There are a couple of notable utilities of note:

-  `even_sphere <even_sphere/README.html>`__ - a utility for generating a text file of evenly-spaced observation locations that can then be used in a perfect model experiment.
-  `obs_error <obs_error/README.html>`__ - modules that specify observation errors based on what is used by ECMWF and NCEP

See the `perfect_model <../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html>`__ program documentation on how to run a model with a set of observations that have only locations, types, and times and have the forward operators compute the observation values.

Adding Your Observations to DART
================================

First, you should understand that DART already supports a tremendous variety of observations. To fully support an observation means that the observation can be converted from its native format to the DART observation sequence format and that the observation forward operator is already implemented. Keep in mind that forward operators are not specific to any one model.

The observation converters are in the *observations/obs_converter* directory, and you should look there for the documentation describing which converters are available.

The forward operators are functionally or logically grouped into Fortran modules in the *observations/forward_operator* directory. DART employs a ‘contractual’ style of programming in that the forward operator requests information from the model, and if the model cannot provide it, the forward operator may request different information in an attempt to collect the information needed to apply the operator. If the model cannot provide any of the required information, the forward operator fails, the DART QC for that observation is set to the appropriate value, and the program continues.


Observation Converter Automation
================================

We have automated parts of the process to streamline the creation of observation converters. The automation tool helps generate observation converters by allowing users to input the name and format their observational data. The tool then produces a draft converter script that should be further customized based on the specific variable and its origin. The automation tool currently supports text, netCDF, HDF, HDF5, and CSV file formats. 

How to Use the Automation Tool
------------------------------

1. **Navigate to the DART/observations/obs_converters directory**: Open a terminal and navigate to the obs_converters directory where the new_converter script is located.

2. **Run the new_converter script**: Use the script to create a new observation converter. Replace `OBSERVATION_NAME` with the name of your observations and `FILE_FORMAT` with the data format you are working with. Currently, only text, netCDF, HDF, HDF5, and CSV formats are supported. See :ref:`Data Sources and Formats <data_sources_and_formats>` if your file type is not supported.

   .. code-block:: sh

      ./new_converter [OBSERVATION_NAME] [FILE_FORMAT]

3. **Review Generated Files**: The script will create a new directory with the necessary files for your observation converter.

   .. code-block:: none

      DART/observations/obs_converters/OBSERVATION_NAME/
          observation_name_to_obs.f90
          observation_name_to_obs.rst
          work/
              quickbuild.sh
              input.nml

4. **Customize the Fortran Code**: Open `observation_name_to_obs.f90` and customize it according to the specifics of your observation data. The template includes boilerplate code to get you started, such as reading latitude, longitude, and variables from netCDF files.

5. **Edit the Documentation**: Update `observation_name_to_obs.rst` with details specific to your new observation converter. This file serves as the documentation for your converter.

6. **Set Up the Work Directory**: The `work` directory contains a `quickbuild.sh` script and an `input.nml` file. Adjust these as necessary for your environment, the nature of your data, and the goals of your analysis.

7. **Compile and Test**: Use the `quickbuild.sh` script to compile your new observation converter. Test it by running conversions on sample data and verifying the output.

   .. code-block:: sh

      cd work
      ./quickbuild.sh

8. **Validate the Converter**: Test conversions and compare the results with expected outputs to ensure your converter works correctly.


Decisions You May Need to Make
------------------------------

Time
~~~~

Time enters into the assimilation system in 3 places: the time of the state vector data (the current model time when this data was produced), the time of each observation, and the assimilation window length. The window length is set by the model-dependent routine ``shortest_time_between_assimilations()``. The internal time-stepping of the model is unrelated to any of these times and is outside the scope of the assimilation system.

The basic time type in DART is a pair of integers: one for the day number and one for the number of seconds. Generally, the low-order models, which aren’t direct geophysical models, use time directly as a sequence of days starting at 0 and incrementing in any appropriate number of seconds or days. The observations assimilated into these systems do not need to use a calendar.

Observations of a real-world system are usually distributed with a year/month/day, hour/min/second timestamp. There are routines in DART to convert back and forth between the (day-number/seconds) format and a variety of (year/month/day) calendars. See `the time manager documentation <../../assimilation_code/modules/utilities/time_manager_mod.html#time_type>`__ for more details on how DART stores time information and the types of available calendars. Some climate models that do long runs (100s or 1000s of years) use a modified calendar for simplicity in computation, e.g., months, which always have 30 days or no leap years. When trying to assimilate real observations into these models, calendar issues may need to be solved.

The smallest resolvable unit of time in DART is a second. To model a system that operates on sub-second time scales, the time can be scaled up by some factor. As long as the observation time, the state data time, and the minimum model advance time are expressed in the same scaled time units, there is no problem.

Error
~~~~~

Observations must specify an associated expected error. Each individual observation stores its own error value so that it can be a constant value for all observations of that type or vary by location, height, magnitude of the observed value, etc. This value is the expected instrument error plus the representativeness error of the model. The model error includes deficiencies in the equations representing the system's processes and errors introduced by representing a continuous system as a series of discrete points. While the instrument error and the representativeness error could be specified separately, they each have the same impact on the assimilation and can be difficult to determine with any real accuracy. For simplicity, in DART (and most current assimilation software), they are combined and specified as a single value.

The instrument error is generally supplied by the instrument maker. Sadly, it is frequently surprisingly difficult to find these values. For the representativeness error, a set of artificial observations could be generated with the `perfect_model_obs <../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html>`__ program, and an assimilation experiment could be run to generate an estimate of the error in the model. In practice, however, most people make an educated guess on the values of the error and then start with a larger-than-expected value and decrease it based on the results of running some test assimilations. For these tests, the namelist for the `outlier threshold <../../assimilation_code/programs/filter/filter.html#Namelist>`__ should be disabled by setting it to -1 (the default value is 3). This value controls whether the observation is rejected because the observed value is too far from the ensemble mean.

If the diagnostics show that the difference between the mean of the forward operators and the observed value is consistently smaller than the specified observation error, then the error is probably too large. A too-large error reduces the impact of an observation on the state. If the specified observation error is too small, the observation will likely be rejected when the outlier threshold is enabled, and the observation will not be assimilated. It is important to look at the output observation sequence files after assimilation to see how many observations were assimilated or rejected, and also at the RMSE (`root mean squared error <http://www.wikipedia.org/wiki/RMSE>`__) versus the total spread. DART includes Matlab diagnostic routines to create these types of plots. The observation RMSE and total spread should be roughly commensurate. The total spread includes contributions from both the ensemble variance and the observational error variance, so it can be adjusted by changing the error values on the incoming observations. There are other ways to adjust the ensemble spread, including `inflation <../../assimilation_code/programs/filter/filter.html#Inflation>`__, so the observation error is not the only factor to consider.

One last recommendation: if possible, the Prior forward operator values should be compared against the observations after several assimilation cycles. If you plot results using the Posterior values, it is always possible for the assimilation to overfit the observations and look good on the diagnostic plots. However, the actual test is to advance the model then and look at how the forecast of the state compares to the observations.

Types
~~~~~

All observations have to have a specific ‘type’. There are namelist controls to turn on and off the assimilation of observations at run-time by type or to only evaluate the forward operator for observation but have no impact on the state. Several of the diagnostics also group observations by type to give aggregate statistics after an assimilation. Generally, types are based on both the observing platform or instrument as well as the kind of observation, e.g., RADIOSONDE_TEMPERATURE, ARGO_SALINITY, etc. Each type is associated with a single underlying generic ‘kind,’ which controls what forward operator code is called inside the model, e.g. QTY_TEMPERATURE, QTY_DENSITY, etc.

See `here <../forward_operators/obs_def_mod.html>`__ for more details on how to use and add new DART types. The DART obs_kind_mod.f90 defines a list of already defined observation kinds, and users can either use existing observation types in ‘obs_def_xxx_mod.f90’ files or define their own.

Locations
~~~~~~~~~

The two most common choices for specifying the location of an observation are the `threed_sphere <../../assimilation_code/location/threed_sphere/location_mod.html>`__
and the `oned <../../assimilation_code/location/oned/location_mod.html>`__ locations. For observations of a real-world system, the 3D Sphere is generally the best choice. For low-order, 1D models, the 1D locations are the most commonly used. The observation locations must match the locations used in the model.

.. _data_sources_and_formats:
Data Sources and Formats
========================

See the various subdirectories here, which generally include information on where the example data was obtained and in what format it is distributed. Most data is available for download from the web. The Data Support Section (DSS) at NSF NCAR has large data repositories, the MADIS data center distributes observations in NetCDF format, and GTS real-time weather data is available from various sources. For new converters, if you can find the format in which the data is distributed, you may be able to adapt one of the existing converters here for your own use. Formats read by the existing converters include NetCDF, HDF, little-r, text, Prepbufr, amongst others.

**See the** :ref:`programs <programs>` **section below for a list of the current converter programs. It might save you from reinventing the wheel.**

If none of the existing converters are suitable for your data or the automation tool does not support your specific data type, start by creating a new subdirectory in the *observations* directory. Copy one of the existing converters with the recursive option (*cp -r*) and adapt it to your needs. Our suggestions for which converter to start from depend on the format of your input observations to be converted:

+-----------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| netCDF                            | Start with the MADIS converters, and in particular, try the convert_madis_profiler.f90 file because it is the most straightforward. Another good option is SST/oi_sst_to_obs.f90 |
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

Contact the `DART development group <mailto:dart@ucar.edu>`__ if you have observations in a different format that you want to convert. We can advise you on how to approach writing the code.

Converting a Series of Observations
===================================

If you are running a series of assimilation steps, you may need a separate observation sequence (obs_seq) file per step. We suggest creating the first few files by hand to check the resulting obs_seq files and then writing scripts (python, shell) to automate the creation of the remainder of the files. The following are some considerations to take into account when creating scripts for a series of obs_seq files.

Looping in Time
---------------

Often, observations are distributed in files containing observations from a particular time period, e.g., a file per day or week. The output obs_seq files need to include observations from the same time period as the assimilation window: how often the assimilation is stopped and how often the model is advanced in time. The conversion process can convert all the observations from an input file into a single output file and, in a subsequent step, break the file into the required time ranges. Alternatively, the conversion process can extract and convert only the observations required for a single output file and loop multiple times over the same input file.

Generally, earth system models use calendar dates, including months, days, years, hours, minutes, and seconds. The `advance_time` program is very useful in adding or subtracting time periods from calendar dates, considering changing months and years, accounting for leap days, etc.

Observation conversion programs usually take one of two strategies for their input and output filenames.

- Have fixed input and output filenames for the converter. Have the script make symbolic links from the actual filenames to the fixed names for the files for each conversion run.
- Have a Fortran namelist variable that sets the input and output filenames for the converter. Have the script generate or edit the namelist file (e.g., with the `sed` stream editor) to set the actual filenames for each conversion run.

Generally, it is a good idea to encode the date information in the output filename so each file is guaranteed to be unique. This can also make it simpler at filter runtime to generate the required input observation sequence filenames using a program like `advance_time`.

Multiple Observation Files
--------------------------

It is common for an assimilation to want to use observations from different sources. Generally, it is easier to convert observations from each source separately and then merge them together with the `obs_sequence_tool`. Creating filenames and directory names that follow a pattern that can be generated with the `advance_time` program makes this easier to do.

The `obs_sequence_tool` can read the input filenames from a separate ASCII file. This makes generating the filenames easy from a script; it can simply concatenate the input filenames echoed to an ASCII file and then run the obs_sequence_tool. The output file can either be set by using `sed` on the namelist, or a fixed output filename can be used, and then the file is renamed after the tool has run.

Conversion Run Time for Large File Counts
-----------------------------------------

If hundreds of files need to be generated and a supercomputer or other multiple-CPU resource is available, batch files, which convert multiple files simultaneously, can save a lot of time. Ensure that each conversion has its own settings and unique filenames. A separate working directory from other conversions running at the same time can simplify scripting needs.

Verification
------------

Observations taken from real-world sources can have missing values, illegal values, missing files, duplicated data, etc. Writing or adapting programs like `obs_info` can be useful for printing out the first and last obs times in a file, counting each obs type, etc. It is easy to find truncated data files, especially for observations that are close to the start/end of a month or year.

If converting a large number of files, it is common for computer system failures to occur at random times. File systems fill up, batch jobs exit early, and power glitches stop programs before they finish. When verifying, look for anomalous observation counts, unexpected first and last times of obs in a file, missing files, files with many fewer bytes than others, and anything else you can think of.

Output Formats
--------------

There are options to write output obs_seq files in binary, which are roughly half the size of ASCII files. However, it greatly increases the effort to examine the contents of a file for problems. Generally, we have used the ASCII format. It is portable between systems of different "endians" (order of bytes in a multi-byte number) and can be browsed much more easily.

