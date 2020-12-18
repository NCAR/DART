##########
CAM-FV README
##########

Contents
========

#. `Overview`_
#. `Setup Scripts`_
#. `Namelist`_
#. `Setup Variations`_
#. `Notes for Continuing an Integration`_
#. `Discussion`_
#. `Files`_
#. `Nitty gritty: Efficiency possibilities`_
#. `References and Acknowledgements`_
#. `Terms of Use`_

Overview
========

The DART system supports data assimilation into the Community Atmosphere Model
(CAM) which is the atmospheric component of the Community Earth System Model
(`CESM <http://www2.cesm.ucar.edu/models>`__). This DART interface is being
used by graduate students, post-graduates, and scientists at universities and
research labs to conduct data assimilation reseearch. Others are using the
products of data assimilation (analyses), which were produced here at NCAR
using CESM+DART, to conduct related research. The variety of research can be
sampled on the DART
`Publications <https://dart.ucar.edu/pages/Publications.html>`__ page.

"CAM" refers to a family of related atmospheric components, which can be built
with two independent main characteristics. CESM labels these as:

**resolution**

   where *resolution* refers to both the horizontal resolution of the grid
   (rather than the vertical resolution) **and** the dynamical core run on the
   specified grid. The dynamical core refers to the fluid dynamical equations
   run on the specified grid.

**compset**

   where *compset* refers to the vertical grid **and** the parameterizations --
   the formulation of the subgridscale physics. These parameterizations 
   encompass the equations describing physical processes such as convection,
   radiation, chemistry.
   
   - The vertical grid is determined by the needs of the chosen
     parameterizations, thus the vertical spacing and the top level of the
     model domain, specified by a variable known as ``ptop``, vary.
   - The combinations of parameterizations and vertical grids are named: CAM3.5,
     CAM5, CAM#, ... WACCM, WACCM#, WACCM-X, CAM-Chem.

`Setup Variations`_ describes the differences in the namelists, build scripts
assimilation setup.

This DART+CAM interface has the following features.

-  Assimilate within the CESM software framework by using the multi-instance
   capability of CESM1.1.1 (and later). This enables assimilation of suitable
   observations into multiple CESM components. The ability to assimilate in the
   previous mode, where DART called 'stand-alone' CAMs when needed, is not being
   actively supported for these CESM versions.
-  Support for the finite-volume (FV) dynamical core.
-  Use any resolution of CAM.
-  Assimilate a variety of observations; to date the observations successfully
   assimilated include the NCEP reanalysis BUFR obs (T,U,V,Q), Global
   Positioning System radio occultation obs, and MOPITT carbon monoxide (when a
   chemistry model is incorporated into CAM-FV). Research has also explored
   assimilating surface observations, cloud liquid water, and aerosols. SABER
   and AURA observations have been assimilated into WACCM.
-  Specify, via namelist entries, the CAM (initial file) variables which will be
   directly affected by the observations, that is, the state vector. This allows
   users to change the model state without recompiling (but other restrictions
   remain).
-  Generate analyses on the CAM grid which have only CAM model error in them,
   rather than another model's.
-  Generate such analyses with as few as 20 ensemble members.

In addition to the standard DART package there are ensembles of initial
condition files at the large file website
http://www.image.ucar.edu/pub/DART/CAM/ that are helpful for interfacing CAM
with DART. In the current (2015) mode, CESM+DART can easily be started from a
single model state, which is perturbed to create an ensemble of the desired
size. A spin-up period is then required to allow the ensemble members to
diverge.

Sample sets of observations, which can be used with CESM+DART assimilations, can
be found at http://www.image.ucar.edu/pub/DART/Obs_sets/ of which the NCEP BUFR
observations are the most widely used.

**Reanalyses**

There have been two large-scale reanalysis efforts using CAM and DART. 
The NCAR Research Archive dataset **"CAM6 Data Assimilation Research Testbed (DART) Reanalysis"**
`DS345.0 <https://rda.ucar.edu/datasets/ds345.0/#!description>`__ contains
just under 120Tb (yes Tb) of data:

   These CAM6 Data Assimilation Research Testbed (DART) Reanalysis data 
   products are designed to facilitate a broad variety of research using 
   NCAR's CESM2 models, ranging from model evaluation to (ensemble) 
   hindcasting, data assimilation experiments, and sensitivity studies. 
   They come from an 80 member ensemble reanalysis of the global 
   troposphere and stratosphere using DART and CAM6 from CESM2.1. 
   The data products represent the actual states of the atmosphere 
   during 2 recent decades at a 1 degree horizontal resolution and 
   6 hourly frequency. Each ensemble member is an equally likely 
   description of the atmosphere, and is also consistent with 
   dynamics and physics of CAM6.
   
   
The NCAR Research Archive dataset **"An Ensemble of Atmospheric Forcing Files from a CAM Reanalysis"**
`DS199.1 <https://rda.ucar.edu/datasets/ds199.1/#!description>`__ contains
about 1.5Tb of data:

   This dataset contains files that are an ensemble of 'coupler history' 
   files from an 80-member reanalysis performed with the Data Assimilation 
   Research Testbed (DART) using the Community Atmosphere Model Version 
   4 with the finite volume core (CAM4 FV) at 1.9 degree by 2.5 degree 
   resolution. The observations assimilated include all those used in 
   the NCEP/NCAR reanalysis (temperature and wind components from 
   radiosondes, aircraft, and satellite drift winds) plus radio 
   occultation observations from the COSMIC satellites starting in late 
   2006. These files are intended to be used as 'DATM stream files' 
   for CESM component sets that require a data atmosphere. Some example 
   stream text files are included to illustrate how to use these data.

**Guidance**

Experience on a variety of machines has shown that it is a very good idea to
make sure your run-time environment has the following:

.. code-block:: bash

   limit stacksize unlimited
   limit datasize unlimited

This page contains the documentation for the DART interface module for the CAM
and WACCM models, using the dynamical cores listed above. This implementation
uses the CAM initial files (**not** restart files) for transferring the model
state to/from the filter. This may change in future versions, but probably only
for CAM-SE. The reasons for this include:

#. The contents of the restart files vary depending on both the model release
   version and the physics packages selected.
#. There is no metadata describing the variables in the restart files. Some
   information can be tracked down in the ``atm.log`` file, but not all of it.
#. The restart files (for non-chemistry model versions) are much larger than
   the initial files (and we need to deal with an ensemble of them).
#. The temperature on the restart files is virtual equivalent potential
   temperature, which requires (at least) surface pressure, specific humidity,
   and sensible temperature to calculate.
#. CAM does not call the initialization routines when restart files are used,
   so fields which are not modified by DART may be inconsistent with fields
   which are.
#. If DART modifies the contents of the ``.r.`` restart file, it might also
   need to modify the contents of the ``.rs.`` restart file, which has similar
   characteristics (1-3 above) to the ``.r.`` file.

The DART interfaces to CAM and many of the other CESM components have been
integrated with the CESM set-up and run scripts.

Setup Scripts
=============

Unlike previous versions of DART-CAM, CESM runs using its normal scripts, then
stops and calls a DART script, which runs a single assimilation step, then
returns to the CESM run script to continue the model advances. See the CESM
interface documentation in ``$DARTROOT/models/CESM`` for more information on
running DART with CESM. Due to the complexity of the CESM software environment,
the versions of CESM which can be used for assimilation are more restricted than
previously. Each supported CESM version has similar, but unique, sets of set-up
scripts and CESM SourceMods. Those generally do not affect the
``cam-fv/model_mod.f90`` interface. Current (April, 2015) set-up scripts are:

-  ``CESM1_2_1_setup_pmo``: sets up a perfect_model_mod experiment, which
   creates synthetic observations from a free model run, based on the user's
   somewhat restricted choice of model, dates, etc. The restrictions are made
   in order to streamline the script, which will shorten the learning curve for
   new users.
-  ``CESM1_2_1_setup_pmo_advanced``: same as ``CESM1_2_1_setup_pmo``, but can
   handle more advanced set-ups: recent dates (non-default forcing files),
   refined-grid CAM-SE, etc.
-  ``CESM1_2_1_setup_hybrid``: streamlined script (see ``CESM1_2_1_setup_pmo``)
   which sets up an ensemble assimilation using CESM's multi-instance
   capability.
-  ``CESM1_2_1_setup_advanced``: like ``CESM1_2_1_setup_pmo_advanced``, but for
   setting up an assimilation.

The DART state vector should include all prognostic variables in the CAM
initial files which cannot be calculated directly from other prognostic
variables. In practice the state vector sometimes contains derived quantities to
enable DART to compute forward operators (expected observation values)
efficiently. The derived quantities are often overwritten when the model runs
the next timestep, so the work DART does to update them is wasted work.

Expected observation values on pressure, scale height, height or model levels
can be requested from ``model_interpolate``. Surface observations can not yet be
interpolated, due to the difference between the model surface and the earth's
surface where the observations are made. Model_interpolate can be queried for
any (non-surface) variable in the state vector (which are variables native to
CAM) plus pressure on height levels. The default state vector is PS, T, U, V, Q,
CLDLIQ, CLDICE and any tracers or chemicals needed for a given study. Variables
which are not in the initial file can be added (see the ``./doc`` directory
but minor modifications to ``model_mod.f90`` and CAM may be necessary.

The 19 public interfaces in ``model_mod`` are standardized for all DART
compliant models. These interfaces allow DART to get the model state and
metadata describing this state, find state variables that are close to a given
location, and do spatial interpolation for a variety of variables required by
observational operators.

Namelist
========

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists start
with an ampersand ``&`` and terminate with a slash ``/``. Character strings that
contain a ``/`` must be enclosed in quotes to prevent them from prematurely
terminating the namelist.

.. code-block:: fortran

   &model_nml
      cam_template_filename               = 'caminput.nc'
      cam_phis_filename                   = 'cam_phis.nc'
      vertical_localization_coord         = 'PRESSURE'
      use_log_vertical_scale              = .false.
      no_normalization_of_scale_heights   = .true.
      no_obs_assim_above_level            = -1,
      model_damping_ends_at_level         = -1,
      state_variables                     = ''
      assimilation_period_days            = 0
      assimilation_period_seconds         = 21600
      suppress_grid_info_in_output        = .false.
      custom_routine_to_generate_ensemble = .true.
      fields_to_perturb                   = ''
      perturbation_amplitude              = 0.0_r8
      using_chemistry                     = .false.
      use_variable_mean_mass              = .false.
      debug_level                         = 0
   /  

The names of the fields to put into the state vector must match the CAM initial
NetCDF file variable names.

+-------------------------------------+-----------------------------------+------------------------------------------+
| Item                                | Type                              | Description                              |
+=====================================+===================================+==========================================+
| cam_template_file                   | character(len=128)                | CAM initial file used to provide         |
|                                     |                                   | configuration information, such as the   |
|                                     |                                   | grid resolution, number of vertical      |
|                                     |                                   | levels, whether fields are staggered or  |
|                                     |                                   | not, etc.                                |
+-------------------------------------+-----------------------------------+------------------------------------------+
| cam_phis                            | character(len=128)                | CAM topography file. Reads the "PHIS"    |
|                                     |                                   | NetCDF variable from this file.          |
|                                     |                                   | Typically this is a CAM History file     |
|                                     |                                   | because this field is not normally found |
|                                     |                                   | in a CAM initial file.                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| vertical_localization_coord         | character(len=128)                | The vertical coordinate to which all     |
|                                     |                                   | vertical locations are converted in      |
|                                     |                                   | model_mod. Valid options are "pressure", |
|                                     |                                   | "height", "scaleheight" or "level".      |
+-------------------------------------+-----------------------------------+------------------------------------------+
| no_normalization_of_scale_heights   | logical                           | If true the scale height is computed as  |
|                                     |                                   | the log of the pressure at the given     |
|                                     |                                   | location. If false the scale height is   |
|                                     |                                   | computed as a ratio of the log of the    |
|                                     |                                   | surface pressure and the log of the      |
|                                     |                                   | pressure aloft. In limited areas of high |
|                                     |                                   | topography the ratio version might be    |
|                                     |                                   | advantageous, and in previous versions   |
|                                     |                                   | of filter this was the default. For      |
|                                     |                                   | global CAM the recommendation is to set  |
|                                     |                                   | this to .true. so the scale height is    |
|                                     |                                   | simply the log of the pressure at any    |
|                                     |                                   | location.                                |
+-------------------------------------+-----------------------------------+------------------------------------------+
| no_obs_assim_above_level            | integer                           | Because the top of the model is highly   |
|                                     |                                   | damped it is recommended to NOT          |
|                                     |                                   | assimilate observations in the top model |
|                                     |                                   | levels. The units here are CAM model     |
|                                     |                                   | level numbers. Set it to equal or below  |
|                                     |                                   | the lowest model level (the highest      |
|                                     |                                   | number) where damping is applied in the  |
|                                     |                                   | model.                                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| model_damping_ends_at_level         | integer                           | Set this to the lowest model level (the  |
|                                     |                                   | highest number) where model damping is   |
|                                     |                                   | applied. Observations below the          |
|                                     |                                   | 'no_obs_assim_above_level' cutoff but    |
|                                     |                                   | close enough to the model top to have an |
|                                     |                                   | impact during the assimilation will have |
|                                     |                                   | their impacts decreased smoothly to 0 at |
|                                     |                                   | this given model level. The assimilation |
|                                     |                                   | should make no changes to the model      |
|                                     |                                   | state above the given level.             |
+-------------------------------------+-----------------------------------+------------------------------------------+
| state_variables                     | character(len=64), dimension(100) | Character string table that includes:    |
|                                     |                                   | Names of fields (NetCDF variable names)  |
|                                     |                                   | to be read into the state vector, the    |
|                                     |                                   | corresponding DART Quantity for that     |
|                                     |                                   | variable, if a bounded quantity the      |
|                                     |                                   | minimum and maximum valid values, and    |
|                                     |                                   | finally the string 'UPDATE' to indicate  |
|                                     |                                   | the updated values should be written     |
|                                     |                                   | back to the output file. 'NOUPDATE' will |
|                                     |                                   | skip writing this field at the end of    |
|                                     |                                   | the assimilation.                        |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_days            | integer                           | Sets the assimilation window width, and  |
|                                     |                                   | should match the model advance time when |
|                                     |                                   | cycling. The scripts distributed with    |
|                                     |                                   | DART always set this to 0 days, 21600    |
|                                     |                                   | seconds (6 hours).                       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_seconds         | integer                           | Sets the assimilation window width, and  |
|                                     |                                   | should match the model advance time when |
|                                     |                                   | cycling. The scripts distributed with    |
|                                     |                                   | DART always set this to 0 days, 21600    |
|                                     |                                   | seconds (6 hours).                       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| suppress_grid_info_in_output        | logical                           | Filter can update fields in existing     |
|                                     |                                   | files or create diagnostic/output files  |
|                                     |                                   | from scratch. By default files created   |
|                                     |                                   | from scratch include a full set of CAM   |
|                                     |                                   | grid information to make the file fully  |
|                                     |                                   | self-contained and plottable. However,   |
|                                     |                                   | to save disk space the grid variables    |
|                                     |                                   | can be suppressed in files created by    |
|                                     |                                   | filter by setting this to true.          |
+-------------------------------------+-----------------------------------+------------------------------------------+
| custom_routine_to_generate_ensemble | logical                           | The default perturbation routine in      |
|                                     |                                   | filter adds gaussian noise equally to    |
|                                     |                                   | all fields in the state vector. It is    |
|                                     |                                   | recommended to set this option to true   |
|                                     |                                   | so code in the model_mod is called       |
|                                     |                                   | instead. This allows only a limited      |
|                                     |                                   | number of fields to be perturbed. For    |
|                                     |                                   | example, only perturbing the temperature |
|                                     |                                   | field T with a small amount of noise and |
|                                     |                                   | then running the model forward for a few |
|                                     |                                   | days is often a recommended way to       |
|                                     |                                   | generate an ensemble from a single       |
|                                     |                                   | state.                                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| fields_to_perturb                   | character(len=32), dimension(100) | If perturbing a single state to generate |
|                                     |                                   | an ensemble, set                         |
|                                     |                                   | 'custom_routine_to_generate_ensemble =   |
|                                     |                                   | .true.' and list list the field(s) to be |
|                                     |                                   | perturbed here.                          |
+-------------------------------------+-----------------------------------+------------------------------------------+
| perturbation_amplitude              | real(r8), dimension(100)          | For each field name in the               |
|                                     |                                   | 'fields_to_perturb' list give the        |
|                                     |                                   | standard deviation for the gaussian      |
|                                     |                                   | noise to add to each field being         |
|                                     |                                   | perturbed.                               |
+-------------------------------------+-----------------------------------+------------------------------------------+
| pert_base_vals                      | real(r8), dimension(100)          | If pert_sd is positive, this the list of |
|                                     |                                   | values to which the field(s) listed in   |
|                                     |                                   | pert_names will be reset if filter is    |
|                                     |                                   | told to create an ensemble from a single |
|                                     |                                   | state vector. Otherwise, it's is the     |
|                                     |                                   | list of values to use for each ensemble  |
|                                     |                                   | member when perturbing the single field  |
|                                     |                                   | named in pert_names. Unused unless       |
|                                     |                                   | pert_names is set and pert_base_vals is  |
|                                     |                                   | not the DART missing value.              |
+-------------------------------------+-----------------------------------+------------------------------------------+
| using_chemistry                     | logical                           | If using CAM-CHEM, set this to .true.    |
+-------------------------------------+-----------------------------------+------------------------------------------+
| using_variable_mean_mass            | logical                           | If using any variant of WACCM with a     |
|                                     |                                   | very high model top, set this to .true.  |
+-------------------------------------+-----------------------------------+------------------------------------------+
| debug_level                         | integer                           | Set this to increasingly larger values   |
|                                     |                                   | to print out more debugging information. |
|                                     |                                   | Note that this can be very verbose. Use  |
|                                     |                                   | with care.                               |
+-------------------------------------+-----------------------------------+------------------------------------------+

+-------------------------------------+-----------------------------------+------------------------------------------+
| Item                                | Type                              | Description                              |
+=====================================+===================================+==========================================+
| cam_template_file                   | character(len=128)                | CAM initial file used to provide         |
|                                     |                                   | configuration information, such as the   |
|                                     |                                   | grid resolution, number of vertical      |
|                                     |                                   | levels, whether fields are staggered or  |
|                                     |                                   | not, etc.                                |
+-------------------------------------+-----------------------------------+------------------------------------------+
| cam_phis                            | character(len=128)                | CAM topography file. Reads the "PHIS"    |
|                                     |                                   | NetCDF variable from this file.          |
|                                     |                                   | Typically this is a CAM History file     |
|                                     |                                   | because this field is not normally found |
|                                     |                                   | in a CAM initial file.                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| vertical_localization_coord         | character(len=128)                | The vertical coordinate to which all     |
|                                     |                                   | vertical locations are converted in      |
|                                     |                                   | model_mod. Valid options are "pressure", |
|                                     |                                   | "height", "scaleheight" or "level".      |
+-------------------------------------+-----------------------------------+------------------------------------------+
| no_normalization_of_scale_heights   | logical                           | If true the scale height is computed as  |
|                                     |                                   | the log of the pressure at the given     |
|                                     |                                   | location. If false the scale height is   |
|                                     |                                   | computed as a ratio of the log of the    |
|                                     |                                   | surface pressure and the log of the      |
|                                     |                                   | pressure aloft. In limited areas of high |
|                                     |                                   | topography the ratio version might be    |
|                                     |                                   | advantageous, and in previous versions   |
|                                     |                                   | of filter this was the default. For      |
|                                     |                                   | global CAM the recommendation is to set  |
|                                     |                                   | this to .true. so the scale height is    |
|                                     |                                   | simply the log of the pressure at any    |
|                                     |                                   | location.                                |
+-------------------------------------+-----------------------------------+------------------------------------------+
| no_obs_assim_above_level            | integer                           | Because the top of the model is highly   |
|                                     |                                   | damped it is recommended to NOT          |
|                                     |                                   | assimilate observations in the top model |
|                                     |                                   | levels. The units here are CAM model     |
|                                     |                                   | level numbers. Set it to equal or below  |
|                                     |                                   | the lowest model level (the highest      |
|                                     |                                   | number) where damping is applied in the  |
|                                     |                                   | model.                                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| model_damping_ends_at_level         | integer                           | Set this to the lowest model level (the  |
|                                     |                                   | highest number) where model damping is   |
|                                     |                                   | applied. Observations below the          |
|                                     |                                   | 'no_obs_assim_above_level' cutoff but    |
|                                     |                                   | close enough to the model top to have an |
|                                     |                                   | impact during the assimilation will have |
|                                     |                                   | their impacts decreased smoothly to 0 at |
|                                     |                                   | this given model level. The assimilation |
|                                     |                                   | should make no changes to the model      |
|                                     |                                   | state above the given level.             |
+-------------------------------------+-----------------------------------+------------------------------------------+
| state_variables                     | character(len=64), dimension(100) | Character string table that includes:    |
|                                     |                                   | Names of fields (NetCDF variable names)  |
|                                     |                                   | to be read into the state vector, the    |
|                                     |                                   | corresponding DART Quantity for that     |
|                                     |                                   | variable, if a bounded quantity the      |
|                                     |                                   | minimum and maximum valid values, and    |
|                                     |                                   | finally the string 'UPDATE' to indicate  |
|                                     |                                   | the updated values should be written     |
|                                     |                                   | back to the output file. 'NOUPDATE' will |
|                                     |                                   | skip writing this field at the end of    |
|                                     |                                   | the assimilation.                        |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_days            | integer                           | Sets the assimilation window width, and  |
|                                     |                                   | should match the model advance time when |
|                                     |                                   | cycling. The scripts distributed with    |
|                                     |                                   | DART always set this to 0 days, 21600    |
|                                     |                                   | seconds (6 hours).                       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_seconds         | integer                           | Sets the assimilation window width, and  |
|                                     |                                   | should match the model advance time when |
|                                     |                                   | cycling. The scripts distributed with    |
|                                     |                                   | DART always set this to 0 days, 21600    |
|                                     |                                   | seconds (6 hours).                       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| suppress_grid_info_in_output        | logical                           | Filter can update fields in existing     |
|                                     |                                   | files or create diagnostic/output files  |
|                                     |                                   | from scratch. By default files created   |
|                                     |                                   | from scratch include a full set of CAM   |
|                                     |                                   | grid information to make the file fully  |
|                                     |                                   | self-contained and plottable. However,   |
|                                     |                                   | to save disk space the grid variables    |
|                                     |                                   | can be suppressed in files created by    |
|                                     |                                   | filter by setting this to true.          |
+-------------------------------------+-----------------------------------+------------------------------------------+
| custom_routine_to_generate_ensemble | logical                           | The default perturbation routine in      |
|                                     |                                   | filter adds gaussian noise equally to    |
|                                     |                                   | all fields in the state vector. It is    |
|                                     |                                   | recommended to set this option to true   |
|                                     |                                   | so code in the model_mod is called       |
|                                     |                                   | instead. This allows only a limited      |
|                                     |                                   | number of fields to be perturbed. For    |
|                                     |                                   | example, only perturbing the temperature |
|                                     |                                   | field T with a small amount of noise and |
|                                     |                                   | then running the model forward for a few |
|                                     |                                   | days is often a recommended way to       |
|                                     |                                   | generate an ensemble from a single       |
|                                     |                                   | state.                                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| fields_to_perturb                   | character(len=32), dimension(100) | If perturbing a single state to generate |
|                                     |                                   | an ensemble, set                         |
|                                     |                                   | 'custom_routine_to_generate_ensemble =   |
|                                     |                                   | .true.' and list list the field(s) to be |
|                                     |                                   | perturbed here.                          |
+-------------------------------------+-----------------------------------+------------------------------------------+
| perturbation_amplitude              | real(r8), dimension(100)          | For each field name in the               |
|                                     |                                   | 'fields_to_perturb' list give the        |
|                                     |                                   | standard deviation for the gaussian      |
|                                     |                                   | noise to add to each field being         |
|                                     |                                   | perturbed.                               |
+-------------------------------------+-----------------------------------+------------------------------------------+
| pert_base_vals                      | real(r8), dimension(100)          | If pert_sd is positive, this the list of |
|                                     |                                   | values to which the field(s) listed in   |
|                                     |                                   | pert_names will be reset if filter is    |
|                                     |                                   | told to create an ensemble from a single |
|                                     |                                   | state vector. Otherwise, it's is the     |
|                                     |                                   | list of values to use for each ensemble  |
|                                     |                                   | member when perturbing the single field  |
|                                     |                                   | named in pert_names. Unused unless       |
|                                     |                                   | pert_names is set and pert_base_vals is  |
|                                     |                                   | not the DART missing value.              |
+-------------------------------------+-----------------------------------+------------------------------------------+
| using_chemistry                     | logical                           | If using CAM-CHEM, set this to .true.    |
+-------------------------------------+-----------------------------------+------------------------------------------+
| using_variable_mean_mass            | logical                           | If using any variant of WACCM with a     |
|                                     |                                   | very high model top, set this to .true.  |
+-------------------------------------+-----------------------------------+------------------------------------------+
| debug_level                         | integer                           | Set this to increasingly larger values   |
|                                     |                                   | to print out more debugging information. |
|                                     |                                   | Note that this can be very verbose. Use  |
|                                     |                                   | with care.                               |
+-------------------------------------+-----------------------------------+------------------------------------------+

Setup Variations
================

The variants of CAM require slight changes to the setup scripts (in
``$DARTROOT/models/cam-fv/shell_scripts``) and in the namelists (in
``$DARTROOT/models/cam-fv/work/input.nml``). From the DART side, assimilations can be
started from a pre-existing ensemble, or an ensemble can be created from a
single initial file before the first assimilation. In addition, there are setup
differences between 'perfect model' runs, which are used to generate synthetic
observations, and assimilation runs. Those differences are extensive enough that
they've been coded into separate `Setup Scripts`_.

Since the CESM compset and resolution, and the initial ensemble source are
essentially independent of each other, changes for each of those may need to be
combined to perform the desired setup.

Perturbed Ensemble
------------------

The default values in ``work/input.nml`` and
``shell_scripts/CESM1_2_1_setup_pmo`` and
``shell_scripts/CESM1_2_1_setup_hybrid`` are set up for a CAM-FV, single
assimilation cycle using the default values as found in ``model_mod.f90`` and
starting from a single model state, which must be perturbed into an ensemble.
The following are suggestions for setting it up for other assimilations.
Namelist variables listed here might be in any namelist within ``input.nml``.

CAM-FV
------

If built with the FV dy-core, the number of model top levels with extra
diffusion in CAM is controlled by ``div24del2flag``. The recommended minium
values of ``highest_state_pressure_Pa`` come from that variable, and
``cutoff*vert_normalization_X``:

.. code-block:: fortran

      2    ("div2") -> 2 levels  -> highest_state_pressure_Pa =  9400. Pa
      4,24 ("del2") -> 3 levels  -> highest_state_pressure_Pa = 10500. Pa

and:

.. code-block:: fortran

      vert_coord          = 'pressure'
      state_num_1d        = 0,
      state_num_2d        = 1,
      state_num_3d        = 6,
      state_names_1d      = ''
      state_names_2d      = 'PS'
      state_names_3d      = 'T', 'US', 'VS', 'Q', 'CLDLIQ', 'CLDICE'
      which_vert_1d       = 0,
      which_vert_2d       = -1,
      which_vert_3d       = 6*1,
      highest_state_pressure_Pa = 9400. or 10500. 

WACCM
-----

WACCM[#][-X] has a much higher top than the CAM versions, which requires the use
of scale height as the vertical coordinate, instead of pressure, during
assimilation. One impact of the high top is that the number of top model levels
with extra diffusion in the FV version is different than in the low-topped
CAM-FV, so the ``div24del2flag`` options lead to the following minimum values
for ``highest_state_pressure_Pa``:

.. code-block:: fortran

   2    ("div2") -> 3 levels  -> highest_state_pressure_Pa = 0.01 Pa
   4,24 ("del2") -> 4 levels  -> highest_state_pressure_Pa = 0.02 Pa

The best choices of ``vert_normalization_scale_height``, ``cutoff``, and
``highest_state_pressure_Pa`` are still being investigated (April, 2015), and
may depend on the observation distribution being assimilated.

WACCM is also typically run with coarser horizontal resolution. There's an
existing 2-degree ensemble, so see the `Continuing after the first cycle`_
section to start from it, instead of a single state. If you use this, ignore any
existing inflation restart file and tell DART to make its own in the first cycle
in ``input.nml``:

.. code-block:: fortran

   inf_initial_from_restart    = .false.,                 .false.,
   inf_sd_initial_from_restart = .false.,                 .false.,

In any case, make the following changes (or similar) to convert from a CAM setup
to a WACCM setup. ``CESM1_2_1_setup_hybrid``:

.. code-block::

   setenv compset     F_2000_WACCM
   setenv resolution  f19_f19  
   setenv refcase     FV1.9x2.5_WACCM4
   setenv refyear     2008
   setenv refmon      12
   setenv refday      20

and the settings within ``input.nml``:

.. code-block::

   vert_normalization_scale_height = 2.5
   vert_coord                = 'log_invP'
   highest_obs_pressure_Pa   = .001,
   highest_state_pressure_Pa = .01,

If built with the SE dy-core (warning; experimental), then 4 levels will have
extra diffusion, and also see the `CAM-SE`_ section.

If there are problems with instability in the WACCM foreasts, try changing some
of the following parameters in either the user_nl_cam section of the setup
script or input.nml.

-  The default div24del2flag in WACCM is 4. Change it in the setup script to

   .. code-block::

      echo " div24del2flag         = 2 "                       >> ${fname}

   which will use the ``cd_core.F90`` in SourceMods, which has doubled diffusion
   in the top layers compared to CAM.

-  Use a smaller dtime (1800 s is the default for 2-degree) in the setup script.
   This can also be changed in the ensemble of ``user_nl_cam_####`` in the
   ``$CASEROOT`` directory.

   .. code-block::

      echo " dtime         = 600 "                             >> ${fname}

-  Increase highest_state_pressure_Pa in input.nml:

   .. code-block::

      div24del2flag = 2    ("div2") -> highest_state_pressure_Pa = 0.1 Pa
      div24del2flag = 4,24 ("del2") -> highest_state_pressure_Pa = 0.2 Pa

-  Use a larger nsplit and/or nspltvrm in the setup script:

   .. code-block::

      echo " nsplit         = 16 "                             >> ${fname}
      echo " nspltvrm       =  4 "                             >> ${fname}

-  Reduce ``inf_damping`` from the default value of ``0.9`` in ``input.nml``:

   .. code-block:: fortran

      inf_damping           = 0.5,                   0,

Notes for Continuing an Integration
===================================

Continuing after the first cycle
--------------------------------

After the first forecast+assimilation cycle, using an ensemble created from a
single file, it is necessary to change to the 'continuing' mode, where CAM will
not perform all of its startup procedures and DART will use the most recent
ensemble. This example applies to an assimiation using prior inflation
(``inf_...= .true.``). If posterior inflation were needed, then the 2nd column
of ``infl_...`` would be set to ``.true..``. Here is an example snippet from
``input.nml``:

.. code-block:: fortran
   
      start_from_restart      = .true.,
      restart_in_file_name    = "filter_ics",
      single_restart_file_in  = .false.,

      inf_initial_from_restart    = .true.,                 .false.,
      inf_sd_initial_from_restart = .true.,                 .false.,

Combining multiple cycles into one job
--------------------------------------

``CESM1_2_1_setup_hybrid`` and ``CESM1_2_1_setup_pmo`` are set up in the default
cycling mode, where each submitted job performs one model advance and one
assimilation, then resubmits the next cycle as a new job. For long series of
cycles, this can result in a lot of time waiting in the queue for short jobs to
run. This can be prevented by using the 'cycles' scripts generated by
``CESM1_2_1_setup_advanced`` (instead of ``CESM1_2_1_setup_hybrid``). This mode
is described in ``$DARTROOT/models/cam-fv/doc/README_cam-fv``.

Discussion
==========

Many CAM initial file variables are already handled in the ``model_mod``. Here
is a list of others, which may be used in the future. Each would need to have a
DART ``*KIND*`` associated with it in ``model_mod``.

.. code-block:: fortran

   Atmos
      CLOUD:       "Cloud fraction" ;
      QCWAT:       "q associated with cloud water" ;
      TCWAT:       "T associated with cloud water" ;
      CWAT:        "Total Grid box averaged Condensate Amount (liquid + ice)" ;
      also? LCWAT

   pbl
      PBLH:        "PBL height" ;
      QPERT:       "Perturbation specific humidity (eddies in PBL)" ;
      TPERT:       "Perturbation temperature (eddies in PBL)" ;

   Surface
      LANDFRAC:    "Fraction of sfc area covered by land" ;
      LANDM:       "Land ocean transition mask: ocean (0), continent (1), transition (0-1)" ;
         also LANDM_COSLAT
      ICEFRAC:     "Fraction of sfc area covered by sea-ice" ;
      SGH:         "Standard deviation of orography" ;
      Z0FAC:       "factor relating z0 to sdv of orography" ;
      TS:          "Surface temperature (radiative)" ;
      TSOCN:       "Ocean tempertare" ;
      TSICE:       "Ice temperature" ;
      TSICERAD:    "Radiatively equivalent ice temperature" ;

   Land/under surface
      SICTHK:      "Sea ice thickness" ;
      SNOWHICE:    "Water equivalent snow depth" ;
      TS1:         "subsoil temperature" ;
      TS2:         "subsoil temperature" ;
      TS3:         "subsoil temperature" ;
      TS4:         "subsoil temperature" ;

   Other fields are not included because they look more CLM oriented.

   Other fields which users may add to the CAM initial files are not listed here.


Files
=====

-  ``model_nml`` in ``input.nml``
-  ``cam_phis.nc`` (CAM surface height file, often CAM's .h0. file in the CESM run
   environment)
-  ``caminput.nc`` (CAM initial file)
-  ``clminput.nc`` (CLM restart file)
-  ``iceinput.nc`` (CICE restart file) by model_mod at the start of each
   assimilation)
-  netCDF output state diagnostics files

Nitty gritty: Efficiency possibilities
======================================

-  index_from_grid (and others?) could be more efficient by calculating and
   globally storing the beginning index of each cfld and/or the size of each
   cfld. Get_state_meta_data too. See ``clm/model_mod.f90``.

-  Global storage of height fields? but need them on staggered grids (only
   sometimes) Probably not; machines going to smaller memory and more
   recalculation.

-  ! Some compilers can't handle passing a section of an array to a
   subroutine/function; I do this in ``nc_write_model_vars(?)`` and/or
   ``write_cam_init(?)``; replace with an exactly sized array?

-  Is the testing of resolution in read_cam_coord overkill in the line that
   checks the size of ``(resol_n - resol_1)*resol``?

-  Replace some do loops with forall (constructs)

-  Subroutine ``write_cam_times(model_time, adv_time)`` is not needed in
   CESM+DART framework? Keep anyway?

-  Remove the code that accommodates old CAM coordinate order (``lon,lev,lat``).

-  Cubed sphere: Convert lon,lat refs into dim1,dim2 in more subroutines.
   get_val_heights is called with (``column_ind,1``) by CAM-SE code, and
   (``lon_ind, lat_ind``) otherwise).

-  ``cam_to_dart_kinds`` and ``dart_to_cam_types`` are dimensioned 300,
   regardless of the number of fields in the state vector and/or *KIND*\ s .

-  Describe:

   - The coordinate orders and translations; CAM initial file, ``model_mod``,
     and ``DART_Diag.nc``.
   - Motivations

     - There need to be 2 sets of arrays for dimensions and dimids;
   
       - one describing the caminput file (``f_...``)
       - and one for the state (``s_...``) (storage in this module).
       - Call them ``f_dim_Nd``, ``f_dimid_Nd``
       - ``s_dim_Nd``, ``s_dimid_Nd``      

-  Change (private only) subroutine argument lists; structures first, regardless
   of in/out then output, and input variables.

-  Change declarations to have dummy argument integers used as dimensions first

-  Implement a ``grid_2d_type``? Convert phis to a ``grid_2d_type``? ps, and
   staggered ps fields could also be this type.

-  Deallocate ``grid_1d_arrays`` using ``end_1d_grid_instance`` in end_model.
   ``end_model`` is called by subroutines ``pert_model_state``,
   ``nc_write_model_vars``; any problem?

-  ISSUE; In ``P[oste]rior_Diag.nc`` ensemble members are written out
   \*between\* the field mean/spread pair and the inflation mean/sd pair. Would
   it make more sense to put members after both pairs? Easy to do?

-  ISSUE?; ``model_interpolate`` assumes that obs with a vertical location have
   2 horizontal locations too. The state vector may have fields for which this
   isn't true, but no obs we've seen so far violate this assumption. It would
   have to be a synthetic/perfect_model obs, like some sort of average or
   parameter value.

-  ISSUE; In convert_vert, if a 2D field has dimensions (lev, lat) then how is
   ``p_surf`` defined? Code would be needed to set the missing dimension to 1,
   or make different calls to ``coord_ind``, etc.

-  ISSUE; The ``QTY_`` list from obs_def_mod must be updated when new fields are
   added to state vector. This could be done by the preprocessor when it inserts
   the code bits corresponding to the lists of observation types, but it
   currently (10/06) does not. Document accordingly.

-  ISSUE: The CCM code (and Hui's packaging) for geopotentials and heights use
   different values of the physical constants than DART's. In one case Shea
   changed g from 9.81 to 9.80616, to get agreement with CCM(?...), so it may be
   important. Also, matching with Hui's tests may require using his values;
   change to DART after verifying?

-  ISSUE: It's possible to figure out the model_version from the NetCDF file
   itself, rather than have that be user-provided (sometimes incorrect and hard
   to debug) meta-data. model_version is also misnamed; it's really the
   ``caminput.nc`` model version. The actual model might be a different
   version(?). The problem with removing it from the namelist is that the
   scripts need it too, so some rewriting there would be needed.

-  ISSUE: ``max_neighbors`` is set to 6, but could be set to 4 for non-refined
   grids. Is there a good mechanism for this? Is it worth the file space
   savings?

-  ISSUE: ``x_planar`` and ``y_planar`` could be reduced in rank, if no longer
   needed for testing and debugging.

-  "Pobs" marks changes for providing expected obs of P break from past
   philosophy; P is not a native CAM variable (but is already calced here)

-  NOVERT marks modifications for fields with no vertical location, i.e. GWD
   parameters.

Terms of Use
============
 
|Copyright| University Corporation for Atmospheric Research
 
Licensed under the `Apache License, Version 2.0
<http://www.apache.org/licenses/LICENSE-2.0>`__. Unless required by applicable
law or agreed to in writing, software distributed under this license is
distributed on an "as is" basis, without warranties or conditions of any kind,
either express or implied.
 
.. |Copyright| unicode:: 0xA9 .. copyright sign

References and Acknowledgements
===============================

-  `CAM homepage <http://www.ccsm.ucar.edu/models/atm-cam/>`__

Ave Arellano did the first work with CAM-Chem, assimilating MOPPITT CO
observations into CAM-Chem. Jerome Barre and Benjamin Gaubert took up the
development work from Ave, and prompted several additions to DART, as well as
``model_mod.f90``.

Nick Pedatella developed the first ``vert_coord = 'log_invP'`` capability to
enable assimilation using WACCM and scale height vertical locations.
