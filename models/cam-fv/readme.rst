Atmospheric Models in CESM 
==========================

Overview
--------

The larger context of the Community Earth System Model and DART interactions
is described in the `CESM readme <../CESM/readme.html>`_
This document focuses on the several `atmospheric models <http://www2.cesm.ucar.edu/models>`__
that have been developed or adapted to run in the CESM environment. 
They are named according to their dynamical core ("dycore").
As of 2021 these include Finite Volume Community Atmosphere Model (CAM-FV), 
Spectral Element (CAM-SE), and MPAS.
The DART system has supported data assimilation into CAM-FV continuously for many years.
An interface to CAM-SE was added to DART in 2022.
An interface to MPAS is being developed (contact us about the current status).

.. |CAM6_Rean| replace:: 1 degree reanalysis wiki
.. _CAM6_Rean: https://github.com/NCAR/DART/wiki/1-degree,-CAM6,-ensemble-reanalysis-for-CESM-experiments-(2011-thru-2019):-DATM,-hindcasts,-model-evaluation

The flexibility of the DART environment has led to its use
by graduate students, post-graduates, and scientists at universities and
research labs to conduct data assimilation research. Others are using the
products of data assimilation (analyses), which were produced at NSF NCAR
using CESM+DART, to conduct related research. 
The latest reanalysis is described in the DART |CAM6_Rean|_

The variety of research can be sampled on the DART  
`Publications <https://dart.ucar.edu/pages/Publications.html>`__ page.

Terminology
~~~~~~~~~~~

The atmospheric component used in CESM is built
with two independent main characteristics. CESM labels these as:

**resolution**

   signifies both the horizontal resolution of the grid
   (not the vertical) **and** the dynamical core run on the specified grid.
   The dynamical core refers to the fluid dynamical equations
   run on the specified grid.
   Examples of resolution (short) names are f19_f19 (~2 degree Finite Volume dycore)
   or ne30np4_gx1v6 (~1 degree Spectral Element dycore).

**compset**

   refers to the vertical grid **and** the parameterizations --
   the formulation of the subgridscale physics -- as well as the combination
   of active, data, or stub model components. These parameterizations 
   consist of the equations describing physical processes such as convection,
   radiation, and chemistry.
   
   - The vertical grid is determined by the needs of the chosen
     parameterizations, thus the vertical spacing and the top level of the
     model domain vary with those choices. 
   - The combinations of parameterizations and vertical grids are named: CAM3.5,
     CAM5, CAM#, ... WACCM, WACCM#, WACCM-X, CAM-Chem.  
   - The compset is specified as described in the `CESM readme <../CESM/readme.html>`_.

**ensemble, multi-instance, and multidriver**

   These are essentially synonyms referring to multiple, closely related models 
   or model states.  "Ensemble" is DART's vocabulary, while "multi-instance"
   is CESM's original term for an ensemble.  
   "Multidriver" is replacing "multi-instance", and refers to the CESM module
   which coordinates the running of all of the model components.
   Similarly, DART ensemble "members" are the same as CESM "instances".

**CASE, CASEROOT**

   Running a DART setup script creates a CESM "CASE" (the name of your experiment)
   in the "CASEROOT" directory (from where jobs will be controlled and launched),
   both of which are defined in the setup script.  There will also be a run directory
   named $CASEROOT in your scratch space (usually), where the fortran executables 
   can also be found ($scratch/$CASEROOT/bld).
   

`Setup Scripts`_ describes how to specify these and other choices 
in the assimilation setup scripts and namelists.

CAM-FV
------

Here are some highlighted features of this DART interface to CAM-FV.

-  Assimilate within the CESM software framework by using the multidriver
   capability of CESM2 (and later). This enables assimilation of suitable
   observations into a variety of CESM components and leverages CESM's
   build, run, and archiving capabilities.
-  Use any horizontal and vertical resolution of CAM-FV.
-  Assimilate a variety of observations.  To date the observations successfully
   assimilated include: 

   * NCEP reanalysis BUFR obs (T,U,V,Q), 
   * Global Positioning System radio occultation observations (refractivity and electron density), 
   * AIRS retrievals (T and Q), 
   * MOPITT (carbon monoxide, when a chemistry model is incorporated into CAM-FV),
   * OCO2 (carbon dioxide), 
   * Aura MLS (T),
   * SABER (T), 
   * GNSS (total electron content, "TEC"),
   * The development of the ability to assimilate RTTOV radiances is nearly complete
     (2021: contact us for the current status).
   * Research has also explored
     assimilating surface observations, cloud liquid water, and aerosols. 

   The Aura MLS, SABER, and GNSS observations have been assimilated into WACCM
   and WACCM-X; "high top" versions of CAM-FV.
-  Specify, via namelist entries, the CAM (initial file) variables which will be
   directly affected by the observations, that is, the state vector. 
-  Generate analyses on the CAM grid which have only CAM model error in them,
   rather than another model's.

Reanalyses
~~~~~~~~~~

There have been two large-scale reanalysis efforts using CAM-FV and DART. 
The **CAM6 Data Assimilation Research Testbed (DART) Reanalysis**
is archived in the NSF NCAR Research Data Archive 
`DS345.0 <https://rda.ucar.edu/datasets/ds345.0/#!description>`__ .
(See the |CAM6_Rean|_ ).
It contains just under 120Tb (yes Tb) of data:

   These CAM6+DART Reanalysis data 
   products are designed to facilitate a broad variety of research using 
   NSF NCAR's CESM2 models, ranging from model evaluation to (ensemble) 
   hindcasting (initial conditions), data assimilation experiments, and sensitivity studies. 
   They come from an 80 member ensemble reanalysis of the global 
   troposphere and stratosphere using CAM6-FV from CESM2.1. 
   The data products represent the actual states of the atmosphere 
   from 2011-2019 at a ~1 degree horizontal resolution and up to 
   6 hourly frequency. Each ensemble member is an equally likely 
   description of the atmosphere, and is also consistent with 
   dynamics and physics of CAM6-FV.
   
   
An earlier, more limited dataset can be found in the 
`**Ensemble of Atmospheric Forcing Files from a CAM4-FV Reanalysis** 
<https://github.com/NCAR/DART/wiki/2-degree-DATM-ensemble-for-CESM-experiments-(1998-thru-2010)>`__
is archived in the NSF NCAR Research Data Archive 
`DS199.1 <https://rda.ucar.edu/datasets/ds199.1/>`__ .
It contains about 1.5Tb of data:

   This dataset contains files that are an ensemble of 'coupler history' 
   files from an 80-member reanalysis performed with the Data Assimilation 
   Research Testbed (DART) using the Community Atmosphere Model Version 
   4 with the finite volume core (CAM4-FV) at 1.9 degree by 2.5 degree 
   resolution. The observations assimilated include all those used in 
   the NCEP/NSF NCAR reanalysis (temperature and wind components from 
   radiosondes, aircraft, and satellite drift winds) plus radio 
   occultation observations from the COSMIC satellites starting in late 
   2006. These files are intended to be used as 'DATM stream files' 
   for CESM component sets that require a data atmosphere. Some example 
   stream text files are included in the RDA to illustrate how to use these data.

..  * CAM4, 2 degree, 2000-2010, `Reanalysis <https://rda.ucar.edu/datasets/ds199.1>`__
..  * files from the old "large file site": http://www.image.ucar.edu/pub/DART/CAM/

Observations
------------

The CAM6+DART Reanalysis used "observation sequence files"
which contain the types of observations in the table below
("T" = temperature, "U" = zonal wind, "V" = meridional wind,
"Q" = specific humidity, "refractivity" = the bending of light by density variations).
These files are available on NSF NCAR's glade file system:
/glade/p/cisl/dares/Observations/NCEP+ACARS+GPS+AIRS/Thinned_x9x10.
Versions of these files, which also have the results of the reanalysis in them,
are available from the RDA ds345.0 linked above.

NCEP
   NCEP's PREPBUFR files (prepqm) in NSF NCAR's Research Data Archive:
   (https://rda.ucar.edu/datasets/ds090.0/)

COSMIC
   This site (http://www.cosmic.ucar.edu/) provides atmospheric refractivity 
   from a variety of satellites (including COSMIC), which receive Global Positioning System 
   radio occultation signals.

AIRS
   Retrievals from `infrared soundings <http://airs.jpl.nasa.gov/>`_  
   from the `AQUA satellite <http://aqua.nasa.gov/>`_
   They are thinned by a factor of 90 to make their density comparable to the radiosonde network.

+----------------------+----------------------------------------+--------------------------+--------+
| Observation or       |                                        |                          | Data   |
| Retrieval            | Platform                               | Distribution             | Source | 
+======================+========================================+==========================+========+
| T, U, V, Q           | Radiosondes from balloons              | mostly land              | NCEP   |
+----------------------+----------------------------------------+--------------------------+--------+
| T, U, V              | ACARS commercial aircraft              | mostly North America     | NCEP   |
+----------------------+----------------------------------------+--------------------------+--------+
| T, U, V              | AIRCRAFT commercial aircraft           | mostly non-North America | NCEP   |
+----------------------+----------------------------------------+--------------------------+--------+
| U, V                 | Cloud drift winds from GOES satellites | midlatitudes and tropics | NCEP   |
+----------------------+----------------------------------------+--------------------------+--------+
| index of refraction  | Global Positioning System receivers    | global                   | COSMIC |
+----------------------+----------------------------------------+--------------------------+--------+
| T, Q                 | AQUA satellite; AIRS instrument        | global                   | AIRS   |
+----------------------+----------------------------------------+--------------------------+--------+
| altimeter            | Radiosondes, bouys                     | global surface           | NCEP   |
+----------------------+----------------------------------------+--------------------------+--------+

Sample sets of observations, which can be used with CAM+DART assimilations, can
be found at http://www.image.ucar.edu/pub/DART/Obs_sets/ of which the NCEP BUFR
observations are the most widely used.

The CAM-FV DART Interface
=========================

The 19 public interface subroutines in ``model_mod.f90`` are standardized for all DART
compliant models. These interfaces allow DART to get the model state and
metadata describing this state, find state variables that are close to a given
location, and do spatial interpolation for a variety of variables required by
observational operators.
Your choices for how the assimilation (not the hindcast) will happen 
are defined in the ``cam-fv/work/input.nml`` file. 
In that file, the ``model_nml`` namelist lets you control the interaction with CAM-FV.
The CAM-FV, which DART will interact with, is defined by the setup scripts,
as described next.

.. _`Setup Scripts`:

Setup Scripts
-------------

Unlike pre-Manhattan versions of DART-CAM, CESM (CAM) runs using its normal scripts, 
then stops and calls a DART script, which does the desired assimilation tasks, 
then returns to the CESM run script for the next model advance. See the CESM
interface documentation in the `CESM readme <../CESM/readme.html>`_
for more general information about
running DART with CESM. Due to the complexity of the CESM software environment,
the versions of CESM which can be used for assimilation are more restricted than
previously. Each supported CESM version has similar, but unique, sets of setup
scripts and CESM `SourceMods`_. Those generally do not affect the
``cam-fv/model_mod.f90`` interface. 

The primary purpose of a setup script is to *set up* a CESM "CASE"
(compset, resolution, etc.), which can be used by DART.
The ability to *use* DART programs is then set up by a second script; ``DART_config``,
which was created by the setup script.
Here is an outline of the scripts, which are currently (2021) in shell_scripts.
They are roughly in order of complexity, which is the order in which
you might want to use them.
The indenting shows which scripts are used by, or associated with, another script.

.. FIXME; code-block with no argument colors random words in the descriptions.
   What's a better format that code-block?  Table doesn't work well because of indenting limits.

.. code-block::

   cesm2_1/                              Directory containing scripts developed for CESM2_1
       spinup_single                     Setup a single instance (member) CAM-FV case to advance a model state 
                                         some months to a desired date.
       setup_hybrid                      Basic script to set up an assimilation case.
          DART_config.template           Modified to create the script which modifies a CESM CASE to do assimilation.
             no_assimilate.csh.template  Modified to create a script which does no assimilation, 
                                         but prepares files for the next model advance.
             assimilate.csh.template     Modified to create the assimilate.csh script
       compress.csh                      Example of compressing assimilation output for efficient archiving.
                                         Can be called by assimilate.csh
       mv_to_campaign.csh                Example of how to use globus to move files to a remote archive.
       setup_advanced                    Like setup_hybrid, but more model and assimilation features can be modified.
                                         It modifies DART_config.template like setup_hybrid does.
       setup_single_from_ens             Set up a single-instance run using initial conditions taken from 
                                         a single instance of a multi-instance CAM hindcast.  Useful for debugging.
       standalone.pbs                    Batch job tests of assimilation with no model advances.
          test_assimilate.csh            A simpler (earlier) form of assimilate.csh.  
    cesm2_0/                             Similar contents to cesm2_1, plus the following.
       obs_seq_tool_series.csh           Script to process a series of obs_seq.final files,
                                         to change any of the properties available to obs_sequence_tool.
       spinup_single_sst.25              Same as cesm2_1/spinup_single, but uses a high resolution SST dataset.
    synth_obs_locs_to_seqs.csh           Take text output from, e.g. even_sphere.m, and create obs_seq.in files
                                         for use in perfect_model_obs.
 
The scripts  in cesm#_# will handle, for that CESM version;

   * all CAM-FV "physics" variants and vertical resolutions.
     For example, CAM5.5, CAM6, ..., WACCM4, WACCM6, WACCM-X, ..., CAM-Chem.
   * all horizontal resolutions of CAM-FV; 1.9x2.5 (f19xf19), 0.9x1.25 (f09xf09), ....

Physics variants of other dycores are handled in other "model" interfaces,
such as models/cam-se.

.. _reanalysis: https://github.com/kdraeder/cesm}{github.com/kdraeder/cesm

.. _SourceMods: 

SourceMods
~~~~~~~~~~

The most recent SourceMods for the CAM6+DART interface can be fetched from
the github `reanalysis`_ repository.
Change to the cesm2_1_forcing_rean branch, which includes a SourceMods tar file.
Unpack that file into the location you specify in the setup script, before building the CASE.

.. The latest (2020) SourceMods are 
   /glade/u/home/raeder/cesm2_1_relsd_m5.6/CAM6+DART_Reanalysis_SourceMods.tgz
   ? Where should these live?
   ? Are there any SourceMods for WACCM(-X)?
   
Namelists
---------

DART assembles the namelists for all of the relevant modules
into a single namelist file; ``models/cam-fv/input.nml``.
This section focuses on ``model_nml``,
but others are referenced, as needed.
Namelists start with an ampersand ``&`` and terminate with a slash ``/``. 
Character strings that contain a ``/`` must be enclosed in quotes to prevent them 
from prematurely terminating the namelist.
Text outside of the &.../ pairs is ignored.

Here's a list of the model_nml variables and default values.
More detailed descriptions follow in a table and subsections.

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

+-------------------------------------+----------------+-------------------------------------------+
| Item                                | Type           | Description                               |
+=====================================+================+===========================================+
| cam_template_filename               | character      | CAM initial file used to provide          |
|                                     | (len=128)      | configuration information, such as the    |
|                                     |                | grid resolution, number of vertical       |
|                                     |                | levels, whether fields are staggered or   |
|                                     |                | not, etc.  Created by the first hindcast. |
+-------------------------------------+----------------+-------------------------------------------+
| cam_phis_filename                   | character      | CAM topography file. Reads the "PHIS"     |
|                                     | (len=128)      | NetCDF variable from this file.           |
|                                     |                | Typically this is a CAM History file      |
|                                     |                | because this field is not normally found  |
|                                     |                | in a CAM initial file. Created by the     |
|                                     |                | first hindcast.                           |
+-------------------------------------+----------------+-------------------------------------------+
| vertical_localization_coord         | character      | The vertical coordinate to which all      |
|                                     | (len=128)      | vertical locations are converted in       |
|                                     |                | model_mod. Valid options are "pressure",  |
|                                     |                | "height", "scaleheight" or "level".       |
+-------------------------------------+----------------+-------------------------------------------+
| use_log_vertical_scale              | logical        | Use the log of the vertical distances     |
|                                     |                | when interpolating.  This is only used    |
|                                     |                | for locations having which_vert =         |
|                                     |                | VERTISPRESSURE. It should be .true. when  |
|                                     |                | vertical_localization_coord =             |
|                                     |                | "scaleheight" or "height".                |
+-------------------------------------+----------------+-------------------------------------------+
| no_normalization_of_scale_heights   | logical        | If true (default), scale height is        |
|                                     |                | computed as the log of the pressure at    |
|                                     |                | the given location.                       |
|                                     |                | Beware: unnormalized scale heights        |
|                                     |                | decrease upward, and may have values < 0. |
|                                     |                | This works because only differences       |
|                                     |                | of scale height are used and              |
|                                     |                | find_enclosing_indices assigns the larger |
|                                     |                | and smaller coordinate values correctly   |
|                                     |                | in the interpolation.                     |
|                                     |                | If false, the scale height is computed    |
|                                     |                | as the log of the ratio of the surface    |
|                                     |                | pressure to the pressure aloft.           |
|                                     |                | In previous versions normalization        |
|                                     |                | was the default.  It is slightly less     |
|                                     |                | efficient.                                |
+-------------------------------------+----------------+-------------------------------------------+
| no_obs_assim_above_level            | integer        | Because the top of the model is highly    |
|                                     |                | damped it is recommended to NOT           |
|                                     |                | assimilate observations in the top model  |
|                                     |                | levels. The units here are CAM model      |
|                                     |                | level numbers. Set it to equal or below   |
|                                     |                | the lowest model level (the highest       |
|                                     |                | number) where damping is applied in the   |
|                                     |                | model.   See `Diffusion`_\ , below.       |
+-------------------------------------+----------------+-------------------------------------------+
| model_damping_ends_at_level         | integer        | Set this to the lowest model level (the   |
|                                     |                | highest number) where model damping is    |
|                                     |                | applied. Observations below the           |
|                                     |                | 'no_obs_assim_above_level' cutoff, but    |
|                                     |                | close enough to the model top to have an  |
|                                     |                | impact during the assimilation, will have |
|                                     |                | their impacts decreased smoothly to 0 at  |
|                                     |                | this given model level. The assimilation  |
|                                     |                | should make no changes to the model       |
|                                     |                | state above the given level.              |
|                                     |                | See `Diffusion`_\ , below.                |
+-------------------------------------+----------------+-------------------------------------------+
| state_variables                     | character      | Character string table that includes:     |
|                                     | (len=64)       | 1. CAM initial file variable names of     |
|                                     | dimension(100) | fields to be read into the state vector,  |
|                                     |                | 2. the corresponding DART QTY (quantity)  |
|                                     |                | 3. if a bounded quantity, the minimum and |
|                                     |                | maximum valid values,                     |
|                                     |                | 4. the string 'UPDATE' indicates that     |
|                                     |                | the updated values should be written      |
|                                     |                | back to the output file. 'NOUPDATE' will  |
|                                     |                | skip writing this field at the end of     |
|                                     |                | the assimilation.                         |
|                                     |                | See `State Variables`_\ , below.          |
+-------------------------------------+----------------+-------------------------------------------+
| assimilation_period_days            | integer        | With assimilation_period_seconds,         |
|                                     |                | sets the assimilation cycle length.       |
|                                     |                | They should match the model advance time. |
|                                     |                | The CAM scripts distributed with          |
|                                     |                | DART set these to 0 days, 21600 seconds   |
|                                     |                | (6 hours).                                |
|                                     |                | They also set the assimilation window     |
|                                     |                | width.                                    |
+-------------------------------------+----------------+-------------------------------------------+
| assimilation_period_seconds         | integer        | See assimilation_period_days              |
+-------------------------------------+----------------+-------------------------------------------+
| suppress_grid_info_in_output        | logical        | Filter can update fields in existing      |
|                                     |                | files or create diagnostic/output files   |
|                                     |                | from scratch. By default files created    |
|                                     |                | from scratch include a full set of CAM    |
|                                     |                | grid information to make the file fully   |
|                                     |                | self-contained and plottable. However,    |
|                                     |                | to save disk space the grid variables     |
|                                     |                | can be suppressed in files created by     |
|                                     |                | filter by setting this to true.           |
+-------------------------------------+----------------+-------------------------------------------+
| custom_routine_to_generate_ensemble | logical        | Use the subroutines in model_mod.f90      |
|                                     |                | to create an ensemble of initial          |
|                                     |                | conditions (with non-0 spread) from a     |
|                                     |                | single CAM initial file.  This is useful  |
|                                     |                | when there is no existing ensemble of     |
|                                     |                | ICs.  See `Perturbed`_\ , below.          |
+-------------------------------------+----------------+-------------------------------------------+
| fields_to_perturb                   | character,     | If perturbing a single state to generate  |
|                                     | (len=32)       | an ensemble, set                          |
|                                     | dimension(100) | 'custom_routine_to_generate_ensemble =    |
|                                     |                | .true.' and list here the DART QTYs of    |
|                                     |                | the field(s) to be perturbed.             |
+-------------------------------------+----------------+-------------------------------------------+
| perturbation_amplitude              | real(r8),      | For each field name in the                |
|                                     | dimension(100) | 'fields_to_perturb' list, give the        |
|                                     |                | standard deviation of the gaussian noise  |
|                                     |                | to add to each field being perturbed.     |
+-------------------------------------+----------------+-------------------------------------------+
| using_chemistry                     | logical        | If using CAM-CHEM, set this to .true.     |
+-------------------------------------+----------------+-------------------------------------------+
| using_variable_mean_mass            | logical        | If using any variant of WACCM (a very     |
|                                     |                | high model top), set this to .true.       |
+-------------------------------------+----------------+-------------------------------------------+
| debug_level                         | integer        | Set this to increasingly larger values    |
|                                     |                | to print out more debugging information.  |
|                                     |                | Note that this can be very verbose. Use   |
|                                     |                | with care.                                |
+-------------------------------------+----------------+-------------------------------------------+

.. _`Setup Variations`:

Setup Variations
----------------

The default values in ``cam-fv/shell_scripts/cesm#_#/setup*`` 
and in the namelists in ``cam-fv/work/input.nml``
are (mostly) set up for a single assimilation cycle of CAM-fV, 
starting from a single model state, which must be perturbed into an ensemble.
The following are suggestions for setting it up for other assimilations.
Namelist variables listed here might be in any namelist within ``input.nml``.

.. _`State variables`:

State Variables
~~~~~~~~~~~~~~~

This implementation of the DART interface module for the CAM and WACCM models
uses the CAM initial files (**not** restart files) for transferring the model
state to and from the ``filter``. 

The DART state vector should include all prognostic variables in the CAM
initial files which cannot be calculated directly from other prognostic
variables. In practice the state vector sometimes contains derived quantities to
enable DART to compute forward operators (expected observation values) efficiently.
The derived quantities are often overwritten when the model runs
the next timestep, so the work DART does to update them is wasted work.
The standard state vector contains the following fields,
as entered into the ``input.nml:model_nml`` namelist.

.. code-block:: fortran
   
   state_variables  = 
         'T',     'QTY_TEMPERATURE',         'NA', 'NA', 'UPDATE'
         'US',    'QTY_U_WIND_COMPONENT',    'NA', 'NA', 'UPDATE'
         'VS',    'QTY_V_WIND_COMPONENT',    'NA', 'NA', 'UPDATE'
         'Q',     'QTY_SPECIFIC_HUMIDITY',   'NA', 'NA', 'UPDATE'
         'CLDLIQ','QTY_CLOUD_LIQUID_WATER',  'NA', 'NA', 'UPDATE'
         'CLDICE','QTY_CLOUD_ICE',           'NA', 'NA', 'UPDATE'
         'PS',    'QTY_SURFACE_PRESSURE',    'NA', 'NA', 'UPDATE'

Any tracers or chemicals ("constituents" in CESM's vocabulary), 
which are needed for a given study and exist in the initial files, 
can be added to ``state_variables``.  
See the list for CAM6, below.
CAM6 variables which are *not* in the initial file can be added to it
if they are in CAM's list of constituents (or "tracers").
Those variables are identified by a ``&IC`` suffix in the "MASTER FIELD LIST"
in an "atm.log..." or "atm_0001.log..." file.
Finally (you're deeply into the weeds here), variables can be added 
to the list of constituents using CAM's ``cnst_add`` function,
which will not be described here.
In all of these cases, minor modifications to ``model_mod.f90`` and CAM may be necessary.

Here is a list of CAM initial file variables, excluding the variables listed as parts
of the most common state vector, above.
Each would need to have a DART ``*QTY*`` associated with it.

Other moisture variables 

  * NUMICE  "cloud ice number  "
  * NUMLIQ  "cloud liquid number  "
  * NUMRAI  "rain number  "
  * NUMSNO  "snow number  "
  * RAINQM  "rain amount  "
  * SNOWQM  "snow amount  "

Aerosols 

  * DMS   "dimethyl sulfide   "
  * H2O2  "H\ :sub:`2`\ O\ :sub:`2`"
  * H2SO4 "H\ :sub:`2`\ SO\ :sub:`4`"
  * SO2   "SO\ :sub:`2`"
  * SOAG  "secondary organic aerosols gas  "

MAM4 modal aerosol scheme variables ("[ ]" means use a single digit.) 

  * bc_a[1,4]   "black carbon, modes 1 and 4  "
  * dst_a[1-3]  "dust, modes 1 through 3"
  * ncl_a[1-3]  "sea salt (NaCl) , modes 1 through 3"
  * num_a[1-4]  "aerosol number density, modes 1 through 4"
  * pom_a[1,4]  "primary-organic aerosols, modes 1 and 4"
  * soa_a[1,2]  "secondary-organic aerosols, modes 1 and 2"
  * so4_a[1-3]  "sulfate (SO\ :sub:`4`) modes 1 through 3"
   
Expected observation values on pressure, scale height, height or model levels
can be requested from ``model_interpolate``. Surface observations can not yet be
interpolated, due to the difference between the model's lowest level (~7 hPa above
the model surface) and the Earth's surface where the observations are made. 
Model_interpolate can be queried for any (non-surface) variable in the state vector 
(which are variables native to CAM) plus pressure on height levels. 

The reasons initial files are used instead of restart files include:

#. The contents of the restart files vary depending on both the model release
   version and the physics packages selected.
#. There is no metadata describing the variables in the restart files. Some
   information can be tracked down in the ``atm.log`` file, but not all of it.
#. The restart files (for non-chemistry model versions) are much larger than
   the initial files (and we need to deal with an ensemble of them).
#. The temperature on the restart files is virtual equivalent potential
   temperature, which requires (at least) surface pressure, specific humidity,
   and sensible temperature to calculate.
#. CAM does not call the initialization routines when a hindcast is started 
   in ''restart'' mode, so fields which are not modified by DART 
   may be inconsistent with fields which are.
#. If DART modifies the contents of the ``.r.`` restart file, it might also
   need to modify the contents of the ``.rs.`` restart file, which has similar
   characteristics (1-3 above) to the ``.r.`` file.
#. There is no need for exact restart performance because filter alters the model state,
   making exact restarts impossible.

Inflation
~~~~~~~~~

Assimilation using CAM and WACCM should generally use one of DART's
adaptive inflation algorithms.  
As of 2021 these are ``inf_flavor`` = 2 (a widely used and tested option)
and flavor 5 (similar to 2, but enhanced by the use of an gamma distribution
instead of a normal distribution).
"Prior" inflation is generally a better choice than "posterior",
so set ``input.nml:filter_nml:``

.. code-block:: fortran

   inf_initial_from_restart    = .true.,   .false.
   inf_sd_initial_from_restart = .true.,   .false.

For the first cycle, if you have inflation restart files,
you should stage those in the $RUNDIR where the other restart files
will be staged, with names which include "dart.rh.cam_output_priorinf_mean"
and "dart.rh.cam_output_priorinf_sd" in them, so that assimilate.csh will find them.
If you don't have restart files, set ``*initial_from_restart`` to .false.
and assimilate.csh will create inflation restart files
using the values in ``inf_initial`` and ``inf_sd_initial``.
You will need to run the assimilation for some days in order to allow the inflation values
to equilibrate with the observation network and model ensemble spread.

.. _Perturbed:

Perturbed Ensemble
~~~~~~~~~~~~~~~~~~

A multidriver configuration of CAM needs an ensemble of initial condition files
for each active component in order to start a hindcast.
The set of files must include, at a minimum, CAM initial files and CLM restart files.
Usually CICE is also active, and other components may be,
which need their own restart files.
If there is no suitable initial ensemble for starting the ensemble hindcast,
one can be generated from a single model state
by linking it into suitably named files 
(see ../CESM/shell_scripts/link_ens_to_single.csh),
running the first ensemble hindcast, 
and then telling DART to perturb each member before the first assimilation.

The default perturbation routine in filter adds gaussian noise equally 
to all fields in the state vector. 
For CAM it is preferable to use the perturbation mechanism
in the cam-fv/model_mod.f90.
This allows the exclusion of fields which are tricky to perturb, 
such as specific humidity. 
The mechanism is controlled by the input.nml:model_nml "perturb" variables.
Typically, ensemble spread is generated from a single state by adding small 
perturbations to only the temperature field "T" and letting the model 
expand the perturbations to other fields and increase the sizes. 
For example,

.. code-block:: fortran

   filter_nml:
      single_file_in               = .false., (Even though your initial ensemble may be linked to a single file)
      perturb_from_single_instance = .true.
      perturbation_amplitude         (ignored, because model_mod defines it)

   model_nml:
      custom_routine_to_generate_ensemble = .true.
      fields_to_perturb                   = 'QTY_TEMPERATURE'
      perturbation_amplitude              = 0.1


Continuing after the first cycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your first hindcast+assimilation cycle uses an ensemble created from a single file,
you will need to change to the 'continuing' mode, 
where CAM will not perform all of its startup procedures 
and DART will use the most recently created ensemble.

.. code-block:: fortran

   ! model_nml:
      custom_routine_to_generate_ensemble = .true.
      fields_to_perturb                   = ''   (Turns off perturbations)
      perturbation_amplitude              = 0.1  (Ignored.  Can change to 0.0_r8 for consistency)

   ! CESM's env_run.xml:
       <entry id="CONTINUE_RUN" value="TRUE">

.. FIXME the ! allow it to be 'lexed' as fortran, but the ' confuses the syntax highlighting.

Combining multiple cycles into one job
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Setup_hybrid`` and ``setup_pmo`` are set up in the default cycling mode,
where each submitted job performs one model advance and one assimilation,
then resubmits the next cycle as a new job. 
For long series of cycles, this can result in a lot of time waiting in the queue 
for short jobs to run. Prevent this by using CESM's multicycling mode.
To request 2 hours to run 8 assimilation cycles, in $CASEROOT run commands:

.. code-block:: csh

 =  ./xmlchange DATA_ASSIMILATION_CYCLES=8
   ./xmlchange --subgroup case.run --id JOB_WALLCLOCK_TIME      --val 2:00:00
   ./xmlchange --subgroup case.run --id USER_REQUESTED_WALLTIME --val 2:00


.. _Diffusion:

Diffusion Near the Model Top
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CAM applies extra diffusion to the top levels of the model.
The number of levels is indirectly controlled 
by the CAM namelist variable ``div24del2flag``. 
It's not productive to assimilate in those levels
because of the distorting effects of the diffusion,
so the cam-fv/model_mod namelist has variables to prevent assimilation there.
``Model_damping_ends_at_level`` can be set to the same value 
that is activated by div24del2flag, or larger.
An alternative way to prevent assimilation in those layers 
is to exclude high observations using ``no_obs_assim_above_level``.
The CAM6 reanalysis mentioned above used this option,
with no_obs_assim_above_level = 5.
If model_damping_ends_at_level is turned on (has a value other than -1)
it's still sensible to exclude high observations using no_obs_assim_above_level.

It's worth considering the vertical localization when setting 
the value of no_obs_assim_above_level.
Observations at one level can affect model variables at other levels
if the localization is broad enough.
The effective vertical localization can be calculated by

.. code-block::

   cutoff * 2 * vert_normalization_{your_vert_coord} 

where ``cutoff`` is the half-width (hence the 2) 
of the horizontal  localization (radians)
and ``vert_normalization_...`` is the conversion from radians 
to the vertical coordinate system you've chosen using ``vertical_localization_coord``.
The resulting number can be compared against CAM's vertical levels
to decide which should be excluded.

Minimum Recommended Values to Control Assimilation Near the Model Top.

.. FIXME Sphinx renders a cell that is 2 rows deep as 1 row deep,
   even if there is text in both rows.
   +-----+
   | one | 
   | two | 
   +-----+
   yields
   ___________
   | one two |
   -----------
   So I can't split model_damping_ends_at_level or Diffusion levels i
   onto 2 lines to make the table narrower.

+---------------+------------------+-----------------------------+------------------------------+
| div24del2flag | Diffusion levels | model_damping_ends_at_level | no_obs_assim_above_level     |
+===============+==================+=============================+==============================+
| CAM: 2        | 2                | 2                           | (2; depends on localization) |
+---------------+------------------+-----------------------------+------------------------------+
| WACCM: 2      | 3                | 3                           | (3; depends on localization) |
+---------------+------------------+-----------------------------+------------------------------+
| CAM  4, 24    | 3                | 3                           | (3; depends on localization) |
+---------------+------------------+-----------------------------+------------------------------+
| WACCM: 4, 24  | 4                | 4                           | (4; depends on localization) |
+---------------+------------------+-----------------------------+------------------------------+

WACCM
~~~~~

WACCM[#][-X] has a much higher top than the CAM versions, 
which requires the use of scale height as the vertical coordinate, 
instead of pressure, during assimilation. 
Another impact of the high top is that the number of top model levels with extra diffusion 
in the FV version is different than in the low-topped CAM-FV, 
so the ``div24del2flag`` options lead to the larger minimum values listed in the table above.

You may need to experiment to find the best choices of DART namelist variables
to use with WACCM, but a good place to start includes

.. code-block:: fortran

   use_log_vertical_scale          = .true.
   use_variable_mean_mass          = .true.
   vertical_localization_coord     = 'SCALEHEIGHT'
   vert_normalization_scale_height = 1.5
   cutoff                          = 0.15
   no_obs_assim_above_level        = 4,
   

In any case, make the following changes (or similar) to convert from a CAM setup
to a WACCM setup in ``setup_hybrid``:

.. code-block:: csh

   setenv compset     FWHIST
   setenv resolution  f19_f19  
   setenv refcase     {the CASE name of the initial condition file(s) (differs from this assimilation)}
   setenv refyear     {\                                           }
   setenv refmon      { >{the date of the initial condition file(s)}
   setenv refday      {/                                           }

If there are problems with instability in the WACCM foreasts, try changing some
of the following parameters in either the setup script or input.nml.

-  The default ``div24del2flag`` in WACCM is 4. 
   Change it in the CAM namelist section of the setup script to

   .. code-block:: csh

      echo " div24del2flag         = 2 "                       >> ${fname}

.. $cesm/components/cam/dynamics/fv/cd_core.F90
   which will use the ``cd_core.F90`` in SourceMods, which has doubled diffusion
   in the top layers compared to CAM.

-  Set a larger ``ATM_NCPL`` in the setup script.  
   The default for WACCM is 144 (per day).
   The default for WACCM-X is 288 (per day).
   It's safest to choose a value which will evenly divide an hour,
   (for WACCM: ATM_NCPL = 168 or 192 ... multiples of 24)
   but evenly dividing the hindcast period might work
   (for a 6 hour hindcast: ATM_NCPL = 148 or 152 ... multiples of 4).
   To convert an existing CASE, try changing the related namelist variables 
   ``$CASEROOT/user_nl_cpl:{component}_cpl_dt`` (component :math:`\neq` "rof")

   .. code-block:: fortran

      user_nl_cpl:
         atm_cpl_dt = 300
         glc_cpl_dt = 300
         ice_cpl_dt = 300
         lnd_cpl_dt = 300
         ocn_cpl_dt = 300
         wav_cpl_dt = 300

-  Increase model_damping_ends_at_level in input.nml

-  Set a larger nsplit and/or nspltvrm in the CAM namelist section
   of the setup script:

   .. code-block:: csh

      echo " nsplit         = 16 "                             >> ${fname}
      echo " nspltvrm       =  4 "                             >> ${fname}

-  Reduce ``inf_damping`` from the default value of ``0.9`` in ``input.nml``:

   .. code-block:: fortran

      inf_damping           = 0.6,                   0,

CAM-SE
------

DART requires more information than what is available in the default output files from CAM-SE.
Set the following options in the CESM ``user_nl_cam`` namelist to have CESM generate
the files required for DART.

    .. code-block:: text

       inithist               = 'ENDOFRUN' 
       se_write_all_corners = .true.


.. Files
   -----

   -  ``model_nml`` in ``input.nml``
   -  ``cam_phis.nc`` (CAM surface height file, often CAM's .h0. file in the CESM run environment)
   -  netCDF output state diagnostics files

Nitty gritty: Efficiency and Issues to Address
----------------------------------------------


.. warning::

   Experience on a variety of machines has shown that it is a very good idea
   to make sure your run-time environment has the following:

   .. code-block:: bash

       limit stacksize unlimited
       limit datasize unlimited

It may be very beneficial to set MPI environment variables to larger values than the defaults
in $CASEROOT/env_mach_specific.xml:

.. code-block:: xml
  
   <environment_variables>
     <env name="MPI_COMM_MAX">16383</env>
     <env name="MPI_GROUP_MAX">1024</env>

Reduce total core hours and queue wait times by finding the minimum number of whole nodes 
on which CAM will run reliably.  Use that number in the setup script for each member of the ensemble.

Reduce core hours wasted by the single tasked creation of the CESM namelists
before each hindcast by:

   * calling case.submit with the --skip-preview-namelists argument
   * replacing the cime/src/drivers/mct/cime_config/buildnml with the one in the `SourceMods`_ tar file.
     
-  ISSUE: Improve this page
    * Add links and references to this document.
    * Publications web page.
    * CAM-chem; link?  More description?

-  ISSUE?; ``model_interpolate`` assumes that obs with a vertical location have
   2 horizontal locations too. The state vector may have fields for which this
   isn't true, but no obs we've seen so far violate this assumption. It would
   have to be a synthetic/perfect_model obs, like some sort of average or
   parameter value.

-  ISSUE: the cam-se variable ``max_neighbors`` is set to 6, but could be set to 4 
   for non-refined grids. Is there a good mechanism for this? Is it worth the file space
   savings?

-  ISSUE: the cam-se variables ``x_planar`` and ``y_planar`` could be reduced in rank, 
   if no longer needed for testing and debugging.

References and Acknowledgements
-------------------------------

-  `CESM homepage <https://www.cesm.ucar.edu/models/cesm1.3/>`__

Ave Arellano did the first work with CAM-Chem, assimilating MOPPITT CO
observations into CAM-Chem. Jerome Barre and Benjamin Gaubert took up the
development work from Ave, and prompted several additions to DART, as well as
``model_mod.f90``.

Nick Pedatella developed the first vertical_localization_coord = 'SCALEHEIGHT'`` capability 
to enable assimilation using WACCM(-X).

Rafael Montuoro designed the first multicoupler in CESM.
