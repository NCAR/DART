ROMS
====

There are several DART users who have working DART interface code
to the Regional Ocean Modeling System (ROMS), as the model is a community ocean
model funded by the Office of Naval Research. Please visit `MyRoms <https://www.myroms.org/>`_
for more information on the model.

The lead developers are at Rutgers and UCLA, but the list of associate
developers is extensive. Please read `ROMS developers <https://www.myroms.org/index.php?page=roms_devs>`_
for more information.

If you are interested in running DART with this model please contact the DART
group at dart@ucar.edu for more information.  We are currently working with
collaborators to optimize the model_mod interface and associated scripting to
run data assimilation experiments with this model. We may be able to put you in
contact with the right people to get a copy of the code.

Overview
--------

This document describes the relationship between ROMS and DART and provides an
overview of how to perform ensemble data assimilation with ROMS to provide ocean
states that are consistent with the information provided by various ocean
observations.

Running ROMS is complicated. It is **strongly** recommended that you become very
familiar with running ROMS before you attempt a ROMS-DART assimilation
experiment. Running DART is complicated. It is **strongly** recommended that you
become very familiar with running DART before you attempt a ROMS-DART
assimilation experiment. Running ROMS-DART takes expertise in both areas.

We recommend working through the :doc:`DART tutorial <../../theory/readme>`
to learn the concepts of ensemble data assimilation and the capabilities of DART.

The ROMS code is not distributed with DART, it can be obtained from the `ROMS website <https://www.myroms.org>`_.
There you will also find instructions on how to compile and run ROMS. DART can
use the 'verification observations' from ROMS (basically the estimate of the
observation at the location and time computed as the model advances) so it
would be worthwhile to become familiar with that capability of ROMS.

DART calls these 'precomputed forward operators'. DART can also use observations
from the `World Ocean Database <https://www.nodc.noaa.gov/OC5/indprod.html>`_ -
WOD. The conversion from the WOD formats to the DART observation sequence format
is accomplished by the converters in the ``DART/observations/obs_converters/WOD``
directory.

The DART forward operators require interpolation from the ROMS terrain-following
and horizontally curvilinear orthogonal coordinates to the observation location.
Please contact us for more information about this interpolation.

Generating an initial ensemble
------------------------------

The ROMS interface provides the ability to create an ensemble of initial ROMS
history files from an initial file by using the
:doc:`/assimilation_code/programs/perturb_single_instance/perturb_single_instance`.
You can specify an ensemble of any size in the ``perturb_single_instance``
namelist in ``input.nml`` and this program will randomly perturb the 
temperature and salinity fields of an initial ROMS history file to generate 
the ensemble.

A note about filenames
----------------------

During the course of an experiment, many files are created. To make them unique,
the *ocean_time* is converted from "seconds since 1900-01-01 00:00:00" to the
equivalent number of DAYS. An *integer* number of days. The intent is to tag the
filename to reflect the valid time of the model state. This could be used as the
DSTART for the next cycle, so it makes sense to me. The confusion comes when
applied to the observation files.

The input observation files for the ROMS 4DVAR
system typically have a DSTART that designates the start of the forecast cycle
and the file must contain observation from DSTART to the end of the forecast.
Makes sense.

The model runs to the end of the forecast, harvesting the verification
observations along the way. So then DART converts all those verification
observations and tags that file ... with the same time tag as all the other
output files ... which reflects the *ocean_time* (converted to days). The input
observation file to ROMS will have a different DSTART time in the filename than
the corresponding verification files. Ugh. You are free to come up with a better
plan.

These are just examples...after all; hopefully good examples.

Procedure
---------

The procedure to perform an assimilation experiment is outlined in the following
steps:

#. Compile ROMS (as per the ROMS instructions).
#. Compile all the DART executables (in the normal fashion).
#. Stage a directory with all the files required to advance an ensemble
   of ROMS models and DART.
#. Modify the run-time controls in ``ocean.in``, ``s4dvar.in`` and
   ``input.nml``. Since ROMS has a *Bin/subsitute* command, it is used to
   replace temporary placeholders with actual values at various parts
   during the process.
#. Advance all the instances of ROMS; each one will produce a restart
   file and a verification observation file.
#. Convert all the verification observation files into a single DART
   observation sequence file with the
   ``convert_roms_obs.f90`` program in ``DART/observations/obs_converters/ROMS/``.
#. Run filter to assimilate the data (DART will read and update the ROMS files
   directly - no conversion is necessary.)
#. Update the control files for ROMS in preparation for the next model
   advance.

Shell scripts
-------------

The ``shell_scripts`` directory has several scripts that are intended to
provide examples. These scripts **WILL** need to be modified to work on
your system and are heavily internally commented. It will be necessary
to read through and understand the scripts. As mentioned before, the
ROMS *Bin/subsitute* command is used to replace temporary placeholders
with actual values at various parts during the process.

+----------------------------------+----------------------------------+
| Script                           | Description                      |
+==================================+==================================+
| ensemble.sh                      | Was written by Hernan Arango to  |
|                                  | run an ensemble of ROMS models.  |
|                                  | It is an appropriate example of  |
|                                  | what is required from the ROMS   |
|                                  | perspective. It does no data     |
|                                  | assimilation.                    |
+----------------------------------+----------------------------------+
| stage_experiment.csh             | prepares a directory for an      |
|                                  | assimilation experiment. The     |
|                                  | idea is basically that           |
|                                  | everything you need should be    |
|                                  | assembled by this script and     |
|                                  | that this should only be run     |
|                                  | ONCE per experiment. After       |
|                                  | everything is staged in the      |
|                                  | experiment directory, another    |
|                                  | script can be run to advance the |
|                                  | model and perform the            |
|                                  | assimilation.                    |
|                                  | *stage_experiment.csh* will also |
|                                  | modify some of the template      |
|                                  | scripts and copy working         |
|                                  | versions into the experiment     |
|                                  | directory. This script may be    |
|                                  | run interactively, i.e. from the |
|                                  | UNIX command line.               |
+----------------------------------+----------------------------------+
| submit_multiple_cycles_lsf.csh   | is an executable script that     |
|                                  | submits a series of dependent    |
|                                  | jobs to an LSF queuing system.   |
|                                  | Each job runs *cycle.csh* in the |
|                                  | experiment directory and only    |
|                                  | runs if the previous dependent   |
|                                  | job completes successfully.      |
+----------------------------------+----------------------------------+
| cycle.csh.template               | is a non-executable template     |
|                                  | that is modified by              |
|                                  | *stage_experiment.csh* and       |
|                                  | results in an exectuable         |
|                                  | *cycle.csh* in the experiment    |
|                                  | directory. *cycle.csh* is        |
|                                  | designed to be run as a batch    |
|                                  | job and advances the ROMS model  |
|                                  | states one-by-one for the        |
|                                  | desired forecast length. The     |
|                                  | assimilation is performed and    |
|                                  | the control information for the  |
|                                  | next ROMS forecast is updated.   |
|                                  | Each model execution and         |
|                                  | *filter* use the same set of MPI |
|                                  | tasks.                           |
+----------------------------------+----------------------------------+
| submit_multiple_jobs_slurm.csh   | is an executable script that     |
|                                  | submits a series of dependent    |
|                                  | jobs to an LSF queuing system.   |
|                                  | It is possible to submit         |
|                                  | **many** jobs the queue, but the |
|                                  | jobs run one-at-a-time. Every    |
|                                  | assimilation cycle is divided    |
|                                  | into two scripts to be able to   |
|                                  | efficiently set the resources    |
|                                  | for each phase.                  |
|                                  | *advance_ensemble.csh* is a job  |
|                                  | array that advances each ROMS    |
|                                  | instance in separate jobs. When  |
|                                  | the entire job array finishes -  |
|                                  | and only if they all finish      |
|                                  | correctly - will the next job    |
|                                  | start to run. *run_filter.csh*   |
|                                  | performs the assimilation and    |
|                                  | prepares the experiment          |
|                                  | directory for another            |
|                                  | assimilation cycle.              |
|                                  | *submit_multiple_jobs_slurm.csh* |
|                                  | may be run from the command line |
|                                  | in the experiment directory.     |
|                                  | Multiple assimilation cycles can |
|                                  | be specified, so it is possible  |
|                                  | to put **many** jobs in the      |
|                                  | queue.                           |
+----------------------------------+----------------------------------+
| advance_ensemble.csh.template    | is a non-executable template     |
|                                  | that is modified by              |
|                                  | *stage_experiment.csh* and       |
|                                  | results in an exectuable         |
|                                  | *advance_ensemble.csh* in the    |
|                                  | experiment directory.            |
|                                  | *advance_ensemble.csh* is        |
|                                  | designed to submit an job array  |
|                                  | to the queueing system           |
|                                  | (PBS,SLURM, or LSF) to advance   |
|                                  | the ensemble members in separate |
|                                  | jobs.                            |
+----------------------------------+----------------------------------+
| run_filter.csh.template          | is a non-executable template     |
|                                  | that is modified by              |
|                                  | *stage_experiment.csh* and       |
|                                  | results in an exectuable         |
|                                  | *run_filter.csh* in the          |
|                                  | experiment directory.            |
|                                  | *run_filter.csh* is very similar |
|                                  | to *cycle.csh* but does not      |
|                                  | advance the ROMS model           |
|                                  | instances.                       |
+----------------------------------+----------------------------------+

The variables from ROMS that are copied into the DART state vector are
controlled by the *input.nml* *model_nml* namelist. See below for the
documentation on the &model_nml entries. The state vector should include all
variables needed to apply the forward observation operators as well as the
prognostic variables important to restart ROMS.

The example *input.nml* *model_nml* demonstrates how to construct the DART state
vector. The following table explains in detail each entry for the *variables*
namelist item:

+----------------+-----------------------------------+
| Variable name  | This is the ROMS variable name as |
|                | it appears in the ROMS netCDF     |
|                | file.                             |
+----------------+-----------------------------------+
| DART QUANTITY  | This is the character string of   |
|                | the corresponding DART QUANTITY.  |
|                | The complete list of possible     |
|                | DART QUANTITY values is available |
|                | in the ``obs_def_mod``            |
|                | that is built by ``preprocess``.  |
+----------------+-----------------------------------+
| minimum        | If the variable is to be updated  |
|                | in the ROMS restart file, this    |
|                | specifies the minimum value. If   |
|                | set to 'NA', there is no minimum  |
|                | value.                            |
+----------------+-----------------------------------+
| maximum        | If the variable is to be updated  |
|                | in the ROMS restart file, this    |
|                | specifies the maximum value. If   |
|                | set to 'NA', there is no maximum  |
|                | value.                            |
+----------------+-----------------------------------+
| update         | The updated variable may or may   |
|                | not be written to the ROMS        |
|                | restart file.                     |
|                | *'UPDATE'*  means the variable in |
|                | the restart file is updated. This |
|                | is case-insensitive.              |
|                | *'NO_COPY_BACK'*  (or anything    |
|                | else) means the variable in the   |
|                | restart file remains unchanged.   |
+----------------+-----------------------------------+

Namelist
--------

This namelist is read from the file *input.nml*. Namelists start with an
ampersand '&' and terminate with a slash '/'. Character strings that
contain a '/' must be enclosed in quotes to prevent them from
prematurely terminating the namelist. The default namelist is presented
below, a more realistic namelist is presented at the end of this
section.

.. code-block:: fortran

   &model_nml
     roms_filename               = 'roms_input.nc'
     assimilation_period_days    = 1
     assimilation_period_seconds = 0
     vert_localization_coord     = 3
     debug                       = 0
     variables                   = ''
   /

+-----------------------+-----------------------+-----------------------+
| Item                  | Type                  | Description           |
+=======================+=======================+=======================+
| roms_filename         | character(len=256)    | This is the name of   |
|                       |                       | the file used to      |
|                       |                       | provide information   |
|                       |                       | about the ROMS        |
|                       |                       | variable dimensions,  |
|                       |                       | etc.                  |
+-----------------------+-----------------------+-----------------------+
| assi                  | integer               | Combined, these       |
| milation_period_days, |                       | specify the width of  |
| assimi                |                       | the assimilation      |
| lation_period_seconds |                       | window. The current   |
|                       |                       | model time is used as |
|                       |                       | the center time of    |
|                       |                       | the assimilation      |
|                       |                       | window. All           |
|                       |                       | observations in the   |
|                       |                       | assimilation window   |
|                       |                       | are assimilated.      |
|                       |                       | BEWARE: if you put    |
|                       |                       | observations that     |
|                       |                       | occur before the      |
|                       |                       | beginning of the      |
|                       |                       | assimilation_period,  |
|                       |                       | DART will error out   |
|                       |                       | because it cannot     |
|                       |                       | move the model 'back  |
|                       |                       | in time' to process   |
|                       |                       | these observations.   |
+-----------------------+-----------------------+-----------------------+
| variables             | character(:, 5)       | A 2D array of         |
|                       |                       | strings, 5 per ROMS   |
|                       |                       | variable to be added  |
|                       |                       | to the dart state     |
|                       |                       | vector.               |
|                       |                       |                       |
|                       |                       | #. ROMS field name -  |
|                       |                       |    must match netCDF  |
|                       |                       |    variable name      |
|                       |                       |    exactly            |
|                       |                       | #. DART QUANTITY -    |
|                       |                       |    must match a valid |
|                       |                       |    DART QTY_xxx       |
|                       |                       |    exactly            |
|                       |                       | #. minimum physical   |
|                       |                       |    value - if none,   |
|                       |                       |    use 'NA'           |
|                       |                       | #. maximum physical   |
|                       |                       |    value - if none,   |
|                       |                       |    use 'NA'           |
|                       |                       | #. case-insensitive   |
|                       |                       |    string describing  |
|                       |                       |    whether to copy    |
|                       |                       |    the updated        |
|                       |                       |    variable into the  |
|                       |                       |    ROMS restart file  |
|                       |                       |    ('UPDATE') or not  |
|                       |                       |    (any other value). |
|                       |                       |    There is generally |
|                       |                       |    no point copying   |
|                       |                       |    diagnostic         |
|                       |                       |    variables into the |
|                       |                       |    restart file. Some |
|                       |                       |    diagnostic         |
|                       |                       |    variables may be   |
|                       |                       |    useful for         |
|                       |                       |    computing forward  |
|                       |                       |    operators,         |
|                       |                       |    however.           |
+-----------------------+-----------------------+-----------------------+
| ve                    | integer               | Vertical coordinate   |
| rt_localization_coord |                       | for vertical          |
|                       |                       | localization.         |
|                       |                       |                       |
|                       |                       | -  1 = model level    |
|                       |                       | -  2 = pressure (in   |
|                       |                       |    pascals)           |
|                       |                       | -  3 = height (in     |
|                       |                       |    meters)            |
|                       |                       | -  4 = scale height   |
|                       |                       |    (unitless)         |
|                       |                       |                       |
|                       |                       | Currently, only 3     |
|                       |                       | (height) is supported |
|                       |                       | for ROMS.             |
+-----------------------+-----------------------+-----------------------+

A more realistic ROMS namelist is presented here, along with one of the
more unusual settings that is generally necessary when running ROMS. The
*use_precomputed_FOs_these_obs_types* variable needs to list the
observation types that are present in the ROMS verification observation
file.

.. code-block:: fortran

   &model_nml
     roms_filename                = 'roms_input.nc'
     assimilation_period_days     = 1
     assimilation_period_seconds  = 0
     vert_localization_coord      = 3
     debug                        = 1
     variables = 'temp',   'QTY_TEMPERATURE',          'NA', 'NA', 'update',
                 'salt',   'QTY_SALINITY',            '0.0', 'NA', 'update',
                 'u',      'QTY_U_CURRENT_COMPONENT',  'NA', 'NA', 'update',
                 'v',      'QTY_V_CURRENT_COMPONENT',  'NA', 'NA', 'update',
                 'zeta',   'QTY_SEA_SURFACE_HEIGHT'    'NA', 'NA', 'update'
   /
   &obs_kind_nml
     evaluate_these_obs_types = ''
     assimilate_these_obs_types =          'SATELLITE_SSH',
                                           'SATELLITE_SSS',
                                           'XBT_TEMPERATURE',
                                           'CTD_TEMPERATURE',
                                           'CTD_SALINITY',
                                           'ARGO_TEMPERATURE',
                                           'ARGO_SALINITY',
                                           'GLIDER_TEMPERATURE',
                                           'GLIDER_SALINITY',
                                           'SATELLITE_BLENDED_SST',
                                           'SATELLITE_MICROWAVE_SST',
                                           'SATELLITE_INFRARED_SST'
     use_precomputed_FOs_these_obs_types = 'SATELLITE_SSH',
                                           'SATELLITE_SSS',
                                           'XBT_TEMPERATURE',
                                           'CTD_TEMPERATURE',
                                           'CTD_SALINITY',
                                           'ARGO_TEMPERATURE',
                                           'ARGO_SALINITY',
                                           'GLIDER_TEMPERATURE',
                                           'GLIDER_SALINITY',
                                           'SATELLITE_BLENDED_SST',
                                           'SATELLITE_MICROWAVE_SST',
                                           'SATELLITE_INFRARED_SST'
   /
