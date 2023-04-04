
WRF/DART Tutorial Materials for the Manhattan Release.
======================================================


Introduction
------------

This document will describe how to get started with your own Weather
Research and Forecasting (WRF) data assimilation experiments using DART
and only covers only the WRF-specific aspects of integrating with DART.
It is not wise to try to run WRF/DART if you have no experience with WRF
and/or no experience with DART.

This tutorial was assembled to be compatible with ~WRF V3.9.1 and the
DART Manhattan release. Other releases of WRF may or may not be
backwards or forwards compatible with this tutorial.

You must already be comfortable running the
`WRF <http://www2.mmm.ucar.edu/wrf/users/download/get_source.html>`__
system (WPS, real_em build of WRF). If not, work through the `WRF model
tutorial <https://www.mmm.ucar.edu/wrf-tutorial-0>`__
first before trying to link WRF and DART together. Check the WRF user
guide or the
`WRFHELP <https://www.mmm.ucar.edu/wrf-user-support-contributor-information>`__
forum for WRF-specific assistance.

If you are new to DART, we recommend that you become familiar with DART
by working through the :doc:`../../../theory/readme` and then
understanding the :doc:`DART getting started <../../../README>` documentation.

before attempting the WRF/DART tutorial as you will find many helpful
resources for learning the base DART configuration.

*We do not claim that this is a “turnkey” or “black box” system.* Be
mentally prepared to invest a reasonable amount of time on the learning
curve. There are many outstanding research issues which have no easy
answers. This is not a one week/grad student/naive user system. Even
after you get the code up and running, you have to be able to interpret
the results, which requires developing specific skills. There are a lot
of ways to alter how the system works – localization, inflation, which
variables and observations are assimilated, the assimilation window
time, the model resolution, etc, etc. This is both good and bad - you
have many ways of improving your results, but you have to take care on
how you leave all the settings of these inputs. Getting a set of scripts
that runs doesn’t mean the system is running well, or producing useful
results. So - if you’re still reading: Let the adventure begin!

This tutorial introduces a “canned” WRF/DART experiment involving an
ensemble of 50 members that will be initialized from GFS initial
conditions at 2017/04/27 00:00 UTC using a domain of the continental
United States. The data included in the tutorial lasts until 2017/04/30
18:00 UTC. During this period, there was a strong rain and wind event
that affected a large portion of the United States, causing record
rains, localized flooding, and numerous tornadoes. For more information
on the physical account of this case, see
`weather.gov <https://www.weather.gov/lot/2017Apr2930_rainfall>`__.

By default, the tutorial case will only cover 12 hours of this event
starting at 2017/04/27 00:00 UTC. The WRF model will be “spun-up” for
six hours to generate a prior distribution. An assimilation of PREPBUFR
observations will then be performed at 06:00 UTC, at which time analysis
files will be generated to begin a new ensemble forecast. The WRF model
will be advanced for 6 hours and a final assimilation cycle will be
performed at 12:00 UTC. This process could then continue in order to
investigate the strong rain and wind event. For what it’s worth, on
NCAR’s *Cheyenne* under the default test configuration for this case, it
can take an hour to complete a forecast/assimilation cycle. Since the
tutorial runs for two cycles, it can take twice as long.

The goals of this tutorial are to demonstrate how WRF/DART works. After
running this tutorial, you will be able to understand the major steps
involved in setting up your own data assimilation (DA) experiments.
However, you will need to do additional work before you can expect to
have a fully functional WRF/DART system, as some of the steps involved
in this tutorial (in particular, the perturbation bank and the
observation sequence files) are provided for you in order to simplify
the process. Furthermore, if you are not running on the UCAR/NCAR
Cheyenne supercomputing system, you will likely need to customize the
assimilation scripts to match the details of your particular system.


.. important ::

  We have provided instructions for the NCAR supercomputer
  Cheyenne, so you may need to tailor these instructions to your system if
  you are not using Cheyenne. These system-specific setup steps may take a
  good deal of effort, especially if you are unfamiliar with details such
  as MPI, NetCDF, etc. Furthermore, even after you get the code up and
  running, you will need to properly interpret your results.


Step 1: Setup
-------------

There are several dependencies for the executables and scripting
components. On Cheyennne, users have reported success building WRF, WPS,
WRFDA, and DART with the default module environment including Intel
compilers, MPT, and netCDF4. In addition, you'll need to load the
`nco <http://nco.sourceforge.net/>`__ and
`ncl <https://www.ncl.ucar.edu/>`__ modules to run the set of scripts
that accompany the tutorial.

There are multiple phases for the setup: building the DART executables,
getting the initial WRF boundary conditions etc., building (or using
existing) WRF executables, and configuring and staging the scripting
needed to perform an experiment.

Build the DART executables.
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have not already, see :doc:`Getting Started <../../../README>` to
download the DART software package. Set an environment variable
*DART_DIR* to point to your base DART directory. How to do this will
depend on which shell you are using.

===== ====================================================
shell command
===== ====================================================
tcsh  ``setenv DART_DIR <path_to_your_dart_installation>``
bash  ``export DART_DIR=<path_to_your_dart_installation>``
===== ====================================================

In either case, you will replace <path_to_your_dart_installation> with
the actual path to your DART installation. If you are using another
shell, refer to your shell-specific documentation on how to set an
environment variable.

Building the DART executables for the tutorial follows the same process
as building any of the DART executables. Configure the ``mkmf.template``
file for your system, configure the ``input.nml`` for the model you want
to compile, and run ``quickbuild.sh`` (which is not necessarily quick,
but it is quicker than doing it by hand) to compile all the programs you
might need for an experiment with that model.

1. It is assumed you have successfully configured the
   ``$DART_DIR/build_templates/mkmf.template`` file for your system. If
   not, you will need to do so now. See the :doc:`Getting Started <../../../README>`
   for more detail, if necessary.

2. [OPTIONAL] Modify the DART code to use 32bit reals. Most WRF/DART
   users run both the WRF model and the DART assimilation code using
   32bit reals. This is not the default for the DART code. Make this
   single code change before building the DART executables to compile
   all reals as 32bit reals.

   Edit ``$DART_DIR/assimilation_code/modules/utilities/types_mod.f90``
   with your favorite editor. Change

   ::

     ! real precision:
     ! TO RUN WITH REDUCED PRECISION REALS (and use correspondingly less memory)
     ! comment OUT the r8 definition below and use the second one:
     integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
     integer, parameter :: r8 = SELECTED_REAL_KIND(12)   ! 8 byte reals
     !integer, parameter :: r8 = r4                      ! alias r8 to r4
  
   to

   ::

       ! real precision:
       ! TO RUN WITH REDUCED PRECISION REALS (and use correspondingly less memory)
       ! comment OUT the r8 definition below and use the second one:
       integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
       ! integer, parameter :: r8 = SELECTED_REAL_KIND(12)   ! 8 byte reals
       integer, parameter :: r8 = r4                      ! alias r8 to r4

3. Copy the tutorial DART namelist from
   ``$DART_DIR/models/wrf/tutorial/template/input.nml.template`` to
   ``$DART_DIR/models/wrf/work/input.nml``.

   ::

      cd $DART_DIR/models/wrf
      cp tutorial/template/input.nml.template work/input.nml

4. Build the WRF/DART executables:

   ::

      cd $DART_DIR/models/wrf/work
      ./quickbuild.sh

   Many executables are built, the following executables are needed for the
   tutorial and will be copied to the right place by the *setup.csh* script
   in a subsequent step:
 
   ::

      advance_time
      fill_inflation_restart
      filter
      obs_diag
      obs_seq_to_netcdf
      obs_sequence_tool
      pert_wrf_bc
      wrf_dart_obs_preprocess

Preparing the experiment directory.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Approximately 100Gb of space is needed to run the tutorial. Create a
"work" directory someplace with a lot of free space. The rest of the
instructions assume you have an environment variable called *BASE_DIR*
that points to this directory.

===== ====================================================
shell command
===== ====================================================
tcsh  ``setenv BASE_DIR <path_to_your_working_directory>``
bash  ``export BASE_DIR=<path_to_your_working_directory>``
===== ====================================================

1. The WRF boundary conditions and perturbations required to make a
   viable ensemble are available in a 15 GB tar file. Put this file in
   your ``$BASE_DIR``. Since this is a large file, we suggest using
   'wget' to download the file directly to your local system:

   ::

       cd $BASE_DIR
       wget http://www.image.ucar.edu/wrfdart/tutorial/wrf_dart_tutorial_23May2018_v3.tar.gz
       tar -xzvf wrf_dart_tutorial_23May2018_v3.tar.gz

   After untarring the file you should see the following directories:
   *icbc, output, perts,* and *template.* The directory names (case
   sensitive) are important, as the scripts rely on these local paths
   and file names.

2. You will need template WRF namelists from the
   ``$DART_DIR/models/wrf/tutorial/template`` directory:

   ::

       cp $DART_DIR/models/wrf/tutorial/template/namelist.input.meso   $BASE_DIR/template/.
       cp $DART_DIR/models/wrf/tutorial/template/namelist.wps.template $BASE_DIR/template/.

3. You will also need scripting to run a WRF/DART experiment. Copy the contents of 
   ``$DART_DIR/models/wrf/shell_scripts`` to the ``$BASE_DIR/scripts`` directory.

   ::

       mkdir $BASE_DIR/scripts
       cp -R $DART_DIR/models/wrf/shell_scripts/* $BASE_DIR/scripts

Build or locate WRF executables.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The
`WRFDA <http://www2.mmm.ucar.edu/wrf/users/wrfda/download/get_source.html>`__
package is needed to generate a set of perturbed initial ensemble member
files and also to generate perturbed boundary condition files. Since the
tutorial provides a perturbation bank for a specific case, it is not
required to actually *run da_wrfvar.exe* but it needs to be in the
``WRF_RUN`` directory for the tutorial.

Build (or locate an appropriate build of) WRF, WPS and WRFDA.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WRF and WRFDA should be built with the "dmpar" option, while WPS can be
built "serial"ly. See the WRF/WRFDA documentation for more information
about building these packages. 

.. note::
	
 For consistency and to avoid errors, you should build WRF, WPS, WRFDA, and DART with the
 same compiler you use for NetCDF. Likewise MPI should use the same compiler.
 You will need the location of the WRF and WRFDA builds to customize the
 *params.csh* script in the next step.

Configure ``$BASE_DIR/scripts/param.csh`` with proper paths, info, etc.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a script that sets variables which will be read by other
WRF/DART scripts. There are some specific parameters for either the
Cheyenne supercomputing system using the
`PBS <https://www.pbsworks.com/>`__ queueing system or the
(decommissioned) Yellowstone system which used the *LSF* queueing
system. If you are not using Cheyenne, you may still want to use this
script to set your queueing-system specific parameters.

.. important::

   All variables that are marked
   ``'set this appropriately #%%%#'`` need to be set. This list is intended
   to provide some guidance on what needs to be set, but it is not an
   exhaustive list.

 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 |     Script variable     |                                                                     Description                                                                     |
 +=========================+=====================================================================================================================================================+
 | module load mpt         | The Environment Modules MPI compiler to use (here the HPE MPI) compiler). Note that on Cheyenne the default compiler is Intel.                      |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | module load nco         | The nco package.                                                                                                                                    |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | module load ncl/6.6.2   | The ncl package.                                                                                                                                    |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | BASE_DIR                | The directory containing icbc, output, perts, etc.                                                                                                  |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | DART_DIR                | The DART directory.                                                                                                                                 |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | WRF_DM_SRC_DIR          | The directory of the WRF dmpar installation.                                                                                                        |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | WPS_SRC_DIR             | The directory of the WPS installation.                                                                                                              |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | VAR_SRC_DIR             | The directory of the WRFDA installation.                                                                                                            |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | GEO_FILES_DIR           | The root directory of the WPS_GEOG files. NOTE: on Cheyenne these are available in the /glade/u/home/wrfhelp/WPS_GEOG directory                     |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | GRIB_DATA_DIR           | The root directory of the GRIB data input into ungrib.exe. For this tutorial the grib files are included, so use ${ICBC_DIR}/grib_data              |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | GRIB_SRC                | The type of GRIB data (e.g. <Vtable.TYPE>) to use with ungrib.exe to copy the appropriate Vtable file. For the tutorial, the value should be 'GFS'. |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | COMPUTER_CHARGE_ACCOUNT | The project account for supercomputing charges. See your supercomputing project administrator for more information.                                 |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
 | EMAIL                   | The e-mail address used by the queueing system to send job summary information.                                                                     |
 +-------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+


Run the *setup.csh* script to create the proper directory structure and
move executables to proper locations.

::

   cd $BASE_DIR/scripts
   ./setup.csh param.csh

So far, your ``$BASE_DIR`` should contain the following directories:

::

   icbc
   obs_diag
   obsproc
   output
   perts
   post
   rundir
   scripts
   template

Your ``$BASE_DIR/rundir`` directory should contain the following:

**executables:**

 
- `advance_time <../../../assimilation_code/programs/advance_time/advance_time.html>`__,
- `fill_inflation_restart <../../../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart.html>`__,
- `filter <../../../assimilation_code/programs/filter/filter.html>`__,
- `obs_diag <../../../assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`__,
- `obs_seq_to_netcdf <../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html>`__,
- `obs_sequence_tool <../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html>`__,
- ``pert_wrf_bc`` (no helper page),
- `wrf_dart_obs_preprocess <../../../models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess.html>`__

**directories:** 

- ``WRFIN`` (empty)
- ``WRFOUT`` (empty)
- ``WRF_RUN`` (wrf executables and support files)


**scripts:** 

- *add_bank_perts.ncl*
- *new_advance_model.csh*

**support data:** 

- *sampling_error_correction_table.nc*

Check to make sure your ``$BASE_DIR/rundir/WRF_RUN`` directory contains:

::

   da_wrfvar.exe
   wrf.exe
   real.exe
   be.dat
   contents of your WRF build run/ directory (support data files for WRF)

.. note::

	
   Be aware that the *setup.csh* script is designed to remove
   ``$BASE_DIR/rundir/WRF_RUN/namelist.input``. Subsequent scripting will
   modify ``$BASE_DIR/template/namlist.input.meso`` to create the
   ``namelist.input`` for the experiment.

For this tutorial, we are providing you with a specified WRF domain. To
make your own, you would need to define your own wps namelist and use
WPS to make your own geogrid files. See the WRF site for help with
building and running those tools as needed. You would also need to get
the appropriate grib files to generate initial and boundary condition
files for the full period you plan to cycle. In this tutorial we have
provided you with geogrid files, a small set of grib files, and a
namelist to generate series of analyses for several days covering a
North American region.

Let's now look inside the ``$BASE_DIR/scripts`` directory. You should
find the following scripts:

+-----------------------+-------------------------------------------------------------------------------------------+
|      Script name      |                                        Description                                        |
+=======================+===========================================================================================+
| add_bank_perts.ncl    | Adds perturbations to each member.                                                        |
+-----------------------+-------------------------------------------------------------------------------------------+
| assim_advance.csh     | Advances 1 WRF ensemble member to the next analysis time.                                 |
+-----------------------+-------------------------------------------------------------------------------------------+
| assimilate.csh        | Runs filter ... i.e. the assimilation.                                                    |
+-----------------------+-------------------------------------------------------------------------------------------+
| diagnostics_obs.csh   | Computes observation-space diagnostics and the model-space mean analysis increment.       |
+-----------------------+-------------------------------------------------------------------------------------------+
| driver.csh            | Primary script for running the cycled analysis system.                                    |
+-----------------------+-------------------------------------------------------------------------------------------+
| first_advance.csh     | Advances 1 WRF ensemble member (on the first time).                                       |
+-----------------------+-------------------------------------------------------------------------------------------+
| gen_pert_bank.csh     | Saves the perturbations generated by WRFDA CV3.                                           |
+-----------------------+-------------------------------------------------------------------------------------------+
| gen_retro_icbc.csh    | Generates the wrfinput and wrfbdy files.                                                  |
+-----------------------+-------------------------------------------------------------------------------------------+
| init_ensemble_var.csh | Creates the perturbed initial conditions from the WRF-VAR system.                         |
+-----------------------+-------------------------------------------------------------------------------------------+
| mean_increment.ncl    | Computes the mean state-space increment, which can be used for plotting.                  |
+-----------------------+-------------------------------------------------------------------------------------------+
| new_advance_model.csh | advances the WRF model after running DART in a cycling context.                           |
+-----------------------+-------------------------------------------------------------------------------------------+
| param.csh             | Contains most of the key settings to run the WRF/DART system.                             |
+-----------------------+-------------------------------------------------------------------------------------------+
| prep_ic.csh           | Prepares the initial conditions for a single ensemble member.                             |
+-----------------------+-------------------------------------------------------------------------------------------+
| real.csh              | Runs the WRF real.exe program.                                                            |
+-----------------------+-------------------------------------------------------------------------------------------+
| setup.csh             | Creates the proper directory structure and place executables/scripts in proper locations. |
+-----------------------+-------------------------------------------------------------------------------------------+



You will need to edit the following scripts to provide the paths to
where you are running the experiment, to connect up files, and to set
desired dates. Search for the string ``'set this appropriately #%%%#'``
for locations that you need to edit.

::

   cd $BASE_DIR/scripts
   grep -r 'set this appropriately #%%%#' .

Other than *param.csh*, which was covered above, make the following
changes:

+--------------------+--------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|      File name     |           Variable / value           |                                                                                                                    Change description                                                                                                                   |
+====================+======================================+=========================================================================================================================================================================================================================================================+
| driver.csh         | datefnl = 2017042712                 | Change to the final target date; here the final date is already set correctly for this tutorial.                                                                                                                                                        |
+--------------------+--------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| gen_retro_icbc.csh | datefnl = 2017042712                 | Set to the final target date of the tutorial.  However, it is possible (not necessary) to create WRF initial/boundary conditions to 2017043000. This is the latest date that files are included in the tutorial.                                        |
+--------------------+--------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| gen_retro_icbc.csh | paramfile = <full path to param.csh> | The full path to param.csh. Change this on the line after the comment. While these two files are in the same directory here, in general it is helpful to have one param.csh for each experiment.                                                        |
+--------------------+--------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| gen_pert_bank.csh  | All changes                          | As the tutorial includes a perturbation bank, you will not need to run this script for the tutorial, so you will not need to change these values. However, you should set appropriate values when you are ready to generate your own perturbation bank. |
+--------------------+--------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Next, move to the ``$BASE_DIR/perts`` directory. Here you will find 100
perturbation files, called a "perturbation bank." For your own case, you
would need to create a perturbation bank of your own. A brief
description for running the script is available inside the comments of
that file. However, again, for this tutorial, this step has already been
run for you. The ``$BASE_DIR/icbc`` directory contains a *geo_em_d01.nc*
file (geo information for our test domain), and grib files that will be
used to generate the initial and boundary condition files. The
``$BASE_DIR/template`` directory should contain namelists for WRF, WPS,
and filter, along with a wrfinput file that matches what will be the
analysis domain. Finally, the ``$BASE_DIR/output`` directory contains
observations within each directory name. Template files will be placed
here once created (done below), and as we get into the cycling the
output will go in these directories.




Step 2: Initial conditions
--------------------------

To get an initial set of ensemble files, depending on the size of your
ensemble and data available to you, you might have options to initialize
the ensemble from, say, a global ensemble set of states. Here, we
develop a set of flow dependent errors by starting with random
perturbations and conducting a short forecast. We will use the WRFDA
random CV option 3 to provide an initial set of random errors, and since
this is already available in the perturbation bank developed in the
setup, we can simply add these to a deterministic GFS state. Further,
lateral boundary uncertainty will come from adding a random perturbation
to the forecast (target) lateral boundary state, such that after the
integration the lateral boundaries have random errors.

First, we need to generate a set of GFS states and boundary conditions
that will be used in the cycling. Use
``$BASE_DIR/scripts/gen_retro_icbc.csh`` to create this set of files,
which will be added to a subdirectory corresponding to the date of the
run in the ``$BASE_DIR/output`` directory. Make sure
*gen_retro_icbc.csh* has the appropriate path to your *param.csh*
script. If the *param.csh* script also has the correct edits for paths
and you have the executables placed in the rundir, etc., then running
*gen_retro_icbc.csh* should execute a series of operations to extract
the grib data, run metgrid, and then twice execute *real.exe* to
generate a pair of WRF files and a boundary file for each analysis time.

::

   cd $BASE_DIR/scripts
   ./gen_retro_icbc.csh


.. note::

  Ignore any ``rm: No match`` errors, as the script attempts to
  delete output files if they already exist, and they will not for the
  first run.

Once the script completes, inside your ``$BASE_DIR/output/2017042700``
directory you should see these files:

::

   wrfbdy_d01_152057_21600_mean
   wrfinput_d01_152057_0_mean
   wrfinput_d01_152057_21600_mean

These filenames include the Gregorian dates for these files, which is
used by the dart software for time schedules. Similar files (with
different dates) should appear in all of the date directories between
the *datea* and *datef* dates set in the *gen_retro_icbc.csh* script.
All directories with later dates will also have an observation sequence
file *obs_seq.out* that contains observations to be assimilated at that
time.

Next, we will execute the script to generate an initial ensemble of
states for the first analysis. For this we run the script
*init_ensemble_var.csh*, which takes two arguments: a date string and
the location of the *param.csh* script.

::

   cd $BASE_DIR/scripts
   ./init_ensemble_var.csh 2017042700 param.csh

This script generates 50 small scripts and submits them to the batch
system. It assumes a PBS batch system and the 'qsub' command for
submitting jobs. If you have a different batch system, edit this script
and look near the end. You will need to modify the lines staring with
#PBS and change 'qsub' to the right command for your system. You might
also want to modify this script to test running a single member first —
just in case you have some debugging to do.

However, be warned that to successfully complete the tutorial, including
running the *driver.csh* script in Step 5, using a smaller ensemble 
(e.g. < 20 members) can lead to spurious updates during the analysis step,
causing the WRF simulation to fail. 

When complete for the full ensemble, you should find 50 new files in the
directory ``output/2017042700/PRIORS`` with names like *prior_d01.0001*,
*prior_d01.0002*, etc... You may receive an e-mail to helpfully inform
you when each ensemble member has finished.


Step 3: Prepare observations [OPTIONAL]
---------------------------------------

For the tutorial exercise, observation sequence files are provided to
enable you to quickly get started running a test WRF/DART system. If you
want to run with the example observations, you can skip to Step
4.

However, observation processing is critical to the success of running
DART and was covered in :doc:`getting started <../../../README>`. In
brief, to add your own observations to WRF/DART you will need to
understand the relationship between observation definitions and
observation sequences, observation types and observation quantities, and
understand how observation converters extract observations from their
native formats into the DART specific format.

The observation sequence files that are provided in this tutorial come
from NCEP BUFR observations from the GDAS system. These observations
contain a wide array of observation types from many platforms within a
single file.

If you wanted to generate your own observation sequence files from
PREPBUFR for an experiment with WRF/DART, you should follow the guidance
on the
`prepbufr <../../../observations/obs_converters/NCEP/prep_bufr/prep_bufr.html>`__
page to build the bufr conversion programs, get observation files for
the dates you plan to build an analysis for, and run the codes to
generate an observation sequence file.

For completeness, we list here how you could generate these observation
sequence files yourself. 

.. important::

   the following steps are **not
   necessary** for the tutorial as the processed PREPBUFR observation
   sequence files have already been provided for you. However, these steps
   are provided in order to help users get started with these observations
   quickly for their own experiments.

To (again, *optionally*) reproduce the observation sequence files in the
*output* directories, you would do the following:

-  Go into your DART prep_bufr observation converter directory and
   install the PREPBUFR utilities as follows:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr
      ./install.sh

   You may need to edit the *install.sh* script to match your compiler
   and system settings.

-  Go to the
   ``$DART_DIR/observations/obs_converters/NCEP/prep_bufr/work/``
   directory and run *quickbuild.sh* to build the DART
   PREPBUFR-to-intermediate-file observation processor:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr/work
      ./quickbuild.sh

-  Download the PREPBUFR observations for your desired time. Go to the
   `NCAR/UCAR Research Data
   Archive <https://rda.ucar.edu/datasets/ds090.0/>`__ page for the
   NCEP/NCAR Global Reanalysis Products. Register on the site, click on
   the "Data Access" tab, and follow either the instructions for
   external users or NCAR internal users.

-  The downloaded *.tar* file will often be COS-blocked. If so, the file
   will appear corrupted if you attempt to untar it without converting
   the data. See the `NCAR COS-block <https://rda.ucar.edu/#!cosb>`__
   page for more information on how to strip the COS-blocking off of
   your downloaded file.

-  Untar the data in your desired directory.

-  In the ``$DART_DIR/observations/obs_converters/NCEP/prep_bufr/work``
   directory, edit the *input.nml* file. This file will control what
   observations will be used for your experiment, so the namelist
   options are worth investigating a bit here. For example, you could
   use the following:

   ::

      &prep_bufr_nml
         obs_window    = 1.0
         obs_window_cw = 1.5
         otype_use     = 120.0, 130.0, 131.0, 132.0, 133.0, 180.0
                         181.0, 182.0, 220.0, 221.0, 230.0, 231.0
                         232.0, 233.0, 242.0, 243.0, 245.0, 246.0
                         252.0, 253.0, 255.0, 280.0, 281.0, 282.0
         qctype_use    = 0,1,2,3,15
         /

   This defines an observation time window of +/- 1.0 hours, while cloud
   motion vectors will be used over a window of +/- 1.5 hours. This will
   use observation types sounding temps (120), aircraft temps (130,131),
   dropsonde temps (132), mdcars aircraft temps, marine temp (180), land
   humidity (181), ship humidity (182), rawinsonde U,V (220), pibal U,V
   (221), Aircraft U,V (230,231,232), cloudsat winds (242,243,245), GOES
   water vapor (246), sat winds (252,253,255), and ship obs (280, 281,
   282). Additionally, it will include observations with specified qc
   types only. See the
   `prepbufr <../../../observations/obs_converters/NCEP/prep_bufr/prep_bufr.html>`__
   page for more available namelist controls.

-  Within the
   ``$DART_DIR/observations/obs_converters/NCEP/prep_bufr/work``
   directory, edit the *prepbufr.csh* file and change *BUFR_dir*,
   *BUFR_idir*, *BUFR_odir*, and *BUFR_in* to match the locations and
   format of the data you downloaded. A little trial and error might be
   necessary to get these set correctly.

-  Copy over the executables from ``../exe``, and run the *prepbufr.csh*
   script for a single day at a time:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr/work
      cp ../exe/\*.x .
      ./prepbufr.csh \<year\> \<month\> \<day\>

-  Your PREPBUFR files have now been converted to an intermediate ASCII
   format. There is another observation converter to take the
   observations from this format and write them into the native DART
   format. Edit the *input.nml* namelist file in the
   *DART_DIR/observations/obs_converters/NCEP/ascii_to_obs/work*
   directory. Here is a basic example:

   ::

      &ncepobs_nml
         year       = 2017,
         month      = 4,
         day        = 27,
         tot_days   = 3,
         max_num    = 800000,
         select_obs = 0,
         ObsBase = '<path to observations>/temp_obs.',
         daily_file = .false.,
         lat1       = 15.0,
         lat2       = 60.0,
         lon1       = 270.0,
         lon2       = 330.0
         /

   Choosing "select_obs = 0" will select all the observations in the
   ASCII file. Set "ObsBase" to the directory you output the files from
   during the last step. If you wish to choose specific observations
   from the ASCII intermediate file or control other program behavior,
   there are many namelist options documented on the
   `create_real_obs <../../../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs.html>`__
   page.

-  It is now time to build *ascii_to_obs* programs. Run the following:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/ascii_to_obs/work
      ./quickbuild.sh

-  Run the *create_real_obs* program to create the DART observation
   sequence files:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/ascii_to_obs/work
      ./create_real_obs

-  The program *create_real_obs* will create observation sequence files
   with one file for each six hour window. For a cycled experiment, the
   typical approach is to put a single set of observations, associated
   with a single analysis step, into a separate directory. For example,
   within the ``output`` directory, we would create directories like
   ``2017042700``, ``2017042706``, ``2017042712``, etc. for 6-hourly
   cycling. Place the observation files in the appropriate directory to
   match the contents in the files (e.g. *obs_seq2017042706*) and rename
   as simply *obs_seq.out* (e.g. ``output/2017042706/obs_seq.out``).

-  It is helpful to also run the
   `wrf_dart_obs_preprocess <../../../models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess.html>`__
   program, which can strip away observations not in the model domain,
   perform superobservations of dense observations, increase observation
   errors near the lateral boundaries, check for surface observations
   far from the model terrain height, and other helpful pre-processing
   steps. These collectively improve system performance and simplify
   interpreting the observation space diagnostics. There are a number of
   namelist options to consider, and you must provide a *wrfinput* file
   for the program to access the analysis domain information.


Step 4: Creating the first set of adaptive inflation files
----------------------------------------------------------

In this section we describe how to create initial adaptive inflation
files. These will be used by DART to control how the ensemble is
inflated during the first assimilation cycle.

It is convenient to create initial inflation files before you start an
experiment. The initial inflation files may be created with
*fill_inflation_restart*, which was built by the *quickbuild.sh* step.
A pair of inflation files is needed for each WRF domain.

Within the ``$BASE_DIR/rundir`` directory, the *input.nml* file has some
settings that control the behavior of *fill_inflation_restart*. Within
this file there is the section:

::

   &fill_inflation_restart_nml
      write_prior_inf = .true.
      prior_inf_mean  = 1.00
      prior_inf_sd    = 0.6

      write_post_inf  = .false.
      post_inf_mean   = 1.00
      post_inf_sd     = 0.6

      input_state_files = 'wrfinput_d01'
      single_file       = .false.
      verbose           = .false.
      /

These settings write a prior inflation file with a inflation mean of 1.0
and a prior inflation standard deviation of 0.6. These are reasonable
defaults to use. The *input_state_files* variable controls which file to
use as a template. You can either modify this namelist value to point to
one of the *wrfinput_d01_XXX* files under ``$BASE_DIR/output/<DATE>``,
for any given date, or you can copy one of the files to this directory.
The actual contents of the file referenced by *input_state_files* do not
matter, as this is only used as a template for the
*fill_inflation_restart* program to write the default inflation values.
Note that the number of files specified by *input_state_files* must
match the number of domains specified in *model_nml:num_domains*, i.e.
the program needs one template for each domain. This is a
comma-separated list of strings in single 'quotes'.

After running the program, the inflation files must then be moved to the
directory expected by the *driver.csh* script.

Run the following commands with the dates for this particular tutorial:

::

   cd $BASE_DIR/rundir
   cp ../output/2017042700/wrfinput_d01_152057_0_mean ./wrfinput_d01
   ./fill_inflation_restart
   mkdir ../output/2017042700/Inflation_input
   mv input_priorinf_*.nc ../output/2017042700/Inflation_input/

Once these files are in the right place, the scripting should take care
of renaming the output from the previous cycle as the input for the next
cycle.




Step 5: Cycled analysis system
------------------------------

While the DART system provides executables to perform individual tasks
necessary for ensemble data assimilation, for large models such as WRF
that are run on a supercomputer queueing system, an additional layer of
scripts is necessary to glue all of the pieces together. A set of
scripts is provided with the tutorial tarball to provide you a starting
point for your own WRF/DART system. You will need to edit these scripts,
perhaps extensively, to run them within your particular computing
environment. If you will run on NCAR's Cheyenne environment, fewer edits
may be needed, but you should familiarize yourself with `running jobs on
Cheyenne <https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/quick-start-cheyenne>`__
if necessary. A single forecast/assimilation cycle of this tutorial can
take an hour on Cheyenne - longer if debug options are enabled or the
shared nodes are busy - shorter if more cores or a higher optimization
level is acceptable.

In this tutorial, we have previously edited the *param.csh* and other
scripts. Throughout the WRF/DART scripts, there are many options to
adjust cycling frequency, domains, ensemble size, etc., which are
available when adapting this set of scripts for your own research. To
become more famililar with this set of scripts and to eventually make
these scripts your own, we advise commenting out all the places the
script submits jobs while debugging, placing an 'exit' in the script at
each job submission step. This way you will be able to understand how
all of the pieces work together.

However, for this tutorial, we will only show you how the major
components work. The next step in our process is the main *driver.csh*
script, which expects a starting date (YYYYMMDDHH) and the full path of
the resource file as command line arguments. In this example (which uses
csh/tcsh syntax), we are also capturing the run-time output into a file
named *run.out* and the entire command will be running in the
background:

::

   cd $BASE_DIR/scripts
   ./driver.csh 2017042706 param.csh >& run.out &

*driver.csh* will - check that the input files are present (wrfinput
files, wrfbdy, observation sequence, and DART restart files), - create a
job script to run *filter* in ``$BASE_DIR/rundir``, - monitor that
expected output from *filter* is created, - submit jobs to advance the
ensemble to the next analysis time, - (simultaneously with the ensemble
advance) compute assimilation diagnostics - archive and clean up - and
continue to cycle until the final analysis time has been reached.



Step 6: Check your results
--------------------------

Once you have run the analysis system, it is time to check if things ran
well or if there are problems that need to be addressed. DART provides
analysis system diagnostics in both state and observation space.

Check to see if the analysis system actually changed the state. You
should find a file in the *$BASE_DIR/output/* directory called
*analysis_increment.nc* which is the change in the ensemble mean state
from the background to the analysis after running *filter*. Use a tool,
such as *ncview*, to look at this file. You should see spatial patterns
that look something like the meteorology of the day. These should be
places where the background (short ensemble forecast) was adjusted based
on the set of observations provided. Please become familiar with the
:doc:`Diagnostics Section <../../../guide/checking-your-assimilation>`
of the DART Documentation.  

The *driver.csh* script also ran the *diagnostics_obs.csh* which runs
the
`obs_diag <../../../assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`__
program to investigate the observation space analysis statistics. You'll
find the results of this in
``$BASE_DIR/output/<DATE>/obs_diag_output.nc``. There are many Matlab
scripts in the ``$DART_DIR/diagnostics/matlab`` directory that help
explore the effectiveness of the assimilation. Look for their examples
in the :doc:`Observation-Space
Diagnostics <../../../guide/matlab-observation-space>`
section.

The additional files enable plotting the time series of recently
assimilated observations once multiple cycles have been run. Be sure to
check that a high percentage (> 90%) of available observations were
assimilated. Low assimilation rates typically point to a problem with
the background analysis, observation quality, and/or observation error
specification which are important to address before using system results
for science.

Additional statistics can be evaluated using the converted final
observation sequence file in netcdf format from the
`obs_seq_to_netcdf <../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html>`__
tool. This file has a name like *obs_epoch_029.nc*, where the number in
the file is largest in the most recent set of observations processed.
There are Matlab tools to explore where and why the observations were
rejected. *plot_obs_netcdf.m* and *link_obs.m* are particularly useful.

If you encounter difficulties setting up, running, or evaluating the
system performance, please consider using the `GitHub
Issue <https://github.com/NCAR/DART/issues>`__ facility or feel free to
contact us at dart(at)ucar(dot)edu.

Agenda from the 22 Jan 2014 tutorial
------------------------------------

-  Introduction (Anderson) - `DART Lab
   materials <../../../guide/DART_LAB/DART_LAB.html>`__
-  WRF/DART basic building blocks (Romine)
   -`slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_building_blocks.pdf>`__
   (some material is outdated)
-  Computing environment support (Collins)
   -`slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_computing_environment.pdf>`__
-  WRF/DART application examples (Romine)
   -`slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_application_examples.pdf>`__
   (some material is outdated)
-  Observation processing (Collins)
   -`slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_observation_processing.pdf>`__
-  DART diagnostics (Hoar) - :doc:`observation diagnostics <../../../guide/matlab-observation-space>`


More Resources
--------------

-  `Check or Submit DART Issues <https://github.com/NCAR/DART/issues>`__
-  `DAReS website <ttp://dart.ucar.edu>`__
-  `Register for
   DART <https://www2.cisl.ucar.edu/software/dart/download>`__
-  `Preparing
   MATLAB <https://dart.ucar.edu/pages/Getting_Started.html#matlab>`__
   to use with DART.
-  `WRF model users page <http://www.mmm.ucar.edu/wrf/users>`__
-  Need help? e-mail dart (at) ucar (dot) edu
