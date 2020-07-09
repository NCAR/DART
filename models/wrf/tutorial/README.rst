+-----------------------+---------------------------------------------------------------------------------------------------+
| |DART project logo|   | Jump to `DART Manhattan Documentation Main Index <../../../docs/html/Manhattan_release.html>`__   |
+-----------------------+---------------------------------------------------------------------------------------------------+

- INTRODUCTION_
- SETUP_
- `INITIAL ENSEMBLE`__
- `PREPARE OBSERVATIONS`_
- CYCLING_
- `CHECK RESULTS`_
- TUTORIAL_

WRF/DART materials for the Manhattan release.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Introduction
~~~~~~~~~~~~

In this document we will describe how to get started with your own
Weather Research and Forecasting (WRF) data assimilation system through
DART. If you are a new user to DART, we recommend you see the DART
`Getting Started <https://dart.ucar.edu/pages/Getting_Started.html>`__
page before attempting this tutorial as you will find many helpful
resources for learning the base DART configuration. This document covers
only the WRF-specific aspects of integrating with DART.

**DISCLAIMER**: We do not claim that this is a turnkey or blackbox system.
Be mentally prepared to invest a reasonable amount of time on the
learning curve. There are many outstanding research issues which have
no easy answers here. This is not a one week/grad student/naive user system.
Even after you get the code up and running, you have to be able to interpret
the results, which requires developing specific skills.
There are a lot of ways to alter how the system works -- localization,
inflation, which variables and observations are assimilated, the assimilation
window time, the model resolution, etc, etc.
This is both good and bad - you have many ways of improving your results,
but you have to take care on how you leave all the settings of these inputs.
Getting a set of scripts that runs doesn't mean the system is running well,
or producing useful results. Let the adventure begin!

This tutorial introduces a "canned" WRF/DART experiment involving an
ensemble of 50 members that will be initialized from GFS initial
conditions at 2017/04/27 00:00 UTC using a domain of the continental
United States. The data included in the tutorial lasts until 2017/04/30
18:00 UTC. During this period, there was a strong rain and wind event
that affected a large portion of the United States, causing record
rains, localized flooding, and numerous tornadoes. For more information
on the physical account of this case, see
`weather.gov <https://www.weather.gov/lot/2017Apr2930_rainfall>`__.

By default, the tutorial case will only cover 12 hours of this event
starting at 2017/04/27 00:00 UTC. The WRF model will be "spun-up" for
six hours to generate a prior distribution. An assimilation of PREPBUFR
observations will then be performed at 06:00 UTC, at which time analysis
files will be generated to begin a new ensemble forecast. The WRF model
will be advanced for 6 hours and a final assimilation cycle will be
performed at 12:00 UTC. This process could then continue in order to
investigate the strong rain and wind event.

The goals of this tutorial are to demonstrate how WRF/DART works. After
running this tutorial, you will be able to understand the major steps
involved in setting up your own data assimilation (DA) experiments.
However, you will need to do additional work before you can expect to
have a fully functional WRF/DART system, as some of the steps involved in
this tutorial (in particular, the perturbation bank and the observation
sequence files) are provided for you in order to simplify the process.
Furthermore, if you are not running on the UCAR/NCAR Cheyenne
supercomputing system, you will likely need to customize the
assimilation scripts to match the details of your particular system.

This tutorial was assembled to be compatible with ~WRF V3.9.1 and the
DART Manhattan release. Other releases of WRF may or may not be
backwards or forwards compatible with this tutorial. Check the WRF
user guide or the
`WRFHELP <http://www2.mmm.ucar.edu/wrf/users/supports/wrfhelp.html>`__
forum for WRF-specific assistance.

*DISCLAIMER*: We have provided instructions for the NCAR supercomputer
  Cheyenne, so you may need to tailor these instructions to your system
  if you are not using Cheyenne. These system-specific setup steps may
  take a good deal of effort, especially if you are unfamiliar with
  details such as MPI, NetCDF, etc. Furthermore, even after you get the
  code up and running, you will need to properly interpret your results.
  There are a lot of ways to alter how the DART system works —
  localization, inflation, which variables and observations are
  assimilated, the assimilation window time, the model resolution, etc,
  etc. This is both good and bad — you have many ways of improving your
  results, but you have to pay special attention to the settings of all
  these inputs. Getting a set of scripts that runs doesn't mean the
  system is running well, or producing useful results. Let's get started!

--------------

.. _SETUP:

Step 1: Setup
~~~~~~~~~~~~~

There are several dependencies for the executables and scripting
components. On Cheyennne, users have reported success building WRF, WPS,
WRFDA, and DART with the default module environment including Intel
compilers, MPT, and netCDF4. In addition, you'll need to load the
`nco <http://nco.sourceforge.net/>`__ and
`ncl <https://www.ncl.ucar.edu/>`__ modules to run the set of scripts
that accompany the tutorial.

If you have not already, see the
`Getting Started <https://dart.ucar.edu/pages/Getting_Started.html>`__
page to download the DART software package. Set an environment variable
*DART_DIR* to point to your base DART directory. How to do this will
depend on which shell you are using. For example, with the *tcsh*
shell, you will use

``
setenv DART_DIR <path_to_your_dart_installation>
``

while for the *bash* shell you will use

.. raw:: html

   <div class="unix">

export DART\_DIR="<path\_to\_your\_dart\_installation>"

.. raw:: html

   </div>

| 
| In either case, you will replace <path\_to\_your\_dart\_installation>
  with the actual path to your DART installation. If you are using
  another shell, refer to your shell-specific documentation on how to
  set an environment variable.

| In the same way, you will need to create a "working" directory and set
  your *BASE\_DIR* variable. Create a work directory someplace with a
  lot of free space (approximately 100 Gb are needed to run this
  tutorial). On most large systems there is a "scratch" filesystem for
  this purpose. For the rest of these instructions we will assume you
  have an environment variable called *BASE\_DIR* that points to this
  directory. For example, for *tcsh*:

.. raw:: html

   <div class="unix">

setenv BASE\_DIR <path\_to\_your\_working\_directory>

.. raw:: html

   </div>

| 
| or *bash*:

.. raw:: html

   <div class="unix">

export BASE\_DIR="<path\_to\_your\_working\_directory>"

.. raw:: html

   </div>

| 

Now that you have your two environment variables setup, download these
additional software packages (if needed):

-  The
   `WRF <http://www2.mmm.ucar.edu/wrf/users/download/get_source.html>`__
   system (WPS, real\_em build of WRF). It is assumed here that you are
   already comfortable running WRF. If not, work through the `WRF model
   tutorial <http://www2.mmm.ucar.edu/wrf/OnLineTutorial/index.htm>`__
   first before trying to link WRF and DART together.
-  The
   `WRFDA <http://www2.mmm.ucar.edu/wrf/users/wrfda/download/get_source.html>`__
   package, which is needed to generate a set of perturbed initial
   ensemble member files and also to generate perturbed boundary
   condition files. (If running this tutorial on NCAR's Cheyenne system
   this step can be skipped.)
-  The tutorial-specific additional files needed to run the examples for
   this tutorial:

   #. In this directory you will need the contents of
      *DART\_DIR/models/wrf/tutorial* from your DART code directory.

      .. raw:: html

         <div class="unix">

      cd *$BASE\_DIR*
      cp -r $DART\_DIR/models/wrf/tutorial .

      .. raw:: html

         </div>

   #. Place `this very large tar
      file <./wrf_dart_tutorial_23May2018_v3.tar.gz>`__ in your
      BASE\_DIR. CAUTION: this is an approximately 15 GB file, so you
      might be better off using 'wget' to download the file directly to
      your local system, e.g.:

      .. raw:: html

         <div class="unix">

      cd *$BASE\_DIR*
      wget http://www.image.ucar.edu/wrfdart/tutorial/wrf\_dart\_tutorial\_23May2018\_v3.tar.gz
      tar -xzvf wrf\_dart\_tutorial\_23May2018\_v3.tar.gz

      .. raw:: html

         </div>

   #. After untarring the file you should see the following directories:
      *icbc, output, perts,* and *template.* The directory names (case
      sensitive) are important, as the scripts rely on these local paths
      and file names.

Build the software packages and copy files into place:

| Copy the contents of *DART\_DIR/models/wrf/shell\_scripts* to the
  *BASE\_DIR/scripts* directory.

.. raw:: html

   <div class="unix">

cd *$BASE\_DIR*
cp -R $DART\_DIR/models/wrf/shell\_scripts ./scripts

.. raw:: html

   </div>

| 

| Copy the contents (three namelist files) of *tutorial/template* to the
  *BASE\_DIR/template* directory.

.. raw:: html

   <div class="unix">

cd *$BASE\_DIR/template*
cp ../tutorial/template/\* .

.. raw:: html

   </div>

| 

| Build the DART executables.

#. Copy the tutorial DART namelist from *template/input.nml.template* to
   *DART\_DIR/models/wrf/work/input.nml*.

   .. raw:: html

      <div class="unix">

   cd *$BASE\_DIR*
   cp template/input.nml.template $DART\_DIR/models/wrf/work/input.nml

   .. raw:: html

      </div>

#. It is assumed you have successfully configured the
   *DART\_DIR/build\_templates/mkmf.template* file for your system. If
   not, you will need to do so now. See the `Getting
   Started <https://dart.ucar.edu/pages/Getting_Started.html>`__ page
   for more detail, if necessary.
#. | Modify the DART code to use single precision reals. Most WRF/DART
     users run both the WRF model and the DART assimilation code using
     single precision floats. This is not the normal default for the
     DART code.
   | Make this code change before building the DART executables to
     compile everything with single precision reals:

   .. raw:: html

      <div class="unix">

   cd *$DART\_DIR/assimilation\_code/modules/utilities*

   .. raw:: html

      </div>

   | 
   | Edit the *types\_mod.f90* file with your favorite editor.
   | (Tip: search "real precision" to find the code block that contains
     the proper lines)
   | Comment out the following line by adding ' ! ' in the first column:

   ::

                 integer, parameter :: r8 = SELECTED_REAL_KIND(12) ! real r8
                 

   | Uncomment the following line by removing the ' ! ' from the first
     column:

   ::

                 !integer, parameter :: r8 = r4 ! alias r8 to r4
                 

#. Build the WRF/DART executables:

   .. raw:: html

      <div class="unix">

   cd *$DART\_DIR/models/wrf/work*
   ./quickbuild.csh

   .. raw:: html

      </div>

| Build (or locate an appropriate build of) WRF, WPS and WRFDA. WRF and
  WRFDA should be built with the "dmpar" option, while WPS can be built
  "serial"ly. See the WRF/WRFDA documentation for more information about
  building these packages. *NOTE*: for consistency and to avoid errors,
  you should build WRF, WPS, WRFDA, and DART with the same compiler you
  use for NetCDF. Likewise MPI should use the same compiler.

| Edit the *param.csh* script in *BASE\_DIR/scripts* with proper paths,
  info, etc. This is a script that sets variables which will be read by
  other WRF/DART scripts. There are some specific parameters for either
  the Cheyenne supercomputing system using the
  `PBS <https://www.pbsworks.com/>`__ queueing system or the older (now
  defunct) Yellowstone system which used
  `LSF <https://www.ibm.com/support/knowledgecenter/en/SSWRJV_10.1.0/lsf_welcome/lsf_welcome.html>`__.
  If you are not using Cheyenne, you may still want to use this script
  to set your queueing-system specific parameters. The following
  environment variables should be changed in the script:

+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Script variable                            | Description                                                                                                                                                                                                                                                                                                   |
+============================================+===============================================================================================================================================================================================================================================================================================================+
| module load mpt                            | The `Environment Modules <http://modules.sourceforge.net/>`__ MPI compiler to use (here the `HPE MPI <https://www.hpe.com/us/en/product-catalog/detail/pip.hpe-performance-software-message-passing-interface.1010144155.html>`__ compiler). Note that on Cheyenne the intel compiler is loaded by default.   |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| module load nco                            | The `nco <http://nco.sourceforge.net/>`__ package.                                                                                                                                                                                                                                                            |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| module load ncl/6.6.2                      | The `ncl <https://www.ncl.ucar.edu/>`__ package.                                                                                                                                                                                                                                                              |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set BASE\_DIR=<BASE DIR>                   | The root *BASE\_DIR* containing *icbc, output, perts,* etc.                                                                                                                                                                                                                                                   |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set DART\_DIR=<DART DIR>                   | The root *DART\_DIR* directory.                                                                                                                                                                                                                                                                               |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set WRF\_DM\_SRC\_DIR=<WRF DIR>            | The root directory of the WRF dmpar installation.                                                                                                                                                                                                                                                             |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set WPS\_SRC\_DIR=<WPS DIR>                | The root directory of the WPS installation.                                                                                                                                                                                                                                                                   |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set VAR\_SRC\_DIR=<WRFDA DIR>              | The root directory of the WRFDA installation.                                                                                                                                                                                                                                                                 |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set GEO\_FILES\_DIR=<WPS\_GEOG DIR>        | The root directory of the `WPS\_GEOG <https://dtcenter.org/wrf-nmm/users/OnLineTutorial/NMM/WPS/index.php>`__ files. NOTE: on Cheyenne these are available in the */glade/u/home/wrfhelp/WPS\_GEOG* directory                                                                                                 |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set GRIB\_DATA\_DIR=<GRIB DIR>             | The root directory of the GRIB data input into *ungrib.exe*. For this tutorial the grib files are included, so use *${ICBC\_DIR}/grib\_data*                                                                                                                                                                  |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set GRIB\_SRC=<Vtable.TYPE>                | Set the type of GRIB data; this will be used by *ungrib.exe* to copy the appropriate Vtable file. For the tutorial, the value should be 'GFS'.                                                                                                                                                                |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set NCAR\_GAU\_ACCOUNT=<project account>   | Set the project account to charge supercomputing hours to. See your supercomputing project administrator for more information.                                                                                                                                                                                |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| set CEMAIL=<your email address>            | Set the e-mail address used by PBS to send you information about when your job completes.                                                                                                                                                                                                                     |
+--------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

| 

| Run the *setup.csh* script to create the proper directory structure
  and move executables to proper locations.

.. raw:: html

   <div class="unix">

cd *$BASE\_DIR/scripts*
./setup.csh param.csh

.. raw:: html

   </div>

| 

So far, your *BASE\_DIR* should contain the following directories:

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
     tutorial

Your *rundir* should contain the following executables:

+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| executables:    | `advance\_time <../../../assimilation_code/programs/advance_time/advance_time.html>`__, `fill\_inflation\_restart <../../../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart.html>`__, `filter <../../../assimilation_code/programs/filter/filter.html>`__, `obs\_diag <../../../assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`__, `obs\_seq\_to\_netcdf <../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html>`__, `obs\_sequence\_tool <../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html>`__, *pert\_wrf\_bc* (no helper page), `wrf\_dart\_obs\_preprocess <../../../models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess.html>`__   |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| directories:    | *WRFIN* (empty), *WRFOUT* (empty), *WRF\_RUN* (wrf executables and support files, except namelist.input)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| scripts:        | *add\_bank\_perts.ncl*, *new\_advance\_model.csh*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| support data:   | *sampling\_error\_correction\_table.nc*                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Check to make sure your *rundir/WRF\_RUN* directory contains:

::

      da_wrfvar.exe
      wrf.exe
      real.exe
      be.dat
      contents of your WRF build run/ directory (support data files for WRF)

For this tutorial, we are providing you with a specified WRF domain. To
make your own, you would need to define your own wps namelist and use
WPS to make your own geogrid files. See the WRF site for help with
building and running those tools as needed. You would also need to get
the appropriate grib files to generate initial and boundary condition
files for the full period you plan to cycle. In this tutorial we have
provided you with geogrid files, a small set of grib files, and a
namelist to generate series of analyses for several days covering a
North American region.

Let's now look inside the *scripts* directory. You should find the
following scripts:

+--------------------------------------+--------------------------------------+
| Script name                          | Description                          |
+======================================+======================================+
| ::                                   | Add perturbations to each member.    |
|                                      |                                      |
|     add_bank_perts.ncl               |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Template for a submitted job to      |
|                                      | advance ensemble members to the next |
|     assim_advance.csh                | analysis time.                       |
+--------------------------------------+--------------------------------------+
| ::                                   | Template for submitted job to        |
|                                      | conduct the assimilation.            |
|     assimilate.csh                   |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Template for submitted job for       |
|                                      | observation specific diagnostics.    |
|     diagnostics_obs.csh              |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Primary script for running the       |
|                                      | cycled analysis system.              |
|     driver.csh                       |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Template for submitted job to        |
|                                      | advance WRF model state (on the      |
|     first_advance.csh                | first time).                         |
+--------------------------------------+--------------------------------------+
| ::                                   | Save the perturbations generated by  |
|                                      | WRFDA CV3.                           |
|     gen_pert_bank.csh                |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Generate the wrfinput and wrfbdy     |
|                                      | files.                               |
|     gen_retro_icbc.csh               |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Create the perturbed initial         |
|                                      | conditions from the WRF-VAR system.  |
|     init_ensemble_var.csh            |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Compute the mean state-space         |
|                                      | increment, which can be used for     |
|     mean_increment.ncl               | plotting.                            |
+--------------------------------------+--------------------------------------+
| ::                                   | Template for submitted job to        |
|                                      | advance the WRF model after running  |
|     new_advance_model.csh            | DART.                                |
+--------------------------------------+--------------------------------------+
| ::                                   | Contains most of the key settings to |
|                                      | run the DART system.                 |
|     param.csh                        |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Template for submitted job to        |
|                                      | prepare the initial conditions.      |
|     prep_ic.csh                      |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Run the WRF real.exe program.        |
|                                      |                                      |
|     real.csh                         |                                      |
+--------------------------------------+--------------------------------------+
| ::                                   | Create the proper directory          |
|                                      | structure and place                  |
|     setup.csh                        | executables/scripts in proper        |
|                                      | locations.                           |
+--------------------------------------+--------------------------------------+

You will need to edit these scripts to provide the paths to where you
are running the experiment, to connect up files, and to set desired
dates. Search for the string ``'set this appropriately #%%%#'`` for
locations that you need to edit.

.. raw:: html

   <div class="unix">

::

    cd $BASE_DIR/scripts

    grep -r 'set this appropriately #%%%#' .

.. raw:: html

   </div>

| Other than *param.csh*, which was covered above, make the following
  changes:

+--------------------------+--------------------------+--------------------------+
| File name                | Variable / value         | Change description       |
+==========================+==========================+==========================+
| *driver.csh*             | ::                       | Change to the final      |
|                          |                          | target date; here the    |
|                          |     set datefnl = 201704 | final date is already    |
|                          | 2712                     | set correctly for this   |
|                          |                          | tutorial.                |
+--------------------------+--------------------------+--------------------------+
| *gen\_retro\_icbc.csh*   | ::                       | This is the final date   |
|                          |                          | to create WRF            |
|                          |     set datefnl = 201704 | initial/boundary         |
|                          | 3000                     | conditions for. This is  |
|                          |                          | set to the last date     |
|                          |                          | that files are included  |
|                          |                          | in the tutorial.         |
+--------------------------+--------------------------+--------------------------+
| *gen\_retro\_icbc.csh*   | ::                       | The full path to         |
|                          |                          | *param.csh*. Change this |
|                          |     set paramfile = <ful | on the next line after   |
|                          | l param.csh path>        | the comment. While these |
|                          |                          | two files are in the     |
|                          |                          | same directory here, in  |
|                          |                          | general it is helpful to |
|                          |                          | have one *param.csh* for |
|                          |                          | each experiment.         |
+--------------------------+--------------------------+--------------------------+
| *gen\_pert\_bank.csh*    | All changes              | As the tutorial includes |
|                          |                          | a perturbation bank, you |
|                          |                          | will not need to run     |
|                          |                          | this script for the      |
|                          |                          | tutorial, so you will    |
|                          |                          | not need to change these |
|                          |                          | values. However, you     |
|                          |                          | should set appropriate   |
|                          |                          | values when you are      |
|                          |                          | ready to generate your   |
|                          |                          | own perturbation bank.   |
+--------------------------+--------------------------+--------------------------+

| 

Next, move to the *perts* directory. Here you will find 100 perturbation
files, called a "perturbation bank." For your own case, you would need
to create a perturbation bank of your own. A brief description for
running the script is available inside the comments of that file.
However, again, for this tutorial, this step has already been run for
you. The *icbc* directory contains a *geo\_em\_d01.nc* file (geo
information for our test domain), and grib files that will be used to
generate the initial and boundary condition files. The *template*
directory should contain namelists for WRF, WPS, and filter, along with
a wrfinput file that matches what will be the analysis domain. Finally,
the *output* directory contains observations within each directory name.
Template files will be placed here once created (done below), and as we
get into the cycling the output will go in these directories.

.. raw:: html

   <div class="top">

[`top <#>`__]

.. raw:: html

   </div>

--------------

.. _`INITIAL ENSEMBLE`:
 
Step 2: Initial conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

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
that will be used in the cycling. Use the script (in the scripts dir)
named *gen\_retro\_icbc.csh* to create this set of files, which will be
added to a subdirectory corresponding to the date of the run under the
"output" directory in *BASE\_DIR*. Make sure *gen\_retro\_icbc.csh* has
the appropriate path to your *param.csh* script. If the *param.csh*
script also has the correct edits for paths and you have the executables
placed in the rundir, etc., then running *gen\_retro\_icbc.csh* should
execute a series of operations to extract the grib data, run metgrid,
and then twice execute *real.exe* to generate a pair of WRF files and a
boundary file for each analysis time.

.. raw:: html

   <div class="unix">

cd *$BASE\_DIR/scripts*
./gen\_retro\_icbc.csh

.. raw:: html

   </div>

| 
| *NOTE:* ignore any *rm: No match* errors, as the script attempts to
  delete output files if they already exist, and they will not for the
  first run.

Once the script completes, inside your *output/2017042700 directory* you
should see these files:

::

       wrfbdy_d01_152057_21600_mean
       wrfinput_d01_152057_0_mean
       wrfinput_d01_152057_21600_mean

These filenames include the Gregorian dates for these files, which is
used by the dart software for time schedules. Similar files (with
different dates) should appear in all of the date directories between
the *datea* and *datef* dates set in the *gen\_retro\_icbc.csh* script.
All directories with later dates will also have an observation sequence
file *obs\_seq.out* that contains observations to be assimilated at that
time.

Next, we will execute the script to generate an initial ensemble of
states for the first analysis. For this we run the script
*init\_ensemble\_var.csh*, which takes two arguments: a date string and
the location of the *param.csh* script.

.. raw:: html

   <div class="unix">

cd *$BASE\_DIR/scripts*
./init\_ensemble\_var.csh 2017042700 param.csh

.. raw:: html

   </div>

This script generates 50 small scripts and submits them to the batch
system. It assumes a PBS batch system and the 'qsub' command for
submitting jobs. If you have a different batch system, edit this script
and look near the end. You will need to modify the lines staring with
#PBS and change 'qsub' to the right command for your system. You might
also want to modify this script to test running a single member first —
just in case you have some debugging to do.

When complete for the full ensemble, you should find 50 new files in the
directory *output/2017042700/PRIORS* with names like *prior\_d01.0001*,
*prior\_d01.0002*, etc... You may receive an e-mail to helpfully inform
you when each ensemble member has finished.

.. raw:: html

   <div class="top">

[`top <#>`__]

.. raw:: html

   </div>

--------------

Step 3: Prepare observations (optional step)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the tutorial exercise, observation sequence files are provided to
enable you to quickly get started running a test WRFDART system.

However, observation processing is critical to the success of running
DART and was covered in the `Getting
Started <https://dart.ucar.edu/pages/Getting_Started.html>`__ page. In
brief, to add your own observations to WRFDART you will need to
understand the relationship between observation definitions and
observation sequences, observation types and observation quantities, and
understand how observation converters extract observations from their
native formats into the DART specific format.

The observation sequence files that are provided in this tutorial come
from NCEP BUFR observations from the GDAS system. These observations
contain a wide array of observation types from many platforms within a
single file.

If you wanted to generate your own observation sequence files from
PREPBUFR for an experiment with WRFDART, you should follow the guidance
on the
`prepbufr <../../../observations/obs_converters/NCEP/prep_bufr/prep_bufr.html>`__
page to build the bufr conversion programs, get observation files for
the dates you plan to build an analysis for, and run the codes to
generate an observation sequence file.

For completeness, we list here how you could generate these observation
sequence files yourself. *IMPORTANT:* the following steps are **not**
necessary for the tutorial as the processed PREPBUFR observation
sequence files have already been provided for you. However, these steps
are provided in order to help users get started with these observations
quickly for their own experiments.

To (again, *optionally*) reproduce the observation sequence files in the
*output* directories, you would do the following:

-  Go into your DART prep\_bufr observation converter directory and
   install the PREPBUFR utilities as follows:

   .. raw:: html

      <div class="unix">

   cd *$DART\_DIR/observations/obs\_converters/NCEP/prep\_bufr*
   ./install.sh

   .. raw:: html

      </div>

   You may need to edit the *install.sh* script to match your compiler
   and system settings.
-  Go to the
   *DART\_DIR/observations/obs\_converters/NCEP/prep\_bufr/work/*
   directory and run *quickbuild.csh* to build the DART
   PREPBUFR-to-intermediate-file observation processor:

   .. raw:: html

      <div class="unix">

   cd *$DART\_DIR/observations/obs\_converters/NCEP/prep\_bufr/work*
   ./quickbuild.csh

   .. raw:: html

      </div>

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
-  In the *DART\_DIR/observations/obs\_converters/NCEP/prep\_bufr/work*
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
   *DART\_DIR/observations/obs\_converters/NCEP/prep\_bufr/work*
   directory, edit the *prepbufr.csh* file and change *BUFR\_dir*,
   *BUFR\_idir*, *BUFR\_odir*, and *BUFR\_in* to match the locations and
   format of the data you downloaded. A little trial and error might be
   necessary to get these set correctly.
-  Copy over the executables from *../exe*, and run the *prepbufr.csh*
   script for a single day at a time:

   .. raw:: html

      <div class="unix">

   cd *$DART\_DIR/observations/obs\_converters/NCEP/prep\_bufr/work*
   cp ../exe/\*.x . ./prepbufr.csh <year> <month> <day>

   .. raw:: html

      </div>

-  Your PREPBUFR files have now been converted to an intermediate ASCII
   format. There is another observation converter to take the
   observations from this format and write them into the native DART
   format. Edit the *input.nml* namelist file in the
   *DART\_DIR/observations/obs\_converters/NCEP/ascii\_to\_obs/work*
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

   Choosing "select\_obs = 0" will select all the observations in the
   ASCII file. Set "ObsBase" to the directory you output the files from
   during the last step. If you wish to choose specific observations
   from the ASCII intermediate file or control other program behavior,
   there are many namelist options documented on the
   `create\_real\_obs <../../../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs.html>`__
   page.

-  It is now time to build *ascii\_to\_obs* programs. Run the following:

   .. raw:: html

      <div class="unix">

   cd *$DART\_DIR/observations/obs\_converters/NCEP/ascii\_to\_obs/work*
   ./quickbuild.csh

   .. raw:: html

      </div>

-  Run the *create\_real\_obs* program to create the DART observation
   sequence files:

   .. raw:: html

      <div class="unix">

   cd *$DART\_DIR/observations/obs\_converters/NCEP/ascii\_to\_obs/work*
   ./create\_real\_obs

   .. raw:: html

      </div>

-  The program *create\_real\_obs* will create observation sequence
   files with one file for each six hour window. For a cycled
   experiment, the typical approach is to put a single set of
   observations, associated with a single analysis step, into a separate
   directory. For example, within the *output* directory, we would
   create directories like *2017042700*, *2017042706*, *2017042712*,
   etc. for 6-hourly cycling. Place the observation files in the
   appropriate directory to match the contents in the files (e.g.
   *obs\_seq2017042706*) and rename as simply *obs\_seq.out* (e.g.
   *output/2017042706/obs\_seq.out*).
-  It is helpful to also run the
   `wrf\_dart\_obs\_preprocess <../../../models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess.html>`__
   program, which can strip away observations not in the model domain,
   perform superobservations of dense observations, increase observation
   errors near the lateral boundaries, check for surface observations
   far from the model terrain height, and other helpful pre-processing
   steps. These collectively improve system performance and simplify
   interpreting the observation space diagnostics. There are a number of
   namelist options to consider, and you must provide a *wrfinput* file
   for the program to access the analysis domain information.

.. raw:: html

   <div class="top">

[`top <#>`__]

.. raw:: html

   </div>

--------------

Step 4: Creating the first set of adaptive inflation files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we describe how to create initial adaptive inflation
files. These will be used by DART to control how the ensemble is
inflated during the first assimilation cycle.

It is convenient to create initial inflation files before you start an
experiment. The initial inflation files may be created with
*fill\_inflation\_restart*, which was built by the *quickbuild.csh*
step. A pair of inflation files is needed for each WRF domain.

Within the *BASE\_DIR/rundir* directory, the *input.nml* file has some
settings that control the behavior of *fill\_inflation\_restart*. Within
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
defaults to use. The *input\_state\_files* variable controls which file
to use as a template. You can either modify this namelist value to point
to one of the *wrfinput\_d01\_XXX* files under
*BASE\_DIR/output/<DATE>*, for any given date, or you can copy one of
the files to this directory. The actual contents of the file referenced
by *input\_state\_files* do not matter, as this is only used as a
template for the *fill\_inflation\_restart* program to write the default
inflation values. Note that the number of files specified by
*input\_state\_files* must match the number of domains specified in
*model\_nml:num\_domains*, i.e. the program needs one template for each
domain. This is a comma-separated list of strings in single 'quotes'.

After running the program, the inflation files must then be moved to the
directory expected by the *driver.csh* script.

Run the following commands with the dates for this particular tutorial:

.. raw:: html

   <div class="unix">

::

    cd $BASE_DIR/rundir

    cp ../output/2017042700/wrfinput_d01_152057_0_mean ./wrfinput_d01

    ./fill_inflation_restart

    mkdir ../output/2017042700/Inflation_input

    mv input_priorinf_*.nc ../output/2017042700/Inflation_input/

.. raw:: html

   </div>

Once these files are in the right place, the scripting should take care
of renaming the output from the previous cycle as the input for the next
cycle.

.. raw:: html

   <div class="top">

[`top <#>`__]

.. raw:: html

   </div>

--------------

Step 5: Cycled analysis system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the DART system provides executables to perform individual tasks
necessary for ensemble data assimilation, for large models such as WRF
that are run on a supercomputer queueing system, an additional layer of
scripts is necessary to glue all of the pieces together. A set of
scripts is provided with the tutorial tarball to provide you a starting
point for your own WRFDART system. You will need to edit these scripts,
perhaps extensively, to run them within your particular computing
environment. If you will run on NCAR's Cheyenne environment, fewer edits
may be needed, but you should familiarize yourself with `running jobs on
Cheyenne <https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/quick-start-cheyenne>`__
if necessary.

In this tutorial, we have previously edited the *param.csh* and other
scripts. Throughout the WRFDART scripts, there are many options to
adjust cycling frequency, domains, ensemble size, etc., which are
available when adapting this set of scripts for your own research. To
become more famililar with this set of scripts and to eventually make
these scripts your own, we advise commenting out all the places the
script submits jobs while debugging, placing an 'exit' in the script at
each job submission step. This way you will be able to understand how
all of the pieces work together.

However, for this tutorial, we will only show you how the major
components work. The next step in our process is the main *driver.csh*
script, which expects a starting date as a command line argument
(YYYYMMDDHH). So you would, for this tutorial, run it as:

.. raw:: html

   <div class="unix">

cd *$BASE\_DIR/scripts*
./driver.csh 2017042706 param.csh >& run.out &

.. raw:: html

   </div>

The script will check that the input files are present (wrfinput files,
wrfbdy, observation sequence, and DART restart files), create a job
script to run filter in rundir, monitor that expected output from filter
is created, then generate job scripts for all of the model advances.
After this completes, the script will check if this is the last analysis
to determine if a new cycle is needed or not. A script is also launched
by the driver to compute some observation space diagnostics and to
convert the final observation sequence file into a netcdf format.

.. raw:: html

   <div class="top">

[`top <#>`__]

.. raw:: html

   </div>

--------------

Step 6: Check your results
~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have run the analysis system, it is time to check if things ran
well or if there are problems that need to be addressed. DART provides
analysis system diagnostics in both state and observation space.

Check to see if the analysis system actually changed the state. You
should find a file in the *output/$date/* directory called
*analysis\_increment.nc* which is the change in the ensemble mean state
from the background to the analysis after running filter. Use a tool,
such as ncview, to look at this file. You should see spatial patterns
that look something like the meteorology of the day. These should be
places where the background (short ensemble forecast) was adjusted based
on the set of observations provided.

You can also use the provided
`obs\_diag <../../../assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`__
program to investigate the observation space analysis statistics. You'll
find the results of this in output/$date/obs\_diag\_output.nc.
Additional statistics can be evaluated using the converted final
observation sequence file in netcdf format from the
`obs\_seq\_to\_netcdf <../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html>`__
tool. This file has a name like *obs\_epoch\_029.nc*, where the number
in the file is largest in the most recent set of observations processed.
The additional files enable plotting the time series of recently
assimilated observations once multiple cycles have been run. Be sure to
check that a high percentage (> 90%) of available observations were
assimilated. Low assimilation rates typically point to a problem with
the background analysis, observation quality, and/or observation error
specification which are important to address before using system results
for science.

If you encounter difficulties setting up, running, or evaluating the
system performance, please contact us at dart(at)ucar(dot)edu.

Agenda from the 22 Jan 2014 tutorial:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Introduction (Anderson) - `DART Lab
   materials <../../../docs/DART_LAB/DART_LAB.html>`__
-  WRF/DART basic building blocks (Romine) -
   `slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_building_blocks.pdf>`__
   (some material is outdated)
-  Computing environment support (Collins) -
   `slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_computing_environment.pdf>`__
-  WRF/DART application examples (Romine) -
   `slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_application_examples.pdf>`__
   (some material is outdated)
-  Observation processing (Collins) -
   `slides <https://www.image.ucar.edu/wrfdart/classic/wrf_workshop_observation_processing.pdf>`__
-  DART diagnostics (Hoar) - `observation
   diagnostics <https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`__,
   `more observation
   diagnostics <https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html>`__

Helpful links
^^^^^^^^^^^^^

-  `DAReS website <http://www.image.ucar.edu/DAReS/DART/>`__
-  `DART Manhattan release <../../../docs/html/index.html>`__
-  `Register for
   DART <https://www2.cisl.ucar.edu/software/dart/download>`__
-  `Preparing
   MATLAB <http://www.image.ucar.edu/DAReS/DART/DART2_Starting.php#matlab>`__
-  `WRF model users page <http://www.mmm.ucar.edu/wrf/users/>`__
-  Need help? e-mail dart (at) ucar (dot) edu

.. |DART project logo| image:: ../../../docs/images/Dartboard7.png
   :height: 70px
