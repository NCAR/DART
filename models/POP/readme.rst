POP
===

Overview
--------

This document describes the DART interface to the Parallel Ocean Program (POP).
It covers the `Development history`_ of the interface with two implementations
of POP:

- the Los Alamos National Laboratory Parallel Ocean Program (LANL POP), and
- the Community Earth System Model Parallel Ocean Program 2
  (CESM POP2; Smith et al. 2010 [1]_).

This document also provides `Detailed instructions for using DART and CESM POP2
on NSF NCAR's supercomputer`_, including information about the availability of
restart files for `Creating an initial ensemble`_ of model states and
`Observation sequence files`_ for assimilation.

Development History
-------------------

When the DART interface to POP was originally developed circa 2009-2010, the
interface worked with both the LANL POP and CESM POP2 implementations of POP.

LANL POP
~~~~~~~~

In years subsequent to the initial development of the DART interface, the
Computer, Computational, and Statistical Sciences Division at LANL transitioned
from using POP as their primary ocean model to using the Model for Prediction
Across Scales-Ocean (MPAS-Ocean). Thus it became difficult for staff in the
Data Assimilation Research Section (DAReS) at NSF NCAR to maintain access to the
`LANL POP <https://climatemodeling.science.energy.gov/projects/climate-ocean-and-sea-ice-modeling-cosim>`_
source code. As a result, LANL POP has been tested using DART's Lanai framework
but has not been tested using DART's Manhattan framework. If you intend to use
LANL POP with DART Manhattan, contact DAReS staff for assistance by emailing
dart@ucar.edu.

CESM POP2
~~~~~~~~~

The NSF NCAR implementation of POP, `CESM POP2
<https://ncar.github.io/POP/doc/build/html/index.html>`_, has been used
extensively with DART throughout multiple generations of NSF NCAR's supercomputer 
(Bluefire, Yellowstone & Cheyenne) and multiple iterations of NSF NCAR's earth
system model (CCSM4, CESM1 and CESM2). CESM POP2 is supported under DART's
Manhattan framework.

For DART's CESM POP2 interface, the CESM Interactive Ensemble facility is used
to manage the ensemble and the Flux Coupler is responsible for stopping POP2 at
the times required to perform an assimilation. CESM runs continuously and all
of the DART routines run at each assimilation time.

Detailed instructions for using DART and CESM POP2 on NSF NCAR's supercomputer
------------------------------------------------------------------------------

If you're using NSF NCAR's supercomputer, you can run the setup scripts after
making minor edits to set details that are specific to your project. The setup
scripts create a CESM case in which POP is configured using a 1° horizontal
grid, and uses the eddy parametrization of  Gent and McWilliams (1990). [2]_
The CICE model is active and atmospheric forcing is provided by the `CAM6 DART
Reanalysis <https://rda.ucar.edu/datasets/ds345.0/>`_.

The filesystem attached to NSF NCAR's supercomputer is known as the Globally
Accessible Data Environment (GLADE). All filepaths on GLADE have the structure:

.. code-block::

   /glade/*

If you aren't using NSF NCAR's supercomputer, take note of when the ``/glade/``
filepath is present in the setup scripts, since this will indicate sections
that you must alter in order to get the scripts to work on your supercomputer.
Additionally, you'll need to generate your own initial condition and
observation sequence files or you'll need to copy these files from GLADE. If
you want to copy these files from GLADE and don't have access, contact DAReS
staff by emailing dart@ucar.edu for assistance.

Summary
-------

To use DART and CESM POP2 on NSF NCAR's supercomputer, you will need to complete
the following steps.

#. Configure the scripts for your specific experiment by editing
   ``DART_params.csh``.
#. Stage your initial ensemble using ``copy_POP_JRA_restarts.py``.
#. Run the appropriate DART setup script to create and build the CESM case.

If the DART setup script runs to completion, it will print instructions to the
screen. Follow these instructions to submit your case.

Shell scripts
-------------

Since CESM requires many third-party modules in order to compile, it is often 
difficult to compile older versions of CESM because the older modules become 
unavailable. You should attempt to use the most recent setup scripts. The
`Discuss CESM bulletin board <https://bb.cgd.ucar.edu/cesm/>`_ specifies which 
releases of CESM are supported.

The setup scripts are stored in:

.. code-block::

   DART/models/POP/shell_scripts

in subdirectories that correspond releases of CESM. For example:

.. code-block::

   DART/models/POP/shell_scripts/cesm2_1

contains scripts that should be used with CESM releases 2.1.0-2.1.3.

copy_POP_JRA_restarts.py
~~~~~~~~~~~~~~~~~~~~~~~~

This script stages an intial ensemble of POP2 restart files by copying files 
from a prior experiment run by *Who Kim*. Thanks Who!

These restart files can be used as an initial ensemble of model
states. The files are kept in a directory on GLADE that is owned by the Climate
and Global Dynamics (CGD) Ocean Section:

.. code-block::

   /glade/campaign/cgd/oce/people/whokim/csm/g210.G_JRA.v14.gx1v7.01

Unless you're already a member of the CGD Ocean Section, you must be granted 
access to this directory by CISL. Use the `Service Desk
<https://servicedesk.ucar.edu/plugins/servlet/desk>`_ to request permission. If
you're unable to get permission, contact DAReS staff for assistance by emailing
dart@ucar.edu.

Filepaths beginning with ``/glade/campaign/*`` can't be accessed from NSF NCAR's 
supercomputer nodes. You must log on to NSF NCAR's data visualization computer to
copy files from ``/glade/campaign/*``.

This python script was created by *Dan Amrhein*. Thanks Dan!

+-------------------------------+-----------------------------------------------------------+
| Script name                   | Description                                               |
+===============================+===========================================================+
| ``copy_POP_JRA_restarts.py``  | This script copies restart files from the                 |
|                               | g210.G_JRA.v14.gx1v7.01 experiment that are saved in      |
|                               | campaign storage. You must be granted access to the CGD   |
|                               | Ocean Section campaign storage directory and be logged on |
|                               | to NSF NCAR's data visualization computer in order to run |
|                               | this script. The assignment of the ``stagedir`` variable  |
|                               | in this script should match the assignment of the         |
|                               | ``stagedir`` variable in ``DART_params.csh``.             |
+-------------------------------+-----------------------------------------------------------+

In order to use this script, log in to NSF NCAR's data visualization computer and
use python to run the script. For example:

.. code-block::

   $ cd DART/models/POP/shell_scripts/cesm2_1
   $ python copy_POP_JRA_restarts.py

DART_params.csh
~~~~~~~~~~~~~~~

This is the essential script you must edit to get your cases to build properly.
While you need to configure this script, you don't need to run this script.
It is run by the setup scripts.

+---------------------+-----------------------------------------------------------+
| Script name         | Description                                               |
+=====================+===========================================================+
| ``DART_params.csh`` | This script contains most, if not all, of the variables   |
|                     | that you need to set in order to build and run cases. You |
|                     | must read this file carefully and configure the variables |
|                     | to match your needs. The assignment of the ``stagedir``   |
|                     | variable in this script should match the assignment of    |
|                     | the ``stagedir`` variable in                              |
|                     | ``copy_POP_JRA_restarts.py``.                             |
+---------------------+-----------------------------------------------------------+

Setup scripts
~~~~~~~~~~~~~

These are the primary scripts used to setup CESM cases in which data
assimilation is enabled in POP2. The only variable that you might need to set
in these scripts is the ``extra_string`` variable. It is appended to the end of
the CESM case name. You can use it to differentiate experiments with the same
configuration.

+------------------------------------+--------------------------------------------+
| Script name                        | Description                                |
+====================================+============================================+
| ``setup_CESM_perfect_model.csh``   | This script creates a CESM case with a     |
|                                    | single model instance in order to run      |
|                                    | DART's ``perfect_model_obs`` program to    |
|                                    | collect observations from the model run.   |
+------------------------------------+--------------------------------------------+
| ``setup_CESM_hybrid_ensemble.csh`` | This script creates a CESM case with       |
|                                    | multiple model instances in order to run   |
|                                    | DART's ``filter`` program to complete      |
|                                    | assimilation.                              |
+------------------------------------+--------------------------------------------+

After configuring your experiment in ``DART_params.csh``, you can setup a case
by running these scripts. For example, to setup an assimilation experiment:

.. code-block::

   $ cd DART/models/POP/shell_scripts/cesm2_1
   $ ./setup_CESM_hybrid_ensemble.csh

If the setup scripts run to completion, they will print instructions that you
can follow to use CESM's case submit tool to begin a model integration.

CESM_DART_config.csh
~~~~~~~~~~~~~~~~~~~~

This script is copied by the setup scripts into the CESM case directory. It 
configures CESM to run DART.

+--------------------------+------------------------------------------------------+
| Script name              | Description                                          |
+==========================+======================================================+
| ``CESM_DART_config.csh`` | This script is copied into the CESM case directory   |
|                          | where it configures CESM to run DART.                |
+--------------------------+------------------------------------------------------+

Runtime scripts
~~~~~~~~~~~~~~~

These scripts are copied into the CESM case directory. They are called by CESM
and contain the logic to run DART's ``perfect_model_obs`` or ``filter``
programs. You shouldn't need to run these scripts directly, unless they exit 
before completion and halt a CESM integration. In this case you may need to run
the script directly to complete an assimilation in order to continue the
integration.

+-----------------------+---------------------------------------------------------+
| Script name           | Description                                             |
+=======================+=========================================================+
| ``perfect_model.csh`` | This script runs ``perfect_model_obs`` to collect       |
|                       | synthetic data in a single-instance CESM case.          |
+-----------------------+---------------------------------------------------------+
| ``assimilate.csh``    | This script runs ``filter`` to perform assimilation in  |
|                       | a multi-instance CESM case.                             |
+-----------------------+---------------------------------------------------------+

Other files needed for assimilation
-----------------------------------

Creating an initial ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Karspeck et al. (2013) [3]_ find that an ensemble of 1 January model states
selected from a multi-decade free-running integration of POP2 can be used as an
initial ensemble.

If you have access to CGD's Ocean Section directory on ``/glade/campaign`` you
can use the `copy_POP_JRA_restarts.py`_ script to stage a collection of POP
restart files from Who Kim's multi-century ``g210.G_JRA.v14.gx1v7.01``
experiment to serve as an initial ensemble. This experiment uses the JRA-55
dataset for atmospheric forcing (Tsujino et al. 2018 [4]_).

Observation sequence files
~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``setup_CESM_hybrid_ensemble.csh`` is used to create an assimilation
experiment, ``DART_params.csh`` configures the experiment to assimilate 
observation sequence files from the World Ocean Database 2013 (WOD13; Boyer et
al. 2013 [5]_).

The WOD13 dataset comprises data from 2005-01-01 to 2016-12-31 and contains the
following observation types:

+--------------------------------------+--------------------------------------+
| FLOAT_SALINITY                       | FLOAT_TEMPERATURE                    |
+--------------------------------------+--------------------------------------+
| DRIFTER_SALINITY                     | DRIFTER_TEMPERATURE                  |
+--------------------------------------+--------------------------------------+
| GLIDER_SALINITY                      | GLIDER_TEMPERATURE                   |
+--------------------------------------+--------------------------------------+
| MOORING_SALINITY                     | MOORING_TEMPERATURE                  |
+--------------------------------------+--------------------------------------+
| BOTTLE_SALINITY                      | BOTTLE_TEMPERATURE                   |
+--------------------------------------+--------------------------------------+
| CTD_SALINITY                         | CTD_TEMPERATURE                      |
+--------------------------------------+--------------------------------------+
| XCTD_SALINITY                        | XCTD_TEMPERATURE                     |
+--------------------------------------+--------------------------------------+
| APB_SALINITY                         | APB_TEMPERATURE                      |
+--------------------------------------+--------------------------------------+
| XBT_TEMPERATURE                      |                                      |
+--------------------------------------+--------------------------------------+

The W0D13 observations have already been converted into DART's observation 
sequence file format by *Fred Castruccio*. Thanks Fred! The files are stored in
the following directory on GLADE:

.. code-block::

   /glade/p/cisl/dares/Observations/WOD13

The subdirectories are formatted in ``YYYYMM`` order.

Observation sequence files converted from the World Ocean Database 2009 (WOD09;
Johnson et al. 2009 [6]_), which comprises data from 1960-01-01 to 2008-12-31,
are also stored in the following directory on GLADE:

.. code-block::

   /glade/p/cisl/dares/Observations/WOD09

These observation sequence files can be assimilated by changing the
``BASEOBSDIR`` variable in ``DART_params.csh``.

DART extracts the following variables from the POP2 restart files and adjusts
them to be consistent with the observations: ``SALT_CUR``, ``TEMP_CUR``,
``UVEL_CUR``, ``VVEL_CUR``, and ``PSURF_CUR``. 

Data atmosphere streams files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The setup scripts configure the CESM case with atmospheric forcing from the 
`CAM6 DART Reanalysis <https://rda.ucar.edu/datasets/ds345.0/>`_. The coupler 
history files from this reanalysis are referenced in
``user_datm.streams*template`` files. These ``user_datm.streams*template``
files are contained in the same directory as the setup scripts and are
configured and  copied into the CESM case directory by the setup scripts.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand, ``&``, and terminate with a slash, ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

The variables and their default values are listed here:

.. code-block:: fortran

   &model_nml
      assimilation_period_days     = -1
      assimilation_period_seconds  = -1
      model_perturbation_amplitude = 0.2
      binary_grid_file_format      = 'big_endian'
      debug                        = 0,
      model_state_variables        = 'SALT_CUR ', 'QTY_SALINITY             ', 'UPDATE',
                                     'TEMP_CUR ', 'QTY_POTENTIAL_TEMPERATURE', 'UPDATE',
                                     'UVEL_CUR ', 'QTY_U_CURRENT_COMPONENT  ', 'UPDATE',
                                     'VVEL_CUR ', 'QTY_V_CURRENT_COMPONENT  ', 'UPDATE',
                                     'PSURF_CUR', 'QTY_SEA_SURFACE_PRESSURE ', 'UPDATE'
   /

This namelist provides control over the assimilation period for the model. All
observations within (+/-) half of the assimilation period are assimilated. The
assimilation period is the minimum amount of time the model can be advanced, and
checks are performed to ensure that the assimilation window is a multiple of the
ocean model dynamical timestep.

+-------------------------------------+-------------------+------------------------------------------------------------+
| Item                                | Type              | Description                                                |
+=====================================+===================+============================================================+
| ``assimilation_period_days``        | integer           | The number of days to advance the model for each           | 
|                                     |                   | assimilation. If both ``assimilation_period_days`` and     |
|                                     |                   | ``assimilation_period_seconds`` are ≤ 0; the value of the  | 
|                                     |                   | POP namelist variables ``restart_freq`` and                |
|                                     |                   | ``restart_freq_opt`` are used to determine the             |
|                                     |                   | assimilation period.                                       |
|                                     |                   |                                                            |
|                                     |                   | *WARNING:* in the CESM framework, the ``restart_freq`` is  |
|                                     |                   | set to a value that is not useful so DART defaults to 1    |
|                                     |                   | day - even if you are using POP in the LANL framework.     |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``assimilation_period_seconds``     | integer           | In addition to ``assimilation_period_days``, the number    |
|                                     |                   | of seconds to advance the model for each assimilation.     |
|                                     |                   | Make sure you read the description of                      |
|                                     |                   | ``assimilation_period_days``.                              |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``model_perturbation_amplitude``    | real(r8)          | Reserved for future use.                                   |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``binary_grid_file_format``         | character(len=32) | The POP grid files are in a binary format. Valid values    |
|                                     |                   | are ``native``, ``big_endian``, or ``little_endian``.      |
|                                     |                   | Modern versions of Fortran allow you to specify the        |
|                                     |                   | endianness of the file you wish to read when they are      |
|                                     |                   | opened as opposed to needing to set a compiler switch or   |
|                                     |                   | environment variable.                                      |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``debug``                           | integer           | The switch to specify the run-time verbosity.              |
|                                     |                   |                                                            |
|                                     |                   | - ``0`` is as quiet as it gets.                            |
|                                     |                   | - ``> 1`` provides more run-time messages.                 |
|                                     |                   | - ``> 5`` provides ALL run-time messages.                  |
|                                     |                   |                                                            |
|                                     |                   | All values above ``0`` will also write a netCDF file of    |
|                                     |                   | the grid information and perform a grid interpolation      |
|                                     |                   | test.                                                      |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``model_state_variables``           | character(:,3)    | Strings that associate POP variables with a DART quantity  |
|                                     |                   | and whether or not to write the updated values to the      |
|                                     |                   | restart files.                                             |
|                                     |                   | These variables will be read from the POP restart          |
|                                     |                   | file and modified by the assimilation. Some (perhaps all)  |
|                                     |                   | will be used by the forward observation operators. If the  |
|                                     |                   | 3rd column is 'UPDATE', the output files will have the     |
|                                     |                   | modified (assimilated,posterior) values. If the 3rd        |
|                                     |                   | column is 'NO_COPY_BACK', that variable will not be        |
|                                     |                   | written to the restart files. **The DART diagnostic files  |
|                                     |                   | will always have the (modified) posterior values.**        |
|                                     |                   | Diagnostic variables that are useful for the calculation   |
|                                     |                   | of the forward observation operator but have no impact on  |
|                                     |                   | the forecast trajectory of the model could have a value of |
|                                     |                   | ``NO_COPY_BACK``.                                          |
+-------------------------------------+-------------------+------------------------------------------------------------+

References
----------

.. [1] Smith, R., and Coauthors, 2010: The Parallel Ocean Program (POP)
       Reference Manual Ocean Component of the Community Climate System Model
       (CCSM) and Community Earth System Model (CESM). NSF National Center for
       Atmospheric Research,
       `http://www.cesm.ucar.edu/ models/cesm1.0/pop2/doc/sci/POPRefManual.pdf <http://www.cesm.ucar.edu/ models/cesm1.0/pop2/doc/sci/POPRefManual.pdf>`_.

.. [2] Gent, P. R., and J. C. McWilliams, 1990: Isopycnal Mixing in Ocean
       Circulation Models. *Journal of Physical Oceanography*, **20**, 150–155,
       `doi:10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2 <https://doi.org/10.1175/1520-0485(1990)020\<0150:IMIOCM\>2.0.CO;2>`_.

.. [3] Karspeck, A., Yeager, S., Danabasoglu, G., Hoar, T. J., Collins, N. S.,
       Raeder, K. D., Anderson, J. L, Tribbia, J. 2013: An ensemble adjustment
       Kalman filter for the CCSM4 ocean component. *Journal of Climate*, **26**, 7392-7413,
       `doi:10.1175/JCLI-D-12-00402.1 <https://doi.org/10.1175/JCLI-D-12-00402.1>`_.

.. [4] Tsujino, H., Urakawa, S., Nakano, H., Small, R. J., Kim, W. M., Yeager,
       S. G., ... Yamazaki, D., 2018: JRA-55 based surface dataset for driving
       ocean-sea-ice models (JRA55-do). *Ocean Modelling*, **130**, 79-139,
       `doi:10.1016/j.ocemod.2018.07.002 <https://doi.org/10.1016/j.ocemod.2018.07.002>`_.

.. [5] Boyer, T.P., J. I. Antonov, O. K. Baranova, C. Coleman, H. E. Garcia,
       A. Grodsky, D. R. Johnson, R. A. Locarnini, A. V. Mishonov, T.D.
       O'Brien, C.R. Paver, J.R. Reagan, D. Seidov, I. V. Smolyar, and M. M.
       Zweng, 2013: World Ocean Database 2013, NOAA Atlas NESDIS 72, S.
       Levitus, Ed., A. Mishonov, Technical Ed.; Silver Spring, MD, 209 pp., `doi:10.7289/V5NZ85MT <http://doi.org/10.7289/V5NZ85MT>`_.

.. [6] Johnson, D.R., T.P. Boyer, H.E. Garcia, R.A. Locarnini, O.K. Baranova,
       and M.M. Zweng,  2009. World Ocean Database 2009 Documentation. Edited
       by Sydney Levitus. NODC Internal Report 20, NOAA Printing Office, Silver
       Spring, MD, 175 pp., http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html.
