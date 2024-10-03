.. _icepack:

Icepack
=======

Overview
--------

DART interface modules for Icepack, the column physics of the sea ice model CICE (`https://github.com/cice-consortium/Icepack <https://github.com/cice-consortium/Icepack>`_). Icepack is maintained by the CICE Consortium.

The column physics package of the sea ice model CICE, “Icepack”, is maintained by the CICE Consortium. A large portion of the physics in sea ice models can be described in a vertical column, without reference to neighboring grid cells. This code includes several options for simulating sea ice thermodynamics, mechanical redistribution (ridging) and associated area and thickness changes. In addition, the model supports a number of tracers, including thickness, enthalpy, ice age, first-year ice area, deformed ice area and volume, melt ponds, and biogeochemistry.

More information about the model and instructions on how to run Icepack can be found in the `Icepack documentation: <https://cice-consortium-icepack.readthedocs.io/en/main/index.html>`_

This model is run as a separate executable from DART, and this means that you must use scripts to alternate the model and DART program execution and allow for the progression of the assimilation through multiple time windows. These scripts will be provided by DART, but they are currently still in progress. 

The assimilation process can be easily executed within a single assimilation window, however. There is a test case available in ``/glade/work/masmith/test_cases/icepack_test`` that contains the necessary input files to run filter, the main program in DART that performs the assimilation. There is a README file in this directory to give more details on the specifics of the test case. If you do not have access to the NSF NCAR Derecho Supercomputer, the reach out to the DAReS team at ``masmith@ucar.edu`` and we will provide you with the test case.

The steps to run this example are as follows:

1.  | ``cd /DART/models/icepack/work``
    | Navigate to the work directory for Icepack

2.  | ``./quickbuild.sh`` (or ``quickbuild.sh nompi`` if you are not building with mpi)
    | Builds all DART executables 

3.  | ``cp -r /glade/work/masmith/test_cases/icepack_test ..``
    | Copy the test directory from the directory stated above to your Icepack directory

4.  | ``cd ../icepack_test``
    | Navigate to the test directory

5.  | In ``work/input.nml``, set ``perturb_from_single_instance = .true.`` in the
      ``&filter_nml``
    | This setting causes filter to perturb a single restart file to generate an
      ensemble

6.  | ``./filter``
    | Runs the assimilation program, resulting in four main output files:
    |    ``analysis_mean.nc`` and ``analysis_sd`` - the mean and standard deviation of the state of all ensemble members after the assimilation
    |    ``obs_seq.final`` - the ensemble members' estimate of the observations.
    |    ``dart_log.out`` - detailed log file for the execution of filter

Namelist
--------

.. code-block:: fortran

&model_nml
    model_perturbation_amplitude = 2e-05
    debug = 1
    model_state_variables = 'aicen', 'QTY_SEAICE_CONCENTR', 'UPDATE', 'vicen',
                            'QTY_SEAICE_VOLUME', 'UPDATE', 'vsnon', 'QTY_SEAICE_SNOWVOLUME',
                            'UPDATE'
    grid_oi = 3
/

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+------------------------------+---------------+---------------------------------+
| Item                         | Type          | Description                     |
+==============================+===============+=================================+
| model_perturbation_amplitude | real(r8)      | Perturbation amplitude          |
+------------------------------+---------------+---------------------------------+
| debug                        | integer       | When set to 0, debug statements |
|                              |               | are not printed. Higher numbers |
|                              |               | mean more debug reporting.      |
+------------------------------+---------------+---------------------------------+
| model_state_variables        | character(*)  | List of model state variables   |
+------------------------------+---------------+---------------------------------+
| grid_oi                      | integer       | Specifies a constant to be used |
|                              |               | as the value for the first      |
|                              |               | dimension in calls to           |
|                              |               | get_dart_vector_index           |
+------------------------------+---------------+---------------------------------+

References
~~~~~~~~~~

.. [1] Hunke E et al. 2018 CICE-Consortium/Icepack version 1.4.1. `doi:10.5281/zenodo.11223808 <https://doi.org/10.5281/zenodo.11223808>`_
