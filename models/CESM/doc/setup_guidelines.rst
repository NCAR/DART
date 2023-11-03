Multi-Component CESM+DART Setup
===============================

CESM+DART setup overview
------------------------

If you found your way to this file without reading more basic DART help files, 
please read those first. :ref:`Getting Started <Welcome page>` is a good place to find pointers to those files. 
Then see :doc:`CESM<../readme>` for an overview of DART's interfaces to CESM.
Finally, see the ../../{your_cesm_component(s)}/readme.html documentation about
the code-level interfaces and namelist values for various CESM component models.
This document gives specific help in setting up a CESM+DART assimilation
for the first time. 

.. Warning::
   The scripts for multi-component assimilation were developed in the context 
   of DART's Lanai release (or earlier) and CESM1.  They won't work in later versions.
   The instructions below should be considered a template for setting up and running
   a multi-component assimilation, to be modified as needed.   
   Some of them reference code that may be found in $DART/models/cam or cam-old.

The overall strategy is to set up an environment where;
   * CESM is set up as a ''B'' component set configuration 
     (''fully coupled'' = active atmosphere, ocean, land, and possibly others)
   * a separate assimilation can be run using each component model interface for which there are observations.

Each CESM hindcast advances all of the active components model states,
which are then used by the several ``filter`` programs.
So you will need to build separate filters in the models/{your_CESM_component_models}/work directories.
You will also need to assemble an initial ensemble of CESM files,
which consists of restart and initial files for all of the active components.
Each filter will read a separate observation sequence file.

Assimilation set-up procedure
-----------------------------

Here is a list of steps to set up an assimilation.  
It assumes you have downloaded DART and learned how to use it with low order models. 
Some of the steps can be skipped if you have a suitable replacement, as noted.

| 

#.  Decide which component(s) you want to use as the assimilating model(s). (The rest of this list assumes that
    you're building a cam-fv assimilation, as an example.  Steps will need to be repeated for your other models.) 
#.  Look in models/{your_models}/shell_scripts to see which CESM versions are supported.
#.  CESM: locate that version on your system, or check it out from http://www.cesm.ucar.edu/models/current.html
#.  Choose the options in $dart/mkmf/mkmf.template that are best for your assimilation. These will not affect the CESM
    build, only filter.
#.  In models/cam-fv/work/input.nml, be sure to include all of your required obs_def_${platform}_mod.f90 file names in
    preprocess_nml:input_files. It's also convenient to modify the rest of input.nml to make it do what you want for the
    first assimilation cycle.   That may include creating spread in the initial ensemble by perturbing it.
    Input.nml will be copied to the $CASEROOT directory and used by assimilate.csh.
    That copy can be modified for whichever cycles will be run next.
#.  Build the DART executables using quickbuild.sh.
#.  Follow the directions in CESM/shell_scripts/\*setup\* to set up the CESM case and integrate DART into it.
    The DART team recommends a tiny ensemble to start with, to more quickly test whether everything is in order.
#.  Choose a start date for your assimilation. Choosing/creating the initial ensemble is a complicated issue.

    -  It's simpler for CAM assimilations. If you don't have an initial state and/or ensemble for this date, build a
       single instance of CESM (Fxxxx compset for cam-fv) and run it from the default Jan 1 start date until 2 weeks
       before your start date. Be sure to set the cam namelist variable inithist = 'ENDOFRUN' during the last stage, 
       so that CAM will write an "initial" file, which DART needs.
    -  For ocean and land assimilations, which cannot spin up as quickly as the atmosphere,
       creating usable initial ensemble is a more complicated process.  See those models' readme files.

#.  In the CESM run directory, create a cam-fv ensemble (virtual in the case of a single instance) 
    by linking files with instance numbers in them 
    to the restart file set (which may have no instance number) using CESM/shell_scripts/link_ens_to_single.csh.
#.  Link the other model's restart file sets into the run directory (also possibly using link_ens_to_single.csh).
#.  After convincing yourself that the CESM+DART framework is working with no_assimilate.csh, activate the assimilation
    by changing CESM's env_run.xml:DATA_ASSIMILATION_SCRIPT to use assimilate.csh.
#.  After the first hindcast+assimilation cycle finishes correctly, change the input.nml, env_run.xml and env_batch.xml
    to do additional cycle(s) without the perturbation of the initial state, and with using the restart files
    just created by the first cycle. You may also want to turn on the st_archive program. 
    Instructions are in setup_hybrid and cam-fv/work/input.nml.
#.  Finally, build a new case with the full ensemble, activate the assimilate.csh script and repeat the previous item.

Output directory
----------------

CESM's short term archiver (case.st_archive) is controlled by its ``env_archive.xml``. 
DART's setup scripts modify that file to archive DART output along with CESM's. 
(See the `<../../../guide/controlling-files-output.html>`_ for a description of DART's output).
DART's output is archived in ``$arch_dir/esp/{hist,rest,logs,...}``, where arch_dir is defined in
``setup_{hybrid,advanced}``, ``hist`` contains all of the state space and observation space output, and ``rest``
contains the inflation restart files.

The cam-XX assimilate.csh script may make a copy of its obs_seq.final files in a scratch space
($scratch/$case/Obs_seqs) which won't be removed by assimilate.csh.

Shell_scripts for building and running multi-component assimilations
--------------------------------------------------------------------

These scripts are outdated relative to Manhattan 
(path names, batch submission, long-term archiver, ...),
but can serve as a template for multi-component assimilations.

 CESM1_1_1_setup_pmo
   * set up, stage, and build a single-instance, B compset configuration of CESM. 
   * The initial state can come from any single member of a reference case.
   * Synthetic observations are harvested from the CESM model states.

 CESM1_1_1_setup_hybrid   
   * Set up, stage, and build an ensemble assimilation 
   * using a B compset configuration of CESM.
   * The initial states come from a single, multi-instance, reference case

 CESM1_1_1_setup_special
   * Same as CESM1_1_1_setup_hybrid, but the initial states for the 5 active models 
   * come from up to 5 sources:
   * The ICs source directories need to be updated.

 CESM1_1_1_setup_initial
   * Same as CESM1_1_1_setup_hybrid, but fewer comments and error checks.

 CESM1_2_1_setup_pmo
   * Same as CESM1_2_1_setup_hybrid, but for _pmo.

 CESM1_2_1_setup_hybrid
   * Same as CESM1_1_1_setup_hybrid, but updated to accommodate CESM's wave and land ice models.
   * (DART has no interfaces for those components).  Somewhat different handling of SourceMods.

 CESM_DART_config
   * Integrates DART into a pre-existing CESM case, either single- or multi-instance.
   * Typically run by or after one of the \_setup\_ scripts.

 perfect_model.csh
   * Run by the CESM $CASE.run batch job, which was created by ...setup_\ **pmo** .
   * Can call the [component]_perfect_model.csh script for each component which will be used for assimilation.  
 {cam,pop,clm}_perfect_model.csh
   * Runs perfect_model_obs_{cam,pop,clm}

 assimilate.csh
   * Run by the CESM $CASE.run batch job, which was created by ...setup_{\ **hybrid,initial,special**\ }.
   * Can call the assimilate.csh script for each component which will be used for assimilation.
   * See [component]_assimilate.csh below (which were derived from 
     $DART/models/[component]/shell_scripts/.../assimilate.csh

 cam_assimilate.csh
   * Sets up and runs filter for CAM and related observations.
   * Uses cam_to_dart and dart_to_cam, which are not used in the Manhattan release and later.

 clm_assimilate.csh 
   * similar to cam_assimilate.csh

 pop_assimilate.csh
   * similar to cam_assimilate.csh

 no_assimilate.csh
   * The script used as a placeholder in the CESM run scripts when a case is set up.

 cam_no_assimilate.csh
   * The CAM no_assimilate script needs to make an initial file available for the next CAM hindcast.

 run_perfect_model_obs.csh
   * Batch script to run perfect_model_obs for POP (only!)

 CLM_convert_restarts.csh
   * Converts 'old' CLM restart files to whatever resolution you like.

 link_ens_to_single.csh
   * Helper script to generate a virtual ensemble from a single instance (member).

 st_archive.sh
   * A CESM archiving script, modified to handle DART output files.

Helpful hints
-------------

You will probably want to use your computer resources efficiently.
In addition to the Tips and Warnings in `<../readme.html>`__,
The DART team recommends:

   + Experiment with a single instance CASE to learn the smallest number of nodes
     on which it will run reliably.  Strange andor nonreproducible errors often
     are the result of giving insufficient memory to the job.
     (node = several to dozens of central processing units which share memory
     in ways that allow very fast communication).  Build the multi-instance case
     using that number of nodes per instance.  This has 2 benefits; it minimizes queue 
     wait times, and it minimizes internode communication, which can increase exponentially 
     with the number of nodes used.
   + Carefully select the output to be saved and the archiving frequency.  
     Output from large ensemble, large model assimilations can quickly fill 
     the available disk space, resulting in an ugly ending to your job, 
     from which it is time consuming to recover; 
     discarding the partial files and keeping the output needed for evaluation
     and restarting the assimilation.
   + Evaluate the output frequently to determine whether it is worthwhile to continue.
     Looking at the model output in its gridded form can be useful, 
     but the DART team has learned that you can do a much more thorough and efficient evaluation 
     in "observation space", using 
     `obs_diag <../../../assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`_ 
     and scripts in "$DART/diagnostics/matlab" described in the 
     `Observation Space <https://dart-documentation.readthedocs.io/en/latest/guide/matlab-observation-space.html>`_


There are, no doubt, things missing from these lists, so don't struggle too long before contacting dart'at'ucar.edu.

