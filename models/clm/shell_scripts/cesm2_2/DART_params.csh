#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Resource file for use when running CESM (CLM specifically) and DART.
# This file has all the configuration items needed and will be copied
# into the CASEROOT directory to be used during an experiment.

# ==============================================================================
# Options defining the experiment:
#
# CASE          The value of "CASE" will be used many ways; directory and file
#               names both locally and (possibly) on the HPSS, and script names;
#               so consider its length and information content.
# compset       Defines the vertical resolution and physics packages to be used.
#               Must be a standard CESM compset; see the CESM documentation.
# resolution    Defines the horizontal resolution and dynamics; see CESM docs.
# cesmtag       The version of the CESM source code to use when building the code.
# num_instances The number of ensemble members.
#
# For list of the pre-defined component sets: ./query_config --compsets
# To create a variant compset, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ==============================================================================

setenv cesmtag        my_cesm_sandbox
setenv resolution     f09_f09_mg17
setenv compset        2000_DATM%GSWP3v1_CLM50%BGC-CROP_SICE_SOCN_MOSART_SGLC_SWAV
setenv num_instances  5

# Since this example was tested while assimilating solar induced fluorescence,
# we are using 'SIF' in the CASE. Assimilating SIF requires the use_SourceMods
# to be TRUE.

if (${num_instances} == 1) then
   setenv CASE clm5_f09_pmo_SIF
else
   setenv CASE clm5_f09_assim_e${num_instances}
endif

# ==============================================================================
# DART SourceMods for different versions of CESM are available as
# part of the DART package at ${dartroot}/models/clm/DART_SourceMods/ 
# where dartroot is the location of your DART installation
#
# SourceMods may be handled in one of two ways. If you have your own GIT clone of
# the repository, you may simply commit your changes to your GIT repo and 
# set use_SourceMods = FALSE . If you prefer to keep your changes separate 
# please put your SourceMods in a directory with 
# the following structure (which is intended to be similar to the structure 
# in the CLM distribution):
#
# ${SourceModDir}/src.clm
#                    |-- biogeochem
#                    |   `-- CNBalanceCheckMod.F90
#                    !-- biogeophys
#                    !   !-- CanopyFluxesMod.F90
#                    !   `-- PhotosynthesisMod.F90
#                    !   `-- SurfaceRadiationMod.F90
#                    `-- cpl/mct/
#                        `-- lnd_import_export.F90
#
# Description of the intent for each file:
#
# biogeochem/CNBalanceCheckMod.F90   Suppress balance checks for first restart step
#
# cpl/mct/lnd_import_export.F90      (deprecated) DS199.1 originally had some slightly
#                                    negative downward radiations that needed to be 
#                                    corrected.
#
# biogeophys/SurfaceRadiationMod.F90 Allows the use of 'PARVEG' in a history file.
#                                    Normally, only 'PARVEGLN' is output.
#
# biogeophys/CanopyFluxesMod.F90,PhotosynthesisMod.F90 calculate SIF
#
setenv dartroot       /glade/work/${USER}/DART
setenv use_SourceMods TRUE
setenv SourceModDir   ${dartroot}/models/clm/DART_SourceMods/cesm2_2_0/SourceMods

# ==============================================================================
# Directories:
# cesmdata     Location of some supporting CESM data files.
# cesmroot     Location of the CESM code base.
# caseroot     Defines the CESM case directory - where the CESM+DART
#              configuration files will be stored.  This should probably not
#              be on a fileystem that is scrubbed.
#              This script WILL DELETE any existing caseroot, so this script,
#              and other useful things should be kept elsewhere.
# rundir       Defines the location of the CESM run directory.  Will need large
#              amounts of disk space, generally on a scratch partition.
# exeroot      Defines the location of the CESM executable directory , where the
#              CESM executables will be built.  Medium amount of space
#              needed, generally on a scratch partition.
# archdir      Defines the location of the CESM short-term archive directories.
#              Files remain here until the long-term archiver moves them to 
#              permanent storage.  Requires large amounts of disk space. Should
#              not be on a scratch partition unless the long-term archiver is 
#              invoked to move these files to permanent storage.

setenv cesmdata         /glade/campaign/cesm/cesmdata/cseg/inputdata
setenv cesmroot         /glade/work/${USER}/CESM/${cesmtag}
setenv caseroot         /glade/work/${USER}/cases/${cesmtag}/${CASE}
setenv cime_output_root /glade/derecho/scratch/${USER}/${cesmtag}/${CASE}
setenv rundir           ${cime_output_root}/run
setenv exeroot          ${cime_output_root}/bld
setenv archdir          ${cime_output_root}/archive

# ==============================================================================
# Set the variables needed for the DART configuration.
# baseobsdir   Part of the directory name containing the observation sequence 
#              files to be used in the assimilation. The observations are presumed
#              to be stored in sub-directories with names built from the year and
#              month. 'baseobsdir' will be inserted into the appropriate scripts.
# ==============================================================================

setenv baseobsdir             /glade/campaign/cisl/dares/glade-p-dares-Oct2023/Observations/land
setenv pmo_input_baseobsdir   /glade/campaign/cisl/dares/glade-p-dares-Oct2023/Observations/land/pmo/input
setenv pmo_output_baseobsdir  /glade/campaign/cisl/dares/glade-p-dares-Oct2023/Observations/land/pmo/output

# ==============================================================================
# configure settings:
#
# refcase    Name of the existing reference case that this run will start from.
# refyear    The specific date/time-of-day in the reference case that this
# refmon     run will start from.  (Also see 'runtime settings' below for
# refday     start_year, start_mon, start_day and start_tod.)
# reftod
#
# stagedir   The directory location of the reference case files.
#
# startdate  The date used as the starting date for the hybrid run.
# ==============================================================================

setenv refcase      clm5.0.06_f09_80
setenv refyear      2011
setenv refmon       01
setenv refday       01
setenv reftod       00000
setenv refdate      ${refyear}-${refmon}-${refday}
setenv reftimestamp ${refyear}-${refmon}-${refday}-${reftod}

setenv stagedir /glade/campaign/cisl/dares/glade-p-dares-Oct2023/RDA_strawman/CESM_ensembles/CLM/CLM5BGC-Crop/ctsm_${reftimestamp}

# In a hybrid configuration, you can set the startdate to whatever you want.
# It does not have to match the reference (although changing the month/day seems bad).
# runtime settings:

setenv start_year    2011
setenv start_month   01
setenv start_day     01
setenv start_tod     00000
setenv startdate     ${start_year}-${start_month}-${start_day}

# ==============================================================================
# OSSE/Perfect Model experiments only.
# If there is an ensemble of CLM states to choose from, which one do you want as
# the truth? There is an argument for picking an instance that will not be part
# of the ensemble used for the assimilation experiment. 

setenv SingleInstanceRefcase FALSE
setenv TRUTHinstance 80

# ==============================================================================
# The forward operators for the flux tower obs REQUIRE that we predict the name of
# of the history file. The history file names of interest are time-tagged with the
# START of the forecast - not the restart time. The obs_def_tower_mod.f90 requires
# the stop_option to be 'nhours', and the stop_n to be accurate.
#
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in each forecast
# resubmit      How many job steps to run on continue runs (should be 0 initially)

setenv stop_option  nhours
setenv stop_n       24
setenv resubmit     0

# clm_dtime     CLM dynamical timestep (in seconds). 1800 is the default
# h1nsteps      is the number of time steps to put in a single CLM .h1. file
#               DART needs to know this and the only time it is known is during
#               this configuration step. Changing the value later has no effect.

@ clm_dtime = 1800
@ h1nsteps = $stop_n * 3600 / $clm_dtime

# ==============================================================================
# Settings for the data atmosphere

setenv stream_year_align 2011
setenv stream_year_first 2011
setenv stream_year_last  2020

# ==============================================================================
# machine-specific commands:

setenv project      P86850054
setenv machine      derecho
setenv nodes_per_instance 2
setenv number_of_threads 1

# ==============================================================================
# The FORCE  options are not optional. You may need to specify full paths
# to alternate locations that support the '-f' option.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
# ==============================================================================
set nonomatch       # suppress "rm" warnings if wildcard does not match anything

set   MOVE = 'mv -v'
set   COPY = 'cp -v --preserve=timestamps'
set   LINK = 'ln -vs'
set REMOVE = 'rm -rf' 

exit 0

