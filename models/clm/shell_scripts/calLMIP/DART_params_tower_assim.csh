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
#               names: both locally and (possibly) on the HPSS, and script names;
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

setenv cesmtag        ctsm5.4.004
setenv resolution     CLM_USRDAT
setenv compset        HIST_DATM%1PT_CLM60%BGC_SICE_SOCN_SROF_SGLC_SWAV_SESP
setenv num_instances  80


if (${num_instances} == 1) then
   setenv CASE clm6_tower
else
   setenv CASE clm6_tower_assim_flux_deluxe2
endif

setenv dartroot       /glade/work/${USER}/DART
setenv use_SourceMods FALSE
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
setenv cesmroot         /glade/derecho/scratch/${USER}/CESM3_source/${cesmtag}
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

#setenv baseobsdir             /glade/work/bmraczka/calLMIP/observations
setenv baseobsdir             /glade/work/bmraczka/calLMIP/observations_flux

# ==============================================================================
# Parameter estimation settings:
#
# estimate_params : TRUE  = enable parameter DA (domain 4 added to DART state).
#                  FALSE = state-only DA; param staging and expansion are skipped.
#                  IMPORTANT: this value controls the shell scripting layer only.
#                  The Fortran layer reads estimate_params independently from
#                  input.nml:model_nml. Both must be set consistently or
#                  assimilate.csh will issue a warning. If they disagree the
#                  shell scripts will skip param handling while the Fortran
#                  executables (or vice versa) will attempt it -- causing errors.
#
# paramdir        : directory containing the original ensemble parameter files.
#                   These originals are NEVER modified. Working copies are staged
#                   to RUNDIR by CLM6_tower_assim at case setup.
#
# param_file_base : base filename for parameter ensemble files. Files are expected
#                   to follow the naming convention: ${param_file_base}XXXX.nc
#                   where XXXX is the zero-padded ensemble member number
#                   (e.g. calLMIP_precalibrated0001.nc ... calLMIP_precalibrated0080.nc).
#
# ==============================================================================

setenv estimate_params  TRUE
setenv paramdir         /glade/work/bmraczka/calLMIP/lhc/paramfiles
setenv param_file_base  calLMIP_precalibrated

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

# This is hardcoded into user_nl_clm
# These are here for reference
setenv refcase      clm6_tower_June
setenv refyear      2010
setenv refmon       06
setenv refday       01
setenv reftod       00000
setenv refdate      ${refyear}-${refmon}-${refday}
setenv reftimestamp ${refyear}-${refmon}-${refday}-${reftod}

setenv stagedir /glade/derecho/scratch/bmraczka/archive/clm6_transient_precal/rest/${reftimestamp}/ctsm_${reftimestamp}


# In a hybrid configuration, you can set the startdate to whatever you want.
# It does not have to match the reference (although changing the month/day seems bad).
# runtime settings:

setenv start_year    2010
setenv start_month   06
setenv start_day     01
setenv start_tod     00000
setenv startdate     ${start_year}-${start_month}-${start_day}

# ==============================================================================
# The forward operators for the flux tower obs REQUIRE that we predict the name of
# of the history file. The history file names of interest are time-tagged with the
# START of the forecast - not the restart time. The obs_def_tower_mod.f90 requires
# the stop_option to be 'nhours', and the stop_n to be accurate.
#
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in each forecast
# resubmit      How many job steps to run on continue runs (should be 0 initially)

#152 days from Jan 1 to Jun 1

setenv stop_option  ndays
setenv stop_n       1
setenv calendar     GREGORIAN
setenv resubmit     1600

# clm_dtime     CLM dynamical timestep (in seconds). 1800 is the default
# h1nsteps      is the number of time steps to put in a single CLM .h1. file
#               DART needs to know this and the only time it is known is during
#               this configuration step. Changing the value later has no effect.

@ clm_dtime = 1800
@ h1nsteps = $stop_n * 3600 / $clm_dtime

# ==============================================================================
# Settings for the data atmosphere

#setenv stream_year_align 1997
#setenv stream_year_first 1997
#setenv stream_year_last  2014


setenv stream_year_align 2010
setenv stream_year_first 2010
setenv stream_year_last  2014

# ==============================================================================
# machine-specific commands:

setenv project      P86850054
setenv machine      derecho
# Enforce only 1 node for site level runs
setenv nodes_per_instance 1
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

