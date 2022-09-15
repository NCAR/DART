#!/bin/csh
#
# Copyright 2020 University Corporation for Atmospheric Research
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License.
# Please view the License at http://www.apache.org/licenses/LICENSE-2.0
#
# ==============================================================================
#
# Resource file for use when running CESM (POP specifically) and DART.
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
# case_string   A unique string used to name a given case.
# resolution    Defines the horizontal resolution and dynamics; see CESM docs.
# cesmtag       The version of the CESM source code to use when building the code.
# cesmtagmajor  The major release version of CESM, used to identify the
#               appropriate CESM directory.
# 
# To create a variant compset, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ==============================================================================

setenv cesmtag        cesm2_1_1
setenv resolution     f09_g17
setenv compset        G
setenv case_string    GNYF
setenv cesmtagmajor   `echo ${cesmtag} | cut -c1-7`

# ==============================================================================
# There are SourceMods that enable POP to recompute the barotropic velocity
# to prevent the barotropic solver from crashing.
# SourceMods may be handled in one of two ways. If you have your own git clone of
# the repository, you may simply commit your changes to your git repo and 
# set use_SourceMods = FALSE . If you prefer to keep your changes separate 
# (as was required under svn), please put your SourceMods in a directory with 
# the following structure (which is intended to be similar to the structure 
# in the CLM distribution):
#
# ${SourceModDir}/src.pop
#                |-- forcing.F90
#                |-- initial.F90
#                |-- overflows.F90
#                `-- restart.F90

setenv use_SourceMods TRUE
setenv SourceModDir   ~/${cesmtag}/SourceMods

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

setenv cesmdata         /glade/p/cesmdata/cseg/inputdata
setenv cesmroot         /glade/p/cesm/releases/$cesmtag

# ==============================================================================
# Set the variables needed for the DART configuration.
# DARTROOT     Location of the root of _your_ DART installation
# BASEOBSDIR   Part of the directory name containing the observation sequence 
#              files to be used in the assimilation. The observations are presumed
#              to be stored in sub-directories with names built from the year and
#              month. 'BASEOBSDIR' will be inserted into the appropriate scripts.
# ==============================================================================

setenv DARTROOT               /glade/work/${USER}/git/DART
setenv BASEOBSDIR             /glade/p/cisl/dares/Observations/WOD13

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

setenv refcase      g210.G_JRA.v14.gx1v7.01
setenv refyear      2010
setenv refmon       01
setenv refday       01
setenv reftod       00000
setenv refdate      ${refyear}-${refmon}-${refday}
setenv reftimestamp ${refyear}-${refmon}-${refday}-${reftod}

setenv stagedir /glade/scratch/${USER}/${refcase}/rest/${reftimestamp}

# In a hybrid configuration, you can set the startdate to whatever you want.
# It does not have to match the reference (although changing the month/day seems bad).
# runtime settings:

setenv start_year    2014
setenv start_month   01
setenv start_day     01
setenv start_tod     00000
setenv startdate     ${start_year}-${start_month}-${start_day}

# ==============================================================================
# Make sure the CESM directories exist.
# VAR is the shell variable name, DIR is the value
# ==============================================================================
 
foreach VAR ( cesmroot DARTROOT stagedir )
   set DIR = `eval echo \${$VAR}`
   if ( ! -d $DIR ) then
      echo "ERROR: directory '$DIR' not found"
      echo " In the setup script check the setting of: $VAR "
      exit -1
   endif
end

# ==============================================================================
# OSSE/Perfect Model experiments only.
# If there is an ensemble of POP states to choose from, which one do you want as
# the truth? There is an argument for picking an instance that will not be part
# of the ensemble used for the assimilation experiment. 

setenv SingleInstanceRefcase 0
setenv TRUTHinstance 23

# ==============================================================================
#
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in each forecast
# resubmit      How many job steps to run on continue runs (should be 0 initially)

setenv resubmit     0
setenv stop_option  ndays
setenv first_STOP_N 3
setenv STOP_N       1
setenv resubmit     0

# ==============================================================================
# Settings for the data atmosphere

setenv stream_year_align 2011
setenv stream_year_first 2011
setenv stream_year_last  2015

setenv short_term_archiver off

# ==============================================================================
# job settings and machine-specific commands:
#
# mach            Computer name
# ACCOUNT         Project that hours will be charged to
# queue           can be changed during a series by changing the ${CASE}.run
# timewall        can be changed during a series by changing the ${CASE}.run
# ==============================================================================

setenv mach         cheyenne
setenv ACCOUNT      P########
setenv queue        regular
setenv timewall     00:15:00

# ==============================================================================
# by setting the values above you should be able to execute this script and
# have it run.  however, for running a real experiment there are still many
# settings below this point - e.g. component namelists, history file options,
# the processor layout, xml file options, etc - that you will almost certainly
# want to change before doing a real science run.
# ==============================================================================

# The CESM compile step takes enough resource that cheyenne requires a wrapper
# If your platform does not have this restriction, set BUILD_WRAPPER to '' 
# setenv BUILD_WRAPPER ''
setenv BUILD_WRAPPER "qcmd -q share -l select=1 -A ${ACCOUNT} --"

# Command to run an MPI executable on your machine
setenv LAUNCHCMD 'mpiexec_mpt dplace -s 1'

# ==============================================================================
#
# The FORCE  options are not optional. You may need to specify full paths
# to alternate locations that support the '-f' option.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
# ==============================================================================
set  nonomatch       # suppress "rm" warnings if wildcard does not match anything

set  MOVE = 'mv -v'
set  COPY = 'cp -v --preserve=timestamps'
set  LINK = 'ln -vfs'
set  REMOVE = 'rm -rf' 

exit 0

