#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#-----------------------------------------------------------------------------
# This script is designed to be run interactively.
# It sets up a "perfect_model" experiment by creating a directory that will
# be used for the duration of the experiment, i.e. 'CENTRALDIR' and
# populating it with the necessary files. 
# After this stage is set - examine the contents of CENTRALDIR to ensure
# the experiment is set up properly. If so, examine 'CENTRALDIR'/run_pmo.sh
# and then execute it if it is correct.
#-----------------------------------------------------------------------------

setenv ORIGINALDIR `pwd`
setenv JOBNAME     jules_perfect
setenv JOBID       $$

#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) started   at "`date`
echo

#----------------------------------------------------------------------
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#
# Set variables containing various directory names where we will GET things:
# DARTDIR      The location of the DART jules model directory
# JULESDIR     The location of the JULES executable
# ENSEMBLEDIR  The location of the initial ensemble of JULES files
# BASEOBSDIR   The directory containing the (empty) observation sequence file.
# CENTRALDIR   The run-time location for the experiment.
#----------------------------------------------------------------------

set nonomatch # suppress "rm" warnings if wildcard does not match anything

# The FORCE option is not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following

switch ("`hostname`")

   case ys*:
      # NCAR "yellowstone"
      set   COPY = 'cp -v --preserve=timestamps'
      set   LINK = 'ln -vs'
      set REMOVE = 'rm -fr'

      set     DARTDIR = /glade/p/work/${USER}/DART/jules/models/jules
      set    JULESDIR = ${DARTDIR}/src/jules-vn4.2/build/bin
      set ENSEMBLEDIR = /glade/p/work/thoar/DART/jules/models/jules/ensembles
      set  BASEOBSDIR = /glade/p/work/thoar/DART/jules/models/jules/work
      set  CENTRALDIR = /glade/scratch/${user}/DART/${JOBNAME}/job_${JOBID}
   breaksw

   default:
      # Bristol "lorax"
      set   COPY = 'cp -v --preserve=timestamps'
      set   LINK = 'ln -vs'
      set REMOVE = 'rm -fr'

      set     DARTDIR = /users/ar15645/DART_JULES_SVN/models/jules
      set    JULESDIR = /users/hydroeng/JULES/jules-vn4.2/build/bin
      set ENSEMBLEDIR = //users/ar15645/coupling_simulations/synthetic_test_case 
      set  BASEOBSDIR = /users/ar15645/DART_JULES_SVN/models/jules/work
      set  CENTRALDIR = /users/ar15645/run_dart_experiment/job_${JOBID}
   breaksw
endsw

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory for the experiment.
# DART literature frequently calls this 'CENTRALDIR'
#----------------------------------------------------------------------

mkdir -p ${CENTRALDIR}
cd ${CENTRALDIR}

set myname = $0          # this is the name of this script

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the JULES executable, control files, and data files.
# The advance_model.csh needs dart_to_jules (but this script does not).
#-----------------------------------------------------------------------------

foreach FILE ( ${DARTDIR}/work/input.nml   \
               ${DARTDIR}/work/jules_to_dart   \
               ${DARTDIR}/work/dart_to_jules   \
               ${DARTDIR}/work/perfect_model_obs   \
               ${DARTDIR}/shell_scripts/advance_model.csh   \
               ${DARTDIR}/shell_scripts/run_pmo.csh    \
               ${JULESDIR}/jules.exe    \
               ${ENSEMBLEDIR}/*nml )
   ${REMOVE} $FILE:t
   ${COPY}   $FILE    . || exit 1
end

chmod 755 run_pmo.csh

# JULES is highly configurable and - as far as I can tell - there are no
# 'standard' filenames. We must read the JULES namelists and make sure the
# required files are copied to the run directory.

set FILESTRING = `grep -A 4 jules_frac ancillaries.nml | grep file | sed -e "s#file=##"`
set TILEFRACTIONS = `echo $FILESTRING | sed -e "s#[']##g"`

if (  -e   ${ENSEMBLEDIR}/${TILEFRACTIONS} ) then
   ${LINK} ${ENSEMBLEDIR}/${TILEFRACTIONS} . || exit 1
else
   echo "ERROR ... no tile fraction/static data file >${ENSEMBLEDIR}/${TILEFRACTIONS}<"
   echo "ERROR ... ancillaries.nml needs it."
   exit -1
endif

set FILESTRING = `grep file drive.nml | sed -e "s#file=##"`
set FORCINGFILE = `echo $FILESTRING | sed -e "s#[',]##g"`

if (  -e   ${ENSEMBLEDIR}/${FORCINGFILE} ) then
   ${LINK} ${ENSEMBLEDIR}/${FORCINGFILE} . || exit 1
else
   echo "ERROR ... no forcing file >${ENSEMBLEDIR}/${FORCINGFILE}<"
   echo "ERROR ... drive.nml needs it."
   exit -1
endif

#-----------------------------------------------------------------------------
# Get the empty observation sequence file ... or die right away.
# This file will dictate the length of the JULES forecast.
#-----------------------------------------------------------------------------

set OBS_FILE = ${BASEOBSDIR}/obs_seq.in

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.in
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#
# DART namelist settings required:
# &perfect_model_obs_nml:  restart_in_file_name    = 'dart_ics'
# &perfect_model_obs_nml:  obs_sequence_in_name    = 'obs_seq.in'
# &perfect_model_obs_nml:  obs_sequence_out_name   = 'obs_seq.perfect'
# &perfect_model_obs_nml:  init_time_days          = -1,
# &perfect_model_obs_nml:  init_time_seconds       = -1,
# &perfect_model_obs_nml:  first_obs_days          = -1,
# &perfect_model_obs_nml:  first_obs_seconds       = -1,
# &perfect_model_obs_nml:  last_obs_days           = -1,
# &perfect_model_obs_nml:  last_obs_seconds        = -1,
# &jules_to_dart_nml:      jules_to_dart_output_file = 'dart_ics'
# &model_nml:              jules_restart_filename    = 'jules_restart.nc'
# &model_nml:              jules_history_filename    = 'jules_history.nc'
#=========================================================================

#sed -e "s#dart_ics#perfect_ics#" < input.nml.original >! input.nml

#=========================================================================
# Block 2: Convert 1 jules restart file to a DART initial conditions file.
# At the end of the block, we have a DART initial condition file  perfect_ics
#=========================================================================

# because pmo does not update the JULES state, we can simply link.
# The advance_model.csh script needs to have filenames similar to what
# naturally come out of a JULES model advance.
# Since we are running jules_to_dart
# model_mod:static_init_model() must have the 'generic' names,
# so we link the same file to two names. 

set RESTART = ${ENSEMBLEDIR}/ensemble.dump.0001.20140101.00000.nc
set  OUTPUT = ${ENSEMBLEDIR}/ensemble.hour.0001.20140101.00000.nc

if (  -e   ${RESTART} ) then
   ${LINK} ${RESTART} .
   ${LINK} ${RESTART} jules_restart.nc
else
   echo "ERROR ... no JULES restart file <${RESTART}>"
   exit -1
endif

if (  -e   ${OUTPUT} ) then
   ${LINK} ${OUTPUT} .
   ${LINK} ${OUTPUT} jules_output.nc
else
   echo "ERROR ... no JULES output file <${OUTPUT}>"
   exit -1
endif

echo "`date` -- BEGIN JULES-TO-DART"

./jules_to_dart

if ($status != 0) then
   echo "ERROR ... DART died in 'jules_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'jules_to_dart' ... ERROR"
   exit -3
endif

echo "`date` -- END JULES-TO-DART"

echo ""
echo "What to do next:"
echo "1) cd ${CENTRALDIR}"
echo "2) examine EVERYTHING"
echo "3) make sure the JULES namlists are correct for this experiment."
echo "4) make sure the DART  namlists are correct for this experiment."
echo "5) After all that - execute (or submit)"
echo "   ${CENTRALDIR}/run_pmo.csh"
echo ""

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
