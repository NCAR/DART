#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Copies and modifies the files responsible for ensemble
# initialization and integration for the special case of
# perfect_model_obs
#
# Requires the number of nodes needed for model integration and
# the path to the COAMPS directory holding the run that we will
# use as the starting point for the truth run.
######

######
# BEGIN DEFINITIONS
######

CP="cp -f"
PERL="perl -i -p -e"

PERF_WRAPPER_LOG=perfect_wrapper_log
PERF_WRAPPER_FILE=perfectwrapper
PERF_LOCK_FILE=perfect_running.lock
PERF_CONCURRENT_RUNS=1
PERF_MEMBER_DIR=perfect

######
# END DEFINITIONS
######

PERF_NODES_NEEDED=$1
COAMPS_PERFECT_DIR=$2

# advance_wrapper.csh
echo "Customizing advance_wrapper.csh for perfect_model_obs..."
${CP} advance_wrapper.csh perfect_wrapper.csh
echo "  Editing log file name..."
${PERL} "s/(\s*set diskfile = \").*(\")/\1.\/${PERF_WRAPPER_LOG}\2/" perfect_wrapper.csh
echo "  Editing wrapper file name..."
${PERL} "s/(set WRAPPERFILE = ).*/\1${PERF_WRAPPER_FILE}/" perfect_wrapper.csh
echo "  Editing lock file name..."
${PERL} "s/(set LOCKFILE = ).*/\1${PERF_LOCK_FILE}/" perfect_wrapper.csh
echo "  Editing qsub command..."
${PERL} "s/(qsub).*/\1 .\/advance_perfect.csh/" perfect_wrapper.csh
echo "  Changing name of called script in log files..."
${PERL} "s/advance_group/advance_perfect/" perfect_wrapper.csh

# advance_group.csh
echo "Customizing advance_group.csh for perfect_model_obs..."
${CP} advance_group.csh advance_perfect.csh
echo "  Editing PBS job name..."
${PERL} "s/(PBS -N).*/\1 coamps_perfect/" advance_perfect.csh
echo "  Editing PBS node count to ${PERF_NODES_NEEDED}..."
${PERL} "s/(PBS.*nodes)=\d+/\1=${PERF_NODES_NEEDED}/" advance_perfect.csh
echo "  Editing wrapper file name..."
${PERL} "s/(WRAPPERFILE = .*\/).*/\1${PERF_WRAPPER_FILE}/" advance_perfect.csh
echo "  Editing for one instance running simultaneously..."
${PERL} "s/(.*MAX_RUNNING)\s*=\s*\d+/\1 = ${PERF_CONCURRENT_RUNS}/" advance_perfect.csh
echo "  Editing log file name..."
${PERL} "s/(mktemp \").*(\.X+)/\1coamps_perfect\2/" advance_perfect.csh
echo "  Editing advance model script name..."
${PERL} "s/advance_model\.csh/advance_perfect_model\.csh/" advance_perfect.csh
echo "  Editing lock file name..."
${PERL} "s/(GROUPLOCK)\s*=\s*.*/\1 = ${PERF_LOCK_FILE}/" advance_perfect.csh

# advance_model.csh
echo "Customizing advance_model.csh for perfect_model_obs..."
${CP} advance_model.csh advance_perfect_model.csh
echo "  Editing model data directory..."
${PERL} "s/(ensdir\s*=).*/\1 ${PERF_MEMBER_DIR}/" advance_perfect_model.csh
echo "  Eliminating -pernode evaluation..."
${PERL} "s/-pernode//g" advance_perfect_model.csh

# initialize_ensemble.sh
echo "Customizing initialize_ensemble.sh for perfect_model_obs..."
${CP} initialize_ensemble.sh initialize_perfect_model.sh
echo "  Removing loop construct..."
${PERL} "s/^(\s*for\b.*)/\#\1/" initialize_perfect_model.sh
${PERL} "s/^(\s*do\b.*)/\#\1/" initialize_perfect_model.sh
${PERL} "s/^(\s*done\b.*)/\#\1/" initialize_perfect_model.sh
echo "  Editing perfect model data source directory..."
${PERL} "s|(MEMBER_DATA_DIR)=.*|\1=${COAMPS_PERFECT_DIR}|" initialize_perfect_model.sh
echo "  Editing perfect model data directory..."
${PERL} "s/(MEMBER_DIR)=.*/\1=${PERF_MEMBER_DIR}/" initialize_perfect_model.sh
echo "  Editing diagnostic output..."
${PERL} "s/ensemble member .*\"/perfect model\"/" initialize_perfect_model.sh

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

