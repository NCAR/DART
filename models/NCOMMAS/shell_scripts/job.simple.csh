#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# job.simple.csh ... Top level script to run a single assimilation experiment.
#
# Unlike the more complex job.csh, this script only processes a single 
# observation file.  Still fairly complex; requires a raft of
# data files and most of them are in hardcoded locations.
#
# This script is designed to be run from the command line (as a single thread)
# and should only take a few seconds to a minute to complete, depending on
# the filesystem performance and data file size.
#
# The script moves the necessary files to the current directory - in DART
# nomenclature, this will be called CENTRALDIR. 
# After everything is confirmed to have been assembled, it is possible
# to edit the data, data.cal, and input.nml files for the specifics of 
# the experiment; as well as allow final configuration of a 'nodelist' file.
#
# Once the 'table is set', all that remains is to start/submit the 
# 'runme_filter' script. That script will spawn 'filter' as a 
# parallel job on the appropriate nodes; each of these tasks will 
# call a separate model_advance.csh when necessary.
#
# The central directory is where the scripts reside and where script and 
# program I/O are expected to happen.

set CENTRALDIR = `pwd`
set experiment = NCOMMASDART
alias submit 'bsub < \!*'

set myname = $0     # this is the name of this script
set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following 
set OSTYPE = `uname -s` 
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -vp'
      setenv   MOVE 'mv -fv'
      breaksw
endsw

echo " "
echo "Running $experiment on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is "`pwd`

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = /fs/image/home/${user}/SVN/DART/models/NCOMMAS
set NCOMMASDIR = /home/coral/${user}/ncommas

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
#-----------------------------------------------------------------------------

# executables
${COPY} ${DARTDIR}/work/filter                     .
${COPY} ${DARTDIR}/work/wakeup_filter              .
${COPY} ${DARTDIR}/work/dart_to_ncommas            .
${COPY} ${DARTDIR}/work/ncommas_to_dart            .

# shell scripts
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh .
${COPY} ${DARTDIR}/shell_scripts/runme_filter.csh  .

# data files
${COPY} ${DARTDIR}/work/obs_seq.out                .
${COPY} ${DARTDIR}/work/filter_ics                 .
${COPY} ${DARTDIR}/work/input.nml                  .

#-----------------------------------------------------------------------------
# Get the ncommas executable, control files, and  data files.
# trying to use the CCSM naming conventions
#-----------------------------------------------------------------------------

${COPY} ${NCOMMASDIR}/ncommas                   .
${COPY} ${NCOMMASDIR}/ncommas_in                .
${COPY} ${NCOMMASDIR}/gx3v5_tavg_contents       .
${COPY} ${NCOMMASDIR}/gx3v5_movie_contents      .
${COPY} ${NCOMMASDIR}/gx3v5_history_contents    .
${COPY} ${NCOMMASDIR}/gx3v5_transport_contents  .

${COPY} ${NCOMMASDIR}/gx3v5_vert_grid           .
${COPY} ${NCOMMASDIR}/horiz_grid_gx3.nc         .
${COPY} ${NCOMMASDIR}/topography_gx3.nc         .

${COPY} ${NCOMMASDIR}/chl_mm_SeaWiFs97-01_20031205.ieeer8   .
${COPY} ${NCOMMASDIR}/sfwf_20040517.ieeer8                  .
${COPY} ${NCOMMASDIR}/shf_20031208.ieeer8                   .
${COPY} ${NCOMMASDIR}/tidal_energy_gx3v5_20081021.ieeer8    .
${COPY} ${NCOMMASDIR}/ts_PHC2_jan_20030806.ieeer8           .

#-----------------------------------------------------------------------------
# Ensure the (output) experiment directory exists
#-----------------------------------------------------------------------------

if (-d ${experiment}) then
   echo "${experiment} already exists"
else
   echo "Making run-time directory ${experiment} ..."
   mkdir -p ${experiment}
endif

#-----------------------------------------------------------------------------
# Early exit just to check that everything moved OK, and the table is set.
# Gives us a chance to edit the local input.nml, data.cal, etc. if needed.
#-----------------------------------------------------------------------------

exit

#-----------------------------------------------------------------------------
# Runs filter which integrates the results of model advances  (async=4).
#-----------------------------------------------------------------------------

submit runme_filter

echo "Finished at "`date`

exit

#-----------------------------------------------------------------------------
# Move the output to storage after filter completes.
# At this point, all the restart,diagnostic files are in the CENTRALDIR
# and need to be moved to the 'experiment permanent' directory.
# We have had problems with some, but not all, files being moved
# correctly, so we are adding bulletproofing to check to ensure the filesystem
# has completed writing the files, etc. Sometimes we get here before
# all the files have finished being written.
#-----------------------------------------------------------------------------

echo "Listing contents of CENTRALDIR before archiving"
ls -l

${MOVE} *.data *.meta              ${experiment}/ncommas
${MOVE} data data.cal              ${experiment}/ncommas
${MOVE} STD*                       ${experiment}/ncommas

${MOVE} filter_restart*            ${experiment}/DART
${MOVE} assim_model_state_ud[1-9]* ${experiment}/DART
${MOVE} assim_model_state_ic[1-9]* ${experiment}/DART
${MOVE} analysis.nc                ${experiment}/DART
${MOVE} preassim.nc                ${experiment}/DART
${MOVE} obs_seq.final              ${experiment}/DART
${MOVE} dart_log.out               ${experiment}/DART

# Good style dictates that you save the scripts so you can see what worked.

${COPY} input.nml                  ${experiment}/DART
${COPY} *.csh                      ${experiment}/DART
${COPY} $myname                    ${experiment}/DART

ls -lrt

exit 0


