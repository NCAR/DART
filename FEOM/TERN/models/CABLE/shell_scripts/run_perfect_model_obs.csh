#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Top level script to generate observations and a TRUE state.
#
# This script only processes a single
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
#-----------------------------------------------------------------------------
#
#BSUB -J perfect
#BSUB -o perfect.%J.log
#BSUB -P P3507xxxx
#BSUB -q economy
#BSUB -n 1
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu
#
#PBS -P xa5
#PBS -N perfect
#PBS -q express
#PBS -l ncpus=16
#PBS -l mem=30gb
#PBS -l walltime=00:10:00
#PBS -j oe

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_HOSTS) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST

else if ($?PBS_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by PBS
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID:r
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST
   setenv MPICMD      "mpirun -np $NCPUS"
   setenv TASKS_PER_NODE $NCPUS

   source /etc/csh.cshrc
   module purge
   module load pbs openmpi nco netcdf


else

   #-------------------------------------------------------------------
   # interactive ...
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     cable
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST

endif

#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $MYHOST"
echo "${JOBNAME} ($JOBID) started   at "`date`
echo

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#----------------------------------------------------------------------

setenv TMPDIR /short/xa5/${user}/CABLE_DART/${JOBNAME}/job_${JOBID}

mkdir -p ${TMPDIR}
cd ${TMPDIR}

set CENTRALDIR = `pwd`
set myname = $0          # this is the name of this script

# some systems don't like the -v option to any of the following

set OSTYPE = `uname -s`
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      setenv   LINK 'ln -sf'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      setenv   LINK 'ln -sf'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -v'
      setenv   MOVE 'mv -fv'
      setenv   LINK 'ln -sfv'
      breaksw
endsw

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set     DARTDIR = /g/data/xa5/${USER}/DART/TERN/models/CABLE
set CABLEEXEDIR = /short/xa5/CABLE-EXE
# set    CABLEDIR = /short/xa5/CABLE-AUX

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/perfect_model_obs          . || exit 1
${COPY} ${DARTDIR}/work/dart_to_cable              . || exit 1
${COPY} ${DARTDIR}/work/cable_to_dart              . || exit 1

${COPY} ${DARTDIR}/shell_scripts/advance_model.csh . || exit 1

${COPY} ${DARTDIR}/observations/obs_seq.in         . || exit 1
${COPY} ${DARTDIR}/work/input.nml                  . || exit 1

# If possible, use the round-robin approach to deal out the tasks.
# This results in better memory management.

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${COPY} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" input.nml.$$ >! input.nml
      ${REMOVE} input.nml.$$
   endif
endif

#-----------------------------------------------------------------------------
# Get the CABLE executable, control files, and data files.
# The cnpipool file must be a local filename
# The cnpepool file must be a local filename
#-----------------------------------------------------------------------------

${LINK} ${CABLEEXEDIR}/cable-mpi                          .  || exit 2
${COPY} ${CABLEEXEDIR}/cable.nml                          .  || exit 2

# FIXME ... remove all the comments from the cable.nml

set STRING=`grep "casafile%cnpipool" cable.nml | sed -e "s#[='\\!]# #g"`
set CNIPOOLFILE=$STRING[$#STRING]
set STRING=`grep "casafile%cnpepool" cable.nml | sed -e "s#[='\\!]# #g"`
set CNEPOOLFILE=$STRING[$#STRING]
set STRING=`grep "filename%type" cable.nml | sed -e "s#[='\\!]# #g"`
set GRIDINFO=$STRING[$#STRING]

${LINK} ${GRIDINFO}                   CABLE_gridinfo.nc         || exit 2
${COPY} ${CABLEEXEDIR}/${CNIPOOLFILE} CABLE_poolcnp_in.0001.csv || exit 2
${COPY} ${CABLEEXEDIR}/restart_in.nc  CABLE_restart.0001.nc     || exit 2

if ( ! -e ${CABLEEXEDIR}/${CNIPOOLFILE} ) then
   echo "ERROR: cable.nml casafile%cnpipool file does not exist"
   echo "ERROR: expected ${CABLEEXEDIR}/${CNIPOOLFILE}"
   exit 2
endif

if ( ! -e ${GRIDINFO} ) then
   echo "ERROR: cable.nml filename%type gridinfo file does not exist"
   echo "ERROR: expected ${GRIDINFO}"
   exit 2
endif

#-----------------------------------------------------------------------------
# Check that everything moved OK, and the table is set.
# Convert a CABLE file 'restart_in.nc' to a DART ics file 'perfect_ics'
#-----------------------------------------------------------------------------

${LINK} CABLE_restart.0001.nc    CABLE_restart.nc

# This is the time to put the desired time in the CABLE restart file
# and make sure the forcing files are for the right year.

# Ensure that input.nml:cable_to_dart_nml:replace_cable_time
# is correct for this context.

echo '1'                       >! ex_commands
echo '/cable_to_dart_nml'      >> ex_commands
echo '/replace_cable_time'     >> ex_commands
echo ':s/\.true\./\.false\./'  >> ex_commands
echo ':wq'                     >> ex_commands

( ex input.nml < ex_commands ) >& /dev/null

./cable_to_dart || exit 3

# Safeguard against having the wrong setting in general

echo '1'                       >! ex_commands
echo '/cable_to_dart_nml'      >> ex_commands
echo '/replace_cable_time'     >> ex_commands
echo ':s/\.true\./\.false\./'  >> ex_commands
echo ':wq'                     >> ex_commands
( ex input.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

${MOVE} dart_ics perfect_ics

#-----------------------------------------------------------------------------
# Run perfect_model_obs ... harvest the observations to populate obs_seq.out
# model_mod expects a generic name // advance_model.csh expects a filename
# with the ensemble member ID tacked on - must provide both.
#-----------------------------------------------------------------------------

# ${LINK} CABLE_restart.nc     CABLE_restart.0001.nc
# ${LINK} CABLE_poolcnp_in.csv CABLE_poolcnp_in.0001.csv

./perfect_model_obs || exit 4

echo "${JOBNAME} ($JOBID) finished at "`date`

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

exit

${MOVE} DART_restart_in.0001.nc    ${EXPERIMENT}/perfect/DART_restart_in.nc
${MOVE} cable.nml                  ${EXPERIMENT}/perfect
${MOVE} obs_seq.out                ${EXPERIMENT}/perfect
${MOVE} True_State.nc              ${EXPERIMENT}/perfect
${MOVE} perfect_restart            ${EXPERIMENT}/perfect

${MOVE} dart_log.out               ${EXPERIMENT}/perfect
${MOVE} dart_log.nml               ${EXPERIMENT}/perfect
# Good style dictates that you save the scripts so you can see what worked.
${COPY} input.nml                  ${experiment}/DART
${COPY} *.csh                      ${experiment}/DART
${COPY} $myname                    ${experiment}/DART

ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

