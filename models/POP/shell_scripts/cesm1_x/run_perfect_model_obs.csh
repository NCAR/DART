#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Top level script to generate observations and a TRUE state.
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
#-----------------------------------------------------------------------------
#
#BXXX -b 18:00
#BSUB -J POP_OSSE
#BSUB -o POP_OSSE.%J.log
#BSUB -q economy
#BSUB -n 16
#BSUB -R "span[ptile=2]"
#BSUB -P 86850054
#BSUB -W 2:00
#BSUB -N -u ${USER}@ucar.edu
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
## PBS is used on the CGD Linux cluster 'bangkok'
## PBS is used on the CGD Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub run_filter
##
## an explanation of the most common directives follows:
## -N     Job name
## -r n   Declare job non-rerunable
## -e <arg>  filename for standard error
## -o <arg>  filename for standard out
## -q <arg>   Queue name (small, medium, long, verylong)
## -l nodes=xx:ppn=2   requests BOTH processors on the node. On both bangkok
##                     and calgary, there is no way to 'share' the processors
##                     on the node with another job, so you might as well use
##                     them both. (ppn == Processors Per Node)
##=============================================================================
#
#PBS -N POP_OSSE
#PBS -r n
#PBS -e POP_OSSE.err
#PBS -o POP_OSSE.log
#PBS -q medium
#PBS -l nodes=8:ppn=2

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv MPI         mpirun.lsf

else if ($?PBS_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by PBS
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST
   setenv MPI         mpirun

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     POP
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv MPI         csh

endif

#----------------------------------------------------------------------
# Just an echo of the job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $HOST"
echo "${JOBNAME} ($JOBID) started      at "`date`
echo

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#----------------------------------------------------------------------

setenv TMPDIR /ptmp/${user}/${JOBNAME}/job_${JOBID}

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

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = /fs/image/home/${user}/SVN/DART/models/POP
set  POPDIR = /ptmp/${user}/POP/osse
set POPFILE = `head -1 ${POPDIR}/rpointer.ocn.20.restart`

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
#-----------------------------------------------------------------------------

# executables
 ${COPY} ${DARTDIR}/work/perfect_model_obs          .
 ${COPY} ${DARTDIR}/work/dart_to_pop                .
 ${COPY} ${DARTDIR}/work/pop_to_dart                .

# shell scripts
 ${COPY} ${DARTDIR}/shell_scripts/cesm1_x/advance_model.csh .

# data files
 ${COPY} ${DARTDIR}/work/input.nml                  .
 ${COPY} ${DARTDIR}/work/obs_seq.in                 .

#-----------------------------------------------------------------------------
# Get the POP executable, control files, and data files.
# trying to use the CCSM naming conventions
#-----------------------------------------------------------------------------

 ${COPY} ${POPDIR}/pop                       .
 ${COPY} ${POPDIR}/pop_in.part1              .
 ${COPY} ${POPDIR}/pop_in.part2              .
 ${COPY} ${POPDIR}/${POPFILE}                pop.r.nc

 ${COPY} ${POPDIR}/gx3v5_tavg_contents       .
 ${COPY} ${POPDIR}/gx3v5_movie_contents      .
 ${COPY} ${POPDIR}/gx3v5_history_contents    .
 ${COPY} ${POPDIR}/gx3v5_transport_contents  .

 ${COPY} ${POPDIR}/vert_grid.gx3v5              .
 ${COPY} ${POPDIR}/horiz_grid.gx3v5.r8ieee.le   .
 ${COPY} ${POPDIR}/topography.gx3v5.i4ieee.le   .

#${COPY} ${POPDIR}/chl_mm_SeaWiFs97-01_20031205.ieeer8   .
#${COPY} ${POPDIR}/sfwf_20040517.ieeer8                  .
#${COPY} ${POPDIR}/shf_20031208.ieeer8                   .
#${COPY} ${POPDIR}/tidal_energy_gx3v5_20081021.ieeer8    .
#${COPY} ${POPDIR}/ts_PHC2_jan_20030806.ieeer8           .

#-----------------------------------------------------------------------------
# Check that everything moved OK, and the table is set.
# Convert the POP restart file to a DART ics file.
#-----------------------------------------------------------------------------

cat pop_in.part1 pop_in.part2 >! pop_in

./pop_to_dart || exit 1

${MOVE} dart_ics perfect_ics

#-----------------------------------------------------------------------------
# Run perfect_model_obs ... harvest the observations to populate obs_seq.out
# This is the 'CENTRALDIR' ... advance_model.csh will expect some things.
#-----------------------------------------------------------------------------

echo "pop.r.nc"       >! rpointer.ocn.1.restart
echo "RESTART_FMT=nc" >> rpointer.ocn.1.restart

./perfect_model_obs || exit 2

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

${MOVE} *.data *.meta         ${experiment}/POP
${MOVE} data data.cal         ${experiment}/POP
${MOVE} STD*                  ${experiment}/POP

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


