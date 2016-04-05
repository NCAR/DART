#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Top level script to run a single assimilation experiment.
#
# Unlike the more complex job.csh, this script only processes a single 
# observation file.  Still fairly complex; requires a raft of
# data files and most of them are in hardcoded locations.
#
# You need to know which of several batch systems you are using.  The most
# common one is LSF.   PBS is also common.  (POE is another but is
# not supported directly by this script.  It is not recommended that you have a
# parallel cluster without a batch system (it schedules which nodes are assigned
# to which processes) but it is possible to run that way -- you have to do
# more work to get the information about which nodes are involved to the 
# parallel tasks -- but anyway, there is a section below that uses ssh and no
# batch.
#
# How to submit this job:
#  1. Look at the #BSUB or #PBS sections below and adjust any of the parameters
#     on your cluster.  Queue names are very system specific; some systems 
#     require wall-clock limits; some require an explicit charge code.
#  2. Submit this script to the queue:
#        LSF:   bsub < job.simple.csh
#        PBS:   qsub job.simple.csh
#       NONE:   job.simple.csh
#
# The script moves the necessary files to the current directory and then
# starts 'filter' as a parallel job on all nodes; each of these tasks will 
# call some a separate model_advance.csh when necessary.
#
# The central directory is where the scripts reside and where script and 
# program I/O are expected to happen.
# 
# 
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system 
# LSF is used on the IBM   Linux cluster 'lightning'
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluevista'
# The queues on lightning and bluevista are supposed to be similar.
#
# the normal way to submit to the queue is:    bsub < job.simple.csh
#
# an explanation of the most common directives follows:
# -J Job name
# -o STDOUT filename
# -e STDERR filename
# -P      account
# -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
# -n number of processors  (really)
# -W hr:mn   max wallclock time (required on some systems)
##=============================================================================
#BSUB -J DARTCAM
#BSUB -o DARTCAM.%J.log
#BSUB -q regular
#BSUB -n 1
#
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system 
## PBS is used on the CGD   Linux cluster 'bangkok'
## PBS is used on the CGD   Linux cluster 'calgary'
##
## the normal way to submit to the queue is:    qsub job.simple.csh
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
##                     them both.  (ppn == Processors Per Node)
##=============================================================================
#PBS -N DARTCAM
#PBS -r n
#PBS -e DARTCAM.err
#PBS -o DARTCAM.log
#PBS -q medium
#PBS -l nodes=2:ppn=2

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder 
# of the script.

if ($?LS_SUBCWD) then

   # LSF has a list of processors already in a variable (LSB_HOSTS)

   set CENTRALDIR = $LS_SUBCWD
   set JOBNAME = $LSB_JOBNAME
   alias submit 'mpirun.lsf \!*'
   
else if ($?PBS_O_WORKDIR) then

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   alias submit 'mpirun \!*'

else if ($?OCOTILLO_NODEFILE) then

   # no queueing system.  Many possible options here.  If you need a list
   # of processors in a file, this makes one. For example ... 
   # setenv NODEFILE  my_favorite_processors.list
   # echo "node1"  > $NODEFILE
   # echo "node5" >> $NODEFILE
   # echo "node7" >> $NODEFILE
   # echo "node3" >> $NODEFILE

   set CENTRALDIR = `pwd`
   set JOBNAME = DARTCAM
   # i think this is what we want, but csh will not let you do multiline
   # executions; this argues for using ksh (line 2 below)...  (and maybe
   # it needs a cd as well?)
   #alias submit 'foreach i ($OCOTILLO_NODEFILE) ; ssh $i csh \!* ; end'
   #alias submit='for i in $OCOTILLO_NODEFILE ; do ssh $i (cd $CENTRALDIR; csh $*) ; done'
   alias submit 'csh \!*'
   
else

   # interactive

   set CENTRALDIR = `pwd`
   set JOBNAME = DARTCAM
   alias submit 'csh \!*'
   
endif

set myname = $0     # this is the name of this script
set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# Set the experiment name.

set experiment = CAM1X

cd ${CENTRALDIR}

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
echo "Running $JOBNAME on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is "`pwd`

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

# where the DART executables are:
set DARTDIR = /home/${user}/DART

# where the CAM executables are, and what the case build name was:
set CAMsrc = /home/${user}/Cam3/cam3.1/models/atm/cam/bld
set CAMcase = T21-O2

# where the CAM data files are:
set CAMdata = /home/${user}/CAMDATA

# where a precomputed set of initial condition files are
set DARTics  = /home/${user}/CAM_init/T21x80/03-01-01/DART_lunes
set CAMics   = /home/${user}/CAM_init/T21x80/03-01-01/CAM/caminput_
set CLMics   = /home/${user}/CAM_init/T21x80/03-01-01/CLM/clminput_

#-----------------------------------------------------------------------------
# Get the DARTCAM executables and scripts
#-----------------------------------------------------------------------------

set DARTCAMDIR = ${DARTDIR}/models/cam
${COPY} ${DARTCAMDIR}/work/filter                     .
${COPY} ${DARTCAMDIR}/work/dart_to_cam                .
${COPY} ${DARTCAMDIR}/work/cam_to_dart                .
${COPY} ${DARTCAMDIR}/shell_scripts/advance_model.csh .
${COPY} ${DARTCAMDIR}/shell_scripts/run-pc.csh        .

${COPY} ${CAMdata}/input.nml   .
${COPY} ${CAMdata}/obs_seq.out .

#-----------------------------------------------------------------------------
# Determine the number of ensemble members from input.nml,
# It may exist in more than one place - we will use the first instance.
# Parse out the filter_nml string and use the next hunk of lines.
# ditto for the advance command
#-----------------------------------------------------------------------------

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set  ADVANCESTRING = `grep -A 42 filter_nml input.nml | grep adv_ens_command`
set  ensemble_size = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`
set        ADV_CMD = `echo  $ADVANCESTRING[3] | sed -e 's#,##' -e 's#"##g'`

set num_ens = $ensemble_size


echo "There are ${num_ens} ensemble members."

#-----------------------------------------------------------------------------
# Get the necessary data files -- this is the hard part.
# This script does not involve 'cold starting' CAM, nor spinning up DART.
# The DARTics directory has one initial conditions file for
# each ensemble member. We need one for each ...
# The input.nml has a restart_in_file_name of 'filter_ic_old'
# which must match the filename here. 
# Because that same namelist has 'single_restart_file_in' as .false.,
# the restart_in_file_name gets an ensemble member number appended to it.
#-----------------------------------------------------------------------------

# This copies the initial conditions for the correct number
# of ensemble members, and renames them with the ensemble number
# as a file extension.  in all cases the extension is 4 digits,
# so use 'printf' to add leading 0s to the ensemble number.


set n = 1
while($n <= ${num_ens})
   set from = `printf "%s.%04d" ${DARTics}/filter_ic $n`
   set to   = `printf "%s.%04d" filter_ic_old        $n`
   echo copying $from to $to
   ${COPY} $from $to
   @ n++ 
end

${COPY} ${CAMdata}/namelistin                      .
${COPY} ${CAMdata}/caminput.nc                     .
${COPY} ${CAMdata}/clminput.nc                     .

#-----------------------------------------------------------------------------
# T21
# The CAMsrc directory is MORE than just the location of the executable.
# There are more support widgets expected in the directory tree.
#-----------------------------------------------------------------------------

set CAMexe = ${CAMsrc}/${CAMcase}


#-----------------------------------------------------------------------------
# Ensure the (output) experiment directory exists
# All the  CAM-related files will get put in ${experiment}/CAM
# All the  CLM-related files will get put in ${experiment}/CLM
# All the DART-related files will get put in ${experiment}/DART
#-----------------------------------------------------------------------------

if (-d ${experiment}) then
   echo "${experiment} already exists"
else
   echo "Making run-time directory ${experiment} ..."
   mkdir -p ${experiment}
endif
mkdir -p ${experiment}/{CLM,CAM,DART}

#-----------------------------------------------------------------------------
# This is where I should check to make sure all the required files exist.
#-----------------------------------------------------------------------------

if (! -e namelistin ) then
   echo "ERROR ... need a namelistin file."
   exit 99
endif

#-----------------------------------------------------------------------------
# get name of file containing PHIS from the CAM namelist.  This will be used by
# static_init_model to read in the PHIS field, which is used for height obs.
#-----------------------------------------------------------------------------
   grep bnd_topo namelistin >! topo_file
   set  STRING = "1,$ s#'##g"
   set ensstring = `sed -e "$STRING" topo_file`
   set topo_name = $ensstring[3]
   cp $topo_name cam_phis.nc

#-----------------------------------------------------------------------------
# Some information about CAM must be made available to advance_model.csh
# 'casemodel' is required (by advance_model.csh) to be in the Central directory
#-----------------------------------------------------------------------------

echo "${experiment} ${CAMexe} ${CAMics} ${CLMics}" >! casemodel

#-----------------------------------------------------------------------------
# Runs filter which integrates the results of model advances  (async=2).
#
# A 20 member ensemble @ T21 can take anywhere between 10-30 minutes.
#-----------------------------------------------------------------------------

submit filter

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

echo "Finished at "`date`

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


${MOVE} cam_out_temp[1-9]*         ${experiment}/CAM
${MOVE} caminput_[1-9]*.nc         ${experiment}/CAM

${MOVE} clminput_[1-9]*.nc         ${experiment}/CLM

${MOVE} filter_ic_old*             ${experiment}/DART
${MOVE} filter_ic_new*             ${experiment}/DART
${MOVE} assim_model_state_ud[1-9]* ${experiment}/DART
${MOVE} assim_model_state_ic[1-9]* ${experiment}/DART
#${MOVE} inflate_ic_new             ${experiment}/DART

${MOVE} Posterior_Diag.nc          ${experiment}/DART
${MOVE} Prior_Diag.nc              ${experiment}/DART
${MOVE} obs_seq.final              ${experiment}/DART
${MOVE} dart_log.out               ${experiment}/DART

${COPY} namelistin                 ${experiment}
${MOVE} namelist                   ${experiment}
${MOVE} casemodel                  ${experiment}

# Good style dictates that you save the scripts so you can see what worked.

${COPY} input.nml                  ${experiment}
${COPY} *.csh                      ${experiment}
${COPY} $myname                    ${experiment}

# CAM leaves a bunch of remnants in your $HOME directory.
# I have not figured out how to use them ... so I clean up.

${REMOVE} ~/lnd.*.rpointer topog_file.nc 

ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

