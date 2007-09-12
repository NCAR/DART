#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

#-----------------------------------------------------------------------------
# job.csh ... Script to run whole assimilation experiment; multiple obs_seq.out files. 
# Resulting series of jobs can take days to run, depending on the numbers of: 
#  > observation sequence files (each is a separate job and may be queued), 
#  > model state variables ( sum of # of fields x # of grid points) 
#  > ensemble members (multiplies the state vector size). 
#  > observations (tens to hundreds of thousands / assim time are common),
#  > processors requested (most efficient machine use is # processors = # ensemble members). 
#
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
#        LSF:   bsub < job_mpi.csh
#        PBS:   qsub job_mpi.csh
#       NONE:   job_mpi.csh
#
# The script moves the necessary files to the current directory and then
# starts 'filter' as a parallel job on all nodes; each of these tasks will 
# call some a separate model_advance.csh when necessary.
#
# This central directory is where the scripts reside and where script and 
# program I/O are expected to happen.
#-----------------------------------------------------------------------------
# 
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system 
# LSF is used on the IBM   Linux cluster 'lightning'
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluevista'
# The queues on lightning and bluevista are supposed to be similar.
#
#
# an explanation of the most common directives follows:
# -J Job name
# -o STDOUT filename
# -e STDERR filename
# -P      account
#    ---------------------------------------------------------------
# coral;
# -q queue        long share standby economy regular premium special
#    min # procs     2     1       1       6      12      24       -
#    max time     5 dy  6 hr   12 hr   12 hr   12 hr   12 hr    5 dy
# blueice
# -q queue        debug share standby economy regular premium special  hold
#    MAX # NODES      2     2      25      25      25      25     100    25   (16 procs/node)
#    max time      30 m  3 hr    2 hr    6 hr    6 hr    6 hr  unlim   6 hr
#    factor         1.0   1.0     0.1     0.5     1.0     1.5     1.0  0.33
#    ---------------------------------------------------------------
# -n number of processors  (it really only needs one; the scripts it creates will 
#                           use more, as specified below)
# -m all possible computational nodes; select your subset here
# CHANGE TO bl... for blueice
# set nodelist = "cr0128en cr0129en cr0130en cr0131en cr0132en cr0133en cr0134en cr0135en cr0136en cr0137en cr0138en cr0140en cr0141en cr0202en cr0201en "
# excluded; cr0139en 
# echo $nodelist

# -W hr:mn   max wallclock time (required on some systems)
# -b [[mm:]dd:]hh:mm    allow job to run only after this time.
##=============================================================================
#BSUB -J job_mpi
#BSUB -o job_mpi.%J.log
#BSUB -q share
#BSUB -W 0:30
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system 
## PBS is used on the CGD   Linux cluster 'bangkok'
## PBS is used on the CGD   Linux cluster 'calgary'
##
##
## an explanation of the most common directives follows:
## -N     Job name
## -r n   Declare job non-rerunable
## -e <arg>  filename for standard error 
## -o <arg>  filename for standard out 
## -q <arg>   Queue name (choose the queue that will allow a single processor job
##                        to execute quickly; share, debug, ...)
## -l nodes=xx:ppn=2   requests BOTH processors on the node. On both bangkok 
##                     and calgary, there is no way to 'share' the processors 
##                     on the node with another job, so you might as well use 
##                     them both.  (ppn == Processors Per Node)
##=============================================================================
#PBS -N job_mpi
#PBS -r n
#PBS -e job_mpi.err
#PBS -o job_mpi.log
#PBS -q debug
##PBS -l nodes=1:ppn=1
##PBS -l ncpus=32

#PBS -S /bin/csh
##PBS -l mem=1GB
##PBS -W group_list=[group name]

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder 
# of the script.

if ($?LS_SUBCWD) then

   # LSF has a list of processors already in a variable (LSB_HOSTS)

   set CENTRALDIR = $LS_SUBCWD
   set JOBNAME = $LSB_JOBNAME
# for multi-thread   
   set run_command = 'mpirun.lsf '
   which $run_command
   if ($status != 0) then
      exit "no mpirun.lsf found"
   endif
# for single-thread (?)
#    set run_command = ' '
# NOT USED on blueice
   alias submit ' bsub < \!* '
   
else if ($?PBS_O_WORKDIR) then

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   set run_command = 'mpirun '
   which $run_command
   if ($status == 0) then
      set run_command = "$run_command -np "
   else
      exit "no mpirun found, exiting"
   endif
   # The number of processors will be added after it's defined below
   set submit = ' qsub '

else if ($?OCOTILLO_NODEFILE) then

   # ocotillo is a 'special case'. It is the only cluster I know of with
   # no queueing system.  You must generate a list of processors in a 
   # file whose name is in $OCOTILLO_NODEFILE.  For example ... 
   # setenv OCOTILLO_NODEFILE  my_favorite_processors
   # echo "node1"  > $OCOTILLO_NODEFILE
   # echo "node5" >> $OCOTILLO_NODEFILE
   # echo "node7" >> $OCOTILLO_NODEFILE
   # echo "node3" >> $OCOTILLO_NODEFILE

   set CENTRALDIR = `pwd`
   set JOBNAME = job_mpi
   # i think this is what we want, but csh will not let you do multiline
   # executions; this argues for using ksh (line 2 below)...  (and maybe
   # it needs a cd as well?)
   #alias submit 'foreach i ($OCOTILLO_NODEFILE) ; ssh $i csh \!* ; end'
   #alias submit='for i in $OCOTILLO_NODEFILE ; do ssh $i (cd $CENTRALDIR; csh $*) ; done'
   set run_command = ' '
   set submit = 'csh '
   
else

   # interactive
   # YOU need to know if you are using the PBS or LSF queuing
   # system ... and set 'submit' accordingly.

   set CENTRALDIR = `pwd`
   set JOBNAME = job_mpi
   set submit = 'csh '
   set run_command = ' '
   
endif

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
      setenv   LINK 'ln -fs'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -vp'
      setenv   MOVE 'mv -fv'
      setenv   LINK 'ln -s'
      breaksw
endsw

echo " "
echo "Running $JOBNAME on host "`hostname`
echo "Initialized at "`date`
echo "CENTRALDIR is " $CENTRALDIR

#========================================================================================
# User set run parameters to change

# Directory where output will be kept (relative to '.')
set resol = T85

set parallel_cam = true

set exp = ${resol}_12


# The "day" of the first obs_seq.out file of the whole experiment
# (if it's not 1, then script changes may be required, esp defining OBS_SEQ)
# Branching a run; easiest to either copy all the desired ICs into the 
# new experiment "previous" output directory and set obs_seq_first = 1,
# OR point to the ICs where they are (DART_ics_1 etc) set obs_seq_first = 
# the current obs_seq to start with, and copy input_n.nml(continuation) 
# into input_#.nml
set obs_seq_first = 1

# number of obs_seq *files* / day  (not timeslots)
set obs_seq_freq = 2

# 'day'/obs_seq.out numbers to assimilate during this job
# First and last obs seq files. 
set obs_seq_1 = 1
set obs_seq_n = 2

# The month of the obs_seq.out files for this run, and 
# the month of the first obs_seq.out file for this experiment.
set mo = 1
set mo_first = 1

# set obs_seq_root = /ptmp/dart/raeder/Obs_sets/Allx12_I/obs_seq2003
# ascii to avoid byteswapping;  but beware of machines filling in extra/missing
# digits, causing locations with latitude > 90 degrees.
set obs_seq_root = /ptmp/dart/Obs_sets/Allx12_I_ascii/obs_seq2003

# DART source code directory trunk, and CAM interface location.
set DARTDIR = /blhome/${user}/J/DART
set DARTCAMDIR =          ${DARTDIR}/models/cam

# The maximum number of processors that will be used by 
# the $exp_#.script jobs spawned by this script.  
# (FV core jobs may use less, depending on the domain decomposition)
# ptile is the number of processors/node on this machine.  
# It has no bearing on whether CAM is MPI or not, as long as filter is MPI.
set max_num_procs = 32
set ptile = 16

# accounting code used for batch jobs (if no accounting needed, you may need
# to remove the -P lines in the script generation sections below.
set proj_num = 12345678

# The queue to which the $exp_#.script scripts will be submitted,
# and the requested time in that queue.
set queue = standby
set wall_clock = 6:00
   
# ICs for obs_seq_first only.  After that the ICs will come from the previous iteration.
# This copies just the initial conditions for the correct number of ensemble members.
# num_lons, num_lats and num_levs are needed for the FV CAM domain decomposition algorithm,
# in order to use larger numbers of processors.
set num_levs  = 26
if ($resol == T21) then
   set DART_ics_1  = /ptmp/dart/CAM_init/T21/03-01-01/DART_MPI
   set CAM_ics_1   = /ptmp/dart/CAM_init/T21/03-01-01/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/T21/03-01-01/CLM/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAMsrc      = /ptmp/dart/CAM/CAMsrc/Cam3/cam3.1/models/atm/cam/bld/T21_3.1-O3
   set num_lons  = 64
   set num_lats  = 32
else if ($resol == T42) then
   # T42
   set DART_ics_1  = /ptmp/dart/CAM_init/T42/03-01-01/DART_MPI
   set CAM_ics_1   = /ptmp/dart/CAM_init/T42/03-01-01/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/T42/03-01-01/CLM/clminput_
   # The CAMsrc directory is MORE than just the location of the executable.
   # There are more support widgets expected in the directory tree.
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAMsrc      = /ptmp/dart/CAM/CAMsrc/Cam3/cam3.1/models/atm/cam/bld/T42_3.1-O3
   set num_lons  = 128
   set num_lats  = 64
else if ($resol == T85) then
   # T85
   set DART_ics_1  = /ptmp/dart/CAM_init/T85/03-01-01/DART_MPI
   set CAM_ics_1   = /ptmp/dart/CAM_init/T85/03-01-01/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/T85/03-01-01/CLM/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAMsrc      = /ptmp/dart/CAM/CAMsrc/Cam3/cam3.1/models/atm/cam/bld/T85_3.1-O3
   set num_lons  = 256
   set num_lats  = 128
else if ($resol == FV4x5) then
   set DART_ics_1  = /ptmp/dart/CAM_init/FV4x5/03-01-01/DART_MPI
   set CAM_ics_1   = /ptmp/dart/CAM_init/FV4x5/03-01-01/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/FV4x5/03-01-01/CLM/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAMsrc      = /ptmp/dart/CAM/CAMsrc/Cam3/cam3.5/models/atm/cam/bld/FV4x5-O2
   set num_lons  = 72
   set num_lats  = 46
else if ($resol == FV2x2.5) then
   set DART_ics_1  = /ptmp/dart/CAM_init/FV2x2.5/03-01-01/DART_MPI
   set CAM_ics_1   = /ptmp/dart/CAM_init/FV2x2.5/03-01-01/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/FV2x2.5/03-01-01/CLM/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAMsrc      = /ptmp/dart/CAM/CAMsrc/Cam3/cam3.5/models/atm/cam/bld/FV2x2.5-O2
   set num_lons  = 144
   set num_lats  = 91
# If another FV resolution is added, then another qualifier is needed in the
# domain decomposition section below.
endif 

if (${parallel_cam} == true) then
   set CAMsrc = ${CAMsrc}-mpi                                                 
endif                                                                         

# e-mail will be sent to this address when obs_seq finishes
set runner = raeder@ucar.edu

# Each obs_seq file gets its own input.nml ... 
# the following variable helps set this up. 

set input = input_

# CHANGE choice of whether forecasts on caminput_##.nc files should be replaced with
#        analyses from filter_ic_new.#### at the end of each obs_seq.  See 'stuff' below

# Choose which sets of restart files will be backed to up a permanent storage place.
# Scripts auto_re2ms*.csh and auto_diag2ms_LSF.csh may need to be modified/replaced on your system.
# Restarts are backed up when
#      if (previous_obs_seq % $save_freq == $mod_save) save to mass store
# where in the obs_seq loop below j = current obs_seq file - 1
set save_freq = 4
set mod_save = 1

# END of run parameters to change
#==========================================================================================

# This is the CENTRAL directory for whole filter job
#    jobs submitted to batch queues from here
#    I/O between filter and advance_model and assim_region goes through here.
#    filter output is put here, then moved to final storage.
cd ${CENTRALDIR}
set myname = $0                              # this is the name of this script
set MASTERLOG = ${CENTRALDIR}/run_job.log    # Set Variable for a 'master' logfile 
${REMOVE} ${MASTERLOG}                       # clean up old links

#----------------------------------------------------------
# try to discover the ensemble size from the input.nml
# this is some gory shell programming ... all to do 'something simple'
grep ens_size input_${obs_seq_first}.nml >! ensstring.$$
set  STRING = "1,$ s#,##g"
set ensstring = `sed -e "$STRING" ensstring.$$`
set num_ens = $ensstring[3]

${REMOVE} ensstring.$$

echo "There are ${num_ens} ensemble members."
# blueice requires the file to exist in order to append to it
# echo "There are ${num_ens} ensemble members."  >> $MASTERLOG
echo "There are ${num_ens} ensemble members."  > $MASTERLOG

#----------------------------------------------------------
# Figure out CAMs domain decomposition and usable number of processors, 
# if it's an FV core.  This information is passed to run-pc.csh via casemodel.
# run-pc.csh uses it in the creation of the CAM namelists.
# User can ignore this 'if' block, unless new resolution is being added
if ($resol == 'FV4x5' || $resol == 'FV2x2.5' ) then
   @ lat_blocks = $num_lats / 3
   if ($lat_blocks >= $max_num_procs) then
      # 1D (lat only) decomposition will work
      set keep_lat_blocks = 0
      set keep_lev_blocks = 0
      set num_procs = $max_num_procs
      echo "Will use $num_procs procs.  Domain decomposed by latitude only "
   else
      # 2D decomposition needed to use available processes
      # lev_blocks = 1 was handled in previous if block
      # This minimizes lev_blocks and maximizes lat_blocks.  Optimal?
      set lev_blocks = 2
      set max_cam_procs = 0
      @ max_lev_blocks = $num_levs / 3
      set done = no
      while ($lev_blocks <= $max_lev_blocks && $done == 'no')
         @ lat_blocks = $num_lats / 3
         @ cam_procs = $lat_blocks * $lev_blocks 
         while ($cam_procs > $max_num_procs && $lat_blocks > 1)
            @ lat_blocks--
            @ cam_procs = $lat_blocks * $lev_blocks
         end
         if ($cam_procs == $max_num_procs) then
            # Good enough; found a combo that uses all available processes
            set done = yes
            set keep_lat_blocks = $lat_blocks
            set keep_lev_blocks = $lev_blocks
         else if ($cam_procs > $max_cam_procs) then
            # test vs previous best number max_CAM_procs
            set max_cam_procs = $cam_procs
            set keep_lat_blocks = $lat_blocks
            set keep_lev_blocks = $lev_blocks
            # not necessarily done
         else
            # Loser; continue looking through possible lev_blocks
         endif
         echo "lev_blocks, lat_blocks, cam_procs = $lev_blocks $lat_blocks $cam_procs "

         @ lev_blocks++
      end

      @ num_procs = $keep_lat_blocks * $keep_lev_blocks
      # CAM likes the following
      set keep_lon_blocks = $keep_lev_blocks 
      echo "Will use $num_procs procs decomposed into "
      echo "    $keep_lat_blocks lat blocks and "
      echo "    $keep_lev_blocks lev blocks for $resol CAM"
   endif
else
   set num_procs = $max_num_procs
   echo "Will use $num_procs procs.  Domain decomposed by latitude only "
endif

# Add information to the run_command, based on User set parameters
if ($?PBS_O_WORKDIR) then
   set run_command = "$run_command $num_procs "
endif

#---------------------------------------------------------
# Get executable programs and scripts from DART-CAM directories

if (! -x filter) then
   ${COPY} ${DARTCAMDIR}/work/filter                     .
   ${COPY} ${DARTCAMDIR}/work/trans_pv_sv                .
   ${COPY} ${DARTCAMDIR}/work/trans_sv_pv                .
   ${COPY} ${DARTCAMDIR}/work/trans_time                 .
   ${COPY} ${DARTCAMDIR}/work/wakeup_filter              .
endif
if (! -e advance_model.csh) then
   ${COPY} ${DARTCAMDIR}/shell_scripts/advance_model.csh     .
   ${COPY} ${DARTCAMDIR}/shell_scripts/run-pc.csh            .
   ${COPY} ${DARTCAMDIR}/shell_scripts/auto_re2ms*.csh       . 
   ${COPY} ${DARTCAMDIR}/shell_scripts/diags.csh             . 
   ${COPY} ${DARTCAMDIR}/shell_scripts/auto_diag2ms_LSF.csh  . 
endif

set days_in_mo = (31 28 31 30 31 30 31 31 30 31 30 31)
# leap years (but year not defined here);   
#    if (($year % 4) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1
# leap year every 4 years except for century marks, but include centuries divisible by 400
#    So, all modern years divisible by 4 are leap years.

#----------------------------------------------------------
echo "exp num_ens obs_seq_1 obs_seq_n obs_seq_first"
echo "$exp $num_ens $obs_seq_1 $obs_seq_n $obs_seq_first"
echo "DART_ics_1 is $DART_ics_1"

echo "exp num_ens obs_seq_1 obs_seq_n obs_seq_first"       >> $MASTERLOG
echo "$exp $num_ens $obs_seq_1 $obs_seq_n $obs_seq_first"  >> $MASTERLOG
echo "DART_ics_1 is $DART_ics_1"                           >> $MASTERLOG

# clean up old CAM inputs that may be laying around

if (-e caminput_1.nc) then
   ${REMOVE} clminput_[1-9]*.nc 
   ${REMOVE} caminput_[1-9]*.nc 
endif

# Remove any possibly stale CAM surface files
rm topog_file.nc

# Ensure the experiment directory exists

if (-d ${exp}) then
   echo "${exp} already exists" >> $MASTERLOG
   echo "${exp} already exists"
else
   echo "Making run-time directory $exp ..." >> $MASTERLOG
   echo "Making run-time directory $exp ..."
   mkdir -p ${exp}
endif

# Subdirectory name root, where output from each obs_seq iteration will be kept.
# obs_diag looks for obs_seq.final files in directories of the form xx_##[#].
# where xx_  = output_root  signifies the month OF THE OBS_SEQ_FIRST.
# and ##[#]# signifies the 2+ digit obs_seq number within this experiment.
set output_root = ${mo_first}_
if (${mo_first} < 10) set output_root = 0${output_root}

#============================================================================================
# Have an overall outer loop over obs_seq.out files

set i = $obs_seq_1
while($i <= $obs_seq_n) ;# start i/obs_seq loop
   echo ' '
   echo ' ' >> $MASTERLOG
   echo '----------------------------------------------------'
   echo '----------------------------------------------------' >> $MASTERLOG
   echo "starting observation sequence file $i at "`date`
   echo "starting observation sequence file $i at "`date` >> $MASTERLOG

#===================================================================================
   # Each iteration of this loop will write out a batch job script named $exp_#.script
   setenv job_i ${exp}_${i}.script

   if ($?LS_SUBCWD) then
      echo "#\!/bin/csh"                   >!  ${job_i}
      echo "##==================================================================" >> ${job_i}
      echo "#BSUB -J ${exp}_${i}"              >> ${job_i}
      echo "#BSUB -o ${exp}_${i}.%J.log"       >> ${job_i}
      echo "#BSUB -P ${proj_num}"              >> ${job_i}
      echo "#BSUB -N -u ${runner} "            >> ${job_i}
      echo "#BSUB -W ${wall_clock}"            >> ${job_i}
         set n_procs = $num_procs
      echo "#BSUB -q ${queue}"                                 >> ${job_i}
      echo "#BSUB -n ${n_procs}"                               >> ${job_i}
# exclusive use of the nodes; still allows > process/node
      echo "#BSUB -x "                                         >> ${job_i}
#      echo '#BSUB -m "'$nodelist '"'           >> ${job_i}
#      if ($?nodelist) echo '#BSUB -m "'$nodelist '"'           >> ${job_i}
      echo '#BSUB -R "span[ptile='$ptile']"'                   >> ${job_i}
# all possible computational nodes; select your subset here

# add this to your other bsub lines

      if ($i > $obs_seq_1) then
         @ previousjobnumber = $i - 1
         set previousjobname = ${exp}_${previousjobnumber}
         echo "#BSUB -w done($previousjobname)" >> ${job_i}
      endif
      echo "##==================================================================" >> ${job_i}

   else if ($?PBS_O_WORKDIR) then
   
      echo "#\!/bin/csh"                                                          >! ${job_i}
      echo "##==================================================================" >> ${job_i}
      echo "#PBS -N ${exp}_${i}"                                                  >> ${job_i}
      echo "#PBS -e ${exp}_${i}.err"                                              >> ${job_i}
      echo "#PBS -o ${exp}_${i}.log"                                              >> ${job_i}
      echo "#PBS -l walltime=${wall_clock}"                                       >> ${job_i}
      echo "#PBS -S /bin/csh"                                                     >> ${job_i}
      echo "#PBS -r n"                                                            >> ${job_i}
      echo "#PBS -m e"                                                            >> ${job_i}
      echo "#PBS -M ${runner}"                                                    >> ${job_i}
      echo "#PBS -q ${queue}"                                                     >> ${job_i}
      echo "#PBS -l ncpus=${num_procs}"                                           >> ${job_i}
#      echo "#PBS -l nodes=${num_nodes}:ppn=${ptile}"                              >> ${job_i}
#      echo "#PBS -l mem=3GB"                                                      >> ${job_i}
      echo "#PBS -W group_list=[your group name]"                                 >> ${job_i}
      if ($i > $obs_seq_1) then
         echo "#PBS -W depend=afterok:$previousjobname"                           >> ${job_i}
      endif
      echo "##==================================================================" >> ${job_i}

   else if ($?OCOTILLO_NODEFILE) then
      echo "#\!/bin/csh"                            >!  ${job_i}
      echo "#BSUB -J ${exp}_${i}"                   >> ${job_i}
      echo "#BSUB -x"                               >> ${job_i}
      echo "#BSUB -n 1"                             >> ${job_i}
      echo '#BSUB -R "span[ptile=1]"'               >> ${job_i}
      echo "#BSUB -W 06:00"                         >> ${job_i}

   else
      echo "#\!/bin/csh"                   >!  ${job_i}
      echo "##==================================================================" >> ${job_i}
   endif

   echo "set myname = "'$0'"     # this is the name of this script"            >> ${job_i}
   echo "set CENTRALDIR =  ${CENTRALDIR} "                                     >> ${job_i}
   echo "cd ${CENTRALDIR}"                                                     >> ${job_i}
   echo "set MASTERLOG = ${CENTRALDIR}/run_job.log"                            >> ${job_i}
   echo 'set start_time = `date  +%s`'                                         >> ${job_i}
   echo ' echo "host is " `hostname` '                                         >> ${job_i}

#===================================================================================

   # Construct directory name of location of restart files
   @ j = $i - 1
   set out_prev = ${output_root}
   if ($j < 10) set out_prev = ${output_root}0
   set out_prev = ${out_prev}$j


   #-----------------------
   # Get filter input files
   # The first one is different than all the rest ...

   echo " " >> ${job_i}
   if ($i == $obs_seq_first) then
      echo "${COPY} ${input}${obs_seq_first}.nml input.nml "         >> ${job_i}
   else
      if (  -e    ${input}n.nml) then
         echo "${COPY}  ${input}n.nml    input.nml "                 >> ${job_i}
      else if (-e ${input}${i}.nml) then
         echo "${COPY}  ${input}${i}.nml input.nml "                 >> ${job_i}
      else
         echo "input_next is MISSING" >> $MASTERLOG  
         echo "input_next is MISSING" >>                
         exit 19 
      endif
   endif

   
   #----------------------------------------------------------------------
   # Get obs_seq file for this assimilation based on date, derived from
   # the obs_seq NUMBER of this iteration ("i").
   # At this time this script assumes that assimilation experiments start
   # at the beginning of the chosen month.
   # As such, we need to ensure any existing obs_seq.out is GONE, which is
   # a little scary.
   #
   # UGLY; need more general calendar capability in here.
   #----------------------------------------------------------------------

   echo " " >> ${job_i}
   echo "${REMOVE} obs_seq.out " >> ${job_i}

   set seq = $i
   if ($seq == 0) then
      set OBS_SEQ = ${CENTRALDIR}/obs_seq.out
   else
      @ month = $mo - 1
      while ($month >= $mo_first)
          @ seq = $seq - $days_in_mo[$month] * $obs_seq_freq
          @ month = $month - 1
      end
      @ month = $month + 1
      if ($month < 10) set month = 0$month

      @ day = (($seq - 1) / $obs_seq_freq) + 1
      if ($day < 10) set day = 0$day

# Kluge to use hour 0 obs_seq file for Test1_n
#      @ hour = (($seq % $obs_seq_freq) - 1) * 24 / $obs_seq_freq 
#      if ($hour < 10) then
# orig      
      if ($obs_seq_freq == 1) then
         set hour = ' '
      else
         @ hour = ($seq % $obs_seq_freq) * 24 / $obs_seq_freq 
         if ($hour == 0) then
            set hour = 24
         else if ($hour < 10) then
            set hour = 0$hour
         endif
      endif

      set OBS_SEQ = ${obs_seq_root}${month}${day}${hour}

      echo "obs_seq, day, hour = $i $day $hour " >> $MASTERLOG
   endif

   if (  -e ${OBS_SEQ} && ! -z ${OBS_SEQ}) then    
      echo "${LINK} $OBS_SEQ  obs_seq.out " >> ${job_i}
   else
      echo "ERROR - no obs_seq.out for $i - looking for:" 
      echo $OBS_SEQ 
      exit 123
   endif

   echo "job- obs_seq $i used is $OBS_SEQ" >> $MASTERLOG
   echo "job- obs_seq $i used is $OBS_SEQ"

   #----------------------------------------------------------------------
   # Get initial conditions for DART, model from 'permanent' storage or
   # from the result of a previous experiment.
   #
   # from_root defines location of DART initial condition(s)
   # cam_init, clm_init are where CAM gets it's initial files (in advance_model.csh)
   #----------------------------------------------------------------------

   if ($i == $obs_seq_first) then 
      # Get 'initial' initial files
      set from_root = ${DART_ics_1} 
      set  cam_init = ${CAM_ics_1}
      set  clm_init = ${CLM_ics_1}
   else
      # Get initial files from result of previous experiment.
      set from_root = `pwd`/$exp/${out_prev}/DART
      if (! -e times) then
         # There was no advance before the first assimilation,
         # so there is not a set of c[al]minput.nc files associated with the filter_ic.####s
         # Use the old ones, which are still valid at this time.
         set  cam_init = ${CAM_ics_1}
         set  clm_init = ${CLM_ics_1}
      else
         set cam_init =  `pwd`/$exp/${out_prev}/CAM/caminput_
         set clm_init =  `pwd`/$exp/${out_prev}/CLM/clminput_
      endif
   endif

   # transmit info to advance_model.csh
   # The second item echoed must be the directory where CAM is kept
   echo "$exp "                                            >! casemodel.$i
   echo "$CAMsrc "                                         >> casemodel.$i
   echo "$cam_init "                                       >> casemodel.$i
   echo "$clm_init"                                        >> casemodel.$i
   echo "$run_command"                                     >> casemodel.$i
   echo "$num_procs $keep_lev_blocks $keep_lat_blocks "    >> casemodel.$i
   # advance_model wants to see a file 'casemodel' and not keep track of which obs_seq it's for
   echo "$REMOVE casemodel"                >> ${job_i}
   echo "$LINK casemodel.$i casemodel "    >> ${job_i}

   # adaptive inflation ic files may (not) exist
   # Should query input.nml to learn whether to get them?
   echo " "                                                      >> ${job_i}
   echo "${REMOVE} *_inf_ic* "                                   >> ${job_i}
   echo "if (-e ${from_root}/prior_inf_ic) \"                    >> ${job_i}
   echo " ${LINK} ${from_root}/prior_inf_ic prior_inf_ic_old "   >> ${job_i}
   echo "if (-e ${from_root}/post_inf_ic)  \"                    >> ${job_i}
   echo " ${LINK} ${from_root}/post_inf_ic  post_inf_ic_old "    >> ${job_i}

#? MPI too?
   # link to filter_ic file(s), so that filter can copy them to a compute node
   echo " " >> ${job_i}
   echo "if (-e ${from_root}/filter_ic.0001) then "            >> ${job_i}
   echo "   set n = 1"                                         >> ${job_i}
   echo '   while ($n ' "<= ${num_ens})"                       >> ${job_i}
   echo "         set from = ${from_root}/filter_ic*[.0]"'$n'  >> ${job_i}
   echo "         ${REMOVE}         filter_ic_old."'$from:e'   >> ${job_i}
   echo "         ${LINK} "'$from'" filter_ic_old."'$from:e'   >> ${job_i}
   echo "         @ n++"                                       >> ${job_i}
   echo "   end"                                               >> ${job_i}
   echo "else if (-e ${from_root}/filter_ic) then "            >> ${job_i}
   echo "   ${REMOVE} filter_ic_old "                          >> ${job_i}
   echo "   ${LINK} ${from_root}/filter_ic filter_ic_old "     >> ${job_i}
   echo "endif "                                               >> ${job_i}  
   echo ' '
   echo 'echo "job- filter_ic_old is/are" '                    >> ${job_i}
   echo "ls -lt filter_ic_old* "                               >> ${job_i}

   #-----------------------------------------------------------------------------
   # link local CAM input files to generic names in CENTRALDIR.  
   # These just provide grid info to filter, not state info.
   # CAM namelist input file; will be augmented with model advance time and
   # domain decomposition info by run-pc.csh
   #-----------------------------------------------------------------------------
   echo " " >> ${job_i}
   echo "${REMOVE} caminput.nc clminput.nc "                 >> ${job_i}
   echo "${LINK} ${CAMsrc}/caminput.nc caminput.nc"          >> ${job_i}
   echo "${LINK} ${CAMsrc}/clminput.nc clminput.nc"          >> ${job_i}
   
   #-----------------------------------------------------------------------------
   # get name of file containing PHIS from the CAM namelist.  This will be used by
   # static_init_model to read in the PHIS field, which is used for height obs.
   #-----------------------------------------------------------------------------
   ${REMOVE} namelistin
   ${LINK} ${CAMsrc}/namelistin namelistin
   sleep 1 
   if (! -e namelistin ) then
      echo "ERROR ... need a namelistin file." >> $MASTERLOG
      echo "ERROR ... need a namelistin file."
      exit 89
   endif
   echo "if (! -e ${exp}/namelistin) ${COPY} namelistin ${exp}/namelistin  "       >> ${job_i}

   # Observations on heights require PHIS, which are on a CAM surface fields file.
   # Only need to do this for obs_seq_1. 
   if (! -e topog_file.nc) then
      grep bnd_topo namelistin >! topo_file
      set  STRING = "1,$ s#'##g"
      set ensstring = `sed -e "$STRING" topo_file`
      set topo_name = $ensstring[3]
      if (-e $topo_name) then
         ${COPY} $topo_name topog_file.nc
      else
         echo "ERROR ... need a topog_file; none in namelistin." >> $MASTERLOG
         echo "ERROR ... need a topog_file; none in namelistin."
         exit 99
      endif
      
      chmod 644 topog_file.nc
      ${REMOVE} topo_file
   endif

   #======================================================================
   # Run the filter in async=2 mode.

   # runs filter, which tells the model to model advance and assimilates obs
   echo " "                                                  >> ${job_i}
   if (${parallel_cam} == true) then
      echo "${COPY} ${CENTRALDIR}/wakeup_filter  .                       " >> ${job_i}
      echo " "                                                             >> ${job_i}
      echo "rm -f  filter_to_model.lock model_to_filter.lock "             >> ${job_i}
      echo "mkfifo filter_to_model.lock model_to_filter.lock "             >> ${job_i}
      echo 'setenv OHOME ${HOME} '                                         >> ${job_i}
      echo 'setenv NHOME ${HOME}/test '                                    >> ${job_i}
      echo 'if (! -d ${NHOME}) mkdir ${NHOME} '                            >> ${job_i}
      echo 'setenv HOME ${NHOME} '                                         >> ${job_i}
      echo "${run_command} ./filter &"                                     >> ${job_i}
      echo 'setenv HOME ${OHOME} '                                         >> ${job_i}
      echo " "                                                             >> ${job_i}
      echo 'while ( -e filter_to_model.lock )          '                   >> ${job_i}
      # read from the fifo file.  this is *not* a busy wait; it puts the
      # job to sleep in the kernel waiting for input.
      echo " "                                                >> ${job_i}
      echo '  set todo = `( echo $< ) < filter_to_model.lock` ' >> ${job_i}
      echo '  echo todo received, value = ${todo}           ' >> ${job_i}
      echo " "                                                >> ${job_i}
      echo '  if ( "${todo}" == "finished" ) then           ' >> ${job_i}     
      echo '    echo finished command received from filter. ' >> ${job_i}
      echo '    echo main script: filter done.              ' >> ${job_i}
     # add this wait to be sure filter task has exited
     # before starting to clean up the files.
     echo '    wait                                            ' >> ${job_i}
     echo '    echo filter finished, removing pipes.           ' >> ${job_i}
     echo "    rm -f filter_to_model.lock model_to_filter.lock " >> ${job_i}

      echo '    break                                       ' >> ${job_i}
      echo " "                                                >> ${job_i}
      echo '  else if ( "${todo}" == "advance" ) then       ' >> ${job_i}
      echo '    echo advance command received from filter.  ' >> ${job_i}
      echo '    echo calling advance_model.csh now:              ' >> ${job_i}
     echo "    ./advance_model.csh 0 $num_ens filter_control00000  ${parallel_cam}" >> ${job_i}
     # do not execute anything here until you have saved
     # the exit status from the advance model script.
     echo '    set advance_status = $status                ' >> ${job_i}
     echo '    echo saved advance_model.csh exit status    ' >> ${job_i}
     echo " "                                                >> ${job_i}
     echo '    echo restarting filter.  this version of wakeup_filter ' >> ${job_i}
     echo '    echo includes restarting the main filter program.      ' >> ${job_i}
     echo "    ${run_command} ./wakeup_filter              " >> ${job_i}
     echo " "                                                >> ${job_i}
     echo '    if ($advance_status != 0) then              ' >> ${job_i}
     echo '       echo "Model advance failed"              ' >> ${job_i}
     echo '       rm -f filter_lock*                       ' >> ${job_i}
     echo '       break                                    ' >> ${job_i}
     echo '    endif                                       ' >> ${job_i}
     echo " "                                                >> ${job_i}
      echo '  else                                          ' >> ${job_i}
      echo " "                                                >> ${job_i}
      echo '    echo main script: unexpected value received.' >> ${job_i}
      echo '    break                                       ' >> ${job_i}
      echo " "                                                >> ${job_i}
      echo '  endif                                         ' >> ${job_i}
      echo " "                                                >> ${job_i}
      echo 'end                                             ' >> ${job_i}     
      echo " "                                                >> ${job_i}
   else
      # Run the filter in async=2 mode.
      # runs filter, which tells the model to model advance and assimilates obs
      echo "${run_command} ./filter "                        >> ${job_i}
   endif

   set KILLCOMMAND = "touch BOMBED; exit"

   #-----------------
   # When filter.f90 finishes an obs_seq.out file, it creates a file called 'go_end_filter' 
   # in the CENTRALDIR.  Under the MPI/async=2 scenario this won't be used.  The successful
   # completion of the $job_# script will tell the (next) $job_[#+1] script to start.


   #-----------------------------------------------
   # Move the output to storage after filter signals completion.
   # At this point, all the restart,diagnostic files are in the CENTRALDIR
   # and need to be moved to the 'experiment permanent' directory.
   # We have had problems with some, but not all, files being moved
   # correctly, so we are adding bulletproofing to check to ensure the filesystem
   # has completed writing the files, etc. Sometimes we get here before
   # all the files have finished being written.
   #-----------------------------------------------
   # This was clean, but did not always work ... sigh ...
   # ${MOVE} clminput_[1-9]*.nc    ${exp}/${output_dir}/CLM
   # ${MOVE} caminput_[1-9]*.nc    ${exp}/${output_dir}/CAM
   #-----------------------------------------------

   echo " " >> ${job_i}
   echo 'echo "Listing contents of CENTRALDIR before archiving at "`date` ' >> ${job_i}
   echo "ls -l "                                                            >> ${job_i}
   echo "Listing contents of CENTRALDIR before archiving at "`date` >> $MASTERLOG
   ls -l >> $MASTERLOG

#-----------------------------------------------------------------------------
# Ensure the (output) experiment directories exist
# All the  CAM-related files will get put in ${exp}/${output_dir}/CAM
# All the  CLM-related files will get put in ${exp}/${output_dir}/CLM
# All the DART-related files will get put in ${exp}/${output_dir}/DART
#-----------------------------------------------------------------------------
   set output_dir = ${output_root}
   if ($i < 10) set output_dir = ${output_root}0
   set output_dir = ${output_dir}$i

   echo " " >> ${job_i}
   echo "mkdir -p ${exp}/${output_dir}/{CAM,CLM,DART} "                               >> ${job_i}

   echo " " >> ${job_i}
   echo "foreach FILE ( Prior_Diag.nc Posterior_Diag.nc obs_seq.final )"              >> ${job_i}
   echo '   if ( -e $FILE && ! -z $FILE) then '                                       >> ${job_i}
   echo "      ${MOVE} "'$FILE'" ${exp}/${output_dir} "                               >> ${job_i}
   echo "      if ( ! "'$status'" == 0 ) then "                                       >> ${job_i}
   echo '         echo "failed moving ${CENTRALDIR}/$FILE" >> $MASTERLOG '            >> ${job_i}
   echo '         echo "failed moving ${CENTRALDIR}/$FILE" '                          >> ${job_i}
   echo "         $KILLCOMMAND "                                                      >> ${job_i}
   echo "      endif           "                                                      >> ${job_i}
   echo "   else               "                                                      >> ${job_i}
   echo '      echo "... THUMP ... ${CENTRALDIR}/$FILE does not exist and should."'   >> ${job_i}
   echo '      echo "directory contents follow" '                                     >> ${job_i}
   echo "      ls -l "                                                                >> ${job_i}
   echo "      $KILLCOMMAND "                                                         >> ${job_i}
   echo "   endif "                                                                   >> ${job_i}
   echo "end "                                                                        >> ${job_i}

   # inflate_diag may or may not exist (when do_obs_inflate=F do_single_ss_inflate=F)
   # so don't die if it's missing.  We'd have to query input.nml to learn if it should exist.

   echo " " >> ${job_i}
   echo "foreach FILE ( prior_inf_diag post_inf_diag ) "                              >> ${job_i}
   echo '   if ( -e ${FILE} && ! -z $FILE) then '                                     >> ${job_i}
   echo "      ${MOVE} "'${FILE}'" ${exp}/${output_dir}  "                            >> ${job_i}
   echo '      if ( ! $status == 0 ) then '                                           >> ${job_i}
   echo '         echo "failed moving ${CENTRALDIR}/${FILE} " >> $MASTERLOG '         >> ${job_i}
   echo '         echo "failed moving ${CENTRALDIR}/${FILE} " '                       >> ${job_i}
   echo "         $KILLCOMMAND "                                                      >> ${job_i}
   echo "      endif "                                                                >> ${job_i}
   echo "   endif "                                                                   >> ${job_i}
   echo "end "                                                                        >> ${job_i}

   # Move the filter restart file(s) to the storage subdirectory
   echo " " >> ${job_i}
   echo 'echo "moving filter_ic_newS to '${exp}/${output_dir}'/DART/filter_icS" '     >> ${job_i}
   echo "if (-e filter_ic_new) then "                                                 >> ${job_i}
   echo "   ${MOVE} filter_ic_new ${exp}/${output_dir}/DART/filter_ic "               >> ${job_i}
   echo "   if (! "'$status'" == 0 ) then "                                           >> ${job_i}
   echo '      echo "failed moving filter_ic_new to '${exp}/${output_dir}'/DART/filter_ic" ' >> ${job_i}
   echo "      $KILLCOMMAND "                                                         >> ${job_i}
   echo "   endif "                                                                   >> ${job_i}
   echo "else if (-e filter_ic_new.0001) then "                                       >> ${job_i}
   echo "   set n = 1 "                                                               >> ${job_i}
   echo "   while("'$n'" <= ${num_ens}) "                                             >> ${job_i}
   echo '        set from = filter_ic_new*[.0]$n'                                     >> ${job_i}
   echo "        set dest = ${exp}/${output_dir}/DART/filter_ic."'$from:e'            >> ${job_i}

   # stuff analyses into CAM initial files using 
   # echo " " >> ${job_i}
   # echo '      echo "$from >\! member "'                                              >> ${job_i}
   # echo '      echo "caminput_${n}.nc >> member "'                                    >> ${job_i}
   # echo "      ./trans_sv_pv "                                                        >> ${job_i}
   # echo '      echo "stuffing analyses from $from into caminput_${n}.nc >> $MASTERLOG ' >> ${job_i}
   # echo '      ls -l caminput_${n}.nc >> $MASTERLOG   '                               >> ${job_i}
   # end stuffing 

   echo " " >> ${job_i}
   echo "        ${MOVE} "'$from $dest'                                               >> ${job_i}
   echo "        if (! "'$status'" == 0 ) then "                                      >> ${job_i}
   echo '           echo "failed moving $from to ${dest} "'                           >> ${job_i}
   echo "           $KILLCOMMAND "                                                    >> ${job_i}
   echo "        endif "                                                              >> ${job_i}
   echo "        @ n++ "                                                              >> ${job_i}
   echo "   end "                                                                     >> ${job_i}
   echo "else "                                                                       >> ${job_i}
   echo '   echo "NO filter_ic_new FOUND" '                                           >> ${job_i}
   echo '   echo "NO filter_ic_new FOUND" '                                           >> ${job_i}
   echo '   echo "NO filter_ic_new FOUND" '                                           >> ${job_i}
   echo '   echo "NO filter_ic_new FOUND" '                                           >> ${job_i}
   echo "   $KILLCOMMAND "                                                            >> ${job_i}
   echo "endif "                                                                      >> ${job_i}

   echo " " >> ${job_i}
   echo "foreach FILE (prior_inf_ic post_inf_ic)  "                                   >> ${job_i}
   echo '   if (-e ${FILE}_new ) then '                                               >> ${job_i}
   echo "      ${MOVE} "'${FILE}_new'" ${exp}/${output_dir}/DART/"'${FILE} '          >> ${job_i}
   echo '      if (! $status == 0 ) then '                                            >> ${job_i}
   echo '         echo "failed moving ${FILE}_new to '${exp}/${output_dir}/DART/'${FILE}s "' >> ${job_i}
   echo "         $KILLCOMMAND "                                                      >> ${job_i}
   echo "      endif "                                                                >> ${job_i}
   echo "   endif "                                                                   >> ${job_i}
   echo "end "                                                                        >> ${job_i}

   echo " " >> ${job_i}
   echo "set n = 1 "                                                                  >> ${job_i}
   echo 'while ($n <= '"${num_ens})    ;# loop over all ensemble members "            >> ${job_i}
   echo '   set CAMINPUT = caminput_${n}.nc  '                                        >> ${job_i}
   echo '   set CLMINPUT = clminput_${n}.nc  '                                        >> ${job_i}

   echo "if (! -e times && ! -e ${exp}/${output_dir}/CAM/caminput_1.nc) then  "       >> ${job_i}
   echo "   # There may have been no advance; "                                       >> ${job_i}
   echo "   # use the CAM_ics_1 as the CAM ics for this time  "                       >> ${job_i}
   echo "   set ens = 1  "                                                            >> ${job_i}
   echo '   while ($ens <= '$num_ens" )  "                                            >> ${job_i}
   echo "      cp ${CAM_ics_1}"'${ens}'".nc ${exp}/${output_dir}/CAM  "               >> ${job_i}
   echo "      cp ${CLM_ics_1}"'${ens}'".nc ${exp}/${output_dir}/CLM  "               >> ${job_i}
   echo "      @ ens++  "                                                             >> ${job_i}
   echo "   end  "                                                                    >> ${job_i}
   echo "endif  "                                                                     >> ${job_i}

   echo " " >> ${job_i}
   echo '   if ( -e $CAMINPUT && ! -z $CAMINPUT) then '                               >> ${job_i}
   echo "      ${MOVE} "'$CAMINPUT'" ${exp}/${output_dir}/CAM  "                      >> ${job_i}
   echo '      if (! $status == 0 ) then '                                            >> ${job_i}
   echo '         echo "failed moving ${CENTRALDIR}/$CAMINPUT " '                     >> ${job_i}
   echo "         $KILLCOMMAND "                                                      >> ${job_i}
   echo "      endif "                                                                >> ${job_i}
   echo "   else "                                                                    >> ${job_i}
   echo '      echo "${CENTRALDIR}/$CAMINPUT does not exist and maybe should." '      >> ${job_i}
# kdr debug; this should be dependent on whether a model advance was needed before
#            the assim.
#   echo "      $KILLCOMMAND "                                                         >> ${job_i}
   echo "   endif "                                                                   >> ${job_i}

   echo " " >> ${job_i}
   echo '   if ( -e $CLMINPUT && ! -z $CLMINPUT) then '                               >> ${job_i}
   echo "      ${MOVE} "'$CLMINPUT'" ${exp}/${output_dir}/CLM "                       >> ${job_i}
   echo '      if (! $status == 0 ) then '                                            >> ${job_i}
   echo '         echo "failed moving ${CENTRALDIR}/$CLMINPUT " '                     >> ${job_i}
   echo "         $KILLCOMMAND "                                                      >> ${job_i}
   echo "      endif "                                                                >> ${job_i}
   echo "   else "                                                                    >> ${job_i}
   echo '      echo "${CENTRALDIR}/$CLMINPUT does not exist and maybe should." '      >> ${job_i}
# kdr debug; this should be dependent on whether a model advance was needed before
#            the assim.
#   echo "      $KILLCOMMAND "                                                         >> ${job_i}
   echo "   endif "                                                                   >> ${job_i}

   echo " " >> ${job_i}
   echo "  @ n++ "                                                                    >> ${job_i}
   echo "end "                                                                        >> ${job_i}


   # test whether it's safe to end this obs_seq_ 
   echo " " >> ${job_i}
   echo "if ((  -e $exp/${output_dir}/DART/filter_ic &&               \"              >> ${job_i}
   echo "     ! -z $exp/${output_dir}/DART/filter_ic   )    ||        \"              >> ${job_i}
   echo "    (  -e $exp/${output_dir}/DART/filter_ic.*${num_ens}  &&  \"              >> ${job_i}
   echo "     ! -z $exp/${output_dir}/DART/filter_ic.*${num_ens}) ) then "            >> ${job_i}
   echo '    echo "it is OK to proceed with next obs_seq at "`date` >> $MASTERLOG '   >> ${job_i}
   echo '    echo "it is OK to proceed with next obs_seq at "`date` '                 >> ${job_i}
   echo 'else '                                                                       >> ${job_i}
   echo '    echo "RETRIEVE filter_ic files from filter temp directory ?" '           >> ${job_i}
   echo '    echo "Then remove temp and cam advance temps" '                          >> ${job_i}
   echo '    echo "RETRIEVE filter_ic files from filter temp directory ?" >> $MASTERLOG ' >> ${job_i}
   echo '    echo "Then remove temp and cam advance temps" >> $MASTERLOG '            >> ${job_i}
   echo '    exit '                                                                   >> ${job_i}
   echo 'endif '                                                                      >> ${job_i}

#? Should eval be outside of echoing?  can the result be passed in then?

   # Compress and archive older output which we won't need immediately
   echo " " >> ${job_i}
   if (($j > 0) && (-e auto_re2ms_LSF.csh)) then
      echo "cd ${exp}/${out_prev} "                                                         >> ${job_i}
      # echo "if ( $j > $obs_seq_first ) then "                                            >> ${job_i}
      if ($j % $save_freq == $mod_save) then
# Not on blueice
#         echo "alias submit ' bsub < \\!* '"                                                >> ${job_i}
#         echo "alias | grep submit "                                                        >> ${job_i}
#         echo "eval submit ../../auto_re2ms_LSF.csh  >>& " ' $MASTERLOG '                   >> ${job_i}
         echo "eval bsub < ../../auto_re2ms_LSF.csh  >>& " ' $MASTERLOG '                   >> ${job_i}
         echo 'echo "Backing up restart '${exp}/${out_prev}' to mass store; >> $MASTERLOG"' >> ${job_i}
         echo 'echo "    in separate batch job"  >> $MASTERLOG  '                           >> ${job_i}
         echo 'echo "Backing up restart '${exp}/${out_prev}' to mass store; "'              >> ${job_i}
         echo 'echo "    in separate batch job"  '                                          >> ${job_i}
      else
         echo "${REMOVE} DART CAM CLM & "                                                   >> ${job_i}   
         echo 'echo "Removing restart from '${exp}/${out_prev} '" >> $MASTERLOG '           >> ${job_i}
         echo 'echo "Removing restart from '${exp}/${out_prev} '"'                          >> ${job_i}
      endif
#     Diagnostics files
      echo "eval bsub < ../../auto_diag2ms_LSF.csh  >>& " ' $MASTERLOG '                     >> ${job_i}
      echo 'echo "Backing up diagnostics '${exp}/${out_prev}' to mass store; >> $MASTERLOG"' >> ${job_i}
      echo 'echo "    in separate batch job"  >> $MASTERLOG  '                               >> ${job_i}
      echo 'echo "Backing up diagnostics '${exp}/${out_prev}' to mass store; "'              >> ${job_i}
      echo 'echo "    in separate batch job"  '                                              >> ${job_i}

      echo "cd ../.. "                                                                      >> ${job_i}
   else
      echo "NO ./auto_re2ms_LSF.csh FOUND, so NO BACKUP OF ${out_prev} " >> ${MASTERLOG}
   endif


# save a representative model advance 
   echo " " >> ${job_i}
   echo "${MOVE} input.nml           ${exp}/${output_dir} "                           >> ${job_i}
   echo "${MOVE} casemodel.$i        ${exp}/${output_dir} "                           >> ${job_i}
# stdout   from filter.f90 
   echo "${MOVE} run_filter.stout    ${exp}/${output_dir} "                           >> ${job_i}
   echo "${MOVE} dart_log.out        ${exp}/${output_dir} "                           >> ${job_i}
   echo "${MOVE} ${job_i}            ${exp}/${output_dir} "                           >> ${job_i}
   echo "${MOVE} cam_out_temp1       ${exp}/${output_dir} "                           >> ${job_i}
   echo "${MOVE} namelist            ${exp}"                                          >> ${job_i}
   echo "${REMOVE} cam_out_temp* *_ud* *_ic[0-9]* *_ic_old* "                         >> ${job_i}
   echo "set nonomatch"                                                               >> ${job_i}
   if ($i != $obs_seq_first) then
      @ j = $i - 1
#     This will fail if >1 file satisfies the existence
      echo "if (-e ${exp}_${j}.*.log) ${MOVE} ${exp}_${j}.*.log  ${exp}/${out_prev} " >> ${job_i}
   endif
   echo 'set end_time = `date  +%s`'                                     >> ${job_i}
   echo '@ length_time = $end_time - $start_time'                        >> ${job_i}
   echo 'echo "duration = $length_time" '                                >> ${job_i}


# Finally, submit the script that we just created.
# Doesn't work on blueice
#   eval submit ${job_i} >! batchsubmit$$
# The beauty of hard-wiring
   bsub < ${job_i} >! batchsubmit$$


if ($?PBS_O_WORKDIR) then
   set previousjobname = `cat batchsubmit$$`
   set FILTERBATCHID = $previousjobname
else if ($?LS_SUBCWD) then
   set STRING = "1,$ s#<##g"
   sed -e "$STRING" batchsubmit$$ >! bill$$
   set STRING = "1,$ s#>##g"
   sed -e "$STRING" bill$$ >! batchsubmit$$
   set STRING = `cat batchsubmit$$`
   set FILTERBATCHID = $STRING[2]
   # set FILTERBATCHID = "none4test"
   ${REMOVE} batchsubmit$$ bill$$
endif

   echo "filter        spawned as job $FILTERBATCHID at "`date`
   echo "filter        spawned as job $FILTERBATCHID at "`date`       >> $MASTERLOG
   # set KILLCOMMAND = "bkill $FILTERBATCHID; touch BOMBED; exit"

   echo "completed iteration $i ($OBS_SEQ) at "`date`
   @ i++
end # end of the huge "i" loop

# Move the stout from this script to the experiment directory (with a meaningful name)
# and make a link to it using the old name, so that each of the $job_i scripts
# can write their stout to the same file.
${MOVE} ${MASTERLOG} ${exp}/run_job_${obs_seq_1}-${obs_seq_n}.log   
${LINK}              ${exp}/run_job_${obs_seq_1}-${obs_seq_n}.log run_job.log

# actual namelist used by CAM for most recent model advance
${REMOVE} ~/lnd.*.rpointer

