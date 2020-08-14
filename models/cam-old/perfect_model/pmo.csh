#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script to run perfect_model_obs to generate synthetic observations
#    from multiple obs_seq.in files. 
#
# The CENTRAL directory needs to have:
#    caminput_0.nc
#    clminput_0.nc
#    input_X.nml   X = the first obs_seq to be processed. (Used in input_X.nml name)
#    perfect_ic
#    cam_phis.nc   (DART as of 1/2008)
# The following will be retrieved from $CAM_src if necessary
#     advance_model.csh
#     perfect_model_obs
#     run-cam.csh
#     trans_pv_sv
#     trans_sv_pv
#     trans_time
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
#        LSF:   bsub < pmo.csh
#        PBS:   qsub pmo.csh
#       NONE:   pmo.csh
#
# The script moves the necessary files to the current directory and then
# starts 'perfect_model_obs' as a single-threaded job; 
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

# -W hr:mn   max wallclock time (required on some systems)
# -b [[mm:]dd:]hh:mm    allow job to run only after this time.
##=============================================================================
#BSUB -J pmo
#BSUB -o pmo.%J.log
#BSUB -P 93300315
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
#PBS -N pmo
#PBS -r n
#PBS -e pmo.err
#PBS -o pmo.log
#PBS -q debug
##PBS -l nodes=1:ppn=1
##PBS -l ncpus=32

#PBS -S /bin/csh
##PBS -l mem=1GB
##PBS -W group_list=[group name]

##=============================================================================
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
#   set run_command = 'mpirun.lsf '
# for single-thread (?)
   set run_command = ' '
   
else if ($?PBS_O_WORKDIR) then

   # PBS has a list of processors in a file whose name is (PBS_NODEFILE)

   set CENTRALDIR = $PBS_O_WORKDIR
   set JOBNAME = $PBS_JOBNAME
   # pmo doesn't run parallel
   #   set run_command = 'mpirun '
   set run_command = ' '

else

   # interactive
   # YOU need to know if you are using the PBS or LSF queuing
   # system ... and set 'submit' accordingly.

   set CENTRALDIR = `pwd`
   set JOBNAME = pmo
   set submit = 'csh '
   set run_command = ' '
   
endif

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
set resol = FV1.9x2.5

set parallel_cam = false

set exp = PMO2


# The "day" of the first obs_seq.in file of the whole experiment
# (if it's not 1, then script changes may be required, esp defining OBS_SEQ)
set obs_seq_first = 1

# Time spacing of obs_seq files
# freq >  0      number of obs_seq *files* / day  (not timeslots)
# freq = 0       set OBS_SEQ = ${CENTRALDIR}/obs_seq.in
# -15 < freq < 0 that number of days from one to the next (3 means 20070101, 20070104, ...)
# freq = -15     look for the first and 15th of each month
# freq < -15     look for the first of each month

set obs_seq_freq = 2

# 'day'/obs_seq.in numbers to assimilate during this job
# First and last obs seq files. 
set obs_seq_1 = 1
set obs_seq_n = 2


# The month of the obs_seq.in files for this run, and 
# the month of the first obs_seq.in file for this experiment.
# This will be a misnomer for the spin-up run that has obs_seq files
# only at the first day of the month; the 'days' will refer to months
# and this mo can be thought of as the year (2001).

set mo = 7
set mo_first = 7

#  set obs_seq_root = /ptmp/dart/raeder/Obs_sets/Allx12_I/obs_seq2003
#  ascii to avoid byteswapping;  but beware of machines filling in extra/missing
#  digits, causing locations with latitude > 90 degrees.
#  set obs_seq_root = ${CENTRALDIR}/obs_seq2007
set obs_seq_root = /ptmp/dart/Obs_sets/Allx12_I_ascii/obs_seq2007

# DART source code directory trunk, and CAM interface location.
set DARTDIR = /blhome/${user}/J-svn/DART
set DARTCAMDIR =          ${DARTDIR}/models/cam

# The maximum number of processors that will be used by 
# the $exp_#.script jobs spawned by this script.  
# (FV core jobs may use less, depending on the domain decomposition)
# ptile is the number of processors/node on this machine.  
# It has no bearing on whether CAM is MPI or not, as long as perfect_model_obs is MPI.
set max_num_procs = 1
set ptile = 1

# accounting code used for batch jobs (if no accounting needed, you may need
# to remove the -P lines in the script generation sections below.
set proj_num = 93300315

# The queue to which the $exp_#.script scripts will be submitted,
# and the requested time in that queue.
set queue = share
set wall_clock = 1:00
   
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
   set CAM_src      = /ptmp/dart/CAM/CAM_src/Cam3/cam3.1/models/atm/cam/bld/T21_3.1-O3
   set CAM_phis  = $CAM_src/cam_phis.nc
   set num_lons  = 64
   set num_lats  = 32
else if ($resol == T42) then
# T42
   set DART_ics_1  = /ptmp/dart/CAM_init/T42/03-01-01/DART_MPI
   set CAM_ics_1   = /ptmp/dart/CAM_init/T42/03-01-01/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/T42/03-01-01/CLM/clminput_
   # The CAM_src directory is MORE than just the location of the executable.
   # There are more support widgets expected in the directory tree.
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAM_src      = /ptmp/dart/CAM/CAM_src/Cam3/cam3.1/models/atm/cam/bld/T42_3.1-O3
   set CAM_phis  = $CAM_src/cam_phis.nc
   set num_lons  = 128
   set num_lats  = 64
else if ($resol == T85) then
# T85
   set DART_ics_1  = /ptmp/dart/CAM_init/T85_cam3.5/Jul_1/DART
   set CAM_ics_1   = /ptmp/dart/CAM_init/T85_cam3.5/Jul_1/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/T85_cam3.5/Jul_1/CLM/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAM_src      = /blhome/raeder/Cam3/cam3.5/models/atm/cam/bld/T85-O3
   # set CAM_src      = /ptmp/dart/CAM/CAM_src/Cam3/cam3.1/models/atm/cam/bld/T85_3.1-O3
   set CAM_phis  = $CAM_src/cam_phis.nc
   set num_lons  = 256
   set num_lats  = 128
else if ($resol == FV4x5) then
   set DART_ics_1  = /ptmp/dart/CAM_init/FV4x5/03-01-01/DART_MPI
   set CAM_ics_1   = /ptmp/dart/CAM_init/FV4x5/03-01-01/CAM/caminput_
   set CLM_ics_1   = /ptmp/dart/CAM_init/FV4x5/03-01-01/CLM/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAM_src      = /ptmp/dart/CAM/CAM_src/Cam3/cam3.5/models/atm/cam/bld/FV4x5-O2
   set CAM_phis  = $CAM_src/cam_phis.nc
   set num_lons  = 72
   set num_lats  = 46
else if ($resol == FV2x2.5) then
# the first ICs should be in c[al]minput_0.nc and perfect_ic
   set DART_ics_1  = ${CENTRALDIR}
   set CAM_ics_1   = ${CENTRALDIR}/caminput_
   set CLM_ics_1   = ${CENTRALDIR}/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAM_src      = /ptmp/dart/CAM/CAM_src/Cam3/cam3.5/models/atm/cam/bld/FV2x2.5-O2
   set CAM_phis  = $CAM_src/cam_phis.nc
   set num_lons  = 144
   set num_lats  = 91
else if ($resol == FV1.9x2.5) then
   set DART_ics_1  = ${CENTRALDIR}
   set CAM_ics_1   = ${CENTRALDIR}/caminput_
   set CLM_ics_1   = ${CENTRALDIR}/clminput_
   # -mpi will be attached to this name if parallel_cam = true; don't add it here
   set CAM_src      = /blhome/raeder/Cam3/cam3.5/models/atm/cam/bld/FV1.9x2.5-O3
   set CAM_phis  = $CAM_src/cam_phis.nc
   set num_lons  = 144
   set num_lats  = 96
# If another FV resolution is added, then another qualifier is needed in the
# domain decomposition section below.
endif 

# e-mail will be sent to this address when obs_seq finishes
set runner = raeder@ucar.edu

# Each obs_seq file gets its own input.nml ... 
# the following variable helps set this up. 

set input = input_

# Choose which sets of restart files will be backed to up a permanent storage place.
# Scripts auto_re2ms*.csh and auto_diag2ms_LSF.csh may need to be modified/replaced on your system.
# Restarts are backed up when
#      if (previous_obs_seq % $save_freq == $mod_save) save to mass store
# where in the obs_seq loop below j = current obs_seq file - 1
# where in the obs_seq loop below j = current obs_seq file - 1
set save_freq = 1
set mod_save = 1

# END of run parameters to change
#==========================================================================================


# This is the CENTRAL directory for whole perfect_model_obs job
#    jobs submitted to batch queues from here
#    I/O between perfect_model_obs and advance_model and assim_region goes through here.
#    perfect_model_obs output is put here, then moved to final storage.
cd ${CENTRALDIR}
set myname = $0                              # this is the name of this script
set MASTERLOG = ${CENTRALDIR}/run_job.log    # Set Variable for a 'master' logfile 
${REMOVE} ${MASTERLOG}                       # clean up old links

#----------------------------------------------------------
set num_ens = 1

echo "There are ${num_ens} ensemble members."
echo "There are ${num_ens} ensemble members."  > $MASTERLOG

#---------------------------------------------------------
# Get executable programs and scripts from DART-CAM directories

if (! -x perfect_model_obs) then
   ${COPY} ${DARTCAMDIR}/work/perfect_model_obs                     .
   ${COPY} ${DARTCAMDIR}/work/trans_pv_sv                .
   ${COPY} ${DARTCAMDIR}/work/trans_sv_pv                .
   ${COPY} ${DARTCAMDIR}/work/trans_time                 .
endif
if (! -e advance_model.csh) then
   ${COPY} ${DARTCAMDIR}/shell_scripts/advance_model.csh     .
   ${COPY} ${DARTCAMDIR}/shell_scripts/run-cam.csh            .
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

# clean up old CAM inputs  and obs_seq.out that may be laying around

if (-e caminput_1.nc) then
   ${REMOVE} clminput_[1-9]*.nc 
   ${REMOVE} caminput_[1-9]*.nc 
endif
if (-e obs_seq.out) ${REMOVE} obs_seq.out

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
# where xx_  = output_root before 2008 signified the month OF THE OBS_SEQ_FIRST.
# Now xx is always 'obs'.
# ##[#] signified the 2+ digit obs_seq number within this experiment before 2008.
# Now it's always 4 digits.
set output_root = obs


#============================================================================================
# Have an overall outer loop over obs_seq.in files

set i = $obs_seq_1
while($i <= $obs_seq_n) ;# start i/obs_seq loop
echo ' '
echo ' ' >> $MASTERLOG
echo '----------------------------------------------------'
echo '----------------------------------------------------' >> $MASTERLOG
echo "starting observation sequence file $i at "`date`
echo "starting observation sequence file $i at "`date` >> $MASTERLOG

#===================================================================================
# Each iteration of this loop will run perfect_model_obs for one obs_seq.in file

# Construct directory name of location of restart files
@ j = $i - 1
set out_prev = `printf "%s_%04d" ${output_root} $j`

set start_time = `date  +%s`
#-----------------------
# Get perfect_model_obs input files
# The first one is different than all the rest ...

if ($i == $obs_seq_first) then
      ${COPY} ${input}${obs_seq_first}.nml input.nml 
else
      if (-e            ${input}n.nml) then
         ${COPY}  ${input}n.nml    input.nml 
      else if (-e ${input}${i}.nml) then
         ${COPY}  ${input}${i}.nml input.nml 
      else
         echo "input_next is MISSING" >> $MASTERLOG  
         echo "input_next is MISSING"
         exit 19 
      endif
endif

   
#----------------------------------------------------------------------
# Get obs_seq file for this assimilation based on date, derived from
# the obs_seq NUMBER of this iteration ("i").
# At this time this script assumes that assimilation experiments start
# at the beginning of the chosen month.
# As such, we need to ensure any existing obs_seq.in is GONE, which is
# a little scary.
#
# UGLY; need more general calendar capability in here.
#----------------------------------------------------------------------

${REMOVE} obs_seq.in 

set seq = $i
if ($seq == 0) then
      set OBS_SEQ = ${CENTRALDIR}/obs_seq.in
else if ($obs_seq_freq < -15) then
      @ month = $i + 1
      set day = 01
      
      set OBS_SEQ = ${obs_seq_root}${month}${day}

else if ($obs_seq_freq == -15) then
      @ month =  ($i / 2) + 1
      @ day   = (($i % 2) * 14) + 1
      if ($month > 12) then
         set month = $month % 12
         set obs_seq_root = ${CENTRALDIR}/obs_seq2002
      endif
      if ($month < 10) set month = 0$month
      if ($day   < 10) set day   = 0$day
      
      set OBS_SEQ = ${obs_seq_root}${month}${day}

      echo "obs_seq, root, month, day = $i $obs_seq_root $month $day " >> $MASTERLOG

else if ($obs_seq_freq > 0) then
      @ month = $mo - 1
      while ($month >= $mo_first)
          @ seq = $seq - $days_in_mo[$month] * $obs_seq_freq
          @ month = $month - 1
      end
      @ month = $month + 1
      if ($month < 10) set month = 0$month

      @ day = (($seq - 1) / $obs_seq_freq) + 1
      if ($day < 10) set day = 0$day

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
      ${LINK} $OBS_SEQ  obs_seq.in 
else
      echo "ERROR - no obs_seq.in for $i - looking for:" 
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
set cam_init =  `pwd`/$exp/${out_prev}/CAM/caminput_
set clm_init =  `pwd`/$exp/${out_prev}/CLM/clminput_
endif

# transmit info to advance_model.csh
# The second item echoed must be the directory where CAM is kept
echo "$exp "                                            >! casemodel.$i
echo "$CAM_src "                                         >> casemodel.$i
echo "$cam_init "                                       >> casemodel.$i
echo "$clm_init"                                        >> casemodel.$i
echo "$parallel_cam"                                    >> casemodel.$i
echo "$run_command"                                     >> casemodel.$i
# advance_model wants to see a file 'casemodel' and not keep track of which obs_seq it's for
${REMOVE} casemodel
${LINK} casemodel.$i casemodel 

# link to perfect_ic_old file(s), so that perfect_model_obs can copy them to a compute node
if (-e ${from_root}/perfect_ic_old) then 
${REMOVE} perfect_ic_old 
   endif 
   ${LINK} ${from_root}/perfect_ic perfect_ic_old 
   echo ' '
   echo "job- perfect_ic_old is/are" 
   ls -lt perfect_ic_old* 

   #-----------------------------------------------------------------------------
   # link local CAM input files to generic names in CENTRALDIR.  
   # These just provide grid info to perfect_model_obs, not state info.
   # CAM namelist input file; will be augmented with model advance time and
   # domain decomposition info by run-cam.csh
   #-----------------------------------------------------------------------------
   
   ${REMOVE} caminput.nc clminput.nc 
   ${LINK} ${CAM_src}/caminput.nc caminput.nc
   ${LINK} ${CAM_src}/clminput.nc clminput.nc
   
   #-----------------------------------------------------------------------------
   # get name of file containing PHIS from the CAM namelist.  This will be used by
   # static_init_model to read in the PHIS field, which is used for height obs.
   #-----------------------------------------------------------------------------
   ${REMOVE} namelistin
   ${LINK} ${CAM_src}/namelistin namelistin
   sleep 1 
   if (! -e namelistin ) then
      echo "ERROR ... need a namelistin file." >> $MASTERLOG
      echo "ERROR ... need a namelistin file."
      exit 89
   endif
   if (! -e ${exp}/namelistin) ${COPY} namelistin ${exp}/namelistin  

   # Observations on heights require PHIS, which are on a CAM surface fields file.
   # Only need to do this for obs_seq_1. 
   rm cam_phis.nc
   if (-e $CAM_phis) then
      ${COPY} $CAM_phis cam_phis.nc
# ?      chmod 744 cam_phis.nc
   else
      echo "ERROR ... need a cam_phis file from CAM h0 history file." >> $MASTERLOG
      echo "ERROR ... need a cam_phis file from CAM h0 history file."
         exit 99
   endif
      
   #======================================================================

   # runs perfect_model_obs, which tells the model to model advance and assimilates obs
   
   # Run the perfect_model_obs in async=2 mode.
   # runs perfect_model_obs, which tells the model to model advance and assimilates obs
   ${run_command} ./perfect_model_obs 

   set KILLCOMMAND = 'touch BOMBED; exit'

   #-----------------------------------------------
   # Move the output to storage after perfect_model_obs signals completion.
   # At this point, all the restart,diagnostic files are in the CENTRALDIR
   # and need to be moved to the 'experiment permanent' directory.
   # We have had problems with some, but not all, files being moved
   # correctly, so we are adding bulletproofing to check to ensure the filesystem
   # has completed writing the files, etc. Sometimes we get here before
   # all the files have finished being written.

   
   echo "Listing contents of CENTRALDIR before archiving at "`date` 
   ls -l 
   echo "Listing contents of CENTRALDIR before archiving at "`date` >> $MASTERLOG
   ls -l >> $MASTERLOG

#-----------------------------------------------------------------------------
# Ensure the (output) experiment directories exist
# All the  CAM-related files will get put in ${exp}/${output_dir}/CAM
# All the  CLM-related files will get put in ${exp}/${output_dir}/CLM
# All the DART-related files will get put in ${exp}/${output_dir}/DART
#-----------------------------------------------------------------------------
   set output_dir = `printf "%s_%04d" ${output_root} $i`

   mkdir -p ${exp}/${output_dir}/{CAM,CLM,DART} 

      if ( -e True_State.nc && ! -z True_State.nc) then 
         ${MOVE} True_State.nc ${exp}/${output_dir} 
         if ( ! $status == 0 ) then 
            echo "failed moving ${CENTRALDIR}/True_State.nc" >> $MASTERLOG 
            echo "failed moving ${CENTRALDIR}/True_State.nc" 
            eval $KILLCOMMAND 
         endif           
      else               
         echo "... THUMP ... ${CENTRALDIR}/True_State.nc does not exist and should"
         echo "directory contents follow" 
         ls -l 
         eval $KILLCOMMAND 
      endif 

      if ( -e obs_seq.out && ! -z obs_seq.out) then 
         set date_name = $OBS_SEQ:t
         ${MOVE} obs_seq.out ${exp}/${output_dir}/syn_$date_name 
         if ( ! $status == 0 ) then 
            echo "failed moving ${CENTRALDIR}/obs_seq.out" >> $MASTERLOG 
            echo "failed moving ${CENTRALDIR}/obs_seq.out" 
            eval $KILLCOMMAND 
         endif           
      else               
         echo "... THUMP ... ${CENTRALDIR}/obs_seq.out does not exist and should"
         echo "directory contents follow" 
         ls -l 
         eval $KILLCOMMAND 
      endif 

   # inflate_diag may or may not exist (when do_obs_inflate=F do_single_ss_inflate=F)
   # so don't die if it's missing.  We'd have to query input.nml to learn if it should exist.

   # Move the perfect_model_obs restart file(s) to the storage subdirectory
   
   echo "moving perfect_ic_new to ${exp}/${output_dir}/DART/perfect_ic" 
   if (-e perfect_ic_new) then 
      ${MOVE} perfect_ic_new ${exp}/${output_dir}/DART/perfect_ic 
      if (! $status == 0 ) then 
         echo "failed moving perfect_ic_new to ${exp}/${output_dir}/DART/perfect_ic" 
         eval $KILLCOMMAND 
      endif 
   else 
      echo "NO perfect_ic_new FOUND" 
      echo "NO perfect_ic_new FOUND" 
      echo "NO perfect_ic_new FOUND" 
      echo "NO perfect_ic_new FOUND" 
      eval $KILLCOMMAND 
   endif 
    
   set n = 1 
   while ($n <= ${num_ens})    ;# loop over all ensemble members 
      set CAMINPUT = caminput_${n}.nc  
      set CLMINPUT = clminput_${n}.nc  
    
      if ( -e $CAMINPUT && ! -z $CAMINPUT) then 
         ${MOVE} $CAMINPUT ${exp}/${output_dir}/CAM  
         if (! $status == 0 ) then 
            echo "failed moving ${CENTRALDIR}/$CAMINPUT " 
            eval $KILLCOMMAND 
         endif 
      else 
         echo "${CENTRALDIR}/$CAMINPUT does not exist and maybe should." 
# kdr debug; this should be dependent on whether a model advance was needed before
#            the assim.
#         eval $KILLCOMMAND 
      endif 
    
      if ( -e $CLMINPUT && ! -z $CLMINPUT) then 
         ${MOVE} $CLMINPUT ${exp}/${output_dir}/CLM 
         if (! $status == 0 ) then 
            echo "failed moving ${CENTRALDIR}/$CLMINPUT " 
            eval $KILLCOMMAND 
         endif 
      else 
         echo "${CENTRALDIR}/$CLMINPUT does not exist and maybe should." 
# kdr debug; this should be dependent on whether a model advance was needed before
#            the assim.
#         eval $KILLCOMMAND 
      endif 
    
     @ n++ 
   end 

   if (! -e times ) then  
      # There may have been no advance; 
      # use the CAM_ics_1 as the CAM ics for this time  
      set ens = 1  
      while ($ens <= "$num_ens" )  
         cp ${CAM_ics_1}${ens}.nc ${exp}/${output_dir}/CAM  
         cp ${CLM_ics_1}${ens}.nc ${exp}/${output_dir}/CLM  
         @ ens++  
      end  
   endif  

   # test whether it's safe to end this obs_seq_ 
    
   if (  -e $exp/${output_dir}/DART/perfect_ic && \
       ! -z $exp/${output_dir}/DART/perfect_ic   ) then
       echo "it is OK to proceed with next obs_seq at "`date` >> $MASTERLOG 
       echo "it is OK to proceed with next obs_seq at "`date` 
   else 
       echo "RETRIEVE perfect_ic files from perfect_model_obs temp directory ?" 
       echo "Then remove temp and cam advance temps" 
       echo "RETRIEVE perfect_ic files from perfect_model_obs temp directory ?" >> $MASTERLOG 
       echo "Then remove temp and cam advance temps" >> $MASTERLOG 
       exit 
   endif 

# save a representative model advance 
    
   ${MOVE} input.nml           ${exp}/${output_dir} 
   ${MOVE} casemodel.$i        ${exp}/${output_dir} 
# stdout   from perfect_model_obs.f90 
   ${MOVE} dart_log.out        ${exp}/${output_dir} 
   ${MOVE} cam_out_temp1       ${exp}/${output_dir} 
   ${MOVE} namelist            ${exp}
   ${REMOVE} cam_out_temp* *_ud* *_ic[0-9]* *_ic_old* 

   set end_time = `date  +%s`
   @ length_time = $end_time - $start_time
   echo "duration = $length_time" 


   echo "completed iteration $i ($OBS_SEQ) at "`date`
   @ i++
end # end of the huge "i" loop

${REMOVE} ~/lnd.*.rpointer

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

