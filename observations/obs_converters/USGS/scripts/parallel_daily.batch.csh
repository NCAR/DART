#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DESCRIPTION:
#
# Generate daily streamflow observation sequence files.
# Submit this script
# to your batch system and it will invoke the 'makedaily.csh'
# script once for each conversion day.
# this one does N conversions in parallel from a command script.
#
#==========================================================================
# PBS directives                qsub paralell_daily.batch
#
# qstat    information on running jobs
# qsub     submit a job
# qdel     kill a job
# qpeek    see output from a running job
#
#PBS -m abe
#PBS -M thoar@ucar.edu         
#PBS -N streamflow         
#PBS -l walltime=00:20:00
#PBS -q regular
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -A P86850054
#
#==========================================================================
# USER SETTINGS HERE

set  start_year = 2016
set start_month = 8
set   start_day = 31

set  end_year = 2016
set end_month = 9
set   end_day = 30 

set DARTDIR = `pwd`
set  obsdir = /glade/scratch/arezoo/dart/matthew/Reforecase/usgs_timeslices_All
set workdir = /glade/scratch/thoar/streamflow

# set things that vary between batch systems here.

if ($?SLURM_JOB_ID) then
  echo 'running SLURM'
  setenv SLURM true
  source /glade/u/apps/opt/slurm_init/init.csh
  setenv RUNCMD "srun --multi-prog"
  setenv JOBNAME $SLURM_JOB_NAME

  # should match the mpiprocs=X setting above
  # todo ... there is an enviroment variable for this ...
  set njobs = 36

else if ($?PBS_NODEFILE) then

  echo 'running PBS'
  setenv PBS true
  setenv MPI_SHEPHERD true
  setenv RUNCMD "mpiexec_mpt launch_cf.sh"
  setenv JOBNAME $PBS_JOBNAME

  set njobs = `cat $PBS_NODEFILE | wc -l`
  echo "njobs from PBS_NODEFILE is $njobs"

else if ($?LSB_HOSTS) then
  echo 'running LSF'
  setenv LSF true
  setenv MP_PGMMODEL mpmd
  setenv RUNCMD "mpirun.lsf -cmdfile"
  setenv JOBNAME $LSB_OUTPUTFILE:ar
  # should match the mpiprocs=X setting above
  # todo ... there is an enviroment variable for this ...
  set njobs = 36

else
  echo 'running without a batch system'
  setenv NOBATCH true
  setenv RUNCMD "csh"
  setenv JOBNAME "streamflow"
  set njobs = 36

endif

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR

mkdir -p $workdir
cd $workdir

# make the converter script happy by making these files consistent
# and link to a working advance_time

ln -sf ${DARTDIR}/../work/advance_time
ln -sf ${DARTDIR}/../work/input.nml
cp  -f ${DARTDIR}/makedaily.csh .

echo job started at `date`

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# the output of advance time with the -g input is:
#   gregorian_day_number  seconds
# use $var[1] to return just the day number

set mon2=`printf %02d $end_month`
set day2=`printf %02d $end_day`
set end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ./advance_time`)

set mon2=`printf %02d $start_month`
set day2=`printf %02d $start_day`
set start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ./advance_time`)

# the output of this call is a string YYYYMMDDHH
# see below for help in how to easily parse this up into words
set curday=`echo ${start_year}${mon2}${day2}00 0 | ./advance_time`

# how many total days are going to be processed (for the loop counter)
# note that the parens below are necessary; otherwise the computation
# does total = end - (start+1), or total = end - start - 1, which is
# not how elementary math usually works on a computer (left to right
# evaluation of ops with equal priority is most common.)
@ totaldays = ( $end_d[1] - $start_d[1] ) + 1

# loop over each day
set d=1
while ( $d <= $totaldays )

   rm -f mycmdfile

   set j=1
   while ( $j <= $njobs && $d <= $totaldays)

      # parse out the parts from a string which is YYYYMMDDHH
      # use cut with the byte option to pull out columns 1-4, 5-6, and 7-8
      # then bc to strip off leading blanks
      set  year=`echo $curday | cut -b1-4`
      set month=`echo $curday | cut -b5-6`
      set   day=`echo $curday | cut -b7-8`
  
      # numeric month/day (no leading 0)
      set nmonth=`echo $month | bc`
      set   nday=`echo $day | bc`
  
      # status/debug - comment in or out as desired.
      echo starting processing for ${year} ${nmonth} ${nday}
  
      # this conversion date
      set tdate = ${year}${month}${day}
  
      # create a command file where each line calls a script with
      # unique arguments, including a unique work directory name.
      echo "csh ./makedaily.csh ${tdate} ${obsdir} ${tdate}" >> mycmdfile
  
      # advance the day; the output is YYYYMMDD00
      set curday=`echo ${year}${month}${day}00 +1d | ./advance_time`
  
      # advance the loop counter
      @ d++
  
      # advance the concurrent job counter
      @ j++
   end

   # the system seems to want the same number of commands in each
   # invocation of the command file as there are cpus on the node.
   # if we aren't running an even multiple of real tasks compared
   # to the cpu count, pad the rest of the script with a call to 'date'

   # avoid echoing the filename by making wc read stdin
   set j=`cat mycmdfile | wc -l`
   while ( $j < $njobs )
     echo "date " >> mycmdfile
     @ j++
   end
 
   echo "running jobs:"
   cat mycmdfile

   # actually launch the jobs here
   $RUNCMD ./mycmdfile >& run_output_${tdate}.txt

end

echo job ended at `date`

exit 0


