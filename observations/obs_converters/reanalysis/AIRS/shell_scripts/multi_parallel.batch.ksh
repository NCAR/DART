#!/bin/ksh -x
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: multi_parallel.lsf 9948 2016-03-03 22:30:57Z nancy $
#
#--------------------------------------------------------------
# DESCRIPTION:
# Convert AIRS T/Q obs from HDF-EOS format to DART format.
# creates daily files from up to 240 individual 6-minute files.
#
# Driver script for the parallel version.  Submit this script
# to your batch system and it will invoke the 'dodaily.sh'
# script once for each conversion day.
#
# this one does N conversions in parallel from a command script.
#
#--------------------------------------------------------------

#==========================================================================
# SLURM directives             sbatch script.csh
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submit a job
# scancel   kill a job
#
#SBATCH --ignore-pbs
#SBATCH --job-name=airsobs1
#SBATCH -n 36
#SBATCH --ntasks-per-node=36
#SBATCH --time=02:10:00
#SBATCH -A P86850054
#SBATCH -p dav
#SBATCH -C casper
#SBATCH -e gpsobs1.%j.err
#SBATCH -o gpsobs1.%j.out
#
#==========================================================================
# PBS directives                qsub script.csh
#
# qstat    information on running jobs
# qsub     submit a job
# qdel     kill a job
# qpeek    see output from a running job
#
#PBS -N airsobs2020         
#PBS -l walltime=02:00:00
#PBS -q regular
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -A P86850054
#
#==========================================================================
# LSF directives                bsub < script.csh
#
# bstat    information on running jobs
# bsub     submit a job
# bdel     kill a job
# bpeek    see output from a running job
#
#BSUB -J airsobs1
#BSUB -o airsobs1.%J.log
#BSUB -q small
#BSUB -n 16
#BSUB -W 0:10:00
#BSUB -P P86850054
#
#==========================================================================

# USER SETTINGS HERE

#--------------------------------------------------------------

# set the first and last days.  can roll over month and year boundaries.
#
# set these WITHOUT leading 0s.  they are numeric values, not strings.
#

TARGETDIR=../AIRS_24_sub9x10_nceperrs

let  start_year=2020
let start_month=1
let   start_day=1

let    end_year=2020
let   end_month=1
let     end_day=31

# should match the mpiprocs=X setting above
let njobs=36

# END USER SETTINGS

# set things that vary between batch systems here.

if [ "$SLURM_JOB_ID" -ne "" ] ; then
  echo running SLURM
  SLURM=true
  RUNCMD="srun --multi-prog"
elif [ "$PBS_NODEFILE" -ne "" ] ;  then
  echo running PBS
  PBS=true
  export MPI_SHEPHERD=true
  RUNCMD="mpiexec_mpt launch_cf.sh"
elif [ "$LSB_HOSTS" -ne "" ] ; then
  echo running LSF
  LSF=true
  export MP_PGMMODEL=mpmd
  RUNCMD="mpirun.lsf -cmdfile"
else
  echo running without a batch system
  NOBATCH=true
  RUNCMD="csh"
fi

TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

echo job started at `date`

if [ ! -f advance_time ]; then
  echo error: advance_time must exist in the current directory
  exit 1
fi

# before running advance_time there must be an input.nml
# with at least a utilities namelist.  it will be overwritten
# in the first real loop by the input.nml.template, so make
# any permanent changes the template and not to input.nml.

rm -f input.nml
echo '&utilities_nml' > input.nml
echo '/' >> input.nml

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# the output of advance time with the -g input is:
#   gregorian_day_number  seconds
# use ${var[0]} to return the day number, ${var[1]} for the seconds

step=1d

mon2c=`printf %02d $end_month`
day2c=`printf %02d $end_day`
hor2c=00
endstring=${end_year}${mon2c}${day2c}${hor2c}
#echo endstring $endstring
set -A end_t `echo $endstring 0 -g | ./advance_time`

#echo end_t $end_t

mon2c=`printf %02d $start_month`
day2c=`printf %02d $start_day`
hor2c=00
startstring=${start_year}${mon2c}${day2c}${hor2c}
#echo startstring $startstring
set -A start_t `echo $startstring 0 -g | ./advance_time`

#echo start_t $start_t

# count up how many steps we are going to be executing

curtime=$startstring
let totalsteps=0

while true ; do

  # compute the current gregorian day here.  
  #  this_t[0] is day, this_t[1] is seconds
  set -A this_t `echo $curtime 0 -g | ./advance_time`

#echo this_t $this_t
  # if the current day is beyond the end day, we are done.  break out of loop.
  if [[ ${this_t[0]} -gt ${end_t[0]} ]] ; then break ; fi

  # advance to the next time and increment the step count
  curtime=`echo $curtime $step | ./advance_time`
  let totalsteps=$totalsteps+1
done


echo This script is going to create $totalsteps DART obs sequence files.
echo Starting at $start_year $start_month $start_day 
echo Ending at $end_year $end_month $end_day 
echo

# the output of this call is a string YYYYMMDDHH
# see below for help in how to easily parse this up into words
  curstep=`echo $startstring         0 | ./advance_time`
 nextstep=`echo $startstring     $step | ./advance_time`

#echo before loop
#echo curstep nextstep $curstep $nextstep

# ok, let's actually do something.  up to now it has
# been computing fiddly time bits.

let s=1
while [[ $s -le $totalsteps ]] ; do

  rm -f mycmdfile

  let j=1
  while [[ $j -le $njobs && $s -le $totalsteps ]] ; do

    # parse out the parts from a string which is YYYYMMDDHH
    # use cut with the byte option to pull out columns 1-4, 5-6, 7-8, and 9-10
    # c = current day, n = next day
  
     cyear=`echo $curstep | cut -b1-4`
    cmonth=`echo $curstep | cut -b5-6`
      cday=`echo $curstep | cut -b7-8`
     chour=`echo $curstep | cut -b9-10`
  
     nyear=`echo $nextstep | cut -b1-4`
    nmonth=`echo $nextstep | cut -b5-6`
      nday=`echo $nextstep | cut -b7-8`
     nhour=`echo $nextstep | cut -b9-10`
  
    # compute the equivalent gregorian days/secs here.
    set -A g `echo ${cyear}${cmonth}${cday}${chour} 0 -g | ./advance_time`
    cgregday=${g[0]}; cgregsec=${g[1]}
  
    # special for dart: step start needs to be +1 second 
    set -A g `echo ${cyear}${cmonth}${cday}${chour} 0 -g | ./advance_time`
    ssgregday=${g[0]}; ssgregsec=${g[1]}
  
    set -A g `echo ${cyear}${cmonth}${cday}${chour} $step -g | ./advance_time`
    segregday=${g[0]}; segregsec=${g[1]}
  
    # status/debug - comment in or out as desired.
#    echo "starting processing for "  ${cyear}  ${cmonth}  ${cday}  ${chour}, gregorian day/sec:  $cgregday  $cgregsec
#    echo 
  
  
    # construct often-used strings below
       ym=${cyear}${cmonth}   ;   ymd=${ym}${cday}      ;   ymdh=${ymd}${chour}
  
#  #echo
#  echo     ym=$ym ; echo   ymd=${ymd}   ; echo   ymdh=${ymdh}
#  echo
  
    # ----------------------------------------
    # start of this merge-specific code
  
  
    # this conversion date and year/month
    tmonth=${cyear}${cmonth}
    tdate=${cyear}${cmonth}${cday}

    # finally, the work this script is actually doing:
    echo converting files for day $cday

    # make a list of the 240 filenames for this day
    rm -f flist
    ls ../Daily_raw_files/$tmonth/$tdate/AIRS.${cyear}.${cmonth}.${cday}.*.L2.RetStd_IR.v6.*.hdf >> flist
#cat flist
  
    # where to do the work.
    workdir=workdir_${tdate}
  
    # make sure the output dir exists
    mkdir -p $TARGETDIR/$tmonth

    # fix up the output filename to have the right date.
    # the converter will write the file to the correct location - no additional move is needed.
    # (do YYYYMMDD first before YYYYMM or this won't work.)
    sed -e "s;OUTDIR;../$TARGETDIR;g" -e "s/YYYYMMDD/$tdate/g" -e "s/YYYYMM/$tmonth/g" input.nml.template > input.nml

    mkdir -p $workdir
    sed -e 's;^;../;' flist > $workdir/flist
    cp input.nml $workdir/
    ( cd $workdir; ln -sf ../convert_airs_L2 . )

    # create a command file where each line calls a script with
    # unique arguments, including a unique work directory name.
    echo "ksh ./dodaily.sh $workdir" >> mycmdfile
    
    # ----------------------------------------

    # advance the step; the output is YYYYMMDDHH
    curstep=$nextstep
    nextstep=`echo $curstep     $step | ./advance_time`
   
#echo bot of loop
#echo curstep nextstep $curstep $nextstep

    # advance the loop counter
    let s=$s+1
   
    # advance the concurrent job counter
    let j=$j+1
  done

  echo running jobs:
  cat mycmdfile

  # the system seems to want the same number of commands in each
  # invocation of the command file as there are cpus on the node.
  # if we aren't running an even multiple of real tasks compared
  # to the cpu count, pad the rest of the script with a call to 'date'

  # avoid echoing the filename by making wc read stdin
  let j=`cat mycmdfile | wc -l`
  while [[ $j -lt $njobs ]] ; do
    echo "sleep 1; date " >> mycmdfile
    let j=$j+1
  done

  # actually launch the jobs here
  $RUNCMD ./mycmdfile 

done


echo job ended at `date`

exit 0

# <next few lines under version control, do not edit>
# $URL: https://subversion.ucar.edu/DAReS/DART/trunk/observations/NCEP/prep_bufr/work/multi_parallel.lsf $
# $Revision: 9948 $
# $Date: 2016-03-03 15:30:57 -0700 (Thu, 03 Mar 2016) $

