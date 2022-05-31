#!/usr/bin/env bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# split the monthly file into "daily" files which start at 12:01Z 
# the previous day and end at 12:00Z on the day that matches the 
# day in the filename.
#
#===============================================================================
# LSF directives                bsub < test_batch.csh
#
#BSUB -J splitobs
#BSUB -W 4:00
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q geyser
#BSUB -P P86850054
#BSUB -n 1
#===============================================================================
# PBS directives                qsub test_batch.csh
#
#PBS -N testpbs         
#PBS -l walltime=00:01:10
#PBS -q economy
#PBS -l select=1:ncpus=1
#PBS -A P86850054

# set the first and last days to be split.  can roll over
# month and year boundaries now!   note that for the first day
# you need the previous day data available.

let start_year=2017
let start_month=1
let start_day=1

let end_year=2017
let end_month=6
let end_day=30

EXEDIR='../work'
DATDIR='../data/iowa'
OUTDIR=/glade/scratch/${USER}/SMAP_obs

# The actual dates in the observation files need to be replaced
# with the similar dates from 2005. Or not.

MODEL_YEAR=2017

# end of things you should have to set in this script

# make sure there is an initial input.nml for advance_time
# input.nml gets overwritten in the subsequent loop.

\cp -f ${EXEDIR}/input.nml.template input.nml || exit -1

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# these outputs from advance time (with the -g flag) are
# 2 integers: gregorian_day_number seconds
# and since we don't set hours, minutes, or seconds, the second
# number is always 0 and uninteresting for us.
mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
end_d=(`echo ${end_year}${mon2}${day2} 0 -g | ${EXEDIR}/advance_time`)
echo "${end_year}${mon2}${day2} is $end_d"

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2} 0 -g | ${EXEDIR}/advance_time`)
echo "${start_year}${mon2}${day2} is $start_d"

# these are a string in the format YYYYMMDDHH
# do them here to prime the loop below which first takes them apart.
currday=(`echo ${start_year}${mon2}${day2}  0 | ${EXEDIR}/advance_time`)
echo $currday

# how many total days are going to be merged (for the loop counter)
# (pull out the first of the 2 numbers which are output from advance_time)
let totaldays=${end_d}-${start_d}+1
echo $totaldays

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  cyear=${currday:0:4}
  cmonth=${currday:4:2}
  cday=${currday:6:2}

  echo ' '
  echo "day $d ... $cyear $cmonth $cday"

  # make sure output dir exists
  OBSDIR=${OUTDIR}/${MODEL_YEAR}${cmonth}
  if [[ ! -d $OBSDIR ]] ; then
     mkdir -p $OBSDIR
  fi

  # compute the equivalent gregorian days here.
  # gp is the day before the day of interest
  # gn is the day of interest
  gp=(`echo ${cyear}${cmonth}${cday}00 -43199s -g | ${EXEDIR}/advance_time`)
  gn=(`echo ${cyear}${cmonth}${cday}00 +43200s -g | ${EXEDIR}/advance_time`)
  gp_day=${gp[0]}
  gp_sec=${gp[1]}
  gn_day=${gn[0]}
  gn_sec=${gn[1]}

  f0=(`echo ${MODEL_YEAR}${cmonth}${cday}00 -1d -g | ${EXEDIR}/advance_time`)
  f1=(`echo ${MODEL_YEAR}${cmonth}${cday}00 +0d -g | ${EXEDIR}/advance_time`)

  fake0=${f0[0]}
  fake1=${f1[0]}

  echo "start = $gp_day $gp_sec  fake day is $fake0"
  echo "end   = $gn_day $gn_sec  fake day is $fake1"

  echo "converting SMAP ${cyear} ${cmonth} ${cday} : gregorian day $greg1"

  today=${cyear}${cmonth}${cday}
  obsdate=${MODEL_YEAR}-${cmonth}-${cday}-00000
  yyyymm=${cyear}${cmonth}

  sed -e "s/DART1D/${gp_day}/g" \
      -e "s/DART1S/${gp_sec}/g" \
      -e "s/DARTND/${gn_day}/g" \
      -e "s/DARTNS/${gn_sec}/g" < ../work/input.nml.template > input.nml

  # do the extract here
  # SMAP_L2_to_obs  reads in a list of files, writes out 'obs_seq.out'
  # obs_sequence_tool  reads in a 'obs_seq.out', writes 'obs_seq.0Z'
  # since we are not using the obs_sequence tool this time, ... just rename
  # file to expected name for the remainder of the processing

  previous=(`echo ${cyear}${cmonth}${cday} -1d | ${EXEDIR}/advance_time`)
  yesterday=${previous:0:8}

  ls -1 ${DATDIR}/${yyyymm}/*_D_${yesterday}* >  file_list.txt
  ls -1 ${DATDIR}/${yyyymm}/*_D_${today}*     >> file_list.txt

  ${EXEDIR}/SMAP_L2_to_obs
#  ${EXEDIR}/obs_sequence_tool
  if [ -f obs_seq.out ]; then
     mv   obs_seq.out obs_seq.0Z
  fi

  # sometimes there are no observations for this time 

  if [ -f obs_seq.0Z ]; then
     # replace the actual date with the model date

     sed -e "s/${gp_day}/${fake0}/" \
         -e "s/${gn_day}/${fake1}/" obs_seq.0Z > obs_seq.${obsdate}

     \mv -v obs_seq.${obsdate} ${OBSDIR}
  else
     # create an and empty observation sequence file
     touch ${OBSDIR}/obs_seq.${obsdate}
  fi

  \rm -f obs_seq.out obs_seq.0Z

  # advance the day; the output is YYYYMMDD00
  prevday=$currday
  currday=(`echo ${cyear}${cmonth}${cday}00 +1d | ${EXEDIR}/advance_time`)

  # advance the loop counter
  let d=d+1
echo d=$d
 
done

exit 0


