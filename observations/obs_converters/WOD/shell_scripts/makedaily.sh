#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# split the monthly file into "daily" files which start at 12:01Z 
# the previous day and end at 12:00Z on the day that matches the 
# day in the filename.
#
#BSUB -J splitobs
#BSUB -W 4:00
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q share
#BSUB -P 93300315
#BSUB -n 1

# set the first and last days to be split.  can roll over
# month and year boundaries now!   note that for the first day
# you need the previous day data available.

let start_year=1998
let start_month=1
let start_day=2

let end_year=1999
let end_month=12
let end_day=31

EXEDIR='../work'
DATDIR='../data'

# end of things you should have to set in this script

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
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ${EXEDIR}/advance_time`)
#echo $end_d

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ${EXEDIR}/advance_time`)
#echo $start_d

# these are a string in the format YYYYMMDDHH
# do them here to prime the loop below which first takes them apart.
currday=(`echo ${start_year}${mon2}${day2}00   0 | ${EXEDIR}/advance_time`)
nextday=(`echo ${start_year}${mon2}${day2}00 +1d | ${EXEDIR}/advance_time`)
prevday=(`echo ${start_year}${mon2}${day2}00 -1d | ${EXEDIR}/advance_time`)

# how many total days are going to be merged (for the loop counter)
# (pull out the first of the 2 numbers which are output from advance_time)
let totaldays=${end_d}-${start_d}+1
#echo $totaldays

# loop over each day
let d=1
while (( d <= totaldays)) ; do
#echo top of loop
#echo $currday $nextday
  # parse out the parts from a string which is YYYYMMDDHH
  # both for the current day and tomorrow
  cyear=${currday:0:4}
  cmonth=${currday:4:2}
  cday=${currday:6:2}
  nyear=${nextday:0:4}
  nmonth=${nextday:4:2}
  nday=${nextday:6:2}
  pyear=${prevday:0:4}
  pmonth=${prevday:4:2}
  pday=${prevday:6:2}
#echo curr $cyear $cmonth $cday
#echo next $nyear $nmonth $nday
#echo prev $pyear $pmonth $pday

  # compute the equivalent gregorian days here.
  g=(`echo ${cyear}${cmonth}${cday}00 -1d -g | ${EXEDIR}/advance_time`)
  greg0=${g[0]}
  g=(`echo ${cyear}${cmonth}${cday}00 0 -g | ${EXEDIR}/advance_time`)
  greg1=${g[0]}
  g=(`echo ${cyear}${cmonth}${cday}00 +1d -g | ${EXEDIR}/advance_time`)
  greg2=${g[0]}
#echo $greg0 $greg1 $greg2

  echo starting WOD obs for ${cyear}${cmonth}${cday} gregorian= $greg1

  # last 12 hrs yesterdays data plus the first 12 hrs of todays
  if [[ ${cmonth} == '01' && ${cday} == '01' ]] ; then
     echo "${DATDIR}/obs_seq${pyear}.wod" > olist
     echo "${DATDIR}/obs_seq${cyear}.wod" >> olist
  else
     echo "${DATDIR}/obs_seq${cyear}.wod" > olist
  fi

  sed -e "s/YYYY/${cyear}/g"    \
      -e "s/MM/${cmonth}/g"     \
      -e "s/PP/${pmonth}/g"     \
      -e "s/DD/${cday}/g"       \
      -e "s/GREG0/${greg0}/g"  \
      -e "s/GREG1/${greg1}/g"  \
      -e "s/GREG2/${greg2}/g"    < ./input.nml.template > input.nml

  # make sure output dir exists
  if [[ ! -d ../${cyear}${cmonth} ]] ; then
     mkdir ../${cyear}${cmonth}
  fi

  # do the extract here
  ${EXEDIR}/obs_sequence_tool
  

  # advance the day; the output is YYYYMMDD00
  prevday=$currday
  currday=$nextday
  nextday=(`echo ${nyear}${nmonth}${nday}00 +1d | ${EXEDIR}/advance_time`)
#echo $currday $nextday $prevday

  # advance the loop counter
  let d=d+1
#echo d=$d
 
done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

