#!/bin/bash
# BLUEFIRE /usr/local/bin/bash
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# split a yearly file into "daily" files which start at 00:00Z
# the previous day and end at 23:59Z on the day that matches the
# day in the filename.

# set the first and last days to be split.  can roll over
# month and year boundaries now!   note that for the first day
# the data from the previous day must be available.

let start_year=2000
let start_month=01
let start_day=01

let end_year=2000
let end_month=12
let end_day=31

# end of things you should have to set in this script

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# make sure there is an initial input.nml for advance_time
cp -f ../work/input.nml.template input.nml || exit 1

# these outputs from advance time (with the -g flag) are
# 2 integers: gregorian_day_number seconds
# and since we don't set hours, minutes, or seconds, the second
# number is always 0 and uninteresting for us.
mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
 echo last day is day $end_d

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
 echo first day is day $start_d

# these are a string in the format YYYYMMDDHH
# do them here to prime the loop below which first takes them apart.
prevday=(`echo ${start_year}${mon2}${day2}00 -1d | ../work/advance_time`)
currday=(`echo ${start_year}${mon2}${day2}00   0 | ../work/advance_time`)
nextday=(`echo ${start_year}${mon2}${day2}00 +1d | ../work/advance_time`)

# how many total days are going to be split (for the loop counter)
# (pull out the first of the 2 numbers which are output from advance_time)
let totaldays=${end_d}-${start_d}+1

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  echo "subsetting $d of $totaldays ..."
  #echo $currday $nextday

  # parse out the parts from a string which is YYYYMMDDHH
  # for yesterday(previous), today(current), and tomorrow(next)
  pyear=${prevday:0:4}
  pmonth=${prevday:4:2}
  pday=${prevday:6:2}

  cyear=${currday:0:4}
  cmonth=${currday:4:2}
  cday=${currday:6:2}

  nyear=${nextday:0:4}
  nmonth=${nextday:4:2}
  nday=${nextday:6:2}

  # compute the equivalent gregorian days here.
  g=(`echo ${cyear}${cmonth}${cday}00 -1d -g | ../work/advance_time`)
  greg0=${g[0]}
  g=(`echo ${cyear}${cmonth}${cday}00   0 -g | ../work/advance_time`)
  greg1=${g[0]}
  g=(`echo ${cyear}${cmonth}${cday}00 +1d -g | ../work/advance_time`)
  greg2=${g[0]}

  echo prev $pyear $pmonth $pday which is gregorian $greg0
  echo curr $cyear $cmonth $cday which is gregorian $greg1
  echo next $nyear $nmonth $nday which is gregorian $greg2

  # I have annual files  ...
  # I'll need to revisit this when I wrap over year boundaries ... TJH

  sed -e "s/YYYY/${cyear}/g"    \
      -e "s/MM/${cmonth}/g"     \
      -e "s/PP/${pmonth}/g"     \
      -e "s/DD/${cday}/g"       \
      -e "s/GREG0/${greg0}/g"  \
      -e "s/GREG1/${greg1}/g"  \
      -e "s/GREG2/${greg2}/g"    < ../work/input.nml.template > input.nml

  # make sure output dir exists
  if [[ ! -d ../${cyear}${cmonth} ]] ; then
     mkdir ../${cyear}${cmonth}
  fi

  # do the extract here
  ../work/obs_sequence_tool

  # advance the day; the output is YYYYMMDD00
  prevday=$currday
  currday=$nextday
  nextday=(`echo ${nyear}${nmonth}${nday}00 +1d | ../work/advance_time`)
  #echo $currday $nextday $prevday

  # advance the loop counter
  let d=d+1

done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

