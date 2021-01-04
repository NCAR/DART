#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# merge the files into "daily" files which start at 03:01Z
# and end at 03:00Z the following day.  (the name of the file
# is the first day.)

# set the first and last days to be merged.  can roll over
# month and year boundaries now!   note that for the end day,
# you need at least the first 40ish files from the following day
# for the merge to have the right data (from 0Z to 3Z) available.

let start_year=2006
let start_month=10
let start_day=1

let end_year=2007
let end_month=1
let end_day=31

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
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ./advance_time`)

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ./advance_time`)

# these are a string in the format YYYYMMDDHH
# do them here to prime the loop below which first takes them apart.
curday=(`echo ${start_year}${mon2}${day2}00 0 | ./advance_time`)
nextday=(`echo ${start_year}${mon2}${day2}00 +1d | ./advance_time`)

# how many total days are going to be merged (for the loop counter)
# (pull out the first of the 2 numbers which are output from advance_time)
let totaldays=${end_d[0]}-${start_d[0]}+1

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  # parse out the parts from a string which is YYYYMMDDHH
  # both for the current day and tomorrow
  cyear=${curday:0:4}
  cmonth=${curday:4:2}
  cday=${curday:6:2}
  nyear=${nextday:0:4}
  nmonth=${nextday:4:2}
  nday=${nextday:6:2}

  # compute the equivalent gregorian days here.
  g=(`echo ${cyear}${cmonth}${cday}00 0 -g | ./advance_time`)
  greg1=${g[0]}
  let greg2=greg1+1

  echo starting AIRS obs merge ${cyear}${cmonth}${cday}
  echo gregorian: $greg

  # all of todays data plus the first 40 of tomorrows
  ls AIRS.${cyear}.${cmonth}.${cday}.*.out > olist
  ls AIRS.${nyear}.${nmonth}.${nday}.0[0123]?.out >> olist

  sed -e "s/YYYY/${cyear}/g"    \
      -e "s/MM/${cmonth}/g"     \
      -e "s/DD/${cday}/g"       \
      -e "s/GREG1/${greg1}/g"  \
      -e "s/GREG2/${greg2}/g"    < ./input.nml.template > input.nml

  # do the merge here
  ./obs_sequence_tool
  

  # advance the day; the output is YYYYMMDD00
  curday=nextday
  nextday=(`echo ${year}${month}${day}00 +1d | ./advance_time`)

  # advance the loop counter
  let d=d+1
 
done

exit 0

