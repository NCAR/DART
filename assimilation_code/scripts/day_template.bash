#!/bin/bash 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# this is a template for a shell script that can loop
# over multiple days and roll over month boundaries
# or even year boundaries.  see the section inside the loop
# for the 'your code goes here' part. look for the string ADDME.
# this script requires the executable 'advance_time' to be
# built and exist in the current directory, and advance_time
# requires an input.nml namelist file.

# it computes the gregorian day if useful, and makes sure
# the days and months are always 2 digits long.

# set the first and last days.  can roll over month and year boundaries.
let start_year=2006
let start_month=1
let start_day=31

let end_year=2006
let end_month=2
let end_day=3

# <ADDME> put more stuff here if you have user settable options

# end of things you should have to set in this script

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# the output of advance time with the -g input is:
#   gregorian_day_number  seconds
# use ${var[0]} to return just the day number

mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ./advance_time`)

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ./advance_time`)

# the output of this call is a string YYYYMMDDHH
# see below for help in how to easily parse this up into words
curday=`echo ${start_year}${mon2}${day2}00 0 | ./advance_time`

# how many total days are going to be processed (for the loop counter)
let totaldays=${end_d[0]}-${start_d[0]}
let totaldays=$totaldays+1

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  # parse out the parts from a string which is YYYYMMDDHH
  # use cut with the byte option to pull out columns 1-4, 5-6, and 7-8
   year=`echo $curday | cut -b1-4`
  month=`echo $curday | cut -b5-6`
    day=`echo $curday | cut -b7-8`

  # compute the equivalent gregorian day here.
  g=(`echo ${year}${month}${day}00 0 -g | ./advance_time`)
  greg=${g[0]}

  # status/debug - comment in or out as desired.
  echo starting processing for ${year} ${month} ${day}
  echo which is gregorian day: $greg


  # <ADDME> your code goes here.
  # use $year, $month, $day, and $greg as needed


  # advance the day; the output is YYYYMMDD00
  curday=`echo ${year}${month}${day}00 +1d | ./advance_time`

  # advance the loop counter
  let d=d+1
 
done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

