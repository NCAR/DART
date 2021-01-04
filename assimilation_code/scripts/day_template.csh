#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

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
set start_year=2006
set start_month=10
set start_day=31

set end_year=2006
set end_month=11
set end_day=3

# <ADDME> put more stuff here if you have user settable options

# end of things you should have to set in this script

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
# not how elementary math is supposed to work.
@ totaldays = ( $end_d[1] - $start_d[1] ) + 1

# loop over each day
set d=1
while ( $d <= $totaldays )

  # parse out the parts from a string which is YYYYMMDDHH
  # use cut with the byte option to pull out columns 1-4, 5-6, and 7-8
  set  year=`echo $curday | cut -b1-4`
  set month=`echo $curday | cut -b5-6`
  set   day=`echo $curday | cut -b7-8`

  # compute the equivalent gregorian day here.
  set g=(`echo ${year}${month}${day}00 0 -g | ./advance_time`)
  set greg=$g[1]

  # status/debug - comment in or out as desired.
  echo starting processing for ${year} ${month} ${day}
  echo which is gregorian day: $greg


  # <ADDME> your code goes here.  
  # use $year, $month, $day, and $greg as needed.


  # advance the day; the output is YYYYMMDD00
  set curday=`echo ${year}${month}${day}00 +1d | ./advance_time`

  # advance the loop counter
  @ d += 1
 
end

exit 0


