#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: day_template.csh 4203 2009-12-17 22:26:35Z thoar $
#
# this script loops over days, calling the GPS convert script
# once per day.  it can roll over month and year boundaries.
#
# this script requires the executable 'advance_time' to be
# built and exist in the current directory, and advance_time
# requires a minimal input.nml namelist file (empty utilities_nml only).

# set the first and last days.  can roll over month and year boundaries.
set start_year=2007
set start_month=12
set start_day=30

set end_year=2008
set end_month=1
set end_day=10

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
  #echo which is gregorian day: $greg

  # use $year, $month, $day, and $greg as needed.
  # month, day have leading 0s if needed so they are always 2 digits

  # THE WORK HAPPENS HERE
  # call the convert script.  in this case the 3 yes's are to
  # download, convert, delete in one go - each day at a time.

  ./cosmic_to_obsseq.csh ${year}${month}${day} ../cosmic yes yes yes


  # advance the day; the output is YYYYMMDD00
  set curday=`echo ${year}${month}${day}00 +1d | ./advance_time`

  # advance the loop counter
  @ d += 1
 
end

exit 0

# <next few lines under version control, do not edit>
# $URL: https://subversion.ucar.edu/DAReS/DART/trunk/shell_scripts/day_template.csh $
# $Revision: 4203 $
# $Date: 2009-12-17 15:26:35 -0700 (Thu, 17 Dec 2009) $

