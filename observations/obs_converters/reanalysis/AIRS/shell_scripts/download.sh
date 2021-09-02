#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: download.sh 10982 2017-02-01 23:43:10Z thoar@ucar.edu $
#
# OBSOLETE: the AIRS obs are no longer on the mass store.

# set the first and last days. can roll over
# month and year boundaries now!
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
mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ./advance_time`)

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ./advance_time`)

curday=(`echo ${start_year}${mon2}${day2}00 0 | ./advance_time`)

# how many total days are going to be converted (for the loop counter)
let totaldays=${end_d[0]}-${start_d[0]}+1

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  # parse out the parts from a string which is YYYYMMDDHH
  year=${curday:0:4}
  month=${curday:4:2}
  day=${curday:6:2}


  echo should be getting ${year}${month}${day}.tar from somewhere
  # PUT SOME OTHER COMMAND HERE IF YOU HAVE A SOURCE FOR
  # THESE OBS.


  # advance the day; the output is YYYYMMDD00
  curday=(`echo ${year}${month}${day}00 +1d | ./advance_time`)

  # advance the loop counter
  let d=d+1
 
done

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/rma_trunk/observations/obs_converters/AIRS/shell_scripts/download.sh $
# $Revision: 10982 $
# $Date: 2017-02-01 16:43:10 -0700 (Wed, 01 Feb 2017) $

