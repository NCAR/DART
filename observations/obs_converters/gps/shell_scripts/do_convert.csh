#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# this script loops over days, calling the GPS convert script
# once per day.  it can roll over month and year boundaries.
#
# this script requires the executable 'advance_time' to be
# built and exist in the current directory, and advance_time
# requires a minimal input.nml namelist file (empty utilities_nml only).

# -------------------

# set the first and last days.  can roll over month and year boundaries.
set start_year=2009
set start_month=9
set start_day=1

set end_year=2009
set end_month=9
set end_day=3


# for each day: download the data or not, convert to daily obs_seq files 
# or not, and delete the data files after conversion or not.
set do_download = 'yes'
set do_convert  = 'yes'
set do_delete   = 'no'


# set the list of satellite data to convert.
# - in the comments below, 'now*' is the current date minus 3-4 months 
#   since the reprocessed datasets lag the realtime data.  'realtime' 
#   is data from up to today's date but with less quality control.
# - dates below are YYYY.DDD  where DDD is day number in the year.
# - only select one of reprocessed or realtime for a particular
#   satellite or you will get duplicate observations.

rm -fr satlist
echo cosmic      >>! satlist  # all 6 COSMIC : 2006.194 - now*
## echo cosmicrt >>! satlist  # COSMIC : realtime
echo sacc        >>! satlist  # Argentinan SAC-C : 2006.068 - now*
## echo saccrt   >>! satlist  # SAC-C : realtime
echo ncofs       >>! satlist  # new Air Force C/NOFS : 2010.335 - now*
## echo ncofsrt  >>! satlist  # C/NOFS : realtime
echo grace       >>! satlist  # Grace-A : 2007.059 - now*
echo tsx         >>! satlist  # German TerraSAR-X : 2008.041 - now*
echo metopa      >>! satlist  # Metop-A/GRAS : 2008.061 - now*
echo champ       >>! satlist  # CHAMP : 2001.139 - 2008.274


# where to download the data and do the conversions, relative to
# this shell_scripts directory.
set datadir = ../gpsro

# end of things you should have to set in this script

# -------------------

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

  # THE WORK HAPPENS HERE:  call the convert script for each day.

  ./gpsro_to_obsseq.csh ${year}${month}${day} $datadir \
                         $do_download $do_convert $do_delete ./satlist


  # advance the day; the output is YYYYMMDD00
  set curday=`echo ${year}${month}${day}00 +1d | ./advance_time`

  # advance the loop counter
  @ d += 1
 
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

