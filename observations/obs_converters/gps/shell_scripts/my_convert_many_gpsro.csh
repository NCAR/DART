#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# this script loops over days, calling the GPS convert script
# once per day.  it can roll over month and year boundaries.
#
# it can download data on demand from the CDAAC web site,
# convert it, and delete it, or any combination of these
# three functions, depending on need.  e.g. if the data is
# already downloaded, it can just convert.  it can download
# and convert and not delete until the results are checked
# and then in a second pass just delete, etc, etc.
#
# this script requires the executable 'advance_time' to be
# built and exist in the current directory, and advance_time
# requires a minimal input.nml namelist file (empty &utilities_nml only).
#
# this script constructs the arguments needed for the script that
# is doing all the "real" work:  gpsro_to_obsseq.csh
# see that script for details of what is involved in doing
# the actual conversion.
# 
# -------------------

# if you want to submit this as a batch job:
#PBS -N gps_obs
#PBS -q caldera
#PBS -j oe
#PBS -A project
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:mpiprocs=1

# -------------------

#BSUB -J gps_obs
#BSUB -q caldera
#BSUB -o conv_gps_%J.out
#BSUB -e conv_gps_%J.err
#BSUB -P project
#BSUB -W 15:00
#BSUB -n 1

# -------------------

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR

# -------------------
# -------------------

# start of things you should have to set in this script

# set the first and last days.  can roll over month and year boundaries.
set start_year=2017
set start_month=8
set start_day=1

set end_year=2017
set end_month=8
set end_day=3

#set end_year=2017
#set end_month=9
#set end_day=30


# for each day: 
#  download the data from the web site or not, 
#  convert to daily obs_seq files or not, and 
#  delete the data files after conversion or not.

set do_download = 'yes'
set do_convert  = 'no'
set do_delete   = 'no'


# set the list of satellite data to convert.
# - in the comments below, 'realtime' is usually data up to the
#   current date with less quality control.
# - 'reprocessed' generally has the highest quality
# - dates are YYYY.DDD  where DDD is day number in the year.
# - only select one of reprocessed or realtime for a particular
#   satellite or you will get duplicate observations if they
#   have overlapping time periods.
#
# WARNING: the available obs are updated frequently as the
# data is reprocessed and new obs arrive.
# check the CDAAC web site for the currently available days.  
#
# log in here: http://cdaac-www.cosmic.ucar.edu/cdaac/login/
# go to 'data center' then 'data access' for the current
# table of satellites and days that have data files.
#
# these dates were current as of May 2016:
#
# name_to_use    date range           what instrument?
# -----------    ----------           ---------------
# champ2016    2001.138 - 2008.279    CHAMP
# cnofs        2010.060 - 2011.365    Air Force C/NOFS
# cnofsrt      2012.001 - 2015.193      "    (realtime)
# cosmic2013   2006.112 - 2014.120    COSMIC, reprocessed
# cosmic       2014.121 - 2015.364      "    (not reprocessed yet)
# cosmicrt     2014.181 - 2016.123      "    (realtime)
# gpsmet       1995.111 - 1997.047    ?
# gpsmetas     1995.237 - 1997.016    ?
# grace        2007.059 - 2015.364    Grace-A
# kompsat5rt   2015.305 - 2016.123    ?
# metopa2016   2007.274 - 2015.365    Metop-A/GRAS, reprocessed 2016
# metopb       2013.032 - 2015.059    Metop-B/GRAS 
# sacc         2006.068 - 2011.215    Argentinan SAC-C
# saccrt       2011.329 - 2013.226      "    (realtime)
# tsx          2008.041 - 2015.333    German TerraSAR-X
#

# which satellites to include:

rm -fr satlist
echo cosmic      >> satlist
echo grace       >> satlist
echo metopa      >> satlist
echo metopb      >> satlist
echo tsx         >> satlist


# make the converter script happy by making these files consistent
# and link to a working advance_time

cp satlist ../work
ln -sf ../work/advance_time
if ( ! -f input.nml ) then
  echo \&utilities_nml > input.nml
  echo / >> input.nml
endif


# where to download the data and do the conversions, relative to
# this shell_scripts directory.  the script below will add YYYYMM
# to the end of this string.

#set datadir = ../gpsro
#set datadir = /glade/p/cisl/dares/Observations/GPS/staged
set datadir = /glade/scratch/nancy/staged

# end of things you should have to set in this script

# -------------------
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

  ./my_gpsro_to_obsseq.csh ${year}${month}${day} $datadir/${year}${month} \
                         $do_download $do_convert $do_delete ./satlist


  # advance the day; the output is YYYYMMDD00
  set curday=`echo ${year}${month}${day}00 +1d | ./advance_time`

  # advance the loop counter
  @ d++
 
end

exit 0


