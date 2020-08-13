#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# madis converter - reads in the madis netcdf hourly files
# and output a DART obs_seq file.  this is a one-for-one process; to do the
# next phase where all obs within particular time windows are selected,
# see the 'windowing.csh' script.
#
# this script requires the executable 'advance_time' to be
# built and exist in the current directory, and advance_time
# requires an input.nml namelist file.


# set the following variable to true if you're planning to pass the 
# start & end times for this script in as command line args.  
# set it to false if you're planning to set the times by editing this file.

set command_line_args = false

# set the type, and first and last times.  can roll over day, month and 
# year boundaries.  hours go from 0 to 23; days 1 to 31; months 1 to 12.
# do not add a preceeding 0 to any single digit values.

if ($command_line_args == 'true') then
  if ($#argv != 9) then
     echo usage: $0 type start_year start_month start_day start_hour end_year end_month end_day end_hour
     exit 1
  endif
  set type        = $argv[1]

  set start_year  = $argv[2]
  set start_month = $argv[3]
  set start_day   = $argv[4]
  set start_hour  = $argv[5]
  
  set end_year    = $argv[6]
  set end_month   = $argv[7]
  set end_day     = $argv[8]
  set end_hour    = $argv[9]
else
  # set this to the specific type:  metar, mesonet, acars, marine, rawin
  set type=metar

  set start_year  = 2005
  set start_month = 8
  set start_day   = 26
  set start_hour  = 0
  
  set end_year    = 2005
  set end_month   = 8
  set end_day     = 26
  set end_hour    = 0
endif

# set this to true to make a single output file for each calendar day.
# set it to false to make one output file per input file.
set daily = true

# locations for input madis netcdf files, output dart obs_seq files.
# assumes directory structure under this location with the naming
# convention: .../YYYYMM/DD/${type}
set src_base_dir = /Volumes/joshua3/romine/data/STEP2009/MADIS
set out_dir = /users/romine/step09/MADIS/obs_sequence/conus

set obs_out = obs_seq.${type}

# end of things you should have to set in this script

# convert the start and stop times to gregorian days, so we can
# compute total number of hours including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial hour while we are doing the total hours calc.

# the output of advance time with the -g input is:
#   gregorian_day_number  seconds
# use $var[1] to return just the day number, $var[2] for secs.

set mon2=`printf %02d $end_month`
set day2=`printf %02d $end_day`
set  hr2=`printf %02d $end_hour`
set end_t=(`echo ${end_year}${mon2}${day2}${hr2} 0 -g | ./advance_time`)

set mon2=`printf %02d $start_month`
set day2=`printf %02d $start_day`
set  hr2=`printf %02d $start_hour`
set start_t=(`echo ${start_year}${mon2}${day2}${hr2} 0 -g | ./advance_time`)

# the output of this call is a string YYYYMMDDHH
# see below for help in how to easily parse this up into words
set curhr=`echo ${start_year}${mon2}${day2}${hr2} 0 | ./advance_time`

# how many total hours are going to be processed? (for the loop counter)
if ( $start_t[2] > $end_t[2]) then
   @ end_t[2] += 86400
   @ end_t[1] -= 1
endif
@ totaldays = ( $end_t[1] - $start_t[1] ) 
@ totalsecs = ( $end_t[2] - $start_t[2] ) 
@ totalhrs = ($totaldays * 24) + ($totalsecs / 3600) + 1
echo days, secs = hrs: $totaldays, $totalsecs = $totalhrs

# clean up before starting loop
rm -f $obs_out

# loop over each hour
set h=1
while ( $h <= $totalhrs )

  # parse out the parts from a string which is YYYYMMDDHH
  # use cut with the byte option to pull out columns 1-4, 5-6, 7-8, 9-10
  set  year=`echo $curhr | cut -b1-4`
  set month=`echo $curhr | cut -b5-6`
  set   day=`echo $curhr | cut -b7-8`
  set  hour=`echo $curhr | cut -b9-10`

  # compute the equivalent gregorian day here.
  set g=(`echo ${year}${month}${day}${hour} 0 -g | ./advance_time`)
  set gregday=$g[1]
  set gregsec=$g[2]

  # status/debug - comment in or out as desired.
  echo starting processing for ${year} ${month} ${day} ${hour}
  #echo which is gregorian day: $gregday, $gregsec

  ##################################################################

  # $year, $month, $day, $hour, and $gregday, $gregsec are available
  # to be used as needed.

  ##################################################################

  # make sure output dir exists.  could make it per month, or per day
  # by adding year, month, day to dirname here.
  if (! -d $out_dir) mkdir -p $out_dir

  # append YYYYMM/DD/type to the base directory
  set src_dir = ${src_base_dir}/${year}${month}/${day}/${type}

  # input and output filenames (without paths)
  set infn  = ${year}${month}${day}_${hour}00

  if ($daily == 'true') then
    set outfn = obs_seq_${type}_${year}${month}${day}
  else
    set outfn = obs_seq_${type}_${year}${month}${day}${hour}
  endif

  # if the input is still zipped, unzip it
  if ( -f $src_dir/${infn}.gz ) gunzip $src_dir/${infn}.gz
  if ( ! -f $src_dir/$infn ) echo input filename $src_dir/$infn not found

  ln -sf $src_dir/$infn  ${type}_input.nc

  # set these to T or F as you wish:
  # the rawinsonde converter reads two logicals from standard input to
  # control whether to ouput significant level winds and/or significant
  # level temperatures in addition to the mandatory level info.
  # the satwind converter reads three logicals to control whether to
  # output IR, VIS, and/or WV band winds.
  # the other converters read nothing from stdin.
  if ( ${type} == 'rawin') then
    set instring = 'T T'
  else if ( ${type} == 'satwnd' ) then
    set instring = 'T T T'
  else
    set instring = ''
  endif

  echo $instring | ./convert_madis_${type} >&! out.convert_madis_${type}

  if ($daily == 'true') then
    if ($hour == 23) mv $obs_out $out_dir/$outfn
  else
    mv $obs_out $out_dir/$outfn
  endif

  ##################################################################

  # advance the hour; the output is YYYYMMDDHH
  set curhr=`echo ${year}${month}${day}${hour} +1h | ./advance_time`

  # advance the loop counter
  @ h += 1
 
end

# finish up any partial days and clean up
if ( -f $obs_out) mv $obs_out $out_dir/$outfn

rm -f ${type}_input.nc

exit 0


