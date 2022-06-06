#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# madis data window program -- reads hourly or daily obs_seq files
# by type, and outputs a single obs_seq file with all obs for that
# window.  the output file can then be input directly to the dart
# filter program, or can be input to the wrf preprocessor, which
# does superobs, excluding obs by region, edge processing for special
# wrf boundary condition handling, etc.
# 
# if no preprocessing is needed, this program may not be required.
# the individual madis files can be combined with the 'combine.csh' 
# script and those files can be input to filter.
#
# this script requires the executable 'advance_time' to be
# built and exist in the current directory, and advance_time
# requires an input.nml namelist file.


# set this to true if you're planning to pass the start & end times
# in as command line args.  set it to false if you're planning to set
# the times by editing this file.

set command_line_args = false

# set the first and last analysis times.  can roll over day, month and 
# year boundaries.  hours go from 0 to 23; days 1 to 31; months 1 to 12.

if ($command_line_args == 'true') then
  if ($#argv != 10) then
     echo usage: $0 start_yr start_mon start_dy start_hr start_min end_yr end_mon end_dy end_hr end_min
     exit 1
  endif
  set start_year  = $argv[1]
  set start_month = $argv[2]
  set start_day   = $argv[3]
  set start_hour  = $argv[4]
  set start_min   = $argv[5]
  
  set end_year    = $argv[6]
  set end_month   = $argv[7]
  set end_day     = $argv[8]
  set end_hour    = $argv[9]
  set end_min     = $argv[10]
else
  set start_year  = 2005
  set start_month = 8
  set start_day   = 26
  set start_hour  = 0
  set start_min   = 0
  
  set end_year    = 2005
  set end_month   = 8
  set end_day     = 26
  set end_hour    = 4
  set end_min     = 0
endif

# <ADDME> put more stuff here if you have user settable options

# set frequency of analysis times.  the value should be specified
# in the format for increments to the advance_time program.  
# these should match the expected assimilation windows
# when running filter.
# the output will be a series of obs files with a name based
# on the analysis time computed from this increment.
# valid strings look like +1h30m, +2h, +20m, +1d, etc. 
set time_window = +1h

# set this to true if the input files are daily.
# set it to false if the input files are hourly.
set daily = true

# set this to true if you want the script to fail if any
# input files are missing.  if it is ok for some files to
# not exist at some times, set this to false and the script
# will continue to run even if no output files are created.
set missing_fatal = true

# set this list to all the types that should be processed
# if you reorder or change this list, you must also update
# the win_xxx and file_win_xxx lines below.
set type_list=(metar mesonet acars marine rawin)

# accept observations only within these time windows
# around the analysis time.  these times must be equal to
# or less than half the value of the time_window above, or
# you will start duplicating observations.  these can all
# be the same times, or you can narrow the windows for certain
# data sources.  note that these windows do not have to be
# symmetric around the analysis time.  the +1sec is for
# hourly analysis times where you don't want to replicate
# observations exactly 30 mins away from the hour.
# the order of these MUST match the type_list array above.
set win_before=( -7m  -7m -15m -30m+1s -30m+1s)
set win_after =( +7m  +7m +15m +30m    +30m)

# locations for input madis netcdf files, output dart obs_seq files
set src_base_dir = ../data
set out_dir      = ../data

# depending on your filename conventions, you may have to look
# further below and change the format string for the filename.
# see the NAMES: tag below.

# otherwise, end of things you should have to set in this script

# the input files as received from MADIS are hourly, but
# each type includes a different time window.  these are
# the actual min and max times for each type.  and note
# that these are +/- the hourly times, not necessarily
# the analysis time.  you should NOT change these values
# unless MADIS changes the times that they include in
# each of their hourly files.
# the order of these MUST match the type array above.
set file_win_before=(-15m   0    0  -15m -30m)
set file_win_after =(+45m +59m +59m +45m +30m)

# convert the start and stop times to gregorian days, so we can
# compute total number of hours including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial hour while we are doing the total hours calc.

if ( ! -x ./advance_time ) then
   echo 'FATAL ERROR:'
   echo 'advance_time program not found in current directory.'
   echo 'should be built by the quickbuild.sh script in the'
   echo 'MADIS/work directory. put a copy here and try again.'
   exit 1
endif

if ( ! -x ./obs_sequence_tool ) then
   echo 'FATAL ERROR:'
   echo 'obs_sequence_tool program not found in current directory.'
   echo 'should be built by the quickbuild.sh script in the'
   echo 'MADIS/work directory. put a copy here and try again.'
   exit 1
endif

# the output of advance time with the -g input is:
#   gregorian_day_number  seconds
# use $var[1] to return just the day number, $var[2] for secs.

set mon2=`printf %02d $end_month`
set day2=`printf %02d $end_day`
set  hr2=`printf %02d $end_hour`
set  mn2=`printf %02d $end_min`
set end_t=(`echo ${end_year}${mon2}${day2}${hr2}${mn2} 0 -g | ./advance_time`)

set mon2=`printf %02d $start_month`
set day2=`printf %02d $start_day`
set  hr2=`printf %02d $start_hour`
set  mn2=`printf %02d $start_min`
set start_t=(`echo ${start_year}${mon2}${day2}${hr2}${mn2} 0 -g | ./advance_time`)

# the output of this call is a string YYYYMMDDHHMM
# see below for help in how to easily parse this up into words
set curtime=`echo ${start_year}${mon2}${day2}${hr2}${mn2} 0 | ./advance_time`

# loop until time > end time
while ( 1 )

  # parse out the parts from a string which is YYYYMMDDHHMM
  # use cut with the byte option to pull out columns 1-4, 5-6, 7-8, 9-10, 11-12
  set  year=`echo $curtime | cut -b1-4`
  set month=`echo $curtime | cut -b5-6`
  set   day=`echo $curtime | cut -b7-8`
  set  hour=`echo $curtime | cut -b9-10`
  set   min=`echo $curtime | cut -b11-12`

  # compute the equivalent gregorian day here.
  set g=(`echo $curtime 0 -g | ./advance_time`)
  set gregday=$g[1]
  set gregsec=$g[2]

  # if the current day is beyond the end day, we are done.  break out of loop.
  # if the current day is equal to end day, check the seconds to see if we quit.
  # otherwise, do the loop again.
  if ($gregday > $end_t[1]) break
  if ($gregday == $end_t[1]  &&  $gregsec > $end_t[2]) break

  # status/debug - comment in or out as desired.
  echo ' '
  echo starting processing for ${year} ${month} ${day} ${hour} ${min}
  #echo which is gregorian day: $gregday, $gregsec

  ##################################################################

  # make sure output dir exists.  could make it per month, or per day
  # by adding year, month, day to dirname here.
  if (! -d $out_dir) mkdir $out_dir

  rm -f obsflist obstemp newobslist

  # madis files are per hour; trim off the minutes here
  set curhr = `echo $curtime | cut -c1-10`00

  # loop over types for this analysis time:

  @ t = 1
  while ($t <= $#type_list) 
    set this_type = $type_list[$t]

    echo $this_type window: $win_before[$t] $win_after[$t]
   
    set  obef=(`echo $curtime $win_before[$t]    | ./advance_time` )
    set  oaft=(`echo $curtime  $win_after[$t]    | ./advance_time` )
    set gobef=(`echo $curtime $win_before[$t] -g | ./advance_time` )
    set goaft=(`echo $curtime  $win_after[$t] -g | ./advance_time` )

    set  wbef=(`echo $curhr $file_win_before[$t]    | ./advance_time` )
    set  waft=(`echo $curhr  $file_win_after[$t]    | ./advance_time` )
    set gwbef=(`echo $curhr $file_win_before[$t] -g | ./advance_time` )
    set gwaft=(`echo $curhr  $file_win_after[$t] -g | ./advance_time` )

    set  fbef=(`echo $curhr 0    | ./advance_time` )
    set  faft=(`echo $curhr 0    | ./advance_time` )
    set gfbef=(`echo $curhr 0 -g | ./advance_time` )
    set gfaft=(`echo $curhr 0 -g | ./advance_time` )

    #echo ' obs:' $obef \( $gobef \) to $oaft \( $goaft \)
    #echo 'file:' $wbef \( $gwbef \) to $waft \( $gwaft \)

    rm -f obstemp

    # wxxx are the file window times, oxxx are the obs times.
    # back wbef until it's earlier than obef, and move waft forward
    # until it's later than oaft.  gxxx are the gregorian equivs
    # of the original times.  and finally, fxxx are the hourly file times.

    set fbefdy = `echo $fbef | cut -c1-8`
    set faftdy = `echo $faft | cut -c1-8`

    # NAMES:  the pattern here must match how you choose how you chose
    # to name the obs_seq files converted from each of the original MADIS files.
    # there are a variety of strings available to help you construct a
    # filename: $this_type is metar, marine, etc.  ${fbef},${faft} is 
    # a time string YYYYMMDDHH.  if you choose to name the obs_seq files
    # something other than the full time string, you can parse them up
    # the same as $curtime above.  if you do change the name scheme,
    # there are 4 total lines where this must be changed.
    if ($daily != 'true') then
      echo $src_base_dir/obs_seq_${this_type}_${fbef}  >>! obstemp
      echo $src_base_dir/obs_seq_${this_type}_${faft}  >>! obstemp
    else
      echo $src_base_dir/obs_seq_${this_type}_${fbefdy}  >>! obstemp
      echo $src_base_dir/obs_seq_${this_type}_${faftdy}  >>! obstemp
    endif

    while ( ( $gwbef[1] >  $gobef[1] )  || \
            ( $gwbef[1] == $gobef[1] && $gwbef[2] > $gobef[2] ) )
      set  wbef=(`echo  $wbef  -1h    | ./advance_time` )
      set gwbef=(`echo  $wbef   0  -g | ./advance_time` )

      set  fbef=(`echo  $fbef  -1h    | ./advance_time` )
      set gfbef=(`echo  $fbef   0  -g | ./advance_time` )

      set fbefdy = `echo $fbef | cut -c1-8`

      # NAMES:
      if ($daily != 'true') then
        echo $src_base_dir/obs_seq_${this_type}_${fbef}  >>! obstemp
      else
        echo $src_base_dir/obs_seq_${this_type}_${fbefdy}  >>! obstemp
      endif

      #echo 'back start:' $wbef \( $gwbef \) to $waft \( $gwaft \)
    end

    while ( ( $gwaft[1] <  $goaft[1] ) || \
            ( $gwaft[1] == $goaft[1] && $gwaft[2] < $goaft[2]) )
      set  waft=(`echo  $waft  +1h    | ./advance_time` )
      set gwaft=(`echo  $waft   0  -g | ./advance_time` )

      set  faft=(`echo  $faft  +1h    | ./advance_time` )
      set gfaft=(`echo  $faft   0  -g | ./advance_time` )

      set faftdy = `echo $faft | cut -c1-8`

      # NAMES:
      if ($daily != 'true') then
        echo $src_base_dir/obs_seq_${this_type}_${faft}  >>! obstemp
      else
        echo $src_base_dir/obs_seq_${this_type}_${faftdy}  >>! obstemp
      endif

      #echo 'fore end: ' $wbef \( $gwbef \) to $waft \( $gwaft \)
    end

    # duplicate hours on files are ok, uniq removes them.
    sort obstemp | uniq > obsflist
    echo 'windowing for this type is using input files:'
    cat obsflist

    # sed the input.nml to set start/stop times
    sed -e "s/BDAY/$gobef[1]/" \
        -e "s/BSEC/$gobef[2]/" \
        -e "s/ADAY/$goaft[1]/" \
        -e "s/ASEC/$goaft[2]/" input.nml.template >! input.nml

    # run obs_seq_tool
    ./obs_sequence_tool 
 
    # do something with the output, and clear the input
    # if you want to change the output filename, do it here 
    set out_name = $out_dir/obs_seq_${this_type}_${curtime}

    # move output to name by type, time, and save a copy
    # in a list for a final merge of all types.  make sure
    # the output file was created before adding it to the
    # list of files to be merged.

    if ( -e obs_seq.out ) then
       mv obs_seq.out $out_name
       echo $out_name >>! newobslist
    else
       echo "$out_name not created; merge was unsuccessful."
       echo ' execution error, or input file(s) do not exist.'
       if ($missing_fatal == 'true') then
          echo 'exiting with fatal error'
          exit 1
       endif
    endif

    rm -f obsflist

    @ t ++
  end
  # end of type loop (metar, marine, etc)

  # now cat the type files together
  sort newobslist | uniq >! obsflist
  echo 'merging all types together is using files: '
  cat obsflist

  # sed the input.nml to unset start/stop times
  sed -e "s/BDAY/-1/" \
      -e "s/BSEC/-1/" \
      -e "s/ADAY/-1/" \
      -e "s/ASEC/-1/" input.nml.template > input.nml

  # run obs_seq_tool one more time to stitch these together
  ./obs_sequence_tool
 
  # do something with the output, and clear the input
  # if you don't like the final hourly filename, change it here.
  mv obs_seq.out $out_dir/obs_seq.${curtime}
  rm -f obsflist


  ##################################################################

  # advance to the next analysis time; the output is YYYYMMDDHHMM
  set curtime=`echo ${curtime} $time_window | ./advance_time`

end

exit 0


