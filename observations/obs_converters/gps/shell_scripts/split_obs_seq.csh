#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#

# Example of a shell script to split a series of obs_sequence files
# into a series of multiple smaller files.  It requires the 'advance_time'
# executable to do the time/calendar computations, and it will almost
# certainly require some customization for use with any input other
# than our NCEP daily obs_seq files.  But it should be a good starting
# point for additional scripting. Thanks to Josh Hacker for contributing it.
#
# Split a series of obs_sequence files into multiple smaller files.
# This depends on knowing how the original files are organized, e.g. are
# the observations one file per day from 0Z to 0Z, or they are one
# file per day but 3Z to 3Z as our NCEP obs files are organized (so
# 6 hour assimilation windows centered on 6,12,18,24Z can be run 
# from a single input file).   The cycle interval below should be
# an exact multiple of the length of the input files.

# set start and end YYYYMMDDHHMM here. 
# this is the midpoint of each time window.
set starth = 201303140010
set   endh = 201303140100

# set time window for each output file.  the syntax of these
# lines must match the valid input formats for the advance_time program
set cycle_interval = 10m      # total interval duration
set first_half     = -5m+1s  # -half interval + 1 second
set last_half      = +5m     # +half interval


# set once and should be able to leave as-is
set input_dir    = ../daily
set output_dir   = ../10min
set advance_exec = ./advance_time
set nml_template = ./input.nml.template
set COPY         = "/bin/cp -f"
set MOVE         = "/bin/mv -f"


# loop from start to end time.
set dtg = $starth
while ( $dtg <= $endh )
  
  echo processing window centered on time $dtg 
  set start_window=(`echo $dtg ${first_half} -g | $advance_exec`)
  set start_day = $start_window[1]
  set start_sec = $start_window[2]
  set end_window=(`echo $dtg ${last_half} -g | $advance_exec`)
  set end_day = $end_window[1]
  set end_sec = $end_window[2]

  # cut off hours
  set ymd = `echo $dtg | cut -c 1-8`

  # get year, month, date for use below
  set y = `echo $dtg | cut -c 1-4`
  set m = `echo $dtg | cut -c 5-6`
  set d = `echo $dtg | cut -c 7-8`

  # set the format of the input and output files
  set input_file  = \'$input_dir/obs_seq.ion.${y}-${m}-${d}\'
  set output_file = \'$output_dir/obs_seq.ion.${dtg}\'

  # decrement string used for input filename by a day if 
  # interval center is 00 h 
  # NOTE: whether you need to do this or not depends on how
  # the rollover of days is managed, also how files are named
  # (0Z vs 24Z).
  set hh = `echo $dtg | cut -c 9-10`
  set mm = `echo $dtg | cut -c 11-12`
  if ( $hh == "00" && $mm == "00") then
    set dtg2 = `echo $dtg -24h | $advance_exec`
    set y2 = `echo $dtg2 | cut -c 1-4`
    set m2 = `echo $dtg2 | cut -c 5-6`
    set d2 = `echo $dtg2 | cut -c 7-8`
    set input_file = \'$input_dir/obs_seq.ion.${y2}-${m2}-${d2}\',$input_file
  endif
    
  echo input will come from file $input_file
  echo output will go into file $output_file

  # change the values in the template file and overwrite the
  # old input.nml - this controls what the execution of the
  # obs_sequence_tool will do.
  sed -e "s/^.*first_obs_days *=.*/first_obs_days\ =\ $start_day\,/" \
      -e "s/^.*first_obs_seconds *=.*/first_obs_seconds\ =\ $start_sec\,/" \
      -e "s/^.*last_obs_days *=.*/last_obs_days\ =\ $end_day\,/" \
      -e "s/^.*last_obs_seconds *=.*/last_obs_seconds\ =\ $end_sec\,/" \
      -e "s;^.*filename_seq *=.*;filename_seq=\ ${input_file},;" \
      -e "s;^.*filename_out *=.*;filename_out=\ ${output_file},;" \
      $nml_template > input.nml      

  # do the splitting here
  ./obs_sequence_tool

  # advance to next time
  set dtg = `echo $dtg ${cycle_interval} | $advance_exec`

end

echo 'Finished'

exit 0


