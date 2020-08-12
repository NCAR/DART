#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# this is a template for a shell script that can loop
# over multiple hours and can roll over day, month, and
# even year boundaries.  see the section inside the loop
# for the 'your code goes here' part. look for the string ADDME.
# this script requires the executable 'advance_time' to be
# built and exist in the current directory, and advance_time
# requires an input.nml namelist file.

# set this to true if you're planning to pass the start & end times
# in as command line args.  set it to false if you're planning to set
# the times by editing this file.

set command_line_args = false

# set the first and last times.  can roll over day, month and year boundaries.
# hours go from 0 to 23; days from 1 to 31, months from 1 to 12.

if ($command_line_args == 'true') then
  if ($#argv != 8) then
     echo usage: $0 start_year start_month start_day start_hour end_year end_month end_day end_hour
     exit 1
  endif
  set start_year=$argv[1]
  set start_month=$argv[2]
  set start_day=$argv[3]
  set start_hour=$argv[4]

  set end_year=$argv[5]
  set end_month=$argv[6]
  set end_day=$argv[7]
  set end_hour=$argv[8]
else
  set start_year=2001
  set start_month=06
  set start_day=01
  set start_hour=0

  set end_year=2001
  set end_month=12
  set end_day=31
  set end_hour=0
endif

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

set mon2=`printf $end_month`
set day2=`printf  $end_day`
set  hr2=`printf  $end_hour`
set end_t=(`echo ${end_year}${mon2}${day2}${hr2} 0 -g | ../work/advance_time`)

set mon2=`printf $start_month`
set day2=`printf  $start_day`
set  hr2=`printf  $start_hour`
set start_t=(`echo ${start_year}${mon2}${day2}${hr2} 0 -g | ../work/advance_time`)

# the output of this call is a string YYYYMMDDHH
# see below for help in how to easily parse this up into words
set curhr=`echo ${start_year}${mon2}${day2}${hr2} 0 | ../work/advance_time`

# how many total hours are going to be processed (for the loop counter)
# note that the parens below are necessary; otherwise the computation
# does total = end - (start+1), or total = end - start - 1, which is
# not how elementary math is supposed to work.
if ( $start_t[2] > $end_t[2]) then
   @ end_t[2] += 86400
   @ end_t[1] -= 1
endif
@ totaldays = ( $end_t[1] - $start_t[1] + 1 )

# loop over each day

set obsindir = "/glade/p/work/yfzhang/observations/modis-tsfc"
set obsoutdir = "/glade/p/work/yfzhang/observations/modis-tsfc/obs_seqs"

set d = 1
while ( $d <= $totaldays )

   # parse out the parts from a string which is YYYYMMDDHH
   # use cut with the byte option to pull out columns 1-4, 5-6, 7-8, 9-10
   set  year=`echo $curhr | cut -b1-4`
   set month=`echo $curhr | cut -b5-6`
   set   day=`echo $curhr | cut -b7-8`
   set  hour=`echo $curhr | cut -b9-10`

   # compute the equivalent gregorian day here.
   set g=(`echo ${year}${month}${day}${hour} 0 -g | ../work/advance_time`)
   set gregday=$g[1]
   set gregsec=$g[2]

   # compute the day of the year here.
   set julian=(`echo ${year}${month}${day}${hour} 0 -j | ../work/advance_time`)
   set doy=$julian[2]

   # status/debug - comment in or out as desired.
   echo starting processing for ${year} ${month} ${day} ${hour}
   echo which is gregorian day: $gregday, $gregsec
   echo which is   day-of-year: $doy

   # <ADDME> your code goes here.
   # use $year, $month, $day, $hour, and $gregday, $gregsec as needed.

   # change the namelist to convert the MODIS file to 
   # a uniquely-named observation sequence file.

   @ doy_add = $doy - 1

   set doy_3d = `printf %03d $doy`
   set filein  = "MOD29E1D.A${year}${doy_3d}.IST.NH.nc"
   set fileout = "obs_seq.${year}-${month}-${day}-00000"
   echo $obsindir/$filein
   echo $obsoutdir/$fileout

   sed -e "/year/ c\   year = $year"  \
       -e "/doy/ c\   doy = $doy_add"  \
       -e "/seaice_input_file/ c\  seaice_input_file = '${obsindir}/$filein'" \
       -e "/obs_out_file/ c\  obs_out_file = '${obsoutdir}/$fileout'" input.nml > temp
   mv temp input.nml

   ../work/modis_ist_to_obs_netcdf

   # advance the time counter; the output is YYYYMMDDHH
   set curhr=`echo ${year}${month}${day}${hour} +1d | ../work/advance_time`

   # advance the loop counter
   @ d += 1

end

exit 0


#%# # example of using sed and lists of obs files to automate
#%# # calling the obs_sequence_tool to split or combine obs_seq files:
#%#
#%# # put a list of filenames into 'obstemp' somehow
#%#
#%# # remove duplicate filenames
#%# sort obstemp | uniq > infilelist
#%# echo 'using input files:'
#%# cat infilelist
#%#
#%# # if the start and stop times are in gregorian format,
#%# # in $start and $stop, use sed to set the input.nml
#%# sed -e "s/BDAY/$start[1]/" \
#%#     -e "s/BSEC/$start[2]/" \
#%#     -e "s/ASEC/$stop[2]/"   \
#%#     -e "s/ASEC/$stop[2]/"    input.nml.template >! input.nml
#%#
#%# # run obs_seq_tool
#%# ./obs_sequence_tool
#%#
#%# # move the output someplace
#%# mv obs_seq.combined obs_seq.$curhr
#%#


