#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Parallel Version - fork up to 32 processes to run at once.
# Invoke this script with 4 args: year, month, start/end day
#
#--------------------------------------------------------------
#
#  This script is used to generate daily (3:01Z to 3:00Z of next day) decoded 
#  NCEP reanalysis PREPBUFR text/ascii data.
#  It should be run only for days within a single month.
#
#  This script requires 2 arguments - start day number and stop day number
#  (it will convert up to and including the stop day)
#
#--------------------------------------------------------------
# USER SET PARAMETERS

# if daily is 'yes', 4 6-hour files will be processed and a single, 1-day
# output file will be created.  if daily is set to anything else, then
# each 6-hour input file will be converted to a single 6-hour output file.
# this script still processes a day's worth of files at a time even if
# daily is 'no' - it just makes 4 individual files per day.

set    daily = yes

# if daily is 'no' and zeroZ is 'yes', input files at 0Z will be translated 
# into output files also marked 0Z.  otherwise, they will be named with the 
# previous day number and 24Z (chose here based on what script will be 
# processing these files next.  the 'create_real_obs' script assumes 
# filenames with the pattern 6Z,12Z,18Z,24Z, so 'no' is right for it.)
# this variable is ignored completely if daily is 'yes'.

set zeroZ = no

# if convert is 'yes', then the big-endian BUFR files will be converted
# to little-endian files before processing. this is needed if you are running
# on a machine that uses Intel chips (e.g. linux clusters, altix, pcs, etc).
# it is not needed for ibm power systems.  any value other than 'yes' will
# skip the convert step.

set  convert = no 

# if block is 'yes', then the cword program will be run to convert an
# unblocked file into a blocked one.  this is not required for recent
# prepbufr files, but older ones may require it.

set block = yes

# starting year, month, day, and ending day.  this script does not allow
# you to do more than a single month at a time, but does handle the last
# day of the month, leap day in feb, and the last day of the year correctly.
# this version of the conversion tool takes up to 3 hours of observations
# from the day *following* the end day, so you must have at least the 6Z
# file from one day beyond the last day for this to finish ok (to get obs
# at exactly 3Z from the next file).

set year     = $argv[1]
set month    = $argv[2]
set beginday = $argv[3]
set endday   = $argv[4]

# directory where the BUFR files are located.  the script assumes the
# files will be located in subdirectories by month, with the names following
# the pattern YYYYMM, and then inside the subdirectories, the files are
# named by the pattern 'prepqmYYMMDDHH'.  for example, if the dir below
# is the default ../data, then the 6Z file for dec 27th, 2010 would be:
#  ../data/201012/prepqm10122706
# but the conventions for names of prepqm files have changed over the years,
# so if the prepqm files do *not* follow this pattern, you will have to edit
# the BUFR_in variable in the script below to match the filenames you have.
# the setting of BUFR_out matches what the 'create_real_obs' script expects
# to have as input, but if you want to generate a different set of names
# you can change it below as well.
# there are several shell variables in the loop you can use to construct
# alternate names:
# 'year' is 4 digits; 'yy' is 2.
# 'mm', 'dd', and 'hh' are 0 padded so they are always 2 digits.
# 'oyear', 'omm', 'odd', 'ohh' are the original date, if the day, month, 
# and/or year have rolled over.

set BUFR_dir = ../data
#set BUFR_dir = /ptmp/dart/Obs_sets/NCEP_bufr

# directory where DART ncep observation programs are located, relative
# to one directory down from the current directory.  (the parallel
# version starts each conversion in a new subdir from here.)

# one deeper than the serial version.
set DART_exec_dir = ../../exe

# END USER SET PARAMETERS
#--------------------------------------------------------------

set days_in_mo = (31 28 31 30 31 30 31 31 30 31 30 31)
# leap years: year 2000 requires that you do even the centuries right
if (($year %   4) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1
if (($year % 100) == 0) @ days_in_mo[2] = $days_in_mo[2] - 1
if (($year % 400) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1



# save original year in case we roll over the year
set oyear = $year
@ oyy = $year % 100

# Loop over days
set day = $beginday
set last = $endday
while ( $day <= $last )

   # clear any old intermediate (text) files
   rm -f temp_obs prepqm.in prepqm.out 

   # save a copy of the original day and month, in case we roll over below,
   # in both single/double digit format and guarenteed 2 digit format.
   set oday = $day
   set odd  = $day
   if ($odd < 10) set odd = 0$odd
   set omonth = $month
   set omm    = $month
   if ($omm < 10) set omm = 0$omm

   # convert 1 days worth (data from 3:01Z to 3:00Z of the next day) 
   # of BUFR files into a single intermediate file if 'daily' is set to true; 
   # into 4 6-hour files otherwise.
   set h = 0
   set next_day = not
   # daily files need to pull out observations at exactly 3Z from the next day.
   # non-daily files should process up to hour 24, one for one with the inputs.
   while ( (($h < 30) && ($daily == yes)) || (($h < 24) && ($daily == no)) )
      @ h  = $h + 6
      @ hh = $h % 24
      @ dd = $day + ($h / 24)
      @ yy = $year % 100
      set mm = $month

      # special handling for the end of the day, month, year
      if ($hh == 0 || ($hh > 0 && $next_day == yes) ) then
         set next_day = yes
         if ($dd > $days_in_mo[$month]) then
            # next month
            # signal that this is the last day to do
            set last = 0
            set dd = 1
            @ mm++
            if ($mm > 12) then
               if ($mm > 12 && $hh == 0) then
                 @ year ++
               endif
               # next year
               set mm = 1
               @ yy = $year % 100
            endif
         endif
      endif

      # format the date for filename construction, guarenteed 2 digits
      if ($yy < 10) set yy = 0$yy
      if ($mm < 10) set mm = 0$mm
      if ($dd < 10) set dd = 0$dd
      if ($hh < 10) set hh = 0$hh
      set ohh = $h   # hour not modulo 24
      if ($h < 10)  set ohh = 0$h


      # the prepqm input files.  match general naming pattern with the
      # data time encoded in the filename.  if the pattern of the filename
      # is different (2 digit year vs 4, extra fixed text in the name, etc)
      # fix the BUFR_in line below to match what you have.  if the file is
      # gzipped, you can leave it and this program will unzip it before
      # processing it.
      set BUFR_in = ${BUFR_dir}/${year}${mm}/prepqm${yy}${mm}${dd}${hh}
      #set BUFR_in = ${BUFR_dir}/${year}${mm}/prepqm${yy}${mm}${dd}${hh}.no_ship_id

      if ( -e ${BUFR_in} ) then
         echo "copying ${BUFR_in} into prepqm.in"
         rm -f prepqm.in
         cp -f ${BUFR_in} prepqm.in
      else if ( -e ${BUFR_in}.gz ) then
         echo "unzipping ${BUFR_in}.gz into prepqm.in"
         rm -f prepqm.in
         gunzip -c -f ${BUFR_in}.gz >! prepqm.in
      else
         echo "MISSING INPUT FILE: cannot find either"
         echo ${BUFR_in}
         echo   or 
         echo ${BUFR_in}.gz
         echo "Script will abort now."
         exit -1
      endif

      # blocking
      if ($block == 'yes') then
         echo "blocking prepqm.in"
         mv -f prepqm.in prepqm.unblocked
         echo 'block' >! in
         echo 'prepqm.unblocked' >> in
         echo 'prepqm.blocked' >> in
         ${DART_exec_dir}/cword.x < in
         mv -f prepqm.blocked prepqm.in
         rm -f prepqm.unblocked in
      endif

      # byte swapping
      if ($convert == 'yes') then
         echo "byteswapping bigendian to littleendian prepqm.in"
         mv -f prepqm.in prepqm.bigendian
         ${DART_exec_dir}/grabbufr.x prepqm.bigendian prepqm.littleendian
         mv -f prepqm.littleendian prepqm.in
         rm -f prepqm.bigendian
      endif

      if ($h == 30) then
         # get any obs exactly at 3Z from the 6Z file of the next day using a
         # modified prepbufr program, since 6 hour assimilation windows 
         # centered on 6Z, 12Z, etc would be 03:01Z-09:00Z, 09:01Z-15:00Z, etc.
         # obs exactly at 03Z need to be part of the same file which spans 
         # 21:01Z-03:00Z.   both of these prepbufr programs also add 24h to 
         # any time after midnight, so the hours in the output files run
         # from 3.016 hours to 27.000 hours.  (the ascii intermediate files
         # do not contain day numbers, but probably should and then the hours
         # can run from a normal 0Z to 23:59Z.)
         ${DART_exec_dir}/prepbufr_03Z.x
      else
         ${DART_exec_dir}/prepbufr.x
      endif

      if ($daily == 'yes') then
         cat prepqm.out >>! temp_obs
         rm -f prepqm.out
      else
         if ($zeroZ == 'yes') then
            # if 0Z, output named with current day and 0Z
            set BUFR_out = ${BUFR_dir}/${year}${mm}/temp_obs.${year}${mm}${dd}${hh}
         else
            # if 0Z, output named with previous day and 24Z
            set BUFR_out = ${BUFR_dir}/${oyear}${omm}/temp_obs.${oyear}${omm}${odd}${ohh}
         endif
         echo "moving output to ${BUFR_out}"
         mv -fv prepqm.out ${BUFR_out}
      endif
   end

   if ($daily == 'yes') then
      # use the original dates without rollover
      set BUFR_out = ${BUFR_dir}/${oyear}${omm}/temp_obs.${oyear}${omm}${odd}
      echo "moving output to ${BUFR_out}"
      mv -fv temp_obs ${BUFR_out}
   endif

   rm -f prepqm.in

   @ day++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

