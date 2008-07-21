#!/bin/csh 
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

#--------------------------------------------------------------
# DESCRIPTION:
#
#  This script is used to generate daily (3:01Z to 3:00Z of next day) decoded 
#  NCEP reanalysis PREPBUFR text/ascii data.
#
#--------------------------------------------------------------

# you may need additional options to submit this job to LSF;
# e.g. -P project charge code, different queue name, etc.

#BSUB -o prepbufr.out
#BSUB -e prepbufr.err
#BSUB -J prepbufr
#BSUB -q regular
#BSUB -W 1:00
#BSUB -P XXXXXXXX
#BSUB -n 1


#
#--------------------------------------------------------------
#
#  This script is used to generate daily (3:01Z to 3:00Z of next day) decoded 
#  NCEP reanalysis PREPBUFR text/ascii data.
#  It should be run only for days within a single month.
#
#--------------------------------------------------------------
# USER SET PARAMETERS

# if daily is 'yes', 4 6-hour files will be processed and a single, 1-day
# output file will be created.  if daily is set to anything else, then
# each 6-hour input file will be converted to a single 6-hour output file.

set    daily = yes

# if convert is 'yes', then the big-endian BUFR files will be converted
# to little-endian files before processing. this is needed if you are running
# on a machine that uses Intel chips (e.g. linux clusters, altix, pcs, etc).
# it is not needed for ibm power system.  any value other than 'yes' will
# skip the convert step.

set  convert = no 

# starting year, month, day, and ending day.  this script does not allow
# you to do more than a single month at a time, but does handle the last
# day of the month, leap day in feb, and the last day of the year correctly.
# this version of the conversion tool takes up to 3 hours of observations
# from the day *following* the end day, so you must have at least the 6Z
# file from one day beyond the last day for this to finish ok.

set year     = 1989
set month    = 1
set beginday = 1
set endday   = 7

# directory where the BUFR files are located.  the script assumes the
# files will be located in subdirectories by month, with the names following
# the pattern YYYYMM, and then inside the subdirectories, the files are
# named by the pattern 'prepqmYYYYMMDDHH'.  for example, if the dir below
# is the default ../data, then the 6Z file for jan 1st, 1989 would be:
#  ../data/198901/prepqm1989010106
# if the prepqm files do *not* follow this pattern, you may have to edit
# the BUFR_file and BUFR_out variables below inside the loop of the script.
# year is 4 digits; month, day, hour are 1 or 2 digits.  yy, mm, dd, hh are
# all truncated or padded to be exactly 2 digits.

set BUFR_dir = ../data

# END USER SET PARAMETERS
#--------------------------------------------------------------

set days_in_mo = (31 28 31 30 31 30 31 31 30 31 30 31)
# leap years: year 2000 requires that you do even the centuries right
if (($year %   4) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1
if (($year % 100) == 0) @ days_in_mo[2] = $days_in_mo[2] - 1
if (($year % 400) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1

rm -f prepqm.out 


# save original year in case we roll over the year
set oyear = $year
@ oyy = $year % 100

# Loop over days
set day = $beginday
set last = $endday
while ( $day <= $last )
   echo '-------------------------------------- '

   # clear any old intermediate (text) BUFR file
   rm -f temp_obs

   # save a copy of the original day and month, in case we roll over below,
   # in both single/double digit format and guarenteed 2 digit format.
   set oday = $day
   set odd  = $day
   if ($odd < 10) set odd = 0$odd
   set omonth = $month
   set omm    = $month
   if ($omm < 10) set omm = 0$omm

   # convert 1 "day"s worth (data from '3z to 3z of the next day) of BUFR files 
   # into a single intermediate file if 'daily' is set to true; 4 6-hour files
   # otherwise.
   set h = 0
   set next_day = not
   # daily files need to pull out observations at exactly 3Z from the next day.
   # non-daily files should just process up to hour 24, one for one with the inputs.
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
      set oh = $h   # hour not modulo 24
      if ($h < 10)  set oh = 0$h


      # link(big endian) or make(little endian) input file 'prepqm' 
      # for prepbufr.x.  if the pattern for the prepqm files is different,
      # fix the BUFR_file below to match.
      set BUFR_file = ${BUFR_dir}/${year}${mm}/prepqm${yy}${mm}${dd}${hh}
      if (! -e ${BUFR_file}) then
         echo "MISSING FILE ${BUFR_file} and aborting"
         exit -1
      endif

      if ($convert == 'yes') then
         echo "converting ${BUFR_file} to littleendian prepqm.in"
         cp -f ${BUFR_file} prepqm.bigendian
         ../exe/grabbufr.x prepqm.bigendian prepqm.littleendian
         mv prepqm.littleendian prepqm.in
         rm prepqm.bigendian
      else
         echo "linking ${BUFR_file} to prepqm.in"
         ln -f ${BUFR_file} prepqm.in
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
         # can run from a normal 0Z to 23:50Z.
         ../exe/prepbufr_03Z.x
      else
         ../exe/prepbufr.x
      endif

      if ($daily == 'yes') then
         cat prepqm.out >>! temp_obs
         rm prepqm.out
      else
         # use the dates that do not roll over at month/year end, with hours
         set BUFR_out = ${BUFR_dir}/${oyear}${omm}/temp_obs.${oyear}${omm}${odd}${oh}
         echo "moving output to ${BUFR_out}"
         mv -v prepqm.out ${BUFR_out}
      endif
   end

   if ($daily == 'yes') then
      # use the original dates without rollover
      set BUFR_out = ${BUFR_dir}/${oyear}${omm}/temp_obs.${oyear}${omm}${odd}
      echo "moving output to ${BUFR_out}"
      mv -v temp_obs ${BUFR_out}
   endif

   @ day++
end

exit 0
