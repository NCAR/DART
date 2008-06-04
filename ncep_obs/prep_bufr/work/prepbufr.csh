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

# Loop over days

set day = $beginday
set last = $endday
while ( $day <= $last )
   echo '-------------------------------------- '

   # clear any old intermediate (text) BUFR file
   rm -f temp_obs

   # convert 1 "day"s worth (data from '3z to 3z of the next day) of BUFR files 
   # into a single intermediate file if 'daily' is set to true; 4 6-hour files
   # otherwise.
   set h = 0
   set next_day = not
   while ($h < 30)
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

      # format the date for filename
      if ($yy < 10) set yy = 0$yy
      if ($mm < 10) set mm = 0$mm
      if ($dd < 10) set dd = 0$dd
      if ($hh < 10) set hh = 0$hh

      # link(big endian) or make(little endian) input file 'prepqm' 
      # for prepbufr.x.  if the pattern for the prepqm files is different,
      # fix the BUFR_file below to match.
      set BUFR_file = ${BUFR_dir}/${year}${mm}/prepqm${yy}${mm}${dd}${hh}
      set BUFR_out  = ${BUFR_dir}/${year}${mm}/temp_obs.${year}${mm}${dd}
      if (! -e ${BUFR_file}) then
         echo "MISSING FILE ${BUFR_file} and aborting"
         exit
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
         # get any obs between 0Z and 3Z from the 6Z file of the next day 
         # using a slightly modified prepbufr program
         ../exe/prepbufr_03Z.x
      else
         ../exe/prepbufr.x
      endif

      if ($daily == 'yes') then
         cat prepqm.out >>! temp_obs
         rm prepqm.out
      else
         echo "moving output to ${BUFR_out}${hh}"
         mv -v prepqm.out ${BUFR_out}${hh}
      endif
   end

   if ($daily == 'yes') then
      set dd = $day
      if ($dd < 10) set dd = 0$dd

      set mm = $month
      if ($mm < 10) set mm = 0$mm

      echo "moving output to ${BUFR_out}"
      mv -v temp_obs ${BUFR_out}
   endif

   @ day++
end

exit
