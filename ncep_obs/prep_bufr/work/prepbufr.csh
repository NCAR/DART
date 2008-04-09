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
#BSUB -q share
#BSUB -W 2:00
#BSUB -P NNNNNNNN
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

# set echo

# Convert from big-endian BUFR files to little-endian for Intel chip systems.
# ('yes' or whatever)
set    daily = yes
set  convert = no 
set     year = 1988
set    month = 12
set beginday = 1
#
# end day (up to and including the last day of the month.  
#  Leap year Februaries are OK.
#  Remember that the prepqm###### file for hour 0 of the first day of the next
#  month is necessary for endday = last day of a month.)
#
set endday = 31

# Location of BUFR files (named prepqmYYYYMMDDHH)
# are assumed to be in subdirectories named YYYYMM of the path listed here.
# Those subdirectory names will be constructed below.
set BUFR_dir = ../data/
set get_year = $year

# END USER SET PARAMETERS
#--------------------------------------------------------------

set days_in_mo = (31 28 31 30 31 30 31 31 30 31 30 31)
# leap years - year 2000 makes this matter that you do the centuries right
if (($year %   4) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1
if (($year % 100) == 0) @ days_in_mo[2] = $days_in_mo[2] - 1
if (($year % 400) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1

rm -f prepqm.out *.err *.out

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
      echo ' '
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
                 @ get_year ++
               endif
               # next year
               set mm = 1
               @ yy = $get_year % 100
            endif
         endif
      endif

      # format the date for filename
      if ($yy < 10) set yy = 0$yy
      if ($mm < 10) set mm = 0$mm
      if ($dd < 10) set dd = 0$dd
      if ($hh < 10) set hh = 0$hh

      # link(big endian) or make(little endian) input file 'prepqm' 
      # for prepbufr.x
      set BUFR_loc = ${BUFR_dir}/${get_year}${mm}
      if (! -e ${BUFR_loc}/prepqm${yy}${mm}${dd}${hh}) then
         echo "MISSING FILE ${BUFR_loc}/prepqm${yy}${mm}${dd}${hh} and aborting"
         exit
      endif

      if ($convert == 'yes') then
         echo "copying bigendian to ${BUFR_loc}/prepqm${yy}${mm}${dd}${hh}"
         cp ${BUFR_loc}/prepqm${yy}${mm}${dd}${hh} prepqm.bigendian
         ls -l prepqm.bigendian
         ../exe/grabbufr.x prepqm.bigendian prepqm.littleendian
         mv prepqm.littleendian prepqm.in
         rm prepqm.bigendian
      else
         echo "linking prepqm.in to ${BUFR_loc}/prepqm${yy}${mm}${dd}${hh}"
         ln -f ${BUFR_loc}/prepqm${yy}${mm}${dd}${hh} prepqm.in
      endif

      if ($h == 30) then
         # scavenge a few stragglers from 6Z of the next day 
         # using a special prepbufr program
         ../exe/prepbufr_03Z.x
      else
         ../exe/prepbufr.x
      endif

      if ($daily == 'yes') then
         cat prepqm.out >>! temp_obs
         rm prepqm.out
      else
         mv -v prepqm.out ${BUFR_dir}/${year}${mm}/temp_obs.${year}${mm}${dd}${hh}
      endif
   end

   if ($daily == 'yes') then
      set dd = $day
      if ($dd < 10) set dd = 0$dd

      set mm = $month
      if ($mm < 10) set mm = 0$mm
   
      mv -v temp_obs   ${BUFR_dir}/${year}${mm}/temp_obs.${year}${mm}${dd}
   endif

   @ day++
end

exit
