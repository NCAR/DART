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

#BSUB -o prepbufr.out
#BSUB -e prepbufr.err
#BSUB -J prepbufr
#BSUB -q regular
#BSUB -P 86850054
#BSUB -n 1

#--------------------------------------------------------------
# USER SET PARAMETERS

# set echo

# Convert from big-endian BUFR files to little-endian for Intel chip systems.
# ('yes' or whatever)
set  convert = yes 
set     year = 1997
set    month = 12      
set beginday = 30
#
# end day (up to and including the last day of the month.  Leap year Februaries are OK.
#          Remember that the prepqm###### file for hour 0 of the first day of the next
#          month is necessary for endday = last day of a month.)
#
set endday = 31

# END USER SET PARAMETERS
#--------------------------------------------------------------

set days_in_mo = (31 28 31 30 31 30 31 31 30 31 30 31)
# leap years 
if (($year % 4) == 0) @ days_in_mo[2] = $days_in_mo[2] + 1

if ( $?LS_SUBCWD ) then
   cd $LS_SUBCWD
endif

rm prepqm.out temp_obs  *.err *.out

# Loop over days

set day = $beginday
set last = $endday
while ( $day <= $last )
   echo '-------------------------------------- '

   # clear any old intermediate (text) BUFR file
   rm temp_obs

   # convert 1 "day"s worth (data from '6z to 6z of the next day) of BUFR files 
   #    into a single intermediate file.
   set h = 0
   set next_day = not
   while ($h < 30)
      echo ' '
      @ h  = $h + 6
      @ hh = $h % 24
      @ dd = $day + $h / 24
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
               # next year
               set mm = 1
               @ yy = ($year + 1) % 100
            endif
         endif
      endif

      # format the date for filename
      if ($yy < 10) set yy = 0$yy
      if ($mm < 10) set mm = 0$mm
      if ($dd < 10) set dd = 0$dd
      if ($hh < 10) set hh = 0$hh

      # link(big endian) or make(little endian) input file 'prepqm' for prepbufr.x
      if (! -e ../data/prepqm${yy}${mm}${dd}${hh}) then
         echo "MISSING FILE ../data/prepqm${yy}${mm}${dd}${hh} and aborting"
         exit
      endif

      if ($convert == 'yes') then
         ln -f -s  ../data/prepqm${yy}${mm}${dd}${hh} prepqm.bigendian
         ../exe/grabbufr.x
         mv prepqm.littleendian prepqm.in
      else
         ln -f -s  ../data/prepqm${yy}${mm}${dd}${hh} prepqm.in
      endif

      if ($h == 30) then
         # scavenge a few stragglers from 6Z of the next day using a special prepbufr program
         ../exe/prepbufr_03Z.x
      else
         ../exe/prepbufr.x
      endif
      cat prepqm.out >>   temp_obs
      rm prepqm.out
   end

   set dd = $day
   if ($dd < 10) set dd = 0$dd

   set mm = $month
   if ($mm < 10) set mm = 0$mm

   mv -v temp_obs   temp_obs.${year}${mm}${dd}

   @ day++
end

exit
