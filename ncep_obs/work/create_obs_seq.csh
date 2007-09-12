#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL: http://subversion.ucar.edu/DAReS/DART/trunk/ncep_obs/prep_bufr/work/prepbufr.csh $
# $Id: prepbufr.csh 3042 2007-08-02 22:03:02Z thoar $
# $Revision: 3042 $
# $Date: 2007-08-02 16:03:02 -0600 (Thu, 02 Aug 2007) $

#--------------------------------------------------------------
# DESCRIPTION:
#    Run create_obs_seq for the timespan described in input.nml
#    If resulting files are ASCII fix the locations of obs at the poles.
#    (porting ASCII to another machine (IBM) allows that machine to add
#    digits to the end of the numbers, resulting in latitudes > 90)
#
#--------------------------------------------------------------

# you might need other options set for LSF on your system; 
#  e.g. -P project_charge_code, a different queue name, etc.

#BSUB -o create_obs_seq.out
#BSUB -e create_obs_seq.err
#BSUB -J create_obs_seq
#BSUB -q share
#BSUB -W 6:00
#BSUB -n 1


if ($?LS_SUBCWD) then
   cd $LS_SUBCWD
endif

touch create_obs_seq.out

set  STRING = "1,$ s#,##g"

grep year input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set year = $ensstring[3]

grep month input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set mo = $ensstring[3]

grep day input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set day1 = $ensstring[3]

grep tot_days input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set tot_days = $ensstring[3]
@ dayn = $day1 + $tot_days - 1

grep write_binary_obs_sequence input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set binary = $ensstring[3]

grep ObsBase input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set ObsBase = $ensstring[3]

grep daily_file input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set daily_file = $ensstring[3]

echo "year, mo, day1, dayn, binary = " $year $mo $day1 $dayn $binary >> create_obs_seq.out
echo ObsBase = $ObsBase                                              >> create_obs_seq.out

./create_real_obs

set n = $day1
while ($n <= $dayn)
   set d = $n
   if ($n < 10) set d = 0$d
   
   if ($daily_file == .true.) then
      set obs_seq = (obs_seq${year}${mo}${d})
   else
      set obs_seq = (obs_seq${year}${mo}${d}12 obs_seq${year}${mo}${d}24)
   endif
   foreach obs ($obs_seq)
      if ($binary == .false.) then
         if (-e $obs) then
            echo "fixing $obs pole locations"                    >> create_obs_seq.out
            sed -e 's/1.57079632679490/1.57079632679488/' $obs   >! fixed_pole
         else
            exit
         endif
         mv fixed_pole $obs
      endif
      mv $obs $Obs_base:h &
   end
   @ n++
end


exit

