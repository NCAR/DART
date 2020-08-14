#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#--------------------------------------------------------------
# DESCRIPTION:
#    Run create_obs_seq for the timespan described in input.nml
#    If resulting files are ASCII fix the locations of obs at the poles.
#    (porting ASCII to another machine (IBM) allows that machine to add
#    digits to the end of the numbers, resulting in latitudes > 90)
#
#--------------------------------------------------------------
#
# you might need other options set for LSF on your system; 
#  e.g. -P project_charge_code, a different queue name, etc.
#
#BSUB -o create_obs_seq.out
#BSUB -e create_obs_seq.err
#BSUB -J create_obs_seq
#BSUB -q share
#BSUB -W 2:00
#BSUB -P NNNNNNNN
#BSUB -n 1

touch create_obs_seq.log

set  STRING = "1,$ s#,##g"

grep year input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set year = $ensstring[3]

grep month input.nml >! ensstring.$$
set ensstring = `sed -e "$STRING" ensstring.$$`
set mo = $ensstring[3]
if ($mo < 10) set mo = 0$mo

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

rm -f ensstring.$$

echo "year, mo, day1, dayn, binary = " $year $mo $day1 $dayn $binary >> create_obs_seq.log
echo ObsBase = $ObsBase                                              >> create_obs_seq.log

./create_real_obs

set n = $day1
while ($n <= $dayn)
   set d = $n
   if ($n < 10) set d = 0$d
   
   set base = obs_seq${year}${mo}${d}
   if ($daily_file == .true.) then
      set obs_seq = ( ${base} )
   else
      # originally, if not doing daily files, the default was 2 files per day.
      # now the default matches the ncep data file blocking; 4 files per day,
      # 6 hours per file.  if you want to go back to 2 files per day, you will
      # have to modify the following line to have only 12 and 24.
      set obs_seq = ( ${base}06 ${base}12 ${base}18 ${base}24 )
   endif
   foreach obs ($obs_seq)
      if ($binary == .false.) then
         if (-e $obs) then
            # on at least one platform (ibm power5), this value is read in
            # and rounded to be 1 bit in the least significant digit larger
            # than 90.000 in degrees, which causes errors in the location module.
            # drop these values just slightly so they read in as 90/-90 exactly.
            echo " fixing $obs pole locations"                     >> create_obs_seq.log
            sed -e 's/ 1.57079632679490/ 1.57079632679488/' \
                -e 's/-1.57079632679490/-1.57079632679488/' $obs   >! fixed_pole${n}
            mv -f fixed_pole${n} $obs
         else
            exit
         endif
      endif
      mv $obs $Obs_base:h &
   end
   @ n++
end

exit 0


