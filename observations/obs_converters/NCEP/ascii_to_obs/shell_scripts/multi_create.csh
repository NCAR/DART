#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# test version - convert multiple days in parallel
#
#--------------------------------------------------------------
# DESCRIPTION:
#   Run create_obs_seq for the timespan described in input.nml
#   Up to a month at a time, in parallel.
#
#--------------------------------------------------------------

# you might need other options set for LSF on your system; 
#  e.g. -P project_charge_code, a different queue name, etc.

#BSUB -J createobs
#BSUB -o createobs.out
#BSUB -e createobs.err
#BSUB -q regular
#BSUB -P xxxxxxxx
#BSUB -W 0:20
#BSUB -n 1

# USER INPUTS - SET HERE
# Set the year, month, and day range here.
# this version of the script doesn't wrap across month boundaries
set     year = 1999
set    month =   12
set startday =    1
set   endday =   31

# dir to move the output files into - will have YYYYMM appended.
set destdir = ./ACARS_24_ascii

# if you want to run more than one of these jobs at a time,
# each must have a different base tempname.
set tempname = tempdir

# END USER INPUTS

# Loop over days, first to convert, then to clean up
set day = $startday
set lastday = $endday

# this version doesn't wrap across months yet.
set mo = `printf %02d $month`

while ( $day <= $lastday )

   set dy = `printf %02d $day`

   if ( ! -d ${tempname}_${dy}) mkdir ${tempname}_${dy}
   cd ${tempname}_${dy}

   cp -f ../create_real_obs .
   sed -e "s/YYYY/${year}/" \
       -e "s/MM/${mo}/"     \
       -e "s/NN/${day}/"    < ../input.nml.template >! input.nml

   echo starting day $day
   ./create_real_obs >! obs.log &

   cd ..

   @ day++
end

wait

echo all obs_seq files created at `date`

set day = $startday
set lastday = $endday

# remove the various bits, and either list or move the new obs_seq
# file to a destination directory.  don't remove the tempdirs so
# you don't miss errors as easily.

set dest = ${destdir}/${year}${mo}
if ( ! -d ${dest}) mkdir ${dest}

while ( $day <= $lastday )
 set dy = `printf %02d $day`

 # mv the obs file to the eventual destination dir here
 mv ${tempname}_${dy}/obs_seq${year}${mo}${dy} ${dest}/
 
 # or, just ls it:
 # ls -l ${tempname}_${dy}/obs_seq${year}${mo}${dy}

 rm ${tempname}_${dy}/input.nml ${tempname}_${dy}/create_real_obs ${tempname}_${dy}/dart_log.nml
 #rmdir ${tempname}_${day} 
 @ day++
end

ls -l ${dest}

exit 0


