#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Convert multiple days of prepbufr data into ascii intermediate files
# in a loop.  
#
# WARNING WARNING WARNING: This script may not be what you want 
# to use!!!  On most batch systems it does NOT in fact run
# these jobs in parallel.  It also only loops over days within
# a single month.  If you aren't running with a batch system this
# script does try to start multiple conversions in the background.
#
# See the more complicated (but functionally correct) 
# multi_parallel.batch script for how to use MPI to run
# multiple serial jobs at the same time.  
#
#--------------------------------------------------------------
# DESCRIPTION:
#
#  This script is used to read PREPBUFR data files 
#  and output ascii/text intermediate files.
#
#--------------------------------------------------------------

# USER SETTINGS HERE

# Set year, month, days for to pass as args to prepbufr.csh
set year     = 2010
set month    =   12
set startday =    1
set endday   =   30

# Subdirectory base name.  If running multiple months at the
# same time, must have different base names.  (could add month
# number onto tempdir names in script below.)

set prefix = tempdir

# END USER SETTINGS

# Loop over days, first to convert, then to clean up
set day = $startday
set lastday = $endday

while ( $day <= $lastday )

   mkdir ${prefix}_${day}
   cd ${prefix}_${day}

   cp -f ../prepbufr.csh .
   cp -f ../input.nml .

   echo starting day $day
   ./prepbufr.csh $year $month $day $day &

   cd ..

   @ day++
end

wait

set day = $startday
set lastday = $endday

# if everything ran ok, the directory should be empty.
# rather than do rm -fr dir, remove what we think is there
# and have it fail, intentionally, if other files still exist.
while ( $day <= $lastday )

 rm ${prefix}_${day}/input.nml ${prefix}_${day}/prepbufr.csh
 rmdir ${prefix}_${day} 

 @ day++
end

exit 0


