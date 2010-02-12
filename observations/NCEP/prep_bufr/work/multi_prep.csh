#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Driver script for the parallel version.  Submit this script
# to your batch system and it will invoke the 'multi_body.csh'
# script once for each conversion day.
#
#--------------------------------------------------------------
# DESCRIPTION:
#
#  This script is used to generate daily (3:01Z to 3:00Z of next day) decoded 
#  NCEP reanalysis PREPBUFR text/ascii data.
#
#--------------------------------------------------------------

#BSUB -o multiprep.out
#BSUB -e multiprep.err
#BSUB -J multiprep
#BSUB -q regular
#BSUB -W 0:15
#BSUB -P xxxxxxxx
#BSUB -n 1

# Set year, month is multi_body.csh, set days only here.

# Only 4 inputs which need to be set in this script.
set year     = 2004
set month    =    1
set startday =    1
set endday   =   30

# Subdirectory base name.  If running multiple months at the
# same time, must have different base names.
set prefix = tempdir

# Loop over days, first to convert, then to clean up
set day = $startday
set lastday = $endday

while ( $day <= $lastday )

   mkdir ${prefix}_${day}
   cd ${prefix}_${day}

   cp -f ../multi_body.csh .
   cp -f ../input.nml .

   echo starting day $day
   ./multi_body.csh $year $month $day $day &

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
 rm ${prefix}_${day}/input.nml ${prefix}_${day}/multi_body.csh
 rmdir ${prefix}_${day} 
 @ day++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

