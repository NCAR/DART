#!/bin/csh 
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# TEST VERSION - try to fork up to 32 processes to run
# at once.
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
# $URL: https://subversion.ucar.edu/DAReS/DART/trunk/ncep_obs/prep_bufr/work/prepbufr.csh $
# $Id: prepbufr.csh 3637 2008-10-30 20:31:21Z nancy $
# $Revision: 3637 $
# $Date: 2008-10-30 14:31:21 -0600 (Thu, 30 Oct 2008) $

