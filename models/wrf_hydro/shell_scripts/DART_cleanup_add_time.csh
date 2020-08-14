#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#-------------------------------------------------------------------------------
# DESCRIPTION:
#
# Since each DART diagnostic file at each timestep has been written to a 
# separate file, it is desirable to consolidate these files. Each DART
# file has a 'time' dimension and coordinate variable, but the data variables
# do not have a 'time' dimension. This script adds the 'time' dimension
# to the variables.
#
# requires 1 argument:
#    $1 - the DATE string in the form YYYYMMDDHHmm
#-------------------------------------------------------------------------------

if ($# != 1) then
   echo "usage: $0 <YYYYMMDDHHmm>"
   exit -1
endif 

set COMPACTDATE = $1

#===============================================================================
# every (non-coordinate) variable in the output directory gets a time dimension

\rm -f file_list.${COMPACTDATE}.txt
\rm -f output/${COMPACTDATE}/*.out.nc
\rm -f output/${COMPACTDATE}/*.out1.nc
\rm -f output/${COMPACTDATE}/*.out2.nc

ls -1 output/${COMPACTDATE}/*.nc >! file_list.${COMPACTDATE}.txt
foreach FILE ( `cat file_list.${COMPACTDATE}.txt` ) 
   set   FBASE = $FILE:r
   set OUTPUT = ${FBASE}.out.nc
   ncecat -O -h -H -u time ${FILE} ${OUTPUT} || exit 2
end

\rm -f file_list.${COMPACTDATE}.txt

exit 0

#===============================================================================
