#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#==========================================================================

set SNAME = $0
set FNAME = $1

if ($#argv != 1) then

   echo " "
   echo "usage: $SNAME:t <some_file_name> "
   echo " "
   echo "Reads the 'ocean_time' variable from the file"
   echo "(which has units of seconds) and prints the number of days."
   echo " "
   exit 1

endif

if (! -e $FNAME) then
   echo "$FNAME does not exist."
   exit 1
endif

# double ocean_time(ocean_time) ;
#        ocean_time:long_name = "time since initialization" ;
#        ocean_time:field = "time, scalar, series" ;
#        ocean_time:units = "seconds since 1900-01-01 00:00:00" ;
#        ocean_time:time_origin = "01-JAN-1900 00:00:00" ;

# NOTE ... bc can handle the 'long' integers that happen when the
# reference time is 1900-01-01

set OCEAN_TIME = `ncdump -v ocean_time ${FNAME} | grep "ocean_time =" | tail -1`
set TIME_SEC = `echo $OCEAN_TIME | grep -oE '[[:digit:]]+'`
set TIME_FLOAT = `echo "scale=6; $TIME_SEC / 86400.0" | bc `
set TIME_INTEGER = `echo "scale=0; $TIME_FLOAT / 1" | bc `

#echo "$OCEAN_TIME"
#echo "The ocean_time is $TIME_SEC seconds or $TIME_INTEGER days"

echo $TIME_INTEGER

exit 0


