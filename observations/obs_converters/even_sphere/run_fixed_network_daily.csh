#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# loop, calling ./create_fixed_network_seq to create separate
# files for each time period.  edit the values below to change
# the dates and intervals.

@ month = 8
@ day   = 1
@ ndays = 31

while($day <= $ndays)
   echo '"set_def.out"' > create_fixed_input
   echo 1 >> create_fixed_input
   echo 2 >> create_fixed_input
   echo "2008 ${month} ${day} 12 0 0" >> create_fixed_input
   echo '0 43200' >> create_fixed_input

   if($day < 10) then
      echo obs_seq2008080${day} >> create_fixed_input
   else
      echo obs_seq200808${day} >> create_fixed_input
   endif

   ./create_fixed_network_seq < create_fixed_input
   @ day++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

