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

# this one makes 2 files/day, 12 hours apart, single time per file.

@ year   = 2008
@ month  = 8
@ day    = 1
@ ndays  = 3
@ nhours = 12


while($day <= $ndays)
   @ hour = 0
   while ($hour < 24)
      echo '"set_def.out"' > create_fixed_input
      echo 1 >> create_fixed_input
      echo 1 >> create_fixed_input
      echo $year $month $day $hour 0 0 >> create_fixed_input
      echo '0 0' >> create_fixed_input
   
      set dstring = `printf %04d%02d%02d%02d $year $month $day $hour`
      echo obs_seq${dstring} >> create_fixed_input
   
      ./create_fixed_network_seq < create_fixed_input

      @ hour += $nhours
   end

   @ day++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

