#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Call ./create_fixed_network_seq to create a separate file for each time period.  
# Edit the values below to change the dates.
# It makes 1 file/day, single time per file, at 12 UTC.
# It assumes that create_fixed_network has any model-specific files it needs in this directory.
# It requires a set_def.out file (usually created by create_obs_sequence).

@ year  = 2008
@ month = 8
@ day   = 1
@ ndays = 31

while($day <= $ndays)
   echo '"set_def.out"' > create_fixed_input
   echo 1 >> create_fixed_input
   echo 2 >> create_fixed_input
   echo "${year} ${month} ${day} 12 0 0" >> create_fixed_input

   echo '0 43200' >> create_fixed_input

   set dstring = `printf %04d%02d%02d $year $month $day`
   echo obs_seq${dstring} >> create_fixed_input


   ./create_fixed_network_seq < create_fixed_input
   @ day++
end

exit 0


