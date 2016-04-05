#!/bin/csh
#
# Given one gcom_restart.nc file with multiple timesteps in it,
# this will create a series of netCDF files - one for each timestep.
# In this case, gcom_restart.nc has 15 timesteps (ncks uses C-style indexing)
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

set member = 1
while ($member <= 15)

   @ timestep = $member - 1

   set fname = `printf gcom_restart_%04d.nc $member`

   \rm -f $fname bob.nc bob.$$.txt

   echo -n "Creating $fname ... "

   # subset out two consecutive timesteps

   ncks -O -d time,$timestep gcom_restart.nc bob.nc

   # overwrite the timestep value with zero
   # so that they are all 'valid' at the same time.
   # Until there is a better way to create an initial
   # ensemble with some spread - this will have to do.
   # It may not be sufficient.

   ncdump bob.nc | sed -e "/ time = /c\ time = 0 ;" > bob.$$.txt

   ncgen -o $fname bob.$$.txt

   echo "done."

   \rm -f bob.nc bob.$$.txt

   @ member++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

