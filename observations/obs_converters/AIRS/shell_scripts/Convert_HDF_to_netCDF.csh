#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download


cd ../data

foreach FILE ( *hdf )

   set BASE = $FILE:r
   set NEWNAME = $BASE.nc

   echo 
   echo "Converting $FILE to"
   echo "           $NEWNAME"
   echo 

   \rm -f bob.nc
   h4tonccf_nc4 $FILE bob.nc || exit 1
   ncatted -a coremetadata,global,d,,, -a StructMetadata_0,global,d,,, bob.nc $NEWNAME

end

\rm -f  bob.nc

