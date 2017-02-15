#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

set SNAME = $0
set clobber

set startdir=`pwd`

foreach project ( `find . -name quickbuild.csh -print` )

   set dir = $project:h
   set FAILURE = 0

   echo 
   echo 
   echo "=================================================================="
   echo "=================================================================="
   echo "Compiling obs converter $dir starting at "`date`
   echo "=================================================================="
   echo "=================================================================="
   echo 
   echo 


   cd $dir
   echo
   echo building in $dir

   ./quickbuild.csh || set FAILURE = 1

   if ( $FAILURE ) then
      echo
      echo "=================================================================="
      echo "ERROR - unsuccessful build in $dir"
      echo 

      switch ( $dir )
   
         case */var/*
            echo "This build expected to fail unless you have the WRF code in-situ."
            echo "=================================================================="
         breaksw
            
         case *AIRS*
            echo "AIRS build is expected to fail due to dependency on hdfeos libs,"
            echo "not required to be part of the standard DART environment."
            echo "=================================================================="
         breaksw
            
         case *quikscat*
            echo "quikscat build is expected to fail due to dependency on mfhdf libs,"
            echo "not required to be part of the standard DART environment."
            echo "=================================================================="
         breaksw
  
         default
            echo " unexpected error"
            echo "=================================================================="
         breaksw
      endsw
  
   else
         
   echo 
   echo 
   echo "=================================================================="
   echo "=================================================================="
   echo "Build of obs converter $dir ended at "`date`
   echo "=================================================================="
   echo "=================================================================="
   echo 
   echo 

   cd $startdir
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

