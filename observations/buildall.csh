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

   cd $dir
   echo
   echo building in $dir

   switch ("$dir")

      case */var/*
         echo "expected to fail unless you have the WRF code in-situ."
      breaksw
         
      case *AIRS*
         ./quickbuild.csh
         echo "AIRS build is expected to fail due to dependency on hdfeos libs,"
         echo "not required to be part of the standard DART environment."
      breaksw
         
      case *quikscat*
         ./quickbuild.csh
         echo "quikscat build is expected to fail due to dependency on mfhdf libs,"
         echo "not required to be part of the standard DART environment."
      breaksw
         
      default:
         ./quickbuild.csh || set FAILURE = 1
      breaksw
         
   endsw

   if ( $FAILURE != 0 ) then
      echo
      echo "ERROR unsuccessful build in $dir"
      echo "ERROR unsuccessful build in $dir"
      echo "ERROR unsuccessful build in $dir"
      echo
#     exit -1
   endif

   cd $startdir
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

