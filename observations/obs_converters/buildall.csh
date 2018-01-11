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

echo 
echo 
echo "=================================================================="
echo "=================================================================="
echo "Compiling NCEP BUFR libs starting at "`date`
echo "=================================================================="
echo "=================================================================="
echo 
echo 

# the NCEP bufr libs are needed by at least one other converter
# (the gps bufr one) and they build differently with their own
# install.sh script.  we do most of our testing with intel and
# gfortran, so try to figure out from our mkmf.template file
# which one is being used and set env vars so the install script
# will use the right one.  this doesn't cover every compiler
# we support but it will get 80% of the cases with 20% of the work.

if ( -f ../../build_templates/mkmf.template ) then
   set fcomp=`grep '^FC' ../../build_templates/mkmf.template | sed -e 's/FC *= *\([A-Za-z][^ ]*\)/\1/' `
   if ( "$fcomp" == "ifort" ) then
      setenv CCOMP intel
      setenv FCOMP intel
      setenv UNDERSCORE add
      echo setting the BUFR lib to build using the intel compilers
   else if ( "$fcomp" == "gfortran" ) then
      setenv CCOMP gnu
      setenv FCOMP gnu
      setenv UNDERSCORE add
      echo setting the BUFR lib to build using the gnu compilers
   else if ( "$fcomp" == "pgf90" ) then
      setenv CCOMP pgi
      setenv FCOMP pgi
      setenv UNDERSCORE add
      echo setting the BUFR lib to build using the pgi compilers
   else
      echo unrecognized compiler in ../../build_templates/mkmf.template
      echo set NCEP BUFR library compiler choice in NCEP/prep_bufr/install.sh
      echo this script will use whatever compiler is selected there
   endif
endif

cd NCEP/prep_bufr

./install.sh

echo 
echo 
echo "=================================================================="
echo "=================================================================="
echo "Build of NCEP BUFR libs ended at "`date`
echo "=================================================================="
echo "=================================================================="
echo 
echo 

cd $startdir

foreach project ( `find . -name quickbuild.csh -print` )

   cd $startdir

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

   echo
   echo
   echo "=================================================================="
   echo "=================================================================="

   if ( $FAILURE ) then
      echo "ERROR - unsuccessful build in $dir at "`date`
      echo 

      switch ( $dir )
   
         case *GSI2DART*
            echo " This build expected to fail on case-insensitive filesystems."
         breaksw
            
         case */var/*
            echo " This build expected to fail unless you have the WRF code in-situ."
         breaksw
            
         case *AIRS*
            echo " AIRS build is expected to fail due to dependency on hdfeos libs,"
            echo " which are not required to be part of the standard DART environment."
         breaksw
            
         case *quikscat*
            echo " quikscat build is expected to fail due to dependency on mfhdf libs,"
            echo " which are not required to be part of the standard DART environment."
         breaksw
  
         default
            echo " unexpected error"
         breaksw
      endsw
   else
      echo "Successful build of obs converter $dir ended at "`date`
   endif

   echo "=================================================================="
   echo "=================================================================="
   echo
   echo
  
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

