#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$


echo 
echo 
echo "=================================================================="
echo "Start of observation converter tests at "`date`
echo "=================================================================="
echo 
echo 

set startdir=`pwd`

set LOGDIR=${startdir}/testing_logs
echo putting build and run logs in:
echo $LOGDIR

mkdir -p ${LOGDIR}
\rm -f ${LOGDIR}/*

echo 
echo 
echo "=================================================================="
echo "Compiling NCEP BUFR libs starting at "`date`
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
   else if ( "$fcomp" == "nagfor" ) then
      setenv CCOMP nag
      setenv FCOMP nag
      setenv UNDERSCORE add
      echo setting the BUFR lib to build using the nag compilers
   else
      echo unrecognized compiler in ../../build_templates/mkmf.template
      echo set NCEP BUFR library compiler choice in NCEP/prep_bufr/install.sh
      echo this script will use whatever compiler is selected there
   endif
endif

cd NCEP/prep_bufr

set FAILURE = 0

( ./install.sh > ${LOGDIR}/buildlog.NCEP.out ) || set FAILURE = 1

echo 
echo 
echo "=================================================================="
echo "Build of NCEP BUFR libs ended at "`date`
if ( $FAILURE ) then
      echo 
      echo "ERROR - build was unsuccessful"
      echo 
endif
echo "=================================================================="
echo 
echo 

cd $startdir

foreach quickb ( `find . -name quickbuild.csh -print` )

   cd $startdir

   # get the working dir name. also, make a project name by stripping off
   # the leading ./ and the /work parts of the dirname, and turning slashes
   # into underscores so we can use the string as part of a log filename.
   set wdir = $quickb:h
   set project = `echo $wdir | sed -e 's;^./;;' -e 's;/[^/]*$;;' -e 's;/;_;g'`

   echo 
   echo 
   echo "=================================================================="
   echo "Compiling obs converter $project starting at "`date`
   echo "=================================================================="
   echo 
   echo 


   cd $wdir
   echo
   echo building in $wdir

   set FAILURE = 0
   ( ./quickbuild.csh > ${LOGDIR}/buildlog.${project}.out ) || set FAILURE = 1

   echo

   if ( $FAILURE ) then
      echo "ERROR - unsuccessful build of $project at "`date`
      echo 

      switch ( $project )
   
         case GSI2DART
            echo " This build expected to fail on case-insensitive filesystems."
         breaksw
            
         case var
            echo " This build expected to fail unless you have the WRF code in-situ."
         breaksw
            
         case AIRS
            echo " AIRS build is expected to fail due to dependency on hdfeos libs,"
            echo " which are not required to be part of the standard DART environment."
         breaksw
            
         case quikscat
            echo " quikscat build is expected to fail due to dependency on mfhdf libs,"
            echo " which are not required to be part of the standard DART environment."
         breaksw
  
         default
            echo " unexpected error"
         breaksw
      endsw
   else
      echo "Successful build of obs converter $project ended at "`date`
      echo 
      echo "Executing converters in directory $wdir"

      \rm -f *.o *.mod
      \rm -f Makefile input.nml.*_default .cppdefs

      foreach TARGET ( mkmf_* )
         set FAILURE = 0
         set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
         echo "Running $PROG"
   
         # for programs which read standard input, put what they need into a prog.in file.
         # if we miss any programs which need input and we don't have a .in file, have it
         # read from /dev/null so it errors out and doesn't just sit there waiting for input
         if ( -f ${PROG}.in ) then
           ( ./$PROG < ${PROG}.in > ${LOGDIR}/runlog.${project}.out ) || set FAILURE = 1
         else
           ( ./$PROG < /dev/null > ${LOGDIR}/runlog.${project}.out ) || set FAILURE = 1
         endif
         if ( $FAILURE ) then
            echo "ERROR - unsuccessful run of $PROG at "`date`
         else
            echo "Successful run of $PROG at "`date`
            \rm -f $PROG
         endif
      end

      echo

   endif

   echo "=================================================================="
   echo
   echo
  
end

echo 
echo 
echo "=================================================================="
echo "End of observation converter tests at "`date`
echo "=================================================================="
echo 
echo 

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

