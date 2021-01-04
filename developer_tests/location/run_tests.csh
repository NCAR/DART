#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# this script builds and  runs the location test code for each of the
# possible location modules.
#
#----------------------------------------------------------------------

# prevent shell warning messages about no files found when trying
# to remove files using wildcards.

set nonomatch

echo
echo
echo "=================================================================="
echo "Start of location module tests at "`date`
echo "=================================================================="
echo
echo

set LOGDIR = `pwd`/testing_logs
mkdir -p $LOGDIR
\rm -f $LOGDIR/*
echo "build and run logs are in: $LOGDIR"

set LOCLIST = ( annulus channel column oned threed \
                threed_cartesian threed_sphere \
                twod twod_annulus twod_sphere )

foreach i ( $LOCLIST )

   echo
   echo
   echo "------------------------------------------------------------------"
   echo "Starting tests of location module $i at "`date`
   echo "------------------------------------------------------------------"
   echo
   echo

   set FAILURE = 0

   cd $i/test

   ./mkmf_location_test
   ( make > $LOGDIR/buildlog.$i.out ) || set FAILURE = 1

   echo
   echo
   if ( $FAILURE ) then
     echo "------------------------------------------------------------------"
     echo "ERROR - unsuccessful build of location module $i at "`date`
     echo "------------------------------------------------------------------"
     echo
     echo
   else

      ls -l location_test
      ( ./location_test < test.in > $LOGDIR/runlog.$i.out ) || set FAILURE = 1

     echo
     echo
     if ( $FAILURE ) then
       echo "------------------------------------------------------------------"
       echo "ERROR - unsuccessful run of location module $i tests at "`date`
       echo "------------------------------------------------------------------"
     else
       echo "------------------------------------------------------------------"
       echo "Tests of location module $i complete at "`date`
       echo "------------------------------------------------------------------"
     endif
     echo
     echo

     \rm -f *.o *.mod input.nml*_default dart_log.* \
            Makefile location_test_file* location_test

   endif

   cd ../..
end

echo
echo
echo "=================================================================="
echo "End of location module tests at "`date`
echo "=================================================================="
echo
echo

exit 0

