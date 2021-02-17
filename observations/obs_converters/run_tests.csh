#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# At present, the GSI2DART converter MUST be compiled with MPI,
# and we need to get the mpi_launch_command to run it.
# Consequently, I am making the interface to this script look similar
# to the developer_tests/run_tests.csh, which is a bit of overkill.
# The mkmf_gsi_to_dart compiles with MPI even without the -mpi flag,
# but I don't want to pass the -mpi flag to all the quickbuild.csh
# scripts in the obs_converters directory ...
#
# usage: [ -mpi | -nompi ] [ -mpicmd name_of_mpi_launch_command ]
#
#----------------------------------------------------------------------

# prevent shell warning messages about no files found when trying
# to remove files using wildcards.
set nonomatch

set usingmpi=no
set MPICMD=""

if ( $#argv > 0 ) then
  if ( "$argv[1]" == "-mpi" ) then
    set usingmpi=yes
  else if ( "$argv[1]" == "-nompi" ) then
    set usingmpi=no
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi ]  [ -mpicmd name_of_mpi_launch_command ]"
    echo " default is to run tests without MPI"
    exit -1
  endif
  shift
endif

if ( $#argv > 1 ) then
  if ( "$argv[1]" == "-mpicmd" ) then
    set MPICMD = "$argv[2]"
  else
    echo "Unrecognized argument to $0: $argv[1]"
    echo "Usage: $0 [ -mpi | -nompi ]  [ -mpicmd name_of_mpi_launch_command ]"
    echo " default is to run tests without MPI"
    exit -1
  endif
  shift
endif

# set the environment variable MPI to anything in order to enable the
# MPI builds and tests.  set the argument to the build scripts so it
# knows which ones to build.
if ( "$usingmpi" == "yes" ) then
  echo "Building with MPI support."
  set QUICKBUILD_ARG='-mpi'
  if ( ! $?MPICMD) then
    set MPICMD='mpirun -n 2'
  endif
  echo "MPI programs will be started with: $MPICMD"
else if ( "$usingmpi" == "no" ) then
  echo "Building WITHOUT MPI support."
  set QUICKBUILD_ARG='-nompi'
  set MPICMD=""
else
  echo "Internal error: unrecognized value of usingmpi; should not happen"
  exit -1
endif

if ( ! $?REMOVE) then
   setenv REMOVE 'rm -f'
endif
if ( ! $?REMOVE_DIR) then
   setenv REMOVE_DIR 'rmdir'
endif
if ( ! $?COPY) then
   setenv COPY 'cp -f'
endif
if ( ! $?MOVE) then
   setenv MOVE 'mv -f'
endif

if ( ! $?host) then
   setenv host `uname -n`
endif

echo 
echo 
echo "=================================================================="
echo "Start of observation converter tests at "`date`
echo "=================================================================="
echo 
echo 
echo "Running observation converter tests on $host"

set startdir=`pwd`
set LOGDIR=${startdir}/testing_logs
mkdir -p ${LOGDIR}
${REMOVE} ${LOGDIR}/*
echo "build and run logs are in: $LOGDIR"

echo 
echo "------------------------------------------------------------------"
echo "Building NCEP BUFR libs starting at "`date`
echo "------------------------------------------------------------------"
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
echo "------------------------------------------------------------------"
echo "Build of NCEP BUFR libs ended at "`date`
if ( $FAILURE ) then
      echo 
      echo "ERROR - build was unsuccessful at "`date`
      echo 
endif
echo "------------------------------------------------------------------"
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
   echo "------------------------------------------------------------------"
   echo "Testing obs converter $project starting at "`date`
   echo "------------------------------------------------------------------"
   echo 
   echo 


   cd $wdir
   echo
   echo "building in $wdir"

   # save original input.nml & obs seq files here
   set SAVEDIR = saveme.test_dart
   mkdir -p ${SAVEDIR}
   if ( -e input.nml ) then 
      ${COPY} input.nml ${SAVEDIR}
   endif
   if ( -e obs_seq.* ) then
      ${COPY} obs_seq.* ${SAVEDIR}
   endif

   # If there is a testing namelist, use it.
   if ( -f input.nml.testing ) then
      ${COPY} input.nml.testing input.nml
   endif

   set FAILURE = 0
   ( ./quickbuild.csh > ${LOGDIR}/buildlog.${project}.out ) || set FAILURE = 1
   echo

   if ( $FAILURE ) then
      echo "ERROR - unsuccessful build of $project at "`date`
      echo 

      switch ( $project )
   
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
      echo "Successful build of obs converter $project"
      echo 
      echo "Executing converters in directory $wdir"

      ${REMOVE} *.o *.mod
      ${REMOVE} Makefile input.nml.*_default .cppdefs

      # @todo FIXME ... can skip running preprocess at this point, SHOULD run the
      # observation converter programs (whatever name) BEFORE running obs_sequence_tool
      # as it is, the obs_sequence_tool is failing because the obs_seq.out has not
      # been created yet.

      foreach TARGET ( mkmf_* )

         if ( $TARGET == "mkmf_preprocess" && (-e preprocess)) goto skip

         set FAILURE = 0
         set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
         echo
         echo "Running $PROG"

         # for programs which read standard input, put what they need into a prog.in file
         # in the tests directory.
         # if we miss any programs which need input and we don't have a .in file, have it
         # read from /dev/null so it errors out and doesn't just sit there waiting for input
         if ( -f ../work/${PROG}.in ) then

           if ( -f using_mpi_for_$PROG ) then
              ( ${MPICMD} ./$PROG < ../work/${PROG}.in > ${LOGDIR}/runlog.${project}.${PROG}.out ) || set FAILURE = 1
           else
              (           ./$PROG < ../work/${PROG}.in > ${LOGDIR}/runlog.${project}.${PROG}.out ) || set FAILURE = 1
           endif

         else

           if ( -f using_mpi_for_$PROG ) then
              ( ${MPICMD} ./$PROG < /dev/null > ${LOGDIR}/runlog.${project}.${PROG}.out ) || set FAILURE = 1
           else
              (           ./$PROG < /dev/null > ${LOGDIR}/runlog.${project}.${PROG}.out ) || set FAILURE = 1
           endif
         endif
         if ( $FAILURE ) then

            switch ( $PROG )

               # These programs rely on the HDF-EOS libraries. Spcecifying them in
               # the mkmf_* file allows them to compile, but unless your run-time
               # environment matches ... the execution will fail.
               case L1_AMSUA_to_netcdf
                  echo
                  echo "If $PROG fails due to 'error while loading shared libraries ...'"
                  echo "make sure your [DY]LD_LIBRARY_PATH is consistent with the library"
                  echo "paths in $TARGET. This may still fail for other reasons."
                  echo
               breaksw
            
               case convert_airs_L2
                  echo
                  echo "If $PROG fails due to 'error while loading shared libraries ...'"
                  echo "make sure your [DY]LD_LIBRARY_PATH is consistent with the library"
                  echo "paths in $TARGET. This may still fail for other reasons."
                  echo
               breaksw
            
               case convert_amsu_L1
                  echo
                  echo "If $PROG fails due to 'error while loading shared libraries ...'"
                  echo "make sure your [DY]LD_LIBRARY_PATH is consistent with the library"
                  echo "paths in $TARGET. This may still fail for other reasons."
                  echo
               breaksw
            
               case convert_L2b
                  echo
                  echo "If $PROG fails due to 'error while loading shared libraries ...'"
                  echo "make sure your [DY]LD_LIBRARY_PATH is consistent with the library"
                  echo "paths in $TARGET. This may still fail for other reasons."
                  echo
               breaksw
            
               default
                  echo "ERROR - unsuccessful run of $PROG at "`date`
               breaksw
            endsw

         else
            echo "Successful run of $PROG"
            ${REMOVE} $PROG
         endif

      skip:
      end

   endif

   echo "Restoring original input.nml and obs_seq files"
   ${MOVE} ${SAVEDIR}/* .
   ${REMOVE_DIR} ${SAVEDIR}

end

echo 
echo 
echo "=================================================================="
echo "End of observation converter tests at "`date`
echo "=================================================================="
echo 
echo 

exit 0

