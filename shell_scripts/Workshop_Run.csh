#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

set SNAME = $0
set clobber

if ($#argv == 1) then

      setenv LOGDIR $1
      echo "Running $SNAME:t and archiving results to $LOGDIR"
  
else
      echo " "
      echo "usage: $SNAME:t destinationdir"
      echo " "
      echo "This script runs the 'workshop_setup.csh' script for a wide range of models"
      echo "and archives the results in 'destinationdir'. If 'destinationdir' does not"
      echo "exist, it is created. If it does exist, any duplicate contents will be overwritten."
      echo " "
      echo "This must be run from the top-level 'DART' directory."
      echo " "
      echo "This is a pretty verbose process, so if you are logging the output,"
      echo "make sure you have plenty of space:"
      echo " "
      echo "./$SNAME:t |& tee DART_test.log"
      echo " "
      echo "can easily result in a 750 Kb log file"
      echo " "
      exit 1
endif

if ( ! -d models/lorenz_96 ) then
   echo "$SNAME:t must be run from the top-level"
   echo "DART directory -- please try again."
   exit 2
else
   set DARTHOME = `pwd`
endif

setenv ARCHIVEDIR ${DARTHOME}/${LOGDIR}
if ( ! -d $ARCHIVEDIR ) then
   echo "Creating $ARCHIVEDIR"
   mkdir -pv $ARCHIVEDIR
endif

echo "The top-level DART directory (DARTHOME) is $DARTHOME"

#----------------------------------------------------------------------
# Run the workshop_setup demo for a wide range of models.
#----------------------------------------------------------------------

@ makenum  = 1
@ modelnum = 101
foreach MODEL ( 9var lorenz_63 lorenz_84 lorenz_96 lorenz_96_2scale \
   forced_lorenz_96 lorenz_04 bgrid_solo pe2lyr )

    echo "----------------------------------------------------------"
    echo "Running $MODEL at "`date`
    echo ""

    cd ${DARTHOME}/models/${MODEL}/work

    \rm -fv ../../../obs_def/obs_def_mod.f90
    \rm -fv ../../../obs_kind/obs_kind_mod.f90
    \rm -fv *.o *.mod Makefile .cppdefs input.nml*default 
    \rm -fv obs_diag filter perfect_model_obs create_fixed_network_seq create_obs_sequence 
    \rm -fv assim_region integrate_model

    \rm -fv input.nml obs_seq.in obs_seq.out obs_seq.final \
            perfect_ics fitler_ics \
            True_State.nc Prior_Diag.nc Posterior_Diag.nc
    cvs update

    csh workshop_setup.csh

    if ( -e obs_seq.final ) then
       ./obs_diag
    endif

    #------------------------------------------------------------------
    # Save the output to a directory for comparison
    #------------------------------------------------------------------

    mkdir -pv ${ARCHIVEDIR}/${MODEL}/work

    mv -v dart_log.out        ${ARCHIVEDIR}/${MODEL}/work
    mv -v True_State.nc       ${ARCHIVEDIR}/${MODEL}/work
    mv -v perfect_restart     ${ARCHIVEDIR}/${MODEL}/work
    mv -v obs_seq.out         ${ARCHIVEDIR}/${MODEL}/work
    mv -v Prior_Diag.nc       ${ARCHIVEDIR}/${MODEL}/work
    mv -v Posterior_Diag.nc   ${ARCHIVEDIR}/${MODEL}/work
    mv -v obs_seq.final       ${ARCHIVEDIR}/${MODEL}/work
    mv -v filter_restart      ${ARCHIVEDIR}/${MODEL}/work
    mv -v assim_tools_restart ${ARCHIVEDIR}/${MODEL}/work

    mv -v *ges_times*dat      ${ARCHIVEDIR}/${MODEL}/work
    mv -v *anl_times*dat      ${ARCHIVEDIR}/${MODEL}/work
    mv -v ObsDiagAtts.m       ${ARCHIVEDIR}/${MODEL}/work

    #------------------------------------------------------------------
    # clean up and move on
    #------------------------------------------------------------------
    \rm -fv *.o *.mod Makefile .cppdefs input.nml*default 
    \rm -fv obs_diag filter perfect_model_obs create_fixed_network_seq create_obs_sequence 
    \rm -fv assim_region integrate_model preprocess

   @ makenum  = $makenum  + 1
   @ modelnum = $modelnum + 1
end

#----------------------------------------------------------------------
# Run the full suite of lorenz_96 tests and compile much 
#----------------------------------------------------------------------

cd ${DARTHOME}

./test_dart.csh

echo ""
echo "Testing complete  at "`date`
echo "-------------------------------------------"

