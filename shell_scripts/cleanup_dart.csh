#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

set SNAME = $0
set clobber

switch ( $#argv )
   case 0:
      # supplying no arguments -- echo usage not
      breaksw
   default:
      echo " "
      echo "usage: $SNAME:t"
      echo " "
      echo "This script compiles 'filter' for a wide range of models and then does"
      echo "relatively extensive tests of the L96 programs with a variety of options."
      echo " "
      echo "This must be run from the top-level 'DART' directory."
      echo " "
      echo "This is a pretty verbose process, so if you are logging the output,"
      echo "make sure you have plenty of space:"
      echo " "
      echo "./$SNAME:t |& tee DART_test.log"
      echo " "
      echo "can easily result in a 750 Kb log file"
      exit 1
      breaksw
endsw

if ( ! -d models/lorenz_96 ) then
   echo "models/lorenz_96 does not exist. $SNAME:t must be run from the top-level"
   echo "DART directory -- please try again."
   exit 2
else
   set DARTHOME = `pwd`
endif

echo "The top-level DART directory (DARTHOME) is $DARTHOME"
echo "----------------------------------------------------------"
echo "Cleaning the lorenz_96 (L96) at "`date`
echo ""

cd ${DARTHOME}/models/lorenz_96/work

# Make sure that all .o, .mod and executables are gone
rm -rfv *.o *.mod assim_region create_fixed_network_seq create_obs_sequence filter
rm -rfv integrate_model perfect_model_obs Makefile vi_script obs_seq.out obs_seq.in
rm -rfv filter_restart.???? go_end_filter dart_log.out ens_manager_ens_file.0001
rm -rfv assim_region.csh advance_ens.csh assim_filter.csh advance_model.csh filter_server.csh 
rm -rfv perfect_ics.spun_up filter_ics.spun_up temp_input set_def.out perfect_ics.10hour filter_ics.10hour
rm -rfv perfect_restart.baseline obs_seq.out.baseline True_State.nc.baseline input.nml
rm -rfv perfect_restart.out_of_core obs_seq.out.out_of_core  True_State.nc.out_of_core
rm -rfv perfect_restart.2 obs_seq.out.2 True_State.nc.2 
rm -rfv perfect_restart.3 obs_seq.out.3 True_State.nc.3 
rm -rfv obs_seq.final.baseline filter_restart.baseline assim_tools_restart.baseline 
rm -rfv Prior_Diag.nc.baseline Posterior_Diag.nc.baseline 
rm -rfv obs_seq.final.in_files filter_restart.in_files assim_tools_restart.in_files 
rm -rfv Prior_Diag.nc.in_files Posterior_Diag.nc.in_files 

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

