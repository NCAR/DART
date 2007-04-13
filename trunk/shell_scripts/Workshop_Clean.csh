#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

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
      echo "This script cleans up after an execution of 'Workshop_Run.csh'"
      echo "It removes all the executables, .o, and .mod files. in all of "
      echo "the same model directories exercised by Workshop_Run.csh."
      echo " "
      echo "This must be run from the top-level 'DART' directory."
      echo " "
      exit 1
      breaksw
endsw

if ( ! -d models/lorenz_96 ) then
   echo "$SNAME:t must be run from the top-level"
   echo "DART directory -- please try again."
   exit 2
else
   set DARTHOME = `pwd`
endif

echo "The top-level DART directory (DARTHOME) is $DARTHOME"

#----------------------------------------------------------------------
# Remove the result of running workshop_setup.csh for each model
#----------------------------------------------------------------------

foreach MODEL ( 9var lorenz_63 lorenz_84 lorenz_96 lorenz_96_2scale \
    forced_lorenz_96 lorenz_04 bgrid_solo pe2lyr )

    cd ${DARTHOME}/models/${MODEL}/work

    \rm -fv *.o 
    \rm -fv *.mod 
    \rm -fv input.nml*default 
    \rm -fv *raw*_times_*.dat 
    \rm -fv ../../../obs_def/obs_def_mod.f90
    \rm -fv ../../../obs_kind/obs_kind_mod.f90
    \rm -fv obs_diag filter perfect_model_obs create_fixed_network_seq .cppdefs Makefile \
            create_obs_sequence preprocess go_end_filter obs_seq.final \
            filter_restart perfect_restart ObsDiagAtts.m assim_tools_restart \
            assim_region integrate_model dart_log.out

end

