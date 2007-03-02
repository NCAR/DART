#!/bin/csh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines updated by version control software, do not edit>
# $Id$
# $Source: /home/thoar/CVS.REPOS/DART/models/PBL_1d/work/workshop_setup.csh,v $
# $Name:  $

#----------------------------------------------------------------------
# Script to manage the compilation of all components for this model;
# executes a known "perfect model" experiment using an existing
# observation sequence file (obs_seq.in) and initial conditions appropriate 
# for both 'perfect_model_obs' (perfect_ics) and 'filter' (filter_ics).
# There are enough initial conditions for 80 ensemble members in filter.
# Use ens_size = 81 and it WILL bomb. Guaranteed.
# The 'input.nml' file controls all facets of this execution.
#
# 'create_obs_sequence' and 'create_fixed_network_sequence' were used to
# create the observation sequence file 'obs_seq.in' - this defines 
# what/where/when we want observations. This script does not run these 
# programs - intentionally. 
#
# 'perfect_model_obs' results in a True_State.nc file that contains 
# the true state, and obs_seq.out - a file that contains the "observations"
# that will be assimilated by 'filter'.
#
# 'filter' results in three files (at least): Prior_Diag.nc - the state 
# of all ensemble members prior to the assimilation (i.e. the forecast), 
# Posterior_Diag.nc - the state of all ensemble members after the 
# assimilation (i.e. the analysis), and obs_seq.final - the ensemble 
# members' estimate of what the observations should have been.
#
# Once 'perfect_model_obs' has advanced the model and harvested the 
# observations for the assimilation experiment, 'filter' may be run 
# over and over by simply changing the namelist parameters in input.nml.
#
# The result of each assimilation can be explored in model-space with
# matlab scripts that directly read the netCDF output, or in observation-space.
# 'obs_diag' is a program that will create observation-space diagnostics
# for any result of 'filter' and results in a couple data files that can
# be explored with yet more matlab scripts.
#
#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the 
# resulting source file is used by all the remaining programs, 
# so this MUST be run first.
#----------------------------------------------------------------------
#
# If you get a ton of compile ERRORS (not warnings) read on ...
#
# Since a lot of the code is inherited from wrf, it comes with a .F 
# extension even though it is free-format. This makes it necessary to
# compile with flags that force interpretation of free-format.
# They also rely on the autopromotion flag ... arghhh ... -r8
# Intel     -free -r8
# gfortran  -ffree-form -fdefault-real-8
# pathscale -freeform -r8
# pgi       -Mfree -Mr8
# absoft    -ffree  (see the mkmf.template for absoft for more on r8)
#----------------------------------------------------------------------

\rm -f preprocess gen_init create_obs_sequence create_fixed_network_seq
\rm -f perfect_model_obs filter obs_diag create_real_network_seq 
\rm -f driver.x merge_obs_seq

csh mkmf_preprocess
make         || exit 1
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90
./preprocess || exit 2

#----------------------------------------------------------------------
echo ""
echo "Building this model generally requires the fortran free-format flag and"
echo "the real*8 override flag be added to the default mkmf.template rules."
echo "If the following compile fails read the comments in the workshop_setup.csh"
echo "script for more help."
echo ""

@ n = 2
foreach TARGET ( mkmf_* )

   switch ( $TARGET )
   case mkmf_preprocess:
      breaksw
   default:
      @ n = $n + 1
      echo
      echo "---------------------------------------------------"
      echo "build number $n is ${TARGET}" 
      csh $TARGET
      make || exit $n
      breaksw
   endsw
end

#----------------------------------------------------------------------
# check for the input data files.  they are large ( > 90Mb ) and
# so are not part of the default distribution.  they need to be 
# downloaded separately and put in $DART/models/PBL_1d/indata.

if ( ! -f ../indata/wrfrt_2006.nc ) then
    echo "Error:"
    echo "This model requires some large data files as input which are"
    echo "not packaged as part of the DART distribution. To run this model"
    echo "contact thoar at ucar dot edu for more information on how"
    echo "to get a copy of the data."
    exit 100
endif

./perfect_model_obs  || exit 20
./filter             || exit 21

#----------------------------------------------------------------------
# The observation-space diagnostics program is not fully developed yet.
# In order to match the bahavior of the other models that use the threed_sphere
# location module, the obs_diag.final  file must exist in a directory.
# We're hardcoding that here. Clearly suboptimal.

if (! -d 06_01) then
   mkdir 06_01
endif

if ( -e 06_01/obs_seq.final ) then
     mv -v 06_01/obs_seq.final 06_01/obs_seq.final.old
endif

\cp -p obs_seq.final 06_01/obs_seq.final

./obs_diag   || exit 99
