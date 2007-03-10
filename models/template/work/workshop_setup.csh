#!/bin/csh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines automatically updated by version control software, do not edit>
# $Id$
# $Source: /home/thoar/CVS.REPOS/DART/models/template/work/workshop_setup.csh,v $

\rm -f preprocess create_obs_sequence create_fixed_network_seq
\rm -f perfect_model_obs filter obs_diag integrate_model
\rm -f merge_obs_seq
\rm -f *.o *.mod

csh mkmf_preprocess
make         || exit 1
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90
./preprocess || exit 2

#----------------------------------------------------------------------

csh mkmf_create_obs_sequence
make         || exit 3
csh mkmf_create_fixed_network_seq
make         || exit 4
csh mkmf_perfect_model_obs
make         || exit 5
csh mkmf_filter
make         || exit 6
csh mkmf_obs_diag
make         || exit 7
csh mkmf_integrate_model
make         || exit 8
csh mkmf_merge_obs_seq
make         || exit 9

#./perfect_model_obs || exit 20
#./filter            || exit 21
#\rm -f go_end_filter
