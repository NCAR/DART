#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next five lines automatically updated by CVS, do not edit>
# $Source$
# $Revision$
# $Date$
# $Author$
# $Name$

csh mkmf_create_fixed_network_seq
make
csh mkmf_create_obs_sequence
make
csh mkmf_filter
make
csh mkmf_perfect_ics
make
csh mkmf_perfect_model_obs
make
csh mkmf_dart_to_MITgcm
make
csh mkmf_MITgcm_to_dart
make



