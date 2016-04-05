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
#
# This script does the mkmf on create_obs_sequence,
# and performs the make.
#
# This is the second of the csh_*.s scripts to execute.
#
csh mkmf_create_obs_sequence
make

