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

# This script does the mkmf on create_fixed_network_seq,
# and performs the make.
#
# This is the first of the csh_*.s scripts to execute.
#
csh mkmf_create_fixed_network_seq
make

