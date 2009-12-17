#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This script builds the makefile for 'perfect_model_obs',
# and the builds 'perfect_model_obs'.
#
# This is the third of the csh_*.s scripts to execute.
#
csh mkmf_perfect_model_obs
make

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

