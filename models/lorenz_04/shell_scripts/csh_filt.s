#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This script builds the makefile for 'filter', and then builds 'filter'.
#
# This is the last of the csh_*.s scripts to execute.

csh mkmf_filter
make

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

