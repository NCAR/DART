#!/bin/csh

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

 set FC = ifort
 # or pgf90, gfortran, xlf90, g95, etc
 set FCFLAGS = "-r8 -pc 64"

 \rm -f *.o

 ${FC} -o stat_test.x ${FCFLAGS} stat_test.f
 ${FC} -o arg_test.x ${FCFLAGS} arg_test.f
 ${FC} -o ../exe/grabbufr.x ${FCFLAGS} grabbufr.f ../lib/bufrlib.a

