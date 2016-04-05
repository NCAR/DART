#!/bin/csh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section, 
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# DESCRIPTION:

 set echo
 set FCFLAGS	= "-r8 -pc 64"
 rm *.o

 pgf90 -c ${FCFLAGS} grabbufr.f
 pgf90 -o ../exe/grabbufr.x ${FCFLAGS} grabbufr.o ../lib/bufrlib.a

