#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

 set FC = ifort
 # or pgf90, gfortran, xlf90, g95, etc
 set FCFLAGS = "-r8 -pc 64"

 \rm -f *.o

 ${FC} -o stat_test.x ${FCFLAGS} stat_test.f
 ${FC} -o arg_test.x ${FCFLAGS} arg_test.f
 ${FC} -o ../exe/grabbufr.x ${FCFLAGS} grabbufr.f ../lib/bufrlib.a

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

