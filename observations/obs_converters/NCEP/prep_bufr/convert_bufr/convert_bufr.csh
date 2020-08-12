#!/bin/csh
#
# This file is not protected by the DART copyright agreement.
# DART $Id$

 set FC = ifort
 # or pgf90, gfortran, xlf90, g95, etc
 set FCFLAGS = "-r8 -pc 64"

 \rm -f *.o

 ${FC} -o stat_test.x ${FCFLAGS} stat_test.f
 ${FC} -o arg_test.x ${FCFLAGS} arg_test.f
 ${FC} -o ../exe/grabbufr.x ${FCFLAGS} grabbufr.f ../lib/bufrlib.a

exit 0


