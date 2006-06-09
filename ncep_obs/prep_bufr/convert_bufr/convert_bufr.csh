#!/bin/csh -f

 set echo
 set FCFLAGS	= "-r8 -pc 64"
 rm *.o

 pgf90 -c ${FCFLAGS} grabbufr.f
 pgf90 -o ../exe/grabbufr.x ${FCFLAGS} grabbufr.o ../lib/bufrlib.a

