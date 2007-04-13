#!/bin/sh 
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

#  ------------------------------------------------------------------------
#  This script will make executables which extract data
#  from ADP BUFR input files, and place the data into a basic text file.
#  prepbufr.x:  used to extract data from prepbufr files
#  ** Make sure the "ar" command location has been set in your path
#  environment variable.  Type "which ar" to check if this is done. **
#  ------------------------------------------------------------------------
 
set -eua
 
#  ------------------------------------------------------------------------
#  CPLAT - platform type (linux, sgi, aix, sun, macosx)
#  ------------------------------------------------------------------------
 
CPLAT=linux
SRC=../src
LIB=../lib
EXE=../exe

#  different platforms use different link name protocols
#  -----------------------------------------------------
 
if [ $CPLAT = sgi ]
then
   openrb=openrb_
   openwb=openwb_
   crdbfr=crdbufr_
   cwrbfr=cwrbufr_
   lenmsg=lenm_
   cc=cc; ff=f77
elif [ $CPLAT = linux ]
then
   openrb=openrb
   openwb=openwb
   crdbfr=crdbufr
   cwrbfr=cwrbufr
   lenmsg=lenm
   cc=cc; ff=pgf90
elif [ $CPLAT = aix ]
then
   openrb=openrb
   openwb=openwb
   crdbfr=crdbufr
   cwrbfr=cwrbufr
   lenmsg=lenm
   cc=cc; ff=f77
elif [ $CPLAT = sun ]
then
   openrb=openrb_
   openwb=openwb_
   crdbfr=crdbufr_
   cwrbfr=cwrbufr_
   lenmsg=lenm_
   cc=cc; ff=f77
elif [ $CPLAT = macosx ]
then
   openrb=openrb
   openwb=openwb
   crdbfr=crdbufr
   cwrbfr=cwrbufr
   lenmsg=lenm
   cc=gcc; ff=gfortran
fi

#  Compile and archive the Bufr Library
#  ------------------------------------

$ff -c $LIB/*.f

ar crv $LIB/bufrlib.a *.o

rm *.o
 
#  Compile the decode programs
#  ---------------------------------------
 
$ff -c $SRC/prepbufr.f

$ff -c $SRC/prepbufr_03Z.f

 
#  link and load the executables
#  -----------------------------

$ff -o $EXE/prepbufr.x prepbufr.o $LIB/bufrlib.a
$ff -o $EXE/prepbufr_03Z.x prepbufr_03Z.o $LIB/bufrlib.a


#  clean up
#  --------

rm -f *.o
