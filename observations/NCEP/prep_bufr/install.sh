#!/bin/sh 
#
# This file is not protected by the DART copyright agreement.
# DART $Id$
#
#  ------------------------------------------------------------------------
#  This script will make executables which extract data
#  from ADP BUFR input files, and place the data into a basic text file.
#  prepbufr.x:  used to extract data from prepbufr files
#
#  If you get an error about 'ar' not being found, make sure the
#  directory containing the 'ar' command is part of your path.
#  'which ar' will print the location if this is done correctly.
# 
#  If you get a link or runtime error about 'bort' being undefined, go into
#  the lib directory and read the README_BUFRLIB file about how to fix it.
#  You can do something like:  cc='cc -DUNDERSCORE' to set the flag for all
#  compiles at once.
#
#  Note that some versions of gfortran give bogus warnings
#  about using goto to go to the end of a do loop when compiling
#  the prepbufr library routines.  This is a mistake in the 
#  gfortran warning logic (the lib code is ok), and can be ignored.
#
#  ------------------------------------------------------------------------

set -eua
 
#  ------------------------------------------------------------------------
#  CPLAT - platform type (linux, sgi, aix, sun, macosx)
#  ------------------------------------------------------------------------
 
# this choice here should really be compiler, i.e. ifort, pgi, gfortran
#CPLAT=macosx
CPLAT=linux
#CPLAT=aix

#  set up the compilers to use
#  -----------------------------------------------------
 
if [ $CPLAT = sgi ]
then
   cc=cc; ff=f77
elif [ $CPLAT = linux ]
then
# possible different compiler choices
#   cc='cc -O'; ff='pgf90 -O'
#   cc='icc -O'; ff='ifort -O'
    cc='gcc -DUNDERSCORE -O'; ff='ifort -O'
elif [ $CPLAT = aix ]
then
   cc='cc -O'; ff='f77 -O'
elif [ $CPLAT = sun ]
then
   cc=cc; ff=f77
elif [ $CPLAT = macosx ]
then
   cc='gcc -DUNDERSCORE'; ff=gfortran
fi


#  Compile and archive the Bufr Library
#  ------------------------------------

echo 'Compiling the Bufr library'

cd lib
$ff -c *.f
$cc -c *.c
ar crv bufrlib.a *.o
rm *.o
 
#  Compile and link the decode programs
#  ---------------------------------------

echo 'Compiling the prepbufr programs'

cd ../src
$ff prepbufr.f     ../lib/bufrlib.a -o ../exe/prepbufr.x 
$ff prepbufr_03Z.f ../lib/bufrlib.a -o ../exe/prepbufr_03Z.x 

 
# If needed, compile the 2 auxiliary conversion programs.

#  Compile the binary format converter program
#  ---------------------------------------
 
echo 'Compiling the grabbufr converter program'
 
cd ../convert_bufr
$ff grabbufr.f ../lib/bufrlib.a -o ../exe/grabbufr.x
 
 
#  Compile the block/unblock converter program
#  ---------------------------------------
 
echo 'Compiling the blk/ublk converter program'
 
cd ../blk_ublk
$ff cwordsh.f ../lib/bufrlib.a -o ../exe/cword.x
 

#  clean up
#  --------

echo 'Finished making executables'
cd ..

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

