#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

#----------------------------------------------------------------------
# compile all programs in the current directory that have a mkmf_xxx file.
#----------------------------------------------------------------------

# this item's name:
set ITEM = "Generate Sampling Error Correction Table"


# ---------------
# shouldn't have to modify this script below here.

\rm -f *.o *.mod Makefile .cppdefs

@ n = 0

foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's/mkmf_//'`

   @ n = $n + 1
   echo
   echo "---------------------------------------------------"
   echo "$ITEM build number $n is $PROG"
   \rm -f $PROG
   csh $TARGET || exit $n
   make        || exit $n

end

echo "Success: All programs compiled."

\rm -f *.o *.mod  input.nml.*_default Makefile .cppdefs

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

