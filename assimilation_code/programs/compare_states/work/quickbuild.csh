#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This script compiles all executables in this directory.

\rm -f *.o *.mod Makefile .cppdefs

set MODEL = "compare_states"

@ n = 0

#----------------------------------------------------------------------
# Build all the single-threaded targets
#----------------------------------------------------------------------

foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   @ n = $n + 1
   echo
   echo "---------------------------------------------------"
   echo "${MODEL} build number ${n} is ${PROG}"
   \rm -f ${PROG}
   csh $TARGET || exit $n
   make        || exit $n

end

# clean up.  comment this out if you want to keep the .o and .mod files around
\rm -f *.o *.mod input.nml.*_default Makefile .cppdefs

echo "Success: All DART programs compiled."

exit 0

