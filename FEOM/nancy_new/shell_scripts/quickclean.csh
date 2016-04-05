#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script cleans out the work dir.


foreach TARGET ( mkmf_* )

  csh $TARGET || exit $n
  make clean  || exit $n

end

\rm -f *.o *.mod 
\rm -f input.nml*_default
\rm -f dart_log.*

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

