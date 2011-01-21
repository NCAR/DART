#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

set SNAME = $0
set clobber

set startdir=`pwd`

foreach dir ( */oned */threed_sphere */*/work */work )
 cd $dir
 echo building in $dir
 echo

 ./quickbuild.csh || echo ERROR unsuccessful build in $dir

 cd $startdir
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

